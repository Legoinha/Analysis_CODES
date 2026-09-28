"""Mass decorrelation of classifier inputs by quantile morphing.

A feature f whose background distribution changes with the candidate mass m is replaced by
f_dc, whose background distribution is the same at every mass:

1. Reference: DATA after the pre-cut in 3.60-4.00 GeV, outside the X(3872) region (3.85-3.89,
   the gap between the training sidebands) and outside +-20 MeV around the Psi(2S).
   This is background only.
2. In 10 MeV mass slices, the background percentiles of f at LEVELS (0.5%, 1.5%, ..., 99.5%).
   Each percentile is fitted versus mass with a polynomial of degree DEGREE in (m - m_X), which
   also gives it under the two peaks by interpolation.
3. For a candidate with mass m, f is mapped piecewise-linearly from the percentiles at m to the
   percentiles at the X(3872) mass: f_dc = Q_ref(F(f | m)). Below the lowest and above the
   highest percentile it is shifted by the difference of those percentiles, so the tails keep
   their order. For background, f_dc has the distribution of background at m_X at every mass;
   at m = m_X the map is the identity, so signal features keep their values.

The fitted coefficients are stored per sample in a ROOT file (one tree map_<feature> with the
level and the coefficients), read here and by ML_tmva/TMVA_decorrelation.h.
"""
from datetime import datetime, timezone

import numpy as np
import uproot


X3872_MASS = 3.87164  # as in plotER/aux/masses.h
PSI2S_MASS = 3.68610
MASS_RANGE = (3.60, 4.00)
EXCLUDED = ((PSI2S_MASS - 0.020, PSI2S_MASS + 0.020), (3.85, 3.89))
SLICE_WIDTH = 0.010
MIN_SLICE_ENTRIES = 500
LEVELS = np.linspace(0.005, 0.995, 100)
DEGREE = 3


def fit_region(bmass):
    bmass = np.asarray(bmass, np.float64)
    keep = (bmass > MASS_RANGE[0]) & (bmass < MASS_RANGE[1])
    for low, high in EXCLUDED:
        keep &= ~((bmass > low) & (bmass < high))
    return keep


def slices(bmass):
    """(low, high) of the 10 MeV slices fully inside the fit region with enough entries."""
    edges = np.arange(MASS_RANGE[0], MASS_RANGE[1] + 1e-9, SLICE_WIDTH)
    result = []
    for low, high in zip(edges[:-1], edges[1:]):
        if any(low < ex_high and high > ex_low for ex_low, ex_high in EXCLUDED):
            continue
        if np.sum((bmass >= low) & (bmass < high)) >= MIN_SLICE_ENTRIES:
            result.append((low, high))
    return result


def build(bmass, features):
    """Percentile curves of every feature. features: {name: values}, same length as bmass."""
    bmass = np.asarray(bmass, np.float64)
    parts = slices(bmass)
    centers = np.array([0.5 * (low + high) for low, high in parts]) - X3872_MASS
    maps = {}
    for name, values in features.items():
        values = np.asarray(values, np.float64)
        quantiles = np.array([np.quantile(values[(bmass >= low) & (bmass < high)], LEVELS) for low, high in parts])
        maps[name] = np.array([np.polyfit(centers, quantiles[:, k], DEGREE) for k in range(len(LEVELS))])
    return maps, parts


class Maps:
    def __init__(self, coefficients, metadata=None):
        self.coefficients = coefficients   # {feature: array (levels, DEGREE + 1), highest power first}
        self.metadata = metadata or {}
        self.reference = {name: self.percentiles(name, np.array([X3872_MASS]))[0] for name in coefficients}

    def percentiles(self, name, bmass):
        """Background percentiles of `name` at each mass, shape (n, levels), non-decreasing."""
        x = np.asarray(bmass, np.float64) - X3872_MASS
        c = self.coefficients[name]
        q = np.zeros((len(x), len(c)))
        for power in range(c.shape[1]):
            q = q * x[:, None] + c[None, :, power]
        return np.maximum.accumulate(q, axis=1)

    def transform(self, name, values, bmass):
        values = np.asarray(values, np.float64)
        q = self.percentiles(name, bmass)
        ref = self.reference[name]
        k = np.sum(values[:, None] >= q, axis=1)          # 0 .. levels
        inner = np.clip(k, 1, q.shape[1] - 1)
        rows = np.arange(len(values))
        q_low, q_high = q[rows, inner - 1], q[rows, inner]
        span = np.where(q_high > q_low, q_high - q_low, 1.0)
        t = np.where(q_high > q_low, (values - q_low) / span, 0.0)
        out = ref[inner - 1] + t * (ref[inner] - ref[inner - 1])
        out = np.where(k == 0, values - q[:, 0] + ref[0], out)
        out = np.where(k == q.shape[1], values - q[:, -1] + ref[-1], out)
        return out.astype(np.float32)

    def features(self):
        return list(self.coefficients)


def save(path, maps, metadata):
    with uproot.recreate(path) as out:
        for name, c in maps.items():
            out[f"map_{name}"] = {"level": LEVELS.astype(np.float64),
                                  **{f"c{DEGREE - p}": c[:, p].astype(np.float64) for p in range(DEGREE + 1)}}
        meta = {
            "reference_mass": X3872_MASS,
            "degree": DEGREE,
            "mass_range": list(MASS_RANGE),
            "excluded": [list(e) for e in EXCLUDED],
            "slice_width": SLICE_WIDTH,
            "levels": len(LEVELS),
            "created_utc": datetime.now(timezone.utc).isoformat(),
            **metadata,
        }
        for key, value in meta.items():
            out[f"metadata/{key}"] = str(value)


def load(path):
    coefficients, metadata = {}, {}
    with uproot.open(path) as f:
        for key in f.keys():
            name = key.split(";")[0]
            if name.startswith("map_"):
                a = f[name].arrays(library="np")
                coefficients[name[4:]] = np.column_stack([a[f"c{p}"] for p in range(DEGREE, -1, -1)])
            elif name.startswith("metadata/"):
                metadata[name.split("/")[-1]] = str(f[name])
    return Maps(coefficients, metadata)


def closure(maps, name, values, bmass, regions):
    """Largest |fraction below reference percentile - level| of the transformed background in
    each (low, high) mass region: 0 means the transform makes the background mass-independent."""
    out = {}
    transformed = maps.transform(name, values, bmass)
    ref = maps.reference[name]
    for label, (low, high) in regions.items():
        inside = (bmass >= low) & (bmass < high)
        if inside.sum() < 100:
            continue
        below = np.array([np.mean(transformed[inside] < r) for r in ref])
        out[label] = float(np.max(np.abs(below - LEVELS)))
    return out
