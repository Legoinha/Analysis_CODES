"""Mass decorrelation of classifier inputs by quantile morphing.

A feature f whose background distribution changes with the candidate mass m is replaced by
f_dc, whose background distribution is the same at every mass:

1. Reference: DATA after the pre-cut in the decorrelation range of the case, outside its excluded
   windows (the signal region and other peaks; ML_selection/settings.py MASS). This is background only.
   X(3872): 3.60-4.00 GeV outside 3.85-3.89 (the gap between the training sidebands) and +-20 MeV around
   the psi(2S).
2. In slices of decorrelation_slice, the background percentiles of f at LEVELS (0.5%, 1.5%, ..., 99.5%).
   Each percentile is fitted versus mass with a polynomial of degree DEGREE in (m - peak), which also
   gives it under the excluded peaks by interpolation.
3. For a candidate with mass m, f is mapped piecewise-linearly from the percentiles at m to the
   percentiles at the peak: f_dc = Q_ref(F(f | m)). Below the lowest and above the highest percentile
   it is shifted by the difference of those percentiles, so the tails keep their order. For
   background, f_dc has the distribution of background at the peak at every mass; at the peak the map
   is the identity, so signal features keep their values.

The fitted coefficients are stored per case in a ROOT file (one tree map_<feature> with the level and
the coefficients, and the reference mass in metadata/), read here and by ML_tmva/TMVA_decorrelation.h.
"""
from datetime import datetime, timezone

import numpy as np
import uproot


MIN_SLICE_ENTRIES = 500
LEVELS = np.linspace(0.005, 0.995, 100)
DEGREE = 3


def fit_region(bmass, c):
    """The DATA of case c (ml_samples.case) the maps are fitted on."""
    bmass = np.asarray(bmass, np.float64)
    keep = (bmass > c["decorrelation_range"][0]) & (bmass < c["decorrelation_range"][1])
    for low, high in c["decorrelation_excluded"]:
        keep &= ~((bmass > low) & (bmass < high))
    return keep


def build(bmass, features, c):
    """Percentile curves of every feature, and the (low, high) of the slices fully inside the fit
    region with enough entries. features: {name: values}, same length as bmass."""
    bmass = np.asarray(bmass, np.float64)
    width = c["decorrelation_slice"]
    edges = np.arange(c["decorrelation_range"][0], c["decorrelation_range"][1] + 1e-9, width)
    parts = [(low, high) for low, high in zip(edges[:-1], edges[1:])
             if not any(low < ex_high and high > ex_low for ex_low, ex_high in c["decorrelation_excluded"])
             and np.sum((bmass >= low) & (bmass < high)) >= MIN_SLICE_ENTRIES]
    centers = np.array([0.5 * (low + high) for low, high in parts]) - c["peak"]
    maps = {}
    for name, values in features.items():
        values = np.asarray(values, np.float64)
        quantiles = np.array([np.quantile(values[(bmass >= low) & (bmass < high)], LEVELS) for low, high in parts])
        maps[name] = np.array([np.polyfit(centers, quantiles[:, k], DEGREE) for k in range(len(LEVELS))])
    return maps, parts


class Maps:
    def __init__(self, coefficients, reference_mass, metadata):
        self.coefficients = coefficients   # {feature: array (levels, DEGREE + 1), highest power first}
        self.reference_mass = reference_mass
        self.metadata = metadata
        self.reference = {name: self.percentiles(name, np.array([reference_mass]))[0] for name in coefficients}

    def percentiles(self, name, bmass):
        """Background percentiles of `name` at each mass, shape (n, levels), non-decreasing."""
        x = np.asarray(bmass, np.float64) - self.reference_mass
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


def save(path, maps, c, metadata):
    with uproot.recreate(path) as out:
        for name, coefficients in maps.items():
            out[f"map_{name}"] = {"level": LEVELS.astype(np.float64),
                                  **{f"c{DEGREE - p}": coefficients[:, p].astype(np.float64) for p in range(DEGREE + 1)}}
        meta = {
            "reference_mass": c["peak"],
            "degree": DEGREE,
            "mass_range": list(c["decorrelation_range"]),
            "excluded": [list(e) for e in c["decorrelation_excluded"]],
            "slice_width": c["decorrelation_slice"],
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
    return Maps(coefficients, float(metadata["reference_mass"]), metadata)


def closure(maps, name, transformed, bmass, regions):
    """Largest |fraction below reference percentile - level| of the transformed background (the _dc
    values) in each (low, high) mass region: 0 means the transform makes the background mass-independent."""
    out = {}
    ref = maps.reference[name]
    for label, (low, high) in regions.items():
        inside = (bmass >= low) & (bmass < high)
        if inside.sum() < 100:
            continue
        below = np.array([np.mean(transformed[inside] < r) for r in ref])
        out[label] = float(np.max(np.abs(below - LEVELS)))
    return out
