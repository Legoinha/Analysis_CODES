"""Model selection shared by ML_pytorch and ML_xgboost, and the study metrics of all three methods.

- the train / validation / test split,
- the four figures of merit (FOM) a scan can maximize: "auc", "sigeff", "pauc", "punzi",
- the overtraining pass/fail test and the working-point efficiencies,
- bootstrap uncertainties of the FOMs and of the FOM gain between two models,
- the mass-sculpting check on DATA never used for training,
- the scores file every training writes (NN, XGB and TMVA alike) and the metrics read from it.

The validation split chooses the model (early stopping, hyperparameters); the test split is only used
for the final report. ML_tmva/TMVA_common.h implements the same definitions in C++. The settings are
in ML_selection/settings.py: the shared ones below, and those of the mass regions of each tree in the
case dict c (ml_samples.case): peak, punzi_half_width, punzi_sideband_start, punzi_sideband_width,
punzi_sides (-1 below the peak, +1 above) and the sculpting_* regions.
"""
import json
import sys
from pathlib import Path

import numpy as np
import uproot
from scipy.special import kolmogorov
from sklearn.metrics import auc as roc_area, roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split

# settings.py is in ML_selection/; batch jobs get it next to this file.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import settings


FOMS = ("auc", "sigeff", "pauc", "punzi")
VALIDATION_FRACTION = 0.2
TEST_FRACTION = 0.2


def fom_settings(c):
    """Every setting a FOM, the working point or the overtraining test of case c depends on. The sides of the
    Punzi sidebands are recorded only when there is one (B+): with both, the key is left out, so the settings
    of the studies made before the sides existed (X(3872)) stay valid."""
    return {
        "bkg_eff_target": settings.BKG_EFF_TARGET,
        "pauc_max_bkg_eff": settings.PAUC_MAX_BKG_EFF,
        "punzi_a": settings.PUNZI_A,
        "punzi_b": settings.PUNZI_B,
        "signal_half_width": c["punzi_half_width"],
        "sideband_start": c["punzi_sideband_start"],
        "sideband_width": c["punzi_sideband_width"],
        "punzi_min_background": settings.PUNZI_MIN_BACKGROUND,
        "threshold_step": c["threshold_step"],
        "working_point": "Punzi maximum on all signal MC and all DATA, B counted in the Punzi sidebands (optimalCUT_X_punzi.C)",
        "ks_pvalue_min": settings.KS_PVALUE_MIN,
        "max_auc_gap": settings.MAX_AUC_GAP,
        "auc_gap_sigma": settings.AUC_GAP_SIGMA,
        "validation_fraction": VALIDATION_FRACTION,
        "test_fraction": TEST_FRACTION,
        **({"sideband_sides": list(c["punzi_sides"])} if len(c["punzi_sides"]) == 1 else {}),
    }


def split_train_validation_test(arrays, stratify, random_state):
    """Stratified random split of a dict of equally long arrays into training, validation and
    test dicts. The same inputs and random_state always give the same split, whatever the
    feature columns are, so different feature sets share their validation events."""
    index = np.arange(len(arrays[stratify]))
    rest, test = train_test_split(index, test_size=TEST_FRACTION, random_state=random_state,
                                  stratify=arrays[stratify])
    train, validation = train_test_split(rest, test_size=VALIDATION_FRACTION / (1.0 - TEST_FRACTION),
                                         random_state=random_state, stratify=arrays[stratify][rest])
    return tuple({name: values[part] for name, values in arrays.items()} for part in (train, validation, test))


def window_background(bmass, c):
    """Background expected in the Punzi signal window from DATA candidates after the pre-ML cut:
    the count in the Punzi sidebands next to the window (on the punzi_sides of the peak), scaled to its
    width. The same estimate as optimalCUT_X_punzi.C without a classifier cut."""
    offset = np.asarray(bmass, dtype=np.float64) - c["peak"]
    start, width, sides = c["punzi_sideband_start"], c["punzi_sideband_width"], c["punzi_sides"]
    in_sidebands = np.any([(side * offset > start) & (side * offset < start + width) for side in sides], axis=0)
    return float(np.sum(in_sidebands)) * c["punzi_half_width"] / width * 2 / len(sides)


def data_weights(bmass, c):
    """Weight of every DATA candidate of a skim: 1 in the training sidebands, kept in full, and the
    data_thinning of the case elsewhere, where make_skims.py keeps 1 in data_thinning (settings.DATA_THINNING)."""
    bmass = np.asarray(bmass, np.float64)
    in_sidebands = np.any([(bmass > low) & (bmass < high) for low, high in c["training_sidebands"]], axis=0)
    return np.where(in_sidebands, 1, c["data_thinning"])


def partial_auc(fpr, tpr):
    """ROC area for fpr from 0 to PAUC_MAX_BKG_EFF, divided by PAUC_MAX_BKG_EFF. The curve
    is closed at the boundary with the linearly interpolated tpr."""
    limit = settings.PAUC_MAX_BKG_EFF
    k = int(np.searchsorted(fpr, limit, side="right"))
    x, y = fpr[:k], tpr[:k]
    if k < len(fpr):
        edge = tpr[k - 1] + (tpr[k] - tpr[k - 1]) * (limit - fpr[k - 1]) / (fpr[k] - fpr[k - 1])
        x, y = np.append(x, limit), np.append(y, edge)
    return float(roc_area(x, y)) / limit


def punzi_curve(y, score, w, bmass, window_bkg, c, thresholds=None):
    """The Punzi FOM of optimalCUT_X_punzi.C on the thresholds k * c["threshold_step"] covering the scores
    (or the given thresholds; score >= threshold): the signal efficiency, and B = the background
    candidates of the Punzi sidebands of case c that pass, scaled to the window. window_bkg is B before
    the cut from all DATA (window_background), so B = window_bkg x the passing fraction of the sample's
    sideband background; with all DATA as the sample this is exactly the count of optimalCUT_X_punzi.C.
    Returns that fraction, the signal efficiency, the thresholds and the FOM, NaN where fewer than
    PUNZI_MIN_BACKGROUND sideband candidates pass (B would rest on a handful of candidates). A coarse
    grid keeps the maximum from landing on a fluctuation of the curve."""
    y, score, w = (np.asarray(a, np.float64) for a in (y, score, w))
    step = c["threshold_step"]
    if thresholds is None:
        thresholds = np.arange(np.floor(score.min() / step), np.ceil(score.max() / step) + 1) * step
    offset = np.asarray(bmass, np.float64) - c["peak"]
    start, width = c["punzi_sideband_start"], c["punzi_sideband_width"]
    signal = y == 1
    sideband = (y == 0) & np.any([(side * offset > start) & (side * offset < start + width) for side in c["punzi_sides"]], axis=0)
    # Weighted fraction of the signal (sideband) candidates with score >= each threshold.
    fractions = []
    for mask, weight in ((signal, w[signal]), (sideband, np.ones(sideband.sum()))):
        order = np.argsort(score[mask])
        cumulative = np.concatenate([[0.0], np.cumsum(weight[order])])
        fractions.append((cumulative[-1] - cumulative[np.searchsorted(score[mask][order], thresholds, side="left")]) / cumulative[-1])
    efficiency, passing = fractions
    background = window_bkg * passing
    sqrt_b = np.sqrt(background)
    a, b = settings.PUNZI_A, settings.PUNZI_B
    s_min = b**2 / 2.0 + a * sqrt_b + b / 2.0 * np.sqrt(b**2 + 4.0 * a * sqrt_b + 4.0 * background)
    fom = np.where(passing * sideband.sum() >= settings.PUNZI_MIN_BACKGROUND, efficiency / s_min, np.nan)
    return passing, efficiency, thresholds, fom


def all_data(scores):
    """Mass and score of every DATA candidate of a scores file: those never used for training
    (sculpt) and the training background. The Punzi sidebands lie inside the training sidebands, so
    their background is counted on both, as optimalCUT_X_punzi.C does on the full DATA."""
    background = scores["train"]["y"] == 0
    return (np.concatenate([scores["sculpt"]["bmass"], scores["train"]["bmass"][background]]).astype(np.float64),
            np.concatenate([scores["sculpt"]["score"], scores["train"]["score"][background]]).astype(np.float64))


def full_sample(scores):
    """y, score, w, bmass of all signal MC (training, validation and test) and all DATA of a scores
    file: the samples optimalCUT_X_punzi.C scans (DATA with unit weights)."""
    parts = [scores[name] for name in ("train", "validation", "test")]
    bmass, score = all_data(scores)
    signal = [p["y"] == 1 for p in parts]
    return (np.concatenate([np.ones(sum(int(s.sum()) for s in signal)), np.zeros(len(score))]),
            np.concatenate([p["score"][s] for p, s in zip(parts, signal)] + [score]).astype(np.float64),
            np.concatenate([p["w"][s] for p, s in zip(parts, signal)] + [np.ones(len(score))]).astype(np.float64),
            np.concatenate([p["bmass"][s] for p, s in zip(parts, signal)] + [bmass]).astype(np.float64))


def figures_of_merit(y, score, w, bmass, window_bkg, c):
    """AUC, signal efficiency at BKG_EFF_TARGET and partial AUC up to PAUC_MAX_BKG_EFF from the
    weighted ROC curve, and the best Punzi FOM with its threshold (punzi_curve)."""
    fpr, tpr, _ = roc_curve(y, score, sample_weight=w, drop_intermediate=False)
    _, _, grid, punzi = punzi_curve(y, score, w, bmass, window_bkg, c)
    best = int(np.nanargmax(punzi))
    return {
        "auc": float(roc_area(fpr, tpr)),
        "sigeff": float(np.max(tpr[fpr <= settings.BKG_EFF_TARGET])),
        "pauc": partial_auc(fpr, tpr),
        "punzi": float(punzi[best]),
        "punzi_threshold": float(grid[best]),
    }


def ks_statistic(sample_a, sample_b, weights_a, weights_b):
    a = np.asarray(sample_a, dtype=np.float64)
    b = np.asarray(sample_b, dtype=np.float64)
    wa = np.asarray(weights_a, dtype=np.float64)
    wb = np.asarray(weights_b, dtype=np.float64)
    order_a = np.argsort(a)
    order_b = np.argsort(b)
    a = a[order_a]
    b = b[order_b]
    cumulative_a = np.concatenate([[0.0], np.cumsum(wa[order_a]) / np.sum(wa)])
    cumulative_b = np.concatenate([[0.0], np.cumsum(wb[order_b]) / np.sum(wb)])
    grid = np.sort(np.unique(np.concatenate([a, b])))
    cdf_a = cumulative_a[np.searchsorted(a, grid, side="right")]
    cdf_b = cumulative_b[np.searchsorted(b, grid, side="right")]
    return float(np.max(np.abs(cdf_a - cdf_b)))


def efficiency(score, w, threshold):
    """Weighted fraction of entries with score >= threshold and its binomial uncertainty."""
    w = np.asarray(w, dtype=np.float64)
    passed = np.asarray(score) >= threshold
    total = np.sum(w)
    eff = np.sum(w[passed]) / total
    variance = ((1.0 - eff) ** 2 * np.sum(w[passed] ** 2) + eff**2 * np.sum(w[~passed] ** 2)) / total**2
    return float(eff), float(np.sqrt(variance))


def evaluate(train, other, window_bkg, c, threshold=None):
    """FOMs on `other` and the overtraining test of `train` against `other`.

    train and other are (y, score, w, bmass) tuples. The working point is the Punzi-optimal
    threshold of `other`, or `threshold` when given; signal and background efficiencies
    at the working point are compared between the two sets. The model passes when the AUC gap is at
    most MAX_AUC_GAP or AUC_GAP_SIGMA times its statistical uncertainty, whichever is larger, and both
    weighted KS p-values (effective sample sizes (sum w)^2 / sum w^2) are at least KS_PVALUE_MIN.
    """
    y_train, score_train, w_train, _ = train
    y_other, score_other, w_other, bmass_other = other
    result = figures_of_merit(y_other, score_other, w_other, bmass_other, window_bkg, c)
    result["train_auc"] = float(roc_auc_score(y_train, score_train, sample_weight=w_train))
    result["auc_gap"] = result["train_auc"] - result["auc"]
    # Hanley-McNeil standard error of each weighted AUC with the effective entries of each class (within
    # 10% of the bootstrap error of the validation AUC); the two AUCs are on disjoint candidates.
    variance = 0.0
    for auc, y, w in ((result["train_auc"], y_train, w_train), (result["auc"], y_other, w_other)):
        w = np.asarray(w, dtype=np.float64)
        n_s, n_b = (np.sum(w[y == label]) ** 2 / np.sum(w[y == label] ** 2) for label in (1, 0))
        q1, q2 = auc / (2.0 - auc), 2.0 * auc**2 / (1.0 + auc)
        variance += (auc * (1.0 - auc) + (n_s - 1.0) * (q1 - auc**2) + (n_b - 1.0) * (q2 - auc**2)) / (n_s * n_b)
    result["auc_gap_error"] = float(np.sqrt(variance))
    result["auc_gap_limit"] = max(settings.MAX_AUC_GAP, settings.AUC_GAP_SIGMA * result["auc_gap_error"])
    result["working_point"] = result["punzi_threshold"] if threshold is None else float(threshold)
    for name, label in (("signal", 1), ("background", 0)):
        in_train = y_train == label
        in_other = y_other == label
        distance = ks_statistic(score_train[in_train], score_other[in_other], w_train[in_train], w_other[in_other])
        n_a, n_b = (np.sum(w) ** 2 / np.sum(w**2) for w in (np.asarray(w_train[in_train], np.float64),
                                                             np.asarray(w_other[in_other], np.float64)))
        eff_train, err_train = efficiency(score_train[in_train], w_train[in_train], result["working_point"])
        eff_other, err_other = efficiency(score_other[in_other], w_other[in_other], result["working_point"])
        result[f"ks_{name}"] = distance
        result[f"ks_{name}_pvalue"] = float(kolmogorov(np.sqrt(n_a * n_b / (n_a + n_b)) * distance))
        result[f"{name}_eff_train"] = eff_train
        result[f"{name}_eff_train_error"] = err_train
        result[f"{name}_eff"] = eff_other
        result[f"{name}_eff_error"] = err_other
        result[f"{name}_eff_pull"] = (eff_train - eff_other) / np.hypot(err_train, err_other)
    result["passed"] = bool(
        result["auc_gap"] <= result["auc_gap_limit"]
        and result["ks_signal_pvalue"] >= settings.KS_PVALUE_MIN
        and result["ks_background_pvalue"] >= settings.KS_PVALUE_MIN
    )
    return result


def print_evaluation(label, result):
    print(f"{label}: AUC {result['auc']:.5f} (train {result['train_auc']:.5f}, gap {result['auc_gap']:+.5f} "
          f"+- {result['auc_gap_error']:.5f}, limit {result['auc_gap_limit']:.5f}), "
          f"sigeff {result['sigeff']:.4f}, pauc {result['pauc']:.4f}, Punzi FOM {result['punzi']:.6f} at threshold {result['punzi_threshold']:.4f}")
    print(f"{label}: KS signal {result['ks_signal']:.4f} (p {result['ks_signal_pvalue']:.3f}), "
          f"KS background {result['ks_background']:.4f} (p {result['ks_background_pvalue']:.3f}), "
          f"overtraining test {'passed' if result['passed'] else 'FAILED'}")
    for name in ("signal", "background"):
        print(f"{label}: {name} efficiency at working point {result['working_point']:.4f}: "
              f"train {result[f'{name}_eff_train']:.4f} +- {result[f'{name}_eff_train_error']:.4f}, "
              f"{label.lower()} {result[f'{name}_eff']:.4f} +- {result[f'{name}_eff_error']:.4f} "
              f"(pull {result[f'{name}_eff_pull']:+.2f})")


# =============================================================================
# Uncertainties
# =============================================================================

N_BOOTSTRAP = 200


def bootstrap_foms(y, score, w, bmass, window_bkg, c, n=N_BOOTSTRAP, seed=1):
    """Standard deviation of each FOM over bootstrap resamples of (y, score, w, bmass)."""
    rng = np.random.default_rng(seed)
    values = {fom: [] for fom in FOMS}
    for _ in range(n):
        pick = rng.integers(0, len(y), len(y))
        result = figures_of_merit(y[pick], score[pick], w[pick], bmass[pick], window_bkg, c)
        for fom in FOMS:
            values[fom].append(result[fom])
    return {fom: float(np.std(v)) for fom, v in values.items()}


def paired_gain(new, old, window_bkg, c, n=N_BOOTSTRAP, seed=1):
    """FOM gain of model `new` over model `old` on the same validation candidates, with its
    bootstrap standard deviation. new and old are dicts of y, idx, w, bmass, score; candidates are
    matched by (y, idx). The same resamples are used for both models, so the common
    fluctuations cancel and the gain is much more precise than either FOM."""
    key_new = np.asarray(new["y"], np.int64) * (1 << 40) + np.asarray(new["idx"], np.int64)
    key_old = np.asarray(old["y"], np.int64) * (1 << 40) + np.asarray(old["idx"], np.int64)
    common, pick_new, pick_old = np.intersect1d(key_new, key_old, return_indices=True)
    y, w, m = new["y"][pick_new], new["w"][pick_new], new["bmass"][pick_new]
    s_new, s_old = new["score"][pick_new], old["score"][pick_old]
    gain = {fom: figures_of_merit(y, s_new, w, m, window_bkg, c)[fom] - figures_of_merit(y, s_old, w, m, window_bkg, c)[fom]
            for fom in FOMS}
    rng = np.random.default_rng(seed)
    samples = {fom: [] for fom in FOMS}
    for _ in range(n):
        pick = rng.integers(0, len(y), len(y))
        a = figures_of_merit(y[pick], s_new[pick], w[pick], m[pick], window_bkg, c)
        b = figures_of_merit(y[pick], s_old[pick], w[pick], m[pick], window_bkg, c)
        for fom in FOMS:
            samples[fom].append(a[fom] - b[fom])
    return {fom: {"gain": float(gain[fom]), "sigma": float(np.std(samples[fom]))} for fom in FOMS}, int(len(common))


# =============================================================================
# Sculpting
# =============================================================================

def sculpting(bmass, score, thresholds, c):
    """Background efficiency versus mass on DATA never used for training (in the training
    sidebands, their held-out part).

    For each threshold: a straight line fitted to the efficiency in sculpting_bin bins over
    sculpting_range outside sculpting_excluded (slope relative to the mean efficiency per 100 MeV,
    and chi2/ndf), and the ratio R of the efficiency measured in the gaps (sculpting_gap off the peak
    to the excluded edges, the primary local test) and in the Punzi sidebands to the line's
    prediction there. R = 1 and a small slope mean the classifier keeps the same fraction of
    background next to the peak as far from it.
    """
    bmass = np.asarray(bmass, np.float64)
    score = np.asarray(score, np.float64)
    peak = c["peak"]
    low_edge, high_edge = c["sculpting_excluded"]
    distance = np.abs(bmass - peak)
    start, width = c["punzi_sideband_start"], c["punzi_sideband_width"]
    regions = {
        "gap": (distance > c["sculpting_gap"]) & (bmass > low_edge) & (bmass < high_edge),
        "punzi_sidebands": np.any([(side * (bmass - peak) > start) & (side * (bmass - peak) < start + width)
                                   for side in c["punzi_sides"]], axis=0),
    }
    width = c["sculpting_bin"]
    edges = np.concatenate([np.arange(c["sculpting_range"][0], low_edge + 1e-9, width),
                            np.arange(high_edge, c["sculpting_range"][1] + 1e-9, width)])
    bins = [(lo, hi) for lo, hi in zip(edges[:-1], edges[1:]) if not (lo >= low_edge and hi <= high_edge)]
    result = {}
    for name, threshold in thresholds.items():
        passed = score >= threshold
        centers, effs, errs = [], [], []
        for lo, hi in bins:
            inside = (bmass >= lo) & (bmass < hi)
            if inside.sum() < 20:
                continue
            eff = passed[inside].sum() / inside.sum()
            centers.append(0.5 * (lo + hi)); effs.append(eff)
            errs.append(np.sqrt(max(eff * (1.0 - eff), 1.0 / inside.sum()) / inside.sum()))
        centers, effs, errs = map(np.asarray, (centers, effs, errs))
        coeffs, cov = np.polyfit(centers - peak, effs, 1, w=1.0 / errs, cov="unscaled")
        chi2 = float(np.sum(((np.polyval(coeffs, centers - peak) - effs) / errs) ** 2))
        mean_eff = float(np.average(effs, weights=1.0 / errs ** 2))
        entry = {
            "threshold": float(threshold),
            "background_efficiency": mean_eff,
            "slope_per_100MeV_relative": float(coeffs[0] * 0.1 / mean_eff),
            "slope_per_100MeV_relative_error": float(np.sqrt(cov[0, 0]) * 0.1 / mean_eff),
            "chi2_ndf": chi2 / max(len(centers) - 2, 1),
            "fit_bins": int(len(centers)),
            "line": [float(v) for v in coeffs],   # efficiency = line[0] * (m - peak) + line[1]
        }
        for region, mask in regions.items():
            prediction = float(np.mean(np.polyval(coeffs, bmass[mask] - peak)))
            gradient = np.stack([np.mean(bmass[mask] - peak), 1.0])
            prediction_error = float(np.sqrt(gradient @ cov @ gradient))
            eff = passed[mask].sum() / mask.sum()
            err = np.sqrt(max(eff * (1.0 - eff), 1.0 / mask.sum()) / mask.sum())
            ratio = eff / prediction
            entry[f"R_{region}"] = float(ratio)
            entry[f"R_{region}_error"] = float(ratio * np.hypot(err / max(eff, 1e-12), prediction_error / prediction))
            entry[f"{region}_entries"] = int(mask.sum())
        result[name] = entry
    return result


# =============================================================================
# Scores file and study metrics
# =============================================================================

def write_scores(path, splits, sculpt, spectators, metadata):
    """Scores of one trained model, the common input of the study metrics.

    splits: {"train"|"validation"|"test": dict of y, w, score, bmass, idx}
    sculpt: dict of bmass, score (DATA never used for training)
    spectators: {key: dict of w, score}
    metadata: dict written as metadata_json (window_background, features, method, ...)
    """
    with uproot.recreate(path) as out:
        for name, part in splits.items():
            out[name] = {k: np.asarray(part[k]) for k in ("y", "w", "score", "bmass", "idx")}
        out["sculpt"] = {k: np.asarray(sculpt[k], np.float64) for k in ("bmass", "score")}
        for key, part in spectators.items():
            out[f"spectator_{key}"] = {k: np.asarray(part[k], np.float64) for k in ("w", "score")}
        out["metadata_json"] = json.dumps(metadata, indent=2, sort_keys=True)


def read_scores(path):
    with uproot.open(path) as f:
        scores = {name.split(";")[0]: f[name].arrays(library="np") for name in f.keys()
                  if f[name].classname.startswith("TTree")}
        scores["metadata"] = json.loads(str(f["metadata_json"]))
    for part in ("train", "validation", "test"):
        scores[part] = {k: (np.asarray(v, np.float64) if k != "idx" else np.asarray(v, np.int64))
                        for k, v in scores[part].items()}
    return scores


def study_metrics(scores, c):
    """Every metric of one trained model of case c, from its scores file.

    working_point: the Punzi maximum on all signal MC and all DATA (punzi_curve on full_sample),
    the threshold and FOM optimalCUT_X_punzi.C finds for the model. validation: the four FOMs, the
    Punzi FOM on the validation split, with bootstrap errors, the overtraining test (train vs
    validation) and the efficiencies at the working point; test: the same on the test split (for the
    final report only); sculpting at the working point and at twice and half its background
    efficiency; spectator efficiencies at the working point.
    """
    # The window background with the current window and sidebands, from all stored DATA.
    window_bkg = window_background(all_data(scores)[0], c)
    passing, signal_eff, grid, fom = punzi_curve(*full_sample(scores), window_bkg, c)
    best = int(np.nanargmax(fom))
    working_point = {"threshold": float(grid[best]), "punzi": float(fom[best]), "signal_eff": float(signal_eff[best]),
                     "background_eff": float(passing[best]), "window_background": float(window_bkg * passing[best])}
    parts = {name: (scores[name]["y"], scores[name]["score"], scores[name]["w"], scores[name]["bmass"])
             for name in ("train", "validation", "test")}
    validation = evaluate(parts["train"], parts["validation"], window_bkg, c, threshold=working_point["threshold"])
    validation["errors"] = bootstrap_foms(*parts["validation"], window_bkg, c)
    test = evaluate(parts["train"], parts["test"], window_bkg, c, threshold=working_point["threshold"])
    # Looser and tighter: twice and half the working-point background efficiency, measured on
    # the untrained DATA of the sculpting fit regions (far more candidates than the validation).
    bmass, score = scores["sculpt"]["bmass"], scores["sculpt"]["score"]
    fit_region = ((bmass > c["sculpting_range"][0]) & (bmass < c["sculpting_range"][1])
                  & ~((bmass > c["sculpting_excluded"][0]) & (bmass < c["sculpting_excluded"][1])))
    # The thinned DATA counted as often as it was thinned (data_weights).
    background = np.repeat(score[fit_region], data_weights(bmass[fit_region], c))
    eff_wp = float(np.mean(background >= validation["working_point"]))
    thresholds = {
        "looser": float(np.quantile(background, 1.0 - min(2.0 * eff_wp, 1.0))),
        "working_point": validation["working_point"],
        "tighter": float(np.quantile(background, 1.0 - 0.5 * eff_wp)),
    }
    spectators = {}
    for name, part in scores.items():
        if name.startswith("spectator_"):
            eff, err = efficiency(part["score"], part["w"], validation["working_point"])
            spectators[name[len("spectator_"):]] = {"efficiency": eff, "error": err}
    return {
        "validation": validation,
        "test": test,
        "sculpting": sculpting(bmass, score, thresholds, c),
        "spectators": spectators,
        "working_point": working_point,
        "window_background": window_bkg,
        "fom_settings": fom_settings(c),
        "metadata": scores["metadata"],
    }
