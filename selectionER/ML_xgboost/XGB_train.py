import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot
import xgboost as xgb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import auc as sklearn_auc, confusion_matrix, roc_curve
import sklearn

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ML_common"))
import ml_fom
import ml_samples


# =============================================================================
# Configuration
# =============================================================================

# Inputs, pre-cuts and training sidebands are in ML_common/ml_samples.py; training reads the
# skims made by selectionER/make_skims.py. Output folders per sample:
OUTPUT_BASES = {
    'PbPb18': Path('xgb_outputs/ntmix_PbPb/pbpb18'),
    'PbPb23': Path('xgb_outputs/ntmix_PbPb/pbpb23'),
    'PbPb24': Path('xgb_outputs/ntmix_PbPb/pbpb24'),
    'ppRef24': Path('xgb_outputs/ntmix_ppRef'),
}
SAMPLE_CONFIGS = {sample: {**ml_samples.FLAT_INPUTS[sample], "output_dir": base} for sample, base in OUTPUT_BASES.items()}
DEFAULT_SAMPLE = 'PbPb23'
DEFAULT_FEATURES = list(ml_samples.FEATURES)
MC_WEIGHT_BRANCH = ml_samples.MC_WEIGHT_BRANCH

# Sample size. The train/validation/test split is set in ML_common/ml_fom.py.
MAX_SIGNAL = 165000
MAX_BACKGROUND = 500000
RANDOM_STATE = 42

# XGBoost training controls.
N_ROUNDS       = 1000
EARLY_STOPPING = 100
VERBOSE_EVAL   = 50

# Default hyperparameters, used without --params-json (a scan summary or a JSON dict).
XGB_PARAMS = {
    "colsample_bytree": 0.8,
    "eta": 0.05,
    "eval_metric": "auc",
    "gamma": 0.0,
    "max_depth": 4,
    "min_child_weight": 1.0,
    "nthread": 4,
    "objective": "binary:logistic",
    "reg_alpha": 0.0,
    "reg_lambda": 1.0,
    "subsample": 0.8,
}
# Output files.
OUTPUT_MODEL = 'xgb_X3872_vs_sideband.json'
OUTPUT_REPORT = 'training_report.root'
OUTPUT_SCORES = 'scores.root'

CURRENT_SAMPLE = DEFAULT_SAMPLE
CURRENT_FOM = None
CURRENT_FEATURES = DEFAULT_FEATURES
OUTPUT_DIR = None


# =============================================================================
# Training Workflow
# =============================================================================

def run_output_dir(sample_name, fom):
    """Output folder of one sample and model-selection FOM, e.g. xgb_outputs/ntmix_PbPb/pbpb23_auc."""
    base = SAMPLE_CONFIGS[sample_name]["output_dir"]
    return base.with_name(f"{base.name}_{fom}")


def configure_sample(sample_name, fom, features=None, output_dir=None):
    global CURRENT_SAMPLE, CURRENT_FOM, CURRENT_FEATURES, OUTPUT_DIR
    CURRENT_SAMPLE = sample_name
    CURRENT_FOM = fom
    CURRENT_FEATURES = list(features or DEFAULT_FEATURES)
    OUTPUT_DIR = Path(output_dir) if output_dir else run_output_dir(sample_name, fom)


def load_params(path):
    """Hyperparameters, boosting rounds and early stopping from a scan summary or a JSON dict."""
    content = json.loads(Path(path).read_text())
    if "best_params" not in content:
        return content, N_ROUNDS, EARLY_STOPPING
    config = content["config"]
    return content["best_params"], config["n_rounds"], config["early_stopping"]


def parse_args():
    parser = argparse.ArgumentParser(description="Train the XGBoost classifier of one ntmix sample.")
    parser.add_argument("--sample", choices=SAMPLE_CONFIGS, default=DEFAULT_SAMPLE)
    parser.add_argument("--fom", choices=ml_fom.FOMS, required=True)
    parser.add_argument("--features", help="Comma-separated inputs (default: ml_samples.FEATURES).")
    parser.add_argument("--params-json", help="Scan summary or JSON dict with the hyperparameters.")
    parser.add_argument("--output-dir", help="Output folder (default: xgb_outputs/.../<sample>_<fom>).")
    parser.add_argument("--max-signal", type=int, default=MAX_SIGNAL)
    parser.add_argument("--max-background", type=int, default=MAX_BACKGROUND)
    return parser.parse_args()


def main(sample_name, fom, features=None, params_json=None, output_dir=None,
         max_signal=MAX_SIGNAL, max_background=MAX_BACKGROUND):
    configure_sample(sample_name, fom, features, output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    params, n_rounds, early_stopping = (load_params(params_json) if params_json
                                        else (XGB_PARAMS, N_ROUNDS, EARLY_STOPPING))

    print("Training sample:", CURRENT_SAMPLE)
    print("Model selection FOM:", CURRENT_FOM)
    print("Features:", ", ".join(CURRENT_FEATURES))
    print("Hyperparameters:", params)
    data = ml_samples.load_training_data(CURRENT_SAMPLE, CURRENT_FEATURES, max_signal, max_background,
                                         RANDOM_STATE, with_extras=True)
    window_bkg = data["window_background"]
    print("Pre-cut of the skims:", data["pre_cut"])
    print("Signal / background candidates used:", data["stats"]["signal_used"]["entries"], "/", data["stats"]["background_used"])
    print("Background in the Punzi signal window before the classifier cut:", window_bkg)

    splits = {name: data[name] for name in ("train", "validation", "test")}
    dmatrix = {name: xgb.DMatrix(part["x"], label=part["y"], weight=part["w"], feature_names=CURRENT_FEATURES)
               for name, part in splits.items()}
    evals_result = {}

    # Early stopping on the validation split; the test split is only used below.
    model = xgb.train(
        params,
        dmatrix["train"],
        num_boost_round=n_rounds,
        evals=[(dmatrix["train"], "train"), (dmatrix["validation"], "validation")],
        early_stopping_rounds=early_stopping,
        evals_result=evals_result,
        verbose_eval=VERBOSE_EVAL,
    )
    print_train_validation_metric_diff(evals_result)

    scores = {name: predict_best(model, dmatrix[name]) for name in splits}
    triplet = lambda name: (splits[name]["y"], scores[name], splits[name]["w"], splits[name]["bmass"])
    validation = ml_fom.evaluate(triplet("train"), triplet("validation"), window_bkg)
    test = ml_fom.evaluate(triplet("train"), triplet("test"), window_bkg, threshold=validation["punzi_threshold"])
    ml_fom.print_evaluation("Validation", validation)
    ml_fom.print_evaluation("Test", test)
    print(f"Selected FOM {CURRENT_FOM} on validation: {validation[CURRENT_FOM]:.6f}")

    predict = lambda x: predict_best(model, xgb.DMatrix(x, feature_names=CURRENT_FEATURES))
    spectators = {key: {"w": part["w"], "score": predict(part["x"])} for key, part in data["spectators"].items()}
    sculpt = {"bmass": data["sculpt"]["bmass"], "score": predict(data["sculpt"]["x"])}

    model_path = OUTPUT_DIR / OUTPUT_MODEL
    model.save_model(model_path)
    print("Saved model to", model_path)

    tr, te = splits["train"], splits["test"]
    save_training_history(evals_result, model.best_iteration)
    save_roc_plot(te["y"], scores["test"], te["w"])
    importance = shap_importance(model, dmatrix["validation"])
    save_shap_importance(importance, CURRENT_FEATURES)
    save_score_distributions(scores["train"], tr["y"], tr["w"], scores["test"], te["y"], te["w"],
                             spectators["spectator"]["score"], spectators["spectator"]["w"])
    save_confusion_matrix_plot(te["y"], scores["test"], validation["punzi_threshold"], te["w"])
    save_training_report(
        model=model,
        importance=importance,
        evals_result=evals_result,
        score_train=scores["train"],
        y_train=tr["y"],
        w_train=tr["w"],
        score_test=scores["test"],
        y_test=te["y"],
        w_test=te["w"],
        score_spec=spectators["spectator"]["score"],
        spec_weights=spectators["spectator"]["w"],
        validation=validation,
        test=test,
        window_bkg=window_bkg,
        pre_cut=data["pre_cut"],
        params=params,
        n_rounds=n_rounds,
        early_stopping=early_stopping,
        max_signal=max_signal,
        max_background=max_background,
        sample_stats=data["stats"],
    )
    ml_fom.write_scores(OUTPUT_DIR / OUTPUT_SCORES,
                        {name: {**part, "score": scores[name]} for name, part in splits.items()},
                        sculpt, spectators, {
                            "method": "xgb", "sample": CURRENT_SAMPLE, "fom": CURRENT_FOM,
                            "features": CURRENT_FEATURES, "window_background": window_bkg,
                            "pre_cut": data["pre_cut"], "params": params, "random_state": RANDOM_STATE,
                            "best_iteration": int(model.best_iteration),
                            "importance": dict(zip(CURRENT_FEATURES, (importance / importance.sum()).tolist())),
                            "importance_method": "SHAP (mean |SHAP|, validation)",
                        })
    print("Saved scores to", OUTPUT_DIR / OUTPUT_SCORES)


def predict_best(model, dmatrix, **kwargs):
    best_iteration = getattr(model, "best_iteration", None)
    if best_iteration is not None:
        kwargs["iteration_range"] = (0, int(best_iteration) + 1)
    return model.predict(dmatrix, **kwargs)


def print_train_validation_metric_diff(evals_result):
    metric = next(iter(evals_result["train"]))
    final_train = evals_result["train"][metric][-1]
    final_validation = evals_result["validation"][metric][-1]
    print(
        f"Final train-validation {metric.upper()} difference "
        f"(round {len(evals_result['train'][metric]) - 1}): {final_train - final_validation:.6f} "
        f"(train={final_train:.6f}, validation={final_validation:.6f})"
    )


# =============================================================================
# Plot Helpers
# =============================================================================

def save_training_history(evals_result, best_iteration):
    metric = next(iter(evals_result["train"]))
    train_values = evals_result["train"][metric]
    validation_values = evals_result["validation"][metric]
    rounds = np.arange(1, len(train_values) + 1)
    best_index = len(train_values) - 1 if best_iteration is None else min(int(best_iteration), len(train_values) - 1)
    best_gap = train_values[best_index] - validation_values[best_index]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(rounds, train_values, label=f"Train {metric.upper()}", lw=2)
    ax.plot(rounds, validation_values, label=f"Validation {metric.upper()}", lw=2)
    ax.set_xlabel("Boosting round")
    ax.set_ylabel(metric.upper())
    ax.set_ylim(0.75, 1.0)
    ax.set_title("Training History")
    ax.grid(alpha=0.3)
    ax.legend(title=f"{metric.upper()} gap at best round {best_index + 1}: {best_gap:.4f}")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "training_history.pdf")
    plt.close(fig)


def save_roc_plot(y_test, score_test, w_test):
    fpr, tpr, _ = roc_curve(y_test, score_test, sample_weight=w_test)
    roc_auc = sklearn_auc(fpr, tpr)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(fpr, tpr, lw=2, label=f"ROC (AUC = {roc_auc:.4f})")
    ax.plot([0, 1], [0, 1], "--", color="gray")
    ax.set_xlabel("Background efficiency")
    ax.set_ylabel("Signal efficiency")
    ax.set_title("ROC Curve")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "roc_curve.pdf")
    plt.close(fig)


def shap_importance(model, dmatrix):
    """Mean |SHAP| of every input on the validation split."""
    return np.mean(np.abs(predict_best(model, dmatrix, pred_contribs=True)[:, :-1]), axis=0)


def save_shap_importance(mean_abs_shap, feature_names):
    shap_percent = 100.0 * mean_abs_shap / np.sum(mean_abs_shap)
    order = np.argsort(shap_percent)[::-1]
    ordered_features = np.array(feature_names)[order]
    ordered_percent = shap_percent[order]
    cumulative_percent = np.cumsum(ordered_percent)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(ordered_features, cumulative_percent, color="tab:green", alpha=0.35, label="Cumulative")
    ax.scatter(ordered_percent, ordered_features, color="black", zorder=3, label="Individual")
    for feature, individual, cumulative in zip(ordered_features, ordered_percent, cumulative_percent):
        ax.text(cumulative + 1.0, feature, f"{cumulative:.1f}%", va="center", fontsize=8)
        ax.text(individual + 1.0, feature, f"{individual:.1f}%", va="center", fontsize=8, color="black")
    ax.invert_yaxis()
    ax.set_xlabel("SHAP importance [%]")
    ax.set_title("Cumulative Feature Importance (SHAP)")
    ax.set_xlim(0.0, 108.0)
    ax.grid(axis="x", alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "shap_importance.pdf")
    plt.close(fig)


def save_score_distributions(score_train, y_train, w_train, score_test, y_test, w_test, score_spec, spec_weights):
    bins = np.linspace(0.0, 1.0, 21)
    bin_centers = 0.5 * (bins[1:] + bins[:-1])
    fig, (ax, ax_diff) = plt.subplots(
        2,
        1,
        figsize=(7, 7),
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05},
    )

    def normalized_weights(weights):
        weights = np.asarray(weights, dtype=np.float64)
        total = np.sum(weights)
        return weights / total if total > 0 else weights

    def normalized_hist(scores, weights):
        hist, _ = np.histogram(scores, bins=bins, weights=normalized_weights(weights))
        return hist

    def percent_difference(test_hist, train_hist):
        diff = np.full_like(train_hist, np.nan, dtype=np.float64)
        valid = train_hist > 0
        diff[valid] = 100.0 * (test_hist[valid] - train_hist[valid]) / train_hist[valid]
        return diff

    sig_train = y_train == 1
    bkg_train = y_train == 0
    sig_test = y_test == 1
    bkg_test = y_test == 0

    ks_signal = ml_fom.ks_statistic(
        score_train[sig_train],
        score_test[sig_test],
        w_train[sig_train],
        w_test[sig_test],
    )
    ks_background = ml_fom.ks_statistic(
        score_train[bkg_train],
        score_test[bkg_test],
        w_train[bkg_train],
        w_test[bkg_test],
    )

    sig_train_hist = normalized_hist(score_train[sig_train], w_train[sig_train])
    bkg_train_hist = normalized_hist(score_train[bkg_train], w_train[bkg_train])
    sig_test_hist = normalized_hist(score_test[sig_test], w_test[sig_test])
    bkg_test_hist = normalized_hist(score_test[bkg_test], w_test[bkg_test])

    ax.stairs(sig_train_hist, bins, lw=2, color="tab:orange", label="Signal train")
    ax.stairs(bkg_train_hist, bins, lw=2, color="tab:blue", label="Background train")
    ax.stairs(sig_test_hist, bins, fill=True, alpha=0.25, color="tab:orange", label="Signal test")
    ax.stairs(bkg_test_hist, bins, fill=True, alpha=0.25, color="tab:blue", label="Background test")
    if score_spec is not None and len(score_spec) > 0:
        normalized_spec_weights = normalized_weights(spec_weights)
        ax.hist(score_spec, bins=bins, weights=normalized_spec_weights, histtype="step", lw=2, linestyle="--", color="gold", label="Psi(2S) spectator")

    ax_diff.axhline(0.0, color="gray", lw=1)
    ax_diff.plot(bin_centers, percent_difference(sig_test_hist, sig_train_hist), marker="o", ms=3, lw=1.5, color="tab:orange", label="Signal")
    ax_diff.plot(bin_centers, percent_difference(bkg_test_hist, bkg_train_hist), marker="o", ms=3, lw=1.5, color="tab:blue", label="Background")

    ax.set_ylabel("Normalized entries")
    ax.set_title("Score Distributions")
    ax.grid(alpha=0.3)
    distribution_legend = ax.legend(fontsize=9, loc="upper left")
    ax.add_artist(distribution_legend)
    ks_handles = [
        Line2D([], [], linestyle="none", label=f"Signal: {ks_signal:.4f}"),
        Line2D([], [], linestyle="none", label=f"Background: {ks_background:.4f}"),
    ]
    ax.legend(
        handles=ks_handles,
        fontsize=9,
        title="Weighted KS",
        title_fontsize=9,
        loc="upper left",
        bbox_to_anchor=(0.48, 1.0),
        handlelength=0,
        handletextpad=0,
    )
    ax_diff.set_xlabel("XGB score")
    ax_diff.set_ylabel("(test-train)/train [%]")
    ax_diff.grid(alpha=0.3)
    ax_diff.legend(fontsize=8, loc="best")
    fig.savefig(OUTPUT_DIR / "score_distributions.pdf")
    plt.close(fig)


def save_confusion_matrix_plot(y_test, score_test, threshold, w_test):
    y_pred = (score_test >= threshold).astype(int)
    cm = confusion_matrix(y_test, y_pred, sample_weight=w_test, labels=[0, 1])

    fig, ax = plt.subplots(figsize=(5, 4))
    im = ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Pred. bkg", "Pred. sig"])
    ax.set_yticklabels(["True bkg", "True sig"])
    ax.set_title(f"Confusion Matrix (thr = {threshold:.3f})")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{cm[i, j]:.3f}", ha="center", va="center", color="black")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "confusion_matrix.pdf")
    plt.close(fig)


# =============================================================================
# Persistent Training Report
# =============================================================================

def save_training_report(
    model,
    importance,
    evals_result,
    score_train,
    y_train,
    w_train,
    score_test,
    y_test,
    w_test,
    score_spec,
    spec_weights,
    validation,
    test,
    window_bkg,
    pre_cut,
    params,
    n_rounds,
    early_stopping,
    max_signal,
    max_background,
    sample_stats,
):
    threshold = validation["punzi_threshold"]
    signal_train = y_train == 1
    signal_test = y_test == 1
    background_train = y_train == 0
    background_test = y_test == 0
    fpr, tpr, roc_thresholds = roc_curve(y_test, score_test, sample_weight=w_test)
    cm = confusion_matrix(
        y_test,
        (score_test >= threshold).astype(np.int8),
        sample_weight=w_test,
        labels=[0, 1],
    )

    mean_abs_shap = importance
    shap_fraction = mean_abs_shap / np.sum(mean_abs_shap)

    bins = np.linspace(0.0, 1.0, 26)
    def normalized_hist(values, weights):
        weights = np.asarray(weights, dtype=np.float64)
        weight_sum = np.sum(weights)
        if weight_sum > 0.0:
            weights = weights / weight_sum
        return np.histogram(values, bins=bins, weights=weights)[0]

    score_histograms = {
        "bin_low": bins[:-1],
        "bin_high": bins[1:],
        "signal_train": normalized_hist(score_train[signal_train], w_train[signal_train]),
        "background_train": normalized_hist(score_train[background_train], w_train[background_train]),
        "signal_test": normalized_hist(score_test[signal_test], w_test[signal_test]),
        "background_test": normalized_hist(score_test[background_test], w_test[background_test]),
    }
    if score_spec is not None:
        score_histograms["psi2s_spectator"] = normalized_hist(score_spec, spec_weights)

    metric_name = next(iter(evals_result["train"]))
    metadata = {
        "format_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "python_version": sys.version,
        "package_versions": {
            "numpy": np.__version__,
            "uproot": uproot.__version__,
            "xgboost": xgb.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "sample": CURRENT_SAMPLE,
        "fom": CURRENT_FOM,
        "fom_settings": ml_fom.settings(),
        "window_background": window_bkg,
        "skims": {key: ml_samples.skim_path(CURRENT_SAMPLE, key) for key in ml_samples.TREES},
        "trees": dict(ml_samples.TREES),
        "features": list(CURRENT_FEATURES),
        "signal_cut": pre_cut,
        "background_sideband_cut": ml_samples.sideband_cut(),
        "mc_weight_branch": MC_WEIGHT_BRANCH,
        "class_weighting": "pThat-weighted signal and unit-weight data, normalized to equal class sums",
        "random_state": RANDOM_STATE,
        "max_signal": max_signal,
        "max_background": max_background,
        "n_rounds": n_rounds,
        "early_stopping": early_stopping,
        "xgb_params": params,
        "model_file": OUTPUT_MODEL,
        "report_file": OUTPUT_REPORT,
        "sample_statistics": sample_stats,
    }
    # validation_*: model selection (train vs validation); test_*: final check (train vs test)
    # at the validation working point.
    metrics = {f"validation_{name}": float(value) for name, value in validation.items()}
    metrics.update({f"test_{name}": float(value) for name, value in test.items()})
    metrics.update({
        "best_iteration": float(getattr(model, "best_iteration", -1)),
        "best_score": float(getattr(model, "best_score", np.nan)),
        "train_entries": float(len(y_train)),
        "test_entries": float(len(y_test)),
        "train_weight_sum": float(np.sum(w_train)),
        "test_weight_sum": float(np.sum(w_test)),
    })

    report_path = OUTPUT_DIR / OUTPUT_REPORT
    with uproot.recreate(report_path) as report:
        report["metadata_json"] = json.dumps(metadata, indent=2, sort_keys=True)
        report["metrics"] = {
            name: np.asarray([value], dtype=np.float64)
            for name, value in metrics.items()
        }
        report["training_history"] = {
            "round": np.arange(1, len(evals_result["train"][metric_name]) + 1, dtype=np.int32),
            f"train_{metric_name}": np.asarray(evals_result["train"][metric_name], dtype=np.float64),
            f"validation_{metric_name}": np.asarray(evals_result["validation"][metric_name], dtype=np.float64),
        }
        report["roc_curve"] = {
            "false_positive_rate": np.asarray(fpr, dtype=np.float64),
            "true_positive_rate": np.asarray(tpr, dtype=np.float64),
            "threshold": np.asarray(roc_thresholds, dtype=np.float64),
        }
        report["confusion_matrix"] = {
            "true_class": np.repeat(np.asarray([0, 1], dtype=np.int8), 2),
            "predicted_class": np.tile(np.asarray([0, 1], dtype=np.int8), 2),
            "weighted_count": np.asarray(cm, dtype=np.float64).reshape(-1),
        }
        report["feature_importance"] = {
            "feature_index": np.arange(len(CURRENT_FEATURES), dtype=np.int32),
            "mean_abs_shap": np.asarray(mean_abs_shap, dtype=np.float64),
            "shap_fraction": np.asarray(shap_fraction, dtype=np.float64),
        }
        report["feature_names_json"] = json.dumps(CURRENT_FEATURES)
        report["score_distributions"] = {
            name: np.asarray(values, dtype=np.float64)
            for name, values in score_histograms.items()
        }
    print("Saved training report to", report_path)


# =============================================================================
# Script Entry Point
# =============================================================================

if __name__ == "__main__":
    args = parse_args()
    main(args.sample, args.fom, args.features.split(",") if args.features else None, args.params_json,
         args.output_dir, args.max_signal, args.max_background)
