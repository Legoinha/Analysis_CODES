import argparse
import json
import os
import re

import numpy as np
import optuna
import xgboost as xgb

import XGB_train as train_cfg
from XGB_train import ml_fom, ml_samples


# =============================================================================
# Configuration
# =============================================================================

# XGB_train.py and ML_common/ml_samples.py are the single source of truth for samples,
# features and cuts; ML_common/ml_fom.py for the split, the figures of merit, and the
# overtraining test. Training data come from the skims of selectionER/make_skims.py.
SAMPLE_CONFIGS = train_cfg.SAMPLE_CONFIGS
DEFAULT_SAMPLE = train_cfg.DEFAULT_SAMPLE
DEFAULT_N_TRIALS = 50
DEFAULT_FEATURES = list(train_cfg.DEFAULT_FEATURES)
DEFAULT_MAX_SIGNAL = 165000
DEFAULT_MAX_BACKGROUND = 500000
DEFAULT_RANDOM_STATE = 42
DEFAULT_N_ROUNDS = 3000
DEFAULT_EARLY_STOPPING = 100
DEFAULT_NTHREAD = 4
OUTPUT_MODEL = "xgb_X3872_vs_sideband.json"
DEFAULT_SUMMARY_SUBDIR = "optuna_summaries"


# =============================================================================
# Data Helpers
# =============================================================================

def load_training_arrays(sample_name, features, max_signal, max_background):
    data = ml_samples.load_training_data(sample_name, features, max_signal, max_background, DEFAULT_RANDOM_STATE)
    print("Pre-cut of the skims:", data["pre_cut"])
    print("Signal / background candidates used:", data["stats"]["signal_used"]["entries"], "/", data["stats"]["background_used"])
    print("Background in the Punzi signal window:", data["window_background"])
    return data


# =============================================================================
# Optuna Objective
# =============================================================================

def build_objective(data, features, fom, n_rounds, early_stopping):
    train, validation, window_bkg = data["train"], data["validation"], data["window_background"]
    x_train, y_train, w_train = train["x"], train["y"], train["w"]
    x_val, y_val, w_val = validation["x"], validation["y"], validation["w"]
    dtrain = xgb.DMatrix(x_train, label=y_train, weight=w_train, feature_names=features)
    dval = xgb.DMatrix(x_val, label=y_val, weight=w_val, feature_names=features)

    def objective(trial):
        params = {
            "objective": "binary:logistic",
            "eval_metric": "auc",
            "max_depth": trial.suggest_int("max_depth", 2, 5),
            "eta": trial.suggest_float("eta", 0.005, 0.12, log=True),
            "min_child_weight": trial.suggest_float("min_child_weight", 0.5, 20.0, log=True),
            "gamma": trial.suggest_float("gamma", 0.0, 8.0),
            "subsample": trial.suggest_float("subsample", 0.55, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.55, 1.0),
            "reg_lambda": trial.suggest_float("reg_lambda", 0.5, 50.0, log=True),
            "reg_alpha": trial.suggest_float("reg_alpha", 0.0, 10.0),
            "nthread": DEFAULT_NTHREAD,
        }

        model = xgb.train(
            params,
            dtrain,
            num_boost_round=n_rounds,
            evals=[(dtrain, "train"), (dval, "validation")],
            early_stopping_rounds=early_stopping,
            verbose_eval=False,
        )

        result = ml_fom.evaluate(
            (y_train, train_cfg.predict_best(model, dtrain), w_train, train["bmass"]),
            (y_val, train_cfg.predict_best(model, dval), w_val, validation["bmass"]),
            window_bkg,
        )
        for name, value in result.items():
            trial.set_user_attr(name, value)
        trial.set_user_attr("best_iteration", int(model.best_iteration) if model.best_iteration is not None else None)
        # Overtraining is a pass/fail requirement: a failed trial never becomes the best one.
        if not result["passed"]:
            raise optuna.TrialPruned()
        return result[fom]

    return objective


# =============================================================================
# Output Helpers
# =============================================================================

def best_params_from_trial(trial):
    return {
        "objective": "binary:logistic",
        "eval_metric": "auc",
        **trial.params,
        "nthread": DEFAULT_NTHREAD,
    }


def build_summary(study, args, features, pre_cut):
    # All trials of a job can fail the overtraining test: the summary then has no best trial.
    completed = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
    best_trial = study.best_trial if completed else None
    rejected = [t for t in study.trials if t.user_attrs.get("passed") is False]
    return {
        "sample": args.sample,
        "fom": args.fom,
        "seed": args.seed,
        "best_value": best_trial.value if best_trial else None,
        "best_params": best_params_from_trial(best_trial) if best_trial else None,
        "best_trial": {
            "number": best_trial.number,
            "params": best_trial.params,
            "user_attrs": best_trial.user_attrs,
        } if best_trial else None,
        "n_trials": len(study.trials),
        # Every trial, in order: the first ones of each job are TPE's random start-up trials.
        "trials": [{"number": t.number, "state": t.state.name, "value": t.value, "params": t.params,
                    "user_attrs": t.user_attrs} for t in study.trials],
        "n_rejected_overtraining": len(rejected),
        "config": {
            "sample_configs": {
                name: {**cfg, "output_dir": str(cfg["output_dir"])}
                for name, cfg in SAMPLE_CONFIGS.items()
            },
            "default_sample": args.sample,
            "output_dir": str(train_cfg.run_output_dir(args.sample, args.fom)),
            "trees": dict(ml_samples.TREES),
            "features": features,
            "scan_tag": args.scan_tag,
            "summary_subdir": args.summary_subdir,
            "pre_cut": pre_cut,
            "mc_weight_branch": ml_samples.MC_WEIGHT_BRANCH,
            "training_sidebands": [list(s) for s in ml_samples.TRAINING_SIDEBANDS],
            "max_signal": args.max_signal,
            "max_background": args.max_background,
            "random_state": DEFAULT_RANDOM_STATE,
            "fom_settings": ml_fom.settings(),
            "n_rounds": args.n_rounds,
            "early_stopping": args.early_stopping,
            "output_model": OUTPUT_MODEL,
        },
    }


def save_summary(summary, summary_path=None):
    if summary_path:
        train_cfg.Path(summary_path).parent.mkdir(parents=True, exist_ok=True)
        train_cfg.Path(summary_path).write_text(json.dumps(summary, indent=2, sort_keys=True))
        print("Saved Optuna summary:", summary_path)
        return summary_path
    output_dir = train_cfg.Path(summary["config"]["output_dir"])
    summary_subdir = summary["config"].get("summary_subdir", "")
    if summary_subdir:
        output_dir = output_dir / summary_subdir
    output_dir.mkdir(parents=True, exist_ok=True)
    tag = summary["config"].get("scan_tag", "")
    tag_part = f"_{sanitize_filename_part(tag)}" if tag else ""
    output_path = output_dir / f"optuna_summary_{summary['sample']}_{summary['fom']}{tag_part}_seed{summary['seed']}.json"
    output_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
    print("Saved Optuna summary:", output_path)
    return output_path


def sanitize_filename_part(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


# =============================================================================
# Application Workflow
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(description="Tune XGB_train.py hyperparameters with Optuna.")
    parser.add_argument("--sample", choices=SAMPLE_CONFIGS, default=DEFAULT_SAMPLE)
    parser.add_argument("--fom", choices=ml_fom.FOMS, required=True,
                        help="Validation figure of merit to maximize.")
    parser.add_argument("--trials", type=int, default=DEFAULT_N_TRIALS)
    parser.add_argument("--max-signal", type=int, default=DEFAULT_MAX_SIGNAL)
    parser.add_argument("--max-background", type=int, default=DEFAULT_MAX_BACKGROUND)
    parser.add_argument("--n-rounds", type=int, default=DEFAULT_N_ROUNDS)
    parser.add_argument("--early-stopping", type=int, default=DEFAULT_EARLY_STOPPING)
    parser.add_argument("--seed", type=int, default=DEFAULT_RANDOM_STATE)
    parser.add_argument("--features", help="Comma-separated inputs (default: ml_samples.FEATURES).")
    parser.add_argument("--summary-path", help="Write the summary JSON here instead of the sample output folder.")
    parser.add_argument(
        "--scan-tag",
        default="",
        help="Optional tag inserted into the summary JSON filename.",
    )
    parser.add_argument(
        "--summary-subdir",
        default=DEFAULT_SUMMARY_SUBDIR,
        help="Optional subdirectory below the sample output directory for summary JSON files.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    batch_quiet = os.environ.get("XGB_BATCH_QUIET", "").lower() in {"1", "true", "yes"}
    if batch_quiet:
        optuna.logging.set_verbosity(optuna.logging.WARNING)

    features = args.features.split(",") if args.features else DEFAULT_FEATURES
    data = load_training_arrays(args.sample, features, args.max_signal, args.max_background)

    sampler = optuna.samplers.TPESampler(seed=args.seed)
    study = optuna.create_study(direction="maximize", sampler=sampler)
    study.optimize(
        build_objective(data, features, args.fom, args.n_rounds, args.early_stopping),
        n_trials=args.trials,
        show_progress_bar=not batch_quiet,
    )

    summary = build_summary(study, args, features, data["pre_cut"])
    save_summary(summary, args.summary_path)
    if summary["best_trial"] is None:
        print("No trial passed the overtraining test:", summary["n_rejected_overtraining"], "of", summary["n_trials"], "rejected")
        return
    best = study.best_trial.user_attrs
    print(f"Best validation {args.fom}:", study.best_value)
    print("Best parameters:", study.best_trial.params)
    print("Best validation AUC:", best["auc"], "train-validation AUC gap:", best["auc_gap"])
    print("Best KS p-values: signal", best["ks_signal_pvalue"], "background", best["ks_background_pvalue"])
    print("Trials rejected by the overtraining test:", summary["n_rejected_overtraining"], "of", summary["n_trials"])


# =============================================================================
# Script Entry Point
# =============================================================================

if __name__ == "__main__":
    main()
