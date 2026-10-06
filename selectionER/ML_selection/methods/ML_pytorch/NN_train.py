import argparse
import copy
import json
import platform
import sys
from datetime import datetime, timezone
import numpy as np
import uproot
import torch
from torch import nn
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import auc as sklearn_auc, confusion_matrix, roc_auc_score, roc_curve
import sklearn

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ML_common"))
import ml_fom
import ml_samples
import ml_plots
from ml_samples import settings


# =============================================================================
# Configuration
# =============================================================================

# Cases, inputs, pre-cuts and training sidebands are in ML_selection/settings.py; training reads the
# skims made by ML_selection/scripts/make_skims.py. Outputs in nn_outputs/<tree>/<year>_<fom>/.
# Features passed through log1p before standardization (heavy-tailed inputs) are settings.LOG_FEATURES,
# by flat-tree name: BmuLeadPt_dc is logged when BmuLeadPt is listed. The transform is stored inside the model.
# Sample size: settings.MAX_SIGNAL and MAX_BACKGROUND. The train/validation/test split is set in ML_common/ml_fom.py.
RANDOM_STATE = 42
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# NN training controls.
MAX_EPOCHS     = 150
EARLY_STOPPING = 15
VERBOSE_EVAL   = 5
NTHREAD        = 4

# NN model hyperparameters.
NN_PARAMS = {
    "activation": "silu",
    "batch_size": 256,
    "dropout": 0.029017660353167773,
    "hidden_layers": 4,
    "hidden_units": 64,
    "learning_rate": 0.0014753881417661754,
    "lr_factor": 0.5,
    "lr_patience": 4,
    "weight_decay": 3.638079234932926e-06,
}
# Output files.
OUTPUT_MODEL = 'nn_X3872_vs_sideband.pt'
OUTPUT_REPORT = 'training_report.root'
OUTPUT_SCORES = 'scores.root'

CURRENT_CASE = None
CURRENT_FOM = None
CURRENT_FEATURES = None
OUTPUT_DIR = None


# =============================================================================
# Model Definition
# =============================================================================

ACTIVATIONS = {
    "relu": nn.ReLU,
    "silu": nn.SiLU,
    "elu": nn.ELU,
    "tanh": nn.Tanh,
    "leaky_relu": nn.LeakyReLU,
}


class X3872Classifier(nn.Module):
    """Fully connected binary classifier.

    Input preprocessing (log1p on selected features, then standardization) is
    stored as buffers, so the saved model takes raw branch values as input.
    forward() returns logits; predict_proba() returns scores in [0, 1].
    """

    def __init__(self, n_features, hidden_layers, hidden_units, dropout, activation):
        super().__init__()
        self.register_buffer("log_mask", torch.zeros(n_features, dtype=torch.bool))
        self.register_buffer("feature_mean", torch.zeros(n_features))
        self.register_buffer("feature_std", torch.ones(n_features))

        layers = []
        width_in = n_features
        for _ in range(hidden_layers):
            layers += [
                nn.Linear(width_in, hidden_units),
                nn.BatchNorm1d(hidden_units),
                ACTIVATIONS[activation](),
                nn.Dropout(dropout),
            ]
            width_in = hidden_units
        layers.append(nn.Linear(width_in, 1))
        self.network = nn.Sequential(*layers)

    def set_preprocessing(self, log_mask, mean, std):
        self.log_mask.copy_(torch.as_tensor(log_mask, dtype=torch.bool))
        self.feature_mean.copy_(torch.as_tensor(mean, dtype=torch.float32))
        self.feature_std.copy_(torch.as_tensor(std, dtype=torch.float32))

    def preprocess(self, x):
        x = torch.where(self.log_mask, torch.log1p(torch.clamp(x, min=0.0)), x)
        return (x - self.feature_mean) / self.feature_std

    def forward(self, x):
        return self.network(self.preprocess(x)).squeeze(-1)

    def predict_proba(self, x):
        return torch.sigmoid(self.forward(x))


def build_model(n_features, params):
    return X3872Classifier(
        n_features=n_features,
        hidden_layers=int(params["hidden_layers"]),
        hidden_units=int(params["hidden_units"]),
        dropout=float(params["dropout"]),
        activation=params["activation"],
    )


def fit_preprocessing(model, x_train, feature_names):
    log_mask = np.array([name.removesuffix("_dc") in settings.LOG_FEATURES for name in feature_names], dtype=bool)
    x = np.asarray(x_train, dtype=np.float64).copy()
    x[:, log_mask] = np.log1p(np.clip(x[:, log_mask], 0.0, None))
    mean = x.mean(axis=0)
    std = x.std(axis=0)
    std[std == 0.0] = 1.0
    model.set_preprocessing(log_mask, mean, std)


def save_model(model, path, params, feature_names, best_epoch):
    torch.save(
        {
            "format_version": 1,
            "state_dict": model.state_dict(),
            "nn_params": dict(params),
            "feature_names": list(feature_names),
            "log_features": list(settings.LOG_FEATURES),
            "best_epoch": int(best_epoch),
            "sample": CURRENT_CASE["name"],
        },
        path,
    )


def load_model(path):
    checkpoint = torch.load(path, map_location=DEVICE, weights_only=True)
    model = build_model(len(checkpoint["feature_names"]), checkpoint["nn_params"])
    model.load_state_dict(checkpoint["state_dict"])
    model.to(DEVICE).eval()
    return model, checkpoint


@torch.no_grad()
def predict_scores(model, x, device=None, batch_size=65536):
    device = device or next(model.parameters()).device
    model.eval()
    x_tensor = torch.as_tensor(np.asarray(x, dtype=np.float32))
    scores = []
    for start in range(0, len(x_tensor), batch_size):
        batch = x_tensor[start:start + batch_size].to(device)
        scores.append(model.predict_proba(batch).cpu())
    return torch.cat(scores).numpy().astype(np.float32)


@torch.no_grad()
def evaluate_loss(model, x_tensor, y_tensor, w_tensor, device, batch_size=65536):
    model.eval()
    loss_sum = 0.0
    weight_sum = 0.0
    for start in range(0, len(x_tensor), batch_size):
        xb = x_tensor[start:start + batch_size].to(device)
        yb = y_tensor[start:start + batch_size].to(device)
        wb = w_tensor[start:start + batch_size].to(device)
        loss = nn.functional.binary_cross_entropy_with_logits(model(xb), yb, reduction="none")
        loss_sum += float(torch.sum(loss * wb))
        weight_sum += float(torch.sum(wb))
    return loss_sum / weight_sum


def train_network(
    x_train, y_train, w_train, x_val, y_val, w_val,
    params, max_epochs, early_stopping, feature_names,
    verbose_eval=0, seed=RANDOM_STATE, epoch_callback=None,
):
    """Train one network with early stopping on the weighted validation AUC.

    Returns (model, history, best_epoch). The returned model carries the
    weights from the best epoch. epoch_callback(epoch, validation_scores) may
    raise to stop training (used for Optuna pruning).
    """
    np.random.seed(seed)
    torch.manual_seed(seed)
    device = DEVICE
    model = build_model(len(feature_names), params)
    fit_preprocessing(model, x_train, feature_names)
    model.to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(params["learning_rate"]),
        weight_decay=float(params["weight_decay"]),
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="max",
        factor=float(params["lr_factor"]),
        patience=int(params["lr_patience"]),
    )

    x_train_t = torch.as_tensor(np.asarray(x_train, dtype=np.float32))
    y_train_t = torch.as_tensor(np.asarray(y_train, dtype=np.float32))
    w_train_t = torch.as_tensor(np.asarray(w_train, dtype=np.float32))
    x_val_t = torch.as_tensor(np.asarray(x_val, dtype=np.float32))
    y_val_t = torch.as_tensor(np.asarray(y_val, dtype=np.float32))
    w_val_t = torch.as_tensor(np.asarray(w_val, dtype=np.float32))

    batch_size = int(params["batch_size"])
    generator = torch.Generator().manual_seed(seed)
    history = {
        "epoch": [],
        "train_loss": [],
        "validation_loss": [],
        "train_auc": [],
        "validation_auc": [],
        "learning_rate": [],
    }
    best_auc = -np.inf
    best_epoch = 0
    best_state = copy.deepcopy(model.state_dict())
    epochs_without_improvement = 0

    for epoch in range(1, max_epochs + 1):
        model.train()
        permutation = torch.randperm(len(x_train_t), generator=generator)
        for start in range(0, len(permutation), batch_size):
            index = permutation[start:start + batch_size]
            if len(index) < 2:
                continue  # BatchNorm needs at least two entries in training mode.
            xb = x_train_t[index].to(device)
            yb = y_train_t[index].to(device)
            wb = w_train_t[index].to(device)
            optimizer.zero_grad(set_to_none=True)
            loss = nn.functional.binary_cross_entropy_with_logits(model(xb), yb, reduction="none")
            loss = torch.sum(loss * wb) / torch.sum(wb)
            loss.backward()
            optimizer.step()

        score_train = predict_scores(model, x_train, device)
        score_val = predict_scores(model, x_val, device)
        train_auc = roc_auc_score(y_train, score_train, sample_weight=w_train)
        val_auc = roc_auc_score(y_val, score_val, sample_weight=w_val)
        history["epoch"].append(epoch)
        history["train_loss"].append(evaluate_loss(model, x_train_t, y_train_t, w_train_t, device))
        history["validation_loss"].append(evaluate_loss(model, x_val_t, y_val_t, w_val_t, device))
        history["train_auc"].append(float(train_auc))
        history["validation_auc"].append(float(val_auc))
        history["learning_rate"].append(float(optimizer.param_groups[0]["lr"]))
        scheduler.step(val_auc)

        if val_auc > best_auc:
            best_auc = val_auc
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            epochs_without_improvement = 0
        else:
            epochs_without_improvement += 1

        if verbose_eval and (epoch == 1 or epoch % verbose_eval == 0):
            print(
                f"[{epoch}]\ttrain-loss:{history['train_loss'][-1]:.5f}"
                f"\tvalidation-loss:{history['validation_loss'][-1]:.5f}"
                f"\ttrain-auc:{train_auc:.5f}\tvalidation-auc:{val_auc:.5f}"
                f"\tlr:{history['learning_rate'][-1]:.2e}"
            )
        if epoch_callback is not None:
            epoch_callback(epoch, score_val)
        if epochs_without_improvement >= early_stopping:
            if verbose_eval:
                print(f"Early stopping at epoch {epoch}; best epoch {best_epoch} (validation AUC {best_auc:.5f})")
            break

    model.load_state_dict(best_state)
    model.eval()
    return model, history, best_epoch


# =============================================================================
# Training Workflow
# =============================================================================

def configure_case(case_name, fom, features=None, output_dir=None):
    """The case, FOM, inputs (default: its production features) and output folder of this run (default:
    nn_outputs/<tree>/<year>_<fom>)."""
    global CURRENT_CASE, CURRENT_FOM, CURRENT_FEATURES, OUTPUT_DIR
    CURRENT_CASE = ml_samples.case(case_name)
    CURRENT_FOM = fom
    CURRENT_FEATURES = list(features or CURRENT_CASE["features"])
    OUTPUT_DIR = Path(output_dir or f"nn_outputs/{CURRENT_CASE['tree']}/{CURRENT_CASE['year']}_{fom}")


def load_params(path):
    """Hyperparameters, maximum epochs, early stopping and training seed from a scan summary
    (the seed of its job, so the refit reproduces the chosen trial) or a JSON dict."""
    content = json.loads(Path(path).read_text())
    if "best_params" not in content:
        return content, MAX_EPOCHS, EARLY_STOPPING, RANDOM_STATE
    config = content["config"]
    return content["best_params"], config["max_epochs"], config["early_stopping"], content["seed"]


def parse_args():
    parser = argparse.ArgumentParser(description="Train the PyTorch NN classifier of one case.")
    parser.add_argument("--case", choices=ml_samples.CASES, required=True, help="<tree>/<year>, e.g. ntKp/PbPb23.")
    parser.add_argument("--fom", choices=ml_fom.FOMS, required=True)
    parser.add_argument("--features", help="Comma-separated inputs (default: the production features of the tree).")
    parser.add_argument("--params-json", help="Scan summary or JSON dict with the hyperparameters.")
    parser.add_argument("--output-dir", help="Output folder (default: nn_outputs/<tree>/<year>_<fom>).")
    parser.add_argument("--max-signal", type=int, default=settings.MAX_SIGNAL)
    parser.add_argument("--max-background", type=int, default=settings.MAX_BACKGROUND)
    return parser.parse_args()


def main(case_name, fom, features=None, params_json=None, output_dir=None,
         max_signal=settings.MAX_SIGNAL, max_background=settings.MAX_BACKGROUND):
    configure_case(case_name, fom, features, output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(NTHREAD)
    params, max_epochs, early_stopping, training_seed = (load_params(params_json) if params_json
                                                         else (NN_PARAMS, MAX_EPOCHS, EARLY_STOPPING, RANDOM_STATE))

    print("Training case:", CURRENT_CASE["name"])
    print("Model selection FOM:", CURRENT_FOM)
    print("Features:", ", ".join(CURRENT_FEATURES))
    print("Hyperparameters:", params)
    print("Device:", DEVICE)
    data = ml_samples.load_training_data(CURRENT_CASE, CURRENT_FEATURES, max_signal, max_background,
                                         RANDOM_STATE, with_extras=True)
    window_bkg = data["window_background"]
    print("Pre-cut of the skims:", data["pre_cut"])
    print("Signal / background candidates used:", data["stats"]["signal_used"]["entries"], "/", data["stats"]["background_used"])
    print("Background in the Punzi signal window before the classifier cut:", window_bkg)
    splits = {name: data[name] for name in ("train", "validation", "test")}
    tr, va, te = splits["train"], splits["validation"], splits["test"]

    # Early stopping on the validation split; the test split is only used below.
    model, history, best_epoch = train_network(
        tr["x"], tr["y"], tr["w"], va["x"], va["y"], va["w"],
        params=params,
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        feature_names=CURRENT_FEATURES,
        verbose_eval=VERBOSE_EVAL,
        seed=training_seed,
    )
    print_train_validation_metric_diff(history, best_epoch)

    scores = {name: predict_scores(model, part["x"]) for name, part in splits.items()}
    parts = {name: (part["y"], scores[name], part["w"], part["bmass"]) for name, part in splits.items()}
    validation = ml_fom.evaluate(parts["train"], parts["validation"], window_bkg, CURRENT_CASE)
    test = ml_fom.evaluate(parts["train"], parts["test"], window_bkg, CURRENT_CASE, threshold=validation["punzi_threshold"])
    ml_fom.print_evaluation("Validation", validation)
    ml_fom.print_evaluation("Test", test)
    print(f"Selected FOM {CURRENT_FOM} on validation: {validation[CURRENT_FOM]:.6f}")

    spectators = {key: {"w": part["w"], "score": predict_scores(model, part["x"])} for key, part in data["spectators"].items()}
    sculpt = {"bmass": data["sculpt"]["bmass"], "score": predict_scores(model, data["sculpt"]["x"])}

    model_path = OUTPUT_DIR / OUTPUT_MODEL
    save_model(model, model_path, params, CURRENT_FEATURES, best_epoch)
    print("Saved model to", model_path)

    importance = permutation_importance(model, va["x"], va["y"], va["w"])
    save_training_history(history, best_epoch)
    save_roc_plot(te["y"], scores["test"], te["w"])
    save_feature_importance(importance, CURRENT_FEATURES)
    save_score_distributions(scores["train"], tr["y"], tr["w"], scores["test"], te["y"], te["w"], spectators)
    save_confusion_matrix_plot(te["y"], scores["test"], validation["punzi_threshold"], te["w"])
    save_training_report(
        model=model,
        history=history,
        best_epoch=best_epoch,
        importance=importance,
        score_train=scores["train"],
        y_train=tr["y"],
        w_train=tr["w"],
        score_test=scores["test"],
        y_test=te["y"],
        w_test=te["w"],
        spectators=spectators,
        validation=validation,
        test=test,
        window_bkg=window_bkg,
        pre_cut=data["pre_cut"],
        params=params,
        max_signal=max_signal,
        max_background=max_background,
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        sample_stats=data["stats"],
    )
    ml_fom.write_scores(OUTPUT_DIR / OUTPUT_SCORES,
                        {name: {**part, "score": scores[name]} for name, part in splits.items()},
                        sculpt, spectators, {
                            "method": "nn", "sample": CURRENT_CASE["name"], "fom": CURRENT_FOM,
                            "features": CURRENT_FEATURES, "window_background": window_bkg,
                            "pre_cut": data["pre_cut"], "params": params, "random_state": RANDOM_STATE,
                            "training_seed": training_seed, "best_epoch": int(best_epoch),
                            "importance": dict(zip(CURRENT_FEATURES, (importance / importance.sum()).tolist())),
                            "importance_method": "permutation (validation AUC drop)",
                        })
    print("Saved scores to", OUTPUT_DIR / OUTPUT_SCORES)


def print_train_validation_metric_diff(history, best_epoch):
    index = best_epoch - 1
    train_auc = history["train_auc"][index]
    validation_auc = history["validation_auc"][index]
    print(
        f"Train-validation AUC difference at best epoch {best_epoch}: "
        f"{train_auc - validation_auc:.6f} (train={train_auc:.6f}, validation={validation_auc:.6f})"
    )


def permutation_importance(model, x_test, y_test, w_test, n_repeats=3):
    """Weighted AUC drop when each feature column is shuffled (on the validation split)."""
    rng = np.random.default_rng(RANDOM_STATE)
    base_auc = roc_auc_score(y_test, predict_scores(model, x_test), sample_weight=w_test)
    drops = np.zeros(x_test.shape[1], dtype=np.float64)
    for column in range(x_test.shape[1]):
        column_drops = []
        for _ in range(n_repeats):
            x_perm = np.array(x_test, copy=True)
            x_perm[:, column] = rng.permutation(x_perm[:, column])
            perm_auc = roc_auc_score(y_test, predict_scores(model, x_perm), sample_weight=w_test)
            column_drops.append(base_auc - perm_auc)
        drops[column] = np.mean(column_drops)
    return np.clip(drops, 0.0, None)


# =============================================================================
# Plot Helpers
# =============================================================================

def save_training_history(history, best_epoch):
    epochs = np.asarray(history["epoch"])
    best_index = best_epoch - 1
    best_gap = history["train_auc"][best_index] - history["validation_auc"][best_index]

    fig, (ax, ax_loss) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    ax.plot(epochs, history["train_auc"], label="Train AUC", lw=2)
    ax.plot(epochs, history["validation_auc"], label="Validation AUC", lw=2)
    ax.axvline(best_epoch, color="gray", ls="--", lw=1)
    ax.set_ylabel("AUC")
    ax.set_ylim(0.75, 1.0)
    ax.set_title("Training History")
    ax.grid(alpha=0.3)
    ax.legend(title=f"AUC gap at best epoch {best_epoch}: {best_gap:.4f}")

    ax_loss.plot(epochs, history["train_loss"], label="Train loss", lw=2)
    ax_loss.plot(epochs, history["validation_loss"], label="Validation loss", lw=2)
    ax_loss.axvline(best_epoch, color="gray", ls="--", lw=1)
    ax_loss.set_xlabel("Epoch")
    ax_loss.set_ylabel("Weighted BCE")
    ax_loss.grid(alpha=0.3)
    ax_loss.legend()
    fig.tight_layout()
    ml_plots.label_sample(fig, CURRENT_CASE["label"])
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
    ml_plots.label_sample(fig, CURRENT_CASE["label"])
    fig.savefig(OUTPUT_DIR / "roc_curve.pdf")
    plt.close(fig)


def save_feature_importance(importance, feature_names):
    total = np.sum(importance)
    percent = 100.0 * importance / total if total > 0 else np.zeros_like(importance)
    order = np.argsort(percent)[::-1]
    ordered_features = np.array(feature_names)[order]
    ordered_percent = percent[order]
    cumulative_percent = np.cumsum(ordered_percent)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(ordered_features, cumulative_percent, color="tab:green", alpha=0.35, label="Cumulative")
    ax.scatter(ordered_percent, ordered_features, color="black", zorder=3, label="Individual")
    for feature, individual, cumulative in zip(ordered_features, ordered_percent, cumulative_percent):
        ax.text(cumulative + 1.0, feature, f"{cumulative:.1f}%", va="center", fontsize=8)
        ax.text(individual + 1.0, feature, f"{individual:.1f}%", va="center", fontsize=8, color="black")
    ax.invert_yaxis()
    ax.set_xlabel("Permutation importance (validation AUC drop) [%]")
    ax.set_title("Cumulative Feature Importance (permutation)")
    ax.set_xlim(0.0, 108.0)
    ax.grid(axis="x", alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    ml_plots.label_sample(fig, CURRENT_CASE["label"])
    fig.savefig(OUTPUT_DIR / "feature_importance.pdf")
    plt.close(fig)


SPECTATOR_STYLES = (("gold", "--"), ("tab:purple", ":"), ("tab:green", "-."))


def save_score_distributions(score_train, y_train, w_train, score_test, y_test, w_test, spectators):
    """Normalized signal and background scores, train and test, with the spectator MC samples of the case,
    and the test-train difference."""
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
    hists, ks = {}, {}
    for name, label in (("signal", 1), ("background", 0)):
        in_train, in_test = y_train == label, y_test == label
        for split, scores, weights, mask in (("train", score_train, w_train, in_train), ("test", score_test, w_test, in_test)):
            hists[name, split] = np.histogram(scores[mask], bins=bins, weights=weights[mask] / np.sum(weights[mask]))[0]
        ks[name] = ml_fom.ks_statistic(score_train[in_train], score_test[in_test], w_train[in_train], w_test[in_test])

    for name, color in (("signal", "tab:orange"), ("background", "tab:blue")):
        ax.stairs(hists[name, "train"], bins, lw=2, color=color, label=f"{name.capitalize()} train")
        ax.stairs(hists[name, "test"], bins, fill=True, alpha=0.25, color=color, label=f"{name.capitalize()} test")
        train = hists[name, "train"]
        difference = np.where(train > 0, 100.0 * (hists[name, "test"] - train) / np.where(train > 0, train, 1.0), np.nan)
        ax_diff.plot(bin_centers, difference, marker="o", ms=3, lw=1.5, color=color, label=name.capitalize())
    for (key, part), (color, style) in zip(spectators.items(), SPECTATOR_STYLES):
        ax.hist(part["score"], bins=bins, weights=part["w"] / np.sum(part["w"]), histtype="step", lw=2, linestyle=style,
                color=color, label=CURRENT_CASE["inputs"][key]["label"])
    ax_diff.axhline(0.0, color="gray", lw=1)

    ax.set_ylabel("Normalized entries")
    ax.set_title("Score Distributions")
    ax.grid(alpha=0.3)
    distribution_legend = ax.legend(fontsize=9, loc="upper left")
    ax.add_artist(distribution_legend)
    ks_handles = [
        Line2D([], [], linestyle="none", label=f"Signal: {ks['signal']:.4f}"),
        Line2D([], [], linestyle="none", label=f"Background: {ks['background']:.4f}"),
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
    ax_diff.set_xlabel("NN score")
    ax_diff.set_ylabel("(test-train)/train [%]")
    ax_diff.grid(alpha=0.3)
    ax_diff.legend(fontsize=8, loc="best")
    ml_plots.label_sample(fig, CURRENT_CASE["label"])
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
    ml_plots.label_sample(fig, CURRENT_CASE["label"])
    fig.savefig(OUTPUT_DIR / "confusion_matrix.pdf")
    plt.close(fig)


# =============================================================================
# Persistent Training Report
# =============================================================================

def save_training_report(
    model,
    history,
    best_epoch,
    importance,
    score_train,
    y_train,
    w_train,
    score_test,
    y_test,
    w_test,
    spectators,
    validation,
    test,
    window_bkg,
    pre_cut,
    params,
    max_signal,
    max_background,
    max_epochs,
    early_stopping,
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

    importance_sum = np.sum(importance)
    importance_fraction = importance / importance_sum if importance_sum > 0.0 else np.zeros_like(importance)

    bins = np.linspace(0.0, 1.0, 26)
    score_histograms = {"bin_low": bins[:-1], "bin_high": bins[1:]}
    for name, scores, weights in (
        ("signal_train", score_train[signal_train], w_train[signal_train]),
        ("background_train", score_train[background_train], w_train[background_train]),
        ("signal_test", score_test[signal_test], w_test[signal_test]),
        ("background_test", score_test[background_test], w_test[background_test]),
    ) + tuple((f"spectator_{key}", part["score"], part["w"]) for key, part in spectators.items()):
        weights = np.asarray(weights, dtype=np.float64)
        score_histograms[name] = np.histogram(scores, bins=bins, weights=weights / np.sum(weights))[0]

    n_parameters = int(sum(p.numel() for p in model.parameters()))
    metadata = {
        "format_version": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "python_version": sys.version,
        "package_versions": {
            "numpy": np.__version__,
            "uproot": uproot.__version__,
            "torch": torch.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "device": str(DEVICE),
        "sample": CURRENT_CASE["name"],
        "fom": CURRENT_FOM,
        "fom_settings": ml_fom.fom_settings(CURRENT_CASE),
        "window_background": window_bkg,
        "skims": {key: part["skim"] for key, part in CURRENT_CASE["inputs"].items()},
        "trees": {key: part["tree"] for key, part in CURRENT_CASE["inputs"].items()},
        "features": list(CURRENT_FEATURES),
        "log_features": list(settings.LOG_FEATURES),
        "signal_cut": pre_cut,
        "background_sideband_cut": CURRENT_CASE["training_sideband_cut"],
        "mc_weight_branch": settings.MC_WEIGHT_BRANCH,
        "class_weighting": "pThat-weighted signal and unit-weight data, normalized to equal class sums",
        "random_state": RANDOM_STATE,
        "max_signal": max_signal,
        "max_background": max_background,
        "max_epochs": max_epochs,
        "early_stopping": early_stopping,
        "nn_params": params,
        "n_parameters": n_parameters,
        "architecture": str(model),
        "model_file": OUTPUT_MODEL,
        "report_file": OUTPUT_REPORT,
        "sample_statistics": sample_stats,
    }
    # validation_*: model selection (train vs validation); test_*: final check (train vs test)
    # at the validation working point.
    metrics = {f"validation_{name}": float(value) for name, value in validation.items()}
    metrics.update({f"test_{name}": float(value) for name, value in test.items()})
    metrics.update({
        "best_epoch": float(best_epoch),
        "train_entries": float(len(y_train)),
        "test_entries": float(len(y_test)),
        "train_weight_sum": float(np.sum(w_train)),
        "test_weight_sum": float(np.sum(w_test)),
        "n_parameters": float(n_parameters),
    })

    report_path = OUTPUT_DIR / OUTPUT_REPORT
    with uproot.recreate(report_path) as report:
        report["metadata_json"] = json.dumps(metadata, indent=2, sort_keys=True)
        report["metrics"] = {
            name: np.asarray([value], dtype=np.float64)
            for name, value in metrics.items()
        }
        report["training_history"] = {
            "epoch": np.asarray(history["epoch"], dtype=np.int32),
            "train_auc": np.asarray(history["train_auc"], dtype=np.float64),
            "validation_auc": np.asarray(history["validation_auc"], dtype=np.float64),
            "train_loss": np.asarray(history["train_loss"], dtype=np.float64),
            "validation_loss": np.asarray(history["validation_loss"], dtype=np.float64),
            "learning_rate": np.asarray(history["learning_rate"], dtype=np.float64),
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
            "auc_drop": np.asarray(importance, dtype=np.float64),
            "importance_fraction": np.asarray(importance_fraction, dtype=np.float64),
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
    main(args.case, args.fom, args.features.split(",") if args.features else None, args.params_json,
         args.output_dir, args.max_signal, args.max_background)
