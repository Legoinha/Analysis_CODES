#!/usr/bin/env python3
"""Correlation matrices of the classifier variables of each case "<tree>/<year>" (its feature-study pool
and production features, ML_selection/settings.py) and Bmass, in the training populations of the case,
read from its skims (after the pre-cut): DATA in the training sidebands (background), and the signal MC
weighted with pThatreweight (signal). Informative only.

Usage, in the LCG view:  python3 correlation_dataMC.py                      every case
                         python3 correlation_dataMC.py ntKp/PbPb23 ntmix/PbPb18
Output: <tree>/<year>/corr_<tree>_<year>_data_sidebands.pdf and corr_<tree>_<year>_mc_signal.pdf
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ML_selection" / "methods" / "ML_common"))
import ml_plots
import ml_samples


def draw(arrays, weights, variables, title, path, label):
    """The weighted correlation matrix of `variables`, with its values written in the cells."""
    x = np.column_stack([np.asarray(arrays[name], np.float64) for name in variables])
    centered = x - np.average(x, axis=0, weights=weights)
    covariance = (centered * weights[:, None]).T @ centered / weights.sum()
    sigma = np.sqrt(np.diag(covariance))
    matrix = covariance / np.outer(sigma, sigma)
    fig, ax = plt.subplots(figsize=(0.6 * len(variables) + 3, 0.55 * len(variables) + 3))
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-1.0, vmax=1.0)
    ax.set_xticks(range(len(variables)))
    ax.set_yticks(range(len(variables)))
    ax.set_xticklabels(variables, rotation=45, ha="right")
    ax.set_yticklabels(variables)
    ax.set_title(title, fontsize=14, fontweight="bold")
    for i in range(len(variables)):
        for j in range(len(variables)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    ml_plots.label_sample(fig, label)
    fig.savefig(path)
    plt.close(fig)
    print("Saved", path)


def main(name):
    c = ml_samples.case(name)
    variables = c["universe"] + ["Bmass"]
    out = HERE / c["tree"] / c["year"]
    out.mkdir(parents=True, exist_ok=True)
    stem = f"corr_{c['tree']}_{c['year']}"
    data = ml_samples.read_branches(c["inputs"]["data"]["skim"], c["inputs"]["data"]["tree"], variables)[1]
    bmass = np.asarray(data["Bmass"], np.float64)
    sideband = np.any([(bmass > low) & (bmass < high) for low, high in c["training_sidebands"]], axis=0)
    data = {key: values[sideband] for key, values in data.items()}
    draw(data, np.ones(int(sideband.sum())), variables, f"{c['label']} DATA training sidebands ({int(sideband.sum())} candidates)",
         out / f"{stem}_data_sidebands.pdf", c["label"])
    mc = ml_samples.read_branches(c["inputs"]["signal"]["skim"], c["inputs"]["signal"]["tree"],
                                  variables + [ml_samples.MC_WEIGHT_BRANCH])[1]
    draw(mc, np.asarray(mc[ml_samples.MC_WEIGHT_BRANCH], np.float64), variables,
         f"{c['label']} {c['inputs']['signal']['label']}, pThat-weighted ({len(mc['Bmass'])} candidates)",
         out / f"{stem}_mc_signal.pdf", c["label"])


if __name__ == "__main__":
    for name in sys.argv[1:] or ml_samples.CASES:
        main(name)
