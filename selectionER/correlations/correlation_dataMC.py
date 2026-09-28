#!/usr/bin/env python3
"""Correlation matrices of the feature-study variables (FEATURE_POOL of
ML_common/ml_samples.py) and Bmass, in the training populations of each sample, after its pre-cut:
DATA in the training sidebands (background), and prompt X(3872) MC weighted with pThatreweight
(signal). Informative only.

Usage, in the LCG view:  python3 correlation_dataMC.py PbPb18 PbPb23
Output: <sample>/corr_<sample>_data_sidebands.pdf and <sample>/corr_<sample>_mc_X3872.pdf
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ML_common"))
import ml_samples

VARIABLES = ml_samples.FEATURE_POOL + ["Bmass"]


def load(sample, key, extra):
    """VARIABLES and `extra` of one flat file after the sample's pre-cut, in double precision."""
    cut = ml_samples.PRE_CUTS[sample]
    _, arrays = ml_samples.read_branches(ml_samples.FLAT_INPUTS[sample][key], ml_samples.TREES[key],
                                         VARIABLES + extra + ml_samples.branches_from_cut(cut))
    keep = ml_samples.cut_mask(arrays, cut)
    return {name: np.asarray(arrays[name][keep], np.float64) for name in VARIABLES + extra}


def correlation(arrays, weights):
    x = np.column_stack([arrays[name] for name in VARIABLES])
    centered = x - np.average(x, axis=0, weights=weights)
    covariance = (centered * weights[:, None]).T @ centered / weights.sum()
    sigma = np.sqrt(np.diag(covariance))
    return covariance / np.outer(sigma, sigma)


def draw(matrix, title, path):
    fig, ax = plt.subplots(figsize=(10, 9))
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-1.0, vmax=1.0)
    ax.set_xticks(range(len(VARIABLES)))
    ax.set_yticks(range(len(VARIABLES)))
    ax.set_xticklabels(VARIABLES, rotation=45, ha="right")
    ax.set_yticklabels(VARIABLES)
    ax.set_title(title, fontsize=14, fontweight="bold")
    for i in range(len(VARIABLES)):
        for j in range(len(VARIABLES)):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=8)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print("Saved", path)


def main(sample):
    out = HERE / sample
    out.mkdir(exist_ok=True)
    data = load(sample, "data", [])
    sideband = ml_samples.in_training_sidebands(data["Bmass"])
    data = {name: values[sideband] for name, values in data.items()}
    n_data = len(data["Bmass"])
    draw(correlation(data, np.ones(n_data)), f"{sample} DATA training sidebands ({n_data} candidates)",
         out / f"corr_{sample}_data_sidebands.pdf")
    mc = load(sample, "signal", [ml_samples.MC_WEIGHT_BRANCH])
    draw(correlation(mc, mc[ml_samples.MC_WEIGHT_BRANCH]),
         f"{sample} prompt X(3872) MC, pThat-weighted ({len(mc['Bmass'])} candidates)",
         out / f"corr_{sample}_mc_X3872.pdf")


if __name__ == "__main__":
    for sample in sys.argv[1:]:
        main(sample)
