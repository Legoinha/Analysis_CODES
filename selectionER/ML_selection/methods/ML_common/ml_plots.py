"""Figures of the ML selection checks.

- decorrelation: made by make_skims.py per case. The background percentiles of every feature
  versus Bmass before and after the mass decorrelation (flat after it), and the signal MC
  distribution before and after (almost unchanged).
- sculpting: made by feature_study.py for every trained candidate, from its scores file. The mass
  spectrum of the untrained DATA before and after the cut, and the background efficiency versus
  Bmass with the line of ml_fom.sculpting, both in its bins, with the mass regions of the
  analysis shaded.
- punzi_mass: also for every trained candidate. The whole DATA mass spectrum after the
  Punzi-optimal cut (page 1) and after the cuts of MASS_CUT_VARIATIONS around it (5% tighter, 5%,
  10%, 15% and 20% looser): an excess only at the optimal cut is a fluctuation of the few candidates
  left at high threshold.
- punzi: also for every trained candidate. The Punzi FOM and the signal and Punzi-sideband
  background efficiencies versus the threshold, on all candidates, validation and training.
Every figure carries its case in a box at its top left corner (label_sample). c is the case dict of
ml_samples.case: its mass regions (ML_selection/settings.py MASS) set the shading and the ranges.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import ml_fom

THRESHOLD_COLORS = {"looser": "tab:blue", "working_point": "tab:red", "tighter": "tab:green"}
MASS_CUT_VARIATIONS = (0.05, -0.05, -0.10, -0.15, -0.20)   # punzi_mass, pages 2-6: threshold + f |threshold|


def label_sample(fig, label):
    """The case (ntmix PbPb23, ntKp ppRef, ...) in a box at the top left corner of the figure, outside the
    axes, so that no page of any study can be read without it."""
    fig.text(0.004, 0.996, label, ha="left", va="top", fontsize=10, fontweight="bold",
             bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "0.5"})


def shade_regions(ax, labels, c):
    """Signal window, training sidebands, the near-peak background of R_gap (from sculpting_gap off the
    peak to the edges of the line fit) and, hatched as they may overlap the training sidebands, the
    Punzi sidebands."""
    m, w, s0 = c["peak"], c["punzi_half_width"], c["punzi_sideband_start"]
    s1, gap = s0 + c["punzi_sideband_width"], c["sculpting_gap"]
    low, high = c["sculpting_excluded"]
    spans = [((m - w, m + w), "tab:red", "signal window")]
    spans += [((low, m - gap), "tab:green", "near-peak background (R_gap)"), ((m + gap, high), "tab:green", None)]
    spans += [(band, "tab:blue", "training sidebands" if i == 0 else None)
              for i, band in enumerate(c["training_sidebands"])]
    for (a, b), color, label in spans:
        ax.axvspan(a, b, color=color, alpha=0.12, lw=0, label=label if labels else None)
    for i, side in enumerate(c["punzi_sides"]):
        a, b = sorted((m + side * s0, m + side * s1))
        ax.axvspan(a, b, facecolor="none", edgecolor="tab:orange", hatch="//", lw=0,
                   label="Punzi sidebands" if labels and i == 0 else None)


def sculpting(bmass, score, result, path, title, c):
    """result: ml_fom.sculpting of these candidates, one entry per threshold. The spectrum and the
    efficiency share the sculpting_bin bins of the line fit."""
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 9), gridspec_kw={"height_ratios": [2, 1.4]})
    width = c["sculpting_bin"]
    edges = np.arange(c["sculpting_range"][0], c["sculpting_range"][1] + 1e-9, width)
    centers = 0.5 * (edges[:-1] + edges[1:])
    in_fit = (centers < c["sculpting_excluded"][0]) | (centers > c["sculpting_excluded"][1])
    shade_regions(top, True, c)
    shade_regions(bottom, False, c)
    # The spectra count the thinned DATA as often as it was thinned; the efficiencies are within bins.
    weights = ml_fom.data_weights(bmass, c)
    before = np.histogram(bmass, edges, weights=weights)[0]
    top.stairs(before, edges, color="black", lw=1.5, label="before the cut")
    highest = before.max()
    total = np.histogram(bmass, edges)[0]
    used = total >= 20   # as the fit of ml_fom.sculpting
    for i, (name, entry) in enumerate(result.items()):
        color, mean = THRESHOLD_COLORS[name], entry["background_efficiency"]
        passed = score >= entry["threshold"]
        if name == "working_point":
            scaled = np.histogram(bmass[passed], edges, weights=weights[passed])[0] / mean
            top.stairs(scaled, edges, color=color, lw=1.2, label=f"after the working-point cut, scaled by 1 / background eff. ({mean:.4f})")
            highest = max(highest, scaled.max())
        eff = np.histogram(bmass[passed], edges)[0][used] / total[used]
        err = np.sqrt(np.maximum(eff * (1 - eff), 1.0 / total[used]) / total[used])
        x = centers[used] + (i - 1) * 0.16 * width
        fit = in_fit[used]
        bottom.errorbar(x[fit], eff[fit] / mean, err[fit] / mean, fmt="o", ms=3, color=color)
        bottom.errorbar(x[~fit], eff[~fit] / mean, err[~fit] / mean, fmt="o", ms=4, mfc="white", color=color,
                        label=f"{name.replace('_', ' ')}: R_gap {entry['R_gap']:.3f} ± {entry['R_gap_error']:.3f}, "
                              f"slope {100 * entry['slope_per_100MeV_relative']:+.1f} ± "
                              f"{100 * entry['slope_per_100MeV_relative_error']:.1f}% / 100 MeV")
        line = np.polyval(entry["line"], edges - c["peak"]) / mean
        bottom.plot(edges, line, color=color, lw=1)
    top.set_ylabel(f"candidates / {1000 * width:g} MeV")
    top.set_ylim(0, 1.8 * highest)
    top.set_title(title + "\n(in the training sidebands only the candidates not used for training)", fontsize=10)
    top.legend(fontsize=10, loc="upper right")
    bottom.axhline(1.0, color="grey", lw=0.8, ls=":")
    bottom.set_ylabel("background eff. / mean")
    bottom.set_xlabel("Bmass [GeV]   (open markers: not in the line fit; the signal window holds signal)")
    bottom.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    label_sample(fig, c["label"])
    fig.savefig(path)
    plt.close(fig)


def punzi(scores, window_background, threshold, path, title, c):
    """scores: a scores file (ml_fom.read_scores); threshold: the working point. The Punzi FOM
    (ml_fom.punzi_curve) on all candidates (which picks the working point), on the validation split
    and on the training split, on one threshold grid; the grey band is the bootstrap +-1 sigma of
    the validation curve, so a training curve inside it differs by statistics only."""
    samples = {"all": ml_fom.full_sample(scores),
               **{name: tuple(scores[name][k] for k in ("y", "score", "w", "bmass")) for name in ("validation", "train")}}
    grid = ml_fom.punzi_curve(*samples["all"], window_background, c)[2]
    y, score, w, bmass = samples["validation"]
    rng = np.random.default_rng(1)
    band = []
    for _ in range(ml_fom.N_BOOTSTRAP):
        pick = rng.integers(0, len(y), len(y))
        band.append(ml_fom.punzi_curve(y[pick], score[pick], w[pick], bmass[pick], window_background, c, grid)[3])
    band = np.array(band)
    filled = np.all(np.isfinite(band), axis=0)
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 8), gridspec_kw={"height_ratios": [1.4, 1]})
    top.fill_between(grid[filled], *np.percentile(band[:, filled], [16, 84], axis=0), color="grey", alpha=0.3, lw=0,
                     label="validation ±1σ")
    for (name, arrays), style, width in zip(samples.items(), ("-", "--", ":"), (2.0, 1.2, 1.2)):
        label = {"train": "training"}.get(name, name)
        background, signal, _, fom = ml_fom.punzi_curve(*arrays, window_background, c, grid)
        top.plot(grid, fom, ls=style, lw=width, marker=".", ms=3, color="black", label=f"Punzi FOM, {label}")
        bottom.plot(grid, signal, ls=style, lw=width, color="tab:red", label=f"signal eff., {label}")
        bottom.plot(grid, background, ls=style, lw=width, color="tab:blue", label=f"background eff., {label}")
    for ax in (top, bottom):
        ax.axvline(threshold, color="grey", ls="--", lw=1, label=f"working point {threshold:.2f}")
    top.set_ylabel(f"e_S / S_min(B), a = {ml_fom.settings.PUNZI_A:g}, b = {ml_fom.settings.PUNZI_B:g}")
    top.set_title(title, fontsize=10)
    top.legend(fontsize=8)
    bottom.set_yscale("log")
    bottom.set_ylim(1e-4, 1.2)
    bottom.set_ylabel("efficiency")
    bottom.set_xlabel("classifier threshold (score >= threshold)")
    bottom.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    label_sample(fig, c["label"])
    fig.savefig(path)
    plt.close(fig)


def punzi_mass(bmass, score, threshold, path, title, c):
    """All DATA after the pre-cut in the decorrelation range in 10 MeV bins, after the cut (top) and
    before it (bottom): page 1 at the Punzi-optimal threshold, then one page per MASS_CUT_VARIATIONS.
    The window background is estimated from the Punzi sidebands as in optimalCUT_X_punzi.C."""
    pages = [(threshold, title)] + [
        (threshold + f * abs(threshold), f"DATA at a cut {100 * abs(f):g}% {'tighter' if f > 0 else 'looser'} "
                                         f"than the Punzi-optimal {threshold:.4f}") for f in MASS_CUT_VARIATIONS]
    with PdfPages(path) as pdf:
        for cut, label in pages:
            fig = mass_after_cut(bmass, score, cut, label, c)
            label_sample(fig, c["label"])
            pdf.savefig(fig)
            plt.close(fig)


def mass_after_cut(bmass, score, threshold, title, c):
    """One page of punzi_mass."""
    m, w, s0, width, sides = c["peak"], c["punzi_half_width"], c["punzi_sideband_start"], c["punzi_sideband_width"], c["punzi_sides"]
    edges = np.arange(c["decorrelation_range"][0], c["decorrelation_range"][1] + 1e-9, 0.010)
    passed = score >= threshold
    weights = ml_fom.data_weights(bmass, c)   # the thinned DATA counted as often as it was thinned
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 9), gridspec_kw={"height_ratios": [2, 1.2]})
    for ax, keep, label in ((top, passed, f"after the cut (score >= {threshold:.4f})"), (bottom, np.ones_like(passed), "before the cut")):
        ax.axvspan(m - w, m + w, color="tab:red", alpha=0.15, lw=0, label="Punzi window")
        for i, side in enumerate(sides):
            ax.axvspan(m + side * s0, m + side * (s0 + width), color="tab:orange", alpha=0.15, lw=0,
                       label="Punzi sidebands" if i == 0 else None)
        for (name, mass), color in zip(c["marks"].items(), ("tab:red", "tab:purple")):
            ax.axvline(mass, color=color, ls="--", lw=0.8, label=name)
        counts = np.histogram(bmass[keep], edges, weights=weights[keep])[0]
        ax.stairs(counts, edges, color="black", lw=1.2, label=label)
        ax.set_ylim(0, 1.45 * max(counts.max(), 1))
        ax.set_ylabel("candidates / 10 MeV")
        ax.legend(fontsize=10, loc="upper right")
    top.set_title(title, fontsize=10)
    bottom.set_xlabel("Bmass [GeV]")
    fig.tight_layout()
    return fig


def decorrelation(data, signal, features, numbers, path, title, c):
    """data: DATA skim arrays (Bmass, f, f_dc); signal: signal MC arrays with pThatreweight;
    numbers: {feature: {"corr_raw", "corr_dc", "closure": {region: value}}}."""
    edges = np.arange(c["decorrelation_range"][0], c["decorrelation_range"][1] + 1e-9, c["decorrelation_slice"])
    centers = 0.5 * (edges[:-1] + edges[1:])
    slices = [(data["Bmass"] >= a) & (data["Bmass"] < b) for a, b in zip(edges[:-1], edges[1:])]
    levels = (0.1, 0.5, 0.9)
    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(9, 0.4 * len(features) + 2))
        y = np.arange(len(features))
        ax.barh(y + 0.2, [numbers[f]["corr_raw"] for f in features], 0.4, color="grey", label="raw")
        ax.barh(y - 0.2, [numbers[f]["corr_dc"] for f in features], 0.4, color="tab:blue", label="decorrelated (_dc)")
        ax.set_yticks(y)
        ax.set_yticklabels(features)
        ax.invert_yaxis()
        ax.axvline(0.0, color="black", lw=0.8)
        ax.set_xlabel("correlation with Bmass, DATA in the fit region")
        ax.set_title(title, fontsize=11)
        ax.legend()
        fig.tight_layout()
        label_sample(fig, c["label"])
        pdf.savefig(fig)
        plt.close(fig)
        for f in features:
            fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5))
            for low, high in c["decorrelation_excluded"]:
                left.axvspan(low, high, color="grey", alpha=0.2, lw=0)
            for level, style in zip(levels, (":", "-", "--")):
                left.plot(centers, [np.quantile(data[f][s], level) for s in slices], color="grey", ls=style,
                          label=f"raw {100 * level:.0f}%")
                left.plot(centers, [np.quantile(data[f + "_dc"][s], level) for s in slices], color="tab:blue", ls=style,
                          label=f"decorrelated {100 * level:.0f}%")
            closure = ", ".join(f"{k} {v:.3f}" for k, v in numbers[f]["closure"].items())
            left.set_title(f"{f}: DATA percentiles vs Bmass (grey bands: not in the fit)\n"
                           f"corr {numbers[f]['corr_raw']:+.3f} → {numbers[f]['corr_dc']:+.3f}; closure {closure}", fontsize=9)
            left.set_xlabel("Bmass [GeV]")
            left.set_ylabel(f)
            left.legend(fontsize=7, ncol=3)
            low, high = np.quantile(signal[f], [0.005, 0.995])
            bins = np.linspace(low, high, 61)
            w = signal["pThatreweight"]
            right.hist(signal[f], bins, weights=w, density=True, histtype="step", color="grey", label="raw")
            right.hist(signal[f + "_dc"], bins, weights=w, density=True, histtype="step", color="tab:blue", label="decorrelated")
            right.set_title(f"{f}: {c['inputs']['signal']['label']}, pThat-weighted", fontsize=9)
            right.set_xlabel(f)
            right.legend(fontsize=8)
            fig.tight_layout()
            label_sample(fig, c["label"])
            pdf.savefig(fig)
            plt.close(fig)


def trimming(comparison, scores_all, scores_trimmed, path, title, c):
    """fullPOOL/trimmed/comparison.pdf: the importance of fullPOOL/all with the trimming rule, the
    metrics of all and trimmed side by side, their DATA mass spectra at their working points, and
    their Punzi curves (comparison: the content of comparison.json)."""
    rule = comparison["rule"]
    runs = {"all": (comparison["all"], scores_all), "trimmed": (comparison["trimmed"], scores_trimmed)}
    colors = {"all": "tab:blue", "trimmed": "tab:red"}
    with PdfPages(path) as pdf:
        ranked = sorted(comparison["importance_of_all"].items(), key=lambda item: -item[1])
        names, shares = [n for n, _ in ranked], np.array([s for _, s in ranked])
        y = np.arange(len(names))
        fig, ax = plt.subplots(figsize=(9, 0.5 * len(names) + 2.5))
        ax.barh(y, shares, color=["tab:green" if n in comparison["kept"] else "lightgrey" for n in names],
                label="individual (green: kept)")
        ax.plot(np.cumsum(shares), y, "o-", color="black", label="cumulative")
        ax.axvline(rule["cumulative"], color="tab:purple", ls="--", label=f"{100 * rule['cumulative']:g}% cumulative")
        ax.set_yticks(y)
        ax.set_yticklabels(names)
        ax.invert_yaxis()
        ax.set_xlim(0.0, 1.05)
        ax.set_xlabel(f"importance share in fullPOOL/all ({comparison['all']['importance_method']})")
        ax.set_title(f"{title}: trimming")
        ax.legend(fontsize=8, loc="lower right")
        fig.tight_layout()
        label_sample(fig, c["label"])
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        entries = [entry for entry, _ in runs.values()]
        panels = [(f"validation {f}", [(e["validation"][f], e["errors"][f]) for e in entries]) for f in ml_fom.FOMS]
        panels += [("Punzi FOM at the working point", [(e["working_point"]["punzi"], 0.0) for e in entries]),
                   ("R_gap at the working point", [(e["sculpting"]["R_gap"], e["sculpting"]["R_gap_error"]) for e in entries])]
        for ax, (label, values) in zip(axes.flat, panels):
            for k, (name, (v, e)) in enumerate(zip(runs, values)):
                ax.errorbar(k, v, e, fmt="o", color=colors[name], ms=7)
            ax.set_xticks(range(len(runs)))
            ax.set_xticklabels([f"{n} ({len(e['features'])})" for n, (e, _) in runs.items()])
            ax.set_xlim(-0.5, len(runs) - 0.5)
            ax.set_title(label, fontsize=10)
            ax.grid(axis="y", alpha=0.3)
        gain = comparison["trimmed_minus_all"]
        fig.suptitle(f"{title}: trimmed minus all (paired): " + ", ".join(
            f"{f} {gain[f]['gain']:+.4f} ± {gain[f]['sigma']:.4f}" for f in ml_fom.FOMS), fontsize=10)
        fig.tight_layout()
        label_sample(fig, c["label"])
        pdf.savefig(fig)
        plt.close(fig)

        edges = np.arange(c["decorrelation_range"][0], c["decorrelation_range"][1] + 1e-9, 0.010)
        fig, ax = plt.subplots(figsize=(9, 6))
        m, w, s0, width, sides = c["peak"], c["punzi_half_width"], c["punzi_sideband_start"], c["punzi_sideband_width"], c["punzi_sides"]
        ax.axvspan(m - w, m + w, color="tab:red", alpha=0.12, lw=0)
        for side in sides:
            ax.axvspan(m + side * s0, m + side * (s0 + width), color="tab:orange", alpha=0.12, lw=0)
        for mass in c["marks"].values():
            ax.axvline(mass, color="tab:purple", ls="--", lw=0.8)
        highest = 0.0
        for name, (entry, scores) in runs.items():
            bmass, score = ml_fom.all_data(scores)
            weights = ml_fom.data_weights(bmass, c)
            passed = score >= entry["working_point"]["threshold"]
            counts = np.histogram(bmass[passed], edges, weights=weights[passed])[0]
            ax.stairs(counts, edges, color=colors[name], lw=1.5, label=name)
            highest = max(highest, counts.max())
        ax.set_ylim(0, 1.45 * highest)
        ax.set_xlabel("Bmass [GeV]")
        ax.set_ylabel("candidates / 10 MeV")
        ax.set_title(f"{title}: DATA at each working point", fontsize=10)
        ax.legend(fontsize=10, loc="upper right")
        fig.tight_layout()
        label_sample(fig, c["label"])
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(9, 5))
        for name, (entry, scores) in runs.items():
            full = ml_fom.full_sample(scores)
            _, _, grid, fom = ml_fom.punzi_curve(*full, ml_fom.window_background(full[3][full[0] == 0], c), c)
            ax.plot(grid, fom, marker=".", ms=3, color=colors[name], label=f"{name}, working point {entry['working_point']['threshold']:.2f}")
            ax.axvline(entry["working_point"]["threshold"], color=colors[name], ls="--", lw=0.8)
        ax.set_xlabel("classifier threshold (score >= threshold)")
        ax.set_ylabel(f"e_S / S_min(B), a = {ml_fom.settings.PUNZI_A:g}, b = {ml_fom.settings.PUNZI_B:g}")
        ax.set_title(f"{title}: Punzi FOM", fontsize=10)
        ax.legend(fontsize=8)
        fig.tight_layout()
        label_sample(fig, c["label"])
        pdf.savefig(fig)
        plt.close(fig)
