"""Figures of the ML selection checks.

- decorrelation: made by make_skims.py per sample. The background percentiles of every feature
  versus Bmass before and after the mass decorrelation (flat after it), and the signal MC
  distribution before and after (almost unchanged).
- sculpting: made by feature_study.py for every trained candidate, from its scores file. The mass
  spectrum of the untrained DATA before and after the cut, and the background efficiency versus
  Bmass with the line of ml_fom.sculpting, with the mass regions of the analysis shaded.
- punzi_mass: also for every trained candidate. The whole DATA mass spectrum after the
  Punzi-optimal cut, with the counts in the Punzi window and sidebands.
- punzi: also for every trained candidate. The Punzi FOM and the signal and Punzi-sideband
  background efficiencies versus the threshold, on all candidates, validation and training.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import ml_decorrelation
import ml_fom
import ml_samples

M_X = ml_fom.X3872_MASS
THRESHOLD_COLORS = {"looser": "tab:blue", "working_point": "tab:red", "tighter": "tab:green"}


def shade_regions(ax, labels):
    """Signal window, Punzi sidebands, gaps and training sidebands around the X(3872)."""
    w, s0, s1 = ml_fom.SIGNAL_HALF_WIDTH, ml_fom.SIDEBAND_START, ml_fom.SIDEBAND_START + ml_fom.SIDEBAND_WIDTH
    low, high = ml_fom.SCULPT_FIT_EXCLUDE
    spans = [((M_X - w, M_X + w), "tab:red", "signal window"),
             ((M_X - s1, M_X - s0), "tab:orange", "Punzi sidebands"), ((M_X + s0, M_X + s1), "tab:orange", None),
             ((low, M_X - ml_fom.SCULPT_GAP[0]), "tab:green", "gaps"), ((M_X + ml_fom.SCULPT_GAP[0], high), "tab:green", None)]
    spans += [(band, "tab:blue", "training sidebands" if i == 0 else None)
              for i, band in enumerate(ml_samples.TRAINING_SIDEBANDS)]
    for (a, b), color, label in spans:
        ax.axvspan(a, b, color=color, alpha=0.12, lw=0, label=label if labels else None)


def sculpting(bmass, score, result, path, title):
    """result: ml_fom.sculpting of these candidates, one entry per threshold."""
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 9), gridspec_kw={"height_ratios": [2, 1.4]})
    fine = np.arange(ml_fom.SCULPT_FIT_RANGE[0], ml_fom.SCULPT_FIT_RANGE[1] + 1e-9, 0.0025)
    coarse = np.arange(ml_fom.SCULPT_FIT_RANGE[0], ml_fom.SCULPT_FIT_RANGE[1] + 1e-9, ml_fom.SCULPT_BIN_WIDTH)
    centers = 0.5 * (coarse[:-1] + coarse[1:])
    in_fit = (centers < ml_fom.SCULPT_FIT_EXCLUDE[0]) | (centers > ml_fom.SCULPT_FIT_EXCLUDE[1])
    shade_regions(top, True)
    shade_regions(bottom, False)
    top.stairs(np.histogram(bmass, fine)[0], fine, color="black", lw=1.5, label="before the cut")
    total = np.histogram(bmass, coarse)[0]
    used = total >= 20   # as the fit of ml_fom.sculpting
    for i, (name, entry) in enumerate(result.items()):
        color, mean = THRESHOLD_COLORS[name], entry["background_efficiency"]
        passed = score >= entry["threshold"]
        if name == "working_point":
            top.stairs(np.histogram(bmass[passed], fine)[0] / mean, fine, color=color, lw=1.2,
                       label=f"after the working-point cut, scaled by 1 / background eff. ({mean:.4f})")
        eff = np.histogram(bmass[passed], coarse)[0][used] / total[used]
        err = np.sqrt(np.maximum(eff * (1 - eff), 1.0 / total[used]) / total[used])
        x = centers[used] + (i - 1) * 0.0008
        fit = in_fit[used]
        bottom.errorbar(x[fit], eff[fit] / mean, err[fit] / mean, fmt="o", ms=3, color=color)
        bottom.errorbar(x[~fit], eff[~fit] / mean, err[~fit] / mean, fmt="o", ms=4, mfc="white", color=color,
                        label=f"{name.replace('_', ' ')}: R_gap {entry['R_gap']:.3f} ± {entry['R_gap_error']:.3f}, "
                              f"slope {100 * entry['slope_per_100MeV_relative']:+.1f} ± "
                              f"{100 * entry['slope_per_100MeV_relative_error']:.1f}% / 100 MeV")
        line = np.polyval(entry["line"], coarse - M_X) / mean
        bottom.plot(coarse, line, color=color, lw=1)
    top.set_ylabel("candidates / 2.5 MeV")
    top.set_title(title + "\n(in the training sidebands only the candidates not used for training)", fontsize=10)
    top.legend(fontsize=7, loc="lower center", ncol=3)
    bottom.axhline(1.0, color="grey", lw=0.8, ls=":")
    bottom.set_ylabel("background eff. / mean")
    bottom.set_xlabel("Bmass [GeV]   (open markers: not in the line fit; the signal window holds X(3872) signal)")
    bottom.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def punzi(scores, window_background, threshold, path, title):
    """scores: a scores file (ml_fom.read_scores); threshold: the working point. The Punzi FOM
    (ml_fom.punzi_curve) on all candidates (which picks the working point), on the validation split
    and on the training split, on one threshold grid; the grey band is the bootstrap +-1 sigma of
    the validation curve, so a training curve inside it differs by statistics only."""
    samples = {"all": ml_fom.full_sample(scores),
               **{name: tuple(scores[name][k] for k in ("y", "score", "w", "bmass")) for name in ("validation", "train")}}
    grid = ml_fom.punzi_curve(*samples["all"], window_background)[2]
    y, score, w, bmass = samples["validation"]
    rng = np.random.default_rng(1)
    band = []
    for _ in range(ml_fom.N_BOOTSTRAP):
        pick = rng.integers(0, len(y), len(y))
        band.append(ml_fom.punzi_curve(y[pick], score[pick], w[pick], bmass[pick], window_background, grid)[3])
    band = np.array(band)
    filled = np.all(np.isfinite(band), axis=0)
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 8), gridspec_kw={"height_ratios": [1.4, 1]})
    top.fill_between(grid[filled], *np.percentile(band[:, filled], [16, 84], axis=0), color="grey", alpha=0.3, lw=0,
                     label="validation ±1σ")
    for (name, sample), style, width in zip(samples.items(), ("-", "--", ":"), (2.0, 1.2, 1.2)):
        label = {"train": "training"}.get(name, name)
        background, signal, _, fom = ml_fom.punzi_curve(*sample, window_background, grid)
        top.plot(grid, fom, ls=style, lw=width, marker=".", ms=3, color="black", label=f"Punzi FOM, {label}")
        bottom.plot(grid, signal, ls=style, lw=width, color="tab:red", label=f"signal eff., {label}")
        bottom.plot(grid, background, ls=style, lw=width, color="tab:blue", label=f"background eff., {label}")
    for ax in (top, bottom):
        ax.axvline(threshold, color="grey", ls="--", lw=1, label=f"working point {threshold:.2f}")
    top.set_ylabel(f"e_S / S_min(B), a = {ml_fom.PUNZI_A:g}, b = {ml_fom.PUNZI_B:g}")
    top.set_title(title, fontsize=10)
    top.legend(fontsize=8)
    bottom.set_yscale("log")
    bottom.set_ylim(1e-4, 1.2)
    bottom.set_ylabel("efficiency")
    bottom.set_xlabel("classifier threshold (score >= threshold)")
    bottom.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def punzi_mass(bmass, score, threshold, path, title):
    """All DATA after the pre-cut, 3.60-4.00 GeV in 10 MeV bins, after the cut (top) and before it
    (bottom). The window background is estimated from the Punzi sidebands as in optimalCUT_X_punzi.C."""
    w, s0, width = ml_fom.SIGNAL_HALF_WIDTH, ml_fom.SIDEBAND_START, ml_fom.SIDEBAND_WIDTH
    edges = np.arange(ml_decorrelation.MASS_RANGE[0], ml_decorrelation.MASS_RANGE[1] + 1e-9, 0.010)
    passed = score >= threshold
    distance = np.abs(bmass - M_X)
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 9), gridspec_kw={"height_ratios": [2, 1.2]})
    for ax, keep, label in ((top, passed, f"after the cut (score >= {threshold:.4f})"), (bottom, np.ones_like(passed), "before the cut")):
        ax.axvspan(M_X - w, M_X + w, color="tab:red", alpha=0.15, lw=0, label="Punzi window")
        for side in (-1, 1):
            ax.axvspan(M_X + side * s0, M_X + side * (s0 + width), color="tab:orange", alpha=0.15, lw=0,
                       label="Punzi sidebands" if side < 0 else None)
        ax.axvline(M_X, color="tab:red", ls="--", lw=0.8, label="X(3872)")
        ax.axvline(ml_decorrelation.PSI2S_MASS, color="tab:purple", ls="--", lw=0.8, label="$\\psi$(2S)")
        ax.stairs(np.histogram(bmass[keep], edges)[0], edges, color="black", lw=1.2, label=label)
        n_window = int(np.sum(keep & (distance < w)))
        n_side = int(np.sum(keep & (distance > s0) & (distance < s0 + width)))
        background = n_side * w / width
        ax.text(0.02, 0.03, f"window: {n_window}, background from sidebands: {background:.1f}\n"
                f"excess: {n_window - background:.1f}, excess / sqrt(background): {(n_window - background) / np.sqrt(background):.2f}",
                transform=ax.transAxes, va="bottom", fontsize=9, bbox={"facecolor": "white", "alpha": 0.8})
        ax.set_ylabel("candidates / 10 MeV")
        ax.legend(fontsize=7, loc="lower right")
    top.set_title(title + "\n(the training-sideband candidates carry their training scores)", fontsize=10)
    bottom.set_xlabel("Bmass [GeV]")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def decorrelation(data, signal, features, numbers, path, title):
    """data: DATA skim arrays (Bmass, f, f_dc); signal: prompt X(3872) MC arrays with weights;
    numbers: {feature: {"corr_raw", "corr_dc", "closure": {region: value}}}."""
    suffix = ml_samples.DECORRELATED_SUFFIX
    edges = np.arange(ml_decorrelation.MASS_RANGE[0], ml_decorrelation.MASS_RANGE[1] + 1e-9, ml_decorrelation.SLICE_WIDTH)
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
        pdf.savefig(fig)
        plt.close(fig)
        for f in features:
            fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5))
            for low, high in ml_decorrelation.EXCLUDED:
                left.axvspan(low, high, color="grey", alpha=0.2, lw=0)
            for level, style in zip(levels, (":", "-", "--")):
                left.plot(centers, [np.quantile(data[f][s], level) for s in slices], color="grey", ls=style,
                          label=f"raw {100 * level:.0f}%")
                left.plot(centers, [np.quantile(data[f + suffix][s], level) for s in slices], color="tab:blue", ls=style,
                          label=f"decorrelated {100 * level:.0f}%")
            closure = ", ".join(f"{k} {v:.3f}" for k, v in numbers[f]["closure"].items())
            left.set_title(f"{f}: DATA percentiles vs Bmass (grey bands: not in the fit)\n"
                           f"corr {numbers[f]['corr_raw']:+.3f} → {numbers[f]['corr_dc']:+.3f}; closure {closure}", fontsize=9)
            left.set_xlabel("Bmass [GeV]")
            left.set_ylabel(f)
            left.legend(fontsize=7, ncol=3)
            low, high = np.quantile(signal[f], [0.005, 0.995])
            bins = np.linspace(low, high, 61)
            w = signal[ml_samples.MC_WEIGHT_BRANCH]
            right.hist(signal[f], bins, weights=w, density=True, histtype="step", color="grey", label="raw")
            right.hist(signal[f + suffix], bins, weights=w, density=True, histtype="step", color="tab:blue", label="decorrelated")
            right.set_title(f"{f}: prompt X(3872) MC, pThat-weighted", fontsize=9)
            right.set_xlabel(f)
            right.legend(fontsize=8)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def trimming(comparison, scores_all, scores_trimmed, path, title):
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
        ax.axvline(rule["minimum"], color="tab:orange", ls="--", label=f"{100 * rule['minimum']:g}% individual")
        ax.axvline(rule["cumulative"], color="tab:purple", ls="--", label=f"{100 * rule['cumulative']:g}% cumulative")
        ax.set_yticks(y)
        ax.set_yticklabels(names)
        ax.invert_yaxis()
        ax.set_xlim(0.0, 1.05)
        ax.set_xlabel(f"importance share in fullPOOL/all ({comparison['all']['importance_method']})")
        ax.set_title(f"{title}: trimming")
        ax.legend(fontsize=8, loc="lower right")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(2, 3, figsize=(14, 8))
        panels = [(f"validation {f}", lambda c, f=f: (c["validation"][f], c["errors"][f])) for f in ml_fom.FOMS]
        panels += [("Punzi FOM at the working point", lambda c: (c["working_point"]["punzi"], 0.0)),
                   ("R_gap at the working point", lambda c: (c["sculpting"]["R_gap"], c["sculpting"]["R_gap_error"]))]
        for ax, (label, value) in zip(axes.flat, panels):
            for k, (name, (entry, _)) in enumerate(runs.items()):
                v, e = value(entry)
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
        pdf.savefig(fig)
        plt.close(fig)

        edges = np.arange(ml_decorrelation.MASS_RANGE[0], ml_decorrelation.MASS_RANGE[1] + 1e-9, 0.010)
        fig, ax = plt.subplots(figsize=(9, 6))
        w, s0, width = ml_fom.SIGNAL_HALF_WIDTH, ml_fom.SIDEBAND_START, ml_fom.SIDEBAND_WIDTH
        ax.axvspan(M_X - w, M_X + w, color="tab:red", alpha=0.12, lw=0)
        for side in (-1, 1):
            ax.axvspan(M_X + side * s0, M_X + side * (s0 + width), color="tab:orange", alpha=0.12, lw=0)
        ax.axvline(ml_decorrelation.PSI2S_MASS, color="tab:purple", ls="--", lw=0.8)
        for name, (entry, scores) in runs.items():
            bmass, score = ml_fom.all_data(scores)
            passed = score >= entry["working_point"]["threshold"]
            distance = np.abs(bmass - M_X)
            n_window = int(np.sum(passed & (distance < w)))
            background = np.sum(passed & (distance > s0) & (distance < s0 + width)) * w / width
            ax.stairs(np.histogram(bmass[passed], edges)[0], edges, color=colors[name], lw=1.5,
                      label=f"{name}: window {n_window}, background {background:.1f}, "
                            f"excess / sqrt(background) {(n_window - background) / np.sqrt(background):.2f}")
        ax.set_xlabel("Bmass [GeV]")
        ax.set_ylabel("candidates / 10 MeV")
        ax.set_title(f"{title}: DATA at each working point", fontsize=10)
        ax.legend(fontsize=8, loc="lower center")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(9, 5))
        for name, (entry, scores) in runs.items():
            sample = ml_fom.full_sample(scores)
            _, _, grid, fom = ml_fom.punzi_curve(*sample, ml_fom.window_background(sample[3][sample[0] == 0]))
            ax.plot(grid, fom, marker=".", ms=3, color=colors[name], label=f"{name}, working point {entry['working_point']['threshold']:.2f}")
            ax.axvline(entry["working_point"]["threshold"], color=colors[name], ls="--", lw=0.8)
        ax.set_xlabel("classifier threshold (score >= threshold)")
        ax.set_ylabel(f"e_S / S_min(B), a = {ml_fom.PUNZI_A:g}, b = {ml_fom.PUNZI_B:g}")
        ax.set_title(f"{title}: Punzi FOM", fontsize=10)
        ax.legend(fontsize=8)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
