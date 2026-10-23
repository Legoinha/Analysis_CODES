# Step 2 (v2): fits, k-scans with bootstrap errors, pulls, universality, closure, transfer tests, figures.
# Data signal is extracted per bin of the studied probability with a binned mass fit whose signal
# shape is seeded from truth-matched MC (data floats the mean and a width scale).
import json, warnings
import numpy as np
from scipy import stats, optimize, special
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C_RAW, C_COR = "#2a78d6", "#eb6834"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
                     "savefig.dpi": 200})
RNG = np.random.default_rng(12345)
NBOOT = 200
R = {}

PARTS = {
    "ppRef_psi2S": dict(m0=3.6861, ndf=6, bkg="poly", label=r"ppRef $\psi(2S)$", dm=0.003),
    "ppRef_X3872": dict(m0=3.8717, ndf=6, bkg="poly", label="ppRef X(3872)", dm=0.003),
    "ppRef_Bp": dict(m0=5.2793, ndf=4, bkg="bplus", label=r"ppRef $B^{+}$", dm=0.008),
    "PbPb23_Bp": dict(m0=5.2793, ndf=4, bkg="bplus", label=r"PbPb 2023 $B^{+}$", dm=0.008),
    "PbPb24_Bp": dict(m0=5.2793, ndf=4, bkg="bplus", label=r"PbPb 2024 data vs 2023 MC $B^{+}$", dm=0.008),
}
PMIN = {"Bchi2Prob": 0.005, "BujvProb": 0.01, "BtktkvProb": 0.0}
NDF_SUB = {"BujvProb": 1, "BtktkvProb": 1}


def load(name):
    z = np.load(f"cache/{name}.npz")
    return ({k[2:]: z[k] for k in z.files if k.startswith("d_")}, {k[2:]: z[k] for k in z.files if k.startswith("m_")})


def edges_for(var):
    return np.array([PMIN[var]] + list(np.arange(0.05, 1.0001, 0.05)))


# ------------------------------------------------------------------ mass model
def dgauss(m, mu, s1, r, f):
    s2 = s1 * r
    return (f * np.exp(-0.5 * ((m - mu) / s1) ** 2) / (s1 * np.sqrt(2 * np.pi))
            + (1 - f) * np.exp(-0.5 * ((m - mu) / s2) ** 2) / (s2 * np.sqrt(2 * np.pi)))


class MassFit:
    def __init__(self, part, mass, win, nb, mc_mass, mc_w):
        self.p = PARTS[part]; self.m0 = self.p["m0"]
        self.edges = np.linspace(win[0], win[1], nb + 1)
        self.c = 0.5 * (self.edges[1:] + self.edges[:-1]); self.bw = self.edges[1] - self.edges[0]
        self.e0, self.ew = 5.14, 0.03
        self.mc_shape(mc_mass, mc_w)
        self.inclusive(mass)

    def mc_shape(self, mm, w):
        e = np.linspace(self.m0 - 0.12 if self.p["bkg"] != "poly" else self.m0 - 0.05,
                        self.m0 + 0.12 if self.p["bkg"] != "poly" else self.m0 + 0.05, 101)
        h = np.histogram(mm, e, weights=w)[0]; h2 = np.histogram(mm, e, weights=w ** 2)[0]
        c = 0.5 * (e[1:] + e[:-1]); bw = e[1] - e[0]
        f = lambda m, n, mu, s1, r, fr: n * bw * dgauss(m, mu, s1, r, fr)
        p0 = [h.sum(), self.m0, 0.005 if self.p["bkg"] == "poly" else 0.02, 2.0, 0.7]
        popt, _ = optimize.curve_fit(f, c, h, p0=p0, sigma=np.sqrt(np.maximum(h2, h2[h2 > 0].min())),
                                     bounds=([0, self.m0 - 0.01, 0.001, 1.05, 0.05], [np.inf, self.m0 + 0.01, 0.08, 6, 0.99]))
        self.mc_mu, self.mc_s1, self.r, self.fr = popt[1:]

    def bkg(self, m, b0, b1, b2):
        x = m - self.m0
        if self.p["bkg"] == "poly":
            return b0 + b1 * x + b2 * x * x
        return b0 * np.exp(b1 * x) + b2 * 0.5 * special.erfc((m - self.e0) / self.ew)

    def _fit(self, h, f, p0, lo, hi):
        sig = np.sqrt(np.maximum(h, 1.0))
        popt, pcov = optimize.curve_fit(f, self.c, h, p0=p0, sigma=sig, bounds=(lo, hi), maxfev=50000)
        for _ in range(2):
            sig = np.sqrt(np.maximum(f(self.c, *popt), 1.0))
            popt, pcov = optimize.curve_fit(f, self.c, h, p0=popt, sigma=sig, bounds=(lo, hi), maxfev=50000)
        return popt, pcov, np.sum((h - f(self.c, *popt)) ** 2 / np.maximum(f(self.c, *popt), 1.0))

    def inclusive(self, mass):
        h = np.histogram(mass, self.edges)[0].astype(float); self.h_incl = h; N = h.sum(); dm = self.p["dm"]
        r, fr = self.r, self.fr
        if self.p["bkg"] == "poly":
            def f(m, ns, mu, sc, b0, b1, b2):
                return ns * self.bw * dgauss(m, mu, self.mc_s1 * sc, r, fr) + self.bkg(m, b0, b1, b2)
            p0 = [0.1 * N, self.mc_mu, 1.1, 0.9 * N / len(h), 0, 0]
            lo = [0, self.mc_mu - dm, 0.7, -np.inf, -np.inf, -np.inf]; hi = [N, self.mc_mu + dm, 1.6, np.inf, np.inf, np.inf]
        else:
            def f(m, ns, mu, sc, b0, b1, b2, e0, ew):
                self.e0, self.ew = e0, ew
                return ns * self.bw * dgauss(m, mu, self.mc_s1 * sc, r, fr) + self.bkg(m, b0, b1, b2)
            p0 = [0.3 * N, self.mc_mu, 1.1, 0.4 * N / len(h), -3.0, 0.3 * N / len(h), 5.14, 0.03]
            lo = [0, self.mc_mu - dm, 0.7, 0, -30, 0, 5.08, 0.005]; hi = [N, self.mc_mu + dm, 1.6, np.inf, 5, np.inf, 5.20, 0.08]
        popt, pcov, chi2 = self._fit(h, f, p0, lo, hi)
        self.incl, self.incl_err = popt, np.sqrt(np.diag(pcov)); self.incl_chi2 = chi2 / (len(h) - len(popt))
        self.mu, self.sc = popt[1], popt[2]; self.s1 = self.mc_s1 * self.sc; self.bpar = popt[3:6]
        if self.p["bkg"] != "poly":
            self.e0, self.ew = popt[6], popt[7]
        self.f_incl = f

    def sigma_eff(self, s1=None):
        s1 = self.s1 if s1 is None else s1
        return np.sqrt(self.fr * s1 ** 2 + (1 - self.fr) * (s1 * self.r) ** 2)

    def yield_in(self, mass, free_shape=True):
        h = np.histogram(mass, self.edges)[0].astype(float); N = h.sum()
        if N < 5:
            return np.nan, np.nan
        frac = N / self.h_incl.sum(); r, fr = self.r, self.fr
        b = self.bpar * (np.array([frac, 1, frac]) if self.p["bkg"] != "poly" else frac)
        blo = [-np.inf] * 3 if self.p["bkg"] == "poly" else [0, -30, 0]
        bhi = [np.inf] * 3 if self.p["bkg"] == "poly" else [np.inf, 5, np.inf]
        if free_shape:
            def f(m, ns, dmu, sc, b0, b1, b2):
                return ns * self.bw * dgauss(m, self.mu + dmu, self.s1 * sc, r, fr) + self.bkg(m, b0, b1, b2)
            p0 = [max(self.incl[0] * frac, 1), 0.0, 1.0, *b]
            lo = [0, -self.p["dm"], 0.75, *blo]; hi = [N, self.p["dm"], 1.35, *bhi]
        else:
            def f(m, ns, b0, b1, b2):
                return ns * self.bw * dgauss(m, self.mu, self.s1, r, fr) + self.bkg(m, b0, b1, b2)
            p0 = [max(self.incl[0] * frac, 1), *b]; lo = [0, *blo]; hi = [N, *bhi]
        try:
            popt, pcov, _ = self._fit(h, f, p0, lo, hi)
            return popt[0], max(np.sqrt(pcov[0, 0]), 1.0)
        except Exception:
            return np.nan, np.nan


def binned_yields(fit, mass, x, edges, free_shape=True):
    y, e = [], []
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        s = (x >= a) & ((x < b) if i < len(edges) - 2 else (x <= b))
        yy, ee = fit.yield_in(mass[s], free_shape); y.append(yy); e.append(ee)
    return np.array(y), np.array(e)


# ------------------------------------------------------------------ k machinery
class MCdist:
    """Weighted MC probabilities. Scaling the MC chi2 by k^2 is monotonic, so histogramming the
    corrected probability in fixed bins equals histogramming the raw one in transformed edges."""

    def __init__(self, p, w):
        o = np.argsort(p); self.p, self.w = p[o], w[o]
        self.cw, self.cw2 = np.r_[0, np.cumsum(self.w)], np.r_[0, np.cumsum(self.w ** 2)]

    def hist(self, edges, k, ndf):
        e = stats.chi2.sf(stats.chi2.isf(edges, ndf) / (k * k), ndf)
        i = np.searchsorted(self.p, e, side="left"); i[-1] = len(self.p)
        return np.diff(self.cw[i]), np.diff(self.cw2[i])

    def resampled(self):
        return MCdist(self.p, self.w * RNG.poisson(1.0, len(self.p)))


def chi2_dm(yd, ed, hm, hm2):
    ok = np.isfinite(yd) & np.isfinite(ed) & (ed > 0)
    norm = yd[ok].sum() / hm[ok].sum()
    return np.sum((yd[ok] - hm[ok] * norm) ** 2 / (ed[ok] ** 2 + hm2[ok] * norm ** 2)), int(ok.sum()) - 1


KS = np.arange(0.90, 1.3001, 0.001)


def best_k(yd, ed, edges, mcd, ndf):
    c = np.array([chi2_dm(yd, ed, *mcd.hist(edges, k, ndf))[0] for k in KS])
    i = int(np.argmin(c)); kb = KS[i]
    lo, hi = max(i - 2, 0), min(i + 3, len(KS))
    if hi - lo == 5:  # local parabola through the 5 grid points around the minimum, clamped to that range
        a, b, _ = np.polyfit(KS[lo:hi], c[lo:hi], 2)
        if a > 0:
            kb = float(np.clip(-b / (2 * a), KS[lo], KS[hi - 1]))
    return float(kb), c


def kfit(yd, ed, edges, p, w, ndf, nboot=NBOOT):
    mcd = MCdist(p, w)
    kb, curve = best_k(yd, ed, edges, mcd, ndf)
    c1, nd = chi2_dm(yd, ed, *mcd.hist(edges, 1.0, ndf))
    cmin, _ = chi2_dm(yd, ed, *mcd.hist(edges, kb, ndf))
    boots = []
    for _ in range(nboot):
        ydb = np.where(np.isfinite(yd), yd + ed * RNG.standard_normal(len(yd)), np.nan)
        boots.append(best_k(ydb, ed, edges, mcd.resampled(), ndf)[0])
    q16, q50, q84 = np.percentile(boots, [16, 50, 84])
    return dict(k=kb, err=float(0.5 * (q84 - q16)), k_boot_median=float(q50), chi2_k1=float(c1), chi2min=float(cmin), ndof=int(nd - 1),
                curve=(curve - curve.min()).tolist())


# ------------------------------------------------------------------ plotting
def band(ax, edges, h, h2, scale, color):
    lo, hi = (h - np.sqrt(h2)) * scale, (h + np.sqrt(h2)) * scale
    ax.fill_between(np.repeat(edges, 2)[1:-1], np.repeat(lo, 2), np.repeat(hi, 2), color=color, alpha=0.18, lw=0)


def draw_compare(ax, axr, edges, yd, ed, mraw, mcor, title, xlabel, lab_cor):
    c = 0.5 * (edges[1:] + edges[:-1]); wbin = np.diff(edges); ok = np.isfinite(yd)
    nd = np.nansum(yd); yn, en = yd / nd / wbin, ed / nd / wbin
    hr, hr2 = mraw; sr = 1.0 / hr[ok].sum() / wbin
    band(ax, edges, hr, hr2, sr, C_RAW); ax.stairs(hr * sr, edges, color=C_RAW, lw=1.6, label="MC raw (±stat)")
    if mcor is not None:
        hc, hc2 = mcor; scr = 1.0 / hc[ok].sum() / wbin
        band(ax, edges, hc, hc2, scr, C_COR); ax.stairs(hc * scr, edges, color=C_COR, lw=1.6, ls="--", label=lab_cor)
    ax.errorbar(c, yn, en, fmt="o", ms=3.5, color=INK, lw=1, label="Data (fitted signal)")
    ax.set_title(title, fontsize=8.5, loc="left", color=INK); ax.set_ylabel("normalized density"); ax.legend(fontsize=7)
    axr.axhline(1, color=INK2, lw=0.8)
    axr.errorbar(c, yn / (hr * sr), en / (hr * sr), fmt="o", ms=3, color=C_RAW, lw=1, label="data / MC raw")
    if mcor is not None:
        axr.errorbar(c + 0.008, yn / (hc * scr), en / (hc * scr), fmt="s", ms=3, color=C_COR, lw=1, mfc="white", label="data / MC corr.")
    rat = np.r_[yn / (hr * sr)] if mcor is None else np.r_[yn / (hr * sr), yn / (hc * scr)]
    rat = rat[np.isfinite(rat)]
    top = min(2.2, max(1.5, np.max(rat) + 0.15)) if len(rat) else 1.5
    axr.set_ylim(min(0.5, 2 - top), top); axr.set_xlabel(xlabel); axr.set_ylabel("data / MC"); axr.legend(fontsize=6.5, ncol=2, loc="lower center")


def fig_compare(items, fname):
    n = len(items)
    fig = plt.figure(figsize=(3.7 * n, 4.5))
    gs = fig.add_gridspec(2, n, height_ratios=[3, 1.3], hspace=0.08, wspace=0.32)
    for j, it in enumerate(items):
        ax = fig.add_subplot(gs[0, j]); axr = fig.add_subplot(gs[1, j], sharex=ax)
        plt.setp(ax.get_xticklabels(), visible=False); draw_compare(ax, axr, *it)
    fig.savefig(f"fig/{fname}.png", bbox_inches="tight"); plt.close(fig)


def mass_figure(fit, part, n, fname):
    fig, ax = plt.subplots(figsize=(5.2, 3.3)); c = fit.c
    ax.errorbar(c, fit.h_incl, np.sqrt(fit.h_incl), fmt="o", ms=2.5, color=INK, lw=0.8, label="Data")
    ax.plot(c, fit.f_incl(c, *fit.incl), color=C_RAW, lw=1.6, label="Fit")
    ax.plot(c, fit.bkg(c, *fit.bpar), color=C_COR, lw=1.4, ls="--", label="Background")
    ax.set_xlabel("mass [GeV]"); ax.set_ylabel(f"candidates / {fit.bw*1000:.1f} MeV")
    ax.set_title(f"{PARTS[part]['label']}: inclusive mass fit ({n:,} candidates)", fontsize=8.5, loc="left")
    ax.text(0.02, 0.96, f"N_sig = {fit.incl[0]:,.0f} ± {fit.incl_err[0]:,.0f}\nσ_eff = {fit.sigma_eff()*1000:.1f} MeV "
            f"(MC {fit.sigma_eff(fit.mc_s1)*1000:.1f})\nχ²/ndf = {fit.incl_chi2:.2f}", transform=ax.transAxes, va="top", fontsize=7.5)
    ax.legend(fontsize=7.5, loc="upper right"); fig.savefig(f"fig/{fname}.png", bbox_inches="tight"); plt.close(fig)


# ------------------------------------------------------------------ studies
def make_fit(part, nb):
    d, m = load(part)
    lo, hi = d["Bmass"].min(), d["Bmass"].max()
    fit = MassFit(part, d["Bmass"], (np.floor(lo * 100) / 100, np.ceil(hi * 100) / 100), nb, m["Bmass"], m["pThatreweight"])
    mass_figure(fit, part, len(d["Bmass"]), f"{part}_mass")
    w = m["pThatreweight"]
    R[part] = dict(n_data=int(len(d["Bmass"])), n_mc=int(len(w)), n_mc_eff=float(w.sum() ** 2 / (w ** 2).sum()),
                   n_sig=float(fit.incl[0]), n_sig_err=float(fit.incl_err[0]), sigma_eff_MeV=float(fit.sigma_eff() * 1000),
                   sigma_eff_mc_MeV=float(fit.sigma_eff(fit.mc_s1) * 1000), mass_chi2ndf=float(fit.incl_chi2))
    return fit, d, m


def study(part, nb, variables, free_shape=True):
    fit, d, m = make_fit(part, nb)
    ndf = PARTS[part]["ndf"]; items = []
    for var in variables:
        nd = ndf if var == "Bchi2Prob" else NDF_SUB[var]; edges = edges_for(var)
        yd, ed = binned_yields(fit, d["Bmass"], d[var], edges, free_shape)
        res = kfit(yd, ed, edges, m[var], m["pThatreweight"], nd)
        mcd = MCdist(m[var], m["pThatreweight"])
        R[part][var] = res
        items.append((edges, yd, ed, mcd.hist(edges, 1.0, nd), mcd.hist(edges, res["k"], nd),
                      f"{PARTS[part]['label']}: {var} (ndf={nd})", var, f"MC corr. k={res['k']:.3f}"))
    fig_compare(items, f"{part}_compare")
    return fit, d, m


def unweighted_check(part, fit, d, m):
    var = "Bchi2Prob"; edges = edges_for(var); ndf = PARTS[part]["ndf"]
    yd, ed = binned_yields(fit, d["Bmass"], d[var], edges)
    R[part]["Bchi2Prob_unweightedMC"] = kfit(yd, ed, edges, m[var], np.ones_like(m[var]), ndf)


def pull_study(part, fit, d, m):
    mu, se = fit.mu, fit.sigma_eff()
    SR = np.abs(d["Bmass"] - mu) < 2.5 * se
    SB = (np.abs(d["Bmass"] - mu) > 4 * se) & (np.abs(d["Bmass"] - mu) < 8 * se)
    sc = (2.5 * se) / (4 * se)
    edges = np.linspace(-6, 6, 121); cen = 0.5 * (edges[1:] + edges[:-1])

    def width(h):
        h = np.clip(h, 0, None); c = np.cumsum(h) / h.sum()
        return 0.5 * (np.interp(0.84, c, cen) - np.interp(0.16, c, cen))

    def arr(a, kind):
        v = np.concatenate([a[f"Btrk1{kind}1"] / a[f"Btrk1{kind}Error1"], a[f"Btrk2{kind}1"] / a[f"Btrk2{kind}Error1"]])
        return v, np.concatenate([a["Btrk1Pt"], a["Btrk2Pt"]]), np.abs(np.concatenate([a["Btrk1Eta"], a["Btrk2Eta"]]))

    classes = [("all", lambda pt, eta: pt > 0), ("pT<0.95", lambda pt, eta: pt < 0.95), ("pT≥0.95", lambda pt, eta: pt >= 0.95),
               ("|η|<1", lambda pt, eta: eta < 1), ("|η|≥1", lambda pt, eta: eta >= 1)]
    out, hists = {}, {}
    for kind in ("Dxy", "Dz"):
        vs, ps, es = arr({k: v[SR] for k, v in d.items()}, kind); vb, pb, eb = arr({k: v[SB] for k, v in d.items()}, kind)
        vm, pm, em = arr(m, kind); wm = np.concatenate([m["pThatreweight"]] * 2)
        out[kind] = {}
        for cname, f in classes:
            a, b, c = f(ps, es), f(pb, eb), f(pm, em)
            hs = np.histogram(vs[a], edges)[0] - sc * np.histogram(vb[b], edges)[0]
            hm = np.histogram(vm[c], edges, weights=wm[c])[0]
            bs = []
            for _ in range(50):
                ia, ib, ic = RNG.integers(0, a.sum(), a.sum()), RNG.integers(0, b.sum(), b.sum()), RNG.integers(0, c.sum(), c.sum())
                bs.append(width(np.histogram(vs[a][ia], edges)[0] - sc * np.histogram(vb[b][ib], edges)[0])
                          / width(np.histogram(vm[c][ic], edges, weights=wm[c][ic])[0]))
            out[kind][cname] = dict(w_data=float(width(hs)), w_mc=float(width(hm)), ratio=float(width(hs) / width(hm)), err=float(np.std(bs)))
            if cname == "all":
                hists[kind] = (hs, hm)
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.3), gridspec_kw=dict(wspace=0.3))
    for ax, kind in zip(axs[:2], ("Dxy", "Dz")):
        hs, hm = hists[kind]
        ax.stairs(hm / hm.sum(), edges, color=C_RAW, lw=1.5, label="MC")
        ax.errorbar(cen, hs / hs.sum(), np.sqrt(np.abs(hs)) / hs.sum(), fmt="o", ms=2.2, color=INK, lw=0.7, label="Data (SB-subtracted)")
        ax.set_yscale("log"); ax.set_ylim(1e-5, 0.1); ax.set_xlabel(f"pion {kind.lower()} / σ({kind.lower()}) w.r.t. PV"); ax.set_ylabel("normalized")
        r = out[kind]["all"]; ax.set_title(f"{kind} pull: width data/MC = {r['ratio']:.3f} ± {r['err']:.3f}", fontsize=8.5, loc="left"); ax.legend(fontsize=7)
    ax = axs[2]; names = [c[0] for c in classes]; x = np.arange(len(names))
    for j, (kind, col, mk) in enumerate((("Dxy", C_RAW, "o"), ("Dz", C_COR, "s"))):
        ax.errorbar(x + (j - 0.5) * 0.18, [out[kind][n]["ratio"] for n in names], [out[kind][n]["err"] for n in names],
                    fmt=mk, color=col, ms=5, lw=1.2, mfc="white" if j else col, label=f"{kind} pull")
    ax.axhline(1, color=INK2, lw=0.8); ax.set_xticks(x); ax.set_xticklabels(names, fontsize=8)
    ax.set_ylabel("pull core width: data / MC"); ax.set_ylim(0.95, 1.2); ax.set_title("by pion class (16–84% half-width)", fontsize=8.5, loc="left"); ax.legend(fontsize=7)
    fig.savefig(f"fig/{part}_pulls.png", bbox_inches="tight"); plt.close(fig)
    R[part]["pulls"] = out
    R[part]["PVnchi2"] = dict(data_SR=float(d["PVnchi2"][SR].mean()), data_SB=float(d["PVnchi2"][SB].mean()),
                              mc=float(np.average(m["PVnchi2"], weights=m["pThatreweight"])))


def universality(part, fit, d, m):
    ndf = PARTS[part]["ndf"]
    edges = np.array([0.005, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0])
    minpt = lambda a: np.minimum(a["Btrk1Pt"], a["Btrk2Pt"]); maxeta = lambda a: np.maximum(np.abs(a["Btrk1Eta"]), np.abs(a["Btrk2Eta"]))
    dmu = lambda a: np.abs(a["Bujmass"] - 3.0969)
    subsets = [("min pT(π) < 0.95 GeV", lambda a: minpt(a) < 0.95), ("min pT(π) ≥ 0.95 GeV", lambda a: minpt(a) >= 0.95),
               ("pT 7.5–10 GeV", lambda a: a["Bpt"] < 10), ("pT 10–15 GeV", lambda a: (a["Bpt"] >= 10) & (a["Bpt"] < 15)),
               ("pT 15–50 GeV", lambda a: a["Bpt"] >= 15), ("N_trk < 30", lambda a: a["nChargedTracks"] < 30),
               ("N_trk 30–45", lambda a: (a["nChargedTracks"] >= 30) & (a["nChargedTracks"] < 45)), ("N_trk ≥ 45", lambda a: a["nChargedTracks"] >= 45),
               ("max |η(π)| < 1", lambda a: maxeta(a) < 1), ("max |η(π)| ≥ 1", lambda a: maxeta(a) >= 1),
               ("|m(μμ) − m(J/ψ)| < 25 MeV", lambda a: dmu(a) < 0.025), ("|m(μμ) − m(J/ψ)| ≥ 25 MeV", lambda a: dmu(a) >= 0.025)]
    out = []
    for nm, f in subsets:
        sd, sm = f(d), f(m)
        yd, ed = binned_yields(fit, d["Bmass"][sd], d["Bchi2Prob"][sd], edges)
        r = kfit(yd, ed, edges, m["Bchi2Prob"][sm], m["pThatreweight"][sm], ndf, nboot=100)
        out.append(dict(subset=nm, k=r["k"], err=r["err"], chi2min=r["chi2min"], ndof=r["ndof"], nsig=float(np.nansum(yd))))
    R[part]["universality"] = out
    kg, eg = R[part]["Bchi2Prob"]["k"], R[part]["Bchi2Prob"]["err"]
    ks, es = np.array([o["k"] for o in out]), np.array([o["err"] for o in out])
    R[part]["universality_chi2"] = dict(chi2=float(np.sum((ks - kg) ** 2 / es ** 2)), n=len(out))
    fig, ax = plt.subplots(figsize=(7.0, 3.8)); y = np.arange(len(out))[::-1]
    ax.axvspan(kg - eg, kg + eg, color=C_COR, alpha=0.18, lw=0); ax.axvline(kg, color=C_COR, lw=1.2, ls="--", label=f"global k = {kg:.3f} ± {eg:.3f}")
    ax.axvline(1.0, color=INK2, lw=0.8)
    ax.errorbar(ks, y, xerr=es, fmt="o", color=C_RAW, ms=4.5, lw=1.2)
    for yy, o in zip(y, out):
        ax.text(1.215, yy, f"{o['k']:.3f} ± {o['err']:.3f}", va="center", fontsize=7, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([o["subset"] for o in out], fontsize=8)
    ax.set_xlabel("k fitted on Bchi2Prob in the subset"); ax.set_xlim(0.95, 1.28)
    ax.set_title(f"{PARTS[part]['label']}: k in subsets (one fit per subset)", fontsize=8.5, loc="left"); ax.legend(fontsize=7.5, loc="lower left")
    fig.savefig(f"fig/{part}_universality.png", bbox_inches="tight"); plt.close(fig)


def effs_from(edges, yd, ed, mraw, mcor, cuts):
    out = []
    for cut in cuts:
        a = edges[:-1] >= cut - 1e-9
        A, B = np.nansum(yd[a]), np.nansum(yd[~a]); sA, sB = np.sqrt(np.nansum(ed[a] ** 2)), np.sqrt(np.nansum(ed[~a] ** 2))
        e_ = A / (A + B); se = np.sqrt((B * sA) ** 2 + (A * sB) ** 2) / (A + B) ** 2
        er, ec = mraw[0][a].sum() / mraw[0].sum(), mcor[0][a].sum() / mcor[0].sum()
        out.append(dict(cut=cut, data=float(e_), data_err=float(se), mc_raw=float(er), mc_cor=float(ec),
                        ratio_raw=float(e_ / er), ratio_cor=float(e_ / ec)))
    return out


def closure(part, fit, d, m):
    ndf = PARTS[part]["ndf"]; var = "Bchi2Prob"; edges = edges_for(var)
    de, do = d["EvtNo"] % 2 == 0, d["EvtNo"] % 2 == 1; me, mo = m["EvtNo"] % 2 == 0, m["EvtNo"] % 2 == 1
    yde, ede = binned_yields(fit, d["Bmass"][de], d[var][de], edges)
    rk = kfit(yde, ede, edges, m[var][me], m["pThatreweight"][me], ndf)
    ydo, edo = binned_yields(fit, d["Bmass"][do], d[var][do], edges)
    mcd = MCdist(m[var][mo], m["pThatreweight"][mo]); hr, hc = mcd.hist(edges, 1.0, ndf), mcd.hist(edges, rk["k"], ndf)
    c_raw, nd_ = chi2_dm(ydo, edo, *hr); c_cor, _ = chi2_dm(ydo, edo, *hc)
    R[part]["closure"] = dict(k_even=rk["k"], k_even_err=rk["err"], chi2_raw=float(c_raw), chi2_cor=float(c_cor), ndof=int(nd_),
                              eff=effs_from(edges, ydo, edo, hr, hc, (0.05, 0.10, 0.20, 0.30)))
    fig_compare([(edges, ydo, edo, hr, hc, f"{PARTS[part]['label']}: ODD events; k from EVEN events", var,
                  f"MC corr. (k_even = {rk['k']:.3f})")], f"{part}_closure")


def transfer(part, nb, k, kerr, src, free_shape=True, edges=None, fit_dm=None):
    fit, d, m = (fit_dm if fit_dm else make_fit(part, nb))
    ndf = PARTS[part]["ndf"]; var = "Bchi2Prob"; edges = edges_for(var) if edges is None else edges
    yd, ed = binned_yields(fit, d["Bmass"], d[var], edges, free_shape)
    mcd = MCdist(m[var], m["pThatreweight"]); hr, hc = mcd.hist(edges, 1.0, ndf), mcd.hist(edges, k, ndf)
    c_raw, nd_ = chi2_dm(yd, ed, *hr); c_cor, _ = chi2_dm(yd, ed, *hc)
    own = kfit(yd, ed, edges, m[var], m["pThatreweight"], ndf)
    R[part]["transfer"] = dict(k_applied=k, k_applied_err=kerr, source=src, chi2_raw=float(c_raw), chi2_cor=float(c_cor), ndof=int(nd_),
                               k_own=own["k"], k_own_err=own["err"], eff=effs_from(edges, yd, ed, hr, hc, (0.05, 0.10, 0.20)))
    fig_compare([(edges, yd, ed, hr, hc, f"{PARTS[part]['label']}: Bchi2Prob (ndf={ndf})", var, f"MC corr. ({src} k = {k:.3f})")],
                f"{part}_transfer")


def kcurve_figure(entries, fname):
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    cols = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6", "#e34948", "#008300"]
    lss = ["-", "--", ":", "-.", (0, (5, 1)), (0, (1, 1)), (0, (3, 1, 1, 1)), (0, (6, 2, 1, 2))]
    for i, (part, var) in enumerate(entries):
        r = R[part][var]
        ax.plot(KS, np.array(r["curve"]) / max(r["ndof"], 1), color=cols[i], ls=lss[i], lw=1.5,
                label=f"{PARTS[part]['label']} {var}: k = {r['k']:.3f} ± {r['err']:.3f}")
    ax.set_xlabel("track-error scale factor k   (MC χ² → k²·χ²)"); ax.set_ylabel("(χ² − χ²_min) / ndf  vs data")
    ax.set_ylim(0, 6); ax.set_xlim(0.95, 1.25); ax.legend(fontsize=6.5, loc="upper center", ncol=1)
    fig.savefig(f"fig/{fname}.png", bbox_inches="tight"); plt.close(fig)


def summary_figure(entries, fname):
    fig, ax = plt.subplots(figsize=(6.8, 3.4)); y = np.arange(len(entries))[::-1]
    for yy, (lab, k, e, col) in zip(y, entries):
        ax.errorbar([k], [yy], xerr=[e], fmt="o", color=col, ms=5, lw=1.4)
        ax.text(1.215, yy, f"{k:.3f} ± {e:.3f}", va="center", fontsize=7.5, color=INK2)
    ax.axvline(1.0, color=INK2, lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels([e[0] for e in entries], fontsize=8); ax.set_xlim(0.97, 1.27)
    ax.set_xlabel("k (data spread / MC spread, relative to the assigned errors)")
    ax.set_title("All measurements of k", fontsize=8.5, loc="left")
    fig.savefig(f"fig/{fname}.png", bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    fit, d, m = study("ppRef_psi2S", 70, ["Bchi2Prob", "BujvProb", "BtktkvProb"])
    unweighted_check("ppRef_psi2S", fit, d, m)
    pull_study("ppRef_psi2S", fit, d, m)
    universality("ppRef_psi2S", fit, d, m)
    closure("ppRef_psi2S", fit, d, m)
    kpsi, epsi = R["ppRef_psi2S"]["Bchi2Prob"]["k"], R["ppRef_psi2S"]["Bchi2Prob"]["err"]

    fitb, db, mb = study("ppRef_Bp", 55, ["Bchi2Prob", "BujvProb"])
    transfer("ppRef_Bp", 55, kpsi, epsi, "ψ(2S)", fit_dm=(fitb, db, mb))
    kbp, ebp = R["ppRef_Bp"]["Bchi2Prob"]["k"], R["ppRef_Bp"]["Bchi2Prob"]["err"]

    transfer("ppRef_X3872", 40, kpsi, epsi, "ψ(2S)", free_shape=False, edges=np.array([0.005, 0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 1.0]))
    for part in ("PbPb23_Bp", "PbPb24_Bp"):
        study(part, 55, ["Bchi2Prob", "BujvProb"])

    kcurve_figure([("ppRef_psi2S", "Bchi2Prob"), ("ppRef_psi2S", "BujvProb"), ("ppRef_psi2S", "BtktkvProb"),
                   ("ppRef_Bp", "Bchi2Prob"), ("ppRef_Bp", "BujvProb"), ("PbPb23_Bp", "Bchi2Prob"), ("PbPb24_Bp", "Bchi2Prob")], "kscan_all")
    P = R["ppRef_psi2S"]
    summary_figure([
        ("ppRef ψ(2S): Bchi2Prob (4 trk + mass c.)", P["Bchi2Prob"]["k"], P["Bchi2Prob"]["err"], C_RAW),
        ("ppRef ψ(2S): Bchi2Prob, unweighted MC", P["Bchi2Prob_unweightedMC"]["k"], P["Bchi2Prob_unweightedMC"]["err"], C_RAW),
        ("ppRef ψ(2S): BujvProb (μμ vertex)", P["BujvProb"]["k"], P["BujvProb"]["err"], C_RAW),
        ("ppRef ψ(2S): BtktkvProb (ππ vertex)", P["BtktkvProb"]["k"], P["BtktkvProb"]["err"], C_RAW),
        ("ppRef ψ(2S): pion dxy pull width", P["pulls"]["Dxy"]["all"]["ratio"], P["pulls"]["Dxy"]["all"]["err"], C_RAW),
        ("ppRef ψ(2S): pion dz pull width", P["pulls"]["Dz"]["all"]["ratio"], P["pulls"]["Dz"]["all"]["err"], C_RAW),
        ("ppRef B+: Bchi2Prob (3 trk + mass c.)", R["ppRef_Bp"]["Bchi2Prob"]["k"], R["ppRef_Bp"]["Bchi2Prob"]["err"], C_COR),
        ("ppRef B+: BujvProb (μμ vertex)", R["ppRef_Bp"]["BujvProb"]["k"], R["ppRef_Bp"]["BujvProb"]["err"], C_COR),
        ("ppRef X(3872): Bchi2Prob (own fit)", R["ppRef_X3872"]["transfer"]["k_own"], R["ppRef_X3872"]["transfer"]["k_own_err"], C_COR),
        ("PbPb 2023 B+: Bchi2Prob", R["PbPb23_Bp"]["Bchi2Prob"]["k"], R["PbPb23_Bp"]["Bchi2Prob"]["err"], "#1baf7a"),
        ("PbPb 2024 data vs 2023 MC B+: Bchi2Prob", R["PbPb24_Bp"]["Bchi2Prob"]["k"], R["PbPb24_Bp"]["Bchi2Prob"]["err"], "#1baf7a"),
    ], "k_summary")
    for part in R:
        for v in list(R[part]):
            if isinstance(R[part][v], dict):
                R[part][v].pop("curve", None)
    json.dump(R, open("results.json", "w"), indent=1)
    print("done")
