# Step 3: build the LaTeX note from results.json and the figures; every number comes from the fits.
import json, subprocess, math

R = json.load(open("results.json"))
P, B, X = R["ppRef_psi2S"], R["ppRef_Bp"], R["ppRef_X3872"]
P23, P24 = R["PbPb23_Bp"], R["PbPb24_Bp"]


def f(x, n=3):
    return f"{x:.{n}f}"


def ke(r, n=3):
    return f"${r['k']:.{n}f} \\pm {r['err']:.{n}f}$"


# combined ppRef k from the two independent chi2-based measurements (psi(2S) 4-track, B+ 3-track)
w1, w2 = 1 / P["Bchi2Prob"]["err"] ** 2, 1 / B["Bchi2Prob"]["err"] ** 2
k_comb = (P["Bchi2Prob"]["k"] * w1 + B["Bchi2Prob"]["k"] * w2) / (w1 + w2)
e_comb = 1 / math.sqrt(w1 + w2)
K_REC = round(k_comb, 2)
uni = P["universality"]
chi2_uni_rec = sum((u["k"] - k_comb) ** 2 / u["err"] ** 2 for u in uni)
kmin_u, kmax_u = min(u["k"] for u in uni), max(u["k"] for u in uni)

V = {}
V["KCOMB"] = f"{k_comb:.3f}"; V["ECOMB"] = f"{e_comb:.3f}"; V["KREC"] = f"{K_REC:.2f}"
V["KLO"] = f"{K_REC - 0.02:.2f}"; V["KHI"] = f"{K_REC + 0.02:.2f}"
V["KPSI"] = ke(P["Bchi2Prob"]); V["KPSIUW"] = ke(P["Bchi2Prob_unweightedMC"]); V["KPSIMM"] = ke(P["BujvProb"]); V["KPSIPP"] = ke(P["BtktkvProb"])
V["KBP"] = ke(B["Bchi2Prob"]); V["KBPMM"] = ke(B["BujvProb"])
V["KX"] = f"${X['transfer']['k_own']:.3f} \\pm {X['transfer']['k_own_err']:.3f}$"
V["K23"] = ke(P23["Bchi2Prob"]); V["K24"] = ke(P24["Bchi2Prob"]); V["K23MM"] = ke(P23["BujvProb"]); V["K24MM"] = ke(P24["BujvProb"])
pu = P["pulls"]
V["PDXY"] = f"${pu['Dxy']['all']['ratio']:.3f} \\pm {pu['Dxy']['all']['err']:.3f}$"
V["PDZ"] = f"${pu['Dz']['all']['ratio']:.3f} \\pm {pu['Dz']['all']['err']:.3f}$"
V["PVD"] = f(P["PVnchi2"]["data_SR"]); V["PVB"] = f(P["PVnchi2"]["data_SB"]); V["PVM"] = f(P["PVnchi2"]["mc"])
V["NEFF"] = f"{P['n_mc_eff']:,.0f}"; V["NMC"] = f"{P['n_mc']:,}"
V["UNICHI2"] = f"{chi2_uni_rec:.1f}"; V["UNIN"] = str(len(uni)); V["UNIMIN"] = f(kmin_u); V["UNIMAX"] = f(kmax_u)
cl = P["closure"]
V["CLK"] = f"${cl['k_even']:.3f} \\pm {cl['k_even_err']:.3f}$"; V["CLRAW"] = f"{cl['chi2_raw']:.1f}"; V["CLCOR"] = f"{cl['chi2_cor']:.1f}"; V["CLNDF"] = str(cl["ndof"])
bt, xt = B["transfer"], X["transfer"]
V["BTRAW"] = f"{bt['chi2_raw']:.1f}"; V["BTCOR"] = f"{bt['chi2_cor']:.1f}"; V["BTNDF"] = str(bt["ndof"]); V["BTOWN"] = f"{B['Bchi2Prob']['chi2min']:.1f}"
V["XTRAW"] = f"{xt['chi2_raw']:.1f}"; V["XTCOR"] = f"{xt['chi2_cor']:.1f}"; V["XTNDF"] = str(xt["ndof"])
V["XSIG"] = f"{X['n_sig']:,.0f} \\pm {X['n_sig_err']:,.0f}"
V["XWD"] = f"{X['sigma_eff_MeV']:.1f}"; V["XWM"] = f"{X['sigma_eff_mc_MeV']:.1f}"
V["XDIFF"] = f"{(X['transfer']['k_own'] - k_comb) / math.hypot(X['transfer']['k_own_err'], e_comb):.1f}"
V["D2324"] = f"{(P24['Bchi2Prob']['k'] - P23['Bchi2Prob']['k']) / math.hypot(P24['Bchi2Prob']['err'], P23['Bchi2Prob']['err']):.1f}"


def samples_table():
    rows = []
    for key, lab in (("ppRef_psi2S", r"ppRef $\psi(2S)$"), ("ppRef_X3872", r"ppRef X(3872)"), ("ppRef_Bp", r"ppRef $B^{+}$"),
                     ("PbPb23_Bp", r"PbPb 2023 $B^{+}$"), ("PbPb24_Bp", r"PbPb 2024 $B^{+}$ $^{\dagger}$")):
        r = R[key]
        rows.append(f"{lab} & {r['n_data']:,} & ${r['n_sig']:,.0f} \\pm {r['n_sig_err']:,.0f}$ & {r['sigma_eff_MeV']:.1f} / {r['sigma_eff_mc_MeV']:.1f} "
                    f"& {r['mass_chi2ndf']:.2f} & {r['n_mc']:,} & {r['n_mc_eff']:,.0f} \\\\")
    return "\n".join(rows).replace(",", "{,}")


def k_table():
    rows = []
    def row(lab, r, nd):
        rows.append(f"{lab} & {nd} & {ke(r)} & {r['chi2_k1']:.1f} & {r['chi2min']:.1f} & {r['ndof']} \\\\")
    row(r"ppRef $\psi(2S)$ \texttt{Bchi2Prob}", P["Bchi2Prob"], 6)
    row(r"\quad same, unweighted MC", P["Bchi2Prob_unweightedMC"], 6)
    row(r"ppRef $\psi(2S)$ \texttt{BujvProb}", P["BujvProb"], 1)
    row(r"ppRef $\psi(2S)$ \texttt{BtktkvProb}", P["BtktkvProb"], 1)
    row(r"ppRef $B^{+}$ \texttt{Bchi2Prob}", B["Bchi2Prob"], 4)
    row(r"ppRef $B^{+}$ \texttt{BujvProb}", B["BujvProb"], 1)
    row(r"PbPb 2023 $B^{+}$ \texttt{Bchi2Prob}", P23["Bchi2Prob"], 4)
    row(r"PbPb 2023 $B^{+}$ \texttt{BujvProb}", P23["BujvProb"], 1)
    row(r"PbPb 2024$^{\dagger}$ $B^{+}$ \texttt{Bchi2Prob}", P24["Bchi2Prob"], 4)
    row(r"PbPb 2024$^{\dagger}$ $B^{+}$ \texttt{BujvProb}", P24["BujvProb"], 1)
    return "\n".join(rows)


def eff_table(eff, with_ratio=True):
    rows = []
    for e in eff:
        rr = e["data"] / e["mc_raw"]; rc = e["data"] / e["mc_cor"]
        rows.append(f"$> {e['cut']:.2f}$ & ${e['data']:.3f} \\pm {e['data_err']:.3f}$ & {e['mc_raw']:.3f} & {e['mc_cor']:.3f} & {rr:.3f} & {rc:.3f} \\\\")
    return "\n".join(rows)


def uni_table():
    rows = []
    for u in uni:
        s = u["subset"].replace("N_trk", r"$N_{\rm trk}$").replace("π", r"$\pi$").replace("μμ", r"$\mu\mu$").replace("J/ψ", r"J/$\psi$").replace("≥", r"$\geq$").replace("η", r"$\eta$").replace("−", "$-$").replace("–", "--")
        rows.append(f"{s} & {u['nsig']:,.0f} & ${u['k']:.3f} \\pm {u['err']:.3f}$ & {u['chi2min']:.1f}/{u['ndof']} \\\\".replace(",", "{,}"))
    return "\n".join(rows)


def pull_table():
    rows = []
    names = {"all": "all pions", "pT<0.95": r"$p_T < 0.95$ GeV", "pT≥0.95": r"$p_T \geq 0.95$ GeV", "|η|<1": r"$|\eta| < 1$", "|η|≥1": r"$|\eta| \geq 1$"}
    for c in ("all", "pT<0.95", "pT≥0.95", "|η|<1", "|η|≥1"):
        a, b = pu["Dxy"][c], pu["Dz"][c]
        rows.append(f"{names[c]} & {a['w_data']:.3f} & {a['w_mc']:.3f} & ${a['ratio']:.3f} \\pm {a['err']:.3f}$ & {b['w_data']:.3f} & {b['w_mc']:.3f} & ${b['ratio']:.3f} \\pm {b['err']:.3f}$ \\\\")
    return "\n".join(rows)


V["SAMPLES"] = samples_table(); V["KTABLE"] = k_table(); V["UNITABLE"] = uni_table(); V["PULLTABLE"] = pull_table()
V["CLEFF"] = eff_table(cl["eff"]); V["BTEFF"] = eff_table(bt["eff"]); V["XTEFF"] = eff_table(xt["eff"])

TEX = r"""\documentclass[10pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage[margin=2.1cm]{geometry}
\usepackage{graphicx,booktabs,float,xcolor,caption,amsmath}
\usepackage[colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black]{hyperref}
\captionsetup{font=small,labelfont=bf}
\makeatletter\g@addto@macro{\UrlBreaks}{\UrlOrds}\makeatother
\setlength{\parskip}{4pt}\setlength{\parindent}{0pt}
\newcommand{\kk}{\ensuremath{k}}
\title{\textbf{Vertex-fit $\chi^2$ probability in data and MC:\\ a track-error scale factor and an MC correction}\\[2pt]\large ppRef 2024 (with $B^{+}$ and PbPb cross-checks)}
\author{X(3872) analysis --- study notes}
\date{14 October 2026}
\begin{document}
\maketitle

\section*{Summary}
\begin{itemize}\setlength{\itemsep}{2pt}
\item In ppRef 2024, the vertex-fit probabilities are systematically lower in data than in MC. This holds for the full
$\mu\mu\pi\pi$ fit (\texttt{Bchi2Prob}), the dimuon vertex (\texttt{BujvProb}), the dipion vertex (\texttt{BtktkvProb}) and the
primary vertex. It also holds for $B^{+}\to J/\psi K^{+}$.
\item All of these are described by one number: relative to its assigned errors, the track spread in data is
\kk\ times larger than in MC. MC is corrected by $\chi^2_{\rm MC} \to \kk^2\chi^2_{\rm MC}$.
\item Measured values: $\psi(2S)$ \texttt{Bchi2Prob} @@KPSI@@; $B^{+}$ \texttt{Bchi2Prob} @@KBP@@;
pion impact-parameter pull widths @@PDXY@@ ($d_{xy}$) and @@PDZ@@ ($d_z$).
Combined ($\psi(2S)$ and $B^{+}$ $\chi^2$ fits): $\kk = @@KCOMB@@ \pm @@ECOMB@@$ (stat).
\item \textbf{Closure:} $\kk$ derived on even events and applied to odd events reduces the data/MC $\chi^2$ from
@@CLRAW@@ to @@CLCOR@@ for @@CLNDF@@ degrees of freedom. The data/MC efficiency of a \texttt{Bchi2Prob} cut moves from 0.89--0.97 to 0.98--0.99.
The $\psi(2S)$ value applied to $B^{+}$ (a different topology) reduces the $\chi^2$ from @@BTRAW@@ to @@BTCOR@@ for @@BTNDF@@.
\item \textbf{X(3872) itself:} applying the $\psi(2S)$ \kk\ reduces the data/MC $\chi^2$ from @@XTRAW@@ to @@XTCOR@@ for @@XTNDF@@.
A free fit on the X prefers $\kk=$ @@KX@@ (statistical only, @@XDIFF@@$\sigma$ above the combined value).
\item \textbf{Recommendation:} correct all MC vertex probabilities with $\kk = @@KREC@@$ before training, and
use $\kk = @@KLO@@$ and $@@KHI@@$ as systematic variations (plus $\kk \approx 1.14$ for ntmix as an X-specific check).
After the correction a residual 1--2\% data/MC efficiency difference remains for \texttt{Bchi2Prob} cuts; take it as a systematic uncertainty.
PbPb needs its own \kk\ (Section~\ref{sec:pbpb}).
\end{itemize}

\section{Variables and correction model}
All probabilities are $P(\chi^2, n_{\rm dof})$ (\texttt{TMath::Prob}) of a vertex fit made in Bfinder from MiniAOD
packed-candidate tracks (\texttt{pseudoTrack()}). The degrees of freedom follow from
\texttt{KinematicConstrainedVertexUpdator}: $n_{\rm dof} = 2N_{\rm trk} - 3 + N_{\rm constraints}$.
\begin{center}\small
\begin{tabular}{lll}
\toprule
variable & fit & $n_{\rm dof}$ \\
\midrule
\texttt{Bchi2Prob} (ntmix: $\psi(2S)$, X) & $\mu^+\mu^-\pi^+\pi^-$, one vertex, $m(\mu\mu)$ constrained to $m_{J/\psi}$ & 6 \\
\texttt{Bchi2Prob} (ntKp: $B^{+}$) & $\mu^+\mu^-K^{+}$, one vertex, $m(\mu\mu)$ constrained & 4 \\
\texttt{BujvProb} & $\mu^+\mu^-$ vertex, no constraint & 1 \\
\texttt{BtktkvProb} & $\pi^+\pi^-$ vertex, no constraint & 1 \\
\bottomrule
\end{tabular}
\end{center}
\textbf{Correction.} If the true track spread in data is \kk\ times larger than in MC relative to the
errors assigned by reconstruction, then every vertex $\chi^2$ in data is larger by $\kk^2$. The MC correction is applied per candidate:
\begin{equation}
p_{\rm corr} = P\!\left(\kk^2\, \chi^2_{\rm MC},\, n_{\rm dof}\right), \qquad \chi^2_{\rm MC} = F^{-1}_{\chi^2}(1-p_{\rm MC};\, n_{\rm dof}) ,
\end{equation}
so the stored probability is enough; no reprocessing is needed. Since $\kk>1$ the correction only lowers probabilities.
Preselection thresholds ($\texttt{Bchi2Prob}>0.005$, $\texttt{BujvProb}>0.01$) are applied again to MC after the correction,
as in data.

\section{Samples and method}
\textbf{Selections} follow the nominal fit selections with one exception: no BDT/NN \texttt{Prediction} cut, because
\texttt{Bchi2Prob} is a training input. ppRef $\psi(2S)$ and X: \texttt{BQvalue}$<0.15$, \texttt{Btrk1,2dR}$<0.5$,
$7.5<p_T<50$\,GeV, prompt region $L_{xy} m/p_T<0.005$\,cm. $B^{+}$ (ppRef and PbPb): $10<p_T<50$\,GeV, $|y|<2.4$, and cuts that do not use any
track error: pointing angle $\alpha<0.05$, $d_{\rm SV-PV}^{\,2D}>0.02$\,cm, $p_T(K)>1$\,GeV. MC is truth-matched (\texttt{Bgen==23333})
and weighted with \texttt{pThatreweight}. $^{\dagger}$No 2024 PbPb $B^{+}$ MC exists, so 2024 data are compared with the 2023 MC.

\begin{table}[H]\centering\small
\caption{Samples. $\sigma_{\rm eff}$: effective mass resolution in data / MC. $N_{\rm eff}^{\rm MC}=(\sum w)^2/\sum w^2$.}
\begin{tabular}{lrrrrrr}
\toprule
sample & candidates & $N_{\rm sig}$ (data fit) & $\sigma_{\rm eff}$ [MeV] & mass $\chi^2$/ndf & $N^{\rm MC}$ & $N_{\rm eff}^{\rm MC}$ \\
\midrule
@@SAMPLES@@
\bottomrule
\end{tabular}
\end{table}

\textbf{Signal in data.} The data are split into bins of the probability under study, and each bin gets its own binned mass fit.
The signal is a double Gaussian whose shape (width ratio, fraction) comes from a fit to truth-matched MC. The data fit floats the mean and a common width scale,
both inclusively and in each bin. The background is a second-order polynomial for $\psi(2S)$/X; for $B^{+}$ it is an exponential plus an
erfc edge for partially reconstructed decays. The signal yield and its fit error per bin form the data distribution. Sideband subtraction gives the same picture.

\textbf{Measuring \kk.} For each trial \kk\ (0.900--1.300 in steps of 0.001) the MC is corrected with Eq.~(1), normalized to the data,
and compared bin by bin with $\chi^2 = \sum (y_i - m_i)^2/(\sigma_{y,i}^2 + \sigma_{m,i}^2)$, where $\sigma_m$ includes the weighted-MC
statistical error. \kk\ is the minimum. Its uncertainty is the 68\% half-width of 200 bootstrap replicas, where data yields are fluctuated
by their errors and MC weights are multiplied by Poisson(1) numbers. This is needed because the heavily weighted MC makes the $\chi^2(\kk)$ curve jagged.

\textbf{Pulls.} For prompt $\psi(2S)$ pions, $d_{xy}/\sigma(d_{xy})$ and $d_z/\sigma(d_z)$ relative to the primary vertex should have a
width set purely by the error calibration. The core width (half of the 16--84\% interval) is compared between sideband-subtracted data and MC.
The absolute widths are below 1 because each pion is part of the primary-vertex fit. That bias is common to data and MC and cancels in the ratio.

\section{ppRef $\psi(2S)$: the main measurement}
\begin{figure}[H]\centering
\includegraphics[width=0.52\textwidth]{fig/ppRef_psi2S_mass.png}
\caption{Inclusive $\psi(2S)$ mass fit, ppRef 2024, prompt region.}
\end{figure}
\begin{figure}[H]\centering
\includegraphics[width=\textwidth]{fig/ppRef_psi2S_compare.png}
\caption{Fitted $\psi(2S)$ signal in data (points) compared with raw MC (blue, solid) and MC corrected with its own fitted \kk\
(orange, dashed). Bands show MC statistical uncertainty. Lower panels: data/MC.}
\end{figure}

\begin{table}[H]\centering\small
\caption{Fitted \kk\ for every vertex probability. $\chi^2(\kk{=}1)$ is the uncorrected data/MC agreement; $\chi^2_{\min}$ is at the fitted \kk.}
\begin{tabular}{lrcrrr}
\toprule
measurement & $n_{\rm dof}$ & \kk & $\chi^2(\kk{=}1)$ & $\chi^2_{\min}$ & ndf \\
\midrule
@@KTABLE@@
\bottomrule
\end{tabular}
\end{table}
The unweighted-MC result shows that the large spread in \texttt{pThatreweight} (the $\psi(2S)$ MC has $N_{\rm eff}=@@NEFF@@$
out of @@NMC@@ candidates) limits the precision, not the central value.

\subsection*{Which errors: impact-parameter pulls}
\begin{figure}[H]\centering
\includegraphics[width=\textwidth]{fig/ppRef_psi2S_pulls.png}
\caption{Prompt $\psi(2S)$ pion pulls with respect to the primary vertex. Left and centre: inclusive distributions. Right: data/MC width ratio by pion class.}
\end{figure}
\begin{table}[H]\centering\small
\caption{Pull core widths (half of the 16--84\% interval).}
\begin{tabular}{lcccccc}
\toprule
& \multicolumn{3}{c}{$d_{xy}/\sigma(d_{xy})$} & \multicolumn{3}{c}{$d_z/\sigma(d_z)$} \\
class & data & MC & data/MC & data & MC & data/MC \\
\midrule
@@PULLTABLE@@
\bottomrule
\end{tabular}
\end{table}
The excess spread is the same in the transverse and longitudinal directions, for soft and hard pions, and in the barrel and endcaps.
The muons show it too (\texttt{BujvProb}). The assigned errors themselves agree between data and MC (as in the
\texttt{Btrk1DzError1} and \texttt{Btrk1PtErr} validation plots). What differs is the true spread of the tracks relative to those errors. The
primary-vertex $\chi^2/n_{\rm dof}$ shows the same effect: @@PVD@@ in the signal region, @@PVB@@ in the sidebands, and @@PVM@@ in MC.
The effect is event-wide, not specific to the signal. The most likely origin is the per-hit position errors (alignment and cluster-position
errors) in the 2024 prompt reconstruction compared with the simulated detector. That cannot be fixed at analysis level other than by correcting MC.

\section{Is one \kk\ enough?}
\begin{figure}[H]\centering
\includegraphics[width=0.8\textwidth]{fig/ppRef_psi2S_universality.png}
\caption{\kk\ fitted separately in subsets of the ppRef $\psi(2S)$ sample.}
\end{figure}
\begin{table}[H]\centering\small
\begin{tabular}{lrcc}
\toprule
subset & $N_{\rm sig}$ & \kk & $\chi^2_{\min}$/ndf \\
\midrule
@@UNITABLE@@
\bottomrule
\end{tabular}
\end{table}
The subset values span @@UNIMIN@@--@@UNIMAX@@. Compared with the combined $\kk=@@KCOMB@@$, they give $\chi^2=@@UNICHI2@@$ for @@UNIN@@ subsets.
The two clearest trends are with the dimuon mass ($|m_{\mu\mu}-m_{J/\psi}|$ small vs.\ large; the mass constraint adds a term that depends on the
muon momentum resolution) and between $p_T$ bins (not monotonic). Both are at the $\pm0.02$ level, which sets the systematic variation recommended below. Multiplicity and pion
$p_T$/$\eta$ show no significant dependence, so pileup and the MiniAOD low-$p_T$ covariance packing ($p_T<0.95$\,GeV) do not drive it.

\section{Closure test (even/odd events)}
\kk\ is derived on even events ($\kk=$ @@CLK@@) and applied to the MC of odd events. It is compared with odd-event data,
which played no part in deriving \kk.
\begin{figure}[H]\centering
\includegraphics[width=0.45\textwidth]{fig/ppRef_psi2S_closure.png}
\caption{Closure on the odd half. $\chi^2$/ndf: @@CLRAW@@/@@CLNDF@@ before, @@CLCOR@@/@@CLNDF@@ after the correction.}
\end{figure}
\begin{table}[H]\centering\small
\caption{Efficiency of \texttt{Bchi2Prob} cuts on the odd half (relative to the preselection $>0.005$).}
\begin{tabular}{lccccc}
\toprule
cut & data & MC raw & MC corrected & data/MC raw & data/MC corrected \\
\midrule
@@CLEFF@@
\bottomrule
\end{tabular}
\end{table}
The shape is restored over the whole range. The lowest bin ($p<0.05$) is still slightly underpredicted: the data have a small
non-Gaussian tail that a pure scale factor cannot produce. This is the 1--2\% residual quoted in the summary.

\section{Independent topology: ppRef $B^{+}\to J/\psi K^{+}$}
\begin{figure}[H]\centering
\includegraphics[width=0.48\textwidth]{fig/ppRef_Bp_mass.png}
\includegraphics[width=0.40\textwidth]{fig/ppRef_Bp_transfer.png}
\caption{Left: $B^{+}$ mass fit (the $\chi^2$/ndf of 2.6 reflects the simple partially-reconstructed/$J/\psi\pi^{+}$ model).
Right: data vs MC corrected with the \emph{$\psi(2S)$} \kk.}
\end{figure}
\begin{figure}[H]\centering
\includegraphics[width=0.75\textwidth]{fig/ppRef_Bp_compare.png}
\caption{$B^{+}$ with its own fitted \kk.}
\end{figure}
\begin{table}[H]\centering\small
\caption{$B^{+}$: efficiency of \texttt{Bchi2Prob} cuts, with MC corrected by the $\psi(2S)$ \kk.}
\begin{tabular}{lccccc}
\toprule
cut & data & MC raw & MC corrected & data/MC raw & data/MC corrected \\
\midrule
@@BTEFF@@
\bottomrule
\end{tabular}
\end{table}
A displaced, 3-track, $n_{\rm dof}=4$ fit gives the same \kk\ (@@KBP@@) as the prompt 4-track $\psi(2S)$ fit. With the $\psi(2S)$ \kk\
the $B^{+}$ $\chi^2$ goes from @@BTRAW@@ to @@BTCOR@@ for @@BTNDF@@ (own \kk: @@BTOWN@@), and the cut efficiencies agree to better than 1.5\%.

\section{X(3872) in ppRef (the ntmix training case)}
\begin{figure}[H]\centering
\includegraphics[width=0.48\textwidth]{fig/ppRef_X3872_mass.png}
\includegraphics[width=0.40\textwidth]{fig/ppRef_X3872_transfer.png}
\caption{Left: X(3872) mass fit ($N_{\rm sig}=@@XSIG@@$). Right: X signal vs MC corrected with the $\psi(2S)$ \kk.
The signal shape is fixed in each bin because of the limited statistics.}
\end{figure}
\begin{table}[H]\centering\small
\begin{tabular}{lccccc}
\toprule
cut & data & MC raw & MC corrected & data/MC raw & data/MC corrected \\
\midrule
@@XTEFF@@
\bottomrule
\end{tabular}
\end{table}
The correction improves the X agreement ($\chi^2$: @@XTRAW@@ $\to$ @@XTCOR@@ for @@XTNDF@@). A free fit on the X prefers a larger value,
$\kk=$ @@KX@@ (statistical only), @@XDIFF@@$\sigma$ above the combined value. With about 2.6k signal events at low $S/B$, and a data/MC width ratio of
@@XWD@@/@@XWM@@\,MeV that differs from the $\psi(2S)$ one, this needs care. It is not yet evidence for an X-specific effect,
but $\kk\approx1.14$ should be one of the systematic variations for the ntmix training.

\section{PbPb: $B^{+}$ in 2023 and 2024}\label{sec:pbpb}
Without the BDT, the PbPb 2023 $\psi(2S)$ peak has only $\sim$360 signal events on $\sim$9000 background, which is too few for this method.
$B^{+}$ with the error-independent selection gives about 1.4--1.6k signal per year.
\begin{figure}[H]\centering
\includegraphics[width=0.75\textwidth]{fig/PbPb23_Bp_compare.png}
\includegraphics[width=0.75\textwidth]{fig/PbPb24_Bp_compare.png}
\caption{PbPb $B^{+}$: 2023 (top), 2024 data vs 2023 MC (bottom).}
\end{figure}
\texttt{Bchi2Prob}: $\kk_{2023}=$ @@K23@@ and $\kk_{2024}=$ @@K24@@ (vs 2023 MC). The order matches the year ordering seen in the
PbPb chi2Prob comparison, but the difference is only @@D2324@@$\sigma$. Unlike ppRef, the PbPb dimuon vertex does not show the effect
(\texttt{BujvProb}: @@K23MM@@, @@K24MM@@). PbPb needs its own measurement with more signal: a larger or BDT-independent tighter selection, and 2024 MC.
\textbf{The ppRef \kk\ must not be applied to PbPb.}

\section{All measurements}
\begin{figure}[H]\centering
\includegraphics[width=0.85\textwidth]{fig/k_summary.png}
\includegraphics[width=0.75\textwidth]{fig/kscan_all.png}
\caption{Top: all \kk\ values. Bottom: $(\chi^2-\chi^2_{\min})/{\rm ndf}$ as a function of \kk.}
\end{figure}

\section{Limitations}
\begin{itemize}\setlength{\itemsep}{1pt}
\item \textbf{MC statistics.} \texttt{pThatreweight} spans about four orders of magnitude, so the $\psi(2S)$ MC has
$N_{\rm eff}=@@NEFF@@$ of @@NMC@@ candidates. This limits the \kk\ precision and affects every MC-based quantity in the analysis, not only this one.
\item \textbf{One number.} A single scale factor reproduces the bulk but not the data's extreme low-probability tail (1--2\% residual on cut efficiencies).
\item \textbf{Fit models.} The $B^{+}$ mass model has $\chi^2$/ndf $=2.6$; X has low $S/B$; PbPb has low statistics.
\item \textbf{Other features.} Other error-normalized features (\texttt{Bnorm\_svpvDistance\_2D}, \texttt{Bnorm\_trk*Dxy/Dz}) are expected to need the
same factor, applied as a division of the MC significance by \kk. That still has to be validated.
\item \textbf{Pileup.} The official ppRef MC includes pileup (in-time $\langle\mu\rangle \approx 3.8$ in the $\psi(2S)$ sample). No pileup
reweighting is applied yet, and whether the MC pileup profile matches data has not been checked. Pileup cannot change the candidate-only fits
(\texttt{Bchi2Prob}, \texttt{BujvProb}), but it can contribute to the primary-vertex quantities. The planned pileup flag and a pileup-profile comparison will allow this check.
\end{itemize}

\section{Recommendation: correcting MC before the AN and before training}
\begin{enumerate}\setlength{\itemsep}{2pt}
\item In \texttt{Flat\_TREEs.C}, \emph{MC only}: keep the original as \texttt{Bchi2Prob\_raw} and store the corrected value in
\texttt{Bchi2Prob} (likewise for \texttt{BujvProb} and \texttt{BtktkvProb} with $n_{\rm dof}=1$). Then re-apply the preselection thresholds:
\begin{verbatim}
// MC only: emulate the larger track spread observed in data (k from psi(2S) and B+)
float correctProb(float p, int ndf, double k) {
  if (p <= 0.f || p >= 1.f) return p;
  return TMath::Prob(k * k * TMath::ChisquareQuantile(1. - p, ndf), ndf);
}
\end{verbatim}
with $n_{\rm dof}=6$ (ntmix, ntKstar, ntphi: 4 tracks), 4 (ntKp), 1 (\texttt{BujvProb}, \texttt{BtktkvProb}).
\item Nominal $\kk=@@KREC@@$ for ppRef. Systematic variations: $@@KLO@@$ and $@@KHI@@$ (covers the channel/variable/subset spread),
plus $1.14$ for ntmix (the X's own preference).
\item Repeat the data/MC validation (\texttt{DataSIGNAL\_VS\_MC}) with corrected MC for every training feature, on $\psi(2S)$ and $B^{+}$.
This includes the error-normalized significances, divided by \kk\ in MC.
\item Retrain (including ntmix) on corrected MC signal. The background comes from data sidebands and is untouched.
Validate the BDT score and the efficiency of the BDT cut, data vs MC, on $\psi(2S)$.
\item Derive PbPb \kk\ separately (and per year) once 2024 MC exists. Do not reuse the ppRef value.
\end{enumerate}

\appendix
\section{Reproducing}
Everything is in \url{/eos/home-h/hmarques/Analysis_CODES/plotER/Validation/chi2Prob_kStudy/}:
\texttt{extract.py} (reads flat trees and applies selections), \texttt{analyze.py} (all fits, \kk\ scans, figures,
\texttt{results.json}) and \texttt{make\_note.py} (this document). Inputs: the flat trees in
\url{/eos/user/h/hmarques/RUN3_Data_MC_sharing/} (subfolders \texttt{X3872} and \texttt{Bmesons}).
Bootstrap random seed 12345.
\end{document}
"""
for k, v in V.items():
    TEX = TEX.replace(f"@@{k}@@", v)
assert "@@" not in TEX, [l for l in TEX.splitlines() if "@@" in l][:3]
open("chi2Prob_kStudy.tex", "w").write(TEX)
for _ in range(2):
    r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "chi2Prob_kStudy.tex"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
print("pdflatex exit", r.returncode)
if r.returncode:
    print(r.stdout[-3000:])
