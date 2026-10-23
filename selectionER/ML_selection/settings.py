"""Settings of the ML selection: the only file to edit by hand.

Every tool reads it: the skims (scripts/make_skims.py), the NN and XGB trainings and the figures of
merit (methods/ML_common), the feature study (scripts/feature_study.py), the scans
(scripts/fresh_ML_scan.sh, scripts/train_best_ML.sh) and correlations/correlation_dataMC.py. TMVA
and optimalCUT_X_punzi.C read its C++ copy, methods/ML_tmva/ml_settings.h, which make_skims.py,
feature_study.py, fresh_ML_scan.sh and train_best_ML.sh write again from this file when they start
(ml_samples.write_cpp_settings).

A case is one candidate tree in one data set, named "<tree>/<year>": ntmix/PbPb23, ntKp/ppRef, ...
Its skims, feature studies and scored samples live in <folder>/<tree>/<year>/. The settings below
are keyed by tree (they hold for every year of it) or, in CASES, by tree and year. Only the cases
listed in CASES exist.
"""

FLAT = "/eos/user/h/hmarques/RUN3_Data_MC_sharing"
# plotER/aux/masses.h
X3872_MASS = 3.87164
PSI2S_MASS = 3.68610
BU_MASS = 5.27941
BD_MASS = 5.27963
BS_MASS = 5.36691

# =============================================================================
# Per tree
# =============================================================================

# The inputs of a case: key -> (tree in the flat file, kind, label). The kind tags the skim and
# scored files. data: the background (DATA in the training sidebands) and every DATA check;
# signal: the MC the classifier is trained on (pThat-weighted); every other key is a spectator,
# an MC sample whose efficiency at the working point is reported.
INPUTS = {
    "ntmix": {
        "data": ("ntmix", "DATA", "DATA"),
        "signal": ("ntmix_X3872", "MC_X3872", "prompt X(3872) MC"),
        "signal_nonprompt": ("ntmix_X3872", "MC_X3872_NONPROMPT", "nonprompt X(3872) MC"),
        "spectator": ("ntmix_PSI2S", "MC_PSI2S", "prompt psi(2S) MC"),
        "spectator_nonprompt": ("ntmix_PSI2S", "MC_PSI2S_NONPROMPT", "nonprompt psi(2S) MC"),
    },
    "ntKp": {"data": ("ntKp", "DATA", "DATA"), "signal": ("ntKp", "MC", "B+ -> J/psi K+ MC")},
    "ntphi": {"data": ("ntphi", "DATA", "DATA"), "signal": ("ntphi", "MC", "Bs -> J/psi phi MC")},
    "ntKstar": {"data": ("ntKstar", "DATA", "DATA"), "signal": ("ntKstar", "MC", "B0 -> J/psi K*0 MC")},
}

# Production classifier inputs, one line per tree (feature_study.py adopt rewrites the line). A name
# ending in _dc is the mass-decorrelated copy of the feature, a branch of the skims; a name ending in
# _log is log10 of the flat branch (ml_samples.add_log_features). The B mesons use one set of inputs, the
# production features and the feature-study pool alike: the pointing angle, the vertex probability, the
# decay-length and track-1 impact-parameter significances, track-1 pT and Delta R, and the J/psi eta.
# Balpha, the decay-length significance and the track pT span 1-4 decades (Balpha: signal 0.0005-2.4
# rad, DATA up to pi): the NN takes them as log10, which also keeps the decorrelation's shifts of the
# tails proportional; XGBoost bins every input by its quantiles, so log10 changes nothing for it.
# Bnorm_trk1Dxy is signed (symmetric), so it stays as it is.
FEATURES = {
    "ntmix": ["Btrk1dR", "Btrk2Pt", "Bmu1pt", "Bmu2pt", "BtrkPtimb", "Bmu1eta", "Bmu2eta", "Bujeta", "Btktketa"],
    "ntKp": ["Balpha_log", "Bchi2Prob", "Bnorm_svpvDistance_2D_log", "Bnorm_trk1Dxy", "Btrk1Pt_log", "Bujeta", "Btrk1dR"],
    "ntphi": ["Balpha_log", "Bchi2Prob", "Bnorm_svpvDistance_2D_log", "Bnorm_trk1Dxy", "Btrk1Pt_log", "Bujeta", "Btrk1dR"],
    "ntKstar": ["Balpha_log", "Bchi2Prob", "Bnorm_svpvDistance_2D_log", "Bnorm_trk1Dxy", "Btrk1Pt_log", "Bujeta", "Btrk1dR"],
}

# Feature study: the features it selects from (fullPOOL/all, then fullPOOL/trimmed with the most
# important ones). The B mesons take their production set (FEATURES) with, for the two-track Bs and B0,
# the pT-leading track pT and the largest track Delta R instead of those of track 1; B+ has one track.
POOL = {
    "ntmix": ["BtrkLeadPt", "BtrkSubPt", "BtrkMaxdR", "BmuLeadPt", "Bpt", "BmuMaxdR"],
    "ntKp": ["Balpha_log", "Bchi2Prob", "Bnorm_svpvDistance_2D_log", "Bnorm_trk1Dxy", "Btrk1Pt_log", "Bujeta", "Btrk1dR"],
    "ntphi": ["Balpha_log", "Bchi2Prob", "Bnorm_svpvDistance_2D_log", "Bnorm_trk1Dxy", "BtrkLeadPt_log", "Bujeta", "BtrkMaxdR"],
    "ntKstar": ["Balpha_log", "Bchi2Prob", "Bnorm_svpvDistance_2D_log", "Bnorm_trk1Dxy", "BtrkLeadPt_log", "Bujeta", "BtrkMaxdR"],
}

# The classifiers of each tree (feature_study.py runs a study for each). No tree uses TMVA at present.
ML_METHODS = {"ntmix": ("xgb", "nn"), "ntKp": ("xgb", "nn"), "ntphi": ("xgb", "nn"), "ntKstar": ("xgb", "nn")}

# Feature study: the forward selection (one input added per step) after fullPOOL/all and trimmed.
FORWARD_SELECTION = {"ntmix": False, "ntKp": False, "ntphi": False, "ntKstar": False}

# Feature study: the final check, the last selection trained again with the other decorrelation setting
# (DECORRELATE).
FINAL_CHECK = {"ntmix": False, "ntKp": False, "ntphi": False, "ntKstar": False}

# Feature study: per method, the inputs it trains on: False the raw features, True their _dc copies
# (methods/ML_common/ml_decorrelation.py). With raw inputs the NN learns the mass (X(3872)
# fullPOOL/all: R_gap 14-26, psi(2S) efficiency 0.1-0.5% at the working point); the final check of
# every study trains the last selection with the other setting.
DECORRELATE = {
    "ntmix": {"xgb": False, "nn": True, "tmva": False},
    "ntKp": {"xgb": False, "nn": True, "tmva": False},
    "ntphi": {"xgb": False, "nn": True, "tmva": False},
    "ntKstar": {"xgb": False, "nn": True, "tmva": False},
}

# Largest max_depth of the XGBoost scans (XGB_optuna.py draws 2 .. this). The B mesons, with far more
# signal MC and sideband DATA than the X(3872), went to the old cap of 5 in most candidates.
# The range of min_child_weight (the smallest hessian sum of a leaf, log-uniform) of the scans. X(3872)
# at pT > 10: 1% of the signal MC carries half of its weight (training weight ~130-180 each, median 0.04),
# so one candidate gives a leaf a hessian up to ~40; with depth up to 5 and min_child_weight up to 20, 96%
# of the trials failed the overtraining test (train - validation AUC gap ~0.05). Depth 3 and a floor of 20
# keep the trees from isolating a few heavy candidates.
XGB_MAX_DEPTH = {"ntmix": 3, "ntKp": 8, "ntphi": 8, "ntKstar": 8}
XGB_MIN_CHILD_WEIGHT = {"ntmix": (20.0, 500.0), "ntKp": (0.5, 20.0), "ntphi": (0.5, 20.0), "ntKstar": (0.5, 20.0)}

# The mass regions of each tree, in GeV:
#   peak                    the signal mass
#   punzi_half_width        Punzi window |Bmass - peak| < punzi_half_width
#   punzi_sideband_start,   its background: the DATA from start to start + width off the peak, on the
#   punzi_sideband_width,   sides of punzi_sides (-1 below the peak, +1 above), scaled to the window
#   punzi_sides             (the background is counted on all DATA, as optimalCUT_X_punzi.C does)
#   training_sidebands      the background candidates of the training
#   decorrelation_range,    the quantile maps of the _dc features are fitted on DATA in the range
#   decorrelation_excluded, outside the excluded windows, in slices of decorrelation_slice
#   decorrelation_slice     (the range is also that of the DATA mass plots)
#   sculpting_range,        the sculpting line is fitted to the background efficiency in bins of
#   sculpting_excluded,     sculpting_bin over the range outside the excluded window; R_gap is the
#   sculpting_gap,          efficiency from sculpting_gap off the peak to the excluded edges over the
#   sculpting_bin           line's prediction
#   marks                   masses drawn on the DATA mass plots
# X(3872): the MC peak holds 68% within 6 MeV, so +-5 MeV keeps ~60% of it and 15-35 MeV lets only
# ~2% of the window signal into the sidebands. Its training sidebands span 20-100 MeV from the peak. B mesons: the MC peak holds 68% within 20 MeV (B+, B0;
# Bs 14 MeV): the window is about 2 sigma. The training sidebands span 65-235 MeV from the peak: wide, for
# a smooth background with many candidates, and close, for a background like the one under the peak
# (after a tight cut the B+ background falls ~10x from the peak to 200-300 MeV) that leaves less room to
# learn the mass (R_gap 1.43 +- 0.12, Bs XGB with 200-300 MeV); the signal MC there (B+ 1.4%, Bs 0.6%,
# B0 3.7%) is negligible next to the DATA before any cut. The Punzi background, counted after the cut
# when little background is left, comes from 100-200 MeV, where the signal MC is 0.27% (B+), 0.08% (Bs)
# and 1.16% (B0, its broad tail). The decorrelation fit and the sculpting line exclude the gap between
# the training sidebands, as for the X(3872); R_gap is measured from sculpting_gap to it. B+: only the
# high side, for the training and for the Punzi background, because the low one holds partially
# reconstructed B decays (lost particles land below the peak).
MASS = {
    "ntmix": {
        "peak": X3872_MASS,
        "punzi_half_width": 0.005, "punzi_sideband_start": 0.015, "punzi_sideband_width": 0.02, "punzi_sides": (-1, 1),
        "training_sidebands": ((X3872_MASS - 0.100, X3872_MASS - 0.020), (X3872_MASS + 0.020, X3872_MASS + 0.100)),
        "decorrelation_range": (3.60, 4.00),
        "decorrelation_excluded": ((PSI2S_MASS - 0.020, PSI2S_MASS + 0.020), (3.85, 3.89)),
        "decorrelation_slice": 0.010,
        "sculpting_range": (3.75, 4.00), "sculpting_excluded": (3.85, 3.89), "sculpting_gap": 0.011, "sculpting_bin": 0.005,
        "marks": {"X(3872)": X3872_MASS, "$\\psi$(2S)": PSI2S_MASS},
    },
    "ntKp": {
        "peak": BU_MASS,
        "punzi_half_width": 0.040, "punzi_sideband_start": 0.10, "punzi_sideband_width": 0.10, "punzi_sides": (1,),
        "training_sidebands": ((BU_MASS + 0.065, BU_MASS + 0.235),),
        "decorrelation_range": (BU_MASS - 0.35, BU_MASS + 0.35),
        "decorrelation_excluded": ((BU_MASS - 0.065, BU_MASS + 0.065),),
        "decorrelation_slice": 0.020,
        "sculpting_range": (BU_MASS - 0.35, BU_MASS + 0.35), "sculpting_excluded": (BU_MASS - 0.065, BU_MASS + 0.065),
        "sculpting_gap": 0.050, "sculpting_bin": 0.010,
        "marks": {"B$^+$": BU_MASS},
    },
    "ntphi": {
        "peak": BS_MASS,
        "punzi_half_width": 0.025, "punzi_sideband_start": 0.10, "punzi_sideband_width": 0.10, "punzi_sides": (-1, 1),
        "training_sidebands": ((BS_MASS - 0.235, BS_MASS - 0.065), (BS_MASS + 0.065, BS_MASS + 0.235)),
        "decorrelation_range": (BS_MASS - 0.35, BS_MASS + 0.35),
        "decorrelation_excluded": ((BS_MASS - 0.065, BS_MASS + 0.065),),
        "decorrelation_slice": 0.020,
        "sculpting_range": (BS_MASS - 0.35, BS_MASS + 0.35), "sculpting_excluded": (BS_MASS - 0.065, BS_MASS + 0.065),
        "sculpting_gap": 0.035, "sculpting_bin": 0.010,
        "marks": {"B$_s^0$": BS_MASS},
    },
    "ntKstar": {
        "peak": BD_MASS,
        "punzi_half_width": 0.040, "punzi_sideband_start": 0.10, "punzi_sideband_width": 0.10, "punzi_sides": (-1, 1),
        "training_sidebands": ((BD_MASS - 0.235, BD_MASS - 0.065), (BD_MASS + 0.065, BD_MASS + 0.235)),
        "decorrelation_range": (BD_MASS - 0.35, BD_MASS + 0.35),
        "decorrelation_excluded": ((BD_MASS - 0.065, BD_MASS + 0.065),),
        "decorrelation_slice": 0.020,
        "sculpting_range": (BD_MASS - 0.35, BD_MASS + 0.35), "sculpting_excluded": (BD_MASS - 0.065, BD_MASS + 0.065),
        "sculpting_gap": 0.050, "sculpting_bin": 0.010,
        "marks": {"B$^0$": BD_MASS},
    },
}

# DATA in the skims: every candidate in the training sidebands (which hold the Punzi sidebands), and 1 in
# DATA_THINNING of the others inside the decorrelation range (the only DATA the selection uses). Training and
# the Punzi background use only the sidebands; the others serve the sculpting check, the mass plots and the
# decorrelation check, which count each of them DATA_THINNING times (ml_fom.data_weights). The decorrelation
# maps are fitted on all of them. The B0 DATA holds ~55M candidates per year in its +-350 MeV window, ~16M of
# them in the sidebands.
DATA_THINNING = {"ntmix": 1, "ntKp": 10, "ntphi": 10, "ntKstar": 10}

# =============================================================================
# Per case
# =============================================================================

# files: the flat file of every input of the tree (INPUTS). PbPb24 has no MC of its own: its DATA is
# paired with the PbPb23 MC. ppRef: only B+ so far (its Bs and B0 MC have no pThatreweight).
# ntmix Run 3: three regions, <year>_<region> (X_REGIONS), with the 2025 pre-cut (X_PRE_CUT_2025, no
# Bchi2Prob cut beyond the flattening's 0.005), raw inputs for every method and looser XGBoost limits (X_XGB_2025):
#   PbPb25            PbPb25 DATA, PbPb25 private MC; inputs Btrk1dR Btrk2Pt Btrk1Pt Bmu1pt
#   <year>_mc2325     PbPb23, PbPb24, PbPb25 DATA, each with the PbPb23 official + PbPb25 private MC merged
#                     (X3872/PbPb2325: pThatreweight 1/(1/w23 + 1/w25) per pThat range, the luminosities
#                     add); inputs BtrkLeadPt BtrkSubPt BtrkMaxdR BmuLeadPt
# pre_cut: ROOT syntax (&& for AND, || for OR), the candidates every training, check and scored-sample
# cut starts from. The skims carry it: after a change, python3 scripts/make_skims.py --tree <tree>
# --year <year>. PbPb18: the legacy flat files hold only Bchi2Prob > 0.10; the cut is kept explicit.
# pt_bins: the pT bins of optimalCUT_X_punzi.C, inside the pT range of the pre-cut (the inclusive
# bin is always added).
# A case may also set any per-tree setting (features, pool, decorrelate, ...) for itself.
X_PbPb18 = f"{FLAT}/X3872/legacy_run2/REflated/flat_ntmix_PbPb18"
X_PbPb25_MC = f"{FLAT}/X3872/PbPb25/flat_ntmix_PbPb25_MC"
X_PbPb25_MC_FILES = {"signal": f"{X_PbPb25_MC}_X3872.root", "signal_nonprompt": f"{X_PbPb25_MC}_X3872_nonPrompt.root",
                     "spectator": f"{X_PbPb25_MC}_PSI2S.root", "spectator_nonprompt": f"{X_PbPb25_MC}_PSI2S_nonPrompt.root"}
X_PbPb2325_MC = f"{FLAT}/X3872/PbPb2325/flat_ntmix_PbPb2325_MC"
X_PbPb2325_MC_FILES = {"signal": f"{X_PbPb2325_MC}_X3872.root", "signal_nonprompt": f"{X_PbPb2325_MC}_X3872_nonPrompt.root",
                       "spectator": f"{X_PbPb2325_MC}_PSI2S.root", "spectator_nonprompt": f"{X_PbPb2325_MC}_PSI2S_nonPrompt.root"}
X_REGIONS = {"pt10y24": ("((Bpt > 10) && (Bpt < 50)) && (abs(By) < 2.4)", [10.0, 15.0, 50.0]),
             "pt15y16": ("((Bpt > 15) && (Bpt < 50)) && (abs(By) < 1.6)", [15.0, 50.0]),
             "pt10y16": ("((Bpt > 10) && (Bpt < 50)) && (abs(By) < 1.6)", [10.0, 15.0, 50.0])}
X_PRE_CUT_2025 = "(BQvalue < 0.15) && (Btrk2dR < 0.3)"
# PbPb25 MC: N_eff 13% of the candidates (PbPb23 official: 3%), so the X(3872) caps are loosened a bit.
X_XGB_2025 = {"xgb_max_depth": 4, "xgb_min_child_weight": (5.0, 200.0)}
X_RAW = {"decorrelate": {"xgb": False, "nn": False, "tmva": False}}
X_INPUTS_PbPb25 = ["Btrk1dR", "Btrk2Pt", "Btrk1Pt", "Bmu1pt"]
X_INPUTS_MERGED = ["BtrkLeadPt", "BtrkSubPt", "BtrkMaxdR", "BmuLeadPt"]
B_PRE_CUT = "((Bpt > 7.5) && (Bpt < 50)) && (abs(By) < 2.4) && (Bchi2Prob > 0.05)"
B_PT_BINS = [7.5, 10.0, 15.0, 20.0, 30.0, 50.0]
CASES = {
    "ntmix": {
        "PbPb18": {
            "files": {"data": f"{X_PbPb18}_DATA.root", "signal": f"{X_PbPb18}_MC_X3872.root",
                      "signal_nonprompt": f"{X_PbPb18}_MC_X3872_nonPrompt.root", "spectator": f"{X_PbPb18}_MC_PSI2S.root",
                      "spectator_nonprompt": f"{X_PbPb18}_MC_PSI2S_nonPrompt.root"},
            "pre_cut": "((Bpt > 15) && (Bpt < 50)) && (abs(By) < 1.6) && (BQvalue < 0.15) && (Bchi2Prob > 0.10)",
            "pt_bins": [15.0, 50.0],
        },
        **{f"PbPb25_{region}": {"files": {"data": f"{FLAT}/X3872/PbPb25/flat_ntmix_PbPb25_DATA.root", **X_PbPb25_MC_FILES},
                                "pre_cut": f"{cut} && {X_PRE_CUT_2025}", "pt_bins": pt_bins,
                                "features": X_INPUTS_PbPb25, "pool": X_INPUTS_PbPb25, **X_RAW, **X_XGB_2025}
           for region, (cut, pt_bins) in X_REGIONS.items()},
        **{f"{year}_mc2325_{region}": {"files": {"data": f"{FLAT}/X3872/{year}/flat_ntmix_{year}_DATA.root", **X_PbPb2325_MC_FILES},
                                       "pre_cut": f"{cut} && {X_PRE_CUT_2025}", "pt_bins": pt_bins,
                                       "features": X_INPUTS_MERGED, "pool": X_INPUTS_MERGED, **X_RAW, **X_XGB_2025}
           for year in ("PbPb23", "PbPb24", "PbPb25") for region, (cut, pt_bins) in X_REGIONS.items()},
        # ppRef, when its private X(3872) MC is forested and flattened (MiniAOD done at T2_US_Vanderbilt,
        # 2026-10-22): the files of X3872/ppRef24/flat_ntmix_ppRef_*.root,
        # with "pt_bins": [7.5, 12.5, 17.5, 22.5, 50.0].
    },
    "ntKp": {
        "PbPb23": {"files": {"data": f"{FLAT}/Bmesons/PbPb23/flat_ntKp_PbPb23_DATA.root",
                             "signal": f"{FLAT}/Bmesons/PbPb23/flat_ntKp_PbPb23_MC.root"},
                   "pre_cut": B_PRE_CUT, "pt_bins": B_PT_BINS},
        "PbPb24": {"files": {"data": f"{FLAT}/Bmesons/PbPb24/flat_ntKp_PbPb24_DATA.root",
                             "signal": f"{FLAT}/Bmesons/PbPb23/flat_ntKp_PbPb23_MC.root"},
                   "pre_cut": B_PRE_CUT, "pt_bins": B_PT_BINS},
        "ppRef": {"files": {"data": f"{FLAT}/Bmesons/ppRef/flat_ntKp_ppRef_DATA.root",
                            "signal": f"{FLAT}/Bmesons/ppRef/flat_ntKp_ppRef_MC.root"},
                  "pre_cut": B_PRE_CUT, "pt_bins": B_PT_BINS},
    },
    "ntphi": {
        "PbPb23": {"files": {"data": f"{FLAT}/Bmesons/PbPb23/flat_ntphi_PbPb23_DATA.root",
                             "signal": f"{FLAT}/Bmesons/PbPb23/flat_ntphi_PbPb23_MC.root"},
                   "pre_cut": B_PRE_CUT, "pt_bins": B_PT_BINS},
        "PbPb24": {"files": {"data": f"{FLAT}/Bmesons/PbPb24/flat_ntphi_PbPb24_DATA.root",
                             "signal": f"{FLAT}/Bmesons/PbPb23/flat_ntphi_PbPb23_MC.root"},
                   "pre_cut": B_PRE_CUT, "pt_bins": B_PT_BINS},
    },
    "ntKstar": {
        "PbPb23": {"files": {"data": f"{FLAT}/Bmesons/PbPb23/flat_ntKstar_PbPb23_DATA.root",
                             "signal": f"{FLAT}/Bmesons/PbPb23/flat_ntKstar_PbPb23_MC.root"},
                   "pre_cut": B_PRE_CUT, "pt_bins": B_PT_BINS},
        "PbPb24": {"files": {"data": f"{FLAT}/Bmesons/PbPb24/flat_ntKstar_PbPb24_DATA.root",
                             "signal": f"{FLAT}/Bmesons/PbPb23/flat_ntKstar_PbPb23_MC.root"},
                   "pre_cut": B_PRE_CUT, "pt_bins": B_PT_BINS},
    },
}

# =============================================================================
# Shared by every case
# =============================================================================

MC_WEIGHT_BRANCH = "pThatreweight"

# NN only: features that get log1p before the standardization, by flat-tree name (X and X_dc alike);
# a name the model does not use is ignored. The B mesons take their _log inputs instead.
LOG_FEATURES = ["Btrk1Pt", "Btrk2Pt", "Bmu1pt", "Bmu2pt", "BtrkLeadPt", "BtrkSubPt", "BmuLeadPt", "Bujpt", "Bpt"]

# Maximum candidates used for training (signal MC, DATA in the training sidebands). MAX_SIGNAL is above
# the signal MC after every pre-cut (X(3872) PbPb23 118k, PbPb18 184k; B+ PbPb23 279k, ppRef 298k): all of
# it is used.
MAX_SIGNAL = 300000
MAX_BACKGROUND = 500000

# Model selection (methods/ML_common/ml_fom.py, methods/ML_tmva/TMVA_common.h). FOM: the validation
# figure of merit a scan maximizes by default (fresh_ML_scan.sh):
#   auc     weighted ROC AUC
#   sigeff  signal efficiency at background efficiency BKG_EFF_TARGET
#   pauc    ROC area for background efficiency 0 to PAUC_MAX_BKG_EFF, divided by it (0 to 1)
#   punzi   best Punzi FOM e_S / S_min(B), with a = PUNZI_A, b = PUNZI_B and B the background in the
#           Punzi window estimated from its sidebands (MASS), as optimalCUT_X_punzi.C does
FOM = "auc"
BKG_EFF_TARGET = 0.02
PAUC_MAX_BKG_EFF = 0.05
PUNZI_A = 5.0
PUNZI_B = 2.0
# The Punzi maximum is searched only among cuts that keep at least this many sideband candidates: with
# fewer, B rests on a handful of candidates and an empty tail wins by chance.
PUNZI_MIN_BACKGROUND = 10.0
# The Punzi FOM is scanned on thresholds k * THRESHOLD_STEP (the same grid in the model selection and in
# optimalCUT_X_punzi.C), so the maximum does not land on a fluctuation of the curve. The B-meson scores
# crowd against 1 (~60% of the signal above 0.999), so they take finer steps; the X(3872) keeps the 0.02 of
# its finished studies.
THRESHOLD_STEP = {"ntmix": 0.02, "ntKp": 0.001, "ntphi": 0.001, "ntKstar": 0.001}
# Overtraining pass/fail (train vs validation): both KS p-values >= KS_PVALUE_MIN and AUC gap <= the
# larger of MAX_AUC_GAP and AUC_GAP_SIGMA times its statistical uncertainty (Hanley-McNeil with effective
# entries). A failing trial is never chosen. The X(3872) PbPb23 signal MC has ~3.4k effective entries of
# 118k (the pThat-5 and pThat-10 campaigns carry ~70% of the weight): a model fitted on the test split
# alone showed a train - validation gap of 0.016 +- 0.012, so a fixed 0.01 rejected every trial.
# AUC_GAP_SIGMA is per tree. X(3872) at pT > 10: ~300 validation signal candidates with weight > 50 carry
# half of its weight (N_eff ~650-950), and the split, the same for every trial, already moves the AUC of a
# single untrained input by up to 0.025 between training and validation; at 2 sigma 95% of the XGB
# trials failed (gap ~3.5 sigma for depth 2 and min_child_weight 500 alike), at 3 sigma 20% pass.
KS_PVALUE_MIN = 0.01
MAX_AUC_GAP = 0.01
AUC_GAP_SIGMA = {"ntmix": 3.0, "ntKp": 2.0, "ntphi": 2.0, "ntKstar": 2.0}
