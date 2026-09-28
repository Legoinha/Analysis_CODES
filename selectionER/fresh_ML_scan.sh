#!/bin/bash
# =============================================================================
# Central control for the ML hyperparameter scans (ML_pytorch, ML_xgboost, ML_tmva).
#
# Edit ONLY the settings block below. The script then
#   1. writes the features, cuts, sidebands and sample sizes into all three methods
#      (ML_common/ml_samples.py for NN and XGB, ML_tmva/TMVA_config.h for TMVA),
#   2. checks them against the real input files: every branch exists, the pre-cuts parse,
#      ROOT and the NumPy translation used by NN/XGB select the same candidates, and the
#      skims (make_skims.py), which every training reads, match the current pre-cut,
#   3. asks for confirmation, removes old scan jobs and deletes the old outputs
#      of the selected methods and samples,
#   4. submits the scans: NN 200, XGB 250, TMVA 384 jobs per sample.
#
# Every scan run is tagged with its model-selection FOM (auc, sigeff, pauc or punzi):
# outputs, logs and scored samples of pbpb23 with FOM auc live in *_outputs/.../pbpb23_auc,
# so runs with different FOMs are kept side by side and cleaned separately.
# The feature study (feature_study.py) uses the same settings; its adopt command writes
# the chosen set into FEATURES below.
#
# Usage, from any directory on lxplus:
#   ./fresh_ML_scan.sh              full chain for all three methods, FOM from the settings
#   ./fresh_ML_scan.sh sigeff       the same with another FOM (auc, sigeff, pauc, punzi)
#   ./fresh_ML_scan.sh tmva         clean and resubmit only the listed methods (nn, xgb, tmva)
#   ./fresh_ML_scan.sh PbPb18       clean and resubmit only the listed samples
#   ./fresh_ML_scan.sh PbPb18 tmva  these can be combined; everything else is left alone
#   ./fresh_ML_scan.sh --sync-only  steps 1 and 2 only; nothing deleted or submitted
# =============================================================================

# ----------------------------- settings --------------------------------------
SAMPLES=(PbPb23 PbPb18)

# Production classifier inputs (all three methods). A name ending in _dc is the
# mass-decorrelated copy of the feature, a branch of the skims.
FEATURES=(Btrk1dR Btrk2Pt Bmu1pt Bmu2pt BtrkPtimb Bmu1eta Bmu2eta Bujeta Btktketa)

# Feature study (feature_study.py): the features it selects from.
FEATURE_POOL=(BtrkLeadPt BtrkMaxdR BtrkSubPt BtrkMindR BtrkMaxAbsEta BmuLeadPt BmuMaxdR)

# Features that get log1p before standardization, NN only. By base name: X and X_dc alike.
LOG_FEATURES=(Btrk2Pt Bmu1pt Bmu2pt BtrkLeadPt BtrkSubPt BmuLeadPt)

# Pre-cut of each sample in ROOT syntax: && for AND, || for OR, case-sensitive branch names.
# The skims carry it: after a change, remake them with python3 make_skims.py --sample <sample>.
PRE_CUT_PbPb23='((Bpt > 15) && (Bpt < 50)) && (abs(By) < 2.4) && (BQvalue < 0.15) && (Bchi2Prob > 0.05)'
PRE_CUT_PbPb18='((Bpt > 15) && (Bpt < 50)) && (abs(By) < 2.4) && (BQvalue < 0.15) && (Bchi2Prob > 0.10)'
# Training background: DATA in these Bmass windows, GeV: 15-50 MeV from m_X = 3.87164 on both sides.
# They overlap the Punzi sidebands (15-35 MeV): the Punzi background is counted on all DATA there,
# the training candidates included, as optimalCUT_X_punzi.C does.
TRAINING_SIDEBANDS="3.82164-3.85664 3.88664-3.92164"

# Maximum candidates used for training (signal MC, sideband DATA).
MAX_SIGNAL=165000
MAX_BACKGROUND=500000

# Model selection, shared by the three methods (ML_common/ml_fom.py, ML_tmva/TMVA_config.h).
# FOM: the validation figure of merit each scan maximizes:
#   auc     weighted ROC AUC
#   sigeff  signal efficiency at background efficiency BKG_EFF_TARGET
#   pauc    ROC area for background efficiency 0 to PAUC_MAX_BKG_EFF, divided by it (0 to 1)
#   punzi   best Punzi FOM, with the window background estimated from the DATA sidebands
#           next to it (same window as optimalCUT_X_punzi.C) times the background efficiency
FOM=auc
BKG_EFF_TARGET=0.02
PAUC_MAX_BKG_EFF=0.05
PUNZI_A=5.0
PUNZI_B=2.0
# Punzi window |Bmass - m_X| < SIGNAL_HALF_WIDTH, and its background from the sidebands
# SIDEBAND_START to SIDEBAND_START + SIDEBAND_WIDTH on both sides, scaled to the window (GeV).
# The signal MC peak holds 68% within 6 MeV: +-5 MeV keeps ~60% of it; 15-35 MeV lets only ~2% of
# the window signal into the background estimate. The Punzi sidebands lie inside the training
# sidebands; their background is counted on all DATA (see TRAINING_SIDEBANDS).
SIGNAL_HALF_WIDTH=0.005
SIDEBAND_START=0.015
SIDEBAND_WIDTH=0.020
# The Punzi maximum is searched only among cuts that keep at least this many background
# candidates (validation split in the model selection, DATA sidebands in optimalCUT_X_punzi.C):
# with fewer, B is estimated from a handful of candidates and an empty tail wins by chance.
PUNZI_MIN_BACKGROUND=10
# The Punzi FOM is scanned on thresholds k * THRESHOLD_STEP (a coarse grid, the same in the model
# selection and optimalCUT_X_punzi.C), so the maximum does not land on a fluctuation of the curve.
THRESHOLD_STEP=0.02
# Overtraining pass/fail (train vs validation): both KS p-values >= KS_PVALUE_MIN and
# AUC gap <= MAX_AUC_GAP. A failing trial is rejected, never chosen as the best one.
KS_PVALUE_MIN=0.05
MAX_AUC_GAP=0.01

# Jobs of earlier scans are found and removed automatically: every queued job
# whose log file lives in one of the three ML condor_logs folders.
# -----------------------------------------------------------------------------

set -eo pipefail
BASE=/eos/user/h/hmarques/Analysis_CODES/selectionER
cd "$BASE"
MODE=full
METHODS=()
CLI_SAMPLES=()
for arg in "$@"; do
  case "$arg" in
    --sync-only) MODE=--sync-only ;;
    nn|xgb|tmva) METHODS+=("$arg") ;;
    auc|sigeff|pauc|punzi) FOM="$arg" ;;
    PbPb18|PbPb23) CLI_SAMPLES+=("$arg") ;;
    *) echo "Unknown option: $arg (use nn, xgb, tmva, auc, sigeff, pauc, punzi, PbPb18, PbPb23, or --sync-only)"; exit 1 ;;
  esac
done
[ ${#METHODS[@]} -eq 0 ] && METHODS=(nn xgb tmva)
[ ${#CLI_SAMPLES[@]} -gt 0 ] && SAMPLES=("${CLI_SAMPLES[@]}")
has() { [[ " ${METHODS[*]} " == *" $1 "* ]]; }

# ---- 1. write the settings into the three methods ----------------------------
echo "=== 1. Writing settings into ML_common, ML_pytorch, ML_xgboost, ML_tmva"
BACKUP_DIR="$BASE/.ML_sync_backups/$(date +%Y%m%d_%H%M%S)"
export BACKUP_DIR PRE_CUT_PbPb23 PRE_CUT_PbPb18 TRAINING_SIDEBANDS MAX_SIGNAL MAX_BACKGROUND
export BKG_EFF_TARGET PAUC_MAX_BKG_EFF PUNZI_A PUNZI_B SIGNAL_HALF_WIDTH SIDEBAND_START SIDEBAND_WIDTH PUNZI_MIN_BACKGROUND THRESHOLD_STEP KS_PVALUE_MIN MAX_AUC_GAP
export FEATURES_LIST="${FEATURES[*]}" LOG_FEATURES_LIST="${LOG_FEATURES[*]}"
export POOL_LIST="${FEATURE_POOL[*]}"
python3 - <<'PY'
import os, re, shutil, sys

feats = os.environ["FEATURES_LIST"].split()
logf = os.environ["LOG_FEATURES_LIST"].split()
pool = os.environ["POOL_LIST"].split()
pre_cuts = {s: os.environ["PRE_CUT_" + s] for s in ("PbPb18", "PbPb23")}
sidebands = [tuple(float(x) for x in w.split("-")) for w in os.environ["TRAINING_SIDEBANDS"].split()]
nsig, nbkg = int(os.environ["MAX_SIGNAL"]), int(os.environ["MAX_BACKGROUND"])
backup = os.environ["BACKUP_DIR"]

for name, items in (("FEATURES", feats), ("FEATURE_POOL", pool)):
    if len(set(items)) != len(items): sys.exit(name + " has duplicates")
for cut in pre_cuts.values():
    if re.search(r"(?<![&|])[&|](?![&|])", cut): sys.exit("Use && and || in cuts, not single & or |: " + cut)
# The same string as ml_samples.sideband_cut().
side = " || ".join("((Bmass > %.5f) && (Bmass < %.5f))" % w for w in sidebands)

def pylist(name, items):
    return name + " = [\n" + "".join("    %r,\n" % i for i in items) + "]"

M, S = re.M, re.M | re.S
fom_names = ("BKG_EFF_TARGET", "PAUC_MAX_BKG_EFF", "PUNZI_A", "PUNZI_B", "SIGNAL_HALF_WIDTH", "SIDEBAND_START",
             "SIDEBAND_WIDTH", "PUNZI_MIN_BACKGROUND", "THRESHOLD_STEP", "KS_PVALUE_MIN", "MAX_AUC_GAP")
fom = {name: repr(float(os.environ[name])) for name in fom_names}
window = ("SIGNAL_HALF_WIDTH", "SIDEBAND_START", "SIDEBAND_WIDTH", "PUNZI_MIN_BACKGROUND", "THRESHOLD_STEP")
edits = {
    "ML_common/ml_fom.py": [(r"^%s = .*$" % n, "%s = %s" % (n, fom[n]), M) for n in fom_names],
    "ML_common/ml_samples.py": [
        (r"^PRE_CUTS = \{\n.*?^\}", "PRE_CUTS = {\n" + "".join("    %r: %r,\n" % kv for kv in pre_cuts.items()) + "}", S),
        (r"^TRAINING_SIDEBANDS = .*$", "TRAINING_SIDEBANDS = (%s)" % ", ".join("(%.5f, %.5f)" % w for w in sidebands), M),
        (r"^FEATURES = \[\n.*?^\]", pylist("FEATURES", feats), S),
        (r"^FEATURE_POOL = \[\n.*?^\]", pylist("FEATURE_POOL", pool), S),
    ],
    "optimalCUT_X_punzi.C": [(r"^static const double %s = .*;$" % n, "static const double %s = %s;" % (n, fom[n]), M)
                             for n in window],
    "ML_pytorch/NN_train.py": [
        (r"^LOG_FEATURES = \[\n.*?^\]", pylist("LOG_FEATURES", logf), S),
        (r"^MAX_SIGNAL = .*$", "MAX_SIGNAL = %d" % nsig, M),
        (r"^MAX_BACKGROUND = .*$", "MAX_BACKGROUND = %d" % nbkg, M),
    ],
    "ML_xgboost/XGB_train.py": [
        (r"^MAX_SIGNAL = .*$", "MAX_SIGNAL = %d" % nsig, M),
        (r"^MAX_BACKGROUND = .*$", "MAX_BACKGROUND = %d" % nbkg, M),
    ],
    "ML_pytorch/NN_optuna.py": [
        (r"^DEFAULT_MAX_SIGNAL = .*$", "DEFAULT_MAX_SIGNAL = %d" % nsig, M),
        (r"^DEFAULT_MAX_BACKGROUND = .*$", "DEFAULT_MAX_BACKGROUND = %d" % nbkg, M),
    ],
    "ML_xgboost/XGB_optuna.py": [
        (r"^DEFAULT_MAX_SIGNAL = .*$", "DEFAULT_MAX_SIGNAL = %d" % nsig, M),
        (r"^DEFAULT_MAX_BACKGROUND = .*$", "DEFAULT_MAX_BACKGROUND = %d" % nbkg, M),
    ],
    "ML_pytorch/submit_optuna.sub": [
        (r"^max_signal = .*$", "max_signal = %d" % nsig, M),
        (r"^max_background = .*$", "max_background = %d" % nbkg, M),
    ],
    "ML_xgboost/submit_optuna.sub": [
        (r"^max_signal = .*$", "max_signal = %d" % nsig, M),
        (r"^max_background = .*$", "max_background = %d" % nbkg, M),
    ],
    "ML_tmva/TMVA_config.h": [
        (r"^static const std::vector<TString> FEATURES = \{\n.*?^\};",
         "static const std::vector<TString> FEATURES = {\n" + ",\n".join('   "%s"' % f for f in feats) + "\n};", S),
        (r"^static const TString BACKGROUND_SIDEBAND_CUT =\n.*?;$",
         'static const TString BACKGROUND_SIDEBAND_CUT =\n   "%s";' % side, S),
    ] + [(r"^static const double %s = .*;$" % n, "static const double %s = %s;" % (n, fom[n]), M) for n in fom_names],
    "ML_tmva/TMVA_optimize.C": [
        (r"Long64_t maxSignal = \d+", "Long64_t maxSignal = %d" % nsig, 0),
        (r"Long64_t maxBackground = \d+", "Long64_t maxBackground = %d" % nbkg, 0),
    ],
    "ML_tmva/TMVA_train.C": [
        (r"Long64_t maxSignal = \d+", "Long64_t maxSignal = %d" % nsig, 0),
        (r"Long64_t maxBackground = \d+", "Long64_t maxBackground = %d" % nbkg, 0),
    ],
}

new_text = {}
for path, rules in edits.items():
    text = open(path).read()
    for pattern, repl, flags in rules:
        text, n = re.subn(pattern, lambda _m, r=repl: r, text, flags=flags)
        if n != 1:
            sys.exit("%s: expected 1 match for %r, found %d. Nothing was written." % (path, pattern, n))
    new_text[path] = text

changed = [p for p, t in new_text.items() if open(p).read() != t]
for path in changed:
    dest = os.path.join(backup, path)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    shutil.copy2(path, dest)
    open(path, "w").write(new_text[path])
print("Changed files:", ", ".join(changed) if changed else "none, all already in sync")
if changed: print("Previous versions saved in", backup)
PY

# ---- 2. check the settings against the input files -------------------------
echo "=== 2. Checking branches, cuts and skims"
if grep -q "release 8" /etc/redhat-release; then
  VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt
else
  VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt
fi
(
  source "$VIEW/setup.sh" >/dev/null 2>&1
  cd ML_common
  SAMPLES_LIST="${SAMPLES[*]}" LOG_FEATURES_LIST="${LOG_FEATURES[*]}" ML_USE_XROOTD=1 PYTHONDONTWRITEBYTECODE=1 python - <<'PY'
import os, sys
import ROOT
import uproot
import ml_samples as cfg

ROOT.gErrorIgnoreLevel = ROOT.kFatal
tmva_text = open("../ML_tmva/TMVA_config.h").read()
needed = cfg.feature_universe() + ["Bmass"]
problems = []
for name in os.environ["LOG_FEATURES_LIST"].split():
    if name not in needed:
        problems.append("LOG_FEATURES: %s is not a model input or pool feature" % name)
for sample in os.environ["SAMPLES_LIST"].split():
    pre_cut = cfg.PRE_CUTS[sample]
    for key, tree_name in cfg.TREES.items():
        path = cfg.FLAT_INPUTS[sample][key]
        if path not in tmva_text:
            problems.append("%s %s: input path differs between ml_samples.py and TMVA_config.h" % (sample, key))
        f = ROOT.TFile.Open(cfg.resolve_input_path(path))
        tree = f.Get(tree_name) if f else None
        if not tree:
            problems.append("%s %s: cannot open tree %s in %s" % (sample, key, tree_name, path)); continue
        branches = {b.GetName() for b in tree.GetListOfBranches()}
        missing = [x for x in needed + ([] if key == "data" else [cfg.MC_WEIGHT_BRANCH]) if x not in branches]
        if missing: problems.append("%s %s: missing branches %s" % (sample, key, missing))
        cuts = [("pre-cut", pre_cut)] + ([("training sidebands", cfg.sideband_cut())] if key == "data" else [])
        counts = {}
        for name, cut in cuts:
            used = [n for n in cfg.branches_from_cut(cut) if n not in branches]
            if used: problems.append("%s %s: %s uses missing branches %s" % (sample, key, name, used)); continue
            if ROOT.TTreeFormula("check", cut, tree).GetNdim() == 0:
                problems.append("%s %s: %s does not parse in ROOT" % (sample, key, name)); continue
            # NN and XGB evaluate the cut with NumPy after turning && and || into & and |;
            # both readings must select the same candidates.
            n_root = tree.GetEntries(cut)
            _, arrays = cfg.read_branches(path, tree_name, list(dict.fromkeys(["Bmass"] + cfg.branches_from_cut(cut))))
            n_numpy = int(cfg.cut_mask(arrays, cut).sum())
            if n_root != n_numpy:
                problems.append("%s %s: %s selects %d entries in ROOT but %d in NumPy" % (sample, key, name, n_root, n_numpy))
            counts[name] = n_root
        entries = tree.GetEntries()
        f.Close()
        # The skim must hold exactly this file after the current pre-cut, with every _dc feature.
        skim = cfg.resolve_input_path(cfg.skim_path(sample, key))
        stale = "skim %s: %%s. Remake it: python3 make_skims.py --sample %s" % (cfg.skim_path(sample, key), sample)
        try:
            with uproot.open(skim) as s:
                meta = {k.split(";")[0].split("/")[-1]: str(s[k]) for k in s.keys() if k.startswith("metadata/")}
                skim_branches = set(s[tree_name].keys())
        except FileNotFoundError:
            problems.append(stale % "missing"); continue
        wanted = needed + [f + cfg.DECORRELATED_SUFFIX for f in cfg.feature_universe()]
        if meta["pre_cut"] != pre_cut:
            problems.append(stale % ("made with pre-cut %r" % meta["pre_cut"]))
        elif meta["source_file"] != path or int(meta["source_entries"]) != entries or int(meta["entries"]) != counts.get("pre-cut"):
            problems.append(stale % "made from another version of the flat file")
        elif [b for b in wanted if b not in skim_branches]:
            problems.append(stale % ("missing branches %s" % [b for b in wanted if b not in skim_branches]))
        print("  ok  %-7s %-20s %d entries; %s" % (sample, key, entries, ", ".join("%s %d" % kv for kv in counts.items())))
if problems:
    print("\nProblems found:"); [print("  " + p) for p in problems]; sys.exit(1)
print("All branches exist, all cuts parse, ROOT and NumPy select the same candidates, and the skims are current.")
PY
)

if [ "$MODE" = "--sync-only" ]; then
  echo "Sync and checks done. Nothing deleted or submitted."; exit 0
fi

# ---- 3. remove old jobs and delete old outputs -----------------------------
echo "=== 3. Cleaning"
module load lxbatch/eossubmit
myschedd out || true

# Only the selected methods and samples are cleaned: their queued jobs, their
# per-sample output folder, their Condor logs, their scored files and the Punzi
# results made from those scored files (optimalCUT/<method>/<sample>).
shopt -s nullglob
FOLDERS=(); TARGETS=()
for m in "${METHODS[@]}"; do
  case $m in nn) dir=ML_pytorch; out=nn_outputs;; xgb) dir=ML_xgboost; out=xgb_outputs;; tmva) dir=ML_tmva; out=tmva_outputs;; esac
  FOLDERS+=("${dir#ML_}")
  for s in "${SAMPLES[@]}"; do
    low=$(echo "$s" | tr 'A-Z' 'a-z')
    TARGETS+=("$dir/$out/ntmix_PbPb/${low}_${FOM}" $dir/condor_logs/*_${s}_${FOM}_* $dir/scored_samples/flat_ntmix_${s}_${FOM}_scored_* "optimalCUT/$m/${s}_${FOM}")
  done
done
has xgb  && TARGETS+=(ML_xgboost/__pycache__)
has tmva && TARGETS+=(ML_tmva/dataset)
EXISTING=(); for d in "${TARGETS[@]}"; do [ -e "$d" ] && EXISTING+=("$d"); done

ML_JOBS="regexp(\"/selectionER/ML_($(IFS='|'; echo "${FOLDERS[*]}"))/condor_logs/[a-z]+_($(IFS='|'; echo "${SAMPLES[*]}"))_${FOM}_\", UserLog)"
OLD_CLUSTERS=$(condor_q -constraint "$ML_JOBS" -af ClusterId 2>/dev/null | sort -u | tr '\n' ' ')
echo "Jobs to remove:"
if [ -n "${OLD_CLUSTERS// }" ]; then
  for c in $OLD_CLUSTERS; do
    echo "  cluster $c: $(condor_q $c -af UserLog 2>/dev/null | head -1 | sed 's#.*/selectionER/##'), $(condor_q $c -af ClusterId 2>/dev/null | wc -l) jobs"
  done
else
  echo "  none"
fi
echo "To delete:"
for d in "${EXISTING[@]}"; do [ -d "$d" ] && du -sh "$d"; done
nfiles=0; for d in "${EXISTING[@]}"; do [ -f "$d" ] && nfiles=$((nfiles + 1)); done
echo "  and $nfiles log or scored files"
[ ${#EXISTING[@]} -eq 0 ] && echo "  nothing"
echo "Methods to submit: ${METHODS[*]}"
echo "Samples to submit: ${SAMPLES[*]}"
echo "Model selection FOM: $FOM"
read -r -p "Remove these jobs, delete these outputs and submit? [y/N] " answer
[ "$answer" = "y" ] || { echo "Aborted. Settings were synced, nothing deleted or submitted."; exit 1; }

if [ -n "${OLD_CLUSTERS// }" ]; then
  condor_rm $OLD_CLUSTERS 2>/dev/null || true
  # Condor keeps writing to the job logs until the jobs have left the queue.
  for i in $(seq 1 60); do
    left=$(condor_q $OLD_CLUSTERS -af ClusterId 2>/dev/null | wc -l)
    [ "$left" -eq 0 ] && break
    echo "  waiting for $left old jobs to leave the queue..."; sleep 10
  done
fi
if [ ${#EXISTING[@]} -gt 0 ]; then
  for attempt in 1 2 3; do
    rm -rf "${EXISTING[@]}" 2>/dev/null && break
    echo "  EOS still busy, retrying in 10 s..."; sleep 10
  done
fi
for d in "${EXISTING[@]}"; do
  [ -e "$d" ] && { echo "Could not delete $d. Delete it by hand and rerun."; exit 1; }
done

for s in "${SAMPLES[@]}"; do
  low=$(echo "$s" | tr 'A-Z' 'a-z')
  has nn   && mkdir -p ML_pytorch/nn_outputs/ntmix_PbPb/${low}_${FOM}/optuna_summaries ML_pytorch/condor_logs
  has xgb  && mkdir -p ML_xgboost/xgb_outputs/ntmix_PbPb/${low}_${FOM}/optuna_summaries ML_xgboost/condor_logs
  has tmva && mkdir -p ML_tmva/tmva_outputs/ntmix_PbPb/${low}_${FOM}/optimization ML_tmva/condor_logs
done

# ---- 4. submit --------------------------------------------------------------
echo "=== 4. Submitting"
for s in "${SAMPLES[@]}"; do
  low=$(echo "$s" | tr 'A-Z' 'a-z')
  has nn   && (cd ML_pytorch && condor_submit -append "sample = $s" -append "fom = $FOM" -append "sample_output_dir = nn_outputs/ntmix_PbPb/${low}_${FOM}" submit_optuna.sub)
  has xgb  && (cd ML_xgboost && condor_submit -append "sample = $s" -append "fom = $FOM" -append "sample_output_dir = xgb_outputs/ntmix_PbPb/${low}_${FOM}" submit_optuna.sub)
  has tmva && (cd ML_tmva    && condor_submit -append "sample = $s" -append "fom = $FOM" -append "sample_output_dir = tmva_outputs/ntmix_PbPb/${low}_${FOM}" submit_optimize.sub)
done
condor_q
