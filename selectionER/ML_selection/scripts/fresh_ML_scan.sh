#!/bin/bash
# =============================================================================
# Central control for the ML hyperparameter scans (ML_pytorch, ML_xgboost, ML_tmva) of one or more
# cases "<tree>/<year>" (ntmix/PbPb23, ntKp/ppRef, ...).
#
# All settings (cases, features, pre-cuts, sidebands, mass windows, sample sizes, FOM settings) are in
# ML_selection/settings.py. The script
#   1. writes their C++ copy for TMVA and optimalCUT_X_punzi.C (methods/ML_tmva/ml_settings.h),
#   2. checks the selected cases against their input files: every branch the case uses exists, the
#      pre-cut and the training sidebands parse, ROOT and the NumPy translation used by NN/XGB select
#      the same candidates, and the skims (make_skims.py), which every training reads, match the
#      current pre-cut and flat files,
#   3. asks for confirmation, removes old scan jobs and deletes the old outputs of the selected
#      methods and cases,
#   4. submits the scans: NN 200, XGB 250, TMVA 384 jobs per case.
#
# Every scan run is tagged with its model-selection FOM (auc, sigeff, pauc or punzi): outputs, logs and
# scored samples of ntKp/PbPb23 with FOM auc live in methods/ML_<method>/<method>_outputs/ntKp/PbPb23_auc,
# methods/ML_<method>/condor_logs/*_ntKp_PbPb23_auc_* and scored_samples/ntKp/PbPb23/*_<method>_auc_*, so
# runs with different FOMs are kept side by side and cleaned separately.
# The feature study (feature_study.py) uses the same settings; its adopt command writes the chosen
# set into the FEATURES of the tree in settings.py.
#
# Usage, from any directory on lxplus:
#   ./fresh_ML_scan.sh                      every method and case, FOM of settings.py
#   ./fresh_ML_scan.sh sigeff               the same with another FOM (auc, sigeff, pauc, punzi)
#   ./fresh_ML_scan.sh tmva                 clean and resubmit only the listed methods (nn, xgb, tmva)
#   ./fresh_ML_scan.sh ntKp/PbPb23          clean and resubmit only the listed cases
#   ./fresh_ML_scan.sh ntKp/PbPb23 tmva     these can be combined; everything else is left alone
#   ./fresh_ML_scan.sh --sync-only          steps 1 and 2 only; nothing deleted or submitted
# =============================================================================

set -eo pipefail
ML=/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_selection
SELECTION=$(dirname "$ML")
cd "$ML/methods"
if grep -q "release 8" /etc/redhat-release; then
  VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt
else
  VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt
fi
# The FOM and the cases of settings.py.
read -r FOM ALL_CASES <<< "$(source "$VIEW/setup.sh" >/dev/null 2>&1; cd ML_common &&
                             PYTHONDONTWRITEBYTECODE=1 python -c 'import ml_samples; print(ml_samples.settings.FOM, " ".join(ml_samples.CASES))')"
MODE=full
METHODS=()
CASES=()
for arg in "$@"; do
  case "$arg" in
    --sync-only) MODE=--sync-only ;;
    nn|xgb|tmva) METHODS+=("$arg") ;;
    auc|sigeff|pauc|punzi) FOM="$arg" ;;
    */*) CASES+=("$arg") ;;
    *) echo "Unknown option: $arg (use nn, xgb, tmva, auc, sigeff, pauc, punzi, a case <tree>/<year> of: $ALL_CASES, or --sync-only)"; exit 1 ;;
  esac
done
[ ${#METHODS[@]} -eq 0 ] && METHODS=(nn xgb tmva)
[ ${#CASES[@]} -eq 0 ] && read -r -a CASES <<< "$ALL_CASES"

# ---- 1. write the C++ copy of the settings ----------------------------------
echo "=== 1. Writing methods/ML_tmva/ml_settings.h from settings.py"
(source "$VIEW/setup.sh" >/dev/null 2>&1; cd ML_common && PYTHONDONTWRITEBYTECODE=1 python -c 'import ml_samples; ml_samples.write_cpp_settings()')

# ---- 2. check the settings against the input files -------------------------
echo "=== 2. Checking branches, cuts and skims of ${CASES[*]}"
(
  source "$VIEW/setup.sh" >/dev/null 2>&1
  cd ML_common
  CASES_LIST="${CASES[*]}" ML_USE_XROOTD=1 PYTHONDONTWRITEBYTECODE=1 python - <<'PY'
import os, sys
import ROOT
import uproot
import ml_samples

ROOT.gErrorIgnoreLevel = ROOT.kFatal
problems = []
for name in os.environ["CASES_LIST"].split():
    c = ml_samples.case(name)
    needed = c["universe"] + ["Bmass"]
    flat = list(dict.fromkeys(f.removesuffix(ml_samples.LOG_SUFFIX) for f in needed))   # <branch>_log: log10 of the branch
    for key, part in c["inputs"].items():
        f = ROOT.TFile.Open(ml_samples.XROOTD + part["file"])
        tree = f.Get(part["tree"])
        branches = {b.GetName() for b in tree.GetListOfBranches()}
        missing = [x for x in flat + ([] if key == "data" else [ml_samples.MC_WEIGHT_BRANCH]) if x not in branches]
        if missing:
            problems.append("%s %s: missing branches %s" % (name, key, missing))
        cuts = [("pre-cut", c["pre_cut"])] + ([("training sidebands", c["training_sideband_cut"])] if key == "data" else [])
        counts = {}
        for label, cut in cuts:
            if ROOT.TTreeFormula("check", cut, tree).GetNdim() == 0:
                problems.append("%s %s: %s does not parse in ROOT" % (name, key, label))
                continue
            # NN and XGB evaluate the cut with NumPy after turning && and || into & and |;
            # both readings must select the same candidates.
            counts[label] = tree.GetEntries(cut)
            names = list(dict.fromkeys(["Bmass"] + [n for n in ml_samples.CUT_NAMES.findall(cut) if n != "abs"]))
            n_numpy = sum(int(ml_samples.cut_mask(arrays, cut).sum()) for arrays in
                          uproot.iterate(f"{ml_samples.XROOTD}{part['file']}:{part['tree']}", names, step_size="500 MB", library="np"))
            if counts[label] != n_numpy:
                problems.append("%s %s: %s selects %d entries in ROOT but %d in NumPy" % (name, key, label, counts[label], n_numpy))
        entries, skim_entries = tree.GetEntries(), tree.GetEntries(part["skim_cut"])
        f.Close()
        # The skim must hold exactly this file after the current pre-cut, with every _dc feature.
        stale = "skim %s: %%s. Remake it: python3 make_skims.py --tree %s --year %s" % (part["skim"], c["tree"], c["year"])
        with uproot.open(ml_samples.XROOTD + part["skim"]) as s:
            meta = {k.split(";")[0].split("/")[-1]: str(s[k]) for k in s.keys() if k.startswith("metadata/")}
            skim_branches = set(s[part["tree"]].keys())
        wanted = needed + [f + ml_samples.DECORRELATED_SUFFIX for f in c["universe"]]
        if meta["pre_cut"] != c["pre_cut"]:
            problems.append(stale % ("made with pre-cut %r" % meta["pre_cut"]))
        elif meta["source_file"] != part["file"] or int(meta["source_entries"]) != entries or int(meta["entries"]) != skim_entries:
            problems.append(stale % "made from another version of the flat file")
        elif int(meta.get("data_thinning", 1)) != (c["data_thinning"] if key == "data" else 1):
            problems.append(stale % ("made with DATA thinning %s" % meta.get("data_thinning", 1)))
        elif [b for b in wanted if b not in skim_branches]:
            problems.append(stale % ("missing branches %s" % [b for b in wanted if b not in skim_branches]))
        print("  ok  %-15s %-20s %d entries; %s" % (name, key, entries, ", ".join("%s %d" % kv for kv in counts.items())))
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

# Only the selected methods and cases are cleaned: their queued jobs, their run output folder,
# their Condor logs, their scored files and the Punzi results made from those scored files.
shopt -s nullglob
FOLDERS=(); TAGS=(); TARGETS=()
for m in "${METHODS[@]}"; do
  case $m in nn) dir=ML_pytorch; out=nn_outputs;; xgb) dir=ML_xgboost; out=xgb_outputs;; tmva) dir=ML_tmva; out=tmva_outputs;; esac
  FOLDERS+=("${dir#ML_}")
  for c in "${CASES[@]}"; do
    tree=${c%/*}; year=${c#*/}
    TARGETS+=("$dir/$out/$tree/${year}_${FOM}" $dir/condor_logs/*_${tree}_${year}_${FOM}_*
              $SELECTION/scored_samples/$tree/$year/flat_${tree}_${year}_${m}_${FOM}_scored_* "$SELECTION/optimalCUT/$tree/$year/${m}_${FOM}")
  done
done
for c in "${CASES[@]}"; do TAGS+=("${c/\//_}"); done
EXISTING=(); for d in "${TARGETS[@]}"; do [ -e "$d" ] && EXISTING+=("$d"); done

ML_JOBS="regexp(\"/ML_selection/methods/ML_($(IFS='|'; echo "${FOLDERS[*]}"))/condor_logs/[a-z]+_($(IFS='|'; echo "${TAGS[*]}"))_${FOM}_\", UserLog)"
OLD_CLUSTERS=$(condor_q -constraint "$ML_JOBS" -af ClusterId 2>/dev/null | sort -u | tr '\n' ' ')
echo "Jobs to remove:"
if [ -n "${OLD_CLUSTERS// }" ]; then
  for cluster in $OLD_CLUSTERS; do
    echo "  cluster $cluster: $(condor_q $cluster -af UserLog 2>/dev/null | head -1 | sed 's#.*/ML_selection/##'), $(condor_q $cluster -af ClusterId 2>/dev/null | wc -l) jobs"
  done
else
  echo "  none"
fi
echo "To delete:"
for d in "${EXISTING[@]}"; do [ -d "$d" ] && du -sh "$d"; done
nfiles=0; for d in "${EXISTING[@]}"; do [ -f "$d" ] && nfiles=$((nfiles + 1)); done
echo "  and $nfiles log or scored files"
echo "Methods to submit: ${METHODS[*]}"
echo "Cases to submit: ${CASES[*]}"
echo "Model selection FOM: $FOM"
read -r -p "Remove these jobs, delete these outputs and submit? [y/N] " answer
[ "$answer" = "y" ] || { echo "Aborted. The C++ settings were written, nothing deleted or submitted."; exit 1; }

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
  # EOS can be busy for a while after the jobs left; the last attempt stops the script if it fails.
  for attempt in 1 2; do
    rm -rf "${EXISTING[@]}" 2>/dev/null && break
    echo "  EOS still busy, retrying in 10 s..."; sleep 10
  done
  rm -rf "${EXISTING[@]}"
fi

# ---- 4. submit --------------------------------------------------------------
echo "=== 4. Submitting"
for c in "${CASES[@]}"; do
  tree=${c%/*}; year=${c#*/}
  for m in "${METHODS[@]}"; do
    case $m in
      nn)   mkdir -p ML_pytorch/nn_outputs/$tree/${year}_${FOM}/optuna_summaries ML_pytorch/condor_logs
            (cd ML_pytorch && condor_submit -append "tree = $tree" -append "year = $year" -append "fom = $FOM" submit_optuna.sub) ;;
      xgb)  mkdir -p ML_xgboost/xgb_outputs/$tree/${year}_${FOM}/optuna_summaries ML_xgboost/condor_logs
            (cd ML_xgboost && condor_submit -append "tree = $tree" -append "year = $year" -append "fom = $FOM" submit_optuna.sub) ;;
      tmva) mkdir -p ML_tmva/tmva_outputs/$tree/${year}_${FOM}/optimization ML_tmva/condor_logs
            (cd ML_tmva && condor_submit -append "tree = $tree" -append "year = $year" -append "fom = $FOM" submit_optimize.sub) ;;
    esac
  done
done
condor_q
