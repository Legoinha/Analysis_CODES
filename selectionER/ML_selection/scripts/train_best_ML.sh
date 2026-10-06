#!/bin/bash
# =============================================================================
# Train the best scan result of every ML method for one case "<tree>/<year>".
#
# For each method (NN, XGBoost, TMVA) the script
#   1. checks the scan: finished jobs, number of results, and that the results were made with the
#      current features, pre-cut, sidebands, sample sizes and FOM settings of ML_selection/settings.py,
#   2. runs the final training locally with the best passing trial (NN, XGB: its scan
#      summary via --params-json; TMVA: ML_USE_OPTIMIZED) and saves the log next to the model,
# and prints a summary of the metrics at the end, and the commands that score the samples and
# optimize the threshold.
#
# Usage, on lxplus (run inside tmux: training takes a while):
#   ./train_best_ML.sh ntKp/PbPb23 auc              all three methods, scans made with FOM auc
#   ./train_best_ML.sh ntKp/PbPb23 auc nn xgb       only the listed methods (nn, xgb, tmva)
# The FOM (auc, sigeff, pauc, punzi) selects the scan run, e.g. <method>_outputs/ntKp/PbPb23_auc.
# =============================================================================
set -o pipefail
ML=/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_selection
cd "$ML/methods"

CASE=$1; FOM=$2; shift 2
TREE=${CASE%/*}; YEAR=${CASE#*/}
METHODS=("$@"); [ ${#METHODS[@]} -eq 0 ] && METHODS=(nn xgb tmva)

if grep -q "release 8" /etc/redhat-release; then
  VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt
else
  VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt
fi
export ML_USE_XROOTD=1 PYTHONDONTWRITEBYTECODE=1

# ---- 1. check the scans ------------------------------------------------------
echo "=== 1. Checking scans for $CASE, FOM $FOM: ${METHODS[*]}"
module load lxbatch/eossubmit >/dev/null 2>&1 || true
unfinished=0
for m in "${METHODS[@]}"; do
  case $m in nn) log_pattern="/ML_pytorch/condor_logs/optuna_";; xgb) log_pattern="/ML_xgboost/condor_logs/optuna_";; tmva) log_pattern="/ML_tmva/condor_logs/optimize_";; esac
  counts=$(condor_q -constraint "regexp(\"${log_pattern}${TREE}_${YEAR}_${FOM}_\", UserLog)" -af JobStatus 2>/dev/null |
           awk '{n[$1]++} END {printf "%d idle, %d running, %d held", n[1], n[2], n[5]}')
  if [ -n "$counts" ] && [ "$counts" != "0 idle, 0 running, 0 held" ]; then
    echo "  $m: scan jobs still in the queue: $counts"; unfinished=1
  fi
done

# The check writes ml_settings.h for TMVA and prints "BEST <method> <summary>" for the NN and XGB
# summaries to train from.
CHECK=$(
  source "$VIEW/setup.sh" >/dev/null 2>&1
  CASE=$CASE FOM=$FOM METHODS="${METHODS[*]}" python - <<'PY' | tee /dev/stderr
import glob, json, os, re, sys
sys.path.insert(0, "ML_common")
import ml_fom, ml_samples
from ml_samples import settings
ml_samples.write_cpp_settings()
c, fom, methods = ml_samples.case(os.environ["CASE"]), os.environ["FOM"], os.environ["METHODS"].split()
run = f"{c['tree']}/{c['year']}_{fom}"
problems = []

def compare(method, where, scan, current):
    for key, now in current.items():
        then = scan.get(key)
        if then is None: continue
        if key.endswith("cut"):
            then, now_cut = (re.sub(r"\s+", "", cut.replace("&&", "&").replace("||", "|")) for cut in (then, now))
            same = then == now_cut
        else:
            same = then == now
        if not same:
            problems.append("%s: scan result %s has %s = %r, current settings have %r" % (method, where, key, scan.get(key), now))

for method, folder in (("nn", "ML_pytorch"), ("xgb", "ML_xgboost")):
    if method not in methods: continue
    files = glob.glob(f"{folder}/{method}_outputs/{run}/optuna_summaries/optuna_summary_{c['tree']}_{c['year']}_{fom}_*.json")
    summaries = [(json.load(open(f)), f) for f in files]
    passing = [s for s in summaries if s[0]["best_value"] is not None]
    if not passing:
        problems.append("%s: %d scan results, none with a trial that passes the overtraining test" % (method, len(files)))
        continue
    best, path = max(passing, key=lambda s: s[0]["best_value"])
    u = best["best_trial"]["user_attrs"]
    print("  %-4s %d results (%d with a passing trial), best validation %s %.5f (validation AUC %.4f, KS p %.3f/%.3f) in %s"
          % (method, len(files), len(passing), fom, best["best_value"], u["auc"], u["ks_signal_pvalue"],
             u["ks_background_pvalue"], os.path.basename(path)))
    cfg = best["config"]
    compare(method, os.path.basename(path), {
        "features": cfg["features"], "pre_cut": cfg["pre_cut"], "training_sidebands": cfg["training_sidebands"],
        "max_signal": cfg["max_signal"], "max_background": cfg["max_background"], "fom_settings": cfg["fom_settings"]}, {
        "features": c["features"], "pre_cut": c["pre_cut"], "training_sidebands": [list(s) for s in c["training_sidebands"]],
        "max_signal": settings.MAX_SIGNAL, "max_background": settings.MAX_BACKGROUND, "fom_settings": ml_fom.fom_settings(c)})
    print("BEST %s %s" % (method, os.path.abspath(path)))

if "tmva" in methods:
    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kError
    # Read the /eos mount directly; the XRootD redirect fails intermittently here.
    ROOT.gEnv.SetValue("TFile.CrossProtocolRedirects", 0)
    current = {"features": ",".join(c["features"]), "signal_cut": c["pre_cut"],
               "background_sideband_cut": c["training_sideband_cut"],
               "max_signal": str(settings.MAX_SIGNAL), "max_background": str(settings.MAX_BACKGROUND)}
    files = glob.glob(f"ML_tmva/tmva_outputs/{run}/optimization/summary_job_*.root")
    best, n_passed, n_points = None, 0, 0
    for f in files:
        rf = ROOT.TFile.Open(f)
        tree = rf.Get("trial")
        for entry in range(tree.GetEntries()):  # one entry per tree-count checkpoint
            tree.GetEntry(entry)
            n_points += 1
            if not tree.passed: continue  # failed the overtraining test
            n_passed += 1
            if best is None or tree.objective > best[0]:
                best = (tree.objective, tree.validation_auc, "%s at %d trees" % (os.path.basename(f), tree.n_trees))
        compare("tmva", os.path.basename(f), {k: rf.Get(k).GetTitle() for k in current}, current)
        rf.Close()
    if best is None:
        problems.append("tmva: none of the %d scan points of %d results passes the overtraining test" % (n_points, len(files)))
    else:
        print("  tmva %d results, %d of %d points pass the overtraining test, best validation %s %.5f "
              "(validation AUC %.4f) in %s" % ((len(files), n_passed, n_points, fom) + best))

if problems:
    print("\nThe scan results do not match the current settings:")
    for p in problems[:10]: print("  " + p)
    if len(problems) > 10: print("  ... and %d more" % (len(problems) - 10))
    print("Run a fresh scan with ./fresh_ML_scan.sh, or restore the settings used by the scan.")
    sys.exit(1)
PY
) || exit 1

if [ $unfinished -eq 1 ]; then
  read -r -p "Some scans are not finished. Train on the results available now? [y/N] " answer
  [ "$answer" = "y" ] || { echo "Stopped. Nothing was changed."; exit 1; }
fi

# ---- 2. train ----------------------------------------------------------------
declare -A STATUS LOGS
for m in "${METHODS[@]}"; do
  case $m in nn) dir=ML_pytorch;; xgb) dir=ML_xgboost;; tmva) dir=ML_tmva;; esac
  out="$dir/${m}_outputs/$TREE/${YEAR}_${FOM}"; mkdir -p "$out"; LOGS[$m]="$out/training_log.txt"
  echo "=== $m: training (log: ${LOGS[$m]})"
  (
    source "$VIEW/setup.sh" >/dev/null 2>&1
    cd "$dir"
    case $m in
      nn)   summary=$(echo "$CHECK" | sed -n "s/^BEST nn //p"); echo "Used summary: $summary"
            python NN_train.py --case "$CASE" --fom "$FOM" --params-json "$summary" ;;
      xgb)  summary=$(echo "$CHECK" | sed -n "s/^BEST xgb //p"); echo "Used summary: $summary"
            python XGB_train.py --case "$CASE" --fom "$FOM" --params-json "$summary" ;;
      tmva) ML_USE_OPTIMIZED=1 root -l -b -q "TMVA_train.C(\"$CASE\",\"$FOM\")" ;;
    esac
  ) 2>&1 | tee "$ML/methods/${LOGS[$m]}"
  if [ $? -eq 0 ]; then STATUS[$m]=ok; else STATUS[$m]=FAILED; fi
done

# ---- summary -----------------------------------------------------------------
echo
echo "=== Summary for $CASE, FOM $FOM"
for m in "${METHODS[@]}"; do
  echo "--- $m: ${STATUS[$m]}"
  grep -hE "^(Used summary|Best optimization objective|Selected FOM|Validation:|Test:)" "${LOGS[$m]}" | sed 's/^/    /'
done
echo
echo "Next: score the samples (into selectionER/scored_samples/$CASE/), then optimize the threshold on them:"
for m in "${METHODS[@]}"; do
  case $m in
    nn)   echo "  (cd $ML/methods/ML_pytorch && ./run_nn_batch.sh NN_apply.py --case $CASE --fom $FOM --kind all)" ;;
    xgb)  echo "  (cd $ML/methods/ML_xgboost && source $VIEW/setup.sh && python XGB_apply.py --case $CASE --fom $FOM --kind all)" ;;
    tmva) echo "  (cd $ML/methods/ML_tmva && root -l -b -q 'TMVA_apply.C(\"$CASE\",\"$FOM\",\"all\")')" ;;
  esac
  echo "  (cd $(dirname "$ML") && root -l -b -q 'optimalCUT_X_punzi.C(\"$CASE\",\"$m\",\"$FOM\")')"
done
