#!/bin/bash
# =============================================================================
# Train the best scan result of every ML method for one sample.
#
# For each method (NN, XGBoost, TMVA) the script
#   1. checks the scan: finished jobs, number of results, and that the results
#      were made with the current features, pre-cut and sidebands from fresh_ML_scan.sh,
#   2. runs the final training locally with the best passing trial (NN, XGB: its scan
#      summary via --params-json; TMVA: ML_USE_OPTIMIZED) and saves the log next to the model,
# and prints a summary of the metrics at the end.
#
# Usage, on lxplus (run inside tmux: training takes a while):
#   ./train_best_ML.sh PbPb23 auc              all three methods, scans made with FOM auc
#   ./train_best_ML.sh PbPb23 auc nn xgb       only the listed methods (nn, xgb, tmva)
# The FOM (auc, sigeff, pauc, punzi) selects the scan run, e.g. *_outputs/ntmix_PbPb/pbpb23_auc.
# =============================================================================
set -o pipefail
BASE=/eos/user/h/hmarques/Analysis_CODES/selectionER
cd "$BASE"

SAMPLE=$1; FOM=$2; shift 2 || true
USAGE="Usage: $0 <PbPb18|PbPb23> <auc|sigeff|pauc|punzi> [nn] [xgb] [tmva]"
case "$FOM" in auc|sigeff|pauc|punzi) ;; *) echo "$USAGE"; exit 1 ;; esac
case "$SAMPLE" in PbPb18|PbPb23) OUT_SUB="ntmix_PbPb/$(echo "$SAMPLE" | tr 'A-Z' 'a-z')_${FOM}" ;; *) echo "$USAGE"; exit 1 ;; esac
METHODS=("$@"); [ ${#METHODS[@]} -eq 0 ] && METHODS=(nn xgb tmva)
for m in "${METHODS[@]}"; do
  case "$m" in nn|xgb|tmva) ;; *) echo "Unknown method: $m (use nn, xgb, tmva)"; exit 1 ;; esac
done

if grep -q "release 8" /etc/redhat-release; then
  VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt
else
  VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt
fi
export ML_USE_XROOTD=1 PYTHONDONTWRITEBYTECODE=1

method_dir()  { case $1 in nn) echo ML_pytorch;; xgb) echo ML_xgboost;; tmva) echo ML_tmva;; esac; }
output_dir()  { case $1 in nn) echo "ML_pytorch/nn_outputs/$OUT_SUB";; xgb) echo "ML_xgboost/xgb_outputs/$OUT_SUB";; tmva) echo "ML_tmva/tmva_outputs/$OUT_SUB";; esac; }
job_prefix()  { case $1 in tmva) echo optimize;; *) echo optuna;; esac; }

# ---- 1. check the scans ------------------------------------------------------
echo "=== 1. Checking scans for $SAMPLE, FOM $FOM: ${METHODS[*]}"
module load lxbatch/eossubmit >/dev/null 2>&1 || true
unfinished=0
for m in "${METHODS[@]}"; do
  log_pattern="/$(method_dir $m)/condor_logs/$(job_prefix $m)_${SAMPLE}_${FOM}_"
  counts=$(condor_q -constraint "regexp(\"$log_pattern\", UserLog)" -af JobStatus 2>/dev/null |
           awk '{n[$1]++} END {printf "%d idle, %d running, %d held", n[1], n[2], n[5]}')
  if [ -n "$counts" ] && [ "$counts" != "0 idle, 0 running, 0 held" ]; then
    echo "  $m: scan jobs still in the queue: $counts"; unfinished=1
  fi
done

# The check prints "BEST <method> <summary>" for the NN and XGB summaries to train from.
CHECK=$(
  source "$VIEW/setup.sh" >/dev/null 2>&1
  SAMPLE=$SAMPLE FOM=$FOM OUT_SUB=$OUT_SUB METHODS="${METHODS[*]}" python - <<'PY' | tee /dev/stderr
import glob, json, os, re, sys
sys.path.insert(0, "ML_common")
import ml_fom, ml_samples
sample, fom, methods = os.environ["SAMPLE"], os.environ["FOM"], os.environ["METHODS"].split()
norm = lambda cut: re.sub(r"\s+", "", cut.replace("&&", "&").replace("||", "|"))
problems = []

def compare(method, where, scan, current):
    for key, now in current.items():
        then = scan.get(key)
        if then is None: continue
        same = norm(then) == norm(now) if key.endswith("cut") else then == now
        if not same:
            problems.append("%s: scan result %s has %s = %r, current settings have %r"
                            % (method, where, key, then, now))

for method, folder, module in (("nn", "ML_pytorch", "NN_train"), ("xgb", "ML_xgboost", "XGB_train")):
    if method not in methods: continue
    sys.path.insert(0, folder)
    cfg = __import__(module)
    files = glob.glob("%s/%s_outputs/%s/optuna_summaries/optuna_summary_%s_%s_*.json"
                      % (folder, method, os.environ["OUT_SUB"], sample, fom))
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
    c = best["config"]
    compare(method, os.path.basename(path), {
        "features": c["features"], "pre_cut": c["pre_cut"], "training_sidebands": c["training_sidebands"],
        "max_signal": c["max_signal"], "max_background": c["max_background"],
        "fom_settings": c["fom_settings"]}, {
        "features": list(ml_samples.FEATURES), "pre_cut": ml_samples.PRE_CUTS[sample],
        "training_sidebands": [list(s) for s in ml_samples.TRAINING_SIDEBANDS],
        "max_signal": cfg.MAX_SIGNAL, "max_background": cfg.MAX_BACKGROUND,
        "fom_settings": ml_fom.settings()})
    print("BEST %s %s" % (method, os.path.abspath(path)))

if "tmva" in methods:
    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kError
    # Read the /eos mount directly; the XRootD redirect fails intermittently here.
    ROOT.gEnv.SetValue("TFile.CrossProtocolRedirects", 0)
    ROOT.gInterpreter.ProcessLine('#include "ML_tmva/TMVA_config.h"')
    T = ROOT.MLTMVA
    train_c = open("ML_tmva/TMVA_train.C").read()
    current = {"features": ",".join(str(f) for f in T.FEATURES), "signal_cut": ml_samples.PRE_CUTS[sample],
               "background_sideband_cut": str(T.BACKGROUND_SIDEBAND_CUT),
               "max_signal": re.search(r"Long64_t maxSignal = (\d+)", train_c).group(1),
               "max_background": re.search(r"Long64_t maxBackground = (\d+)", train_c).group(1)}
    files = glob.glob("ML_tmva/tmva_outputs/%s/optimization/summary_job_*.root" % os.environ["OUT_SUB"])
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
    if not files:
        problems.append("tmva: no scan results found")
    elif best is None:
        problems.append("tmva: none of the %d scan points passes the overtraining test" % n_points)
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
best_summary() { echo "$CHECK" | sed -n "s/^BEST $1 //p"; }

if [ $unfinished -eq 1 ]; then
  read -r -p "Some scans are not finished. Train on the results available now? [y/N] " answer
  [ "$answer" = "y" ] || { echo "Stopped. Nothing was changed."; exit 1; }
fi

# ---- 2. train ----------------------------------------------------------------
declare -A STATUS
for m in "${METHODS[@]}"; do
  out="$(output_dir $m)"; mkdir -p "$out"; log="$out/training_log.txt"
  echo "=== $m: training (log: $log)"
  (
    source "$VIEW/setup.sh" >/dev/null 2>&1
    cd "$(method_dir $m)"
    case $m in
      nn)   echo "Used summary: $(best_summary nn)"
            python NN_train.py --sample "$SAMPLE" --fom "$FOM" --params-json "$(best_summary nn)" ;;
      xgb)  echo "Used summary: $(best_summary xgb)"
            python XGB_train.py --sample "$SAMPLE" --fom "$FOM" --params-json "$(best_summary xgb)" ;;
      tmva) ML_USE_OPTIMIZED=1 root -l -b -q "TMVA_train.C(\"$SAMPLE\",\"$FOM\")" ;;
    esac
  ) 2>&1 | tee "$BASE/$log"
  if [ $? -eq 0 ]; then STATUS[$m]=ok; else STATUS[$m]=FAILED; fi
done

# ---- summary -----------------------------------------------------------------
echo
echo "=== Summary for $SAMPLE, FOM $FOM"
for m in "${METHODS[@]}"; do
  log="$(output_dir $m)/training_log.txt"
  echo "--- $m: ${STATUS[$m]}"
  grep -hE "^(Used summary|Best optimization objective|Selected FOM|Validation:|Test:)" "$log" | sed 's/^/    /'
done
echo
echo "Next: score the samples, then optimize the threshold on the scored samples:"
for m in "${METHODS[@]}"; do
  case $m in
    nn)   echo "  (cd ML_pytorch && ./run_nn_batch.sh NN_apply.py --sample $SAMPLE --fom $FOM --kind all)" ;;
    xgb)  echo "  (cd ML_xgboost && source $VIEW/setup.sh && python XGB_apply.py --sample $SAMPLE --fom $FOM --kind all)" ;;
    tmva) echo "  (cd ML_tmva && root -l -b -q 'TMVA_apply.C(\"$SAMPLE\",\"$FOM\",\"all\")')" ;;
  esac
  echo "  root -l -b -q 'optimalCUT_X_punzi.C(\"${SAMPLE/ppRef24/ppRef}\",\"$m\",\"$FOM\")'"
done
