# TMVA BDTG workflow

> **Central settings** are in `ML_selection/settings.py` for all three methods; this folder reads their
> C++ copy `ml_settings.h`, which the scripts of `ML_selection/scripts/` write again when they start
> (`python3 -c 'import ml_samples; ml_samples.write_cpp_settings()'` from `../ML_common` by hand).
> `../../scripts/fresh_ML_scan.sh` checks the inputs, cleans and submits the scans of all three methods;
> after the scans, `../../scripts/train_best_ML.sh ntKp/PbPb23 auc` trains the best passing trial. The
> feature study (`../../scripts/feature_study.py`) drives the same macros.

The ROOT-native reference classifier parallel to `ML_xgboost` and `ML_pytorch`: a gradient-boosted
decision tree with the signal MC of a case as signal and its DATA mass sidebands as background. The
spectator MC samples of the X(3872) cases (nonprompt X(3872), prompt and nonprompt ψ(2S)) are only
scored. Every macro takes the case `"<tree>/<year>"` and the FOM as its first arguments (or `ML_CASE`
and `ML_FOM`).

## Files

- `TMVA_config.h`: TMVA settings (method, scan grid, random-search ranges) and `LoadRun`, the case at run
  time: skims (with `ML_VARIANT`, through XRootD with `ML_USE_XROOTD=1`), inputs (`ML_FEATURES`, comma-
  separated, `_dc` names allowed, or the production features), pre-cut of the skims, output folder
  (`ML_OUTPUT_DIR`, or `tmva_outputs/<tree>/<year>_<fom>`).
- `ml_settings.h`: written from `settings.py`, do not edit.
- `TMVA_common.h`: forest evaluation tree by tree, weighted ROC/AUC, the FOMs and the overtraining test
  with the definitions of `ML_common/ml_fom.py`.
- `TMVA_decorrelation.h`: the `_dc` transform for scoring, from the maps of the skims.
- `TMVA_train.C`: final training, plots, XML model, ROOT report and `scores.root`.
- `TMVA_optimize.C`: one hyperparameter scan job.
- `TMVA_apply.C`: the XML model applied to the flat files of the case, with `Prediction` and the metadata.
- `run_tmva.sh`: EOS-resident launcher required by CERN EosSubmit.
- `submit_optimize.sub`, `submit_train.sub`, `submit_score.sub`: scan (384 jobs), final training, scoring.

No Python environment is required (the ROOT of the LCG view). Signal MC uses `pThatreweight`, DATA unit
weights; `NormMode=EqualNumEvents` gives both classes equal total importance.

## Small local test

```bash
cd /eos/user/h/hmarques/Analysis_CODES/selectionER/ML_selection/methods/ML_tmva
root -l -b -q 'TMVA_train.C("ntKp/PbPb23","auc",1000,1000,50)'
ML_FEATURES=BtrkLeadPt_dc,BtrkMaxdR_dc ML_OUTPUT_DIR=/tmp/$USER/tmva_test root -l -b -q 'TMVA_train.C("ntmix/PbPb18","auc",3000,3000,200)'
root -l -b -q 'TMVA_optimize.C("ntKp/PbPb23","auc",0,1000,1000)'
```

1,000 signal and 1,000 background candidates and 50 trees give a complete test model, report and
`scores.root` (the input of the feature-study metrics).

## Condor submission from EOS

In every new lxplus shell, select CERN's EOS-aware schedd before any `condor_submit`:

```bash
module load lxbatch/eossubmit
myschedd out
```

## Hyperparameter scan

`../../scripts/fresh_ML_scan.sh tmva ntKp/PbPb23` submits it. By hand:

```bash
mkdir -p condor_logs tmva_outputs/ntKp/PbPb23_auc/optimization
condor_submit -append 'tree = ntKp' -append 'year = PbPb23' -append 'fom = auc' submit_optimize.sub
```

TMVA's test sample holds 40% of the events and is split into the validation and test parts at random.
Each of the 384 jobs trains one combination of the four settings below with 2000 trees, then evaluates
the validation FOMs and the overtraining test at every tree count of the first line on the same forest
(a boosted forest with 2000 trees contains every smaller one): 384 x 14 = 5376 grid points.

```text
NTrees (checkpoints): 25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000
MaxDepth:             3, 4, 5, 6
Shrinkage:            0.001, 0.002, 0.003, 0.005, 0.0075, 0.01
BaggedSampleFraction: 0.4, 0.5, 0.6, 0.7
MinNodeSize:          1%, 1.5%, 2.5%, 5%
nCuts:                20
```

This is grid v2, with the tree counts of v3. The first scan (v1: 500 to 5000 trees, depth 1 to 4,
shrinkage 0.005 to 0.05, bagging 0.6 to 0.9, min node 0.25% to 1.5%) put the X(3872) PbPb23 optimum at
the edge of four axes: 500 trees, depth 4, shrinkage 0.005, min node 1.5%, bagging 0.6. If you change the
grid in `TMVA_config.h`, set `njobs` in `submit_optimize.sub` to MaxDepth x Shrinkage x
BaggedSampleFraction x MinNodeSize values.

With `ML_RANDOM_SEARCH=1` (the feature study), job j instead draws its four settings from the `RANDOM_*`
ranges of `TMVA_config.h` with seed j + 1 (depth 2–10, shrinkage 0.002–0.5 and min node 0.5–10%
log-uniform, bagging 0.3–1.0) and still evaluates every tree count. `ML_SUMMARY_PATH` sets the summary
file. Each job writes one ROOT file under `optimization/`; its `trial` tree has one entry per tree count,
with all FOMs, the KS p-values, the AUC gap, `passed`, and `objective` (the FOM of the run). The final
training uses the best `objective` among the points that pass.

## Final training

```bash
ML_USE_OPTIMIZED=1 root -l -b -q 'TMVA_train.C("ntKp/PbPb23","auc")'
condor_submit -append 'tree = ntKp' -append 'year = PbPb23' -append 'fom = auc' submit_train.sub
```

`ML_OPTIMIZATION_DIR` points it at another folder of summaries (the feature-study refits). Without
`ML_USE_OPTIMIZED=1` (or with `-append 'use_optimized = 0'`), `TMVA_train.C` uses its default arguments
(`NTrees=2000`, `MaxDepth=2`, `Shrinkage=0.02`, `BaggedSampleFraction=0.7`, `MinNodeSize=1.0%`,
`nCuts=20`), which predate grid v2.

`training_report.root` stores the configuration, sample counts, validation and test FOMs, overtraining
tests and working-point efficiencies, ROC curve, confusion matrix, score and train-test-difference
histograms (with the spectators), permutation importance, staged AUC history, and the boost weights.
TMVA has no SHAP values, and ROOT discards the separation-gain statistics when it rebuilds a BDTG from
XML, so the importance is the validation AUC drop when one input is shuffled.

## Scoring

```bash
root -l -b -q 'TMVA_apply.C("ntKp/PbPb23","auc","all")'      # or one kind: DATA, MC, MC_X3872, ...
mkdir -p /eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/ntKp/PbPb23
condor_submit -append 'tree = ntKp' -append 'year = PbPb23' -append 'fom = auc' -queue 'kind in (DATA, MC)' submit_score.sub
```

The outputs go to `selectionER/scored_samples/<tree>/<year>/flat_<tree>_<year>_tmva_<fom>_scored_<KIND>.root`
(`ML_SCORED_DIR` overrides the folder). The train and score jobs read the scan summaries, the model and
its report from the `/eos` mount of the worker (condor on EOS cannot transfer a folder in). Every original
branch is kept and the response added as `Prediction`, in [-1, 1], with the `metadata/` directory of all
three methods (see `selectionER/README.md`); `_dc` inputs are computed from their base branch and `Bmass`
with `TMVA_decorrelation.h`.

## Threshold optimization

From `selectionER/`: `root -l -b -q 'optimalCUT_X_punzi.C("ntKp/PbPb23","tmva","auc")'`, results in
`optimalCUT/ntKp/PbPb23/tmva_auc/`.
