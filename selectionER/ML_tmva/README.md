# ntmix TMVA BDTG workflow

> **Central settings.** Features, pre-cuts, sidebands, and sample sizes are
> set in `selectionER/fresh_ML_scan.sh` for all three ML methods. Edit them
> there, not in this folder. `../fresh_ML_scan.sh --sync-only` writes and checks
> them; `../fresh_ML_scan.sh` also cleans all ML outputs and submits the scans.
> Training reads the skims of `../make_skims.py`; scoring reads the flat files.
> After the scans finish, `../train_best_ML.sh PbPb23 punzi` checks the scan results
> and trains all three methods for that sample with their best passing trial.
> Add `nn`, `xgb`, or `tmva` to train only some. The feature study
> (`../feature_study.py`, see `selectionER/README.md`) drives the same scripts.

This is the ROOT-native reference classifier parallel to `ML_xgboost`.
It trains a TMVA gradient-boosted decision tree with prompt X(3872) MC as
signal and DATA mass sidebands as background. Psi(2S) is only a spectator.

The trained model is applied to DATA, prompt and nonprompt X(3872), and prompt
and nonprompt Psi(2S). Paths, trees, branches, sample names, and score kinds are
case-sensitive.

## Files

- `TMVA_common.h`: helpers shared by training and scan: tree-by-tree forest
  evaluation and exact weighted ROC/AUC.
- `TMVA_config.h`: samples, trees, features, training sidebands, skim and output paths, scan settings.
- `TMVA_decorrelation.h`: the `_dc` transform for scoring, from the maps of the skims.
- `TMVA_train.C`: final training, plots, XML model, and ROOT report.
- `TMVA_optimize.C`: one hyperparameter scan trial.
- `TMVA_apply.C`: apply the XML model and add `Prediction` and the metadata.
- `run_tmva.sh`: EOS-resident launcher required by CERN EosSubmit.
- `submit_optimize.sub`: 384 independent scan jobs, 12 tree counts each.
- `submit_train.sub`: one final training job.
- `submit_score.sub`: five independent scoring jobs.

No Python environment is required. The code uses the ROOT 6.30 TMVA installation
available on lxplus and AlmaLinux 9 batch workers.

## Shared physics configuration

Signal MC uses `pThatreweight`; DATA background uses
unit event weights. `NormMode=EqualNumEvents` gives signal and background equal
total importance in TMVA after weighting.

Configured samples are `PbPb18`, `PbPb23`, `PbPb24`, and `ppRef24`.

## Small local test

```bash
cd /eos/user/h/hmarques/Analysis_CODES/selectionER/ML_tmva

root -l -b -q 'TMVA_train.C("PbPb23","punzi",1000,1000,50)'
```

This uses 1,000 selected signal candidates, 1,000 selected background
candidates, and 50 trees. It creates a complete test model, the report and
`scores.root` (the input of the feature-study metrics). `ML_FEATURES` (comma-separated,
`_dc` names allowed) overrides the inputs and `ML_OUTPUT_DIR` the output folder:

```bash
ML_FEATURES=BtrkLeadPt_dc,BtrkMaxdR_dc ML_OUTPUT_DIR=/tmp/$USER/tmva_test root -l -b -q 'TMVA_train.C("PbPb18","auc",3000,3000,200)'
```

Score DATA with the test model:

```bash
root -l -b -q 'TMVA_apply.C("PbPb23","punzi","DATA")'
```

Score all five inputs locally:

```bash
root -l -b -q 'TMVA_apply.C("PbPb23","punzi","all")'
```

## Condor submission from EOS

This repository and its log directories are on EOS. In every new lxplus shell,
select CERN's EOS-aware schedd before running any `condor_submit` command:

```bash
module load lxbatch/eossubmit
myschedd out
```

`myschedd out` prints the schedd selected from the EosSubmit pool. The exact
hostname can vary. The module setting remains active for the current shell.

## Hyperparameter optimization

Test one scan point locally with small samples:

```bash
root -l -b -q 'TMVA_optimize.C("PbPb23","punzi",0,1000,1000)'
```

Production scan:

```bash
mkdir -p condor_logs
condor_submit -append 'sample = PbPb23' -append 'fom = punzi' -append 'sample_output_dir = tmva_outputs/ntmix_PbPb/pbpb23_punzi' submit_optimize.sub
condor_submit -append 'sample = PbPb18' -append 'fom = punzi' -append 'sample_output_dir = tmva_outputs/ntmix_PbPb/pbpb18_punzi' submit_optimize.sub
```

TMVA's test sample holds 40% of the events and is split into the validation and
test parts at random. Each of the 384 jobs trains one combination of the four settings
below with 1000 trees. It then evaluates the validation FOMs and the overtraining
test (train vs validation) at every tree count in the
first line on the same forest: a boosted forest with 1000 trees contains every
smaller one. So each job gives 13 results, and the scan covers
384 x 12 = 4608 grid points with 384 jobs.

```text
NTrees (checkpoints): 25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000
MaxDepth:             3, 4, 5, 6
Shrinkage:            0.001, 0.002, 0.003, 0.005, 0.0075, 0.01
BaggedSampleFraction: 0.4, 0.5, 0.6, 0.7
MinNodeSize:          1%, 1.5%, 2.5%, 5%
nCuts:                20
```

This is grid v2, with the tree counts of v3 (up to 1000). The first scan (v1: 500 to 5000 trees, depth 1 to 4,
shrinkage 0.005 to 0.05, bagging 0.6 to 0.9, min node 0.25% to 1.5%) put the
PbPb23 optimum at the edge of four axes: 500 trees, depth 4, shrinkage 0.005,
min node 1.5%, bagging 0.6. Depth 1 to 2 and shrinkage 0.015 and above never
came close.

The grid is set in `TMVA_config.h`. If you change it, set `njobs` in
`submit_optimize.sub` to MaxDepth x Shrinkage x BaggedSampleFraction x
MinNodeSize values. A job number outside the grid stops with an error.

With `ML_RANDOM_SEARCH=1` (the feature study), job j instead draws its four settings
from the `RANDOM_*` ranges of `TMVA_config.h` with seed j + 1 (depth 2–8, shrinkage
0.002–0.2 and min node 0.5–10% log-uniform, bagging 0.3–1.0) and still evaluates
every tree count. `ML_SUMMARY_PATH` sets the summary file.

Each job writes one ROOT file under the sample's `optimization/` directory. Its
`trial` tree has one entry per tree count, with all three FOMs, the KS p-values,
the AUC gap, `passed`, and `objective` (the FOM of the run). Final training uses
the best `objective` among the points that pass. Model selection (train/validation/test split, the FOM a scan maximizes, the overtraining
pass/fail test) is shared by the three methods and described in `selectionER/README.md`. Jobs read
the skims, so almost all their time is training.

## Final training

Use the best available optimization summary locally:

```bash
ML_USE_OPTIMIZED=1 root -l -b -q 'TMVA_train.C("PbPb23","punzi")'
```

`ML_OPTIMIZATION_DIR` points it at another folder of summaries (the feature-study
refits). Without `ML_USE_OPTIMIZED=1`, `TMVA_train.C` uses its fixed default arguments:
`NTrees=2000`, `MaxDepth=2`, `Shrinkage=0.02`,
`BaggedSampleFraction=0.7`, `MinNodeSize=1.0%`, and `nCuts=20`. They predate
grid v2.

Submit final training after all optimization jobs have finished:

```bash
mkdir -p condor_logs
condor_submit -append 'sample = PbPb23' -append 'fom = punzi' -append 'sample_output_dir = tmva_outputs/ntmix_PbPb/pbpb23_punzi' submit_train.sub
```

`submit_train.sub` uses the best completed scan point automatically. To train
with the fixed defaults instead:

```bash
condor_submit -append 'sample = PbPb23' -append 'fom = punzi' -append 'sample_output_dir = tmva_outputs/ntmix_PbPb/pbpb23_punzi' -append 'use_optimized = 0' submit_train.sub
```

## Final-training outputs

`tmva_training.root` is the native TMVA output. `training_report.root` stores
the configuration metadata, selected sample counts, validation and test FOMs,
overtraining tests and working-point efficiencies, weighted
signal and background KS, selection-threshold metric, ROC curve, balanced
confusion matrix, score and train-test-difference histograms, BDT variable
importance, staged AUC history, and per-tree boost weights.

TMVA does not provide XGBoost-style SHAP values. ROOT 6.30 also discards the
training-only separation-gain statistics when its Factory reconstructs BDTG
from XML, which made the old importance output `NaN`. The
`feature_importance.pdf` plot and `feature_importance` report histogram now
contain finite split-frequency importance percentages, ordered and displayed
like the XGBoost cumulative-importance plot.

`training_history.pdf` now shows train and validation AUC versus boosting round,
evaluated at up to 100 checkpoints on all training and validation events, on the
same 0.75 to 1 range as the XGB and NN histories. The raw gradient-boost weights (constant for BDTG) remain
available in the report's `boost_weights` tree.

## Scoring submission

After final training succeeds:

```bash
mkdir -p condor_logs
condor_submit -append 'sample = PbPb23' -append 'fom = punzi' -append 'sample_output_dir = tmva_outputs/ntmix_PbPb/pbpb23_punzi' submit_score.sub
```

The train and score jobs read the scan summaries, the model and its report from the `/eos`
mount of the worker (condor on EOS cannot transfer a folder in), so `sample_output_dir` must
match `sample` and `fom`.

The five outputs are written under `scored_samples/`:

```text
flat_ntmix_<system>_<fom>_scored_DATA.root
flat_ntmix_<system>_<fom>_scored_MC_X3872.root
flat_ntmix_<system>_<fom>_scored_MC_X3872_NONPROMPT.root
flat_ntmix_<system>_<fom>_scored_MC_PSI2S.root
flat_ntmix_<system>_<fom>_scored_MC_PSI2S_NONPROMPT.root
```

Every original branch is preserved and the TMVA response is added as
`Prediction`, in [-1, 1]. Every scored file also has a `metadata/` directory of
`TObjString`s, with the same keys for all three ML methods: `classifier`,
`sample`, `kind`, `tree`, `input_file`, `model`, `features`, `pre_cut` (the
pre-cut of the training that wrote the model, read from its
`training_report.root`), `mc_weight`, `score_branch`, `score_min`, `score_max`,
and `decorrelation_maps` when the model has `_dc` inputs. Those are computed from
their base branch and `Bmass` with `TMVA_decorrelation.h`.

## Threshold optimization

From `selectionER/`:

```bash
root -l -b -q 'optimalCUT_X_punzi.C("PbPb23","tmva","punzi")'
```

It reads the scored DATA and prompt X(3872) MC from `scored_samples/`, takes
the pre-cut and the [-1, 1] score range from their metadata, and writes its
results to `optimalCUT/tmva/PbPb23_punzi/`.
