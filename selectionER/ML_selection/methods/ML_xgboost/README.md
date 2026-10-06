# XGBoost workflow

> **Central settings** (cases, features, pre-cuts, sidebands, mass windows, sample sizes, FOM settings)
> are in `ML_selection/settings.py` for all three methods. `../../scripts/fresh_ML_scan.sh` checks the
> inputs, cleans and submits the scans; after the scans, `../../scripts/train_best_ML.sh ntKp/PbPb23 auc`
> trains the best passing trial. Training reads the skims of `../../scripts/make_skims.py`; scoring reads
> the flat files. The feature study (`../../scripts/feature_study.py`) drives the same scripts.

The classifier is trained with the signal MC of a case as signal and its DATA mass sidebands as
background, and applied to every input of the case (X(3872): DATA, prompt and nonprompt X(3872) and
ψ(2S) MC; B mesons: DATA and MC). Every script takes `--case <tree>/<year>` and `--fom`. Paths, tree,
branch and variable names are case-sensitive and used directly: uproot, NumPy, XGBoost or the filesystem
fail when an input is wrong.

## Files

- `XGB_train.py`: final training (outputs in `xgb_outputs/<tree>/<year>_<fom>/`).
- `XGB_optuna.py`: hyperparameter scan.
- `XGB_apply.py`: the `Prediction` branch and the metadata added to every flat file of the case.
- `run_optuna_batch.sh`: sets up the LCG view and runs a Python script on Condor.
- `submit_optuna.sub`: 250 Optuna Condor processes per case.

Signal MC uses `pThatreweight`, DATA unit weights; both classes are normalized to equal total weight
before the train/validation/test split, and early stopping uses the validation split. The saved model
holds the trees up to the best iteration. Model selection (split, FOMs, overtraining test) is shared by
the three methods and described in `selectionER/README.md`.

## Environment on lxplus

Work from this directory (output paths are relative to it), in the LCG view of the node:

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh   # EL9
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt/setup.sh   # EL8
module load lxbatch/eossubmit; myschedd out                                  # before any condor_submit
```

## Optuna

```bash
python XGB_optuna.py --case ntKp/PbPb23 --fom auc --trials 5 --n-rounds 200 --early-stopping 5 \
  --max-signal 1000 --max-background 1000 --seed 0          # small local test
mkdir -p condor_logs xgb_outputs/ntKp/PbPb23_auc/optuna_summaries
condor_submit -append 'tree = ntKp' -append 'year = PbPb23' -append 'fom = auc' submit_optuna.sub
```

Every job writes one summary JSON (`optuna_summary_<tree>_<year>_<fom>_<cluster>_seed<N>.json`); a job
whose trials all fail the overtraining test records `best_value` null.

## Final training

```bash
python XGB_train.py --case ntKp/PbPb23 --fom auc --params-json xgb_outputs/ntKp/PbPb23_auc/optuna_summaries/<best>.json
```

Without `--params-json` the defaults `XGB_PARAMS` are used. `--features` overrides the inputs,
`--output-dir` the output folder. The training also writes `scores.root`, the input of the feature-study
metrics. `training_report.root` stores the configuration, sample statistics, AUC and weighted KS metrics,
best iteration, ROC curve, training history, confusion matrix, SHAP importance, and score distributions.

## Scoring

```bash
python XGB_apply.py --case ntKp/PbPb23 --fom auc --kind all     # or one kind: DATA, MC, MC_X3872, ...
```

The scorer reads each input in `250 MB` chunks, keeps all original branches, adds `Prediction` in
[0, 1], and writes `selectionER/scored_samples/<tree>/<year>/flat_<tree>_<year>_xgb_<fom>_scored_<KIND>.root`
(`--scored-dir` overrides the folder), with the `metadata/` directory of all three methods (see
`selectionER/README.md`). `_dc` inputs are computed from their base branch and `Bmass` with the maps of
the training skims.

## Threshold optimization

From `selectionER/`: `root -l -b -q 'optimalCUT_X_punzi.C("ntKp/PbPb23","xgb","auc")'`, results in
`optimalCUT/ntKp/PbPb23/xgb_auc/`.
