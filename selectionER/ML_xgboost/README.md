# ntmix XGBoost workflow

> **Central settings.** Features, pre-cuts, sidebands, and sample sizes are
> set in `selectionER/fresh_ML_scan.sh` for all three ML methods. Edit them
> there, not in this folder. `../fresh_ML_scan.sh --sync-only` writes and checks
> them; `../fresh_ML_scan.sh` also cleans all ML outputs and submits the scans.
> Training reads the skims of `../make_skims.py`; scoring reads the flat files.
> After the scans finish, `../train_best_ML.sh PbPb23 punzi` checks the scan results
> and trains all three methods for that sample with their best passing trial.
> Add `nn`, `xgb`, or `tmva` to train only some. The feature study
> (`../feature_study.py`, see `selectionER/README.md`) drives the same scripts.

The classifier is trained with prompt X(3872) MC as signal and DATA mass
sidebands as background. Psi(2S) is not used as a training class.

The trained model is applied to:

- DATA
- prompt X(3872) MC
- nonprompt X(3872) MC
- prompt Psi(2S) MC
- nonprompt Psi(2S) MC

Paths, tree names, branch names, sample names, and variable names are
case-sensitive. The scripts use them directly and let uproot, NumPy, XGBoost,
or the filesystem fail when an input is wrong.

## Files

- `XGB_train.py`: sample configuration and final training.
- `XGB_optuna.py`: hyperparameter scan.
- `XGB_apply.py`: add the `Prediction` branch and the metadata to the five samples.
- `run_optuna_batch.sh`: sets up the LCG view and runs a Python script (e.g. `XGB_optuna.py`) on Condor.
- `submit_optuna.sub`: 250 Optuna Condor processes per sample.

Final training and scoring run locally.

## Physics configuration

Signal MC uses `pThatreweight` as its event weight. DATA background uses unit
weights. Signal and background weights are normalized to equal total class
weight before the train/validation/test split. Early stopping uses the validation
split.

Model selection (train/validation/test split, the FOM a scan maximizes, the overtraining
pass/fail test) is shared by the three methods and described in `selectionER/README.md`.

All input paths are defined once in `FLAT_INPUTS` of `../ML_common/ml_samples.py`,
together with the pre-cuts, training sidebands and production features. PbPb24
paths follow the expected naming convention and can be used when its MC files
become available.

Configured samples are `PbPb18`, `PbPb23`, `PbPb24`, and `ppRef24`. The PbPb18
inputs are the reflated legacy Run-2 samples under:

```text
/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/
```

## Environment on lxplus

Work from this directory; output paths are relative to it. Source the LCG view
that matches the node (`cat /etc/redhat-release`), once per shell:

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh   # EL9
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt/setup.sh   # EL8
```

This repository and its log directories are on EOS. Before any Condor
submission in a new lxplus shell, select CERN's EOS-aware schedd:

```bash
module load lxbatch/eossubmit
myschedd out
```

`myschedd out` prints the schedd selected from the EosSubmit pool. The exact
hostname can vary.

## Optuna

Small local test:

```bash
python XGB_optuna.py \
  --sample PbPb23 \
  --fom punzi \
  --trials 5 \
  --n-rounds 200 \
  --early-stopping 5 \
  --max-signal 1000 \
  --max-background 1000 \
  --summary-subdir optuna_summaries \
  --seed 0
```

The full Condor scan is submitted by `../fresh_ML_scan.sh xgb`. By hand:

```bash
module load lxbatch/eossubmit
mkdir -p condor_logs
condor_submit -append 'sample = PbPb23' -append 'fom = punzi' -append 'sample_output_dir = xgb_outputs/ntmix_PbPb/pbpb23_punzi' submit_optuna.sub
condor_submit -append 'sample = PbPb18' -append 'fom = punzi' -append 'sample_output_dir = xgb_outputs/ntmix_PbPb/pbpb18_punzi' submit_optuna.sub
```

Every job writes one summary JSON; a job whose trials all fail the overtraining
test records `best_value` null. `../train_best_ML.sh` picks the best summary.

## Final training

```bash
python XGB_train.py --sample PbPb23 --fom punzi --params-json xgb_outputs/ntmix_PbPb/pbpb23_punzi/optuna_summaries/<best>.json
```

Without `--params-json` the defaults `XGB_PARAMS` are used. `--features` overrides
the inputs, `--output-dir` the output folder. The training also writes
`scores.root`, the input of the feature-study metrics (`selectionER/README.md`).

`training_report.root` stores the configuration metadata, pThat/sample
statistics, AUC and weighted KS metrics, best iteration, ROC curve, training
history, confusion matrix, SHAP importance, and score distributions.

## Scoring

Score one file:

```bash
python XGB_apply.py --sample PbPb23 --fom punzi --kind DATA
python XGB_apply.py --sample PbPb23 --fom punzi --kind MC_X3872
python XGB_apply.py --sample PbPb23 --fom punzi --kind MC_X3872_NONPROMPT
python XGB_apply.py --sample PbPb23 --fom punzi --kind MC_PSI2S
python XGB_apply.py --sample PbPb23 --fom punzi --kind MC_PSI2S_NONPROMPT
```

Score all five:

```bash
python XGB_apply.py --sample PbPb23 --fom punzi --kind all
```

The scorer reads each input in `250 MB` chunks, preserves all original
branches, and adds `Prediction` in [0, 1]. Outputs are written under
`scored_samples/` as:

```text
flat_ntmix_<system>_<fom>_scored_DATA.root
flat_ntmix_<system>_<fom>_scored_MC_X3872.root
flat_ntmix_<system>_<fom>_scored_MC_X3872_NONPROMPT.root
flat_ntmix_<system>_<fom>_scored_MC_PSI2S.root
flat_ntmix_<system>_<fom>_scored_MC_PSI2S_NONPROMPT.root
```

Every scored file also has a `metadata/` directory of `TObjString`s, with the
same keys for all three ML methods: `classifier`, `sample`, `kind`, `tree`,
`input_file`, `model`, `features`, `pre_cut` (the pre-cut of the training
that wrote the model, read from its `training_report.root`), `mc_weight`,
`score_branch`, `score_min`, `score_max`, and `decorrelation_maps` when the
model has `_dc` inputs. Those are computed from their base branch and `Bmass`
with the maps of the training skims.

## Threshold optimization

From `selectionER/`:

```bash
root -l -b -q 'optimalCUT_X_punzi.C("PbPb23","xgb","punzi")'
```

It reads the scored DATA and prompt X(3872) MC from `scored_samples/`, takes
the pre-cut and score range from their metadata, and writes its results to
`optimalCUT/xgb/PbPb23_punzi/`.
