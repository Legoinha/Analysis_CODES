# PyTorch NN workflow

> **Central settings** (cases, features, `LOG_FEATURES`, pre-cuts, sidebands, mass windows, sample sizes,
> FOM settings) are in `ML_selection/settings.py` for all three methods. `../../scripts/fresh_ML_scan.sh`
> checks the inputs, cleans and submits the scans; after the scans,
> `../../scripts/train_best_ML.sh ntKp/PbPb23 auc nn` trains the best passing trial. Training reads the
> skims of `../../scripts/make_skims.py`; scoring reads the flat files. The feature study
> (`../../scripts/feature_study.py`) drives the same scripts.

A fully connected network that separates the signal MC of a case from its DATA mass sidebands, the NN
counterpart of `ML_xgboost` and `ML_tmva` with the same samples, features, cuts and
train/validation/test split. The output is a `Prediction` branch in [0, 1] added to every input of the
case. Every script takes `--case <tree>/<year>` and `--fom`; only the scan runs on Condor.

## Files

| File | Role |
|---|---|
| `NN_train.py` | Network definition and final training (outputs in `nn_outputs/<tree>/<year>_<fom>/`) |
| `NN_optuna.py` | Hyperparameter scan |
| `NN_apply.py` | Adds the `Prediction` branch and the metadata to every flat file of the case |
| `run_nn_batch.sh` | Sets up the LCG view and runs Python |
| `submit_optuna.sub` | 200 Optuna Condor jobs per case |

## Environment

Work from this directory (output paths are relative to it). `run_nn_batch.sh` sources an LCG view with
PyTorch, uproot, Optuna, scikit-learn and matplotlib, and reads inputs through XRootD
(`ML_USE_XROOTD=1`; set `ML_USE_XROOTD=0` for the `/eos` mount). Its default view is for EL9; on an EL8
node set `export NN_LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt` once per shell.
Before any Condor submission: `module load lxbatch/eossubmit; myschedd out`. Run long local steps
(training, scoring DATA) in `tmux`.

## Scan

```bash
mkdir -p condor_logs nn_outputs/ntKp/PbPb23_auc/optuna_summaries
condor_submit -append 'tree = ntKp' -append 'year = PbPb23' -append 'fom = auc' submit_optuna.sub
./run_nn_batch.sh NN_optuna.py --case ntKp/PbPb23 --fom auc --trials 5 --max-epochs 10 --early-stopping 3 \
  --max-signal 3000 --max-background 3000 --seed 0        # small local test
```

Each of the 200 jobs runs 100 trials with its own seed (no new trial after 20 hours, inside the
`tomorrow` flavour) and writes one summary,
`nn_outputs/<tree>/<year>_<fom>/optuna_summaries/optuna_summary_<tree>_<year>_<fom>_<cluster>_seed<N>.json`.
The objective is the validation FOM of `--fom`; trials failing the overtraining test are rejected (see
`selectionER/README.md`), and a median pruner stops trials whose per-epoch validation FOM falls behind.
Scanned space:

```text
hidden_layers   1 - 8
hidden_units    16, 32, 64, 128, 256
activation      relu, silu, elu, tanh
dropout         0.0 - 0.4
learning_rate   1e-4 - 3e-2 (log)
weight_decay    1e-6 - 1e-2 (log)
batch_size      256, 512, 1024, 2048, 4096
```

## Final training

```bash
./run_nn_batch.sh NN_train.py --case ntKp/PbPb23 --fom auc --params-json <best summary> 2>&1 | tee nn_outputs/ntKp/PbPb23_auc/training_log.txt
```

`--params-json` takes the hyperparameters, the epoch settings and the training seed of the scan job, so
the training reproduces the chosen trial; without it the defaults `NN_PARAMS` are used. `--features`
overrides the inputs and `--output-dir` the output folder. Check that the validation FOM matches the best
trial, that the test overtraining test passes (train and test efficiencies at the working point agree,
pulls of order 1), and the PDFs (`score_distributions.pdf`, with the spectators of the case;
`training_history.pdf`). Outputs: `nn_X3872_vs_sideband.pt` (weights, preprocessing, hyperparameters,
features, best epoch), `training_report.root`, `scores.root` (the input of the feature-study metrics),
`roc_curve.pdf`, `feature_importance.pdf` (validation AUC drop when one input is shuffled),
`confusion_matrix.pdf`.

## Scoring

```bash
./run_nn_batch.sh NN_apply.py --case ntKp/PbPb23 --fom auc --kind all     # or one kind: DATA, MC, MC_X3872, ...
```

The scorer keeps all original branches, adds `Prediction`, and writes
`selectionER/scored_samples/<tree>/<year>/flat_<tree>_<year>_nn_<fom>_scored_<KIND>.root`
(`--scored-dir` overrides the folder), with the `metadata/` directory of all three methods. The NN score
scale differs from the BDT score, so a BDT cut value does not carry over. Then, from `selectionER/`:
`root -l -b -q 'optimalCUT_X_punzi.C("ntKp/PbPb23","nn","auc")'`.

## Network

```text
raw features -> log1p(settings.LOG_FEATURES) -> standardize
            -> [Linear -> BatchNorm -> activation -> Dropout] x hidden_layers
            -> Linear -> logit
```

The log transform (by base name, `X` and `X_dc` alike) and the standardization are fitted on the
training split and stored inside the model, so the saved model takes raw branch values; the score is
`sigmoid(logit)`. Signal MC uses `pThatreweight`, DATA unit weights, both classes normalized to the same
total weight; the weights enter the loss, and every AUC and KS value is weighted. Training uses AdamW,
lowers the learning rate when the validation AUC plateaus, and stops early on the weighted validation
AUC, keeping the weights of the best epoch. `NN_train.load_model(path)` loads a trained model. CUDA is
used when available (add `request_gpus = 1` and a CUDA LCG view in `NN_LCG_VIEW`).
