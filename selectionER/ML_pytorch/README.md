# ntmix PyTorch NN workflow

A fully connected neural network that separates X(3872) signal from
combinatorial background. It is the NN counterpart of `ML_xgboost` and
`ML_tmva`, and uses the same samples, features, cuts, and train/validation/test split.

- **Signal:** prompt X(3872) MC.
- **Background:** DATA in the B-mass sidebands.
- **Output:** a `Prediction` branch in [0, 1] added to DATA and to the four MC
  samples (prompt and nonprompt X(3872), prompt and nonprompt Psi(2S)).

This README takes you through the whole chain, from the hyperparameter scan to
the scored ROOT files. Follow the steps in order.

---

## The chain at a glance

Run each sample (for example PbPb23, then PbPb18) through steps 3 to 5.

| Step | What | Where it runs | Command |
|---|---|---|---|
| 0 | Set up the environment | lxplus shell | see below |
| 1 | Optuna hyperparameter scan | Condor, 200 jobs per sample | `condor_submit ... submit_optuna.sub` |
| 2 | Check the scan finished | lxplus | log checks and ranking script |
| 3 | Choose the best summary | lxplus | ranking script, or `../train_best_ML.sh` |
| 4 | Final training and metrics | lxplus | `NN_train.py --sample X --params-json <summary>` |
| 5 | Score the five samples | lxplus | `NN_apply.py --sample X --kind all` |
| 6 | Re-optimize the cut on `Prediction` | ROOT | `optimalCUT_X_punzi.C("PbPb23","nn","punzi")` |

**Features, cuts, and sample sizes are set centrally** in
`selectionER/fresh_ML_scan.sh`, shared with `ML_xgboost` and `ML_tmva`, and written
by it into `../ML_common/ml_samples.py` (`LOG_FEATURES` and sample sizes into
`NN_train.py`). Do not edit them there: the script overwrites them. Training reads
the skims of `../make_skims.py`. Run
`../fresh_ML_scan.sh --sync-only` to write and check new settings, or
`../fresh_ML_scan.sh` to also clean all outputs and submit the scans of all
three methods (this replaces step 1 below). After the scans finish,
`../train_best_ML.sh PbPb23 punzi` checks the results and runs steps 2 to 4 for all
three methods; add `nn` to train only the network. The feature study
(`../feature_study.py`) drives the same scripts.

Only step 1 runs on Condor. Steps 4 and 5 run locally.

---

## Step 0: environment

Always work from this directory. Output paths in the scripts are relative to
it.

```bash
cd /eos/user/h/hmarques/Analysis_CODES/selectionER/ML_pytorch
```

No virtual environment is needed. `run_nn_batch.sh` sources an LCG view that
provides PyTorch, uproot, Optuna, scikit-learn, and matplotlib, then runs
Python. Every command below goes through it.

The default view is for EL9 (AlmaLinux 9) nodes. Check your node with
`cat /etc/redhat-release`. On an EL8 node, select the EL8 view first, once per
shell:

```bash
export NN_LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el8-gcc11-opt
```

`run_nn_batch.sh` reads inputs through XRootD (`ML_USE_XROOTD=1`). Set
`ML_USE_XROOTD=0` to read through the `/eos` mount instead.

Before any Condor submission, select CERN's EOS-aware scheduler, once per
shell:

```bash
module load lxbatch/eossubmit
myschedd out
```

Long local steps (training, scoring DATA) should run inside `tmux` or `screen`,
so they survive a dropped connection.

---

## Step 1: Optuna hyperparameter scan (Condor)

Each sample gets 200 Condor jobs. Each job runs 100 trials with its own random
seed and writes one JSON summary of its best trial.

```bash
module load lxbatch/eossubmit
mkdir -p condor_logs \
         nn_outputs/ntmix_PbPb/pbpb23_punzi/optuna_summaries \
         nn_outputs/ntmix_PbPb/pbpb18_punzi/optuna_summaries

condor_submit -append 'sample = PbPb23' -append 'fom = punzi' -append 'sample_output_dir = nn_outputs/ntmix_PbPb/pbpb23_punzi' submit_optuna.sub
condor_submit -append 'sample = PbPb18' -append 'fom = punzi' -append 'sample_output_dir = nn_outputs/ntmix_PbPb/pbpb18_punzi' submit_optuna.sub
```

Follow the jobs with `condor_q`. A job stops starting new trials after 20
hours, so the scan finishes within the 24-hour `tomorrow` flavour.

The summaries land in:

```text
nn_outputs/ntmix_PbPb/pbpb23_punzi/optuna_summaries/optuna_summary_PbPb23_punzi_<ClusterId>_seed<N>.json
```

Scan settings are in `submit_optuna.sub`: 100 trials per job, up to 150 epochs,
early stopping after 15 epochs without improvement, and at most 150k signal and
500k background candidates.

<details>
<summary>What the scan optimizes</summary>

The objective is the validation FOM chosen with `--fom`; trials failing the
overtraining test are rejected (see `selectionER/README.md`). A median pruner also
stops trials whose per-epoch validation FOM falls behind earlier trials. Pruned
trials are expected and are not failures.

Scanned space:

```text
hidden_layers   1 - 5
hidden_units    16, 32, 64, 128, 256
activation      relu, silu, elu, tanh
dropout         0.0 - 0.4
learning_rate   1e-4 - 1e-2 (log)
weight_decay    1e-6 - 1e-2 (log)
batch_size      256, 512, 1024, 2048, 4096
```

</details>

---

## Step 2: check that the scan finished

**Count finished jobs.** Replace `<ClusterId>` with the cluster number of your
submission. Every job should report return value 0.

```bash
grep -c "Job terminated" condor_logs/optuna_PbPb23_<ClusterId>.log
grep -h "return value" condor_logs/optuna_PbPb23_<ClusterId>.log | sort | uniq -c
```

**List seeds without a summary.** A few missing seeds are harmless. Those jobs
are usually still running or were evicted. Check them with `condor_q`.

```bash
for i in $(seq 0 199); do
  [ -f nn_outputs/ntmix_PbPb/pbpb23_punzi/optuna_summaries/optuna_summary_PbPb23_punzi_<ClusterId>_seed$i.json ] || echo -n "$i "
done; echo
```

**Look for misplaced summaries.** Occasionally a JSON file ends up in a nested
`optuna_summaries/nn_outputs/...` folder. `train_best_ML.sh` does not look there, so
move any such file up one level into `optuna_summaries/`:

```bash
find nn_outputs -path '*optuna_summaries/nn_outputs*' -name '*.json'
```

**Rank the results.** This plain `python3` snippet needs no LCG view:

```bash
python3 - <<'PY'
import glob, json, statistics
sample, folder = "PbPb23", "nn_outputs/ntmix_PbPb/pbpb23_punzi/optuna_summaries"
rows = []
for path in glob.glob(folder + "/*.json"):
    s = json.load(open(path))
    if s["best_value"] is None: continue  # every trial of this job failed the overtraining test
    u = s["best_trial"]["user_attrs"]; p = s["best_params"]
    rows.append((s["best_value"], s["seed"], u["auc"], u["train_auc"], min(u["ks_signal_pvalue"], u["ks_background_pvalue"]), p, path))
rows.sort(key=lambda r: -r[0])
values = [r[0] for r in rows]
print(f"{sample}: {len(rows)} summaries with a passing trial, best {values[0]:.5f}, median {statistics.median(values):.5f}")
for v, seed, val, train, ks, p, path in rows[:5]:
    print(f"{v:.5f} seed{seed:<4} val_auc={val:.4f} train_auc={train:.4f} min KS p={ks:.3f} "
          f"{p['hidden_layers']}x{p['hidden_units']} {p['activation']} lr={p['learning_rate']:.1e} bs={p['batch_size']} {path}")
PY
```

A healthy scan has a flat top: the best few results differ only in the fourth
decimal, and train and validation AUC are close.

---

## Step 3: choose the best summary

The first line of the ranking is the summary to train from; `../train_best_ML.sh`
makes the same choice and also checks that the scan used the current features,
pre-cut, sidebands and sample sizes.

---

## Step 4: final training and metrics

```bash
./run_nn_batch.sh NN_train.py --sample PbPb23 --fom punzi --params-json <best summary> 2>&1 | tee nn_outputs/ntmix_PbPb/pbpb23_punzi/training_log.txt
```

`--params-json` takes the hyperparameters, the epoch settings and the training
seed of the scan job, so the training reproduces the chosen trial. Without it the
defaults `NN_PARAMS` are used. `--features` overrides the inputs and `--output-dir`
the output folder. The training also writes `scores.root`, the input of the
feature-study metrics (`selectionER/README.md`).

This takes a few minutes on CPU. The last lines of the log give the result:

The `Validation:` and `Test:` lines at the end of the log give the FOMs, the
overtraining test, and the working-point efficiencies.

**What to check:**

- **The validation FOM** should match the best Optuna trial.
- **The test overtraining test** should pass, and the train and test signal and
  background efficiencies at the working point should agree within their
  uncertainties (pull of order 1).
- **The PDFs.** Look at `score_distributions.pdf` for train/test agreement and
  at `training_history.pdf` for a smooth loss curve.

Outputs in `nn_outputs/ntmix_PbPb/pbpb23_punzi/`:

| File | Content |
|---|---|
| `nn_X3872_vs_sideband.pt` | Model weights, preprocessing, hyperparameters, feature list, best epoch |
| `training_report.root` | Configuration, sample statistics, AUC and KS metrics, ROC curve, per-epoch history, confusion matrix, feature importance, score distributions |
| `training_log.txt` | The console output captured by `tee` |
| `training_history.pdf` | Loss, AUC, and learning rate per epoch |
| `roc_curve.pdf` | ROC curve |
| `score_distributions.pdf` | Train and test score distributions for signal and background |
| `feature_importance.pdf` | Permutation importance: the drop in test AUC when one feature is shuffled |
| `confusion_matrix.pdf` | Confusion matrix on the test split at the working point (Punzi-optimal threshold of the validation split) |

The model file stores its own hyperparameters, features and preprocessing.

---

## Step 5: score the samples

```bash
./run_nn_batch.sh NN_apply.py --sample PbPb23 --fom punzi --kind all
```

This loads `nn_X3872_vs_sideband.pt` from the sample's output folder and scores
all five samples. The DATA file is several GB, so this step takes the longest.
To score one sample at a time, replace `all` with `DATA`, `MC_X3872`,
`MC_X3872_NONPROMPT`, `MC_PSI2S`, or `MC_PSI2S_NONPROMPT`.

The scorer keeps all original branches, adds `Prediction`, and writes to
`scored_samples/`:

```text
flat_ntmix_PbPb23_punzi_scored_DATA.root
flat_ntmix_PbPb23_punzi_scored_MC_X3872.root
flat_ntmix_PbPb23_punzi_scored_MC_X3872_NONPROMPT.root
flat_ntmix_PbPb23_punzi_scored_MC_PSI2S.root
flat_ntmix_PbPb23_punzi_scored_MC_PSI2S_NONPROMPT.root
```

The names match the XGBoost and TMVA outputs, so downstream code only needs to
point at this folder.

Every scored file also has a `metadata/` directory of `TObjString`s, with the
same keys for all three ML methods: `classifier`, `sample`, `kind`, `tree`,
`input_file`, `model`, `features`, `pre_cut` (the pre-cut of the training
that wrote the model, read from its `training_report.root`), `mc_weight`,
`score_branch`, `score_min`, `score_max`, and `decorrelation_maps` when the
model has `_dc` inputs. Those are computed from their base branch and `Bmass`
with the maps of the training skims.

---

## Step 6: next sample, then downstream

Repeat steps 3 to 5 for the next sample:

```bash
./run_nn_batch.sh NN_train.py --sample PbPb18 --fom punzi --params-json <best PbPb18 summary> 2>&1 | tee nn_outputs/ntmix_PbPb/pbpb18_punzi/training_log.txt
./run_nn_batch.sh NN_apply.py --sample PbPb18 --fom punzi --kind all
```

Then re-optimize the working point on `Prediction`, from `selectionER/`:

```bash
root -l -b -q 'optimalCUT_X_punzi.C("PbPb23","nn","punzi")'
```

The macro reads the scored DATA and prompt X(3872) MC from `scored_samples/`,
takes the pre-cut and score range from their metadata, and writes its results
to `optimalCUT/nn/PbPb23_punzi/`. The NN score scale differs from the BDT score, so a
BDT cut value does not carry over.

---

## Reference results

The first scan of September 2026 (clusters 341026 and 341027) gave the
following. It used the previous configuration: 4 features (`Btrk1dR`,
`Btrk2Pt`, `Btrk1Pt`, `Bchi2Prob`), the pre-cut Bchi2Prob > 0.1 instead of
Btrk2dR < 0.35, and sidebands at 3.75 to 3.80 and 3.95 to 4.00 GeV. Results with
the current configuration will differ. Replace this table after the new scan.

| Sample | Summaries | Best score | Median score | Best test AUC | Best trial |
|---|---|---|---|---|---|
| PbPb23 | 196 of 200 | 0.89645 | 0.89499 | 0.9021 | seed 85: 5 x 128, relu, dropout 0.034, lr 4.2e-3, wd 1.5e-6, batch 512 |
| PbPb18 | 198 of 200 | 0.86279 | 0.86101 | 0.8729 | seed 117: 5 x 64, silu, dropout 0.020, lr 5.2e-3, wd 3.5e-5, batch 512 |

No trial failed. About two thirds of the trials were pruned, which is normal.

---

## Known issues

- **Run from this directory.** Output folders are relative paths. Running from
  elsewhere writes outputs to the wrong place or fails to find the model.
- **EL9 view on an EL8 node.** The default LCG view fails on EL8. Set
  `NN_LCG_VIEW` as in step 0.

---

## Small local test

Before a full scan, check that the code runs end to end on a few thousand
candidates. These commands overwrite the sample's model and plots, so rerun the
full training afterwards.

```bash
./run_nn_batch.sh NN_train.py --sample PbPb23 --fom punzi --max-signal 5000 --max-background 5000 --max-epochs 12

./run_nn_batch.sh NN_optuna.py --sample PbPb23 --fom punzi --trials 5 --max-epochs 10 --early-stopping 3 \
  --max-signal 3000 --max-background 3000 --summary-subdir optuna_summaries --seed 0
```

---

## Background

### Files

| File | Role |
|---|---|
| `NN_train.py` | Sample configuration, network definition, final training |
| `NN_optuna.py` | Hyperparameter scan |
| `NN_apply.py` | Adds the `Prediction` branch and the metadata to the five samples |
| `run_nn_batch.sh` | Sets up the LCG view and runs Python |
| `submit_optuna.sub` | 200 Optuna Condor jobs per sample |

### Samples

All input paths are defined once in `FLAT_INPUTS` of `../ML_common/ml_samples.py`.
Configured samples are `PbPb18`, `PbPb23`, `PbPb24`, and `ppRef24`. Paths, tree
names, branch names, and variable names are case-sensitive.

### Physics configuration

Identical to `ML_xgboost`: same samples, trees, features, pre-cuts, sideband
definition, subsampling, and 60/20/20 stratified train/validation/test split, all
from `../ML_common/ml_samples.py`.

- **Features:** `ml_samples.FEATURES` (production), or `--features`. The
  `LOG_FEATURES` of `NN_train.py` get a log transform, by base name (`X` and `X_dc`
  alike).
- **Pre-cuts:** 15 < Bpt < 50 GeV, |By| < 2.4, BQvalue < 0.15, and Bchi2Prob > 0.05
  (PbPb23) or > 0.10 (PbPb18).
- **Background sidebands:** 15–50 MeV from m_X on both sides, 3.82164–3.85664 and 3.88664–3.92164 GeV.
- **Sample sizes:** at most 165k signal and 500k background candidates.
- **Weights:** signal MC uses `pThatreweight`, DATA uses unit weights. Both
  classes are normalized to the same total weight before the split. The weights
  enter the loss, and every AUC and KS value is weighted.

### Network

```text
raw features -> log1p(LOG_FEATURES) -> standardize
            -> [Linear -> BatchNorm -> activation -> Dropout] x hidden_layers
            -> Linear -> logit
```

The log transform and the standardization are fitted on the training split and
stored inside the model, so the saved model takes raw branch values. The score
is `sigmoid(logit)`.

Training uses AdamW, lowers the learning rate when the validation AUC plateaus,
and stops early on the weighted validation AUC, keeping the weights of the best
epoch. The test split is only used for the final report.

A trained model can be loaded in Python with `NN_train.load_model(path)`.

### GPU

The code uses CUDA automatically when it is available. The default view and
submit files are CPU-only, which is enough for four input features. To use a
GPU, add `request_gpus = 1` to a submit file and point `NN_LCG_VIEW` at a
CUDA-enabled LCG view.
