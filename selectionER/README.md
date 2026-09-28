# selectionER: ML selection and threshold optimization

Three classifiers (`ML_pytorch` NN, `ML_xgboost` XGB, `ML_tmva` BDTG) separate
prompt X(3872) MC (signal) from DATA in the B-mass sidebands (background). They
share samples, features, pre-cuts, sample sizes, and the model-selection rules
below. Each method folder has its own README with the method-specific commands.

## The chain

```bash
python3 make_skims.py --sample PbPb23                                 # once per pre-cut: skims and decorrelation maps
./feature_study.py run                                                 # optional: choose the inputs (below)
./fresh_ML_scan.sh [auc|sigeff|pauc|punzi] [nn xgb tmva] [PbPb18 PbPb23]   # sync settings, check inputs, clean, submit scans
./train_best_ML.sh PbPb23 auc [nn xgb tmva]                           # check scans, train the best trial
# score: the commands are printed at the end of train_best_ML.sh
root -l -b -q 'optimalCUT_X_punzi.C("PbPb23","xgb","auc")'            # Punzi threshold on the scored samples
```

All settings (features, pre-cuts, training sidebands, feature-study pool, sample sizes,
FOM settings) are edited only in the settings block of `fresh_ML_scan.sh`. Its step 1
writes them into `ML_common/ml_samples.py` (NN and XGB), `ML_tmva/TMVA_config.h`,
`ML_common/ml_fom.py` and `optimalCUT_X_punzi.C`; step 2 checks every input tree
(branches exist, cuts parse, ROOT and the NumPy translation used by NN and XGB select
the same candidates) and that the skims match the current pre-cut and flat files.

## Skims and mass decorrelation (`make_skims.py`, `ML_common/ml_decorrelation.py`)

Every training and scan reads the skims in `skims/<sample>/`: the flat trees after the
sample's pre-cut, with all branches, plus `<feature>_dc` for every feature a model may
use (production `FEATURES`, anchor and pool). Scoring (`*_apply`) reads the full flat
files and computes the `_dc` values itself.

`<feature>_dc` is the feature with its background dependence on the candidate mass
removed by quantile morphing. On DATA in 3.60–4.00 GeV outside the X(3872) region
(3.85–3.89) and ±20 MeV around the ψ(2S), the background percentiles (0.5%, 1.5%, …,
99.5%) are measured in 10 MeV slices and fitted versus mass (cubic in m − m_X). A
candidate's value is then mapped from the percentiles at its mass to those at m_X.
Background then has the same distribution at every mass, and at m_X the map is the
identity. `make_skims.py` checks this on the skims it writes:
- `decorrelation_check_<sample>.json` holds the feature–mass correlation before and after,
  and the closure (fraction below each percentile minus the level) in the fit region and in
  the X gaps the fit never saw;
- `decorrelation_check_<sample>.pdf` shows the correlations of all features, then per feature
  the background 10/50/90% percentiles vs Bmass before and after (flat after), and the signal
  MC distribution before and after (almost unchanged).

`python3 make_skims.py --sample PbPb18 --checks-only` redoes the check on the existing skims. The maps (`decorrelation_<sample>.root`) are shared by Python
and C++ (`ML_tmva/TMVA_decorrelation.h`). A model with `_dc` inputs records them in
the `decorrelation_maps` metadata of its scored files.

## Runs are tagged with their FOM

The FOM a scan maximizes (`auc`, `sigeff`, `pauc`, `punzi`) is part of the run name, so runs
with different FOMs are kept side by side:

```text
ML_<method>/<method>_outputs/ntmix_PbPb/pbpb23_punzi/          scans, model, training report
ML_<method>/condor_logs/optuna_PbPb23_punzi_<cluster>_<job>.*  logs (optimize_* for TMVA)
ML_<method>/scored_samples/flat_ntmix_PbPb23_punzi_scored_<KIND>.root
optimalCUT/<nn|xgb|tmva>/PbPb23_punzi/                        Punzi results
```

## Model selection (`ML_common/ml_fom.py`, `ML_tmva/TMVA_common.h`)

**Split.** 60% training, 20% validation, 20% test, stratified by class. The
validation split chooses the model: early stopping (NN, XGB), the number of trees
(TMVA), and the hyperparameters of every scan. The test split is used once, in the
final training report. For TMVA, TMVA's test sample (40%) is split into the
validation and test parts at random (fixed seed).

**Figures of merit**, all on the weighted validation ROC curve:

- `auc`: the ROC area.
- `sigeff`: the signal efficiency at background efficiency `BKG_EFF_TARGET` (0.02,
  close to the background efficiency of the Punzi working point).
- `pauc`: the ROC area for background efficiency from 0 to `PAUC_MAX_BKG_EFF` (0.05),
  divided by that width (0 to 1). It covers the region around the working point but
  averages over many ROC points, so it fluctuates less than `sigeff`.
- `punzi`: the best `e_S / S_min(B)` over thresholds, with the Punzi `a`, `b` of the
  settings, computed as `optimalCUT_X_punzi.C` computes it. `e_S` is the signal efficiency;
  `B` is the background left in the window after the cut: the DATA candidates in the Punzi
  sidebands that pass it, scaled to the window width. In the scans and step comparisons it is
  evaluated on the validation split, where B is the whole-DATA sideband count before the cut
  times the passing fraction of the split's sideband candidates.

  **Window and sidebands** (settings `SIGNAL_HALF_WIDTH`, `SIDEBAND_START`, `SIDEBAND_WIDTH`):
  |m − m_X| < 5 MeV, background from 15–35 MeV on both sides scaled by 10/40. The signal MC
  peak holds 68% within 6 MeV, so the window keeps about 60% of it, and only about 2% of the
  window signal falls in the sidebands. The sidebands lie inside the training sidebands
  (15–50 MeV from m_X), so their background is counted on all DATA, the training candidates
  included, exactly as `optimalCUT_X_punzi.C` counts it on the full DATA. The maximum is taken only over cuts that keep at least
  `PUNZI_MIN_BACKGROUND` (10) background candidates. Otherwise a cut that leaves no
  validation background, and so gets B = 0, would win by chance: that happened for TMVA at
  1.7% signal efficiency. `optimalCUT_X_punzi.C` applies the same window, sidebands and
  minimum (to the DATA sidebands).

  **Threshold grid.** The Punzi FOM is scanned on thresholds k × `THRESHOLD_STEP` (0.02), the
  same grid in the model selection (Python and TMVA) and in `optimalCUT_X_punzi.C`. Scanning
  every candidate's score would let the maximum land on an upward fluctuation of the curve.
  `refit/punzi.pdf` shows the validation curve with its bootstrap ±1σ band: a training curve
  inside the band differs from it by statistics only.

**Overtraining is pass/fail, not a penalty.** A model passes when the
train-validation AUC gap is at most `MAX_AUC_GAP` and the weighted KS p-values of
the signal and background scores (train vs validation, effective sample sizes
`(sum w)^2 / sum w^2`) are both at least `KS_PVALUE_MIN`. Failing trials are
rejected (Optuna: pruned; TMVA: `passed = 0`) and are never chosen. A scan job in
which every trial fails still writes its summary, with `best_value` null.

**Working point.** The working point is the Punzi maximum on all signal MC (training, validation
and test) and all DATA of the model, the threshold and value `optimalCUT_X_punzi.C` finds for it:
the same formula, window, sidebands, threshold grid and 10-candidate minimum. The training report compares the signal and background
efficiencies at that threshold between training and test, with binomial
uncertainties for weighted events and their pull: overtraining in the high-score
tail, which biases the efficiency of the analysis cut, shows up there even when the
global AUC gap is small.

The `Validation:` and `Test:` lines of the training log and the `metrics` of
`training_report.root` (`validation_*`, `test_*`) hold all these numbers.

**Scores file.** Every training also writes `scores.root`, identical for the three
methods (`ml_fom.write_scores`):
- `train`, `validation`, `test`: y, w, score, bmass and idx (the skim entry);
- `sculpt`: the DATA never used for training;
- `spectator_<key>`: nonprompt X(3872), prompt and nonprompt ψ(2S) MC;
- `metadata_json`.

`ml_fom.study_metrics` computes everything the feature study compares from it:
- the four FOMs with bootstrap errors, the overtraining test, and the working-point
  efficiencies;
- the test split at the validation working point;
- the spectator efficiencies;
- the sculpting check.

**Sculpting check.** The background efficiency versus Bmass on the untrained DATA, at
the working point and at twice and half its background efficiency:
- a straight line fitted in 5 MeV bins over 3.75–3.85 and 3.89–4.00 GeV (slope per
  100 MeV relative to the mean, and χ²/ndf);
- `R_gap`: the efficiency measured from 11 MeV to the training sidebands, divided by
  the line's prediction there. This is the primary test.
- `R_punzi_sidebands`: the same ratio in the Punzi sidebands (10–17 MeV), which hold a little signal.

R = 1 and a flat line mean the cut keeps the same fraction of background next to the
peak as far from it. In the feature study every candidate also gets `refit/sculpting.pdf`
(`ML_common/ml_plots.py`): the untrained DATA mass spectrum before and after the working-point
cut, and the background efficiency vs Bmass at the three cuts with the fitted line and the
mass regions shaded. `refit/mass_punzi_cut.pdf` shows the whole DATA spectrum (3.60–4.00 GeV, 10 MeV bins)
before and after the Punzi-optimal cut, with the X(3872) and ψ(2S) masses, the Punzi window and
sidebands, and the window count against the sideband estimate. `refit/punzi.pdf` shows the Punzi FOM and the
signal and background efficiencies versus the threshold (validation and training), with the
working point.

## Feature study (`feature_study.py`)

One study per method (xgb, nn, tmva) and sample (PbPb23, PbPb18), selecting from the 7
features of `FEATURE_POOL`. AUC is the central FOM; sigeff, pauc and the Punzi FOM are
recorded for every model and replayed.
- **Candidates.** All inputs are decorrelated (`_dc`). Each candidate is tuned by 100 condor
  jobs of 5 Optuna trials (TMVA: 100 random draws, every tree count evaluated); the best
  passing trial is refitted once and measured with `ml_fom.study_metrics`. Every refit also
  measures the importance of its inputs on the validation split: XGB mean |SHAP|, NN and
  TMVA the AUC drop when the input is shuffled.
- **fullPOOL**, the tentative selections:
  - `all`: all pool features at once. `production` (the production `FEATURES`, raw) is
    trained alongside as the reference.
  - `trimmed`: the inputs of `all` in decreasing importance up to the one that brings the
    cumulative share to 90%, without those below 5% on their own; done once.
    `comparison.json` and `comparison.pdf` compare it with `all`: the trimming, the FOMs with
    their paired differences, the working point, sculpting, spectator efficiencies, both
    DATA mass spectra and both Punzi curves.
- **Forward selection** (`FORWARD_SELECTION`). Step 0 is the two most important inputs of
  `all`; step k tries the step k−1 winner plus each remaining pool feature, up to the number
  of inputs of `trimmed`. The step winner has the highest validation AUC among the candidates
  whose refit passes the overtraining test and does not sculpt the mass (R_gap at the working
  point within 3σ of 1, `SCULPT_SIGMA`); its gain over the previous winner is a paired
  bootstrap on the shared validation candidates.
  - **Plateau:** the first step whose gain is below 2σ and stays below for 2 steps.
  - **Recommended set:** the smallest path set within 1σ (paired) of the best one.
  - **Replay:** for each FOM, where it would have chosen another candidate.
- **Final:** the last selection (the recommended set, or `trimmed` without the forward
  selection) without decorrelation.

The pool holds no feature that is a function of others: BtrkPtimb = (BtrkLeadPt −
BtrkSubPt)/(BtrkLeadPt + BtrkSubPt) was removed, because after decorrelation (a separate
mass-dependent map per feature) the relation between the three depends on Bmass, and an NN
learned the mass from it.

**A new setup starts every study again.** When the pool, the production inputs, the
decorrelation, the budget, the trimming rule, the forward flag, the pre-cut, the skims or the
training sidebands change, the next pass removes the study's jobs and results. A change of the
FOM settings or of `SCULPT_SIGMA` only measures the candidates again from their stored scores.

```bash
./feature_study.py run                                                  # in tmux: every study, a pass every 10 min
./feature_study.py advance --method xgb --sample PbPb18 --dry-run      # write the next submit files only
./feature_study.py status
./feature_study.py replay
./feature_study.py replot                                               # redraw every plot from the stored results
./feature_study.py adopt --method xgb --sample PbPb23 --selection trimmed   # all | trimmed | forward -> FEATURES
./feature_study.py candidate --method xgb --sample PbPb18 --features BtrkLeadPt,BtrkMaxdR   # one candidate, locally
./feature_study.py candidate --method xgb --sample PbPb18 --features BtrkLeadPt,BtrkMaxdR --condor   # the same on condor
```

The outputs live in `feature_study/<method>_<sample>/`:
- `fullPOOL/{all,trimmed,production}/` and `step_<k>/<candidate>/`, `final/no_decorrelation/`:
  `scan/`, `best.json`, `refit/`, `metrics.json`;
- `step_<k>/result.json`: the ranking, the winner and the paired gains;
- `summary.json`, `summary.pdf`: the selections and references side by side, the forward
  path per FOM and the sculpting of its winners.

Every candidate's `refit/` keeps the model, `scores.root` and the plots, about 10 MB (TMVA's
duplicate `tmva_training.root` is removed); scan jobs keep only their summaries. A scan job
that ends without output is resubmitted once; a candidate then goes on with the jobs it has,
or is marked `no_passing_trial`.

## Scored samples and metadata

Every scored file carries a `metadata/` directory of `TObjString`s: `classifier`,
`sample`, `fom`, `kind`, `tree`, `input_file`, `model`, `features`, `pre_cut` (the
pre-cut of the training that wrote the model), `mc_weight`, `score_branch`,
`score_min`, `score_max`, and `decorrelation_maps` for models with `_dc` inputs.
`optimalCUT_X_punzi.C` takes its pre-cut, trees, weight, and score range from there.
