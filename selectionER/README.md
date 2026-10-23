# selectionER: ML selection and threshold optimization

Three classifiers (NN, XGBoost, TMVA BDTG; every tree uses the first two at present, `ML_METHODS`) separate signal MC from DATA in the mass sidebands, for
every **case**: one candidate tree in one data set, named `<tree>/<year>`.

| tree | signal | years | MC of PbPb24 | pre-cut | sidebands | feature-study pool | forward selection, final check |
|---|---|---|---|---|---|---|---|
| `ntmix` | prompt X(3872) → J/ψ ρ | PbPb18 | — | 15 < pT < 50, \|y\| < 1.6, BQ < 0.15, Bχ²Prob > 0.10 | both | BtrkLeadPt BtrkSubPt BtrkMaxdR BmuLeadPt Bpt BmuMaxdR | off |
| `ntmix` | prompt X(3872) → J/ψ ρ | Run 3, `<year>_<region>` with region `pt10y24`, `pt15y16`, `pt10y16` (pT > 10 or 15, \|y\| < 2.4 or 1.6): `PbPb25` (own MC), `PbPb23/24/25_mc2325` (PbPb23 official + PbPb25 private MC merged) | — | pT < 50, BQ < 0.15, Btrk2dR < 0.3 | both | PbPb25: Btrk1dR Btrk2Pt Btrk1Pt Bmu1pt; mc2325: BtrkLeadPt BtrkSubPt BtrkMaxdR BmuLeadPt (also the production features, raw for every method) | off |
| `ntKp` | B+ → J/ψ K+ | PbPb23, PbPb24, ppRef | PbPb23 | 7.5 < pT < 50, \|y\| < 2.4, Bχ²Prob > 0.05 | high only | Balpha_log Bchi2Prob Bnorm_svpvDistance_2D_log Bnorm_trk1Dxy Btrk1Pt_log Bujeta Btrk1dR | off |
| `ntphi` | Bs → J/ψ φ | PbPb23, PbPb24 | PbPb23 | 7.5 < pT < 50, \|y\| < 2.4, Bχ²Prob > 0.05 | both | as ntKp | off |
| `ntKstar` | B0 → J/ψ K*0 | PbPb23, PbPb24 | PbPb23 | 7.5 < pT < 50, \|y\| < 2.4, Bχ²Prob > 0.05 | both | as ntKp | off |

The B mesons use one set of inputs for the production features and the feature-study pool (fullPOOL/all,
then fullPOOL/trimmed with the most important ones), with XGBoost and the NN only (`ML_METHODS`; no TMVA).
B+ takes its training background and its Punzi background from the high sideband only: the low one holds
partially reconstructed B decays (a lost particle puts the mass below the peak). A feature
`<branch>_log` is log10 of the flat branch, computed where the `_dc` copies are (skims, scoring): Balpha,
the decay-length significance and the track pT span 1–4 decades (Balpha: signal 0.0005–2.4 rad), which
the NN needs as log10; XGBoost bins its inputs by quantiles, so log10 changes nothing for it. The B-meson
pre-cut tightens the Bχ²Prob > 0.005 of the flattening to 0.05. The X(3872) cases also report the efficiency of the spectator MC
samples (nonprompt X(3872), prompt and nonprompt ψ(2S)); the B mesons have none. Only the NN trains on
the mass-decorrelated features by default (`DECORRELATE`).

**Every setting is in [`ML_selection/settings.py`](ML_selection/settings.py)**: the cases and their flat
files, pre-cuts and pT bins, the production features and pool of each tree, the mass regions of each
tree (peak, Punzi window and sidebands, training sidebands, decorrelation and sculpting regions), the
sample sizes and the model-selection settings. A new case (ntmix/ppRef when its MC is ready, ppRef Bs
and B0, ...) is one more entry of `CASES`. An entry can also set any per-tree setting (`features`, `pool`,
`decorrelate`, ...) for itself (the Run 3 ntmix cases). TMVA and `optimalCUT_X_punzi.C` read the C++ copy
`ML_selection/methods/ML_tmva/ml_settings.h`, which `make_skims.py`, `feature_study.py`,
`fresh_ML_scan.sh` and `train_best_ML.sh` write again from `settings.py` when they start.

## Layout

```text
selectionER/
├── optimalCUT_X_punzi.C            Punzi threshold on the scored samples (standalone)
├── optimalCUT/<tree>/<year>/<method>_<fom>/        its results
├── scored_samples/<tree>/<year>/   the final samples: flat_<tree>_<year>_<method>_<fom>_scored_<KIND>.root
├── correlations/                   correlation_dataMC.py, and <tree>/<year>/ signal and background matrices
└── ML_selection/
    ├── settings.py                 all settings
    ├── scripts/                    make_skims.py, feature_study.py, fresh_ML_scan.sh, train_best_ML.sh
    ├── skims/<tree>/<year>/        skims, decorrelation maps and their check
    ├── feature_study/<tree>/<year>/<method>/       feature studies
    └── methods/                    the code the scripts drive
        ├── ML_common/              samples and training data, FOMs, decorrelation, plots
        ├── ML_pytorch/  ML_xgboost/  ML_tmva/      trainers, scans, scoring; outputs in
        │                           <method>_outputs/<tree>/<year>_<fom>/, condor_logs/
```

A test variant (`ML_VARIANT`, see `methods/ML_common/ml_samples.py`: `2018MCall`, `dR035`, `XJpsi`,
`small` for ntmix) has its skims and studies in `<tree>/<year>_<variant>/`; the C++
copy holds its cases as `<tree>/<year>:<variant>`.

## The chain

From `ML_selection/scripts/`, in the LCG view:

```bash
python3 make_skims.py --tree ntKp --year PbPb23                # once per pre-cut: skims and decorrelation maps
./feature_study.py run --tree ntKp                              # optional: choose the inputs (below)
./fresh_ML_scan.sh [auc|sigeff|pauc|punzi] [nn xgb tmva] [ntKp/PbPb23 ...]   # check inputs, clean, submit scans
./train_best_ML.sh ntKp/PbPb23 auc [nn xgb tmva]               # check scans, train the best trial
# score: the commands are printed at the end of train_best_ML.sh, into scored_samples/ntKp/PbPb23/
root -l -b -q 'optimalCUT_X_punzi.C("ntKp/PbPb23","xgb","auc")'   # from selectionER/: Punzi threshold
python3 ../../correlations/correlation_dataMC.py ntKp/PbPb23     # correlation matrices, from the skims
```

`fresh_ML_scan.sh` writes `ml_settings.h`, then checks every input tree of the selected cases
(branches exist, cuts parse, ROOT and the NumPy translation used by NN and XGB select the same
candidates) and that the skims match the current pre-cut and flat files, before it cleans and submits.

**PbPb24** has no MC of its own: its cases pair the 2024 DATA with the PbPb23 MC (and the PbPb23
pre-cut). Its background, decorrelation maps and Punzi sidebands come from the 2024 DATA.

## Skims and mass decorrelation (`scripts/make_skims.py`, `methods/ML_common/ml_decorrelation.py`)

Every training and scan reads the skims in `skims/<tree>/<year>/`: the flat trees after the case's
pre-cut (DATA also inside the decorrelation range of the tree, ±350 MeV around a B peak: the only DATA
the selection uses), with the branches the trainings read (every feature a model of the case may use,
from its production features and pool, Bmass, the pre-cut branches and, for MC, `pThatreweight`), plus
`<feature>_dc` for every feature. The DATA skim keeps every candidate of the training sidebands (which
hold the Punzi sidebands), the only DATA the training and the Punzi background use, and 1 in
`DATA_THINNING` of the others (B mesons 10, X(3872) 1): those serve the sculpting check, the mass plots
and the decorrelation check, which count each of them `DATA_THINNING` times (`ml_fom.data_weights`). The
decorrelation maps are fitted on all of the DATA before the thinning. The skims exist so that the hundreds
of scan jobs of a candidate read only the candidates and branches they need, with the `_dc` features
computed once: the B0 DATA flat file holds 264M candidates in 35 GB (both K/π mass assignments in a
±150 MeV K*0 window), of which ~90M pass its pre-cut and ~55M lie in the mass window, ~16M of them in the
sidebands, so its DATA skim keeps ~20M (ppRef B+: 4.2M in the window, 1.0M kept); about 2/3 of the
B-meson DATA has pT < 10 GeV. Scoring (`*_apply`) reads the full flat files and computes the `_dc` values
itself. The skims are written in chunks, on local disk first. The X(3872) skims made before 2026-10-04
still hold every flat branch.

`<feature>_dc` is the feature with its background dependence on the candidate mass removed by quantile
morphing. On DATA in the decorrelation range of the tree outside its excluded windows (X(3872):
3.60–4.00 GeV outside 3.85–3.89 and ±20 MeV around the ψ(2S); B mesons: ±350 MeV around the peak
outside ±65 MeV, the gap between the training sidebands), the background percentiles (0.5%, 1.5%, …, 99.5%) are measured in mass slices (10 MeV;
B mesons 20 MeV) and fitted versus mass (cubic in m − peak). A candidate's value is then mapped from the
percentiles at its mass to those at the peak. Background then has the same distribution at every mass,
and at the peak the map is the identity. `make_skims.py` checks this on the skims it writes:
- `decorrelation_check_<tree>_<year>.json` holds the feature–mass correlation before and after, and the
  closure (fraction below each percentile minus the level) in the fit region and in the gaps next to
  the peak the fit never saw;
- `decorrelation_check_<tree>_<year>.pdf` shows the correlations of all features, then per feature the
  background 10/50/90% percentiles vs Bmass before and after (flat after), and the signal MC
  distribution before and after (almost unchanged).

`python3 make_skims.py --tree ntmix --year PbPb18 --checks-only` redoes the check on the existing skims.
The maps (`decorrelation_<tree>_<year>.root`) are shared by Python and C++
(`methods/ML_tmva/TMVA_decorrelation.h`). A model with `_dc` inputs records them in the
`decorrelation_maps` metadata of its scored files.

## Runs are tagged with their FOM

The FOM a scan maximizes (`auc`, `sigeff`, `pauc`, `punzi`) is part of the run name, so runs with
different FOMs are kept side by side:

```text
ML_selection/methods/ML_<method>/<method>_outputs/ntKp/PbPb23_punzi/        scans, model, training report
ML_selection/methods/ML_<method>/condor_logs/optuna_ntKp_PbPb23_punzi_<cluster>_<job>.*   logs (optimize_* for TMVA)
scored_samples/ntKp/PbPb23/flat_ntKp_PbPb23_<method>_punzi_scored_<KIND>.root
optimalCUT/ntKp/PbPb23/<method>_punzi/                                        Punzi results
```

## Model selection (`methods/ML_common/ml_fom.py`, `methods/ML_tmva/TMVA_common.h`)

**Split.** 60% training, 20% validation, 20% test, stratified by class. The validation split chooses
the model: early stopping (NN, XGB), the number of trees (TMVA), and the hyperparameters of every
scan. The test split is used once, in the final training report. For TMVA, TMVA's test sample (40%)
is split into the validation and test parts at random (fixed seed).

**Figures of merit**, all on the weighted validation ROC curve:

- `auc`: the ROC area.
- `sigeff`: the signal efficiency at background efficiency `BKG_EFF_TARGET` (0.02, close to the
  background efficiency of the X(3872) Punzi working point).
- `pauc`: the ROC area for background efficiency from 0 to `PAUC_MAX_BKG_EFF` (0.05), divided by that
  width (0 to 1). It covers the region around the working point but averages over many ROC points,
  so it fluctuates less than `sigeff`.
- `punzi`: the best `e_S / S_min(B)` over thresholds, with the Punzi `a`, `b` of the settings, computed
  as `optimalCUT_X_punzi.C` computes it. `e_S` is the signal efficiency; `B` is the background left in
  the window after the cut: the DATA candidates in the Punzi sidebands that pass, scaled to the window
  width. In the scans and step comparisons it is evaluated on the validation split, where B is the
  whole-DATA sideband count before the cut times the passing fraction of the split's sideband
  candidates.

  **Window and sidebands** (`MASS` of each tree in `settings.py`). X(3872): |m − m_X| < 5 MeV,
  background from 15–35 MeV on both sides scaled by 10/40; the signal MC peak holds 68% within 6 MeV,
  so the window keeps about 60% of it, and only about 2% of the window signal falls in the sidebands.
  B mesons: |m − m_B| < 40 MeV (Bs 25 MeV, about 2σ of the MC peak), background from 100–200 MeV on
  both sides (B+: the high side only, scaled by 80/100). The training
  sidebands span 20–100 MeV from the peak for the X(3872) and 65–235 MeV for the B mesons; the Punzi background is counted on all DATA,
  the training candidates included, exactly as
  `optimalCUT_X_punzi.C` counts it on the full DATA. The maximum is taken only over cuts that keep at
  least `PUNZI_MIN_BACKGROUND` (10) background candidates. Otherwise a cut that leaves no validation
  background, and so gets B = 0, would win by chance: that happened for TMVA at 1.7% signal
  efficiency. `optimalCUT_X_punzi.C` applies the same window, sidebands and minimum (to the DATA
  sidebands).

  **Threshold grid.** The Punzi FOM is scanned on thresholds k × `THRESHOLD_STEP` (per tree: X(3872) 0.02, B mesons 0.001), the same grid
  in the model selection (Python and TMVA) and in `optimalCUT_X_punzi.C`. Scanning every candidate's
  score would let the maximum land on an upward fluctuation of the curve. `refit/punzi.pdf` shows the
  validation curve with its bootstrap ±1σ band: a training curve inside the band differs from it by
  statistics only.

**Overtraining is pass/fail, not a penalty.** A model passes when the train-validation AUC gap is at
most `MAX_AUC_GAP` or `AUC_GAP_SIGMA` (per tree: X(3872) 3, B mesons 2) times its statistical uncertainty, whichever is larger
(Hanley-McNeil errors of the two AUCs with the effective entries of each class, within 10% of the
bootstrap), and the weighted KS p-values of the signal and background scores (train vs validation,
effective sample sizes `(sum w)^2 / sum w^2`) are both at least `KS_PVALUE_MIN`. The uncertainty matters
for the X(3872) PbPb23 signal MC (X → J/ψ ρ): about 3.4k effective entries of 118k after the pre-cut
(the pThat-5 and pThat-10 campaigns carry about 70% of the weight); a model fitted on the test split
alone already showed a train-validation gap of 0.012-0.016 ± 0.011, so a fixed 0.01 rejected every NN
trial of a two-input model. The B-meson MC has more (7.5–60 GeV: B+ PbPb23 37k of 279k, ppRef 39k of 298k; Bs 19k of 156k; B0 7k of 109k;
nearly the same in 7.5–50 GeV). Failing
trials are rejected (Optuna: pruned; TMVA: `passed = 0`) and are never chosen. A scan job in which every
trial fails still writes its summary, with `best_value` null.

**Working point.** The working point is the Punzi maximum on all signal MC (training, validation and
test) and all DATA of the model, the threshold and value `optimalCUT_X_punzi.C` finds for it: the same
formula, window, sidebands, threshold grid and 10-candidate minimum. The training report compares the
signal and background efficiencies at that threshold between training and test, with binomial
uncertainties for weighted events and their pull: overtraining in the high-score tail, which biases
the efficiency of the analysis cut, shows up there even when the global AUC gap is small.

The `Validation:` and `Test:` lines of the training log and the `metrics` of `training_report.root`
(`validation_*`, `test_*`) hold all these numbers.

**Scores file.** Every training also writes `scores.root`, identical for the three methods
(`ml_fom.write_scores`):
- `train`, `validation`, `test`: y, w, score, bmass and idx (the skim entry);
- `sculpt`: the DATA never used for training;
- `spectator_<key>`: the spectator MC samples of the case (X(3872) only);
- `metadata_json`.

`ml_fom.study_metrics` computes everything the feature study compares from it:
- the four FOMs with bootstrap errors, the overtraining test, and the working-point efficiencies;
- the test split at the validation working point;
- the spectator efficiencies;
- the sculpting check.

**Sculpting check.** The background efficiency versus Bmass on the untrained DATA, at the working
point and at twice and half its background efficiency:
- a straight line fitted in bins of `sculpting_bin` over `sculpting_range` outside `sculpting_excluded`
  (X(3872): 5 MeV bins over 3.75–3.85 and 3.89–4.00 GeV; B mesons: 10 MeV bins from 150 to 350 MeV off
  the peak) (slope per 100 MeV relative to the mean, and χ²/ndf);
- `R_gap`: the efficiency measured from `sculpting_gap` off the peak (X(3872) 11 MeV, B+ and B0 60–70 MeV,
  Bs 40 MeV) to the excluded edges, divided by the line's prediction there. This is the primary test.
- `R_punzi_sidebands`: the same ratio in the Punzi sidebands.

R = 1 and a flat line mean the cut keeps the same fraction of background next to the peak as far from
it. For the B+ the partially reconstructed decays below the peak make the efficiency versus mass less
linear on the low side. In the feature study every candidate also gets `refit/sculpting.pdf`
(`methods/ML_common/ml_plots.py`): the untrained DATA mass spectrum before and after the working-point
cut, and the background efficiency vs Bmass at the three cuts with the fitted line and the mass regions
shaded. Every figure (study, trainers, TMVA, `optimalCUT_X_punzi.C`) carries its case in a box in a
corner. `refit/mass_punzi_cut.pdf` shows the whole DATA spectrum (the decorrelation range, 10 MeV bins)
before and after the Punzi-optimal cut, with the signal (and ψ(2S)) masses, the Punzi window and
sidebands, and the window count against the sideband estimate; pages 2–6 repeat it at cuts 5% tighter,
and 5%, 10%, 15% and 20% looser (threshold + f·|threshold|, f = +0.05, −0.05, −0.10, −0.15, −0.20), to
tell a real excess from a fluctuation of the few candidates left at high threshold. `refit/punzi.pdf`
shows the Punzi FOM and the signal and background efficiencies versus the threshold (validation and
training), with the working point.

## Feature study (`scripts/feature_study.py`)

One study per method (xgb, nn, tmva) and case, selecting from the pool of its tree (`POOL` in
`settings.py`), in `feature_study/<tree>/<year>/<method>/`. The central FOM (`CENTRAL_FOM` = auc, or
`--fom auc|sigeff|pauc|punzi`) tunes every candidate and chooses every forward step; all four FOMs are
recorded for every model and replayed. A study refuses another `--fom` than its own.
- **Test variants** (X(3872) only). `ML_VARIANT=<variant>` (`methods/ML_common/ml_samples.py`, passed to
  every job) runs a test round from the same code with its own skims (`skims/ntmix/<year>_<variant>/`,
  made by `ML_VARIANT=<variant> python3 make_skims.py`) and studies (`feature_study/ntmix/<year>_<variant>/`):
  `2018MCall` trains PbPb23 and PbPb24 on the PbPb18 MC, `dR035` adds Btrk1dR < 0.35 and Btrk2dR < 0.35
  to every pre-cut, `XJpsi` adds Bujeta, Bujpt, By and Bpt to the pool with |By| < 2.4, `small` uses the
  pool BtrkLeadPt, BtrkSubPt, BtrkMaxdR without the forward selection.
- **Candidates.** `DECORRELATE` sets the inputs per tree and method: the NN trains on the `_dc` features
  (with raw inputs it learned the mass: R_gap 14-26 and 0.1-0.5% ψ(2S) efficiency for the X(3872)
  `fullPOOL/all`), XGB and TMVA on the raw features, whose mass dependence the sculpting metrics watch.
  Each candidate is tuned by 100 condor jobs of 1 Optuna trial each, 100 independent random draws
  (TMVA: 100 random draws, every tree count evaluated); the best passing trial is refitted once and
  measured with `ml_fom.study_metrics`. Every refit also measures the importance of its inputs on the
  validation split: XGB mean |SHAP|, NN and TMVA the AUC drop when the input is shuffled.
- **fullPOOL**, the tentative selections:
  - `all`: all pool features at once. `production` (the production `FEATURES` of the tree, raw) is
    trained alongside as the reference.
  - `trimmed`: the inputs of `all` in decreasing importance up to the one that brings the cumulative
    share to 90%; done once. `comparison.json` and `comparison.pdf` compare it with `all`: the trimming,
    the FOMs with their paired differences, the working point, sculpting, spectator efficiencies, both
    DATA mass spectra and both Punzi curves.
- **Forward selection** (`FORWARD_SELECTION` of the tree: off for every tree at present). Step 0 is
  the two most important inputs of `all`; step k tries the step k−1 winner plus each remaining pool
  feature, for at most `MAX_ADDED` = 6 steps. The step winner has the highest validation central FOM
  among the candidates whose refit passes the overtraining test and does not sculpt the mass (R_gap at
  the working point within 3σ of 1, `SCULPT_SIGMA`); its gain over the previous winner is a paired
  bootstrap on the shared validation candidates.
  - **Plateau:** the first step whose gain is below 2σ and stays below for 2 steps.
  - **Recommended set:** the smallest path set within 1σ (paired) of the best one.
  - **Replay:** for each FOM, where it would have chosen another candidate.
- **Final** (`FINAL_CHECK` of the tree: off for every tree at present): the last selection (the
  recommended set, or `trimmed` without the forward selection or when no step-0 candidate passes) with
  the other decorrelation setting: XGB and TMVA with the `_dc` inputs (`final/decorrelation/`), the NN
  with the raw ones (`final/no_decorrelation/`).

The pool holds no feature that is a function of others: BtrkPtimb = (BtrkLeadPt − BtrkSubPt)/(BtrkLeadPt
+ BtrkSubPt) was removed, because after decorrelation (a separate mass-dependent map per feature) the
relation between the three depends on Bmass, and an NN learned the mass from it.

**A new setup starts every study again.** When the pool, the production inputs, the decorrelation, the
budget, the forward flag, the pre-cut, the skims, the training sidebands or the overtraining test
change, the next pass removes the study's jobs and results. A new trimming rule keeps `fullPOOL/all`
and `fullPOOL/production`; when it changes the trimmed inputs, `trimmed`, the forward path, `final` and
the summary start again. A change of the FOM settings (window, sidebands, ...) or of `SCULPT_SIGMA` only
measures the candidates again from their stored scores.

```bash
./feature_study.py run --tree ntKp --tree ntphi --tree ntKstar          # in tmux: those studies, a pass every 10 min
./feature_study.py advance --method xgb --tree ntKp --year ppRef --dry-run   # write the next submit files only
./feature_study.py status --year PbPb23
./feature_study.py replay --tree ntmix
./feature_study.py replot --tree ntmix                                       # redraw every plot from the stored results
./feature_study.py adopt --method xgb --tree ntKp --year PbPb23 --selection trimmed   # all | trimmed | forward -> FEATURES of the tree
./feature_study.py candidate --method xgb --tree ntmix --year PbPb18 --features BtrkLeadPt,BtrkMaxdR   # one candidate, locally
./feature_study.py candidate --method xgb --tree ntmix --year PbPb18 --features BtrkLeadPt,BtrkMaxdR --condor   # the same on condor
```

The outputs live in `feature_study/<tree>/<year>/<method>/`:
- `fullPOOL/{all,trimmed,production}/` and `step_<k>/<candidate>/`, `final/decorrelation/`: `scan/`,
  `best.json`, `refit/`, `metrics.json`;
- `step_<k>/result.json`: the ranking, the winner and the paired gains;
- `summary.json`, `summary.pdf`: the selections and references side by side, the forward path per FOM
  and the sculpting of its winners.

Every candidate's `refit/` keeps the model, `scores.root` and the plots, about 10 MB for X(3872) (TMVA's
duplicate `tmva_training.root` is removed); `scores.root` holds every DATA candidate of the skim, so a B0
refit takes about 0.3 GB. Once the candidate is measured, the results of all its scan jobs go into
`scan/trials.json`, only the best job's summary (which the refit reads) is kept, and the job logs and
condor log are removed (`./feature_study.py compact` does the same for older studies); submit files are
removed once submitted. The jobs cannot drop their stdout: the EOS-submit schedds refuse `output = /dev/null`.
A scan job that ends without output is resubmitted once; a candidate then goes on with the jobs it has,
or is marked `no_passing_trial`.

## Scored samples and metadata

The apply scripts write `scored_samples/<tree>/<year>/flat_<tree>_<year>_<method>_<fom>_scored_<KIND>.root`:
the flat tree with the `Prediction` branch, for every input of the case (X(3872): DATA, MC_X3872,
MC_X3872_NONPROMPT, MC_PSI2S, MC_PSI2S_NONPROMPT; B mesons: DATA, MC). Every scored file carries a
`metadata/` directory of `TObjString`s: `classifier`, `sample` (the case), `fom`, `kind`, `tree`,
`input_file`, `model`, `features`, `pre_cut` (the pre-cut of the training that wrote the model),
`mc_weight`, `score_branch`, `score_min`, `score_max`, and `decorrelation_maps` for models with `_dc`
inputs. `optimalCUT_X_punzi.C` takes its pre-cut, trees, weight, and score range from there, and the
Punzi window, sidebands and pT bins of the case from `ml_settings.h`.
