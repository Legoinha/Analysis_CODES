# combinedANA: simultaneous fit of the PbPb years

`roofitCombined.C` fits the PbPb years of one tree together, bin by bin: one `RooSimultaneous` over a `year`
category. It works for the five fit trees, `ntmix_X3872`, `ntmix_PSI2S`, `ntKp`, `ntphi` and `ntKstar`, with
the model of `../roofitB.h` in every year:

| tree | signal | background |
|---|---|---|
| `ntmix_X3872`, `ntmix_PSI2S` | double Gaussian | 2nd-order Chebychev |
| `ntphi`, `ntKstar` | double Gaussian | 3rd-order Chebychev |
| `ntKp` | double Gaussian | exponential + erfc (partially reconstructed B) |

Only the yields belong to each year: signal, combinatorial background, partially reconstructed B (ntKp) and the
MC normalization. Every other parameter is one object for all years:

- the signal shape (sigmas and fractions), fitted once on the MC of all years and then locked;
- the peak position `mean` and the data/MC resolution scale factor `scale`, floating in the data fit;
- the combinatorial background shape (Chebychev coefficients or exponential slope);
- the shape of the partially reconstructed B of `ntKp` (`m_nonprompt_scale`, `m_nonprompt_shift`).

The single-year fit, `../roofitB.C`, stays the fit of ppRef (one year per tree) and the check of each PbPb year on
its own. Both write the same files in the same place (Products).

## Inputs and selection

`doCombinedRoofit.sh <TREE>` takes, for every year, the scored samples of selectionER:

```text
selectionER/scored_samples/<case>/<year>/flat_<case>_<year>_<METHOD>_<FOM>_scored_DATA.root
selectionER/scored_samples/<case>/<year>/flat_<case>_<year>_<METHOD>_<FOM>_scored_<MC kind>.root
```

(`<case>` is `ntmix` for both X(3872) and psi(2S), whose MC kinds are `MC_X3872` and `MC_PSI2S`; the B mesons
use their tree and `MC`.) The selection of a year is the `pre_cut` of its scored DATA metadata and
`<score_branch> > threshold`. The threshold is the Punzi optimum of the inclusive bin,
`selectionER/optimalCUT/<case>/<year>/<METHOD>_<FOM>/punzi_summary.root`, so `optimalCUT_X_punzi.C` must have
run for every year first. X(3872) and psi(2S) take the threshold of the ntmix case: the same cut.

The driver settings (top of `doCombinedRoofit.sh`): `YEARS` per tree (PbPb23 and PbPb24; PbPb18 can be added to
ntmix), `METHOD` and `FOM` of the scored samples, the inclusive and pT bins, and `SYST`.

```bash
cd fitER/combinedANA
bash doCombinedRoofit.sh ntKp
bash doCombinedRoofit.sh ntmix_X3872
bash doCombinedRoofit.sh ntmix_PSI2S
```

The macro itself:

```cpp
roofitCombined(TREE, FULL, YEARS, INPUTDATA, INPUTMC, THRESHOLDS, VAR, VAR_BINS, SYST)
```

`FULL=1` is the inclusive fit (validation), `FULL=0` the binned fit (Acc x Eff); `SYST=1` adds the fit-model
variations (Systematics) to the binned fit. With a single year in `YEARS` the combined fit is a single-year fit
with the products below, which tests the chain before the other years exist.

## The fit of one bin

1. The MC of all years, fitted simultaneously in the signal window (pThat-weighted, `scale` fixed to 1, one shape
   and a normalization per year), sets the signal shape.
2. With the shape locked, the data of all years are fitted simultaneously over the mass range of `roofitB.C`:
   per-year yields, and one `mean`, one `scale` and one set of background shapes for all years.

## Systematics (`SYST=1`, binned fit)

Every variation of `GetBackgroundSystematicModels` and `GetSignalSystematicModels` (`../aux/uti.h`) is a whole
combined fit again, with the variation in every year at once and, again, everything but the yields common:

- background variations (2nd, 3rd, 4th-order Chebychev) replace the nominal background of every year; the MC fit
  is unchanged;
- signal variations (1gauss, 3gauss, gauss_cb) refit the MC of every year with the varied shape, then the data;
  `fixed` keeps the (shared) mean at its MC value.

For every year and bin, the deviation of the year's signal yield from the nominal one, |Y_var - Y_nom| / Y_nom, is
recorded; the background and signal systematics are the largest deviation of each family, the total their sum
in quadrature, as in `roofitB.C`. Each year gets its tables, its graph and
`fit_totalSystematic_leadingUnc_<TREE>_<YEAR>_<VAR>.root`, the resultER input. Because a variation is applied to
all years together, the deviations of the years are fully correlated by construction: a result summed over years
adds them linearly, not in quadrature. Sharing the shapes across the years is not varied.

## Significance (X(3872) and B0: `TreeSetup::significance`)

- On the plot of each year: S/sqrt(S+B) in +-2 sigma_eff around the mean (`GetSignificance`), computed on the
  per-year export of the combined fit.
- On the add-up plot and in `hCombinedSignificance`: the significance of the signal of all years together, from
  the likelihood ratio of the simultaneous fit against no signal in any year (all yields fixed to 0, everything
  else refitted), q0 = 2 (NLL_0 - NLL). With k years whose yields are bounded at 0, q0 follows the chi-bar-squared
  mixture sum_j C(k,j)/2^k chi2_j, which gives the p-value and Z_LR. Past q0 ~ 1400 the p-value underflows a
  double; Z_LR is then sqrt(q0), the asymptotic value for one bounded yield.

## Products

Everything goes to the folders of `roofitB.C`. **Each year gets the files of a single-year fit, under the same
names, in `fitER/ROOTfiles/<YEAR>/`, replacing the ones there**: validation, Acc x Eff and resultER read one place,
whichever fit wrote it last. To be sure the files of a PbPb year are those of the combined fit, run the combined
fit again; its files carry `fitStrategy = combined` (a single-year fit writes no such key).

```text
fitER/ROOTfiles/<YEAR>/nominalFitModel_<TREE>_<YEAR>.root          FULL=1   validation
fitER/ROOTfiles/<YEAR>/fitResults_<TREE>_<VAR>_<YEAR>.root         FULL=0   Acc x Eff (data)
fitER/ROOTfiles/<YEAR>/mcFitResults_<TREE>_<VAR>_<YEAR>.root       FULL=0   Acc x Eff (MC closure)
fitER/ROOTfiles/<YEAR>/systematicFILES/fit_totalSystematic_leadingUnc_<TREE>_<YEAR>_<VAR>.root   SYST=1, resultER
fitER/ROOTfiles/<COMBO>/combinedFitModel_<TREE>_<COMBO>.root              FULL=1
fitER/ROOTfiles/<COMBO>/combinedFitResults_<TREE>_<VAR>_<COMBO>.root      FULL=0

fitER/results/<YEAR>/<TREE>/<VAR>/data_<YEAR>_<VAR>_<low>_<high>[_<variation>]_<TREE>.pdf   (and mc_)
fitER/results/<YEAR>/tables/systematics_tables_<TREE>_<YEAR>_<VAR>.tex/.pdf
fitER/results/<YEAR>/Graphs/general_systematics_plot_<TREE>_<YEAR>_<VAR>.pdf
fitER/results/<COMBO>/<TREE>/<VAR>/data_<COMBO>_<VAR>_<low>_<high>_<TREE>.pdf   add-up plots (below)
fitER/results/<COMBO>/Graphs/{raw_yield,chi2,resolution,scale_variation}_<TREE>_<COMBO>_<VAR>.pdf
```

`<COMBO>` is the years joined by `_` (`PbPb23_PbPb24`).

**Add-up plot.** No fit is made to the summed data: the fit is simultaneous, and each year's model is normalized to
that year's yields. The add-up plot shows the data of all years summed, the combined fit (the sum of the yearly
fits, red, drawn above everything), the signal of every year (dashed: PbPb23 orange, PbPb24 yellowish) and the
background of all years ("Bkg.": combinatorial and partially reconstructed B, dashed blue); the pulls and the chi2/ndf compare the summed data with the
combined fit, with the floating parameters of the combined fit.

The per-year files have the keys and object names of `roofitB.C`: the year's data as `data`, `model<N>_`,
`nsig<N>_`, `nbkg<N>_`, `nbkg_part_r<N>_`, `mean<N>_`, `sigma1<N>_`, `sig_doubleG<N>_`, `bkg<N>_`, `erfc<N>`,
`scale`, the snapshots `nominalPars_bin<N>` (the floating parameters of the combined fit, shared ones included)
and `mcFitPars_bin<N>`, `hPt`, `hPtMC`, and the metadata keys `inputData`, `inputMC`, `selectionCut`,
`variableName`, `analysisVarBins`, `massMin`, `massMax`, `nMassBins`. They add `fitStrategy` (`combined`),
`combinedYears`, `sharedParameters` (`all but the yields`), `combinedFile` and `scoreThreshold`. So the validation and the Acc x Eff
correction run on each year independently, exactly as on a single-year fit.

The combined file holds `ws_combined` (the `RooSimultaneous` of every bin, the combined data with its `year`
category, the snapshots `nominalPars_bin<N>`), `hPt_<year>` of every year, `hPtSum` (the summed yield, its error
from the covariance of the combined fit), `hCombinedSignificance` (X(3872)) and the inputs and selection of every
year (`inputData_<year>`, `inputMC_<year>`, `selectionCut_<year>`).

## Downstream

Validation (`plotER/Validation/run_DataSIGNAL_VS_MC.sh <TREE> <YEAR>`), Acc x Eff (`effER/`) and resultER read
`fitER/ROOTfiles/<YEAR>/` for every year, ppRef and PbPb alike.
