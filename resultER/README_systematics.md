# Systematic uncertainties of the X(3872) / psi(2S) ratio

Run the macros from `resultER/`, after the nominal and varied `effER` corrected-yield products and the `fitER`
fit-systematic products have been regenerated.

```bash
root -l -b -q 'ntmix_UNCpropagator.C("ppRef","Bpt")'
root -l -b -q 'ntmix_UNCpropagator.C("ppRef","nChargedTracks")'
root -l -b -q 'Ratio_CorrectedYields.C("ppRef","Bpt")'
```

The nominal observable in every analysis bin is

\[
R_{\rm nom} = Y^{\rm corr}_{X(3872),\,\rm 2D\ raw\ sPlot}
              /Y^{\rm corr}_{\psi(2S),\,\rm 2D\ raw\ sPlot}.
\]

The current sources are handled as follows.

- **Fit model:** `fitER` supplies one leading relative uncertainty per particle and bin. The source names and signed
  shifts are not available in that artifact, so the X(3872) and psi(2S) values are treated as independent:
  \(u_{\rm fit}^{R}=\sqrt{(u_{\rm fit}^{X})^2+(u_{\rm fit}^{\psi})^2}\).
- **Data-MC discrepancy in Acc x Eff:** the same reweighting case is applied to both corrected yields. For every named
  case, `resultER` calculates \(R_k=Y^X_k/Y^\psi_k\), then uses
  \(\max_k |R_k/R_{\rm nom}-1|\) in that bin. Common shifts therefore cancel and opposite shifts add naturally.
- **Acc x Eff method:** map-dimension and reading-method cases are paired in the same way. The largest ratio shift in
  each bin defines the source.
- **Total:** the fit-ratio, data-MC-ratio, and Acc x Eff method-ratio sources are added in quadrature.

`effER/accXeff_COMPARISONS.C` stores the corrected yield of each case in its comparison ROOT file. It also stores the
case list, reference case, tree, system, analysis variable, and varied observable as metadata. The ratio-level output
stores the nominal and varied ratios, every per-case percentage shift, each source total, the overall total, and the
correlation treatment. This makes the selected envelope traceable without reopening the individual corrected-yield
files.

The final uncertainty file is

```text
output_ntmix/systematicFILES/ntmix_ratioUnc_X3872_over_PSI2S_<YEAR>_<VAR>.root
```

`Ratio_CorrectedYields.C` uses `hTotalUncPercent` from this file for the inclusive, prompt, and nonprompt systematic
boxes. Prompt and nonprompt fraction errors are currently propagated into the point error together with the inclusive
statistical error; no separate systematic source for those fractions is defined yet.

This implementation is per year. A future combined-year result needs source-by-source correlations across years
before combining the yearly ratios. The stored efficiency cases and their signed ratio shifts prepare that extension,
but the fit-model products will first need source-resolved variations instead of only a leading magnitude. No
cross-year correlation assumption is made here.
