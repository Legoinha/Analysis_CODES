# Efficiency Chain Runbook

## 1. Build ACCxEFF Maps

Build the efficiency maps once per particle and map choice.

```bash
root -l -b -q 'accXeff_2D.C("ntmix_X3872","ppRef","raw")'
root -l -b -q 'accXeff_2D.C("ntmix_X3872","ppRef","usePw")'
root -l -b -q 'accXeff_2D.C("ntmix_X3872","ppRef","useXw")'

root -l -b -q 'accXeff_2D.C("ntmix_PSI2S","ppRef","raw")'
root -l -b -q 'accXeff_2D.C("ntmix_PSI2S","ppRef","usePw")'
root -l -b -q 'accXeff_2D.C("ntmix_PSI2S","ppRef","useXw")'
```

Map meaning:
- `usePw`: nominal map, using the Psi2S Prediction weight.
- `useXw`: X(3872) Prediction-weight map variation.
- `raw`: no-reweight map variation.


















## 2. Method Variations ### X(3872) ### Psi2S
For method variations, keep the map fixed to `usePw`. This varies only the correction method: `1D`, `2D`, and `splot`.

```bash
root -l -b -q 'Apply_EffxAcc.C("ntmix_X3872","ppRef","Bpt","all","usePw")'
root -l -b -q 'Apply_EffxAcc.C("ntmix_X3872","ppRef","nChargedTracks","all","usePw")'
root -l -b -q 'Apply_EffxAcc.C("ntmix_PSI2S","ppRef","Bpt","all","usePw")'
root -l -b -q 'Apply_EffxAcc.C("ntmix_PSI2S","ppRef","nChargedTracks","all","usePw")'

root -l -b -q 'Compare_methods.C("ntmix_X3872","ppRef","Bpt")'
root -l -b -q 'Compare_methods.C("ntmix_X3872","ppRef","nChargedTracks")'
root -l -b -q 'Compare_methods.C("ntmix_PSI2S","ppRef","Bpt")'
root -l -b -q 'Compare_methods.C("ntmix_PSI2S","ppRef","nChargedTracks")'
```

## 3. Map Variations ### X(3872) ### Psi2S
For map variations, keep the method fixed to `splot`. This varies only the map: `raw`, `usePw`, and `useXw`.

```bash
root -l -b -q 'Apply_EffxAcc.C("ntmix_X3872","ppRef","Bpt","splot","all")'
root -l -b -q 'Apply_EffxAcc.C("ntmix_X3872","ppRef","nChargedTracks","splot","all")'
root -l -b -q 'Apply_EffxAcc.C("ntmix_PSI2S","ppRef","Bpt","splot","all")'
root -l -b -q 'Apply_EffxAcc.C("ntmix_PSI2S","ppRef","nChargedTracks","splot","all")'

root -l -b -q 'Compare_MapWeights.C("ntmix_X3872","ppRef","Bpt")'
root -l -b -q 'Compare_MapWeights.C("ntmix_X3872","ppRef","nChargedTracks")'
root -l -b -q 'Compare_MapWeights.C("ntmix_PSI2S","ppRef","Bpt")'
root -l -b -q 'Compare_MapWeights.C("ntmix_PSI2S","ppRef","nChargedTracks")'
```






## 4. Method Closure Test ### X(3872) ### Psi2S
For the method closure test, use reconstructed signal MC with the nominal `usePw` map fixed. The `1D`, `2D`, and `splot` correction factors are compared directly, with `splot` taken as the nominal method. Each command produces only the final comparison plot tagged with `_CLOSURE_`.

```bash
root -l -b -q 'Closure_methods.C("ntmix_X3872","ppRef","Bpt")'
root -l -b -q 'Closure_methods.C("ntmix_PSI2S","ppRef","Bpt")'
```

The closure test reuses the existing ACCxEFF maps. Rebuild the maps only when the map construction itself changes.



## Important

Do not mix method and map variations for systematic studies.

- Method variation: use `Apply_EffxAcc.C(...,"all","usePw")`, then `Compare_methods.C(...)` once at the end.
- Map variation: use `Apply_EffxAcc.C(...,"splot","all")`, then `Compare_MapWeights.C(...)` once at the end.

Avoid using `Apply_EffxAcc.C(...,"all","all")` for the systematic studies, because that produces both method and map variations together and makes it easier to mix the two effects.