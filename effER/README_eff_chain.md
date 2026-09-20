# Acc x Eff production chain

Run these commands from **effER/**. Supported names are exact and
case-sensitive:

- trees: **ntmix_X3872**, **ntmix_PSI2S**
- systems: **ppRef**, **PbPb23**
- map dimensions: **0D**, **1D**, **2D**
- efficiency weights used here: **raw**, **Bpt**, and **Score** for PbPb23
- reading methods: **sPlot**, **mWindow**

There are no aliases, fallback samples, or name normalization.

Only binned fitER products enter this chain:

~~~text
../fitER/ROOTfiles/<SYSTEM>/fitResults_<TREE>_Bpt_<SYSTEM>.root
../fitER/ROOTfiles/<SYSTEM>/fitResults_<TREE>_<VAR>_<SYSTEM>.root
../fitER/ROOTfiles/<SYSTEM>/mcFitResults_<TREE>_<VAR>_<SYSTEM>.root
~~~

The first file supplies **inputMC**, the reconstructed selection, and the Bpt
analysis bins to map production. There is no fallback MC path: binned products
made with an older fitER workflow must be regenerated with the current
**roofitB.C** first. The second file supplies fitted data yields, stored data,
snapshots, and models to map reading. The third is the MC-only closure input.

Validation supplies these exact files:

~~~text
../plotER/Validation/WEIGHTS/ntmix_<SYSTEM>_X3872_weight.root
../plotER/Validation/WEIGHTS/ntmix_<SYSTEM>_PSI2S_weight.root
~~~

The map code reads **hWeight_<WEIGHT>** and its stored
**weightExpression_<WEIGHT>**. A validation weight modifies only the selected
reconstructed numerator. Generator and acceptance counts are never
variable-reweighted. Every generated, accepted, and selected MC count uses
**pThatreweight**.

## 1. Produce the maps

~~~cpp
accXeff_MAPS(TREE, SYSTEM, DIMENSION, WEIGHT)
~~~

- **2D**: fine pT-versus-|y| map from MC counts; ROOT and PDF output.
- **1D**: pT axis of those counts, integrated over the configured |y| range;
  ROOT output only.
- **0D**: one direct **Nsel/Nacc** value in every Bpt analysis bin; ROOT output
  only.

Produce the three raw maps and the available reweighted 2D maps:

~~~bash
# ppRef X(3872)
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","ppRef","0D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","ppRef","1D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","ppRef","2D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","ppRef","2D","Bpt")'

# ppRef psi(2S)
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","ppRef","0D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","ppRef","1D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","ppRef","2D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","ppRef","2D","Bpt")'

# PbPb23 X(3872)
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","PbPb23","0D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","PbPb23","1D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","PbPb23","2D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","PbPb23","2D","Bpt")'
root -l -b -q 'accXeff_MAPS.C("ntmix_X3872","PbPb23","2D","Score")'

# PbPb23 psi(2S)
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","PbPb23","0D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","PbPb23","1D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","PbPb23","2D","raw")'
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","PbPb23","2D","Bpt")'
root -l -b -q 'accXeff_MAPS.C("ntmix_PSI2S","PbPb23","2D","Score")'
~~~

~~~text
output/<SYSTEM>/ROOTs/<TREE>_<SYSTEM><DIMENSION>map_ACCxEFF_<WEIGHT>.root
~~~

Only 2D calls write PDFs, under **output/<SYSTEM>/2Dmaps/**.

## 2. Read the maps in each analysis bin

~~~cpp
accXeff_READ(TREE, SYSTEM, VAR, METHOD, DIMENSION, WEIGHT)
~~~

**sPlot** reloads **nominalPars_binN** from the requested binned fit, floats the
signal and background yields for the sPlot covariance, and computes:

~~~text
corrected yield = sum_i signal_sWeight_i / AccEff_i
~~~

**mWindow** averages **1/AccEff** over all selected data candidates in the
resonance window, then multiplies it by the binned fitted signal yield. The
windows are X(3872) +/-20 MeV and psi(2S) +/-15 MeV.

Produce sPlot on every available map, plus the raw 2D mass-window case:

~~~bash
# ppRef X(3872)
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","Bpt","sPlot","0D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","Bpt","sPlot","1D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","Bpt","sPlot","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","Bpt","mWindow","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","Bpt","sPlot","2D","Bpt")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","nChargedTracks","sPlot","1D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","nChargedTracks","sPlot","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","nChargedTracks","mWindow","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","ppRef","nChargedTracks","sPlot","2D","Bpt")'

# ppRef psi(2S)
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","Bpt","sPlot","0D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","Bpt","sPlot","1D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","Bpt","sPlot","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","Bpt","mWindow","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","Bpt","sPlot","2D","Bpt")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","nChargedTracks","sPlot","1D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","nChargedTracks","sPlot","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","nChargedTracks","mWindow","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","ppRef","nChargedTracks","sPlot","2D","Bpt")'

# PbPb23 X(3872)
root -l -b -q 'accXeff_READ.C("ntmix_X3872","PbPb23","Bpt","sPlot","0D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","PbPb23","Bpt","sPlot","1D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","PbPb23","Bpt","sPlot","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","PbPb23","Bpt","mWindow","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","PbPb23","Bpt","sPlot","2D","Bpt")'
root -l -b -q 'accXeff_READ.C("ntmix_X3872","PbPb23","Bpt","sPlot","2D","Score")'

# PbPb23 psi(2S)
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","PbPb23","Bpt","sPlot","0D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","PbPb23","Bpt","sPlot","1D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","PbPb23","Bpt","sPlot","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","PbPb23","Bpt","mWindow","2D","raw")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","PbPb23","Bpt","sPlot","2D","Bpt")'
root -l -b -q 'accXeff_READ.C("ntmix_PSI2S","PbPb23","Bpt","sPlot","2D","Score")'
~~~

The raw 2D sPlot case is nominal, but it uses the same code path as every other
combination.

~~~text
output/<SYSTEM>/ROOTs/<TREE>_<SYSTEM>_<VAR>_<DIMENSION>_<WEIGHT>_<METHOD>_CorrectedYields.root
~~~

## 3. Produce the two systematic comparisons

The comparison receives explicit **{DIMENSION, WEIGHT, METHOD}** triples:

~~~cpp
accXeff_COMPARISONS(TREE, SYSTEM, VAR, CASES, TAG, SAVE_ROOT)
~~~

The raw 2D sPlot case is the reference. A PDF is always written. With
**SAVE_ROOT=true**, the uncertainty ROOT file is written too. With
**SAVE_ROOT=false**, any existing systematic ROOT artifact remains untouched,
so an exploratory plot cannot replace a resultER input.

Only these two production comparisons are required:

~~~bash
# ppRef X(3872)
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_X3872","ppRef","Bpt",{{"0D","raw","sPlot"},{"1D","raw","sPlot"},{"2D","raw","sPlot"},{"2D","raw","mWindow"}},"METHODS_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_X3872","ppRef","Bpt",{{"2D","raw","sPlot"},{"2D","Bpt","sPlot"}},"DATA_MC_AGREEMENT_reweight_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_X3872","ppRef","nChargedTracks",{{"1D","raw","sPlot"},{"2D","raw","sPlot"},{"2D","raw","mWindow"}},"METHODS_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_X3872","ppRef","nChargedTracks",{{"2D","raw","sPlot"},{"2D","Bpt","sPlot"}},"DATA_MC_AGREEMENT_reweight_comparison",true)'

# ppRef psi(2S)
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_PSI2S","ppRef","Bpt",{{"0D","raw","sPlot"},{"1D","raw","sPlot"},{"2D","raw","sPlot"},{"2D","raw","mWindow"}},"METHODS_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_PSI2S","ppRef","Bpt",{{"2D","raw","sPlot"},{"2D","Bpt","sPlot"}},"DATA_MC_AGREEMENT_reweight_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_PSI2S","ppRef","nChargedTracks",{{"1D","raw","sPlot"},{"2D","raw","sPlot"},{"2D","raw","mWindow"}},"METHODS_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_PSI2S","ppRef","nChargedTracks",{{"2D","raw","sPlot"},{"2D","Bpt","sPlot"}},"DATA_MC_AGREEMENT_reweight_comparison",true)'

# PbPb23 X(3872)
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_X3872","PbPb23","Bpt",{{"0D","raw","sPlot"},{"1D","raw","sPlot"},{"2D","raw","sPlot"},{"2D","raw","mWindow"}},"METHODS_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_X3872","PbPb23","Bpt",{{"2D","raw","sPlot"},{"2D","Score","sPlot"},{"2D","Bpt","sPlot"}},"DATA_MC_AGREEMENT_reweight_comparison",true)'

# PbPb23 psi(2S)
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_PSI2S","PbPb23","Bpt",{{"0D","raw","sPlot"},{"1D","raw","sPlot"},{"2D","raw","sPlot"},{"2D","raw","mWindow"}},"METHODS_comparison",true)'
root -l -b -q 'accXeff_COMPARISONS.C("ntmix_PSI2S","PbPb23","Bpt",{{"2D","raw","sPlot"},{"2D","Score","sPlot"},{"2D","Bpt","sPlot"}},"DATA_MC_AGREEMENT_reweight_comparison",true)'
~~~

~~~text
output/<SYSTEM>/systematicFILES/METHODS_comparison_<TREE>_<SYSTEM>_<VAR>.pdf
output/<SYSTEM>/ROOTs/METHODS_comparison_<TREE>_<SYSTEM>_<VAR>.root
output/<SYSTEM>/systematicFILES/DATA_MC_AGREEMENT_reweight_comparison_<TREE>_<SYSTEM>_<VAR>.pdf
output/<SYSTEM>/ROOTs/DATA_MC_AGREEMENT_reweight_comparison_<TREE>_<SYSTEM>_<VAR>.root
~~~

resultER reads only these exact tagged ROOT filenames. Other tags do not enter
uncertainty propagation.

## 4. Run the MC-only closure

**Closure_methods.C** never opens data. It reads the raw maps and the fitER
**mcFitResults** companion. All candidates in each stored **ntmix_X3872** or
**ntmix_PSI2S** MC dataset are selected and generator matched already.

The required ppRef and PbPb23 **mcFitResults** companions already exist. Run
the desired closure calls individually:

~~~bash
c~~~

The single closure PDF contains:

- raw 0D using all selected signal MC
- raw 1D using all selected signal MC
- raw 2D using all selected signal MC
- raw 2D using the fitted-MC-mean mass window
- raw 2D using the saved fitER MC model and sPlot

All five cases use pThat weights. No validation-derived Bpt, Score, or other
variable weight is used in closure.

~~~text
output/<SYSTEM>/closure/closure_<TREE>_<SYSTEM>_<VAR>.pdf
output/<SYSTEM>/ROOTs/closure_<TREE>_<SYSTEM>_<VAR>.root
~~~

## 5. resultER contract

~~~text
../effER/output/<SYSTEM>/ROOTs/ntmix_X3872_<SYSTEM>_<VAR>_2D_raw_sPlot_CorrectedYields.root
../effER/output/<SYSTEM>/ROOTs/ntmix_PSI2S_<SYSTEM>_<VAR>_2D_raw_sPlot_CorrectedYields.root
~~~

**ntmix_UNCpropagator.C** consumes the two exact comparison ROOT files from
section 3 as the independent Acc x Eff method and data-MC discrepancy sources.
