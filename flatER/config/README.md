Store file lists under:

- `flatER/filelists/DATA/`
- `flatER/filelists/MC/`

Each `.txt` file must contain one ROOT file path or ROOT filename wildcard per line.

`run_flat.sh` interface:

```bash
bash flatER/run_flat.sh SYSTEM KIND TREE PARTICLE PVSNP
```

Arguments:

- `SYSTEM`: for example `ppRef`, `PbPb23`, `PbPb24`, `PbPb25`
- `KIND`  : `DATA` or `MC`
- `TREE`  : `ntmix`, `ntphi`, `ntKp`, `ntKstar`

- `PARTICLE`: only used for `ntmix` MC, for example `"_X3872"` or `"_PSI2S"`
- `PVSNP`   : only used for `ntmix` MC, `""` for prompt or `"_nonPrompt"` for nonprompt

Examples:

```bash
bash run_flat.sh ppRef DATA ntmix "" ""
bash run_flat.sh ppRef DATA ntKp "" ""

bash run_flat.sh ppRef MC ntKp "" ""
bash run_flat.sh ppRef MC ntmix _X3872 ""
bash run_flat.sh ppRef MC ntmix _PSI2S _nonPrompt
```

PbPb23 B-meson MC commands, run from `Analysis_CODES`:

```bash
bash flatER/run_flat.sh PbPb23 MC ntKp "" ""
bash flatER/run_flat.sh PbPb23 MC ntphi "" ""
bash flatER/run_flat.sh PbPb23 MC ntKstar "" ""
```

Matching rules:

- `DATA`: only the filename keywords `DATA` and `SYSTEM` are used
- `MC` with `TREE=ntmix` and `PVSNP=""`: this is the default prompt case, and the filename must match `MC`, `SYSTEM`, `TREE`, and `PARTICLE`, while not containing `nonprompt`
- `MC` with `TREE=ntmix` and `PVSNP="_nonPrompt"`: the filename must match `MC`, `SYSTEM`, `TREE`, and `PARTICLE`, and must contain `nonprompt`
- `MC` with `TREE!=ntmix`: the filename must match `MC`, `SYSTEM`, and `TREE`

Important details:

- flattened outputs are written under `/eos/user/h/hmarques/RUN3_Data_MC_sharing`
- `TREE=ntmix` writes to `X3872/<SYSTEM>`, with `ppRef` mapped to `X3872/ppRef24`
- `TREE!=ntmix` writes to `Bmesons/<SYSTEM>`
- `run_flat.sh` does not take a filelist path directly
- it looks inside `flatER/filelists/DATA/` or `flatER/filelists/MC/`
- your `.txt` files must be placed in the right subfolder, with names that match the requested case

## MC pThat normalization

MC flattening reads `config/mc_normalization.csv` and adds a `Double_t` branch named
`pThatreweight` to both the flattened reconstructed tree and `ntGen`. Every candidate
from a given source campaign receives the same campaign weight:

```text
pThatreweight = xsec_pb * filter_eff / n_gen
```

The value is stored, not automatically multiplied into the other branches. This is
intentional: a ROOT `TTree` cannot permanently attach a different implicit weight to
each entry, and changing the physical feature values would be incorrect. Downstream
histograms or fits should therefore use `pThatreweight` explicitly.

This implements the requested direct campaign normalization; it does not impose
exclusive generator-pThat intervals. If the productions are truly inclusive
`pThat > threshold` samples, the overlap/stitching convention should be confirmed
before a physics result is made from their sum.

The normalization table columns are:

```text
system,tree,particle,promptness,pthat,path_pattern,xsec_pb,filter_eff,n_gen
```

`path_pattern` identifies the production campaign and is independent of whether it
is an official or private sample. Before creating an output, the flattener validates
every MC file against exactly one row and checks that all `pThat`/`phat` tokens in its
path agree with the table. Missing, ambiguous, duplicated, or non-physical entries
stop the job instead of assigning a default weight.

The current registry uses ppRef and PbPb23 production-table cross sections/filter
efficiencies and actual full-campaign input event counts for `n_gen`. Nominal
requested counts must not be substituted for the count in the files being
flattened. Campaigns that are not present on EOS are omitted until they can be
counted. To enable PbPb24 or another system later, add its rows to the same CSV; no
change to the flattening logic is required. Data flattening neither reads this table
nor creates the `pThatreweight` branch.














## Reweighting comparison plots

After producing the final merged MC file, `run_flat.sh` automatically compares the
reconstructed `Bpt` distribution before and after applying `pThatreweight`. Both
histograms are normalized to unit area and the y-axis is logarithmic. The style follows
`plotER/plot_dataMC.C`: a blue hatched
unweighted distribution and an orange reweighted line on a 600x600 canvas.

The comparison presentation is fixed; there is no plot-mode argument. The logarithmic
y-axis range is fixed to 10^-5 through 10^-1 for direct comparison across samples.

The Bpt comparison uses 100 bins from 0 to 55 GeV/c. Each comparison is saved only
as a PDF under:

```text
flatER/reweighting_comparisons/
├── ntmix/    # prompt/nonprompt Psi(2S) and X(3872)
└── Bmeson/   # B+, Bs, and B0 channels
```

Filenames include the system, particle/promptness, variable, and weight branch, so
parallel jobs for different samples do not collide.

The plotting macro is standalone: if the flattened ROOT file already exists, no new
flattening is needed. Run these commands from `Analysis_CODES`.

### ppRef X(3872) and Psi(2S)

```bash
# Prompt Psi(2S)
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root","ntmix_PSI2S","ntmix","ppRef","_PSI2S","","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'

# Nonprompt Psi(2S)
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S_nonPrompt.root","ntmix_PSI2S","ntmix","ppRef","_PSI2S","_nonPrompt","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'

# Prompt X(3872)
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root","ntmix_X3872","ntmix","ppRef","_X3872","","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'

# Nonprompt X(3872)
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872_nonPrompt.root","ntmix_X3872","ntmix","ppRef","_X3872","_nonPrompt","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'
```

### PbPb23 X(3872) and Psi(2S)

```bash
# Prompt Psi(2S)
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S.root","ntmix_PSI2S","ntmix","PbPb23","_PSI2S","","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'

# Nonprompt Psi(2S)
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S_nonPrompt.root","ntmix_PSI2S","ntmix","PbPb23","_PSI2S","_nonPrompt","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'

# Prompt X(3872)
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872.root","ntmix_X3872","ntmix","PbPb23","_X3872","","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'

# Nonprompt X(3872)
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872_nonPrompt.root","ntmix_X3872","ntmix","PbPb23","_X3872","_nonPrompt","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'
```

### ppRef B mesons

```bash
# B+
root -l -b -q 'flatER/PlotReweightComparison.C("/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKp_ppRef_MC.root","ntKp","ntKp","ppRef","","","flatER/reweighting_comparisons","Bpt","pThatreweight",100)'
```

The current ppRef `ntKstar` (B0) and `ntphi` (Bs) flat files do not contain the
`pThatreweight` branch, so the comparison cannot run on them yet. Add their commands
after weighted flat outputs become available.

Future centrality or multiplicity comparisons can use the same function by passing
their variable and weight branch. Automatic runner production currently remains
restricted to `Bpt` with `pThatreweight`.

For example, with:

```bash
bash run_flat.sh PbPb24 DATA ntmix "" ""
```

the script will pick `.txt` files in `flatER/filelists/DATA/` whose names contain:

- `DATA`
- `PbPb24`

So names like these are good:

```text
DATA_PbPb24_00.txt
DATA_PbPb24_01.txt
...
```

Workflow:

1. `run_flat.sh` scans `filelists/DATA` or `filelists/MC`
2. it keeps only the `.txt` files matching the requested case
3. it runs `Flat_TREEs.C` once per matched list, using `_0`, `_1`, `_2`, ... as `NUN`
4. for MC, it validates each input campaign and assigns `pThatreweight`
5. it writes chunk outputs like `flat_ntmix_ppRef_DATA_0.root` in `flatER/` and keeps them at that creation path while merging
6. it always rebuilds a self-contained final file with `hadd`, even for one chunk
7. it removes the temporary chunk `.root` files
8. for MC, it saves the unweighted/reweighted `Bpt` comparison PDF and ROOT file
