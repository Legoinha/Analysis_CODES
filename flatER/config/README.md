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

MC flattening reads `config/mc_normalization.csv` and adds two branches to both the
flattened reconstructed tree and `ntGen`: `pthat`, the generator pThat copied from the
forest, and `pThatreweight`, the event weight.

The MC productions are inclusive: the sample `pThat-X` contains every event with
pThat > X. The samples therefore overlap, and an event with generator pThat `p` can
come from every sample with threshold X_j < p. The merged sample has luminosity
`L_1 + ... + L_k` there, so each event is weighted by

```text
pThatreweight(p) = 1 / sum_{X_j < p} L_j,     L_j = n_gen_j / (xsec_pb_j * filter_eff_j)
```

The sum runs over all pThat rows of the same system, tree, particle and promptness.
Events with p between the two lowest thresholds keep the plain weight of the lowest
sample; above each further threshold the weight drops because more samples cover it.
This is the combined-sample form of the legacy recipe in
`Bfinder/Bfinder/weighPthat/weighPurePthat.C`, (sigma_k - sigma_k+1) / N_k, and agrees
with it in expectation. Adding the samples with the per-campaign weight
`xsec * filter_eff / n_gen` instead would count pThat > 10 twice, pThat > 15 three
times, and so on.

This needs forests made with the Bfinder version that stores `pthat` in `ntmix` and
`ntGen` (from the `GenEventInfoProduct` binning value). Older forests have no
`pthat` branch.

The value is stored, not automatically multiplied into the other branches.
Downstream histograms or fits should use `pThatreweight` explicitly.

The normalization table columns are:

```text
system,tree,particle,promptness,pthat,path_pattern,xsec_pb,filter_eff
```

`pthat` is the generation threshold X of the campaign. `n_gen` is not in the table:
`Flat_TREEs.C` assigns every input file to its campaign by `path_pattern` (exactly one
row of the group must match its path) and counts `n_gen` as the number of
Bfinder/ntGen events in the files of that campaign. The log prints the files and
`n_gen` found per campaign. A campaign without input files has `n_gen = 0` and adds no
luminosity. Because the weight sums over all campaigns of a group, `run_flat.sh` joins
all matched MC file lists into one flattening run. To enable PbPb24 or another system,
add its rows to the same CSV. Data flattening neither reads this table nor creates the
two branches.

The flattener stops with an error when an MC input file has no `pthat` branch or a
branch that is not `Float_t`. Without this check ROOT would only print an error and
leave `pthat = 0`, giving every event the weight of the lowest pThat campaign.

## Reweighting comparison plots

After producing the final merged MC file, `run_flat.sh` automatically compares the
reconstructed `Bpt` and the generator `pthat` distributions before and after applying
`pThatreweight`. The reweighted `pthat` spectrum must be smooth across the campaign
thresholds (5, 10, 15, 30, 50); a step there points to a wrong `xsec_pb`,
`filter_eff`, or missing input files of one campaign. Both
histograms are normalized to unit area and the y-axis is logarithmic. The style follows
`plotER/plot_dataMC.C`: a blue hatched
unweighted distribution and an orange reweighted line on a 600x600 canvas.

The comparison presentation is fixed; there is no plot-mode argument. The logarithmic
y-axis range of the Bpt plots is fixed to 10^-5 through 10^-1 for direct comparison
across samples.

The Bpt comparison uses 100 bins from 0 to 55 GeV/c; the pthat comparison uses 100
bins from 0 to 150 GeV/c with the y-axis up to 1. Each comparison is saved only
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
their variable and weight branch. The runner produces the `Bpt` and `pthat` plots.

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
3. it runs `Flat_TREEs.C` once per matched list, using `_0`, `_1`, `_2`, ... as `NUN`;
   for MC, all matched lists are joined into one run (`_0`)
4. for MC, it counts `n_gen` per pThat campaign and assigns `pThatreweight` per event
   from the generator `pthat`
5. it writes chunk outputs like `flat_ntmix_ppRef_DATA_0.root` in `flatER/` and keeps them at that creation path while merging
6. it always rebuilds a self-contained final file with `hadd`, even for one chunk
7. it removes the temporary chunk `.root` files
8. for MC, it saves the unweighted/reweighted `Bpt` and `pthat` comparison PDFs
