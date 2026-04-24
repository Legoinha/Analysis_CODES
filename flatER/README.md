Store file lists under:

- `flatER/filelists/DATA/`
- `flatER/filelists/MC/`

Each `.txt` file must contain one ROOT file path per line.

`run_flat.sh` interface:

```bash
bash flatER/run_flat.sh SYSTEM KIND TREE PARTICLE PVSNP
```

Arguments:

- `SYSTEM`: for example `ppRef`, `PbPb23`, `PbPb24`, `PbPb25`
- `KIND`  : `DATA` or `MC`
- `TREE`  : `ntmix`, `ntphi`, `ntKp`, `ntKstar`

- `PARTICLE`: only used for `ntmix` MC, for example `"_X3872"` or `"_Psi2S"`
- `PVSNP`   : only used for `ntmix` MC, `""` for prompt or `"_nonPrompt"` for nonprompt

Examples:

```bash
bash run_flat.sh ppRef DATA ntmix "" ""
bash run_flat.sh ppRef DATA ntKp "" ""

bash run_flat.sh ppRef MC ntKp "" ""
bash run_flat.sh ppRef MC ntmix _X3872 ""
bash run_flat.sh ppRef MC ntmix _PSI2S _nonPrompt
```

Matching rules:

- `DATA`: only the filename keywords `DATA` and `SYSTEM` are used
- `MC` with `TREE=ntmix` and `PVSNP=""`: this is the default prompt case, and the filename must match `MC`, `SYSTEM`, `TREE`, and `PARTICLE`, while not containing `nonprompt`
- `MC` with `TREE=ntmix` and `PVSNP="_nonPrompt"`: the filename must match `MC`, `SYSTEM`, `TREE`, and `PARTICLE`, and must contain `nonprompt`
- `MC` with `TREE!=ntmix`: the filename must match `MC`, `SYSTEM`, and `TREE`

Workflow:

1. `run_flat.sh` scans `filelists/DATA` or `filelists/MC`
2. it keeps only the `.txt` files matching the requested case
3. it runs `Flat_TREEs.C` once per matched list, using `_0`, `_1`, `_2`, ... as `NUN`
4. it writes chunk outputs like `flat_ntmix_ppRef_DATA_0.root`
5. it merges them with `hadd` into the final file
6. it removes the temporary chunk `.root` files