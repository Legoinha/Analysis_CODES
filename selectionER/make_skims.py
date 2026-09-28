#!/usr/bin/env python3
"""Build the training skims and the mass-decorrelation maps of one or more samples.

For each sample (default: every sample with a pre-cut in ML_common/ml_samples.py):
1. every flat tree (DATA, prompt and nonprompt X(3872) and Psi(2S) MC) after the sample's
   pre-cut, applied by ROOT, so the cut has exactly its ROOT meaning;
2. the decorrelation maps, fitted on the DATA skim (ML_common/ml_decorrelation.py);
3. every skim with all flat-tree branches plus <feature>_dc for every feature of
   ml_samples.feature_universe(), and a metadata/ directory (sample, kind, pre-cut, source);
4. the decorrelation check, from the written skims: correlation of every feature with Bmass
   before and after, and the closure (largest |CDF - level| of the transformed background) in
   the fit region and in the mass gaps next to the X(3872) the fit never saw.

Output in selectionER/skims/<sample>/: skim_<sample>_<KIND>.root, decorrelation_<sample>.root,
decorrelation_check_<sample>.pdf (percentiles vs Bmass per feature) and .json (the numbers).

    python3 make_skims.py                 # all samples
    python3 make_skims.py --sample PbPb18
    python3 make_skims.py --sample PbPb18 --checks-only    # step 4 on the existing skims
"""
import json
import argparse
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot
import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parent / "ML_common"))
import ml_decorrelation
import ml_plots
import ml_samples

ROOT.gErrorIgnoreLevel = ROOT.kWarning


def skim_tree(sample, key, work):
    """Pre-cut applied by ROOT into a local temporary file; returns its path and entry counts."""
    source = ml_samples.FLAT_INPUTS[sample][key]
    tree_name = ml_samples.TREES[key]
    local = os.path.join(work, f"{key}.root")
    infile = ROOT.TFile.Open(source)
    tree = infile.Get(tree_name)
    outfile = ROOT.TFile(local, "RECREATE")
    skim = tree.CopyTree(ml_samples.PRE_CUTS[sample])
    counts = (int(tree.GetEntries()), int(skim.GetEntries()))
    skim.Write()
    outfile.Close()
    infile.Close()
    return local, counts


def make_skims(sample):
    """Steps 1 to 3."""
    universe = ml_samples.feature_universe()
    out_dir = Path(ml_samples.skim_path(sample, "data")).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as work:
        skims = {}
        for key in ml_samples.TREES:
            skims[key] = skim_tree(sample, key, work)
            print(f"  {key:20s} {skims[key][1][1]:>9d} of {skims[key][1][0]:>10d} entries pass the pre-cut")

        data = uproot.open(skims["data"][0])[ml_samples.TREES["data"]].arrays(universe + ["Bmass"], library="np")
        fit = ml_decorrelation.fit_region(data["Bmass"])
        coefficients, parts = ml_decorrelation.build(data["Bmass"][fit], {f: data[f][fit] for f in universe})
        ml_decorrelation.save(ml_samples.decorrelation_path(sample), coefficients, {
            "sample": sample, "pre_cut": ml_samples.PRE_CUTS[sample],
            "source": ml_samples.FLAT_INPUTS[sample]["data"], "slices": len(parts),
        })
        maps = ml_decorrelation.load(ml_samples.decorrelation_path(sample))
        print(f"  decorrelation maps: {len(parts)} slices of {1000 * ml_decorrelation.SLICE_WIDTH:.0f} MeV, "
              f"degree {ml_decorrelation.DEGREE}")

        for key, (local, (source_entries, entries)) in skims.items():
            with uproot.open(local) as f:
                arrays = f[ml_samples.TREES[key]].arrays(library="np")
            for name in universe:
                arrays[name + ml_samples.DECORRELATED_SUFFIX] = maps.transform(name, arrays[name], arrays["Bmass"])
            path = ml_samples.skim_path(sample, key)
            with uproot.recreate(path) as out:
                out[ml_samples.TREES[key]] = arrays
                for meta_key, value in {
                    "sample": sample, "kind": ml_samples.KINDS[key], "pre_cut": ml_samples.PRE_CUTS[sample],
                    "source_file": ml_samples.FLAT_INPUTS[sample][key], "source_tree": ml_samples.TREES[key],
                    "source_entries": source_entries, "entries": entries,
                    "decorrelated_features": ",".join(universe),
                    "decorrelation_maps": ml_samples.decorrelation_path(sample),
                    "created_utc": datetime.now(timezone.utc).isoformat(),
                }.items():
                    out[f"metadata/{meta_key}"] = str(value)
            print(f"  wrote {path} ({os.path.getsize(path) / 1e6:.0f} MB)")


def check_decorrelation(sample):
    """Step 4, from the skims and maps on disk: prints the numbers and writes the .json and .pdf."""
    universe = ml_samples.feature_universe()
    dc = [f + ml_samples.DECORRELATED_SUFFIX for f in universe]
    maps = ml_decorrelation.load(ml_samples.decorrelation_path(sample))
    data = ml_samples.read_skim(sample, "data", universe + dc + ["Bmass"])
    signal = ml_samples.read_skim(sample, "signal", universe + dc + [ml_samples.MC_WEIGHT_BRANCH])
    fit = ml_decorrelation.fit_region(data["Bmass"])
    m = ml_decorrelation.X3872_MASS
    regions = {"fit region": (ml_decorrelation.MASS_RANGE[0], ml_decorrelation.MASS_RANGE[1]),
               "X gap below": (3.85, m - 0.011), "X gap above": (m + 0.011, 3.89)}
    print(f"  decorrelation check: correlation with Bmass (fit region) raw -> _dc, "
          f"and largest |CDF - level| of the transformed background per region")
    numbers = {}
    for f in universe:
        closure = {}
        for label, (low, high) in regions.items():
            mask = (data["Bmass"] >= low) & (data["Bmass"] < high) & (fit if label == "fit region" else True)
            closure.update(ml_decorrelation.closure(maps, f, data[f][mask], data["Bmass"][mask], {label: (low, high)}))
        numbers[f] = {"corr_raw": float(np.corrcoef(data[f][fit], data["Bmass"][fit])[0, 1]),
                      "corr_dc": float(np.corrcoef(data[f + ml_samples.DECORRELATED_SUFFIX][fit], data["Bmass"][fit])[0, 1]),
                      "closure": closure}
        print(f"    {f:15s} {numbers[f]['corr_raw']:+.3f} -> {numbers[f]['corr_dc']:+.3f}   "
              + "  ".join(f"{k} {v:.3f}" for k, v in closure.items()))
    out = Path(ml_samples.decorrelation_path(sample)).parent
    (out / f"decorrelation_check_{sample}.json").write_text(json.dumps(numbers, indent=2))
    ml_plots.decorrelation(data, signal, universe, numbers, out / f"decorrelation_check_{sample}.pdf",
                           f"{sample}: mass decorrelation of the classifier inputs")
    print(f"  wrote {out}/decorrelation_check_{sample}.pdf and .json")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", action="append", choices=ml_samples.PRE_CUTS)
    parser.add_argument("--checks-only", action="store_true", help="Only step 4, on the existing skims.")
    args = parser.parse_args()
    for sample in args.sample or list(ml_samples.PRE_CUTS):
        print(f"=== {sample}: pre-cut {ml_samples.PRE_CUTS[sample]}")
        if not args.checks_only:
            make_skims(sample)
        check_decorrelation(sample)


if __name__ == "__main__":
    main()
