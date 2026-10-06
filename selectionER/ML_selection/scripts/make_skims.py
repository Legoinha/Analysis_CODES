#!/usr/bin/env python3
"""Build the training skims and the mass-decorrelation maps of one or more cases.

For each case "<tree>/<year>" (default: every case of ML_selection/settings.py; --tree and --year
narrow them):
1. every flat tree of the case (DATA, signal MC and the spectator MC samples) after the case's
   pre-cut (DATA also inside the decorrelation range of the tree, the only DATA the selection uses:
   the training and Punzi sidebands, the sculpting fit and the mass plots lie inside it), applied by
   ROOT, so the cut has exactly its ROOT meaning, with the branches the trainings
   read: the features of the case's universe (pool and production features), Bmass, the branches of the
   pre-cut and, for MC, pThatreweight (ROOT then reads only these of the ~80 flat branches: the B0 DATA
   flat file holds 264M candidates in 35 GB);
2. the decorrelation maps, fitted on all this DATA (methods/ML_common/ml_decorrelation.py);
3. every skim with these branches plus <feature>_dc for every feature of the universe, and a metadata/
   directory (sample, kind, pre-cut, source); the DATA skim keeps every candidate of the training
   sidebands and 1 in DATA_THINNING of the others (settings.py: B mesons 10, X(3872) 1);
4. the decorrelation check, from the written skims (the thinned DATA counted DATA_THINNING times):
   correlation of every feature with Bmass before and after, and the closure (largest |CDF - level| of
   the transformed background) in the fit region and in the gaps next to the peak the fit never saw.

Output in ML_selection/skims/<tree>/<year>/ (<year>_<variant> with ML_VARIANT, see ml_samples.py):
skim_<tree>_<year>_<KIND>.root, decorrelation_<tree>_<year>.root,
decorrelation_check_<tree>_<year>.pdf (percentiles vs Bmass per feature) and .json (the numbers).
Every file is written on local disk and then copied to EOS: uproot writing straight to the EOS fuse
mount left corrupt baskets.

    python3 make_skims.py                                   # all cases
    python3 make_skims.py --tree ntKp                       # every year of ntKp
    python3 make_skims.py --tree ntmix --year PbPb18
    python3 make_skims.py --tree ntmix --year PbPb18 --checks-only    # step 4 on the existing skims
"""
import argparse
import json
import os
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot
import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "methods" / "ML_common"))
import ml_decorrelation
import ml_fom
import ml_plots
import ml_samples

ROOT.gErrorIgnoreLevel = ROOT.kWarning


def make_skims(c):
    """Steps 1 to 3."""
    universe = c["universe"]
    flat = list(dict.fromkeys(f.removesuffix(ml_samples.LOG_SUFFIX) for f in universe))
    Path(c["skim_dir"]).mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as work:
        counts = {}
        cut_branches = [name for name in ml_samples.CUT_NAMES.findall(c["pre_cut"]) if name != "abs"]
        for key, part in c["inputs"].items():
            infile = ROOT.TFile.Open(part["file"])
            tree = infile.Get(part["tree"])
            tree.SetBranchStatus("*", 0)
            for name in dict.fromkeys(flat + ["Bmass"] + cut_branches + ([] if key == "data" else [ml_samples.MC_WEIGHT_BRANCH])):
                tree.SetBranchStatus(name, 1)
            outfile = ROOT.TFile(os.path.join(work, f"{key}.root"), "RECREATE")
            skim = tree.CopyTree(part["skim_cut"])
            counts[key] = (int(tree.GetEntries()), int(skim.GetEntries()))
            skim.Write()
            outfile.Close()
            infile.Close()
            print(f"  {key:20s} {counts[key][1]:>9d} of {counts[key][0]:>10d} entries pass {part['skim_cut']}")

        data = uproot.open(os.path.join(work, "data.root"))[c["inputs"]["data"]["tree"]].arrays(flat + ["Bmass"], library="np")
        ml_samples.add_log_features(data, universe)
        fit = ml_decorrelation.fit_region(data["Bmass"], c)
        coefficients, parts = ml_decorrelation.build(data["Bmass"][fit], {f: data[f][fit] for f in universe}, c)
        del data
        local_maps = os.path.join(work, "decorrelation.root")
        ml_decorrelation.save(local_maps, coefficients, c, {
            "sample": c["name"], "pre_cut": c["pre_cut"], "source": c["inputs"]["data"]["file"], "slices": len(parts),
        })
        shutil.copy(local_maps, c["maps"] + ".new")
        os.replace(c["maps"] + ".new", c["maps"])
        maps = ml_decorrelation.load(local_maps)
        print(f"  decorrelation maps: {len(parts)} slices of {1000 * c['decorrelation_slice']:.0f} MeV, "
              f"degree {ml_decorrelation.DEGREE}")

        for key, part in c["inputs"].items():
            local_skim = os.path.join(work, f"skim_{key}.root")
            with uproot.open(os.path.join(work, f"{key}.root")) as f, uproot.recreate(local_skim) as out:
                # Metadata first: written after a tree past 2 GB, uproot moves the streamer record to the
                # end of the file and overwrites its last bytes, which ROOT then cannot read.
                for meta_key, value in {
                    "sample": c["name"], "kind": part["kind"], "pre_cut": c["pre_cut"], "skim_cut": part["skim_cut"],
                    "source_file": part["file"], "source_tree": part["tree"],
                    "source_entries": counts[key][0], "entries": counts[key][1],
                    "decorrelated_features": ",".join(universe),
                    "decorrelation_maps": c["maps"],
                    "data_thinning": c["data_thinning"] if key == "data" else 1,
                    "created_utc": datetime.now(timezone.utc).isoformat(),
                }.items():
                    out[f"metadata/{meta_key}"] = str(value)
                # In chunks: the B0 DATA holds ~55M candidates. DATA: every candidate of the training
                # sidebands, 1 in data_thinning of the others (ml_fom.data_weights).
                entry = 0
                for chunk, arrays in enumerate(f[part["tree"]].iterate(step_size="200 MB", library="np")):
                    if key == "data":
                        index = entry + np.arange(len(arrays["Bmass"]))
                        entry += len(index)
                        keep = (ml_fom.data_weights(arrays["Bmass"], c) == 1) | (index % c["data_thinning"] == 0)
                        arrays = {name: values[keep] for name, values in arrays.items()}
                    ml_samples.add_log_features(arrays, universe)
                    for name in universe:
                        arrays[name + ml_samples.DECORRELATED_SUFFIX] = maps.transform(name, arrays[name], arrays["Bmass"])
                    if chunk == 0:
                        out.mktree(part["tree"], {name: array.dtype for name, array in arrays.items()})
                    out[part["tree"]].extend(arrays)
            shutil.copy(local_skim, part["skim"] + ".new")
            os.replace(part["skim"] + ".new", part["skim"])
            print(f"  wrote {part['skim']} ({os.path.getsize(part['skim']) / 1e6:.0f} MB)")


def check_decorrelation(c):
    """Step 4, from the skims and maps on disk: prints the numbers and writes the .json and .pdf."""
    universe = c["universe"]
    dc = [f + ml_samples.DECORRELATED_SUFFIX for f in universe]
    maps = ml_decorrelation.load(c["maps"])
    data = ml_samples.read_branches(c["inputs"]["data"]["skim"], c["inputs"]["data"]["tree"], universe + dc + ["Bmass"])[1]
    # The thinned DATA counted as often as it was thinned.
    expand = np.repeat(np.arange(len(data["Bmass"])), ml_fom.data_weights(data["Bmass"], c))
    data = {name: values[expand] for name, values in data.items()}
    signal = ml_samples.read_branches(c["inputs"]["signal"]["skim"], c["inputs"]["signal"]["tree"],
                                      universe + dc + [ml_samples.MC_WEIGHT_BRANCH])[1]
    fit = ml_decorrelation.fit_region(data["Bmass"], c)
    m, gap = c["peak"], c["sculpting_gap"]
    regions = {"fit region": c["decorrelation_range"],
               "gap below": (c["sculpting_excluded"][0], m - gap), "gap above": (m + gap, c["sculpting_excluded"][1])}
    print("  decorrelation check: correlation with Bmass (fit region) raw -> _dc, "
          "and largest |CDF - level| of the transformed background per region")
    numbers = {}
    for f in universe:
        closure = {}
        for label, (low, high) in regions.items():
            mask = (data["Bmass"] >= low) & (data["Bmass"] < high) & (fit if label == "fit region" else True)
            closure.update(ml_decorrelation.closure(maps, f, data[f + ml_samples.DECORRELATED_SUFFIX][mask], data["Bmass"][mask],
                                                    {label: (low, high)}))
        numbers[f] = {"corr_raw": float(np.corrcoef(data[f][fit], data["Bmass"][fit])[0, 1]),
                      "corr_dc": float(np.corrcoef(data[f + ml_samples.DECORRELATED_SUFFIX][fit], data["Bmass"][fit])[0, 1]),
                      "closure": closure}
        print(f"    {f:22s} {numbers[f]['corr_raw']:+.3f} -> {numbers[f]['corr_dc']:+.3f}   "
              + "  ".join(f"{k} {v:.3f}" for k, v in closure.items()))
    stem = f"{c['skim_dir']}/decorrelation_check_{c['tree']}_{c['year']}"
    Path(stem + ".json").write_text(json.dumps(numbers, indent=2))
    ml_plots.decorrelation(data, signal, universe, numbers, stem + ".pdf",
                           f"{c['label']}: mass decorrelation of the classifier inputs", c)
    print(f"  wrote {stem}.pdf and .json")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tree", action="append", help="Only these trees (ntmix, ntKp, ntphi, ntKstar).")
    parser.add_argument("--year", action="append", help="Only these years (PbPb18, PbPb23, PbPb24, ppRef).")
    parser.add_argument("--checks-only", action="store_true", help="Only step 4, on the existing skims.")
    args = parser.parse_args()
    ml_samples.write_cpp_settings()
    for name in ml_samples.CASES:
        c = ml_samples.case(name)
        if (args.tree and c["tree"] not in args.tree) or (args.year and c["year"] not in args.year):
            continue
        print(f"=== {name}: pre-cut {c['pre_cut']}")
        if not args.checks_only:
            make_skims(c)
        check_decorrelation(c)


if __name__ == "__main__":
    main()
