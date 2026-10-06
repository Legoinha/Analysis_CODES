import argparse
import json
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
import uproot
import xgboost as xgb

import XGB_train as train_cfg
from XGB_train import ml_samples
import ml_decorrelation


SCRIPT_DIR = Path(__file__).resolve().parent
SCORE_BRANCH = "Prediction"
# selectionER/scored_samples/<tree>/<year>/flat_<tree>_<year>_xgb_<fom>_scored_<KIND>.root
SCORED_SAMPLE_DIR = Path(__file__).resolve().parents[3] / "scored_samples"
CLASSIFIER = "xgb"
KINDS = sorted({kind for inputs in ml_samples.settings.INPUTS.values() for _, kind, _ in inputs.values()})


def score_sample(model, part, output_path, step_size, metadata, maps):
    """One flat file of the case with the XGB score added; a _log input is computed from its flat branch
    (ml_samples.add_log_features), a _dc input from it and Bmass with the decorrelation maps of the
    training skims. The scored tree keeps only the flat branches."""
    feature_names = model.feature_names
    total_entries = 0
    print(f"Scoring {part['label']}: {part['file']}")
    # Written on local disk, then copied: uproot writing straight to the EOS fuse mount left corrupt baskets.
    local = Path(tempfile.mkdtemp()) / Path(output_path).name
    with uproot.open(ml_samples.XROOTD + part["file"]) as input_root:
        tree = input_root[part["tree"]]
        with uproot.recreate(local) as output_root:
            # Metadata first: written after a tree past 2 GB (PbPb23 DATA), uproot moves the
            # streamer record to the end of the file and overwrites its last bytes, which ROOT
            # then cannot read.
            for key, value in {**metadata, "kind": part["kind"], "tree": part["tree"], "input_file": part["file"]}.items():
                output_root[f"metadata/{key}"] = value

            for arrays in tree.iterate(step_size=step_size, library="np"):
                arrays = dict(arrays)
                inputs = dict(arrays)
                ml_samples.add_log_features(inputs, [name.removesuffix("_dc") for name in feature_names])
                features = np.column_stack([
                    maps.transform(name.removesuffix("_dc"), inputs[name.removesuffix("_dc")], arrays["Bmass"])
                    if name.endswith("_dc") else np.asarray(inputs[name], dtype=np.float32)
                    for name in feature_names
                ])
                arrays[SCORE_BRANCH] = model.predict(xgb.DMatrix(features, feature_names=feature_names)).astype(np.float32)
                if total_entries == 0:
                    output_root.mktree(part["tree"], {name: array.dtype for name, array in arrays.items()})
                output_root[part["tree"]].extend(arrays)
                total_entries += len(arrays[SCORE_BRANCH])
                print(f"  processed entries: {total_entries}")

    shutil.copy(local, f"{output_path}.new")
    os.replace(f"{output_path}.new", output_path)
    shutil.rmtree(local.parent)
    print(f"  saved: {output_path}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=ml_samples.CASES, required=True, help="<tree>/<year>, e.g. ntKp/PbPb23.")
    parser.add_argument("--fom", choices=train_cfg.ml_fom.FOMS, required=True)
    parser.add_argument("--kind", choices=["all"] + KINDS, default="all")
    parser.add_argument("--step-size", default="250 MB")
    parser.add_argument("--output-dir", help="Folder of the trained model (default: the run output folder).")
    parser.add_argument("--scored-dir", help="Folder of the scored files (default: scored_samples/<tree>/<year>).")
    return parser.parse_args()


def main():
    args = parse_args()
    train_cfg.configure_case(args.case, args.fom, output_dir=args.output_dir)
    c = train_cfg.CURRENT_CASE
    scored_dir = Path(args.scored_dir or SCORED_SAMPLE_DIR / c["tree"] / c["year"])
    scored_dir.mkdir(parents=True, exist_ok=True)

    # The saved model holds the trees up to the best iteration of its early stopping (XGB_train.py).
    model_path = train_cfg.OUTPUT_DIR / train_cfg.OUTPUT_MODEL
    model = xgb.Booster()
    model.load_model(model_path)

    print("Scoring case:", c["name"])
    print("Loaded model:", model_path)
    print("Using features:", ", ".join(model.feature_names))

    # The pre-ML cut is read from the report of the training that wrote this model.
    with uproot.open(train_cfg.OUTPUT_DIR / train_cfg.OUTPUT_REPORT) as report_file:
        report = json.loads(report_file["metadata_json"])
    metadata = {
        "classifier": CLASSIFIER,
        "sample": c["name"],
        "fom": train_cfg.CURRENT_FOM,
        "model": str(model_path) if model_path.is_absolute() else f"ML_selection/methods/{SCRIPT_DIR.name}/{model_path}",
        "features": ",".join(model.feature_names),
        "pre_cut": report["signal_cut"],
        "mc_weight": report["mc_weight_branch"],
        "score_branch": SCORE_BRANCH,
        "score_min": "0",
        "score_max": "1",
    }

    maps = None
    if any(name.endswith("_dc") for name in model.feature_names):
        metadata["decorrelation_maps"] = c["maps"]
        maps = ml_decorrelation.load(c["maps"])
    for part in c["inputs"].values():
        if args.kind in ("all", part["kind"]):
            score_sample(model, part,
                         scored_dir / f"flat_{c['tree']}_{c['year']}_{CLASSIFIER}_{args.fom}_scored_{part['kind']}.root",
                         args.step_size, metadata, maps)


if __name__ == "__main__":
    main()
