import argparse
import json
from pathlib import Path

import numpy as np
import torch
import uproot

import NN_train as train_cfg
from NN_train import ml_samples
import ml_decorrelation


SCRIPT_DIR = Path(__file__).resolve().parent
SCORE_BRANCH = "Prediction"
SCORED_SAMPLE_SYSTEMS = {
    "ppRef24": "ppRef",
    "PbPb18": "PbPb18",
    "PbPb23": "PbPb23",
    "PbPb24": "PbPb24",
}
SCORED_SAMPLE_DIR = SCRIPT_DIR / "scored_samples"
CLASSIFIER = "nn"


def samples_to_score(output_kind):
    config = train_cfg.SAMPLE_CONFIGS[train_cfg.CURRENT_SAMPLE]
    samples = [
        {"label": key, "input_file": config[key], "tree_name": ml_samples.TREES[key], "output_kind": ml_samples.KINDS[key]}
        for key in ml_samples.TREES
    ]
    if output_kind == "all":
        return samples
    return [sample for sample in samples if sample["output_kind"] == output_kind]


def score_sample(model, feature_names, sample, step_size, metadata, maps):
    system_name = SCORED_SAMPLE_SYSTEMS[train_cfg.CURRENT_SAMPLE]
    output_path = SCORED_SAMPLE_DIR / (
        f"flat_ntmix_{system_name}_{train_cfg.CURRENT_FOM}_scored_{sample['output_kind']}.root"
    )
    total_entries = 0

    print(f"Scoring {sample['label']}: {sample['input_file']}")
    with uproot.open(ml_samples.resolve_input_path(sample["input_file"])) as input_root:
        tree = input_root[sample["tree_name"]]
        with uproot.recreate(output_path) as output_root:
            first_chunk = True
            for arrays in tree.iterate(step_size=step_size, library="np"):
                arrays = dict(arrays)
                features = ml_samples.feature_matrix(arrays, feature_names, maps)
                arrays[SCORE_BRANCH] = train_cfg.predict_scores(model, features)

                if first_chunk:
                    branch_types = {
                        name: array.dtype
                        for name, array in arrays.items()
                    }
                    output_root.mktree(sample["tree_name"], branch_types)
                    first_chunk = False
                output_root[sample["tree_name"]].extend(arrays)
                total_entries += len(arrays[SCORE_BRANCH])
                print(f"  processed entries: {total_entries}")

            sample_metadata = {
                **metadata,
                "kind": sample["output_kind"],
                "tree": sample["tree_name"],
                "input_file": sample["input_file"],
            }
            for key, value in sample_metadata.items():
                output_root[f"metadata/{key}"] = value

    print(f"  saved: {output_path}")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample",
        choices=train_cfg.SAMPLE_CONFIGS,
        default=train_cfg.DEFAULT_SAMPLE,
    )
    parser.add_argument("--fom", choices=train_cfg.ml_fom.FOMS, required=True)
    parser.add_argument(
        "--kind",
        choices=(
            "all",
            "DATA",
            "MC_X3872",
            "MC_X3872_NONPROMPT",
            "MC_PSI2S",
            "MC_PSI2S_NONPROMPT",
        ),
        default="all",
    )
    parser.add_argument("--step-size", default="250 MB")
    parser.add_argument("--output-dir", help="Folder of the trained model (default: the run output folder).")
    parser.add_argument("--scored-dir", default=str(SCORED_SAMPLE_DIR), help="Folder of the scored files.")
    return parser.parse_args()


def main():
    global SCORED_SAMPLE_DIR
    args = parse_args()
    train_cfg.configure_sample(args.sample, args.fom, output_dir=args.output_dir)
    SCORED_SAMPLE_DIR = Path(args.scored_dir)
    torch.set_num_threads(train_cfg.NTHREAD)
    SCORED_SAMPLE_DIR.mkdir(parents=True, exist_ok=True)

    model_path = train_cfg.OUTPUT_DIR / train_cfg.OUTPUT_MODEL
    model, checkpoint = train_cfg.load_model(model_path)
    feature_names = checkpoint["feature_names"]

    print("Scoring sample:", train_cfg.CURRENT_SAMPLE)
    print("Loaded model:", model_path)
    print("Model trained on sample:", checkpoint["sample"])
    print("Using model from epoch:", checkpoint["best_epoch"])
    print("Using features:", ", ".join(feature_names))

    # The pre-ML cut is read from the report of the training that wrote this model.
    with uproot.open(train_cfg.OUTPUT_DIR / train_cfg.OUTPUT_REPORT) as report_file:
        report = json.loads(report_file["metadata_json"])
    metadata = {
        "classifier": CLASSIFIER,
        "sample": train_cfg.CURRENT_SAMPLE,
        "fom": train_cfg.CURRENT_FOM,
        "model": f"{SCRIPT_DIR.name}/{model_path}",
        "features": ",".join(feature_names),
        "pre_cut": report["signal_cut"],
        "mc_weight": report["mc_weight_branch"],
        "score_branch": SCORE_BRANCH,
        "score_min": "0",
        "score_max": "1",
    }

    maps = None
    if any(ml_samples.is_decorrelated(name) for name in feature_names):
        metadata["decorrelation_maps"] = ml_samples.decorrelation_path(train_cfg.CURRENT_SAMPLE)
        maps = ml_decorrelation.load(metadata["decorrelation_maps"])
    for sample in samples_to_score(args.kind):
        score_sample(model, feature_names, sample, args.step_size, metadata, maps)


if __name__ == "__main__":
    main()
