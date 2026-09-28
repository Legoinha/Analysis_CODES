"""Samples, pre-cuts, skims, and the loading of training data shared by ML_pytorch and ML_xgboost.

- FLAT_INPUTS: the flattened files of every sample (single source for the Python methods;
  ML_tmva/TMVA_config.h holds the same paths in C++).
- PRE_CUTS, TRAINING_SIDEBANDS, FEATURE_POOL: written by selectionER/fresh_ML_scan.sh.
- Skims: selectionER/make_skims.py writes, per sample, the flat trees after the pre-cut, plus the
  mass-decorrelated copy <feature>_dc of every feature it may be trained on. All trainings and
  scans read the skims; scoring (the *_apply scripts) reads the full flat files.
- load_training_data: signal MC and DATA sidebands from the skims, subsampled, balanced, and
  split into training / validation / test, with the DATA never used for training (for the
  sculpting check) and the spectator MC samples (for their efficiencies).
"""
import os
import re
import sys
import time

import numpy as np
import uproot

import ml_fom


FLAT_INPUTS = {
    'PbPb18': {
        'data': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_DATA.root',
        'signal': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_MC_X3872.root',
        'signal_nonprompt': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_MC_X3872_nonPrompt.root',
        'spectator': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_MC_PSI2S.root',
        'spectator_nonprompt': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_MC_PSI2S_nonPrompt.root',
    },
    'PbPb23': {
        'data': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root',
        'signal': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872.root',
        'signal_nonprompt': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872_nonPrompt.root',
        'spectator': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S.root',
        'spectator_nonprompt': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S_nonPrompt.root',
    },
    'PbPb24': {
        'data': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root',
        'signal': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_X3872.root',
        'signal_nonprompt': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_X3872_nonPrompt.root',
        'spectator': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S.root',
        'spectator_nonprompt': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S_nonPrompt.root',
    },
    'ppRef24': {
        'data': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root',
        'signal': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root',
        'signal_nonprompt': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872_nonPrompt.root',
        'spectator': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root',
        'spectator_nonprompt': '/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S_nonPrompt.root',
    },
}
TREES = {
    "data": "ntmix",
    "signal": "ntmix_X3872",
    "signal_nonprompt": "ntmix_X3872",
    "spectator": "ntmix_PSI2S",
    "spectator_nonprompt": "ntmix_PSI2S",
}
# Output kinds of the scored and skimmed files.
KINDS = {
    "data": "DATA",
    "signal": "MC_X3872",
    "signal_nonprompt": "MC_X3872_NONPROMPT",
    "spectator": "MC_PSI2S",
    "spectator_nonprompt": "MC_PSI2S_NONPROMPT",
}
SPECTATORS = ("signal_nonprompt", "spectator", "spectator_nonprompt")
MC_WEIGHT_BRANCH = "pThatreweight"

# Settings written by fresh_ML_scan.sh.
PRE_CUTS = {
    'PbPb18': '((Bpt > 15) && (Bpt < 50)) && (abs(By) < 2.4) && (BQvalue < 0.15) && (Bchi2Prob > 0.10)',
    'PbPb23': '((Bpt > 15) && (Bpt < 50)) && (abs(By) < 2.4) && (BQvalue < 0.15) && (Bchi2Prob > 0.05)',
}
TRAINING_SIDEBANDS = ((3.82164, 3.85664), (3.88664, 3.92164))
# Production classifier inputs (a name ending in _dc is the mass-decorrelated feature).
FEATURES = [
    'Btrk1dR',
    'Btrk2Pt',
    'Bmu1pt',
    'Bmu2pt',
    'BtrkPtimb',
    'Bmu1eta',
    'Bmu2eta',
    'Bujeta',
    'Btktketa',
]
# Feature study: the features it selects from.
FEATURE_POOL = [
    'BtrkLeadPt',
    'BtrkMaxdR',
    'BtrkSubPt',
    'BtrkMindR',
    'BtrkMaxAbsEta',
    'BmuLeadPt',
    'BmuMaxdR',
]

SKIM_DIR = "/eos/user/h/hmarques/Analysis_CODES/selectionER/skims"
DECORRELATED_SUFFIX = "_dc"


def feature_universe():
    """Every flat-tree feature a model may use: the skims carry each of them and its _dc copy."""
    return list(dict.fromkeys(FEATURE_POOL + [base_feature(f) for f in FEATURES]))


def skim_path(sample, key):
    return f"{SKIM_DIR}/{sample}/skim_{sample}_{KINDS[key]}.root"


def decorrelation_path(sample):
    return f"{SKIM_DIR}/{sample}/decorrelation_{sample}.root"


def base_feature(name):
    """Name of the flat-tree branch a feature is made from: BtrkMaxdR_dc -> BtrkMaxdR."""
    return name[:-len(DECORRELATED_SUFFIX)] if name.endswith(DECORRELATED_SUFFIX) else name


def is_decorrelated(name):
    return name.endswith(DECORRELATED_SUFFIX)


# =============================================================================
# Reading
# =============================================================================

def resolve_input_path(file_path):
    path = str(file_path)
    if os.environ.get("ML_USE_XROOTD") == "1" and path.startswith("/eos/"):
        return f"root://eosuser.cern.ch/{path}"
    return path


# XRootD reads can stall when many batch jobs read the same EOS file at once.
# uproot gives up after 30 s by default, so allow a longer wait and retry.
XROOTD_TIMEOUT_S = 600
READ_ATTEMPTS = 4


def read_branches(file_path, tree_name, branch_names):
    for attempt in range(1, READ_ATTEMPTS + 1):
        try:
            root_file = uproot.open(resolve_input_path(file_path), timeout=XROOTD_TIMEOUT_S)
            tree = root_file[tree_name]
            result = int(tree.num_entries), tree.arrays(list(dict.fromkeys(branch_names)), library="np")
        except (OSError, TimeoutError) as error:
            if attempt == READ_ATTEMPTS:
                raise
            wait = 60 * attempt
            print(f"Read of {file_path} failed ({error!r}); retry {attempt} of {READ_ATTEMPTS - 1} in {wait} s",
                  file=sys.stderr, flush=True)
            time.sleep(wait)
            continue
        # All data is in memory now. A failed close on an overloaded EOS server
        # must not throw it away.
        try:
            root_file.close()
        except (OSError, TimeoutError) as error:
            print(f"Ignoring failed close of {file_path}: {error!r}", file=sys.stderr, flush=True)
        return result


def read_skim(sample, key, branches):
    return read_branches(skim_path(sample, key), TREES[key], branches)[1]


def skim_metadata(sample, key):
    with uproot.open(resolve_input_path(skim_path(sample, key))) as skim:
        return {name.split(";")[0].split("/")[-1]: str(skim[name])
                for name in skim.keys() if name.startswith("metadata/")}


def root_cut_to_numpy_expr(cut):
    expr = cut.strip()
    if not expr or expr == "1":
        return "True"
    return expr.replace("&&", "&").replace("||", "|")


def branches_from_cut(cut):
    if not cut or cut.strip() == "1":
        return []
    names = re.findall(r"\b[A-Za-z_][A-Za-z0-9_]*\b", cut)
    reserved = {"True", "False", "and", "or", "not", "abs"}
    return [name for name in names if name not in reserved]


def cut_mask(arrays, cut):
    # Compare in double precision like ROOT's TTreeFormula: with float32 branches NumPy
    # would compare in float32 and treat entries at a float32-rounded cut value differently.
    namespace = {name: np.asarray(arrays[name], dtype=np.float64) for name in branches_from_cut(cut)}
    namespace["abs"] = np.abs
    mask = eval(root_cut_to_numpy_expr(cut), {"__builtins__": {}}, namespace)
    if np.isscalar(mask):
        return np.full(len(arrays["Bmass"]), bool(mask), dtype=bool)
    return np.asarray(mask, dtype=bool)


def sideband_cut():
    """TRAINING_SIDEBANDS as a ROOT cut string."""
    return " || ".join(f"((Bmass > {low:.5f}) && (Bmass < {high:.5f}))" for low, high in TRAINING_SIDEBANDS)


def in_training_sidebands(bmass):
    bmass = np.asarray(bmass, dtype=np.float64)
    return np.any([(bmass > low) & (bmass < high) for low, high in TRAINING_SIDEBANDS], axis=0)


# =============================================================================
# Training data
# =============================================================================

def subsample_indices(n, max_entries, random_state):
    if max_entries is None or n <= max_entries:
        return np.arange(n)
    return np.sort(np.random.default_rng(random_state).choice(n, size=max_entries, replace=False))


def make_balanced_weights(signal_weights, background_weights):
    signal_weights = np.asarray(signal_weights, dtype=np.float64)
    background_weights = np.asarray(background_weights, dtype=np.float64)
    target_class_sum = 0.5 * (len(signal_weights) + len(background_weights))
    balanced_signal = signal_weights * (target_class_sum / np.sum(signal_weights))
    balanced_background = background_weights * (target_class_sum / np.sum(background_weights))
    return balanced_signal.astype(np.float32), balanced_background.astype(np.float32)


def columns(arrays, features):
    return np.column_stack([np.asarray(arrays[name], dtype=np.float32) for name in features])


def feature_matrix(arrays, names, maps=None):
    """Classifier inputs from flat-tree arrays: a name ending in _dc is computed from its flat
    branch and Bmass with the decorrelation maps (ml_decorrelation.load(decorrelation_path(sample)))."""
    return np.column_stack([
        maps.transform(base_feature(name), arrays[base_feature(name)], arrays["Bmass"]) if is_decorrelated(name)
        else np.asarray(arrays[name], dtype=np.float32)
        for name in names
    ])


def weight_stats(weights):
    weights = np.asarray(weights, dtype=np.float64)
    return {
        "entries": int(len(weights)),
        "weight_sum": float(np.sum(weights)),
        "effective_entries": float(ml_fom.effective_entries(weights)),
        "weight_min": float(np.min(weights)),
        "weight_max": float(np.max(weights)),
    }


def load_training_data(sample, features, max_signal, max_background, random_state, with_extras=False):
    """Training data of one sample from its skims.

    Signal: prompt X(3872) MC, weighted with pThatreweight. Background: DATA in the training
    sidebands, unit weights. Both are subsampled to max_signal / max_background, balanced to
    equal class sums, and split into training / validation / test (ml_fom). Each split is a
    dict of x, y, w, bmass and idx (the skim entry, which with y identifies a candidate).

    with_extras also returns "sculpt" (the DATA never used for training: every candidate
    outside the training sidebands plus the sideband candidates outside the training split)
    and "spectators" (nonprompt X(3872), prompt and nonprompt Psi(2S) MC, x and w).
    """
    features = list(features)
    signal = read_skim(sample, "signal", features + ["Bmass", MC_WEIGHT_BRANCH])
    data = read_skim(sample, "data", features + ["Bmass"])
    sideband = in_training_sidebands(data["Bmass"])

    sig_idx = subsample_indices(len(signal["Bmass"]), max_signal, random_state)
    bkg_idx = subsample_indices(int(sideband.sum()), max_background, random_state)
    bkg_idx = np.flatnonzero(sideband)[bkg_idx]
    w_sig, w_bkg = make_balanced_weights(signal[MC_WEIGHT_BRANCH][sig_idx], np.ones(len(bkg_idx)))
    arrays = {
        "x": np.vstack([columns(signal, features)[sig_idx], columns(data, features)[bkg_idx]]),
        "y": np.concatenate([np.ones(len(sig_idx), np.float32), np.zeros(len(bkg_idx), np.float32)]),
        "w": np.concatenate([w_sig, w_bkg]),
        "bmass": np.concatenate([signal["Bmass"][sig_idx], data["Bmass"][bkg_idx]]).astype(np.float64),
        "idx": np.concatenate([sig_idx, bkg_idx]).astype(np.int64),
    }
    train, validation, test = ml_fom.split_train_validation_test(arrays, "y", random_state)

    result = {
        "train": train,
        "validation": validation,
        "test": test,
        "window_background": ml_fom.window_background(data["Bmass"]),
        "pre_cut": skim_metadata(sample, "data")["pre_cut"],
        "stats": {
            "signal_selected": weight_stats(signal[MC_WEIGHT_BRANCH]),
            "signal_used": weight_stats(signal[MC_WEIGHT_BRANCH][sig_idx]),
            "background_selected_sidebands": int(sideband.sum()),
            "background_used": int(len(bkg_idx)),
            "data_entries": int(len(data["Bmass"])),
        },
    }
    if with_extras:
        trained = np.zeros(len(data["Bmass"]), dtype=bool)
        trained[train["idx"][train["y"] == 0]] = True
        result["sculpt"] = {"x": columns(data, features)[~trained], "bmass": data["Bmass"][~trained].astype(np.float64)}
        result["spectators"] = {}
        for key in SPECTATORS:
            arrays = read_skim(sample, key, features + [MC_WEIGHT_BRANCH])
            result["spectators"][key] = {"x": columns(arrays, features), "w": np.asarray(arrays[MC_WEIGHT_BRANCH], np.float64)}
    return result
