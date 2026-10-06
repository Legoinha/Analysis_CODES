"""Cases, skims and the training data shared by ML_pytorch, ML_xgboost, the feature study and the skims.

A case is one candidate tree in one data set, "<tree>/<year>" (ML_selection/settings.py, CASES).
case(name) gathers everything about it: its settings, the test variant applied, and its paths.
- VARIANT (environment ML_VARIANT): a test variant of the cases of one tree, with its own skims and
  studies (skims/<tree>/<year>_<variant>/, feature_study/<tree>/<year>_<variant>/); everything else as
  standard (VARIANT_TREES):
    2018MCall  ntmix: PbPb23 and PbPb24 DATA with the PbPb18 MC (PbPb18 itself is the standard)
    dR035      ntmix: both pions within dR < 0.35 of the J/psi, in every pre-cut
    XJpsi      ntmix: pool + Bujeta, Bujpt, By, Bpt, with |By| < 2.4
    small      ntmix: pool BtrkLeadPt, BtrkSubPt, BtrkMaxdR, without the forward selection
- Skims: ML_selection/scripts/make_skims.py writes, per case, the flat trees after the pre-cut, plus the
  mass-decorrelated copy <feature>_dc of every feature a model may use. A feature <branch>_log is
  log10 of the flat branch (add_log_features), and its _dc copy is decorrelated in log10. All trainings
  and scans read the skims; scoring (the *_apply scripts) reads the full flat files.
- load_training_data: signal MC and DATA sidebands from the skims, subsampled, balanced, and split
  into training / validation / test, with the DATA never used for training (for the sculpting check)
  and the spectator MC samples (for their efficiencies).
- write_cpp_settings: the C++ copy of settings.py for TMVA and optimalCUT_X_punzi.C.
"""
import copy
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import uproot

# settings.py is in ML_selection/; batch jobs get it next to this file.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import settings
import ml_fom


ML_SELECTION = "/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_selection"
VARIANT = os.environ.get("ML_VARIANT", "")
# Batch jobs (ML_USE_XROOTD=1) read the /eos/ files through XRootD.
XROOTD = "root://eosuser.cern.ch/" if os.environ.get("ML_USE_XROOTD") == "1" else ""
DECORRELATED_SUFFIX = "_dc"
LOG_SUFFIX = "_log"
MC_WEIGHT_BRANCH = settings.MC_WEIGHT_BRANCH
# Branch names in a ROOT cut (abs is the only function the cuts use).
CUT_NAMES = re.compile(r"\b[A-Za-z_]\w*")
# The tree of every test variant; with a variant, CASES holds only the cases it applies to.
VARIANT_TREES = {"2018MCall": "ntmix", "dR035": "ntmix", "XJpsi": "ntmix", "small": "ntmix"}
VARIANT_CASES = {variant: [f"{tree}/{year}" for tree, years in settings.CASES.items() for year in years
                           if tree == VARIANT_TREES.get(variant, tree) and not (variant == "2018MCall" and year == "PbPb18")]
                 for variant in [""] + list(VARIANT_TREES)}
CASES = VARIANT_CASES[VARIANT]


def case(name, variant=VARIANT):
    """Settings and paths of one case, "<tree>/<year>", in a test variant (default ML_VARIANT):
    settings.MASS and settings.CASES entries, the per-tree features, pool, forward selection, final check,
    decorrelation, DATA thinning, XGBoost depth cap and Punzi threshold step, and
      inputs                  per input key (settings.INPUTS): its flat tree, kind, label, flat file, skim and
                              skim_cut (the pre-cut; DATA also inside the decorrelation range, the only DATA
                              the selection uses; make_skims.py then thins it, settings.DATA_THINNING)
      spectators              the input keys other than data and signal
      universe                every feature a model may use: the skims carry each and its _dc copy (and, for
                              a <branch>_log, the flat branch)
      training_sideband_cut   the training sidebands as a ROOT cut
      skim_dir, maps          the skims and the decorrelation maps"""
    tree, year = name.split("/")
    c = copy.deepcopy({**settings.MASS[tree], **settings.CASES[tree][year]})
    c.update({"name": name, "tree": tree, "year": year, "label": f"{tree} {year}",
              "features": list(settings.FEATURES[tree]), "pool": list(settings.POOL[tree]),
              "forward_selection": settings.FORWARD_SELECTION[tree], "final_check": settings.FINAL_CHECK[tree],
              "decorrelate": dict(settings.DECORRELATE[tree]), "data_thinning": settings.DATA_THINNING[tree],
              "xgb_max_depth": settings.XGB_MAX_DEPTH[tree], "threshold_step": settings.THRESHOLD_STEP[tree]})
    if variant == "2018MCall":
        c["files"].update({key: path for key, path in settings.CASES["ntmix"]["PbPb18"]["files"].items() if key != "data"})
    elif variant == "dR035":
        c["pre_cut"] += " && (Btrk1dR < 0.35) && (Btrk2dR < 0.35)"
    elif variant == "XJpsi":
        c["pool"] += ["Bujeta", "Bujpt", "By", "Bpt"]
        c["pre_cut"] = c["pre_cut"].replace("abs(By) < 1.6", "abs(By) < 2.4")
    elif variant == "small":
        c["pool"] = ["BtrkLeadPt", "BtrkSubPt", "BtrkMaxdR"]
        c["forward_selection"] = False
    c["skim_dir"] = f"{ML_SELECTION}/skims/{tree}/{year}" + (f"_{variant}" if variant else "")
    c["maps"] = f"{c['skim_dir']}/decorrelation_{tree}_{year}.root"
    low, high = c["decorrelation_range"]
    c["inputs"] = {key: {"tree": flat_tree, "kind": kind, "label": label, "file": c["files"][key],
                         "skim": f"{c['skim_dir']}/skim_{tree}_{year}_{kind}.root",
                         "skim_cut": f"({c['pre_cut']}) && (Bmass >= {low:.5f}) && (Bmass <= {high:.5f})" if key == "data" else c["pre_cut"]}
                   for key, (flat_tree, kind, label) in settings.INPUTS[tree].items()}
    c["spectators"] = [key for key in c["inputs"] if key not in ("data", "signal")]
    c["universe"] = list(dict.fromkeys(c["pool"] + [f.removesuffix(DECORRELATED_SUFFIX) for f in c["features"]]))
    c["training_sideband_cut"] = " || ".join(f"((Bmass > {low:.5f}) && (Bmass < {high:.5f}))"
                                             for low, high in c["training_sidebands"])
    return c


# =============================================================================
# Reading
# =============================================================================

# XRootD reads can stall when many batch jobs read the same EOS file at once.
# uproot gives up after 30 s by default, so allow a longer wait and retry.
XROOTD_TIMEOUT_S = 600
READ_ATTEMPTS = 4


def add_log_features(arrays, names):
    """arrays (flat branches, NumPy) gets <branch>_log = log10(<branch>) for every such name."""
    for name in names:
        if name.endswith(LOG_SUFFIX):
            arrays[name] = np.log10(arrays[name.removesuffix(LOG_SUFFIX)]).astype(np.float32)


def read_branches(file_path, tree_name, branch_names):
    """(entries, {branch: array}) of one tree of a /eos/ file."""
    for attempt in range(1, READ_ATTEMPTS + 1):
        try:
            root_file = uproot.open(XROOTD + str(file_path), timeout=XROOTD_TIMEOUT_S)
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


def cut_mask(arrays, cut):
    """The entries passing a ROOT cut (&&, ||, abs), evaluated with NumPy in double precision like ROOT's
    TTreeFormula: with float32 branches NumPy would compare in float32 and treat entries at a
    float32-rounded cut value differently."""
    namespace = {name: np.asarray(arrays[name], dtype=np.float64) for name in CUT_NAMES.findall(cut) if name != "abs"}
    namespace["abs"] = np.abs
    return np.asarray(eval(cut.replace("&&", "&").replace("||", "|"), {"__builtins__": {}}, namespace), dtype=bool)


# =============================================================================
# Training data
# =============================================================================

def weight_stats(weights):
    weights = np.asarray(weights, dtype=np.float64)
    return {
        "entries": int(len(weights)),
        "weight_sum": float(np.sum(weights)),
        "effective_entries": float(np.sum(weights) ** 2 / np.sum(weights ** 2)),
        "weight_min": float(np.min(weights)),
        "weight_max": float(np.max(weights)),
    }


def load_training_data(c, features, max_signal, max_background, random_state, with_extras=False):
    """Training data of one case (case()) from its skims.

    Signal: the signal MC, weighted with pThatreweight. Background: DATA in the training sidebands,
    unit weights. Both are subsampled to max_signal / max_background, balanced to equal class sums, and
    split into training / validation / test (ml_fom). Each split is a dict of x, y, w, bmass and idx (the
    skim entry, which with y identifies a candidate).

    with_extras also returns "sculpt" (the DATA never used for training: every candidate outside the
    training sidebands plus the sideband candidates outside the training split) and "spectators" (x and
    w of every spectator MC sample of the case).
    """
    features = list(features)
    signal = read_branches(c["inputs"]["signal"]["skim"], c["inputs"]["signal"]["tree"], features + ["Bmass", MC_WEIGHT_BRANCH])[1]
    data = read_branches(c["inputs"]["data"]["skim"], c["inputs"]["data"]["tree"], features + ["Bmass"])[1]
    bmass = np.asarray(data["Bmass"], dtype=np.float64)
    sideband = np.any([(bmass > low) & (bmass < high) for low, high in c["training_sidebands"]], axis=0)

    n_signal, n_sideband = len(signal["Bmass"]), int(sideband.sum())
    sig_idx = (np.arange(n_signal) if n_signal <= max_signal else
               np.sort(np.random.default_rng(random_state).choice(n_signal, size=max_signal, replace=False)))
    bkg_idx = np.flatnonzero(sideband)[np.arange(n_sideband) if n_sideband <= max_background else
                                       np.sort(np.random.default_rng(random_state).choice(n_sideband, size=max_background, replace=False))]
    # Both classes sum to half the candidates.
    target_class_sum = 0.5 * (len(sig_idx) + len(bkg_idx))
    w_sig = np.asarray(signal[MC_WEIGHT_BRANCH][sig_idx], dtype=np.float64)
    w_sig = (w_sig * (target_class_sum / np.sum(w_sig))).astype(np.float32)
    w_bkg = np.full(len(bkg_idx), target_class_sum / len(bkg_idx), dtype=np.float32)
    arrays = {
        "x": np.vstack([np.column_stack([np.asarray(signal[f], np.float32) for f in features])[sig_idx],
                        np.column_stack([np.asarray(data[f], np.float32) for f in features])[bkg_idx]]),
        "y": np.concatenate([np.ones(len(sig_idx), np.float32), np.zeros(len(bkg_idx), np.float32)]),
        "w": np.concatenate([w_sig, w_bkg]),
        "bmass": np.concatenate([signal["Bmass"][sig_idx], data["Bmass"][bkg_idx]]).astype(np.float64),
        "idx": np.concatenate([sig_idx, bkg_idx]).astype(np.int64),
    }
    train, validation, test = ml_fom.split_train_validation_test(arrays, "y", random_state)
    with uproot.open(XROOTD + c["inputs"]["data"]["skim"]) as skim:
        pre_cut = str(skim["metadata/pre_cut"])

    result = {
        "train": train,
        "validation": validation,
        "test": test,
        "window_background": ml_fom.window_background(bmass, c),
        "pre_cut": pre_cut,
        "stats": {
            "signal_selected": weight_stats(signal[MC_WEIGHT_BRANCH]),
            "signal_used": weight_stats(signal[MC_WEIGHT_BRANCH][sig_idx]),
            "background_selected_sidebands": n_sideband,
            "background_used": int(len(bkg_idx)),
            "data_entries": int(len(bmass)),
        },
    }
    if with_extras:
        trained = np.zeros(len(bmass), dtype=bool)
        trained[train["idx"][train["y"] == 0]] = True
        result["sculpt"] = {"x": np.column_stack([np.asarray(data[f], np.float32)[~trained] for f in features]),
                            "bmass": bmass[~trained]}
        result["spectators"] = {}
        for key in c["spectators"]:
            spectator = read_branches(c["inputs"][key]["skim"], c["inputs"][key]["tree"], features + [MC_WEIGHT_BRANCH])[1]
            result["spectators"][key] = {"x": np.column_stack([np.asarray(spectator[f], np.float32) for f in features]),
                                         "w": np.asarray(spectator[MC_WEIGHT_BRANCH], np.float64)}
    return result


# =============================================================================
# C++ copy of the settings
# =============================================================================

def write_cpp_settings():
    """methods/ML_tmva/ml_settings.h, the settings TMVA and optimalCUT_X_punzi.C read, written from
    settings.py: every case, and every case of every test variant under "<tree>/<year>:<variant>". The
    file is written only when its content changes."""
    cases = []
    for variant, names in VARIANT_CASES.items():
        for name in names:
            c = case(name, variant)
            inputs = ",\n".join(f'       {{"{key}", {{"{part["tree"]}", "{part["kind"]}", "{part["label"]}", "{part["file"]}"}}}}'
                                for key, part in c["inputs"].items())
            spectators = ", ".join(f'"{key}"' for key in c["spectators"])
            features = ", ".join(f'"{f}"' for f in c["features"])
            pt_bins = ", ".join(repr(float(b)) for b in c["pt_bins"])
            sides = ", ".join(str(s) for s in c["punzi_sides"])
            cases.append(
                f'   {{"{name}{":" + variant if variant else ""}",\n'
                f'    {{"{c["tree"]}", "{c["year"]}",\n'
                f'     {{\n{inputs}\n     }},\n'
                f'     {{{spectators}}},\n'
                f'     {{{features}}},\n'
                f'     "{c["training_sideband_cut"]}",\n'
                f'     {c["peak"]!r}, {c["punzi_half_width"]!r}, {c["punzi_sideband_start"]!r}, {c["punzi_sideband_width"]!r},\n'
                f'     {{{sides}}},\n'
                f'     {{{pt_bins}}},\n'
                f'     "{c["skim_dir"]}",\n'
                f'     {c["threshold_step"]!r}}}}}')
    text = f"""// Written from ML_selection/settings.py by ml_samples.write_cpp_settings(): edit settings.py, not this file.
// The settings of TMVA (TMVA_config.h) and of optimalCUT_X_punzi.C.
#ifndef ML_SETTINGS_H
#define ML_SETTINGS_H

#include <map>
#include <vector>

#include "TString.h"

namespace MLSettings {{

struct Input {{
   TString tree;    // in the flat file and the skim
   TString kind;    // tag of the skim and scored files
   TString label;
   TString file;    // flat file
}};

// One case, "<tree>/<year>" (in a test variant "<tree>/<year>:<variant>"): its inputs (data, signal
// and spectators), production features, training sidebands, Punzi window and sidebands around the peak
// (GeV) and the sides of these (-1 below the peak, +1 above), pT bins of optimalCUT_X_punzi.C, skim
// folder, and step of the Punzi threshold grid.
struct Case {{
   TString tree;
   TString year;
   std::map<TString, Input> inputs;
   std::vector<TString> spectators;
   std::vector<TString> features;
   TString trainingSidebandCut;
   double peak;
   double punziHalfWidth;
   double punziSidebandStart;
   double punziSidebandWidth;
   std::vector<int> punziSides;
   std::vector<double> ptBins;
   TString skimDir;
   double thresholdStep;
}};

static const std::map<TString, Case> CASES = {{
{(","+chr(10)).join(cases)}
}};

static const TString MC_WEIGHT_BRANCH = "{settings.MC_WEIGHT_BRANCH}";
static const Long64_t MAX_SIGNAL = {settings.MAX_SIGNAL};
static const Long64_t MAX_BACKGROUND = {settings.MAX_BACKGROUND};
static const double BKG_EFF_TARGET = {settings.BKG_EFF_TARGET!r};
static const double PAUC_MAX_BKG_EFF = {settings.PAUC_MAX_BKG_EFF!r};
static const double PUNZI_A = {settings.PUNZI_A!r};
static const double PUNZI_B = {settings.PUNZI_B!r};
static const double PUNZI_MIN_BACKGROUND = {settings.PUNZI_MIN_BACKGROUND!r};
static const double KS_PVALUE_MIN = {settings.KS_PVALUE_MIN!r};
static const double MAX_AUC_GAP = {settings.MAX_AUC_GAP!r};
static const double AUC_GAP_SIGMA = {settings.AUC_GAP_SIGMA!r};

}} // namespace MLSettings

#endif
"""
    path = Path(__file__).resolve().parents[1] / "ML_tmva" / "ml_settings.h"
    if not path.exists() or path.read_text() != text:
        path.write_text(text)
