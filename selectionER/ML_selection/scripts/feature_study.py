#!/usr/bin/env python3
"""Feature selection for the classifiers of every case (ML_selection/settings.py).

A study is one method (xgb, nn, tmva) on one case "<tree>/<year>" (ntmix/PbPb23, ntKp/ppRef, ...). The
central figure of merit (CENTRAL_FOM = auc, or --fom: auc, sigeff, pauc or punzi) tunes every candidate
and chooses every forward step; all four FOMs are recorded for every model and replayed afterwards. The
studies live in ML_selection/feature_study/<tree>/<year>/<method>/ (--study-dir for another folder);
with ML_VARIANT=<variant> (methods/ML_common/ml_samples.py) a test variant's studies live in
<tree>/<year>_<variant>/ and read the skims of skims/<tree>/<year>_<variant>/.

A candidate is one set of inputs, tuned by N_JOBS condor jobs of N_TRIALS Optuna trials each (a
TMVA job is one random draw from the RANDOM_* ranges of TMVA_config.h, with every tree count
evaluated). Its best passing trial is refitted once, and ml_fom.study_metrics measures the refit
from its scores file, with the importance of every input (XGB: mean |SHAP|, NN and TMVA: the
validation AUC drop when the input is shuffled).

fullPOOL, the tentative selections
  all          every feature of the pool of the tree at once; the production inputs (settings.FEATURES
  production   of the tree, raw) are trained alongside as the reference, unless they are fullPOOL/all
               itself (the same inputs, raw for the method): then fullPOOL/all is the reference
  trimmed      the inputs of all in decreasing importance up to the one that reaches
               IMPORTANCE_CUMULATIVE of it; comparison.json and comparison.pdf compare it with all
Forward selection (settings.FORWARD_SELECTION of the tree; off for the B mesons)
  step 0       the ANCHOR_SIZE most important inputs of all
  step k       the step k-1 winner plus each remaining pool feature, one candidate each, up to
               MAX_ADDED steps. The winner has the highest validation central FOM among
               the candidates whose refit passes the overtraining test and does not sculpt the mass
               (R_gap within SCULPT_SIGMA of 1); its gain over the previous winner is a paired
               bootstrap on the shared validation candidates.
  plateau      the first step whose gain is below GAIN_SIGMA sigma and stays so for PATIENCE steps
  recommended  the smallest path set within WITHIN_SIGMA paired sigma of the best path set
  replay       for each FOM: the step winners it would have chosen, its plateau and recommended set
Final (settings.FINAL_CHECK of the tree; off for the B mesons): the last selection (recommended, or
trimmed without the forward selection or when no step-0 candidate passes) with the other decorrelation
setting. settings.DECORRELATE picks the baseline
per tree and method: the NN trains on the _dc features (with raw inputs it learns the mass from the
track pT and dR) and ends with the raw check (final/no_decorrelation); XGB and TMVA train on the raw
features and end with the decorrelated one (final/decorrelation).

When the setup changes (pool, production inputs, decorrelation, budget, forward selection, pre-cut,
skims, training sidebands, overtraining test), a study removes its jobs and results and starts again. A change of the
trimming rule alone keeps fullPOOL/all and fullPOOL/production: when the trimmed inputs change,
fullPOOL/trimmed, the forward path, the final check and the summary start again. A change of the FOM settings or of SCULPT_SIGMA only measures its candidates again from
their stored scores and closes the steps again.

Layout: feature_study/<tree>/<year>/<method>/
  study.json                  settings
  fullPOOL/{all,trimmed,production}/   candidates: candidate.json, scan/, best.json, refit/, metrics.json
                              refit/: model, training_report.root, scores.root, the trainer's plots
                              (ROC, scores, importance, history), sculpting.pdf, mass_punzi_cut.pdf, punzi.pdf
  fullPOOL/trimmed/comparison.json, comparison.pdf
  step_<k>/<name>/, step_<k>/result.json      the forward path: candidates, ranking, winner, gains
  final/decorrelation/        (final/no_decorrelation/ when DECORRELATE is True for the method)
  summary.json, summary.pdf   selections, references, forward path and replay

Usage, in the LCG view, from ML_selection/scripts/. --method, --tree and --year narrow a command to
some of the studies (every method of every case by default); --campaign FILE takes them from a file
("method tree/year" per line). --fom auc|sigeff|pauc|punzi picks the central FOM (default CENTRAL_FOM);
a study keeps the FOM it was started with.
  ./feature_study.py advance   [--method xgb] [--tree ntKp] [--year PbPb23] [--dry-run]   one pass: collect, submit
  ./feature_study.py run       [--tree ntKp] [--dry-run]   in tmux: a pass every POLL_MINUTES
  ./feature_study.py status    [--method xgb] [--year PbPb23]
  ./feature_study.py candidate --method xgb --tree ntmix --year PbPb18 --features BtrkLeadPt,BtrkMaxdR [--condor]
                               one candidate outside any study, in feature_study/<tree>/<year>/local/
                               (locally with --jobs 2 --trials 3 by default; --close-scan: stop waiting
                               for scan jobs)
  ./feature_study.py replay    [--method xgb] [--tree ntmix]            summary.json, printed
  ./feature_study.py replot    [--method xgb] [--tree ntmix]            redraw all plots from stored results
  ./feature_study.py adopt     --method xgb --tree ntKp --year PbPb23 --selection trimmed   FEATURES of the tree in settings.py
  ./feature_study.py compact   [--tree ntmix]   per-job scan results into scan/trials.json, job logs removed
                               (done automatically after every measurement)
--study-dir NAME works on ML_selection/NAME instead of feature_study (a round kept for reference); it
must hold the setup of the current skims, or the next advance starts its studies again, and its central
FOM must be the --fom one.
--dry-run writes the submit files of what would be submitted, and submits nothing.
"""
import argparse
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ML_selection/, with condor paths in the /eos/user/<initial>/<user> form, like the other submit files.
SELECTION_DIR = Path(re.sub(r"^/eos/home-(\w)/", r"/eos/user/\1/", str(Path(__file__).resolve().parents[1])))
sys.path.insert(0, str(SELECTION_DIR / "methods" / "ML_common"))
import ml_fom
import ml_plots
import ml_samples
from ml_samples import settings


# =============================================================================
# Settings
# =============================================================================

# The pool, the production inputs, the forward selection and the decorrelation of each tree are in
# settings.py (POOL, FEATURES, FORWARD_SELECTION, DECORRELATE).
CENTRAL_FOM = "auc"            # default of --fom: auc, sigeff, pauc or punzi (ml_fom.FOMS)
IMPORTANCE_CUMULATIVE = 0.90 # fullPOOL/trimmed: the inputs of fullPOOL/all in decreasing importance up to
                             # the one reaching 90% of it
ANCHOR_SIZE = 2              # forward selection: from the ANCHOR_SIZE most important inputs of
                             # fullPOOL/all, one input more per step
MAX_ADDED = 6                # hard cap on the steps (inputs added one at a time) after the anchor
N_JOBS = 100                 # scan jobs per candidate
N_TRIALS = 1                 # Optuna trials per job (xgb, nn): 100 jobs x 1 trial, 100 independent
                             # random draws per candidate, like TMVA
GAIN_SIGMA = 2.0             # a step is significant when its paired gain is >= GAIN_SIGMA sigma
PATIENCE = 2                 # non-significant steps in a row that make the plateau
WITHIN_SIGMA = 1.0           # recommended: smallest path set within this many paired sigma of the best
SCULPT_SIGMA = 3.0           # a candidate can win a step only if its R_gap (ml_fom.sculpting, working
                             # point) is within this many sigma of 1: a model that learns the mass
                             # changes the background next to the peak (R in the Punzi sidebands
                             # holds some real signal and is only reported)
MEASURE_WORKERS = 8          # candidates measured in parallel when the FOM settings change
MAX_RESUBMITS = 1            # a job that ends without output is resubmitted this many times
STALE_HOURS = 72             # a job still without output or end event after this is taken as lost
POLL_MINUTES = 10            # run: time between two passes
STUDY_DIR = SELECTION_DIR / "feature_study"
CONDOR_SETUP = "module load lxbatch/eossubmit >/dev/null 2>&1; "

COMMON_FILES = [SELECTION_DIR / "methods" / "ML_common" / name
                for name in ("ml_fom.py", "ml_samples.py", "ml_decorrelation.py", "ml_plots.py")] + [SELECTION_DIR / "settings.py"]
METHODS = {
    "xgb": {
        "dir": SELECTION_DIR / "methods" / "ML_xgboost", "executable": "run_optuna_batch.sh",
        "scan": "XGB_optuna.py", "train": "XGB_train.py", "summary": "summary_{seed}.json", "scan_args": [],
        "cpus": 4, "memory": "8 GB", "flavour": "workday",
    },
    "nn": {
        "dir": SELECTION_DIR / "methods" / "ML_pytorch", "executable": "run_nn_batch.sh",
        "scan": "NN_optuna.py", "train": "NN_train.py", "summary": "summary_{seed}.json",
        "scan_args": ["--timeout-hours", "6"],   # the job still writes its summary inside "workday"
        "cpus": 4, "memory": "8 GB", "flavour": "workday",
    },
    "tmva": {
        "dir": SELECTION_DIR / "methods" / "ML_tmva", "executable": "run_tmva.sh",
        "scan": "TMVA_optimize.C", "train": "TMVA_train.C", "summary": "summary_job_{seed}.root", "scan_args": [],
        "headers": ["TMVA_config.h", "TMVA_common.h", "TMVA_decorrelation.h", "ml_settings.h"],
        "cpus": 1, "memory": "10 GB", "flavour": "tomorrow",   # TMVA holds every sideband candidate: B0 ~7.5 GB
    },
}


# =============================================================================
# Small helpers
# =============================================================================

def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def to_json(value):
    if isinstance(value, dict):
        return {str(k): to_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_json(v) for v in value]
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.ndarray):
        return to_json(value.tolist())
    return value


def read_json(path, default=None):
    path = Path(path)
    return json.loads(path.read_text()) if path.exists() else default


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(to_json(value), indent=2, sort_keys=True))


def terminated_jobs(log_path):
    """(cluster, proc) of every job whose end (005 terminated, 009 aborted) is in a condor user log."""
    if not Path(log_path).exists():
        return set()
    pattern = re.compile(r"^(?:005|009) \((\d+)\.(\d+)\.\d+\)", re.M)
    return {(int(c), int(p)) for c, p in pattern.findall(Path(log_path).read_text(errors="replace"))}


def condor_submit(submit_file):
    out = subprocess.run(["bash", "-lc", CONDOR_SETUP + f"condor_submit {submit_file}"],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True).stdout
    print(out.strip())
    return int(re.search(r"submitted to cluster (\d+)", out).group(1))


def job_states(keys, done, submissions, log_path, closed=False):
    """Keys (scan seeds) by the state of their jobs: done (output written and no job
    of the key still running), lost (ended or stale without output, resubmits used up, or the scan
    closed), retry (ended without output, may be resubmitted), pending, unsubmitted.
    submissions: [{"cluster", "procs": {proc: key}, "submitted_utc"}]."""
    ended = terminated_jobs(log_path)
    tries = {}
    for sub in submissions:
        for proc, key in sub["procs"].items():
            age = (datetime.now(timezone.utc) - datetime.fromisoformat(sub["submitted_utc"])).total_seconds() / 3600.0
            finished = (sub["cluster"], int(proc)) in ended or age > STALE_HOURS
            tries.setdefault(key, []).append(finished)
    state = {"done": [], "lost": [], "retry": [], "pending": [], "unsubmitted": []}
    for key in keys:
        if key not in tries:
            state["done" if done(key) else "lost" if closed else "unsubmitted"].append(key)
        elif not all(tries[key]):
            state["pending"].append(key)
        elif done(key):
            state["done"].append(key)
        elif closed or len(tries[key]) > MAX_RESUBMITS:
            state["lost"].append(key)
        else:
            state["retry"].append(key)
    return state


def record_cluster(cluster, jobs, submissions_path):
    """Appends a submitted cluster to the submissions file of each of its candidates, with the
    candidate's procs: jobs are (candidate, key) pairs in proc order."""
    by_candidate = {}
    for proc, (cand, key) in enumerate(jobs):
        by_candidate.setdefault(cand.path, (cand, {}))[1][str(proc)] = key
    for cand, procs in by_candidate.values():
        submissions = read_json(submissions_path(cand), [])
        submissions.append({"cluster": cluster, "procs": procs, "submitted_utc": now()})
        write_json(submissions_path(cand), submissions)


# =============================================================================
# Candidates
# =============================================================================

def measure_again(path):
    """Metrics and plots of a measured candidate from its stored scores (a multiprocessing task)."""
    cand = Candidate(path)
    cand.measure(cand.metrics["best"], cand.metrics["scan"])


class Candidate:
    """One feature set: scan jobs, best trial, refit and metrics, all in one folder."""

    def __init__(self, path):
        self.path = Path(path)
        self.info = read_json(self.path / "candidate.json")
        self.method, self.sample = self.info["method"], self.info["sample"]
        self.case = ml_samples.case(self.sample)
        self.name, self.features = self.info["name"], self.info["features"]
        self.scan_dir = self.path / "scan"

    @staticmethod
    def create(path, method, sample, name, base, decorrelate, role, added=None,
               jobs=N_JOBS, trials=N_TRIALS, max_signal=None, max_background=None):
        path = Path(path)
        if not (path / "candidate.json").exists():
            write_json(path / "candidate.json", {
                "name": name, "method": method, "sample": sample, "role": role, "added": added,
                "base_features": list(base), "decorrelated": decorrelate,
                "features": [f + ml_samples.DECORRELATED_SUFFIX for f in base] if decorrelate else list(base),
                "fom": CENTRAL_FOM,
                "jobs": jobs, "trials": trials, "max_signal": max_signal, "max_background": max_background,
                "created_utc": now(),
            })
        for sub in ("scan/logs", "logs"):
            (path / sub).mkdir(parents=True, exist_ok=True)
        return Candidate(path)

    # ---- scan -----------------------------------------------------------------
    def scan_status(self):
        """Seeds by the state of their scan jobs (job_states; done: summary written)."""
        return job_states(range(self.info["jobs"]), lambda seed: (self.scan_dir / METHODS[self.method]["summary"].format(seed=seed)).exists(),
                          read_json(self.scan_dir / "submissions.json", []), self.scan_dir / "condor.log",
                          (self.scan_dir / "closed.json").exists())

    def close_scan(self):
        """Stops waiting for the scan: its jobs still in condor are removed, and the candidate goes
        on with the jobs that have written their summary."""
        pending = self.scan_status()["pending"]
        procs = [f"{sub['cluster']}.{proc}" for sub in read_json(self.scan_dir / "submissions.json", [])
                 for proc, seed in sub["procs"].items() if seed in pending]
        write_json(self.scan_dir / "closed.json", {"closed_utc": now(), "removed": procs})
        if procs:
            subprocess.run(["bash", "-lc", CONDOR_SETUP + "condor_rm " + " ".join(procs)], check=True)
        print(f"  scan closed: {len(procs)} running jobs removed")

    def best_trial(self):
        """Best passing trial over the finished scan jobs (None if none passed), and totals."""
        best, totals = None, {"jobs_with_output": 0, "trials": 0, "rejected_overtraining": 0}
        for seed in range(self.info["jobs"]):
            path = self.scan_dir / METHODS[self.method]["summary"].format(seed=seed)
            if not path.exists():
                continue
            totals["jobs_with_output"] += 1
            if self.method == "tmva":
                with uproot.open(path) as f:
                    rows = f["trial"].arrays(library="np")
                totals["trials"] += len(rows["objective"])
                totals["rejected_overtraining"] += int(np.sum(rows["passed"] == 0))
                for i in np.flatnonzero(rows["passed"] == 1):
                    if best is None or rows["objective"][i] > best["value"]:
                        best = {"job": seed, "summary": path.name, "value": float(rows["objective"][i]),
                                "params": {k: to_json(rows[k][i]) for k in (
                                    "n_trees", "max_depth", "shrinkage", "bagged_sample_fraction",
                                    "min_node_size_percent", "n_cuts")}}
            else:
                summary = read_json(path)
                totals["trials"] += summary["n_trials"]
                totals["rejected_overtraining"] += summary["n_rejected_overtraining"]
                if summary["best_value"] is not None and (best is None or summary["best_value"] > best["value"]):
                    best = {"job": seed, "summary": path.name, "value": summary["best_value"],
                            "params": summary["best_params"]}
        return best, totals

    # ---- refit ----------------------------------------------------------------
    def refit_status(self):
        """done (scores written), retry, lost, pending or unsubmitted."""
        if (self.path / "refit" / "scores.root").exists():
            return "done"
        submissions = read_json(self.path / "refit_submissions.json", [])
        if not submissions:
            return "unsubmitted"
        ended = terminated_jobs(self.path / "refit.log")
        if not all((s["cluster"], s["proc"]) in ended or
                   (datetime.now(timezone.utc) - datetime.fromisoformat(s["submitted_utc"])).total_seconds() / 3600.0 > STALE_HOURS
                   for s in submissions):
            return "pending"
        return "retry" if len(submissions) <= MAX_RESUBMITS else "lost"

    # ---- metrics --------------------------------------------------------------
    @property
    def metrics(self):
        return read_json(self.path / "metrics.json")

    def measure(self, best, totals):
        """metrics.json and the plots from the refit's scores file. TMVA's own output
        (tmva_training.root) only repeats what scores.root holds, so it is removed."""
        refit = self.path / "refit"
        scores = ml_fom.read_scores(refit / "scores.root")
        metrics = ml_fom.study_metrics(scores, self.case)
        status = "ok" if metrics["validation"]["passed"] else "refit_overtrained"
        write_json(self.path / "metrics.json", {**self.record(status, best, totals), "metrics": metrics})
        self.plot(scores, metrics)
        (refit / "tmva_training.root").unlink(missing_ok=True)

    def compact_scan(self):
        """After the measurement only the best scan job is needed (the refit reads its summary): the results
        of every job go into scan/trials.json, and the other summaries, the job logs and the condor log are
        removed. Done once per candidate."""
        if (self.scan_dir / "trials.json").exists():
            return
        best = read_json(self.path / "best.json")
        trials = {}
        for seed in range(self.info["jobs"]):
            path = self.scan_dir / METHODS[self.method]["summary"].format(seed=seed)
            if not path.exists():
                continue
            if self.method == "tmva":
                with uproot.open(path) as f:
                    trials[seed] = {k: v.tolist() for k, v in f["trial"].arrays(library="np").items()}
            else:
                summary = read_json(path)
                trials[seed] = {k: summary[k] for k in ("best_value", "best_params", "n_trials", "trials")}
            if best is None or path.name != best["summary"]:
                path.unlink()
        write_json(self.scan_dir / "trials.json", trials)
        shutil.rmtree(self.scan_dir / "logs", ignore_errors=True)
        (self.scan_dir / "condor.log").unlink(missing_ok=True)

    def plot(self, scores, metrics):
        """refit/sculpting.pdf, mass_punzi_cut.pdf and punzi.pdf (ML_common/ml_plots.py)."""
        refit = self.path / "refit"
        ml_plots.sculpting(scores["sculpt"]["bmass"], scores["sculpt"]["score"], metrics["sculpting"],
                           refit / "sculpting.pdf", "Untrained DATA", self.case)
        v, wp = metrics["validation"], metrics["working_point"]
        ml_plots.punzi_mass(*ml_fom.all_data(scores), v["working_point"], refit / "mass_punzi_cut.pdf",
                            f"DATA at the Punzi-optimal cut (signal eff. {wp['signal_eff']:.3f}, "
                            f"background eff. {wp['background_eff']:.4f})", self.case)
        ml_plots.punzi(scores, metrics["window_background"], v["working_point"], refit / "punzi.pdf",
                       f"Punzi FOM {wp['punzi']:.5f}", self.case)

    def record(self, status, best, totals):
        return {"name": self.name, "features": self.features, "base_features": self.info["base_features"],
                "added": self.info["added"], "status": status, "best": best, "scan": totals, "measured_utc": now()}

    def sculpting_pull(self):
        s = self.metrics["metrics"]["sculpting"]["working_point"]
        return abs(s["R_gap"] - 1.0) / s["R_gap_error"]

    def eligible(self):
        """Measured, passes the overtraining test, and does not sculpt the mass."""
        m = self.metrics
        return m is not None and m["status"] == "ok" and self.sculpting_pull() <= SCULPT_SIGMA

    def value(self, fom):
        return self.metrics["metrics"]["validation"][fom]

    def error(self, fom):
        return self.metrics["metrics"]["validation"]["errors"][fom]

    def describe(self):
        m = self.metrics
        if m is not None:
            if m["status"] != "ok":
                return m["status"]
            text = f"{CENTRAL_FOM} {self.value(CENTRAL_FOM):.4g} +- {self.error(CENTRAL_FOM):.2g}"
            return text if self.eligible() else text + f", sculpting (R_gap {self.sculpting_pull():.1f} sigma from 1)"
        scan = self.scan_status()
        text = f"scan {len(scan['done'])}/{self.info['jobs']}"
        if scan["lost"]:
            text += f" ({len(scan['lost'])} lost)"
        if not (scan["pending"] or scan["unsubmitted"] or scan["retry"]):
            text += ", refit " + self.refit_status()
        return text


# =============================================================================
# Commands for the scan and refit jobs
# =============================================================================

def scan_command(cand, seed, local):
    """(argv, extra environment) of one scan job. A local job gets absolute paths and the
    candidate's values; a batch job the queue variables of submit_scans."""
    m = METHODS[cand.method]
    if local:
        script, summary = str(m["dir"] / m["scan"]), str(cand.scan_dir / m["summary"].format(seed=seed))
        seed_text, trials, features = str(seed), str(cand.info["trials"]), ",".join(cand.features)
    else:
        script, summary = m["scan"], m["summary"].format(seed="$(seed)")
        seed_text, trials, features = "$(seed)", "$(trials)", "$(features)"
    if cand.method == "tmva":
        env = {"ML_CASE": cand.sample, "ML_FOM": cand.info["fom"], "ML_JOB": seed_text, "ML_RANDOM_SEARCH": "1",
               "ML_FEATURES": features, "ML_SUMMARY_PATH": summary}
        if not local:
            return ["-l", "-b", "-q", script], env
        call = script
        if cand.info["max_signal"]:
            call = f'{script}("{cand.sample}","{cand.info["fom"]}",{seed},{cand.info["max_signal"]},{cand.info["max_background"]})'
        return ["root", "-l", "-b", "-q", call], env
    argv = [script, "--case", cand.sample, "--fom", cand.info["fom"], "--trials", trials, "--seed", seed_text,
            "--features", features, "--summary-path", summary] + m["scan_args"]
    if local and cand.info["max_signal"]:
        argv += ["--max-signal", str(cand.info["max_signal"]), "--max-background", str(cand.info["max_background"])]
    return ([sys.executable] + argv if local else argv), {}


def refit_command(cand, best, local):
    m = METHODS[cand.method]
    script = str(m["dir"] / m["train"]) if local else m["train"]
    refit = str(cand.path / "refit") if local else "refit"
    features = ",".join(cand.features) if local else "$(features)"
    if cand.method == "tmva":
        env = {"ML_CASE": cand.sample, "ML_FOM": cand.info["fom"], "ML_USE_OPTIMIZED": "1",
               # Batch: only the best job's summary is transferred (condor on EOS cannot transfer a
               # folder); the best passing point in it is the best of the whole scan.
               "ML_OPTIMIZATION_DIR": str(cand.scan_dir) if local else ".",
               "ML_OUTPUT_DIR": refit, "ML_FEATURES": features}
        if not local:
            return ["-l", "-b", "-q", script], env
        call = script
        if cand.info["max_signal"]:
            call = f'{script}("{cand.sample}","{cand.info["fom"]}",{cand.info["max_signal"]},{cand.info["max_background"]})'
        return ["root", "-l", "-b", "-q", call], env
    params = str(cand.scan_dir / best["summary"]) if local else "$(params)"
    argv = [script, "--case", cand.sample, "--fom", cand.info["fom"], "--features", features,
            "--params-json", params, "--output-dir", refit]
    if local and cand.info["max_signal"]:
        argv += ["--max-signal", str(cand.info["max_signal"]), "--max-background", str(cand.info["max_background"])]
    return ([sys.executable] + argv if local else argv), {}


def submit_file(method, argv, env, inputs, outputs, initialdir, log, output, error):
    """Submit description without its queue line; batch jobs read the skims through XRootD."""
    m = METHODS[method]
    env = {**env, "ML_USE_XROOTD": "1", "ML_VARIANT": ml_samples.VARIANT}
    lines = [
        "universe = vanilla",
        f"executable = {m['dir'] / m['executable']}",
        f"initialdir = {initialdir}",
        "arguments = " + " ".join(argv),
        'environment = "' + " ".join(f"{k}={v}" for k, v in env.items()) + '"',
        f"transfer_input_files = {', '.join(inputs)}",
        f"transfer_output_files = {outputs}",
        f"output = {output}",
        f"error = {error}",
        f"log = {log}",
        "getenv = False",
        "should_transfer_files = YES",
        "when_to_transfer_output = ON_EXIT",
        f"request_cpus = {m['cpus']}",
        f"request_memory = {m['memory']}",
        "request_disk = 4 GB",
        'requirements = (OpSysAndVer == "AlmaLinux9")',
        f'+JobFlavour = "{m["flavour"]}"',
        # Held jobs are removed after an hour; the driver then resubmits them once.
        "periodic_remove = (JobStatus == 5) && ((time() - EnteredCurrentStatus) > 3600)",
    ]
    return "\n".join(lines) + "\n"


def method_inputs(method, script):
    """The files a job of `script` needs: with the TMVA headers, or with the trainer that the scan
    imports and ML_common."""
    m = METHODS[method]
    if method == "tmva":
        return [str(m["dir"] / script)] + [str(m["dir"] / h) for h in m["headers"]]
    scripts = [script] + ([m["train"]] if script != m["train"] else [])
    return [str(m["dir"] / s) for s in scripts] + [str(p) for p in COMMON_FILES]


def submit_scans(folder, jobs, dry_run):
    """One cluster for (candidate, seed) pairs of one method; returns the number of jobs."""
    if not jobs:
        return 0
    method = jobs[0][0].method
    tag = datetime.now().strftime("%Y%m%d_%H%M%S")
    items = folder / f"scan_{tag}.jobs"
    # features last: condor gives the last variable the rest of the line, commas included.
    # TMVA jobs have no trials (one random draw each).
    trials = (lambda c: "") if method == "tmva" else (lambda c: f"{c.info['trials']} ")
    items.write_text("".join(f"{c.path} {seed} {trials(c)}{','.join(c.features)}\n" for c, seed in jobs))
    argv, env = scan_command(jobs[0][0], 0, local=False)
    m = METHODS[method]
    # The job logs are kept only until the candidate is measured (Candidate.compact_scan). EOS-submit schedds
    # refuse output = /dev/null (every absolute path must be in /eos).
    text = submit_file(method, argv, env, method_inputs(method, m["scan"]), m["summary"].format(seed="$(seed)"),
                       "$(cand_dir)/scan", "condor.log", "logs/scan_$(seed).out", "logs/scan_$(seed).err")
    text += f"queue cand_dir, seed, {'' if method == 'tmva' else 'trials, '}features from {items}\n"
    sub = folder / f"scan_{tag}.sub"
    sub.write_text(text)
    print(f"  {len(jobs)} scan jobs: {sub}")
    if dry_run:
        return len(jobs)
    record_cluster(condor_submit(sub), jobs, lambda cand: cand.scan_dir / "submissions.json")
    sub.unlink()
    items.unlink()
    return len(jobs)


def submit_refits(folder, refits, dry_run):
    """One cluster for the refits of one method: (candidate, best trial) pairs."""
    if not refits:
        return 0
    method = refits[0][0].method
    m = METHODS[method]
    tag = datetime.now().strftime("%Y%m%d_%H%M%S")
    items = folder / f"refit_{tag}.jobs"
    # Every refit gets the summary of its best scan job; NN and XGB read it with --params-json.
    items.write_text("".join(f"{c.path} {b['summary']} {','.join(c.features)}\n" for c, b in refits))
    argv, env = refit_command(refits[0][0], refits[0][1], local=False)
    text = submit_file(method, argv, env, method_inputs(method, m["train"]) + ["scan/$(params)"],
                       "refit", "$(cand_dir)", "refit.log", "logs/refit.out", "logs/refit.err")
    text += f"queue cand_dir, params, features from {items}\n"
    sub = folder / f"refit_{tag}.sub"
    sub.write_text(text)
    print(f"  {len(refits)} refits: {sub}")
    if dry_run:
        return len(refits)
    cluster = condor_submit(sub)
    for proc, (cand, _) in enumerate(refits):
        submissions = read_json(cand.path / "refit_submissions.json", [])
        submissions.append({"cluster": cluster, "proc": proc, "submitted_utc": now()})
        write_json(cand.path / "refit_submissions.json", submissions)
    sub.unlink()
    items.unlink()
    return len(refits)


def progress(candidates, folder, dry_run):
    """Moves every candidate as far as it can go now; True when all are measured."""
    scans, refits = [], []
    for cand in candidates:
        if cand.metrics is not None:
            continue
        scan = cand.scan_status()
        if scan["unsubmitted"] or scan["retry"]:
            scans += [(cand, seed) for seed in scan["unsubmitted"] + scan["retry"]]
            continue
        if scan["pending"]:
            continue
        best, totals = cand.best_trial()
        if best is None:
            write_json(cand.path / "metrics.json", cand.record("no_passing_trial", None, totals))
            continue
        write_json(cand.path / "best.json", {**best, "scan": totals})
        state = cand.refit_status()
        if state == "done":
            print(f"  measuring {cand.path}")
            cand.measure(best, totals)
            cand.compact_scan()
        elif state in ("unsubmitted", "retry"):
            refits.append((cand, best))
        elif state == "lost":
            write_json(cand.path / "metrics.json", cand.record("refit_failed", best, totals))
    submit_scans(folder, scans, dry_run)
    submit_refits(folder, refits, dry_run)
    return all(c.metrics is not None for c in candidates)


# =============================================================================
# Studies
# =============================================================================

def trim(importance):
    """The inputs kept from fullPOOL/all: in decreasing importance up to the one that brings the
    cumulative share to IMPORTANCE_CUMULATIVE."""
    kept, cumulative = [], 0.0
    for name, share in sorted(importance.items(), key=lambda item: -item[1]):
        kept.append(name)
        cumulative += share
        if cumulative >= IMPORTANCE_CUMULATIVE:
            break
    return kept


def remove_jobs(paths):
    """condor_rm of every cluster recorded under the given folders."""
    clusters = sorted({s["cluster"] for p in paths for f in p.glob("**/*submissions.json") for s in read_json(f)})
    if clusters:
        # Finished clusters are no longer in the queue, so condor_rm reports them as not found.
        subprocess.run(["bash", "-lc", CONDOR_SETUP + "condor_rm " + " ".join(map(str, clusters))])


def brief(cand):
    """The numbers of one candidate that the comparisons and summaries show."""
    m = cand.metrics
    entry = {"path": str(cand.path), "features": cand.info["base_features"], "decorrelated": cand.info["decorrelated"],
             "status": m["status"]}
    if "metrics" in m:
        v = m["metrics"]["validation"]
        entry.update({
            "validation": {f: v[f] for f in ml_fom.FOMS}, "errors": v["errors"], "auc_gap": v["auc_gap"],
            "ks_pvalues": [v["ks_signal_pvalue"], v["ks_background_pvalue"]], "passed": v["passed"],
            "test_auc": m["metrics"]["test"]["auc"], "working_point": m["metrics"]["working_point"],
            "sculpting": m["metrics"]["sculpting"]["working_point"],
            "spectators": {k: s["efficiency"] for k, s in m["metrics"]["spectators"].items()},
            "importance": m["metrics"]["metadata"]["importance"],
            "importance_method": m["metrics"]["metadata"]["importance_method"],
        })
    return entry


class Study:
    def __init__(self, method, sample):
        self.method, self.sample = method, sample
        self.case = ml_samples.case(sample)
        variant = f"_{ml_samples.VARIANT}" if ml_samples.VARIANT else ""
        self.path = STUDY_DIR / self.case["tree"] / (self.case["year"] + variant) / method
        self.settings = read_json(self.path / "study.json")
        self._paired = None
        if self.settings is not None and self.settings["central_fom"] != CENTRAL_FOM:
            sys.exit(f"{self.path} is a {self.settings['central_fom']} study: use --fom {self.settings['central_fom']}")

    def setup(self):
        """What defines the candidates: when any of it changes, the study starts again from nothing."""
        c = self.case
        with uproot.open(ml_samples.XROOTD + c["inputs"]["data"]["skim"]) as skim:
            pre_cut, created = str(skim["metadata/pre_cut"]), str(skim["metadata/created_utc"])
        return to_json({
            "pool": c["pool"], "production_features": c["features"], "decorrelate": c["decorrelate"][self.method],
            "jobs": N_JOBS, "trials": N_TRIALS, "forward_selection": c["forward_selection"], "anchor_size": ANCHOR_SIZE,
            "max_added": MAX_ADDED, "importance_cumulative": IMPORTANCE_CUMULATIVE,
            "pre_cut": pre_cut, "skims_created_utc": created,
            "training_sidebands": c["training_sidebands"],
            # The test every scan trial must pass: trials judged by another rule are not comparable.
            "overtraining": {"ks_pvalue_min": settings.KS_PVALUE_MIN, "max_auc_gap": settings.MAX_AUC_GAP,
                             "auc_gap_sigma": settings.AUC_GAP_SIGMA},
        })

    def init(self):
        setup = self.setup()
        if self.settings is not None and self.settings.get("setup") != setup:
            trimming = ("importance_cumulative", "importance_min")
            if ({k: v for k, v in self.settings["setup"].items() if k not in trimming}
                    == {k: v for k, v in setup.items() if k not in trimming}):
                self.retrim(setup)
            else:
                print(f"  new setup: the study starts again, {self.path} is removed")
                self.remove()
        if self.settings is None:
            self.settings = {"method": self.method, "sample": self.sample, "central_fom": CENTRAL_FOM, "setup": setup,
                             "sculpt_sigma": SCULPT_SIGMA, "fom_settings": ml_fom.fom_settings(self.case), "created_utc": now()}
            write_json(self.path / "study.json", self.settings)
            print(f"Created {self.path}")
        return self

    def remove(self):
        """Its jobs still in condor and all its results."""
        remove_jobs([self.path])
        shutil.rmtree(self.path)
        self.settings, self._paired = None, None

    def retrim(self, setup):
        """A new trimming rule: when it changes the inputs of fullPOOL/trimmed, the jobs and results of
        fullPOOL/trimmed, the forward path, the final check and the summary are removed;
        fullPOOL/all and fullPOOL/production are kept."""
        folder = self.path / "fullPOOL"
        old = Candidate(folder / "trimmed").info if (folder / "trimmed" / "candidate.json").exists() else None
        measured = (folder / "all" / "metrics.json").exists() and Candidate(folder / "all").metrics["status"] == "ok"
        if old is not None and not (measured and trim(self.importance()) == old["base_features"]):
            print(f"  new trimming rule: fullPOOL/trimmed {old['base_features']}, the forward path and the "
                  f"final check start again")
            stale = [folder / "trimmed", self.path / "final"] + sorted(self.path.glob("step_*"))
            remove_jobs(stale)
            for path in stale:
                shutil.rmtree(path, ignore_errors=True)
            for name in ("summary.json", "summary.pdf", "paired_cache.json"):
                (self.path / name).unlink(missing_ok=True)
            self._paired = None
        self.settings["setup"] = setup
        write_json(self.path / "study.json", self.settings)

    def candidate(self, folder, name, features, decorrelate, role, added=None):
        s = self.settings["setup"]
        return Candidate.create(folder / name, self.method, self.sample, name, features, decorrelate, role, added,
                                s["jobs"], s["trials"])

    # ---- fullPOOL -----------------------------------------------------------------
    def importance(self):
        """Importance shares of fullPOOL/all, by feature name."""
        shares = Candidate(self.path / "fullPOOL" / "all").metrics["metrics"]["metadata"]["importance"]
        return {name.removesuffix(ml_samples.DECORRELATED_SUFFIX): share for name, share in shares.items()}

    def trimmed(self):
        return self.candidate(self.path / "fullPOOL", "trimmed", trim(self.importance()), self.settings["setup"]["decorrelate"],
                              "fullpool")

    def compare_trimmed(self):
        """fullPOOL/trimmed/comparison.json and comparison.pdf: what the trimming dropped, and all the
        metrics of fullPOOL/all and fullPOOL/trimmed side by side, with their paired differences."""
        folder = self.path / "fullPOOL"
        full, trimmed = Candidate(folder / "all"), Candidate(folder / "trimmed")
        importance = self.importance()
        comparison = {
            "rule": {"cumulative": IMPORTANCE_CUMULATIVE},
            "importance_of_all": importance, "kept": trimmed.info["base_features"],
            "dropped": [f for f in importance if f not in trimmed.info["base_features"]],
            "all": brief(full), "trimmed": brief(trimmed), "trimmed_minus_all": self.paired(trimmed.path, full.path),
        }
        write_json(folder / "trimmed" / "comparison.json", comparison)
        ml_plots.trimming(comparison, ml_fom.read_scores(full.path / "refit" / "scores.root"),
                          ml_fom.read_scores(trimmed.path / "refit" / "scores.root"),
                          folder / "trimmed" / "comparison.pdf", f"{self.method} {self.case['label']}", self.case)
        print(f"  fullPOOL: trimmed to {comparison['kept']}, dropped {comparison['dropped']}")

    # ---- forward selection -------------------------------------------------------
    def results(self):
        results, k = [], 0
        while (self.path / f"step_{k}" / "result.json").exists():
            results.append(read_json(self.path / f"step_{k}" / "result.json"))
            k += 1
        return results

    def path_complete(self, results):
        """The path ends after MAX_ADDED steps (or when the pool is used up or no candidate of a step
        is eligible)."""
        if not results:
            return False
        last = results[-1]
        return (last["winner"] is None or len(results) > MAX_ADDED
                or not [f for f in self.settings["setup"]["pool"] if f not in last["features"]])

    def step_candidates(self, k, results):
        s = self.settings["setup"]
        if k == 0:
            ranked = sorted(self.importance().items(), key=lambda item: -item[1])
            specs = [("anchor", [name for name, _ in ranked[:s["anchor_size"]]], None)]
        else:
            base = results[k - 1]["features"]
            specs = [(f"add_{f}", base + [f], f) for f in s["pool"] if f not in base]
        return [self.candidate(self.path / f"step_{k}", name, features, s["decorrelate"], "path", added)
                for name, features, added in specs]

    def close_step(self, k, candidates, results):
        eligible = sorted([c for c in candidates if c.eligible()], key=lambda c: -c.value(CENTRAL_FOM))
        status = lambda c: "sculpting" if c.metrics["status"] == "ok" and not c.eligible() else c.metrics["status"]
        ranking = [{"name": c.name, "added": c.info["added"], "status": status(c),
                    "R_gap_pull": c.sculpting_pull() if c.metrics["status"] == "ok" else None,
                    "values": {f: c.value(f) for f in ml_fom.FOMS} if c.eligible() else None,
                    "errors": c.metrics["metrics"]["validation"]["errors"] if c.eligible() else None}
                   for c in eligible + [c for c in candidates if not c.eligible()]]
        result = {"step": k, "ranking": ranking, "winner": None, "features": None, "closed_utc": now()}
        if eligible:
            winner = eligible[0]
            result.update({"winner": winner.name, "features": winner.info["base_features"],
                           "added": winner.info["added"], "path": str(winner.path.relative_to(self.path)),
                           "versus_winner": {c.name: self.paired(winner.path, c.path) for c in eligible[1:]}})
            if k > 0:
                result["gain"] = self.paired(winner.path, self.path / results[k - 1]["path"])
            result["within_sigma"] = [name for name, g in result["versus_winner"].items()
                                      if g[CENTRAL_FOM]["gain"] <= WITHIN_SIGMA * g[CENTRAL_FOM]["sigma"]]
        write_json(self.path / f"step_{k}" / "result.json", result)
        text = f"{result['winner']} ({winner.value(CENTRAL_FOM):.4g})" if eligible else "none eligible, path ends"
        print(f"  step {k} closed: winner {text}")

    # ---- paired gains, cached per study -----------------------------------------
    def paired(self, new_path, old_path):
        """Paired validation gain of the refit in new_path over the one in old_path, every FOM. The
        cache holds the time stamps of both scores files, so a retrained candidate is recomputed."""
        cache_path = self.path / "paired_cache.json"
        if self._paired is None:
            self._paired = read_json(cache_path, {})
        new, old = Candidate(new_path), Candidate(old_path)
        key = f"{new.path.relative_to(self.path)}|{old.path.relative_to(self.path)}"
        stamps = [os.path.getmtime(c.path / "refit" / "scores.root") for c in (new, old)]
        entry = self._paired.get(key, {})
        fom_settings = ml_fom.fom_settings(self.case)
        if (entry.get("scores_mtime"), entry.get("fom_settings")) != (stamps, fom_settings):
            scores = ml_fom.read_scores(new.path / "refit" / "scores.root")
            gains, n_common = ml_fom.paired_gain(scores["validation"], ml_fom.read_scores(old.path / "refit" / "scores.root")["validation"],
                                                 ml_fom.window_background(ml_fom.all_data(scores)[0], self.case), self.case)
            self._paired[key] = {**gains, "n_common": n_common, "scores_mtime": stamps, "fom_settings": fom_settings}
            write_json(cache_path, self._paired)
        return self._paired[key]

    # ---- plateau, recommended set, replay ---------------------------------------
    def analyse(self, results, fom):
        path = [r for r in results if r["winner"]]
        values = [Candidate(self.path / r["path"]).value(fom) for r in path]
        gains = [None] + [self.paired(self.path / r["path"], self.path / path[k - 1]["path"])[fom]
                          for k, r in enumerate(path) if k > 0]
        significant = [None] + [g["gain"] >= GAIN_SIGMA * g["sigma"] for g in gains[1:]]
        plateau = next((k for k in range(1, len(path))
                        if not any(significant[j] for j in range(k, min(k + PATIENCE, len(path))))), None)
        best = int(np.argmax(values))
        distance = [self.paired(self.path / path[best]["path"], self.path / r["path"])[fom] if k != best
                    else {"gain": 0.0, "sigma": 0.0} for k, r in enumerate(path)]
        recommended = next(k for k, d in enumerate(distance) if d["gain"] <= WITHIN_SIGMA * d["sigma"])
        # The candidate this FOM would have picked at each step, among those the CENTRAL_FOM path offered.
        choices = []
        for r in path:
            ranked = [c for c in r["ranking"] if c["values"]]
            pick = max(ranked, key=lambda c: c["values"][fom])["name"]
            choices.append({"step": r["step"], "path_winner": r["winner"], "choice": pick,
                            "diverges": pick != r["winner"],
                            "winner_minus_choice": r["versus_winner"][pick][fom] if pick != r["winner"] else None})
        return {
            "fom": fom,
            "steps": [{"step": r["step"], "features": r["features"], "added": r.get("added"), "value": values[k],
                       "error": Candidate(self.path / r["path"]).error(fom), "gain": gains[k],
                       "significant": significant[k], "distance_to_best": distance[k]} for k, r in enumerate(path)],
            "plateau_step": plateau,
            "stop_features": path[(plateau - 1) if plateau is not None else -1]["features"],
            "best_step": best,
            "recommended_step": recommended,
            "recommended_features": path[recommended]["features"],
            "first_divergence": next((c["step"] for c in choices if c["diverges"]), None),
            "choices": choices,
        }

    # ---- summary ---------------------------------------------------------------
    def forward_path(self, results):
        """The step winners of the forward selection: none without it, or when no candidate of
        step 0 is eligible."""
        return [r for r in results if r["winner"]] if self.settings["setup"]["forward_selection"] else []

    def reference(self, results):
        """The last tentative selection: the recommended set of the forward selection, or
        fullPOOL/trimmed when there is no forward path."""
        path = self.forward_path(results)
        if not path:
            return self.path / "fullPOOL" / "trimmed"
        return self.path / path[self.analyse(results, CENTRAL_FOM)["recommended_step"]]["path"]

    def production_name(self):
        """The candidate of the production inputs: fullPOOL/all when they are its inputs, raw for the method
        (B+ XGB), otherwise fullPOOL/production."""
        s = self.settings["setup"]
        return "all" if sorted(s["production_features"]) == sorted(s["pool"]) and not s["decorrelate"] else "production"

    def summarize(self, results):
        """summary.json and summary.pdf: the tentative selections (fullPOOL/all, fullPOOL/trimmed and,
        with the forward selection, its recommended set), the references (production inputs, the
        last selection with the other decorrelation setting), and the forward path with its replay per FOM."""
        folder = self.path / "fullPOOL"
        reference = self.reference(results)
        final = "no_decorrelation" if self.settings["setup"]["decorrelate"] else "decorrelation"
        cands = {"all": Candidate(folder / "all"), "trimmed": Candidate(folder / "trimmed")}
        summary = {"method": self.method, "sample": self.sample, "central_fom": CENTRAL_FOM, "updated_utc": now(),
                   "trimming": read_json(folder / "trimmed" / "comparison.json")}
        if self.forward_path(results):
            summary["analysis"] = {fom: self.analyse(results, fom) for fom in ml_fom.FOMS}
            cands["forward"] = Candidate(reference)
        summary["selections"] = {name: cand.info["base_features"] for name, cand in cands.items()}
        cands["production"] = Candidate(folder / self.production_name())
        if self.case["final_check"]:
            cands[final] = Candidate(self.path / "final" / final)
        summary["candidates"] = {name: {**brief(cand), "reference_minus_this": self.paired(reference, cand.path)
                                        if cand.path != reference and "metrics" in cand.metrics else None}
                                 for name, cand in cands.items()}
        summary["reference"] = str(reference.relative_to(self.path))
        write_json(self.path / "summary.json", summary)
        self.plot_summary(results, summary)
        return summary

    def plot_summary(self, results, summary):
        """summary.pdf: the selections and references side by side; with the forward selection, its
        path in every FOM, the paired gain of every step, and the sculpting of the step winners."""
        with PdfPages(self.path / "summary.pdf") as pdf:
            measured = {n: c for n, c in summary["candidates"].items() if "validation" in c}
            y = np.arange(len(measured))
            fig, (left, middle, right) = plt.subplots(1, 3, figsize=(15, 4 + 0.3 * len(measured)))
            left.errorbar([c["validation"][CENTRAL_FOM] for c in measured.values()], y,
                          xerr=[c["errors"][CENTRAL_FOM] for c in measured.values()], fmt="o", color="black")
            left.set_xlabel(f"validation {CENTRAL_FOM}")
            middle.plot([c["working_point"]["punzi"] for c in measured.values()], y, "o", color="black")
            middle.set_xlabel("Punzi FOM at the working point")
            right.errorbar([c["sculpting"]["R_gap"] for c in measured.values()], y,
                           xerr=[c["sculpting"]["R_gap_error"] for c in measured.values()], fmt="o", color="tab:green")
            right.axvline(1.0, color="grey", lw=0.8)
            right.set_xlabel("R_gap at the working point")
            labels = [f"{n} ({len(c['features'])})" for n, c in measured.items()]
            for ax in (left, middle, right):
                ax.set_yticks(y)
                ax.set_yticklabels(labels, fontsize=9)
                ax.invert_yaxis()
                ax.grid(axis="x", alpha=0.3)
            fig.suptitle(f"{self.method} {self.case['label']}: selections and references (number of inputs)")
            fig.tight_layout()
            ml_plots.label_sample(fig, self.case["label"])
            pdf.savefig(fig)
            plt.close(fig)
            if "analysis" not in summary:
                return
            path = [r for r in results if r["winner"]]
            labels = ["anchor"] + ["+" + r["added"] for r in path[1:]]
            x = np.arange(len(path))
            winners = [Candidate(self.path / r["path"]).metrics["metrics"] for r in path]
            for page in ("values", "gains"):
                fig, axes = plt.subplots(2, 2, figsize=(13, 9))
                for ax, fom in zip(axes.flat, ml_fom.FOMS):
                    a = summary["analysis"][fom]
                    if page == "values":
                        for k, r in enumerate(path):
                            others = [c["values"][fom] for c in r["ranking"] if c["values"] and c["name"] != r["winner"]]
                            ax.plot([k] * len(others), others, "o", color="lightgrey", ms=4)
                        ax.errorbar(x, [s["value"] for s in a["steps"]], [s["error"] for s in a["steps"]], fmt="o-",
                                    color="black", ms=4, label=f"step winner (chosen by {CENTRAL_FOM})")
                        rec = a["recommended_step"]
                        ax.plot(rec, a["steps"][rec]["value"], "*", color="tab:red", ms=14, label="recommended")
                        ax.set_ylabel(f"validation {fom}")
                    else:
                        for k, g in zip(x[1:], [s["gain"] for s in a["steps"][1:]]):
                            filled = g["gain"] >= GAIN_SIGMA * g["sigma"]
                            ax.errorbar(k, g["gain"], g["sigma"], fmt="o", color="black", mfc="black" if filled else "white")
                        ax.axhline(0.0, color="grey", lw=0.8, label="no gain")
                        ax.set_ylabel(f"paired gain in {fom} (open: < {GAIN_SIGMA:g} sigma)")
                    if a["plateau_step"] is not None:
                        ax.axvline(a["plateau_step"], color="tab:orange", ls="--", label="plateau")
                    ax.set_xticks(x)
                    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
                    ax.set_title(fom)
                    ax.legend(fontsize=7)
                fig.suptitle(f"{self.method} {self.case['label']}: " + ("validation FOM per step" if page == "values"
                             else "gain of each step winner over the previous one"))
                fig.tight_layout()
                ml_plots.label_sample(fig, self.case["label"])
                pdf.savefig(fig)
                plt.close(fig)
            fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5))
            wp = [m["sculpting"]["working_point"] for m in winners]
            start, width = self.case["punzi_sideband_start"], self.case["punzi_sideband_width"]
            sidebands = f"{1000 * start:g}-{1000 * (start + width):g} MeV"
            for key, color, name in (("R_gap", "tab:green", "gaps"), ("R_punzi_sidebands", "tab:orange", sidebands)):
                left.errorbar(x, [s[key] for s in wp], [s[key + "_error"] for s in wp], fmt="o-", color=color, label=name)
            left.axhline(1.0, color="grey", lw=0.8)
            left.set_ylabel("background eff. / line prediction, working point")
            right.errorbar(x, [100 * s["slope_per_100MeV_relative"] for s in wp],
                           [100 * s["slope_per_100MeV_relative_error"] for s in wp], fmt="o-", color="black")
            right.axhline(0.0, color="grey", lw=0.8)
            right.set_ylabel("background eff. slope [% / 100 MeV], working point")
            for ax in (left, right):
                ax.set_xticks(x)
                ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
            left.legend(fontsize=8)
            fig.suptitle(f"{self.method} {self.case['label']}: sculpting of the step winners")
            fig.tight_layout()
            ml_plots.label_sample(fig, self.case["label"])
            pdf.savefig(fig)
            plt.close(fig)

    # ---- metric settings changed during the study ------------------------------------
    def revise(self):
        """Brings a study to the current SCULPT_SIGMA and FOM settings (ml_fom.fom_settings), which change
        only metrics: every measured candidate is measured again from its stored scores (nothing is
        trained again), the path steps are closed again, and from the first step whose winner
        changes on, the later steps and the final check are removed, so that advance rebuilds them.
        The settings are written last: an interrupted revision is simply done again."""
        current = {"sculpt_sigma": SCULPT_SIGMA, "fom_settings": ml_fom.fom_settings(self.case)}
        changed = {k: v for k, v in current.items() if self.settings.get(k) != v}
        if not changed:
            return
        print(f"  revising the study to {changed}")
        old = self.results()
        self.settings.update(changed)
        measured = [(info.parent, read_json(info)) for info in sorted(self.path.glob("*/*/metrics.json"))]
        stale = [path for path, m in measured if "metrics" in m and m["metrics"].get("fom_settings") != current["fom_settings"]]
        if stale:
            print(f"  measuring {len(stale)} candidates again with the current FOM settings")
            with multiprocessing.Pool(MEASURE_WORKERS) as pool:
                pool.map(measure_again, stale)
        results, keep = [], len(old)
        for k, r in enumerate(old):
            if self.path_complete(results):
                keep = k
                break
            candidates = self.step_candidates(k, results)
            if not all(c.metrics is not None for c in candidates):
                keep = k + 1
                break
            self.close_step(k, candidates, results)
            results = self.results()[:k + 1]
            if results[k]["winner"] != r["winner"]:
                keep = k + 1
                break
        for step in self.path.glob("step_*"):
            if int(step.name.split("_")[1]) >= keep:
                shutil.rmtree(step)
        for k in range(len(results), keep):   # a step kept open: its candidates are not all measured
            (self.path / f"step_{k}" / "result.json").unlink(missing_ok=True)
        if [r["winner"] for r in results] != [r["winner"] for r in old]:
            if (self.path / "final").exists():
                shutil.rmtree(self.path / "final")
        for path in (self.path / "summary.json", self.path / "summary.pdf",
                     self.path / "fullPOOL" / "trimmed" / "comparison.json", self.path / "fullPOOL" / "trimmed" / "comparison.pdf"):
            path.unlink(missing_ok=True)
        write_json(self.path / "study.json", self.settings)

    # ---- one pass ---------------------------------------------------------------
    def advance(self, dry_run=False):
        """Collects what finished and submits what can start; returns True when the study is done."""
        self.init()
        print(f"== {self.method} {self.sample}")
        self.revise()
        s, folder = self.settings["setup"], self.path / "fullPOOL"
        full = self.candidate(folder, "all", s["pool"], s["decorrelate"], "fullpool")
        production = (full if self.production_name() == "all" else
                      self.candidate(folder, "production", s["production_features"], False, "production"))
        production_done = progress([production], folder, dry_run)
        if not progress([full], folder, dry_run):
            return False
        if "metrics" not in full.metrics:
            # No trial of fullPOOL/all passes the overtraining test: nothing to trim or to start from.
            if production_done and not (self.path / "summary.json").exists():
                write_json(self.path / "summary.json", {
                    "method": self.method, "sample": self.sample, "updated_utc": now(),
                    "stopped": f"fullPOOL/all: {full.metrics['status']}",
                    "candidates": {"all": brief(full), "production": brief(production)}})
                print(f"  stopped: fullPOOL/all has no model ({full.metrics['status']})")
            return production_done
        trimmed = self.trimmed()
        if not progress([trimmed], folder, dry_run):
            return False
        if not (folder / "trimmed" / "comparison.json").exists():
            self.compare_trimmed()
        results = self.results()
        while s["forward_selection"] and not self.path_complete(results):
            k = len(results)
            candidates = self.step_candidates(k, results)
            if not progress(candidates, self.path / f"step_{k}", dry_run):
                return False
            self.close_step(k, candidates, results)
            results = self.results()
        if self.case["final_check"]:
            reference = Candidate(self.reference(results))
            final_name = "no_decorrelation" if s["decorrelate"] else "decorrelation"
            final = self.candidate(self.path / "final", final_name, reference.info["base_features"], not s["decorrelate"], final_name)
            if not progress([final], self.path / "final", dry_run):
                return False
        if not production_done:
            return False
        if not (self.path / "summary.json").exists():
            summary = self.summarize(results)
            print(f"  done: selections {summary['selections']}")
        return True

    def status(self):
        print(f"== {self.method} {self.sample}: {self.path}")
        if self.settings is None:
            print("  not started")
            return
        for path in sorted((self.path / "fullPOOL").glob("*/candidate.json")):
            print(f"    fullPOOL/{path.parent.name} {Candidate(path.parent).info['base_features']}: "
                  f"{Candidate(path.parent).describe()}")
        results = self.results()
        for r in results:
            text = f"{r['winner']} {r['ranking'][0]['values'][CENTRAL_FOM]:.4g}" if r["winner"] else "no winner"
            if r.get("gain"):
                g = r["gain"][CENTRAL_FOM]
                text += f", gain {g['gain']:+.2g} +- {g['sigma']:.2g}"
            print(f"  step {r['step']}: {text}")
        for path in sorted((self.path / f"step_{len(results)}").glob("*/candidate.json")) + sorted((self.path / "final").glob("*/candidate.json")):
            print(f"    {path.parent.parent.name}/{path.parent.name}: {Candidate(path.parent).describe()}")
        summary = read_json(self.path / "summary.json")
        if summary:
            print(f"  finished: {summary.get('stopped') or summary['selections']}")


# =============================================================================
# Command line
# =============================================================================

def one_candidate(args):
    """One candidate end to end, outside any study: args.jobs scan jobs, the refit of the best
    passing trial, and the metrics. Locally the jobs run one after the other. With --condor the
    command does one pass like advance (submit the scan, then the refit, then measure), so it is
    rerun until the metrics are printed; condor jobs always use the full samples."""
    base = args.features.split(",")
    name = args.name or "set_" + "_".join(base)
    case_name = f"{args.tree[0]}/{args.year[0]}"
    variant = f"_{ml_samples.VARIANT}" if ml_samples.VARIANT else ""
    folder = Path(args.output_dir or STUDY_DIR / args.tree[0] / (args.year[0] + variant) / "local" / args.method / name)
    jobs = args.jobs or (N_JOBS if args.condor else 2)
    trials = args.trials or (N_TRIALS if args.condor else 3)
    cand = Candidate.create(folder, args.method, case_name, name, base, args.decorrelate, "local",
                            None, jobs, trials, args.max_signal, args.max_background)
    print(f"Candidate {name}: {cand.features}, {cand.info['jobs']} jobs x {cand.info['trials']} trials, in {folder}")
    if args.close_scan:
        cand.close_scan()
    if args.condor:
        if not progress([cand], cand.path.parent, args.dry_run):
            print(f"  {cand.describe()}. Rerun this command to collect and go on.")
            return
    else:
        # The scan jobs one after the other, then the refit of the best passing trial.
        for seed in range(cand.info["jobs"]):
            if not (cand.scan_dir / METHODS[cand.method]["summary"].format(seed=seed)).exists():
                argv, env = scan_command(cand, seed, local=True)
                log = cand.scan_dir / "logs" / f"scan_{seed}.log"
                print("  $", " ".join(argv), f"(log: {log})")
                with open(log, "w") as out:
                    subprocess.run(argv, cwd=cand.scan_dir, env={**os.environ, **env}, stdout=out, stderr=subprocess.STDOUT, check=True)
        best, totals = cand.best_trial()
        print(f"Scan: {totals}; best {best and best['value']}")
        if best is None:
            write_json(cand.path / "metrics.json", cand.record("no_passing_trial", None, totals))
        else:
            write_json(cand.path / "best.json", {**best, "scan": totals})
            if not (cand.path / "refit" / "scores.root").exists():
                argv, env = refit_command(cand, best, local=True)
                log = cand.path / "logs" / "refit.log"
                print("  $", " ".join(argv), f"(log: {log})")
                with open(log, "w") as out:
                    subprocess.run(argv, cwd=cand.path, env={**os.environ, **env}, stdout=out, stderr=subprocess.STDOUT, check=True)
            cand.measure(best, totals)
    m = cand.metrics
    print(f"Scan: {m['scan']}")
    if "metrics" not in m:
        print(f"Status {m['status']}")   # no_passing_trial or refit_failed
        return
    v = m["metrics"]["validation"]
    print(f"Status {m['status']}: " + ", ".join(f"{f} {v[f]:.4f} +- {v['errors'][f]:.4f}" for f in ml_fom.FOMS))
    for label, s in m["metrics"]["sculpting"].items():
        print(f"  sculpting {label}: R_gap {s['R_gap']:.3f} +- {s['R_gap_error']:.3f}, "
              f"slope {s['slope_per_100MeV_relative']:+.3f} +- {s['slope_per_100MeV_relative_error']:.3f} per 100 MeV")
    print(f"Metrics: {cand.path / 'metrics.json'}")


def adopt(study, selection):
    """The chosen selection of a study becomes the production FEATURES of its tree in settings.py."""
    summary = read_json(study.path / "summary.json")
    base = summary["selections"][selection]
    features = [f + ml_samples.DECORRELATED_SUFFIX for f in base] if study.settings["setup"]["decorrelate"] else base
    path = SELECTION_DIR / "settings.py"
    text = path.read_text()
    block = re.search(r"^FEATURES = \{\n.*?^\}", text, flags=re.M | re.S)
    line = f'    "{study.case["tree"]}": {json.dumps(features)},'
    new_block = re.sub(rf'^    "{study.case["tree"]}": \[.*\],$', line, block.group(0), flags=re.M)
    path.write_text(text[:block.start()] + new_block + text[block.end():])
    print(f"FEATURES of {study.case['tree']} in {path} set to {features} ({selection} of {study.path}):\n{line}")
    print("Next: ./fresh_ML_scan.sh --sync-only, then the scans as usual.")


def studies_from(args, existing_only):
    """The campaign file's studies, or every method of its tree (settings.ML_METHODS) x case, narrowed by
    --method, --tree and --year."""
    if args.campaign:
        lines = [l.split("#")[0].split() for l in Path(args.campaign).read_text().splitlines()]
        return [Study(*l) for l in lines if l]
    studies = [Study(m, name) for m in ([args.method] if args.method else METHODS) for name in ml_samples.CASES
               if m in settings.ML_METHODS[name.split("/")[0]]
               and (not args.tree or name.split("/")[0] in args.tree) and (not args.year or name.split("/")[1] in args.year)]
    return [s for s in studies if s.settings is not None] if existing_only else studies


def main():
    global CENTRAL_FOM, STUDY_DIR
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("advance", "run", "status", "candidate", "replay", "replot", "adopt", "compact"))
    parser.add_argument("--method", choices=METHODS)
    parser.add_argument("--tree", action="append", help="Only these trees (ntmix, ntKp, ntphi, ntKstar); repeatable.")
    parser.add_argument("--year", action="append", help="Only these years (PbPb18, PbPb23, PbPb24, ppRef); repeatable.")
    parser.add_argument("--campaign", help="File with one 'method tree/year' per line.")
    parser.add_argument("--fom", choices=ml_fom.FOMS, default=CENTRAL_FOM,
                        help=f"Central FOM: tunes the candidates and chooses the forward steps (default {CENTRAL_FOM}).")
    parser.add_argument("--study-dir", help="Folder of the studies, in ML_selection/ (default feature_study): "
                                            "a finished round renamed for reference is advanced or replotted there.")
    parser.add_argument("--dry-run", action="store_true", help="Write the submit files, submit nothing.")
    parser.add_argument("--features", help="candidate: comma-separated base features.")
    parser.add_argument("--selection", choices=("all", "trimmed", "forward"), help="adopt: which selection.")
    parser.add_argument("--name", help="candidate: folder name.")
    parser.add_argument("--output-dir", help="candidate: folder (default feature_study/<tree>/<year>/local/<method>/<name>).")
    parser.add_argument("--condor", action="store_true", help="candidate: run the scan and refit as condor jobs.")
    parser.add_argument("--close-scan", action="store_true",
                        help="candidate: stop waiting for the last scan jobs (removes them) and go on with the rest.")
    parser.add_argument("--jobs", type=int, help=f"candidate: scan jobs (default 2 locally, {N_JOBS} on condor).")
    parser.add_argument("--trials", type=int, help=f"candidate: trials per job, xgb and nn (default 3 locally, {N_TRIALS} on condor).")
    parser.add_argument("--max-signal", type=int, help="candidate: fewer signal candidates for a quick test.")
    parser.add_argument("--max-background", type=int, help="candidate: fewer background candidates.")
    parser.add_argument("--decorrelate", action="store_true", help="candidate: the _dc features (default: raw).")
    args = parser.parse_args()
    CENTRAL_FOM = args.fom
    STUDY_DIR = SELECTION_DIR / (args.study_dir or STUDY_DIR.name)
    if args.command in ("candidate", "adopt") and not (args.method and args.tree and args.year):
        parser.error(f"{args.command} needs --method, --tree and --year")
    if (args.command, args.selection) == ("adopt", None) or (args.command, args.features) == ("candidate", None):
        parser.error("adopt needs --selection, candidate needs --features")
    if args.command == "candidate" and args.condor and (args.max_signal or args.max_background):
        parser.error("condor jobs use the full samples; --max-signal and --max-background are local only")
    ml_samples.write_cpp_settings()
    if args.command == "candidate":
        return one_candidate(args)
    studies = studies_from(args, existing_only=args.command in ("status", "replay", "replot", "adopt", "compact"))
    if args.command == "advance":
        for study in studies:
            study.advance(args.dry_run)
    elif args.command == "run":
        while not all([study.advance(args.dry_run) for study in studies]) and not args.dry_run:
            print(f"-- next pass in {POLL_MINUTES} min ({now()})", flush=True)
            time.sleep(POLL_MINUTES * 60)
    elif args.command == "status":
        if not studies:
            print(f"No study in {STUDY_DIR} yet: ./feature_study.py advance (or run)")
        for study in studies:
            study.status()
    elif args.command == "replay":
        for study in studies:
            # A finished study from its summary, otherwise the forward path so far.
            summary = read_json(study.path / "summary.json")
            analysis = summary.get("analysis", {}) if summary else {
                fom: study.analyse(study.results(), fom) for fom in ml_fom.FOMS if study.forward_path(study.results())}
            for fom, a in analysis.items():
                print(f"{study.method} {study.sample} {fom:6s}: plateau step {a['plateau_step']}, recommended step "
                      f"{a['recommended_step']} {a['recommended_features']}, first divergence {a['first_divergence']}")
    elif args.command == "replot":
        # After a change in methods/ML_common/ml_plots.py: every measured candidate's plots, from its stored
        # scores and metrics (nothing is trained or measured again), and the study summaries.
        for study in studies:
            for info in sorted(study.path.glob("*/*/metrics.json")):
                cand = Candidate(info.parent)
                if "metrics" in cand.metrics:
                    cand.plot(ml_fom.read_scores(cand.path / "refit" / "scores.root"), cand.metrics["metrics"])
            if (study.path / "fullPOOL" / "trimmed" / "comparison.json").exists():
                study.compare_trimmed()
            if (study.path / "summary.json").exists():
                study.plot_summary(study.results(), read_json(study.path / "summary.json"))
            print(f"{study.method} {study.sample}: plots redrawn")
    elif args.command == "adopt":
        adopt(studies[0], args.selection)
    elif args.command == "compact":
        # Studies made before the compaction existed: every measured candidate's scan folder, and the
        # submit files of finished passes.
        for study in studies:
            for info in sorted(study.path.glob("**/metrics.json")):
                Candidate(info.parent).compact_scan()
            for path in list(study.path.glob("**/*.sub")) + list(study.path.glob("**/*.jobs")):
                path.unlink()
            print(f"{study.method} {study.sample}: compacted")


if __name__ == "__main__":
    main()
