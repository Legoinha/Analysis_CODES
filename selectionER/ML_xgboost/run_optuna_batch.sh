#!/bin/bash
# EOS-resident launcher for the XGBoost workflow in condor jobs.
# Usage: ./run_optuna_batch.sh XGB_optuna.py --sample PbPb23 --fom auc
LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt
source "${LCG_VIEW}/setup.sh"
set -eo pipefail
export ML_USE_XROOTD=1
export XGB_BATCH_QUIET=1
export PYTHONDONTWRITEBYTECODE=1

exec python "$@"
