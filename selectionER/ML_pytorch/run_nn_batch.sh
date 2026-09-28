#!/bin/bash
# EOS-resident launcher for the PyTorch NN workflow.
# Usage: ./run_nn_batch.sh NN_train.py --sample PbPb23
# Override the LCG view (e.g. on an EL8 node) with NN_LCG_VIEW.

LCG_VIEW=${NN_LCG_VIEW:-/cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt}
# Source before enabling errexit: the LCG setup runs `gnuplot --version | cut`,
# which fails on nodes without libreadline.so.8 and would abort the launcher.
source "${LCG_VIEW}/setup.sh"
set -eo pipefail
export ML_USE_XROOTD=${ML_USE_XROOTD:-1}
export NN_BATCH_QUIET=${NN_BATCH_QUIET:-1}
export PYTHONDONTWRITEBYTECODE=1

exec python "$@"
