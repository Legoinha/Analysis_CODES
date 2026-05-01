#!/bin/bash
set -euo pipefail

cd /eos/user/h/hmarques/Analysis_CODES/selectionER
source .envs/xgb_train/bin/activate

RUN_ID="${1:-0}"
SEED_BASE="${SEED_BASE:-42}"
SEED=$((SEED_BASE + RUN_ID))

TRIALS=400
N_ROUNDS=5000
EARLY_STOPPING=100
MAX_BACKGROUND=500000
SCAN_TAG="${SCAN_TAG:-}"
SUMMARY_SUBDIR="${SUMMARY_SUBDIR:-optuna_summaries}"
DROP_FEATURES="${DROP_FEATURES:-}"

EXTRA_ARGS=()
if [[ -n "${SCAN_TAG}" ]]; then
  EXTRA_ARGS+=(--scan-tag "${SCAN_TAG}")
fi
if [[ -n "${SUMMARY_SUBDIR}" ]]; then
  EXTRA_ARGS+=(--summary-subdir "${SUMMARY_SUBDIR}")
fi
if [[ -n "${DROP_FEATURES}" ]]; then
  IFS=',' read -ra FEATURES_TO_DROP <<< "${DROP_FEATURES}"
  for FEATURE in "${FEATURES_TO_DROP[@]}"; do
    EXTRA_ARGS+=(--drop-feature "${FEATURE}")
  done
fi

python XGB_optuna.py \
  --sample ppRef24 \
  --trials "${TRIALS}" \
  --n-rounds "${N_ROUNDS}" \
  --early-stopping "${EARLY_STOPPING}" \
  --max-background "${MAX_BACKGROUND}" \
  --seed "${SEED}" \
  "${EXTRA_ARGS[@]}"
