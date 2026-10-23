#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_prompt_plus_NP_Validation.sh <X3872|PSI2S> <system>
#
# Required beforehand: the inclusive and B-enriched fits, followed by
# ../run_nonPrompt.sh to produce the measured inclusive nonprompt fraction.
# Every comparison applies pThatreweight to every MC entry. Prompt and
# nonprompt comparisons use complementary lxy selections read from the
# fraction metadata; the inclusive comparison uses their measured mixture.

PARTICLE="${1:?Missing particle: X3872 or PSI2S}"
SYSTEM="${2:?Missing system}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NONPROMPT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
VALIDATION_DIR="$(cd "${SCRIPT_DIR}/../../Validation" && pwd)"
BASE="/eos/user/h/hmarques/Analysis_CODES"

case "$PARTICLE" in
  X3872)
    TREE="ntmix_X3872"
    FRACTION_HIST="hX3872_dataDriven_nonprompt_fraction_inclusive"
    ;;
  PSI2S|psi2s)
    TREE="ntmix_PSI2S"
    PARTICLE="PSI2S"
    FRACTION_HIST="hPsi2S_dataDriven_nonprompt_fraction_inclusive"
    ;;
esac

NONPROMPT_MODEL="${BASE}/fitER/ROOTfiles/${SYSTEM}_nonPrompt/nominalFitModel_${TREE}_${SYSTEM}_nonPrompt.root"
FRACTION="${NONPROMPT_DIR}/${SYSTEM}/nonprompt_fraction_inclusive.root"
NONPROMPT_OUTPUT="${SCRIPT_DIR}/${SYSTEM}/${TREE}/nonprompt"
PROMPT_OUTPUT="${SCRIPT_DIR}/${SYSTEM}/${TREE}/prompt"
INCLUSIVE_OUTPUT="${SCRIPT_DIR}/${SYSTEM}/${TREE}/inclusive"

# Keep the generic inclusive-data/prompt-MC validation and its weight payload.
bash "${VALIDATION_DIR}/run_DataSIGNAL_VS_MC.sh" "$TREE" "$SYSTEM"

# B-enriched data and nonprompt MC: both use lxy > cut.
bash "${VALIDATION_DIR}/run_DataSIGNAL_VS_MC.sh" \
  "$TREE" "${SYSTEM}_nonPrompt" "" "" 0 \
  "$NONPROMPT_OUTPUT" "${SYSTEM} B-enriched Signal" "MC np" 0 \
  nonprompt "$FRACTION" "" "" 0

# B-suppressed data and inclusive MC: both use lxy < cut.
bash "${VALIDATION_DIR}/run_DataSIGNAL_VS_MC.sh" \
  "$TREE" "$SYSTEM" "" "" 0 \
  "$PROMPT_OUTPUT" "${SYSTEM} B-suppressed Signal" "MC" 0 \
  prompt "$FRACTION" "$NONPROMPT_MODEL" "$FRACTION_HIST" 0

# Inclusive data and the measured prompt/nonprompt MC mixture: no lxy cut.
bash "${VALIDATION_DIR}/run_DataSIGNAL_VS_MC.sh" \
  "$TREE" "$SYSTEM" "" "" 0 \
  "$INCLUSIVE_OUTPUT" "${SYSTEM} Inclusive Signal" "MC" 0 \
  "" "$FRACTION" "$NONPROMPT_MODEL" "$FRACTION_HIST" 0
