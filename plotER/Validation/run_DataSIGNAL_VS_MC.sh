#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_DataSIGNAL_VS_MC.sh <tree> <system> [weight_to_apply] [weight_particle] [show_background]
#
# The fit model is fitER/ROOTfiles/<system>/nominalFitModel_<tree>_<system>.root: the single-year fit (roofitB.C)
# or, for a PbPb year, the per-year file of the combined fit (combinedANA), whichever ran last.
# Input files, the fit selection, and the fitted-variable range come from the fit model metadata.
# With no weight_to_apply, the nominal run writes all hWeight_<variable> histograms.
# Providing one weight variable implicitly requests a reweighted validation.
#
# Examples:
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef
#   bash run_DataSIGNAL_VS_MC.sh ntmix_PSI2S ppRef
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef Bpt
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef Bpt PSI2S
#   bash run_DataSIGNAL_VS_MC.sh ntKp ppRef '' '' 0
#   bash run_DataSIGNAL_VS_MC.sh ntKp PbPb23
#   bash run_DataSIGNAL_VS_MC.sh ntKp PbPb24

TREE="${1:?Missing required tree argument. See usage above.}"
SYSTEM="${2:?Missing required system argument. See usage above.}"
WEIGHT_TO_APPLY="${3:-}"
WEIGHT_PARTICLE="${4:-}"
SHOW_BACKGROUND="${5:-1}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BASE="/eos/user/h/hmarques/Analysis_CODES"

cleanup_aclic() {
  rm -f \
    DataSIGNAL_VS_MC_C.so \
    DataSIGNAL_VS_MC_C.d \
    DataSIGNAL_VS_MC_C_ACLiC_dict_rdict.pcm \
    DataSIGNAL_VS_MC_C_ACLiC_dict.cxx \
    DataSIGNAL_VS_MC_C_ACLiC_linkdef.h \
    DataSIGNAL_VS_MC_C_ACLiC_map
}
trap cleanup_aclic EXIT

case "$TREE" in
  ntmix|ntmix_X3872)
    TREE="ntmix_X3872"
    DATA_TREE="ntmix"
    PARTICLE="X3872"
    WEIGHT_TREE="ntmix"
    MASS_AXIS_TITLE="m_{J/#psi #pi^{-} #pi^{+}} [GeV/c^{2}]"
    ;;
  ntmix_psi2s|ntmix_PSI2S)
    TREE="ntmix_PSI2S"
    DATA_TREE="ntmix"
    PARTICLE="PSI2S"
    WEIGHT_TREE="ntmix"
    MASS_AXIS_TITLE="m_{J/#psi #pi^{-} #pi^{+}} [GeV/c^{2}]"
    ;;
  ntphi)
    DATA_TREE="ntphi"
    PARTICLE="Bs"
    WEIGHT_TREE="ntphi"
    MASS_AXIS_TITLE="m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]"
    ;;
  ntKp)
    DATA_TREE="ntKp"
    PARTICLE="Bp"
    WEIGHT_TREE="ntKp"
    MASS_AXIS_TITLE="m_{J/#psi K^{+}} [GeV/c^{2}]"
    ;;
  ntKstar)
    DATA_TREE="ntKstar"
    PARTICLE="B0"
    WEIGHT_TREE="ntKstar"
    MASS_AXIS_TITLE="m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]"
    ;;
esac

MODEL="${BASE}/fitER/ROOTfiles/${SYSTEM}/nominalFitModel_${TREE}_${SYSTEM}.root"

if [[ -n "$WEIGHT_TO_APPLY" ]]; then
  RUN_TYPE="reweighted"
  WEIGHT_PARTICLE="${WEIGHT_PARTICLE:-$PARTICLE}"
else
  RUN_TYPE="nominal (write all comparison weights)"
  WEIGHT_PARTICLE="$PARTICLE"
fi
WEIGHT_FILE="WEIGHTS/${WEIGHT_TREE}_${SYSTEM}_${WEIGHT_PARTICLE}_weight.root"

echo "Running DataSIGNAL_VS_MC.C with:"
echo "  TREE        = $TREE"
echo "  SYSTEM      = $SYSTEM"
echo "  MODEL       = $MODEL"
echo "  RUN_TYPE    = $RUN_TYPE"
echo "  WEIGHT_TO_APPLY = ${WEIGHT_TO_APPLY:-<none>}"
echo "  WEIGHT_PARTICLE = $WEIGHT_PARTICLE"
echo "  WEIGHT_FILE = $WEIGHT_FILE"
echo "  SHOW_BACKGROUND = $SHOW_BACKGROUND"
echo "  OUTPUT_DIR  = Compare_${SYSTEM}/${TREE}"

root -l -b -q "DataSIGNAL_VS_MC.C(\"${MODEL}\",\"${TREE}\",\"${DATA_TREE}\",\"${MASS_AXIS_TITLE}\",\"${WEIGHT_FILE}\",\"${WEIGHT_TO_APPLY}\",\"${WEIGHT_PARTICLE}\",\"${SYSTEM}\",${SHOW_BACKGROUND})"
