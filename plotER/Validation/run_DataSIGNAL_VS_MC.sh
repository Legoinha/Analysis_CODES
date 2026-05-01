#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_DataSIGNAL_VS_MC.sh [tree] [system] [custom_cut] [data] [mc] [model]
# Examples:
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef
#   bash run_DataSIGNAL_VS_MC.sh ntmix_PSI2S ppRef
#   bash run_DataSIGNAL_VS_MC.sh ntphi ppRef "Bnorm_svpvDistance_2D > 4"

TREE="${1:-ntphi}"
SYSTEM="${2:-ppRef}"
CUT_DEFAULT="1"
BASE="/eos/user/h/hmarques/Analysis_CODES"

# Auto defaults per tree/system (override with args 4/5/6)
case "$TREE" in
  ntmix|ntmix_X3872)
    TREE="ntmix_X3872"
    DATA_DEFAULT="${BASE}/selectionER/scored_samples/flat_ntmix_${SYSTEM}_scored_DATA.root"
    MC_DEFAULT="${BASE}/selectionER/scored_samples/flat_ntmix_${SYSTEM}_scored_MC_X3872.root"
    CUT_DEFAULT="xgb_score > 0.7 && abs(By) < 1.2 && Bpt > 10 && BQvalue < 0.1"
    ;;
  ntmix_psi2s|ntmix_PSI2S)
    TREE="ntmix_PSI2S"
    DATA_DEFAULT="${BASE}/selectionER/scored_samples/flat_ntmix_${SYSTEM}_scored_DATA.root"
    MC_DEFAULT="${BASE}/selectionER/scored_samples/flat_ntmix_${SYSTEM}_scored_MC_PSI2S.root"
    CUT_DEFAULT="1"
    ;;
  ntphi)
    DATA_DEFAULT="${BASE}/flatER/Bmeson/flat_ntphi_${SYSTEM}_DATA.root"
    MC_DEFAULT="${BASE}/flatER/Bmeson/flat_ntphi_${SYSTEM}_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKp)
    DATA_DEFAULT="${BASE}/flatER/Bmeson/flat_ntKp_${SYSTEM}_DATA.root"
    MC_DEFAULT="${BASE}/flatER/Bmeson/flat_ntKp_${SYSTEM}_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKstar)
    DATA_DEFAULT="${BASE}/flatER/Bmeson/flat_ntKstar_${SYSTEM}_DATA.root"
    MC_DEFAULT="${BASE}/flatER/Bmeson/flat_ntKstar_${SYSTEM}_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  *)
    echo "Unknown tree: $TREE"
    echo "Use one of: ntmix_X3872, ntmix_PSI2S, ntphi, ntKp, ntKstar"
    exit 1
    ;;
esac

CUT="${3:-$CUT_DEFAULT}"

DATA="${4:-$DATA_DEFAULT}"
MC="${5:-$MC_DEFAULT}"
MODEL="${6:-${BASE}/fitER/ROOTfiles/${SYSTEM}/nominalFitModel_${TREE}_${SYSTEM}.root}"
if [[ ! -f "$MODEL" && -f "${BASE}/fitER/ROOTfiles/nominalFitModel_${TREE}_${SYSTEM}.root" ]]; then
  MODEL="${BASE}/fitER/ROOTfiles/nominalFitModel_${TREE}_${SYSTEM}.root"
fi

echo "Running DataSIGNAL_VS_MC.C with:"
echo "  TREE  = $TREE"
echo "  SYSTEM= $SYSTEM"
echo "  CUT   = $CUT"
echo "  DATA  = $DATA"
echo "  MC    = $MC"
echo "  MODEL = $MODEL"

root -l -b -q "DataSIGNAL_VS_MC.C(\"${DATA}\",\"${MC}\",\"${MODEL}\",\"${CUT}\",\"${TREE}\")"
