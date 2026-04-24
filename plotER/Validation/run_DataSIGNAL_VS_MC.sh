#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_DataSIGNAL_VS_MC.sh [tree] [custom_cut] [data] [mc] [model]
# Examples:
#   bash run_DataSIGNAL_VS_MC.sh ntphi
#   bash run_DataSIGNAL_VS_MC.sh ntmix
#   bash run_DataSIGNAL_VS_MC.sh ntmix_psi2s
#   bash run_DataSIGNAL_VS_MC.sh ntphi "Bnorm_svpvDistance_2D > 4"

TREE="${1:-ntphi}"
CUT_DEFAULT="1"

# Auto defaults per tree (override with args 3/4)
case "$TREE" in
  ntmix)
    DATA_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_DATA_wScore.root"
    MC_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_MC_wScore_X3872.root"
    CUT_DEFAULT="xgb_score > 0.55 && BQvalue<0.15"
    ;;
  ntmix_psi2s)
    DATA_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_DATA_wScore.root"
    MC_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_MC_wScore_psi2s.root"
    CUT_DEFAULT="xgb_score > 0.55 && BQvalue<0.15"
    ;;
  ntphi)
    DATA_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntphi_ppRef_DATA.root"
    MC_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntphi_ppRef_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKp)
    DATA_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKp_ppRef_DATA.root"
    MC_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKp_ppRef_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKstar)
    DATA_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_DATA.root"
    MC_DEFAULT="/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_MC.root"
    CUT_DEFAULT="Bnorm_svpvDistance_2D > 4"
    ;;
  *)
    echo "Unknown tree: $TREE"
    echo "Use one of: ntmix, ntmix_psi2s, ntphi, ntKp, ntKstar"
    exit 1
    ;;
esac

CUT="${2:-$CUT_DEFAULT}"

DATA="${3:-$DATA_DEFAULT}"
MC="${4:-$MC_DEFAULT}"
MODEL="${5:-/eos/user/h/hmarques/Analysis_CODES/fitER/ROOTfiles/nominalFitModel_${TREE}_ppRef.root}"

echo "Running DataSIGNAL_VS_MC.C with:"
echo "  TREE  = $TREE"
echo "  CUT   = $CUT"
echo "  DATA  = $DATA"
echo "  MC    = $MC"
echo "  MODEL = $MODEL"

root -l -b -q "DataSIGNAL_VS_MC.C(\"${DATA}\",\"${MC}\",\"${MODEL}\",\"${CUT}\",\"${TREE}\")"
