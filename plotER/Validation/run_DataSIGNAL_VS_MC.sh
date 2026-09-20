#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_DataSIGNAL_VS_MC.sh [tree] [system] [custom_cut] [weight_to_apply] [weight_particle]
#
# With no weight_to_apply, the nominal run writes all hWeight_<variable> histograms.
# Providing one weight variable implicitly requests a reweighted validation.
#
# Examples:
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef
#   bash run_DataSIGNAL_VS_MC.sh ntmix_PSI2S ppRef

#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' Bpt
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' Bchi2Prob
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' Bpt PSI2S


#   bash run_DataSIGNAL_VS_MC.sh ntKp ppRef "Bnorm_svpvDistance_2D > 4 && Bpt > 7.5 && Bpt <= 60"

TREE="${1:-ntKp}"
SYSTEM="${2:-ppRef}"
CUSTOM_CUT="${3:-}"
WEIGHT_TO_APPLY="${4:-}"
WEIGHT_PARTICLE="${5:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CUTs="1"
BASE="/eos/user/h/hmarques/Analysis_CODES"
SAMPLE_BASE="/eos/user/h/hmarques/RUN3_Data_MC_sharing"

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
    if [[ "$SYSTEM" == "PbPb23" ]]; then
      DATA="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_DATA.root"
      MC="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_MC_X3872.root"
      CUTs="(Bpt > 15 && Bpt < 50) && (abs(By) < 1.6) && (BQvalue < 0.15) && Btrk2dR <= 0.25 && Score > 0.85"
    else
      DATA="${SAMPLE_BASE}/X3872/ppRef24/flat_ntmix_ppRef_DATA.root"
      MC="${SAMPLE_BASE}/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root"
      CUTs="Btrk1dR < 0.5 && Btrk2dR < 0.5 && BQvalue < 0.15  && (Bpt > 7.5  && Bpt < 50)"
    fi
    ;;
  ntmix_psi2s|ntmix_PSI2S)
    TREE="ntmix_PSI2S"
    DATA_TREE="ntmix"
    PARTICLE="PSI2S"
    WEIGHT_TREE="ntmix"
    MASS_AXIS_TITLE="m_{J/#psi #pi^{-} #pi^{+}} [GeV/c^{2}]"
    if [[ "$SYSTEM" == "PbPb23" ]]; then
      DATA="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_DATA.root"
      MC="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_MC_PSI2S.root"
      CUTs="(Bpt > 15 && Bpt < 50) && (abs(By) < 1.6) && (BQvalue < 0.15) && Btrk2dR <= 0.25 && Score > 0.85"
    else
      DATA="${SAMPLE_BASE}/X3872/ppRef24/flat_ntmix_ppRef_DATA.root"
      MC="${SAMPLE_BASE}/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root"
      CUTs="Btrk1dR < 0.5 && Btrk2dR < 0.5 && BQvalue < 0.15  && (Bpt > 7.5  && Bpt < 50)"
    fi
    ;;
  ntphi)
    DATA_TREE="ntphi"
    PARTICLE="Bs"
    WEIGHT_TREE="ntphi"
    MASS_AXIS_TITLE="m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]"
    DATA="${SAMPLE_BASE}/Bmesons/ppRef/flat_ntphi_ppRef_DATA.root"
    MC="${SAMPLE_BASE}/Bmesons/ppRef/flat_ntphi_ppRef_MC.root"
    CUTs="Bnorm_svpvDistance_2D > 4 && Bpt > 7.5"
    ;;
  ntKp)
    DATA_TREE="ntKp"
    PARTICLE="Bp"
    WEIGHT_TREE="ntKp"
    MASS_AXIS_TITLE="m_{J/#psi K^{+}} [GeV/c^{2}]"
    DATA="${SAMPLE_BASE}/Bmesons/ppRef/flat_ntKp_ppRef_DATA.root"
    MC="${SAMPLE_BASE}/Bmesons/ppRef/flat_ntKp_ppRef_MC.root"
    CUTs="Bnorm_svpvDistance_2D > 4 && Bpt > 7.5 && Bpt <= 60"
    ;;
  ntKstar)
    DATA_TREE="ntKstar"
    PARTICLE="B0"
    WEIGHT_TREE="ntKstar"
    MASS_AXIS_TITLE="m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]"
    DATA="${SAMPLE_BASE}/Bmesons/ppRef/flat_ntKstar_ppRef_DATA.root"
    MC="${SAMPLE_BASE}/Bmesons/ppRef/flat_ntKstar_ppRef_MC.root"
    CUTs="Bnorm_svpvDistance_2D > 4 && Bpt > 7.5"
    ;;
esac

CUT="${CUSTOM_CUT:-$CUTs}"
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
echo "  CUT         = $CUT"
echo "  DATA        = $DATA"
echo "  MC          = $MC"
echo "  MODEL       = $MODEL"
echo "  RUN_TYPE    = $RUN_TYPE"
echo "  WEIGHT_TO_APPLY = ${WEIGHT_TO_APPLY:-<none>}"
echo "  WEIGHT_PARTICLE = $WEIGHT_PARTICLE"
echo "  WEIGHT_FILE = $WEIGHT_FILE"
echo "  OUTPUT_DIR  = Compare_${SYSTEM}/${TREE}"

root -l -b -q "DataSIGNAL_VS_MC.C(\"${DATA}\",\"${MC}\",\"${MODEL}\",\"${CUT}\",\"${TREE}\",\"${DATA_TREE}\",\"${MASS_AXIS_TITLE}\",\"${WEIGHT_FILE}\",\"${WEIGHT_TO_APPLY}\",\"${WEIGHT_PARTICLE}\",\"${SYSTEM}\")"
