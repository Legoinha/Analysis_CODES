#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash run_DataSIGNAL_VS_MC.sh [tree] [system] [custom_cut] [mode] [reweight_variable] [weight_particle]
#
# mode:
#   nominal   Run the nominal comparison and write sPlot + ML-discrepancy weights.
#   reweight  Re-run the comparison with the existing ML-discrepancy weight.
#
# Examples:
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef
#   bash run_DataSIGNAL_VS_MC.sh ntmix_PSI2S ppRef

#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' reweight
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' reweight Bchi2Prob
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' reweight Btrk1PtErr,Bchi2Prob
#   bash run_DataSIGNAL_VS_MC.sh ntmix_X3872 ppRef '' reweight Prediction PSI2S


#   bash run_DataSIGNAL_VS_MC.sh ntphi ppRef "Bnorm_svpvDistance_2D > 4"

TREE="${1:-ntphi}"
SYSTEM="${2:-ppRef}"
CUSTOM_CUT="${3:-}"
MODE="${4:-nominal}"
REWEIGHT_VARIABLE="${5:-Bpt}"
WEIGHT_PARTICLE="${6:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CUTs="1"
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
    DATA="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root"
    MC="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root"
    CUTs="Btrk1dR < 0.5 && Btrk2dR < 0.5 && BQvalue < 0.15"
    ;;
  ntmix_psi2s|ntmix_PSI2S)
    TREE="ntmix_PSI2S"
    DATA_TREE="ntmix"
    PARTICLE="PSI2S"
    WEIGHT_TREE="ntmix"
    MASS_AXIS_TITLE="m_{J/#psi #pi^{-} #pi^{+}} [GeV/c^{2}]"
    DATA="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root"
    MC="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root"
    CUTs="Btrk1dR < 0.5 && Btrk2dR < 0.5 && BQvalue < 0.15"
    ;;
  ntphi)
    DATA_TREE="ntphi"
    PARTICLE="Bs"
    WEIGHT_TREE="ntphi"
    MASS_AXIS_TITLE="m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]"
    DATA="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/Data_2024ppRef_Bs.root"
    MC="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/MC_2024ppRef_Bs.root"
    CUTs="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKp)
    DATA_TREE="ntKp"
    PARTICLE="Bp"
    WEIGHT_TREE="ntKp"
    MASS_AXIS_TITLE="m_{J/#psi K^{+}} [GeV/c^{2}]"
    DATA="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/Data_2024ppRef_Bu.root"
    MC="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/MC_2024ppRef_Bu.root"
    #DATA="./../../../RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKp_ppRef_DATA.root"
    #MC="./../../../RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKp_ppRef_MC.root"
    CUTs="Bnorm_svpvDistance_2D > 4"
    ;;
  ntKstar)
    DATA_TREE="ntKstar"
    PARTICLE="B0"
    WEIGHT_TREE="ntKstar"
    MASS_AXIS_TITLE="m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]"
    DATA="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/Data_2024ppRef_B0.root"
    MC="/eos/user/c/ctorresc/BmesonsHIN/PreXGBFiles/MC_2024ppRef_B0.root"
    CUTs="Bnorm_svpvDistance_2D > 4"
    ;;
esac

CUT="${CUSTOM_CUT:-$CUTs}"
WEIGHT_PARTICLE="${WEIGHT_PARTICLE:-$PARTICLE}"
MODEL="${BASE}/fitER/ROOTfiles/${SYSTEM}/nominalFitModel_${TREE}_${SYSTEM}.root"

MODE_LC="${MODE,,}"
case "$MODE_LC" in
  nominal|raw|0|false|no|"")
    REWEIGHT_MC=0
    ;;
  reweight|reweighted|rw|1|true|yes)
    REWEIGHT_MC=1
    ;;
esac

WEIGHT_FILE="WEIGHTS/${WEIGHT_TREE}_${SYSTEM}_${WEIGHT_PARTICLE}_weight.root"

echo "Running DataSIGNAL_VS_MC.C with:"
echo "  TREE        = $TREE"
echo "  SYSTEM      = $SYSTEM"
echo "  CUT         = $CUT"
echo "  DATA        = $DATA"
echo "  MC          = $MC"
echo "  MODEL       = $MODEL"
echo "  MODE        = $MODE_LC"
echo "  REW_VAR     = $REWEIGHT_VARIABLE"
echo "  WEIGHT_PARTICLE = $WEIGHT_PARTICLE"
echo "  WEIGHT_FILE = $WEIGHT_FILE"

root -l -b -q "DataSIGNAL_VS_MC.C(\"${DATA}\",\"${MC}\",\"${MODEL}\",\"${CUT}\",\"${TREE}\",\"${DATA_TREE}\",\"${MASS_AXIS_TITLE}\",${REWEIGHT_MC},\"${WEIGHT_FILE}\",\"${REWEIGHT_VARIABLE}\",\"${WEIGHT_PARTICLE}\")"
