DOANALYSISPbPb_FULL_X=0
DOANALYSISPbPb_BINNED_PT_X=1
DOANALYSISPbPb_BINNED_Y_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="ppRef"



#Data and MC Samples
MC_X="/eos/user/h/hmarques/X3872/DATA_sharing/dataMC_2024/flat_ntmix_ppRef_MC.root"
Data_X="/eos/user/h/hmarques/X3872/DATA_sharing/dataMC_2024/flat_ntmix_ppRef_DATA.root"    ###Rectandular cuts
#Data_X="/eos/user/h/hmarques/Analysis_CODES/DATA_Signal_X_PSI2S_TrM.root"            ###From Prince ML
#Data and MC Samples

## SELECTION CUTs go here 
CUTs="BQvalue<0.15 && Btrk1dR < 0.55 && Btrk2dR < 0.55"
#CUTs="1"


mkdir -p ROOTfiles/
OutputFile_X_FULL="ROOTfiles/yields_X_full"
OutputFile_X_BINNED_Y="ROOTfiles/yields_X_binned_y"
OutputFile_X_BINNED_PT="ROOTfiles/yields_X_binned_pt"
OutputFile_X_BINNED_MULT="ROOTfiles/yields_X_binned_Mult"

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString outputfile = "", TString outplotf = "", TString ExtraMCsample = "", TString SYSTEM = "ppRef"){
#

if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix\", \
                      1, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$OutputFile_X_FULL\", \
                      \"results/X/Bpt\", \
                      \" \", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_X  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$OutputFile_X_BINNED_PT\", \
                      \"results/X/Bpt\", \
                      \" \", \
                      \"$syst\")"
fi

#if [ $DOANALYSISPbPb_BINNED_Y_X  -eq 1  ]; then
#root -b -q 'roofitB.C('\"ntmix\"','0','\"$Data_X\"','\"$MC_X\"','\"By\"'   ,'\"$CUTs\"','\"$OutputFile_X_BINNED_Y\"'   ,'\"results/X/By\"'  , '\" \"', '\"$syst\"')'
#fi

#if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
#root -b -q 'roofitB.C('\"ntmix\"','0','\"$Data_X\"','\"$MC_X\"','\"nMult\"','\"$CUTs\"','\"$OutputFile_X_BINNED_MULT\"','\"results/X/nMult\"', '\" \"', '\"$syst\"')'
#fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so