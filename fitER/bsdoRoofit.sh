DOANALYSISPbPb_FULL_BS=0
DOANALYSISPbPb_BINNED_PT_BS=0
DOANALYSISPbPb_BINNED_Y_BS=0
DOANALYSISPbPb_BINNED_MULT_BS=0

#Data and MC Samples
MC_Bs="/lstore/cms/lekai/Bmeson/MC/ppRef/MC_Bs.root"
Data_Bs="/lstore/cms/lekai/Bmeson/data/data_Bs.root"
#Data and MC Samples

## NEW CUTS ? here 
CUTs="1"

##
syst = "ppRef"

mkdir -p ROOTfiles/
OutputFile_BS_FULL="ROOTfiles/yields_Bs_full"
OutputFile_BS_BINNED_Y="ROOTfiles/yields_Bs_binned_y"
OutputFile_BS_BINNED_PT="ROOTfiles/yields_Bs_binned_pt"
OutputFile_BS_BINNED_MULT="ROOTfiles/yields_Bs_binned_Mult"

#The Function to be called:
#
#void roofitB(TString tree = "ntphi", int full = 0, TString inputdata = "", TString inputmc = "", TString varExp = "", TString cut = "", TString outputfile = "", TString outplotf = "", TString jpsiFile = ""){
#

if [ $DOANALYSISPbPb_FULL_BS  -eq 1  ]; then
root -b  -q 'roofitB.C+('\"ntphi\"','1','\"$Data_Bs\"','\"$MC_Bs\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BS_FULL\"','\"results/Bs/Bpt\"')'
fi

if [ $DOANALYSISPbPb_BINNED_PT_BS  -eq 1  ]; then
root -b  -q 'roofitB.C+('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"Bpt\"','\"$CUTPbPb\"','\"$OutputFile_BS_BINNED_PT\"','\"results/Bs/Bpt\"')'
fi

if [ $DOANALYSISPbPb_BINNED_Y_BS  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"By\"','\"$CUTPbPb\"','\"$OutputFile_BS_BINNED_Y\"','\"results/Bs/By\"')'
fi

if [ $DOANALYSISPbPb_BINNED_MULT_BS  -eq 1  ]; then
root -b  -q 'roofitB.C('\"ntphi\"','0','\"$Data_Bs\"','\"$MC_Bs\"','\"nMult\"','\"$CUTPbPb\"','\"$OutputFile_BS_BINNED_MULT\"','\"results/Bs/nMult\"')'
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so