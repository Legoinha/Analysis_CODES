DOANALYSISPbPb_FULL_Bp=0
DOANALYSISPbPb_BINNED_PT_Bp=1
DOANALYSISPbPb_BINNED_Y_Bp=0
DOANALYSISPbPb_BINNED_MULT_Bp=0

#Data and MC Samples
Data_Bp="/eos/user/h/hmarques/Analysis_CODES/flatER/flat_ntKp_ppRef_DATA.root"
MC_Bp="/eos/user/h/hmarques/Analysis_CODES/flatER/flat_ntKp_ppRef_MC.root"
#Data and MC Samples

## NEW CUTS ? here 
CUTs="Bnorm_svpvDistance_2D > 4"

##
syst="ppRef"

mkdir -p ROOTfiles/

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString OUTPLOTF = "", TString ExtraMCsample = "", TString SYSTEM = "ppRef"){
#


if [ $DOANALYSISPbPb_FULL_Bp  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntKp\", \
                      1, \
                      \"$Data_Bp\", \
                      \"$MC_Bp\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"results/Bp/Bpt\", \
                      \" \", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_Bp  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntKp\",\
                      0, \
                      \"$Data_Bp\", \
                      \"$MC_Bp\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"results/Bp/Bpt\", \
                      \" \", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
