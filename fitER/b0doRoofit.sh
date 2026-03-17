DOANALYSISPbPb_FULL_B0=0
DOANALYSISPbPb_BINNED_PT_B0=1
DOANALYSISPbPb_BINNED_Y_B0=0
DOANALYSISPbPb_BINNED_MULT_B0=0

#Data and MC Samples
Data_B0="/eos/user/h/hmarques/Analysis_CODES/flatER/flat_ntKstar_ppRef_DATA.root"
MC_B0="/eos/user/h/hmarques/Analysis_CODES/flatER/flat_ntKstar_ppRef_MC.root"
#Data and MC Samples

## CUTS (SELECTION ?) here 
CUTs="Bnorm_svpvDistance_2D > 4"

##
syst="ppRef"

mkdir -p ROOTfiles/

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString OUTPLOTF = "", TString ExtraMCsample = "", TString SYSTEM = "ppRef"){
#

if [ $DOANALYSISPbPb_FULL_B0  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntKstar\", \
                      1, \
                      \"$Data_B0\", \
                      \"$MC_B0\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"results/B0/Bpt\", \
                      \" \", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_B0  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntKstar\",\
                      0, \
                      \"$Data_B0\", \
                      \"$MC_B0\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"results/B0/Bpt\", \
                      \" \", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
