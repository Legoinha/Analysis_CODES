DOANALYSISPbPb_FULL_X=1
DOANALYSISPbPb_BINNED_PT_X=0
DOANALYSISPbPb_BINNED_Y_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="PbPb23"

#Data and MC Samples
MC_X="./../../RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC.root"
Data_X="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/pbpb23_24_comb.root"  
#Data and MC Samples

## SELECTION CUTs go here 
#CUTs="BQvalue<0.150 && Btrk1dR < 0.55 && Btrk2dR < 0.55 "
#CUTs="BQvalue<0.300 && Btrk1dR < 0.55 && Btrk2dR < 0.55 && Bpt > 10 && Bpt < 50 && abs(By)<1.2"  ## RUN1
#CUTs="xgb_score > 0.55 && BQvalue<0.13"
CUTs="abs(By) < 1.2 && Bpt > 10 && BQvalue < 0.2 && CentBin > 20 && Btrk1dR < .25 && Btrk2dR < .25 && BtrkPtimb > 0.15"
#CUTs="1"  #"((Bpt > 5 && Bpt < 7.5) && abs(By) > 1.4) ||  (Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4)"


mkdir -p ROOTfiles/

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString OUTPLOTF = "", TString ExtraMCsample = "", TString SYSTEM = "ppRef"){
#

if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix\", \
                      1, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
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
                      \"results/X/Bpt\", \
                      \" \", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"nSelectedChargedTracks\", \
                      \"$CUTs\", \
                      \"results/X/nMult\", \
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