DOANALYSISPbPb_FULL_X=1
DOANALYSISPbPb_BINNED_PT_X=0
DOANALYSISPbPb_BINNED_Y_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="ppRef"

#Data and MC Samples
#MC_X="/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_ppRef_scored_MC_X3872.root"
#Data_X="/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_ppRef_scored_DATA.root"
MC_X="/eos/user/h/hmarques/Analysis_CODES/fitER/MC_Signal_pp_AANN.root"
Data_X="/eos/user/h/hmarques/Analysis_CODES/fitER/DATA_Signal_pp_AANN.root"
#Data and MC Samples

## SELECTION CUTs go here 
#CUTs="BQvalue<0.150 && Btrk1dR < 0.55 && Btrk2dR < 0.55 "
#CUTs="Btrk1dR < 0.55 && Btrk2dR < 0.55 && Bpt < 50"  ## RUN1 && Bpt > 10 && Bpt < 50 && abs(By)<1.2
#CUTs="xgb_score > 0.65 && abs(By) < 1.2 && Bpt > 10 && BQvalue < 0.1"  ##
CUTs="BQvalue < 0.1"  #"((Bpt > 5 && Bpt < 7.5) && abs(By) > 1.4) ||  (Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4)"


mkdir -p ROOTfiles/

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString SYSTEM = "ppRef"){
#

if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\", \
                      1, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"nSelectedChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

#if [ $DOANALYSISPbPb_BINNED_Y_X  -eq 1  ]; then
#root -b -q 'roofitB.C('\"ntmix\"','0','\"$Data_X\"','\"$MC_X\"','\"By\"'   ,'\"$CUTs\"','\"$OutputFile_X_BINNED_Y\"'   ,'\"results/X/By\"'  , '\" \"', '\"$syst\"')'
#fi

#if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
#root -b -q 'roofitB.C('\"ntmix\"','0','\"$Data_X\"','\"$MC_X\"','\"nMult\"','\"$CUTs\"','\"$OutputFile_X_BINNED_MULT\"','\"results/X/nMult\"', '\" \"', '\"$syst\"')'
#fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
