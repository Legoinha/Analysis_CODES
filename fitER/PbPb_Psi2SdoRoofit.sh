DOANALYSISPbPb_FULL_X=1
DOANALYSISPbPb_BINNED_PT_X=1
DOANALYSISPbPb_BINNED_Y_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="PbPb23"

#Data and MC Samples
MC_X="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_MC_PSI2S.root"
Data_X="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_DATA.root"
#Data and MC Samples

## SELECTION CUTs go here 
CUTs="(Bpt > 15 && Bpt < 50) && (abs(By) < 1.6) && (BQvalue < 0.15) && Btrk2dR <= 0.25 && Score > 0.85"

#CUTs="1"  #"((Bpt > 5 && Bpt < 7.5) && abs(By) > 1.4) ||  (Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4)"


mkdir -p "ROOTfiles/$syst" "results/$syst"

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString SYSTEM = "ppRef", std::vector<double> VAR_BINS = {}){
#

if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_PSI2S\", \
                      1, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{15,50})"
fi

if [ $DOANALYSISPbPb_BINNED_PT_X  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_PSI2S\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{15,50})"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_PSI2S\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"nChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{0,15,30,50,100})"
fi

#if [ $DOANALYSISPbPb_BINNED_Y_X  -eq 1  ]; then
#root -b -q 'roofitB.C('\"ntmix\"','0','\"$Data_X\"','\"$MC_X\"','\"By\"'   ,'\"$CUTs\"','\"$OutputFile_X_BINNED_Y\"'   ,'\"results/X/By\"'  , '\" \"', '\"$syst\"')'
#fi

#if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
#root -b -q 'roofitB.C('\"ntmix\"','0','\"$Data_X\"','\"$MC_X\"','\"nMult\"','\"$CUTs\"','\"$OutputFile_X_BINNED_MULT\"','\"results/X/nMult\"', '\" \"', '\"$syst\"')'
#fi
