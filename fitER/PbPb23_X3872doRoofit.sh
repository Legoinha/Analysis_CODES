DOANALYSISPbPb_FULL_X=1
DOANALYSISPbPb_BINNED_PT_X=0
DOANALYSISPbPb_BINNED_Y_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="PbPb23"

#Data and MC Samples
MC_X="/eos/user/h/hmarques/Analysis_CODES/X_pb23_v36_fid3_7v7_rw0_xgb_v1/flat_ntmix_PbPb23_MC_X3872.root"
Data_X="/eos/user/h/hmarques/Analysis_CODES/X_pb23_v36_fid3_7v7_rw0_xgb_v1/flat_ntmix_PbPb23_DATA.root"
#Data and MC Samples

## SELECTION CUTs go here 
CUTs="(Bpt > 15 && Bpt < 50) && BQvalue < 0.15 && (abs(By) < 1.6) && (Btrk2dR < 0.35) && (Btrk1dR < 0.35) && Prediction > 0.86"  

#CUTs="1"  #"((Bpt > 5 && Bpt < 7.5) && abs(By) > 1.4) ||  (Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4)"


mkdir -p "ROOTfiles/$syst" "results/$syst"

#The Function to be called:
#
#void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString SYSTEM = "ppRef", std::vector<double> VAR_BINS = {}){
#

if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_X3872\", \
                      1, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{15,50})"
fi

if [ $DOANALYSISPbPb_BINNED_PT_X  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_X3872\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{15,50})"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_X3872\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"nChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{0,15,30,50,100})"
fi
