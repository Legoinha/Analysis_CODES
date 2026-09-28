DOANALYSISPbPb_FULL_PSI=1
DOANALYSISPbPb_BINNED_PT_PSI=0
DOANALYSISPbPb_BINNED_MULT_PSI=0

##
syst="PbPb23_nonPrompt"

#Data and MC Samples
MC_PSI="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_MC_PSI2S_nonPrompt.root"
Data_PSI="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_DATA.root"
#Data and MC Samples

## SELECTION CUTs go here
## Same PbPb23 baseline as PbPb_Psi2SdoRoofit.sh + the non-prompt (displaced) requirement
CUTs="(Bpt > 15 && Bpt < 50) && (abs(By) < 1.6) && (BQvalue < 0.15) && Btrk2dR <= 0.25 && Score > 0.85 && BLxy*(Bmass/Bpt)>0.01"

mkdir -p "ROOTfiles/$syst" "results/$syst"



if [ $DOANALYSISPbPb_FULL_PSI  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_PSI2S\", \
                      1, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{15,50})"
fi

if [ $DOANALYSISPbPb_BINNED_PT_PSI  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_PSI2S\",\
                      0, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{15,20,50})"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_PSI  -eq 1  ]; then
root -b -q "roofitB.C(\"ntmix_PSI2S\",\
                      0, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"nChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{0,3800,8000})"
fi

