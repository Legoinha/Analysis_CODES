DOANALYSISPbPb_FULL_X=1
DOANALYSISPbPb_BINNED_PT_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="PbPb23_nonPrompt"

#Data and MC Samples
MC_X="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_MC_X3872_nonPrompt.root"
Data_X="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored/flat_ntmix_PbPb23_DATA.root"
#Data and MC Samples

## SELECTION CUTs go here
## Same PbPb23 baseline as PbPb_X3872doRoofit.sh + the non-prompt (displaced) requirement
CUTs="(Bpt > 15 && Bpt < 50) && (abs(By) < 1.6) && (BQvalue < 0.15) && Btrk2dR <= 0.25 && Score > 0.85 && BLxy*(Bmass/Bpt)>0.01"


mkdir -p "ROOTfiles/$syst" "results/$syst"



if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\", \
                      1, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{15,50})"
fi

if [ $DOANALYSISPbPb_BINNED_PT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{15,20,50})"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$Data_X\", \
                      \"$MC_X\", \
                      \"nChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{0,3800,8000})"
fi

rm -f roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
