DOANALYSISPbPb_FULL_PSI=1
DOANALYSISPbPb_BINNED_PT_PSI=0
DOANALYSISPbPb_BINNED_Y_PSI=0
DOANALYSISPbPb_BINNED_MULT_PSI=0

##
syst="ppRef"

#Data and MC Samples
#MC_PSI="/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_ppRef_scored_MC_PSI2S.root"
#Data_PSI="/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_ppRef_scored_DATA.root"
MC_PSI="/eos/user/h/hmarques/Analysis_CODES/fitER/MC_Signal_pp_AANN.root"
Data_PSI="/eos/user/h/hmarques/Analysis_CODES/fitER/DATA_Signal_pp_AANN.root"
#Data and MC Samples

## SELECTION CUTs go here
#CUTs="xgb_score > 0.65 && abs(By) < 1.2 && Bpt > 10 && BQvalue < 0.1"
CUTs="BQvalue < 0.1"

mkdir -p ROOTfiles/

if [ $DOANALYSISPbPb_FULL_PSI -eq 1 ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\", \
                      1, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_PSI -eq 1 ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\",\
                      0, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_PSI -eq 1 ]; then
root -b -q "roofitB.C++(\"ntmix_PSI2S\",\
                      0, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"nSelectedChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
