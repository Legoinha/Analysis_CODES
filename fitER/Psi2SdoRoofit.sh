DOANALYSISPbPb_FULL_PSI=0
DOANALYSISPbPb_BINNED_PT_PSI=1
DOANALYSISPbPb_BINNED_Y_PSI=0
DOANALYSISPbPb_BINNED_MULT_PSI=0

##
syst="ppRef"

#Data and MC Samples
MC_PSI="/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_MC_wScore_psi2s.root"
Data_PSI="/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_DATA_wScore.root"
#Data and MC Samples

## SELECTION CUTs go here
CUTs="xgb_score > 0.55 && BQvalue<0.10"

mkdir -p ROOTfiles/

if [ $DOANALYSISPbPb_FULL_PSI -eq 1 ]; then
root -b -q "roofitB.C++(\"ntmix_psi2s\", \
                      1, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_PSI -eq 1 ]; then
root -b -q "roofitB.C(\"ntmix_psi2s\",\
                      0, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_PSI -eq 1 ]; then
root -b -q "roofitB.C(\"ntmix_psi2s\",\
                      0, \
                      \"$Data_PSI\", \
                      \"$MC_PSI\", \
                      \"nSelectedChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
