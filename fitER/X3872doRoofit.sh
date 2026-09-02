DOANALYSISPbPb_FULL_X=1
DOANALYSISPbPb_BINNED_PT_X=0
DOANALYSISPbPb_BINNED_Y_X=0
DOANALYSISPbPb_BINNED_MULT_X=0

##
syst="ppRef"

#DATA and MC Samples
#MC="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/MC_with_score.root"
#DATA="/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pp24_v3_fid2_4v1_xgb_v1/DATA_with_score.root"
MC="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root"
DATA="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root"

## SELECTION CUTs go here
CUTs_INC="Btrk1dR < 0.5 && Btrk2dR < 0.5 && BQvalue < 0.15"
CUTs="((Bpt > 7.5  && Bpt < 12.5 && Prediction > 0.24) || (Bpt > 12.5 && Bpt < 17.5 && Prediction > 0.38) || (Bpt > 17.5 && Bpt < 22.5 && Prediction > 0.44) || (Bpt > 22.5 && Bpt < 50 && Prediction > 0.10)) && BQvalue < 0.15  "



mkdir -p ROOTfiles/

# The Function to be called:
# void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString SYSTEM = "ppRef"){


if [ $DOANALYSISPbPb_FULL_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\", \
                      1, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs_INC\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_PT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

if [ $DOANALYSISPbPb_BINNED_MULT_X  -eq 1  ]; then
root -b -q "roofitB.C++(\"ntmix_X3872\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"nChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\")"
fi

rm roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
