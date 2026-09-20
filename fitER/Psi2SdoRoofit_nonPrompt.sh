syst="ppRef_nonPrompt"


MC="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S_nonPrompt.root"
DATA="/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root"

## SELECTION CUTs go here
CUTs_INC="Btrk1dR < 0.5 && Btrk2dR < 0.5 && BQvalue < 0.15 && BLxy*(Bmass/Bpt)>0.03 "
CUTs=" BQvalue < 0.15 && Btrk1dR < .5 && Btrk2dR < .5 && BLxy*(Bmass/Bpt)>0.03 "



mkdir -p ROOTfiles/

root -b -q "roofitB.C(\"ntmix_PSI2S\", \
                      1, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs_INC\", \
                      \"$syst\", \
                      std::vector<double>{7.5,50})"

root -b -q "roofitB.C(\"ntmix_PSI2S\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"Bpt\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{7.5,12.5,17.5,22.5,50})"

root -b -q "roofitB.C(\"ntmix_PSI2S\",\
                      0, \
                      \"$DATA\", \
                      \"$MC\", \
                      \"nChargedTracks\", \
                      \"$CUTs\", \
                      \"$syst\", \
                      std::vector<double>{0,15,30,50,100})"

rm -f roofitB_C.d roofitB_C_ACLiC_dict_rdict.pcm roofitB_C.so
