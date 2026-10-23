# Run from inside plotER/nonPrompt_STUDY_lxy -- the fit paths are relative to it.

DOSHAPES=1                # prompt/nonprompt lxy MC distributions, reads the inclusive-fit selection metadata
DOFRACTION_INCLUSIVE=1    # pT-inclusive point only, needs just the FULL (nominalFitModel) fits
DOFRACTION_BINNED=1       # adds binned pT + nChargedTracks, needs those fits to have converged

SYSTEMS="ppRef" ## just include other systems without commas : "ppRef PbPb23"

for syst in $SYSTEMS; do

	if [ $DOSHAPES -eq 1 ]; then
	root -l -b -q "plot_nonPrompt_lxyShapes.C(\"$syst\")"
	fi

	if [ $DOFRACTION_INCLUSIVE -eq 1 ]; then
	root -l -b -q "plot_nonPrompt_fraction.C(\"$syst\", true)"
	fi

	if [ $DOFRACTION_BINNED -eq 1 ]; then
	root -l -b -q "plot_nonPrompt_fraction.C(\"$syst\", false)"
	fi

done
