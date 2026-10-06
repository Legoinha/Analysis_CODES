#!/usr/bin/env bash
set -euo pipefail

# Combined (simultaneous) fit of the PbPb years of one tree, roofitCombined.C:
#
#   bash doCombinedRoofit.sh <TREE>        TREE: ntmix_X3872 ntmix_PSI2S ntKp ntphi ntKstar
#
# Inputs of every year: the scored samples of selectionER,
#   scored_samples/<case>/<year>/flat_<case>_<year>_<METHOD>_<FOM>_scored_{DATA,<MC kind>}.root,
# whose metadata give the pre-cut and the score branch, and the score threshold: the Punzi optimum of the
# inclusive bin, optimalCUT/<case>/<year>/<METHOD>_<FOM>/punzi_summary.root (selectionER/optimalCUT_X_punzi.C).
# X(3872) and psi(2S) both take the threshold of the ntmix case, so they are selected with the same cut.
# Products (README.md): every year gets the files of a single-year fit, under the same names, in
# fitER/ROOTfiles/<YEAR>/ (replacing the ones there), which validation, Acc x Eff and resultER read.

TREE="${1:?Missing tree: ntmix_X3872 ntmix_PSI2S ntKp ntphi ntKstar}"

DO_FULL=1             # inclusive fit: nominalFitModel_<TREE>_<YEAR>.root (validation)
DO_BINNED_PT=1        # pT-binned fit: fitResults_/mcFitResults_<TREE>_Bpt_<YEAR>.root (Acc x Eff)
SYST=1                # model variations of the binned fit: systematicFILES/ (resultER)
METHOD="nn"
FOM="auc"
SELECTION="/eos/user/h/hmarques/Analysis_CODES/selectionER"

case "$TREE" in
  ntmix_X3872|ntmix_PSI2S)
    CASE="ntmix"
    MC_KIND="MC_${TREE#ntmix_}"
    YEARS=(PbPb23 PbPb24)
    FULL_BINS="15,50"
    PT_BINS="15,50"
    ;;
  ntKp|ntphi|ntKstar)
    CASE="$TREE"
    MC_KIND="MC"
    YEARS=(PbPb23 PbPb24)
    FULL_BINS="7.5,50"
    PT_BINS="7.5,10,15,20,30,50"
    ;;
esac

cd "$(dirname "${BASH_SOURCE[0]}")"
trap 'rm -f roofitCombined_C.so roofitCombined_C.d roofitCombined_C_ACLiC_dict_rdict.pcm' EXIT

YEARS_V=""
DATA_V=""
MC_V=""
THRESHOLDS_V=""
for YEAR in "${YEARS[@]}"; do
  SCORED="${SELECTION}/scored_samples/${CASE}/${YEAR}/flat_${CASE}_${YEAR}_${METHOD}_${FOM}_scored"
  PUNZI="${SELECTION}/optimalCUT/${CASE}/${YEAR}/${METHOD}_${FOM}/punzi_summary.root"
  THRESHOLD=$(python3 -c "import ROOT; f = ROOT.TFile.Open('${PUNZI}'); print(next(e.optimalThreshold for e in f.Get('punziResults') if e.inclusive))" | tail -1)
  echo "${TREE} ${YEAR}: score threshold ${THRESHOLD} (${PUNZI})"
  YEARS_V+="${YEARS_V:+,}\"${YEAR}\""
  DATA_V+="${DATA_V:+,}\"${SCORED}_DATA.root\""
  MC_V+="${MC_V:+,}\"${SCORED}_${MC_KIND}.root\""
  THRESHOLDS_V+="${THRESHOLDS_V:+,}${THRESHOLD}"
done

# void roofitCombined(TString TREE, int FULL, std::vector<TString> YEARS, std::vector<TString> INPUTDATA,
#                     std::vector<TString> INPUTMC, std::vector<double> THRESHOLDS, TString VAR,
#                     std::vector<double> VAR_BINS, int SYST)

if [ $DO_FULL -eq 1 ]; then
root -l -b -q "roofitCombined.C+(\"${TREE}\", \
                      1, \
                      std::vector<TString>{${YEARS_V}}, \
                      std::vector<TString>{${DATA_V}}, \
                      std::vector<TString>{${MC_V}}, \
                      std::vector<double>{${THRESHOLDS_V}}, \
                      \"Bpt\", \
                      std::vector<double>{${FULL_BINS}}, \
                      0)"
fi

if [ $DO_BINNED_PT -eq 1 ]; then
root -l -b -q "roofitCombined.C+(\"${TREE}\", \
                      0, \
                      std::vector<TString>{${YEARS_V}}, \
                      std::vector<TString>{${DATA_V}}, \
                      std::vector<TString>{${MC_V}}, \
                      std::vector<double>{${THRESHOLDS_V}}, \
                      \"Bpt\", \
                      std::vector<double>{${PT_BINS}}, \
                      ${SYST})"
fi
