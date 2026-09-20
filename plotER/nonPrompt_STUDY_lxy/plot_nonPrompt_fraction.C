// Data-driven nonprompt fraction with the B-enriched method:
//
//     f_nonprompt = N(lxy cut, data fit) / [ eff_lxy(nonprompt MC) * N(inclusive, data fit) ]
//
// The binning is read back from the fit histograms, so it always matches the fits.
//
// inclusiveOnly = true  -> only the pT-inclusive point, which needs nothing but the
//                          FULL (nominalFitModel) fits. This is the case to use while
//                          the binned fits do not converge.
// inclusiveOnly = false -> adds the binned pT and nChargedTracks fits.
//
//   root -b -q 'plot_nonPrompt_fraction.C("PbPb23", true)'
//   root -b -q 'plot_nonPrompt_fraction.C("ppRef", false)'

#include <fstream>
#include <iomanip>
#include <vector>

#include <TFile.h>
#include <TH1D.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "aux_nonPrompt.h"

void plot_nonPrompt_fraction(TString systemNAME = "ppRef", bool inclusiveOnly = false, bool drawBinnedPoints = true)
{
	gStyle->SetOptStat(0);
	const NonPromptSetup s = NonPromptSetupFor(systemNAME);
	gSystem->mkdir(s.outDir, true);

	TFile fNonPromptPsi(s.nonPromptPsi);
	TFile fNonPromptX(s.nonPromptX);
	TTree* tNonPromptPsi = (TTree*) fNonPromptPsi.Get("ntmix_PSI2S");
	TTree* tNonPromptX   = (TTree*) fNonPromptX.Get("ntmix_X3872");

	std::vector<TString> blockName;
	std::vector<TH1D*> effPsi, effX, fracPsi, fracX;

	// PbPb23 prompt Bpt fits have one inclusive [15,50] bin; other systems retain
	// their separate FULL prompt fit when the Bpt result is genuinely binned.
	const TString promptPsiInclusive = s.system == "PbPb23"
		? Form("%s/fitResults_ntmix_PSI2S_Bpt_%s.root", s.fitDir.Data(), s.system.Data())
		: Form("%s/nominalFitModel_ntmix_PSI2S_%s.root", s.fitDir.Data(), s.system.Data());
	const TString promptXInclusive = s.system == "PbPb23"
		? Form("%s/fitResults_ntmix_X3872_Bpt_%s.root", s.fitDir.Data(), s.system.Data())
		: Form("%s/nominalFitModel_ntmix_X3872_%s.root", s.fitDir.Data(), s.system.Data());
	TH1D* hPsiInclFitFull = LoadHpt(promptPsiInclusive, "hPsi2S_inclusiveFitYield_inclusive");
	TH1D* hXInclFitFull   = LoadHpt(promptXInclusive, "hX3872_inclusiveFitYield_inclusive");
	TH1D* hPsiBenrFitFull = LoadHpt(Form("%s/nominalFitModel_ntmix_PSI2S_%s_nonPrompt.root", s.fitDirNonPrompt.Data(), s.system.Data()), "hPsi2S_BenrichedFitYield_inclusive");
	TH1D* hXBenrFitFull   = LoadHpt(Form("%s/nominalFitModel_ntmix_X3872_%s_nonPrompt.root", s.fitDirNonPrompt.Data(), s.system.Data()), "hX3872_BenrichedFitYield_inclusive");

	printf("\n[%s] psi(2S) lxy efficiency, inclusive (cut %g cm)\n", s.system.Data(), s.lxyCutPsi);
	TH1D* hPsiEffFull = LxyPassFraction(tNonPromptPsi, "hPsi2S_nonprompt_lxy_fraction_inclusive", hPsiInclFitFull, "Bpt", s, s.lxyCutPsi);
	printf("\n[%s] X(3872) lxy efficiency, inclusive (cut %g cm)\n", s.system.Data(), s.lxyCutX);
	TH1D* hXEffFull   = LxyPassFraction(tNonPromptX, "hX3872_nonprompt_lxy_fraction_inclusive", hXInclFitFull, "Bpt", s, s.lxyCutX);

	TH1D* hPsiFracFull = EmptyLike(hPsiInclFitFull, "hPsi2S_dataDriven_nonprompt_fraction_inclusive", ";p_{T} [GeV/c];f_{nonprompt}");
	TH1D* hXFracFull   = EmptyLike(hXInclFitFull, "hX3872_dataDriven_nonprompt_fraction_inclusive", ";p_{T} [GeV/c];f_{nonprompt}");
	FillNonPromptFraction(*hPsiFracFull, hPsiInclFitFull, hPsiBenrFitFull, hPsiEffFull);
	FillNonPromptFraction(*hXFracFull, hXInclFitFull, hXBenrFitFull, hXEffFull);

	blockName.push_back("inclusive");
	effPsi.push_back(hPsiEffFull);   effX.push_back(hXEffFull);
	fracPsi.push_back(hPsiFracFull); fracX.push_back(hXFracFull);

	// ------------------------------------------------------- binned pT and multiplicity
	TH1D *hPsiFrac = nullptr, *hXFrac = nullptr, *hPsiFracMult = nullptr, *hXFracMult = nullptr;
	if (!inclusiveOnly) {
		TH1D* hPsiInclFit     = LoadHpt(Form("%s/fitResults_ntmix_PSI2S_Bpt_%s.root", s.fitDir.Data(), s.system.Data()), "hPsi2S_inclusiveFitYield_Bpt");
		TH1D* hXInclFit       = LoadHpt(Form("%s/fitResults_ntmix_X3872_Bpt_%s.root", s.fitDir.Data(), s.system.Data()), "hX3872_inclusiveFitYield_Bpt");
		TH1D* hPsiBenrFit     = LoadHpt(Form("%s/fitResults_ntmix_PSI2S_Bpt_%s_nonPrompt.root", s.fitDirNonPrompt.Data(), s.system.Data()), "hPsi2S_BenrichedFitYield_Bpt");
		TH1D* hXBenrFit       = LoadHpt(Form("%s/fitResults_ntmix_X3872_Bpt_%s_nonPrompt.root", s.fitDirNonPrompt.Data(), s.system.Data()), "hX3872_BenrichedFitYield_Bpt");
		TH1D* hPsiInclFitMult = LoadHpt(Form("%s/fitResults_ntmix_PSI2S_nChargedTracks_%s.root", s.fitDir.Data(), s.system.Data()), "hPsi2S_inclusiveFitYield_nChargedTracks");
		TH1D* hXInclFitMult   = LoadHpt(Form("%s/fitResults_ntmix_X3872_nChargedTracks_%s.root", s.fitDir.Data(), s.system.Data()), "hX3872_inclusiveFitYield_nChargedTracks");
		TH1D* hPsiBenrFitMult = LoadHpt(Form("%s/fitResults_ntmix_PSI2S_nChargedTracks_%s_nonPrompt.root", s.fitDirNonPrompt.Data(), s.system.Data()), "hPsi2S_BenrichedFitYield_nChargedTracks");
		TH1D* hXBenrFitMult   = LoadHpt(Form("%s/fitResults_ntmix_X3872_nChargedTracks_%s_nonPrompt.root", s.fitDirNonPrompt.Data(), s.system.Data()), "hX3872_BenrichedFitYield_nChargedTracks");

		printf("\n[%s] psi(2S) lxy efficiency vs Bpt\n", s.system.Data());
		TH1D* hPsiEff     = LxyPassFraction(tNonPromptPsi, "hPsi2S_nonprompt_lxy_fraction_Bpt", hPsiInclFit, "Bpt", s, s.lxyCutPsi);
		printf("\n[%s] X(3872) lxy efficiency vs Bpt\n", s.system.Data());
		TH1D* hXEff       = LxyPassFraction(tNonPromptX, "hX3872_nonprompt_lxy_fraction_Bpt", hXInclFit, "Bpt", s, s.lxyCutX);
		printf("\n[%s] psi(2S) lxy efficiency vs nChargedTracks\n", s.system.Data());
		TH1D* hPsiEffMult = LxyPassFraction(tNonPromptPsi, "hPsi2S_nonprompt_lxy_fraction_nChargedTracks", hPsiInclFitMult, "nChargedTracks", s, s.lxyCutPsi);
		printf("\n[%s] X(3872) lxy efficiency vs nChargedTracks\n", s.system.Data());
		TH1D* hXEffMult   = LxyPassFraction(tNonPromptX, "hX3872_nonprompt_lxy_fraction_nChargedTracks", hXInclFitMult, "nChargedTracks", s, s.lxyCutX);

		hPsiFrac     = EmptyLike(hPsiInclFit, "hPsi2S_dataDriven_nonprompt_fraction_Bpt", ";p_{T} [GeV/c];f_{nonprompt}");
		hXFrac       = EmptyLike(hXInclFit, "hX3872_dataDriven_nonprompt_fraction_Bpt", ";p_{T} [GeV/c];f_{nonprompt}");
		hPsiFracMult = EmptyLike(hPsiInclFitMult, "hPsi2S_dataDriven_nonprompt_fraction_nChargedTracks", ";N_{trk};f_{nonprompt}");
		hXFracMult   = EmptyLike(hXInclFitMult, "hX3872_dataDriven_nonprompt_fraction_nChargedTracks", ";N_{trk};f_{nonprompt}");
		FillNonPromptFraction(*hPsiFrac, hPsiInclFit, hPsiBenrFit, hPsiEff);
		FillNonPromptFraction(*hXFrac, hXInclFit, hXBenrFit, hXEff);
		FillNonPromptFraction(*hPsiFracMult, hPsiInclFitMult, hPsiBenrFitMult, hPsiEffMult);
		FillNonPromptFraction(*hXFracMult, hXInclFitMult, hXBenrFitMult, hXEffMult);

		blockName.push_back("Bpt");
		effPsi.push_back(hPsiEff);   effX.push_back(hXEff);
		fracPsi.push_back(hPsiFrac); fracX.push_back(hXFrac);
		blockName.push_back("nChargedTracks");
		effPsi.push_back(hPsiEffMult);   effX.push_back(hXEffMult);
		fracPsi.push_back(hPsiFracMult); fracX.push_back(hXFracMult);
	}

	// ----------------------------------------------------------------------- bookkeeping
	const TString tag = inclusiveOnly ? "_inclusive" : "";
	std::ofstream out(Form("%s/nonprompt_fraction%s.txt", s.outDir.Data(), tag.Data()));
	out << std::fixed << std::setprecision(6);
	out << "System: " << s.system.Data() << "\n";
	out << "Base cut: " << s.baseCut.Data() << "\n";
	out << "lxy expression: " << s.lxyExpr.Data() << "\n";
	out << "lxy cut psi(2S) [cm]: " << s.lxyCutPsi << "\n";
	out << "lxy cut X(3872) [cm]: " << s.lxyCutX << "\n";
	out << "MC weight: pThatreweight\n";

	printf("\n[%s] data-driven nonprompt fractions\n", s.system.Data());
	for (size_t b = 0; b < blockName.size(); ++b) {
		out << "\n# " << blockName[b].Data() << "\n";
		out << "low high psi2S_effLxy psi2S_fnp psi2S_unc X3872_effLxy X3872_fnp X3872_unc\n";
		printf("  -- %s --\n", blockName[b].Data());
		for (int i = 1; i <= fracPsi[b]->GetNbinsX(); ++i) {
			const double low = fracPsi[b]->GetXaxis()->GetBinLowEdge(i);
			const double high = fracPsi[b]->GetXaxis()->GetBinUpEdge(i);
			out << low << " " << high << " "
			    << effPsi[b]->GetBinContent(i) << " "
			    << fracPsi[b]->GetBinContent(i) << " " << fracPsi[b]->GetBinError(i) << " "
			    << effX[b]->GetBinContent(i) << " "
			    << fracX[b]->GetBinContent(i) << " " << fracX[b]->GetBinError(i) << "\n";
			printf("  [%g, %g]: psi2S = %.4f +/- %.4f, X3872 = %.4f +/- %.4f\n", low, high,
			       fracPsi[b]->GetBinContent(i), fracPsi[b]->GetBinError(i),
			       fracX[b]->GetBinContent(i), fracX[b]->GetBinError(i));
		}
	}
	out.close();

	TFile fOut(Form("%s/nonprompt_fraction%s.root", s.outDir.Data(), tag.Data()), "RECREATE");
	for (size_t b = 0; b < blockName.size(); ++b) {
		effPsi[b]->Write();
		effX[b]->Write();
		fracPsi[b]->Write();
		fracX[b]->Write();
	}
	TObjString baseCutObj(s.baseCut);
	TObjString lxyExprObj(s.lxyExpr);
	TParameter<double> lxyCutPsiPar("lxyCutPsi2S_cm", s.lxyCutPsi);
	TParameter<double> lxyCutXPar("lxyCutX3872_cm", s.lxyCutX);
	baseCutObj.Write("baseCut");
	lxyExprObj.Write("lxyExpression");
	lxyCutPsiPar.Write();
	lxyCutXPar.Write();
	fOut.Close();

	// ---------------------------------------------------------------------------- plots
	if (inclusiveOnly) {
		DrawFractionCanvas(*hPsiFracFull, *hXFracFull, *hPsiFracFull, *hXFracFull,
		                   "p_{T} [GeV/c]", Form("%s/nonprompt_fraction_X_PSI2S_inclusive.pdf", s.outDir.Data()),
		                   s, false);
	} else {
		DrawFractionCanvas(*hPsiFrac, *hXFrac, *hPsiFracFull, *hXFracFull,
		                   "p_{T} [GeV/c]", Form("%s/nonprompt_fraction_X_PSI2S.pdf", s.outDir.Data()),
		                   s, drawBinnedPoints);
		DrawFractionCanvas(*hPsiFracMult, *hXFracMult, *hPsiFracFull, *hXFracFull,
		                   "N_{trk}", Form("%s/nonprompt_fraction_X_PSI2S_nChargedTracks.pdf", s.outDir.Data()),
		                   s, true);
	}
}
