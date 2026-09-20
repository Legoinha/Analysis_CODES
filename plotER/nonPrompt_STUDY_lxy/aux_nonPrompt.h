#pragma once

#include <vector>

#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>


// Everything that changes between ppRef and PbPb23 lives here and only here.
// The lxy cuts MUST stay identical to the ones in the fitER nonPrompt .sh files,
// otherwise the MC efficiency and the B-enriched yields do not describe the same sample.
struct NonPromptSetup {
	TString system;
	TString promptPsi;
	TString promptX;
	TString nonPromptPsi;
	TString nonPromptX;
	TString baseCut;          // the complete rectangular selection, pT window included
	TString lxyExpr;          // displacement variable used to enrich the B component
	TString lxyLabel;         // how that variable is written on the axes
	double  lxyCutPsi;        // [cm]  -> fitER/<...>Psi2SdoRoofit_nonPrompt.sh
	double  lxyCutX;          // [cm]  -> fitER/<...>X3872doRoofit_nonPrompt.sh
	TString ptLabel;          // caption only, never part of a selection
	int     lxyNonPromptBins; // axis of the nonprompt lxy panel; the prompt panel is
	double  lxyNonPromptMin;  // -0.2 -- 0.2 in both systems and stays hardcoded
	double  lxyNonPromptMax;
	TString collisionText;
	TString fitDir;           // inclusive (prompt+nonprompt) fits
	TString fitDirNonPrompt;  // B-enriched (lxy cut) fits
	TString outDir;
};

inline NonPromptSetup NonPromptSetupFor(TString systemNAME)
{
	NonPromptSetup s;
	s.system          = systemNAME;
	s.fitDir          = Form("../../fitER/ROOTfiles/%s", systemNAME.Data());
	s.fitDirNonPrompt = Form("../../fitER/ROOTfiles/%s_nonPrompt", systemNAME.Data());
	s.outDir          = systemNAME;
	s.lxyExpr         = "BLxy*(Bmass/Bpt)";
	s.lxyLabel        = "l_{xy}";

	if (systemNAME == "PbPb23") {
		// the PbPb23 selection needs the XGBoost Score, so the scored samples are used
		const TString dir = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/output/selected/X_pb23_v19_fid13_9v9_rw0_xgb_v1/root_scored";
		s.promptPsi     = dir + "/flat_ntmix_PbPb23_MC_PSI2S.root";
		s.promptX       = dir + "/flat_ntmix_PbPb23_MC_X3872.root";
		s.nonPromptPsi  = dir + "/flat_ntmix_PbPb23_MC_PSI2S_nonPrompt.root";
		s.nonPromptX    = dir + "/flat_ntmix_PbPb23_MC_X3872_nonPrompt.root";
		s.baseCut       = "(Bpt > 15 && Bpt < 50) && (abs(By) < 1.6) && (BQvalue < 0.15) && Btrk2dR <= 0.25 && Score > 0.85";
		s.lxyCutPsi     = 0.01;
		s.lxyCutX       = 0.01;
		s.ptLabel       = "15 < p_{T} < 50 GeV/c";
		s.lxyNonPromptBins = 50;
		s.lxyNonPromptMin  = -0.02;
		s.lxyNonPromptMax  = 0.5;
		s.collisionText = "PbPb #sqrt{s_{NN}}=5.36 TeV";
	} else {
		const TString dir = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24";
		s.promptPsi     = dir + "/flat_ntmix_ppRef_MC_PSI2S.root";
		s.promptX       = dir + "/flat_ntmix_ppRef_MC_X3872.root";
		s.nonPromptPsi  = dir + "/flat_ntmix_ppRef_MC_PSI2S_nonPrompt.root";
		s.nonPromptX    = dir + "/flat_ntmix_ppRef_MC_X3872_nonPrompt.root";
		s.baseCut       = "BQvalue < 0.15 && Btrk1dR < .5 && Btrk2dR < .5 && (Bpt > 7.5 && Bpt < 50)";
		s.lxyCutPsi     = 0.03;
		s.lxyCutX       = 0.03;
		s.ptLabel       = "7.5 < p_{T} < 50 GeV/c";
		s.lxyNonPromptBins = 50;
		s.lxyNonPromptMin  = -0.2;
		s.lxyNonPromptMax  = .5;
		s.collisionText = "pp #sqrt{s}=5.36 TeV";
	}
	return s;
}


// pThatreweight is present in every X3872/psi(2S) MC sample of both systems.
inline double WeightedEntries(TTree* tree, const TString& selection)
{
	tree->SetEstimate(tree->GetEntries() + 1);
	const Long64_t n = tree->Draw("pThatreweight", selection, "goff");
	const double* values = tree->GetV1();
	double sum = 0.;
	for (Long64_t i = 0; i < n; ++i) sum += values[i];
	return sum;
}

// The binning is always inherited from the fit output, so the MC efficiency and
// the fitted yields can never end up on different bin edges.
inline TH1D* LoadHpt(TString path, TString cloneName)
{
	TFile* f = TFile::Open(path, "READ");
	TH1D* out = (TH1D*) f->Get("hPt")->Clone(cloneName);
	out->SetDirectory(nullptr);
	f->Close();
	return out;
}

// A fresh histogram on ref's binning. Built from the bin edges rather than cloned:
// the fit's hPt carries alphanumeric bin labels (the mean of the variable in each
// bin, set in roofitB.C), and a Clone drags them along and turns every axis drawn
// here into a labelled one instead of a numeric one.
inline TH1D* EmptyLike(const TH1D* ref, TString name, TString title)
{
	const TAxis* ax = ref->GetXaxis();
	const int nbins = ax->GetNbins();
	std::vector<double> edges(nbins + 1);
	for (int i = 1; i <= nbins; ++i) edges[i - 1] = ax->GetBinLowEdge(i);
	edges[nbins] = ax->GetBinUpEdge(nbins);
	TH1D* h = new TH1D(name, title, nbins, edges.data());
	h->SetDirectory(nullptr);
	h->SetStats(0);
	return h;
}

// Fraction of nonprompt MC surviving the lxy cut, bin by bin of ref.
inline TH1D* LxyPassFraction(TTree* tree, TString name, const TH1D* ref, TString varExpr,
                             const NonPromptSetup& s, double lxyCut)
{
	TH1D* h = EmptyLike(ref, name, Form(";%s;N(%s > %g cm) / N", varExpr.Data(), s.lxyLabel.Data(), lxyCut));
	for (int i = 1; i <= h->GetNbinsX(); ++i) {
		const double lo = h->GetXaxis()->GetBinLowEdge(i);
		const double hi = h->GetXaxis()->GetBinUpEdge(i);
		const TString binCut = (i == h->GetNbinsX())
			? Form("%s >= %.8f && %s <= %.8f", varExpr.Data(), lo, varExpr.Data(), hi)
			: Form("%s >= %.8f && %s < %.8f",  varExpr.Data(), lo, varExpr.Data(), hi);
		const TString den = Form("(%s) && (%s)", s.baseCut.Data(), binCut.Data());
		const TString num = Form("(%s) && (%s > %.6f)", den.Data(), s.lxyExpr.Data(), lxyCut);
		const double d = WeightedEntries(tree, den);
		const double n = WeightedEntries(tree, num);
		const double f = (d > 0.) ? n / d : 0.;
		h->SetBinContent(i, f);
		h->SetBinError(i, (d > 0.) ? TMath::Sqrt(f * (1. - f) / d) : 0.);
		printf("  %s [%g, %g]: %.4f (%g / %g)\n", varExpr.Data(), lo, hi, f, n, d);
	}
	return h;
}

// f_nonprompt = N_Benriched / ( eff_lxy^MC * N_inclusive ), fully propagated.
inline void FillNonPromptFraction(TH1D& hOut, TH1D* hIncl, TH1D* hBenr, TH1D* hBenrNonPromptMC)
{
	for (int i = 1; i <= hOut.GetNbinsX(); ++i) {
		const double incl = hIncl->GetBinContent(i);
		const double inclErr = hIncl->GetBinError(i);
		const double benr = hBenr->GetBinContent(i);
		const double benrErr = hBenr->GetBinError(i);
		const double fmc = hBenrNonPromptMC->GetBinContent(i);
		const double fmcErr = hBenrNonPromptMC->GetBinError(i);
		if (incl <= 0. || fmc <= 0.) continue;
		const double nonprompt = benr / (fmc * incl);
		const double dB = 1. / (fmc * incl);
		const double dI = -benr / (fmc * incl * incl);
		const double dF = -benr / (fmc * fmc * incl);
		hOut.SetBinContent(i, nonprompt);
		hOut.SetBinError(i, TMath::Sqrt(dB * dB * benrErr * benrErr +
		                                dI * dI * inclErr * inclErr +
		                                dF * dF * fmcErr * fmcErr));
	}
}

inline TLine* DrawInclusiveBand(TH1D& hInclusive, Color_t color, double xLow, double xHigh)
{
	const double y = hInclusive.GetBinContent(1);
	const double e = hInclusive.GetBinError(1);
	TBox* band = new TBox(xLow, y - e, xHigh, y + e);
	band->SetFillColorAlpha(color, 0.14);
	band->SetLineColor(color);
	band->SetLineStyle(2);
	band->SetLineWidth(1);
	band->Draw("same");
	TLine* line = new TLine(xLow, y, xHigh, y);
	line->SetLineColor(color);
	line->SetLineStyle(2);
	line->SetLineWidth(3);
	line->Draw("same");
	return line;
}

inline void DrawFractionCanvas(TH1D& hPsi, TH1D& hX, TH1D& hPsiIncl, TH1D& hXIncl,
                               TString xTitle, TString outPdf,
                               const NonPromptSetup& s, bool drawBinnedPoints)
{
	hPsi.SetLineColor(kOrange - 2);   hPsi.SetMarkerColor(kOrange - 2);
	hPsi.SetMarkerStyle(20);          hPsi.SetMarkerSize(1.0);
	hPsi.SetLineWidth(3);             hPsi.SetStats(0);
	hX.SetLineColor(kOrange - 3);     hX.SetMarkerColor(kOrange - 3);
	hX.SetMarkerStyle(21);            hX.SetMarkerSize(1.0);
	hX.SetLineWidth(3);               hX.SetStats(0);

	double yMax = 0.;
	if (drawBinnedPoints) {
		for (int i = 1; i <= hPsi.GetNbinsX(); ++i) {
			yMax = TMath::Max(yMax, hPsi.GetBinContent(i) + hPsi.GetBinError(i));
			yMax = TMath::Max(yMax, hX.GetBinContent(i) + hX.GetBinError(i));
		}
	}
	yMax = TMath::Max(yMax, hPsiIncl.GetBinContent(1) + hPsiIncl.GetBinError(1));
	yMax = TMath::Max(yMax, hXIncl.GetBinContent(1) + hXIncl.GetBinError(1));
	yMax = TMath::Max(1.25, 1.35 * yMax);

	const double xLow = hPsi.GetXaxis()->GetXmin();
	const double xHigh = hPsi.GetXaxis()->GetXmax();
	TH1D* frame = EmptyLike(&hPsi, "hFractionFrame", Form(";%s;f_{nonprompt}", xTitle.Data()));
	frame->SetMinimum(0.);
	frame->SetMaximum(yMax);
	frame->GetXaxis()->SetTitleSize(0.045);
	frame->GetYaxis()->SetTitleSize(0.050);
	frame->GetYaxis()->SetTitleOffset(1.15);

	TCanvas c("cFraction", "", 700, 700);
	c.SetLeftMargin(0.15);
	c.SetRightMargin(0.05);
	c.SetTopMargin(0.07);
	c.SetBottomMargin(0.12);
	frame->Draw("AXIS");
	TLine* psiInclusiveLine = DrawInclusiveBand(hPsiIncl, kOrange - 2, xLow, xHigh);
	TLine* xInclusiveLine = DrawInclusiveBand(hXIncl, kOrange - 3, xLow, xHigh);
	if (drawBinnedPoints) {
		hPsi.Draw("E1 SAME");
		hX.Draw("E1 SAME");
	}
	frame->Draw("AXIS SAME");

	TLegend leg(0.58, 0.68, 0.93, 0.82);
	leg.SetBorderSize(0);
	leg.SetFillStyle(0);
	leg.SetTextFont(42);
	leg.SetTextSize(0.032);
	if (drawBinnedPoints) {
		leg.AddEntry(&hPsi, "#psi(2S)", "lep");
		leg.AddEntry(&hX, "X(3872)", "lep");
		leg.AddEntry((TObject*)0, "Dashed bands: inclusive", "");
	} else {
		leg.AddEntry(psiInclusiveLine, "#psi(2S) inclusive", "l");
		leg.AddEntry(xInclusiveLine, "X(3872) inclusive", "l");
	}
	leg.Draw();

	TLatex text;
	text.SetNDC();
	text.SetTextFont(42);
	text.SetTextSize(0.042);
	text.SetTextAlign(11);
	text.DrawLatex(0.16, 0.95, "#bf{CMS} #it{Preliminary}");
	text.SetTextAlign(31);
	text.SetTextSize(0.035);
	text.DrawLatex(0.95, 0.95, s.collisionText);
	text.SetTextSize(0.040);
	text.DrawLatex(0.93, 0.87, "#bf{B-enriched method}");
	text.SetTextSize(0.035);
	text.DrawLatex(0.93, 0.63, s.ptLabel);
	c.SaveAs(outPdf);
}
