// Prompt and nonprompt lxy shapes from MC only -- needs no fit output.
//
//   root -b -q 'plot_nonPrompt_lxyShapes.C("ppRef")'
//   root -b -q 'plot_nonPrompt_lxyShapes.C("PbPb23")'

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include "aux_nonPrompt.h"

void plot_nonPrompt_lxyShapes(TString systemNAME = "ppRef")
{
	gStyle->SetOptStat(0);
	const NonPromptSetup s = NonPromptSetupFor(systemNAME);
	gSystem->mkdir(s.outDir, true);

	// the boolean selection must sit inside the parentheses, otherwise "*" binds
	// tighter than "&&" and the pThat weight is silently dropped
	const TString plotCut = Form("(%s) * pThatreweight", s.baseCut.Data());

	TFile fPromptPsi(s.promptPsi);
	TFile fPromptX(s.promptX);
	TFile fNonPromptPsi(s.nonPromptPsi);
	TFile fNonPromptX(s.nonPromptX);
	TTree* tPromptPsi    = (TTree*) fPromptPsi.Get("ntmix_PSI2S");
	TTree* tPromptX      = (TTree*) fPromptX.Get("ntmix_X3872");
	TTree* tNonPromptPsi = (TTree*) fNonPromptPsi.Get("ntmix_PSI2S");
	TTree* tNonPromptX   = (TTree*) fNonPromptX.Get("ntmix_X3872");

	TH1F hPromptPsi("hPromptPsi", Form(";%s [cm];Events", s.lxyLabel.Data()), 50, s.lxyPromptMin, s.lxyPromptMax);
	TH1F hPromptX("hPromptX", Form(";%s [cm];Events", s.lxyLabel.Data()), 50, s.lxyPromptMin, s.lxyPromptMax);
	TH1F hNonPromptPsi("hNonPromptPsi", Form(";%s [cm];Events", s.lxyLabel.Data()), 50, s.lxyNonPromptMin, s.lxyNonPromptMax);
	TH1F hNonPromptX("hNonPromptX", Form(";%s [cm];Events", s.lxyLabel.Data()), 50, s.lxyNonPromptMin, s.lxyNonPromptMax);

	tPromptPsi   ->Draw(Form("%s >> hPromptPsi", s.lxyExpr.Data()), plotCut, "goff");
	tPromptX     ->Draw(Form("%s >> hPromptX", s.lxyExpr.Data()), plotCut, "goff");
	tNonPromptPsi->Draw(Form("%s >> hNonPromptPsi", s.lxyExpr.Data()), plotCut, "goff");
	tNonPromptX  ->Draw(Form("%s >> hNonPromptX", s.lxyExpr.Data()), plotCut, "goff");

	hPromptPsi.Scale(1. / hPromptPsi.Integral());
	hPromptX.Scale(1. / hPromptX.Integral());
	hNonPromptPsi.Scale(1. / hNonPromptPsi.Integral());
	hNonPromptX.Scale(1. / hNonPromptX.Integral());

	auto lxyPassPercentage = [&](TTree* tree, double cut) {
		const WeightedCount total = WeightedEntries(tree, s.baseCut);
		const TString retainedCut = Form("(%s) && (%s > %.8f)",
		                                 s.baseCut.Data(), s.lxyExpr.Data(), cut);
		const WeightedCount retained = WeightedEntries(tree, retainedCut);
		return 100. * retained.sumW / total.sumW;
	};
	const double promptPsiRetained = 100. - lxyPassPercentage(tPromptPsi, s.lxyCut);
	const double promptXRetained = 100. - lxyPassPercentage(tPromptX, s.lxyCut);
	const double nonPromptPsiRetained = lxyPassPercentage(tNonPromptPsi, s.lxyCut);
	const double nonPromptXRetained = lxyPassPercentage(tNonPromptX, s.lxyCut);

	hPromptPsi.SetLineColor(kOrange - 2);     hPromptPsi.SetLineWidth(3);     hPromptPsi.SetFillStyle(0);
	hNonPromptPsi.SetLineColor(kOrange - 2);  hNonPromptPsi.SetLineWidth(3);  hNonPromptPsi.SetFillStyle(0);
	hPromptX.SetLineColor(kOrange - 3);       hPromptX.SetLineWidth(3);       hPromptX.SetFillStyle(0);
	hNonPromptX.SetLineColor(kOrange - 3);    hNonPromptX.SetLineWidth(3);    hNonPromptX.SetFillStyle(0);

	// One panel per production mechanism. The common lxy working point is shown
	// by the black dashed line.
	for (int isNonPrompt = 0; isNonPrompt < 2; ++isNonPrompt) {
		TH1F& hPsi = isNonPrompt ? hNonPromptPsi : hPromptPsi;
		TH1F& hX   = isNonPrompt ? hNonPromptX : hPromptX;
		const double yMin = 1.e-5;
		const double yMax = 5. * TMath::Max(hPsi.GetMaximum(), hX.GetMaximum());
		hPsi.SetMinimum(yMin);
		hPsi.SetMaximum(yMax);
		hPsi.GetXaxis()->SetTitleSize(0.045);
		hPsi.GetYaxis()->SetTitleSize(0.050);
		hPsi.GetYaxis()->SetTitleOffset(1.15);

		TCanvas c("c", "", 700, 700);
		c.SetLogy();
		c.SetLeftMargin(0.15);
		c.SetRightMargin(0.05);
		c.SetTopMargin(0.07);
		c.SetBottomMargin(0.12);
		hPsi.Draw("HIST");
		hX.Draw("HIST SAME");

		TH1F hPsiShaded(hPsi);
		TH1F hXShaded(hX);
		hPsiShaded.SetName(Form("hPsiShaded_%d", isNonPrompt));
		hXShaded.SetName(Form("hXShaded_%d", isNonPrompt));
		for (int i = 1; i <= hPsiShaded.GetNbinsX(); ++i) {
			const bool outsideRetainedRegion = isNonPrompt
			    ? hPsiShaded.GetBinCenter(i) < s.lxyCut
			    : hPsiShaded.GetBinCenter(i) >= s.lxyCut;
			if (outsideRetainedRegion) {
				hPsiShaded.SetBinContent(i, 0.);
				hXShaded.SetBinContent(i, 0.);
			}
		}
		hPsiShaded.SetFillColor(kOrange - 2);
		hPsiShaded.SetFillStyle(3002);
		hPsiShaded.SetLineWidth(0);
		hXShaded.SetFillColor(kOrange - 3);
		hXShaded.SetFillStyle(3002);
		hXShaded.SetLineWidth(0);
		hPsiShaded.Draw("HIST SAME");
		hXShaded.Draw("HIST SAME");
		hPsi.Draw("HIST SAME");
		hX.Draw("HIST SAME");

		TLine lxyCut(s.lxyCut, yMin, s.lxyCut, yMax);
		lxyCut.SetLineColor(kBlack);
		lxyCut.SetLineStyle(2);
		lxyCut.SetLineWidth(2);
		lxyCut.Draw("SAME");

		const double legendX1 = isNonPrompt ? 0.48 : 0.18;
		const double legendX2 = isNonPrompt ? 0.93 : 0.62;
		TLegend leg(legendX1, 0.61, legendX2, 0.87);
		leg.SetBorderSize(0);
		leg.SetFillStyle(0);
		leg.SetTextFont(42);
		leg.SetTextSize(0.028);
		const double psiRetained = isNonPrompt ? nonPromptPsiRetained : promptPsiRetained;
		const double xRetained = isNonPrompt ? nonPromptXRetained : promptXRetained;
		const TString psiLabel = Form("#psi(2S), retained %.2f%%", psiRetained);
		const TString xLabel = Form("X(3872), retained %.2f%%", xRetained);
		leg.SetHeader(Form("#bf{MC %s, %s}",
		                   isNonPrompt ? "non Prompt" : "Prompt", s.ptLabel.Data()),
		              "L");
		leg.AddEntry(&hPsi, psiLabel, "l");
		leg.AddEntry(&hX, xLabel, "l");
		leg.AddEntry(&lxyCut,
		             Form("%s %s %g cm", s.lxyLabel.Data(), isNonPrompt ? ">" : "<", s.lxyCut),
		             "l");
		leg.Draw();

		TLatex text;
		text.SetNDC();
		text.SetTextFont(42);
		text.SetTextSize(0.042);
		text.SetTextAlign(11);
		text.DrawLatex(0.16, 0.95, "#bf{CMS} #it{Simulation}");
		c.SaveAs(Form("%s/lxy_%s_X_PSI2S.pdf", s.outDir.Data(), isNonPrompt ? "nonPrompt" : "prompt"));
	}

	TFile out(Form("%s/lxy_shapes.root", s.outDir.Data()), "RECREATE");
	hPromptPsi.Write();
	hPromptX.Write();
	hNonPromptPsi.Write();
	hNonPromptX.Write();
	out.Close();
}
