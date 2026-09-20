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

	TH1F hPromptPsi("hPromptPsi", Form(";%s [cm];Events", s.lxyLabel.Data()), 50, -0.02, 0.02);
	TH1F hPromptX("hPromptX", Form(";%s [cm];Events", s.lxyLabel.Data()), 50, -0.02, 0.02);
	TH1F hNonPromptPsi("hNonPromptPsi", Form(";%s [cm];Events", s.lxyLabel.Data()), s.lxyNonPromptBins, s.lxyNonPromptMin, s.lxyNonPromptMax);
	TH1F hNonPromptX("hNonPromptX", Form(";%s [cm];Events", s.lxyLabel.Data()), s.lxyNonPromptBins, s.lxyNonPromptMin, s.lxyNonPromptMax);

	tPromptPsi   ->Draw(Form("%s >> hPromptPsi", s.lxyExpr.Data()), plotCut, "goff");
	tPromptX     ->Draw(Form("%s >> hPromptX", s.lxyExpr.Data()), plotCut, "goff");
	tNonPromptPsi->Draw(Form("%s >> hNonPromptPsi", s.lxyExpr.Data()), plotCut, "goff");
	tNonPromptX  ->Draw(Form("%s >> hNonPromptX", s.lxyExpr.Data()), plotCut, "goff");

	hPromptPsi.Scale(1. / hPromptPsi.Integral());
	hPromptX.Scale(1. / hPromptX.Integral());
	hNonPromptPsi.Scale(1. / hNonPromptPsi.Integral());
	hNonPromptX.Scale(1. / hNonPromptX.Integral());

	hPromptPsi.SetLineColor(kOrange - 2);     hPromptPsi.SetLineWidth(3);     hPromptPsi.SetFillStyle(0);
	hNonPromptPsi.SetLineColor(kOrange - 2);  hNonPromptPsi.SetLineWidth(3);  hNonPromptPsi.SetFillStyle(0);
	hPromptX.SetLineColor(kOrange - 3);       hPromptX.SetLineWidth(3);       hPromptX.SetFillStyle(0);
	hNonPromptX.SetLineColor(kOrange - 3);    hNonPromptX.SetLineWidth(3);    hNonPromptX.SetFillStyle(0);

	// one panel per production mechanism; the two lxy working points are drawn in
	// the colour of the state they are applied to
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

		TLine lPsi(s.lxyCutPsi, yMin, s.lxyCutPsi, yMax);
		lPsi.SetLineColor(kOrange - 2);
		lPsi.SetLineStyle(2);
		lPsi.SetLineWidth(2);
		lPsi.Draw("SAME");
		TLine lX(s.lxyCutX, yMin, s.lxyCutX, yMax);
		lX.SetLineColor(kOrange - 3);
		lX.SetLineStyle(2);
		lX.SetLineWidth(2);
		lX.Draw("SAME");

		TLegend leg(0.62, 0.68, 0.93, 0.82);
		leg.SetBorderSize(0);
		leg.SetFillStyle(0);
		leg.SetTextFont(42);
		leg.SetTextSize(0.035);
		leg.AddEntry(&hPsi, Form("#psi(2S), cut %g cm", s.lxyCutPsi), "l");
		leg.AddEntry(&hX, Form("X(3872), cut %g cm", s.lxyCutX), "l");
		leg.Draw();

		TLatex text;
		text.SetNDC();
		text.SetTextFont(42);
		text.SetTextSize(0.042);
		text.SetTextAlign(11);
		text.DrawLatex(0.16, 0.95, "#bf{CMS} #it{Simulation}");
		text.SetTextAlign(31);
		text.SetTextSize(0.035);
		text.DrawLatex(0.95, 0.95, s.collisionText);
		text.SetTextSize(0.040);
		text.DrawLatex(0.93, 0.87, isNonPrompt ? "#bf{Nonprompt}" : "#bf{Prompt}");
		text.SetTextSize(0.035);
		text.DrawLatex(0.93, 0.65, s.ptLabel);
		c.SaveAs(Form("%s/lxy_%s_X_PSI2S.pdf", s.outDir.Data(), isNonPrompt ? "nonPrompt" : "prompt"));
	}

	TFile out(Form("%s/lxy_shapes.root", s.outDir.Data()), "RECREATE");
	hPromptPsi.Write();
	hPromptX.Write();
	hNonPromptPsi.Write();
	hNonPromptX.Write();
	out.Close();
}
