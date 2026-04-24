#pragma once

#include "aux/uti.h"
#include "RooWorkspace.h"
#include "TString.h"

using namespace RooFit;
using namespace std;

Int_t _count=0;

RooFitResult *fit(TString variation, TString pdf, TString tree, TCanvas* c, TCanvas* cMC, RooDataSet* ds, RooDataSet* dsMC, RooRealVar* mass, float binmin, float binmax, RooWorkspace& w, TString which_var, int NBIN){
	double init_mean = Bs_MASS;
	if (tree == "ntKp") init_mean = Bu_MASS;
	else if (tree == "ntKstar") init_mean = Bd_MASS;
	else if (tree == "ntmix") init_mean = X3872_MASS;
	else if (tree == "ntmix_psi2s") init_mean = PSI2S_MASS;

	RooRealVar mean(Form("mean%d_%s", _count, pdf.Data()), "", init_mean, init_mean - 0.01, init_mean + 0.01);
	RooRealVar sigma1(Form("sigma1%d_%s", _count, pdf.Data()), "", 0.01, 0.001, 0.1);
	RooRealVar sigma2(Form("sigma2%d_%s", _count, pdf.Data()), "", 0.005, 0.001, 0.1);
	RooRealVar sigma3(Form("sigma3%d_%s", _count, pdf.Data()), "", 0.01, 0.001, 0.03);
	RooRealVar sigma4cb(Form("sigma4cb%d_%s", _count, pdf.Data()), "", 0.005, 0.001, 0.05);
	RooRealVar alpha(Form("alpha%d_%s", _count, pdf.Data()), "", 4., 0, 15);
	RooRealVar n(Form("n_%d_%s", _count, pdf.Data()), "", 10, -100, 200);
	RooRealVar* scale = new RooRealVar("scale", "scale", 1, 0.2, 2);
	RooProduct scaled_sigma1(Form("scaled_sigma1%d_%s", _count, pdf.Data()), "scaled_sigma1", RooArgList(*scale, sigma1));
	RooProduct scaled_sigma2(Form("scaled_sigma2%d_%s", _count, pdf.Data()), "scaled_sigma2", RooArgList(*scale, sigma2));
	RooProduct scaled_sigma3(Form("scaled_sigma3%d_%s", _count, pdf.Data()), "scaled_sigma3", RooArgList(*scale, sigma3));
	RooProduct scaled_sigmacb(Form("scaled_sigmacb%d_%s", _count, pdf.Data()), "scaled_sigmacb", RooArgList(*scale, sigma4cb));
	RooGaussian sig1(Form("sig1%d_%s", _count, pdf.Data()), "", *mass, mean, scaled_sigma1);
	RooGaussian sig2(Form("sig2%d_%s", _count, pdf.Data()), "", *mass, mean, scaled_sigma2);
	RooGaussian sig3(Form("sig3%d_%s", _count, pdf.Data()), "", *mass, mean, scaled_sigma3);
	RooCBShape CB(Form("CB%d_%s", _count, pdf.Data()), "", *mass, mean, scaled_sigmacb, alpha, n);
	RooRealVar sig1frac(Form("sig1frac%d_%s", _count, pdf.Data()), "", 0.5, 0.01, 1);
	RooRealVar sig2frac(Form("sig2frac%d_%s", _count, pdf.Data()), "", 0.5, 0.01, 1);

	RooAddPdf* sig = nullptr;
	if ((variation == "" && pdf == "") || variation == "background" || (variation == "signal" && pdf == "fixed")) sig = new RooAddPdf(Form("sig_doubleG%d_%s", _count, pdf.Data()), "", RooArgList(sig1, sig2), sig1frac);
	if (variation == "signal" && pdf == "1gauss") sig = new RooAddPdf(Form("sig_Gaussian%d_%s", _count, pdf.Data()), "", RooArgList(sig1), RooArgList(), true);
	if (variation == "signal" && pdf == "3gauss") sig = new RooAddPdf(Form("sig_tripleG%d_%s", _count, pdf.Data()), "", RooArgList(sig1, sig2, sig3), RooArgList(sig1frac, sig2frac), true);
	if (variation == "signal" && pdf == "gauss_cb") sig = new RooAddPdf(Form("sig_gaussCB%d_%s", _count, pdf.Data()), "", RooArgList(sig1, CB), sig1frac);

	RooRealVar* nsigMC = new RooRealVar(Form("nsigMC%d_%s", _count, pdf.Data()), "", dsMC->sumEntries(), 0.9 * dsMC->sumEntries(), 1.1 * dsMC->sumEntries());
	RooAddPdf* modelMC = new RooAddPdf(Form("modelMC%d_%s", _count, pdf.Data()), "", RooArgList(*sig), RooArgList(*nsigMC));
	scale->setConstant();

	mass->setRange("signal", init_mean - 0.07, init_mean + 0.07);
	if (tree == "ntKp") mass->setRange("signal", init_mean - 0.1, init_mean + 0.1);
	else if (tree == "ntmix") mass->setRange("signal", init_mean - 0.035, init_mean + 0.035);
	else if (tree == "ntmix_psi2s") mass->setRange("signal", init_mean - 0.03, init_mean + 0.03);

	RooFitResult* fitResultMC = modelMC->fitTo(*dsMC, Save(), Extended(), Range("signal"));
	w.import(*nsigMC);

	cMC->Clear();
	cMC->cd();
	TPad* pMC1 = new TPad(Form("pMC1_%d", _count), Form("pMC1_%d", _count), 0., 0., 1., 1.);
	pMC1->SetBorderMode(1);
	pMC1->SetFrameBorderMode(0);
	pMC1->SetBorderSize(2);
	pMC1->SetBottomMargin(0.22);
	pMC1->SetLeftMargin(0.14);
	pMC1->SetRightMargin(0.04);
	pMC1->Draw();

	double xMinPlot = mass->getMin();
	double xMaxPlot = mass->getMax();
	if (tree == "ntphi") { xMinPlot = 5.2; xMaxPlot = 5.55; }
	else if (tree == "ntKp") { xMinPlot = 5.05; xMaxPlot = 5.55; }
	else if (tree == "ntKstar") { xMinPlot = 5.05; xMaxPlot = 5.55; }

	TString xTtile_decayC = "";
	if (tree == "ntmix" || tree == "ntmix_psi2s") xTtile_decayC = "m_{J/#psi #pi^{+} #pi^{-}} [GeV/c^{2}]";
	else if (tree == "ntphi") xTtile_decayC = "m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]";
	else if (tree == "ntKstar") xTtile_decayC = "m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]";
	else if (tree == "ntKp") xTtile_decayC = "m_{J/#psi K^{+}} [GeV/c^{2}]";

	RooPlot* frameMC = mass->frame(Range(xMinPlot, xMaxPlot));
	frameMC->SetTitle("");
	frameMC->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})", (mass->getMax() - mass->getMin()) / NBIN * 1000));
	frameMC->GetYaxis()->SetTitleOffset(2.15);
	frameMC->GetYaxis()->SetTitleSize(0.035);
	frameMC->GetXaxis()->SetTitle(xTtile_decayC);
	frameMC->GetXaxis()->SetTitleSize(0.030);
	frameMC->GetXaxis()->SetTitleOffset(1.10);
	frameMC->GetXaxis()->CenterTitle();
	frameMC->GetYaxis()->SetTitleFont(42);
	frameMC->GetXaxis()->SetLabelFont(42);
	frameMC->GetYaxis()->SetLabelFont(42);
	frameMC->GetXaxis()->SetLabelOffset(0.012);
	frameMC->GetXaxis()->SetLabelSize(0.031);
	frameMC->GetXaxis()->SetTickLength(0.035);
	frameMC->GetYaxis()->SetTitleSize(0.027);
	frameMC->GetYaxis()->SetLabelSize(0.027);
	frameMC->SetStats(0);

	pMC1->cd();
	dsMC->plotOn(frameMC, Name(Form("dsMC%d_%s", _count, pdf.Data())), Binning(NBIN), MarkerSize(0.5), MarkerStyle(8), LineColor(1), LineWidth(1));
	modelMC->plotOn(frameMC, Name(Form("sigMC%d_%s", _count, pdf.Data())), Range("signal"), NormRange("signal"), Normalization(nsigMC->getVal(), RooAbsReal::NumEvent), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(1));
	modelMC->plotOn(frameMC, Name(Form("modelMCcurve%d_%s", _count, pdf.Data())), DrawOption("L"), LineWidth(0));
	modelMC->paramOn(frameMC, Layout(0.15, 0.3, 0.85), Format("NEU", AutoPrecision(2)));
	frameMC->getAttFill()->SetFillStyle(0);
	frameMC->getAttLine()->SetLineWidth(0);
	frameMC->Draw();
	cMC->RedrawAxis();

	TLegend* legMC = new TLegend(0.6, 0.78, 0.92, 0.90, NULL, "brNDC");
	legMC->SetBorderSize(0);
	legMC->SetTextSize(0.035);
	legMC->SetTextFont(42);
	legMC->SetFillStyle(0);
	legMC->AddEntry(frameMC->findObject(Form("dsMC%d_%s", _count, pdf.Data())), "Signal MC", "lp");
	TString ident_par = "";
	if (tree == "ntphi") ident_par = "B^{0}_{s}";
	else if (tree == "ntKp") ident_par = "B^{+}";
	else if (tree == "ntKstar") ident_par = "B^{0}";
	else if (tree == "ntmix") ident_par = "X(3872)";
	else if (tree == "ntmix_psi2s") ident_par = "#psi(2S)";
	legMC->AddEntry(frameMC->findObject(Form("sigMC%d_%s", _count, pdf.Data())), Form("%s, %s ", ident_par.Data(), pdf.Data()), "f");
	legMC->Draw();

	double n_signal_initial = ds->sumEntries(TString::Format("abs(Bmass-%g)<%g", init_mean, (tree == "ntmix" || tree == "ntmix_psi2s") ? 0.005 : 0.05));
	RooRealVar nsig(Form("nsig%d_%s", _count, pdf.Data()), "", n_signal_initial * 0.4, n_signal_initial * 0.2, n_signal_initial * 1.5);

	RooRealVar nbkg(Form("nbkg%d_%s", _count, pdf.Data()), "", ds->sumEntries() * 0.7, ds->sumEntries() * 0.1, ds->sumEntries());
	RooRealVar a0(Form("a0%d_%s", _count, pdf.Data()), "", -0.35, -2, 2);
	RooRealVar a1(Form("a1%d_%s", _count, pdf.Data()), "", -0.05, -2, 2);
	RooRealVar a2(Form("a2%d_%s", _count, pdf.Data()), "", 0.01, -2, 2);
	RooRealVar a3(Form("a3%d_%s", _count, pdf.Data()), "", 0, -2, 2);
	RooChebychev bkg_2nd(Form("bkg%d_%s", _count, pdf.Data()), "", *mass, RooArgList(a0, a1));
	RooChebychev bkg_3rd(Form("bkg%d_%s", _count, pdf.Data()), "", *mass, RooArgSet(a0, a1, a2));
	RooChebychev bkg_4th(Form("bkg%d_%s", _count, pdf.Data()), "", *mass, RooArgSet(a0, a1, a2, a3));
	RooRealVar lambda(Form("lambda%d_%s", _count, pdf.Data()), "lambda", -0.5, -5., 0.1);
	RooExponential bkg(Form("bkg%d_%s", _count, pdf.Data()), "", *mass, lambda);

	RooRealVar nbkg_part_r(Form("nbkg_part_r%d_%s", _count, pdf.Data()), "", 6000, 250, 1e4);
	RooRealVar* m_nonprompt_scale = new RooRealVar(Form("m_nonprompt_scale%d_%s", _count, ""), "m_nonprompt_scale", 0.01, 0.001, 0.1);
	RooRealVar* m_nonprompt_shift = new RooRealVar(Form("m_nonprompt_shift%d_%s", _count, ""), "m_nonprompt_shift", 5.15, 5.1, 5.2);
	RooGenericPdf* erfc = new RooGenericPdf(Form("erfc%d", _count), "0.5*TMath::Erfc((@0-@2)/@1)", RooArgList(*mass, *m_nonprompt_scale, *m_nonprompt_shift));

	RooAddPdf* model = nullptr;
	if (tree == "ntmix" || tree == "ntmix_psi2s") {
		if ((variation == "" && pdf == "") || variation == "signal") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg_3rd), RooArgList(nsig, nbkg));
		if (variation == "background" && pdf == "4th") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg_4th), RooArgList(nsig, nbkg));
	}
	if (tree == "ntphi") {
		if ((variation == "" && pdf == "") || variation == "signal") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg_3rd), RooArgList(nsig, nbkg));
		if (variation == "background" && pdf == "2nd") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg_2nd), RooArgList(nsig, nbkg));
		if (variation == "background" && pdf == "4th") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg_4th), RooArgList(nsig, nbkg));
	}
	if (tree == "ntKstar") {
		if ((variation == "" && pdf == "") || variation == "signal") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg_3rd), RooArgList(nsig, nbkg));
		if (variation == "background" && pdf == "2nd") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg_2nd), RooArgList(nsig, nbkg));
		if (variation == "background" && pdf == "3rd") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg_3rd), RooArgList(nsig, nbkg));
	}
	if (tree == "ntKp") {
		if ((variation == "" && pdf == "")) model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg, *erfc), RooArgList(nsig, nbkg, nbkg_part_r));
		if (pdf == "mass_range") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg), RooArgList(nsig, nbkg));
		if (variation == "background" && pdf == "2nd") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(bkg_2nd, *sig, *erfc), RooArgList(nbkg, nsig, nbkg_part_r));
		if (variation == "background" && pdf == "3rd") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(bkg_3rd, *sig, *erfc), RooArgList(nbkg, nsig, nbkg_part_r));
		if (variation == "signal" && pdf == "1gauss") model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(bkg, sig1, *erfc), RooArgList(nbkg, nsig, nbkg_part_r));
		if (variation == "signal" && (pdf == "3gauss" || pdf == "fixed" || pdf == "gauss_cb")) model = new RooAddPdf(Form("model%d_%s", _count, pdf.Data()), "", RooArgList(*sig, bkg, *erfc), RooArgList(nsig, nbkg, nbkg_part_r));
	}

	scale->setConstant(false);
	sigma1.setConstant();
	if (pdf != "1gauss") {
		sigma2.setConstant();
		sig1frac.setConstant();
	}
	if (variation == "signal" && pdf == "3gauss") {
		sigma3.setConstant();
		sig2frac.setConstant();
	}
	if (variation == "signal" && pdf == "gauss_cb") {
		sigma4cb.setConstant();
		n.setConstant();
		alpha.setConstant();
	}
	if (variation == "signal" && pdf == "fixed") mean.setConstant();

	TString fitRange = (pdf == "mass_range") ? "m_rangeB" : "all";
	RooFitResult* fitResult = model->fitTo(*ds, Save(), Extended(kTRUE), Range(fitRange));
	fitResult->Print("v");
	w.import(*model);

	c->cd();
	RooPlot* frame = mass->frame(Range(xMinPlot, xMaxPlot));
	frame->SetStats(0);
	frame->SetTitle("");
	frame->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})", (mass->getMax() - mass->getMin()) / NBIN * 1000));
	frame->GetYaxis()->SetTitleOffset(2.);
	frame->GetXaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleSize(0.035);
	frame->GetXaxis()->SetTitleSize(0.);
	frame->GetXaxis()->SetTitleFont(0);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelSize(0);
	frame->GetYaxis()->SetLabelSize(0.035);
	if (tree == "ntmix" || tree == "ntmix_psi2s" && true) { frame->GetYaxis()->SetRangeUser(0.0, GetNtmixDataYmax(binmin, binmax));}	
	TPad* p1 = new TPad("p1", "p1", 0., 0.22, 1., 1);
	p1->SetBorderMode(1);
	p1->SetFrameBorderMode(0);
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.01);
	p1->SetLeftMargin(0.14);
	p1->SetRightMargin(0.04);
	p1->Draw();
	TPad* p2 = new TPad("p2", "p2", 0., 0., 1., 0.22);
	p2->SetTopMargin(0.0);
	p2->SetBottomMargin(0.34);
	p2->SetLeftMargin(0.14);
	p2->SetRightMargin(0.04);
	p2->SetBorderMode(0);
	p2->SetBorderSize(2);
	p2->SetFrameBorderMode(0);
	p2->SetTicks(1, 1);
	p2->Draw();

	p1->cd();
	ds->plotOn(frame, Name(Form("ds_cut%d", _count)), Binning(NBIN), MarkerSize(0.5), MarkerStyle(8), MarkerColor(1), LineColor(1), LineWidth(1));
	model->plotOn(frame, Name(Form("model%d_%s", _count, pdf.Data())), Range(fitRange), NormRange(fitRange), Precision(1e-6), DrawOption("L"), LineColor(2), LineWidth(1));
	model->plotOn(frame, Name(Form("sig%d_%s", _count, pdf.Data())), Components(*sig), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(1));
	if (tree == "ntKp") model->plotOn(frame, RooFit::Name(Form("erfc%d_%s", _count, "")), Components(*erfc), Range(fitRange), NormRange(fitRange), LineColor(kGreen + 3), LineStyle(9), LineWidth(2), DrawOption("L"));
	model->plotOn(frame, Name(Form("bkg%d_%s", _count, pdf.Data())), Components(bkg), Range(fitRange), NormRange(fitRange), Precision(1e-6), DrawOption("L"), LineStyle(7), LineColor(4), LineWidth(1));
	double chi2Ndf = frame->chiSquare(Form("model%d_%s", _count, pdf.Data()), Form("ds_cut%d", _count), fitResult->floatParsFinal().getSize());
	if (!std::isfinite(chi2Ndf) || chi2Ndf < 0) chi2Ndf = -1.0;
	RooRealVar chi2Var(Form("chi2_data_norm%d_%s", _count, pdf.Data()), "", chi2Ndf);
	w.import(chi2Var);
	frame->Draw();

	TLegend* leg = new TLegend(0.68, 0.60, 0.92, 0.90, NULL, "brNDC");
	if (tree != "ntKp" && tree != "ntmix") leg = new TLegend(0.68, 0.66, 0.92, 0.90, NULL, "brNDC");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);
	leg->AddEntry(frame->findObject(Form("ds_cut%d", _count)), "Data", "LEP");
	leg->AddEntry(frame->findObject(Form("model%d_%s", _count, pdf.Data())), "Fit Model", "l");
	leg->AddEntry(frame->findObject(Form("bkg%d_%s", _count, pdf.Data())), " Comb. Bkg.", "l");
	if (tree == "ntKp") {
		leg->AddEntry(frame->findObject(Form("sig%d_%s", _count, pdf.Data())), " B^{+} #rightarrow J/#psi K^{+}", "f");
		leg->AddEntry(frame->findObject(Form("erfc%d_%s", _count, pdf.Data())), " B #rightarrow J/#psi X", "l");
	} else if (tree == "ntmix") {
		leg->AddEntry(frame->findObject(Form("sig%d_%s", _count, pdf.Data())), " X(3872) #rightarrow J/#psi #pi^{+} #pi^{-}", "f");
	} else if (tree == "ntmix_psi2s") {
		leg->AddEntry(frame->findObject(Form("sig%d_%s", _count, pdf.Data())), " #psi(2S) #rightarrow J/#psi #pi^{+} #pi^{-}", "f");
	} else if (tree == "ntphi") {
		leg->AddEntry(frame->findObject(Form("sig%d_%s", _count, pdf.Data())), " B_{s}^{0} #rightarrow J/#psi K^{+}K^{-}", "f");
	} else if (tree == "ntKstar") {
		leg->AddEntry(frame->findObject(Form("sig%d_%s", _count, pdf.Data())), " B^{0} #rightarrow J/#psi K^{*}", "f");
	}
	leg->Draw();

	p2->cd();
	RooHist* pull_hist = frame->pullHist(Form("ds_cut%d", _count), Form("model%d_%s", _count, pdf.Data()));
	pull_hist->SetMarkerSize(0.5);
	RooPlot* pull_plot = mass->frame(Range(xMinPlot, xMaxPlot));
	pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist), "XP");
	pull_plot->SetTitle("");
	pull_plot->SetXTitle(xTtile_decayC.Data());
	pull_plot->SetYTitle("Pull");
	pull_plot->GetYaxis()->SetTitleFont(42);
	pull_plot->GetYaxis()->SetTitleSize(0.13);
	pull_plot->GetYaxis()->CenterTitle(kTRUE);
	pull_plot->GetYaxis()->SetLabelOffset(0.01);
	pull_plot->GetYaxis()->SetLabelFont(42);
	pull_plot->GetYaxis()->SetLabelSize(0.12);
	pull_plot->GetYaxis()->SetNdivisions(305);
	pull_plot->GetYaxis()->SetTitleOffset(0.55);
	pull_plot->GetYaxis()->SetRangeUser(-3.5, 3.5);
	pull_plot->GetXaxis()->SetTitleSize(0.13);
	pull_plot->GetXaxis()->SetTitleOffset(1.0);
	pull_plot->GetXaxis()->CenterTitle();
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelOffset(0.01);
	pull_plot->GetXaxis()->SetLabelSize(0.14);
	pull_plot->GetXaxis()->SetTickLength(0.16);
	pull_plot->Draw();
	TLine* line_ref = new TLine(xMinPlot, 0., xMaxPlot, 0.);
	line_ref->SetLineStyle(1);
	line_ref->SetLineColor(2);
	line_ref->SetLineWidth(1);
	line_ref->Draw("same");

	cout << "\n-------------------------------------------------------------------------------------- \n" << endl;
	cout << "Signal Yield Y_s = " << nsig.getVal() << "     yield Error = " << nsig.getError();
	cout << "\n-------------------------------------------------------------------------------------- \n" << endl;

	p1->cd();
	return fitResult;
}
// END OF MAIN FITTING FUNCTION

















void validate_fit(
  RooWorkspace* w, 
  TString pdf, 
  TString tree, 
  TString variable, 
  int full, 
  float binmin, 
  float binmax, 
  string Path)
{
	std::cout << "Performing Check on Fit" << std::endl;
	RooRealVar mass = *(w->var("mass"));
	RooAbsPdf* model  = w->pdf(Form("model%d_%s",_count,pdf.Data()));
	//RooDataSet* data = (RooDataSet*) w->data("data");

	//model->fitTo(*data);
	vector<RooRealVar> params;
	params.push_back(*(w->var(Form("nsig%d_%s",_count,pdf.Data()))));
	//params.push_back(*(w->var(Form("mean%d_%s",_count,pdf.Data()))));

	double n_signal_init = params[0].getVal();
	double n_signal_error_init = params[0].getError();
	int params_size = params.size();

	RooMCStudy* mcstudy = new RooMCStudy(*model, mass,  Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));
	mcstudy->generateAndFit(5000);

	cout << "DONE Generate and Fit " << endl;

	TString XName[2] = {"P(Y)","mean"};
	vector<RooPlot*> framesPull, framesParam, framesError;

	for(int i = 0; i < params_size; ++i)
	{
		framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(50),FrameRange(-5,5)));
		framesPull[i]->SetTitle("");
		framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(50)));
		framesParam[i]->SetTitle("");
		framesError.push_back(mcstudy->plotError(params.at(i),FrameBins(50)));
		framesError[i]->SetTitle("");
	}

	vector<TGraph*> h1;
	vector<TGraph*> h2;
	vector<TGraph*> h3;

	for(int i = 0; i < params_size; ++i){
		h1.push_back(static_cast<TGraph*>(framesPull.at(i)->getObject(0)));
		h2.push_back(static_cast<TGraph*>(framesParam.at(i)->getObject(0)));
		h3.push_back(static_cast<TGraph*>(framesError.at(i)->getObject(0)));
	}

	TCanvas* c_pull = new TCanvas("pulls", "pulls", 700, 700);
	TCanvas* c_params = new TCanvas("params", "params", 700, 700);
	TCanvas* c_errors = new TCanvas("errors", "errors", 700, 700);
	gStyle->SetOptFit(0111);
	gPad->SetLeftMargin(0.15);
	gStyle->SetStatX(0.95);		//Stat box x position (top right hand corner)	
	gStyle->SetStatY(0.9); 		//Stat box y position 	
	gStyle->SetStatW(0.1);	 		//Stat box width as fraction of pad size	0.05	
	gStyle->SetStatFont(62);  		//Stat box font
	gStyle->SetStatFontSize(0);
	gStyle->SetStatBorderSize(0);
	gStyle->SetStatColor(0);
	gStyle->SetStatStyle(0);		//Stat box fill style hollow

	for(int i = 0; i < params_size; ++i){

		n_signal_init = params[i].getVal();
		n_signal_error_init = params[i].getError();

		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Fit("gaus","","",-5,5);
		
		c_pull->Update();
		h1[i]->GetFunction("gaus")->SetLineColor(kCyan+1);
		h1[i]->GetFunction("gaus")->SetLineWidth(1);
		h1[i]->GetFunction("gaus")->SetFillStyle(3019);
		h1[i]->GetFunction("gaus")->SetFillColor(kCyan+1);
		h1[i]->GetXaxis()->SetTitle(Form("%s",XName[i].Data()));
		h1[i]->GetYaxis()->SetTitle("Toy MCs");
		h1[i]->GetXaxis()->SetRangeUser(-5, 5);
		h1[i]->Draw("APsame");
	
		c_errors->cd();
		h3[i]->SetTitle("");
		h3[i]->Fit("gaus","","",n_signal_error_init*0.5,n_signal_error_init*1.5);
		
		c_errors->Update();
		h3[i]->GetFunction("gaus")->SetLineColor(4);
		h3[i]->GetFunction("gaus")->SetLineWidth(1);
		h3[i]->GetXaxis()->SetTitle(Form("%s Error",XName[i].Data()));
		h3[i]->GetYaxis()->SetTitle("");
		h3[i]->Draw("APsame");

		c_params->cd();
		h2[i]->SetTitle("");
		h2[i]->Fit("gaus","","",n_signal_init * 0.5,n_signal_init * 1.5);

		c_params->Update();
		h2[i]->GetFunction("gaus")->SetLineColor(4);
		h2[i]->GetFunction("gaus")->SetLineWidth(1);
		h2[i]->GetXaxis()->SetTitle(Form("%s Mean",XName[i].Data()));
		h2[i]->GetYaxis()->SetTitle("");
		h2[i]->Draw("APsame");

		c_pull->SaveAs(Form("%s/pull_signal_%s_%0.1f_%0.1f_%d_%s.pdf",Path.c_str(),variable.Data(),binmin,binmax,i,tree.Data()));
		c_params->SaveAs(Form("%s/param_signal_%s_%0.1f_%0.1f_%d_%s.pdf",Path.c_str(),variable.Data(),binmin,binmax,i,tree.Data()));
		c_errors->SaveAs(Form("%s/error_signal_%s_%0.1f_%0.1f_%d_%s.pdf",Path.c_str(),variable.Data(),binmin,binmax,i,tree.Data()));
	}
}
