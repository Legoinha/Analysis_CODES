#pragma once

// Combined (simultaneous) fit of several PbPb years: the model of one year, the fit of one analysis bin and the
// drawings. roofitCombined.C runs the bins and writes the products. The model is the one of ../roofitB.h.

#include "../aux/uti.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "TTreeFormula.h"
#include "TMatrixDSym.h"
#include "Math/QuantFuncMathCore.h"
#include <map>
#include <memory>

using namespace RooFit;
using namespace std;

// Everything of the fit and its drawing that depends on the tree, with the values of roofitB.C and roofitB.h. A tree
// that is not here breaks the lookup (std::map::at).
struct TreeSetup {
	TString dataTree;           // tree of the DATA file
	double peak;                // signal mass (plotER/aux/masses.h)
	double massMin, massMax;    // fitted mass range
	int nMassBins;              // mass bins of the plots and of the chi2
	double signalHalfWidth;     // MC fit window around the peak
	double peakCountHalfWidth;  // DATA counted around the peak for the initial signal yield
	TString background;         // nominal background: 2nd, 3rd (Chebychev), expo (exponential)
	bool partialReco;           // erfc of the partially reconstructed B
	bool significance;          // the plots show the significance of the signal (per year and of the combined fit)
	double plotMin, plotMax;    // drawn mass range
	TString massTitle;
	TString decay;              // legend of the signal
};

inline const TreeSetup& GetTreeSetup(const TString& tree)
{
	static const std::map<TString, TreeSetup> setups = {
		{"ntmix_X3872", {"ntmix", X3872_MASS, 3.8, 4.0, 20, 0.035, 0.005, "2nd", false, true, 3.8, 4.0,
		                 "m_{J/#psi #pi^{+} #pi^{-}} [GeV/c^{2}]", " X(3872) #rightarrow J/#psi #pi^{+} #pi^{-}"}},
		{"ntmix_PSI2S", {"ntmix", PSI2S_MASS, 3.6, 3.8, 20, 0.03, 0.005, "2nd", false, false, 3.6, 3.8,
		                 "m_{J/#psi #pi^{+} #pi^{-}} [GeV/c^{2}]", " #psi(2S) #rightarrow J/#psi #pi^{+} #pi^{-}"}},
		{"ntKp", {"ntKp", Bu_MASS, 5.0, 5.8, 40, 0.1, 0.05, "expo", true, false, 5.05, 5.55,
		          "m_{J/#psi K^{+}} [GeV/c^{2}]", " B^{+} #rightarrow J/#psi K^{+}"}},
		{"ntKstar", {"ntKstar", Bd_MASS, 5.0, 5.8, 40, 0.07, 0.05, "3rd", false, true, 5.05, 5.55,
		             "m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]", " B^{0} #rightarrow J/#psi K^{*}"}},
		{"ntphi", {"ntphi", Bs_MASS, 5.0, 5.8, 40, 0.07, 0.05, "3rd", false, false, 5.2, 5.55,
		           "m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]", " B_{s}^{0} #rightarrow J/#psi K^{+}K^{-}"}},
	};
	return setups.at(tree);
}

// The model of one year in one analysis bin, for the nominal fit or one model variation. pars holds every parameter
// in the order the builder creates it: the order is the same for every year and for the per-year export, so two
// builds of the same model are matched parameter by parameter. shapePars are the signal-shape parameters the MC
// fit determines and the data fit keeps constant.
struct YearModel {
	std::vector<RooRealVar*> pars;
	std::vector<RooRealVar*> shapePars;
	RooRealVar *mean = nullptr, *scale = nullptr, *sigma1 = nullptr, *sigma2 = nullptr, *sig1frac = nullptr;
	RooRealVar *nsig = nullptr, *nbkg = nullptr, *nbkgPart = nullptr, *nsigMC = nullptr;
	RooAbsPdf *sig = nullptr, *bkg = nullptr, *erfc = nullptr, *model = nullptr, *modelMC = nullptr;
};

// Builds the model of one year as roofitB.h builds it: the double-Gaussian signal (or the signal variation SIGCODE:
// 1gauss, 3gauss, gauss_cb, fixed), the nominal background of the tree (TreeSetup; or the background variation
// BKGCODE: 2nd, 3rd, 4th) and, for ntKp, the partially reconstructed B (erfc). Only the yields belong to the year
// (nsig, nbkg, nbkg_part_r and the MC normalization nsigMC): their names are those of roofitB.h (nsig1_, ...) with
// "_<YEAR>" appended. Every other parameter (mean, scale, MC signal shape, background and erfc shapes) carries no
// year and is created once per POOL, so all years of one fit hold the same object. The pdfs belong to the year.
// YEAR = "" builds the per-year export, whose names are exactly those of roofitB.h.
YearModel BuildYearModel(const TString& tree, int bin, const TString& sigCode, const TString& bkgCode, const TString& year,
                         std::map<TString, RooRealVar*>& pool, RooRealVar* mass, RooDataSet* data, RooDataSet* mc)
{
	YearModel m;
	const TreeSetup& setup = GetTreeSetup(tree);
	const TString code = sigCode + bkgCode;
	const char* c = code.Data();
	auto name = [&](TString base, bool perYear) {
		if (year != "" && perYear) base += (base.EndsWith("_") ? "" : "_") + year;
		return base;
	};
	auto par = [&](TString base, bool perYear, double value, double low, double high) {
		const TString n = name(base, perYear);
		if (!pool.count(n)) pool[n] = new RooRealVar(n, n, value, low, high);
		m.pars.push_back(pool[n]);
		return pool[n];
	};

	const double peak = setup.peak;
	m.mean = par(Form("mean%d_%s", bin, c), false, peak, peak - 0.01, peak + 0.01);
	m.scale = par("scale", false, 1, 0.9, 1.15);
	m.sigma1 = par(Form("sigma1%d_%s", bin, c), false, 0.01, 0.001, 0.1);
	m.shapePars.push_back(m.sigma1);
	RooAbsPdf* sig1 = new RooGaussian(name(Form("sig1%d_%s", bin, c), true), "", *mass, *m.mean,
		*new RooProduct(name(Form("scaled_sigma1%d_%s", bin, c), true), "", RooArgList(*m.scale, *m.sigma1)));
	if (sigCode == "1gauss") m.sig = sig1;
	else if (sigCode == "gauss_cb") {
		RooRealVar* sigmaCB = par(Form("sigma4cb%d_%s", bin, c), false, 0.005, 0.001, 0.05);
		RooRealVar* alpha = par(Form("alpha%d_%s", bin, c), false, 4., 0., 15.);
		RooRealVar* n = par(Form("n_%d_%s", bin, c), false, 10., -100., 200.);
		m.sig1frac = par(Form("sig1frac%d_%s", bin, c), false, 0.5, 0.01, 1.);
		m.shapePars.insert(m.shapePars.end(), {sigmaCB, alpha, n, m.sig1frac});
		RooAbsPdf* cb = new RooCBShape(name(Form("CB%d_%s", bin, c), true), "", *mass, *m.mean,
			*new RooProduct(name(Form("scaled_sigmacb%d_%s", bin, c), true), "", RooArgList(*m.scale, *sigmaCB)), *alpha, *n);
		m.sig = new RooAddPdf(name(Form("sig_gaussCB%d_%s", bin, c), true), "", RooArgList(*sig1, *cb), RooArgList(*m.sig1frac));
	} else {
		m.sigma2 = par(Form("sigma2%d_%s", bin, c), false, 0.005, 0.001, 0.1);
		m.sig1frac = par(Form("sig1frac%d_%s", bin, c), false, 0.5, 0.01, 1.);
		m.shapePars.insert(m.shapePars.end(), {m.sigma2, m.sig1frac});
		RooAbsPdf* sig2 = new RooGaussian(name(Form("sig2%d_%s", bin, c), true), "", *mass, *m.mean,
			*new RooProduct(name(Form("scaled_sigma2%d_%s", bin, c), true), "", RooArgList(*m.scale, *m.sigma2)));
		if (sigCode == "3gauss") {
			RooRealVar* sigma3 = par(Form("sigma3%d_%s", bin, c), false, 0.01, 0.001, 0.03);
			RooRealVar* sig2frac = par(Form("sig2frac%d_%s", bin, c), false, 0.5, 0.01, 1.);
			m.shapePars.insert(m.shapePars.end(), {sigma3, sig2frac});
			RooAbsPdf* sig3 = new RooGaussian(name(Form("sig3%d_%s", bin, c), true), "", *mass, *m.mean,
				*new RooProduct(name(Form("scaled_sigma3%d_%s", bin, c), true), "", RooArgList(*m.scale, *sigma3)));
			m.sig = new RooAddPdf(name(Form("sig_tripleG%d_%s", bin, c), true), "", RooArgList(*sig1, *sig2, *sig3),
			                      RooArgList(*m.sig1frac, *sig2frac), true);
		} else {
			m.sig = new RooAddPdf(name(Form("sig_doubleG%d_%s", bin, c), true), "", RooArgList(*sig1, *sig2), RooArgList(*m.sig1frac));
		}
	}

	m.nsigMC = par(Form("nsigMC%d_%s", bin, c), true, mc->sumEntries(), 0.9 * mc->sumEntries(), 1.1 * mc->sumEntries());
	m.modelMC = new RooAddPdf(name(Form("modelMC%d_%s", bin, c), true), "", RooArgList(*m.sig), RooArgList(*m.nsigMC));

	const double nData = data->sumEntries();
	const double nPeak = data->sumEntries(Form("abs(Bmass-%g)<%g", peak, setup.peakCountHalfWidth));
	m.nsig = par(Form("nsig%d_%s", bin, c), true, std::max(0.0, 0.4 * nPeak), 0.0, std::max(10.0, 2.0 * nPeak));
	m.nbkg = par(Form("nbkg%d_%s", bin, c), true, 0.7 * nData, 0.1 * nData, nData);

	const TString bkgType = bkgCode == "" ? setup.background : bkgCode;
	if (bkgType == "expo") {
		RooRealVar* lambda = par(Form("lambda%d_%s", bin, c), false, -0.5, -5., 0.1);
		m.bkg = new RooExponential(name(Form("bkg%d_%s", bin, c), true), "", *mass, *lambda);
	} else {
		const double init[] = {-0.35, -0.05, 0.01, 0.0};
		RooArgList coefficients;
		for (int k = 0; k < std::atoi(bkgType.Data()); ++k) coefficients.add(*par(Form("a%d%d_%s", k, bin, c), false, init[k], -2., 2.));
		m.bkg = new RooChebychev(name(Form("bkg%d_%s", bin, c), true), "", *mass, coefficients);
	}

	if (setup.partialReco) {
		m.nbkgPart = par(Form("nbkg_part_r%d_%s", bin, c), true, std::max(1.0, 0.05 * nData), 0.0, nData);
		RooRealVar* erfcScale = par(Form("m_nonprompt_scale%d_%s", bin, c), false, 0.01, 0.001, 0.1);
		RooRealVar* erfcShift = par(Form("m_nonprompt_shift%d_%s", bin, c), false, 5.15, 5.1, 5.2);
		m.erfc = new RooGenericPdf(name(code == "" ? Form("erfc%d", bin) : Form("erfc%d_%s", bin, c), true),
		                           "0.5*TMath::Erfc((@0-@2)/@1)", RooArgList(*mass, *erfcScale, *erfcShift));
		m.model = new RooAddPdf(name(Form("model%d_%s", bin, c), true), "", RooArgList(*m.sig, *m.bkg, *m.erfc),
		                        RooArgList(*m.nsig, *m.nbkg, *m.nbkgPart));
	} else {
		m.model = new RooAddPdf(name(Form("model%d_%s", bin, c), true), "", RooArgList(*m.sig, *m.bkg), RooArgList(*m.nsig, *m.nbkg));
	}
	return m;
}

// The per-year export of a fitted year: the same model built with the names of roofitB.h, every parameter at the
// value, error and constant flag its counterpart has in VALUES (a snapshot of the combined fit). It is written for
// validation and Acc x Eff (nominal model) and drawn on the plots of the year, which so read as those of roofitB.C.
YearModel ExportYearModel(const TString& tree, int bin, const TString& sigCode, const TString& bkgCode, const YearModel& fitted,
                          const RooArgSet& values, RooRealVar* mass, RooDataSet* data, RooDataSet* mc)
{
	std::map<TString, RooRealVar*> pool;
	YearModel exported = BuildYearModel(tree, bin, sigCode, bkgCode, "", pool, mass, data, mc);
	for (size_t k = 0; k < exported.pars.size(); ++k) {
		const RooRealVar* v = static_cast<const RooRealVar*>(values.find(fitted.pars[k]->GetName()));
		exported.pars[k]->setVal(v->getVal());
		exported.pars[k]->setError(v->getError());
		exported.pars[k]->setConstant(v->isConstant());
	}
	return exported;
}

// The candidates of TREENAME in FILE that pass CUT inside the mass window, as a RooDataSet of VARS, weighted by
// pThatreweight when VARS holds it (MC). Only the branches the cut and VARS read are copied: a scored DATA tree
// holds ~80 branches and several GB.
RooDataSet* ReadSample(const TString& file, const TString& treeName, const TString& cut, const RooArgSet& vars,
                       const TString& name, double massMin, double massMax)
{
	const TString fullCut = Form("(%s) && (Bmass > %f) && (Bmass < %f)", cut.Data(), massMin, massMax);
	cout << "Reading " << treeName << " from " << file << endl << "  selection: " << fullCut << endl;
	TFile* input = TFile::Open(file);
	TTree* tree = input->Get<TTree>(treeName);
	std::vector<TString> branches;
	for (const RooAbsArg* v : vars) branches.push_back(v->GetName());
	{
		TTreeFormula formula("readSampleCut", fullCut, tree);
		for (int i = 0; i < formula.GetNcodes(); ++i) branches.push_back(formula.GetLeaf(i)->GetBranch()->GetName());
	}
	tree->SetBranchStatus("*", 0);
	for (const TString& b : branches) tree->SetBranchStatus(b, 1);
	gROOT->cd();
	TTree* selected = tree->CopyTree(fullCut);
	RooDataSet* ds = vars.find("pThatreweight")
		? new RooDataSet(name, name, selected, vars, "", "pThatreweight")
		: new RooDataSet(name, name, selected, vars);
	input->Close();
	cout << "  " << ds->numEntries() << " candidates (sum of weights " << ds->sumEntries() << ")" << endl;
	return ds;
}

// Draws one MC signal fit as roofitB.h draws it.
void DrawMCFit(TCanvas* cMC, const TString& tree, RooRealVar* mass, RooDataSet* mc, const YearModel& m, int nMassBins)
{
	static int nDrawn = 0;
	++nDrawn;
	const TreeSetup& setup = GetTreeSetup(tree);
	const double xMin = setup.plotMin, xMax = setup.plotMax;
	const TString xTitle = setup.massTitle;
	const int signalColor = (tree == "ntmix_PSI2S") ? kOrange - 2 : kOrange - 3;

	cMC->Clear();
	cMC->cd();
	TPad* pad = new TPad(Form("pMC1_%d", nDrawn), "", 0., 0., 1., 1.);
	pad->SetBorderMode(1);
	pad->SetFrameBorderMode(0);
	pad->SetBorderSize(2);
	pad->SetBottomMargin(0.22);
	pad->SetLeftMargin(0.14);
	pad->SetRightMargin(0.04);
	pad->Draw();

	RooPlot* frame = mass->frame(Range(xMin, xMax));
	frame->SetTitle("");
	frame->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})", (mass->getMax() - mass->getMin()) / nMassBins * 1000));
	frame->GetYaxis()->SetTitleOffset(2.15);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetTitleSize(0.027);
	frame->GetYaxis()->SetLabelSize(0.027);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetTitle(xTitle);
	frame->GetXaxis()->SetTitleSize(0.030);
	frame->GetXaxis()->SetTitleOffset(1.10);
	frame->GetXaxis()->CenterTitle();
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelOffset(0.012);
	frame->GetXaxis()->SetLabelSize(0.031);
	frame->GetXaxis()->SetTickLength(0.035);
	frame->SetStats(0);

	pad->cd();
	const TString dataName = Form("dsMC_%d", nDrawn), sigName = Form("sigMC_%d", nDrawn);
	mc->plotOn(frame, Name(dataName), Binning(nMassBins), MarkerSize(0.5), MarkerStyle(8), LineColor(1), LineWidth(1));
	m.modelMC->plotOn(frame, Name(sigName), Range("signal"), NormRange("signal"), Normalization(m.nsigMC->getVal(), RooAbsReal::NumEvent),
	                  DrawOption("LF"), FillStyle(3002), FillColor(signalColor), LineStyle(7), LineColor(signalColor), LineWidth(1));
	m.modelMC->plotOn(frame, DrawOption("L"), LineWidth(0));
	m.modelMC->paramOn(frame, Layout(0.18, 0.48, 0.78), Format("NEU", AutoPrecision(2)));
	frame->getAttFill()->SetFillStyle(0);
	frame->getAttLine()->SetLineWidth(0);
	frame->Draw();
	TLatex* particle = new TLatex(0.2, 0.85, FitParticleLabel(tree, true));
	setupLABELS(particle, 0.060, true);
	cMC->RedrawAxis();

	TLegend* leg = new TLegend(0.6, 0.78, 0.92, 0.90, NULL, "brNDC");
	leg->SetBorderSize(0);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(frame->findObject(dataName), "Signal MC", "lp");
	leg->AddEntry(frame->findObject(sigName), "Signal PDF", "f");
	leg->Draw();
}

// The two pads of a data fit (roofitB.h): the fit above, the pulls below. Returns the frame of the fit, P1 current.
RooPlot* MakeFitFrame(TCanvas* c, const TreeSetup& setup, RooRealVar* mass, int nMassBins, TPad*& p1, TPad*& p2)
{
	c->Clear();
	c->cd();
	RooPlot* frame = mass->frame(Range(setup.plotMin, setup.plotMax));
	frame->SetStats(0);
	frame->SetTitle("");
	frame->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})", (mass->getMax() - mass->getMin()) / nMassBins * 1000));
	frame->GetYaxis()->SetTitleOffset(2.);
	frame->GetYaxis()->SetTitleSize(0.035);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetLabelSize(0.035);
	frame->GetXaxis()->SetTitleOffset(1.2);
	frame->GetXaxis()->SetTitleSize(0.);
	frame->GetXaxis()->SetTitleFont(0);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelSize(0);
	p1 = new TPad("p1", "p1", 0., 0.22, 1., 1);
	p1->SetBorderMode(1);
	p1->SetFrameBorderMode(0);
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.01);
	p1->SetLeftMargin(0.14);
	p1->SetRightMargin(0.04);
	p1->Draw();
	p2 = new TPad("p2", "p2", 0., 0., 1., 0.22);
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
	return frame;
}

// The pulls of the data DATANAME against the curve CURVENAME of FRAME, in the lower pad P2.
void DrawPulls(TPad* p2, RooPlot* frame, const TreeSetup& setup, RooRealVar* mass, const TString& dataName, const TString& curveName)
{
	p2->cd();
	RooHist* pulls = frame->pullHist(dataName, curveName);
	pulls->SetMarkerSize(0.5);
	RooPlot* pullFrame = mass->frame(Range(setup.plotMin, setup.plotMax));
	pullFrame->addPlotable(static_cast<RooPlotable*>(pulls), "XP");
	pullFrame->SetTitle("");
	pullFrame->SetXTitle(setup.massTitle);
	pullFrame->SetYTitle("Pull");
	pullFrame->GetYaxis()->SetTitleFont(42);
	pullFrame->GetYaxis()->SetTitleSize(0.13);
	pullFrame->GetYaxis()->CenterTitle(kTRUE);
	pullFrame->GetYaxis()->SetLabelOffset(0.01);
	pullFrame->GetYaxis()->SetLabelFont(42);
	pullFrame->GetYaxis()->SetLabelSize(0.12);
	pullFrame->GetYaxis()->SetNdivisions(305);
	pullFrame->GetYaxis()->SetTitleOffset(0.55);
	pullFrame->GetYaxis()->SetRangeUser(-3.5, 3.5);
	pullFrame->GetXaxis()->SetTitleSize(0.13);
	pullFrame->GetXaxis()->SetTitleOffset(1.0);
	pullFrame->GetXaxis()->CenterTitle();
	pullFrame->GetXaxis()->SetLabelFont(42);
	pullFrame->GetXaxis()->SetLabelOffset(0.01);
	pullFrame->GetXaxis()->SetLabelSize(0.14);
	pullFrame->GetXaxis()->SetTickLength(0.16);
	pullFrame->Draw();
	TLine* zero = new TLine(setup.plotMin, 0., setup.plotMax, 0.);
	zero->SetLineColor(2);
	zero->SetLineWidth(1);
	zero->Draw("same");
}

// Draws one data fit as roofitB.h draws it: data, total model, signal (filled), combinatorial background,
// partially reconstructed B (ntKp, when ERFC is not empty), the fitted parameters and the pulls. Leaves the upper
// pad current and returns the chi2/ndf of the total model with NFLOAT parameters.
double DrawDataFit(TCanvas* c, const TString& tree, RooRealVar* mass, RooDataSet* data, RooAbsPdf* model,
                   const RooArgSet& sig, const RooArgSet& bkg, const RooArgSet& erfc, int nFloat, int nMassBins)
{
	static int nDrawn = 0;
	++nDrawn;
	const TreeSetup& setup = GetTreeSetup(tree);
	const int signalColor = (tree == "ntmix_PSI2S") ? kOrange - 2 : kOrange - 3;
	TPad *p1, *p2;
	RooPlot* frame = MakeFitFrame(c, setup, mass, nMassBins, p1, p2);

	const TString dataName = Form("ds_cut_%d", nDrawn), modelName = Form("model_%d", nDrawn);
	const TString sigName = Form("sig_%d", nDrawn), bkgName = Form("bkg_%d", nDrawn), erfcName = Form("erfc_%d", nDrawn);
	data->plotOn(frame, Name(dataName), Binning(nMassBins), MarkerSize(0.5), MarkerStyle(8), MarkerColor(1), LineColor(1), LineWidth(1));
	model->plotOn(frame, Name(modelName), Range("all"), NormRange("all"), Precision(1e-6), DrawOption("L"), LineColor(2), LineWidth(1));
	model->plotOn(frame, Name(sigName), Components(sig), DrawOption("LF"), FillStyle(3002), FillColor(signalColor), LineStyle(7),
	              LineColor(signalColor), LineWidth(1));
	if (erfc.getSize() > 0) model->plotOn(frame, Name(erfcName), Components(erfc), Range("all"), NormRange("all"), LineColor(kGreen + 3),
	                                      LineStyle(9), LineWidth(2), DrawOption("L"));
	model->plotOn(frame, Name(bkgName), Components(bkg), Range("all"), NormRange("all"), Precision(1e-6), DrawOption("L"), LineStyle(7),
	              LineColor(4), LineWidth(1));
	double chi2 = frame->chiSquare(modelName, dataName, nFloat);
	if (!std::isfinite(chi2) || chi2 < 0) chi2 = -1.0;
	model->paramOn(frame, Layout(0.18, 0.48, 0.78), Format("NEU", AutoPrecision(2)));
	frame->getAttFill()->SetFillStyle(0);  // the box of the parameters, the last item of the frame
	frame->getAttLine()->SetLineWidth(0);
	frame->Draw();

	TLegend* leg = new TLegend(0.67, 0.60, 0.91, 0.90, NULL, "brNDC");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);
	leg->AddEntry(frame->findObject(dataName), "Data", "LEP");
	leg->AddEntry(frame->findObject(modelName), "Fit Model", "l");
	leg->AddEntry(frame->findObject(bkgName), " Comb. Bkg.", "l");
	leg->AddEntry(frame->findObject(sigName), setup.decay, "f");
	if (erfc.getSize() > 0) leg->AddEntry(frame->findObject(erfcName), " B #rightarrow J/#psi X", "l");
	leg->Draw();

	DrawPulls(p2, frame, setup, mass, dataName, modelName);
	p1->cd();
	return chi2;
}

// The add-up plot of a combined fit. No fit is made to the summed data: the fit is simultaneous, and each year's
// model is normalized to the yields of that year. The plot shows the data of all years summed, the combined fit (the
// sum of the yearly fits, red, drawn above everything), the signal of every year (dashed, PbPb23 orange, PbPb24
// yellowish, ...) and the background of all years (combinatorial and partially reconstructed B, dashed blue), with the
// line styles of the plot of each year, and the pulls of the summed data against the combined fit. Leaves the upper pad current and returns the chi2/ndf of the combined fit with NFLOAT
// parameters (the floating parameters of the combined fit).
double DrawAddUpFit(TCanvas* c, const TString& tree, RooRealVar* mass, RooDataSet* data, const std::vector<RooAbsPdf*>& yearModels,
                    const std::vector<RooAbsPdf*>& yearSignals, const RooArgSet& backgrounds, const std::vector<TString>& years,
                    int nFloat, int nMassBins)
{
	static int nDrawn = 0;
	++nDrawn;
	const TreeSetup& setup = GetTreeSetup(tree);
	const int signalColors[] = {kOrange + 7, kOrange - 2, kRed - 7, kMagenta - 7};
	TPad *p1, *p2;
	RooPlot* frame = MakeFitFrame(c, setup, mass, nMassBins, p1, p2);

	RooArgList models;
	for (RooAbsPdf* m : yearModels) models.add(*m);
	RooAddPdf* sumModel = new RooAddPdf(Form("sumOfYearlyFits_%d", nDrawn), "", models);
	const TString dataName = Form("ds_sum_%d", nDrawn), sumName = Form("sum_%d", nDrawn), bkgName = Form("bkg_sum_%d", nDrawn);
	data->plotOn(frame, Name(dataName), Binning(nMassBins), MarkerSize(0.5), MarkerStyle(8), MarkerColor(1), LineColor(1), LineWidth(1));
	sumModel->plotOn(frame, Name(bkgName), Components(backgrounds), Range("all"), NormRange("all"), Precision(1e-6), DrawOption("L"),
	                 LineStyle(7), LineColor(4), LineWidth(2));
	for (size_t k = 0; k < years.size(); ++k) {
		sumModel->plotOn(frame, Name(Form("sig_%d_%zu", nDrawn, k)), Components(*yearSignals[k]), Range("all"), NormRange("all"),
		                 Precision(1e-6), DrawOption("L"), LineColor(signalColors[k % 4]), LineStyle(7), LineWidth(2));
	}
	sumModel->plotOn(frame, Name(sumName), Range("all"), NormRange("all"), Precision(1e-6), DrawOption("L"), LineColor(2), LineWidth(2));
	double chi2 = frame->chiSquare(sumName, dataName, nFloat);
	if (!std::isfinite(chi2) || chi2 < 0) chi2 = -1.0;
	frame->Draw();

	TLegend* leg = new TLegend(0.60, 0.75 - 0.05 * years.size(), 0.91, 0.90, NULL, "brNDC");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);
	leg->AddEntry(frame->findObject(dataName), "Data, all years", "LEP");
	leg->AddEntry(frame->findObject(sumName), "Combined fit", "l");
	for (size_t k = 0; k < years.size(); ++k) leg->AddEntry(frame->findObject(Form("sig_%d_%zu", nDrawn, k)), years[k] + "' signal", "l");
	leg->AddEntry(frame->findObject(bkgName), "Bkg.", "l");
	leg->Draw();

	DrawPulls(p2, frame, setup, mass, dataName, sumName);
	p1->cd();
	return chi2;
}

// The labels roofitB.C draws on a data fit, in the current pad: particle, chi2/ndf, analysis bin, signal yield, the
// model variation (when not empty) and the significance (when not empty).
void DrawFitLabels(const TString& tree, const TString& var, double low, double high, double chi2, double yield,
                   double yieldError, const TString& variation, const TString& significance)
{
	setupLABELS(new TLatex(0.2, 0.85, FitParticleLabel(tree, true)), 0.060, true);
	setupLABELS(new TLatex(0.68, 0.55, Form("#chi^{2}/ndf = %.2f", chi2)));
	setupLABELS(new TLatex(0.68, 0.50, Form("%.1f #leq %s < %.1f", low, FitVarLabel(var).Data(), high)));
	setupLABELS(new TLatex(0.68, 0.45, Form("Y_{s} = %.0f #pm %.0f", yield, yieldError)));
	if (variation != "") setupLABELS(new TLatex(0.68, 0.40, Form("(%s)", variation.Data())), 0.030);
	if (significance != "") setupLABELS(new TLatex(0.68, 0.35, significance), 0.035);
}

// One simultaneous fit of one analysis bin: the signal shape from the MC of all years (one simultaneous MC fit, one
// shape and per-year normalizations), then the data of all years with the shape fixed and everything but the yields
// common to all years.
struct CombinedFit {
	std::vector<YearModel> years;
	RooSimultaneous* sim = nullptr;
	RooSimultaneous* simMC = nullptr;
	RooFitResult* result = nullptr;
	RooFitResult* resultMC = nullptr;
	RooArgSet* mcValues = nullptr;    // every parameter right after the MC fit
	RooArgSet* dataValues = nullptr;  // every parameter right after the data fit
};

CombinedFit FitCombined(const TString& tree, int bin, const TString& sigCode, const TString& bkgCode,
                        const std::vector<TString>& years, RooRealVar* mass,
                        RooCategory* cat, const std::vector<RooDataSet*>& data, const std::vector<RooDataSet*>& mc,
                        RooDataSet* combData, RooDataSet* combMC)
{
	CombinedFit f;
	std::map<TString, RooRealVar*> pool;
	const TString code = sigCode + bkgCode;
	f.sim = new RooSimultaneous(Form("simPdf%d_%s", bin, code.Data()), "", *cat);
	f.simMC = new RooSimultaneous(Form("simPdfMC%d_%s", bin, code.Data()), "", *cat);
	RooArgSet allPars;
	for (size_t y = 0; y < years.size(); ++y) {
		f.years.push_back(BuildYearModel(tree, bin, sigCode, bkgCode, years[y], pool, mass, data[y], mc[y]));
		f.sim->addPdf(*f.years[y].model, years[y]);
		f.simMC->addPdf(*f.years[y].modelMC, years[y]);
		f.years[y].scale->setConstant(true);
		for (RooRealVar* p : f.years[y].pars) allPars.add(*p, true);
	}

	f.resultMC = f.simMC->fitTo(*combMC, Save(), Extended(), Range("signal"), SumW2Error(true));
	f.mcValues = static_cast<RooArgSet*>(allPars.snapshot());

	for (YearModel& m : f.years) {
		for (RooRealVar* p : m.shapePars) p->setConstant(true);
		if (sigCode == "fixed") m.mean->setConstant(true);
		m.scale->setConstant(false);
	}
	f.result = f.sim->fitTo(*combData, Save(), Extended(true), Range("all"));
	f.result->Print("v");
	f.dataValues = static_cast<RooArgSet*>(allPars.snapshot());
	return f;
}

// One canvas with one graph per series over the analysis bins (X at the mean of the bin, horizontal bars to its
// edges), saved to PATH.
void DrawBinGraphs(const TString& path, const TString& xTitle, const TString& yTitle, const std::vector<double>& x,
                   const std::vector<double>& exLow, const std::vector<double>& exHigh, const std::vector<TString>& labels,
                   const std::vector<std::vector<double>>& y, const std::vector<std::vector<double>>& ey)
{
	const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2};
	const int markers[] = {20, 21, 22, 23, 33, 34, 29};
	TCanvas c("cBinGraphs", "", 700, 600);
	c.SetLeftMargin(0.14);
	TMultiGraph* mg = new TMultiGraph();
	TLegend* leg = new TLegend(0.60, 0.88 - 0.05 * y.size(), 0.92, 0.90);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextFont(42);
	leg->SetTextSize(0.035);
	double yMax = 0.0;
	for (size_t s = 0; s < y.size(); ++s) {
		TGraphAsymmErrors* g = new TGraphAsymmErrors(x.size(), x.data(), y[s].data(), exLow.data(), exHigh.data(), ey[s].data(), ey[s].data());
		g->SetLineColor(colors[s % 7]);
		g->SetMarkerColor(colors[s % 7]);
		g->SetMarkerStyle(markers[s % 7]);
		g->SetLineWidth(2);
		mg->Add(g, "PE");
		leg->AddEntry(g, labels[s], "lep");
		for (size_t i = 0; i < x.size(); ++i) yMax = std::max(yMax, y[s][i] + ey[s][i]);
	}
	mg->SetTitle(";" + xTitle + ";" + yTitle);
	mg->Draw("A");
	mg->GetYaxis()->SetRangeUser(0.0, yMax > 0.0 ? 1.4 * yMax : 1.0);
	mg->GetYaxis()->SetTitleOffset(1.6);
	leg->Draw();
	c.SaveAs(path);
}
