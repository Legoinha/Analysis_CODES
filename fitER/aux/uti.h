#pragma once

#include "TAxis.h"
#include "TSystem.h"
#include "TLine.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooHist.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include <RooBifurGauss.h>
#include <RooCmdArg.h>
#include <TLegend.h>
#include <TLatex.h>
#include <fstream>
#include <string>
#include <iomanip>
#include <TCanvas.h>
#include <TPad.h>
#include "RooMCStudy.h"
#include "../../plotER/aux/masses.h"
#include <vector>
#include <array>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <utility>
#include <cmath>
#include <TFile.h>
#include <TDirectoryFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TMathText.h>
#include <TBox.h>
#include <TCut.h>
#include <TColor.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TFitResult.h>
#include "TMultiGraph.h"
#include <stdio.h>

using namespace std;



inline double GetSignificance(
	RooWorkspace* ws,
	int count,
	RooRealVar* mass,
	double nSigma = 2.0)
{
	if (!ws || !mass) return -1.0;

	RooRealVar* nsigVar = ws->var(Form("nsig%d_%s", count, ""));
	RooRealVar* nbkgVar = ws->var(Form("nbkg%d_%s", count, ""));
	RooRealVar* scaleVar = ws->var("scale");
	RooRealVar* sigma1Var = ws->var(Form("sigma1%d_%s", count, ""));
	RooRealVar* sigma2Var = ws->var(Form("sigma2%d_%s", count, ""));
	RooRealVar* fracVar = ws->var(Form("sig1frac%d_%s", count, ""));

	if (!nsigVar || !nbkgVar || !scaleVar || !sigma1Var || !sigma2Var || !fracVar) return -1.0;

	const double signalYield = nsigVar->getVal();
	const double bkgTotalYield = nbkgVar->getVal();
	const double scale = scaleVar->getVal();
	const double sigma1 = sigma1Var->getVal();
	const double sigma2 = sigma2Var->getVal();
	const double fracSigma1 = fracVar->getVal();

	RooAbsPdf* bkgPdf = ws->pdf(Form("bkg%d_%s", count, ""));
	RooRealVar* meanVar = ws->var(Form("mean%d_%s", count, ""));
	if (!meanVar) return -1.0;

	const double s1 = scale * sigma1;
	const double s2 = scale * sigma2;
	const double sigmaEff = std::sqrt(fracSigma1 * s1 * s1 + (1.0 - fracSigma1) * s2 * s2);

	if (!bkgPdf || sigmaEff <= 0.0 || nSigma <= 0.0) return -1.0;

	const double sigLow = meanVar->getVal() - nSigma * sigmaEff;
	const double sigHigh = meanVar->getVal() + nSigma * sigmaEff;
	mass->setRange("rangeSIG", sigLow, sigHigh);

	auto* intFrac = bkgPdf->createIntegral(*mass, RooFit::NormSet(*mass), RooFit::Range("rangeSIG"));
	const double fracInSignal = intFrac ? intFrac->getVal() : 0.0;
	const double bkgInWindow = fracInSignal * bkgTotalYield;
	const double denom = signalYield + bkgInWindow;
	const double signif = (denom > 0.0) ? signalYield / std::sqrt(denom) : -1.0;

	cout << "Signal window: [" << sigLow << ", " << sigHigh << "] (±" << nSigma << "σeff)\n";
	cout << "Background in Sig. Region " << fracInSignal << " ====> " << bkgInWindow << "\n";
	cout << "Signal Yield: " << signalYield << "\n";
	cout << "Significance: " << signif << "\n";

	return signif;
}

inline double GetChi2Ndf(
	RooDataSet* ds,
	RooAbsPdf* pdf,
	RooRealVar* mass,
	int nBins,
	const char* rangeName = "",
	RooFitResult* fitResult = nullptr,
	RooWorkspace* ws = nullptr,
	const char* wsVarName = "")
{

	const bool hasRange = (rangeName && rangeName[0] != '\0');
	RooPlot* chi2Frame = mass->frame();
	RooDataSet* dsForChi2 = ds;
	if (hasRange) {dsForChi2 = static_cast<RooDataSet*>(ds->reduce(RooFit::CutRange(rangeName)));}

	dsForChi2->plotOn(chi2Frame, RooFit::Name("chi2_data_tmp"), RooFit::Binning(nBins));
	if (hasRange) { pdf->plotOn(chi2Frame, RooFit::Name("chi2_pdf_tmp"), RooFit::Range(rangeName), RooFit::NormRange(rangeName));}
	else { pdf->plotOn(chi2Frame, RooFit::Name("chi2_pdf_tmp"), RooFit::Range(""), RooFit::NormRange(""));}

	int nPar = 0;
	if (fitResult) {
		nPar = fitResult->floatParsFinal().getSize();
	} else {
		RooArgSet* pars = pdf->getParameters(*ds);
		if (pars) {
			RooAbsCollection* floatPars = pars->selectByAttrib("Constant", kFALSE);
			if (floatPars) {
				nPar = floatPars->getSize();
				delete floatPars;
			}
			delete pars;
		}
	}

	double chi2Ndf = chi2Frame->chiSquare("chi2_pdf_tmp", "chi2_data_tmp", nPar);
	if (!std::isfinite(chi2Ndf) || chi2Ndf < 0) chi2Ndf = -1.0;

	if (ws && wsVarName && wsVarName[0] != '\0') {
		RooRealVar chi2Var(wsVarName, "", chi2Ndf);
		ws->import(chi2Var);
	}

	delete chi2Frame;
	if (dsForChi2 != ds) delete dsForChi2;
	return chi2Ndf;
}

inline void setupLABELS(TLatex* latexTEXT, double tSize = 0.035, bool DrawText = true){
	latexTEXT->SetNDC();
	latexTEXT->SetTextAlign(13);
	latexTEXT->SetTextFont(42);
	latexTEXT->SetTextSize(tSize);
	latexTEXT->SetLineWidth(2);
	if (DrawText) latexTEXT->Draw();
}

inline void DrawCmsHeader(
	TPad* pad,
	TString COLsystem = "",
	const TString& leftText = "#bf{CMS} #it{Preliminary}",
	float textSize = 0.045,
	float yOffset = 0.30)
{
	if (!pad) return;
	TString rightText = "";
	if (COLsystem=="ppRef") rightText = "pp #sqrt{s}=5.36 TeV, (L=455.7 pb^{-1})" ;
	else if (COLsystem=="PbPb23") rightText = "PbPb #sqrt{s_{NN}}=5.36 TeV, (L=1.72 nb^{-1})" ;


	pad->cd();
	const float l = pad->GetLeftMargin();
	const float t = pad->GetTopMargin();
	const float r = pad->GetRightMargin();
	const float y = 1.f - t + yOffset * t;

	TLatex latex;
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack);
	latex.SetTextFont(42);
	latex.SetTextSize(textSize);
	latex.SetTextAlign(11);
	latex.DrawLatex(l+0.04, y, leftText);
	latex.SetTextAlign(31);
	latex.SetTextSize(textSize-0.015);
	latex.DrawLatex(1.f - r +0.06, y, rightText);
}







struct SystVariationConfig {
	std::string code;
	std::string label;
};

inline std::vector<SystVariationConfig> GetBackgroundSystematicModels(const TString& tree)
{
	if (tree == "ntphi")  {//3rd-order Chebyshev
		return {{"2nd", "2nd-order Chebyshev"}, {"4th", "4th-order Chebyshev"}};
	}
	if (tree == "ntKstar"){//3rd-order Chebyshev
		return {{"2nd", "2nd-order Chebyshev"}, {"3rd", "3rd-order Chebyshev"}};
	}
	if (tree == "ntKp")  {//Exponential
		return {{"2nd", "2nd-order Chebyshev + erfc"}, {"3rd", "3rd-order Chebyshev + erfc"}};
	}
	if (tree == "ntmix") {//3rd-order Chebyshev
		return {{"4th", "4th-order Chebyshev"}};
	}
	return {};
}

inline std::vector<SystVariationConfig> GetSignalSystematicModels(const TString& tree)
{

	if (tree == "ntphi")  {// Double Gaussian
		return {{"1gauss", "Gaussian"}, {"fixed", "Fixed mean"}};
	}
	if (tree == "ntKstar"){// Double Gaussian
		return {{"3gauss", "Triple Gaussian"}, {"gauss_cb", "Gaussian + Crystal Ball"}, {"fixed", "Fixed mean"}};
	}
	if (tree == "ntKp")  {// Double Gaussian
		return {{"3gauss", "Triple Gaussian"}, {"gauss_cb", "Gaussian + Crystal Ball"}, {"fixed", "Fixed mean"}};
	}
	if (tree == "ntmix") {// Double Gaussian
		return {{"1gauss", "Gaussian"}, {"3gauss", "Triple Gaussian"}};
	}
	return {};
}

inline std::string GetSystematicColumnLabel(const TString& var, double lowEdge, double highEdge)
{
	std::ostringstream clabel;
	if (var == "Bpt") { clabel << lowEdge << "$<p_T<$" << highEdge   ;} 
  else if (var == "By") { clabel << lowEdge << "$<|y|<$" << highEdge ;} 
  else if (var == "nSelectedChargedTracks") { clabel << lowEdge << "$<nTrks<$" << highEdge ;} 
  else if (var == "CentBin") { clabel << lowEdge << "$<Cent<$" << highEdge ;} 
  else { clabel << lowEdge << "<" << var.Data() << "<" << highEdge  ;}
	return clabel.str();
}

inline std::vector<std::string> GetGeneralSystematicLabels()
{
	return {"Background PDF", "Signal PDF", "Total Systematic"};
}

inline std::vector<std::vector<double> > BuildGeneralSystematicNumbers(
	const std::vector<std::vector<double> >& general_syst,
	const std::vector<std::vector<double> >& stat_error)
{
	(void)stat_error;
	return general_syst;
}









inline void latex_table_block(
	std::ostream& out,
	int n_col,
	int n_lin,
	const std::vector<std::string>& col_name,
	const std::vector<std::string>& labels,
	const std::vector<std::vector<double> >& numbers)
{
	std::string col = "c";
	for (int i=1; i<n_col; i++) col += "|c";

	out << "\\begin{center}" << std::endl;
	out << "\\small" << std::endl;
	out << "\\begin{tabular}{" + col + "}" << std::endl;
	out << "\\toprule" << std::endl;

	for (int c=0; c<n_col-1; c++) out << col_name[c] << " & ";
	out << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;

	for (int i=1; i<n_lin; i++) {
		out << labels[i-1] << " & ";
		for (int c=1; c<n_col-1; c++) out << numbers[c-1][i-1] << " \\% & ";
		out << numbers[n_col-2][i-1] << " \\% \\\\" << std::endl;
	}

	out << "\\bottomrule" << std::endl;
	out << "\\end{tabular}" << std::endl;
	out << "\\end{center}" << std::endl;
}

inline void latex_tables_document(
	std::string filename,
	const std::vector<int>& n_cols,
	const std::vector<int>& n_lines,
	const std::vector<std::vector<std::string> >& col_names,
	const std::vector<std::vector<std::string> >& labels,
	const std::vector<std::vector<std::vector<double> > >& numbers)
{
	std::ofstream file(filename + ".tex");
	std::ofstream file_check(filename + "_check.tex");
	file << std::fixed << std::setprecision(2);
	file_check << std::fixed << std::setprecision(2);

	file_check << "\\documentclass{article}" << std::endl;
	file_check << "\\usepackage{geometry}" << std::endl;
	file_check << "\\usepackage{booktabs}" << std::endl;
	file_check << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;
	file_check << "\\begin{document}" << std::endl;

	for (size_t i = 0; i < n_cols.size(); ++i) {
		latex_table_block(file, n_cols[i], n_lines[i], col_names[i], labels[i], numbers[i]);
		latex_table_block(file_check, n_cols[i], n_lines[i], col_names[i], labels[i], numbers[i]);
		if (i + 1 < n_cols.size()) {
			file << std::endl << std::endl;
			file_check << "\\vspace{0.5cm}" << std::endl;
		}
	}

	file_check << "\\end{document}" << std::endl;
	file.close();
	file_check.close();

	std::string pdfCmd = "pdflatex -interaction=nonstopmode -output-directory=./files " + filename + "_check.tex";
	(void)system(pdfCmd.c_str());
	std::remove((filename + "_check.aux").c_str());
	std::remove((filename + "_check.log").c_str());
}

inline void WriteSystematicsTablesDocument(
	const std::string& filename,
	const std::vector<std::string>& col_name_signal,
	const std::vector<std::string>& col_name_back,
	const std::vector<std::string>& col_name_general,
	const std::vector<std::string>& labels_signal,
	const std::vector<std::string>& labels_back,
	const std::vector<std::string>& labels_general,
	const std::vector<std::vector<double> >& signal_numbers,
	const std::vector<std::vector<double> >& back_numbers,
	const std::vector<std::vector<double> >& general_numbers)
{
	std::vector<int> table_n_cols = {
		static_cast<int>(col_name_signal.size()),
		static_cast<int>(col_name_back.size()),
		static_cast<int>(col_name_general.size())
	};
	std::vector<int> table_n_lines = {
		static_cast<int>(1 + labels_signal.size()),
		static_cast<int>(1 + labels_back.size()),
		static_cast<int>(1 + labels_general.size())
	};
	std::vector<std::vector<std::string> > table_col_names = {col_name_signal, col_name_back, col_name_general};
	std::vector<std::vector<std::string> > table_labels = {labels_signal, labels_back, labels_general};
	std::vector<std::vector<std::vector<double> > > table_numbers = {signal_numbers, back_numbers, general_numbers};

	latex_tables_document(filename, table_n_cols, table_n_lines, table_col_names, table_labels, table_numbers);
}
