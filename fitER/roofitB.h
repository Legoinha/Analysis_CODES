#pragma once

#include "aux/uti.h"
#include "RooWorkspace.h"
#include "TString.h"

using namespace RooFit;
using namespace std;

Int_t _count=0;

RooFitResult *fit(TString variation, TString pdf,TString tree, TCanvas* c, TCanvas* cMC, RooDataSet* ds, RooDataSet* dsMC, RooRealVar* mass, float binmin, float binmax, RooWorkspace& w, TString which_var, int NBIN){
	
	// set initial masses (means)
	RooDataSet* dsMC_spec = nullptr;
	double init_mean = Bs_MASS;
	double init_mean_PSI2S = PSI2S_MASS;
	if(tree == "ntKp")        {init_mean = Bu_MASS;}
    else if(tree == "ntKstar"){init_mean = Bd_MASS;}
    else if(tree == "ntmix")  {
		dsMC_spec = static_cast<RooDataSet*>(dsMC->reduce("isX3872==0"));
		dsMC = static_cast<RooDataSet*>(dsMC->reduce("isX3872==1"));
		init_mean = X3872_MASS;
	}

	// Signal-shape parameters (used for both MC and data fits)
	RooRealVar mean(Form("mean%d_%s",_count,pdf.Data()),"",init_mean,init_mean-0.01,init_mean+0.01) ;
	RooRealVar sigma1(Form("sigma1%d_%s",_count,pdf.Data()),"",0.01, 0.001, 0.1) ;
	RooRealVar sigma2(Form("sigma2%d_%s",_count,pdf.Data()),"",0.005, 0.001, 0.1) ;
	RooRealVar sigma3(Form("sigma3%d_%s",_count,pdf.Data()),"",0.01, 0.001, 0.03) ;
	RooRealVar sigma4cb(Form("sigma4cb%d_%s",_count,pdf.Data()),"",0.005,0.001,0.05) ;
	RooRealVar alpha(Form("alpha%d_%s",_count,pdf.Data()),"",4.,0,15);
	RooRealVar n(Form("n_%d_%s", _count, pdf.Data()),"",10,-100,200);
	RooRealVar* scale = new RooRealVar("scale","scale",1,0.2,2);
	RooProduct scaled_sigma1(Form("scaled_sigma1%d_%s",_count,pdf.Data()),"scaled_sigma1", RooArgList(*scale,sigma1));
	RooProduct scaled_sigma2(Form("scaled_sigma2%d_%s",_count,pdf.Data()),"scaled_sigma2", RooArgList(*scale,sigma2));
	RooProduct scaled_sigma3(Form("scaled_sigma3%d_%s",_count,pdf.Data()),"scaled_sigma3", RooArgList(*scale,sigma3));
	RooProduct scaled_sigmacb(Form("scaled_sigmacb%d_%s",_count,pdf.Data()),"scaled_sigmacb", RooArgList(*scale,sigma4cb));
	RooGaussian sig1(Form("sig1%d_%s",_count,pdf.Data()),"",*mass,mean,scaled_sigma1);  
	RooGaussian sig2(Form("sig2%d_%s",_count, pdf.Data()),"",*mass,mean,scaled_sigma2);  
	RooGaussian sig3(Form("sig3%d_%s",_count, pdf.Data()),"",*mass,mean,scaled_sigma3);  
	RooCBShape  CB(Form("CB%d_%s",_count, pdf.Data()),"",*mass,mean,scaled_sigmacb, alpha, n);
	RooRealVar sig1frac(Form("sig1frac%d_%s",_count, pdf.Data()),"", 0.5, 0.01, 1);
	RooRealVar sig2frac(Form("sig2frac%d_%s",_count, pdf.Data()),"", 0.5, 0.01, 1);  

	RooAddPdf* sig_doubleG = new RooAddPdf(Form("sig_doubleG%d_%s",_count,pdf.Data()),"",RooArgList(sig1,sig2),sig1frac);
	RooAddPdf* sig_Gaussian = new RooAddPdf(Form("sig_Gaussian%d_%s",_count,pdf.Data()), "", RooArgList(sig1), RooArgList(), true);
	RooAddPdf* sig_tripleG = new RooAddPdf(Form("sig_tripleG%d_%s",_count,pdf.Data()), "", RooArgList(sig1, sig2, sig3), RooArgList(sig1frac, sig2frac), true);
	RooAddPdf* sig_gaussCB = new RooAddPdf(Form("sig_gaussCB%d_%s",_count,pdf.Data()), "", RooArgList(sig1, CB), sig1frac);
	RooAddPdf* sig = nullptr;
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed" )) sig = sig_doubleG;
	if(variation=="signal" && pdf=="1gauss")   sig = sig_Gaussian;
	if(variation=="signal" && pdf=="3gauss")   sig = sig_tripleG;
	if(variation=="signal" && pdf=="gauss_cb") sig = sig_gaussCB;

	// Special cases (Psi2s or B0 wt) signal-shape parameters
	RooRealVar mean_spec(Form("mean_spec%d_%s",_count,pdf.Data()),"",init_mean_PSI2S,init_mean_PSI2S-0.01,init_mean_PSI2S+0.01) ;
	RooRealVar sigma1_spec(Form("sigma1_spec%d_%s",_count,pdf.Data()),"",0.005,0.0005,0.01);
	RooRealVar sigma2_spec(Form("sigma2_spec%d_%s",_count,pdf.Data()),"",0.001,0.0005,0.01);
	RooRealVar sigma3_spec(Form("sigma3_spec%d_%s",_count,pdf.Data()),"",0.005,0.0005,0.01);
	RooRealVar* scale_spec = new RooRealVar("scale_spec","scale_spec",1,0.1,2);
	RooProduct scaled_sigma1_spec(Form("scaled_sigma1_spec%d_%s",_count,pdf.Data()),"scaled_sigma1_spec", RooArgList(*scale_spec,sigma1_spec));
	RooProduct scaled_sigma2_spec(Form("scaled_sigma2_spec%d_%s",_count,pdf.Data()),"scaled_sigma2_spec", RooArgList(*scale_spec,sigma2_spec));
	RooProduct scaled_sigma3_spec(Form("scaled_sigma3_spec%d_%s",_count,pdf.Data()),"scaled_sigma3_spec", RooArgList(*scale_spec,sigma3_spec));
	RooGaussian sig1_spec(Form("sig1_spec%d_%s",_count,pdf.Data()),"",*mass,mean_spec,scaled_sigma1_spec);  
	RooGaussian sig2_spec(Form("sig2_spec%d_%s",_count,pdf.Data()),"",*mass,mean_spec,scaled_sigma2_spec);  
	RooGaussian sig3_spec(Form("sig3_spec%d_%s",_count,pdf.Data()),"",*mass,mean_spec,scaled_sigma3_spec);  
	RooRealVar sig1frac_spec(Form("sig1frac_spec%d_%s",_count,pdf.Data()),"", 0.5, 0.01, 1); 
	RooRealVar sig2frac_spec(Form("sig2frac_spec%d_%s",_count,pdf.Data()),"", 0.5, 0.01, 1);
	// Keep the same PDF concept as X(3872), but always build a dedicated psi(2S) shape.
	RooAddPdf* sig_spec_doubleG = new RooAddPdf(Form("sig_spec_doubleG%d_%s",_count,pdf.Data()),"",RooArgList(sig1_spec,sig2_spec),RooArgList(sig1frac_spec), true);
	RooAddPdf* sig_spec_Gaussian = new RooAddPdf(Form("sig_spec_Gaussian%d_%s",_count,pdf.Data()), "", RooArgList(sig1_spec), RooArgList(), true);
	RooAddPdf* sig_spec_tripleG = new RooAddPdf(Form("sig_spec_tripleG%d_%s",_count,pdf.Data()), "", RooArgList(sig1_spec, sig2_spec, sig3_spec), RooArgList(sig1frac_spec, sig2frac_spec), true);
	RooAddPdf* sig_spec = nullptr;
	if((variation=="" && pdf=="") || variation=="background" || (variation=="signal" && pdf=="fixed" )) sig_spec = sig_spec_doubleG;
	if(variation=="signal" && pdf=="1gauss")   sig_spec = sig_spec_Gaussian;
	if(variation=="signal" && pdf=="3gauss")   sig_spec = sig_spec_tripleG;

	// MC fit MODEL
	RooRealVar* nsigMC = nsigMC = new RooRealVar(Form("nsigMC%d_%s",_count, pdf.Data()),"", dsMC->sumEntries(), 0.9*dsMC->sumEntries(), 1.1 * dsMC->sumEntries());
	RooRealVar* nsig_specMC = nullptr;
	if(tree == "ntmix" ){nsig_specMC = new RooRealVar(Form("nsig_specMC%d_%s",_count,pdf.Data()), "", dsMC_spec->sumEntries(), 0.9*dsMC_spec->sumEntries(), 1.1 * dsMC_spec->sumEntries());}
	
	RooAddPdf* modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sig),RooArgList(*nsigMC));
	RooAddPdf* modelMC_spec = nullptr;
	if (tree == "ntmix") modelMC_spec = new RooAddPdf(Form("modelMC_spec%d_%s",_count, pdf.Data()),"",RooArgList(*sig_spec), RooArgList(*nsig_specMC));
	scale->setConstant();
	scale_spec->setConstant();
	// MC fit MODEL

	//////////ROOFIT ROOFIT ROOFIT MC MC MC MC MC
	// focus the MC fit to the signal region to prevent statistical flutuations
	mass->setRange("signal",init_mean-0.07, init_mean+0.07);  
	if(tree =="ntKp"){ mass->setRange("signal",init_mean-0.09, init_mean+0.09); }
	else if(tree=="ntmix"){ 
		mass->setRange("signal",init_mean-0.035, init_mean+0.035);
		mass->setRange("signal_psi2s",init_mean_PSI2S-0.03, init_mean_PSI2S+0.03);
	}
	
	RooFitResult * fitResultMC = modelMC->fitTo(*dsMC, Save(), Extended(), Range("signal"));
	w.import(*nsigMC);
	RooFitResult * fitResultMC_spec = nullptr;
	if (tree == "ntmix") {
		fitResultMC_spec = modelMC_spec->fitTo(*dsMC_spec, Save(), Range("signal_psi2s"));
		w.import(*nsig_specMC);
	}
	//////////ROOFIT ROOFIT ROOFIT MC MC MC MC MC

	// Get chi2/ndf for MC fits
	GetChi2Ndf(dsMC, modelMC, mass, NBIN, "signal", fitResultMC, &w, Form("chi2MC_norm%d_%s", _count, pdf.Data()));
	if (tree == "ntmix"){GetChi2Ndf(dsMC_spec, modelMC_spec, mass, NBIN, "signal_psi2s", fitResultMC_spec, &w, Form("chi2MC_PSI2S_norm%d_%s", _count, pdf.Data()));}
	// Get chi2/ndf for MC fits



	cMC->Clear();
	cMC->cd();
	TPad *pMC1 = new TPad(Form("pMC1_%d", _count),Form("pMC1_%d", _count), 0., 0.22, 1., 1);
	TPad *pMC2 = new TPad(Form("pMC2_%d", _count),Form("pMC2_%d", _count), 0., 0., 1., 0.22);
	pMC1->SetBorderMode(1);
	pMC1->SetFrameBorderMode(0);
	pMC1->SetBorderSize(2);
	pMC1->SetBottomMargin(0.01);
	pMC1->SetLeftMargin(0.14);
	pMC1->SetRightMargin(0.04);
	pMC2->SetTopMargin(0.0);
	pMC2->SetBottomMargin(0.34);
	pMC2->SetLeftMargin(0.14);
	pMC2->SetRightMargin(0.04);
	pMC2->SetBorderMode(0);
	pMC2->SetBorderSize(2);
	pMC2->SetFrameBorderMode(0);
	pMC2->SetTicks(1,1);
	pMC1->Draw();
	pMC2->Draw();
	double xMinPlot = 5.05;
	double xMaxPlot = 5.55;
	if(tree=="ntphi") { xMinPlot = 5.2; xMaxPlot = 5.55; }
	if(tree=="ntmix") { xMinPlot = 3.6;  xMaxPlot = 4.0; }
	RooPlot* frameMC = mass->frame(Range(xMinPlot, xMaxPlot));
	frameMC->SetTitle("");
	frameMC->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(mass->getMax()-mass->getMin())/NBIN*1000));    
	frameMC->GetYaxis()->SetTitleOffset(2.);
	frameMC->GetXaxis()->SetTitleOffset(1.2);
	frameMC->GetYaxis()->SetTitleSize(0.035);
	frameMC->GetXaxis()->SetTitleSize(0.);
	frameMC->GetXaxis()->SetTitleFont(0);
	frameMC->GetYaxis()->SetTitleFont(42);
	frameMC->GetXaxis()->SetLabelFont(42);
	frameMC->GetYaxis()->SetLabelFont(42);
	frameMC->GetXaxis()->SetLabelSize(0);
	frameMC->GetYaxis()->SetLabelSize(0.035);
	frameMC->SetStats(0);
	TString xTtile_decayC = "";
	if(tree == "ntmix"){        xTtile_decayC = "m_{J/#psi #pi^{+} #pi^{-}} [GeV/c^{2}]";}
	else if(tree == "ntphi"){   xTtile_decayC = "m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]";}
	else if (tree == "ntKstar"){xTtile_decayC = "m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]";}
	else if (tree == "ntKp" ){  xTtile_decayC = "m_{J/#psi K^{+}} [GeV/c^{2}]";}

	//for pull
	RooDataSet* dsMC_plot = dsMC;
	RooAbsPdf* modelMC_pull = modelMC;
	if (tree=="ntmix" && dsMC_spec && modelMC_spec) {
		dsMC_plot = static_cast<RooDataSet*>(dsMC->Clone(Form("dsMC_combined%d_%s", _count, pdf.Data())));
		dsMC_plot->append(*dsMC_spec);
		modelMC_pull = new RooAddPdf(Form("modelMC_combined%d_%s", _count, pdf.Data()), "", RooArgList(*sig, *sig_spec), RooArgList(*nsigMC, *nsig_specMC));
	}

	//PLOT MC FIT
	pMC1->cd();
	dsMC_plot->plotOn(frameMC, Name(Form("dsMC%d_%s", _count, pdf.Data())), Binning(NBIN), MarkerSize(0.5), MarkerStyle(8), LineColor(1), LineWidth(1));
	modelMC->plotOn(frameMC, Name(Form("sigMC%d_%s", _count, pdf.Data())), Range("signal"), NormRange("signal"), Normalization(nsigMC->getVal(), RooAbsReal::NumEvent), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(1));
	if (modelMC_spec){ modelMC_spec->plotOn(frameMC, Name(Form("sigSpecMC%d_%s", _count, pdf.Data())), Range("signal_psi2s"), NormRange("signal_psi2s"), Normalization(nsig_specMC->getVal(), RooAbsReal::NumEvent), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-2), LineStyle(7), LineColor(kOrange-2), LineWidth(1));}
	modelMC_pull->plotOn(frameMC, Name(Form("modelMCcurve%d_%s", _count, pdf.Data())), DrawOption("L"), LineWidth(0));

	// PARAMETERS
	modelMC->paramOn(frameMC,Layout(0.15,0.3,0.85), Format("NEU",AutoPrecision(2)));
	if (modelMC_spec) modelMC_spec->paramOn(frameMC, Layout(0.15,0.3,0.55), Format("NEU",AutoPrecision(2)));
	//frameMC->getAttText()->SetTextSize(0.03);
	frameMC->getAttFill()->SetFillStyle(0);
	frameMC->getAttLine()->SetLineWidth(0);
	frameMC->Draw();
	cMC->RedrawAxis();

	TLegend *legMC = new TLegend(0.6,0.78,0.92,0.90,NULL,"brNDC"); 
	if(tree == "ntmix") legMC = new TLegend(0.6,0.72,0.92,0.90,NULL,"brNDC");
	legMC->SetBorderSize(0);
	legMC->SetTextSize(0.035);     
	legMC->SetTextFont(42);
	legMC->SetFillStyle(0);
	legMC->AddEntry(frameMC->findObject(Form("dsMC%d_%s", _count, pdf.Data())), "Signal MC","lp");
	TString ident_par = "";
	if(tree=="ntphi"){    ident_par = "B^{0}_{s}";}
	else if(tree=="ntKp"){ident_par = "B^{+}";}
	else if(tree =="ntKstar"){
		ident_par = "B^{0}";
		//legMC->AddEntry(frameMC->findObject(sigSpecMC->GetName()),"WT PDF","f");
	}
	else if(tree=="ntmix"){
		ident_par = "X(3872)";
		legMC->AddEntry(frameMC->findObject(Form("sigSpecMC%d_%s", _count, pdf.Data())), Form("#psi(2S), %s ", pdf.Data()),"f");
	}
	legMC->AddEntry(frameMC->findObject(Form("sigMC%d_%s", _count, pdf.Data())),Form("%s, %s ",ident_par.Data(),pdf.Data()),"f");
  	legMC->Draw();

	// Draw pulls in bottom pad
	pMC2->cd();
	RooHist* pull_hist_mc = frameMC->pullHist(Form("dsMC%d_%s", _count, pdf.Data()), Form("modelMCcurve%d_%s", _count, pdf.Data()));
	pull_hist_mc->SetMarkerSize(0.6);
	pull_hist_mc->SetMarkerStyle(20);
	pull_hist_mc->SetMarkerColor(1);
	pull_hist_mc->SetLineColor(1);
	if (tree == "ntmix" && pull_hist_mc) {
		for (int ip = 0; ip < pull_hist_mc->GetN(); ++ip) {
			double x = 0., y = 0.;
			pull_hist_mc->GetPoint(ip, x, y);
			bool inXsignal = (x >= mass->getMin("signal")  && x <= mass->getMax("signal"));
			bool inPsi2s   = (x >= mass->getMin("signal_psi2s") && x <= mass->getMax("signal_psi2s"));
			if (!inXsignal && !inPsi2s) {
				pull_hist_mc->SetPoint(ip, x, 1e6);
				pull_hist_mc->SetPointError(ip, 0., 0., 0., 0.);
			}
		}
	}

	if (pull_hist_mc) {
		RooPlot* pull_plot_mc = mass->frame(Range(xMinPlot, xMaxPlot));
		pull_plot_mc->GetXaxis()->SetRangeUser(xMinPlot, xMaxPlot);
		pull_plot_mc->addPlotable(static_cast<RooPlotable*>(pull_hist_mc),"XP");
		pull_plot_mc->SetTitle("");
		pull_plot_mc->SetXTitle(xTtile_decayC);
		pull_plot_mc->SetYTitle("Pull");
		pull_plot_mc->GetYaxis()->SetTitleFont(42);
		pull_plot_mc->GetYaxis()->SetTitleSize(0.11);
		pull_plot_mc->GetYaxis()->CenterTitle(kFALSE);
		pull_plot_mc->GetYaxis()->SetLabelOffset(0.01);
		pull_plot_mc->GetYaxis()->SetLabelFont(42);
		pull_plot_mc->GetYaxis()->SetLabelSize(0.10);
		pull_plot_mc->GetYaxis()->SetNdivisions(305);
		pull_plot_mc->GetYaxis()->SetTitleOffset(0.55);
		pull_plot_mc->GetYaxis()->SetRangeUser(-3.5, 3.5);
		pull_plot_mc->GetXaxis()->SetTitleSize(0.13);
		pull_plot_mc->GetXaxis()->SetTitleOffset(1.0);
		pull_plot_mc->GetXaxis()->CenterTitle();
		pull_plot_mc->GetXaxis()->SetLabelFont(42);
		pull_plot_mc->GetXaxis()->SetLabelOffset(0.01);
		pull_plot_mc->GetXaxis()->SetLabelSize(0.14);
		pull_plot_mc->GetXaxis()->SetTickLength(0.16);
		pull_plot_mc->Draw();

		// Draw horizontal reference line in pull plot
		TLine *line_ref_mc_bg = new TLine(xMinPlot, 0., xMaxPlot, 0.);
		line_ref_mc_bg->SetLineStyle(2);
		line_ref_mc_bg->SetLineColor(kGray+2);
		line_ref_mc_bg->SetLineWidth(1);
		line_ref_mc_bg->Draw("same");
		TLine *line_ref_mc_sig = new TLine(mass->getMin("signal"), 0., mass->getMax("signal"), 0.);
		line_ref_mc_sig->SetLineStyle(1);
		line_ref_mc_sig->SetLineColor(kOrange-3);
		line_ref_mc_sig->SetLineWidth(2);
		line_ref_mc_sig->Draw("same");
		if (tree == "ntmix") {
			TLine *line_ref_mc_spec = new TLine(mass->getMin("signal_psi2s"), 0., mass->getMax("signal_psi2s"), 0.);
			line_ref_mc_spec->SetLineStyle(1);
			line_ref_mc_spec->SetLineColor(kOrange-2);
			line_ref_mc_spec->SetLineWidth(2);
			line_ref_mc_spec->Draw("same");
		}
	}





	/////////////////////////////////////////////////////////////////////	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////	/////////////////////////////////////////////////////////////////////

	//FIT THE DATA FIT THE DATA FIT THE DATA

	////////// SIGNAL FUNCTIONS
	double n_signal_initial = ds->sumEntries(TString::Format("abs(Bmass-%g)<0.05",init_mean)) ;
	double n_signal_initial_spec=0;
	if (tree == "ntmix"){
		n_signal_initial = ds->sumEntries(TString::Format("abs(Bmass-%g)<0.005",init_mean)) ;
		n_signal_initial_spec = ds->sumEntries(TString::Format("abs(Bmass-%g)<0.006",init_mean_PSI2S));
	}
	RooRealVar nsig(Form("nsig%d_%s",_count,pdf.Data()),"",n_signal_initial*0.4,n_signal_initial*0.2,n_signal_initial*1.5);
	RooRealVar c1(Form("c1%d_%s",_count,pdf.Data()),"",1.,0.,5.);

    // special cases 
	// PSI2S 
	RooRealVar* nsig_spec = nullptr;
    if(tree == "ntmix" ){
		nsig_spec = new RooRealVar(Form("nsig_spec%d_%s",_count,pdf.Data()),"",n_signal_initial_spec*0.3,n_signal_initial_spec*0.05,n_signal_initial_spec*1.2);
		sigma1_spec.setConstant();
		if(pdf!="1gauss"){
			sigma2_spec.setConstant();
			sig1frac_spec.setConstant();
		}
		if(variation=="signal" && pdf=="3gauss"){
			sigma3_spec.setConstant();
			sig2frac_spec.setConstant();
		}
		if(variation=="signal" && pdf=="fixed") mean_spec.setConstant();
    }
	////////// SIGNAL FUNCTIONS

	//////////  BACKGROUND FUNCTIONS
	RooRealVar nbkg(Form("nbkg%d_%s",_count,pdf.Data()),"",ds->sumEntries() * 0.7,ds->sumEntries()*0.1,ds->sumEntries());
	RooRealVar a0(Form("a0%d_%s",_count,pdf.Data()),"",-0.35,-2,2);
	RooRealVar a1(Form("a1%d_%s",_count,pdf.Data()),"",-0.05,-2,2);
	RooRealVar a2(Form("a2%d_%s",_count,pdf.Data()),"",0.01,-2,2);
	RooRealVar a3(Form("a3%d_%s",_count,pdf.Data()),"",0,-2,2);
	// // //
	//RooRealVar a0(Form("a0%d_%s",_count,pdf.Data()),"",-0.77,-1,1); 
	//RooRealVar a1(Form("a1%d_%s",_count,pdf.Data()),"",-0.15,-1,1);
	//RooRealVar a2(Form("a2%d_%s",_count,pdf.Data()),"",0.23,-2,1);
	//RooRealVar a3(Form("a3%d_%s",_count,pdf.Data()),"",-0.1,-1,1);
	// // //
 	RooChebychev bkg_2nd(Form("bkg%d_%s",_count,pdf.Data()), "", *mass, RooArgList(a0,a1));
	RooChebychev bkg_3rd(Form("bkg%d_%s",_count,pdf.Data()), "", *mass,RooArgSet(a0,a1,a2));
	RooChebychev bkg_4th(Form("bkg%d_%s",_count,pdf.Data()), "", *mass,RooArgSet(a0,a1,a2,a3));
	RooRealVar lambda(Form("lambda%d_%s", _count,pdf.Data()), "lambda",-0.5, -5., 0.1);
	RooExponential bkg(  Form("bkg%d_%s",_count,pdf.Data()),"",*mass,lambda);

	// special cases
    // B+ Part. Reconstructed Background 
	RooRealVar nbkg_part_r(Form("nbkg_part_r%d_%s",_count,pdf.Data()),"",6000,50,1e4);
	RooRealVar* m_nonprompt_scale=0;
	RooRealVar* m_nonprompt_shift=0;
	m_nonprompt_scale = new RooRealVar(Form("m_nonprompt_scale%d_%s",_count,""), "m_nonprompt_scale",0.01, 0.001, 0.1);
	m_nonprompt_shift = new RooRealVar(Form("m_nonprompt_shift%d_%s",_count,""), "m_nonprompt_shift", 5.15, 5.1, 5.2);
	RooGenericPdf* erfc = new RooGenericPdf( Form("erfc%d",_count), "0.5*TMath::Erfc((@0-@2)/@1)", RooArgList(*mass, *m_nonprompt_scale, *m_nonprompt_shift));

	//////////  BACKGROUND FUNCTIONS

	//////////////// MODELs MODELs MODELs MODELs MODELs MODELs MODELs MODELs MODELs 
	RooAddPdf* model = nullptr;

	/////////////////X X X X X X X X X X X     nominal: DG + 3rO
	if(tree == "ntmix"){
		if((variation=="" && pdf=="")|| variation=="signal")  model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig_spec,*sig,bkg_3rd),RooArgList(*nsig_spec,nsig,nbkg));
		if(variation=="background" && pdf=="4th") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig_spec,*sig,bkg_4th),RooArgList(*nsig_spec,nsig,nbkg));
	}
	/////////////////X X X X X X X X X X X

	/////////////////Bs Bs Bs Bs Bs Bs Bs Bs     nominal: DG + 3rO (for now)
	if(tree == "ntphi"){
		if((variation=="" && pdf=="") || variation=="signal") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_3rd),RooArgList(nsig,nbkg));
		if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_2nd),RooArgList(nsig,nbkg));
		if(variation=="background" && pdf=="4th") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_4th),RooArgList(nsig,nbkg));
	}
	/////////////////Bs Bs Bs Bs Bs Bs Bs Bs

	/////////////////B0 B0 B0 B0 B0 B0 B0 B0     nominal: DG + 3rO       // MISSING WT component fit!
	if(tree == "ntKstar"){
		if((variation=="" && pdf=="") || variation=="signal") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_3rd),RooArgList(nsig,nbkg));
		if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_2nd),RooArgList(nsig,nbkg));
		if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_3rd),RooArgList(nsig,nbkg));
}
	/////////////////B0 B0 B0 B0 B0 B0 B0 B0 

	/////////////////BP BP BP BP BP BP BP BP     nominal: DG + exp (+ Erfc)
	if(tree == "ntKp"){
		if((variation=="" && pdf=="")) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg,*erfc),RooArgList(nsig,nbkg,nbkg_part_r));
		if(pdf=="mass_range"){         model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg)      ,RooArgList(nsig,nbkg));}
		if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg_2nd,*sig,*erfc),RooArgList(nbkg,nsig,nbkg_part_r));
		if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg_3rd,*sig,*erfc),RooArgList(nbkg,nsig,nbkg_part_r));
		if(variation=="signal" && pdf=="1gauss" ) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg,sig1,*erfc),RooArgList(nbkg,nsig,nbkg_part_r));
		if((variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb" ))) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig, bkg, *erfc),RooArgList(nsig, nbkg, nbkg_part_r));
	}
	/////////////////BP BP BP BP BP BP BP BP

	//////////////// MODELs MODELs MODELs MODELs MODELs MODELs MODELs MODELs 

	//////////////// set MC learned parameters CONSTANT 
	scale->setConstant(false);	
	scale_spec->setConstant(false);	
	sigma1.setConstant();
	if(pdf!="1gauss"){
		sigma2.setConstant();
		sig1frac.setConstant();
	}
	if(variation=="signal" && pdf=="3gauss"){
		sigma3.setConstant();
		sig2frac.setConstant();  
	}
	if(variation=="signal" && pdf=="gauss_cb"){
		sigma4cb.setConstant();
		n.setConstant();
		alpha.setConstant();
	}
	if(variation=="signal" && pdf=="fixed") mean.setConstant();
	//////////////// SET PARAMETERS FROM MC FITS

	////// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT

	TString fitRange = (pdf == "mass_range") ? "m_rangeB" : "all";
	RooFitResult* fitResult = model->fitTo(*ds, Save(), Extended(kTRUE), Range(fitRange));
	fitResult->Print("v"); 
	w.import(*model);
	
	// Get chi2/ndf for data fit
	TString chi2DataVarName = Form("chi2_data_norm%d_%s", _count, pdf.Data());
	GetChi2Ndf(ds, model, mass, NBIN, "", fitResult, &w, chi2DataVarName.Data());
	// Get chi2/ndf for data fits

	////// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT









	///// PLOT DATA FIT PLOT DATA FIT PLOT DATA FIT PLOT DATA FIT PLOT DATA FIT PLOT DATA FIT 
	c->cd();
	RooPlot* frame = mass->frame(Range(xMinPlot, xMaxPlot));
	//model->paramOn(frame,Layout(0.2, 0.45, 0.5), Format("NEU",AutoPrecision(2)));
	//frame->getAttText()->SetTextSize(0.035);
	//frame->getAttFill()->SetFillStyle(0);
	//frame->getAttLine()->SetLineWidth(0);
	frame->SetStats(0);
	frame->SetTitle("");
	frame->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(mass->getMax()-mass->getMin())/NBIN*1000));
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
	TPad *p1 = new TPad("p1","p1",0.,0.22,1.,1);
	p1->SetBorderMode(1); 
	p1->SetFrameBorderMode(0); 
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.01);
	p1->SetLeftMargin(0.14);
	p1->SetRightMargin(0.04);
	p1->Draw();
	TPad *p2 = new TPad("p2","p2",0.,0.,1.,0.22); 
	p2->SetTopMargin(0.0);
	p2->SetBottomMargin(0.34);
	p2->SetLeftMargin(0.14);
	p2->SetRightMargin(0.04);
	p2->SetBorderMode(0);
	p2->SetBorderSize(2); 
	p2->SetFrameBorderMode(0); 
	p2->SetTicks(1,1); 
	p2->Draw();

	p1->cd();
	ds->plotOn(frame, Name(Form("ds_cut%d", _count)), Binning(NBIN), MarkerSize(0.5), MarkerStyle(8), MarkerColor(1), LineColor(1), LineWidth(1)); 
	model->plotOn(frame, Name(Form("model%d_%s", _count, pdf.Data())), Range(fitRange), NormRange(fitRange), Precision(1e-6), DrawOption("L"), LineColor(2), LineWidth(1));
	model->plotOn(frame, Name(Form("sig%d_%s", _count, pdf.Data())),  Components(*sig), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(1)); 
	if(tree== "ntmix"){model->plotOn(frame, Name(Form("SIG_spec%d_%s", _count, pdf.Data())), Components(*sig_spec), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-2), LineStyle(7), LineColor(kOrange-2), LineWidth(1)); }
	else if(tree== "ntKp"){ model->plotOn(frame, RooFit::Name(Form("erfc%d_%s",_count,"")) , Components(*erfc), Range(fitRange), NormRange(fitRange), LineColor(kGreen+3), LineStyle(9), LineWidth(2), DrawOption("L"));}
	model->plotOn(frame, Name(Form("bkg%d_%s",_count,pdf.Data())) ,  Components(bkg), Range(fitRange), NormRange(fitRange), Precision(1e-6),  DrawOption("L"), LineStyle(7), LineColor(4), LineWidth(1));

	frame->Draw();

	TLegend *leg = new TLegend(0.68,0.60,0.92,0.90,NULL,"brNDC"); 
	if(tree != "ntKp" && tree != "ntmix") leg = new TLegend(0.68,0.66,0.92,0.90,NULL,"brNDC");
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);
	leg->AddEntry(frame->findObject(Form("ds_cut%d", _count)), Form("Data"),"LEP");
	leg->AddEntry(frame->findObject(Form("model%d_%s",_count,pdf.Data())),"Fit Model","l");
	leg->AddEntry(frame->findObject(Form("bkg%d_%s",_count,pdf.Data()))," Comb. Bkg.","l");
	if(tree== "ntKp"){
		leg->AddEntry(frame->findObject(Form("sig%d_%s",_count,pdf.Data()))," B^{+} #rightarrow J/#psi K^{+}","f");
		leg->AddEntry(frame->findObject(Form("erfc%d_%s",_count,pdf.Data()))," B #rightarrow J/#psi X","l");
	}
	else if(tree== "ntmix"){
		leg->AddEntry(frame->findObject(Form("SIG_spec%d_%s",_count,pdf.Data()))," #psi(2S) #rightarrow J/#psi #pi^{+} #pi^{-}","l");
		leg->AddEntry(frame->findObject(Form("sig%d_%s",_count,pdf.Data()))," X(3872) #rightarrow J/#psi #pi^{+} #pi^{-}","f");
	}
	else if(tree== "ntphi"){
		leg->AddEntry(frame->findObject(Form("sig%d_%s",_count,pdf.Data()))," B_{s}^{0} #rightarrow J/#psi K^{+}K^{-}","f");
	}
	else if(tree== "ntKstar"){
		leg->AddEntry(frame->findObject(Form("sig%d_%s",_count,pdf.Data()))," B^{0} #rightarrow J/#psi K^{*}","f");
	}
	leg->Draw();
	
	p2->cd();
	RooHist* pull_hist = frame->pullHist(Form("ds_cut%d",_count),Form("model%d_%s",_count,pdf.Data()));
	pull_hist->SetMarkerSize(0.5);
	RooPlot* pull_plot = mass->frame(Range(xMinPlot, xMaxPlot));
	pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"XP");
	pull_plot->SetTitle("");
	pull_plot->SetXTitle(xTtile_decayC.Data());
	pull_plot->SetYTitle("Pull");
	pull_plot->GetYaxis()->SetTitleFont(42);  
	pull_plot->GetYaxis()->SetTitleSize(0.11);
	pull_plot->GetYaxis()->CenterTitle(kFALSE);
	pull_plot->GetYaxis()->SetLabelOffset(0.01);
	pull_plot->GetYaxis()->SetLabelFont(42);
	pull_plot->GetYaxis()->SetLabelSize(0.10);
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
	// draw horizontal zero line on top of the pull plot
	TLine *line_ref = new TLine(xMinPlot, 0., xMaxPlot, 0.);
	line_ref->SetLineStyle(1);
	line_ref->SetLineColor(2);
	line_ref->SetLineWidth(1);
	line_ref->Draw("same");

	///// PLOT DATA FIT PLOT DATA FIT  PLOT DATA FIT  PLOT DATA FIT  PLOT DATA FIT  PLOT DATA FIT 

	cout << "\n-------------------------------------------------------------------------------------- \n" << endl;
	cout << "Signal Yield Y_s = " << nsig.getVal() << "     yield Error = " << nsig.getError() ;
	if(tree=="ntmix"){
	cout << "PSI(2S) Signal Yield Y_s = " << nsig_spec->getVal() << "     yield Error = " << nsig_spec->getError() ;
	}
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
