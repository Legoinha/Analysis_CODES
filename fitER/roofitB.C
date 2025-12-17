#include "roofitB.h"
#include "aux/CMS_lumi.C"
#include "TSystem.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include <TMath.h>
#include <string>
#include <sstream>
#include <TGraph.h>
#include <stdio.h>

#include "../plotER/aux/parameters.h"  
#include "../plotER/aux/ACCSEL.h"

void read_samples(RooWorkspace& w, vector<TString> label, TString fName, TString treeName, TString sample, TString system="ppRef", TString DOselCUTS="1");
std::pair<int, std::vector<double>> defineBinning(const TString& var, const TString& tree, int full);

// PDF VARIATION FOR SYST STUDIES
int syst_study=0;

// VALIDATION STUDIES
int val=0;

void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString outputfile = "", TString OUTPLOTF = "", TString ExtraMCsample = "", TString SYSTEM = "ppRef"){

	//Setup the working area
	gSystem->mkdir(Form("./%s/validation",OUTPLOTF.Data()),true); 
	gSystem->mkdir("./results/tables",true); 
	gSystem->mkdir("./results/Graphs", true); 
	gSystem->mkdir(Form("%s", OUTPLOTF.Data()),true); 
	
	gStyle->SetTextSize(0.05);
	gStyle->SetTextFont(42);
	gStyle->SetPadRightMargin(cRightMargin);
	gStyle->SetPadLeftMargin(cLeftMargin);
	gStyle->SetPadTopMargin(cTopMargin);
	gStyle->SetPadBottomMargin(cBottomMargin);
	gStyle->SetPadBottomMargin(0.45);
	gStyle->SetTitleX(.0f);

	// BINING DEFINITION
	int _nBins;
	std::vector<double> _varBINS;
	std::tie(_nBins, _varBINS) = defineBinning(VAR, TREE, FULL);

	// PRINT WHAT IS ABOUT TO HAPPEN //
	cout << endl << endl;
	cout << "Variable "<< VAR << endl;
	cout << "number of bins: " << _nBins << endl;	
	cout << TREE << " (" << _nBins << ") " << " BINS: ";
	for(int t=0; t < _nBins+1 ;t++){cout <<"__"<<  _varBINS[t]<<"__";}
	cout << endl;
	//cout << endl << endl;
	// PRINT WHAT IS ABOUT TO HAPPEN //

	TString NEWcuts = Form("%s",CUT.Data()); 
	cout << "seldata= "<< NEWcuts << endl;
	cout << endl << endl;

	//DEFINE VARIABLES
	RooRealVar* mass = nullptr;
	if (TREE == "ntmix"){ minhisto = minhisto_X, maxhisto = maxhisto_X;} 
	else { minhisto = minhisto_B, maxhisto = maxhisto_B; }
	mass = new RooRealVar("Bmass", "Bmass", minhisto, maxhisto);
	mass->setRange("m_rangeB", 5.2 , 5.5);        //set a range to be used if pdf = mass_rangeB
	mass->setRange("all", minhisto, maxhisto);
	RooRealVar* pt    = new RooRealVar("Bpt","Bpt",0,200);
	RooRealVar* y     = new RooRealVar("By","By",-2.4, 2.4);
	RooRealVar* nMult = new RooRealVar("nMult","nMult",0,100);

	RooWorkspace* ws = new RooWorkspace("ws");
	ws->import(*mass);
	ws->import(*y);
	ws->import(*pt);
	ws->import(*nMult);

	//DATA and MC SAMPLES
	vector<TString>   ANA_vars = {"Bpt", "By"};
	read_samples(*ws, ANA_vars, INPUTDATA.Data(), TREE.Data(), "data", SYSTEM.Data(), NEWcuts);
	read_samples(*ws, ANA_vars, INPUTMC.Data()  , TREE.Data(),   "mc", SYSTEM.Data(), NEWcuts);

	RooDataSet* data = (RooDataSet*) ws->data("data");
	RooDataSet* mc   = (RooDataSet*) ws->data("mc");
	cout << "Total DATA entries: " << data->sumEntries() << "\n"; 
	cout << "Total   MC entries: " << mc  ->sumEntries() << "\n";

	RooDataSet* mc_spec = nullptr;
	if ("ntmix"==TREE){
		read_samples(*ws, ANA_vars, ExtraMCsample.Data(), TREE.Data(), "mc_spec", SYSTEM.Data(), NEWcuts);
		mc_spec = (RooDataSet*) ws->data("mc_spec");		
		//mc_spec = (RooDataSet*) mc_spec->reduce(NEWcuts);
		cout << "Total MC (PSI2S) entries: " << mc_spec->sumEntries() << "\n";
		mc->append(*mc_spec);
        cout << "Total MC (merged) entries: " << mc->sumEntries() << "\n";
	}
	//DATA and MC SAMPLES

	//MODELS for syst studies
	vector<string> background = {"2nd", "mass_range"};
	vector<string> signal     = {"3gauss", "fixed", "gauss_cb"};
	vector<vector<double>> background_syst;
	vector<vector<double>> signal_syst;
	vector<vector<double>> general_syst;
	vector<vector<double>> back_syst_rel_values;
	vector<vector<double>> sig_syst_rel_values;
	vector<vector<double>> stat_error;
	//MODELS for syst studies

	// FIT INFO
	double yield_vec[_nBins];
	double yield_Stat_unc[_nBins];
	double yield_vec_systerr_low[_nBins];
	double yield_vec_systerr_high[_nBins];
	double scale_vec[_nBins];
	double scale_vec_unc[_nBins];
	double resol_vec_unc[_nBins];
	double var_mean_av[_nBins];
	double hori_av_low[_nBins];
	double hori_av_high[_nBins];
	double chi2_vec[_nBins];
	double chi2MC_vec[_nBins];
	double chi2_vec_sig[signal.size()][_nBins];
	double chi2_vec_back[background.size()][_nBins];
	double chi2MC_vec_sig[signal.size()][_nBins];
	double chi2MC_vec_back[background.size()][_nBins];
	double resol_vec[_nBins];
	// FIT INFO

	TH1D* hPt = new TH1D("hPt","",_nBins,_varBINS.data());  

	//Fit the J/psi pi MC sample (shapes of J/psi pi peak is determined inclusively)
	if(TREE=="ntKp"){

		//PDF MODELS PDF MODELS PDF MODELS
		//inclusive MC jpsipi Model
		RooRealVar* m_jpsipi_fraction2 = 0;
		RooRealVar* m_jpsipi_mean1 = 0;
		RooRealVar* m_jpsipi_sigma1l = 0;
		RooRealVar* m_jpsipi_sigma1r = 0;
		m_jpsipi_fraction2 = new RooRealVar("m_jpsipi_fraction2","m_jpsipi_fraction2",0.4,0.0,0.8);
		m_jpsipi_mean1 = new RooRealVar("m_jpsipi_mean1","m_jpsipi_mean1",5.35, 5.3, 5.5);
		RooRealVar m_jpsipi_sigma2l("m_jpsipi_sigma2l","m_jpsipi_sigma2l",0.05,0.020,0.500);
		RooRealVar m_jpsipi_sigma2r("m_jpsipi_sigma2r","m_jpsipi_sigma2r",0.02,0.0050,0.500);
		m_jpsipi_sigma1l = new RooRealVar("m_jpsipi_sigma1l","m_jpsipi_sigma1l",0.05,0.010,0.150);
		m_jpsipi_sigma1r = new RooRealVar("m_jpsipi_sigma1r","m_jpsipi_sigma1r",0.17,0.010,0.350);
		RooBifurGauss m_jpsipi_gaussian2("m_jpsipi_gaussian2", "m_jpsipi_gaussian2", *mass, *m_jpsipi_mean1, m_jpsipi_sigma2l, m_jpsipi_sigma2r);
		RooBifurGauss m_jpsipi_gaussian1("m_jpsipi_gaussian1", "m_jpsipi_gaussian1", *mass, *m_jpsipi_mean1, *m_jpsipi_sigma1l, *m_jpsipi_sigma1r);
		RooAddPdf* jpsipi = new RooAddPdf("jpsipi", "jpsipi", RooArgList(m_jpsipi_gaussian2, m_jpsipi_gaussian1), RooArgList(*m_jpsipi_fraction2));
		//inclusive MC jpsipi Model
		//PDF MODELS PDF MODELS PDF MODELS

		// PREPARE DATA SETS
		vector<TString> jpsi_vars = {"By", "Bpt", "Bgen"};
		read_samples(*ws, jpsi_vars, ExtraMCsample.Data(), "ntKp", "jpsinp", SYSTEM.Data());
		RooDataSet* full_data_MC = (RooDataSet*) ws->data("jpsinp");

		// FORM PEAKING Background BINS
		RooDataSet* fullds_JPSI_shape_fix = (RooDataSet*)full_data_MC->reduce("Bgen == 23335");
		// FORM PEAKING Background BINS
		// PREPARE DATA SETS

		//[START] FIX SHAPE (J/Psi pi) 
		RooRealVar n_jpsipi_ext("n_jpsipi_ext", "n_jpsipi_ext", 1000 , 0., (fullds_JPSI_shape_fix->sumEntries())*2);
		RooExtendPdf jpsipi_ext("jpsipi_ext", "extended jpsipi", *jpsipi, n_jpsipi_ext);
		// FIT
		mass->setRange("bjpsipi", 5.2, 5.9);
		auto jpsipi_result = jpsipi_ext.fitTo(*fullds_JPSI_shape_fix, Range("bjpsipi"), Save(), Extended());
		// FIT
		//plot_mcfit(*ws, &jpsipi_ext, fullds_JPSI_shape_fix, Form("./results/BP/InclusiveMC_JPsipi_fit%s.pdf"), NormRange("bjpsipi"), DrawOption("LF"), FillStyle(3008), FillColor(kMagenta+1), LineStyle(1), LineColor(kMagenta+1), LineWidth(1)); 
		ws->import(*jpsipi);
		fix_parameters(*ws, "jpsipi" );
		//[END] FIX SHAPE (J/Psi pi) 
	}

	//BIN ANALYSIS START
	for(int i=0;i<_nBins;i++){
		_count++;
		TCanvas* c  = new TCanvas(Form("c%d"  ,_count),"",700,700);
		TCanvas* cMC= new TCanvas(Form("cMC%d",_count),"",700,700);

		double b_width = _varBINS[i+1]-_varBINS[i];

		//BINNING OF THE DATA
		RooDataSet* ith_DATA_bin = (RooDataSet*) data->reduce(Form("(abs(%s)>=%f && abs(%s)<=%f)",VAR.Data(),_varBINS[i],VAR.Data(),_varBINS[i+1]));
		RooDataSet* ith_MC_bin   = (RooDataSet*) mc  ->reduce(Form("(abs(%s)>=%f && abs(%s)<=%f)",VAR.Data(),_varBINS[i],VAR.Data(),_varBINS[i+1]));
		ith_DATA_bin->SetName(Form("data%d", _count));
		ith_MC_bin  ->SetName(Form("mc%d"  , _count));
		cout << "Total DATA entries in "<< i << "-th BIN: " << ith_DATA_bin->numEntries() << "\n";
		cout << "Total   MC entries in "<< i << "-th BIN: " << ith_MC_bin  ->numEntries() << "\n";
		RooRealVar * MC_candidates = new RooRealVar(Form("MC_candidates_%d",_count),"MC_candidates", ith_MC_bin->sumEntries());
		ws->import(*MC_candidates);
		//BINNING OF THE DATA

		if(VAR == "Bpt"){var_mean_av[i] = ith_DATA_bin->mean(*pt);}     	
		else if(VAR == "By"){
    		double sumAbs = 0.0;
			// Loop over the dataset and compute the sum of absolute values
			for (int iy = 0; iy < ith_DATA_bin->numEntries(); iy++) {
				RooRealVar* y_abs = (RooRealVar*) ith_DATA_bin->get(iy)->find("By");
				double abs_y = TMath::Abs(y_abs->getVal());
				sumAbs += abs_y;
			}
			var_mean_av[i] = sumAbs / ith_DATA_bin->numEntries();
		}
		else if(VAR == "nMult"){var_mean_av[i] = ith_DATA_bin->mean(*nMult);}
		
		// for the bin range in the histograms
		hori_av_low[i] = var_mean_av[i]-_varBINS[i];
		hori_av_high[i] = _varBINS[i+1]-var_mean_av[i];

		////////// FITFITFITFITFITFITFITFITFITFITFITFIT //////////

		cout << "Starting the fiting function for " << TREE.Data() << " " << VAR.Data() << " ["<< _varBINS[i] << ", " << _varBINS[i+1] << "]" << endl;
		RooFitResult* f_results = fit("", "", TREE, c, cMC, ith_DATA_bin, ith_MC_bin, mass, _varBINS[i], _varBINS[i+1], *ws, VAR.Data(),nbinsmasshisto);

		////////// FITFITFITFITFITFITFITFITFITFITFITFIT //////////
		


		// Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results 
		// Yield
		RooRealVar* bkgYield = static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index(Form("nbkg%d_%s",_count,""))));
		RooRealVar* fitYield = static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index(Form("nsig%d_%s",_count,""))));
		yield_vec[i]=fitYield->getVal();
		yield_Stat_unc[i] = fitYield->getError();	
		vector<double> stat_un;
		stat_un.push_back((double) yield_Stat_unc[i]/yield_vec[i]*100);

		//Resolution
		RooRealVar* width_scale = static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index("scale")));
		RooRealVar* sigma1 = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sigma1%d_", _count))));
		RooRealVar* sigma2 = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sigma2%d_", _count))));
		RooRealVar* weight = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sig1frac%d_", _count))));
		scale_vec[i]     = width_scale->getVal();
		scale_vec_unc[i] = width_scale->getError();
		resol_vec[i]     = sqrt(weight->getVal() * pow(sigma1->getVal(), 2) + (1 - weight->getVal()) * pow(sigma2->getVal(), 2)) * scale_vec[i] ;
		resol_vec_unc[i] = (scale_vec_unc[i] / scale_vec[i]) * resol_vec[i] ;

		//chi2
		TH1D* h = new TH1D(Form("h%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
		ith_DATA_bin->fillHistogram(h, *mass);
		RooDataHist* dh = new RooDataHist(Form("dh%d",_count),"",*mass,Import(*h));
		RooAbsPdf* model = (RooAbsPdf*)ws->pdf(Form("model%d_%s",_count,""));
		RooChi2Var chi2(Form("chi2%d",_count),"chi2",*model,*dh);
		chi2_vec[i] = chi2.getVal()/(nbinsmasshisto - f_results->floatParsFinal().getSize()); //normalised chi square

		cout << "Normalised Chi2 value (NdF=" << f_results->floatParsFinal().getSize() << "): " << chi2_vec[i] << endl;
		cout << "Probability of Chi2 value is " << TMath::Prob(chi2.getVal(), (nbinsmasshisto - f_results->floatParsFinal().getSize()) ) << endl;

		TH1D* hMC = new TH1D(Form("hMC%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
		ith_MC_bin->fillHistogram(hMC, *mass);
		RooDataHist* dhMC = new RooDataHist(Form("dhMC%d",_count),"",*mass,Import(*hMC));
		RooAbsPdf* modelMC = (RooAbsPdf*)ws->pdf(Form("modelMC%d_%s",_count,""));
		RooChi2Var chi2MC(Form("chi2MC%d",_count),"chi2MC",*modelMC,*dhMC);
		RooRealVar * fMC_params = new RooRealVar("fMC_params", "", ws->var(Form("ndfMC_%d_%s", _count, ""))->getVal() );
		chi2MC_vec[i] = chi2MC.getVal()/(nbinsmasshisto - fMC_params->getVal()); //normalised chi square
		// Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results 

		//////////////////////////// FILL HISTOGRAMS
		
		hPt->SetBinContent(i+1,yield_vec[i]/b_width);
		hPt->SetBinError(i+1,yield_Stat_unc[i]/b_width);
    	hPt->GetXaxis()->SetBinLabel(i+1, Form("%f_results", var_mean_av[i]));

		//////////////////////////// FILL HISTOGRAMS

		////////////////////////////////////////////////////////// LABELS IN PLOTS
	  	
		TLatex* tex_y = new TLatex(0.5,0.5,"");
		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.035);
		tex_y->SetLineWidth(2);

		TLatex* tex_yCUT = new TLatex(0.5,0.5,"");
		tex_yCUT->SetNDC();
		tex_yCUT->SetTextFont(42);
		tex_yCUT->SetTextSize(0.035);
		tex_yCUT->SetLineWidth(2);

		TLatex* chi_square = new TLatex(0.64, 0.55, Form("#chi^{2}/ndf = %.2f",chi2_vec[i]));
		chi_square->SetNDC();
		chi_square->SetTextFont(42);
		chi_square->SetTextSize(0.035);
		chi_square->SetLineWidth(2);
		chi_square->Draw();

		TLatex* varBIN = new TLatex(0.64, 0.50, Form("%d < p_{T} < %d [GeV/c]", (int)_varBINS[i], (int)_varBINS[i+1]));
		varBIN->SetNDC();
		varBIN->SetTextFont(42);
		varBIN->SetTextSize(0.035);
		varBIN->SetLineWidth(2);
		varBIN->Draw();

		if(_nBins == 1){ //inclusive bin case
			tex_yCUT->SetText(0.21, 0.68, "p_{T} < 10 GeV/c : 1.5 < |y| < 2.4" );
			//tex_yCUT->Draw();
			tex_y->SetText(0.21, 0.61, "p_{T} > 10 GeV/c : |y| < 2.4" );
			//tex_y->Draw();
		}
		else {
			if(_varBINS[i] >= 10 && _nBins != 1) {
				tex_y->SetText(0.21, 0.68, "|y| < 2.4");
				//tex_y->Draw();
			} else if( _varBINS[i+1] <= 10){
				tex_yCUT->SetText(0.21, 0.68, "1.5 < |y| < 2.4");
				//tex_yCUT->Draw();
			} 
		}	

		//SIGNIFICANCE CALCULATION
		TLatex *Signf = new TLatex(0.5,0.5,"");
		if (TREE == "ntmix" || TREE == "ntKstar" && VAR == "Bpt" && _nBins !=1){
			RooAbsPdf* bkg_pdf_ = ws->pdf(Form("bkg%d_%s", _count, "")) ;
			mass->setRange("sigX", X3872_MASS - 0.03, X3872_MASS + 0.03);
			RooAbsReal* intFrac = bkg_pdf_->createIntegral(*mass, NormSet(*mass), Range("sigX"));
			double fracInSignal = intFrac->getVal() ; 
			cout << "Background in Sig. Region " << fracInSignal << " ====> " << fracInSignal*bkgYield->getVal() << std::endl;
			cout << "Signal Yield: " << yield_vec[i] << endl;
			cout << "Significance: " << yield_vec[i] / sqrt(yield_vec[i] + fracInSignal*bkgYield->getVal()) << endl;
			Signf->SetText(0.64,0.35,Form("Significance = %.2f", yield_vec[i] / sqrt(yield_vec[i] + fracInSignal*bkgYield->getVal()) ) );
			Signf->SetNDC();
			Signf->SetTextFont(42);
			Signf->SetTextSize(0.035);
			Signf->SetLineWidth(2);
			Signf->Draw();
		}

		///////////////////////////////////////////////////////// /LABELS IN PLOTS


		CMS_lumi(c,19011,0);  //CMS PRELIMINARY + etc
		c->Update();

		if(VAR == "By"){
			c->SaveAs(  Form("%s/data_%s_%s_%0.1f_%0.1f_",OUTPLOTF.Data(),SYSTEM.Data(),Form("abs(%s)",VAR.Data()), (float)_varBINS[i],(float)_varBINS[i+1])+TREE+".pdf");
			cMC->SaveAs(Form("%s/mc_%s_%s_%0.1f_%0.1f_"  ,OUTPLOTF.Data(),SYSTEM.Data(),Form("abs(%s)",VAR.Data()), (float)_varBINS[i],(float)_varBINS[i+1])+TREE+".pdf");
		}else{
			c->SaveAs(  Form("%s/data_%s_%s_%i_%i_",OUTPLOTF.Data(),SYSTEM.Data(),VAR.Data(),(int)_varBINS[i],(int)_varBINS[i+1])+TREE+".pdf");
			cMC->SaveAs(Form("%s/mc_%s_%s_%i_%i_"  ,OUTPLOTF.Data(),SYSTEM.Data(),VAR.Data(),(int)_varBINS[i],(int)_varBINS[i+1])+TREE+".pdf");
		}

		vector<double> back_variation; 
		vector<double> back_unc; 
		vector<double> signal_variation; 
		vector<double> signal_unc; 
		vector<double> general_unc; 
		double max_signal=0.; 
		double max_back=0.; 
		TLatex* chi_back = new TLatex(0.5,0.5,"");  //for model systematics
		TLatex* chi_sig  = new TLatex(0.5,0.5,"");

		if(syst_study==1){
			//BACKGROUND MODEL SYSTEMATIC STUDY
			for(int j=0; j < static_cast<int>(background.size()); j++)
			{
				RooFitResult* f_back = fit("background", background[j], TREE, c, cMC, ith_DATA_bin, ith_MC_bin, mass, _varBINS[i], _varBINS[i+1], *ws, VAR, nbinsmasshisto);
				RooAbsPdf* model_back = (RooAbsPdf*)ws->pdf(Form("model%d_%s",_count,background[j].c_str()));
				TString chi2_fitRange = (background[j] == "mass_range") ? "m_range" : "all";
				RooChi2Var chi2_back("chi2_back","chi2_back",*model_back,*dh, Range(chi2_fitRange));
				chi2_vec_back[j][i] = chi2_back.getVal()/(nbinsmasshisto - f_back->floatParsFinal().getSize()); //normalised chi square

				RooAbsPdf* modelMC_back = (RooAbsPdf*)ws->pdf(Form("modelMC%d_%s",_count,background[j].c_str()));
				RooChi2Var chi2MC_back("chi2MC_back","chi2MC_back",*modelMC_back,*dhMC);
				RooRealVar * fMC_back_params = new RooRealVar("fMC_back_params", "", ws->var(Form("ndfMC_%d_%s", _count, background[j].c_str()))->getVal() );
				chi2MC_vec_back[j][i] = chi2MC_back.getVal()/(nbinsmasshisto - fMC_back_params->getVal());

				RooRealVar* fitYield_b_sys = static_cast<RooRealVar*>(f_back->floatParsFinal().at(f_back->floatParsFinal().index(Form("nsig%d_%s",_count, background[j].c_str()))));
				chi_square->SetText(0.21,0.54,Form("#chi^{2}/ndf = %.2f ",chi2_vec_back[j][i]));
				chi_square->Draw();

				CMS_lumi(c,19011,0);
				c->Update();

				if(VAR == "By"){c->SaveAs(Form("%s/data_%s_%s_%0.1f_%0.1f_%s_", OUTPLOTF.Data(), SYSTEM.Data(), Form("abs(%s)",VAR.Data()),(float)_varBINS[i],(float)_varBINS[i+1],background[j].c_str())+TREE+ ".pdf");}
				else { c->SaveAs(Form("%s/data_%s_%s_%i_%i_%s_", OUTPLOTF.Data(), SYSTEM.Data(), VAR.Data(),(int)_varBINS[i],(int)_varBINS[i+1],background[j].c_str())+TREE+".pdf");}

				RooRealVar* fitYield_back = static_cast<RooRealVar*>(f_back->floatParsFinal().at(f_back->floatParsFinal().index(Form("nsig%d_%s",_count,background[j].c_str()))));
				back_variation.push_back(fitYield_back->getVal());
				back_unc.push_back(abs(((yield_vec[i]-fitYield_back->getVal())/yield_vec[i])*100));
				if(abs(((yield_vec[i]-fitYield_back->getVal())/yield_vec[i])*100)>max_back) max_back=abs(((yield_vec[i]-fitYield_back->getVal())/yield_vec[i])*100);
			}
			general_unc.push_back(max_back);
			background_syst.push_back(back_variation);

			//SIGNAL MODEL SYSTEMATIC STUDY
			for(int j=0; j< static_cast<int>(signal.size()); j++)
			{
				RooFitResult* f_signal = fit("signal", signal[j], TREE, c, cMC, ith_DATA_bin, ith_MC_bin, mass, _varBINS[i], _varBINS[i+1], *ws, VAR, nbinsmasshisto);
				RooAbsPdf* model_sig = (RooAbsPdf*)ws->pdf(Form("model%d_%s",_count,signal[j].c_str()));
				RooChi2Var chi2_sig("chi2_sig","chi2_sig",*model_sig,*dh);
				chi2_vec_sig[j][i] = chi2_sig.getVal()/(nbinsmasshisto - f_signal->floatParsFinal().getSize()); //normalised chi square
				RooAbsPdf* modelMC_signal = (RooAbsPdf*)ws->pdf(Form("modelMC%d_%s",_count,signal[j].c_str()));
				RooChi2Var chi2MC_signal("chi2MC_signal","chi2MC_signal", *modelMC_signal, *dhMC);
				RooRealVar * fMC_signal_params = new RooRealVar("fMC_signal_params", "", ws->var(Form("ndfMC_%d_%s", _count, signal[j].c_str()))->getVal() );
				chi2MC_vec_sig[j][i] = chi2MC_signal.getVal()/(nbinsmasshisto - fMC_signal_params->getVal());
                
				RooRealVar* fitYield_b_sig = static_cast<RooRealVar*>(f_signal->floatParsFinal().at(f_signal->floatParsFinal().index(Form("nsig%d_%s",_count, signal[j].c_str()))));
				chi_square->SetText(0.21, 0.54, Form("#chi^{2}/ndf = %.2f ", chi2_vec_sig[j][i]));
				chi_square->Draw();

				CMS_lumi(c,19011,0);
				c->Update();
				
				if (signal[j] != "fixed") {
					if(VAR == "By"){ cMC->SaveAs(Form("%s/mc_%s_%s_%0.1f_%0.1f_%s_",OUTPLOTF.Data(),SYSTEM.Data(),Form("abs(%s)",VAR.Data()), (float)_varBINS[i], (float)_varBINS[i+1],signal[j].c_str())+TREE+".pdf");} 
					else { cMC->SaveAs(Form("%s/mc_%s_%s_%i_%i_%s_",OUTPLOTF.Data(),SYSTEM.Data(),VAR.Data(), (int)_varBINS[i], (int)_varBINS[i+1],signal[j].c_str() )+TREE+".pdf");}
				}
				if(VAR == "By"){ c->SaveAs(Form("%s/data_%s_%s_%0.1f_%0.1f_%s_",OUTPLOTF.Data(),SYSTEM.Data(),Form("abs(%s)",VAR.Data()),(float)_varBINS[i],(float)_varBINS[i+1],signal[j].c_str() )+TREE+".pdf");}
				else{ c->SaveAs(Form("%s/data_%s_%s_%i_%i_%s_",OUTPLOTF.Data(),SYSTEM.Data(),VAR.Data(),(int)_varBINS[i],(int)_varBINS[i+1],signal[j].c_str() )+TREE+".pdf");}
				
				RooRealVar* fitYield_signal = static_cast<RooRealVar*>(f_signal->floatParsFinal().at(f_signal->floatParsFinal().index(Form("nsig%d_%s",_count,signal[j].c_str()))));
				signal_variation.push_back(fitYield_signal->getVal());
				signal_unc.push_back(abs(((yield_vec[i]-fitYield_signal->getVal())/yield_vec[i])*100));
				if(abs(((yield_vec[i]-fitYield_signal->getVal())/yield_vec[i])*100)>max_signal) max_signal=abs(((yield_vec[i]-fitYield_signal->getVal())/yield_vec[i])*100);
			}

			general_unc.push_back(max_signal);
			signal_syst.push_back(signal_variation);
			general_unc.push_back(sqrt(max_back*max_back+max_signal*max_signal));
			stat_error.push_back(stat_un);
			back_syst_rel_values.push_back(back_unc);
			sig_syst_rel_values.push_back(signal_unc);			
			general_syst.push_back(general_unc);
			yield_vec_systerr_low[i] = general_unc[2] / 100 * yield_vec[i];
			yield_vec_systerr_high[i] = general_unc[2] / 100 * yield_vec[i];
		}

		//VALIDATION STUDIES
		if (val==1){
			string Path_val=Form("./%s/validation",OUTPLOTF.Data());
			validate_fit(ws, "", TREE, VAR, FULL, _varBINS[i], _varBINS[i+1],Path_val);
		}
		//VALIDATION STUDIES
	}
	//BIN ANALYSIS END
	//BIN ANALYSIS END
	
	TFile* outf = new TFile(Form("%s.root",outputfile.Data()),"recreate");
	outf->cd();
	hPt->Write();	
	outf->Close();

	if(syst_study==1 && FULL==0){
		vector<string> labels_back = {"1st Poly", "2nd Poly", "mass range", "jpsipi/JpsiK" };
		vector<string> col_name_back;
		vector<string> labels_signal = {"Triple Gaussian", "Fixed Mean", "CB+Gaussian", "Double CB"};
		vector<string> labels_general = {"Background", "Signal", "Total"};
		vector<string> labels_general_stat = {"Statistical Uncertainty"};
		vector<string> col_name_general;
		vector<string> col_name_signal;
		vector<string> col_name_general_stat;
		string name;
		col_name_back.push_back("Background Model");
		col_name_signal.push_back("Signal Model");
		col_name_general.push_back("Systematic Source");
		col_name_general_stat.push_back(" ");
		if(VAR=="Bpt"){name="$<p_T<$";} 
		else if(VAR=="By"){name="$<y<$";} 
		else if(VAR=="nMult"){name="$<nTrks<$";}
		for(int i=0;i<_nBins;i++){
			ostringstream clabel;
			clabel<<_varBINS[i]<<name<<_varBINS[i+1];
			string label1 = clabel.str();
			col_name_back.push_back(label1);
			col_name_signal.push_back(label1);
			col_name_general.push_back(label1);
			col_name_general_stat.push_back(label1);
		}

		latex_table("./files/background_systematics_table_" +string(VAR.Data())+"_"+string(TREE.Data()), _nBins+1, (int)(1+background.size()), col_name_back,labels_back,back_syst_rel_values, "Background PDF Systematic Errors");
		latex_table("./files/signal_systematics_table_"     +string(VAR.Data())+"_"+string(TREE.Data()), _nBins+1, (int)(1+signal.size()), col_name_signal, labels_signal,sig_syst_rel_values, "Signal PDF Systematic Errors");
		latex_table("./files/general_systematics_table_"    +string(VAR.Data())+"_"+string(TREE.Data()), _nBins+1, 4 , col_name_general, labels_general, general_syst, "Overall PDF Variation Systematic Errors");	
		latex_table("./files/Statistical_error_table_"      +string(VAR.Data())+"_"+string(TREE.Data()), _nBins+1, 2 , col_name_general_stat, labels_general_stat, stat_error, "Statistical error");	
		vector<string> tabeltype ={"background_systematics_table_", "signal_systematics_table_", "general_systematics_table_", "Statistical_error_table_"};
		vector<string> filetype ={"_check.aux", "_check.log", "_check.pdf"};
		for (int i=0;i<(int)(tabeltype.size());i++){
			for (int j=0;j<(int)(filetype.size());j++){
				rename((tabeltype[i]+string (VAR.Data())+"_"+string (TREE.Data())+filetype[j]).c_str(), ("./files/"+tabeltype[i]+string (VAR.Data())+"_"+string (TREE.Data())+filetype[j]).c_str());
			}
		}

		double zero[_nBins];
		for (int i=0;i<_nBins;i++){zero[i]=0.;}
		double low_high_b[_nBins];
		//These are only used to plot the systematics (PDF Var) in the same plot (it helps visualizing them!)
		Double_t x[_nBins];
		for (int i=0;i<_nBins;i++){ 
			x[i] = (_varBINS[i]+_varBINS[i+1])/2 ;
			low_high_b[i] = _varBINS[i+1] - x[i] ;
		}
		//These are only used to plot the systematics (PDF Var) in the same plot (it helps visualizing them!)

		TGraph *binning= new TGraphAsymmErrors (_nBins,x,zero,low_high_b,low_high_b,zero,zero);
		binning->SetMarkerColorAlpha(kBlack, 0); //transparent
		binning->SetLineWidth(6);

		TMultiGraph* m_back_sig= new TMultiGraph(); //to be used latter to acomodate both sig and back
		TLegend *legsigback = new TLegend(0.75,0.71,0.89,0.89, NULL, "brNDC");
		legsigback->SetBorderSize(0);
		legsigback->SetTextSize(0.035);
		legsigback->SetTextFont(42);
		legsigback->SetFillStyle(0);
		m_back_sig->Add(binning);
		m_back_sig->GetXaxis()->SetTitle("p_{T}");
		m_back_sig->GetYaxis()->SetTitle("Systematic Uncertainty(%)");

		TCanvas* c_back= new TCanvas("c_back","",700,700);
		TLegend *legback = new TLegend(0.8,0.75,0.89,0.89,NULL,"brNDC");
		legback->SetBorderSize(0);
		legback->SetTextSize(0.035);
		legback->SetTextFont(42);
		legback->SetFillStyle(0);
		TMultiGraph* m_back= new TMultiGraph();
		const char* backlabel[4]={"Linear", "2nd Poly", "mass range", "J/#psi#pi^{+}/J/#psiK^{+}"};
		double y_max_back=0;
		for (int j=0;j<(int)(background.size());j++){
			Double_t y[_nBins];
			for (int i=0;i<_nBins;i++){
				y[i]=	back_syst_rel_values[i][j];
				if (y[i]>y_max_back){y_max_back=y[i];}
			}
			TGraph *g_back= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
			g_back->SetMarkerColor(j+1);
			g_back->SetMarkerStyle(22);
			m_back->Add(g_back);
			m_back_sig->Add(g_back);
			legback->AddEntry(g_back,backlabel[j],"p");
			legsigback->AddEntry(g_back,backlabel[j],"p");
		}
		m_back->Add(binning);
		m_back->GetXaxis()->SetTitle("p_{T}");
		m_back->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
		m_back->GetYaxis()->SetRangeUser(0, y_max_back*1.5);
		m_back->Draw("AE1*");
		legback->Draw();
		c_back->SaveAs(Form("./results/tables/background_systematics_plot_%s_%s.pdf",TREE.Data(),VAR.Data())); 

		TCanvas* c_sig= new TCanvas("c_sig","",700,700);
		TLegend *legsig = new TLegend(0.75,0.8,0.89,0.89,NULL,"brNDC");
		legsig->SetBorderSize(0);
		legsig->SetTextSize(0.035);
		legsig->SetTextFont(42);
		legsig->SetFillStyle(0);
		TMultiGraph* m_sig= new TMultiGraph();
		const char* siglabel[3]={"Triple Gaussian", "Fixed Mean", "CB+Gaussian"};
		double y_max_sig=0;
		for (int j=0;j<(int)(signal.size());j++){
			Double_t y[_nBins];
			for (int i=0;i<_nBins;i++){
				y[i]=	sig_syst_rel_values[i][j];
				if (y[i]>y_max_sig){y_max_sig=y[i];}
			}
			TGraph *g_sig= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
			g_sig->SetMarkerColor(j+5);
			g_sig->SetMarkerStyle(21);
			m_sig->Add(g_sig);
			m_back_sig->Add(g_sig);
			legsig->AddEntry(g_sig,siglabel[j],"p");
			legsigback->AddEntry(g_sig,siglabel[j],"p");
		}
		m_sig->Add(binning);
		m_sig->GetXaxis()->SetTitle("p_{T}");
		m_sig->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
		m_sig->GetYaxis()->SetRangeUser(0, y_max_sig*1.5);
		m_sig->Draw("AE1*");
		legsig->Draw();
		c_sig->SaveAs(Form("./results/tables/signal_systematics_plot_%s_%s.pdf",TREE.Data(),VAR.Data())); 

		TCanvas *c_sig_back= new TCanvas("c_sig_back","",700,700);
		m_back_sig->GetYaxis()->SetRangeUser(0, 4);
		
		if(VAR == "By"){
			m_back_sig->GetXaxis()->SetTitle("Rapidity (y)");
			m_back_sig->GetYaxis()->SetTitle("dY_{S}/dy");
			m_back_sig->GetXaxis()->SetLimits(0,2.4);
		}
		else if(VAR == "Bpt"){
			m_back_sig->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			m_back_sig->GetYaxis()->SetTitle("dY_{S}/dp_{T}");
			if (TREE == "ntKp"){ m_back_sig->GetXaxis()->SetLimits(0 ,65); }
			if (TREE == "ntphi"){ m_back_sig->GetXaxis()->SetLimits(0 ,65); }
		}
		else if(VAR == "nMult"){
			m_back_sig->GetXaxis()->SetTitle("Multiplicity (Mult)");
			m_back_sig->GetYaxis()->SetTitle("dY_{S}/dMult");
			m_back_sig->GetXaxis()->SetLimits(0, 110);
		}
		m_back_sig->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
		m_back_sig->Draw("AE1*");
		legsigback->Draw();
		c_sig_back->SaveAs(Form("./results/tables/background_signal_systematics_plot_%s_%s.pdf",TREE.Data(),VAR.Data())); 

		TCanvas* c_gen= new TCanvas("c_gen","",700,700);
		TLegend *legen = new TLegend(0.8,0.77,0.89,0.89,NULL,"brNDC");
		legen->SetBorderSize(0);
		legen->SetTextSize(0.035);
		legen->SetTextFont(42);
		legen->SetFillStyle(0);	
		TMultiGraph* m_gen= new TMultiGraph();
		const char* genlabel[3]={"Background", "Signal", "Total"};
		double y_max_gen=0;
		for (int j=0;j<(int)(labels_general.size());j++){
			Double_t y[_nBins];
			for (int i=0;i<_nBins;i++){
				y[i]=general_syst[i][j];
				if (y[i]>y_max_gen){y_max_gen=y[i];}
			}
			TGraph *g_gen= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
			g_gen->SetMarkerColor(j+1);
			g_gen->SetMarkerStyle(21);
			m_gen->Add(g_gen);
			legen->AddEntry(g_gen,genlabel[j],"p");
		}
		Double_t y[_nBins];
		for (int i=0;i<_nBins;i++){y[i]=	stat_error[i][0];}
		TGraph *g_gen= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
		m_gen->Add(binning);
		g_gen->SetMarkerColor(9);
		g_gen->SetMarkerStyle(21);
		m_gen->Add(g_gen);
		legen->AddEntry(g_gen,"Statistical","p");
		m_gen->GetXaxis()->SetTitle("p_{T}");  //VAR.Data()
		m_gen->GetYaxis()->SetTitle("Total Uncertainty(%)");
		m_gen->GetYaxis()->SetRangeUser(0, y_max_gen*1.5);
		m_gen->Draw("AE1*");
		legen->Draw();
		c_gen->SaveAs(Form("./results/tables/general_systematics_plot_%s_%s.pdf",TREE.Data(),VAR.Data())); 
	}

	//Differential plot part
	TCanvas c_diff;
	TMultiGraph* mg = new TMultiGraph();
	TLegend *leg_d = new TLegend(0.7,0.7,0.9,0.9);
	TGraphAsymmErrors* gr_staterr = new TGraphAsymmErrors(_nBins,var_mean_av,yield_vec,hori_av_low,hori_av_high,yield_Stat_unc,yield_Stat_unc);
	gr_staterr->SetLineColor(1); 
	gr_staterr->SetName("Y_stat");
	mg->Add(gr_staterr, "stat");
	leg_d->AddEntry(gr_staterr, "Statistical Uncertainty", "e");
	if(syst_study==1){
		TGraphAsymmErrors* gr_systerr = new TGraphAsymmErrors(_nBins, var_mean_av, yield_vec, nullptr, nullptr, yield_vec_systerr_low, yield_vec_systerr_high);
		gr_systerr->SetLineColor(2);
	 	gr_staterr->SetName("Y_syst");
		mg->Add(gr_systerr,"syst");
		leg_d->AddEntry(gr_systerr, "Systematic Uncertainty", "e");
	}
	if(VAR == "By"){
		mg->GetXaxis()->SetTitle("Rapidity (y)");
		mg->GetYaxis()->SetTitle("dY_{S}/dy");
		mg->GetXaxis()->SetLimits(0,2.4);
	} else if(VAR == "Bpt"){
		mg->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		mg->GetYaxis()->SetTitle("dY_{S}/dp_{T}");
		mg->GetXaxis()->SetLimits(0 ,80);
	} else if(VAR == "nMult"){
		mg->GetXaxis()->SetTitle("Multiplicity (Mult)");
		mg->GetYaxis()->SetTitle("dY_{S}/dMult");
		mg->GetXaxis()->SetLimits(0, 110);
	}
	mg->Draw("ap");
	leg_d->SetBorderSize(0);
	leg_d->SetFillStyle(0);
	leg_d->SetTextSize(0);
	leg_d->Draw();
	const char* pathc =Form("./results/Graphs/raw_yield_%s_%s.pdf", TREE.Data(), VAR.Data());
	c_diff.SaveAs(pathc);
	//Differential plot part ends

	//Scale part starts
	double scale_max = 0;
	for(int i = 0; i < _nBins; i++){if(scale_vec[i] > scale_max){scale_max = scale_vec[i];}}

	TCanvas c_par;
	TMultiGraph* mg_par = new TMultiGraph();
	TGraphAsymmErrors* gr_scale = new TGraphAsymmErrors(_nBins,var_mean_av,scale_vec,hori_av_low,hori_av_high,scale_vec_unc,scale_vec_unc);
	gr_scale->SetLineColor(1); 

	if(VAR == "By"){
		mg_par->GetXaxis()->SetTitle("Rapidity (y)");
		mg_par->GetXaxis()->SetLimits(0 ,2.4);
	} else if(VAR == "Bpt"){
		mg_par->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		mg_par->GetXaxis()->SetLimits(0 ,80); 
	} else if(VAR == "nMult"){
		mg_par->GetXaxis()->SetTitle("Multiplicity (Mult)");
		mg_par->GetXaxis()->SetLimits(0, 110);
	}
	mg_par->GetYaxis()->SetTitle("Scale factor");
	mg_par->Add(gr_scale);
	mg_par->GetYaxis()->SetRangeUser(0,scale_max*1.4);
	mg_par->Draw("ap");
	const char* pathc_par =Form("./results/Graphs/scale_variation_%s_%s.pdf",TREE.Data(),VAR.Data()); 
	c_par.SaveAs(pathc_par);
	//Scale part ends

	//Resolution plot part
	TCanvas c_resol;
	TMultiGraph* mg_resol = new TMultiGraph();
	TGraphAsymmErrors* gr_resol = new TGraphAsymmErrors(_nBins, var_mean_av, resol_vec, hori_av_low, hori_av_high, resol_vec_unc, resol_vec_unc);
	gr_resol->SetLineColor(1); 
	mg_resol->GetYaxis()->SetTitle("Resolution");
	if(VAR == "By"){
		mg_resol->GetXaxis()->SetTitle("Rapidity (y)");
		mg_resol->GetXaxis()->SetLimits(0 ,2.4);
	} else if(VAR == "Bpt"){
		mg_resol->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		mg_resol->GetXaxis()->SetLimits(0, 100); 
	} else if(VAR == "nMult"){
		mg_resol->GetXaxis()->SetTitle("Multiplicity (Mult)");
		mg_resol->GetXaxis()->SetLimits(0, 110);
	}
	mg_resol->GetYaxis()->SetRangeUser(0,0.2);
	mg_resol->Add(gr_resol);
	mg_resol->Draw("ap");
	const char* pathc_resol =Form("./results/Graphs/resolution_%s_%s.pdf",TREE.Data(),VAR.Data()); 
	c_resol.SaveAs(pathc_resol);
	//Resolution plot part ends

	//Chi2 plot part 
	double chi2_max = 0;
	for(int i = 0; i < _nBins; i++){
		if(chi2_vec[i] > chi2_max){chi2_max = chi2_vec[i];}
		if(chi2MC_vec[i] > chi2_max){chi2_max = chi2MC_vec[i];}
	}
	TCanvas c_chi2;
	TMultiGraph* mg_chi2 = new TMultiGraph();
	TLegend *leg_chi2 = new TLegend(0.7,0.8,0.9,0.9);
	TGraphAsymmErrors* gr_chi2   = new TGraphAsymmErrors(_nBins,var_mean_av,chi2_vec,hori_av_low,hori_av_high,nullptr,nullptr);
	gr_chi2->SetLineColor(1); 
	TGraphAsymmErrors* grMC_chi2 = new TGraphAsymmErrors(_nBins,var_mean_av,chi2MC_vec,hori_av_low,hori_av_high,nullptr,nullptr);
	grMC_chi2->SetLineColor(2); 

	if(VAR == "By"){
		mg_chi2->GetXaxis()->SetTitle("Rapidity (y)");
		mg_chi2->GetXaxis()->SetLimits(0 ,2.4);
	} else if(VAR == "Bpt"){
		mg_chi2->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		mg_chi2->GetXaxis()->SetLimits(0 ,80); 
	} else if(VAR == "nMult"){
		mg_chi2->GetXaxis()->SetTitle("Multiplicity (Mult)");
		mg_chi2->GetXaxis()->SetLimits(0, 110);
	}
	mg_chi2->GetYaxis()->SetTitle("#chi^{2}/NDF");
	mg_chi2->GetYaxis()->SetRangeUser(0.0, chi2_max*1.5);
	mg_chi2->Add(gr_chi2);
	mg_chi2->Add(grMC_chi2);
	mg_chi2->Draw("ap");
	leg_chi2->AddEntry(gr_chi2, "Data", "e");
	leg_chi2->AddEntry(grMC_chi2, "MC", "e");
	leg_chi2->SetBorderSize(0);
	leg_chi2->SetFillStyle(0);
	leg_chi2->SetTextSize(0);
	leg_chi2->Draw();

	const char* pathc_chi2 =Form("./results/Graphs/chi2_%s_%s.pdf",TREE.Data(),VAR.Data()); 
	c_chi2.SaveAs(pathc_chi2);
	// Nominal part ONLY

	if(syst_study==1){
		//chi2 plot part (sigsum) starts
		TCanvas c_chi2_sigsum;
		TMultiGraph* mg_chi2_sigsum = new TMultiGraph();
		TLegend *leg_chi2_sigsum = new TLegend(0.7,0.8,0.9,0.9);
		double chi2_max_sigsum = 0;

		for(int j=0; j<static_cast<int>(signal.size()); j++){
			for(int i = 0; i < _nBins; i++){if(chi2_vec_sig[j][i] > chi2_max_sigsum){chi2_max_sigsum = chi2_vec_sig[j][i];}}
			TGraphAsymmErrors* gr_chi2_sigsum = new TGraphAsymmErrors(_nBins,var_mean_av,chi2_vec_sig[j],hori_av_low,hori_av_high,nullptr,nullptr);
			gr_chi2_sigsum->SetLineColor(j+2);
			mg_chi2_sigsum->Add(gr_chi2_sigsum);
			leg_chi2_sigsum->AddEntry(gr_chi2_sigsum, Form("%s",signal[j].c_str()), "e");
		}
		if(VAR == "By"){
			mg_chi2_sigsum->GetXaxis()->SetTitle("Rapidity (y)");
			mg_chi2_sigsum->GetXaxis()->SetLimits(0 ,2.4);
		} else if(VAR == "Bpt"){
			mg_chi2_sigsum->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			mg_chi2_sigsum->GetXaxis()->SetLimits(0 ,80);
		} else if(VAR == "nMult"){
			mg_chi2_sigsum->GetXaxis()->SetTitle("Multiplicity (Mult)");
			mg_chi2_sigsum->GetXaxis()->SetLimits(0, 110);
		}
		mg_chi2_sigsum->GetYaxis()->SetTitle("#chi^{2}/NDF");
		mg_chi2_sigsum->Add(gr_chi2);
		mg_chi2_sigsum->GetYaxis()->SetRangeUser(0.0, chi2_max_sigsum*1.4);
		mg_chi2_sigsum->Draw("ap");
		leg_chi2_sigsum->AddEntry(gr_chi2, "Nominal", "e");
		leg_chi2_sigsum->SetFillStyle(0);
		leg_chi2_sigsum->SetTextSize(0);
		leg_chi2_sigsum->Draw();

		const char* pathc_chi2_sigsum =Form("./results/Graphs/chi2_%s_%s_signal_summary.pdf",TREE.Data(),VAR.Data()); 
		c_chi2_sigsum.SaveAs(pathc_chi2_sigsum);
		//chi2 plot part (sigsum) ends

		//chi2 plot part (backsum) starts
		TCanvas c_chi2_backsum;
		TMultiGraph* mg_chi2_backsum = new TMultiGraph();
		TLegend *leg_chi2_backsum = new TLegend(0.7,0.8,0.9,0.9);

		double chi2_max_backsum = 0;

		for(int j=0; j<static_cast<int>(signal.size()); j++){
			for(int i = 0; i < _nBins; i++){if(chi2_vec_back[j][i] > chi2_max_backsum){chi2_max_backsum = chi2_vec_back[j][i];}}
			TGraphAsymmErrors* gr_chi2_backsum = new TGraphAsymmErrors(_nBins,var_mean_av,chi2_vec_back[j],hori_av_low,hori_av_high,nullptr,nullptr);
			gr_chi2_backsum->SetLineColor(j+2);
			mg_chi2_backsum->Add(gr_chi2_backsum);
			leg_chi2_backsum->AddEntry(gr_chi2_backsum, Form("%s",background[j].c_str()), "e");
		}

		if(VAR == "By"){
			mg_chi2_backsum->GetXaxis()->SetTitle("Rapidity (y)");
			mg_chi2_backsum->GetXaxis()->SetLimits(0,2.4);
		} else if(VAR == "Bpt"){
			mg_chi2_backsum->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			mg_chi2_backsum->GetXaxis()->SetLimits(0 ,80);
		} else if(VAR == "nMult"){
			mg_chi2_backsum->GetXaxis()->SetTitle("Multiplicity (Mult)");
			mg_chi2_backsum->GetXaxis()->SetLimits(0, 110);
		}
		mg_chi2_backsum->GetYaxis()->SetTitle("#chi^{2}/NDF");
		mg_chi2_backsum->Add(gr_chi2);
		mg_chi2_backsum->GetYaxis()->SetRangeUser(0.0, chi2_max_backsum*1.4);
		mg_chi2_backsum->Draw("ap");
		leg_chi2_backsum->AddEntry(gr_chi2, "Nominal", "e");
		leg_chi2_backsum->SetFillStyle(0);
		leg_chi2_backsum->SetTextSize(0);
		leg_chi2_backsum->Draw();

		const char* pathc_chi2_backsum =Form("./results/Graphs/chi2_%s_%s_background_summary.pdf",TREE.Data(),VAR.Data()); 
		c_chi2_backsum.SaveAs(pathc_chi2_backsum);
		//chi2 plot part (backsum) ends
	}
	//Chi2 plot part ends
}









void read_samples(RooWorkspace& w, vector<TString> label, TString fName, TString treeName, TString sample, TString colsys, TString DOselCUTS){
	
	cout << "Reading " <<  colsys.Data() << " " << treeName.Data() << " " << sample.Data() << " samples" << endl;
	cout << "at: " << fName.Data() << endl;

	TFile* fin = new TFile(fName);
	TTree* tIn = (TTree*) fin->Get(treeName);
	//std::cout << "Tree entries: " << tIn->GetEntries() << std::endl;

	// Ensure data is well selected
	TString SELcuts    = SELcuts_ppRef;
	TString trgmatches = TRGmatching;
	TString ACCcuts    = ACCcuts_ppRef;
	if (treeName == "ntKp") {
		SELcuts = SELcuts_ppRef_Bu;
		ACCcuts = ACCcuts_ppRef_Bu;
	}
	if (colsys == "pbpb" || colsys == "PbPb") {
		SELcuts = SELcuts_PbPb;
		ACCcuts = ACCcuts_PbPb;
		if (treeName == "ntKp") {
			SELcuts = SELcuts_PbPb_Bu;
			ACCcuts = ACCcuts_PbPb_Bu;
		}
		trgmatches = "1";
	}

	TString fullCut = Form("(%s) && (%s) && (%s) && "
						"((Bpt < 10 && abs(By) > 1.5) || (Bpt > 10)) && "
						"(%s) &&"
						"(Bmass > %f) && (Bmass < %f)",
						ACCcuts.Data(), SELcuts.Data(), trgmatches.Data(),
						DOselCUTS.Data(),
						minhisto, maxhisto);

	TDirectory* savedDir = gDirectory;
	gROOT->cd();                 
	TTree* t1 = tIn->CopyTree(fullCut);
	savedDir->cd();              
	//std::cout << "Skimmed entries: " << t1->GetEntries() << std::endl;

	if(true){
		// Create a canvas to draw the mass histogram
        TCanvas *canvas = new TCanvas("canvas", "Bmass Distribution", 600, 600);
        canvas->SetLeftMargin(0.15); // or try 0.18 for more space

        // Define histogram parameters
        double hist_Xlow  = minhisto;  // Minimum Bmass
        double hist_Xhigh = maxhisto;   // Maximum Bmass
        double bin_length_MEV = (hist_Xhigh - hist_Xlow)*1000 / nbinsmasshisto;

        // Create a histogram for Bmass
        TH1F *hist_Bmass = new TH1F("hist_Bmass", Form("; m_{J/#Psi #pi^{+} #pi^{-}} ; Entries / %.1f MeV", bin_length_MEV), nbinsmasshisto, hist_Xlow, hist_Xhigh);
		t1->Draw("Bmass >> hist_Bmass");
		hist_Bmass->SetLineColor(kBlack);
		hist_Bmass->SetLineWidth(1);
		hist_Bmass->SetFillColor(kBlack);    
		hist_Bmass->SetFillStyle(3017); 
		hist_Bmass->Draw("");
		gPad->Update();

		// Save the canvas as an image
		canvas->SaveAs(Form("%s_%s_Bmass.pdf", sample.Data(), treeName.Data()));
		// Clean up
		delete hist_Bmass;
		delete canvas;
	}

	// Get the variables	
	RooArgList arg_list("arg_list");
	arg_list.add(*(w.var("Bmass")));
	for(auto lab : label){arg_list.add(*(w.var(lab)));}
	RooDataSet* data_s = new RooDataSet(sample, sample, t1, arg_list);
	w.import(*data_s);
}

std::pair<int, std::vector<double>>
defineBinning(const TString& var, const TString& tree, int full)
{
    // Number of bins
    int nBins = 1;

    if (var == "Bpt" && full == 0) {
        if (tree == "ntmix") nBins = N_pt_Bins_X;
        else                 nBins = N_pt_Bins_B;
    } else if (var == "By")   {nBins = N_y_Bins_X;
    } else if (var == "nMult"){nBins = N_mult_Bins_X;
    } else if (var == "Cent") {nBins = N_cent_Bins_X;}

    std::vector<double> varBINS;
    varBINS.resize(nBins + 1);

    if(var=="Bpt"){
        if(full == 1){
            if(tree == "ntmix"){
                varBINS[0] = ptbinsvec_X.front();
                varBINS[1] = ptbinsvec_X.back();
            }else{
                varBINS[0] = ptbinsvec_B.front();
                varBINS[1] = ptbinsvec_B.back();
            }
        }else{
            if(tree=="ntmix"){for(int c = 0; c <= nBins; ++c){varBINS[c] = ptbinsvec_X[c];}} 
			else{for(int c = 0; c <= nBins; ++c){varBINS[c] = ptbinsvec_B[c];}}
        }
    }
	else if (var == "By")   {for(int c = 0; c <= nBins; ++c){varBINS[c] = ybinsvec[c];}}
	else if (var == "nMult"){for(int c = 0; c <= nBins; ++c){varBINS[c] = nmbinsvec[c];}} 
	else if (var == "Cent") {for(int c = 0; c <= nBins; ++c){varBINS[c] = centbinsvec[c];}}

    return { nBins, varBINS };
}

