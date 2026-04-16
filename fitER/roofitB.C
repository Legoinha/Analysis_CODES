#include "roofitB.h"
#include "TSystem.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include <TMath.h>
#include <string>
#include <sstream>
#include <TGraph.h>
#include <stdio.h>

#include "../plotER/aux/parameters.h"   //

void read_samples(RooWorkspace& w, vector<TString> label, TString fName, TString treeName, TString sample, TString system="ppRef", TString DOselCUTS="1");
std::pair<int, std::vector<double>> defineBinning(const TString& var, const TString& tree, int full);

// PDF VARIATION FOR SYST STUDIES
int syst_study=1;

// VALIDATION STUDIES
int val=0;

void roofitB(TString TREE = "ntphi", int FULL = 0, TString INPUTDATA = "", TString INPUTMC = "", TString VAR = "", TString CUT = "", TString OUTPLOTF = "", TString ExtraMCsample = "", TString SYSTEM = "ppRef"){

	//Setup the working area
	gSystem->mkdir(Form("./%s/validation",OUTPLOTF.Data()),true); 
	gSystem->mkdir("./files",true);
	gSystem->mkdir("./results/tables",true); 
	gSystem->mkdir("./results/Graphs", true); 
	gSystem->mkdir(Form("%s", OUTPLOTF.Data()),true); 
	
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

	TString SELcuts = Form("%s",CUT.Data()); 
	cout << "SELECTION cuts = "<< SELcuts << endl;
	cout << endl << endl;

	//DEFINE VARIABLES
	RooRealVar* mass = nullptr;
	if (TREE == "ntmix"){ minhisto = minhisto_X, maxhisto = maxhisto_X;} 
	else { minhisto = minhisto_B, maxhisto = maxhisto_B; }
	mass = new RooRealVar("Bmass", "Bmass", minhisto, maxhisto);
	mass->setRange("m_rangeB", 5.2 , 5.5);        //set a range to be used if pdf = mass_rangeB
	mass->setRange("all", minhisto, maxhisto);
	RooRealVar* pt    = new RooRealVar("Bpt","Bpt",0,300);
	RooRealVar* y     = new RooRealVar("By","By",-2.4, 2.4);
	RooRealVar* nSelectedChargedTracks = new RooRealVar("nSelectedChargedTracks","nSelectedChargedTracks",0,2000000000);
	RooRealVar* CentBin = new RooRealVar("CentBin","CentBin",0,100);
	RooRealVar* isX3872 = new RooRealVar("isX3872","isX3872",0 ,1); // 1 = X(3872) , 0 = PSI2S

	RooWorkspace* ws = new RooWorkspace("ws");
	ws->import(*mass);
	ws->import(*y);
	ws->import(*pt);
	ws->import(*nSelectedChargedTracks);
	ws->import(*CentBin);
	ws->import(*isX3872); // <- only needed for X3872 analysis

	//DATA and MC SAMPLES
	vector<TString>   ANA_vars = {"Bpt", "By", "CentBin", "nSelectedChargedTracks"};
	read_samples(*ws, ANA_vars, INPUTDATA.Data(), TREE.Data(), "data", SYSTEM.Data(), SELcuts);
	read_samples(*ws, ANA_vars, INPUTMC.Data()  , TREE.Data(),   "mc", SYSTEM.Data(), SELcuts);
	RooDataSet* data = (RooDataSet*) ws->data("data");
	RooDataSet* mc   = (RooDataSet*) ws->data("mc");
	//DATA and MC SAMPLES


	//MODELS for syst studies
	vector<SystVariationConfig> background = GetBackgroundSystematicModels(TREE);
	vector<SystVariationConfig> signal = GetSignalSystematicModels(TREE);
	vector<vector<double>> background_syst;
	vector<vector<double>> signal_syst;
	vector<vector<double>> general_syst;
	vector<vector<double>> back_syst_rel_values;
	vector<vector<double>> sig_syst_rel_values;
	vector<vector<double>> stat_error;
	vector<vector<double>> general_syst_2S;
	vector<vector<double>> back_syst_rel_values_2S;
	vector<vector<double>> sig_syst_rel_values_2S;
	vector<vector<double>> stat_error_2S;
	//MODELS for syst studies

	// FIT INFO
	double yield_vec[_nBins];
	double yield_Stat_unc[_nBins];
	double yield_vec_systerr_low[_nBins];
	double yield_vec_systerr_high[_nBins];
	double yield_vec_2S[_nBins];
	double yield_Stat_unc_2S[_nBins];
	double yield_vec_2S_systerr_low[_nBins];
	double yield_vec_2S_systerr_high[_nBins];
	double scale_vec[_nBins];
	double scale_vec_unc[_nBins];
	double scale_vec_2S[_nBins];
	double scale_vec_2S_unc[_nBins];
	double resol_vec_unc[_nBins];
	double resol_vec_2S_unc[_nBins];
	double var_mean_av[_nBins];
	double hori_av_low[_nBins];
	double hori_av_high[_nBins];
	double chi2_vec[_nBins];
	double chi2MC_vec[_nBins];
	double chi2MC_PSI2S_vec[_nBins];
	vector<vector<double>> chi2_vec_sig(signal.size(), vector<double>(_nBins, -1.0));
	vector<vector<double>> chi2_vec_back(background.size(), vector<double>(_nBins, -1.0));
	double resol_vec[_nBins];
	double resol_vec_2S[_nBins];
	// FIT INFO

	TH1D* hPt = new TH1D("hPt","",_nBins,_varBINS.data());  
	TH1D* hPt_2S = new TH1D("hPt_2S","",_nBins,_varBINS.data());
	TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_varBINS.data());
	TH1D* hPtMC_2S = new TH1D("hPtMC_2S","",_nBins,_varBINS.data());

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
		else if(VAR == "nSelectedChargedTracks"){var_mean_av[i] = ith_DATA_bin->mean(*nSelectedChargedTracks);}
		
		// for the bin range in the histograms
		hori_av_low[i]  = var_mean_av[i]-_varBINS[i];
		hori_av_high[i] = _varBINS[i+1]-var_mean_av[i];

		////////// FITFITFITFITFITFITFITFITFITFITFITFIT //////////

		cout << "Starting the fiting function for " << TREE.Data() << " " << VAR.Data() << " ["<< _varBINS[i] << ", " << _varBINS[i+1] << "]" << endl;
		RooFitResult* f_results = fit("", "", TREE, c, cMC, ith_DATA_bin, ith_MC_bin, mass, _varBINS[i], _varBINS[i+1], *ws, VAR.Data(),nbinsmasshisto);
		ws->saveSnapshot(Form("nominalPars_bin%d", _count), f_results->floatParsFinal(), true);

		////////// FITFITFITFITFITFITFITFITFITFITFITFIT //////////


		// Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results 
		// Yield
			RooRealVar* fitYield = static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index(Form("nsig%d_%s",_count,""))));

			yield_vec[i]=fitYield->getVal();
			yield_Stat_unc[i] = fitYield->getError();	
			vector<double> stat_un;
			stat_un.push_back((double) yield_Stat_unc[i]/yield_vec[i]*100);
			vector<double> stat_un_2S;
			if (TREE == "ntmix") {
				RooRealVar* fitYield_spec = static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index(Form("nsig_spec%d_%s",_count,""))));
				yield_vec_2S[i] = fitYield_spec->getVal();
				yield_Stat_unc_2S[i] = fitYield_spec->getError();
				stat_un_2S.push_back((double) yield_Stat_unc_2S[i]/yield_vec_2S[i]*100);
			}

		//Resolution
		RooRealVar* width_scale = static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index("scale")));
		RooRealVar* sigma1 = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sigma1%d_", _count))));
		RooRealVar* sigma2 = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sigma2%d_", _count))));
		RooRealVar* weight = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sig1frac%d_", _count))));
			scale_vec[i]     = width_scale->getVal();
			scale_vec_unc[i] = width_scale->getError();
			resol_vec[i]     = sqrt(weight->getVal() * pow(sigma1->getVal(), 2) + (1 - weight->getVal()) * pow(sigma2->getVal(), 2)) * scale_vec[i] ;
			resol_vec_unc[i] = (scale_vec_unc[i] / scale_vec[i]) * resol_vec[i] ;
			if (TREE == "ntmix") {
				RooRealVar* width_scale_spec = static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index("scale_spec")));
				RooRealVar* sigma1_spec = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sigma1_spec%d_", _count))));
				RooRealVar* sigma2_spec = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sigma2_spec%d_", _count))));
				RooRealVar* weight_spec = static_cast<RooRealVar*>(f_results->constPars().at(f_results->constPars().index(Form("sig1frac_spec%d_", _count))));
				scale_vec_2S[i]     = width_scale_spec->getVal();
				scale_vec_2S_unc[i] = width_scale_spec->getError();
				resol_vec_2S[i]     = sqrt(weight_spec->getVal() * pow(sigma1_spec->getVal(), 2) + (1 - weight_spec->getVal()) * pow(sigma2_spec->getVal(), 2)) * scale_vec_2S[i];
				resol_vec_2S_unc[i] = (scale_vec_2S_unc[i] / scale_vec_2S[i]) * resol_vec_2S[i];
			}

		//chi2 (computed in fit() and stored in workspace)
		RooRealVar* chi2_data_norm = ws->var(Form("chi2_data_norm%d_%s", _count, ""));
		RooRealVar* chi2MC_norm  = ws->var(Form("chi2MC_norm%d_%s",  _count, ""));
		RooRealVar* chi2MC_PSI2S_norm = ws->var(Form("chi2MC_PSI2S_norm%d_%s", _count, ""));
		chi2_vec[i]   = (chi2_data_norm ? chi2_data_norm->getVal() : -1.0);
		chi2MC_vec[i] = (chi2MC_norm  ? chi2MC_norm->getVal()  : -1.0);
		chi2MC_PSI2S_vec[i] = (chi2MC_PSI2S_norm ? chi2MC_PSI2S_norm->getVal() : -1.0);
		// Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results // Get Fit results 

		//////////////////////////// FILL YIELD HISTOGRAMS
		hPt->SetBinContent(i+1,yield_vec[i]/b_width);
		hPt->SetBinError(i+1,yield_Stat_unc[i]/b_width);
    	hPt->GetXaxis()->SetBinLabel(i+1, Form("%.2f", var_mean_av[i]));
		RooRealVar* nsigMC = ws->var(Form("nsigMC%d_%s",_count,""));
		hPtMC->SetBinContent(i+1,nsigMC->getVal()/b_width);
		hPtMC->SetBinError(i+1,nsigMC->getError()/b_width);
		hPtMC->GetXaxis()->SetBinLabel(i+1, Form("%.2f", var_mean_av[i]));
		if (TREE == "ntmix") {
			hPt_2S->SetBinContent(i+1,static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index(Form("nsig_spec%d_%s",_count,""))))->getVal()  /b_width);
			hPt_2S->SetBinError(i+1  ,static_cast<RooRealVar*>(f_results->floatParsFinal().at(f_results->floatParsFinal().index(Form("nsig_spec%d_%s",_count,""))))->getError()/b_width);
			hPt_2S->GetXaxis()->SetBinLabel(i+1, Form("%.2f", var_mean_av[i]));
			RooRealVar* nsig_specMC = ws->var(Form("nsig_specMC%d_%s",_count,""));
			hPtMC_2S->SetBinContent(i+1,nsig_specMC->getVal()/b_width);
			hPtMC_2S->SetBinError(i+1,nsig_specMC->getError()/b_width);
			hPtMC_2S->GetXaxis()->SetBinLabel(i+1, Form("%.2f", var_mean_av[i]));
		}
		//////////////////////////// FILL YIELD HISTOGRAMS

		////////////////////////////////////////////////////////// LABELS IN PLOTS
		// print fit meson NAME
		if (TREE == "ntKp" || TREE == "ntKstar" || TREE == "ntphi") {
			TString mesonLabel = "";
			if (TREE == "ntKp") mesonLabel = "#bf{B^{+}}";
			else if (TREE == "ntKstar") mesonLabel = "#bf{B^{0}}";
			else if (TREE == "ntphi") mesonLabel = "#bf{B_{s}^{0}}";
			TLatex* mesonName = new TLatex(0.18, 0.8, mesonLabel);
			setupLABELS(mesonName, 0.060, true);
		}

		// print fit xi2
		TLatex* chi_square = new TLatex(0.68, 0.55, Form("#chi^{2}/ndf = %.2f",chi2_vec[i]));
		setupLABELS(chi_square);

		// print var bin
		TString varLabel;
		if (VAR == "Bpt"){varLabel = "p_{T} [GeV/c]";}
		else if (VAR == "By"){varLabel = "|y|";}
		else if (VAR == "nSelectedChargedTracks"){varLabel = "n_{ch}";}
		else if (VAR == "CentBin"){varLabel = "Centrality (%)";}
		TLatex* varBIN = new TLatex(0.68, 0.50, Form("%d < %s < %d", (int)_varBINS[i], varLabel.Data(), (int)_varBINS[i+1]));
		setupLABELS(varBIN);

		// print Nsig in plot
		TLatex* N_signal = new TLatex(0.68, 0.45, Form("Y_{s} = %.0f #pm %.0f", yield_vec[i], yield_Stat_unc[i]));
		setupLABELS(N_signal);
		TLatex* N_signal_psi = nullptr;
		if (TREE == "ntmix") {
			N_signal_psi = new TLatex(0.68, 0.4, Form("Y_{s}^{#psi(2s)} = %.0f #pm %.0f", yield_vec_2S[i], yield_Stat_unc_2S[i]));
			setupLABELS(N_signal_psi);
		}

		TLatex* variationLabel = new TLatex(0.68, (TREE == "ntmix") ? 0.35 : 0.40, "");
		setupLABELS(variationLabel, 0.030, false);

		if( (TREE == "ntmix") || (TREE == "ntKstar") ){		//SIGNIFICANCE
			double signif = GetSignificance( ws, _count, mass, 2.0);
			TLatex *Signf = new TLatex(0.68, 0.3, Form("S = %.2f", signif));
			setupLABELS(Signf);
		}

		DrawCmsHeader(c, SYSTEM);
		c->Update();
		///////////////////////////////////////////////////////// /LABELS IN PLOTS


		// Save
		if(VAR == "By"){
			c->SaveAs(  Form("%s/data_%s_%s_%0.1f_%0.1f_",OUTPLOTF.Data(),SYSTEM.Data(),Form("abs(%s)",VAR.Data()), (float)_varBINS[i],(float)_varBINS[i+1])+TREE+".pdf");
			cMC->SaveAs(Form("%s/mc_%s_%s_%0.1f_%0.1f_"  ,OUTPLOTF.Data(),SYSTEM.Data(),Form("abs(%s)",VAR.Data()), (float)_varBINS[i],(float)_varBINS[i+1])+TREE+".pdf");
		}else{
			c->SaveAs(  Form("%s/data_%s_%s_%i_%i_",OUTPLOTF.Data(),SYSTEM.Data(),VAR.Data(),(int)_varBINS[i],(int)_varBINS[i+1])+TREE+".pdf");
			cMC->SaveAs(Form("%s/mc_%s_%s_%i_%i_"  ,OUTPLOTF.Data(),SYSTEM.Data(),VAR.Data(),(int)_varBINS[i],(int)_varBINS[i+1])+TREE+".pdf");
		}
		// Save

		vector<double> back_variation; 
		vector<double> back_unc; 
		vector<double> back_variation_2S;
		vector<double> back_unc_2S;
		vector<double> signal_variation; 
		vector<double> signal_unc; 
		vector<double> signal_variation_2S;
		vector<double> signal_unc_2S;
		vector<double> general_unc; 
		vector<double> general_unc_2S;
		double max_signal=0.; 
		double max_signal_2S=0.;
		double max_back=0.; 
		double max_back_2S=0.;

		if(syst_study==1 && FULL==0){
			//BACKGROUND MODEL SYSTEMATIC STUDY
			for(int j=0; j < static_cast<int>(background.size()); j++)
			{
				RooFitResult* f_back = fit("background", background[j].code.c_str(), TREE, c, cMC, ith_DATA_bin, ith_MC_bin, mass, _varBINS[i], _varBINS[i+1], *ws, VAR, nbinsmasshisto);
				RooRealVar* chi2_data_norm_back = ws->var(Form("chi2_data_norm%d_%s", _count, background[j].code.c_str()));
				RooRealVar* fitYield_back = static_cast<RooRealVar*>(f_back->floatParsFinal().at(f_back->floatParsFinal().index(Form("nsig%d_%s",_count,background[j].code.c_str()))));
				RooRealVar* fitYield_back_2S = nullptr;
				if (TREE == "ntmix") {
					fitYield_back_2S = static_cast<RooRealVar*>(f_back->floatParsFinal().at(f_back->floatParsFinal().index(Form("nsig_spec%d_%s",_count,background[j].code.c_str()))));
				}
				chi2_vec_back[j][i] = (chi2_data_norm_back ? chi2_data_norm_back->getVal() : -1.0);
				chi_square->SetText(0.68,0.55,Form("#chi^{2}/ndf = %.2f ",chi2_vec_back[j][i]));
				chi_square->Draw();
				varBIN->Draw();
				N_signal->SetTitle(Form("Y_{s} = %.0f #pm %.0f", fitYield_back->getVal(), fitYield_back->getError()));
				N_signal->Draw();
				if (N_signal_psi && fitYield_back_2S) {
					N_signal_psi->SetTitle(Form("Y_{s}^{#psi(2s)} = %.0f #pm %.0f", fitYield_back_2S->getVal(), fitYield_back_2S->getError()));
					N_signal_psi->Draw();
				}
				variationLabel->SetTitle(Form("(%s)", background[j].label.c_str()));
				variationLabel->Draw();
				c->Update();

				if(VAR == "By"){c->SaveAs(Form("%s/data_%s_%s_%0.1f_%0.1f_%s_", OUTPLOTF.Data(), SYSTEM.Data(), Form("abs(%s)",VAR.Data()),(float)_varBINS[i],(float)_varBINS[i+1],background[j].code.c_str())+TREE+ ".pdf");}
				else { c->SaveAs(Form("%s/data_%s_%s_%i_%i_%s_", OUTPLOTF.Data(), SYSTEM.Data(), VAR.Data(),(int)_varBINS[i],(int)_varBINS[i+1],background[j].code.c_str())+TREE+".pdf");}

				back_variation.push_back(fitYield_back->getVal());
				double back_rel_unc = (yield_vec[i] != 0.) ? abs(((yield_vec[i]-fitYield_back->getVal())/yield_vec[i])*100.) : 0.;
				back_unc.push_back(back_rel_unc);
				if(back_rel_unc > max_back) max_back = back_rel_unc;
				if (TREE == "ntmix" && fitYield_back_2S) {
					back_variation_2S.push_back(fitYield_back_2S->getVal());
					double back_rel_unc_2S = (yield_vec_2S[i] != 0.) ? abs(((yield_vec_2S[i]-fitYield_back_2S->getVal())/yield_vec_2S[i])*100.) : 0.;
					back_unc_2S.push_back(back_rel_unc_2S);
					if(back_rel_unc_2S > max_back_2S) max_back_2S = back_rel_unc_2S;
				}
			}
			general_unc.push_back(max_back);
			background_syst.push_back(back_variation);

			//SIGNAL MODEL SYSTEMATIC STUDY
			for(int j=0; j< static_cast<int>(signal.size()); j++)
			{
				RooFitResult* f_signal = fit("signal", signal[j].code.c_str(), TREE, c, cMC, ith_DATA_bin, ith_MC_bin, mass, _varBINS[i], _varBINS[i+1], *ws, VAR, nbinsmasshisto);
				RooRealVar* chi2_data_norm_sig = ws->var(Form("chi2_data_norm%d_%s", _count, signal[j].code.c_str()));
				RooRealVar* fitYield_signal = static_cast<RooRealVar*>(f_signal->floatParsFinal().at(f_signal->floatParsFinal().index(Form("nsig%d_%s",_count,signal[j].code.c_str()))));
				RooRealVar* fitYield_signal_2S = nullptr;
				if (TREE == "ntmix") {
					fitYield_signal_2S = static_cast<RooRealVar*>(f_signal->floatParsFinal().at(f_signal->floatParsFinal().index(Form("nsig_spec%d_%s",_count,signal[j].code.c_str()))));
				}
					chi2_vec_sig[j][i] = (chi2_data_norm_sig ? chi2_data_norm_sig->getVal() : -1.0);
					chi_square->SetText(0.68, 0.55, Form("#chi^{2}/ndf = %.2f ", chi2_vec_sig[j][i]));
					chi_square->Draw();
					varBIN->Draw();
					N_signal->SetTitle(Form("Y_{s} = %.0f #pm %.0f", fitYield_signal->getVal(), fitYield_signal->getError()));
					N_signal->Draw();
					if (N_signal_psi && fitYield_signal_2S) {
						N_signal_psi->SetTitle(Form("Y_{s}^{#psi(2s)} = %.0f #pm %.0f", fitYield_signal_2S->getVal(), fitYield_signal_2S->getError()));
						N_signal_psi->Draw();
					}
					variationLabel->SetTitle(Form("(%s)", signal[j].label.c_str()));
					variationLabel->Draw();
					c->Update();
				
				if (signal[j].code != "fixed") {
					if(VAR == "By"){ cMC->SaveAs(Form("%s/mc_%s_%s_%0.1f_%0.1f_%s_",OUTPLOTF.Data(),SYSTEM.Data(),Form("abs(%s)",VAR.Data()), (float)_varBINS[i], (float)_varBINS[i+1],signal[j].code.c_str())+TREE+".pdf");} 
					else { cMC->SaveAs(Form("%s/mc_%s_%s_%i_%i_%s_",OUTPLOTF.Data(),SYSTEM.Data(),VAR.Data(), (int)_varBINS[i], (int)_varBINS[i+1],signal[j].code.c_str() )+TREE+".pdf");}
				}
				if(VAR == "By"){ c->SaveAs(Form("%s/data_%s_%s_%0.1f_%0.1f_%s_",OUTPLOTF.Data(),SYSTEM.Data(),Form("abs(%s)",VAR.Data()),(float)_varBINS[i],(float)_varBINS[i+1],signal[j].code.c_str() )+TREE+".pdf");}
				else{ c->SaveAs(Form("%s/data_%s_%s_%i_%i_%s_",OUTPLOTF.Data(),SYSTEM.Data(),VAR.Data(),(int)_varBINS[i],(int)_varBINS[i+1],signal[j].code.c_str() )+TREE+".pdf");}
				
				signal_variation.push_back(fitYield_signal->getVal());
				double signal_rel_unc = (yield_vec[i] != 0.) ? abs(((yield_vec[i]-fitYield_signal->getVal())/yield_vec[i])*100.) : 0.;
				signal_unc.push_back(signal_rel_unc);
				if(signal_rel_unc > max_signal) max_signal = signal_rel_unc;
				if (TREE == "ntmix" && fitYield_signal_2S) {
					signal_variation_2S.push_back(fitYield_signal_2S->getVal());
					double signal_rel_unc_2S = (yield_vec_2S[i] != 0.) ? abs(((yield_vec_2S[i]-fitYield_signal_2S->getVal())/yield_vec_2S[i])*100.) : 0.;
					signal_unc_2S.push_back(signal_rel_unc_2S);
					if(signal_rel_unc_2S > max_signal_2S) max_signal_2S = signal_rel_unc_2S;
				}
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
			if (TREE == "ntmix") {
				general_unc_2S.push_back(max_back_2S);
				general_unc_2S.push_back(max_signal_2S);
				general_unc_2S.push_back(sqrt(max_back_2S*max_back_2S+max_signal_2S*max_signal_2S));
				stat_error_2S.push_back(stat_un_2S);
				back_syst_rel_values_2S.push_back(back_unc_2S);
				sig_syst_rel_values_2S.push_back(signal_unc_2S);
				general_syst_2S.push_back(general_unc_2S);
				yield_vec_2S_systerr_low[i] = general_unc_2S[2] / 100 * yield_vec_2S[i];
				yield_vec_2S_systerr_high[i] = general_unc_2S[2] / 100 * yield_vec_2S[i];
			}
		}

		//VALIDATION STUDIES
		if (val==1 && syst_study==0){
			string Path_val=Form("./%s/validation",OUTPLOTF.Data());
			validate_fit(ws, "", TREE, VAR, FULL, _varBINS[i], _varBINS[i+1],Path_val);
		}
		//VALIDATION STUDIES
	}
	//BIN ANALYSIS END
	//BIN ANALYSIS END
	
	// Save yields in ROOT file
	TFile* outf = new TFile(Form("ROOTfiles/yields_%s_%s_%s.root",TREE.Data(),VAR.Data(), SYSTEM.Data()),"recreate");
	outf->cd();
	hPt->Write();	
	hPtMC->Write();
	if (TREE == "ntmix"){
		hPt_2S->Write();
		hPtMC_2S->Write();
	}
	outf->Close();

	// for sPlot purpose
	if (FULL==1){ 
		TFile* nominalModelOut = new TFile(Form("ROOTfiles/nominalFitModel_%s_%s.root", TREE.Data(), SYSTEM.Data()), "recreate");
		nominalModelOut->cd();
		ws->Write("ws_nominal");
		nominalModelOut->Close();
	}

	//systematic tables!
	if(syst_study==1 && FULL==0){
		vector<string> labels_back;
		vector<string> col_name_back;
		vector<string> labels_signal;
		vector<string> labels_general = GetGeneralSystematicLabels();
		vector<string> col_name_general;
		vector<string> col_name_signal;
		col_name_back.push_back("Background Model");
		col_name_signal.push_back("Signal Model");
		col_name_general.push_back("Systematic Source");
		for (const auto& variation : background) labels_back.push_back(variation.label);
		for (const auto& variation : signal) labels_signal.push_back(variation.label);
		for(int i=0;i<_nBins;i++){
			string label1 = GetSystematicColumnLabel(VAR, _varBINS[i], _varBINS[i+1]);
			col_name_back.push_back(label1);
			col_name_signal.push_back(label1);
			col_name_general.push_back(label1);
		}
		WriteSystematicsTablesDocument(
			"./files/systematics_tables_" + string(VAR.Data()) + "_" + string(TREE.Data()),
			col_name_signal, col_name_back, col_name_general, labels_signal, labels_back, labels_general, sig_syst_rel_values,
			back_syst_rel_values, BuildGeneralSystematicNumbers(general_syst, stat_error)
		);
		if (TREE == "ntmix"){
			WriteSystematicsTablesDocument(
				"./files/systematics_tables_" + string(VAR.Data()) + "_" + string(TREE.Data()) + "_psi2s",
				col_name_signal, col_name_back, col_name_general, labels_signal, labels_back, labels_general,
				sig_syst_rel_values_2S, back_syst_rel_values_2S, BuildGeneralSystematicNumbers(general_syst_2S, stat_error_2S)
			);
		}

		double zero[_nBins];
		for (int i=0;i<_nBins;i++){zero[i]=0.;}
		double low_high_b[_nBins];
		Double_t x[_nBins];
		for (int i=0;i<_nBins;i++){
			x[i] = (_varBINS[i]+_varBINS[i+1])/2;
			low_high_b[i] = _varBINS[i+1] - x[i];
		}

		TGraph *binning= new TGraphAsymmErrors (_nBins,x,zero,low_high_b,low_high_b,zero,zero);
		binning->SetMarkerColorAlpha(kBlack, 0);
		binning->SetLineWidth(6);

		TMultiGraph* m_back_sig= new TMultiGraph();
		TLegend *legsigback = new TLegend(0.75,0.71,0.89,0.89,NULL,"brNDC");
		legsigback->SetBorderSize(0);
		legsigback->SetTextSize(0.035);
		legsigback->SetTextFont(42);
		legsigback->SetFillStyle(0);
		m_back_sig->Add(binning);

		TCanvas* c_back= new TCanvas("c_back","",700,700);
		TLegend *legback = new TLegend(0.8,0.75,0.89,0.89,NULL,"brNDC");
		legback->SetBorderSize(0);
		legback->SetTextSize(0.035);
		legback->SetTextFont(42);
		legback->SetFillStyle(0);
		TMultiGraph* m_back= new TMultiGraph();
		double y_max_back=0;
		for (int j=0;j<(int)(background.size());j++){
			Double_t y[_nBins];
			for (int i=0;i<_nBins;i++){
				y[i] = back_syst_rel_values[i][j];
				if (y[i]>y_max_back){y_max_back=y[i];}
			}
			TGraphAsymmErrors *g_back= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
			g_back->SetMarkerColor(j+1);
			g_back->SetLineColor(j+1);
			g_back->SetMarkerStyle(22);
			g_back->SetMarkerSize(1.2);
			g_back->SetLineWidth(2);
			m_back->Add(g_back,"P");
			m_back_sig->Add(g_back,"P");
			legback->AddEntry(g_back,background[j].label.c_str(),"p");
			legsigback->AddEntry(g_back,background[j].label.c_str(),"p");
		}
		if (y_max_back <= 0.) y_max_back = 1.;
		m_back->Add(binning);
		m_back->GetXaxis()->SetTitle("p_{T}");
		m_back->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
		m_back->GetYaxis()->SetRangeUser(0, y_max_back*1.5);
		m_back->Draw("AE1");
		if(VAR == "By"){
			m_back->GetXaxis()->SetTitle("Rapidity (y)");
			m_back->GetXaxis()->SetLimits(0,2.4);
		}
		else if(VAR == "Bpt"){
			m_back->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			if (TREE != "ntmix"){ m_back->GetXaxis()->SetLimits(0,65); }
		}
		else if(VAR == "nSelectedChargedTracks"){
			m_back->GetXaxis()->SetTitle("Multiplicity (Mult)");
			m_back->GetXaxis()->SetLimits(0,110);
		}
		legback->Draw();
		c_back->SaveAs(Form("./results/tables/background_systematics_plot_%s_%s.pdf",TREE.Data(),VAR.Data()));

		TCanvas* c_sig= new TCanvas("c_sig","",700,700);
		TLegend *legsig = new TLegend(0.75,0.8,0.89,0.89,NULL,"brNDC");
		legsig->SetBorderSize(0);
		legsig->SetTextSize(0.035);
		legsig->SetTextFont(42);
		legsig->SetFillStyle(0);
		TMultiGraph* m_sig= new TMultiGraph();
		double y_max_sig=0;
		for (int j=0;j<(int)(signal.size());j++){
			Double_t y[_nBins];
			for (int i=0;i<_nBins;i++){
				y[i] = sig_syst_rel_values[i][j];
				if (y[i]>y_max_sig){y_max_sig=y[i];}
			}
			TGraphAsymmErrors *g_sig= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
			g_sig->SetMarkerColor(j+5);
			g_sig->SetLineColor(j+5);
			g_sig->SetMarkerStyle(21);
			g_sig->SetMarkerSize(1.2);
			g_sig->SetLineWidth(2);
			m_sig->Add(g_sig,"P");
			m_back_sig->Add(g_sig,"P");
			legsig->AddEntry(g_sig,signal[j].label.c_str(),"p");
			legsigback->AddEntry(g_sig,signal[j].label.c_str(),"p");
		}
		if (y_max_sig <= 0.) y_max_sig = 1.;
		m_sig->Add(binning);
		m_sig->GetXaxis()->SetTitle("p_{T}");
		m_sig->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
		m_sig->GetYaxis()->SetRangeUser(0, y_max_sig*1.5);
		m_sig->Draw("AE1");
		if(VAR == "By"){
			m_sig->GetXaxis()->SetTitle("Rapidity (y)");
			m_sig->GetXaxis()->SetLimits(0,2.4);
		}
		else if(VAR == "Bpt"){
			m_sig->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			m_sig->GetXaxis()->SetLimits(0,80);
		}
		else if(VAR == "nSelectedChargedTracks"){
			m_sig->GetXaxis()->SetTitle("Multiplicity (Mult)");
			m_sig->GetXaxis()->SetLimits(0,110);
		}
		legsig->Draw();
		c_sig->SaveAs(Form("./results/tables/signal_systematics_plot_%s_%s.pdf",TREE.Data(),VAR.Data()));

		TCanvas *c_sig_back= new TCanvas("c_sig_back","",700,700);
		double y_max_combined = (y_max_back > y_max_sig) ? y_max_back : y_max_sig;
		if (y_max_combined <= 0.) y_max_combined = 1.;
		m_back_sig->GetYaxis()->SetRangeUser(0, y_max_combined*1.5);
		if(VAR == "By"){
			m_back_sig->GetXaxis()->SetTitle("Rapidity (y)");
			m_back_sig->GetXaxis()->SetLimits(0,2.4);
		}
		else if(VAR == "Bpt"){
			m_back_sig->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			if (TREE != "ntmix"){ m_back_sig->GetXaxis()->SetLimits(0,65); }
		}
		else if(VAR == "nSelectedChargedTracks"){
			m_back_sig->GetXaxis()->SetTitle("Multiplicity (Mult)");
			m_back_sig->GetXaxis()->SetLimits(0,110);
		}
		m_back_sig->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
		m_back_sig->Draw("AE1");
		legsigback->Draw();
		c_sig_back->SaveAs(Form("./results/tables/background_signal_systematics_plot_%s_%s.pdf",TREE.Data(),VAR.Data()));

		TCanvas* c_gen= new TCanvas("c_gen","",700,700);
		TLegend *legen = new TLegend(0.8,0.77,0.89,0.89,NULL,"brNDC");
		legen->SetBorderSize(0);
		legen->SetTextSize(0.035);
		legen->SetTextFont(42);
		legen->SetFillStyle(0);
		TMultiGraph* m_gen= new TMultiGraph();
		const char* genlabel[3]={"Background PDF", "Signal PDF", "Total Systematic"};
		double y_max_gen=0;
		for (int j=0;j<3;j++){
			Double_t y[_nBins];
			for (int i=0;i<_nBins;i++){
				y[i] = general_syst[i][j];
				if (y[i]>y_max_gen){y_max_gen=y[i];}
			}
			TGraphAsymmErrors *g_gen= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
			g_gen->SetMarkerColor(j+1);
			g_gen->SetLineColor(j+1);
			g_gen->SetMarkerStyle(21);
			g_gen->SetMarkerSize(1.2);
			g_gen->SetLineWidth(2);
			m_gen->Add(g_gen,"P");
			legen->AddEntry(g_gen,genlabel[j],"p");
		}
		Double_t y[_nBins];
		for (int i=0;i<_nBins;i++){y[i]=stat_error[i][0]; if (y[i]>y_max_gen){y_max_gen=y[i];}}
		TGraphAsymmErrors *g_gen= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
		m_gen->Add(binning);
		g_gen->SetMarkerColor(9);
		g_gen->SetLineColor(9);
		g_gen->SetMarkerStyle(21);
		g_gen->SetMarkerSize(1.2);
		g_gen->SetLineWidth(2);
		m_gen->Add(g_gen,"P");
		legen->AddEntry(g_gen,"Statistical","p");
		if (y_max_gen <= 0.) y_max_gen = 1.;
		m_gen->GetXaxis()->SetTitle("p_{T}");
		m_gen->GetYaxis()->SetTitle("Total Uncertainty(%)");
		m_gen->GetYaxis()->SetRangeUser(0, y_max_gen*1.5);
		m_gen->Draw("AE1");
		if(VAR == "By"){
			m_gen->GetXaxis()->SetTitle("Rapidity (y)");
			m_gen->GetXaxis()->SetLimits(0,2.4);
		}
		else if(VAR == "Bpt"){
			m_gen->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			m_gen->GetXaxis()->SetLimits(0,80);
		}
		else if(VAR == "nSelectedChargedTracks"){
			m_gen->GetXaxis()->SetTitle("Multiplicity (Mult)");
			m_gen->GetXaxis()->SetLimits(0,110);
		}
		legen->Draw();
		c_gen->SaveAs(Form("./results/tables/general_systematics_plot_%s_%s.pdf",TREE.Data(),VAR.Data()));

		if (TREE == "ntmix") {
			TCanvas* c_back_2S= new TCanvas("c_back_2S","",700,700);
			TLegend *legback_2S = new TLegend(0.8,0.75,0.89,0.89,NULL,"brNDC");
			legback_2S->SetBorderSize(0);
			legback_2S->SetTextSize(0.035);
			legback_2S->SetTextFont(42);
			legback_2S->SetFillStyle(0);
			TMultiGraph* m_back_2S= new TMultiGraph();
			double y_max_back_2S=0;
			for (int j=0;j<(int)(background.size());j++){
				Double_t ypsi[_nBins];
				for (int i=0;i<_nBins;i++){
					ypsi[i] = back_syst_rel_values_2S[i][j];
					if (ypsi[i]>y_max_back_2S){y_max_back_2S=ypsi[i];}
				}
				TGraphAsymmErrors *g_back_2S= new TGraphAsymmErrors (_nBins,x,ypsi,zero,zero,zero,zero);
				g_back_2S->SetMarkerColor(j+1);
				g_back_2S->SetLineColor(j+1);
				g_back_2S->SetMarkerStyle(22);
				g_back_2S->SetMarkerSize(1.2);
				g_back_2S->SetLineWidth(2);
				m_back_2S->Add(g_back_2S,"P");
				legback_2S->AddEntry(g_back_2S,background[j].label.c_str(),"p");
			}
			if (y_max_back_2S <= 0.) y_max_back_2S = 1.;
			m_back_2S->Add(binning);
			m_back_2S->GetXaxis()->SetTitle("p_{T}");
			m_back_2S->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
			m_back_2S->GetYaxis()->SetRangeUser(0, y_max_back_2S*1.5);
			m_back_2S->Draw("AE1");
			if(VAR == "By"){
				m_back_2S->GetXaxis()->SetTitle("Rapidity (y)");
				m_back_2S->GetXaxis()->SetLimits(0,2.4);
			}
			else if(VAR == "Bpt"){
				m_back_2S->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			}
			else if(VAR == "nSelectedChargedTracks"){
				m_back_2S->GetXaxis()->SetTitle("Multiplicity (Mult)");
				m_back_2S->GetXaxis()->SetLimits(0,110);
			}
			legback_2S->Draw();
			c_back_2S->SaveAs(Form("./results/tables/background_systematics_plot_%s_%s_psi2s.pdf",TREE.Data(),VAR.Data()));

			TCanvas* c_sig_2S= new TCanvas("c_sig_2S","",700,700);
			TLegend *legsig_2S = new TLegend(0.75,0.8,0.89,0.89,NULL,"brNDC");
			legsig_2S->SetBorderSize(0);
			legsig_2S->SetTextSize(0.035);
			legsig_2S->SetTextFont(42);
			legsig_2S->SetFillStyle(0);
			TMultiGraph* m_sig_2S= new TMultiGraph();
			double y_max_sig_2S=0;
			for (int j=0;j<(int)(signal.size());j++){
				Double_t ypsi[_nBins];
				for (int i=0;i<_nBins;i++){
					ypsi[i] = sig_syst_rel_values_2S[i][j];
					if (ypsi[i]>y_max_sig_2S){y_max_sig_2S=ypsi[i];}
				}
				TGraphAsymmErrors *g_sig_2S= new TGraphAsymmErrors (_nBins,x,ypsi,zero,zero,zero,zero);
				g_sig_2S->SetMarkerColor(j+5);
				g_sig_2S->SetLineColor(j+5);
				g_sig_2S->SetMarkerStyle(21);
				g_sig_2S->SetMarkerSize(1.2);
				g_sig_2S->SetLineWidth(2);
				m_sig_2S->Add(g_sig_2S,"P");
				legsig_2S->AddEntry(g_sig_2S,signal[j].label.c_str(),"p");
			}
			if (y_max_sig_2S <= 0.) y_max_sig_2S = 1.;
			m_sig_2S->Add(binning);
			m_sig_2S->GetXaxis()->SetTitle("p_{T}");
			m_sig_2S->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
			m_sig_2S->GetYaxis()->SetRangeUser(0, y_max_sig_2S*1.5);
			m_sig_2S->Draw("AE1");
			if(VAR == "By"){
				m_sig_2S->GetXaxis()->SetTitle("Rapidity (y)");
				m_sig_2S->GetXaxis()->SetLimits(0,2.4);
			}
			else if(VAR == "Bpt"){
				m_sig_2S->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
				m_sig_2S->GetXaxis()->SetLimits(0,80);
			}
			else if(VAR == "nSelectedChargedTracks"){
				m_sig_2S->GetXaxis()->SetTitle("Multiplicity (Mult)");
				m_sig_2S->GetXaxis()->SetLimits(0,110);
			}
			legsig_2S->Draw();
			c_sig_2S->SaveAs(Form("./results/tables/signal_systematics_plot_%s_%s_psi2s.pdf",TREE.Data(),VAR.Data()));

			TCanvas *c_sig_back_2S= new TCanvas("c_sig_back_2S","",700,700);
			TMultiGraph* m_back_sig_2S= new TMultiGraph();
			TLegend *legsigback_2S = new TLegend(0.75,0.71,0.89,0.89,NULL,"brNDC");
			legsigback_2S->SetBorderSize(0);
			legsigback_2S->SetTextSize(0.035);
			legsigback_2S->SetTextFont(42);
			legsigback_2S->SetFillStyle(0);
			m_back_sig_2S->Add(binning);
			double y_max_combined_2S = 0.;
			for (int j=0;j<(int)(background.size());j++){
				Double_t ypsi[_nBins];
				for (int i=0;i<_nBins;i++){
					ypsi[i] = back_syst_rel_values_2S[i][j];
					if (ypsi[i]>y_max_combined_2S){y_max_combined_2S=ypsi[i];}
				}
				TGraphAsymmErrors *g_back_2S= new TGraphAsymmErrors (_nBins,x,ypsi,zero,zero,zero,zero);
				g_back_2S->SetMarkerColor(j+1);
				g_back_2S->SetLineColor(j+1);
				g_back_2S->SetMarkerStyle(22);
				g_back_2S->SetMarkerSize(1.2);
				g_back_2S->SetLineWidth(2);
				m_back_sig_2S->Add(g_back_2S,"P");
				legsigback_2S->AddEntry(g_back_2S,background[j].label.c_str(),"p");
			}
			for (int j=0;j<(int)(signal.size());j++){
				Double_t ypsi[_nBins];
				for (int i=0;i<_nBins;i++){
					ypsi[i] = sig_syst_rel_values_2S[i][j];
					if (ypsi[i]>y_max_combined_2S){y_max_combined_2S=ypsi[i];}
				}
				TGraphAsymmErrors *g_sig_2S= new TGraphAsymmErrors (_nBins,x,ypsi,zero,zero,zero,zero);
				g_sig_2S->SetMarkerColor(j+5);
				g_sig_2S->SetLineColor(j+5);
				g_sig_2S->SetMarkerStyle(21);
				g_sig_2S->SetMarkerSize(1.2);
				g_sig_2S->SetLineWidth(2);
				m_back_sig_2S->Add(g_sig_2S,"P");
				legsigback_2S->AddEntry(g_sig_2S,signal[j].label.c_str(),"p");
			}
			if (y_max_combined_2S <= 0.) y_max_combined_2S = 1.;
			m_back_sig_2S->GetYaxis()->SetRangeUser(0, y_max_combined_2S*1.5);
			if(VAR == "By"){
				m_back_sig_2S->GetXaxis()->SetTitle("Rapidity (y)");
				m_back_sig_2S->GetXaxis()->SetLimits(0,2.4);
			}
			else if(VAR == "Bpt"){
				m_back_sig_2S->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			}
			else if(VAR == "nSelectedChargedTracks"){
				m_back_sig_2S->GetXaxis()->SetTitle("Multiplicity (Mult)");
				m_back_sig_2S->GetXaxis()->SetLimits(0,110);
			}
			m_back_sig_2S->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
			m_back_sig_2S->Draw("AE1");
			legsigback_2S->Draw();
			c_sig_back_2S->SaveAs(Form("./results/tables/background_signal_systematics_plot_%s_%s_psi2s.pdf",TREE.Data(),VAR.Data()));

			TCanvas* c_gen_2S= new TCanvas("c_gen_2S","",700,700);
			TLegend *legen_2S = new TLegend(0.8,0.77,0.89,0.89,NULL,"brNDC");
			legen_2S->SetBorderSize(0);
			legen_2S->SetTextSize(0.035);
			legen_2S->SetTextFont(42);
			legen_2S->SetFillStyle(0);
			TMultiGraph* m_gen_2S= new TMultiGraph();
			double y_max_gen_2S=0;
			for (int j=0;j<3;j++){
				Double_t ypsi[_nBins];
				for (int i=0;i<_nBins;i++){
					ypsi[i] = general_syst_2S[i][j];
					if (ypsi[i]>y_max_gen_2S){y_max_gen_2S=ypsi[i];}
				}
				TGraphAsymmErrors *g_gen_2S= new TGraphAsymmErrors (_nBins,x,ypsi,zero,zero,zero,zero);
				g_gen_2S->SetMarkerColor(j+1);
				g_gen_2S->SetLineColor(j+1);
				g_gen_2S->SetMarkerStyle(21);
				g_gen_2S->SetMarkerSize(1.2);
				g_gen_2S->SetLineWidth(2);
				m_gen_2S->Add(g_gen_2S,"P");
				legen_2S->AddEntry(g_gen_2S,genlabel[j],"p");
			}
			Double_t ystatpsi[_nBins];
			for (int i=0;i<_nBins;i++){ystatpsi[i]=stat_error_2S[i][0]; if (ystatpsi[i]>y_max_gen_2S){y_max_gen_2S=ystatpsi[i];}}
			TGraphAsymmErrors *g_gen_2S= new TGraphAsymmErrors (_nBins,x,ystatpsi,zero,zero,zero,zero);
			m_gen_2S->Add(binning);
			g_gen_2S->SetMarkerColor(9);
			g_gen_2S->SetLineColor(9);
			g_gen_2S->SetMarkerStyle(21);
			g_gen_2S->SetMarkerSize(1.2);
			g_gen_2S->SetLineWidth(2);
			m_gen_2S->Add(g_gen_2S,"P");
			legen_2S->AddEntry(g_gen_2S,"Statistical","p");
			if (y_max_gen_2S <= 0.) y_max_gen_2S = 1.;
			m_gen_2S->GetXaxis()->SetTitle("p_{T}");
			m_gen_2S->GetYaxis()->SetTitle("Total Uncertainty(%)");
			m_gen_2S->GetYaxis()->SetRangeUser(0, y_max_gen_2S*1.5);
			m_gen_2S->Draw("AE1");
			if(VAR == "By"){
				m_gen_2S->GetXaxis()->SetTitle("Rapidity (y)");
				m_gen_2S->GetXaxis()->SetLimits(0,2.4);
			}
			else if(VAR == "Bpt"){
				m_gen_2S->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
				m_gen_2S->GetXaxis()->SetLimits(0,80);
			}
			else if(VAR == "nSelectedChargedTracks"){
				m_gen_2S->GetXaxis()->SetTitle("Multiplicity (Mult)");
				m_gen_2S->GetXaxis()->SetLimits(0,110);
			}
			legen_2S->Draw();
			c_gen_2S->SaveAs(Form("./results/tables/general_systematics_plot_%s_%s_psi2s.pdf",TREE.Data(),VAR.Data()));
		}
	}


	//Differential plot part
	TCanvas c_diff;
	TMultiGraph* mg = new TMultiGraph();
	TLegend *leg_d = new TLegend(0.7,0.7,0.9,0.9);
	TGraphAsymmErrors* gr_staterr = new TGraphAsymmErrors(_nBins,var_mean_av,yield_vec,hori_av_low,hori_av_high,yield_Stat_unc,yield_Stat_unc);
	gr_staterr->SetLineColor(1);
	gr_staterr->SetMarkerColor(1);
	gr_staterr->SetMarkerStyle(20);
	gr_staterr->SetLineWidth(2);
	gr_staterr->SetName("Y_stat");
	mg->Add(gr_staterr, "PE");
	leg_d->AddEntry(gr_staterr, "Statistical Uncertainty", "e");
	if(syst_study==1){
		TGraphAsymmErrors* gr_systerr = new TGraphAsymmErrors(_nBins, var_mean_av, yield_vec, nullptr, nullptr, yield_vec_systerr_low, yield_vec_systerr_high);
		gr_systerr->SetLineColor(2);
		gr_systerr->SetMarkerColor(2);
		gr_systerr->SetLineWidth(2);
	 	gr_systerr->SetName("Y_syst");
		mg->Add(gr_systerr,"E");
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
	} else if(VAR == "nSelectedChargedTracks"){
		mg->GetXaxis()->SetTitle("Multiplicity (Mult)");
		mg->GetYaxis()->SetTitle("dY_{S}/dMult");
		mg->GetXaxis()->SetLimits(0, 110);
	}
	mg->Draw("ap");
	leg_d->SetBorderSize(0);
	leg_d->SetFillStyle(0);
	leg_d->SetTextFont(42);
	leg_d->SetTextSize(0.035);
	leg_d->Draw();
	const char* pathc =Form("./results/Graphs/raw_yield_%s_%s.pdf", TREE.Data(), VAR.Data());
	c_diff.SaveAs(pathc);

	if (TREE == "ntmix") {
		TCanvas c_diff_2S;
		TMultiGraph* mg_2S = new TMultiGraph();
		TLegend *leg_d_2S = new TLegend(0.68,0.78,0.90,0.89,NULL,"brNDC");
		TGraphAsymmErrors* gr_staterr_2S = new TGraphAsymmErrors(_nBins,var_mean_av,yield_vec_2S,hori_av_low,hori_av_high,yield_Stat_unc_2S,yield_Stat_unc_2S);
		gr_staterr_2S->SetLineColor(1);
		gr_staterr_2S->SetMarkerColor(1);
		gr_staterr_2S->SetMarkerStyle(20);
		gr_staterr_2S->SetLineWidth(2);
		mg_2S->Add(gr_staterr_2S, "PE");
		leg_d_2S->AddEntry(gr_staterr_2S, "Statistical Uncertainty", "e");
		if(syst_study==1){
			TGraphAsymmErrors* gr_systerr_2S = new TGraphAsymmErrors(_nBins, var_mean_av, yield_vec_2S, nullptr, nullptr, yield_vec_2S_systerr_low, yield_vec_2S_systerr_high);
			gr_systerr_2S->SetLineColor(2);
			gr_systerr_2S->SetMarkerColor(2);
			gr_systerr_2S->SetLineWidth(2);
			mg_2S->Add(gr_systerr_2S,"E");
			leg_d_2S->AddEntry(gr_systerr_2S, "Systematic Uncertainty", "e");
		}
		if(VAR == "By"){
			mg_2S->GetXaxis()->SetTitle("Rapidity (y)");
			mg_2S->GetYaxis()->SetTitle("dY_{S}/dy");
			mg_2S->GetXaxis()->SetLimits(0,2.4);
		} else if(VAR == "Bpt"){
			mg_2S->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			mg_2S->GetYaxis()->SetTitle("dY_{S}/dp_{T}");
			mg_2S->GetXaxis()->SetLimits(0 ,80);
		} else if(VAR == "nSelectedChargedTracks"){
			mg_2S->GetXaxis()->SetTitle("Multiplicity (Mult)");
			mg_2S->GetYaxis()->SetTitle("dY_{S}/dMult");
			mg_2S->GetXaxis()->SetLimits(0, 110);
		}
		mg_2S->Draw("AP");
		leg_d_2S->SetBorderSize(0);
		leg_d_2S->SetFillStyle(0);
		leg_d_2S->SetTextFont(42);
		leg_d_2S->SetTextSize(0.035);
		leg_d_2S->Draw();
		const char* pathc_2S =Form("./results/Graphs/raw_yield_%s_psi2s_%s.pdf", TREE.Data(), VAR.Data());
		c_diff_2S.SaveAs(pathc_2S);
	}
	//Differential plot part ends

	//Scale part starts
	double scale_max = 0;
	for(int i = 0; i < _nBins; i++){
		if(scale_vec[i] > scale_max){scale_max = scale_vec[i];}
		if(TREE == "ntmix" && scale_vec_2S[i] > scale_max){scale_max = scale_vec_2S[i];}
	}

	TCanvas c_par;
	TMultiGraph* mg_par = new TMultiGraph();
	TGraphAsymmErrors* gr_scale = new TGraphAsymmErrors(_nBins,var_mean_av,scale_vec,hori_av_low,hori_av_high,scale_vec_unc,scale_vec_unc);
	gr_scale->SetLineColor(TREE == "ntmix" ? kOrange-3 : 1);
	gr_scale->SetMarkerColor(TREE == "ntmix" ? kOrange-3 : 1);
	gr_scale->SetMarkerStyle(20);
	gr_scale->SetLineWidth(2);
	TLegend *leg_scale = nullptr;
	if (TREE == "ntmix") {
		leg_scale = new TLegend(0.66,0.78,0.90,0.89,NULL,"brNDC");
		leg_scale->SetBorderSize(0);
		leg_scale->SetFillStyle(0);
		leg_scale->SetTextFont(42);
		leg_scale->SetTextSize(0.035);
	}

	if(VAR == "By"){
		mg_par->GetXaxis()->SetTitle("Rapidity (y)");
		mg_par->GetXaxis()->SetLimits(0 ,2.4);
	} else if(VAR == "Bpt"){
		mg_par->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		mg_par->GetXaxis()->SetLimits(0 ,80); 
	} else if(VAR == "nSelectedChargedTracks"){
		mg_par->GetXaxis()->SetTitle("Multiplicity (Mult)");
		mg_par->GetXaxis()->SetLimits(0, 110);
	}
	mg_par->GetYaxis()->SetTitle("Scale factor");
	mg_par->Add(gr_scale,"PE");
	if (TREE == "ntmix") {
		TGraphAsymmErrors* gr_scale_2S = new TGraphAsymmErrors(_nBins,var_mean_av,scale_vec_2S,hori_av_low,hori_av_high,scale_vec_2S_unc,scale_vec_2S_unc);
		gr_scale_2S->SetLineColor(kOrange-2);
		gr_scale_2S->SetMarkerColor(kOrange-2);
		gr_scale_2S->SetMarkerStyle(21);
		gr_scale_2S->SetLineWidth(2);
		mg_par->Add(gr_scale_2S,"PE");
		leg_scale->AddEntry(gr_scale, "X(3872)", "lp");
		leg_scale->AddEntry(gr_scale_2S, "#psi(2S)", "lp");
	}
	mg_par->GetYaxis()->SetRangeUser(0,scale_max*1.4);
	mg_par->Draw("AP");
	if (leg_scale) leg_scale->Draw();
	const char* pathc_par =Form("./results/Graphs/scale_variation_%s_%s.pdf",TREE.Data(),VAR.Data()); 
	c_par.SaveAs(pathc_par);
	//Scale part ends

	//Resolution plot part
	double resol_max = 0;
	for(int i = 0; i < _nBins; i++){
		if(resol_vec[i] > resol_max){resol_max = resol_vec[i];}
		if(TREE == "ntmix" && resol_vec_2S[i] > resol_max){resol_max = resol_vec_2S[i];}
	}
	TCanvas c_resol;
	TMultiGraph* mg_resol = new TMultiGraph();
	TGraphAsymmErrors* gr_resol = new TGraphAsymmErrors(_nBins, var_mean_av, resol_vec, hori_av_low, hori_av_high, resol_vec_unc, resol_vec_unc);
	gr_resol->SetLineColor(TREE == "ntmix" ? kOrange-3 : 1);
	gr_resol->SetMarkerColor(TREE == "ntmix" ? kOrange-3 : 1);
	gr_resol->SetMarkerStyle(20);
	gr_resol->SetLineWidth(2);
	TLegend *leg_resol = nullptr;
	if (TREE == "ntmix") {
		leg_resol = new TLegend(0.66,0.78,0.90,0.89,NULL,"brNDC");
		leg_resol->SetBorderSize(0);
		leg_resol->SetFillStyle(0);
		leg_resol->SetTextFont(42);
		leg_resol->SetTextSize(0.035);
	}
	mg_resol->GetYaxis()->SetTitle("Resolution");
	if(VAR == "By"){
		mg_resol->GetXaxis()->SetTitle("Rapidity (y)");
		mg_resol->GetXaxis()->SetLimits(0 ,2.4);
	} else if(VAR == "Bpt"){
		mg_resol->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		mg_resol->GetXaxis()->SetLimits(0, 100); 
	} else if(VAR == "nSelectedChargedTracks"){
		mg_resol->GetXaxis()->SetTitle("Multiplicity (Mult)");
		mg_resol->GetXaxis()->SetLimits(0, 110);
	}
	mg_resol->GetYaxis()->SetRangeUser(0,(resol_max > 0. ? resol_max*1.4 : 0.2));
	mg_resol->Add(gr_resol,"PE");
	if (TREE == "ntmix") {
		TGraphAsymmErrors* gr_resol_2S = new TGraphAsymmErrors(_nBins, var_mean_av, resol_vec_2S, hori_av_low, hori_av_high, resol_vec_2S_unc, resol_vec_2S_unc);
		gr_resol_2S->SetLineColor(kOrange-2);
		gr_resol_2S->SetMarkerColor(kOrange-2);
		gr_resol_2S->SetMarkerStyle(21);
		gr_resol_2S->SetLineWidth(2);
		mg_resol->Add(gr_resol_2S,"PE");
		leg_resol->AddEntry(gr_resol, "X(3872)", "lp");
		leg_resol->AddEntry(gr_resol_2S, "#psi(2S)", "lp");
	}
	mg_resol->Draw("AP");
	if (leg_resol) leg_resol->Draw();
	const char* pathc_resol =Form("./results/Graphs/resolution_%s_%s.pdf",TREE.Data(),VAR.Data()); 
	c_resol.SaveAs(pathc_resol);
	//Resolution plot part ends

	//Chi2 plot part 
	double chi2_max = 0;
	for(int i = 0; i < _nBins; i++){
		if(chi2_vec[i] > chi2_max){chi2_max = chi2_vec[i];}
		if(chi2MC_vec[i] > chi2_max){chi2_max = chi2MC_vec[i];}
		if(TREE == "ntmix" && chi2MC_PSI2S_vec[i] > chi2_max){chi2_max = chi2MC_PSI2S_vec[i];}
	}
	TCanvas c_chi2;
	TMultiGraph* mg_chi2 = new TMultiGraph();
	TGraphAsymmErrors* gr_chi2   = new TGraphAsymmErrors(_nBins,var_mean_av,chi2_vec,hori_av_low,hori_av_high,nullptr,nullptr);
	gr_chi2->SetLineColor(kRed+1);
	gr_chi2->SetMarkerColor(kRed+1);
	gr_chi2->SetMarkerStyle(20);
	gr_chi2->SetLineWidth(2);
	TGraphAsymmErrors* grMC_chi2 = new TGraphAsymmErrors(_nBins,var_mean_av,chi2MC_vec,hori_av_low,hori_av_high,nullptr,nullptr);
	grMC_chi2->SetLineColor(kOrange-3);
	grMC_chi2->SetMarkerColor(kOrange-3);
	grMC_chi2->SetMarkerStyle(20);
	grMC_chi2->SetLineWidth(2);
	TGraphAsymmErrors* grMC_PSI2S_chi2 = nullptr;
	if (TREE == "ntmix") {
		grMC_PSI2S_chi2 = new TGraphAsymmErrors(_nBins,var_mean_av,chi2MC_PSI2S_vec,hori_av_low,hori_av_high,nullptr,nullptr);
		grMC_PSI2S_chi2->SetLineColor(kOrange-2);
		grMC_PSI2S_chi2->SetMarkerColor(kOrange-2);
		grMC_PSI2S_chi2->SetMarkerStyle(21);
		grMC_PSI2S_chi2->SetLineWidth(2);
	}

	if(VAR == "By"){
		mg_chi2->GetXaxis()->SetTitle("Rapidity (y)");
		mg_chi2->GetXaxis()->SetLimits(0 ,2.4);
	} else if(VAR == "Bpt"){
		mg_chi2->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		mg_chi2->GetXaxis()->SetLimits(0 ,80); 
	} else if(VAR == "nSelectedChargedTracks"){
		mg_chi2->GetXaxis()->SetTitle("Multiplicity (Mult)");
		mg_chi2->GetXaxis()->SetLimits(0, 110);
	}
	mg_chi2->GetYaxis()->SetTitle("#chi^{2}/NDF");
	mg_chi2->GetYaxis()->SetRangeUser(0.0, chi2_max*1.5);
	mg_chi2->Add(gr_chi2);
	mg_chi2->Add(grMC_chi2);
	if (TREE == "ntmix") mg_chi2->Add(grMC_PSI2S_chi2);
	mg_chi2->Draw("ap");
	TString particleLabel = "Particle";
	TLegend *leg_chi2 = new TLegend(0.68,0.78,0.92,0.90,NULL,"brNDC");
	if (TREE == "ntmix"){
		particleLabel = "X(3872)";
		leg_chi2 = new TLegend(0.68,0.72,0.92,0.90,NULL,"brNDC");
	}
	else if (TREE == "ntphi") particleLabel = "B_{s}^{0}";
	else if (TREE == "ntKp") particleLabel = "B^{+}";
	else if (TREE == "ntKstar") particleLabel = "B^{0}";
	leg_chi2->AddEntry(gr_chi2, "Data fit", "lp");
	leg_chi2->AddEntry(grMC_chi2, Form("%s MC fit", particleLabel.Data()), "lp");
	if (TREE == "ntmix") {
		leg_chi2->AddEntry(grMC_PSI2S_chi2, "#psi(2S) MC fit", "lp");
	}
	leg_chi2->SetBorderSize(0);
	leg_chi2->SetFillStyle(0);
	leg_chi2->SetTextFont(42);
	leg_chi2->SetTextSize(0.035);
	leg_chi2->Draw();

	const char* pathc_chi2 =Form("./results/Graphs/chi2_%s_%s.pdf",TREE.Data(),VAR.Data()); 
	c_chi2.SaveAs(pathc_chi2);
	// Nominal part ONLY

	if(syst_study==1){
		//chi2 plot part (sigsum) starts
		TCanvas c_chi2_sigsum;
		TMultiGraph* mg_chi2_sigsum = new TMultiGraph();
		TLegend *leg_chi2_sigsum = new TLegend(0.68,0.78,0.92,0.90,NULL,"brNDC"); 
		if(TREE == "ntmix") leg_chi2_sigsum = new TLegend(0.68,0.72,0.92,0.90,NULL,"brNDC");
		double chi2_max_sigsum = 0;

		for(int j=0; j<static_cast<int>(signal.size()); j++){
			for(int i = 0; i < _nBins; i++){if(chi2_vec_sig[j][i] > chi2_max_sigsum){chi2_max_sigsum = chi2_vec_sig[j][i];}}
			TGraphAsymmErrors* gr_chi2_sigsum = new TGraphAsymmErrors(_nBins,var_mean_av,chi2_vec_sig[j].data(),hori_av_low,hori_av_high,nullptr,nullptr);
			gr_chi2_sigsum->SetLineColor(j+2);
			mg_chi2_sigsum->Add(gr_chi2_sigsum);
			leg_chi2_sigsum->AddEntry(gr_chi2_sigsum, signal[j].label.c_str(), "e");
		}
		if(VAR == "By"){
			mg_chi2_sigsum->GetXaxis()->SetTitle("Rapidity (y)");
			mg_chi2_sigsum->GetXaxis()->SetLimits(0 ,2.4);
		} else if(VAR == "Bpt"){
			mg_chi2_sigsum->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			mg_chi2_sigsum->GetXaxis()->SetLimits(0 ,80);
		} else if(VAR == "nSelectedChargedTracks"){
			mg_chi2_sigsum->GetXaxis()->SetTitle("Multiplicity (Mult)");
			mg_chi2_sigsum->GetXaxis()->SetLimits(0, 110);
		}
		mg_chi2_sigsum->GetYaxis()->SetTitle("#chi^{2}/NDF");
		mg_chi2_sigsum->Add(gr_chi2);
		mg_chi2_sigsum->GetYaxis()->SetRangeUser(0.0, chi2_max_sigsum*1.4);
		mg_chi2_sigsum->Draw("ap");
		leg_chi2_sigsum->AddEntry(gr_chi2, "Nominal", "e");
		leg_chi2_sigsum->SetBorderSize(0);
		leg_chi2_sigsum->SetFillStyle(0);
		leg_chi2_sigsum->SetTextFont(42);
		leg_chi2_sigsum->SetTextSize(0.035);
		leg_chi2_sigsum->Draw();

		const char* pathc_chi2_sigsum =Form("./results/Graphs/chi2_%s_%s_signal_summary.pdf",TREE.Data(),VAR.Data()); 
		c_chi2_sigsum.SaveAs(pathc_chi2_sigsum);
		//chi2 plot part (sigsum) ends

		//chi2 plot part (backsum) starts
		TCanvas c_chi2_backsum;
		TMultiGraph* mg_chi2_backsum = new TMultiGraph();
		TLegend *leg_chi2_backsum = new TLegend(0.68,0.78,0.92,0.90,NULL,"brNDC"); 
		if(TREE == "ntmix") leg_chi2_backsum = new TLegend(0.68,0.72,0.92,0.90,NULL,"brNDC");
		double chi2_max_backsum = 0;

		for(int j=0; j<static_cast<int>(background.size()); j++){
			for(int i = 0; i < _nBins; i++){if(chi2_vec_back[j][i] > chi2_max_backsum){chi2_max_backsum = chi2_vec_back[j][i];}}
			TGraphAsymmErrors* gr_chi2_backsum = new TGraphAsymmErrors(_nBins,var_mean_av,chi2_vec_back[j].data(),hori_av_low,hori_av_high,nullptr,nullptr);
			gr_chi2_backsum->SetLineColor(j+2);
			mg_chi2_backsum->Add(gr_chi2_backsum);
			leg_chi2_backsum->AddEntry(gr_chi2_backsum, background[j].label.c_str(), "e");
		}

		if(VAR == "By"){
			mg_chi2_backsum->GetXaxis()->SetTitle("Rapidity (y)");
			mg_chi2_backsum->GetXaxis()->SetLimits(0,2.4);
		} else if(VAR == "Bpt"){
			mg_chi2_backsum->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
			mg_chi2_backsum->GetXaxis()->SetLimits(0 ,80);
		} else if(VAR == "nSelectedChargedTracks"){
			mg_chi2_backsum->GetXaxis()->SetTitle("Multiplicity (Mult)");
			mg_chi2_backsum->GetXaxis()->SetLimits(0, 110);
		}
		mg_chi2_backsum->GetYaxis()->SetTitle("#chi^{2}/NDF");
		mg_chi2_backsum->Add(gr_chi2);
		mg_chi2_backsum->GetYaxis()->SetRangeUser(0.0, chi2_max_backsum*1.4);
		mg_chi2_backsum->Draw("ap");
		leg_chi2_backsum->AddEntry(gr_chi2, "Nominal", "e");
		leg_chi2_backsum->SetBorderSize(0);
		leg_chi2_backsum->SetFillStyle(0);
		leg_chi2_backsum->SetTextFont(42);
		leg_chi2_backsum->SetTextSize(0.035);
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

	TString fullCut = Form("(%s) && (Bmass > %f) && (Bmass < %f)", DOselCUTS.Data(), minhisto, maxhisto);

	TDirectory* savedDir = gDirectory;
	gROOT->cd();                 
	TTree* t1 = tIn->CopyTree(fullCut);
	savedDir->cd();              

	if(true){ // Create a canvas to draw the mass histogram
        TCanvas *canvas = new TCanvas("canvas", "Bmass Distribution", 600, 600);
        canvas->SetLeftMargin(0.15); // or try 0.18 for more space

        // Define histogram parameters
        double hist_Xlow  = minhisto;  // Minimum Bmass
        double hist_Xhigh = maxhisto;   // Maximum Bmass
        double bin_length_MEV = (hist_Xhigh - hist_Xlow)*1000 / nbinsmasshisto;

        // Create a histogram for Bmass
		TH1F *hist_Bmass = new TH1F("hist_Bmass", Form("; m_{J/#psi #pi^{-} #pi^{+}} ; Entries / %.1f MeV", bin_length_MEV), nbinsmasshisto, hist_Xlow, hist_Xhigh);
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
	if(treeName=="ntmix" && sample=="mc"){arg_list.add(*(w.var("isX3872")));}
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
    } else if (var == "nSelectedChargedTracks"){nBins = N_mult_Bins_X;
	} else if (var == "Cent" || var == "CentBin") {nBins = N_cent_Bins_X;}

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
	else if (var == "nSelectedChargedTracks"){for(int c = 0; c <= nBins; ++c){varBINS[c] = nmbinsvec[c];}} 
	else if (var == "Cent" || var == "CentBin") {for(int c = 0; c <= nBins; ++c){varBINS[c] = centbinsvec[c];}}

    return { nBins, varBINS };
}
