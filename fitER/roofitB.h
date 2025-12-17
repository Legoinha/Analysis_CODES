#include "TAxis.h"
#include "aux/uti.h"
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
#include <fstream>
#include <string>
#include <iomanip>
#include "RooMCStudy.h"
//#include <RooMinuit.h>
#include "../plotER/aux/masses.h"
#include <vector>
#include <array>

using namespace RooFit;
using namespace std;

void plot_jpsifit(RooWorkspace& w, int nbin_hist, TString pdf, RooAbsPdf* model, RooDataSet* ds, TString plotName,  bool with_sig);
void fix_parameters(RooWorkspace& w, TString pdfName, bool release=false);
void fit_jpsinp (RooWorkspace& w, int nbin_hist, TString pdf, float bin_i, float bin_f, TString Var);
template<typename... Targs>
void plot_mcfit(int counter, TString pdf_info, TCanvas* cMC, RooWorkspace& w, TString ANvar, RooDataSet* dsMC, TString tree, int numBin, float binmin , float binmax);

Int_t _count=0;

RooFitResult *fit(TString variation, TString pdf,TString tree, TCanvas* c, TCanvas* cMC, RooDataSet* ds, RooDataSet* dsMC, RooRealVar* mass, float binmin, float binmax, RooWorkspace& w, TString which_var, int NBIN){
	
	// give better initial values to the masses
	double init_mean = Bs_MASS;
	double init_mean_PSI2S = 0.;
	if(tree == "ntKp")        {init_mean = Bu_MASS;}
    else if(tree == "ntKstar"){init_mean = Bd_MASS;}
    else if(tree == "ntmix")  {init_mean = X3872_MASS;
							   init_mean_PSI2S = PSI2S_MASS;}
	double subt_window = 0.01; // mass window around the Bmeson mass to perform the fit

	// Parmeters for SIGNAL PDFs in MC
	RooRealVar meanMC(Form("meanMC%d_%s",_count,pdf.Data()),"",init_mean,init_mean-subt_window,init_mean+subt_window) ;
	RooRealVar sigma1MC(Form("sigma1MC%d_%s",_count,pdf.Data()),"",0.01, 0.001, 0.1) ;
	RooRealVar sigma2MC(Form("sigma2MC%d_%s",_count, pdf.Data()),"",0.005, 0.001, 0.02) ;
	RooRealVar sigma3MC(Form("sigma3MC%d_%s",_count, pdf.Data()),"",0.01, 0.001, 0.03) ;
	RooRealVar sigma4cbMC(Form("sigma4cbMC%d_%s",_count, pdf.Data()),"",0.005,0.001,0.05) ;
	RooRealVar alphaMC(Form("alphaMC%d_%s",_count,pdf.Data()),"",4.,0,15);
	RooRealVar nMC(Form("nMC_%d_%s", _count, pdf.Data()),"",10,-100,200);
	RooRealVar* scale = new RooRealVar("scale","scale",1,0.2,2);
	RooProduct scaled_sigma1MC(Form("scaled_sigma1MC%d_%s",_count,pdf.Data()),"scaled_sigma1MC", RooArgList(*scale,sigma1MC));
	RooProduct scaled_sigma2MC(Form("scaled_sigma2MC%d_%s",_count,pdf.Data()),"scaled_sigma2MC", RooArgList(*scale,sigma2MC));
	RooProduct scaled_sigma3MC(Form("scaled_sigma3MC%d_%s",_count,pdf.Data()),"scaled_sigma3MC", RooArgList(*scale,sigma3MC));
	RooProduct scaled_sigma4cbMC(Form("scaled_sigma4cbMC%d_%s",_count,pdf.Data()),"scaled_sigma4cbMC", RooArgList(*scale,sigma4cbMC));
	RooGaussian sig1MC(Form("sig1MC%d_%s",_count,pdf.Data()),"",*mass,meanMC,scaled_sigma1MC);  
	RooGaussian sig2MC(Form("sig2MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma2MC);  
	RooGaussian sig3MC(Form("sig3MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma3MC);  
	RooCBShape  CBMC(Form("CBMC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma4cbMC, alphaMC, nMC);
	RooRealVar sig1fracMC(Form("sig1fracMC%d_%s",_count, pdf.Data()),"", 0.5, 0.01, 1);
	RooRealVar sig2fracMC(Form("sig2fracMC%d_%s",_count, pdf.Data()),"", 0.5, 0.01, 1);  
	// Parameters for SIGNAL PDFs in MC
	
	RooAddPdf* sigMC = nullptr;
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed" )) sigMC = new RooAddPdf(Form("sigMC%d_%s",_count,pdf.Data()),"",RooArgList(sig1MC,sig2MC),sig1fracMC);
	if(variation=="signal" && pdf=="3gauss") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, sig2MC, sig3MC), RooArgList(sig1fracMC, sig2fracMC), true);
	if(variation=="signal" && pdf=="gauss_cb") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, CBMC), sig1fracMC);

	// VARS for PSI2S signal MC
	RooRealVar mean_specMC(  Form("mean_specMC%d_%s",_count,pdf.Data()),"",init_mean_PSI2S,init_mean_PSI2S-subt_window,init_mean_PSI2S+subt_window) ;
	RooRealVar sigma1_specMC(Form("sigma1_specMC%d_%s",_count,pdf.Data()),"",0.005,0.0005,0.01);
	RooRealVar sigma2_specMC(Form("sigma2_specMC%d_%s",_count,pdf.Data()),"",0.001,0.0005,0.01);
	RooRealVar sigma3_specMC(Form("sigma3_specMC%d_%s",_count,pdf.Data()),"",0.005,0.0005,0.01);

	RooRealVar* scale_specMC = new RooRealVar("scale_specMC","scale_specMC",1,0.1,2);
	RooProduct scaled_sigma1_specMC(Form("scaled_sigma1_specMC%d_%s",_count,pdf.Data()),"scaled_sigma1_specMC", RooArgList(*scale_specMC,sigma1_specMC));
	RooProduct scaled_sigma2_specMC(Form("scaled_sigma2_specMC%d_%s",_count,pdf.Data()),"scaled_sigma2_specMC", RooArgList(*scale_specMC,sigma2_specMC));
	RooProduct scaled_sigma3_specMC(Form("scaled_sigma3_specMC%d_%s",_count,pdf.Data()),"scaled_sigma3_specMC", RooArgList(*scale_specMC,sigma3_specMC));
	RooGaussian sig1_specMC(Form("sig1_specMC%d_%s",_count,pdf.Data()),"",*mass,mean_specMC,scaled_sigma1_specMC);  
	RooGaussian sig2_specMC(Form("sig2_specMC%d_%s",_count,pdf.Data()),"",*mass,mean_specMC,scaled_sigma2_specMC);  
	RooGaussian sig3_specMC(Form("sig3_specMC%d_%s",_count,pdf.Data()),"",*mass,mean_specMC,scaled_sigma3_specMC);  
	RooRealVar sig1frac_specMC(Form("sig1frac_specMC%d_%s",_count,pdf.Data()),"", 0.5, 0.01, 1); 
	RooRealVar sig2frac_specMC(Form("sig2frac_specMC%d_%s",_count,pdf.Data()),"", 0.5, 0.01, 1);
	RooAddPdf* sig_specMC = new RooAddPdf(Form("sig_specMC%d_%s",_count,pdf.Data()),"",RooArgList(sig1_specMC,sig2_specMC),RooArgList(sig1frac_specMC), true);

	RooRealVar* nsig_specMC=nullptr;
	RooRealVar* nsigMC=nullptr;
    if(tree == "ntmix"){
		float nsigMC_0      = (dsMC->reduce("Bmass > 3.8"))->sumEntries();  //X3872 MC signal
		float nsig_specMC_0 = (dsMC->reduce("Bmass < 3.8"))->sumEntries();  //Psi2s MC signal
		nsig_specMC = new RooRealVar(Form("nsig_specMC%d_%s",_count,pdf.Data()),"",nsig_specMC_0, 0.9*nsig_specMC_0, 1.2 * nsig_specMC_0);
		nsigMC      = new RooRealVar(Form("nsigMC%d_%s",_count, pdf.Data()),"",nsigMC_0, 0.9*nsigMC_0, 1.2 * nsigMC_0);
		w.import(*nsig_specMC);
	}else{
		nsigMC = new RooRealVar(Form("nsigMC%d_%s",_count, pdf.Data()),"",dsMC->sumEntries(), 0.9*dsMC->sumEntries(), 1.2 * dsMC->sumEntries());
	}
	w.import(*nsigMC);

	RooAddPdf* modelMC=nullptr;
	if(tree == "ntmix"){
		modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC,*sig_specMC),RooArgList(*nsigMC,*nsig_specMC));
	}else{
		if((variation=="signal" && (pdf=="gauss_cb"|| pdf=="3gauss"|| pdf=="fixed"))||variation=="background") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(*nsigMC));
		if(variation =="" && pdf=="") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(*nsigMC));
	}

//////////ROOFIT ROOFIT ROOFIT  MC MC MC MC MC

	mass->setRange("signal",init_mean-0.032, init_mean+0.032);  //focus the MC fit to the signal region to prevent statistical flutuations
	if(tree=="ntmix"){mass->setRange("signal_psi2s",init_mean_PSI2S-0.022, init_mean_PSI2S+0.022);}
	scale->setConstant();
	scale_specMC->setConstant();

	RooFitResult * fitResultMC = nullptr;
	if(tree=="ntmix"){fitResultMC = modelMC->fitTo(*dsMC,Save(), Range("signal,signal_psi2s"));}
	else{fitResultMC = modelMC->fitTo(*dsMC, Save(), Range("signal"));}
	w.import(*modelMC);
	RooRealVar * fMCchi2_params = new RooRealVar(Form("ndfMC_%d_%s",_count, pdf.Data()), "", fitResultMC->floatParsFinal().getSize() );
	w.import(*fMCchi2_params);

	cout << "Signal MC FIT Done" << endl;
	cout << "PLOTing IT !!" << endl;

	if ((variation=="signal" || variation=="") && pdf!="fixed"){
		plot_mcfit(_count, pdf.Data(), cMC, w, which_var, dsMC, tree, NBIN, binmin, binmax);
	}

	///////////ROOFIT ROOFIT ROOFIT MC MC MC MC MC.  



    // FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp
    if(tree == "ntKp" && variation=="" && pdf==""){ 

        // DEFINE MODEL to fit the non prompt background
        //inclusive MC signal Model
        RooRealVar* meannp = 0;
        RooRealVar* sigma1np = 0;
        RooProduct* sigma2np;
        RooRealVar* ratio_sigma12np = 0;
        RooRealVar* cofs_b_np = 0;
        meannp = new RooRealVar(Form("meannp%d_%s",_count,""),"meannp",5.281,5.24,5.32);
        sigma1np = new RooRealVar(Form("sigma1np%d_%s",_count,""),"sigma1np",0.02,0.01,0.03);
        ratio_sigma12np = new RooRealVar(Form("ratio_sigma12np%d_%s",_count,""),"ratio_sigma12np", 2.4, 0.1, 10);
        sigma2np = new RooProduct(Form("sigma2np%d_%s",_count,""), "sigma2np", RooArgList(*sigma1np, *ratio_sigma12np));
        cofs_b_np = new RooRealVar(Form("cofs_b_np%d_%s",_count,""), "cofs_b_np", 0.5, 0., 1.);
        RooGaussian* signal1_b_np = new RooGaussian(Form("signal1_b_np%d_%s",_count,""),"signal_gauss1_b_np",*mass,*meannp,*sigma1np);
        RooGaussian* signal2_b_np = new RooGaussian(Form("signal2_b_np%d_%s",_count,""),"signal_gauss2_b_np",*mass,*meannp,*sigma2np); 
        RooAddPdf* signalnp = new RooAddPdf(Form("signalnp%d_%s",_count,""), "signalnp", RooArgList(*signal1_b_np,*signal2_b_np),*cofs_b_np);
        w.import(*signalnp);
        //inclusive MC signal Model

        // MC Part. Reconstructed Background Model
        RooRealVar* m_nonprompt_scale=0;
        RooRealVar* m_nonprompt_shift=0;
        m_nonprompt_scale = new RooRealVar(Form("m_nonprompt_scale%d_%s",_count,""), "m_nonprompt_scale",0.01, 0.005, 1);
        m_nonprompt_shift = new RooRealVar(Form("m_nonprompt_shift%d_%s",_count,""), "m_nonprompt_shift", 5.15, 5.1, 5.2);
        RooGenericPdf* erfc = new RooGenericPdf(Form("erfc%d_%s",_count,""), "erfc", Form("TMath::Erfc((Bmass-m_nonprompt_shift%d_%s)/m_nonprompt_scale%d_%s)",_count,"",_count,""), RooArgList(*mass, *m_nonprompt_scale, *m_nonprompt_shift));
        // MC Part. Reconstructed Background Model

        // MC Combinatorial Background Model
        //RooRealVar* alpha_np; 
        RooRealVar* alpha_np; 
        alpha_np = new RooRealVar(Form("alpha_np%d_%s", _count,""), "alpha_np",-0.6 , -10.0 , -0.1);
        RooExponential COMB_jpsi(Form("COMB_jpsi%d_%s",_count,""), "COMB_jpsi", *mass, *alpha_np);
            
        if ( (which_var == "Bpt") && ((int)binmin == 50 && (int)binmax == 60)){
            cout << "not enough data for bin " << (int)binmin << "_" << (int)binmax << " use last parameters instead" << endl;
            alpha_np->setVal(w.var(Form("alpha_np%d_%s", _count-1,""))->getVal());
            alpha_np->setConstant();
            m_nonprompt_scale->setVal(w.var(Form("m_nonprompt_scale%d_%s",_count-1,""))->getVal());
            m_nonprompt_scale->setConstant();
            m_nonprompt_shift->setVal(w.var(Form("m_nonprompt_shift%d_%s",_count-1,""))->getVal());
            m_nonprompt_shift->setConstant();
        } 
        // MC Combinatorial Background Model

        RooRealVar jpsinp_fraction(Form("jpsinp_fraction%d_%s",_count,""), "fraction", 0.35, 0.05, 1);
        RooAddPdf* m_jpsinp_cont = new RooAddPdf(Form("m_jpsinp_cont%d_%s",_count,""), "model for jpsi nonprompt bg", RooArgList(COMB_jpsi, *erfc), RooArgList(jpsinp_fraction));
        w.import(*m_jpsinp_cont);
        // DEFINE MODEL to fit the non prompt background

        fit_jpsinp(w,  NBIN, pdf, binmin, binmax, which_var.Data());
        }
    // FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp

    // FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT
    // FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT FIT MC WT



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//FIT THE DATA FIT THE DATA FIT THE DATA
	scale->setConstant(false);	
	scale_specMC->setConstant(false);	
	RooAddPdf* model = nullptr;

	c->cd();
	RooPlot* frame = mass->frame();
	TPad *p1 = new TPad("p1","p1",0.,0.215,1.,1);
	p1->SetBorderMode(1); 
	p1->SetFrameBorderMode(0); 
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.10);
	p1->Draw(); 

	TPad *p2 = new TPad("p2","p2",0.,0.,1.,0.24); 
	p2->SetTopMargin(0.);    
	p2->SetBorderMode(0);
	p2->SetBorderSize(2); 
	p2->SetFrameBorderMode(0); 
	p2->SetTicks(1,1); 
	p2->Draw();

	double n_signal_initial = ds->sumEntries(TString::Format("abs(Bmass-%g)<0.005",init_mean)) ;
	cout << "n_signal_initial " << n_signal_initial << endl; 

	double n_signal_initial_spec=0;
	if (tree == "ntmix"){
		n_signal_initial_spec = ds->sumEntries(TString::Format("abs(Bmass-%g)<0.005",init_mean_PSI2S));
		cout << "n_signal_initial_spec " << n_signal_initial_spec << endl; 
	}

	////////// SIGNAL FUNCTIONS
	RooRealVar nsig(Form("nsig%d_%s",_count,pdf.Data()),"",n_signal_initial*0.4,n_signal_initial*0.2,n_signal_initial);
	RooRealVar mean(Form("mean%d_%s",_count,pdf.Data()),"",meanMC.getVal(),meanMC.getVal()-0.03,meanMC.getVal()+0.03) ;
	RooRealVar sigma1(Form("sigma1%d_%s",_count,pdf.Data()),"",sigma1MC.getVal(),0.001,0.1) ;
	RooRealVar sigma2(Form("sigma2%d_%s",_count,pdf.Data()),"",sigma2MC.getVal(),0.001,0.1) ;
	RooRealVar sigma3(Form("sigma3%d_%s",_count,pdf.Data()),"",sigma3MC.getVal(),0.001,0.1) ;
	RooRealVar sigma4cb(Form("sigma4cb%d_%s",_count,pdf.Data()),"",sigma4cbMC.getVal(),0.001,0.1) ;
	RooRealVar alpha(Form("alpha%d_%s",_count,pdf.Data()),"",alphaMC.getVal(),0,50);
	RooRealVar n(Form("n_%d_%s", _count, pdf.Data()),"",nMC.getVal(),0,500);
	RooProduct scaled_sigma1(Form("scaled_sigma1%d_%s",_count,pdf.Data()),"scaled_sigma1", RooArgList(*scale,sigma1));
	RooProduct scaled_sigma2(Form("scaled_sigma2%d_%s",_count,pdf.Data()),"scaled_sigma2", RooArgList(*scale,sigma2));
	RooProduct scaled_sigma3(Form("scaled_sigma3%d_%s",_count,pdf.Data()),"scaled_sigma3", RooArgList(*scale,sigma3));
	RooProduct scaled_sigmacb(Form("scaled_sigmacb%d_%s",_count,pdf.Data()),"scaled_sigmacb", RooArgList(*scale,sigma4cb));
	RooGaussian sig1(Form("sig1%d_%s",_count,pdf.Data()),"",*mass,mean,scaled_sigma1);  
	RooGaussian sig2(Form("sig2%d_%s",_count,pdf.Data()),"",*mass,mean,scaled_sigma2);  
	RooGaussian sig3(Form("sig3%d_%s",_count,pdf.Data()),"",*mass,mean,scaled_sigma3);  
	RooCBShape  CB(Form("CB%d_%s",_count, pdf.Data()),""   ,*mass,mean,scaled_sigmacb, alpha, n);
	RooRealVar c1(Form("c1%d_%s",_count,pdf.Data()),"",1.,0.,5.);
	RooRealVar sig1frac(Form("sig1frac%d_%s",_count,pdf.Data()),"",sig1fracMC.getVal(),0.,1.);
	RooRealVar sig2frac(Form("sig2frac%d_%s",_count,pdf.Data()),"",sig2fracMC.getVal(),0.,1.);

	RooAddPdf* sig = nullptr;
	if(variation=="signal" && pdf=="3gauss")   sig = new RooAddPdf(Form("sig%d_%s",_count,pdf.Data()),"",RooArgList(sig1, sig2, sig3), RooArgList(sig1frac, sig2frac), true);
	if(variation=="signal" && pdf=="gauss_cb") sig = new RooAddPdf(Form("sig%d_%s",_count,pdf.Data()),"",RooArgList(sig1, CB), sig1frac);
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed")) sig = new RooAddPdf(Form("sig%d_%s",_count,pdf.Data()),"",RooArgList(sig1,sig2),sig1frac);

    // special cases 
	// (PSI2S)
    RooAddPdf* SIG_spec   = nullptr;
	RooRealVar* nsig_spec = nullptr;
	RooRealVar* mean_spec      = nullptr;
	RooRealVar* sigma1_spec    = nullptr;
	RooRealVar* sigma2_spec    = nullptr;
	RooRealVar* sigma3_spec    = nullptr;
	RooProduct* scaled_sigma1_spec = nullptr;
	RooProduct* scaled_sigma2_spec = nullptr;
	RooProduct* scaled_sigma3_spec = nullptr;
	RooGaussian* sig1_spec     = nullptr;
	RooGaussian* sig2_spec     = nullptr; 
	RooGaussian* sig3_spec     = nullptr; 
	RooRealVar* sig1frac_spec  = nullptr;
	RooRealVar* sig2frac_spec  = nullptr;
    if(tree == "ntmix" ){
		nsig_spec = new RooRealVar(Form("nsig_spec%d_%s",_count,pdf.Data()),"",n_signal_initial_spec*0.8,n_signal_initial_spec*0.4,n_signal_initial_spec);
        mean_spec = new RooRealVar(Form("mean_spec%d_%s",_count,pdf.Data()),"",mean_specMC.getVal(),mean_specMC.getVal()-0.03,mean_specMC.getVal()+0.03) ;
	    sigma1_spec = new RooRealVar(Form("sigma1_spec%d_%s",_count,pdf.Data()),"",sigma1_specMC.getVal(),0.001,0.1) ;
	    sigma2_spec = new RooRealVar(Form("sigma2_spec%d_%s",_count,pdf.Data()),"",sigma2_specMC.getVal(),0.001,0.1) ;
	    sigma3_spec = new RooRealVar(Form("sigma3_spec%d_%s",_count,pdf.Data()),"",sigma3_specMC.getVal(),0.001,0.1) ;
        scaled_sigma1_spec = new RooProduct(Form("scaled_sigma1_spec%d_%s",_count,pdf.Data()),"scaled_sigma1_spec", RooArgList(*scale_specMC,*sigma1_spec));
        scaled_sigma2_spec = new RooProduct(Form("scaled_sigma2_spec%d_%s",_count,pdf.Data()),"scaled_sigma2_spec", RooArgList(*scale_specMC,*sigma2_spec));
		scaled_sigma3_spec = new RooProduct(Form("scaled_sigma3_spec%d_%s",_count,pdf.Data()),"scaled_sigma3_spec", RooArgList(*scale_specMC,*sigma3_spec));
        sig1_spec = new RooGaussian(Form("sig1_spec%d_%s",_count,pdf.Data()),"",*mass,*mean_spec,*scaled_sigma1_spec);  
	    sig2_spec = new RooGaussian(Form("sig2_spec%d_%s",_count,pdf.Data()),"",*mass,*mean_spec,*scaled_sigma2_spec);  
		sig3_spec = new RooGaussian(Form("sig3_spec%d_%s",_count,pdf.Data()),"",*mass,*mean_spec,*scaled_sigma3_spec);  
        sig1frac_spec = new RooRealVar(Form("sig1frac_spec%d_%s",_count,pdf.Data()),"",sig1frac_specMC.getVal(),0.,1.);
		sig2frac_spec = new RooRealVar(Form("sig2frac_spec%d_%s",_count,pdf.Data()),"",sig2frac_specMC.getVal(),0.,1.);

		//add variations in the future
        SIG_spec = new RooAddPdf(Form("SIG_spec%d_%s",_count,pdf.Data()),"",RooArgList(*sig1_spec,*sig2_spec,*sig3_spec),RooArgList(*sig1frac_spec, *sig2frac_spec), true);
		sigma1_spec->setConstant();
		sigma2_spec->setConstant();
		sigma3_spec->setConstant();
		sig1frac_spec->setConstant();
		sig2frac_spec->setConstant();
    }
	////////// SIGNAL FUNCTIONS

	//////////  BACKGROUND FUNCTIONS
	RooRealVar nbkg(Form("nbkg%d_%s",_count,pdf.Data()),"",ds->sumEntries() * 0.6,ds->sumEntries()*0.9,ds->sumEntries());
	RooRealVar a0(Form("a0%d_%s",_count,pdf.Data()),"",-0.35,-2,2);
	RooRealVar a1(Form("a1%d_%s",_count,pdf.Data()),"",-0.05,-2,2);
	RooRealVar a2(Form("a2%d_%s",_count,pdf.Data()),"",0.01,-2,2);
	RooRealVar a3(Form("a3%d_%s",_count,pdf.Data()),"",0,-2,2);
	RooChebychev bkg_2nd(Form("bkg%d_%s",_count,pdf.Data()), "", *mass, RooArgList(a0,a1));
	RooChebychev bkg_3rd(Form("bkg%d_%s",_count,pdf.Data()), "", *mass,RooArgSet(a0,a1,a2));
	RooChebychev bkg_4th(Form("bkg%d_%s",_count,pdf.Data()), "", *mass,RooArgSet(a0,a1,a2,a3));
	RooRealVar lambda(Form("lambda%d_%s", _count,pdf.Data()), "lambda",-1.5, -5., 1.);
	RooExponential bkg(Form("bkg%d_%s",_count,pdf.Data()),"",*mass,lambda);

	// B+ PEAKING AND PART. RECONSTRUCTED BACKGROUNDS
	RooRealVar nbkg_part_r(Form("nbkg_part_r%d_%s",_count,pdf.Data()),"",100,0,1e5);
	RooProduct *nbkg_peaking = nullptr;
	RooAbsPdf* jpsipi = nullptr ;   
	RooAbsPdf* erfc = nullptr ;   
	if(tree== "ntKp"){ 
		jpsipi = w.pdf("jpsipi");   
		erfc = w.pdf(Form("erfc%d_%s",_count,""));   
		RooRealVar* jpsipi_to_signal_ratio= w.var("jpsipi_to_signal_ratio");
		nbkg_peaking = new RooProduct(Form("nbkg_peaking%d_%s",_count,pdf.Data()), "number of jpsi pi with fixed ratio to n_signal", RooArgList(nsig, *jpsipi_to_signal_ratio));
	}
	// B+ PEAKING AND PART. RECONSTRUCTED BACKGROUNDS
	//////////  BACKGROUND FUNCTIONS

	//////////////// MODEL MODEL MODEL MODEL
	/////////////////X X X X X X X X X X X
	if(tree == "ntmix"){
		if(variation=="" && pdf=="")              model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*SIG_spec,*sig,bkg_4th),RooArgList(*nsig_spec,nsig,nbkg));
		if(variation=="signal" && pdf=="1gauss" ) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*SIG_spec,sig1,bkg_4th),RooArgList(*nsig_spec,nsig,nbkg));
		if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*SIG_spec,*sig,bkg_3rd),RooArgList(*nsig_spec,nsig,nbkg));
		if(variation=="signal" && (pdf=="3gauss"|| pdf=="fixed")) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*SIG_spec,*sig,bkg_4th),RooArgList(*nsig_spec,nsig,nbkg));
	}

	/////////////////Bs Bs Bs Bs Bs Bs Bs Bs
	if(tree == "ntphi"){
		if((variation=="" && pdf=="") || (pdf=="mass_range")) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg),RooArgList(nsig,nbkg));
		if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_2nd),RooArgList(nsig,nbkg));
		if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_3rd),RooArgList(nsig,nbkg));
		if(variation=="signal" && pdf=="1gauss" ) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(sig1,bkg),RooArgList(nsig,nbkg));
		if(variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb" )) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig, bkg),RooArgList(nsig, nbkg));
	}

	/////////////////BP BP BP BP BP BP BP BP
	if(tree == "ntKp"){
		if((variation=="" && pdf=="")) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg,*erfc,*jpsipi),RooArgList(nsig,nbkg,nbkg_part_r,*nbkg_peaking));
		if(pdf=="mass_range"){         model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg,*jpsipi)      ,RooArgList(nsig,nbkg,*nbkg_peaking));}
		if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg_2nd,*sig,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
		if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg_3rd,*sig,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
		if(variation=="signal" && pdf=="1gauss" ) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg,sig1,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
		if((variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb" ))) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig, bkg, *erfc, *jpsipi),RooArgList(nsig, nbkg, nbkg_part_r,*nbkg_peaking));
	}

	/////////////////B0 B0 BO B0 B0 B0 B0 B0

	//////////////// MODEL MODEL MODEL MODEL

	//////////////// SET PARAMETERS FROM MC FITS
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

	TString fitRange = (pdf == "mass_range") ? "m_range" : "all";
	//RooFitResult* fitResult = model->fitTo(*ds, Save(), Minos(), Extended(kTRUE), Range(fitRange));
	RooFitResult* fitResult = model->fitTo(*ds, Save(), Extended(kTRUE), Range(fitRange));
	fitResult->Print("v"); 
	w.import(*model);

	////// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT

  p1->cd();
  ds->plotOn(frame, Name(Form("ds_cut%d", _count)), Binning(NBIN), MarkerSize(0.5), MarkerStyle(8), MarkerColor(1), LineColor(1), LineWidth(1)); 
  model->plotOn(frame, Name(Form("model%d_%s", _count, pdf.Data())), Range(fitRange), Precision(1e-6), DrawOption("L"), LineColor(2), LineWidth(1));
  model->plotOn(frame, Name(Form("sig%d_%s", _count, pdf.Data())),  Components(*sig), Range(fitRange), Precision(1e-6), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(1)); 
  if(tree== "ntKp"){
	//TString option = (pdf == "mass_range")? "L" : "LF";
    //RooCmdArg drawRange = (pdf == "mass_range")? Range(fitRange) : RooCmdArg();
    model->plotOn(frame, RooFit::Name(Form("erfc%d_%s",_count,"")) , Components(*erfc), Range(fitRange),  NormRange(fitRange), LineColor(kGreen+3), LineStyle(9), LineWidth(2), DrawOption("L"));
	model->plotOn(frame, RooFit::Name("B->J/#psi #pi"), Components(*jpsipi), NormRange(fitRange), DrawOption("LF"), FillColor(kMagenta+1), LineStyle(1), LineColor(kMagenta+1), LineWidth(1)); 
}
	if(tree== "ntmix"){
  		model->plotOn(frame, Name(Form("SIG_spec%d_%s", _count, pdf.Data())),  Components(*SIG_spec), Range(fitRange), Precision(1e-6), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-2), LineStyle(7), LineColor(kOrange-2), LineWidth(1)); 
	}
   	model->plotOn(frame, Name(Form("bkg%d_%s",_count,pdf.Data())) ,  Components(bkg), Range(fitRange), Precision(1e-6),  DrawOption("L"), LineStyle(7), LineColor(4), LineWidth(1));

	//model->paramOn(frame,Layout(0.2, 0.45, 0.5), Format("NEU",AutoPrecision(2)));
	//frame->getAttText()->SetTextSize(0.035);
	//frame->getAttFill()->SetFillStyle(0);
	//frame->getAttLine()->SetLineWidth(0);
	frame->SetTitle("");
	frame->SetXTitle("");
	frame->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(mass->getMax()-mass->getMin())/NBIN*1000));
	frame->GetYaxis()->SetTitleOffset(1.7);
	frame->GetYaxis()->SetTitleSize(0.04);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->SetStats(0);
	double plot_min = 5.25; 
	if(tree=="ntKp") plot_min = 5.05; 
	frame->GetXaxis()->SetRangeUser(plot_min,5.5); 
	frame->GetXaxis()->SetNdivisions(-50205);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleFont(42);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetLabelOffset(0.01);
	frame->Draw();

	TLegend *leg = new TLegend(0.63,0.6,0.83,0.9,NULL,"brNDC"); 
	leg->SetBorderSize(0);
	leg->SetTextSize(0.035);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(frame->findObject(Form("ds_cut%d", _count)), Form("Data"),"LEP");
	leg->AddEntry(frame->findObject(Form("model%d_%s",_count,pdf.Data())),"Fit Model","l");
	leg->AddEntry(frame->findObject(Form("bkg%d_%s",_count,pdf.Data()))," Comb. Bkg.","l");
	if(tree== "ntKp"){
		leg->AddEntry(frame->findObject(Form("sig%d_%s",_count,pdf.Data()))," B^{+} #rightarrow J/#psi #K^{+}","f");
		leg->AddEntry(frame->findObject("B->J/#psi #pi")," B^{+} #rightarrow J/#psi #pi^{+}","f");
		leg->AddEntry(frame->findObject(Form("erfc%d_%s",_count,pdf.Data()))," B #rightarrow J/#psi X","l");
	}
	else if(tree== "ntmix"){
		leg->AddEntry(frame->findObject(Form("SIG_spec%d_%s",_count,pdf.Data()))," #psi(2S) #rightarrow J/#psi #pi^{+} #pi^{-}","f");
		leg->AddEntry(frame->findObject(Form("sig%d_%s",_count,pdf.Data()))," X(3872) #rightarrow J/#psi #pi^{+} #pi^{-}","f");
	}
	else if(tree== "ntphi"){
		leg->AddEntry(frame->findObject(Form("sig%d_%s",_count,pdf.Data()))," B_{s}^{0} #rightarrow J/#psi K^{+}K^{-}","f");
	}

	leg->Draw();
	
	p2->cd();
	RooHist* pull_hist = frame->pullHist(Form("ds_cut%d",_count),Form("model%d_%s",_count,pdf.Data()));
	pull_hist->SetMarkerSize(0.5);
	RooPlot* pull_plot = mass->frame();
	(pull_plot->GetXaxis())->SetRangeUser(mass->getMin(), mass->getMax()); //maxhisto
	pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"XP");
	pull_plot->SetTitle("");
	if(tree=="ntKp")pull_plot->SetXTitle("m_{J/#psiK^{+}} [GeV/c^{2}]");
	if(tree=="ntphi")pull_plot->SetXTitle("m_{J/#psiK^{+}K^{-}} [GeV/c^{2}]");
	if(tree=="ntmix")pull_plot->SetXTitle("m_{J/#psiK^{+}#pi^{-}} [GeV/c^{2}]");
	pull_plot->SetYTitle("Pull");
	pull_plot->GetYaxis()->SetTitleFont(42);  
	pull_plot->GetYaxis()->SetTitleSize(0.15);
	pull_plot->GetYaxis()->CenterTitle();
	pull_plot->GetYaxis()->SetLabelOffset(0.01);
	pull_plot->GetYaxis()->SetLabelSize(0.14);
	pull_plot->GetYaxis()->SetNdivisions(305);
	pull_plot->GetYaxis()->SetTitleOffset(0.45);
	pull_plot->GetXaxis()->SetTitleSize(0.15);
	pull_plot->GetXaxis()->SetTitleOffset(1.2);
	pull_plot->GetXaxis()->CenterTitle();
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelOffset(0.01);
	pull_plot->GetXaxis()->SetLabelSize(0.14);
	pull_plot->GetXaxis()->SetTickLength(0.16);
	pull_plot->GetXaxis()->SetNdivisions(-50205);
	pull_plot->Draw();

	// draw horizontal zero line on top of the pull plot
	TLine *line_ref = new TLine(mass->getMin(), 0., mass->getMax(), 0.);
	line_ref->SetLineStyle(1);
	line_ref->SetLineColor(2);
	line_ref->SetLineWidth(1);
	line_ref->Draw("same");

	cout << "-----------------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "Signal Yield = " << nsig.getVal() << "     yield Error = " << nsig.getError() << "     yield Error Up = " << nsig.getAsymErrorHi() << "     yieldPrintErrDown = " << nsig.getAsymErrorLo() << endl;
	cout << "-----------------------------------------------------------------------------------------------------------------------------" << endl;
	if(tree=="ntmix"){
		cout << "PSI(2S) Signal Yield = " << nsig_spec->getVal() << "     yield Error = " << nsig_spec->getError() << "     yield Error Up = " << nsig_spec->getAsymErrorHi() << "     yieldPrintErrDown = " << nsig_spec->getAsymErrorLo() << endl;
		cout << "-----------------------------------------------------------------------------------------------------------------------------" << endl;
	}

	p1->cd();
	return fitResult;
} // END OF MAIN FITTING FUNCTION











// FIT TO BINED PAR_R_BKG
void fit_jpsinp(RooWorkspace& w, int nbin_hist, TString pdf, float bin_i, float bin_f, TString Var) {

	RooDataSet* d_s = (RooDataSet*) w.data("jpsinp");
	cout << Var.Data() << " bins study of PART_R_BKG" << endl;  

	if(Var == "Bpt" ){d_s = (RooDataSet*) d_s->reduce(Form("(Bpt > %d && Bpt < %d)", (int) bin_i , (int) bin_f) );}
	else if(Var == "By") { d_s = (RooDataSet*) d_s->reduce(Form("abs(By)>%f && abs(By)< %f", bin_i , bin_f) );}
	else if(Var == "nMult"){ cout << "for now MULTIPLICITY uses the INCLUSIVE bin";}
  	
	// Get rid of B+ at gen level
	RooDataSet* ds_cont = (RooDataSet*) d_s->reduce("Bgen != 23333 && Bgen != 23335 && Bgen > 5000");
	// Signal MC inclusive
	RooDataSet* ds_sig = (RooDataSet*) d_s->reduce("Bgen == 23333");

	// Create the necessary folders and define paths
	TString path_to_save = Form("./results/BP/PAR_R_bkg/");
  	gSystem->mkdir(path_to_save, true);

	TString jpsi_fit_plot = "" ;
	TString incSIG_fit_plot = "" ;
	TString jpsi_plot_with_sig = "" ;
	if(Var == "Bpt" ){
		jpsi_fit_plot = path_to_save + TString::Format("/np_fit_%i-%i_%s.pdf", (int) bin_i, (int) bin_f, Var.Data() );
		incSIG_fit_plot = path_to_save + TString::Format("/InclusiveMC_Signal_fit_%i-%i_%s.pdf",(int) bin_i,(int) bin_f, Var.Data() );
  		jpsi_plot_with_sig = path_to_save + TString::Format("/np_fit_signal_%i-%i_%s.pdf",(int) bin_i,(int) bin_f, Var.Data() );
	} else if (Var == "By" ){
		jpsi_fit_plot = path_to_save + TString::Format("/np_fit_%0.1f-%0.1f_%s.pdf", bin_i, bin_f, Var.Data() );
		incSIG_fit_plot = path_to_save + TString::Format("/InclusiveMC_Signal_fit_%0.1f-%0.1f_%s.pdf", bin_i, bin_f, Var.Data() );
  		jpsi_plot_with_sig = path_to_save + TString::Format("/np_fit_signal_%0.1f-%0.1f_%s.pdf", bin_i, bin_f, Var.Data() );
	} else if(Var == "nMult" ){
		jpsi_fit_plot = path_to_save + TString::Format("/np_fit_%s.pdf", Var.Data());
		incSIG_fit_plot = path_to_save + TString::Format("/InclusiveMC_Signal_fit_%s.pdf", Var.Data());
  		jpsi_plot_with_sig = path_to_save + TString::Format("/np_fit_signal_%s.pdf", Var.Data());
	}

	//[START] FIX SHAPE (NP background)
	RooRealVar Bmass = *(w.var("Bmass"));
	RooAbsPdf* m_jpsinp_cont = w.pdf(Form("m_jpsinp_cont%d_%s",_count, pdf.Data()));
	// FIT
	auto cont_result = m_jpsinp_cont->fitTo(*ds_cont, Save());
	// FIT
	plot_jpsifit(w, nbin_hist, pdf, m_jpsinp_cont, ds_cont, jpsi_fit_plot, false);
	fix_parameters(w, Form("m_jpsinp_cont%d_%s",_count,pdf.Data()));
	//[END] FIX SHAPE (NP background) 

	//[START] FIT inclusiveMC SIGNAL 
	RooAbsPdf* signalnp = w.pdf(Form("signalnp%d_%s",_count,pdf.Data()));  
  	RooRealVar n_signal_np(Form("n_signal_np%d_%s",_count,pdf.Data()), "n_signal_np", 1000, 0., 150000); 
	RooExtendPdf signal_ext(Form("signal_ext%d_%s",_count,pdf.Data()), "extended signal pdf", *signalnp, n_signal_np);
	Bmass.setRange("bmc", 5.15, 5.4);
	// FIT
	auto signal_result = signal_ext.fitTo(*ds_sig, Range("bmc"), Save(), Extended());
	// FIT
	//plot_mcfit(w, &signal_ext, ds_sig, incSIG_fit_plot,  Range("bmc"), LineColor(kRed),LineStyle(1), LineWidth(2));		
	fix_parameters(w, Form("signalnp%d_%s",_count,pdf.Data()));	
	//[END] FIT inclusiveMC SIGNAL 

	//FIX RATIO
		// Fix the ratio of jpsipi to signal
		RooRealVar jpsipi_to_signal_ratio("jpsipi_to_signal_ratio", "jpsipi_to_signal_ratio",0.05, 0, 1);
		jpsipi_to_signal_ratio.setVal(0.0384);   		// from PDG 	
  		jpsipi_to_signal_ratio.setConstant();
		w.import(jpsipi_to_signal_ratio);
	//FIX RATIO

	// TEST THE FIT TEplotNameST THE FIT
	// Import from WorkSpace and define variables
	RooAbsPdf* erfc = w.pdf(Form("erfc%d_%s",_count,pdf.Data()));
	RooAbsPdf* COMB_jpsi = w.pdf(Form("COMB_jpsi%d_%s",_count,pdf.Data()));  
	RooAbsPdf* jpsipi = w.pdf("jpsipi");    
  	RooRealVar* sigma1_np = w.var(Form("sigma1np%d_%s",_count,pdf.Data()));      
	RooRealVar n_cont(Form("n_cont_np%d_%s",_count,pdf.Data()), "n_cont_np", 1000, 0., (d_s->sumEntries())*2);
	RooRealVar n_erfc(Form("n_nonprompt%d_%s",_count,pdf.Data()), "n_nonprompt", 1000, 0., (d_s->sumEntries())*2);
	RooProduct* n_jpsipi = new RooProduct(Form("n_jpsipi_by_signal%d_%s",_count,pdf.Data()), "number of jpsis pi with fixed ratio to n_signal", RooArgList(n_signal_np, jpsipi_to_signal_ratio));
	RooRealVar* alpha_comb_m =  w.var(Form("alpha_np%d_%s", _count,""));

	// Unfix the signal and back to better describe each pT bin peak
  	sigma1_np->setConstant(false);
	alpha_comb_m->setConstant(false);

	// BUILD TOTAL PDF
	RooAddPdf* model_inclusive = new RooAddPdf(Form("model_inclusive%d_%s",_count,pdf.Data()), "NP with B+", RooArgList(*signalnp, *jpsipi, *erfc, *COMB_jpsi), RooArgList(n_signal_np, *n_jpsipi, n_erfc, n_cont));
    // FITFITFIT JUST TO CHECK 
    model_inclusive->fitTo(*d_s, Save(), Extended(), NumCPU(4));
	// Plot
    plot_jpsifit(w, nbin_hist, pdf, model_inclusive, d_s, jpsi_plot_with_sig, true);
	// TEST THE FIT TEST THE FIT
}
// FIT TO BINED PAR_R_BKG

void plot_jpsifit(RooWorkspace& w, int nbin_hist, TString pdf, RooAbsPdf* model, RooDataSet* ds, TString plotName, bool with_sig) {
  
  TCanvas* can_np= new TCanvas("can_mc","",600,600);
  can_np->cd();
  TPad *p1 = new TPad("p1","p1",0.,0.2,1.,0.99);
  p1->SetBorderMode(1); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.10);
  p1->Draw(); 
  p1->cd();
  
  RooRealVar Bmass = *(w.var("Bmass"));
  Bmass.setRange("bmass", 5.0, 6.0);
  RooPlot* massframe = Bmass.frame(Title(" "));
  massframe->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(Bmass.getMax()-Bmass.getMin())/100*1000));
  massframe->GetXaxis()->SetTitle("m_{J/#psi #pi^{+}} [GeV/c^{2}]");
  massframe->GetXaxis()->CenterTitle();
  massframe->GetYaxis()->SetTitleOffset(1.5);
  massframe->GetXaxis()->SetTitleOffset(1.2);
  massframe->GetYaxis()->SetTitleSize(0.035);
  massframe->GetYaxis()->SetTitleFont(42);
  massframe->GetXaxis()->SetLabelFont(42);
  massframe->GetYaxis()->SetLabelFont(42);
  massframe->GetXaxis()->SetLabelSize(0.035);
  massframe->GetYaxis()->SetLabelSize(0.035);	
  massframe->GetXaxis()->SetRangeUser(5.0,5.6);
  ds->plotOn(massframe, RooFit::Name("NP"),Binning(nbin_hist), MarkerSize(0.5),MarkerStyle(8),LineColor(1),LineWidth(1));
  model->plotOn(massframe, RooFit::Name("NP Fit"), NormRange("bmass"),LineColor(kRed), LineStyle(1), LineWidth(1));
  model->plotOn(massframe, RooFit::Name("par"),Components(Form("erfc%d_%s",_count,pdf.Data())), NormRange("bmass"), LineColor(kGreen+3), LineStyle(9), LineWidth(2), DrawOption("L"));
  model->plotOn(massframe, RooFit::Name("COMB_jpsi"),Components(Form("COMB_jpsi%d_%s",_count,pdf.Data())), NormRange("bmass"),LineColor(kBlue), LineWidth(1), LineStyle(kDashed));
  if (with_sig) {
   model->plotOn(massframe, RooFit::Name("signal"),Components(Form("signalnp%d_%s",_count,pdf.Data())), NormRange("bmass"), LineColor(kOrange-3), LineStyle(1), LineWidth(1), FillStyle(3002),FillColor(kOrange-3), VLines(), DrawOption("LF"));
   model->plotOn(massframe, RooFit::Name("B->J/#psi #pi"),Components("jpsipi"), NormRange("bmass"),DrawOption("LF"), FillColor(kMagenta+1), LineStyle(1), LineColor(kMagenta+1), LineWidth(1)); 
  }

  model->paramOn(massframe,  Layout(0.6, 0.95, 0.65),Format("NEU", AutoPrecision(1)));
  massframe->getAttText()->SetTextSize(0.035);
  massframe->getAttFill()->SetFillStyle(0);
  massframe->getAttLine()->SetLineWidth(0);
  massframe->Draw();

  TLatex txt;
  TLegend *leg = new TLegend(0.75,0.71,0.89,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);     
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(massframe->findObject("NP"), " MC", "ple");
  leg->AddEntry(massframe->findObject("par"), " Part. Rec. Bkg.", "l");
  leg->AddEntry(massframe->findObject("NP Fit"),"Part. + Comb Bkg.","l");
  // compare yields with gen particles
  if (with_sig) {
    leg->AddEntry(massframe->findObject("signal"), "signal", "f");
    leg->AddEntry(massframe->findObject("COMB_jpsi"), " Comb. Bkg.", "l");
    leg->AddEntry(massframe->findObject("B->J/#psi #pi"), "B->J/#psi #pi", "f");
  				}
  leg->Draw();
  can_np->SaveAs(plotName);
				}



template<typename... Targs>

void plot_mcfit(int counter, TString pdf_info, TCanvas* cMC, RooWorkspace& w, TString ANvar, RooDataSet* dsMC, TString tree, int numBin, float binmin, float binmax){
	
	TString pdf_function = "";
	if(pdf_info == ""){pdf_function = "Double Gaussian";}
	else{pdf_function = pdf_info.Data();}

	cMC->Clear();
	cMC->cd();
	TPad *pad1 = new TPad(Form("pMC1_%d", counter),Form("pMC1_%d", counter), 0., 0.2, 1., 1);
	pad1->SetBorderMode(1);
	pad1->SetFrameBorderMode(0);
	pad1->SetBorderSize(2);
	pad1->SetBottomMargin(0.10);
	pad1->Draw();

	auto* Bmass = w.var("Bmass");
	RooPlot* frame = Bmass->frame(RooFit::Title(" "));
	frame->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(Bmass->getMax()-Bmass->getMin())/numBin*1000));    
	frame->GetYaxis()->SetTitleOffset(2.);
	frame->GetXaxis()->SetTitleOffset(1.2);
	frame->GetYaxis()->SetTitleSize(0.035);
	frame->GetXaxis()->SetTitleSize(0.035);
	frame->GetXaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelSize(0.035);
	frame->GetYaxis()->SetLabelSize(0.035);
	frame->GetXaxis()->SetRangeUser(Bmass->getMin(), Bmass->getMax());
	frame->SetStats(0);
	if(tree == "ntmix"){     frame->SetXTitle("m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}]");}
	else if(tree == "ntphi"){frame->SetXTitle("m_{J/#psi K^{+}K^{-}} [GeV/c^{2}]");}
	else if (tree == "ntKstar"){}
	else if (tree == "ntKp" ){frame->SetXTitle("m_{J/#psi K^{+}} [GeV/c^{2}]");}
	frame->GetXaxis()->CenterTitle();

	//get the model MC
	RooAbsPdf* MODEL   = w.pdf(Form("modelMC%d_%s",counter, pdf_info.Data()));
	RooAbsPdf* sigMC     = w.pdf(Form("sigMC%d_%s",counter,pdf_info.Data()));
	RooAbsPdf* sigSpecMC = nullptr;
	if (tree == "ntmix" ){sigSpecMC = w.pdf(Form("sig_specMC%d_%s",counter,pdf_info.Data()));}

	//PLOT MC FIT
	pad1->cd();
	dsMC->plotOn(frame,Name(MODEL->GetName()),Binning(numBin),MarkerSize(0.5),MarkerStyle(8),LineColor(1),LineWidth(1));
	if (tree=="ntmix"){
	MODEL->plotOn(frame, Components(*sigSpecMC), Name(sigSpecMC->GetName()), Range("signal_psi2s"), NormRange("signal_psi2s"), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-2), LineStyle(7), LineColor(kOrange-2), LineWidth(1));
	}
	MODEL->plotOn(frame, Components(*sigMC)    , Name(sigMC->GetName())    , Range("signal")      , NormRange("signal"),       DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(1));

	//MODEL->paramOn(frame,Layout(0.2, 0.5, 0.70), Format("NEU",AutoPrecision(2)));
	//frame->getAttText()->SetTextSize(0.035);
	//frame->getAttFill()->SetFillStyle(0);
	//frame->getAttLine()->SetLineWidth(0);

	frame->Draw();
	cMC->RedrawAxis();

	TLegend *legMC = new TLegend(0.23,0.72,0.4,0.89,NULL,"brNDC"); 
	legMC->SetBorderSize(0);
	legMC->SetTextSize(0.035);     
	legMC->SetTextFont(42);
	legMC->SetFillStyle(0);
	legMC->AddEntry(frame->findObject(MODEL->GetName()), "Signal MC","lp");
	TString ident_par = "";
	if(tree=="ntphi"){    ident_par = "B^{0}_{s}";}
	else if(tree=="ntKp"){ident_par = "B^{+}";}
	else if(tree =="ntKstar"){
		ident_par = "B^{0}";
		legMC->AddEntry(frame->findObject(sigSpecMC->GetName()),"WT PDF","f");
	}
	else if(tree=="ntmix"){
		ident_par = "X(3872)";
		legMC->AddEntry(frame->findObject(sigSpecMC->GetName()),"#psi(2S), Double Gaussian PDF","f");
	}
	legMC->AddEntry(frame->findObject(sigMC->GetName()),Form("%s, %s PDF",ident_par.Data(),pdf_function.Data()),"f");
  	legMC->Draw();

	TLatex* tex_VAR_bins = new TLatex(0.5,0.5,"");
	if     (ANvar == "Bpt")  {tex_VAR_bins = new TLatex(0.235,0.68, Form("%d < p_{T} < %d GeV", (int) binmin, (int) binmax));}
	else if(ANvar == "By")   {tex_VAR_bins = new TLatex(0.235,0.68, Form("%0.1f < |y| < %0.1f GeV", binmin, binmax));}
	else if(ANvar == "nMult"){tex_VAR_bins = new TLatex(0.235,0.68, Form("%d < Mult < %d GeV",(int)binmin, (int)binmax));}
	tex_VAR_bins->SetNDC();
	tex_VAR_bins->SetTextFont(42);
	tex_VAR_bins->SetTextSize(0.035);
	tex_VAR_bins->SetLineWidth(1);
	tex_VAR_bins->Draw();
}

/* Fix or release the parameters for a given PDF */
void fix_parameters(RooWorkspace& w, TString pdfName, bool release) {
  RooAbsPdf* pdf = w.pdf(pdfName);
  RooAbsData* ds = w.data("jpsinp");
  RooArgSet* par_set = pdf->getParameters(*ds);
  auto itr = par_set->createIterator();
  bool toFix = ! release;
  std::string fix_or_float = (toFix)? "fix " : "float ";
  std::cout << fix_or_float << "parameters:";
  for (auto i = 0; i < par_set->getSize(); ++i) {
    RooRealVar* var = (RooRealVar*) itr->Next();
    var->setConstant(true);
    TString name = var->GetName();
    std::cout << name << ", ";
    var->setConstant(toFix);
  }
  std::cout << "\n";
}

void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, 
		std::vector<std::vector<double> > numbers, std::string caption){
	
	std::ofstream file_check;
	std::ofstream file;

	//Begin Document                                                                                    
	file.open(filename + ".tex");
	file_check.open(filename + "_check.tex");

	file_check << "\\documentclass{article}" << std::endl;
	//file << "\\usepackage[utf8]{inputenc}" << std::endl;     
	file_check << "\\usepackage{rotating}" << std::endl;                                                                                   
	// file_check << "\\usepackage{cancel}" << std::endl;
	file_check << "\\usepackage{geometry}" << std::endl;
	file_check << "\\usepackage{booktabs}" << std::endl;
	file_check << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;

	file_check << "\\begin{document}" << std::endl;
	// Create table                                                                                                                                                                                                                                                
	std::string col="c";
	for(int i=1; i<n_col; i++)
		col+="|c";

		file_check << "\\begin{sidewaystable}"<< std::endl;
		file_check << "\\begin{tabular}{"+col+"}" << std::endl;
		file_check << "\\toprule" << std::endl;
		file << "\\begin{tabular}{"+col+"}" << std::endl;
		file << "\\toprule" << std::endl;

	for(int c=0; c<n_col-1; c++){
		file << col_name[c] << " & ";
		file_check << col_name[c] << " & ";
	}

	file << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;
	file_check << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;

	for(int i=1; i<n_lin; i++)
	{
		file << labels[i-1] << " & ";
		file_check << labels[i-1] << " & ";

		for(int c=1; c<n_col-1; c++){
			file << /*std::setprecision(3)<<*/  numbers[c-1][i-1]<< " \\% & ";
			file_check << /*std::setprecision(3)<<*/  numbers[c-1][i-1]<< " \\% & ";
									}
		file << /*std::setprecision(3)<<*/  numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl;
		file_check << /*std::setprecision(3)<<*/  numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl; 
	}

	file << "\\bottomrule" << std::endl;
	file_check << "\\bottomrule" << std::endl;

	//End Table                                                                                                                                    
	file << "\\end{tabular}" << std::endl;
	file_check << "\\end{tabular}" << std::endl;
	file_check << "\\caption{"+caption+"}" << std::endl;
	file_check << "\\end{sidewaystable}"<< std::endl;

	//file_check << "\\end{table}" << std::endl;
	//End document                                                                                                                                 

	file_check << "\\end{document}" << std::endl;

	file.close();
	file_check.close();
	system(("pdflatex " + filename+ "_check.tex").c_str());
	//system(("open " + filename + "_check.pdf").c_str());
}

void validate_fit(RooWorkspace* w, TString pdf, TString tree, TString variable, int full, float binmin, float binmax, string Path)
{
	std::cout << "Performing Check on Fit" << std::endl;
	RooRealVar Bmass = *(w->var("Bmass"));
	RooAbsPdf* model  = w->pdf(Form("model%d_%s",_count,pdf.Data()));
	//RooDataSet* data = (RooDataSet*) w->data("data");

	//model->fitTo(*data);
	vector<RooRealVar> params;
	params.push_back(*(w->var(Form("nsig%d_%s",_count,pdf.Data()))));
	//params.push_back(*(w->var(Form("mean%d_%s",_count,pdf.Data()))));

	double n_signal_init = params[0].getVal();
	double n_signal_error_init = params[0].getError();
	int params_size = params.size();

	RooMCStudy* mcstudy = new RooMCStudy(*model, Bmass,  Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));
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