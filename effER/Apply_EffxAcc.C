#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include <vector>
#include <iostream>

#include "aux/uti.h"
#include "../plotER/aux/masses.h"

void Apply_EffxAcc(
	TString treename = "ntmix",
	TString SYSTEM = "ppRef",
	TString VAR = "Bpt"
) {

    cout << "Applying EffxAcc correction for " << treename.Data() << " in " << SYSTEM.Data() << " for variable " << VAR.Data() << endl;

	TString dataFilePath   = GetDataEffPath(treename, SYSTEM);
	TString accEffFilePath = Form("./output/ROOTs/%s_%s2Dmap_ACCxEFF.root"       , treename.Data(), SYSTEM.Data());
	TString yieldsFilePath = Form("./../fitER/ROOTfiles/yields_%s_%s_%s.root"   , treename.Data(), VAR.Data()    , SYSTEM.Data());
	TString outputFilePath = Form("output/ROOTs/%s_%s_%s_CorrectedYields.root"   , treename.Data(), SYSTEM.Data() , VAR.Data());
	
    gSystem->Exec("mkdir -p output/");
	gSystem->Exec("mkdir -p output/ROOTs/");
	gStyle->SetOptStat(0);

	// Load ACCxEFF map
	TFile* fAccEff = TFile::Open(accEffFilePath, "READ");
	TH2D* hACCxEFF = (TH2D*)fAccEff->Get("hACCxEFF");

	// Load data tree
	TFile* fData = TFile::Open(dataFilePath, "READ");
	TString dataTreeName = treename;
	if (treename == "ntmix_psi2s") dataTreeName = "ntmix";
	TTree* tReco = (TTree*)fData->Get(dataTreeName);
	float Bpt = 0.f, By = 0.f, Bmass = 0.f;
    int nMult = 0;
	tReco->SetBranchAddress("Bpt", &Bpt);
	tReco->SetBranchAddress("By" , &By) ;
	tReco->SetBranchAddress("Bmass", &Bmass);
	tReco->SetBranchAddress("nSelectedChargedTracks", &nMult);
    
	// Load yields for binning and corrected yield 
	TFile* fYield = TFile::Open(yieldsFilePath,"READ");
	TH1D* hYield = (TH1D*)fYield->Get("hPt");

	// Determine bins
	std::vector<double> Yield_varBin;
	if (hYield) {
		int nBins = hYield->GetNbinsX();
		Yield_varBin.reserve(nBins + 1);
		for (int i = 1; i <= nBins; ++i){Yield_varBin.push_back(hYield->GetXaxis()->GetBinLowEdge(i));}
		Yield_varBin.push_back(hYield->GetXaxis()->GetBinUpEdge(nBins));
	}
    const int nBins = (int)Yield_varBin.size() - 1;

	// Masses for signal region definition
	double mass_main  = 0.0;
	if (treename == "ntmix")       mass_main = X3872_MASS;
	else if (treename == "ntmix_psi2s") mass_main = PSI2S_MASS;
	else if (treename == "ntphi")   mass_main = Bs_MASS;
	else if (treename == "ntKp")    mass_main = Bu_MASS;
	else if (treename == "ntKstar") mass_main = Bd_MASS;

	// Accumulators for 1/(AccxEff) and its propagated uncertainty
	std::vector<double> sum_INV_EffxAcc(nBins, 0.0), sum_INV_EffxAcc_Err2(nBins, 0.0);
	std::vector<int> cntMain(nBins, 0);

	// Loop reco candidates
	const Long64_t nReco = tReco->GetEntries();
	for (Long64_t i = 0; i < nReco; ++i) {
		tReco->GetEntry(i);
		double absBy = std::abs(By);
		double varVal = 0.0;
		if (VAR == "Bpt") varVal = Bpt;
		else if (VAR == "By") varVal = absBy;
		else if (VAR == "nMult") varVal = nMult;

		int bin = -1;
		for (size_t b = 0; b + 1 < Yield_varBin.size(); ++b) {
			if (varVal >= Yield_varBin[b] && varVal < Yield_varBin[b + 1]) {
				bin = (int)b;
				break;
			}
		}
		if (bin < 0 && varVal == Yield_varBin.back()) bin = (int)Yield_varBin.size() - 2;  // extreme case where value matches upper edge of last bin
		if (bin < 0){
			//std::cout << "[Apply_EffxAcc] WARNING: Candidate with " << VAR << "=" << varVal << " outside of bins. Skipping." << std::endl;
            continue;
        }

		int hbin = hACCxEFF->FindBin(Bpt, absBy);
		double accEff = hACCxEFF->GetBinContent(hbin);
		double accEffUnc = hACCxEFF->GetBinError(hbin);
		if (accEff <= 0.0){ 
            std::cout << "[Apply_EffxAcc] WARNING: Candidate with Bpt=" << Bpt << " and absBy=" << absBy << " has non-positive AccxEff. Skipping." << std::endl;
            continue;
        }

		if (std::abs(Bmass - mass_main) <= 0.05 ) { // consider only candidates in signal region
			double invEff = 1.0 / accEff;
			double invEffErr = accEffUnc / (accEff * accEff);
			sum_INV_EffxAcc[bin] += invEff;
			sum_INV_EffxAcc_Err2[bin] += invEffErr * invEffErr;
			cntMain[bin] ++;
		}
	}

	// Build correction factor histogram (mean 1/AccxEff derived from mean AccxEff)
	TH1D* hAvg_Inv_EffxAcc = new TH1D("hAvg_Inv_EffxAcc", ";p_{T} [GeV];<#frac{1}{Acc#timesEff}>", nBins, Yield_varBin.data());
	TH1D* hYieldCorr = new TH1D("hYieldCorr", ";p_{T} [GeV];Corrected Yield", nBins, Yield_varBin.data());
	hYieldCorr->SetStats(0);
	hAvg_Inv_EffxAcc->SetStats(0);

	for(int i = 0; i < nBins; ++i){
		double avg_Inv_EffAcc = (cntMain[i] > 0) ? (sum_INV_EffxAcc[i] / cntMain[i]) : 0.0;
		double avg_Inv_ErrEffAcc = (cntMain[i] > 0) ? (std::sqrt(sum_INV_EffxAcc_Err2[i]) / cntMain[i]) : 0.0;
		hAvg_Inv_EffxAcc->SetBinContent(i + 1, avg_Inv_EffAcc);
		hAvg_Inv_EffxAcc->SetBinError(  i + 1, avg_Inv_ErrEffAcc);
		std::cout << "[Apply_EffxAcc] Bin " << i << " [" << Yield_varBin[i] << "," << Yield_varBin[i + 1] << "] "
				  << "Correction Factor=" << avg_Inv_EffAcc << "+-" << avg_Inv_ErrEffAcc << " (n=" << cntMain[i] << ")" << std::endl;
	
	// Correct yields (combine fit and eff uncertainties in quadrature)
	double width = Yield_varBin[i + 1] - Yield_varBin[i];
	double rawYield = (hYield) ? (hYield->GetBinContent(i + 1) * width) : 0.0;
	double rawYieldErr = (hYield) ? (hYield->GetBinError(i + 1) * width) : 0.0;
	double yieldCorr = rawYield * avg_Inv_EffAcc;
	double yieldCorr_Unc = std::sqrt((rawYieldErr * avg_Inv_EffAcc) * (rawYieldErr * avg_Inv_EffAcc) +
								 (rawYield * avg_Inv_ErrEffAcc) * (rawYield * avg_Inv_ErrEffAcc));
	hYieldCorr->SetBinContent(i + 1, yieldCorr);
	hYieldCorr->SetBinError(  i + 1, yieldCorr_Unc);
	std::cout << "[Apply_EffxAcc] Corrected yield bin " << i << " [" << Yield_varBin[i] << "," << Yield_varBin[i + 1] << "] "
			  << "yieldCorr=" << yieldCorr << " +- " << yieldCorr_Unc << std::endl;
    }

	// Save output
	TFile* fout = new TFile(outputFilePath, "RECREATE");
	hAvg_Inv_EffxAcc->Write();
	hYieldCorr->Write();
	if (hYield) hYield->Write("hYieldRaw");
	fout->Close();

	// Save <1/ea> plot
	TCanvas *cCorr = new TCanvas("cCorr", "<1/ea>", 700, 600);
	cCorr->SetLeftMargin(0.15);
	hAvg_Inv_EffxAcc->GetYaxis()->SetTitleOffset(1.6);
	hAvg_Inv_EffxAcc->Draw("E1");
	cCorr->SaveAs(Form("output/%s_%s_%s_AvgInvEffxAcc.pdf", treename.Data(), SYSTEM.Data(), VAR.Data()));
	delete cCorr;

	// Save corrected yields plot
	TCanvas *cYield = new TCanvas("cYield", "Corrected yields", 700, 600);
	cYield->SetLeftMargin(0.15);
	hYieldCorr->GetYaxis()->SetTitleOffset(1.6);
	hYieldCorr->Draw("E1");
	cYield->SaveAs(Form("output/%s_%s_%s_CorrectedYields.pdf", treename.Data(), SYSTEM.Data(), VAR.Data()));
	delete cYield;

	fAccEff->Close();
	fData->Close();
	fYield->Close();
}
