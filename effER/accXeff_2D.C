#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TCut.h"
#include "TString.h"
#include "TObjArray.h"
#include "TSystem.h"
#include <iostream>
#include <vector>

#include "aux/uti.h"

// Make a 2D ACC ratio map: (passes ACC) / (total) in (Bpt, By)
// Inputs:
//  - treename:  tree name (e.g., "ntmix", "ntKp")
//  - SYSTEM:    "ppRef" or one of PbPb23/PbPb24/PbPb25 (PbPb all use same ACC)


void accXeff_2D(TString treename = "ntmix", TString SYSTEM = "ppRef")
{

    // Setting Up 
    TString path_to_MC = GetMCEffPath(treename, SYSTEM);
    TString path_to_Gen = GetGenEffPath(treename, SYSTEM);

    gSystem->Exec("mkdir -p output/");
    gSystem->Exec("mkdir -p output/ROOTs/");


    ///////////////////// BINNING /////////////////////

    // variable pT binning
    std::vector<double> ptBins = { 0,5,7.5,10,15,20, 25, 30, 35, 40, 50};
    const int nPtBins = ptBins.size() - 1;

    // variable y binning (optional)
    std::vector<double> yBins = {0.0, 0.5, 1, 1.4, 1.8, 2, 2.4};
    const int nYBins = yBins.size() - 1;

    ///////////////////// BINNING /////////////////////



    // Open file and tree
    TFile *finReco = TFile::Open(path_to_MC, "READ");
    TFile *finGen  = TFile::Open(path_to_Gen, "READ");
    TTree *tree_reco = (TTree*)finReco->Get(treename);
    TTree *tree_gen  = (TTree*)finGen->Get("ntGen");

    // 2D hists (create AFTER opening TFile to avoid directory issues)
    TH2D *hDen_ACC = new TH2D("hDen_ACC",";p_{T};|y|;denom", nPtBins, ptBins.data(), nYBins, yBins.data());
    TH2D *hNum_ACC = new TH2D("hNum_ACC",";p_{T};|y|;numer", nPtBins, ptBins.data(), nYBins, yBins.data());
    hDen_ACC->Sumw2();
    hNum_ACC->Sumw2();

    std::cout << "2D_ACC_EFF: \nreco input=" << path_to_MC << "\ngen input=" << path_to_Gen << "\ntree=" << treename << "\nSYSTEM=" << SYSTEM << std::endl;


    // ------------------------- ACC ------------------------- //
    std::cout << "\n\n--- ACC calculation --- \n" << std::endl;

    // ACC cuts
    // Muons
    TString MUONacc = "("
                      "(Gmu1pt >= 3.5 && abs(Gmu1eta) < 1.2)"
                      " || (Gmu1pt >= (5.47 - 1.89 * abs(Gmu1eta)) && abs(Gmu1eta) >= 1.2 && abs(Gmu1eta) < 2.1)"
                      " || (Gmu1pt >= 1.5 && abs(Gmu1eta) >= 2.1 && abs(Gmu1eta) < 2.4)"
                      ") && ("
                      "(Gmu2pt >= 3.5 && abs(Gmu2eta) < 1.2)"
                      " || (Gmu2pt >= (5.47 - 1.89 * abs(Gmu2eta)) && abs(Gmu2eta) >= 1.2 && abs(Gmu2eta) < 2.1)"
                      " || (Gmu2pt >= 1.5 && abs(Gmu2eta) >= 2.1 && abs(Gmu2eta) < 2.4)"
                      ")";

    //Tracks (depend on channel and system!)
    TString TRK_y_ACC = "(abs(Gtk1eta) < 2.4";
    if(treename != "ntKp"){TRK_y_ACC += " && abs(Gtk2eta) < 2.4";}
    TRK_y_ACC += ")";
    float pTmin_TRK;
    SYSTEM.Contains("PbPb") ? pTmin_TRK = 0.9: pTmin_TRK = 0.5;
    TString TRK_pT_ACC = Form("Gtk1pt >= %.1f", pTmin_TRK);
    if(treename != "ntKp"){TRK_pT_ACC += Form(" && Gtk2pt >= %.1f", pTmin_TRK);}  
    TString TRK_acc = "(" + TRK_y_ACC + " && " + TRK_pT_ACC + ")";

    //Full ACC cut (muons + tracks)
    TString ACCcut = "(" + MUONacc + " && " + TRK_acc + ")";
    std::cout << "ACC cuts: " << "(" + MUONacc + " && " + TRK_acc + ")" << "\n" <<std::endl;

    TString GENcut = "1";
    std::cout << "No ACC cut: "<<  tree_gen->GetEntries(GENcut) <<" Gen signals"       << std::endl;
    std::cout << "w/ ACC cut: "<<  tree_gen->GetEntries(ACCcut) <<" surviving ACC cut" << std::endl;

    tree_gen->Draw("abs(Gy):Gpt>>hDen_ACC", GENcut, "goff");  // Fill Denominator 
    tree_gen->Draw("abs(Gy):Gpt>>hNum_ACC", ACCcut, "goff");  // Fill Numerator 

    // Compute ratio with binomial errors
    TH2D *hACC = (TH2D*)hNum_ACC->Clone("hACC");
    hACC->SetTitle("Acc; p_{T} [GeV]; y; ACC fraction");
    hACC->Divide(hNum_ACC, hDen_ACC, 1.0, 1.0, "B");

    // Style and draw
//    gStyle->SetOptStat(0)  ;
    TCanvas *cACC = new TCanvas("cACC", "ACC Ratio 2D", 900, 700);
    cACC->SetRightMargin(0.15);
    hACC->SetMinimum(0.0);
    hACC->SetMaximum(1.0);
    hACC->SetContour(50);
    hACC->Draw("COLZ")  ;

    // Save outputs
    TString Acc_out  = "output/" + treename + "_" + SYSTEM + "2Dmap_ACC.pdf";
    TString rootName_ACC = "output/ROOTs/" + treename + "_" + SYSTEM + "2Dmap_ACC.root";
    cACC->SaveAs(Acc_out);

    TFile *fout_ACC = new TFile(rootName_ACC, "RECREATE");
    hDen_ACC->Write();
    hNum_ACC->Write();
    hACC->Write();
    cACC->Write();
    fout_ACC->Close();




    // ------------------------- EFF ------------------------- //

    // Quality cuts are already in RECO level tree
    // just apply the FID region & Selection cuts
    std::cout << "\n\n --- EFF calculation --- \n " << std::endl;

    // Selection cuts 
    TString SELcuts = "Bnorm_svpvDistance_2D > 4";
    if(treename == "ntmix" || treename == "ntmix_psi2s"){SELcuts = "xgb_score > 0.55 && BQvalue<0.10";}

    TString SelectionEFF = "(" + SELcuts + ")" ; // apply selection + fiducial cuts at reco level for EFF calculation
    std::cout << "EFF Selection cuts: " << SelectionEFF << "\n" << std::endl;
    std::cout << "W/ quality + selection cuts: " <<  tree_reco->GetEntries(SelectionEFF) << " surviving SEL cuts" << std::endl;

    // 2D Histograms (DENOMINATOR is ACC 2D histogram ==> assume perfect gen-matching)
    TH2D *hNum_EFF = new TH2D("hNum_EFF", ";p_{T} [GeV];|y|;numer", nPtBins, ptBins.data(), nYBins, yBins.data());
    hNum_EFF->Sumw2();
    TH2D *hDen_EFF = (TH2D*)hNum_ACC->Clone("hDen_EFF"); // Denominator is the ACC histogram (gen-level)
    hDen_EFF->SetTitle("Denominator for EFF; p_{T} [GeV];|y|;denom");

    tree_reco->Draw("abs(By):Bpt>>hNum_EFF", SelectionEFF, "goff");  // Fill Numerator 

    // Compute ratio with binomial errors
    TH2D *hEFF = (TH2D*)hNum_EFF->Clone("hEFF");
    hEFF->SetTitle("Eff; p_{T} [GeV]; y; EFF fraction");
    hEFF->Divide(hNum_EFF, hDen_EFF, 1.0, 1.0, "B");

    // Style and draw
//    gStyle->SetOptStat(0)  ;
    TCanvas *cEFF = new TCanvas("cEFF", "EFF Ratio 2D", 900, 700);
    cEFF->SetRightMargin(0.15);
    hEFF->SetMinimum(0.0);
    hEFF->SetMaximum(1.0);
    hEFF->SetContour(50);
    hEFF->Draw("COLZ")  ;

    // Save outputs
    TString Eff_out  = "output/" + treename + "_" + SYSTEM + "2Dmap_EFF.pdf";
    TString rootName_EFF = "output/ROOTs/" + treename + "_" + SYSTEM + "2Dmap_EFF.root";
    cEFF->SaveAs(Eff_out);

    TFile *fout_EFF = new TFile(rootName_EFF, "RECREATE");
    hDen_EFF->Write();
    hNum_EFF->Write();
    hEFF->Write();
    cEFF->Write();
    fout_EFF->Close();










    // ------------------------- ACC x EFF ------------------------- //
    TH2D *hACCxEFF = (TH2D*)hACC->Clone("hACCxEFF");
    hACCxEFF->SetTitle("Acc #times Eff; p_{T} [GeV]; |y|; ACC#timesEFF");
    hACCxEFF->Multiply(hEFF);

    TCanvas *cACCxEFF = new TCanvas("cACCxEFF", "Acc x Eff 2D", 900, 700);
    cACCxEFF->SetRightMargin(0.15);
    hACCxEFF->SetMinimum(0.0);
    hACCxEFF->SetMaximum(1.0);
    hACCxEFF->SetContour(50);
    hACCxEFF->Draw("COLZ");

    TString AccEff_out  = "output/" + treename + "_" + SYSTEM + "2Dmap_ACCxEFF.pdf";
    TString rootName_AccEff = "output/ROOTs/" + treename + "_" + SYSTEM + "2Dmap_ACCxEFF.root";
    cACCxEFF->SaveAs(AccEff_out);

    TFile *fout_AccEff = new TFile(rootName_AccEff, "RECREATE");
    hACCxEFF->Write();
    cACCxEFF->Write();
    fout_AccEff->Close();

    finReco->Close();
    finGen->Close();

}
