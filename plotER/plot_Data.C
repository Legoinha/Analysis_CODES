#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLegend.h>
#include <TChain.h>
#include <TSystem.h>
#include <iostream>

#include "aux/ACCSEL.h"
#include "aux/parameters.h"

bool isPbPb = false; // Set to true for PbPb, false for pp
TString whichTREE = "ntpi"; 

void plot_Data(){
    // Create a TChain and add all files from the directory
    TChain chain(Form("Bfinder/%s", whichTREE.Data()));

    if (isPbPb){
        //PbPb data
        chain.Add("/lstore/cms/hlegoinha/DATA_X3872_PbPb.root");    
    } else {//ppRef data
        
        // X3872 pp reference data
        //chain.Add("/lstore/cms/lekai/X/data/DATA_pp_XPSI.root");

        // Bmesons pp reference data
        chain.Add("/eos/user/h/hmarques/PPRef2024/Bmesons/PPRefDoubleMuon0/PPRefDoubleMuon0/251212_080026/0000/*");
        chain.Add("/eos/user/h/hmarques/PPRef2024/Bmesons/PPRefDoubleMuon1/PPRefDoubleMuon1/251212_080030/0000/*");
        chain.Add("/eos/user/h/hmarques/PPRef2024/Bmesons/PPRefDoubleMuon2/PPRefDoubleMuon2/251212_000818/0000/*");
        chain.Add("/eos/user/h/hmarques/PPRef2024/Bmesons/PPRefDoubleMuon3/PPRefDoubleMuon3/251212_000822/0000/*");
    }
    std::cout << "DATA entries: " << chain.GetEntries() << std::endl;

    // Create a canvas to draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Bmass Distribution", 600, 600);
    canvas->SetLeftMargin(0.15); // or try 0.18 for more space

    // Define histogram parameters
    double hist_Xlow = 5;       // Minimum Bmass
    double hist_Xhigh = 6;      // Maximum Bmass
    if (whichTREE == "ntmix"){ hist_Xlow = 3.6; hist_Xhigh = 4.0; }
    else if (whichTREE == "ntpi"){ hist_Xlow = 5.5; hist_Xhigh = 6.9; }
    int nbinsmasshisto = 100;    
    double bin_length_MEV = (hist_Xhigh - hist_Xlow)*1000 / nbinsmasshisto;

    /*
        TString EXTRAcuts = "BQvalueuj <  0.12"
                            "&& Btrk1Pt > 1.4"
                            "&& Bpt > 15"
                            "&& Bnorm_trk1Dxy < 7.5"
                            "&& Bnorm_svpvDistance_2D < 2"
                            "&& Btrk1dR < 0.3"
                            "&& Btktkmass < 0.8"
                            "&& Bchi2cl > 0.5"; 
    */

    //TString EXTRAcuts = "BQvalueuj <  0.15";

    TString EXTRAcuts = "Bnorm_svpvDistance_2D > 4";
    if(whichTREE=="ntpi"){EXTRAcuts="Bnorm_svpvDistance_2D > 1"
                                    "&& Bpt > 15"
                                    "&& cos(Balpha)>0.9 "
                                    "&& Bchi2cl > 0.05"
                                    "&& Btrk1Pt < 2.5";}

    // Create a histogram for Bmass
    TH1F *hist_Bmass = new TH1F("hist_Bmass", Form("; m_{J/#Psi #pi^{+} #pi^{-}} ; Entries / %.1f MeV", bin_length_MEV), nbinsmasshisto, hist_Xlow, hist_Xhigh);

    if (isPbPb){
        chain.Draw("Bmass >> hist_Bmass", Form(" %s && %s && %s ",
                                                ACCcuts_PbPb.Data(),
                                                EXTRAcuts.Data(),
                                                SELcuts_PbPb.Data()));
    } else { 
        if (TString(chain.GetName()).Contains("ntKp") || TString(chain.GetName()).Contains("ntpi") ) {
            chain.Draw("Bmass >> hist_Bmass", Form("%s && %s && %s && %s " ,
                                                    ACCcuts_ppRef_Bu.Data(),
                                                    TRGmatching.Data()     ,
                                                    SELcuts_ppRef_Bu.Data(),
                                                    EXTRAcuts.Data()))     ;
        } else {
            chain.Draw("Bmass >> hist_Bmass", Form("%s && %s && %s && %s ",
                                                    ACCcuts_ppRef.Data()  ,
                                                    TRGmatching.Data()    ,
                                                    SELcuts_ppRef.Data()  ,
                                                    EXTRAcuts.Data()))    ;
        }
    }
    
    // Customize the histogram
    hist_Bmass->SetLineColor(kBlack);
    hist_Bmass->SetLineWidth(1);
    hist_Bmass->SetFillColor(kBlack);    
    hist_Bmass->SetFillStyle(3017); 
    hist_Bmass->Draw("");
    gPad->Update();

    TPaveStats *st = (TPaveStats*)hist_Bmass->GetListOfFunctions()->FindObject("stats");
    if (st) {
        st->SetTextColor(kBlack);
        st->SetX1NDC(0.7);
        st->SetX2NDC(0.9);
        st->SetY1NDC(0.8);
        st->SetY2NDC(0.9);
        st->Draw();
    }
    gPad->Update();

    // Save the canvas as an image
    if(isPbPb){ canvas->SaveAs(Form("DATA_PbPb_%s_Bmass.pdf" , whichTREE.Data()));}
    else{       canvas->SaveAs(Form("DATA_ppRef_%s_Bmass.pdf", whichTREE.Data()));}

    // Clean up
    delete hist_Bmass;
    delete canvas;
}

int main() {
    plot_Data();
    return 0;
}