#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLegend.h>
#include <TChain.h>
#include <TSystem.h>
#include <iostream>

#include "aux/parameters.h"

void plot_Data(TString TREE ="ntphi", TString systemNAME = "ppRef"){
    // Create a TChain and add all files from the directory
    TChain chain(Form("%s", TREE.Data()));

    if(systemNAME.Contains("PbPb")){//PbPb data
        chain.Add("/lstore/cms/hlegoinha/DATA_X3872_PbPb.root");    
    } else {//ppRef data
        if (TREE == "ntmix"){       chain.Add("/eos/user/h/hmarques/Analysis_CODES/Flattening/flat_ntmix_ppRef_DATA.root");}
        else if (TREE == "ntKp") {  chain.Add("/eos/user/h/hmarques/Analysis_CODES/Flattening/flat_ntKp_ppRef_DATA.root");}
        else if (TREE == "ntphi"){  chain.Add("/eos/user/h/hmarques/Analysis_CODES/Flattening/flat_ntphi_ppRef_DATA.root");}
        else if (TREE == "ntKstar"){chain.Add("/eos/user/h/hmarques/Analysis_CODES/Flattening/flat_ntKstar_ppRef_DATA.root");}
    }
    std::cout << "DATA entries: " << chain.GetEntries() << std::endl;

    // Create a canvas to draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Bmass Distribution", 600, 600);
    canvas->SetLeftMargin(0.15); // or try 0.18 for more space

    // Define histogram parameters
    double hist_Xlow = 5;      // Minimum Bmass
    double hist_Xhigh = 5.8;   // Maximum Bmass
    if (TREE == "ntmix"){hist_Xlow = 3.6; hist_Xhigh = 4.0;}
    int nbinsmasshisto = 100;    
    double bin_length_MEV = (hist_Xhigh - hist_Xlow)*1000 / nbinsmasshisto;

    TString SELECTIONcuts = "Bnorm_svpvDistance_2D > 4 "
                            "&& Bpt > 5";
    if (TREE == "ntmix") {
        SELECTIONcuts = "BQvalue <  0.15 "
                        "&& Bpt > 5 "
                        "&& Btrk1dR < .55 "
                        "&& Btrk2dR < .55 ";
    }
    //SELECTIONcuts = "1";
    cout << "Applying selection cuts: " << SELECTIONcuts.Data() << std::endl;

    // Create an histogram for Bmass
    TH1F *hist_Bmass = new TH1F("hist_Bmass", Form("; m_{J/#Psi #pi^{+} #pi^{-}} ; Entries / %.1f MeV", bin_length_MEV), nbinsmasshisto, hist_Xlow, hist_Xhigh);
    chain.Draw("Bmass >> hist_Bmass", Form("%s ",SELECTIONcuts.Data()));
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
    if(systemNAME.Contains("PbPb")){ canvas->SaveAs(Form("DATA_PbPb_%s_Bmass.pdf" , TREE.Data()));}
    else{       canvas->SaveAs(Form("DATA_ppRef_%s_Bmass.pdf", TREE.Data()));}

    // Clean up
    delete hist_Bmass;
    delete canvas;
}

