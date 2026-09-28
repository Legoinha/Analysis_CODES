#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TChain.h>
#include <TSystem.h>
#include <TStyle.h>
#include <iostream>

#include "aux/parameters.h"
#include "aux/masses.h"

void plot_Data(TString TREE ="ntmix", TString systemNAME = "PbPb23"){
    gStyle->SetOptStat(0);

    // Create a TChain and add all files from the directory
    TChain chain(Form("%s", TREE.Data()));

    if (systemNAME.Contains("PbPb23")) { // PbPb23 data
        if (TREE == "ntmix") {
            chain.Add("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root");
            //chain.Add("/eos/user/h/hmarques/Analysis_CODES/selectionER/ML_tmva/scored_samples/flat_ntmix_PbPb23_scored_DATA.root");
        } else {
            chain.Add(Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/PbPb23/flat_%s_PbPb23_DATA.root", TREE.Data()));
        }
    }  else if (systemNAME.Contains("PbPb18")) {
        if (TREE == "ntmix") {
            chain.Add("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_DATA.root");
        } 
    } else { // ppRef data
        if (TREE == "ntmix"){     chain.Add("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root");}
        else if (TREE == "ntKp") {  chain.Add("./../../RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKp_ppRef_DATA.root");}
        else if (TREE == "ntphi"){  chain.Add("./../../RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntphi_ppRef_DATA.root");}
        else if (TREE == "ntKstar"){chain.Add("./../../RUN3_Data_MC_sharing/Bmesons/ppRef/flat_ntKstar_ppRef_DATA.root");}
    }

    // Create a canvas to draw the histogram
    TCanvas *canvas = new TCanvas("canvas", "Bmass Distribution", 600, 600);
    canvas->SetLeftMargin(0.15); // or try 0.18 for more space

    // Define histogram parameters
    double hist_Xlow = 5;     // Minimum Bmass
    double hist_Xhigh = 6;  // Maximum Bmass
    if (TREE == "ntmix"){hist_Xlow = 3.6; hist_Xhigh = 4.0;}
    int nbinsmasshisto = 40;    
    double bin_length_MEV = (hist_Xhigh - hist_Xlow)*1000 / nbinsmasshisto;

    TString SELECTIONcuts = "1";

    if (TREE == "ntmix") {
        if (systemNAME.Contains("PbPb23")) {
            SELECTIONcuts = "(Bpt > 15 && Bpt < 50) && (abs(By) < 2.4) && (BQvalue < 0.15) && Bchi2Prob > 0.1";
        } else if (systemNAME.Contains("PbPb18")) {
            SELECTIONcuts = "(Bpt > 15 && Bpt < 50) && (abs(By) < 2.4) && (BQvalue < 0.15) && (Bchi2Prob > 0.1)";  
        } else {
            SELECTIONcuts = "(Bpt > 7.5 && Bpt < 50) && BQvalue < 0.15 ";
            //SELECTIONcuts = "(Bpt > 7.5 && Bpt < 50) && Prediction > 0.58 && BQvalue < 0.15 ";
        }
    }

    else if (TREE != "ntmix"){ SELECTIONcuts = "Bnorm_svpvDistance_2D > 4 && Bpt > 7.5"; } // Bnorm_svpvDistance_2D > 4 && Bpt > 7.5 && BtrkPtimb < 0.2 Bnorm_svpvDistance_2D > 4 &&

    //KSTAR_MASS 0.89594

    //SELECTIONcuts = "1";
    cout << "Applying selection cuts: " << SELECTIONcuts.Data() << std::endl;
    //std::cout << "DATA entries (after cuts): " << chain.GetEntries(SELECTIONcuts.Data()) << std::endl;

    TString Xlabel;
    if (TREE == "ntmix")       {Xlabel = "m_{J/#Psi #pi^{+} #pi^{-}} (GeV)";} 
    else if (TREE == "ntphi")  {Xlabel = "m_{J/#Psi K^{+} K^{-}} (GeV)";    }
    else if (TREE == "ntKp")   {Xlabel = "m_{J/#Psi K^{+}} (GeV)";          }
    else if (TREE == "ntKstar"){Xlabel = "m_{J/#Psi K^{+} #pi^{-}} (GeV)";  }

    // Create an histogram for Bmass
    TH1F *hist_Bmass = new TH1F("hist_Bmass", Form("; %s ; Entries / %.1f MeV", Xlabel.Data(), bin_length_MEV), nbinsmasshisto, hist_Xlow, hist_Xhigh);

    chain.Draw("Bmass >> hist_Bmass", Form("%s ",SELECTIONcuts.Data()));
    // Customize the histogram
    hist_Bmass->SetLineColor(kBlack);
    hist_Bmass->SetLineWidth(1);
    hist_Bmass->SetFillColor(kBlack);    
    hist_Bmass->SetFillStyle(3017); 
    hist_Bmass->SetMinimum(0.0);
    hist_Bmass->Draw("");

    // Print system and number of entries in top right corner
    int nentries = hist_Bmass->GetEntries();
    TLatex* sys_entries = new TLatex(0.2, 0.88, Form("%s, N = %d", systemNAME.Data(), nentries));
	sys_entries->SetNDC();
	sys_entries->SetTextAlign(13);
	sys_entries->SetTextFont(42);
	sys_entries->SetTextSize(0.035);
	sys_entries->SetLineWidth(2);
	sys_entries->Draw();

    if (true && TREE == "ntmix") {
        const double yMin = 0.0;
        const double yMax = hist_Bmass->GetMaximum() * 1.02;

        TLine *lineX = new TLine(X3872_MASS, yMin, X3872_MASS, yMax);
        lineX->SetLineStyle(2);
        lineX->SetLineColor(kOrange-3);
        lineX->SetLineWidth(2);
        lineX->Draw("same");

        TLine *linePSI = new TLine(PSI2S_MASS, yMin, PSI2S_MASS, yMax);
        linePSI->SetLineStyle(2);
        linePSI->SetLineColor(kOrange-2);
        linePSI->SetLineWidth(2);
        linePSI->Draw("same");
    }

    gPad->Update();

    // Save the canvas as an image
    canvas->SaveAs(Form("DATA_%s_%s_Bmass.pdf", systemNAME.Data(), TREE.Data()));

    // Clean up
    delete hist_Bmass;
    delete canvas;
}
