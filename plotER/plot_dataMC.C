#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLegend.h>
#include <iostream>
#include <TStyle.h>

#include "aux/parameters.h"
#include "aux/masses.h"

void plot_dataMC(TString TREE ="ntmix", TString systemNAME = "PbPb24")
{

    //VARIABLES to be plotted
    //VARIABLES
    const char * variables[] = {"Bmass", "abs(BLxy)",  "BsvpvDistance_2D",  "abs(By)",  "BtktkvProb",    "Bpt",     "BQvalue",   "Bcos_dtheta", "BtrkPtimb", "Bchi2Prob", "Btrk1dR", "Btrk1Pt", "Btrk2Pt", "Bnorm_svpvDistance_2D",   };
    const double ranges[][2] = {{3.6,4},    {0,0.1},             {0,0.25},    {0,2.4},         {0,1},   {0,50},     {0.0,0.2},         {0.9,1},    {0,0.45},     {0.0,1},     {0,1},    {0.,5},    {0.,5},                  {0,20},   };    
    //const char * variables[] = {"xgb_score"};
    //const double ranges[][2] = {{0,1}};
    
    //VARIABLES
    //VARIABLES
    //VARIABLES

    TChain chain(TREE.Data());

    ////// OPEN FILES (MC AND DATA) //////
    TTree *tree_MC = nullptr;
    TString path_to_file = "";
    TString path_to_MC   = "";

    if (systemNAME.Contains("PbPb24")){ //PbPb system
        if (TREE == "ntmix"){
            chain.Add(   "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root");
            path_to_MC = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC.root";
        }
    }
    else if (systemNAME.Contains("PbPb23")){ //PbPb system
        if (TREE == "ntmix"){
            chain.Add(   "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root");
            path_to_MC = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC.root";
        }
    }
    else{ //ppRef system
        if (TREE == "ntmix"){     //X3872
            chain.Add( "/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_DATA_wScore.root");
            path_to_MC="/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_MC_wScore.root";
        }
        else if (TREE == "ntphi"){ //B0s
            chain.Add("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntphi_ppRef_DATA.root"); 
            path_to_MC = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntphi_ppRef_MC.root";
        }
        else if (TREE == "ntKp"){  //B+
            chain.Add("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKp_ppRef_DATA.root"); 
            path_to_MC = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKp_ppRef_MC.root";
        }
        else if (TREE == "ntKstar"){  //B0
            chain.Add("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_DATA.root"); 
            path_to_MC = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_MC.root";
        }
    }

    std::cout << "DATA entries: " << chain.GetEntries()    << std::endl;
    TFile::Open(path_to_MC.Data())->GetObject(TREE.Data(), tree_MC  );
    if (TREE == "ntmix" && false) {
        std::cout << "X3872 MC entries: " << tree_MC->GetEntries("isX3872==1") << std::endl;
        std::cout << "PSI2S MC entries: " << tree_MC->GetEntries("isX3872==0") << std::endl;
    } 
    else {std::cout << " MC entries: " << tree_MC->GetEntries() << std::endl;}

    int nVars = sizeof(variables)/sizeof(variables[0]);
    for (int i = 0; i < nVars; ++i){
        TString var = variables[i];

        // Create a canvas to draw the histograms
        TCanvas *canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetLeftMargin(0.15);
        canvas->SetTopMargin(0.05);
        canvas->SetRightMargin(0.05);
        // Hide stats boxes globally
        gStyle->SetOptStat(0);

        int nbinsVARhistos = 100;
        double hist_Xhigh      = ranges[i][1];
        double hist_Xlow       = ranges[i][0];
        double bin_length_MEV  = (hist_Xhigh - hist_Xlow) / nbinsVARhistos;
        
        TString Xlabel ;
        if (var == "Bmass")   {
            if (TREE == "ntKp")        {Xlabel = "m_{J/#Psi K^{+}} [GeV/c^{2}]";}
            else if (TREE == "ntKstar"){Xlabel = "m_{J/#Psi K^{+} #pi^{-}} [GeV/c^{2}]";}
            else if (TREE == "ntphi")  {Xlabel = "m_{J/#Psi K^{+} K^{-}} [GeV/c^{2}]";}
            else if (TREE == "ntmix")  {Xlabel = "m_{J/#Psi #pi^{+} #pi^{-}} [GeV/c^{2}]";}
        }
        else if (var == "Bpt"){Xlabel = "p_{T} [GeV/c]";}
        else {                 Xlabel = var.Data();}

        // Create histograms
        TH1F *hist_SIG = new TH1F("hist_SIG" , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh); 
        TH1F *hist_BKG = new TH1F("hist_BKG" , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh);
        TH1F *hist_spec = new TH1F("hist_spec", Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh);
               
        TString sideband = "1"; 
        if (true){ // BKG from DATA sideband
            if(TREE == "ntmix"){                            sideband = "(Bmass < 3.63 || Bmass > 3.95 || (Bmass > 3.75 && Bmass < 3.8))";}
            else if (TREE == "ntphi" || TREE == "ntKstar") {sideband = "(Bmass > 5.55)";}
            else if (TREE == "ntKp")                       {sideband = "((Bmass > 5.15 && Bmass < 5.20) || (Bmass > 5.40 && Bmass < 5.55))";}
        }

        TString ANYsel = "1"; // customize if needed

        if (TREE == "ntmix") {
            tree_MC->Draw(Form("%s >> hist_SIG" , var.Data()), Form(" %s && isX3872==1 ", ANYsel.Data()));
            tree_MC->Draw(Form("%s >> hist_spec", var.Data()), Form(" %s && isX3872==0 ", ANYsel.Data()));
        } else {
            tree_MC->Draw(Form("%s >> hist_SIG", var.Data()), Form(" %s ", ANYsel.Data()));
        }
        chain.Draw(Form("%s >> hist_BKG", var.Data()), Form(" %s && %s ", sideband.Data(), ANYsel.Data()));  

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

        // Customize Histograms
        hist_SIG->SetLineColor(kOrange-3);
        hist_SIG->SetLineWidth(3);
        hist_SIG->SetFillStyle(0);
        hist_SIG->Scale(1.0 / hist_SIG->Integral());
        if(TREE == "ntmix"){
            hist_spec->SetLineColor(kOrange-2);
            hist_spec->SetLineWidth(3);
            hist_spec->SetFillStyle(0);
            hist_spec->Scale(1.0 / hist_spec->Integral());
        }
        hist_BKG->Scale(1.0 / hist_BKG->Integral());
        hist_BKG->SetLineColor(kBlue);
        hist_BKG->SetFillColor(kBlue);     
        hist_BKG->SetFillStyle(3358); 
        // Customize the Histograms

        if(1){// set the y-axis maximum if needed
            Double_t max_val = TMath::Max(hist_BKG->GetMaximum(), hist_SIG->GetMaximum()) * 1.1;
            hist_SIG->SetMaximum(max_val * 1.1);    // Increase the max range to give some space
            hist_BKG->SetMaximum(max_val * 1.1);
        }

        // Force all plots to start at zero on the vertical axis
        hist_SIG->SetMinimum(0.0);
        hist_BKG->SetMinimum(0.0);
        if (TREE == "ntmix") hist_spec->SetMinimum(0.0);

        // Draw the histograms (background first, signal on top)
        hist_BKG->Draw("HIST");
        hist_SIG->Draw("HIST SAME");
        if (TREE == "ntmix") hist_spec->Draw("HIST SAME");
        gPad->Update();

        // Add legend instead of stats boxes
        TLegend *leg = new TLegend(0.68, 0.78, 0.88, 0.94, systemNAME.Data(), "brNDC");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);
        if (TREE == "ntmix") {
            leg->AddEntry(hist_SIG, "MC Sig: X(3872)", "l");
            leg->AddEntry(hist_spec, "MC Sig: #Psi(2S)", "l");
        } else {
            leg->AddEntry(hist_SIG, "MC Sig", "l");
        }
        leg->AddEntry(hist_BKG, "Data sideband", "f");
        leg->Draw();

        // Save the canvas as an image
        canvas->SaveAs(Form("./presel_STUDY_vars/%s_%s_%s.pdf", TREE.Data(), systemNAME.Data() , var.Data()));

        // Clean up
        delete hist_SIG;
        delete hist_BKG;
        delete hist_spec;
        delete canvas;
    }
}

int main() {
    plot_dataMC();
    return 0;
}