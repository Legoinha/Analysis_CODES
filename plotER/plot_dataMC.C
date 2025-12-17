#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLegend.h>
#include <iostream>
#include <TStyle.h>

#include "aux/ACCSEL.h"
#include "aux/parameters.h"
#include "aux/masses.h"

void plot_dataMC(TString PAR_NAME = "X3872",  TString TREE ="ntmix", TString systemNAME = "ppRef")
{

//VARIABLES
//VARIABLES
const char * variables[] = {"Balpha", "BQvalueuj", "Btktkmass", "Bcos_dtheta", "BtrkPtimb", "Bchi2cl", "Btrk1dR", "Btrk1Pt", "Bnorm_svpvDistance_2D", "Bnorm_svpvDistance", "Bnorm_trk1Dxy"};
const double ranges[][2] = {{0,3.15},    {0.,1.5},    {0.,1.5},         {0,1},     {0,0.9},  {0.05,1},     {0,1},  {0.5, 5},                   {0,5},                {0,10},       {-11,11}};
//VARIABLES
//VARIABLES

/////////////////////////////////  ///////////////////////////  ////////////////

TString PRESELcutlevel = ""; // "_RAW", "_ACC", "_SEL", """, 

/////////////////////////////////  ////////////////////////////////  ///////////

    TChain chain(Form("Bfinder/%s", TREE.Data()));

    ////// SELECT PARTICLE
    TString path_to_file = "";
    TString path_to_MC   = "";
    if (systemNAME.Contains("pbpb") || systemNAME.Contains("PbPb")){
        if (PAR_NAME == "X3872"){
            chain.Add("/lstore/cms/hlegoinha/DATA_X3872_PbPb.root");
            path_to_MC   = "/lstore/cms/lekai/X/MC/MC_X3872_pbpb.root";
        }
        else if (PAR_NAME == "PSI2S"){
            chain.Add("/lstore/cms/hlegoinha/DATA_X3872_PbPb.root");
            path_to_MC   = "";
        }
    }else{
        if (PAR_NAME == "X3872"){
            chain.Add("/lstore/cms/lekai/X/data/DATA_pp_XPSI.root");
            path_to_MC   = "/lstore/cms/lekai/X/MC/MC_X3872.root";
        } else if (PAR_NAME == "PSI2S"){
            chain.Add("/lstore/cms/lekai/X/data/DATA_pp_XPSI.root");
            path_to_MC   = "/lstore/cms/lekai/X/MC/MC_PSI2S.root";
        }
        else{
            chain.Add("/eos/user/h/hmarques/PPRef2024/Bmesons/PPRefDoubleMuon0/PPRefDoubleMuon0/251212_080026/0000/*");
            chain.Add("/eos/user/h/hmarques/PPRef2024/Bmesons/PPRefDoubleMuon1/PPRefDoubleMuon1/251212_080030/0000/*");
            chain.Add("/eos/user/h/hmarques/PPRef2024/Bmesons/PPRefDoubleMuon2/PPRefDoubleMuon2/251212_000818/0000/*");
            chain.Add("/eos/user/h/hmarques/PPRef2024/Bmesons/PPRefDoubleMuon3/PPRefDoubleMuon3/251212_000822/0000/*");
            if      (TREE == "ntKp" )  {path_to_MC = "/lstore/cms/lekai/Bmeson/MC/ppRef/MC_Bu.root";}
            else if (TREE == "ntKstar"){path_to_MC = "/lstore/cms/lekai/Bmeson/MC/ppRef/MC_Bd.root";}
            else if (TREE == "ntphi")  {path_to_MC = "/lstore/cms/lekai/Bmeson/MC/ppRef/MC_Bs.root";}
        }
    }

    // Get the trees from the file
    TTree *tree_MC   = nullptr;
    TFile::Open(path_to_MC.Data())->GetObject(TREE.Data(), tree_MC  );
    std::cout << "DATA entries: " << chain.GetEntries()    << std::endl;
    std::cout << "  MC entries: " << tree_MC->GetEntries() << std::endl;

    int nVars = sizeof(variables)/sizeof(variables[0]);
    for (int i = 0; i < nVars; ++i){
        TString var = variables[i];

        // Create a canvas to draw the histograms
        TCanvas *canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetLeftMargin(0.15);
        canvas->SetTopMargin(0.05);
        canvas->SetRightMargin(0.05);

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
        TH1F *hist_SIG = new TH1F("hist_SIG", Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh); 
        TH1F *hist_BKG = new TH1F("hist_BKG", Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh);
        TH1F *hist     = new TH1F("hist"    , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh);
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        //SELECT THE acc + presel CUT 
        
        TString SELcuts    = SELcuts_ppRef.Data(); //ppRef
        TString trgmatches = TRGmatching.Data();   //TRG matching only in ppRef
        TString ACCcuts    = ACCcuts_ppRef.Data(); //ppRef
        if (TREE == "ntKp"){
            SELcuts = SELcuts_ppRef_Bu.Data();
            ACCcuts = ACCcuts_ppRef_Bu.Data();}
        if (path_to_file.Contains("pbpb") || path_to_file.Contains("PbPb")){ 
            SELcuts = SELcuts_PbPb.Data();
            ACCcuts = ACCcuts_PbPb.Data();
            if (TREE == "ntKp"){
                SELcuts = SELcuts_PbPb_Bu.Data();
                ACCcuts = ACCcuts_PbPb_Bu.Data();}
            trgmatches = "1";
        }
        TString cut = "";
        if      (PRESELcutlevel == "_RAW"){cut = "1";}                                                                         //RAW (inside fid reg only)
        else if (PRESELcutlevel == "_ACC"){cut = Form(" %s "            , ACCcuts.Data());}                                    //ACC
        else if (PRESELcutlevel == "_SEL"){cut = Form(" %s && %s  "     , ACCcuts.Data(), SELcuts.Data());}                    //SEL
        else if (PRESELcutlevel == ""    ){cut = Form(" %s && %s && %s ", ACCcuts.Data(), SELcuts.Data(), trgmatches.Data());}                                                                                              

        TString sepcCASES = "1";
        
        TString Resonances = "1";
        if      (PAR_NAME == "Bs"){Resonances = Form("abs(Btktkmass-%f)<0.015", PHI_MASS);}
        else if (PAR_NAME == "Bd"){Resonances = Form("abs(Btktkmass-%f)<0.15" , KSTAR_MASS);}
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        TString sideband = "";
        if(PAR_NAME == "X3872"    || PAR_NAME == "PSI2S"){sideband = "(Bmass < 3.60 || Bmass > 3.95 || (Bmass > 3.75 && Bmass < 3.8))";}
        else if (PAR_NAME == "Bs" || PAR_NAME == "Bd")   {sideband = "(Bmass > 5.55 || Bmass < 5.15  )";}
        else if (PAR_NAME == "Bu")                       {sideband = "(Bmass > 5.55 )";}

        tree_MC->Draw(Form("%s >> hist_SIG", var.Data()), Form(" %s && %s && %s", isMCsignal.Data(), cut.Data(), sepcCASES.Data()));
        chain.Draw(Form("%s >> hist_BKG", var.Data()), Form(" %s && %s && %s", sideband.Data()  , cut.Data(), sepcCASES.Data()));  // BKG from DATA sideband

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

        // Customize the Histograms
        hist->SetLineColor(kBlack);
        hist->SetLineWidth(2);
        hist_SIG->SetLineColor(kOrange+7);
        hist_SIG->SetFillColor(kOrange+7);    
        //hist_SIG->SetFillStyle(3001); 
        //hist_SIG->SetLineWidth(2);
        //hist_SIG->SetLineStyle(2);

        //hist_SIG_WT->SetLineColor(kOrange);

        hist_BKG->SetLineColor(kBlue);
        hist_BKG->SetFillColor(kBlue);     
        hist_BKG->SetFillStyle(3358);
        //hist_BKG->SetLineStyle(2);
        //hist_BKG->SetLineWidth(2);

        //hist_BKG->SetMarkerStyle(20);// Circle marker
        //hist_BKG->SetMarkerSize(.8); // Bigger dots
        // Customize the Histograms

        //double nEntries = hist->GetEntries();
        if (hist_SIG->Integral()) {
            hist_SIG->Scale(1.0 / hist_SIG->Integral());
            hist_BKG->Scale(1.0 / hist_BKG->Integral());
            //hist->Scale(1.0 / nEntries);
        }

        if(1){// set the y-axis maximum if needed
            Double_t max_val = TMath::Max(hist_BKG->GetMaximum(), hist_SIG->GetMaximum()) * 1.1;
            hist_SIG->SetMaximum(max_val * 1.1);  // Increase the max range to give some space
            hist_BKG->SetMaximum(max_val * 1.1);
        }

        // Draw the histograms
        hist->SetStats(0);
        hist_SIG->Draw("HIST");
        hist_BKG->Draw("HIST SAMES");
        hist->Draw("HIST SAME");
        gPad->Update();

        // Move and color the stat boxes
        TPaveStats *st_bkg = (TPaveStats*)hist_BKG->GetListOfFunctions()->FindObject("stats");
        if (st_bkg) {
            st_bkg->SetTextColor(kBlue);
            st_bkg->SetLineColor(kBlue); 
            st_bkg->SetX1NDC(0.75);
            st_bkg->SetX2NDC(0.95);
            st_bkg->SetY1NDC(0.85);
            st_bkg->SetY2NDC(0.95);
            st_bkg->Draw();
        }
        TPaveStats *st_sig = (TPaveStats*)hist_SIG->GetListOfFunctions()->FindObject("stats");
        if (st_sig) {
            st_sig->SetTextColor(kOrange+7);
            st_sig->SetLineColor(kOrange+7);
            st_sig->SetX1NDC(0.75);
            st_sig->SetX2NDC(0.95);
            st_sig->SetY1NDC(0.75);
            st_sig->SetY2NDC(0.85);
            st_sig->Draw();
        }
        
        // LATEX text
        if(0){
            double Nsignal = hist_SIG->GetEntries();
            double Nbkg = hist_BKG->GetEntries();
            double significance = (Nbkg > 0) ? Nsignal / sqrt(Nbkg) : 0;
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.022);
            latex.SetTextColor(kOrange+7); // Same as hist_SIG
            latex.DrawLatex(0.18, 0.82, Form("N_{sig} = %.0f", Nsignal));
            latex.SetTextColor(kBlue);     // Same as hist_BKG
            latex.DrawLatex(0.18, 0.85, Form("N_{bkg} = %.0f", Nbkg));
        }

        // Add a legend
        auto legend = new TLegend(0.15, 0.7, 0.25, 0.9);
        legend->AddEntry(hist_SIG, "MC SIG", "l");
        legend->AddEntry(hist_BKG, "MC BKG", "l");
        //legend->Draw();

        // Save the canvas as an image
        if(path_to_file.Contains("pbpb") || path_to_file.Contains("PbPb")){systemNAME = "PbPb";}
        canvas->SaveAs(Form("./presel_STUDY_vars/%s_%s_%s_%s.pdf",PAR_NAME.Data(), systemNAME.Data() , var.Data(),  PRESELcutlevel.Data()));

        // Clean up
        delete hist_SIG;
        delete hist_BKG;
        delete hist;
        delete canvas;
    }
}

int main() {
    plot_dataMC();
    return 0;
}