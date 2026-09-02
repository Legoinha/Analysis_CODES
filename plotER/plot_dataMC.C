#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLegend.h>
#include <TLatex.h>
#include <iostream>
#include <TStyle.h>
#include <TSystem.h>
#include "aux/parameters.h"
#include "aux/masses.h"



//// TO RUN

// root -l -b -q plot_dataMC.C'("ntmix","ppRef")'
// root -l -b -q plot_dataMC.C'("ntKp","ppRef")'

//// TO RUN



TString getPlotParticleLabel(TString treeName){
    if (treeName == "X3872") return "X(3872)";
    if (treeName == "Psi2S") return "#Psi(2S)";
    if (treeName == "ntKp") return "B^{+}";
    if (treeName == "ntKstar") return "B^{0}";
    if (treeName == "ntphi") return "B_{s}^{0}";
    return treeName;
}

bool hasVariableForDraw(TTree *tree, TString var)
{
    if (!tree) return false;
    if (tree->GetBranch(var)) return true;
    if (var == "abs(By)") return tree->GetBranch("By") != nullptr;
    return false;
}

void plot_dataMC(TString TREE ="ntmix_X3872", TString systemNAME = "ppRef")
{
    gSystem->Exec("mkdir -p ./presel_STUDY_vars/");

    //VARIABLES
    //VARIABLES
    //VARIABLES
    struct PlotVariable {
        const char* expression;
        double xmin;
        double xmax;
    };
    const PlotVariable variables[] = {
        {"PVx",                       -0.1,  0.1},
        {"PVy",                       -0.1,  0.1},
        {"PVz",                      -30.0, 30.0},
        {"PVnchi2",                    0.0,  1.0},
        {"nChargedTracks",             0.0, 200.0},
        {"nChargedTracks_LOOSE",       0.0, 200.0},
        {"nChargedTracks_TIGHT",       0.0, 200.0},
        {"CentBin",                    0.0, 200.0},
        {"Bmass",                      3.6,  4.0},
        {"Bpt",                        0.0, 50.0},
        {"abs(By)",                    0.0,  2.4},
        {"Bchi2Prob",                  0.0,  1.0},
        {"Btrk1dR",                    0.0,  1.5},
        {"Btrk2dR",                    0.0,  1.5},
        {"BtrkPtimb",                  0.0,  1.0},
        {"Btktkpt",                    0.0, 10.0},
        {"Bujmass",                    2.9,  3.3},
        {"BujvProb",                   0.0,  1.0},
        {"Bnorm_svpvDistance_2D",      0.0, 20.0},
        {"BsvpvDistance_2D",           0.0,  0.25},
        {"BsvpvDisErr_2D",             0.0,  0.05},
        {"BQvalue",                    0.0,  0.6},
        {"Bnorm_trk1Dxy",             -5.0,  5.0},
        {"Bnorm_trk2Dxy",             -5.0,  5.0},
        {"Balpha",                     0.0,  3.2},
        {"Bdtheta",                   -3.2,  3.2},
        {"Bcos_dtheta",               -1.0,  1.0},
        {"Btktkmass",                  0.0,  2.0},
        {"Btrk1Pt",                    0.0, 10.0},
        {"Btrk2Pt",                    0.0, 10.0},
        {"Btrk1Eta",                  -2.4,  2.4},
        {"Btrk2Eta",                  -2.4,  2.4},
        {"Btrk1Phi",                  -3.2,  3.2},
        {"Btrk2Phi",                  -3.2,  3.2},
        {"Btrk1PtErr",                 0.0,  0.1},
        {"Btrk2PtErr",                 0.0,  0.1},
        {"BtktkvProb",                 0.0,  1.0},
        {"BLxy",                      -0.1,  0.1},
        {"BvtxX",                     -0.05, 0.05},
        {"BvtxY",                     -0.05, 0.05},
        {"Bmu1pt",                     0.0, 25.0},
        {"Bmu2pt",                     0.0, 25.0},
        {"Bmu1eta",                   -2.4,  2.4},
        {"Bmu2eta",                   -2.4,  2.4},
        {"Bmu1phi",                   -3.2,  3.2},
        {"Bmu2phi",                   -3.2,  3.2},
        {"Bujpt",                      0.0, 50.0},
        {"Bujeta",                    -2.4,  2.4},
        {"Bujphi",                    -3.2,  3.2},
        {"Bujlxy",                    -0.1,  0.1},
        {"BdiTrackFitValid",          -0.5,  1.5},
        {"Btrk1Dz1",                  -0.5,  0.5},
        {"Btrk2Dz1",                  -0.5,  0.5},
        {"Btrk1DzError1",              0.0,  0.1},
        {"Btrk2DzError1",              0.0,  0.1},
        {"Btrk1Dxy1",                 -0.2,  0.2},
        {"Btrk2Dxy1",                 -0.2,  0.2},
        {"Btrk1DxyError1",             0.0,  0.1},
        {"Btrk2DxyError1",             0.0,  0.1},
        {"Btktketa",                  -2.4,  2.4},
        {"Btktkphi",                  -3.2,  3.2},
        {"Btktky",                    -2.4,  2.4},
        {"Bdoubletpt",                 0.0, 50.0},
        {"Bdoubleteta",               -2.4,  2.4},
        {"Bdoubletphi",               -3.2,  3.2},
        {"Bdoublety",                 -2.4,  2.4},
        {"Bnorm_trk1Dz",              -5.0,  5.0},
        {"Bnorm_trk2Dz",              -5.0,  5.0}
    };
    //const char * variables[] = {"BQvalue"};
    //const double ranges[][2] = {{0,0.6}};
    //const char * variables[] = {"Btrk1dR","Btrk2dR",};
    //const double ranges[][2] = {{0,1.5}, {0,1.5}};
    
    //VARIABLES
    //VARIABLES
    //VARIABLES

    TString dataTreeName = TREE;
    TString mcTreeName = TREE;
    TString mcTreeNameSpec = TREE;

    ////// OPEN FILES (MC AND DATA) //////
    ////// OPEN FILES (MC AND DATA) //////
    TTree *tree_MC = nullptr;
    TTree *tree_MC_spec = nullptr;
    TString path_to_data = "";
    TString path_to_MC   = "";
    TString path_to_MC_spec = "";
    TString baseDir = "";
    if (TREE == "ntmix_X3872") {
        dataTreeName = "ntmix";
        if (systemNAME.Contains("PbPb23")) {
            path_to_MC      = Form("/eos/user/k/kprince/X3872_PbPb/MC_X3872_PbPb_AANN.root");
            path_to_MC_spec = Form("/eos/user/k/kprince/X3872_PbPb/MC_PSI2S_PbPb_AANN.root");
            path_to_data    = Form("/eos/user/k/kprince/X3872_PbPb/DATA_PbPb_AANN.root");
            //path_to_MC = Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872.root");
            //path_to_MC_spec = Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S.root");
            //path_to_data = Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root");
        } 
        else if (systemNAME.Contains("PbPb24")) {
            path_to_MC      = Form("/eos/user/k/kprince/X3872_PbPb/MC_X3872_24b_PbPb_AANN.root");
            path_to_MC_spec = Form("/eos/user/k/kprince/X3872_PbPb/MC_PSI2S_24b_PbPb_AANN.root");
            path_to_data    = Form("/eos/user/k/kprince/X3872_PbPb/DATA_24b_PbPb_AANN.root");
        } else {
            //path_to_data    = Form("/eos/user/k/kprince/X3872_pp_new/DATA_pp_VAANN.root");
            //path_to_MC_spec = Form("/eos/user/k/kprince/X3872_pp_new/MC_PSI2S_pp_VAANN.root");
            //path_to_MC      = Form("/eos/user/k/kprince/X3872_pp_new/MC_X3872_pp_VAANN.root");
            path_to_MC = Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root");
            path_to_MC_spec = Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root");
            path_to_data = Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root");
        }
        mcTreeNameSpec = "ntmix_PSI2S";
    } else {
        path_to_MC   = Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/%s/flat_%s_%s_MC.root",systemNAME.Data(), TREE.Data(), systemNAME.Data());
        path_to_data = Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/%s/flat_%s_%s_DATA.root",systemNAME.Data(), dataTreeName.Data(), systemNAME.Data());
    }
    TChain chain(dataTreeName.Data());
    chain.Add(path_to_data);
    ////// OPEN FILES (MC AND DATA) //////
    ////// OPEN FILES (MC AND DATA) //////

    std::cout << "DATA entries: " << chain.GetEntries()    << std::endl;
    TFile::Open(path_to_MC.Data())->GetObject(mcTreeName.Data(), tree_MC  );
    std::cout << " MC entries: " << tree_MC->GetEntries() << std::endl;
    if (path_to_MC_spec != "") {
        TFile::Open(path_to_MC_spec.Data())->GetObject(mcTreeNameSpec.Data(), tree_MC_spec);
        std::cout << " MC spec entries: " << tree_MC_spec->GetEntries() << std::endl;
    }

    int nVars = sizeof(variables)/sizeof(variables[0]);
    for (int i = 0; i < nVars; ++i){
        TString var = variables[i].expression;
        if (!hasVariableForDraw(&chain, var) || !hasVariableForDraw(tree_MC, var) ||
            (tree_MC_spec && !hasVariableForDraw(tree_MC_spec, var))) {
            std::cout << "[plot_dataMC] Skipping missing variable: " << var << std::endl;
            continue;
        }

        // Create a canvas to draw the histograms
        TCanvas *canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetLeftMargin(0.15);
        canvas->SetTopMargin(0.05);
        canvas->SetRightMargin(0.05);
        // Hide stats boxes globally
        gStyle->SetOptStat(0);

        int nbinsVARhistos = 100;
        double hist_Xhigh      = variables[i].xmax;
        double hist_Xlow       = variables[i].xmin;
        if (var == "Bmass") {
            if (TREE == "ntmix_X3872") {hist_Xlow = 3.6; hist_Xhigh = 4.0;}
            else {hist_Xlow = 5.05; hist_Xhigh = 5.8;}
        }
        double bin_length_MEV  = (hist_Xhigh - hist_Xlow) / nbinsVARhistos;
        
        TString Xlabel ;
        if (var == "Bmass")   {
            if (TREE == "ntKp")        {Xlabel = "m_{J/#Psi K^{+}} [GeV/c^{2}]";}
            else if (TREE == "ntKstar"){Xlabel = "m_{J/#Psi K^{+} #pi^{-}} [GeV/c^{2}]";}
            else if (TREE == "ntphi")  {Xlabel = "m_{J/#Psi K^{+} K^{-}} [GeV/c^{2}]";}
            else if (TREE == "ntmix_X3872")  {Xlabel = "m_{J/#Psi #pi^{+} #pi^{-}} [GeV/c^{2}]";}
        }
        else if (var == "Bpt"){Xlabel = "p_{T} [GeV/c]";}
        else {                 Xlabel = var.Data();}

        // Create histograms
        TH1F *hist_SIG = new TH1F("hist_SIG" , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh); 
        TH1F *hist_BKG = new TH1F("hist_BKG" , Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh);
        TH1F *hist_spec = new TH1F("hist_spec", Form("; %s; Entries / %.3f ", Xlabel.Data(), bin_length_MEV) , nbinsVARhistos, hist_Xlow ,hist_Xhigh);
               
        TString sideband = "1";
        if(TREE == "ntmix_X3872"){ sideband = "(((Bmass > 3.95) & (Bmass < 4.00)) || ((Bmass > 3.75) & (Bmass < 3.80)))";}
        else {sideband = "(Bmass > 5.55)";}

        TString ANYsel = "1"; // Prediction > 0.59
        TString ANA_region = "Bpt > 7.5 && Bpt < 50"; // Bpt > 10 && abs(By) < 1.6
        tree_MC->Draw(Form("%s >> hist_SIG", var.Data()), Form(" %s && %s", ANYsel.Data(), ANA_region.Data()));
        chain.Draw(Form("%s >> hist_BKG", var.Data()), Form(" %s && %s && %s", sideband.Data(), ANYsel.Data(), ANA_region.Data()));

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

        // Customize Histograms

        hist_SIG->SetLineColor(kOrange-3);
        hist_SIG->SetLineWidth(3);
        hist_SIG->SetFillStyle(0);
        if (hist_SIG->Integral() > 0) hist_SIG->Scale(1.0 / hist_SIG->Integral(""));
        hist_SIG->SetMinimum(0.0);

        if (tree_MC_spec) {
            tree_MC_spec->Draw(Form("%s >> hist_spec", var.Data()), Form(" %s && %s", ANYsel.Data(), ANA_region.Data()));
            hist_spec->SetLineWidth(3);
            hist_spec->SetFillStyle(0);
            if (hist_spec->Integral() > 0) hist_spec->Scale(1.0 / hist_spec->Integral(""));
            hist_spec->SetMinimum(0.0);
            hist_spec->SetLineColor(kOrange-2);
        }

        if (hist_BKG->Integral() > 0) hist_BKG->Scale(1.0 / hist_BKG->Integral(""));
        hist_BKG->SetLineColor(kBlue);
        hist_BKG->SetFillColor(kBlue);     
        hist_BKG->SetFillStyle(3358); 
        hist_BKG->SetMinimum(0.0);

        if(1){// set the y-axis maximum if needed
            Double_t max_val = TMath::Max(hist_BKG->GetMaximum(), hist_SIG->GetMaximum()) * 1.1;
            if (tree_MC_spec) max_val = TMath::Max(max_val, hist_spec->GetMaximum() * 1.1);
            hist_SIG->SetMaximum(max_val * 1.1);    // Increase the max range to give some space
            hist_BKG->SetMaximum(max_val * 1.1);
            if (tree_MC_spec) hist_spec->SetMaximum(max_val * 1.1);
        }
        // Customize the Histograms

        // Draw the histograms (background first, signal on top)
        hist_BKG->Draw("HIST");
        hist_SIG->Draw("HIST SAME");
        if (tree_MC_spec) hist_spec->Draw("HIST SAME");
        gPad->Update();

        // Add legend instead of stats boxes
        TString sidebandLatex = sideband;
        if (ANYsel == "1") { ANYsel = ""; }
        else {var += "_SELECTED";}
        sidebandLatex.ReplaceAll("B", "");
        // Split sideband into two lines only when a || operator is present.
        TString sidebandLine1 = sidebandLatex;
        TString sidebandLine2 = "";
        Ssiz_t splitPos = sidebandLatex.First('|');
        if (splitPos != kNPOS) {
            sidebandLine1.Remove(splitPos);
            sidebandLine2 = sidebandLatex(splitPos, sidebandLatex.Length() - splitPos);
        }
        
        TLegend *leg = new TLegend(0.18, 0.62, 0.4, 0.93, NULL, "brNDC");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->SetHeader(Form("#bf{%s}, %s", systemNAME.Data(), ANYsel.Data()));
        if (TREE == "ntmix_X3872") leg->AddEntry(hist_SIG, "X(3872) MC", "l");
        else leg->AddEntry(hist_SIG, Form("%s MC", getPlotParticleLabel(TREE).Data()), "l");
        if (tree_MC_spec) {
            leg->AddEntry(hist_spec, "#Psi(2S) MC", "l");
        }
        leg->AddEntry(hist_BKG, "Data sideband", "f");
        leg->AddEntry((TObject*)0, sidebandLine1.Data(), "");
        if (!sidebandLine2.IsNull()) leg->AddEntry((TObject*)0, sidebandLine2.Data(), "");
        leg->Draw();

        // Save the canvas as an image
        canvas->SaveAs(Form("./presel_STUDY_vars/%s_%s_%s.pdf", TREE.Data(), systemNAME.Data() , var.Data()));

        // Clean up
        delete hist_SIG;
        delete hist_BKG;
        delete hist_spec;
        delete leg;
        delete canvas;
    }
}

int main() {
    plot_dataMC();
    return 0;
}
