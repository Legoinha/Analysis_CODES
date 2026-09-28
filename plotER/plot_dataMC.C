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

// root -l -b -q 'plot_dataMC.C("ntmix_X3872","ppRef")'
// root -l -b -q 'plot_dataMC.C("ntmix_X3872","PbPb23")'
// root -l -b -q 'plot_dataMC.C("ntmix_X3872","PbPb18")'
// root -l -b -q 'plot_dataMC.C("ntKp","PbPb23")'    // B+
// root -l -b -q 'plot_dataMC.C("ntKstar","PbPb23")' // B0
// root -l -b -q 'plot_dataMC.C("ntphi","PbPb23")'   // Bs

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

TString mcSelection(TTree *tree, const TString &selection)
{
    if (!tree->GetBranch("pThatreweight")) return selection;
    return Form("(%s) * pThatreweight", selection.Data());
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
        {"nChargedTracks",             0.0, 10000.0},
        {"CentBin",                    0.0, 100.0},
        {"Bmass",                      3.6,  4.0},
        {"Bpt",                        0.0, 50.0},
        {"abs(By)",                    0.0,  2.4},
        {"Bchi2Prob",                  0.0,  1.0},
        {"Btrk1dR",                    0.0,  .5},
        {"Btrk2dR",                    0.0,  .5},
        {"BtrkPtimb",                  0.0,  1.0},
        {"Btktkpt",                    0.0, 10.0},
        {"Bujmass",                    2.9,  3.3},
        {"BujvProb",                   0.0,  1.0},
        {"Bnorm_svpvDistance_2D",      0.0,  10.0},
        {"BsvpvDistance_2D",           0.0,  0.25},
        {"BsvpvDisErr_2D",             0.0,  0.05},
        {"BQvalue",                    0.0,  0.6},
        {"Bnorm_trk1Dxy",             -7.0,  7.0},
        {"Bnorm_trk2Dxy",             -5.0,  5.0},
        {"Balpha",                     0.0,  3.2},
        {"Bdtheta",                    0,  3.2},
        {"Bcos_dtheta",               -1.0,  1.0},
        {"Btktkmass",                  0.3,  1},
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
        {"Btrk1Dxy1",                 -0.05,  0.05},
        {"Btrk2Dxy1",                 -0.05,  0.05},
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
        {"Bnorm_trk2Dz",              -5.0,  5.0},
        {"BtrkLeadPt",                  0.0, 10.0},
        {"BtrkSubPt",                   0.0, 10.0},
        {"BtrkLeadPtFrac",              0.0,  1.0},
        {"BtrkSubPtFrac",               0.0,  1.0},
        {"BtrkMaxdR",                   0.0,  0.5},
        {"BtrkMindR",                   0.0,  0.5},
        {"BtrkMaxAbsEta",               0.0,  2.4},
        {"BmuLeadPt",                   0.0, 25.0},
        {"BmuSubPt",                    0.0, 25.0},
        {"BmuMaxdR",                    0.0,  1.0},
        {"BmuMindR",                    0.0,  1.0},
        {"BmuMaxAbsEta",                0.0,  2.4},
        {"Prediction",                 0.0,  1.0}
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
    const TString sharingDir = "/eos/user/h/hmarques/RUN3_Data_MC_sharing";
    if (TREE == "ntmix_X3872") {
        TString sampleDir;
        TString sampleTag;
        dataTreeName = "ntmix";
        mcTreeName = "ntmix_X3872";
        mcTreeNameSpec = "ntmix_PSI2S";

        if (systemNAME == "ppRef") {
            sampleDir = "ppRef24";
            sampleTag = "ppRef";
            const TString sampleBase = Form("%s/X3872/%s", sharingDir.Data(), sampleDir.Data());
            path_to_data = Form("%s/flat_ntmix_%s_DATA.root", sampleBase.Data(), sampleTag.Data());
            path_to_MC = Form("%s/flat_ntmix_%s_MC_X3872.root", sampleBase.Data(), sampleTag.Data());
            path_to_MC_spec = Form("%s/flat_ntmix_%s_MC_PSI2S.root", sampleBase.Data(), sampleTag.Data());
        } else if (systemNAME == "PbPb23") {
            const TString sampleBase = Form("%s/X3872/PbPb23", sharingDir.Data());
            path_to_data = sampleBase + "/flat_ntmix_PbPb23_DATA.root";
            path_to_MC = sampleBase + "/flat_ntmix_PbPb23_MC_X3872.root";
            path_to_MC_spec = sampleBase + "/flat_ntmix_PbPb23_MC_PSI2S.root";
        } else if (systemNAME == "PbPb18") {
            const TString sampleBase = Form("%s/X3872/legacy_run2/REflated", sharingDir.Data());
            path_to_data = sampleBase + "/flat_ntmix_PbPb18_DATA.root";
            path_to_MC = sampleBase + "/flat_ntmix_PbPb18_MC_X3872.root";
            path_to_MC_spec = sampleBase + "/flat_ntmix_PbPb18_MC_PSI2S.root";
        } else if (systemNAME == "PbPb24") {
            sampleDir = systemNAME;
            sampleTag = systemNAME;
            const TString sampleBase = Form("%s/X3872/%s", sharingDir.Data(), sampleDir.Data());
            path_to_data = Form("%s/flat_ntmix_%s_DATA.root", sampleBase.Data(), sampleTag.Data());
            path_to_MC = Form("%s/flat_ntmix_%s_MC_X3872.root", sampleBase.Data(), sampleTag.Data());
            path_to_MC_spec = Form("%s/flat_ntmix_%s_MC_PSI2S.root", sampleBase.Data(), sampleTag.Data());
        } else {
            std::cerr << "[plot_dataMC] Unsupported X(3872) system: " << systemNAME << std::endl;
            return;
        }
    } else if (TREE == "ntKp" || TREE == "ntKstar" || TREE == "ntphi") {
        
        const TString sampleBase = Form("%s/Bmesons/%s", sharingDir.Data(), systemNAME.Data());
        path_to_data = Form("%s/flat_%s_%s_DATA.root", sampleBase.Data(), TREE.Data(), systemNAME.Data());
        path_to_MC = Form("%s/flat_%s_%s_MC.root", sampleBase.Data(), TREE.Data(), systemNAME.Data());
    } 

    std::cout << "[plot_dataMC] DATA: " << path_to_data << " (" << dataTreeName << ")" << std::endl;
    std::cout << "[plot_dataMC] MC:   " << path_to_MC << " (" << mcTreeName << ")" << std::endl;
    if (!path_to_MC_spec.IsNull()) {
        std::cout << "[plot_dataMC] MC2:  " << path_to_MC_spec << " (" << mcTreeNameSpec << ")" << std::endl;
    }

    TChain chain(dataTreeName.Data());
    chain.Add(path_to_data.Data());
    TFile *file_MC = TFile::Open(path_to_MC.Data(), "READ");
    file_MC->GetObject(mcTreeName.Data(), tree_MC);
    TFile *file_MC_spec = nullptr;
    if (!path_to_MC_spec.IsNull()) {
        file_MC_spec = TFile::Open(path_to_MC_spec.Data(), "READ");
        file_MC_spec->GetObject(mcTreeNameSpec.Data(), tree_MC_spec);
    }
    ////// OPEN FILES (MC AND DATA) //////
    ////// OPEN FILES (MC AND DATA) //////

    std::cout << "DATA entries: " << chain.GetEntries() << std::endl;
    std::cout << " MC entries: " << tree_MC->GetEntries() << std::endl;
    if (tree_MC_spec) {
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
            else if (TREE == "ntmix_X3872") {Xlabel = "m_{J/#Psi #pi^{+} #pi^{-}} [GeV/c^{2}]";}
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

        TString ANYsel = "1";
        TString ANA_region = "Bpt > 7.5 && Bpt < 50";
        if (TREE == "ntmix_X3872" && (systemNAME == "ppRef")) {
            ANA_region = "(Bpt > 7.5 && Bpt < 50) && "
                         "Btrk1dR < 0.5 && Btrk2dR < 0.5 &&  "
                         "BLxy*(Bmass/Bpt) < 0.05"; // enforce prompt component
        }
        if (TREE == "ntmix_X3872" && systemNAME == "PbPb23") {
            ANA_region = "(Bpt > 15 && Bpt < 50) && (abs(By) < 2.4) && "
                         "(BQvalue < 0.15) && (Bchi2Prob > 0.05)";
        }
        else if (TREE == "ntmix_X3872" && systemNAME == "PbPb18") {
            ANA_region = "(Bpt > 15 && Bpt < 50) && (abs(By) < 2.4) && "
                         "(BQvalue < 0.15) && (Bchi2Prob > 0.1)";
        }
        const TString mcCut = Form("%s && %s", ANYsel.Data(), ANA_region.Data());
        tree_MC->Draw(Form("%s >> hist_SIG", var.Data()), mcSelection(tree_MC, mcCut));
        chain.Draw(Form("%s >> hist_BKG", var.Data()), Form(" %s && %s && %s", sideband.Data(), ANYsel.Data(), ANA_region.Data()));
        const double selectedSignalEntries = hist_SIG->GetEntries();
        const double selectedDataEntries = hist_BKG->GetEntries();

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

        // Customize Histograms

        hist_SIG->SetLineColor(kOrange-3);
        hist_SIG->SetLineWidth(3);
        hist_SIG->SetFillStyle(0);
        if (hist_SIG->Integral() > 0) hist_SIG->Scale(1.0 / hist_SIG->Integral(""));
        hist_SIG->SetMinimum(0.0);

        if (tree_MC_spec) {
            tree_MC_spec->Draw(Form("%s >> hist_spec", var.Data()),
                               mcSelection(tree_MC_spec, mcCut));
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
        if (TREE == "ntmix_X3872" && var == "Bmass") {
            leg->AddEntry(hist_SIG, Form("X(3872) MC, N_{sig} = %.0f", selectedSignalEntries), "l");
        }
        else if (TREE == "ntmix_X3872") leg->AddEntry(hist_SIG, "X(3872) MC", "l");
        else leg->AddEntry(hist_SIG, Form("%s MC", getPlotParticleLabel(TREE).Data()), "l");
        if (tree_MC_spec) {
            if (var == "Bmass") {
                leg->AddEntry(hist_spec, Form("#Psi(2S) MC, N_{sig} = %.0f", hist_spec->GetEntries()), "l");
            }
            else leg->AddEntry(hist_spec, "#Psi(2S) MC", "l");
        }
        if (var == "Bmass") {
            leg->AddEntry(hist_BKG, Form("Data sideband, N_{data} = %.0f", selectedDataEntries), "f");
        }
        else leg->AddEntry(hist_BKG, "Data sideband", "f");
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

    file_MC->Close();
    delete file_MC;
    if (file_MC_spec) {
        file_MC_spec->Close();
        delete file_MC_spec;
    }
}

int main() {
    plot_dataMC();
    return 0;
}
