#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>

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

TString safePlotTag(TString var)
{
    var.ReplaceAll("abs(", "abs_");
    var.ReplaceAll(")", "");
    var.ReplaceAll("(", "");
    var.ReplaceAll("/", "_over_");
    var.ReplaceAll("*", "_times_");
    var.ReplaceAll("+", "_plus_");
    var.ReplaceAll("-", "_minus_");
    var.ReplaceAll(" ", "");
    return var;
}

void plot2D_dataMC(TString TREE = "ntmix_PSI2S", TString systemNAME = "ppRef", TString cut = "1")
{
    TString dataPath;
    TString mcPath;
    TString dataTree;
    TString mcTree;
    TString outTag = TREE;
    TString mcLabel = "MC";
    const bool isNtmixSignal = (TREE == "ntmix_X3872" || TREE == "ntmix_PSI2S");

    const TString sharingDir = "/eos/user/h/hmarques/RUN3_Data_MC_sharing";

    if (isNtmixSignal) {
        TString sampleDir;
        TString sampleTag;
        if (systemNAME == "ppRef" || systemNAME == "ppRef24") {
            sampleDir = "ppRef24";
            sampleTag = "ppRef";
        } else if (systemNAME == "PbPb23" || systemNAME == "PbPb24") {
            sampleDir = systemNAME;
            sampleTag = systemNAME;
        } else {
            std::cerr << "[plot2D_dataMC] Unsupported X(3872)/Psi(2S) system: "
                      << systemNAME << std::endl;
            return;
        }

        const TString particleTag = TREE == "ntmix_X3872" ? "X3872" : "PSI2S";
        const TString sampleBase = Form("%s/X3872/%s", sharingDir.Data(), sampleDir.Data());
        dataTree = "ntmix";
        mcTree = TREE;
        mcLabel = TREE == "ntmix_X3872" ? "X(3872) MC" : "#psi(2S) MC";
        dataPath = Form("%s/flat_ntmix_%s_DATA.root", sampleBase.Data(), sampleTag.Data());
        mcPath = Form("%s/flat_ntmix_%s_MC_%s.root",
                      sampleBase.Data(), sampleTag.Data(), particleTag.Data());
    } else if (TREE == "ntKp" || TREE == "ntKstar" || TREE == "ntphi") {
        if (systemNAME != "ppRef" && systemNAME != "PbPb23" && systemNAME != "PbPb24") {
            std::cerr << "[plot2D_dataMC] Unsupported B-meson system: " << systemNAME
                      << ". Use ppRef, PbPb23, or PbPb24." << std::endl;
            return;
        }

        const TString sampleBase = Form("%s/Bmesons/%s", sharingDir.Data(), systemNAME.Data());
        dataTree = TREE;
        mcTree = TREE;
        dataPath = Form("%s/flat_%s_%s_DATA.root",
                        sampleBase.Data(), TREE.Data(), systemNAME.Data());
        mcPath = Form("%s/flat_%s_%s_MC.root",
                      sampleBase.Data(), TREE.Data(), systemNAME.Data());
    } else {
        std::cerr << "[plot2D_dataMC] Unsupported tree: " << TREE
                  << ". Use ntmix_X3872, ntmix_PSI2S, ntKp, ntKstar, or ntphi."
                  << std::endl;
        return;
    }

    const TString outDir = Form("2Dplots/%s_%s", outTag.Data(), systemNAME.Data());
    gSystem->mkdir(outDir, true);

    std::cout << "[plot2D_dataMC] TREE   = " << TREE << std::endl;
    std::cout << "[plot2D_dataMC] SYSTEM = " << systemNAME << std::endl;
    std::cout << "[plot2D_dataMC] CUT    = " << cut << std::endl;
    std::cout << "[plot2D_dataMC] DATA   = " << dataPath << " (" << dataTree << ")" << std::endl;
    std::cout << "[plot2D_dataMC] MC     = " << mcPath << " (" << mcTree << ")" << std::endl;

    TFile* dataFile = TFile::Open(dataPath.Data(), "READ");
    if (!dataFile || dataFile->IsZombie()) {
        std::cerr << "[plot2D_dataMC] Cannot open DATA file: " << dataPath << std::endl;
        if (dataFile) { dataFile->Close(); delete dataFile; }
        return;
    }

    TFile* mcFile = TFile::Open(mcPath.Data(), "READ");
    if (!mcFile || mcFile->IsZombie()) {
        std::cerr << "[plot2D_dataMC] Cannot open MC file: " << mcPath << std::endl;
        if (mcFile) { mcFile->Close(); delete mcFile; }
        dataFile->Close();
        delete dataFile;
        return;
    }

    TTree* treeData = nullptr;
    TTree* treeMC = nullptr;
    dataFile->GetObject(dataTree.Data(), treeData);
    mcFile->GetObject(mcTree.Data(), treeMC);
    if (!treeData || !treeMC) {
        if (!treeData) {
            std::cerr << "[plot2D_dataMC] Missing DATA tree " << dataTree
                      << " in " << dataPath << std::endl;
        }
        if (!treeMC) {
            std::cerr << "[plot2D_dataMC] Missing MC tree " << mcTree
                      << " in " << mcPath << std::endl;
        }
        dataFile->Close();
        mcFile->Close();
        delete dataFile;
        delete mcFile;
        return;
    }

    struct MapDef {
        TString yVar;
        TString xVar;
        TString yLabel;
        TString xLabel;
        TString tag;
        int xBins;
        double xMin;
        double xMax;
        int yBins;
        double yMin;
        double yMax;
    };

    double bMassMin = isNtmixSignal ? 3.6 : 5.0;
    double bMassMax = isNtmixSignal ? 4.0 : 6.0;
    if (TREE == "ntmix_PSI2S") {
        bMassMin = 3.6;
        bMassMax = 3.8;
    } else if (TREE == "ntmix_X3872") {
        bMassMin = 3.8;
        bMassMax = 4.0;
    }

    double tktkMassMin = 0.0;
    double tktkMassMax = 1.2;
    if (TREE == "ntphi") {
        tktkMassMin = 1.00;
        tktkMassMax = 1.04;
    } else if (TREE == "ntKstar") {
        tktkMassMin = 0.70;
        tktkMassMax = 1.20;
    } else if (isNtmixSignal) {
        tktkMassMin = 0.20;
        tktkMassMax = 1.20;
    }

    double qValueMin = 0.0;
    double qValueMax = 0.6;
    if (TREE == "ntphi" || TREE == "ntKstar") {
        qValueMin = 0.9;
        qValueMax = 1.6;
    }

    const int nMaps = 17;
    MapDef maps[nMaps] = {
        // Original compact set
        {"abs(By)", "Bpt", "|y|", "p_{T} [GeV/c]", "absY_vs_Bpt", 60, 0.0, 60.0, 48, 0.0, 2.4},
        {"Bpt", "Bmass", "p_{T} [GeV/c]", "B mass [GeV/c^{2}]", "Bpt_vs_Bmass", 80, bMassMin, bMassMax, 60, 0.0, 60.0},
        {"abs(By)", "Bmass", "|y|", "B mass [GeV/c^{2}]", "absY_vs_Bmass", 80, bMassMin, bMassMax, 48, 0.0, 2.4},
        {"Btktkmass", "Bmass", "m_{trk trk} [GeV/c^{2}]", "B mass [GeV/c^{2}]", "Btktkmass_vs_Bmass", 80, bMassMin, bMassMax, 80, tktkMassMin, tktkMassMax},
        {"Bchi2Prob", "Bmass", "Bchi2Prob", "B mass [GeV/c^{2}]", "Bchi2Prob_vs_Bmass", 80, bMassMin, bMassMax, 80, 0.0, 1.0},
        {"Bujmass", "Bmass", "m_{#mu#mu} [GeV/c^{2}]", "B mass [GeV/c^{2}]", "Bujmass_vs_Bmass", 80, bMassMin, bMassMax, 80, 2.8, 3.4},
        {"BtrkPtimb", "Bmass", "BtrkPtimb", "B mass [GeV/c^{2}]", "BtrkPtimb_vs_Bmass", 80, bMassMin, bMassMax, 80, 0.0, 1.0},
        {"Bnorm_svpvDistance_2D", "Bmass", "Bnorm_svpvDistance_2D", "B mass [GeV/c^{2}]", "Bnorm_svpvDistance_2D_vs_Bmass", 80, bMassMin, bMassMax, 80, 0.0, 20.0},
        {"Btktkmass", "Bnorm_svpvDistance_2D", "m_{trk trk} [GeV/c^{2}]", "Bnorm_svpvDistance_2D", "Btktkmass_vs_Bnorm_svpvDistance_2D", 80, 0.0, 20.0, 80, tktkMassMin, tktkMassMax},
        {"BtrkPtimb", "Bnorm_svpvDistance_2D", "BtrkPtimb", "Bnorm_svpvDistance_2D", "BtrkPtimb_vs_Bnorm_svpvDistance_2D", 80, 0.0, 20.0, 80, 0.0, 1.0},
        {"Btktkmass", "BQvalue", "m_{trk trk} [GeV/c^{2}]", "BQvalue", "Btktkmass_vs_BQvalue", 80, qValueMin, qValueMax, 80, tktkMassMin, tktkMassMax},
        {"BtrkPtimb", "BQvalue", "BtrkPtimb", "BQvalue", "BtrkPtimb_vs_BQvalue", 80, qValueMin, qValueMax, 80, 0.0, 1.0},

        // Btktkmass correlation set
        {"Btktkmass", "BtrkPtimb", "m_{trk trk} [GeV/c^{2}]", "BtrkPtimb", "Btktkmass_vs_BtrkPtimb", 80, 0.0, 1.0, 80, tktkMassMin, tktkMassMax},
        {"Btktkmass", "abs(By)", "m_{trk trk} [GeV/c^{2}]", "|y|", "Btktkmass_vs_absY", 48, 0.0, 2.4, 80, tktkMassMin, tktkMassMax},
        {"Btktkmass", "Bpt", "m_{trk trk} [GeV/c^{2}]", "p_{T} [GeV/c]", "Btktkmass_vs_Bpt", 60, 0.0, 60.0, 80, tktkMassMin, tktkMassMax},
        {"Btktkmass", "Bchi2Prob", "m_{trk trk} [GeV/c^{2}]", "Bchi2Prob", "Btktkmass_vs_Bchi2Prob", 80, 0.0, 1.0, 80, tktkMassMin, tktkMassMax},
        {"Btktkmass", "Btrk1Pt", "m_{trk trk} [GeV/c^{2}]", "Btrk1 p_{T} [GeV/c]", "Btktkmass_vs_Btrk1Pt", 60, 0.0, 20.0, 80, tktkMassMin, tktkMassMax}
    };

    TString massWindowCut = "";
    if (TREE == "ntmix_PSI2S") {
        massWindowCut = "(Bmass >= 3.6 && Bmass <= 3.8)";
    } else if (TREE == "ntmix_X3872") {
        massWindowCut = "(Bmass >= 3.8 && Bmass <= 4.0)";
    }
    TString drawCut = cut;
    if (!massWindowCut.IsNull()) {
        drawCut = Form("(%s) && %s", cut.Data(), massWindowCut.Data());
    }

    TCanvas canvas("canvas_2d", "", 850, 720);
    canvas.SetLeftMargin(0.12);
    canvas.SetRightMargin(0.16);
    canvas.SetBottomMargin(0.12);
    gStyle->SetOptStat(0);

    for (int i = 0; i < nMaps; ++i) {
        const MapDef &map = maps[i];
        if (!hasVariableForDraw(treeData, map.yVar) || !hasVariableForDraw(treeData, map.xVar) ||
            !hasVariableForDraw(treeMC, map.yVar) || !hasVariableForDraw(treeMC, map.xVar)) {
            std::cout << "[plot2D_dataMC] Skipping missing map: " << map.yVar << " vs " << map.xVar << std::endl;
            continue;
        }

        TH2F hData(Form("hData_%d", i), Form(";%s;%s;Entries", map.xLabel.Data(), map.yLabel.Data()), map.xBins, map.xMin, map.xMax, map.yBins, map.yMin, map.yMax);
        treeData->Draw(Form("%s:%s>>%s", map.yVar.Data(), map.xVar.Data(), hData.GetName()), drawCut, "goff");
        hData.SetDirectory(nullptr);
        hData.Draw("COLZ");

        TLatex latex;
        latex.SetNDC();
        latex.SetTextFont(42);
        latex.SetTextSize(0.035);
        latex.DrawLatex(0.15, 0.91, Form("%s DATA", systemNAME.Data()));
        canvas.SaveAs(Form("%s/%s_%s_DATA_%s.pdf", outDir.Data(), outTag.Data(), systemNAME.Data(), map.tag.Data()));

        TH2F hMC(Form("hMC_%d", i), Form(";%s;%s;Entries", map.xLabel.Data(), map.yLabel.Data()), map.xBins, map.xMin, map.xMax, map.yBins, map.yMin, map.yMax);
        treeMC->Draw(Form("%s:%s>>%s", map.yVar.Data(), map.xVar.Data(), hMC.GetName()),
                 mcSelection(treeMC, drawCut), "goff");
        hMC.SetDirectory(nullptr);
        hMC.Draw("COLZ");
        latex.DrawLatex(0.15, 0.91, Form("%s %s", systemNAME.Data(), mcLabel.Data()));
        canvas.SaveAs(Form("%s/%s_%s_MC_%s.pdf", outDir.Data(), outTag.Data(), systemNAME.Data(), map.tag.Data()));
    }

    dataFile->Close();
    mcFile->Close();
    delete dataFile;
    delete mcFile;
}

int main()
{
    plot2D_dataMC();
    return 0;
}
