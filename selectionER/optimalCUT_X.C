#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

TString scoredSampleSystem(TString system) {
    TString key = system;
    key.ToLower();

    if (key == "ppref" || key == "ppref24") return "ppRef";
    if (key == "pbpb23") return "PbPb23";
    if (key == "pbpb24") return "PbPb24";
    if (key == "pbpb") return "PbPb";
    return "";
}

void optimalCUT_X(TString system = "ppRef") {
    gStyle->SetOptStat(0);

    TString sampleSystem = scoredSampleSystem(system);
    if (sampleSystem.IsNull()) {
        std::cerr << "Unknown system '" << system << "'. Use ppRef, ppRef24, PbPb23, PbPb24, or PbPb." << std::endl;
        return;
    }

    TString sampleDir = "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples";
    TString dataPath = Form("%s/flat_ntmix_%s_scored_DATA.root", sampleDir.Data(), sampleSystem.Data());
    TString mcXPath = Form("%s/flat_ntmix_%s_scored_MC_X3872.root", sampleDir.Data(), sampleSystem.Data());

    std::cout << "Reading scored data sample: " << dataPath << std::endl;
    std::cout << "Reading scored X(3872) MC sample: " << mcXPath << std::endl;

    TFile* fileData = TFile::Open(dataPath);
    TFile* fileX = TFile::Open(mcXPath);
    if (!fileData || fileData->IsZombie() || !fileX || fileX->IsZombie()) {
        std::cerr << "Could not open one of the scored input files." << std::endl;
        return;
    }

    TTree *data = nullptr, *mcX = nullptr;
    fileData->GetObject("ntmix", data);
    fileX->GetObject("ntmix_X3872", mcX);
    if (!data || !mcX) {
        std::cerr << "Could not find required trees: data=ntmix, MC=ntmix_X3872." << std::endl;
        return;
    }

    TString sideband = "((Bmass > 3.92 && Bmass < 3.95) || (Bmass > 3.82 && Bmass < 3.84))";
    double Fs_X = 414 / 36074.00;
    double bkgScale = (3500) / 6722.00;
    double bestThr = 0., bestFom = -1., ymax = 0.;

    TH1F htmp("htmp", "", 1, 0, 1);
    TGraph gX;
    int ip = 0;

    TString preCut = "abs(By) < 1.6 && Bpt > 10 && BQvalue < 0.1";
    if (sampleSystem.Contains("PbPb")) preCut += " && CentBin > 10";

    for (double thr = 0.; thr <= 1.0001; thr += 0.01) {
        TString sel = Form("(%s) && (xgb_score > %.3f)", preCut.Data(), thr);
        TString selData = Form("(%s) && (%s) && (xgb_score > %.3f)", sideband.Data(), preCut.Data(), thr);

        htmp.Reset();
        double sx = mcX->Project("htmp", "xgb_score", sel);
        htmp.Reset();
        double b = data->Project("htmp", "xgb_score", selData) ;

        double denX = Fs_X * sx + bkgScale * b;
        double fomX = (denX > 0.) ? Fs_X * sx / TMath::Sqrt(denX) : 0.;
        if (!std::isfinite(fomX)) fomX = 0.;
        fomX = std::round(fomX * 10.) / 10.;

        gX.SetPoint(ip, thr, fomX);
        if (fomX > bestFom) {
            bestFom = fomX;
            bestThr = thr;
        }
        if (fomX > ymax) ymax = fomX;
        // Print each iteration with Bkg*fb and Signal*fs
        std::cout << Form("thr = %.3f, FOM_X = %.2f, Bkg_DATA = %.2f, Signal_MC = %.2f", thr, fomX, b, sx) << std::endl;
        ip++;
    }

    TCanvas c("c", "", 800, 600);
    TH1F* frame = c.DrawFrame(0., 0., 1., 1.2 * ymax);
    frame->SetTitle(Form("%s; xgb_score; FOM", sampleSystem.Data()));

    gX.SetMarkerStyle(20);
    gX.Draw("LP SAME");

    TLine bestLine(bestThr, 0., bestThr, bestFom);
    bestLine.SetLineStyle(2);
    bestLine.SetLineWidth(2);
    bestLine.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);
    label.DrawLatex(0.16, 0.86, Form("Best threshold = %.2f", bestThr));
    label.DrawLatex(0.16, 0.81, Form("FOM = %.1f", bestFom));

    c.SaveAs(Form("bdt_optimization_%s.pdf", sampleSystem.Data()));

    std::cout << Form("%s: best X(3872) thr. = %.2f, FOM = %.1f", system.Data(), bestThr, bestFom) << std::endl;
}

int main(int argc, char** argv) {
    TString system = "ppRef";
    if (argc > 1) system = argv[1];
    optimalCUT_X(system);
    return 0;
}
