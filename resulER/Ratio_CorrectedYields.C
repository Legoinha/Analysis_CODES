#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include <iostream>

#include "../plotER/aux/parameters.h"

// Usage:
// root -l -b -q 'Ratio_CorrectedYields.C("ntmix","ppRef","Bpt","X3872","Psi2S")'

void Ratio_CorrectedYields(
    TString treenameN = "ntmix",
    TString treenameD = "ntmix",
    TString SYSTEM = "ppRef",
    TString VAR    = "Bpt",
    TString numTag = "_X3872",
    TString denTag = "_Psi2S"
) {

    TString XLabel = "";
    if (VAR == "Bpt") XLabel = "p_{T} [GeV/c]";
    else if (VAR == "By") XLabel = "|y|";   

    gSystem->Exec("mkdir -p output/");
    gStyle->SetOptStat(0);

    if (treenameN != "ntmix") {
        numTag = "";
        denTag = "";
    }
    TString numFile = Form("./../effER/output/ROOTs/%s_%s_%s_CorrectedYields%s.root", treenameN.Data(), SYSTEM.Data(), VAR.Data(), numTag.Data());
    TString denFile = Form("./../effER/output/ROOTs/%s_%s_%s_CorrectedYields%s.root", treenameD.Data(), SYSTEM.Data(), VAR.Data(), denTag.Data());

    TFile *fNum = TFile::Open(numFile, "READ");
    TFile *fDen = TFile::Open(denFile, "READ");
    TH1D *hNum = (TH1D*)fNum->Get("hYieldCorr");
    TH1D *hDen = (TH1D*)fDen->Get("hYieldCorr");

    TH1D *hRatio = (TH1D*)hNum->Clone("hRatio");
    hRatio->SetDirectory(nullptr);
    hRatio->SetTitle(Form("; %s; Ratio", XLabel.Data()));
    if (treenameN == "ntphi" && treenameD == "ntKp") hRatio->SetTitle(Form("; %s; fs/fu", XLabel.Data()));
    hRatio->Divide(hDen);
    if (treenameN == "ntphi" && treenameD == "ntKp") hRatio->Scale( (1/BR_B0s_jpsiphi_mumuKK) / (1/BR_Bp_jpsiK_mumuK ));

    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);

    if (treenameN == "ntmix") { hRatio->GetYaxis()->SetRangeUser(0.01, .10);}

    TCanvas *c = new TCanvas("cRatio", "ratio", 700, 600);
    c->SetLeftMargin(0.15);
    hRatio->GetYaxis()->SetTitleOffset(1.6);
    hRatio->Draw("E1");

    TString outPdf = Form("output/%s_OVER_%s_%s_%s_Ratio.pdf", treenameN.Data(), treenameD.Data(), SYSTEM.Data(), VAR.Data());
    c->SaveAs(outPdf);

    delete c;
    fNum->Close();
    fDen->Close();
}
