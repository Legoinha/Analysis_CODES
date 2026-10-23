#ifndef VALIDATION_AUX_H
#define VALIDATION_AUX_H

#include <TString.h>
#include <TH1D.h>
#include <TMath.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

struct VarCfgSignal {
    TString expr;
    TString title;
    int nbins;
    double xmin;
    double xmax;
    bool absVal;
    TString tag = "";
};

struct AgreementMetrics {
    double ksDistance = -1.0;
    double ksPValue = -1.0;
    double chi2 = 0.0;
    int ndf = 0;
    double chi2PValue = -1.0;
    int binsUsed = 0;
};

struct ReweightInput {
    TString expr;
    TString tag;
    TH1D* hist = nullptr;
};

static constexpr int kNBins = 10;


static TString validationWeightExpression(const VarCfgSignal& var)
{
    return var.absVal ? Form("abs(%s)", var.expr.Data()) : var.expr;
}

static TString validationVariableTag(const VarCfgSignal& var)
{
    return var.tag.IsNull() ? var.expr : var.tag;
}

static double lookupWeight1D(const TH1D* hist, double value)
{
    const int bin = hist->GetXaxis()->FindBin(value);
    return hist->GetBinContent(bin);
}


static std::vector<VarCfgSignal> getSignalVars(TString treeName, TString systemTag,
                                                TString fitVariable,
                                                double fitVariableMin, double fitVariableMax)
{
    std::vector<VarCfgSignal> vars = {
        {"By", ";|y|;", kNBins, 0.0, 1.6, true},
        {"nChargedTracks", ";N_{trk}^{corr};", kNBins, 0.0, 100.0, false},
        {"CentBin", ";Centrality (%);", kNBins, 0.0, 100.0, false},
        {"Bchi2Prob", ";Bchi2Prob;", kNBins, 0.0, 1.0, false},
        {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, .5, false},
        {"Btrk2dR", ";Btrk2dR;", kNBins, 0.0, .5, false},
        {"BtrkPtimb", ";BtrkPtimb;", kNBins, 0.0, 1.0, false},
        {"Btktkpt", ";Btktkpt;", kNBins, 0.0, 5.0, false},
        {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.95, 3.25, false},
        {"BujvProb", ";BujvProb;", kNBins, 0.0, 1.0, false},
        {"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 35.0, false},
        {"BsvpvDistance_2D", ";BsvpvDistance_2D;", kNBins, 0.0, 1.0, false},
        {"BsvpvDisErr_2D", ";BsvpvDisErr_2D;", kNBins, 0.0, 0.025, false},
        {"BQvalue", ";BQvalue;", kNBins, 0.0, 0.3, false},
        {"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -5.0, 5.0, false},
        {"Bnorm_trk2Dxy", ";Bnorm_trk2Dxy;", kNBins, -5.0, 5.0, false},
        {"Balpha", ";Balpha;", kNBins, 0.0, 3.0, false},
        {"Bcos_dtheta", ";Bcos_dtheta;", kNBins, -1.0, 1.0, false},
        {"Btktkmass", ";Btktkmass;", kNBins, 0.0, 2.0, false},
        {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.0, 5.0, false},
        {"Btrk2Pt", ";Btrk2Pt;", kNBins, 0.0, 5.0, false},
        {"BtrkLeadPt", ";Leading-track p_{T};", kNBins, 0.0, 5.0, false},
        {"BtrkSubPt", ";Subleading-track p_{T};", kNBins, 0.0, 2.5, false},
        {"BtrkMaxdR", ";max #DeltaR(track,B);", kNBins, 0.0, 0.5, false},
        {"BtrkMindR", ";min #DeltaR(track,B);", kNBins, 0.0, 0.5, false},
        {"BtrkMaxAbsEta", ";max |#eta(track)|;", kNBins, 0.0, 2.4, false},
        {"Btrk1Eta", ";Btrk1Eta;", kNBins, -2.4, 2.4, false},
        {"Btrk2Eta", ";Btrk2Eta;", kNBins, -2.4, 2.4, false},
        {"Btrk1Phi", ";Btrk1Phi;", kNBins, -3.2, 3.2, false},
        {"Btrk2Phi", ";Btrk2Phi;", kNBins, -3.2, 3.2, false},
        {"Btrk1PtErr", ";Btrk1PtErr;", kNBins, 0.0, 0.1, false},
        {"Btrk2PtErr", ";Btrk2PtErr;", kNBins, 0.0, 0.1, false},
        {"BtktkvProb", ";BtktkvProb;", kNBins, 0.0, 1.0, false},
        {"BLxy*(Bmass/Bpt)", ";l_{xy} [cm];", kNBins, -0.2, 0.5, false, "lxy"},
        {"Bmu1pt", ";Bmu1pt;", kNBins, 0.0, 20.0, false},
        {"Bmu2pt", ";Bmu2pt;", kNBins, 0.0, 20.0, false},
        {"BmuLeadPt", ";Leading-muon p_{T};", kNBins, 0.0, 20.0, false},
        {"BmuSubPt", ";Subleading-muon p_{T};", kNBins, 0.0, 15.0, false},
        {"BmuMaxdR", ";max #DeltaR(#mu,J/#psi);", kNBins, 0.0, 1.5, false},
        {"BmuMindR", ";min #DeltaR(#mu,J/#psi);", kNBins, 0.0, 1.5, false},
        {"BmuMaxAbsEta", ";max |#eta(#mu)|;", kNBins, 0.0, 2.4, false},
        {"Bmu1eta", ";Bmu1eta;", kNBins, -2.4, 2.4, false},
        {"Bmu2eta", ";Bmu2eta;", kNBins, -2.4, 2.4, false},
        {"Bmu1phi", ";Bmu1phi;", kNBins, -3.2, 3.2, false},
        {"Bmu2phi", ";Bmu2phi;", kNBins, -3.2, 3.2, false},
        {"Bujpt", ";Bujpt;", kNBins, 0.0, 35.0, false},
        {"Bujeta", ";Bujeta;", kNBins, -2.4, 2.4, false},
        {"Bujphi", ";Bujphi;", kNBins, -3.2, 3.2, false},
        {"Bujlxy", ";Bujlxy;", kNBins, -0.1, 0.1, false},
        {"Btrk1Dz1", ";Btrk1Dz1;", kNBins, -0.5, 0.5, false},
        {"Btrk2Dz1", ";Btrk2Dz1;", kNBins, -0.5, 0.5, false},
        {"Btrk1DzError1", ";Btrk1DzError1;", kNBins, 0.0, 0.1, false},
        {"Btrk2DzError1", ";Btrk2DzError1;", kNBins, 0.0, 0.1, false},
        {"Btrk1Dxy1", ";Btrk1Dxy1;", kNBins, -0.06, 0.06, false},
        {"Btrk2Dxy1", ";Btrk2Dxy1;", kNBins, -0.06, 0.06, false},
        {"Btrk1DxyError1", ";Btrk1DxyError1;", kNBins, 0.0, 0.05, false},
        {"Btrk2DxyError1", ";Btrk2DxyError1;", kNBins, 0.0, 0.05, false},
        {"Btktketa", ";Btktketa;", kNBins, -2.4, 2.4, false},
        {"Btktkphi", ";Btktkphi;", kNBins, -3.2, 3.2, false},
        {"Btktky", ";Btktky;", kNBins, -2.4, 2.4, false},
        {"Bdoubletpt", ";Bdoubletpt;", kNBins, 0.0, 10.0, false},
        {"Bdoubleteta", ";Bdoubleteta;", kNBins, -2.4, 2.4, false},
        {"Bdoubletphi", ";Bdoubletphi;", kNBins, -3.2, 3.2, false},
        {"Bdoublety", ";Bdoublety;", kNBins, -2.4, 2.4, false},
        {"Bnorm_trk1Dz", ";Bnorm_trk1Dz;", kNBins, -5.0, 5.0, false},
        {"Bnorm_trk2Dz", ";Bnorm_trk2Dz;", kNBins, -5.0, 5.0, false}
    };

    if (treeName == "ntKp") {
        vars = {
            {"By", ";|y|;", kNBins, 0.0, 2.4, true},
            {"nChargedTracks", ";N_{trk}^{corr};", kNBins, 0.0, 100.0, false},
            {"CentBin", ";Centrality (%);", kNBins, 0.0, 100.0, false},
            {"Prediction", ";Classifier score;", kNBins, 0.0, 1.0, false},
            {"Bchi2Prob", ";Bchi2Prob;", kNBins, 0.0, 1.0, false},
            {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, 1.5, false},
            {"Btrk1Dz1", ";Btrk1Dz1;", kNBins, -0.5, 0.5, false},
            {"Bnorm_trk1Dz", ";Bnorm_trk1Dz;", kNBins, -20.0, 20.0, false},
            {"Bmu1phi", ";Bmu1phi;", kNBins, -3.2, 3.2, false},
            {"Btrk1Eta", ";Btrk1Eta;", kNBins, -2.4, 2.4, false},
            {"Btrk1Phi", ";Btrk1Phi;", kNBins, -3.2, 3.2, false},
            {"Bmu1pt", ";Bmu1pt;", kNBins, 0.0, 20.0, false},
            {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.9, 3.3, false},
            {"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 50.0, false},
            {"BsvpvDistance_2D", ";BsvpvDistance_2D;", kNBins, 0.0, 1.0, false},
            {"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -25.0, 25.0, false},
            {"Balpha", ";Balpha;", kNBins, 0.0, 1.0, false},
            {"Bcos_dtheta", ";Bcos_dtheta;", kNBins, 0.995, 1.0, false},
            {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.0, 15.0, false},
            {"Bujpt", ";Bujpt;", kNBins, 0.0, 35.0, false},
            {"Bujeta", ";Bujeta;", kNBins, -2.4, 2.4, false},
            {"Bujphi", ";Bujphi;", kNBins, -3.2, 3.2, false},
            {"Bujlxy", ";Bujlxy;", kNBins, 0.0, 0.5, false}
        };
    }

    if (treeName.BeginsWith("ntmix")) {
        for (auto& v : vars) {
            if (v.expr == "BmuMaxdR") { v.xmin = 0.0; v.xmax = 1.0; }
            if (v.expr == "BmuMindR") { v.xmin = 0.0; v.xmax = 0.8; }
            if (v.expr == "BQvalue") { v.xmin = 0.0; v.xmax = 0.2; }
            if (v.expr == "BtrkMaxdR") { v.xmin = 0.0; v.xmax = 0.8; }
            if (v.expr == "BsvpvDistance_2D") { v.xmin = 0.0; v.xmax = 1.0; }
        }
    }

    if (treeName == "ntmix_X3872") {
        for (auto& v : vars) {
            if (v.expr == "Btktkmass") {
                v.xmin = 0.55;
                v.xmax = 0.80;
            }
            if (v.expr == "Bmass") { v.xmin = 3.6; v.xmax = 4.0; }
        }

    } else if (treeName == "ntmix_PSI2S") {
        for (auto& v : vars) {
            if (v.expr == "Btktkmass") { v.xmin = 0.42; v.xmax = 0.61; }
            if (v.expr == "Bmass") { v.xmin = 3.6; v.xmax = 4.0; }
        }
    }

    if (systemTag.BeginsWith("PbPb")) {
        for (auto& v : vars) {
            if (v.expr == "By" && treeName.BeginsWith("ntmix")) {
                v.xmin = 0.0;
                v.xmax = 1.6;
            }
            if (v.expr == "nChargedTracks") {
                v.xmin = 0.0;
                v.xmax = 10000.0;
            }
        }
    }

    const bool isBMeson = treeName == "ntKp" || treeName == "ntphi" || treeName == "ntKstar";
    if (isBMeson) {
        for (auto& v : vars) {
            if (v.expr == "Bchi2Prob") { v.xmin = 0.0; v.xmax = 1.0; }
            if (v.expr == "Bdoubletpt") { v.xmin = 0.0; v.xmax = 10.0; }
            if (v.expr == "Bnorm_svpvDistance_2D") { v.xmin = 0.0; v.xmax = 50.0; }
            if (v.expr == "BsvpvDistance_2D") { v.xmin = 0.0; v.xmax = 1.0; }
            if (v.expr == "Bnorm_trk1Dxy" || v.expr == "Bnorm_trk2Dxy") {
                v.xmin = -25.0; v.xmax = 25.0;
            }
            if (v.expr == "Bnorm_trk1Dz" || v.expr == "Bnorm_trk2Dz") {
                v.xmin = -20.0; v.xmax = 20.0;
            }
            if (v.expr == "Btrk1dR" || v.expr == "Btrk2dR") { v.xmin = 0.0; v.xmax = 1.5; }
            if (v.expr == "Btrk1Pt" || v.expr == "Btrk2Pt") { v.xmin = 0.0; v.xmax = 15.0; }
            if (v.expr == "BtrkLeadPt") { v.xmin = 0.0; v.xmax = 15.0; }
            if (v.expr == "BtrkSubPt") { v.xmin = 0.0; v.xmax = 2.5; }
            if (v.expr == "BtrkMaxdR" || v.expr == "BtrkMindR") { v.xmin = 0.0; v.xmax = 1.5; }
            if (v.expr == "Bujlxy") { v.xmin = 0.0; v.xmax = 0.5; }
        }
    }
    vars.erase(std::remove_if(vars.begin(), vars.end(), [&](const VarCfgSignal& var) {
        return var.expr == fitVariable || validationVariableTag(var) == fitVariable;
    }), vars.end());
    const TString fitAxisTitle = fitVariable == "Bpt"
        ? ";p_{T} [GeV/c];" : Form(";%s;", fitVariable.Data());
    vars.insert(vars.begin(), {fitVariable, fitAxisTitle, kNBins,
                              fitVariableMin, fitVariableMax, fitVariable == "By"});

    return vars;
}

static TH1D* makeWeightHist(const TH1D* hData, const TH1D* hMC, const TString& name)
{
    TH1D* hWeight = (TH1D*)hData->Clone(name);
    hWeight->Reset("ICES");
    hWeight->SetTitle(hData->GetTitle());
    hWeight->GetXaxis()->SetTitle(hData->GetXaxis()->GetTitle());
    hWeight->GetYaxis()->SetTitle("Data / MC");
    for (int i = 0; i <= hWeight->GetNbinsX() + 1; ++i) {
        const double data = hData->GetBinContent(i);
        const double mc = hMC->GetBinContent(i);
        const double dataErr = hData->GetBinError(i);
        const double mcErr = hMC->GetBinError(i);
        // No MC candidate falls in an empty MC bin, so its weight is never applied: it is 0, whatever the data.
        if (mc == 0.0) {
            if (data != 0.0) {
                std::cerr << "[validation] WARNING: " << name << " bin " << i << ": MC is empty while data is "
                          << data << " +- " << dataErr << "; weight set to 0" << std::endl;
            }
            hWeight->SetBinContent(i, 0.0);
            hWeight->SetBinError(i, 0.0);
            continue;
        }
        const double ratio = data / mc;
        const double ratioErr = std::sqrt(std::pow(dataErr / mc, 2) +
                                          std::pow(data * mcErr / (mc * mc), 2));
        hWeight->SetBinContent(i, ratio);
        hWeight->SetBinError(i, ratioErr);
    }
    return hWeight;
}

static AgreementMetrics computeAgreementMetrics1D(const TH1D* hData, const TH1D* hMC)
{
    AgreementMetrics metrics;
    const double dataIntegral = hData->Integral();
    const double mcIntegral = hMC->Integral();
    if (dataIntegral > 0.0 && mcIntegral > 0.0) {
        metrics.ksPValue = hData->KolmogorovTest(hMC);
        metrics.ksDistance = hData->KolmogorovTest(hMC, "M");
    }

    for (int i = 1; i <= hData->GetNbinsX(); ++i) {
        const double data = hData->GetBinContent(i);
        const double mc = hMC->GetBinContent(i);
        const double dataErr = hData->GetBinError(i);
        const double mcErr = hMC->GetBinError(i);
        const double variance = dataErr * dataErr + mcErr * mcErr;
        if (!(variance > 0.0)) continue;

        const double diff = data - mc;
        metrics.chi2 += diff * diff / variance;
        ++metrics.binsUsed;
    }

    metrics.ndf = metrics.binsUsed - 1;
    if (metrics.ndf > 0) metrics.chi2PValue = TMath::Prob(metrics.chi2, metrics.ndf);
    return metrics;
}

#endif
