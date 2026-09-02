#ifndef VALIDATION_AUX_H
#define VALIDATION_AUX_H

#include <TString.h>
#include <TH1D.h>
#include <TMath.h>
#include <vector>
#include <cmath>

struct VarCfgSignal {
    TString expr;
    TString title;
    int nbins;
    double xmin;
    double xmax;
    bool absVal;
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

static constexpr int kNBins = 15;

static double lookupWeight1D(const TH1D* hist, double value)
{
    const int bin = hist->GetXaxis()->FindBin(value);
    if (bin < 1 || bin > hist->GetNbinsX()) return 1.0;
    return hist->GetBinContent(bin);
}


static std::vector<VarCfgSignal> getSignalVars(TString treeName)
{
    std::vector<VarCfgSignal> vars = {
        {"PVnchi2", ";PVnchi2;", kNBins, 0.0, 1.0, false},
        {"nChargedTracks", ";nChargedTracks;", kNBins, 0.0, 200.0, false},
        {"Bpt", ";p_{T} [GeV/c];", kNBins, 7.5, 50.0, false},
        {"By", ";|y|;", kNBins, 0.0, 2.4, true},
        {"Prediction", ";Prediction;", kNBins, 0.55, 1.0, false},
        {"Bchi2Prob", ";Bchi2Prob;", kNBins, 0.0, 1.0, false},
        {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, .5, false},
        {"Btrk2dR", ";Btrk2dR;", kNBins, 0.0, .5, false},
        {"BtrkPtimb", ";BtrkPtimb;", kNBins, 0.0, 1.0, false},
        {"Btktkpt", ";Btktkpt;", kNBins, 0.0, 10.0, false},
        {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.9, 3.3, false},
        {"BujvProb", ";BujvProb;", kNBins, 0.0, 1.0, false},
        {"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 15.0, false},
        {"BsvpvDistance_2D", ";BsvpvDistance_2D;", kNBins, 0.0, 0.25, false},
        {"BsvpvDisErr_2D", ";BsvpvDisErr_2D;", kNBins, 0.0, 0.025, false},
        {"BQvalue", ";BQvalue;", kNBins, 0.0, 0.6, false},
        {"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -5.0, 5.0, false},
        {"Bnorm_trk2Dxy", ";Bnorm_trk2Dxy;", kNBins, -5.0, 5.0, false},
        {"Balpha", ";Balpha;", kNBins, 0.0, 3.2, false},
        {"Bcos_dtheta", ";Bcos_dtheta;", kNBins, -1.0, 1.0, false},
        {"Btktkmass", ";Btktkmass;", kNBins, 0.0, 2.0, false},
        {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.0, 5.0, false},
        {"Btrk2Pt", ";Btrk2Pt;", kNBins, 0.0, 5.0, false},
        {"Btrk1Eta", ";Btrk1Eta;", kNBins, -2.4, 2.4, false},
        {"Btrk2Eta", ";Btrk2Eta;", kNBins, -2.4, 2.4, false},
        {"Btrk1Phi", ";Btrk1Phi;", kNBins, -3.2, 3.2, false},
        {"Btrk2Phi", ";Btrk2Phi;", kNBins, -3.2, 3.2, false},
        {"Btrk1PtErr", ";Btrk1PtErr;", kNBins, 0.0, 0.1, false},
        {"Btrk2PtErr", ";Btrk2PtErr;", kNBins, 0.0, 0.1, false},
        {"BtktkvProb", ";BtktkvProb;", kNBins, 0.0, 1.0, false},
        {"BLxy", ";BLxy;", kNBins, -0.1, 0.1, false},
        {"Bmu1pt", ";Bmu1pt;", kNBins, 0.0, 25.0, false},
        {"Bmu2pt", ";Bmu2pt;", kNBins, 0.0, 25.0, false},
        {"Bmu1eta", ";Bmu1eta;", kNBins, -2.4, 2.4, false},
        {"Bmu2eta", ";Bmu2eta;", kNBins, -2.4, 2.4, false},
        {"Bmu1phi", ";Bmu1phi;", kNBins, -3.2, 3.2, false},
        {"Bmu2phi", ";Bmu2phi;", kNBins, -3.2, 3.2, false},
        {"Bujpt", ";Bujpt;", kNBins, 0.0, 50.0, false},
        {"Bujeta", ";Bujeta;", kNBins, -2.4, 2.4, false},
        {"Bujphi", ";Bujphi;", kNBins, -3.2, 3.2, false},
        {"Bujlxy", ";Bujlxy;", kNBins, -0.1, 0.1, false},
        {"Btrk1Dz1", ";Btrk1Dz1;", kNBins, -0.5, 0.5, false},
        {"Btrk2Dz1", ";Btrk2Dz1;", kNBins, -0.5, 0.5, false},
        {"Btrk1DzError1", ";Btrk1DzError1;", kNBins, 0.0, 0.1, false},
        {"Btrk2DzError1", ";Btrk2DzError1;", kNBins, 0.0, 0.1, false},
        {"Btrk1Dxy1", ";Btrk1Dxy1;", kNBins, -0.06, 0.06, false},
        {"Btrk2Dxy1", ";Btrk2Dxy1;", kNBins, -0.06, 0.06, false},
        {"Btrk1DxyError1", ";Btrk1DxyError1;", kNBins, 0.0, 0.1, false},
        {"Btrk2DxyError1", ";Btrk2DxyError1;", kNBins, 0.0, 0.1, false},
        {"Btktketa", ";Btktketa;", kNBins, -2.4, 2.4, false},
        {"Btktkphi", ";Btktkphi;", kNBins, -3.2, 3.2, false},
        {"Btktky", ";Btktky;", kNBins, -2.4, 2.4, false},
        {"Bdoubletpt", ";Bdoubletpt;", kNBins, 0.0, 15.0, false},
        {"Bdoubleteta", ";Bdoubleteta;", kNBins, -2.4, 2.4, false},
        {"Bdoubletphi", ";Bdoubletphi;", kNBins, -3.2, 3.2, false},
        {"Bdoublety", ";Bdoublety;", kNBins, -2.4, 2.4, false},
        {"Bnorm_trk1Dz", ";Bnorm_trk1Dz;", kNBins, -5.0, 5.0, false},
        {"Bnorm_trk2Dz", ";Bnorm_trk2Dz;", kNBins, -5.0, 5.0, false}
    };

    if (treeName == "ntphi") {
        vars = {
            {"Bpt", ";p_{T} [GeV/c];", kNBins, 7.5, 50.0, false},
            {"By", ";|y|;", kNBins, 0.0, 2.4, true},
            {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, 0.5, false},
            {"Btrk2dR", ";Btrk2dR;", kNBins, 0.0, 0.5, false},
            {"BtrkPtimb", ";BtrkPtimb;", kNBins, 0.0, 1.0, false},
            {"Btktkpt", ";Btktkpt;", kNBins, 0.0, 10.0, false},
            {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.9, 3.3, false},
            {"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 15.0, false},
            {"BQvalue", ";BQvalue;", kNBins, 0.0, 0.6, false},
            {"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -5.0, 5.0, false},
            {"Bnorm_trk2Dxy", ";Bnorm_trk2Dxy;", kNBins, -5.0, 5.0, false},
            {"Balpha", ";Balpha;", kNBins, 0.0, 3.2, false},
            {"Bcos_dtheta", ";Bcos_dtheta;", kNBins, -1.0, 1.0, false},
            {"Btktkmass", ";Btktkmass;", kNBins, 0.0, 2.0, false},
            {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.0, 5.0, false},
            {"Btrk2Pt", ";Btrk2Pt;", kNBins, 0.0, 5.0, false}
        };
    } else if (treeName == "ntKp") {
        vars = {
            {"Bpt", ";p_{T} [GeV/c];", kNBins, 7.5, 50.0, false},
            {"By", ";|y|;", kNBins, 0.0, 2.4, true},
            {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, 1.5, false},
            {"BtrkPtimb", ";BtrkPtimb;", kNBins, 0.0, 1.0, false},
            {"Btktkpt", ";Btktkpt;", kNBins, 0.0, 10.0, false},
            {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.9, 3.3, false},
            {"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 15.0, false},
            {"BsvpvDistance_2D", ";BsvpvDistance_2D;", kNBins, 0.0, 0.25, false},
            {"BQvalue", ";BQvalue;", kNBins, 0.0, 0.6, false},
            {"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -5.0, 5.0, false},
            {"Balpha", ";Balpha;", kNBins, 0.0, 3.2, false},
            {"Bcos_dtheta", ";Bcos_dtheta;", kNBins, -1.0, 1.0, false},
            {"Btktkmass", ";Btktkmass;", kNBins, 0.0, 2.0, false},
            {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.0, 10.0, false},
            {"BtktkvProb", ";BtktkvProb;", kNBins, 0.0, 1.0, false}
        };
    } else if (treeName == "ntKstar") {
        vars = {
            {"Bpt", ";p_{T} [GeV/c];", kNBins, 7.5, 50.0, false},
            {"By", ";|y|;", kNBins, 0.0, 2.4, true},
            {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, 0.5, false},
            {"Btrk2dR", ";Btrk2dR;", kNBins, 0.0, 0.5, false},
            {"BtrkPtimb", ";BtrkPtimb;", kNBins, 0.0, 1.0, false},
            {"Btktkpt", ";Btktkpt;", kNBins, 0.0, 10.0, false},
            {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.9, 3.3, false},
            {"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 15.0, false},
            {"BsvpvDistance_2D", ";BsvpvDistance_2D;", kNBins, 0.0, 0.25, false},
            {"BQvalue", ";BQvalue;", kNBins, 0.0, 0.6, false},
            {"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -5.0, 5.0, false},
            {"Bnorm_trk2Dxy", ";Bnorm_trk2Dxy;", kNBins, -5.0, 5.0, false},
            {"Balpha", ";Balpha;", kNBins, 0.0, 3.2, false},
            {"Bcos_dtheta", ";Bcos_dtheta;", kNBins, -1.0, 1.0, false},
            {"Btktkmass", ";Btktkmass;", kNBins, 0.0, 2.0, false},
            {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.0, 5.0, false},
            {"Btrk2Pt", ";Btrk2Pt;", kNBins, 0.0, 5.0, false},
            {"BtktkvProb", ";BtktkvProb;", kNBins, 0.0, 1.0, false}
        };
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
            if (v.expr == "Btktkmass") { v.xmin = 0.4; v.xmax = 0.65; }
            if (v.expr == "Bmass") { v.xmin = 3.6; v.xmax = 4.0; }
        }
    }
    return vars;
}

static TH1D* makeWeightHist(const TH1D* hData, const TH1D* hMC, const TString& name)
{
    TH1D* hWeight = (TH1D*)hData->Clone(name);
    hWeight->Reset("ICES");
    hWeight->SetTitle(hData->GetTitle());
    hWeight->GetXaxis()->SetTitle(hData->GetXaxis()->GetTitle());
    hWeight->GetYaxis()->SetTitle("Data / MC");
    for (int i = 1; i <= hWeight->GetNbinsX(); ++i) {
        const double data = hData->GetBinContent(i);
        const double mc = hMC->GetBinContent(i);
        const double dataErr = hData->GetBinError(i);
        const double mcErr = hMC->GetBinError(i);
        double ratio = 1.0;
        double ratioErr = 0.0;
        if (mc != 0.0) {
            ratio = data / mc;
            const double relData = (data != 0.0) ? dataErr / std::abs(data) : 0.0;
            const double relMC = mcErr / std::abs(mc);
            ratioErr = std::abs(ratio) * sqrt(relData * relData + relMC * relMC);
        }
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
