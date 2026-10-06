#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

// PbPb18 (5.02 TeV) vs PbPb23 (5.36 TeV) MC shapes of X(3872) and Psi(2S), pThat-weighted.


namespace {
struct Variable {
    const char *branch;
    const char *axis;
    int bins;
    double min;
    double max;
};

const std::vector<Variable> kKinematicVars = {
    {"Bpt",     "p_{T}(X) [GeV/c]",             50, 5., 55.},
    {"Btrk1Pt", "p_{T}(#pi_{1}) [GeV/c]",       50, 0., 8.},
    {"Bmu1pt",  "p_{T}(#mu_{1}) [GeV/c]",       50, 0., 30.},
    {"Btrk1dR", "#DeltaR(#pi_{1}, J/#psi)",      50, 0., 0.5}
};

const std::vector<Variable> kQualityVars = {
    {"Btrk1PtErr",       "#sigma(p_{T}(#pi_{1})) [GeV/c]", 50, 0., 0.1},
    {"Btrk1DxyError1",   "#sigma(d_{xy}(#pi_{1})) [cm]",   50, 0., 0.02},
    {"Bchi2Prob",        "vertex #chi^{2} probability",    50, 0., 1.},
    {"BtrkPtimb",        "track p_{T} imbalance",          50, 0., 1.}
};

// Feature-study pool (selectionER/fresh_ML_scan.sh FEATURE_POOL).
const std::vector<Variable> kPoolVars = {
    {"BtrkLeadPt", "leading p_{T}(#pi) [GeV/c]",           50, 0., 10.},
    {"BtrkMaxdR",  "max #DeltaR(#pi, J/#psi)",             50, 0., 0.5},
    {"BtrkSubPt",  "subleading p_{T}(#pi) [GeV/c]",        50, 0., 6.},
    {"BmuLeadPt",  "leading p_{T}(#mu) [GeV/c]",           50, 0., 40.},
    {"BmuMaxdR",   "max #DeltaR(#mu, J/#psi)",             50, 0., 0.6}
};

const std::vector<Variable> kJpsiDipionVars = {
    {"Bujpt",   "p_{T}(J/#psi) [GeV/c]", 50, 0., 50.},
    {"Bujeta",  "#eta(J/#psi)",          50, -2.5, 2.5},
    {"Bujlxy",  "L_{xy}(J/#psi) [cm]",   50, -0.1, 0.1},
    {"Btktky",  "y(#pi#pi)",             50, -2.5, 2.5}
};

// Production classifier inputs (selectionER/fresh_ML_scan.sh FEATURES).
const std::vector<Variable> kProductionVars = {
    {"Btrk1dR",   "#DeltaR(#pi_{1}, J/#psi)",  50, 0., 0.5},
    {"Btrk2Pt",   "p_{T}(#pi_{2}) [GeV/c]",    50, 0., 10.},
    {"Bmu1pt",    "p_{T}(#mu_{1}) [GeV/c]",    50, 0., 40.},
    {"Bmu2pt",    "p_{T}(#mu_{2}) [GeV/c]",    50, 0., 40.},
    {"BtrkPtimb", "track p_{T} imbalance",     50, 0., 1.},
    {"Bmu1eta",   "#eta(#mu_{1})",             50, -2.5, 2.5},
    {"Bmu2eta",   "#eta(#mu_{2})",             50, -2.5, 2.5},
    {"Bujeta",    "#eta(J/#psi)",              50, -2.5, 2.5},
    {"Btktketa",  "#eta(#pi#pi)",              50, -2.5, 2.5}
};

void DrawComparison(TTree *first, TTree *second, const TString &tag, const TString &title,
                    const TString &firstLabel, const TString &secondLabel, const TString &cut,
                    const TString &outputName, const std::vector<Variable> &vars)
{
    const TString weightedCut = "pThatreweight*(" + cut + ")";

    const int columns = vars.size() <= 4 ? 2 : 3;
    const int rows = (vars.size() + columns - 1) / columns;
    TCanvas canvas("c_" + tag, "", 500 * columns, 400 * rows);
    canvas.Divide(columns, rows);
    gStyle->SetOptStat(0);
    std::vector<std::unique_ptr<TH1D>> histograms;
    std::vector<std::unique_ptr<TLegend>> legends;
    gROOT->cd();

    for (size_t i = 0; i < vars.size(); ++i) {
        canvas.cd(i + 1);
        gPad->SetTicks();
        const auto &v = vars[i];
        auto hFirst = std::make_unique<TH1D>(Form("hFirst_%s_%zu", tag.Data(), i),
                                             Form(";%s;Normalized entries", v.axis), v.bins, v.min, v.max);
        auto hSecond = std::make_unique<TH1D>(Form("hSecond_%s_%zu", tag.Data(), i),
                                              Form(";%s;Normalized entries", v.axis), v.bins, v.min, v.max);
        hFirst->Sumw2();
        hSecond->Sumw2();
        first->Draw(Form("%s>>%s", v.branch, hFirst->GetName()), weightedCut, "goff");
        second->Draw(Form("%s>>%s", v.branch, hSecond->GetName()), weightedCut, "goff");
        hFirst->Scale(1. / hFirst->Integral());
        hSecond->Scale(1. / hSecond->Integral());

        hFirst->SetLineColor(kBlue + 1);
        hSecond->SetLineColor(kRed + 1);
        hFirst->SetLineWidth(2);
        hSecond->SetLineWidth(2);
        hFirst->SetMaximum(1.25 * std::max(hFirst->GetMaximum(), hSecond->GetMaximum()));
        hFirst->Draw("HIST");
        hSecond->Draw("HIST SAME");

        TLatex label;
        label.SetNDC();
        label.SetTextSize(0.045);
        label.DrawLatex(0.15, 0.92, title);
        label.SetTextSize(0.032);
        label.DrawLatex(0.15, 0.86, cut);

        if (i == 0) {
            auto legend = std::make_unique<TLegend>(0.48, 0.68, 0.88, 0.88);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.030);
            legend->AddEntry(hFirst.get(), firstLabel, "l");
            legend->AddEntry(hSecond.get(), secondLabel, "l");
            legend->Draw();
            legends.push_back(std::move(legend));
        }
        histograms.push_back(std::move(hFirst));
        histograms.push_back(std::move(hSecond));
    }

    const TString outputDir = "flatER/reweighting_comparisons/ntmix";
    gSystem->mkdir(outputDir, true);
    canvas.Print(outputDir + "/" + outputName);
}

void ComparePbPbYears(const TString &particle, const TString &particleLabel, const TString &cut)
{
    const TString base = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872";
    const TString file18 = base + "/legacy_run2/REflated/flat_ntmix_PbPb18_MC_" + particle + ".root";
    const TString file23 = base + "/PbPb23/flat_ntmix_PbPb23_MC_" + particle + ".root";
    std::unique_ptr<TFile> run18(TFile::Open(file18));
    std::unique_ptr<TFile> run23(TFile::Open(file23));
    TTree *tree18 = run18->Get<TTree>("ntmix_" + particle);
    TTree *tree23 = run23->Get<TTree>("ntmix_" + particle);

    const TString outputName = "PbPb18_vs_PbPb23_" + particle + "_reweighted.pdf";
    DrawComparison(tree18, tree23, "PbPbYears_kinematic_" + particle,
                   particleLabel + ", PbPb MC",
                   "2018 (5.02 TeV)", "2023 (5.36 TeV)",
                   cut, outputName + "(", kKinematicVars);
    DrawComparison(tree18, tree23, "PbPbYears_quality_" + particle,
                   particleLabel + ", PbPb MC",
                   "2018 (5.02 TeV)", "2023 (5.36 TeV)",
                   cut, outputName, kQualityVars);
    DrawComparison(tree18, tree23, "PbPbYears_jpsiDipion_" + particle,
                   particleLabel + ", PbPb MC, J/#psi and dipion",
                   "2018 (5.02 TeV)", "2023 (5.36 TeV)",
                   cut, outputName, kJpsiDipionVars);
    DrawComparison(tree18, tree23, "PbPbYears_pool_" + particle,
                   particleLabel + ", PbPb MC, feature-study pool",
                   "2018 (5.02 TeV)", "2023 (5.36 TeV)",
                   cut, outputName, kPoolVars);
    DrawComparison(tree18, tree23, "PbPbYears_production_" + particle,
                   particleLabel + ", PbPb MC, production inputs",
                   "2018 (5.02 TeV)", "2023 (5.36 TeV)",
                   cut, outputName + ")", kProductionVars);
}

} // namespace

// The legacy PbPb18 flat files already contain only candidates with Bchi2Prob > 0.10.
void ComparePbPbMCYears(
    TString cut = "Bpt > 15 && Bpt < 50 && abs(By) < 2.4 && BQvalue < 0.15 && Bchi2Prob > 0.10")
{
    ComparePbPbYears("X3872", "X(3872)", cut);
    ComparePbPbYears("PSI2S", "#psi(2S)", cut);
}
