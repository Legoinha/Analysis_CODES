#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {
struct Variable {
    const char *branch;
    const char *axis;
    int bins;
    double min;
    double max;
};

struct Sample {
    Sample(const TString &inputFile, const TString &legendLabel, Color_t lineColor)
        : fileName(inputFile), label(legendLabel), color(lineColor) {}

    TString fileName;
    TString label;
    Color_t color;
    std::unique_ptr<TFile> file;
    TTree *tree = nullptr;
};

// Keep these definitions synchronized with ComparePbPbMCYears.C.
const std::vector<Variable> kKinematicVars = {
    {"Bpt",     "p_{T}(X) [GeV/c]",             50, 5., 55.},
    {"Btrk1Pt", "p_{T}(#pi_{1}) [GeV/c]",       50, 0., 8.},
    {"Bmu1pt",  "p_{T}(#mu_{1}) [GeV/c]",       50, 0., 30.},
    {"Btrk1dR", "#DeltaR(#pi_{1}, J/#psi)",      50, 0., 0.5}
};

const std::vector<Variable> kQualityVars = {
    {"Btrk1PtErr",     "#sigma(p_{T}(#pi_{1})) [GeV/c]", 50, 0., 0.1},
    {"Btrk1DxyError1", "#sigma(d_{xy}(#pi_{1})) [cm]",   50, 0., 0.02},
    {"Bchi2Prob",      "vertex #chi^{2} probability",    50, 0., 1.},
    {"BtrkPtimb",      "track p_{T} imbalance",          50, 0., 1.}
};

// Feature-study pool (selectionER/fresh_ML_scan.sh FEATURE_POOL).
const std::vector<Variable> kPoolVars = {
    {"BtrkLeadPt", "leading p_{T}(#pi) [GeV/c]",    50, 0., 10.},
    {"BtrkMaxdR",  "max #DeltaR(#pi, J/#psi)",      50, 0., 0.5},
    {"BtrkSubPt",  "subleading p_{T}(#pi) [GeV/c]", 50, 0., 6.},
    {"BmuLeadPt",  "leading p_{T}(#mu) [GeV/c]",    50, 0., 40.},
    {"BmuMaxdR",   "max #DeltaR(#mu, J/#psi)",      50, 0., 0.6}
};

const std::vector<Variable> kJpsiDipionVars = {
    {"Bujpt",   "p_{T}(J/#psi) [GeV/c]", 50, 0., 50.},
    {"Bujeta",  "#eta(J/#psi)",          50, -2.5, 2.5},
    {"Bujlxy",  "L_{xy}(J/#psi) [cm]",   50, -0.1, 0.1},
    {"Btktky",  "y(#pi#pi)",             50, -2.5, 2.5}
};

// Production classifier inputs (selectionER/fresh_ML_scan.sh FEATURES).
const std::vector<Variable> kProductionVars = {
    {"Btrk1dR",   "#DeltaR(#pi_{1}, J/#psi)", 50, 0., 0.5},
    {"Btrk2Pt",   "p_{T}(#pi_{2}) [GeV/c]",   50, 0., 10.},
    {"Bmu1pt",    "p_{T}(#mu_{1}) [GeV/c]",   50, 0., 40.},
    {"Bmu2pt",    "p_{T}(#mu_{2}) [GeV/c]",   50, 0., 40.},
    {"BtrkPtimb", "track p_{T} imbalance",    50, 0., 1.},
    {"Bmu1eta",   "#eta(#mu_{1})",            50, -2.5, 2.5},
    {"Bmu2eta",   "#eta(#mu_{2})",            50, -2.5, 2.5},
    {"Bujeta",    "#eta(J/#psi)",             50, -2.5, 2.5},
    {"Btktketa",  "#eta(#pi#pi)",             50, -2.5, 2.5}
};

void DrawComparison(const std::vector<std::unique_ptr<Sample>> &samples,
                    const TString &tag, const TString &title, const TString &cut,
                    const TString &outputPath, const std::vector<Variable> &vars)
{
    const int columns = vars.size() <= 4 ? 2 : 3;
    const int rows = (vars.size() + columns - 1) / columns;
    TCanvas canvas("c_" + tag, "", 500 * columns, 400 * rows);
    canvas.Divide(columns, rows);
    gStyle->SetOptStat(0);
    gROOT->cd();

    std::vector<std::unique_ptr<TH1D>> histograms;
    histograms.reserve(vars.size() * samples.size());
    std::unique_ptr<TLegend> legend;

    for (size_t i = 0; i < vars.size(); ++i) {
        canvas.cd(i + 1);
        gPad->SetTicks();
        const auto &variable = vars[i];
        double maximum = 0.;
        std::vector<TH1D *> variableHistograms;

        for (size_t j = 0; j < samples.size(); ++j) {
            const auto &sample = *samples[j];
            auto histogram = std::make_unique<TH1D>(
                Form("h_%s_%zu_%zu", tag.Data(), i, j),
                Form(";%s;Normalized entries", variable.axis),
                variable.bins, variable.min, variable.max);
            histogram->Sumw2();

            sample.tree->Draw(Form("%s>>%s", variable.branch, histogram->GetName()), cut, "goff");
            histogram->Scale(1. / histogram->Integral());
            histogram->SetLineColor(sample.color);
            histogram->SetLineWidth(2);
            maximum = std::max(maximum, histogram->GetMaximum());
            variableHistograms.push_back(histogram.get());
            histograms.push_back(std::move(histogram));
        }

        variableHistograms.front()->SetMaximum(1.25 * maximum);
        variableHistograms.front()->Draw("HIST");
        for (size_t j = 1; j < variableHistograms.size(); ++j)
            variableHistograms[j]->Draw("HIST SAME");

        TLatex label;
        label.SetNDC();
        label.SetTextSize(0.045);
        label.DrawLatex(0.15, 0.92, title);
        label.SetTextSize(0.030);
        label.DrawLatex(0.15, 0.86, cut.IsNull() ? "No additional selection" : cut);

        if (i == 0) {
            legend = std::make_unique<TLegend>(0.52, 0.66, 0.88, 0.88);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            for (size_t j = 0; j < samples.size(); ++j)
                legend->AddEntry(variableHistograms[j], samples[j]->label, "l");
            legend->Draw();
        }
    }

    canvas.Print(outputPath);
}
} // namespace

// Run from the repository root with:
// root -l -b -q 'flatER/ComparePbPbDataYears.C()'
//
// Run 3 PbPb data years and the 2024 ppRef. Unlike the MC comparison, data are not multiplied by
// pThatreweight. Each distribution is normalized to unit area after applying the common selection
// (default: the pt10y24 region and pre-cut of selectionER/ML_selection/settings.py).
void ComparePbPbDataYears(
    TString cut = "Bpt > 10 && Bpt < 50 && abs(By) < 2.4 && BQvalue < 0.15 && Btrk2dR < 0.3")
{
    const TString base = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872";
    std::vector<std::unique_ptr<Sample>> samples;
    samples.push_back(std::make_unique<Sample>(
        base + "/PbPb23/flat_ntmix_PbPb23_DATA.root", "PbPb 2023", kRed + 1));
    samples.push_back(std::make_unique<Sample>(
        base + "/PbPb24/flat_ntmix_PbPb24_DATA.root", "PbPb 2024", kGreen + 2));
    samples.push_back(std::make_unique<Sample>(
        base + "/PbPb25/flat_ntmix_PbPb25_DATA.root", "PbPb 2025", kBlue + 1));
    samples.push_back(std::make_unique<Sample>(
        base + "/ppRef24/flat_ntmix_ppRef_DATA.root", "ppRef 2024", kBlack));

    for (auto &sample : samples) {
        sample->file.reset(TFile::Open(sample->fileName));
        sample->tree = sample->file->Get<TTree>("ntmix");
        std::cout << sample->label << ": " << sample->tree->GetEntries() << " entries in " << sample->fileName << std::endl;
    }

    const TString outputDir = "flatER/data_year_comparisons/ntmix";
    gSystem->mkdir(outputDir, true);

    const TString output = outputDir + "/PbPb23_vs_PbPb24_vs_PbPb25_vs_ppRef_DATA.pdf";
    DrawComparison(samples, "dataYears_kinematic", "Run 3 data", cut,
                   output + "(", kKinematicVars);
    DrawComparison(samples, "dataYears_quality", "Run 3 data, track quality", cut,
                   output, kQualityVars);
    DrawComparison(samples, "dataYears_jpsiDipion", "Run 3 data, J/#psi and dipion", cut,
                   output, kJpsiDipionVars);
    DrawComparison(samples, "dataYears_pool", "Run 3 data, feature-study pool", cut,
                   output, kPoolVars);
    DrawComparison(samples, "dataYears_production", "Run 3 data, production inputs", cut,
                   output + ")", kProductionVars);

    std::cout << "Saved Run 3 data-year comparison: " << output << std::endl;
}
