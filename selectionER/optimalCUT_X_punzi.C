#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include "../plotER/aux/masses.h"

// Punzi optimization of the classifier threshold for one system, one classifier, and the
// FOM its hyperparameters were selected with (auc, sigeff, pauc, punzi; see fresh_ML_scan.sh).
//
// Inputs are the scored samples of the classifier, read from its fixed folder
//   <ML folder>/scored_samples/flat_ntmix_<system>_<fom>_scored_{DATA,MC_X3872}.root
// Every scored file carries a metadata/ directory written by the apply step: tree name,
// score branch and range, MC weight, and the pre-ML cut the classifier was trained with.
// That pre-ML cut is the pre-cut of the optimization.
//
// Outputs go to optimalCUT/<classifier>/<system>_<fom>/:
//   punzi_scan_<bin>.pdf   FOM versus threshold, one per pT bin
//   punzi_summary.tex      table of the optimal thresholds
//   punzi_summary.root     metadata/ and the punziResults tree, one entry per bin,
//                          including the full selection string for downstream code
//
// to run, from selectionER/:
//   root -l -b -q 'optimalCUT_X_punzi.C("PbPb23","xgb","punzi")'
// classifier is nn, xgb or tmva, as in fresh_ML_scan.sh and train_best_ML.sh.
// The signal window and sidebands below are written by fresh_ML_scan.sh.
//
// Samples scored outside this framework have no metadata/ directory: set MANUAL_INPUT
// below and give the files, trees, pre-cut, and score by hand. The classifier and fom
// arguments then only name the output folder, e.g.
//   root -l -b -q 'optimalCUT_X_punzi.C("PbPb23","external","manual")'

static const std::map<TString, TString> CLASSIFIER_DIRS = {
    {"nn", "ML_pytorch"},
    {"xgb", "ML_xgboost"},
    {"tmva", "ML_tmva"},
};

// Manual input, used instead of the scored-sample metadata when MANUAL_INPUT is true.
static const bool MANUAL_INPUT = false;
static const TString MANUAL_DATA_FILE = "/eos/user/h/hmarques/path/to/flat_ntmix_PbPb23_scored_DATA.root";
static const TString MANUAL_MC_FILE = "/eos/user/h/hmarques/path/to/flat_ntmix_PbPb23_scored_MC_X3872.root";
static const TString MANUAL_DATA_TREE = "ntmix";
static const TString MANUAL_MC_TREE = "ntmix_X3872";
static const TString MANUAL_PRE_CUT = "((Bpt > 15) && (Bpt < 50)) && (abs(By) < 2.4) && (BQvalue < 0.15)";
static const TString MANUAL_SCORE_BRANCH = "Prediction";
static const double MANUAL_SCORE_MIN = 0.0;
static const double MANUAL_SCORE_MAX = 1.0;
static const TString MANUAL_MC_WEIGHT = "pThatreweight";
static const TString MANUAL_MODEL = "external classifier";

// Analysis pT bins. They must lie inside the pT range of the classifier pre-cut.
static const std::map<TString, std::vector<double>> PT_BINS = {
    {"ppRef", {7.5, 12.5, 17.5, 22.5, 50.0}},
    {"PbPb18", {15.0, 50.0}},
    {"PbPb23", {15.0, 50.0}},
};

static const double THRESHOLD_STEP = 0.02;

// Signal window and sidebands around the X(3872) mass, in GeV. The MC peak has a core
// sigma of about 5 MeV; the sidebands stay clear of the training sidebands, which start
// 18 MeV above and 22 MeV below the peak.
static const double SIGNAL_HALF_WIDTH = 0.005;
static const double SIDEBAND_START = 0.015;
static const double SIDEBAND_WIDTH = 0.02;
// The Punzi maximum only among cuts keeping at least this many DATA sideband candidates.
static const double PUNZI_MIN_BACKGROUND = 10.0;

struct PunziBin {
    double low;
    double high;
    bool inclusive;
    TString cut;
    TString plotLabel;
    TString fileTag;
};

struct PunziResult {
    PunziBin bin;
    TString preCut;
    TString selection;
    double bestThreshold = 0.0;
    double bestRatio = -1.0;
    double bestSignalEff = 0.0;
    double bestBkg = 0.0;
};

// Thresholds scanned on the classifier score: min + k * step for k = 0 .. nSteps - 1.
struct ScoreGrid {
    double min;
    double max;
    int nSteps;
    double threshold(int k) const { return min + k * (max - min) / nSteps; }
};

// The inputs shared by every pT bin.
struct PunziInputs {
    TString system;
    TString classifier;
    TString fom;
    TString model;
    TString dataPath;
    TString mcPath;
    TTree* data;
    TTree* mc;
    TString scoreBranch;
    TString mcWeight;
    TString classifierPreCut;
    TString sidebandCut;
    ScoreGrid grid;
};

TString metadata(TFile* file, const char* key)
{
    return file->Get<TObjString>(Form("metadata/%s", key))->GetString();
}

double punziSmin(double backgroundYield, double a, double b)
{
    const double sqrtB = TMath::Sqrt(backgroundYield);
    return b * b / 2.0
           + a * sqrtB
           + b / 2.0 * TMath::Sqrt(b * b + 4.0 * a * sqrtB + 4.0 * backgroundYield);
}

// Weighted number of candidates passing `selection` with a score above threshold k,
// for every threshold of the grid, filled in a single pass over the tree.
std::vector<double> countsAboveThresholds(TTree* tree, const TString& scoreBranch,
                                          const TString& weight, const TString& selection,
                                          const ScoreGrid& grid)
{
    TH1D histogram("punziScoreHistogram", "", grid.nSteps, grid.min, grid.max);
    tree->Project("punziScoreHistogram", scoreBranch, Form("(%s) * (%s)", weight.Data(), selection.Data()));

    std::vector<double> counts(grid.nSteps);
    for (int k = 0; k < grid.nSteps; ++k) counts[k] = histogram.Integral(k + 1, grid.nSteps + 1);
    return counts;
}

PunziResult optimizeBin(const PunziInputs& in, const PunziBin& bin, const TString& outDir, double a, double b)
{
    PunziResult result;
    result.bin = bin;
    result.preCut = Form("(%s) && (%s)", in.classifierPreCut.Data(), bin.cut.Data());

    const std::vector<double> signal = countsAboveThresholds(in.mc, in.scoreBranch, in.mcWeight, result.preCut, in.grid);
    const std::vector<double> sideband = countsAboveThresholds(
        in.data, in.scoreBranch, "1", Form("(%s) && (%s)", result.preCut.Data(), in.sidebandCut.Data()), in.grid);
    const double sidebandToSignalWindow = SIGNAL_HALF_WIDTH / SIDEBAND_WIDTH;

    TGraph graph;
    for (int k = 0; k < in.grid.nSteps; ++k) {
        const double sigEff = signal[k] / signal[0];
        if (!(sigEff > 0.0)) continue;

        const double bkg = sideband[k] * sidebandToSignalWindow;
        const double fom = sigEff / punziSmin(bkg, a, b);
        if (!std::isfinite(fom) || sideband[k] < PUNZI_MIN_BACKGROUND) continue;

        graph.SetPoint(graph.GetN(), in.grid.threshold(k), fom);
        if (fom > result.bestRatio) {
            result.bestRatio = fom;
            result.bestThreshold = in.grid.threshold(k);
            result.bestSignalEff = sigEff;
            result.bestBkg = bkg;
        }
    }
    result.selection = Form("(%s) && (%s > %.3f)", result.preCut.Data(), in.scoreBranch.Data(), result.bestThreshold);

    const double yMin = 0.0;
    const double yMax = 1.2 * result.bestRatio;

    TCanvas c(Form("c_%s", bin.fileTag.Data()), "", 800, 600);
    c.SetLeftMargin(0.14);
    TH1F* frame = c.DrawFrame(in.grid.min, yMin, in.grid.max, yMax);
    frame->SetTitle(Form(" ; %s (%s); FOM", in.scoreBranch.Data(), in.classifier.Data()));
    graph.SetMarkerStyle(20);
    graph.Draw("LP SAME");

    TLine bestLine(result.bestThreshold, yMin, result.bestThreshold, result.bestRatio);
    bestLine.SetLineStyle(2);
    bestLine.SetLineWidth(2);
    bestLine.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.035);
    label.DrawLatex(0.16, 0.86, Form("%s, %s (%s)", in.system.Data(), in.classifier.Data(), in.fom.Data()));
    label.DrawLatex(0.16, 0.81, bin.plotLabel);
    label.DrawLatex(0.16, 0.76, Form("Best threshold = %.2f", result.bestThreshold));

    c.SaveAs(Form("%s/punzi_scan_%s.pdf", outDir.Data(), bin.fileTag.Data()));
    std::cout << Form("%s %s %s: best threshold = %.2f, FOM = %.6f, signal eff. = %.4f, bkg = %.1f",
                      in.system.Data(), in.classifier.Data(), bin.plotLabel.Data(),
                      result.bestThreshold, result.bestRatio, result.bestSignalEff, result.bestBkg)
              << std::endl;

    return result;
}

void writeSummaryTable(const TString& outDir, const PunziInputs& in, const std::vector<PunziResult>& results)
{
    std::ofstream out(Form("%s/punzi_summary.tex", outDir.Data()));
    out << std::fixed << std::setprecision(4);
    out << "\\documentclass{article}\n";
    out << "\\usepackage{geometry}\n";
    out << "\\usepackage{booktabs}\n";
    out << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm}\n";
    out << "\\begin{document}\n";
    out << "\\begin{center}\n";
    out << "\\small\n";
    out << in.system << ", " << in.classifier << " (" << in.fom << ")\\\\[4pt]\n";
    out << "\\begin{tabular}{c|c|c|c|c}\n";
    out << "\\toprule\n";
    out << "Bin ($p_{T}$ [GeV/c]) & Best threshold & Punzi FOM & Signal eff. & Bkg. estimate \\\\ \\midrule\n";
    for (const auto& result : results) {
        out << Form("%.1f--%.1f%s", result.bin.low, result.bin.high, result.bin.inclusive ? " (incl.)" : "") << " & "
            << result.bestThreshold << " & "
            << result.bestRatio << " & "
            << result.bestSignalEff << " & "
            << result.bestBkg << " \\\\\n";
    }
    out << "\\bottomrule\n";
    out << "\\end{tabular}\n";
    out << "\\end{center}\n";
    out << "\\end{document}\n";
}

void writeRootSummary(const TString& outDir, const PunziInputs& in,
                      const std::vector<PunziResult>& results, double a, double b)
{
    TFile output(Form("%s/punzi_summary.root", outDir.Data()), "RECREATE");

    TDirectory* meta = output.mkdir("metadata");
    auto put = [meta](const char* key, const TString& value) {
        TObjString text(value);
        meta->WriteTObject(&text, key);
    };
    put("input", MANUAL_INPUT ? "manual (MANUAL_INPUT in optimalCUT_X_punzi.C)" : "scored-sample metadata");
    put("classifier", in.classifier);
    put("sample", in.system);
    put("fom", in.fom);
    put("model", in.model);
    put("data_file", in.dataPath);
    put("mc_file", in.mcPath);
    put("data_tree", in.data->GetName());
    put("mc_tree", in.mc->GetName());
    put("score_branch", in.scoreBranch);
    put("mc_weight", in.mcWeight);
    put("pre_cut", in.classifierPreCut);
    put("sideband_cut", in.sidebandCut);
    put("signal_window", Form("|Bmass - %.5f| < %.4f", X3872_MASS, SIGNAL_HALF_WIDTH));
    TParameter<double> aParameter("a", a);
    TParameter<double> bParameter("b", b);
    meta->WriteTObject(&aParameter);
    meta->WriteTObject(&bParameter);

    output.cd();
    TTree summary("punziResults", "Per-bin Punzi cut-optimization results");

    TString binCut;
    TString preCut;
    TString selection;
    double ptLow = 0.0;
    double ptHigh = 0.0;
    Bool_t inclusive = false;
    double optimalThreshold = 0.0;
    double punziFOM = 0.0;
    double signalEfficiency = 0.0;
    double backgroundEstimate = 0.0;

    summary.Branch("binCut", &binCut);
    summary.Branch("preCut", &preCut);
    summary.Branch("selection", &selection);
    summary.Branch("ptLow", &ptLow);
    summary.Branch("ptHigh", &ptHigh);
    summary.Branch("inclusive", &inclusive);
    summary.Branch("optimalThreshold", &optimalThreshold);
    summary.Branch("punziFOM", &punziFOM);
    summary.Branch("signalEfficiency", &signalEfficiency);
    summary.Branch("backgroundEstimate", &backgroundEstimate);

    for (const auto& result : results) {
        binCut = result.bin.cut;
        preCut = result.preCut;
        selection = result.selection;
        ptLow = result.bin.low;
        ptHigh = result.bin.high;
        inclusive = result.bin.inclusive;
        optimalThreshold = result.bestThreshold;
        punziFOM = result.bestRatio;
        signalEfficiency = result.bestSignalEff;
        backgroundEstimate = result.bestBkg;
        summary.Fill();
    }

    summary.Write();
    output.Close();
}

// a is the discovery Z-score; b is the power/exclusion Z-score.
// The defaults target 5-sigma discovery with the approximate 2-sigma exclusion convention.
void optimalCUT_X_punzi(TString system = "PbPb23", TString classifier = "xgb", TString fom = "punzi",
                        double a = 5.0, double b = 2.0)
{
    gStyle->SetOptStat(0);

    const TString macroDir = gSystem->DirName(__FILE__);

    PunziInputs in;
    in.system = system;
    in.classifier = classifier;
    in.fom = fom;
    TFile* fileData = nullptr;
    TFile* fileX = nullptr;
    if (MANUAL_INPUT) {
        std::cout << "Manual input (MANUAL_INPUT = true): scored-sample metadata is not read" << std::endl;
        in.dataPath = MANUAL_DATA_FILE;
        in.mcPath = MANUAL_MC_FILE;
        fileData = TFile::Open(in.dataPath);
        fileX = TFile::Open(in.mcPath);
        in.data = fileData->Get<TTree>(MANUAL_DATA_TREE);
        in.mc = fileX->Get<TTree>(MANUAL_MC_TREE);
        in.model = MANUAL_MODEL;
        in.scoreBranch = MANUAL_SCORE_BRANCH;
        in.mcWeight = MANUAL_MC_WEIGHT;
        in.classifierPreCut = MANUAL_PRE_CUT;
        in.grid.min = MANUAL_SCORE_MIN;
        in.grid.max = MANUAL_SCORE_MAX;
    } else {
        const TString scoredDir = Form("%s/%s/scored_samples", macroDir.Data(), CLASSIFIER_DIRS.at(classifier).Data());
        in.dataPath = Form("%s/flat_ntmix_%s_%s_scored_DATA.root", scoredDir.Data(), system.Data(), fom.Data());
        in.mcPath = Form("%s/flat_ntmix_%s_%s_scored_MC_X3872.root", scoredDir.Data(), system.Data(), fom.Data());
        fileData = TFile::Open(in.dataPath);
        fileX = TFile::Open(in.mcPath);
        in.data = fileData->Get<TTree>(metadata(fileData, "tree"));
        in.mc = fileX->Get<TTree>(metadata(fileX, "tree"));
        in.model = metadata(fileData, "model");
        in.scoreBranch = metadata(fileData, "score_branch");
        in.mcWeight = metadata(fileX, "mc_weight");
        in.classifierPreCut = metadata(fileData, "pre_cut");
        in.grid.min = metadata(fileData, "score_min").Atof();
        in.grid.max = metadata(fileData, "score_max").Atof();
    }
    in.grid.nSteps = TMath::Nint((in.grid.max - in.grid.min) / THRESHOLD_STEP);
    std::cout << "Reading " << system << " data sample: " << in.dataPath << std::endl;
    std::cout << "Reading prompt X(3872) MC sample: " << in.mcPath << std::endl;

    const double mass = X3872_MASS;
    in.sidebandCut = Form("(Bmass > %.6f && Bmass < %.6f) || (Bmass > %.6f && Bmass < %.6f)",
                          mass - SIDEBAND_START - SIDEBAND_WIDTH, mass - SIDEBAND_START,
                          mass + SIDEBAND_START, mass + SIDEBAND_START + SIDEBAND_WIDTH);

    const std::vector<double>& ptBins = PT_BINS.at(system);
    std::vector<PunziBin> bins;
    const double pMin = ptBins.front();
    const double pMax = ptBins.back();
    TString inclusiveFileTag = Form("Bpt_%.1f_%.1f_inclusive", pMin, pMax);
    inclusiveFileTag.ReplaceAll(".", "p");
    bins.push_back({pMin, pMax, true,
                    Form("Bpt > %.8f && Bpt < %.8f", pMin, pMax),
                    Form("%.1f < p_{T} [GeV/c] < %.1f (incl.)", pMin, pMax),
                    inclusiveFileTag});

    for (size_t i = 0; i + 1 < ptBins.size(); ++i) {
        const double low = ptBins[i];
        const double high = ptBins[i + 1];
        const char* lowerCut = (i == 0) ? ">" : ">=";
        const char* lowerLabel = (i == 0) ? "<" : "#leq";
        TString fileTag = Form("Bpt_%.1f_%.1f", low, high);
        fileTag.ReplaceAll(".", "p");
        bins.push_back({low, high, false,
                Form("Bpt %s %.8f && Bpt < %.8f", lowerCut, low, high),
                Form("%.1f %s p_{T} [GeV/c] < %.1f", low, lowerLabel, high),
                fileTag});
    }

    const TString outDir = Form("%s/optimalCUT/%s/%s_%s", macroDir.Data(), classifier.Data(), system.Data(), fom.Data());
    gSystem->mkdir(outDir, true);

    std::cout << "Classifier: " << classifier << ", model " << in.model << std::endl;
    std::cout << Form("Using Punzi FOM = signal_eff / S_min(B), with a = %.1f and b = %.1f", a, b) << std::endl;
    std::cout << "Signal efficiency is weighted with " << in.mcWeight << std::endl;
    std::cout << Form("Scanning %s from %.2f to %.2f in steps of %.2f", in.scoreBranch.Data(), in.grid.min,
                      in.grid.max, THRESHOLD_STEP) << std::endl;
    std::cout << Form("Using X3872 mass = %.5f GeV, signal window = +-%.1f MeV, sidebands = %.1f MeV wide starting at +-%.1f MeV",
                      X3872_MASS, 1000 * SIGNAL_HALF_WIDTH, 1000 * SIDEBAND_WIDTH, 1000 * SIDEBAND_START) << std::endl;
    std::cout << "Using classifier pre-cut: " << in.classifierPreCut << std::endl;

    std::vector<PunziResult> results;
    for (const auto& bin : bins) {
        results.push_back(optimizeBin(in, bin, outDir, a, b));
    }
    writeSummaryTable(outDir, in, results);
    writeRootSummary(outDir, in, results, a, b);

    fileData->Close();
    fileX->Close();
}
