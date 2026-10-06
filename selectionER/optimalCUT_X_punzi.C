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
#include "TPaveText.h"
#include "TObjString.h"
#include "TParameter.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

// The cases and their mass windows, pT bins and Punzi settings, written from ML_selection/settings.py.
#include "ML_selection/methods/ML_tmva/ml_settings.h"

// Punzi optimization of the classifier threshold for one case ("<tree>/<year>": ntmix/PbPb23,
// ntKp/ppRef, ...), one classifier, and the FOM its hyperparameters were selected with (auc, sigeff,
// pauc, punzi; see ML_selection/scripts/fresh_ML_scan.sh).
//
// Inputs are the scored samples of the classifier, read from the folder of the case
//   scored_samples/<tree>/<year>/flat_<tree>_<year>_<classifier>_<fom>_scored_{DATA,<signal kind>}.root
// Every scored file carries a metadata/ directory written by the apply step: tree name,
// score branch and range, MC weight, and the pre-ML cut the classifier was trained with.
// That pre-ML cut is the pre-cut of the optimization.
//
// Outputs go to optimalCUT/<tree>/<year>/<classifier>_<fom>/:
//   punzi_scan_<bin>.pdf   FOM versus threshold, one per pT bin
//   punzi_summary.tex      table of the optimal thresholds
//   punzi_summary.root     metadata/ and the punziResults tree, one entry per bin,
//                          including the full selection string for downstream code
//
// to run, from selectionER/:
//   root -l -b -q 'optimalCUT_X_punzi.C("ntKp/PbPb23","xgb","punzi")'
// classifier is nn, xgb or tmva, as in fresh_ML_scan.sh and train_best_ML.sh. The Punzi window and
// sidebands around the peak, the pT bins and the threshold grid are those of the case in
// ML_selection/settings.py.
//
// Samples scored outside this framework have no metadata/ directory: set MANUAL_INPUT
// below and give the files, trees, pre-cut, and score by hand. The classifier and fom
// arguments then only name the output folder, e.g.
//   root -l -b -q 'optimalCUT_X_punzi.C("ntmix/PbPb23","external","manual")'

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

// Thresholds scanned on the classifier score: min + k * (max - min) / nSteps for k = 0 .. nSteps - 1.
struct ScoreGrid {
    double min;
    double max;
    int nSteps;
};

// The inputs shared by every pT bin.
struct PunziInputs {
    TString system;
    MLSettings::Case settings;
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
    const double sidebandToSignalWindow = in.settings.punziHalfWidth / in.settings.punziSidebandWidth * 2 / in.settings.punziSides.size();

    TGraph graph;
    for (int k = 0; k < in.grid.nSteps; ++k) {
        const double sigEff = signal[k] / signal[0];
        if (!(sigEff > 0.0)) continue;

        const double threshold = in.grid.min + k * (in.grid.max - in.grid.min) / in.grid.nSteps;
        const double bkg = sideband[k] * sidebandToSignalWindow;
        const double sqrtB = TMath::Sqrt(bkg);
        const double fom = sigEff / (b * b / 2.0 + a * sqrtB + b / 2.0 * TMath::Sqrt(b * b + 4.0 * a * sqrtB + 4.0 * bkg));
        if (!std::isfinite(fom) || sideband[k] < MLSettings::PUNZI_MIN_BACKGROUND) continue;

        graph.SetPoint(graph.GetN(), threshold, fom);
        if (fom > result.bestRatio) {
            result.bestRatio = fom;
            result.bestThreshold = threshold;
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
    label.DrawLatex(0.16, 0.76, Form("Best threshold = %.3f", result.bestThreshold));

    // The sample name in a box at the top right corner, as on every figure of the ML selection.
    TPaveText sampleBox(0.85, 0.935, 0.995, 0.995, "NDC");
    sampleBox.SetFillColor(kWhite);
    sampleBox.SetLineColor(kGray + 1);
    sampleBox.SetBorderSize(1);
    sampleBox.SetTextFont(62);
    sampleBox.SetTextSize(0.035);
    sampleBox.AddText(in.system);
    sampleBox.Draw();

    c.SaveAs(Form("%s/punzi_scan_%s.pdf", outDir.Data(), bin.fileTag.Data()));
    std::cout << Form("%s %s %s: best threshold = %.3f, FOM = %.6f, signal eff. = %.4f, bkg = %.1f",
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
    put("sample", in.settings.tree + "/" + in.settings.year);
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
    put("signal_window", Form("|Bmass - %.5f| < %.4f", in.settings.peak, in.settings.punziHalfWidth));
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

// a is the discovery Z-score; b is the power/exclusion Z-score (settings.py PUNZI_A and PUNZI_B:
// 5-sigma discovery with the approximate 2-sigma exclusion convention).
void optimalCUT_X_punzi(TString caseName = "ntmix/PbPb23", TString classifier = "xgb", TString fom = "punzi",
                        double a = MLSettings::PUNZI_A, double b = MLSettings::PUNZI_B)
{
    gStyle->SetOptStat(0);

    const TString macroDir = gSystem->DirName(__FILE__);

    PunziInputs in;
    in.settings = MLSettings::CASES.at(caseName);
    const MLSettings::Case &c = in.settings;
    in.system = c.tree + " " + c.year;
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
        const TString stem = Form("%s/scored_samples/%s/%s/flat_%s_%s_%s_%s_scored_", macroDir.Data(), c.tree.Data(),
                                  c.year.Data(), c.tree.Data(), c.year.Data(), classifier.Data(), fom.Data());
        in.dataPath = stem + c.inputs.at("data").kind + ".root";
        in.mcPath = stem + c.inputs.at("signal").kind + ".root";
        fileData = TFile::Open(in.dataPath);
        fileX = TFile::Open(in.mcPath);
        in.data = fileData->Get<TTree>(fileData->Get<TObjString>("metadata/tree")->GetString());
        in.mc = fileX->Get<TTree>(fileX->Get<TObjString>("metadata/tree")->GetString());
        in.model = fileData->Get<TObjString>("metadata/model")->GetString();
        in.scoreBranch = fileData->Get<TObjString>("metadata/score_branch")->GetString();
        in.mcWeight = fileX->Get<TObjString>("metadata/mc_weight")->GetString();
        in.classifierPreCut = fileData->Get<TObjString>("metadata/pre_cut")->GetString();
        in.grid.min = fileData->Get<TObjString>("metadata/score_min")->GetString().Atof();
        in.grid.max = fileData->Get<TObjString>("metadata/score_max")->GetString().Atof();
    }
    in.grid.nSteps = TMath::Nint((in.grid.max - in.grid.min) / c.thresholdStep);
    std::cout << "Reading " << in.system << " data sample: " << in.dataPath << std::endl;
    std::cout << "Reading " << c.inputs.at("signal").label << " sample: " << in.mcPath << std::endl;

    const double mass = c.peak;
    const double start = c.punziSidebandStart;
    const double width = c.punziSidebandWidth;
    for (int side : c.punziSides)
        in.sidebandCut += Form("%s(Bmass > %.6f && Bmass < %.6f)", in.sidebandCut.IsNull() ? "" : " || ",
                               side < 0 ? mass - start - width : mass + start, side < 0 ? mass - start : mass + start + width);

    const std::vector<double>& ptBins = c.ptBins;
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

    const TString outDir = Form("%s/optimalCUT/%s/%s/%s_%s", macroDir.Data(), c.tree.Data(), c.year.Data(), classifier.Data(), fom.Data());
    gSystem->mkdir(outDir, true);

    std::cout << "Classifier: " << classifier << ", model " << in.model << std::endl;
    std::cout << Form("Using Punzi FOM = signal_eff / S_min(B), with a = %.1f and b = %.1f", a, b) << std::endl;
    std::cout << "Signal efficiency is weighted with " << in.mcWeight << std::endl;
    std::cout << Form("Scanning %s from %.2f to %.2f in steps of %g", in.scoreBranch.Data(), in.grid.min,
                      in.grid.max, c.thresholdStep) << std::endl;
    std::cout << Form("Using peak mass = %.5f GeV, signal window = +-%.1f MeV, sidebands: %s",
                      mass, 1000 * c.punziHalfWidth, in.sidebandCut.Data()) << std::endl;
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
