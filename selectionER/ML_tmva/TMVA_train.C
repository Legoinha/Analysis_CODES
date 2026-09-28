#include "TMVA_config.h"
#include "TMVA_common.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TCut.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMarker.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TTree.h"
#include "TTreeFormula.h"

#include "TMVA/Config.h"
#include "TMVA/DataLoader.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/DecisionTreeNode.h"
#include "TMVA/Factory.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

namespace {

struct TrainingHistory {
   std::vector<double> rounds;
   std::vector<double> trainAUC;
   std::vector<double> testAUC;
};

void LoadBestParameters(const TString &outputDir,
                        int &nTrees,
                        int &maxDepth,
                        double &shrinkage,
                        double &baggedSampleFraction,
                        double &minNodeSize,
                        int &nCuts)
{
   // ML_OPTIMIZATION_DIR (feature_study.py refits) or the optimization folder of the run.
   const TString directoryPath =
      gSystem->Getenv("ML_OPTIMIZATION_DIR") ? TString(gSystem->Getenv("ML_OPTIMIZATION_DIR")) : outputDir + "/optimization";
   TSystemDirectory directory("optimization", directoryPath);
   TIter next(directory.GetListOfFiles());
   double bestObjective = -1.0e30;
   while (auto *systemFile = dynamic_cast<TSystemFile *>(next())) {
      const TString name = systemFile->GetName();
      if (!name.BeginsWith("summary_job_") || !name.EndsWith(".root")) continue;
      TFile summary(directoryPath + "/" + name, "READ");
      TTree *trial = summary.Get<TTree>("trial");
      int trialTrees = 0;
      int trialDepth = 0;
      int trialCuts = 0;
      double trialShrinkage = 0.0;
      double trialBagging = 0.0;
      double trialMinNode = 0.0;
      double trialObjective = 0.0;
      int trialPassed = 0;
      trial->SetBranchAddress("n_trees", &trialTrees);
      trial->SetBranchAddress("max_depth", &trialDepth);
      trial->SetBranchAddress("shrinkage", &trialShrinkage);
      trial->SetBranchAddress("bagged_sample_fraction", &trialBagging);
      trial->SetBranchAddress("min_node_size_percent", &trialMinNode);
      trial->SetBranchAddress("n_cuts", &trialCuts);
      trial->SetBranchAddress("objective", &trialObjective);
      trial->SetBranchAddress("passed", &trialPassed);
      // One entry per tree-count checkpoint of this scan job; only those passing the
      // overtraining test can be chosen.
      for (Long64_t entry = 0; entry < trial->GetEntries(); ++entry) {
         trial->GetEntry(entry);
         if (trialPassed && trialObjective > bestObjective) {
            bestObjective = trialObjective;
            nTrees = trialTrees;
            maxDepth = trialDepth;
            shrinkage = trialShrinkage;
            baggedSampleFraction = trialBagging;
            minNodeSize = trialMinNode;
            nCuts = trialCuts;
         }
      }
   }
   std::cout << "Best optimization objective: " << bestObjective << std::endl;
   std::cout << "Best parameters: NTrees=" << nTrees << " MaxDepth=" << maxDepth << " Shrinkage=" << shrinkage
             << " BaggedSampleFraction=" << baggedSampleFraction << " MinNodeSize=" << minNodeSize
             << "% nCuts=" << nCuts << std::endl;
}

TrainingHistory BuildTrainingHistory(const TMVA::MethodBDT &method,
                                     const std::vector<ScoredEvent> &training,
                                     const std::vector<ScoredEvent> &testing)
{
   // All training and validation events: a subsample of the pThat-weighted signal fluctuates too
   // much, so its AUC gap would not be the one the overtraining test measures.
   std::vector<ScoredEvent> trainSample = training;
   std::vector<ScoredEvent> testSample = testing;
   for (auto &event : trainSample) event.score = 0.0;
   for (auto &event : testSample) event.score = 0.0;

   TrainingHistory history;
   const auto &forest = method.GetForest();
   const std::size_t checkpointStep = std::max<std::size_t>(1, forest.size() / 100);
   for (std::size_t i = 0; i < forest.size(); ++i) {
      for (auto &event : trainSample) event.score += EvaluateTree(forest[i], event.features);
      for (auto &event : testSample) event.score += EvaluateTree(forest[i], event.features);
      const bool checkpoint = i == 0 || (i + 1) % checkpointStep == 0 || i + 1 == forest.size();
      if (!checkpoint) continue;
      history.rounds.push_back(i + 1);
      history.trainAUC.push_back(BuildROC(trainSample).second);
      history.testAUC.push_back(BuildROC(testSample).second);
   }
   return history;
}

// Permutation importance: the validation AUC drop when one input is shuffled among the events,
// in percent of the sum over the inputs (a drop below zero counts as zero).
std::vector<double> PermutationImportance(const TMVA::MethodBDT &method, const std::vector<ScoredEvent> &validation)
{
   const auto &forest = method.GetForest();
   auto auc = [&forest](std::vector<ScoredEvent> events) {
      for (auto &event : events) {
         event.score = 0.0;
         for (const auto *tree : forest) event.score += EvaluateTree(tree, event.features);
      }
      return BuildROC(events).second;
   };
   const double reference = auc(validation);
   TRandom3 random(42);
   std::vector<double> importance(validation.front().features.size());
   for (std::size_t j = 0; j < importance.size(); ++j) {
      std::vector<ScoredEvent> shuffled = validation;
      for (std::size_t i = shuffled.size() - 1; i > 0; --i)
         std::swap(shuffled[i].features[j], shuffled[random.Integer(i + 1)].features[j]);
      importance[j] = std::max(0.0, reference - auc(shuffled));
   }
   const double total = std::accumulate(importance.begin(), importance.end(), 0.0);
   for (double &value : importance) value = 100.0 * value / total;
   return importance;
}

struct SkimScores {
   std::vector<double> score;
   std::vector<double> weight;   // pThatreweight for MC, 1 for DATA
   std::vector<double> bmass;
};

// Every entry of one skim scored with the trained model (the skims already carry the pre-cut).
SkimScores ScoreSkim(const TString &sample, const TString &kind, const TString &treeName, const TString &modelPath)
{
   const std::vector<TString> features = MLTMVA::Features();
   std::unique_ptr<TFile> input(TFile::Open(MLTMVA::SkimPath(sample, kind), "READ"));
   TTree *tree = input->Get<TTree>(treeName);
   std::vector<float> values(features.size());
   float bmass = 0.f;
   float entryNumber = 0.f;
   TMVA::Reader reader("!Color:Silent");
   for (std::size_t i = 0; i < features.size(); ++i) {
      reader.AddVariable(features[i], &values[i]);
      tree->SetBranchAddress(features[i], &values[i]);
   }
   reader.AddSpectator("Bmass", &bmass);
   reader.AddSpectator("skim_entry := Entry$", &entryNumber);
   tree->SetBranchAddress("Bmass", &bmass);
   reader.BookMVA(MLTMVA::METHOD_NAME, modelPath);
   TLeaf *weightLeaf = kind == "DATA" ? nullptr : tree->GetLeaf(MLTMVA::MC_WEIGHT_BRANCH);
   SkimScores scores;
   for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
      tree->GetEntry(entry);
      scores.score.push_back(reader.EvaluateMVA(MLTMVA::METHOD_NAME));
      scores.weight.push_back(weightLeaf ? weightLeaf->GetValue() : 1.0);
      scores.bmass.push_back(bmass);
   }
   return scores;
}

// The common scores file of ML_common/ml_fom.py (write_scores): train / validation / test
// with y (1 signal), w, score, bmass and idx (the skim entry); sculpt, the DATA never used
// for training; spectator_<key>, the MC spectators; metadata_json.
void WriteScores(const TString &path,
                 const std::vector<std::pair<TString, const std::vector<ScoredEvent> *>> &splits,
                 const SkimScores &data,
                 const std::vector<ScoredEvent> &training,
                 const std::vector<std::pair<TString, SkimScores>> &spectators,
                 const TString &metadata)
{
   TFile output(path, "RECREATE");
   for (const auto &[name, events] : splits) {
      float y = 0.f;
      double w = 0.0, score = 0.0, bmass = 0.0;
      Long64_t idx = 0;
      TTree tree(name, name);
      tree.Branch("y", &y);
      tree.Branch("w", &w);
      tree.Branch("score", &score);
      tree.Branch("bmass", &bmass);
      tree.Branch("idx", &idx);
      for (const auto &event : *events) {
         y = event.classID == 0 ? 1.f : 0.f;
         w = event.weight;
         score = event.score;
         bmass = event.bmass;
         idx = event.entry;
         tree.Fill();
      }
      tree.Write();
   }
   std::vector<bool> trained(data.score.size(), false);
   for (const auto &event : training)
      if (event.classID == 1) trained[event.entry] = true;
   double score = 0.0, bmass = 0.0, w = 0.0;
   TTree sculpt("sculpt", "DATA never used for training");
   sculpt.Branch("bmass", &bmass);
   sculpt.Branch("score", &score);
   for (std::size_t i = 0; i < data.score.size(); ++i) {
      if (trained[i]) continue;
      bmass = data.bmass[i];
      score = data.score[i];
      sculpt.Fill();
   }
   sculpt.Write();
   for (const auto &[key, spectator] : spectators) {
      TTree tree("spectator_" + key, key);
      tree.Branch("w", &w);
      tree.Branch("score", &score);
      for (std::size_t i = 0; i < spectator.score.size(); ++i) {
         w = spectator.weight[i];
         score = spectator.score[i];
         tree.Fill();
      }
      tree.Write();
   }
   TObjString(metadata).Write("metadata_json");
}

void Normalize(TH1D &histogram)
{
   const double integral = histogram.Integral();
   if (integral > 0.0) histogram.Scale(1.0 / integral);
}

void FillPercentDifference(TH1D &difference, const TH1D &test, const TH1D &training)
{
   for (int bin = 1; bin <= difference.GetNbinsX(); ++bin) {
      const double reference = training.GetBinContent(bin);
      difference.SetBinContent(
         bin, reference > 0.0 ? 100.0 * (test.GetBinContent(bin) - reference) / reference : 0.0);
   }
}

void SaveOutputs(const TString &sample,
                 const TString &fom,
                 const TString &preCut,
                 const TString &outputDir,
                 const TString &methodOptions,
                 const std::vector<ScoredEvent> &training,
                 const std::vector<ScoredEvent> &testing,
                 const Evaluation &validation,
                 const Evaluation &test,
                 double windowBackground,
                 const SkimScores &spectator,
                 const std::vector<double> &importance,
                 const TrainingHistory &history,
                 const std::vector<double> &boostWeights,
                 Long64_t signalEntries,
                 Long64_t backgroundEntries,
                 int nTrees,
                 int maxDepth,
                 double shrinkage,
                 double baggedSampleFraction,
                 double minNodeSize,
                 int nCuts)
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleFont(42, "XYZ");
   gStyle->SetLabelFont(42, "XYZ");
   gStyle->SetTextFont(42);
   gStyle->SetTitleSize(0.050, "XYZ");
   gStyle->SetLabelSize(0.043, "XYZ");
   gStyle->SetLineWidth(1);
   const auto testROC = BuildROC(testing);
   double ksSignal = test.ksSignal;
   double ksBackground = test.ksBackground;
   const double threshold = validation.punziThreshold;

   TGraph rocGraph(testROC.first.size());
   for (std::size_t i = 0; i < testROC.first.size(); ++i)
      rocGraph.SetPoint(i, testROC.first[i].fpr, testROC.first[i].tpr);
   rocGraph.SetTitle("ROC Curve;Background efficiency;Signal efficiency");
   rocGraph.SetLineWidth(2);
   rocGraph.SetLineColor(kBlue + 1);
   rocGraph.SetMinimum(0.0);
   rocGraph.SetMaximum(1.0);
   TCanvas rocCanvas("rocCanvas", "ROC", 600, 600);
   rocCanvas.SetLeftMargin(0.14);
   rocCanvas.SetRightMargin(0.04);
   rocCanvas.SetBottomMargin(0.12);
   rocCanvas.SetTopMargin(0.09);
   rocCanvas.SetGrid();
   rocGraph.Draw("AL");
   rocGraph.GetXaxis()->SetLimits(0.0, 1.0);
   TLine diagonal(0.0, 0.0, 1.0, 1.0);
   diagonal.SetLineColor(kGray + 1);
   diagonal.SetLineStyle(2);
   diagonal.Draw();
   rocGraph.Draw("L SAME");
   TLegend rocLegend(0.18, 0.78, 0.55, 0.88);
   rocLegend.SetBorderSize(0);
   rocLegend.SetFillStyle(0);
   rocLegend.AddEntry(&rocGraph, Form("ROC (AUC = %.4f)", testROC.second), "l");
   rocLegend.Draw();
   rocCanvas.SaveAs(outputDir + "/roc_curve.pdf");

   constexpr int scoreBins = 20;
   TH1D signalTrain("signal_train", "Score Distributions;;Normalized entries", scoreBins, -1.0, 1.0);
   TH1D backgroundTrain("background_train", "", scoreBins, -1.0, 1.0);
   TH1D signalTest("signal_test", "", scoreBins, -1.0, 1.0);
   TH1D backgroundTest("background_test", "", scoreBins, -1.0, 1.0);
   TH1D psi2s("psi2s_spectator", "", scoreBins, -1.0, 1.0);
   signalTrain.Sumw2();
   backgroundTrain.Sumw2();
   signalTest.Sumw2();
   backgroundTest.Sumw2();
   psi2s.Sumw2();
   for (const auto &event : training)
      (event.classID == 0 ? signalTrain : backgroundTrain).Fill(event.score, event.weight);
   for (const auto &event : testing)
      (event.classID == 0 ? signalTest : backgroundTest).Fill(event.score, event.weight);
   for (std::size_t i = 0; i < spectator.score.size(); ++i) psi2s.Fill(spectator.score[i], spectator.weight[i]);
   Normalize(signalTrain);
   Normalize(backgroundTrain);
   Normalize(signalTest);
   Normalize(backgroundTest);
   Normalize(psi2s);

   const int signalColor = kOrange + 7;
   const int backgroundColor = kAzure + 2;
   signalTrain.SetLineColor(signalColor);
   signalTrain.SetLineWidth(2);
   backgroundTrain.SetLineColor(backgroundColor);
   backgroundTrain.SetLineWidth(2);
   signalTest.SetLineColor(signalColor);
   signalTest.SetFillColorAlpha(signalColor, 0.25);
   backgroundTest.SetLineColor(backgroundColor);
   backgroundTest.SetFillColorAlpha(backgroundColor, 0.25);
   psi2s.SetLineColor(kOrange - 3);
   psi2s.SetLineStyle(2);
   psi2s.SetLineWidth(2);

   TH1D signalDifference("signal_test_train_difference",
                         ";TMVA BDTG score;(test-train)/train [%]", scoreBins, -1.0, 1.0);
   TH1D backgroundDifference("background_test_train_difference", "", scoreBins, -1.0, 1.0);
   FillPercentDifference(signalDifference, signalTest, signalTrain);
   FillPercentDifference(backgroundDifference, backgroundTest, backgroundTrain);
   signalDifference.SetLineColor(signalColor);
   signalDifference.SetMarkerColor(signalColor);
   signalDifference.SetMarkerStyle(20);
   signalDifference.SetMarkerSize(0.65);
   backgroundDifference.SetLineColor(backgroundColor);
   backgroundDifference.SetMarkerColor(backgroundColor);
   backgroundDifference.SetMarkerStyle(20);
   backgroundDifference.SetMarkerSize(0.65);

   TCanvas scoreCanvas("scoreCanvas", "Scores", 700, 700);
   TPad scoreTop("scoreTop", "score distributions", 0.0, 0.30, 1.0, 1.0);
   TPad scoreBottom("scoreBottom", "test minus train", 0.0, 0.0, 1.0, 0.30);
   scoreTop.SetLeftMargin(0.13);
   scoreTop.SetRightMargin(0.04);
   scoreTop.SetTopMargin(0.10);
   scoreTop.SetBottomMargin(0.02);
   scoreBottom.SetLeftMargin(0.13);
   scoreBottom.SetRightMargin(0.04);
   scoreBottom.SetTopMargin(0.03);
   scoreBottom.SetBottomMargin(0.34);
   scoreTop.SetGrid();
   scoreBottom.SetGrid();
   scoreTop.Draw();
   scoreBottom.Draw();

   scoreTop.cd();
   signalTrain.SetMaximum(1.25 * std::max({signalTrain.GetMaximum(), backgroundTrain.GetMaximum(),
                                          signalTest.GetMaximum(), backgroundTest.GetMaximum(), psi2s.GetMaximum()}));
   signalTrain.SetMinimum(0.0);
   signalTrain.GetXaxis()->SetLabelSize(0.0);
   signalTrain.Draw("HIST");
   signalTest.Draw("HIST SAME");
   backgroundTest.Draw("HIST SAME");
   backgroundTrain.Draw("HIST SAME");
   signalTrain.Draw("HIST SAME");
   psi2s.Draw("HIST SAME");
   TLegend scoreLegend(0.15, 0.60, 0.48, 0.88);
   scoreLegend.SetBorderSize(0);
   scoreLegend.SetFillStyle(0);
   scoreLegend.AddEntry(&signalTrain, "Signal train", "l");
   scoreLegend.AddEntry(&backgroundTrain, "Background train", "l");
   scoreLegend.AddEntry(&signalTest, "Signal test", "f");
   scoreLegend.AddEntry(&backgroundTest, "Background test", "f");
   scoreLegend.AddEntry(&psi2s, "#psi(2S) spectator", "l");
   scoreLegend.Draw();
   TLegend ksLegend(0.55, 0.68, 0.88, 0.88);
   ksLegend.SetBorderSize(0);
   ksLegend.SetFillStyle(0);
   ksLegend.SetHeader("Weighted KS", "C");
   ksLegend.AddEntry(static_cast<TObject *>(nullptr), Form("Signal: %.4f", ksSignal), "");
   ksLegend.AddEntry(static_cast<TObject *>(nullptr), Form("Background: %.4f", ksBackground), "");
   ksLegend.Draw();

   scoreBottom.cd();
   double maximumDifference = 1.0;
   for (int bin = 1; bin <= scoreBins; ++bin)
      maximumDifference = std::max(
         maximumDifference,
         std::max(std::abs(signalDifference.GetBinContent(bin)),
                  std::abs(backgroundDifference.GetBinContent(bin))));
   maximumDifference = std::min(250.0, 1.2 * maximumDifference);
   signalDifference.SetMinimum(-maximumDifference);
   signalDifference.SetMaximum(maximumDifference);
   signalDifference.GetXaxis()->SetTitleSize(0.12);
   signalDifference.GetXaxis()->SetLabelSize(0.10);
   signalDifference.GetYaxis()->SetTitleSize(0.10);
   signalDifference.GetYaxis()->SetLabelSize(0.08);
   signalDifference.GetYaxis()->SetTitleOffset(0.55);
   signalDifference.GetYaxis()->SetNdivisions(505);
   signalDifference.Draw("LP");
   backgroundDifference.Draw("LP SAME");
   TLine differenceZero(-1.0, 0.0, 1.0, 0.0);
   differenceZero.SetLineColor(kGray + 1);
   differenceZero.Draw();
   signalDifference.Draw("LP SAME");
   backgroundDifference.Draw("LP SAME");
   TLegend differenceLegend(0.68, 0.70, 0.94, 0.94);
   differenceLegend.SetBorderSize(0);
   differenceLegend.SetFillStyle(0);
   differenceLegend.AddEntry(&signalDifference, "Signal", "lp");
   differenceLegend.AddEntry(&backgroundDifference, "Background", "lp");
   differenceLegend.Draw();
   scoreCanvas.SaveAs(outputDir + "/score_distributions.pdf");

   TH2D confusion("confusion_matrix", Form("Confusion Matrix (thr = %.3f)", threshold),
                  2, 0, 2, 2, 0, 2);
   confusion.GetXaxis()->SetBinLabel(1, "Pred. bkg");
   confusion.GetXaxis()->SetBinLabel(2, "Pred. sig");
   confusion.GetYaxis()->SetBinLabel(1, "True sig");
   confusion.GetYaxis()->SetBinLabel(2, "True bkg");
   for (const auto &event : testing) {
      const bool predictedSignal = event.score >= threshold;
      const int predictedBin = predictedSignal ? 2 : 1;
      const int truthBin = event.classID == 0 ? 1 : 2;
      confusion.Fill(predictedBin - 0.5, truthBin - 0.5, event.weight);
   }

   double paletteStops[] = {0.0, 1.0};
   double paletteRed[] = {0.97, 0.03};
   double paletteGreen[] = {0.98, 0.25};
   double paletteBlue[] = {1.00, 0.60};
   TColor::CreateGradientColorTable(2, paletteStops, paletteRed, paletteGreen, paletteBlue, 100);
   gStyle->SetNumberContours(100);
   gStyle->SetPaintTextFormat(".3f");
   confusion.SetMarkerSize(1.35);
   confusion.GetXaxis()->SetLabelSize(0.052);
   confusion.GetYaxis()->SetLabelSize(0.052);
   TCanvas confusionCanvas("confusionCanvas", "Confusion", 500, 400);
   confusionCanvas.SetLeftMargin(0.18);
   confusionCanvas.SetRightMargin(0.14);
   confusionCanvas.SetBottomMargin(0.14);
   confusionCanvas.SetTopMargin(0.12);
   confusion.Draw("COLZ");
   TLatex confusionText;
   confusionText.SetTextAlign(22);
   confusionText.SetTextSize(0.040);
   confusionText.SetTextColor(kBlack);
   for (int xBin = 1; xBin <= 2; ++xBin)
      for (int yBin = 1; yBin <= 2; ++yBin)
         confusionText.DrawLatex(xBin - 0.5, yBin - 0.5,
                                 Form("%.3f", confusion.GetBinContent(xBin, yBin)));
   confusionCanvas.SaveAs(outputDir + "/confusion_matrix.pdf");

   const std::vector<TString> features = MLTMVA::Features();
   TH1D importanceHistogram(
      "feature_importance", "Permutation importance (validation AUC drop);Variable;Importance [%]",
      features.size(), 0, features.size());
   for (std::size_t i = 0; i < features.size(); ++i) {
      importanceHistogram.GetXaxis()->SetBinLabel(i + 1, features[i]);
      importanceHistogram.SetBinContent(i + 1, importance[i]);
   }

   std::vector<std::size_t> importanceOrder(importance.size());
   std::iota(importanceOrder.begin(), importanceOrder.end(), 0);
   std::sort(importanceOrder.begin(), importanceOrder.end(),
             [&importance](std::size_t a, std::size_t b) { return importance[a] > importance[b]; });
   std::vector<double> cumulativeImportance(importance.size(), 0.0);
   double runningImportance = 0.0;
   for (std::size_t i = 0; i < importanceOrder.size(); ++i) {
      runningImportance += importance[importanceOrder[i]];
      cumulativeImportance[i] = runningImportance;
   }

   TCanvas importanceCanvas("importanceCanvas", "Importance", 700, 500);
   importanceCanvas.SetLeftMargin(0.20);
   importanceCanvas.SetRightMargin(0.04);
   importanceCanvas.SetBottomMargin(0.13);
   importanceCanvas.SetTopMargin(0.11);
   importanceCanvas.SetGridx();
   TH2D importanceFrame("importance_frame",
                        "Cumulative Feature Importance (permutation);Permutation importance (validation AUC drop) [%];",
                        100, 0.0, 116.0, importance.size(), 0.0, importance.size());
   std::vector<TBox> importanceBoxes;
   std::vector<TMarker> importanceMarkers;
   importanceBoxes.reserve(importance.size());
   importanceMarkers.reserve(importance.size());
   for (std::size_t i = 0; i < importanceOrder.size(); ++i) {
      const int yBin = importance.size() - i;
      importanceFrame.GetYaxis()->SetBinLabel(yBin, features[importanceOrder[i]]);
   }
   importanceFrame.GetYaxis()->SetLabelSize(0.045);
   importanceFrame.Draw("AXIS");
   TLatex importanceText;
   importanceText.SetTextSize(0.025);
   importanceText.SetTextAlign(12);
   for (std::size_t i = 0; i < importanceOrder.size(); ++i) {
      const double y = importance.size() - i - 0.5;
      const double individual = importance[importanceOrder[i]];
      const double cumulative = cumulativeImportance[i];
      importanceBoxes.emplace_back(0.0, y - 0.30, cumulative, y + 0.30);
      importanceBoxes.back().SetFillColorAlpha(kGreen + 2, 0.35);
      importanceBoxes.back().SetLineColor(kGreen + 2);
      importanceBoxes.back().Draw();
      importanceMarkers.emplace_back(individual, y, 20);
      importanceMarkers.back().SetMarkerColor(kBlack);
      importanceMarkers.back().SetMarkerSize(0.9);
      importanceMarkers.back().Draw();
      importanceText.SetTextColor(kGreen + 3);
      importanceText.DrawLatex(cumulative + 1.0, y + 0.14,
                               Form("%.1f%%", cumulative));
      importanceText.SetTextColor(kBlack);
      importanceText.DrawLatex(individual + 1.0, y - 0.14,
                               Form("%.1f%%", individual));
   }
   TLegend importanceLegend(0.67, 0.75, 0.93, 0.88);
   importanceLegend.SetBorderSize(0);
   importanceLegend.SetFillStyle(0);
   if (!importanceBoxes.empty()) importanceLegend.AddEntry(&importanceBoxes.front(), "Cumulative", "f");
   if (!importanceMarkers.empty()) importanceLegend.AddEntry(&importanceMarkers.front(), "Individual", "p");
   importanceLegend.Draw();
   importanceCanvas.SaveAs(outputDir + "/feature_importance.pdf");

   TGraph trainHistory(history.rounds.size(), history.rounds.data(), history.trainAUC.data());
   TGraph testHistory(history.rounds.size(), history.rounds.data(), history.testAUC.data());
   trainHistory.SetTitle("Training History;Boosting round;AUC");
   trainHistory.SetLineColor(kBlue + 1);
   trainHistory.SetLineWidth(2);
   testHistory.SetLineColor(kOrange + 7);
   testHistory.SetLineWidth(2);
   trainHistory.SetMinimum(0.75);   // the same range as the XGB and NN training histories
   trainHistory.SetMaximum(1.0);
   TCanvas historyCanvas("historyCanvas", "Training history", 700, 500);
   historyCanvas.SetLeftMargin(0.13);
   historyCanvas.SetRightMargin(0.04);
   historyCanvas.SetBottomMargin(0.13);
   historyCanvas.SetTopMargin(0.10);
   historyCanvas.SetGrid();
   trainHistory.Draw("AL");
   testHistory.Draw("L SAME");
   TLegend historyLegend(0.48, 0.18, 0.92, 0.35);
   historyLegend.SetBorderSize(0);
   historyLegend.SetFillStyle(0);
   historyLegend.SetHeader(
      Form("AUC gap at round %.0f: %.4f", history.rounds.back(),
           history.trainAUC.back() - history.testAUC.back()));
   historyLegend.AddEntry(&trainHistory, "Train AUC", "l");
   historyLegend.AddEntry(&testHistory, "Validation AUC", "l");
   historyLegend.Draw();
   historyCanvas.SaveAs(outputDir + "/training_history.pdf");

   TFile report(outputDir + "/training_report.root", "RECREATE");
   TNamed(TString("sample"), sample).Write();
   TNamed(TString("signal_file"), MLTMVA::SkimPath(sample, "MC_X3872")).Write();
   TNamed(TString("background_file"), MLTMVA::SkimPath(sample, "DATA")).Write();
   TNamed(TString("spectator_file"), MLTMVA::SkimPath(sample, "MC_PSI2S")).Write();
   TNamed(TString("signal_tree"), MLTMVA::SIGNAL_TREE).Write();
   TNamed(TString("background_tree"), MLTMVA::DATA_TREE).Write();
   TNamed(TString("spectator_tree"), MLTMVA::SPECTATOR_TREE).Write();
   TNamed(TString("signal_cut"), preCut).Write();
   TNamed(TString("background_cut"), MLTMVA::BACKGROUND_SIDEBAND_CUT).Write();
   TNamed(TString("decorrelation_maps"), MLTMVA::DecorrelationPath(sample)).Write();
   TNamed(TString("mc_weight_branch"), MLTMVA::MC_WEIGHT_BRANCH).Write();
   TNamed(TString("method"), MLTMVA::METHOD_NAME).Write();
   TNamed(TString("method_options"), methodOptions).Write();
   TNamed(TString("host"), TString(gSystem->HostName())).Write();
   TNamed(TString("root_version"), TString(gROOT->GetVersion())).Write();
   TNamed(TString("features"), MLTMVA::FeatureList(features)).Write();

   TNamed(TString("fom"), fom).Write();
   TNamed(TString("window_background"), TString::Format("%.6g", windowBackground)).Write();
   // validation_*: model selection (train vs validation); test_*: final check (train vs test)
   // at the validation working point.
   double trainSignalWeight = WeightSum(training, 0);
   double trainBackgroundWeight = WeightSum(training, 1);
   double testSignalWeight = WeightSum(testing, 0);
   double testBackgroundWeight = WeightSum(testing, 1);
   TTree metrics("metrics", "Validation and test metrics");
   Evaluation validationCopy = validation;
   Evaluation testCopy = test;
   Int_t validationPassed = validation.passed;
   Int_t testPassed = test.passed;
   for (auto &[prefix, evaluation] : std::vector<std::pair<TString, Evaluation *>>{{"validation", &validationCopy},
                                                                              {"test", &testCopy}}) {
      metrics.Branch(prefix + "_auc", &evaluation->auc);
      metrics.Branch(prefix + "_train_auc", &evaluation->trainAUC);
      metrics.Branch(prefix + "_auc_gap", &evaluation->aucGap);
      metrics.Branch(prefix + "_sigeff", &evaluation->sigeff);
      metrics.Branch(prefix + "_pauc", &evaluation->pauc);
      metrics.Branch(prefix + "_punzi", &evaluation->punzi);
      metrics.Branch(prefix + "_punzi_threshold", &evaluation->punziThreshold);
      metrics.Branch(prefix + "_working_point", &evaluation->workingPoint);
      metrics.Branch(prefix + "_ks_signal", &evaluation->ksSignal);
      metrics.Branch(prefix + "_ks_signal_pvalue", &evaluation->ksSignalPValue);
      metrics.Branch(prefix + "_ks_background", &evaluation->ksBackground);
      metrics.Branch(prefix + "_ks_background_pvalue", &evaluation->ksBackgroundPValue);
      metrics.Branch(prefix + "_signal_eff_train", &evaluation->signalEffTrain);
      metrics.Branch(prefix + "_signal_eff_train_error", &evaluation->signalEffTrainError);
      metrics.Branch(prefix + "_signal_eff", &evaluation->signalEff);
      metrics.Branch(prefix + "_signal_eff_error", &evaluation->signalEffError);
      metrics.Branch(prefix + "_background_eff_train", &evaluation->backgroundEffTrain);
      metrics.Branch(prefix + "_background_eff_train_error", &evaluation->backgroundEffTrainError);
      metrics.Branch(prefix + "_background_eff", &evaluation->backgroundEff);
      metrics.Branch(prefix + "_background_eff_error", &evaluation->backgroundEffError);
   }
   metrics.Branch("validation_passed", &validationPassed);
   metrics.Branch("test_passed", &testPassed);
   metrics.Branch("train_signal_weight_sum", &trainSignalWeight);
   metrics.Branch("train_background_weight_sum", &trainBackgroundWeight);
   metrics.Branch("test_signal_weight_sum", &testSignalWeight);
   metrics.Branch("test_background_weight_sum", &testBackgroundWeight);
   metrics.Fill();
   metrics.Write();

   Long64_t spectatorEntries = spectator.score.size();
   Long64_t trainSignalEntries = std::count_if(training.begin(), training.end(),
                                               [](const auto &event) { return event.classID == 0; });
   Long64_t trainBackgroundEntries = training.size() - trainSignalEntries;
   Long64_t testSignalEntries = std::count_if(testing.begin(), testing.end(),
                                              [](const auto &event) { return event.classID == 0; });
   Long64_t testBackgroundEntries = testing.size() - testSignalEntries;
   TTree sampleStatistics("sample_statistics", "Selected sample statistics");
   sampleStatistics.Branch("selected_signal_entries", &signalEntries);
   sampleStatistics.Branch("selected_background_entries", &backgroundEntries);
   sampleStatistics.Branch("train_signal_entries", &trainSignalEntries);
   sampleStatistics.Branch("train_background_entries", &trainBackgroundEntries);
   sampleStatistics.Branch("test_signal_entries", &testSignalEntries);
   sampleStatistics.Branch("test_background_entries", &testBackgroundEntries);
   sampleStatistics.Branch("spectator_entries", &spectatorEntries);
   sampleStatistics.Fill();
   sampleStatistics.Write();

   TTree parameters("training_parameters", "TMVA BDTG parameters");
   parameters.Branch("n_trees", &nTrees);
   parameters.Branch("max_depth", &maxDepth);
   parameters.Branch("shrinkage", &shrinkage);
   parameters.Branch("bagged_sample_fraction", &baggedSampleFraction);
   parameters.Branch("min_node_size_percent", &minNodeSize);
   parameters.Branch("n_cuts", &nCuts);
   parameters.Fill();
   parameters.Write();

   double falsePositiveRate = 0.0;
   double truePositiveRate = 0.0;
   TTree rocTree("roc_curve", "Weighted ROC curve");
   rocTree.Branch("false_positive_rate", &falsePositiveRate);
   rocTree.Branch("true_positive_rate", &truePositiveRate);
   for (const auto &point : testROC.first) {
      falsePositiveRate = point.fpr;
      truePositiveRate = point.tpr;
      rocTree.Fill();
   }
   rocTree.Write();

   double historyRound = 0.0;
   double historyTrainAUC = 0.0;
   double historyTestAUC = 0.0;
   TTree historyTree("training_history", "Staged train and test AUC");
   historyTree.Branch("boosting_round", &historyRound);
   historyTree.Branch("train_auc", &historyTrainAUC);
   historyTree.Branch("validation_auc", &historyTestAUC);
   for (std::size_t i = 0; i < history.rounds.size(); ++i) {
      historyRound = history.rounds[i];
      historyTrainAUC = history.trainAUC[i];
      historyTestAUC = history.testAUC[i];
      historyTree.Fill();
   }
   historyTree.Write();

   int treeIndex = 0;
   double boostWeight = 0.0;
   TTree boostWeightTree("boost_weights", "TMVA boost weight for each tree");
   boostWeightTree.Branch("tree", &treeIndex);
   boostWeightTree.Branch("boost_weight", &boostWeight);
   for (std::size_t i = 0; i < boostWeights.size(); ++i) {
      treeIndex = i + 1;
      boostWeight = boostWeights[i];
      boostWeightTree.Fill();
   }
   boostWeightTree.Write();

   signalTrain.Write();
   backgroundTrain.Write();
   signalTest.Write();
   backgroundTest.Write();
   signalDifference.Write();
   backgroundDifference.Write();
   psi2s.Write();
   confusion.Write();
   importanceHistogram.Write();
   report.Close();

   PrintEvaluation("Validation", validation);
   PrintEvaluation("Test", test);
   std::cout << "Selected FOM " << fom << " on validation: " << validation.Value(fom) << std::endl;
   std::cout << "Saved model: " << MLTMVA::ModelPath(outputDir) << std::endl;
   std::cout << "Saved report: " << outputDir + "/training_report.root" << std::endl;
}

} // namespace

void TMVA_train(TString sample = "",
                TString fom = "",
                Long64_t maxSignal = 165000,
                Long64_t maxBackground = 500000,
                int nTrees = 2000,
                int maxDepth = 2,
                double shrinkage = 0.02,
                double baggedSampleFraction = 0.7,
                double minNodeSize = 1.0,
                int nCuts = 20)
{
   sample = MLTMVA::SelectedSample(sample);
   fom = MLTMVA::SelectedFom(fom);
   const auto &config = MLTMVA::GetSampleConfig(sample);
   const TString outputDir = MLTMVA::OutputDir(config, fom);
   if (gSystem->Getenv("ML_USE_OPTIMIZED") && TString(gSystem->Getenv("ML_USE_OPTIMIZED")) == "1")
      LoadBestParameters(outputDir, nTrees, maxDepth, shrinkage, baggedSampleFraction, minNodeSize, nCuts);
   gSystem->mkdir(outputDir, true);
   gSystem->mkdir(outputDir + "/dataset/weights", true);
   const std::vector<TString> features = MLTMVA::Features();
   const TString preCut = MLTMVA::SkimPreCut(sample);

   TFile *signalFile = TFile::Open(MLTMVA::SkimPath(sample, "MC_X3872"), "READ");
   TFile *backgroundFile = TFile::Open(MLTMVA::SkimPath(sample, "DATA"), "READ");
   TTree *signalTree = signalFile->Get<TTree>(MLTMVA::SIGNAL_TREE);
   TTree *backgroundTree = backgroundFile->Get<TTree>(MLTMVA::DATA_TREE);
   const Long64_t selectedSignal = signalTree->GetEntries();
   const Long64_t selectedBackground = backgroundTree->GetEntries(MLTMVA::BACKGROUND_SIDEBAND_CUT);
   const double windowBackground = WindowBackground(backgroundTree);
   const Long64_t usedSignal = std::min(maxSignal, selectedSignal);
   const Long64_t usedBackground = std::min(maxBackground, selectedBackground);
   // TMVA's test sample holds the validation and the test halves.
   const double heldOut = MLTMVA::VALIDATION_FRACTION + MLTMVA::TEST_FRACTION;
   const Long64_t testSignal = std::max<Long64_t>(2, std::llround(heldOut * usedSignal));
   const Long64_t testBackground = std::max<Long64_t>(2, std::llround(heldOut * usedBackground));
   const Long64_t trainSignal = usedSignal - testSignal;
   const Long64_t trainBackground = usedBackground - testBackground;

   std::cout << "Training sample: " << sample << std::endl;
   std::cout << "Model selection FOM: " << fom << std::endl;
   std::cout << "Signal: " << signalFile->GetName() << std::endl;
   std::cout << "Background: " << backgroundFile->GetName() << std::endl;
   std::cout << "Pre-cut of the skims: " << preCut << std::endl;
   std::cout << "Background sidebands: " << MLTMVA::BACKGROUND_SIDEBAND_CUT << std::endl;
   std::cout << "Features: " << MLTMVA::FeatureList(features) << std::endl;
   std::cout << "Selected signal/background: " << selectedSignal << "/" << selectedBackground << std::endl;
   std::cout << "Background in the Punzi signal window before the classifier cut: " << windowBackground << std::endl;

   TMVA::Tools::Instance();
   TMVA::gConfig().GetIONames().fWeightFileDirPrefix = outputDir;
   TMVA::gConfig().GetIONames().fWeightFileDir = "weights";
   TFile *tmvaOutput = TFile::Open(outputDir + "/tmva_training.root", "RECREATE");
   auto *factory = new TMVA::Factory("TMVAClassification", tmvaOutput,
                                     "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification");
   auto *loader = new TMVA::DataLoader("dataset");
   for (const auto &feature : features) loader->AddVariable(feature, 'F');
   AddSkimSpectators(*loader);
   loader->AddSignalTree(signalTree, 1.0);
   loader->AddBackgroundTree(backgroundTree, 1.0);
   loader->SetSignalWeightExpression(MLTMVA::MC_WEIGHT_BRANCH);
   loader->SetBackgroundWeightExpression("1.0");
   const TString splitOptions = Form(
      "nTrain_Signal=%lld:nTest_Signal=%lld:nTrain_Background=%lld:nTest_Background=%lld:"
      "SplitMode=Random:SplitSeed=42:NormMode=EqualNumEvents:!V",
      trainSignal, testSignal, trainBackground, testBackground);
   loader->PrepareTrainingAndTestTree(TCut(""), TCut(MLTMVA::BACKGROUND_SIDEBAND_CUT), splitOptions);

   const TString methodOptions = Form(
      "!H:!V:NTrees=%d:MinNodeSize=%.6g%%:MaxDepth=%d:BoostType=Grad:Shrinkage=%.8g:"
      "UseBaggedBoost:BaggedSampleFraction=%.8g:nCuts=%d:PruneMethod=NoPruning",
      nTrees, minNodeSize, maxDepth, shrinkage, baggedSampleFraction, nCuts);
   factory->BookMethod(loader, TMVA::Types::kBDT, MLTMVA::METHOD_NAME, methodOptions);
   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();

   auto *method = dynamic_cast<TMVA::MethodBDT *>(
      factory->GetMethod("dataset", MLTMVA::METHOD_NAME));
   const std::vector<double> boostWeights = method->GetBoostWeights();
   tmvaOutput->Write();

   std::vector<ScoredEvent> training =
      ReadScoredEvents(tmvaOutput->Get<TTree>("dataset/TrainTree"));
   auto [validation, testing] = SplitValidationTest(ReadScoredEvents(tmvaOutput->Get<TTree>("dataset/TestTree")));
   BalanceClassWeights(training);
   BalanceClassWeights(validation);
   BalanceClassWeights(testing);
   const TrainingHistory history = BuildTrainingHistory(*method, training, validation);
   const std::vector<double> importance = PermutationImportance(*method, validation);
   const Evaluation validationResult = Evaluate(training, validation, windowBackground);
   const Evaluation testResult = Evaluate(training, testing, windowBackground, validationResult.punziThreshold);

   delete factory;
   delete loader;
   tmvaOutput->Close();
   delete tmvaOutput;

   const TString modelPath = MLTMVA::ModelPath(outputDir);
   const SkimScores data = ScoreSkim(sample, "DATA", MLTMVA::DATA_TREE, modelPath);
   const std::vector<std::pair<TString, SkimScores>> spectators = {
      {"signal_nonprompt", ScoreSkim(sample, "MC_X3872_NONPROMPT", MLTMVA::SIGNAL_TREE, modelPath)},
      {"spectator", ScoreSkim(sample, "MC_PSI2S", MLTMVA::SPECTATOR_TREE, modelPath)},
      {"spectator_nonprompt", ScoreSkim(sample, "MC_PSI2S_NONPROMPT", MLTMVA::SPECTATOR_TREE, modelPath)},
   };
   SaveOutputs(sample, fom, preCut, outputDir, methodOptions, training, testing, validationResult, testResult,
               windowBackground, spectators[1].second, importance, history, boostWeights, selectedSignal,
               selectedBackground, nTrees, maxDepth, shrinkage, baggedSampleFraction, minNodeSize, nCuts);

   TString featureJson, importanceJson;
   for (std::size_t i = 0; i < features.size(); ++i) {
      featureJson += (i ? ", \"" : "\"") + features[i] + "\"";
      importanceJson += TString::Format("%s\"%s\": %.10g", i ? ", " : "", features[i].Data(), importance[i] / 100.0);
   }
   const TString metadata = TString::Format(
      "{\"method\": \"tmva\", \"sample\": \"%s\", \"fom\": \"%s\", \"features\": [%s], "
      "\"window_background\": %.10g, \"pre_cut\": \"%s\", \"random_state\": 42, "
      "\"importance\": {%s}, \"importance_method\": \"permutation (validation AUC drop)\", "
      "\"params\": {\"n_trees\": %d, \"max_depth\": %d, \"shrinkage\": %.10g, "
      "\"bagged_sample_fraction\": %.10g, \"min_node_size_percent\": %.10g, \"n_cuts\": %d}, "
      "\"method_options\": \"%s\"}",
      sample.Data(), fom.Data(), featureJson.Data(), windowBackground, preCut.Data(), importanceJson.Data(), nTrees, maxDepth, shrinkage,
      baggedSampleFraction, minNodeSize, nCuts, methodOptions.Data());
   WriteScores(outputDir + "/scores.root", {{"train", &training}, {"validation", &validation}, {"test", &testing}},
               data, training, spectators, metadata);
   std::cout << "Saved scores: " << outputDir + "/scores.root" << std::endl;
}
