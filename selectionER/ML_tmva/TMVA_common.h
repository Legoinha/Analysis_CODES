#ifndef ML_TMVA_COMMON_H
#define ML_TMVA_COMMON_H

// Helpers shared by TMVA_train.C and TMVA_optimize.C: read TMVA's scored
// train/test trees, evaluate the BDT forest tree by tree, compute weighted
// ROC curves, and evaluate the figures of merit and the overtraining test with
// the same definitions as ML_common/ml_fom.py.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <tuple>
#include <utility>
#include <vector>

#include "TLeaf.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/DecisionTree.h"
#include "TMVA/DecisionTreeNode.h"
#include "TMVA/MethodBDT.h"

#include "TMVA_config.h"

namespace {

struct ScoredEvent {
   double score;
   double weight;
   int classID;   // TMVA: 0 signal, 1 background
   std::vector<double> features;
   double bmass;
   Long64_t entry;   // entry in the skim of its class
};

// Spectators every TMVA training books, so each event can be traced back to its skim entry.
inline void AddSkimSpectators(TMVA::DataLoader &loader)
{
   loader.AddSpectator("Bmass", 'F');
   loader.AddSpectator("skim_entry := Entry$", 'F');
}

std::vector<ScoredEvent> ReadScoredEvents(TTree *tree)
{
   std::vector<ScoredEvent> events;
   events.reserve(tree->GetEntries());
   TLeaf *scoreLeaf = tree->GetLeaf(MLTMVA::METHOD_NAME);
   TLeaf *weightLeaf = tree->GetLeaf("weight");
   TLeaf *classLeaf = tree->GetLeaf("classID");
   TLeaf *massLeaf = tree->GetLeaf("Bmass");
   TLeaf *entryLeaf = tree->GetLeaf("skim_entry");
   const std::vector<TString> names = MLTMVA::Features();
   std::vector<TLeaf *> featureLeaves(names.size());
   for (std::size_t i = 0; i < featureLeaves.size(); ++i)
      featureLeaves[i] = tree->GetLeaf(names[i]);
   for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
      tree->GetEntry(entry);
      std::vector<double> features(featureLeaves.size());
      for (std::size_t i = 0; i < features.size(); ++i) features[i] = featureLeaves[i]->GetValue();
      events.push_back({scoreLeaf->GetValue(), weightLeaf->GetValue(), static_cast<int>(classLeaf->GetValue()),
                        features, massLeaf->GetValue(), std::llround(entryLeaf->GetValue())});
   }
   return events;
}

// TMVA's test sample divided at random into the validation and test parts. The fixed
// seed gives TMVA_optimize.C and TMVA_train.C the same validation events.
std::pair<std::vector<ScoredEvent>, std::vector<ScoredEvent>> SplitValidationTest(
   const std::vector<ScoredEvent> &events)
{
   const double validationShare =
      MLTMVA::VALIDATION_FRACTION / (MLTMVA::VALIDATION_FRACTION + MLTMVA::TEST_FRACTION);
   TRandom3 random(42);
   std::pair<std::vector<ScoredEvent>, std::vector<ScoredEvent>> parts;
   for (const auto &event : events)
      (random.Rndm() < validationShare ? parts.first : parts.second).push_back(event);
   return parts;
}

double WeightSum(const std::vector<ScoredEvent> &events, int classID)
{
   double sum = 0.0;
   for (const auto &event : events)
      if (event.classID == classID) sum += event.weight;
   return sum;
}

void BalanceClassWeights(std::vector<ScoredEvent> &events)
{
   const double signalSum = WeightSum(events, 0);
   const double backgroundSum = WeightSum(events, 1);
   const double targetClassSum = 0.5 * events.size();
   if (signalSum <= 0.0 || backgroundSum <= 0.0) return;
   for (auto &event : events)
      event.weight *= event.classID == 0 ? targetClassSum / signalSum : targetClassSum / backgroundSum;
}

struct RocPoint {
   double fpr;
   double tpr;
   double threshold;   // score >= threshold
};

std::pair<std::vector<RocPoint>, double> BuildROC(const std::vector<ScoredEvent> &events)
{
   std::vector<ScoredEvent> ordered = events;
   std::sort(ordered.begin(), ordered.end(), [](const ScoredEvent &a, const ScoredEvent &b) {
      return a.score > b.score;
   });
   const double totalSignal = WeightSum(ordered, 0);
   const double totalBackground = WeightSum(ordered, 1);
   std::vector<RocPoint> roc = {{0.0, 0.0, std::numeric_limits<double>::infinity()}};
   if (totalSignal <= 0.0 || totalBackground <= 0.0) return {roc, 0.0};

   double signal = 0.0;
   double background = 0.0;
   double auc = 0.0;
   std::size_t begin = 0;
   while (begin < ordered.size()) {
      std::size_t end = begin;
      while (end < ordered.size() && ordered[end].score == ordered[begin].score) {
         if (ordered[end].classID == 0)
            signal += ordered[end].weight;
         else
            background += ordered[end].weight;
         ++end;
      }
      const double fpr = background / totalBackground;
      const double tpr = signal / totalSignal;
      const auto &previous = roc.back();
      auc += 0.5 * (tpr + previous.tpr) * (fpr - previous.fpr);
      roc.push_back({fpr, tpr, ordered[begin].score});
      begin = end;
   }
   return {roc, auc};
}

double WeightedKS(const std::vector<ScoredEvent> &first, const std::vector<ScoredEvent> &second, int classID)
{
   std::vector<std::pair<double, double>> a;
   std::vector<std::pair<double, double>> b;
   for (const auto &event : first)
      if (event.classID == classID) a.push_back({event.score, event.weight});
   for (const auto &event : second)
      if (event.classID == classID) b.push_back({event.score, event.weight});
   std::sort(a.begin(), a.end());
   std::sort(b.begin(), b.end());
   const double sumA = std::accumulate(a.begin(), a.end(), 0.0,
                                      [](double sum, const auto &entry) { return sum + entry.second; });
   const double sumB = std::accumulate(b.begin(), b.end(), 0.0,
                                      [](double sum, const auto &entry) { return sum + entry.second; });
   if (sumA <= 0.0 || sumB <= 0.0) return 0.0;
   std::size_t ia = 0;
   std::size_t ib = 0;
   double cumulativeA = 0.0;
   double cumulativeB = 0.0;
   double distance = 0.0;
   while (ia < a.size() || ib < b.size()) {
      const double value = ib == b.size() || (ia < a.size() && a[ia].first <= b[ib].first)
                            ? a[ia].first
                            : b[ib].first;
      while (ia < a.size() && a[ia].first <= value) cumulativeA += a[ia++].second / sumA;
      while (ib < b.size() && b[ib].first <= value) cumulativeB += b[ib++].second / sumB;
      distance = std::max(distance, std::abs(cumulativeA - cumulativeB));
   }
   return distance;
}

// (sum w)^2 / sum w^2 of one class.
double EffectiveEntries(const std::vector<ScoredEvent> &events, int classID)
{
   double sum = 0.0;
   double sum2 = 0.0;
   for (const auto &event : events)
      if (event.classID == classID) {
         sum += event.weight;
         sum2 += event.weight * event.weight;
      }
   return sum * sum / sum2;
}

// Asymptotic KS p-value with the effective sample sizes.
double KSPValue(double distance, const std::vector<ScoredEvent> &first, const std::vector<ScoredEvent> &second,
                int classID)
{
   const double nFirst = EffectiveEntries(first, classID);
   const double nSecond = EffectiveEntries(second, classID);
   return TMath::KolmogorovProb(std::sqrt(nFirst * nSecond / (nFirst + nSecond)) * distance);
}

// Weighted fraction of one class with score >= threshold, and its binomial uncertainty.
std::pair<double, double> Efficiency(const std::vector<ScoredEvent> &events, int classID, double threshold)
{
   double total = 0.0;
   double passed = 0.0;
   double passed2 = 0.0;
   double failed2 = 0.0;
   for (const auto &event : events) {
      if (event.classID != classID) continue;
      total += event.weight;
      if (event.score >= threshold) {
         passed += event.weight;
         passed2 += event.weight * event.weight;
      } else {
         failed2 += event.weight * event.weight;
      }
   }
   const double efficiency = passed / total;
   const double variance =
      ((1.0 - efficiency) * (1.0 - efficiency) * passed2 + efficiency * efficiency * failed2) / (total * total);
   return {efficiency, std::sqrt(variance)};
}

double PunziSmin(double background)
{
   const double sqrtB = std::sqrt(background);
   const double a = MLTMVA::PUNZI_A;
   const double b = MLTMVA::PUNZI_B;
   return b * b / 2.0 + a * sqrtB + b / 2.0 * std::sqrt(b * b + 4.0 * a * sqrtB + 4.0 * background);
}

// Background expected in the Punzi signal window in the DATA skim (after the pre-ML cut), from
// the sidebands next to it; the same estimate as optimalCUT_X_punzi.C without a classifier cut.
double WindowBackground(TTree *data)
{
   const double m = MLTMVA::X3872_MASS;
   const double start = MLTMVA::SIDEBAND_START;
   const double width = MLTMVA::SIDEBAND_WIDTH;
   const TString sidebands = Form("(Bmass > %.6f && Bmass < %.6f) || (Bmass > %.6f && Bmass < %.6f)",
                                  m - start - width, m - start, m + start, m + start + width);
   const Long64_t count = data->GetEntries(sidebands);
   return count * MLTMVA::SIGNAL_HALF_WIDTH / width;
}

// FOMs on one sample and the overtraining test of the training sample against it.
struct Evaluation {
   double auc = 0.0;
   double trainAUC = 0.0;
   double aucGap = 0.0;
   double sigeff = 0.0;
   double pauc = 0.0;
   double punzi = 0.0;
   double punziThreshold = 0.0;
   double workingPoint = 0.0;
   double ksSignal = 0.0;
   double ksSignalPValue = 0.0;
   double ksBackground = 0.0;
   double ksBackgroundPValue = 0.0;
   double signalEffTrain = 0.0, signalEffTrainError = 0.0, signalEff = 0.0, signalEffError = 0.0;
   double backgroundEffTrain = 0.0, backgroundEffTrainError = 0.0, backgroundEff = 0.0, backgroundEffError = 0.0;
   bool passed = false;

   double Value(const TString &fom) const
   {
      return std::map<TString, double>{{"auc", auc}, {"sigeff", sigeff}, {"pauc", pauc}, {"punzi", punzi}}.at(fom);
   }
};

// The working point is the Punzi-optimal threshold of `other`, or `threshold` when given.
Evaluation Evaluate(const std::vector<ScoredEvent> &training, const std::vector<ScoredEvent> &other,
                    double windowBackground, double threshold = std::numeric_limits<double>::quiet_NaN())
{
   Evaluation result;
   const auto roc = BuildROC(other);
   result.auc = roc.second;
   result.trainAUC = BuildROC(training).second;
   result.aucGap = result.trainAUC - result.auc;
   result.punzi = -1.0;
   // Partial AUC up to PAUC_MAX_BKG_EFF, closed with the linearly interpolated tpr at the boundary.
   const double paucMax = MLTMVA::PAUC_MAX_BKG_EFF;
   for (std::size_t i = 1; i < roc.first.size() && roc.first[i - 1].fpr < paucMax; ++i) {
      const auto &previous = roc.first[i - 1];
      const auto &point = roc.first[i];
      const double fpr = std::min(point.fpr, paucMax);
      const double tpr = point.fpr <= paucMax
                            ? point.tpr
                            : previous.tpr + (point.tpr - previous.tpr) * (paucMax - previous.fpr) / (point.fpr - previous.fpr);
      result.pauc += 0.5 * (tpr + previous.tpr) * (fpr - previous.fpr) / paucMax;
   }
   for (const auto &point : roc.first)
      if (point.fpr <= MLTMVA::BKG_EFF_TARGET) result.sigeff = std::max(result.sigeff, point.tpr);
   // Punzi as optimalCUT_X_punzi.C, on the thresholds k * THRESHOLD_STEP: B = the background
   // candidates of the Punzi sidebands that pass, scaled to the window (windowBackground x their
   // passing fraction), among cuts keeping at least PUNZI_MIN_BACKGROUND of them.
   std::vector<double> sideband;
   for (const auto &event : other) {
      const double distance = std::abs(event.bmass - MLTMVA::X3872_MASS);
      if (event.classID == 1 && distance > MLTMVA::SIDEBAND_START && distance < MLTMVA::SIDEBAND_START + MLTMVA::SIDEBAND_WIDTH)
         sideband.push_back(event.score);
   }
   std::sort(sideband.begin(), sideband.end());
   const auto [lowest, highest] = std::minmax_element(other.begin(), other.end(),
                                                      [](const auto &a, const auto &b) { return a.score < b.score; });
   for (double k = std::floor(lowest->score / MLTMVA::THRESHOLD_STEP); k <= std::ceil(highest->score / MLTMVA::THRESHOLD_STEP); ++k) {
      const double threshold = k * MLTMVA::THRESHOLD_STEP;
      const double tpr = Efficiency(other, 0, threshold).first;
      const double passing = sideband.end() - std::lower_bound(sideband.begin(), sideband.end(), threshold);
      const double punzi = tpr / PunziSmin(windowBackground * passing / sideband.size());
      if (passing >= MLTMVA::PUNZI_MIN_BACKGROUND && punzi > result.punzi) {
         result.punzi = punzi;
         result.punziThreshold = threshold;
      }
   }
   result.workingPoint = std::isnan(threshold) ? result.punziThreshold : threshold;
   result.ksSignal = WeightedKS(training, other, 0);
   result.ksBackground = WeightedKS(training, other, 1);
   result.ksSignalPValue = KSPValue(result.ksSignal, training, other, 0);
   result.ksBackgroundPValue = KSPValue(result.ksBackground, training, other, 1);
   std::tie(result.signalEffTrain, result.signalEffTrainError) = Efficiency(training, 0, result.workingPoint);
   std::tie(result.signalEff, result.signalEffError) = Efficiency(other, 0, result.workingPoint);
   std::tie(result.backgroundEffTrain, result.backgroundEffTrainError) = Efficiency(training, 1, result.workingPoint);
   std::tie(result.backgroundEff, result.backgroundEffError) = Efficiency(other, 1, result.workingPoint);
   result.passed = result.aucGap <= MLTMVA::MAX_AUC_GAP && result.ksSignalPValue >= MLTMVA::KS_PVALUE_MIN &&
                   result.ksBackgroundPValue >= MLTMVA::KS_PVALUE_MIN;
   return result;
}

void PrintEvaluation(const TString &label, const Evaluation &result)
{
   std::cout << Form("%s: AUC %.5f (train %.5f, gap %+.5f), sigeff %.4f, pauc %.4f, Punzi FOM %.6f at threshold %.4f",
                     label.Data(), result.auc, result.trainAUC, result.aucGap, result.sigeff, result.pauc, result.punzi,
                     result.punziThreshold) << std::endl;
   std::cout << Form("%s: KS signal %.4f (p %.3f), KS background %.4f (p %.3f), overtraining test %s",
                     label.Data(), result.ksSignal, result.ksSignalPValue, result.ksBackground,
                     result.ksBackgroundPValue, result.passed ? "passed" : "FAILED") << std::endl;
   std::cout << Form("%s: signal efficiency at working point %.4f: train %.4f +- %.4f, %s %.4f +- %.4f",
                     label.Data(), result.workingPoint, result.signalEffTrain, result.signalEffTrainError,
                     label.Data(), result.signalEff, result.signalEffError) << std::endl;
   std::cout << Form("%s: background efficiency at working point %.4f: train %.4f +- %.4f, %s %.4f +- %.4f",
                     label.Data(), result.workingPoint, result.backgroundEffTrain, result.backgroundEffTrainError,
                     label.Data(), result.backgroundEff, result.backgroundEffError) << std::endl;
}

double EvaluateTree(const TMVA::DecisionTree *tree, const std::vector<double> &features)
{
   const TMVA::DecisionTreeNode *node = tree->GetRoot();
   while (node && node->GetSelector() >= 0) {
      const std::size_t selector = node->GetSelector();
      bool goesRight = features[selector] >= node->GetCutValue();
      if (!node->GetCutType()) goesRight = !goesRight;
      node = goesRight ? node->GetRight() : node->GetLeft();
   }
   return node ? node->GetResponse() : 0.0;
}

// Evaluation of train against validation after the first n trees, for every n in
// `checkpoints`. A gradient-boosted forest with N trees contains every smaller model:
// the score after n trees is the sum of the first n tree responses. The BDTG output is
// a monotonic function of that sum, so every ROC-based quantity is the same.
std::vector<std::pair<int, Evaluation>> EvaluateStaged(const TMVA::MethodBDT &method,
                                                       std::vector<ScoredEvent> training,
                                                       std::vector<ScoredEvent> validation,
                                                       const std::vector<int> &checkpoints,
                                                       double windowBackground)
{
   for (auto &event : training) event.score = 0.0;
   for (auto &event : validation) event.score = 0.0;
   std::vector<std::pair<int, Evaluation>> result;
   const auto &forest = method.GetForest();
   std::size_t next = 0;
   for (std::size_t i = 0; i < forest.size() && next < checkpoints.size(); ++i) {
      for (auto &event : training) event.score += EvaluateTree(forest[i], event.features);
      for (auto &event : validation) event.score += EvaluateTree(forest[i], event.features);
      if (static_cast<int>(i + 1) != checkpoints[next]) continue;
      result.push_back({checkpoints[next], Evaluate(training, validation, windowBackground)});
      ++next;
   }
   return result;
}

} // namespace

#endif
