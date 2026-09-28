#include "TMVA_config.h"
#include "TMVA_common.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "TCut.h"
#include "TFile.h"
#include "TMemFile.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TTree.h"

#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/MethodBDT.h"
#include "TMVA/Tools.h"

void TMVA_optimize(TString sample = "",
                   TString fom = "",
                   int job = -1,
                   Long64_t maxSignal = 165000,
                   Long64_t maxBackground = 500000)
{
   sample = MLTMVA::SelectedSample(sample);
   fom = MLTMVA::SelectedFom(fom);
   if (job < 0 && gSystem->Getenv("ML_JOB")) job = std::atoi(gSystem->Getenv("ML_JOB"));
   if (job < 0) job = 0;
   const auto &config = MLTMVA::GetSampleConfig(sample);
   // Each job trains the largest tree count once and evaluates every smaller
   // tree count in OPTIMIZE_N_TREES from the same forest. The job number
   // therefore only indexes the other four hyperparameters.
   const int gridSize = MLTMVA::OPTIMIZE_MAX_DEPTH.size() * MLTMVA::OPTIMIZE_SHRINKAGE.size() *
                        MLTMVA::OPTIMIZE_BAGGED_SAMPLE_FRACTION.size() *
                        MLTMVA::OPTIMIZE_MIN_NODE_SIZE_PERCENT.size();
   // ML_RANDOM_SEARCH=1: the job number seeds a random draw from the RANDOM_* ranges instead.
   const bool randomSearch = gSystem->Getenv("ML_RANDOM_SEARCH") && TString(gSystem->Getenv("ML_RANDOM_SEARCH")) == "1";
   if (!randomSearch && job >= gridSize) {
      std::cerr << "Job " << job << " is outside the scan grid of " << gridSize << " points." << std::endl;
      gSystem->Exit(1);
   }
   const std::vector<int> &checkpoints = MLTMVA::OPTIMIZE_N_TREES;
   int nTrees = checkpoints.back();
   int maxDepth, nCuts = MLTMVA::OPTIMIZE_N_CUTS;
   double shrinkage, baggedSampleFraction, minNodeSize;
   if (randomSearch) {
      TRandom3 random(job + 1);
      auto logUniform = [&random](double low, double high) { return low * std::pow(high / low, random.Rndm()); };
      maxDepth = MLTMVA::RANDOM_MAX_DEPTH_MIN + random.Integer(MLTMVA::RANDOM_MAX_DEPTH_MAX - MLTMVA::RANDOM_MAX_DEPTH_MIN + 1);
      shrinkage = logUniform(MLTMVA::RANDOM_SHRINKAGE_MIN, MLTMVA::RANDOM_SHRINKAGE_MAX);
      baggedSampleFraction = MLTMVA::RANDOM_BAGGING_MIN + (MLTMVA::RANDOM_BAGGING_MAX - MLTMVA::RANDOM_BAGGING_MIN) * random.Rndm();
      minNodeSize = logUniform(MLTMVA::RANDOM_MIN_NODE_PERCENT_MIN, MLTMVA::RANDOM_MIN_NODE_PERCENT_MAX);
   } else {
      int index = job;
      maxDepth = MLTMVA::OPTIMIZE_MAX_DEPTH[index % MLTMVA::OPTIMIZE_MAX_DEPTH.size()];
      index /= MLTMVA::OPTIMIZE_MAX_DEPTH.size();
      shrinkage = MLTMVA::OPTIMIZE_SHRINKAGE[index % MLTMVA::OPTIMIZE_SHRINKAGE.size()];
      index /= MLTMVA::OPTIMIZE_SHRINKAGE.size();
      baggedSampleFraction =
         MLTMVA::OPTIMIZE_BAGGED_SAMPLE_FRACTION[index % MLTMVA::OPTIMIZE_BAGGED_SAMPLE_FRACTION.size()];
      index /= MLTMVA::OPTIMIZE_BAGGED_SAMPLE_FRACTION.size();
      minNodeSize = MLTMVA::OPTIMIZE_MIN_NODE_SIZE_PERCENT[index % MLTMVA::OPTIMIZE_MIN_NODE_SIZE_PERCENT.size()];
   }
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

   TMVA::Tools::Instance();
   TMemFile output("tmva_optimization", "RECREATE");
   auto *factory = new TMVA::Factory(
      "TMVAOptimization", &output,
      "!V:Silent:!Color:!DrawProgressBar:!ModelPersistence:AnalysisType=Classification");
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

   auto *method = dynamic_cast<TMVA::MethodBDT *>(factory->GetMethod("dataset", MLTMVA::METHOD_NAME));
   output.Write();
   std::vector<ScoredEvent> training = ReadScoredEvents(output.Get<TTree>("dataset/TrainTree"));
   std::vector<ScoredEvent> validation =
      SplitValidationTest(ReadScoredEvents(output.Get<TTree>("dataset/TestTree"))).first;
   BalanceClassWeights(training);
   BalanceClassWeights(validation);
   const auto staged = EvaluateStaged(*method, training, validation, checkpoints, windowBackground);
   if (staged.size() != checkpoints.size()) {
      std::cerr << "Forest has " << method->GetForest().size() << " trees, fewer than the checkpoints." << std::endl;
      gSystem->Exit(1);
   }
   delete factory;
   delete loader;

   // ML_SUMMARY_PATH (feature_study.py) or the optimization folder of the run.
   TString summaryPath = gSystem->Getenv("ML_SUMMARY_PATH") ? gSystem->Getenv("ML_SUMMARY_PATH") : "";
   if (summaryPath.IsNull()) {
      const TString summaryDirectory = MLTMVA::OutputDir(config, fom) + "/optimization";
      gSystem->mkdir(summaryDirectory, true);
      summaryPath = Form("%s/summary_job_%d.root", summaryDirectory.Data(), job);
   }
   TFile summary(summaryPath, "RECREATE");
   Evaluation result;
   int passed = 0;
   double objective = 0.0;
   TTree trial("trial", "TMVA BDTG optimization trial, train vs validation");
   trial.Branch("job", &job);
   trial.Branch("n_trees", &nTrees);
   trial.Branch("max_depth", &maxDepth);
   trial.Branch("shrinkage", &shrinkage);
   trial.Branch("bagged_sample_fraction", &baggedSampleFraction);
   trial.Branch("min_node_size_percent", &minNodeSize);
   trial.Branch("n_cuts", &nCuts);
   trial.Branch("train_auc", &result.trainAUC);
   trial.Branch("validation_auc", &result.auc);
   trial.Branch("auc_gap", &result.aucGap);
   trial.Branch("sigeff", &result.sigeff);
   trial.Branch("pauc", &result.pauc);
   trial.Branch("punzi", &result.punzi);
   trial.Branch("ks_signal_pvalue", &result.ksSignalPValue);
   trial.Branch("ks_background_pvalue", &result.ksBackgroundPValue);
   trial.Branch("signal_eff_train", &result.signalEffTrain);
   trial.Branch("signal_eff", &result.signalEff);
   trial.Branch("background_eff_train", &result.backgroundEffTrain);
   trial.Branch("background_eff", &result.backgroundEff);
   trial.Branch("passed", &passed);
   trial.Branch("objective", &objective);
   std::cout << "Optimization job: " << job << ", FOM " << fom << std::endl;
   std::cout << "NTrees   train AUC   val. AUC   sigeff    pauc     Punzi FOM   KS p sig/bkg   passed" << std::endl;
   for (const auto &[checkpoint, evaluation] : staged) {
      nTrees = checkpoint;
      result = evaluation;
      passed = result.passed;
      objective = result.Value(fom);
      trial.Fill();
      std::cout << Form("%6d   %9.6f   %8.6f   %7.4f   %6.4f   %9.6f   %5.3f/%5.3f    %d", nTrees, result.trainAUC,
                        result.auc, result.sigeff, result.pauc, result.punzi, result.ksSignalPValue,
                        result.ksBackgroundPValue, passed) << std::endl;
   }
   trial.Write();
   TNamed(TString("sample"), sample).Write();
   TNamed(TString("fom"), fom).Write();
   TNamed(TString("method_options"), methodOptions).Write();
   TNamed(TString("features"), MLTMVA::FeatureList(features)).Write();
   TNamed(TString("signal_cut"), preCut).Write();
   TNamed(TString("random_search"), TString(randomSearch ? "1" : "0")).Write();
   TNamed(TString("background_sideband_cut"), MLTMVA::BACKGROUND_SIDEBAND_CUT).Write();
   TNamed(TString("max_signal"), TString::Format("%lld", maxSignal)).Write();
   TNamed(TString("max_background"), TString::Format("%lld", maxBackground)).Write();
   TNamed(TString("window_background"), TString::Format("%.6g", windowBackground)).Write();
   summary.Close();

   std::cout << "Options: " << methodOptions << std::endl;
}
