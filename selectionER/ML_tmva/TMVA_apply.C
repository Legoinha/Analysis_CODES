#include "TMVA_config.h"
#include "TMVA_decorrelation.h"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

namespace {

struct ScoreTarget {
   TString label;
   TString inputFile;
   TString treeName;
   TString outputKind;
};

// Key-value pairs written to the metadata/ directory of every scored file,
// with the same keys as ML_pytorch/NN_apply.py and ML_xgboost/XGB_apply.py.
using Metadata = std::vector<std::pair<TString, TString>>;

Metadata ModelMetadata(const TString &sample, const TString &fom, const TString &outputDir)
{
   // The pre-ML cut and the inputs are read from the report of the training that wrote this model.
   TFile report(outputDir + "/training_report.root", "READ");
   return {
      {"classifier", "tmva"},
      {"sample", sample},
      {"fom", fom},
      {"model", gSystem->IsAbsoluteFileName(outputDir) ? MLTMVA::ModelPath(outputDir) : "ML_tmva/" + MLTMVA::ModelPath(outputDir)},
      {"features", report.Get<TNamed>("features")->GetTitle()},
      {"pre_cut", report.Get<TNamed>("signal_cut")->GetTitle()},
      {"mc_weight", report.Get<TNamed>("mc_weight_branch")->GetTitle()},
      {"score_branch", MLTMVA::SCORE_BRANCH},
      {"score_min", "-1"},
      {"score_max", "1"},
   };
}

std::vector<TString> ModelFeatures(const TString &outputDir)
{
   TFile report(outputDir + "/training_report.root", "READ");
   std::vector<TString> features;
   std::unique_ptr<TObjArray> parts(TString(report.Get<TNamed>("features")->GetTitle()).Tokenize(","));
   for (auto *part : *parts) features.push_back(static_cast<TObjString *>(part)->GetString());
   return features;
}

std::vector<ScoreTarget> Targets(const MLTMVA::SampleConfig &config)
{
   return {
      {"DATA", config.data, MLTMVA::DATA_TREE, "DATA"},
      {"prompt X3872 MC", config.signal, MLTMVA::SIGNAL_TREE, "MC_X3872"},
      {"nonprompt X3872 MC", config.signalNonprompt, MLTMVA::SIGNAL_TREE, "MC_X3872_NONPROMPT"},
      {"prompt Psi2S MC", config.spectator, MLTMVA::SPECTATOR_TREE, "MC_PSI2S"},
      {"nonprompt Psi2S MC", config.spectatorNonprompt, MLTMVA::SPECTATOR_TREE, "MC_PSI2S_NONPROMPT"},
   };
}

// A _dc input is computed from its base branch and Bmass with the maps of the training skims
// (TMVA_decorrelation.h); the scored tree keeps only the flat branches.
void ScoreFile(const MLTMVA::SampleConfig &config, const TString &fom, const TString &outputDir,
               const ScoreTarget &target, Metadata metadata, const std::vector<TString> &features,
               const MLTMVA::Decorrelation *decorrelation)
{
   std::vector<float> values(features.size());
   float bmass = 0.0f;
   float skimEntry = 0.0f;
   TMVA::Reader reader("!Color:!Silent");
   for (std::size_t i = 0; i < features.size(); ++i) reader.AddVariable(features[i], &values[i]);
   reader.AddSpectator("Bmass", &bmass);
   reader.AddSpectator("skim_entry := Entry$", &skimEntry);
   reader.BookMVA(MLTMVA::METHOD_NAME, MLTMVA::ModelPath(outputDir));

   TFile *input = TFile::Open(target.inputFile, "READ");
   TTree *inputTree = input->Get<TTree>(target.treeName);
   std::map<TString, float> raw;
   for (const auto &feature : features) raw[MLTMVA::BaseFeature(feature)] = 0.0f;
   for (auto &[name, value] : raw) inputTree->SetBranchAddress(name, &value);
   inputTree->SetBranchAddress("Bmass", &bmass);

   const TString outputDirectory = "scored_samples";
   gSystem->mkdir(outputDirectory, true);
   const TString outputPath = outputDirectory + "/flat_ntmix_" + config.system + "_" + fom +
                              "_scored_" + target.outputKind + ".root";
   TFile output(outputPath, "RECREATE");
   TTree *outputTree = inputTree->CloneTree(0);
   float prediction = 0.0f;
   outputTree->Branch(MLTMVA::SCORE_BRANCH, &prediction);

   std::cout << "Scoring " << target.label << ": " << target.inputFile << std::endl;
   for (Long64_t entry = 0; entry < inputTree->GetEntries(); ++entry) {
      inputTree->GetEntry(entry);
      for (std::size_t i = 0; i < features.size(); ++i) {
         const float value = raw[MLTMVA::BaseFeature(features[i])];
         values[i] = features[i].EndsWith(MLTMVA::DECORRELATED_SUFFIX)
                        ? decorrelation->Transform(MLTMVA::BaseFeature(features[i]), value, bmass)
                        : value;
      }
      prediction = reader.EvaluateMVA(MLTMVA::METHOD_NAME);
      outputTree->Fill();
      if ((entry + 1) % 1000000 == 0)
         std::cout << "  processed entries: " << entry + 1 << std::endl;
   }
   outputTree->Write();

   metadata.push_back({"kind", target.outputKind});
   metadata.push_back({"tree", target.treeName});
   metadata.push_back({"input_file", target.inputFile});
   TDirectory *metadataDirectory = output.mkdir("metadata");
   for (const auto &[key, value] : metadata) {
      TObjString text(value);
      metadataDirectory->WriteTObject(&text, key);
   }
   output.Close();
   std::cout << "  saved: " << outputPath << std::endl;
}

} // namespace

void TMVA_apply(TString sample = "", TString fom = "", TString kind = "")
{
   sample = MLTMVA::SelectedSample(sample);
   fom = MLTMVA::SelectedFom(fom);
   kind = MLTMVA::SelectedKind(kind);
   const auto &config = MLTMVA::GetSampleConfig(sample);
   const TString outputDir = MLTMVA::OutputDir(config, fom);
   Metadata metadata = ModelMetadata(sample, fom, outputDir);
   const std::vector<TString> features = ModelFeatures(outputDir);
   std::unique_ptr<MLTMVA::Decorrelation> decorrelation;
   for (const auto &feature : features)
      if (feature.EndsWith(MLTMVA::DECORRELATED_SUFFIX) && !decorrelation) {
         decorrelation = std::make_unique<MLTMVA::Decorrelation>(MLTMVA::DecorrelationPath(sample));
         metadata.push_back({"decorrelation_maps", MLTMVA::DecorrelationPath(sample)});
      }
   for (const auto &target : Targets(config))
      if (kind == "all" || kind == target.outputKind)
         ScoreFile(config, fom, outputDir, target, metadata, features, decorrelation.get());
}
