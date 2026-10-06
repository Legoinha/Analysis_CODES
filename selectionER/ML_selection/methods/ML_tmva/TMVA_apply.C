#include "TMVA_config.h"
#include "TMVA_decorrelation.h"

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TNamed.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TTree.h"

#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

// Scores the flat files of case caseName ("<tree>/<year>", default ML_CASE) with the model trained with
// FOM fom (default ML_FOM): every input of the case, or only the one of kind (DATA, MC, MC_X3872, ...;
// default ML_KIND, else all). Each scored file is the flat tree plus the Prediction branch, with a
// metadata/ directory of the same keys as ML_pytorch/NN_apply.py and ML_xgboost/XGB_apply.py. Output:
// ML_SCORED_DIR, or selectionER/scored_samples/<tree>/<year>/flat_<tree>_<year>_tmva_<fom>_scored_<KIND>.root.
// A <branch>_log input is log10 of its flat branch, a _dc input is computed from its base input and Bmass
// with the maps of the training skims (TMVA_decorrelation.h); the scored tree keeps only the flat branches.
void TMVA_apply(TString caseName = "", TString fom = "", TString kind = "")
{
   const MLTMVA::Run run = MLTMVA::LoadRun(caseName, fom);
   const MLSettings::Case &c = run.settings;
   if (kind.IsNull()) kind = gSystem->Getenv("ML_KIND");
   if (kind.IsNull()) kind = "all";
   const TString scoredDir = gSystem->Getenv("ML_SCORED_DIR") ? TString(gSystem->Getenv("ML_SCORED_DIR"))
      : "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/" + c.tree + "/" + c.year;
   gSystem->mkdir(scoredDir, true);

   // The inputs and the pre-ML cut are read from the report of the training that wrote this model.
   TFile report(run.outputDir + "/training_report.root", "READ");
   std::vector<TString> features;
   std::unique_ptr<TObjArray> parts(TString(report.Get<TNamed>("features")->GetTitle()).Tokenize(","));
   for (auto *part : *parts) features.push_back(static_cast<TObjString *>(part)->GetString());
   std::vector<std::pair<TString, TString>> metadata = {
      {"classifier", "tmva"},
      {"sample", run.name},
      {"fom", run.fom},
      {"model", gSystem->IsAbsoluteFileName(run.model) ? run.model : "ML_selection/methods/ML_tmva/" + run.model},
      {"features", report.Get<TNamed>("features")->GetTitle()},
      {"pre_cut", report.Get<TNamed>("signal_cut")->GetTitle()},
      {"mc_weight", report.Get<TNamed>("mc_weight_branch")->GetTitle()},
      {"score_branch", MLTMVA::SCORE_BRANCH},
      {"score_min", "-1"},
      {"score_max", "1"},
   };
   std::unique_ptr<MLTMVA::Decorrelation> decorrelation;
   for (const auto &feature : features)
      if (feature.EndsWith("_dc") && !decorrelation) {
         decorrelation = std::make_unique<MLTMVA::Decorrelation>(run.maps);
         metadata.push_back({"decorrelation_maps", run.maps});
      }

   for (const auto &[key, input] : c.inputs) {
      if (kind != "all" && kind != input.kind) continue;
      std::vector<float> values(features.size());
      float bmass = 0.0f;
      float skimEntry = 0.0f;
      TMVA::Reader reader("!Color:!Silent");
      for (std::size_t i = 0; i < features.size(); ++i) reader.AddVariable(features[i], &values[i]);
      reader.AddSpectator("Bmass", &bmass);
      reader.AddSpectator("skim_entry := Entry$", &skimEntry);
      reader.BookMVA(MLTMVA::METHOD_NAME, run.model);

      TFile *inputFile = TFile::Open(run.xrootd + input.file, "READ");
      TTree *inputTree = inputFile->Get<TTree>(input.tree);
      // Per input: its name without _dc (the key of its decorrelation map) and its flat branch.
      std::vector<TString> inputs, branches;
      for (const auto &feature : features) {
         inputs.push_back(TString(feature).ReplaceAll("_dc", ""));
         branches.push_back(inputs.back().EndsWith("_log") ? TString(inputs.back()(0, inputs.back().Length() - 4)) : inputs.back());
      }
      std::map<TString, float> raw;
      for (const auto &branch : branches) raw[branch] = 0.0f;
      for (auto &[name, value] : raw) inputTree->SetBranchAddress(name, &value);
      inputTree->SetBranchAddress("Bmass", &bmass);

      const TString outputPath = scoredDir + "/flat_" + c.tree + "_" + c.year + "_tmva_" + run.fom + "_scored_" + input.kind + ".root";
      TFile output(outputPath, "RECREATE");
      TTree *outputTree = inputTree->CloneTree(0);
      float prediction = 0.0f;
      outputTree->Branch(MLTMVA::SCORE_BRANCH, &prediction);

      std::cout << "Scoring " << input.label << ": " << input.file << std::endl;
      for (Long64_t entry = 0; entry < inputTree->GetEntries(); ++entry) {
         inputTree->GetEntry(entry);
         for (std::size_t i = 0; i < features.size(); ++i) {
            const double value = branches[i] != inputs[i] ? std::log10(raw[branches[i]]) : raw[branches[i]];
            values[i] = features[i].EndsWith("_dc") ? decorrelation->Transform(inputs[i], value, bmass) : value;
         }
         prediction = reader.EvaluateMVA(MLTMVA::METHOD_NAME);
         outputTree->Fill();
         if ((entry + 1) % 1000000 == 0)
            std::cout << "  processed entries: " << entry + 1 << std::endl;
      }
      outputTree->Write();

      TDirectory *metadataDirectory = output.mkdir("metadata");
      auto fileMetadata = metadata;
      fileMetadata.insert(fileMetadata.end(), {{"kind", input.kind}, {"tree", input.tree}, {"input_file", input.file}});
      for (const auto &[name, value] : fileMetadata) {
         TObjString text(value);
         metadataDirectory->WriteTObject(&text, name);
      }
      output.Close();
      inputFile->Close();
      std::cout << "  saved: " << outputPath << std::endl;
   }
}
