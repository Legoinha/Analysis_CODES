#ifndef ML_TMVA_CONFIG_H
#define ML_TMVA_CONFIG_H

#include <map>
#include <memory>
#include <vector>

#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TString.h"

// The cases, features, cuts, mass windows and model-selection settings, written from
// ML_selection/settings.py by ml_samples.write_cpp_settings() (every Python script of the chain
// writes it again when it starts).
#include "ml_settings.h"

namespace MLTMVA {

static const TString METHOD_NAME = "BDTG";
static const TString SCORE_BRANCH = "Prediction";
// The validation half of TMVA's test sample chooses the model; the other half is only used for the report.
static const std::vector<TString> FOMS = {"auc", "sigeff", "pauc", "punzi"};
static const double VALIDATION_FRACTION = 0.2;
static const double TEST_FRACTION = 0.2;

// Grid v2 (Sep 2026), moved toward the PbPb23 v1 optimum (500 trees, depth 4,
// shrinkage 0.005, bagging 0.6, min node 1.5%, all at the v1 grid edges).
// Tree counts evaluated inside every scan job. Each job trains the largest one
// and reads the smaller ones off the same forest, so they add no jobs.
// Grid v3 of the tree counts (Sep 2026): the v1 and v2 optima sat at few trees, so the counts
// above 1000 were dropped and the low end starts at 25. Oct 2026: back up to 2000, since 13-18% of the
// best B-meson trials sat at 1000 (each job trains 2000 trees, about twice as long).
static const std::vector<int> OPTIMIZE_N_TREES = {25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000};
static const std::vector<int> OPTIMIZE_MAX_DEPTH = {3, 4, 5, 6};
static const std::vector<double> OPTIMIZE_SHRINKAGE = {0.001, 0.002, 0.003, 0.005, 0.0075, 0.01};
static const std::vector<double> OPTIMIZE_BAGGED_SAMPLE_FRACTION = {0.4, 0.5, 0.6, 0.7};
static const std::vector<double> OPTIMIZE_MIN_NODE_SIZE_PERCENT = {1.0, 1.5, 2.5, 5.0};
static const int OPTIMIZE_N_CUTS = 20;
// Random search (ML_RANDOM_SEARCH=1, used by ML_selection/scripts/feature_study.py): job j draws its
// settings from these ranges with seed j + 1, and evaluates every tree count above.
// Oct 2026: depth up to 10 and shrinkage up to 0.5, where ~20-26% of the best B-meson trials sat at the
// old edges (8 and 0.2).
static const int RANDOM_MAX_DEPTH_MIN = 2, RANDOM_MAX_DEPTH_MAX = 10;
// Shrinkage x trees sets how far the forest learns: 0.002 x 1000 trees still reaches the
// v1 optimum (0.005 x 500), below that every draw would underfit.
static const double RANDOM_SHRINKAGE_MIN = 0.002, RANDOM_SHRINKAGE_MAX = 0.5;          // log-uniform
static const double RANDOM_BAGGING_MIN = 0.3, RANDOM_BAGGING_MAX = 1.0;
static const double RANDOM_MIN_NODE_PERCENT_MIN = 0.5, RANDOM_MIN_NODE_PERCENT_MAX = 10.0; // log-uniform
// Scan jobs = MAX_DEPTH x SHRINKAGE x BAGGED_SAMPLE_FRACTION x MIN_NODE_SIZE_PERCENT = 4 x 6 x 4 x 4 = 384.
// Keep njobs in submit_optimize.sub equal to this product.

// One case "<tree>/<year>" (ML_CASE when not given) and model-selection FOM (ML_FOM) at run time:
// its settings (ml_settings.h, those of the test variant ML_VARIANT when set); the skim of every input
// and the decorrelation maps, through XRootD in batch jobs (ML_USE_XROOTD=1); the inputs of the
// run (ML_FEATURES, comma-separated, or the production features of the tree); the pre-cut the skims
// were made with; and the output folder (ML_OUTPUT_DIR, or tmva_outputs/<tree>/<year>_<fom>).
struct Run {
   TString name;
   TString fom;
   MLSettings::Case settings;
   TString xrootd;   // prefix of every /eos/ path
   std::map<TString, TString> skims;
   TString maps;
   std::vector<TString> features;
   TString featureList;   // comma-separated
   TString preCut;
   TString outputDir;
   TString model;   // the BDTG weights in the output folder
};

inline Run LoadRun(TString name, TString fom)
{
   Run run;
   run.name = name.IsNull() ? TString(gSystem->Getenv("ML_CASE")) : name;
   run.fom = fom.IsNull() ? TString(gSystem->Getenv("ML_FOM")) : fom;
   const TString variant = gSystem->Getenv("ML_VARIANT");
   run.settings = MLSettings::CASES.at(variant.IsNull() ? run.name : run.name + ":" + variant);
   run.xrootd = TString(gSystem->Getenv("ML_USE_XROOTD")) == "1" ? "root://eosuser.cern.ch/" : "";
   const TString folder = run.xrootd + run.settings.skimDir + "/";
   const TString stem = run.settings.tree + "_" + run.settings.year;
   for (const auto &[key, input] : run.settings.inputs) run.skims[key] = folder + "skim_" + stem + "_" + input.kind + ".root";
   run.maps = folder + "decorrelation_" + stem + ".root";
   run.features = run.settings.features;
   if (gSystem->Getenv("ML_FEATURES")) {
      run.features.clear();
      std::unique_ptr<TObjArray> parts(TString(gSystem->Getenv("ML_FEATURES")).Tokenize(","));
      for (auto *part : *parts) run.features.push_back(static_cast<TObjString *>(part)->GetString());
   }
   for (const auto &feature : run.features) run.featureList += (run.featureList.IsNull() ? "" : ",") + feature;
   std::unique_ptr<TFile> skim(TFile::Open(run.skims.at("data"), "READ"));
   run.preCut = skim->Get<TObjString>("metadata/pre_cut")->GetString();
   run.outputDir = gSystem->Getenv("ML_OUTPUT_DIR") ? TString(gSystem->Getenv("ML_OUTPUT_DIR"))
                                                    : "tmva_outputs/" + run.settings.tree + "/" + run.settings.year + "_" + run.fom;
   run.model = run.outputDir + "/dataset/weights/TMVAClassification_BDTG.weights.xml";
   return run;
}

} // namespace MLTMVA

#endif
