#ifndef ML_TMVA_CONFIG_H
#define ML_TMVA_CONFIG_H

#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TString.h"

namespace MLTMVA {

struct SampleConfig {
   TString data;
   TString signal;
   TString signalNonprompt;
   TString spectator;
   TString spectatorNonprompt;
   TString outputDir;
   TString system;
};

static const std::map<std::string, SampleConfig> SAMPLE_CONFIGS = {
   {"PbPb18",
    {"root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_DATA.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_MC_X3872.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_MC_X3872_nonPrompt.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_MC_PSI2S.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_MC_PSI2S_nonPrompt.root",
     "tmva_outputs/ntmix_PbPb/pbpb18",
     "PbPb18"}},
   {"PbPb23",
    {"root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_X3872_nonPrompt.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC_PSI2S_nonPrompt.root",
     "tmva_outputs/ntmix_PbPb/pbpb23",
     "PbPb23"}},
   {"PbPb24",
    {"root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_DATA.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_X3872.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_X3872_nonPrompt.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb24/flat_ntmix_PbPb24_MC_PSI2S_nonPrompt.root",
     "tmva_outputs/ntmix_PbPb/pbpb24",
     "PbPb24"}},
   {"ppRef24",
    {"root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_DATA.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_X3872_nonPrompt.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S.root",
     "root://eosuser.cern.ch//eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24/flat_ntmix_ppRef_MC_PSI2S_nonPrompt.root",
     "tmva_outputs/ntmix_ppRef",
     "ppRef"}},
};

// Production inputs; ML_FEATURES (comma-separated) overrides them for one run. A name ending
// in _dc is the mass-decorrelated feature (branch of the skims, see selectionER/make_skims.py).
static const std::vector<TString> FEATURES = {
   "Btrk1dR",
   "Btrk2Pt",
   "Bmu1pt",
   "Bmu2pt",
   "BtrkPtimb",
   "Bmu1eta",
   "Bmu2eta",
   "Bujeta",
   "Btktketa"
};

static const TString DATA_TREE = "ntmix";
static const TString SIGNAL_TREE = "ntmix_X3872";
static const TString SPECTATOR_TREE = "ntmix_PSI2S";
static const TString MC_WEIGHT_BRANCH = "pThatreweight";
static const TString BACKGROUND_SIDEBAND_CUT =
   "((Bmass > 3.82164) && (Bmass < 3.85664)) || ((Bmass > 3.88664) && (Bmass < 3.92164))";
// Training reads the skims of selectionER/make_skims.py: the flat trees after each sample's
// pre-cut (ML_common/ml_samples.py), with the _dc features. The pre-cut is read from them.
static const TString SKIM_DIR = "/eos/user/h/hmarques/Analysis_CODES/selectionER/skims";
static const TString DECORRELATED_SUFFIX = "_dc";
static const TString METHOD_NAME = "BDTG";
static const TString SCORE_BRANCH = "Prediction";

// Model selection, as in ML_common/ml_fom.py (written by fresh_ML_scan.sh). The validation
// half of TMVA's test sample chooses the model; the other half is only used for the report.
static const std::vector<TString> FOMS = {"auc", "sigeff", "pauc", "punzi"};
static const double BKG_EFF_TARGET = 0.02;
static const double PAUC_MAX_BKG_EFF = 0.05;
static const double PUNZI_A = 5.0;
static const double PUNZI_B = 2.0;
static const double SIGNAL_HALF_WIDTH = 0.005;
static const double SIDEBAND_START = 0.015;
static const double SIDEBAND_WIDTH = 0.02;
// The Punzi maximum only among cuts keeping at least this many background candidates.
static const double PUNZI_MIN_BACKGROUND = 10.0;
// The Punzi FOM on thresholds k * THRESHOLD_STEP, the grid of optimalCUT_X_punzi.C.
static const double THRESHOLD_STEP = 0.02;
static const double KS_PVALUE_MIN = 0.05;
static const double MAX_AUC_GAP = 0.01;
static const double X3872_MASS = 3.87164;  // as in plotER/aux/masses.h
static const double VALIDATION_FRACTION = 0.2;
static const double TEST_FRACTION = 0.2;

// Grid v2 (Sep 2026), moved toward the PbPb23 v1 optimum (500 trees, depth 4,
// shrinkage 0.005, bagging 0.6, min node 1.5%, all at the v1 grid edges).
// Tree counts evaluated inside every scan job. Each job trains the largest one
// and reads the smaller ones off the same forest, so they add no jobs.
// Grid v3 of the tree counts (Sep 2026): the v1 and v2 optima sat at few trees, so the counts
// above 1000 are gone (each job trains 1000 trees instead of 3000, about 3x faster) and the
// low end starts at 25.
static const std::vector<int> OPTIMIZE_N_TREES = {25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 800, 1000};
static const std::vector<int> OPTIMIZE_MAX_DEPTH = {3, 4, 5, 6};
static const std::vector<double> OPTIMIZE_SHRINKAGE = {0.001, 0.002, 0.003, 0.005, 0.0075, 0.01};
static const std::vector<double> OPTIMIZE_BAGGED_SAMPLE_FRACTION = {0.4, 0.5, 0.6, 0.7};
static const std::vector<double> OPTIMIZE_MIN_NODE_SIZE_PERCENT = {1.0, 1.5, 2.5, 5.0};
static const int OPTIMIZE_N_CUTS = 20;
// Random search (ML_RANDOM_SEARCH=1, used by selectionER/feature_study.py): job j draws its
// settings from these ranges with seed j + 1, and evaluates every tree count above.
static const int RANDOM_MAX_DEPTH_MIN = 2, RANDOM_MAX_DEPTH_MAX = 8;
// Shrinkage x trees sets how far the forest learns: 0.002 x 1000 trees still reaches the
// v1 optimum (0.005 x 500), below that every draw would underfit.
static const double RANDOM_SHRINKAGE_MIN = 0.002, RANDOM_SHRINKAGE_MAX = 0.2;          // log-uniform
static const double RANDOM_BAGGING_MIN = 0.3, RANDOM_BAGGING_MAX = 1.0;
static const double RANDOM_MIN_NODE_PERCENT_MIN = 0.5, RANDOM_MIN_NODE_PERCENT_MAX = 10.0; // log-uniform
// Scan jobs = MAX_DEPTH x SHRINKAGE x BAGGED_SAMPLE_FRACTION x MIN_NODE_SIZE_PERCENT = 4 x 6 x 4 x 4 = 384.
// Keep njobs in submit_optimize.sub equal to this product.

inline const SampleConfig &GetSampleConfig(const TString &sample)
{
   return SAMPLE_CONFIGS.at(sample.Data());
}

inline TString SelectedSample(TString sample)
{
   if (sample.IsNull()) sample = gSystem->Getenv("ML_SAMPLE");
   if (sample.IsNull()) sample = "PbPb23";
   return sample;
}

inline TString SelectedFom(TString fom)
{
   if (fom.IsNull()) fom = gSystem->Getenv("ML_FOM");
   if (std::find(FOMS.begin(), FOMS.end(), fom) == FOMS.end())
      throw std::invalid_argument("The FOM must be auc, sigeff, pauc or punzi, not '" + std::string(fom.Data()) + "'");
   return fom;
}

inline TString SelectedKind(TString kind)
{
   if (kind.IsNull()) kind = gSystem->Getenv("ML_KIND");
   if (kind.IsNull()) kind = "all";
   return kind;
}

// Output folder of one sample and model-selection FOM, e.g. tmva_outputs/ntmix_PbPb/pbpb23_auc;
// ML_OUTPUT_DIR overrides it for one run.
inline TString OutputDir(const SampleConfig &config, const TString &fom)
{
   if (gSystem->Getenv("ML_OUTPUT_DIR")) return gSystem->Getenv("ML_OUTPUT_DIR");
   return config.outputDir + "_" + fom;
}

inline std::vector<TString> Features()
{
   if (!gSystem->Getenv("ML_FEATURES")) return FEATURES;
   std::vector<TString> features;
   std::unique_ptr<TObjArray> parts(TString(gSystem->Getenv("ML_FEATURES")).Tokenize(","));
   for (auto *part : *parts) features.push_back(static_cast<TObjString *>(part)->GetString());
   return features;
}

inline TString FeatureList(const std::vector<TString> &features)
{
   TString list;
   for (const auto &feature : features) list += (list.IsNull() ? "" : ",") + feature;
   return list;
}

// ML_USE_XROOTD=1 (batch jobs) reads /eos/ files through XRootD, as ml_samples.resolve_input_path.
inline TString ResolveInputPath(const TString &path)
{
   if (gSystem->Getenv("ML_USE_XROOTD") && TString(gSystem->Getenv("ML_USE_XROOTD")) == "1" && path.BeginsWith("/eos/"))
      return "root://eosuser.cern.ch/" + path;
   return path;
}

// Skim of one sample and kind (DATA, MC_X3872, MC_X3872_NONPROMPT, MC_PSI2S, MC_PSI2S_NONPROMPT).
inline TString SkimPath(const TString &sample, const TString &kind)
{
   return ResolveInputPath(SKIM_DIR + "/" + sample + "/skim_" + sample + "_" + kind + ".root");
}

inline TString DecorrelationPath(const TString &sample)
{
   return ResolveInputPath(SKIM_DIR + "/" + sample + "/decorrelation_" + sample + ".root");
}

inline TString BaseFeature(const TString &name)
{
   return name.EndsWith(DECORRELATED_SUFFIX) ? TString(name(0, name.Length() - DECORRELATED_SUFFIX.Length())) : name;
}

inline TString SkimPreCut(const TString &sample)
{
   std::unique_ptr<TFile> skim(TFile::Open(SkimPath(sample, "DATA"), "READ"));
   return skim->Get<TObjString>("metadata/pre_cut")->GetString();
}

inline TString ModelPath(const TString &outputDir)
{
   return outputDir + "/dataset/weights/TMVAClassification_BDTG.weights.xml";
}

} // namespace MLTMVA

#endif
