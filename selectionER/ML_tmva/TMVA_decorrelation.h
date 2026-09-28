#ifndef ML_TMVA_DECORRELATION_H
#define ML_TMVA_DECORRELATION_H

// Mass decorrelation of classifier inputs by quantile morphing, the C++ twin of
// ML_common/ml_decorrelation.py: it reads the same maps (one tree map_<feature> with the
// percentile level and the polynomial coefficients c<degree>..c0 in (m - m_X)) and applies the
// same piecewise-linear map from the percentiles at the candidate mass to those at m_X.

#include <algorithm>
#include <map>
#include <memory>
#include <vector>

#include "TFile.h"
#include "TKey.h"
#include "TObjString.h"
#include "TTree.h"

namespace MLTMVA {

class Decorrelation {
public:
   explicit Decorrelation(const TString &path)
   {
      std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
      referenceMass_ = file->Get<TObjString>("metadata/reference_mass")->GetString().Atof();
      degree_ = file->Get<TObjString>("metadata/degree")->GetString().Atoi();
      for (auto *object : *file->GetListOfKeys()) {
         const TString name = object->GetName();
         if (!name.BeginsWith("map_")) continue;
         TTree *tree = file->Get<TTree>(name);
         std::vector<double> c(degree_ + 1);
         for (int p = 0; p <= degree_; ++p) tree->SetBranchAddress(Form("c%d", p), &c[p]);
         auto &levels = coefficients_[name(4, name.Length() - 4)];
         for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);
            levels.push_back(c);
         }
      }
      for (const auto &[name, levels] : coefficients_) reference_[name] = Percentiles(name, referenceMass_);
   }

   // Non-decreasing background percentiles of one feature at mass m.
   std::vector<double> Percentiles(const TString &name, double bmass) const
   {
      const double x = bmass - referenceMass_;
      const auto &levels = coefficients_.at(name);
      std::vector<double> q(levels.size());
      for (std::size_t k = 0; k < levels.size(); ++k) {
         double value = 0.0;
         for (int p = degree_; p >= 0; --p) value = value * x + levels[k][p];
         q[k] = k == 0 ? value : std::max(value, q[k - 1]);
      }
      return q;
   }

   double Transform(const TString &name, double value, double bmass) const
   {
      const std::vector<double> q = Percentiles(name, bmass);
      const std::vector<double> &ref = reference_.at(name);
      const std::size_t n = q.size();
      const std::size_t k = std::upper_bound(q.begin(), q.end(), value) - q.begin();   // entries <= value
      if (k == 0) return value - q[0] + ref[0];
      if (k == n) return value - q[n - 1] + ref[n - 1];
      const double span = q[k] - q[k - 1];
      const double t = span > 0.0 ? (value - q[k - 1]) / span : 0.0;
      return ref[k - 1] + t * (ref[k] - ref[k - 1]);
   }

private:
   double referenceMass_ = 0.0;
   int degree_ = 0;
   std::map<TString, std::vector<std::vector<double>>> coefficients_;   // c0 .. c<degree> per level
   std::map<TString, std::vector<double>> reference_;
};

} // namespace MLTMVA

#endif
