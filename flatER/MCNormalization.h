#ifndef FLATER_MCNORMALIZATION_H
#define FLATER_MCNORMALIZATION_H

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace MCNormalization {

struct Context {
    std::string system;
    std::string tree;
    std::string particle;
    std::string promptness;
};

struct Entry {
    std::string system;
    std::string tree;
    std::string particle;
    std::string promptness;
    int pthat = -1;
    std::string pathPattern;
    double xsecPb = 0.;
    double filterEfficiency = 0.;
    long long nGenerated = 0;   // counted from Bfinder/ntGen of the flattened input files

    double Weight() const
    {
        return xsecPb * filterEfficiency / static_cast<double>(nGenerated);
    }
};

inline std::string Trim(const std::string &value)
{
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

inline std::string Lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

inline std::vector<std::string> SplitCsvLine(const std::string &line)
{
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(Trim(field));
    return fields;
}

inline Context BuildContext(const std::string &system,
                            const std::string &tree,
                            const std::string &particle,
                            const std::string &promptness)
{
    Context context;
    context.system = Lower(Trim(system));
    context.tree = Lower(Trim(tree));

    if (context.tree == "ntmix") {
        const std::string particleLower = Lower(particle);
        if (particleLower.find("psi2s") != std::string::npos) context.particle = "psi2s";
        else if (particleLower.find("x3872") != std::string::npos) context.particle = "x3872";
        else throw std::runtime_error("Unknown ntmix MC particle tag: '" + particle + "'");

        context.promptness = Lower(promptness).find("nonprompt") != std::string::npos
                           ? "nonprompt" : "prompt";
    } else if (context.tree == "ntkp") {
        context.particle = "bplus";
        context.promptness = "inclusive";
    } else if (context.tree == "ntphi") {
        context.particle = "bs";
        context.promptness = "inclusive";
    } else if (context.tree == "ntkstar") {
        context.particle = "bzero";
        context.promptness = "inclusive";
    } else {
        throw std::runtime_error("No MC normalization particle mapping for tree '" + tree + "'");
    }
    return context;
}

class Registry {
public:
    Registry() = default;

    explicit Registry(const std::string &fileName)
    {
        Load(fileName);
    }

    void Load(const std::string &fileName)
    {
        entries_.clear();
        std::ifstream input(fileName);
        if (!input.is_open()) {
            throw std::runtime_error("Cannot open MC normalization table: " + fileName);
        }

        std::string line;
        int lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;

            const auto fields = SplitCsvLine(line);
            if (Lower(fields.empty() ? "" : fields[0]) == "system") continue;
            if (fields.size() != 8) {
                throw std::runtime_error("Expected 8 columns in " + fileName + ":" +
                                         std::to_string(lineNumber));
            }

            Entry entry;
            try {
                entry.system = Lower(fields[0]);
                entry.tree = Lower(fields[1]);
                entry.particle = Lower(fields[2]);
                entry.promptness = Lower(fields[3]);
                entry.pthat = std::stoi(fields[4]);
                entry.pathPattern = Lower(fields[5]);
                entry.xsecPb = std::stod(fields[6]);
                entry.filterEfficiency = std::stod(fields[7]);
            } catch (const std::exception &error) {
                throw std::runtime_error("Invalid value in " + fileName + ":" +
                                         std::to_string(lineNumber) + " (" + error.what() + ")");
            }

            if (entry.system.empty() || entry.tree.empty() || entry.particle.empty() ||
                entry.promptness.empty() || entry.pathPattern.empty() || entry.pthat < 0 ||
                entry.xsecPb <= 0. || entry.filterEfficiency <= 0. ||
                entry.filterEfficiency > 1.) {
                throw std::runtime_error("Non-physical or empty value in " + fileName + ":" +
                                         std::to_string(lineNumber));
            }

            entries_.push_back(entry);
        }

        if (entries_.empty()) {
            throw std::runtime_error("MC normalization table contains no entries: " + fileName);
        }
    }

    // All rows of one MC group (system, tree, particle, promptness), sorted by pThat threshold.
    std::vector<Entry> Group(const Context &context) const
    {
        std::vector<Entry> group;
        for (const auto &entry : entries_)
            if (entry.system == context.system && entry.tree == context.tree &&
                entry.particle == context.particle && entry.promptness == context.promptness)
                group.push_back(entry);
        std::sort(group.begin(), group.end(),
                  [](const Entry &a, const Entry &b) { return a.pthat < b.pthat; });
        return group;
    }

private:
    std::vector<Entry> entries_;
};

// The pThat sample of one input file: the single group row whose path_pattern is in its path.
inline std::size_t SampleIndex(const std::vector<Entry> &group, const std::string &fileName)
{
    const std::string pathLower = Lower(fileName);
    std::vector<std::size_t> matches;
    for (std::size_t k = 0; k < group.size(); ++k)
        if (pathLower.find(group[k].pathPattern) != std::string::npos) matches.push_back(k);
    if (matches.size() != 1) {
        throw std::runtime_error("Expected exactly one MC normalization row for '" + fileName +
                                 "', found " + std::to_string(matches.size()));
    }
    return matches.front();
}

// Weight for merged inclusive "pThat > X" samples. An event with generator pThat p
// can come from every sample j with X_j < p, so the merged sample has luminosity
// sum_{X_j < p} L_j there, with L_j = n_gen_j / (xsec_j * filter_eff_j):
//     w(p) = 1 / sum_{X_j < p} L_j
// A sample without input files has n_gen = 0 and adds no luminosity.
class PthatWeight {
public:
    PthatWeight() = default;
    explicit PthatWeight(const std::vector<Entry> &group)
    {
        double luminosity = 0.;
        for (const auto &entry : group) {
            luminosity += 1. / entry.Weight();
            thresholds_.push_back(entry.pthat);
            cumulativeLuminosity_.push_back(luminosity);
        }
    }

    double operator()(float pthat) const
    {
        std::size_t k = 0;
        while (k + 1 < thresholds_.size() && pthat > thresholds_[k + 1]) ++k;
        return 1. / cumulativeLuminosity_[k];
    }

    std::vector<int> thresholds_;
    std::vector<double> cumulativeLuminosity_;
};

} // namespace MCNormalization

#endif
