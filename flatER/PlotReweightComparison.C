#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

struct SampleDescription {
    TString fileTag;
    TString plotLabel;
    TString group;
};

SampleDescription DescribeSample(const TString &sourceTree,
                                 const TString &particle,
                                 const TString &promptness)
{
    if (sourceTree == "ntmix") {
        TString particleTag = particle;
        particleTag.ReplaceAll("_", "");
        TString particleLabel = particleTag;
        if (particleTag.EqualTo("PSI2S", TString::kIgnoreCase)) {
            particleTag = "PSI2S";
            particleLabel = "#Psi(2S)";
        } else if (particleTag.EqualTo("X3872", TString::kIgnoreCase)) {
            particleTag = "X3872";
            particleLabel = "X(3872)";
        }

        const bool isNonprompt = promptness.Contains("nonprompt", TString::kIgnoreCase);
        const TString production = isNonprompt ? "nonprompt" : "prompt";
        return {particleTag + "_" + production,
                particleLabel + " " + production,
                "ntmix"};
    }

    if (sourceTree == "ntKp") return {"Bplus", "B^{+}", "Bmeson"};
    if (sourceTree == "ntphi") return {"Bs", "B_{s}^{0}", "Bmeson"};
    if (sourceTree == "ntKstar") return {"Bzero", "B^{0}", "Bmeson"};
    throw std::runtime_error("No comparison-plot label for source tree: " +
                             std::string(sourceTree.Data()));
}

TString SafeFileToken(TString value)
{
    value.ReplaceAll("abs(", "abs_");
    value.ReplaceAll(")", "");
    value.ReplaceAll("/", "_");
    value.ReplaceAll(" ", "_");
    return value;
}

TString VariableAxisLabel(const TString &variable)
{
    if (variable == "Bpt") return "p_{T} [GeV/c]";
    if (variable == "pthat") return "Generator #hat{p}_{T} [GeV/c]";
    if (variable == "CentBin") return "Centrality bin";
    if (variable.Contains("nChargedTracks")) return "Track multiplicity";
    return variable;
}

} // namespace

void PlotReweightComparison(TString inputFile = "",
                            TString flatTreeName = "",
                            TString sourceTree  = "ntmix",
                            TString systemName  = "ppRef",
                            TString particle   = "",
                            TString promptness = "",
                            TString outputBase = "reweighting_comparisons",
                            TString variable = "Bpt",
                            TString weightBranch = "pThatreweight",
                            Int_t nBins = 100)
{
    const bool isPthat = (variable == "pthat");
    const Double_t xMin = 0.;
    const Double_t xMax = isPthat ? 150. : 55.;
    const Double_t yMax = isPthat ? 1. : 1.e-1;

    std::unique_ptr<TFile> input(TFile::Open(inputFile, "READ"));
    if (!input || input->IsZombie()) {
        throw std::runtime_error("Cannot open flattened MC file: " +
                                 std::string(inputFile.Data()));
    }

    TTree *tree = nullptr;
    input->GetObject(flatTreeName, tree);
    if (!tree) {
        throw std::runtime_error("Cannot find tree '" + std::string(flatTreeName.Data()) +
                                 "' in " + std::string(inputFile.Data()));
    }
    if (!tree->GetBranch(variable)) {
        throw std::runtime_error("Cannot find comparison variable branch: " +
                                 std::string(variable.Data()));
    }
    if (!tree->GetBranch(weightBranch)) {
        throw std::runtime_error("Cannot find reweighting branch: " +
                                 std::string(weightBranch.Data()));
    }
    if (nBins <= 0) {
        throw std::runtime_error("Invalid comparison histogram binning");
    }

    const SampleDescription sample = DescribeSample(sourceTree, particle, promptness);
    const TString outputDirectory = outputBase + "/" + sample.group;
    if (gSystem->mkdir(outputDirectory, true) != 0 &&
        gSystem->AccessPathName(outputDirectory)) {
        throw std::runtime_error("Cannot create comparison output directory: " +
                                 std::string(outputDirectory.Data()));
    }

    gROOT->cd();
    gStyle->SetOptStat(0);

    const Double_t binWidth = (xMax - xMin) / nBins;
    const TString yAxisTitle = Form("Normalized entries / %.3f", binWidth);
    const TString title = Form(";%s;%s", VariableAxisLabel(variable).Data(),
                               yAxisTitle.Data());
    TH1D unweighted("unweighted", title, nBins, xMin, xMax);
    TH1D reweighted("reweighted", title, nBins, xMin, xMax);
    unweighted.Sumw2();
    reweighted.Sumw2();

    const Long64_t expectedEntries = tree->GetEntries();
    const Long64_t unweightedDrawn =
        tree->Draw(Form("%s>>unweighted", variable.Data()), "", "goff");
    const Long64_t reweightedDrawn =
        tree->Draw(Form("%s>>reweighted", variable.Data()), weightBranch, "goff");
    if (unweightedDrawn != expectedEntries || reweightedDrawn != expectedEntries) {
        throw std::runtime_error(
            "Comparison plot read only " + std::to_string(unweightedDrawn) +
            " unweighted and " + std::to_string(reweightedDrawn) + " weighted entries out of " +
            std::to_string(expectedEntries));
    }

    if (unweighted.Integral() <= 0. || reweighted.Integral() <= 0.) {
        throw std::runtime_error("Empty unweighted or reweighted comparison histogram");
    }
    unweighted.Scale(1. / unweighted.Integral());
    reweighted.Scale(1. / reweighted.Integral());

    // Match the visual language used by plotER/plot_dataMC.C.
    unweighted.SetLineColor(kBlue);
    unweighted.SetFillColor(kBlue);
    unweighted.SetFillStyle(3358);
    reweighted.SetLineColor(kOrange - 2);
    reweighted.SetLineWidth(3);
    reweighted.SetFillStyle(0);
    unweighted.SetMinimum(1.e-5);
    reweighted.SetMinimum(1.e-5);
    unweighted.SetMaximum(yMax);
    reweighted.SetMaximum(yMax);

    TCanvas canvas("canvas", "", 600, 600);
    canvas.SetLeftMargin(0.15);
    canvas.SetTopMargin(0.05);
    canvas.SetRightMargin(0.05);
    canvas.SetLogy();
    unweighted.Draw("HIST");
    reweighted.Draw("HIST SAME");

    TLegend legend(0.18, 0.70, 0.52, 0.93, nullptr, "brNDC");
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.035);
    legend.SetHeader(Form("#bf{%s, %s}", systemName.Data(), sample.plotLabel.Data()));
    legend.AddEntry(&unweighted, "Unweighted MC", "f");
    legend.AddEntry(&reweighted, "pThat-reweighted MC", "l");
    legend.Draw();
    canvas.Update();

    const TString outputStem = outputDirectory + "/" + systemName + "_" +
                               sample.fileTag + "_" + SafeFileToken(variable) + "_" +
                               SafeFileToken(weightBranch) + "_normalized_log";
    const TString pdfFile = outputStem + ".pdf";
    canvas.SaveAs(pdfFile);

    std::cout << "Saved unweighted/reweighted comparison:" << std::endl;
    std::cout << "  " << pdfFile << std::endl;
}
