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
    return {"Bzero", "B^{0}", "Bmeson"};
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
    TTree *tree = input->Get<TTree>(flatTreeName);

    const SampleDescription sample = DescribeSample(sourceTree, particle, promptness);
    const TString outputDirectory = outputBase + "/" + sample.group;
    gSystem->mkdir(outputDirectory, true);

    gROOT->cd();
    gStyle->SetOptStat(0);

    const Double_t binWidth = (xMax - xMin) / nBins;
    const TString yAxisTitle = Form("Normalized entries / %.3f", binWidth);
    const TString axisLabel = variable == "Bpt" ? "p_{T} [GeV/c]"
                            : variable == "pthat" ? "Generator #hat{p}_{T} [GeV/c]"
                            : variable == "CentBin" ? "Centrality bin"
                            : variable.Contains("nChargedTracks") ? "Track multiplicity" : variable.Data();
    const TString title = Form(";%s;%s", axisLabel.Data(), yAxisTitle.Data());
    TH1D unweighted("unweighted", title, nBins, xMin, xMax);
    TH1D reweighted("reweighted", title, nBins, xMin, xMax);
    unweighted.Sumw2();
    reweighted.Sumw2();

    tree->Draw(Form("%s>>unweighted", variable.Data()), "", "goff");
    tree->Draw(Form("%s>>reweighted", variable.Data()), weightBranch, "goff");
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
                               sample.fileTag + "_" + variable + "_" + weightBranch + "_normalized_log";
    const TString pdfFile = outputStem + ".pdf";
    canvas.SaveAs(pdfFile);

    std::cout << "Saved unweighted/reweighted comparison:" << std::endl;
    std::cout << "  " << pdfFile << std::endl;
}
