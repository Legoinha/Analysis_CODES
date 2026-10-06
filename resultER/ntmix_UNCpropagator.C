#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../fitER/aux/uti.h"
#include "aux.h"

// Build the uncertainty of the final X(3872)/psi(2S) corrected-yield ratio.
// Efficiency variations are paired by their stored case name and propagated
// through the ratio before the leading variation is selected.
//
// root -b -q 'ntmix_UNCpropagator.C("ppRef","Bpt")'
// root -b -q 'ntmix_UNCpropagator.C("ppRef","nChargedTracks")'
// root -b -q 'ntmix_UNCpropagator.C("PbPb23","Bpt")'

static TH1D* LoadUncHist(TString path, TString histName, TString cloneName)
{
    TFile* f = TFile::Open(path, "READ");
    TH1D* h = nullptr;
    f->GetObject(histName, h);
    TH1D* out = (TH1D*)h->Clone(cloneName);
    out->SetDirectory(nullptr);
    f->Close();
    delete f;
    return out;
}

struct EffVariationSet {
    TString referenceCase;
    std::vector<TString> cases;
    std::map<std::string, TH1D*> correctedYields;
};

struct RatioSourceUncertainty {
    TH1D* leading = nullptr;
    std::vector<TH1D*> variedRatios;
    std::vector<TH1D*> deviations;
    TString referenceCase;
    std::vector<TString> cases;
};

static EffVariationSet LoadEffVariationSet(TString treename, TString system,
                                           TString var, TString comparisonTag)
{
    const TString path = Form(
        "../effER/output/%s/%s/ROOTs/%s_%s_%s_%s.root",
        treename.Data(), system.Data(), comparisonTag.Data(),
        treename.Data(), system.Data(), var.Data());
    TFile input(path, "READ");

    TNamed* storedReference = nullptr;
    input.GetObject("referenceCase", storedReference);

    EffVariationSet result;
    result.referenceCase = storedReference->GetTitle();
    for (int index = 0; ; ++index) {
        TNamed* storedCase = nullptr;
        input.GetObject(Form("case%d", index), storedCase);
        if (!storedCase) break;

        const TString caseName = storedCase->GetTitle();
        TH1D* storedYield = nullptr;
        input.GetObject(Form("hYieldCorr_%s", caseName.Data()), storedYield);
        TH1D* correctedYield = static_cast<TH1D*>(storedYield->Clone(Form(
            "hYieldCorr_%s_%s", treename.Data(), caseName.Data())));
        correctedYield->SetDirectory(nullptr);
        result.cases.push_back(caseName);
        result.correctedYields[caseName.Data()] = correctedYield;
    }
    input.Close();
    return result;
}

static void DeleteEffVariationSet(EffVariationSet& variations)
{
    for (auto& item : variations.correctedYields) delete item.second;
}

static RatioSourceUncertainty BuildPairedRatioUncertainty(
    TString system, TString var, TString comparisonTag, TString sourceName,
    TString sourceTitle, TH1D* nominalRatio)
{
    EffVariationSet x = LoadEffVariationSet(
        "ntmix_X3872", system, var, comparisonTag);
    EffVariationSet psi = LoadEffVariationSet(
        "ntmix_PSI2S", system, var, comparisonTag);

    RatioSourceUncertainty result;
    result.referenceCase = x.referenceCase;
    result.cases = x.cases;
    result.leading = static_cast<TH1D*>(nominalRatio->Clone(
        Form("h%sUncPercent", sourceName.Data())));
    result.leading->SetDirectory(nullptr);
    result.leading->Reset("ICES");
    result.leading->SetTitle(Form(";%s;%s uncertainty (%%)",
                                  RatioAxisTitle(var).Data(), sourceTitle.Data()));

    for (const TString& caseName : x.cases) {
        if (caseName == x.referenceCase) continue;

        TH1D* variedRatio = static_cast<TH1D*>(
            x.correctedYields.at(caseName.Data())->Clone(Form(
                "h%sRatio_%s", sourceName.Data(), caseName.Data())));
        variedRatio->SetDirectory(nullptr);
        variedRatio->SetTitle(Form(";%s;X(3872) / #psi(2S)",
                                   RatioAxisTitle(var).Data()));
        variedRatio->Divide(psi.correctedYields.at(caseName.Data()));

        TH1D* deviation = static_cast<TH1D*>(nominalRatio->Clone(Form(
            "h%sDeviationPercent_%s", sourceName.Data(), caseName.Data())));
        deviation->SetDirectory(nullptr);
        deviation->Reset("ICES");
        deviation->SetTitle(Form(";%s;|R_{var}/R_{nom}-1| (%%)",
                                 RatioAxisTitle(var).Data()));

        for (int bin = 1; bin <= nominalRatio->GetNbinsX(); ++bin) {
            const double nominal = nominalRatio->GetBinContent(bin);
            const double varied = variedRatio->GetBinContent(bin);
            const double shift = std::abs(varied / nominal - 1.0) * 100.0;
            deviation->SetBinContent(bin, shift);
            result.leading->SetBinContent(
                bin, std::max(result.leading->GetBinContent(bin), shift));
        }
        result.variedRatios.push_back(variedRatio);
        result.deviations.push_back(deviation);
    }

    DeleteEffVariationSet(x);
    DeleteEffVariationSet(psi);
    return result;
}

static void DeleteRatioSourceUncertainty(RatioSourceUncertainty& source)
{
    delete source.leading;
    for (TH1D* h : source.variedRatios) delete h;
    for (TH1D* h : source.deviations) delete h;
}

void ntmix_UNCpropagator(
    TString SYSTEM = "ppRef",
    TString VAR = "Bpt"
) {
    if (VAR == "nMult") VAR = "nChargedTracks";
    gSystem->mkdir("output_ntmix", true);
    gSystem->mkdir("output_ntmix/systematicFILES", true);
    gStyle->SetOptStat(0);

    TH1D* hFitX = LoadUncHist(Form(
        "../fitER/ROOTfiles/%s/systematicFILES/fit_totalSystematic_leadingUnc_ntmix_X3872_%s_%s.root",
        SYSTEM.Data(), SYSTEM.Data(), VAR.Data()),
        "hLeadingTotalUncPercent", "hFitModelX3872UncPercent");
    TH1D* hFitPsi = LoadUncHist(Form(
        "../fitER/ROOTfiles/%s/systematicFILES/fit_totalSystematic_leadingUnc_ntmix_PSI2S_%s_%s.root",
        SYSTEM.Data(), SYSTEM.Data(), VAR.Data()),
        "hLeadingTotalUncPercent", "hFitModelPSI2SUncPercent");

    TH1D* hNominalRatio = BuildNominalRatio(SYSTEM, VAR);
    hNominalRatio->SetName("hNominalRatio");

    RatioSourceUncertainty dataMC = BuildPairedRatioUncertainty(
        SYSTEM, VAR, "DATA_MC_AGREEMENT_reweight_comparison",
        "DataMC", "Data-MC discrepancy", hNominalRatio);
    RatioSourceUncertainty method = BuildPairedRatioUncertainty(
        SYSTEM, VAR, "METHODS_comparison",
        "AccEffMethod", "Acc#timesEff method", hNominalRatio);

    TH1D* hFitRatio = static_cast<TH1D*>(hFitX->Clone("hFitModelUncPercent"));
    hFitRatio->SetDirectory(nullptr);
    hFitRatio->Reset("ICES");
    hFitRatio->SetTitle(Form(";%s;Fit-model uncertainty in ratio (%%)",
                             RatioAxisTitle(VAR).Data()));

    TH1D* hTotal = static_cast<TH1D*>(hFitRatio->Clone("hTotalUncPercent"));
    hTotal->SetDirectory(nullptr);
    hTotal->SetTitle(Form(";%s;Total systematic uncertainty (%%)",
                          RatioAxisTitle(VAR).Data()));

    std::vector<std::vector<double>> tableNumbers(hTotal->GetNbinsX());
    for (int bin = 1; bin <= hTotal->GetNbinsX(); ++bin) {
        const double fitX = hFitX->GetBinContent(bin);
        const double fitPsi = hFitPsi->GetBinContent(bin);
        const double fitRatio = std::hypot(fitX, fitPsi);
        const double mapRatio = dataMC.leading->GetBinContent(bin);
        const double methodRatio = method.leading->GetBinContent(bin);
        const double total = std::sqrt(fitRatio * fitRatio
                                     + mapRatio * mapRatio
                                     + methodRatio * methodRatio);
        hFitRatio->SetBinContent(bin, fitRatio);
        hTotal->SetBinContent(bin, total);
        tableNumbers[bin - 1] = {
            fitX, fitPsi, fitRatio, mapRatio, methodRatio, total};
    }

    const TString stem = Form(
        "ntmix_ratioUnc_X3872_over_PSI2S_%s_%s",
        SYSTEM.Data(), VAR.Data());
    TFile* output = new TFile(Form(
        "output_ntmix/systematicFILES/%s.root", stem.Data()), "RECREATE");
    hNominalRatio->Write();
    hFitX->Write();
    hFitPsi->Write();
    hFitRatio->Write();
    dataMC.leading->Write("hDataMCUncPercent");
    method.leading->Write("hAccEffMethodUncPercent");
    hTotal->Write();
    for (TH1D* h : dataMC.variedRatios) h->Write();
    for (TH1D* h : dataMC.deviations) h->Write();
    for (TH1D* h : method.variedRatios) h->Write();
    for (TH1D* h : method.deviations) h->Write();

    TNamed("uncertaintyTarget", "corrected-yield ratio X3872/PSI2S").Write();
    TNamed("system", SYSTEM.Data()).Write();
    TNamed("analysisVariable", VAR.Data()).Write();
    TNamed("numeratorTree", "ntmix_X3872").Write();
    TNamed("denominatorTree", "ntmix_PSI2S").Write();
    TNamed("fitModelTreatment", "independent X3872 and PSI2S uncertainties").Write();
    TNamed("dataMCTreatment", "paired corrected-yield ratio variation").Write();
    TNamed("accEffMethodTreatment", "paired corrected-yield ratio variation").Write();
    TNamed("dataMCReferenceCase", dataMC.referenceCase.Data()).Write();
    TNamed("accEffMethodReferenceCase", method.referenceCase.Data()).Write();
    for (std::size_t index = 0; index < dataMC.cases.size(); ++index) {
        TNamed(Form("dataMCCase%zu", index), dataMC.cases[index].Data()).Write();
    }
    for (std::size_t index = 0; index < method.cases.size(); ++index) {
        TNamed(Form("accEffMethodCase%zu", index), method.cases[index].Data()).Write();
    }
    output->Close();

    std::vector<std::string> columnNames = {"Systematic source"};
    for (int bin = 1; bin <= hTotal->GetNbinsX(); ++bin) {
        columnNames.push_back(GetSystematicColumnLabel(
            VAR, hTotal->GetXaxis()->GetBinLowEdge(bin),
            hTotal->GetXaxis()->GetBinUpEdge(bin)));
    }
    const std::vector<std::string> rowLabels = {
        "Fit model X(3872)",
        "Fit model $\\psi(2S)$",
        "Fit model ratio (independent)",
        "Data-MC discrepancy (paired ratio)",
        "Acc$\\times$Eff method (paired ratio)",
        "Total"
    };
    latex_tables_document(
        Form("output_ntmix/systematicFILES/%s_table", stem.Data()),
        {static_cast<int>(columnNames.size())},
        {static_cast<int>(rowLabels.size() + 1)},
        {columnNames}, {rowLabels}, {tableNumbers});

    std::cout << "[ntmix_UNCpropagator] wrote ratio-level systematics to "
              << "output_ntmix/systematicFILES/" << stem << ".root"
              << std::endl;

    delete hFitX;
    delete hFitPsi;
    delete hFitRatio;
    delete hNominalRatio;
    delete hTotal;
    DeleteRatioSourceUncertainty(dataMC);
    DeleteRatioSourceUncertainty(method);
    delete output;
}
