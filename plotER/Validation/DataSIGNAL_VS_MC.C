#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTreeFormula.h>
#include <TParameter.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TPRegexp.h>
#include <TLatex.h>

#include <RooRealVar.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooWorkspace.h>
#include <RooAbsPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooStats/SPlot.h>

#include <iostream>
#include <vector>
#include <cmath>
#include <memory>
#include <algorithm>
#include <stdexcept>

#include "aux.h"
#include "../../fitER/aux/uti.h"

using namespace RooFit;

struct ValidationHist {
    VarCfgSignal var;
    TString tag;
    TString baseVar;
    TH1D* backgroundSideband;
    TH1D* backgroundSplot;
    TH1D* sideband;
    TH1D* splot;
    TH1D* mc;
};

static void fillMCHistFromTree(TTree* tree, TH1D* hist, const TString& expr, const TString& cut,
                               const ReweightInput* reweightInput)
{
    TTreeFormula exprFormula("exprFormula", expr.Data(), tree);
    TTreeFormula cutFormula("cutFormula", cut.Data(), tree);
    TTreeFormula pThatFormula("pThatFormula", "pThatreweight", tree);
    std::unique_ptr<TTreeFormula> reweightFormula;
    if (reweightInput) {
        reweightFormula.reset(new TTreeFormula("reweightFormula", reweightInput->expr.Data(), tree));
    }

    Int_t currentTree = -1;
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->LoadTree(i);
        tree->GetEntry(i);
        if (tree->GetTreeNumber() != currentTree) {
            currentTree = tree->GetTreeNumber();
            exprFormula.UpdateFormulaLeaves();
            cutFormula.UpdateFormulaLeaves();
            pThatFormula.UpdateFormulaLeaves();
            if (reweightFormula) reweightFormula->UpdateFormulaLeaves();
        }

        exprFormula.GetNdata();
        cutFormula.GetNdata();
        pThatFormula.GetNdata();
        if (reweightFormula) reweightFormula->GetNdata();
        if (cutFormula.EvalInstance() == 0.0) continue;

        double weight = pThatFormula.EvalInstance();
        if (reweightFormula) {
            weight *= lookupWeight1D(reweightInput->hist, reweightFormula->EvalInstance());
        }
        hist->Fill(exprFormula.EvalInstance(), weight);
    }
}

void DataSIGNAL_VS_MC(
    TString modelPath,
    TString treeName,
    TString dataTree,
    TString massAxisTitle,
    TString weightPath,
    TString weightToApply,
    TString selectedWeightParticleTag,
    TString systemTag,
    bool showBackgroundProjections = true)
{
    std::cout << "[validation] Fit model: " << modelPath << std::endl;
    const bool REWEIGHT_MC = !weightToApply.IsNull();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    const bool isNtKp = (treeName == "ntKp");

    TString compareBaseDir = Form("Compare_%s", systemTag.Data());
    TString treeOutDir = Form("%s/%s", compareBaseDir.Data(), treeName.Data());
    TString outDir = treeOutDir;
    if (REWEIGHT_MC) outDir = Form("%s/reweightMC_comparison", treeOutDir.Data());
    gSystem->mkdir(compareBaseDir, true);
    gSystem->mkdir(treeOutDir, true);
    gSystem->mkdir(outDir, true);

    TFile* fModel = TFile::Open(modelPath.Data(), "READ");
    const TString dataPath =
        static_cast<TObjString*>(fModel->Get("inputData"))->GetString();
    const TString mcPath =
        static_cast<TObjString*>(fModel->Get("inputMC"))->GetString();
    const TString fitSelection =
        static_cast<TObjString*>(fModel->Get("selectionCut"))->GetString();
    const TString fitVariable =
        static_cast<TObjString*>(fModel->Get("variableName"))->GetString();
    TH1D* analysisVarBins =
        static_cast<TH1D*>(fModel->Get("analysisVarBins"));
    const double fitVariableMin = analysisVarBins->GetXaxis()->GetBinLowEdge(1);
    const double fitVariableMax =
        analysisVarBins->GetXaxis()->GetBinUpEdge(analysisVarBins->GetNbinsX());
    const TString baseCut = Form(
        "(%s) && (abs(%s)>=%0.17g && abs(%s)<%0.17g)",
        fitSelection.Data(), fitVariable.Data(), fitVariableMin,
        fitVariable.Data(), fitVariableMax);

    std::cout << "[validation] Fit metadata input DATA: " << dataPath << std::endl;
    std::cout << "[validation] Fit metadata input MC:   " << mcPath << std::endl;
    std::cout << "[validation] Fit metadata selection: " << baseCut << std::endl;
    std::cout << "[validation] Fit variable range: " << fitVariable << " in ["
              << fitVariableMin << ", " << fitVariableMax << ")" << std::endl;

    TFile* fData = TFile::Open(dataPath.Data(), "READ");
    TFile* fMC = TFile::Open(mcPath.Data(), "READ");
    auto closeInputFiles = [&]() {
        fData->Close(); delete fData; fData = nullptr;
        fMC->Close(); delete fMC; fMC = nullptr;
        fModel->Close(); delete fModel; fModel = nullptr;
    };


    TTree* tData = nullptr;
    TTree* tMC = nullptr;
    fData->GetObject(dataTree.Data(), tData);
    fMC->GetObject(treeName.Data(), tMC);
    const Long64_t invalidPThatWeights = tMC->Draw(
        "pThatreweight",
        "pThatreweight <= 0 || !TMath::Finite(pThatreweight)",
        "goff");
    if (invalidPThatWeights > 0) {
        throw std::runtime_error(Form("[validation] Found %lld non-finite or nonpositive pThatreweight entries in %s",
                                      invalidPThatWeights, mcPath.Data()));
    }
    const double minPThatWeight = tMC->GetMinimum("pThatreweight");
    const double maxPThatWeight = tMC->GetMaximum("pThatreweight");
    std::cout << "[validation] Applying pThatreweight to every selected MC entry; range = ["
              << minPThatWeight << ", " << maxPThatWeight << "]" << std::endl;

    RooWorkspace* ws = nullptr;
    fModel->GetObject("ws_nominal", ws);

    RooRealVar* meanVar = ws->var("mean1_");
    RooRealVar* sigma1 = ws->var("sigma11_");
    RooRealVar* sigma2 = ws->var("sigma21_");
    RooRealVar* sig1frac = ws->var("sig1frac1_");
    RooRealVar* scale = ws->var("scale");
    RooAbsPdf* model = ws->pdf("model1_");
    RooRealVar* nsig = ws->var("nsig1_");
    RooRealVar* nbkg = ws->var("nbkg1_");
    RooRealVar* nbkgPartR = isNtKp ? ws->var("nbkg_part_r1_") : nullptr;

    TParameter<int>* nMassBinsPar = nullptr;
    TParameter<double>* massMinPar = nullptr;
    TParameter<double>* massMaxPar = nullptr;
    fModel->GetObject("nMassBins", nMassBinsPar);
    fModel->GetObject("massMin", massMinPar);
    fModel->GetObject("massMax", massMaxPar);

    const int nMassBins = nMassBinsPar->GetVal();
    const double massMin = massMinPar->GetVal();
    const double massMax = massMaxPar->GetVal();
    double massPlotMin = massMin;
    double massPlotMax = massMax;
    if (treeName == "ntphi") {
        massPlotMin = 5.20;
        massPlotMax = 5.55;
    } else if (treeName == "ntKp" || treeName == "ntKstar") {
        massPlotMin = 5.05;
        massPlotMax = 5.55;
    }

    ReweightInput reweightInput;
    TString reweightTag = "";
    if (!REWEIGHT_MC) {
        std::cout << "Running nominal validation; every available comparison weight will be saved to "
                  << weightPath << std::endl;
    }

    const double signalNSigma = 2;
    const double sidebandInNSigma = 4.0;
    const double sidebandOutNSigma = 8.0;
    const double sigWeight = sig1frac->getVal();
    const double mean = meanVar->getVal();
    const double sigma = sqrt(sigWeight * pow(sigma1->getVal(), 2) + (1.0 - sigWeight) * pow(sigma2->getVal(), 2)) * scale->getVal();
    const double sigLo = mean - signalNSigma * sigma;
    const double sigHi = mean + signalNSigma * sigma;
    const double sbLLo = mean - sidebandOutNSigma * sigma;
    const double sbLHi = mean - sidebandInNSigma * sigma;
    const double sbRLo = mean + sidebandInNSigma * sigma;
    const double sbRHi = mean + sidebandOutNSigma * sigma;
    const double sidebandWidth = isNtKp ? (sbRHi - sbRLo) : ((sbLHi - sbLLo) + (sbRHi - sbRLo));
    const double alpha = (sidebandWidth > 0.0) ? (sigHi - sigLo) / sidebandWidth : 0.0;

    std::cout << "[sideband] " << treeName << ": mean=" << mean << ", effSigma=" << sigma << std::endl;
    std::cout << "  Signal region (+/-" << signalNSigma << "sigma): [" << sigLo << ", " << sigHi << "]" << std::endl;
    std::cout << "  Sideband region (" << sidebandInNSigma << "-" << sidebandOutNSigma << "sigma): ["
              << sbLLo << ", " << sbLHi << "] U [" << sbRLo << ", " << sbRHi << "]" << std::endl;

    TH1D* hMassWin = new TH1D("hMassWin_tmp", Form(";%s;", massAxisTitle.Data()), nMassBins, massMin, massMax);
    const double massBinWidthMeV = 1000.0 * (massMax - massMin) / nMassBins;
    hMassWin->GetYaxis()->SetTitle(Form("Entries / %.4g MeV/c^{2}", massBinWidthMeV));
    hMassWin->GetXaxis()->SetRangeUser(massPlotMin, massPlotMax);
    tData->Draw("Bmass>>hMassWin_tmp", baseCut, "goff");

    TCanvas* cWin = new TCanvas(Form("cWin_%s", treeName.Data()), "", 760, 650);
    cWin->SetLeftMargin(0.14);
    hMassWin->SetLineColor(kBlack);
    hMassWin->SetMarkerStyle(20);
    hMassWin->SetMarkerSize(0.7);
    hMassWin->Draw("E");

    const double yMaxWin = hMassWin->GetMaximum() * 1.15;
    hMassWin->SetMinimum(0.0);
    hMassWin->SetMaximum(yMaxWin);

    TBox* bSig = new TBox(sigLo, 0.0, sigHi, yMaxWin);
    bSig->SetFillColorAlpha(kRed + 1, 0.20);
    bSig->SetLineColor(kRed + 1);
    bSig->Draw("SAME");

    TBox* bSBL = nullptr;
    if (!isNtKp) {
        bSBL = new TBox(sbLLo, 0.0, sbLHi, yMaxWin);
        bSBL->SetFillColorAlpha(kBlue + 1, 0.18);
        bSBL->SetLineColor(kBlue + 1);
        bSBL->Draw("SAME");
    }

    TBox* bSBR = new TBox(sbRLo, 0.0, sbRHi, yMaxWin);
    bSBR->SetFillColorAlpha(kBlue + 1, 0.18);
    bSBR->SetLineColor(kBlue + 1);
    bSBR->Draw("SAME");
    hMassWin->Draw("E SAME");

    TLegend* legWin = new TLegend(0.62, 0.70, 0.90, 0.90);
    legWin->SetBorderSize(0);
    legWin->SetFillStyle(0);
    legWin->AddEntry(hMassWin, "Data", "lep");
    legWin->AddEntry(bSig, Form("Signal region (#pm%g#sigma)", signalNSigma), "f");
    legWin->AddEntry(isNtKp ? bSBR : bSBL, Form("Sideband region (%g-%g#sigma)", sidebandInNSigma, sidebandOutNSigma), "f");
    TLatex winLabel;
    winLabel.SetNDC();
    winLabel.SetTextFont(42);
    winLabel.SetTextSize(0.060);
    winLabel.DrawLatex(0.18, 0.84, FitParticleLabel(treeName, true));
    legWin->Draw();
    cWin->SaveAs(Form("%s/mass_windows_%s_%s.pdf", treeOutDir.Data(), treeName.Data(), systemTag.Data()));

    delete legWin;
    delete bSig;
    delete bSBL;
    delete bSBR;
    delete cWin;
    delete hMassWin;

    auto vars = getSignalVars(treeName, systemTag, fitVariable,
                              fitVariableMin, fitVariableMax);

    std::vector<VarCfgSignal> availableVars;
    for (const auto& var : vars) {
        const bool missingInData = !tData->GetBranch(var.expr);
        const bool missingInMC = !tMC->GetBranch(var.expr);
        if (missingInData || missingInMC) {
            std::cerr << "[validation] WARNING: skipping comparison " << var.expr << "; branch missing from"
                      << (missingInData ? " DATA" : "")
                      << (missingInMC ? " MC" : "") << "." << std::endl;
            continue;
        }
        availableVars.push_back(var);
    }
    vars.swap(availableVars);

    // The classifier score is drawn from its cut in the fit selection ("Prediction > <threshold>") up to 1.
    for (auto& var : vars) {
        if (var.expr == "Prediction") {
            std::unique_ptr<TObjArray> threshold(TPRegexp("Prediction\\s*>=?\\s*([0-9.eE+-]+)").MatchS(baseCut));
            var.xmin = static_cast<TObjString*>(threshold->At(1))->GetString().Atof();
            var.xmax = 1.0;
        }
    }

    if (REWEIGHT_MC) {
        const auto reweightCfg = std::find_if(vars.begin(), vars.end(), [&](const VarCfgSignal& var) {
            return weightToApply == var.expr;
        });
        reweightTag = reweightCfg->expr;
        const TString reweightFolderTag = Form("%s__weight%s", reweightTag.Data(),
                                                selectedWeightParticleTag.Data());
        outDir = Form("%s/reweightMC_comparison/%s", treeOutDir.Data(),
                      reweightFolderTag.Data());
        gSystem->mkdir(outDir, true);

        std::cout << "Using hWeight_" << reweightTag << " from the "
                  << selectedWeightParticleTag << " all-comparisons weight file" << std::endl;
        std::unique_ptr<TFile> fWeight(TFile::Open(Form("file:%s", weightPath.Data()), "READ"));
        const TString histName = Form("hWeight_%s", reweightTag.Data());
        TH1D* sourceWeight = nullptr;
        fWeight->GetObject(histName.Data(), sourceWeight);
        TH1D* clonedWeight = (TH1D*)sourceWeight->Clone(
            Form("hWeight_runtime_%s", reweightTag.Data()));
        clonedWeight->SetDirectory(nullptr);
        reweightInput.expr = validationWeightExpression(*reweightCfg);
        reweightInput.tag = reweightTag;
        reweightInput.hist = clonedWeight;
        std::cout << "Running validation with " << histName
                  << " from " << weightPath << std::endl;
    }

    std::vector<ValidationHist> hists;
    TString cutSig = Form("(%s) && (Bmass>%f && Bmass<%f)", baseCut.Data(), sigLo, sigHi);
    TString cutSB = isNtKp ? Form("(%s) && (Bmass>%f && Bmass<%f)", baseCut.Data(), sbRLo, sbRHi)
                           : Form("(%s) && ((Bmass>%f && Bmass<%f) || (Bmass>%f && Bmass<%f))", baseCut.Data(), sbLLo, sbLHi, sbRLo, sbRHi);
    TString cutMC = Form("(%s)", baseCut.Data());

    for (const auto& v : vars) {
        TString expr = v.absVal ? Form("abs(%s)", v.expr.Data()) : v.expr;
        TString tag = v.expr;

        TH1D* hDataSR = new TH1D(Form("hDataSR_%s_%s", tag.Data(), treeName.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hDataSB = new TH1D(Form("hDataSB_%s_%s", tag.Data(), treeName.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hSideband = new TH1D(Form("hSideband_%s_%s", tag.Data(), treeName.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hSPlot = new TH1D(Form("hSPlot_%s_%s", tag.Data(), treeName.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hBackgroundSplot = new TH1D(Form("hBackgroundSplot_%s_%s", tag.Data(), treeName.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hMC = new TH1D(Form("hMC_%s_%s", tag.Data(), treeName.Data()), v.title, v.nbins, v.xmin, v.xmax);
        hDataSR->Sumw2();
        hDataSB->Sumw2();
        hSideband->Sumw2();
        hSPlot->Sumw2();
        hBackgroundSplot->Sumw2();
        hMC->Sumw2();

        tData->Draw(Form("%s>>%s", expr.Data(), hDataSR->GetName()), cutSig, "goff");
        tData->Draw(Form("%s>>%s", expr.Data(), hDataSB->GetName()), cutSB, "goff");
        TH1D* hBackgroundSideband = (TH1D*)hDataSB->Clone(Form("hBackgroundSideband_%s_%s", tag.Data(), treeName.Data()));
        hSideband->Add(hDataSR);
        hSideband->Add(hDataSB, -alpha);

        fillMCHistFromTree(tMC, hMC, expr, cutMC,
                           REWEIGHT_MC ? &reweightInput : nullptr);

        if (hBackgroundSideband->Integral() > 0) hBackgroundSideband->Scale(1.0 / hBackgroundSideband->Integral());
        if (hSideband->Integral() > 0) hSideband->Scale(1.0 / hSideband->Integral());
        if (hMC->Integral() > 0) hMC->Scale(1.0 / hMC->Integral());

        hists.push_back({v, tag, v.expr, hBackgroundSideband, hBackgroundSplot, hSideband, hSPlot, hMC});
        delete hDataSR;
        delete hDataSB;
    }

    RooRealVar Bmass("Bmass", "Bmass", massMin, massMax);
    RooArgSet obs(Bmass);
    std::vector<std::unique_ptr<RooRealVar>> extraObs;
    auto addObsIfAvailable = [&](const TString& name) {
        if (name != "Bmass" && !obs.find(name.Data())) {
            extraObs.emplace_back(new RooRealVar(name, name, -1e6, 1e6));
            obs.add(*extraObs.back());
        }
    };
    for (const auto& h : hists) {
        addObsIfAvailable(h.baseVar);
    }

    // Every DATA branch the selection reads (pre-cut, score) joins the observables: the RooDataSet applies the cut on them.
    TTreeFormula selectionFormula("selectionFormula", baseCut, tData);
    for (int i = 0; i < selectionFormula.GetNcodes(); ++i) {
        addObsIfAvailable(selectionFormula.GetLeaf(i)->GetBranch()->GetName());
    }
    TString dataCut = Form("(%s) && (Bmass>%f && Bmass<%f)", baseCut.Data(), massMin, massMax);
    RooDataSet data("data", "data", tData, obs, dataCut.Data());
    std::vector<RooRealVar*> fitYields = {nsig, nbkg};
    if (isNtKp) fitYields.push_back(nbkgPartR);
    const double yieldMax = std::max(1.0, static_cast<double>(data.numEntries()));
    for (RooRealVar* y : fitYields) y->setRange(0.0, yieldMax);

    RooArgSet* allPars = model->getParameters(data);
    std::unique_ptr<TIterator> it(allPars->createIterator());
    TObject* obj = nullptr;
    while ((obj = it->Next())) {
        RooRealVar* v = dynamic_cast<RooRealVar*>(obj);
        v->setConstant(true);
    }
    nsig->setConstant(false);
    nbkg->setConstant(false);
    if (isNtKp) nbkgPartR->setConstant(false);

    RooFitResult* fitRes = model->fitTo(data, Extended(true), Save(true), PrintLevel(-1));
    std::cout << "[sPlot] yield-only fit: status=" << fitRes->status()
              << ", covQual=" << fitRes->covQual()
              << ", edm=" << fitRes->edm() << std::endl;
    if (fitRes->status() != 0 || fitRes->covQual() < 2) {
        std::cerr << "[sPlot] Yield-only fit is not reliable; refusing to build sWeights."
                  << std::endl;
        for (auto& h : hists) {
            delete h.backgroundSideband;
            delete h.backgroundSplot;
            delete h.sideband;
            delete h.splot;
            delete h.mc;
        }
        delete reweightInput.hist;
        delete fitRes;
        delete allPars;
        closeInputFiles();
        return;
    }
    RooArgList splotYields;
    splotYields.add(*nsig);
    if (isNtKp) splotYields.add(*nbkgPartR);
    splotYields.add(*nbkg);
    RooStats::SPlot sData("sData", "An SPlot", data, model, splotYields);

    double swMin = 1e9, swMax = -1e9, swMean = 0.0;
    int negWeights = 0, totalWeights = 0;
    TString swVarName = Form("%s_sw", nsig->GetName());
    TString backgroundSwVarName = Form("%s_sw", nbkg->GetName());
    TString partialRecoSwVarName = isNtKp ? Form("%s_sw", nbkgPartR->GetName()) : "";
    for (int i = 0; i < data.numEntries(); ++i) {
        const RooArgSet* row = data.get(i);
        double w = row->getRealValue(swVarName.Data());
        if (w < 0) negWeights++;
        swMin = std::min(swMin, w);
        swMax = std::max(swMax, w);
        swMean += w;
        totalWeights++;
    }
    if (totalWeights > 0) swMean /= totalWeights;
    std::cout << "[sPlot] sWeight stats: min=" << swMin << ", max=" << swMax << ", mean=" << swMean
              << ", negative=" << negWeights << "/" << totalWeights << std::endl;

    TCanvas* cMass = new TCanvas(Form("cMass_%s", treeName.Data()), "", 760, 680);
    cMass->SetLeftMargin(0.14);
    cMass->SetRightMargin(0.04);
    cMass->SetBottomMargin(0.13);
    Bmass.setBins(nMassBins);
    RooPlot* frame = Bmass.frame(Range(massPlotMin, massPlotMax));
    frame->SetTitle("");
    frame->SetStats(0);
    frame->GetXaxis()->SetTitle(massAxisTitle);
    frame->GetXaxis()->SetTitleSize(0.030);
    frame->GetXaxis()->SetTitleOffset(1.25);
    frame->GetXaxis()->CenterTitle();
    frame->GetXaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetLabelFont(42);
    frame->GetXaxis()->SetLabelOffset(0.012);
    frame->GetXaxis()->SetLabelSize(0.031);
    frame->GetXaxis()->SetTickLength(0.035);
    frame->GetYaxis()->SetTitle("Events");
    frame->GetYaxis()->SetTitleOffset(1.65);
    frame->GetYaxis()->SetTitleSize(0.035);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetLabelFont(42);
    frame->GetYaxis()->SetLabelSize(0.035);
    frame->SetMinimum(0.0);
    const int signalColor = (treeName == "ntmix_PSI2S") ? kOrange - 2 : kOrange - 3;
    data.plotOn(frame, Name("data_splot_fit"), Binning(nMassBins), MarkerSize(0.5), MarkerStyle(8), MarkerColor(kBlack), LineColor(kBlack), LineWidth(1));
    std::unique_ptr<RooArgSet> components(model->getComponents());
    RooAbsPdf* sigPdf = dynamic_cast<RooAbsPdf*>(components->find("sig_doubleG1_"));
    RooAbsPdf* bkgPdf = dynamic_cast<RooAbsPdf*>(components->find("bkg1_"));
    RooAbsPdf* partPdf = nullptr;
    if (isNtKp) partPdf = dynamic_cast<RooAbsPdf*>(components->find("erfc1"));

    model->plotOn(frame, Name("signal_splot_fit"), Components(*sigPdf), DrawOption("LF"), FillStyle(3002), FillColor(signalColor), LineStyle(7), LineColor(signalColor), LineWidth(1), Precision(1e-6));
    if (isNtKp) model->plotOn(frame, Name("partial_reco_splot_fit"), Components(*partPdf), DrawOption("L"), LineStyle(9), LineColor(kGreen + 3), LineWidth(2), Precision(1e-6));
    model->plotOn(frame, Name("model_splot_fit"), Precision(1e-6), DrawOption("L"), LineColor(kRed), LineWidth(1));
    model->plotOn(frame, Name("background_splot_fit"), Components(*bkgPdf), Precision(1e-6), DrawOption("L"), LineStyle(7), LineColor(kBlue + 1), LineWidth(2));
    frame->getAttFill()->SetFillStyle(0);
    frame->Draw();

    TLegend* legMass = new TLegend(0.62, 0.66, 0.91, 0.90);
    legMass->SetBorderSize(0);
    legMass->SetFillStyle(0);
    legMass->SetTextFont(42);
    legMass->SetTextSize(0.035);
    legMass->AddEntry(frame->findObject("data_splot_fit"), "Data", "lep");
    legMass->AddEntry(frame->findObject("model_splot_fit"), "Fit Model", "l");
    legMass->AddEntry(frame->findObject("background_splot_fit"), "Comb. Bkg.", "l");
    legMass->AddEntry(frame->findObject("signal_splot_fit"), FitParticleLabel(treeName, true), "f");
    if (isNtKp) legMass->AddEntry(frame->findObject("partial_reco_splot_fit"), "Partial reco.", "l");
    legMass->Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(42);
    label.SetTextSize(0.036);
    label.DrawLatex(0.16, 0.86, FitParticleLabel(treeName, true));
    label.SetTextSize(0.030);
    label.DrawLatex(0.16, 0.80, Form("N_{sig} = %.1f #pm %.1f", nsig->getVal(), nsig->getError()));
    label.DrawLatex(0.16, 0.75, Form("N_{bkg} = %.1f #pm %.1f", nbkg->getVal(), nbkg->getError()));
    if (isNtKp) label.DrawLatex(0.16, 0.70, Form("N_{part} = %.1f #pm %.1f", nbkgPartR->getVal(), nbkgPartR->getError()));

    cMass->SaveAs(Form("%s/massFit_splot_%s_%s.pdf", treeOutDir.Data(), treeName.Data(), systemTag.Data()));
    delete legMass;
    delete frame;
    delete cMass;

    for (auto& h : hists) {
        for (int i = 0; i < data.numEntries(); ++i) {
            const RooArgSet* row = data.get(i);
            double val = row->getRealValue(h.baseVar.Data());
            if (h.var.absVal) val = std::abs(val);
            double w = row->getRealValue(swVarName.Data());
            h.splot->Fill(val, w);
            double backgroundWeight = row->getRealValue(backgroundSwVarName.Data());
            if (isNtKp) backgroundWeight += row->getRealValue(partialRecoSwVarName.Data());
            h.backgroundSplot->Fill(val, backgroundWeight);
        }
        if (h.splot->Integral() > 0) h.splot->Scale(1.0 / h.splot->Integral());
        if (h.backgroundSplot->Integral() > 0) h.backgroundSplot->Scale(1.0 / h.backgroundSplot->Integral());
    }

    TFile* fWeights = nullptr;
    TString weightOutputPath = weightPath;
    if (!REWEIGHT_MC) {
        gSystem->mkdir("WEIGHTS", true);
        fWeights = TFile::Open(Form("file:%s", weightOutputPath.Data()), "RECREATE");
        TString weightVariables;
        for (const auto& h : hists) {
            if (!weightVariables.IsNull()) weightVariables += ",";
            weightVariables += h.tag;
        }
        fWeights->cd();
        TNamed("weightFileContract", "All available nominal sPlot/MC comparison weights; select hWeight_<variable>").Write();
        TNamed("weightVariables", weightVariables.Data()).Write();
        TNamed("fitModel", modelPath.Data()).Write();
    }

    for (auto& h : hists) {

        h.mc->SetLineColor(kOrange + 7);
        h.mc->SetLineWidth(2);
        h.sideband->SetLineColor(kBlue + 1);
        h.sideband->SetMarkerColor(kBlue + 1);
        h.sideband->SetMarkerStyle(20);
        h.sideband->SetLineWidth(2);
        h.splot->SetLineColor(kRed + 1);
        h.splot->SetMarkerColor(kRed + 1);
        h.splot->SetMarkerStyle(24);
        h.splot->SetLineWidth(2);
        h.backgroundSideband->SetLineColor(kGray + 1);
        h.backgroundSideband->SetMarkerColor(kGray + 1);
        h.backgroundSideband->SetMarkerStyle(25);
        h.backgroundSideband->SetLineStyle(1);
        h.backgroundSideband->SetLineWidth(2);
        h.backgroundSplot->SetLineColor(kGray + 1);
        h.backgroundSplot->SetMarkerColor(kGray + 1);
        h.backgroundSplot->SetMarkerStyle(26);
        h.backgroundSplot->SetLineStyle(2);
        h.backgroundSplot->SetLineWidth(2);

        TH1D* hWeight = nullptr;
        if (!REWEIGHT_MC) {
            hWeight = makeWeightHist(h.splot, h.mc, Form("hWeight_%s", h.tag.Data()));
            fWeights->cd();
            hWeight->Write();
            TNamed(Form("weightExpression_%s", h.tag.Data()),
                   validationWeightExpression(h.var).Data()).Write();
        }

        TH1D* hMCBand = (TH1D*)h.mc->Clone(Form("hMCBand_cmp_%s", h.tag.Data()));
        hMCBand->SetFillColorAlpha(kOrange + 7, 0.30);
        hMCBand->SetFillStyle(1001);
        hMCBand->SetLineColor(kOrange + 7);
        hMCBand->SetLineWidth(1);
        hMCBand->SetMarkerSize(0);

        const AgreementMetrics agreement = computeAgreementMetrics1D(h.splot, h.mc);
        const double reducedChi2 = (agreement.ndf > 0) ? agreement.chi2 / agreement.ndf : -1.0;

        TH1D* hRatioSP = makeWeightHist(h.splot, h.mc, Form("hDataOverMCSP_%s", h.tag.Data()));
        hRatioSP->SetLineColor(kRed + 1);
        hRatioSP->SetMarkerColor(kRed + 1);
        hRatioSP->SetMarkerStyle(20);
        hRatioSP->SetMarkerSize(0.9);
        hRatioSP->SetLineWidth(2);

        double ymax = std::max(h.mc->GetMaximum(), std::max(h.sideband->GetMaximum(), h.splot->GetMaximum()));
        if (showBackgroundProjections) {
            ymax = std::max(ymax, h.backgroundSideband->GetMaximum());
            ymax = std::max(ymax, h.backgroundSplot->GetMaximum());
        }
        h.sideband->SetMaximum(1.35 * ymax);
        h.sideband->SetMinimum(0.0);

        double rMin = hRatioSP->GetMinimum();
        double rMax = hRatioSP->GetMaximum();
        if (!(rMax > rMin)) {
            rMin = 0.0;
            rMax = 2.0;
        }
        const double ratioYmin = (rMin < 0.0) ? 1.2 * rMin : 0.0;
        const double ratioYmax = std::max(2.0, 1.2 * rMax);

        TCanvas* c = new TCanvas(Form("c_cmp_%s", h.tag.Data()), "", 760, 650);
        TPad* pTop = new TPad(Form("pTop_%s", h.tag.Data()), "", 0.0, 0.30, 1.0, 1.0);
        TPad* pBot = new TPad(Form("pBot_%s", h.tag.Data()), "", 0.0, 0.00, 1.0, 0.30);
        pTop->SetBottomMargin(0.01);
        pTop->SetLeftMargin(0.14);
        pTop->SetRightMargin(0.04);
        pBot->SetTopMargin(0.01);
        pBot->SetBottomMargin(0.33);
        pBot->SetLeftMargin(0.14);
        pBot->SetRightMargin(0.04);
        pTop->Draw();
        pBot->Draw();

        pTop->cd();
        h.sideband->GetXaxis()->SetLabelSize(0.0);
        h.sideband->GetXaxis()->SetTitleSize(0.0);
        h.sideband->GetXaxis()->SetTickLength(0.0);
        h.sideband->Draw("E");
        if (showBackgroundProjections) {
            h.backgroundSideband->Draw("E SAME");
            h.backgroundSplot->Draw("E SAME");
        }
        h.splot->Draw("E SAME");
        hMCBand->Draw("E2 SAME");
        h.mc->Draw("HIST SAME");
        if (showBackgroundProjections) {
            h.backgroundSideband->Draw("E SAME");
            h.backgroundSplot->Draw("E SAME");
        }
        h.sideband->Draw("E SAME");
        h.splot->Draw("E SAME");

        TLegend* signalLeg = new TLegend(0.68, showBackgroundProjections ? 0.65 : 0.64, 0.99, 0.91);
        signalLeg->SetBorderSize(0);
        signalLeg->SetFillStyle(0);
        signalLeg->SetTextFont(42);
        signalLeg->SetTextSize(showBackgroundProjections ? 0.045 : 0.055);
        signalLeg->SetHeader("#bf{Signal}", "L");
        signalLeg->AddEntry(h.mc, REWEIGHT_MC ? "Reweighted MC" : "MC", "l");
        signalLeg->AddEntry(h.sideband, "Sideband sub.", "lep");
        signalLeg->AddEntry(h.splot, "sPlot", "lep");
        TLegend* backgroundLeg = nullptr;
        if (showBackgroundProjections) {
            backgroundLeg = new TLegend(0.68, 0.47, 0.99, 0.64);
            backgroundLeg->SetBorderSize(0);
            backgroundLeg->SetFillStyle(0);
            backgroundLeg->SetTextFont(42);
            backgroundLeg->SetTextSize(0.045);
            backgroundLeg->SetHeader("#bf{Background}", "L");
            backgroundLeg->AddEntry(h.backgroundSideband, "Sideband", "lep");
            backgroundLeg->AddEntry(h.backgroundSplot, "sPlot", "lep");
        }

        TLatex cmpLabel;
        cmpLabel.SetNDC();
        cmpLabel.SetTextFont(42);
        cmpLabel.SetTextSize(0.060);
        cmpLabel.DrawLatex(0.18, 0.84, FitParticleLabel(treeName, true));
        const TString mcAgreementLabel = REWEIGHT_MC ? "Rew. MC" : "MC";
        cmpLabel.SetTextSize(0.040);
        cmpLabel.DrawLatex(0.18, 0.76, Form("#bf{sPlot-%s (%s)}", mcAgreementLabel.Data(), systemTag.Data()));
        if (agreement.ksPValue >= 0.0) {
            cmpLabel.DrawLatex(0.18, 0.70, Form("#bf{KS=%.3f, p=%.3g}",
                                                agreement.ksDistance, agreement.ksPValue));
        } else {
            cmpLabel.DrawLatex(0.18, 0.70, "#bf{KS=n/a}");
        }
        if (agreement.ndf > 0) {
            cmpLabel.DrawLatex(0.18, 0.64, Form("#bf{#chi^{2}/ndf=%.2f, p=%.3g (ndf=%d)}",
                                                reducedChi2, agreement.chi2PValue, agreement.ndf));
        } else {
            cmpLabel.DrawLatex(0.18, 0.64, "#bf{#chi^{2}/ndf=n/a}");
        }
        signalLeg->Draw();
        if (backgroundLeg) backgroundLeg->Draw();

        pBot->cd();
        hRatioSP->SetTitle("");
        hRatioSP->GetYaxis()->SetTitle("sPlot / MC");
        hRatioSP->GetYaxis()->SetTitleSize(0.09);
        hRatioSP->GetYaxis()->SetLabelSize(0.08);
        hRatioSP->GetYaxis()->SetTitleOffset(0.7);
        hRatioSP->GetYaxis()->SetNdivisions(304);
        hRatioSP->GetXaxis()->SetTitleSize(0.11);
        hRatioSP->GetXaxis()->SetLabelSize(0.10);
        hRatioSP->GetXaxis()->SetTitleOffset(1.1);
        hRatioSP->GetYaxis()->SetRangeUser(ratioYmin, ratioYmax);
        hRatioSP->Draw("E1 P");

        TLine* l0 = new TLine(hRatioSP->GetXaxis()->GetXmin(), 1.0, hRatioSP->GetXaxis()->GetXmax(), 1.0);
        l0->SetLineStyle(2);
        l0->SetLineWidth(2);
        l0->SetLineColor(kOrange + 7);
        l0->Draw("SAME");

        c->SaveAs(Form("%s/%s_%s_%s_compare%s.pdf", outDir.Data(), treeName.Data(), systemTag.Data(), h.tag.Data(),
                       REWEIGHT_MC ? Form("_reweighted_%s__weight%s", reweightTag.Data(), selectedWeightParticleTag.Data()) : ""));

        delete l0;
        delete signalLeg;
        delete backgroundLeg;
        delete pTop;
        delete pBot;
        delete c;
        delete hRatioSP;
        delete hWeight;
        delete hMCBand;
    }

    if (!REWEIGHT_MC) {
        fWeights->Close();
        delete fWeights;
    }
    for (auto& h : hists) {
        delete h.backgroundSideband;
        delete h.backgroundSplot;
        delete h.sideband;
        delete h.splot;
        delete h.mc;
    }
    delete reweightInput.hist;
    delete fitRes;
    delete allPars;

    closeInputFiles();

    if (REWEIGHT_MC) {
        std::cout << "Done. Outputs: " << outDir << "/*.pdf" << std::endl;
    } else {
        std::cout << "Done. Outputs: " << treeOutDir << "/*.pdf, "
                  << weightOutputPath << std::endl;
    }
}
