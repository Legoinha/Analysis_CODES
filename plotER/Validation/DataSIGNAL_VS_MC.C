#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TError.h>

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

using namespace RooFit;

struct VarCfgSignal {
    TString expr;
    TString title;
    int nbins;
    double xmin;
    double xmax;
    bool absVal;
};

static constexpr int kNBins = 15;

static bool isMixFamily(TString treeName)
{
    return (treeName == "ntmix" || treeName == "ntmix_psi2s");
}

static TString dataTreeName(TString treeName)
{
    return (treeName == "ntmix_psi2s") ? "ntmix" : treeName;
}

static TString makeTag(TString expr) {
    expr.ReplaceAll("(", "");
    expr.ReplaceAll(")", "");
    expr.ReplaceAll("/", "_");
    expr.ReplaceAll("*", "x");
    expr.ReplaceAll("+", "plus");
    expr.ReplaceAll("-", "minus");
    expr.ReplaceAll(" ", "");
    return expr;
}

static TString baseVarFromExpr(TString expr)
{
    expr.ReplaceAll(" ", "");
    if (expr.BeginsWith("abs(") && expr.EndsWith(")")) {
        expr = expr(4, expr.Length() - 5);
    }
    return expr;
}

static std::vector<VarCfgSignal> getSignalVars(TString treeName) {
    std::vector<VarCfgSignal> vars = {
        {"Bpt", ";p_{T} [GeV/c];", kNBins, 5, 50, false},
        {"By", ";|y|;", kNBins, 0, 2.4, true},
        {"BtrkPtimb", ";BtrkPtimb;", kNBins, 0.0, 1.0, false},
        {"Bchi2Prob", ";Bchi2Prob;", kNBins, 0.0, 1.0, false},
        {"Btrk1dR", ";Btrk1dR;", kNBins, 0.0, 0.6, false},
        {"Btrk2dR", ";Btrk2dR;", kNBins, 0.0, 0.6, false},
        {"Btrk1Pt", ";Btrk1Pt;", kNBins, 0.0, 6.0, false},
        {"Btrk2Pt", ";Btrk2Pt;", kNBins, 0.0, 6.0, false},
        {"Bnorm_svpvDistance_2D", ";Bnorm_svpvDistance_2D;", kNBins, 0.0, 20.0, false},
        {"Bnorm_trk1Dxy", ";Bnorm_trk1Dxy;", kNBins, -50.0, 50.0, false},
        {"Balpha", ";Balpha;", kNBins, 0.0, 3.14, false},
        {"BQvalue", ";BQvalue;", kNBins, 0.0, .2, false},
        {"BtktkvProb", ";BtktkvProb;", kNBins, 0.0, 1.0, false},
        {"Btktkmass", ";Btktkmass;", kNBins, 0.5, 1., false},
        {"Bcos_dtheta", ";Bcos_dtheta;", kNBins, 0.95, 1.0, false},
        {"BLxy", ";|BLxy|;", kNBins, 0.0, 0.5, true},
        {"BsvpvDistance_2D", ";BsvpvDistance_2D;", kNBins, 0.0, 0.25, false},
        {"Bujmass", ";Bujmass [GeV/c^{2}];", kNBins, 2.9, 3.25, false}
    };
    if (isMixFamily(treeName)) vars.insert(vars.begin() + 3, {"xgb_score", ";xgb_score;", kNBins, 0., 1.0, false});
    return vars;
}

static bool hasVarInTree(TTree* t, const TString& var)
{
    if (!t) return false;
    return (t->GetBranch(var.Data()) != nullptr) || (t->GetLeaf(var.Data()) != nullptr);
}

static std::vector<VarCfgSignal> getAvailableSignalVars(TString treeName, TTree* tData, TTree* tMC)
{
    std::vector<VarCfgSignal> out;
    for (const auto& v : getSignalVars(treeName)) {
        TString baseVar = baseVarFromExpr(v.expr);
        bool okData = hasVarInTree(tData, baseVar);
        bool okMC   = hasVarInTree(tMC, baseVar);
        if (!okData || !okMC) {
            std::cout << "[VarSkip] Skipping " << baseVar << " for " << treeName
                      << " (in data=" << okData << ", in mc=" << okMC << ")" << std::endl;
            continue;
        }
        out.push_back(v);
    }
    return out;
}

static TString resolveModelPath(TString modelPath, TString treeName)
{
    std::vector<TString> tries;
    tries.push_back(modelPath);
    tries.push_back(Form("../fitER/ROOTfiles/nominalFitModel_%s_ppRef.root", treeName.Data()));
    tries.push_back(Form("../../fitER/ROOTfiles/nominalFitModel_%s_ppRef.root", treeName.Data()));
    tries.push_back(Form("/eos/user/h/hmarques/Analysis_CODES/fitER/ROOTfiles/nominalFitModel_%s_ppRef.root", treeName.Data()));

    for (const auto& p : tries) {
        if (p.IsNull() || TString(p).Length() == 0) continue;
        if (!gSystem->AccessPathName(p.Data())) return p;
    }
    return "";
}

struct SigmaInfo {
    double mean;
    double sigma;
    bool valid;
};

struct MassWindows {
    double sigLo;
    double sigHi;
    double sbLLo;
    double sbLHi;
    double sbRLo;
    double sbRHi;
    bool valid;
};

static MassWindows windowsFromSigma(const SigmaInfo& s, double signalNSigma, double sidebandInNSigma, double sidebandOutNSigma)
{
    MassWindows w = {0, 0, 0, 0, 0, 0, false};
    if (!s.valid || !(s.sigma > 0.0)) return w;
    w.sigLo = s.mean - signalNSigma * s.sigma;
    w.sigHi = s.mean + signalNSigma * s.sigma;
    w.sbLLo = s.mean - sidebandOutNSigma * s.sigma;
    w.sbLHi = s.mean - sidebandInNSigma * s.sigma;
    w.sbRLo = s.mean + sidebandInNSigma * s.sigma;
    w.sbRHi = s.mean + sidebandOutNSigma * s.sigma;
    w.valid = true;
    return w;
}

static SigmaInfo extractSigmaFromModel(TString modelPath, TString treeName)
{
    SigmaInfo info = {0.0, 0.0, false};
    TString resolvedPath = resolveModelPath(modelPath, treeName);
    TFile* fModel = TFile::Open(resolvedPath, "READ");
    RooWorkspace* ws = (RooWorkspace*)fModel->Get("ws_nominal");

    RooRealVar* meanVar  = ws->var("mean1_");
    RooRealVar* sigma1   = ws->var("sigma11_");
    RooRealVar* sigma2   = ws->var("sigma21_");
    RooRealVar* sig1frac = ws->var("sig1frac1_");
    RooRealVar* scale    = ws->var("scale");
    
    if (!meanVar || !sigma1 || !sigma2 || !sig1frac || !scale) {
        std::cout << "[extractSigmaFromModel] Could not extract "
                  << "signal"
                  << " parameters from workspace" << std::endl;
        fModel->Close();
        return info;
    }
    
    // Calculate effective sigma as in roofitB.C
    double w = sig1frac->getVal();
    double effSigma = sqrt(w * pow(sigma1->getVal(), 2) + (1 - w) * pow(sigma2->getVal(), 2)) * scale->getVal();
    info.mean = meanVar->getVal();
    info.sigma = effSigma;
    info.valid = true;
    
    std::cout << "[extractSigmaFromModel] " << treeName
              << ": mean=" << info.mean
              << ", effSigma=" << info.sigma << std::endl;
    
    fModel->Close();
    return info;
}

static TString particleLabel(TString treeName)
{
    if (treeName == "ntmix") return "#bf{X(3872)}";
    if (treeName == "ntmix_psi2s") return "#bf{#psi(2S)}";
    if (treeName == "ntKp") return "#bf{B^{+}}";
    if (treeName == "ntKstar") return "#bf{B^{0}}";
    if (treeName == "ntphi") return "#bf{B_{s}^{0}}";
    return Form("#bf{%s}", treeName.Data());
}

static TString massFinalStateAxisTitle(TString treeName)
{
    if (isMixFamily(treeName)) return "m_{J/#psi #pi^{-} #pi^{+}} [GeV/c^{2}]";
    if (treeName == "ntKp")    return "m_{J/#psi K^{+}} [GeV/c^{2}]";
    if (treeName == "ntphi")   return "m_{J/#psi K^{+} K^{-}} [GeV/c^{2}]";
    if (treeName == "ntKstar") return "m_{J/#psi #pi^{+} K^{-}} [GeV/c^{2}]";
    return "m [GeV/c^{2}]";
}

static TH1D* makeMCUncBand(const TH1D* hMC, const TString& name)
{
    if (!hMC) return nullptr;
    TH1D* hBand = (TH1D*)hMC->Clone(name);
    hBand->SetFillColorAlpha(kOrange + 7, 0.30);
    hBand->SetFillStyle(1001);
    hBand->SetLineColor(kOrange + 7);
    hBand->SetLineWidth(1);
    hBand->SetMarkerSize(0);
    return hBand;
}

static TFile* openLocalRootFile(const TString& path, const char* mode)
{
    TString localUri = Form("file:%s", path.Data());
    return TFile::Open(localUri, mode);
}

static void ensureSPlotYieldRange(RooRealVar* y)
{
    if (!y) return;
    double ymin = y->getMin();
    double ymax = y->getMax();

    if (ymin > 0.0) ymin = 0.0;
    if (ymax < 1.0) ymax = 1.0;
    if (!(ymax > ymin)) ymax = ymin + 1.0;

    y->setRange(ymin, ymax);
}

void run_sideband_method(
    TString dataPath,
    TString mcPath,
    TString modelPath,
    TString baseCut,
    TString treeName)
{
    gSystem->mkdir("./BKGsub/", true);
    gSystem->mkdir(Form("./BKGsub/%s", treeName.Data()), true);
    TFile* fData = TFile::Open(dataPath, "READ");
    TFile* fMC   = TFile::Open(mcPath, "READ");
    TTree* tData = nullptr;
    TTree* tMC   = nullptr;
    fData->GetObject(dataTreeName(treeName), tData);
    fMC->GetObject(treeName, tMC);

    // Region definition in sigma units
    const double signalNSigma = (treeName == "ntmix_psi2s") ? 3. : 2.;
    const double sidebandInNSigma = 4.;
    const double sidebandOutNSigma = 8.0;

    TString sigLabelPlain     = Form("Signal region (+/-%gsigma)", signalNSigma);
    TString sbLabelPlain      = Form("Sideband region (%g-%gsigma)", sidebandInNSigma, sidebandOutNSigma);
    TString sigLabelRoot      = Form("Signal region (#pm%g#sigma)", signalNSigma);
    TString sbLabelRoot       = Form("Sideband region (%g-%g#sigma)", sidebandInNSigma, sidebandOutNSigma);

    SigmaInfo sigInfo = extractSigmaFromModel(modelPath, treeName);
    MassWindows w = windowsFromSigma(sigInfo, signalNSigma, sidebandInNSigma, sidebandOutNSigma);
    if (!w.valid) {
        Error("run_sideband_method", "Could not retrieve sigma-based regions for %s", treeName.Data());
        fData->Close();
        fMC->Close();
        return;
    }

    const double sigLo = w.sigLo;
    const double sigHi = w.sigHi;
    const double sbLLo = w.sbLLo;
    const double sbLHi = w.sbLHi;
    const double sbRLo = w.sbRLo;
    const double sbRHi = w.sbRHi;

    std::cout << "[run_sideband_method] Using sigma-based regions for " << treeName << std::endl;
    std::cout << "  " << sigLabelPlain << ": [" << sigLo << ", " << sigHi << "]" << std::endl;
    std::cout << "  " << sbLabelPlain << ": [" << sbLLo << ", " << sbLHi << "] U [" << sbRLo << ", " << sbRHi << "]" << std::endl;

    const double wSig = sigHi - sigLo;
    const bool isNtKp = (treeName == "ntKp");
    const double wSB  = isNtKp ? (sbRHi - sbRLo) : ((sbLHi - sbLLo) + (sbRHi - sbRLo));
    const double alpha = (wSB > 0.0 ? wSig / wSB : 0.0);

    TString cutSig = Form("(%s) && (Bmass>%f && Bmass<%f)", baseCut.Data(), sigLo, sigHi);
    TString cutSB;
    if (!isNtKp) {
        cutSB = Form("(%s) && ((Bmass>%f && Bmass<%f) || (Bmass>%f && Bmass<%f))",
                     baseCut.Data(), sbLLo, sbLHi, sbRLo, sbRHi);
    } else {
        cutSB = Form("(%s) && (Bmass>%f && Bmass<%f)", baseCut.Data(), sbRLo, sbRHi);
    }
    TString cutMCSig = Form("(%s)", baseCut.Data());

    // Diagnostic plot: show mass windows used for sideband subtraction
    double massMinPlot = isMixFamily(treeName) ? 3.6 : 5.0;
    double massMaxPlot = isMixFamily(treeName) ? 4.0 : 5.8;
    const int nMassBins = 80;
    TH1D* hMassWin = new TH1D("hMassWin_tmp", Form(";%s;", massFinalStateAxisTitle(treeName).Data()), nMassBins, massMinPlot, massMaxPlot);
    const double massBinWidthMeV = 1000.0 * (massMaxPlot - massMinPlot) / nMassBins;
    hMassWin->GetYaxis()->SetTitle(Form("Entries / %.4g MeV/c^{2}", massBinWidthMeV));
    tData->Draw("Bmass>>hMassWin_tmp", baseCut, "goff");

    TCanvas* cWin = new TCanvas(Form("cWin_%s", treeName.Data()), "", 760, 650);
    cWin->SetLeftMargin(0.14);
    hMassWin->SetLineColor(kBlack);
    hMassWin->SetMarkerStyle(20);
    hMassWin->SetMarkerSize(0.7);
    hMassWin->Draw("E");

    double yMaxWin = hMassWin->GetMaximum() * 1.15;
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
    legWin->SetHeader(particleLabel(treeName), "C");
    legWin->AddEntry(hMassWin, "Data", "lep");
    legWin->AddEntry(bSig, sigLabelRoot, "f");
    legWin->AddEntry(isNtKp ? bSBR : bSBL, sbLabelRoot, "f");
    legWin->Draw();

    cWin->SaveAs(Form("./BKGsub/%s/mass_windows_%s.pdf", treeName.Data(), treeName.Data()));

    delete legWin;
    delete bSig;
    delete bSBL;
    delete bSBR;
    delete cWin;
    delete hMassWin;

    TFile* fout = openLocalRootFile(Form("./bkgSub_%s.root", treeName.Data()), "RECREATE");

    auto vars = getAvailableSignalVars(treeName, tData, tMC);
    for (const auto& v : vars) {
        TString expr = v.expr;
        if (v.absVal) expr = Form("abs(%s)", v.expr.Data());
        TString tag = makeTag(expr);

        TH1D* hDataSR = new TH1D(Form("hDataSR_%s", tag.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hDataSB = new TH1D(Form("hDataSB_%s", tag.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hMCsig  = new TH1D(Form("hMCsig_%s",  tag.Data()), v.title, v.nbins, v.xmin, v.xmax);
        hMCsig->Sumw2();
        tData->Draw(Form("%s>>%s", expr.Data(), hDataSR->GetName()), cutSig, "goff");
        tData->Draw(Form("%s>>%s", expr.Data(), hDataSB->GetName()), cutSB,  "goff");
        tMC->Draw(Form("%s>>%s",   expr.Data(), hMCsig->GetName()),  cutMCSig, "goff");

        TH1D* hDataSub = (TH1D*)hDataSR->Clone(Form("hDataSub_%s", tag.Data()));
        hDataSub->Add(hDataSB, -alpha);
        TH1D* hDataSBShape = (TH1D*)hDataSB->Clone(Form("hDataSBShape_%s", tag.Data()));

        if (hDataSub->Integral() > 0) hDataSub->Scale(1.0 / hDataSub->Integral());
        if (hDataSBShape->Integral() > 0) hDataSBShape->Scale(1.0 / hDataSBShape->Integral());
        if (hMCsig->Integral() > 0)   hMCsig->Scale(1.0 / hMCsig->Integral());

        hDataSub->SetLineColor(kBlack);
        hDataSub->SetMarkerColor(kBlack);
        hDataSub->SetMarkerStyle(20);
        hDataSub->SetLineWidth(2);
        hDataSBShape->SetLineColorAlpha(kGray + 1, 0.45);
        hDataSBShape->SetMarkerColorAlpha(kGray + 1, 0.45);
        hDataSBShape->SetMarkerStyle(21);
        hDataSBShape->SetLineWidth(1);
        hMCsig->SetLineColor(kOrange + 7);
        hMCsig->SetLineWidth(2);
        TH1D* hMCsigBand = makeMCUncBand(hMCsig, Form("hMCsigBand_%s", tag.Data()));

        double ymax = std::max(hDataSub->GetMaximum(), hMCsig->GetMaximum());
        ymax = std::max(ymax, hDataSBShape->GetMaximum());
        hDataSub->SetMinimum(0.0);
        hDataSub->SetMaximum(1.35 * ymax);

        TCanvas* c = new TCanvas(Form("c_sb_%s", tag.Data()), "", 760, 650);
        c->SetLeftMargin(0.14);
        hDataSub->Draw("E");
        hDataSBShape->Draw("E SAME");
        hMCsigBand->Draw("E2 SAME");
        hMCsig->Draw("HIST SAME");
        hDataSub->Draw("E SAME");

        TLegend* leg = new TLegend(0.65, 0.74, 0.90, 0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetHeader(particleLabel(treeName), "C");
        leg->AddEntry(hDataSub, "Sideband Subtraction (sig)", "lep");
        leg->AddEntry(hDataSBShape, "Sideband bkg", "lep");
        leg->AddEntry(hMCsig, "MC", "l");
        leg->Draw();

        c->SaveAs(Form("./BKGsub/%s/%s_bkgSub.pdf", treeName.Data(), tag.Data()));

        fout->cd();
        hDataSR->Write();
        hDataSB->Write();
        hDataSub->Write();
        hMCsig->Write();

        delete leg;
        delete c;
        delete hDataSR;
        delete hDataSB;
        delete hDataSub;
        delete hDataSBShape;
        delete hMCsigBand;
        delete hMCsig;

    }

    fout->Close();
    fData->Close();
    fMC->Close();
}

void run_splot_method(
    TString dataPath,
    TString mcPath,
    TString modelPath,
    TString baseCut,
    TString treeName)
{
    gSystem->mkdir("./SPLOT/", true);
    gSystem->mkdir(Form("./SPLOT/%s", treeName.Data()), true);

    TString modelPathResolved = resolveModelPath(modelPath, treeName);
    if (modelPathResolved.IsNull() || modelPathResolved.Length() == 0) {
        Error("run_splot_method", "Could not resolve model path for %s", treeName.Data());
        return;
    }

    TFile* fData = TFile::Open(dataPath, "READ");
    TFile* fMC   = TFile::Open(mcPath, "READ");
    TFile* fModel = TFile::Open(modelPathResolved, "READ");
    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie() || !fModel || fModel->IsZombie()) {
        Error("run_splot_method", "Could not open data/mc/model files for %s", treeName.Data());
        if (fData) fData->Close();
        if (fMC) fMC->Close();
        if (fModel) fModel->Close();
        return;
    }
    TTree* tData = nullptr;
    TTree* tMC   = nullptr;
    fData->GetObject(dataTreeName(treeName), tData);
    fMC->GetObject(treeName, tMC);
    RooWorkspace* ws = (RooWorkspace*)fModel->Get("ws_nominal");
    if (!tData || !tMC || !ws) {
        Error("run_splot_method", "Missing tree/workspace for %s", treeName.Data());
        fData->Close();
        fMC->Close();
        fModel->Close();
        return;
    }
    double massMin = 5.0;
    double massMax = 5.8;
    if (treeName == "ntmix") {
        massMin = 3.8;
        massMax = 4.0;
    } else if (treeName == "ntmix_psi2s") {
        massMin = 3.6;
        massMax = 3.8;
    }

    auto vars = getAvailableSignalVars(treeName, tData, tMC);

    RooRealVar Bmass("Bmass", "Bmass", massMin, massMax);
    RooArgSet obs(Bmass);
    std::vector<std::unique_ptr<RooRealVar>> extraObs;
    for (const auto& v : vars) {
        TString baseVar = baseVarFromExpr(v.expr);
        if (baseVar == "Bmass") continue;
        bool alreadyThere = false;
        for (const auto& rv : extraObs) {
            if (baseVar == rv->GetName()) {
                alreadyThere = true;
                break;
            }
        }
        if (alreadyThere) continue;
        extraObs.emplace_back(new RooRealVar(baseVar, baseVar, -1e6, 1e6));
        obs.add(*extraObs.back());
    }

    TString dataCut = Form("(%s) && (Bmass>%f && Bmass<%f)", baseCut.Data(), massMin, massMax);
    RooDataSet data("data", "data", tData, obs, dataCut.Data());
    RooAbsPdf* model = ws->pdf("model1_");
    RooRealVar* nsig = ws->var("nsig1_");
    RooRealVar* nbkg = ws->var("nbkg1_");
    RooRealVar* nbkgPartR = ws->var("nbkg_part_r1_");
    if (!model || !nsig || !nbkg) {
        Error("run_splot_method", "Missing model/yield params in workspace for %s", treeName.Data());
        fData->Close();
        fMC->Close();
        fModel->Close();
        return;
    }

    ensureSPlotYieldRange(nsig);
    ensureSPlotYieldRange(nbkg);
    ensureSPlotYieldRange(nbkgPartR);

    RooArgSet* allPars = model->getParameters(data);
    std::unique_ptr<TIterator> it(allPars->createIterator());
    TObject* obj = nullptr;
    while ((obj = it->Next())) {
        RooRealVar* v = dynamic_cast<RooRealVar*>(obj);
        if (!v) continue;
        v->setConstant(true);
    }
    nsig->setConstant(false);
    nbkg->setConstant(false);
    if (nbkgPartR) nbkgPartR->setConstant(false);

    RooFitResult* fitRes = model->fitTo(data, Extended(true), Save(true), PrintLevel(-1));

    RooArgList splotYields;
    splotYields.add(*nsig);
    if (nbkgPartR) splotYields.add(*nbkgPartR);
    splotYields.add(*nbkg);
    RooStats::SPlot sData("sData", "An SPlot", data, model, splotYields);
    if (!std::isfinite(sData.GetYieldFromSWeight(nsig->GetName()))) {
        Error("run_splot_method", "sPlot failed for %s (non-finite signal yield from sWeights).", treeName.Data());
        if (fitRes) delete fitRes;
        fData->Close();
        fMC->Close();
        fModel->Close();
        return;
    }

    // Diagnostic: check sWeight statistics
    double swMin = 1e9, swMax = -1e9, swMean = 0;
    int negWeights = 0, totalWeights = 0;
    TString swVarName = Form("%s_sw", nsig->GetName());
    for (int i = 0; i < data.numEntries(); ++i) {
        const RooArgSet* row = data.get(i);
        if (!row) continue;
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

    TCanvas* cMass = new TCanvas(Form("cMass_%s", treeName.Data()), "", 760, 650);
    RooPlot* frame = Bmass.frame();
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(massFinalStateAxisTitle(treeName));
    frame->SetMinimum(0.0);
    data.plotOn(frame);
    model->plotOn(frame);
    if (fitRes) model->paramOn(frame, Layout(0.55, 0.99, 0.97));
    frame->Draw();

    cMass->SaveAs(Form("./SPLOT/%s/massFit_splot_%s.pdf", treeName.Data(), treeName.Data()));

    TFile* fout = openLocalRootFile(Form("./splot_%s.root", treeName.Data()), "RECREATE");
    TString mcCut = Form("(%s)", baseCut.Data());

    for (const auto& v : vars) {
        TString exprPlot = v.absVal ? Form("abs(%s)", v.expr.Data()) : v.expr;
        TString tag = makeTag(exprPlot);
        TString baseVar = baseVarFromExpr(v.expr);

        TH1D* hDataSPlot = new TH1D(Form("hDataSPlot_%s", tag.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hDataBkgSPlot = new TH1D(Form("hDataBkgSPlot_%s", tag.Data()), v.title, v.nbins, v.xmin, v.xmax);
        TH1D* hMCsig     = new TH1D(Form("hMCsig_%s", tag.Data()), v.title, v.nbins, v.xmin, v.xmax);
        hMCsig->Sumw2();

        for (int i = 0; i < data.numEntries(); ++i) {
            const RooArgSet* row = data.get(i);
            if (!row) continue;
            double val = row->getRealValue(baseVar.Data());
            if (v.absVal) val = std::abs(val);
            double w = row->getRealValue(Form("%s_sw", nsig->GetName()));
            hDataSPlot->Fill(val, w);
            double wb = row->getRealValue(Form("%s_sw", nbkg->GetName()));
            if (nbkgPartR) wb += row->getRealValue(Form("%s_sw", nbkgPartR->GetName()));
            hDataBkgSPlot->Fill(val, wb);
        }

        TString mcExpr = v.absVal ? Form("abs(%s)", v.expr.Data()) : v.expr;
        tMC->Draw(Form("%s>>%s", mcExpr.Data(), hMCsig->GetName()), mcCut, "goff");

        if (hDataSPlot->Integral() > 0) hDataSPlot->Scale(1.0 / hDataSPlot->Integral());
        if (hDataBkgSPlot->Integral() > 0) hDataBkgSPlot->Scale(1.0 / hDataBkgSPlot->Integral());
        if (hMCsig->Integral() > 0)     hMCsig->Scale(1.0 / hMCsig->Integral());

        hDataSPlot->SetLineColor(kBlue + 1);
        hDataSPlot->SetMarkerColor(kBlue + 1);
        hDataSPlot->SetMarkerStyle(24);
        hDataSPlot->SetLineWidth(2);
        hMCsig->SetLineColor(kOrange + 7);
        hMCsig->SetLineWidth(2);
        TH1D* hMCsigBand = makeMCUncBand(hMCsig, Form("hMCsigBand_sp_%s", tag.Data()));

        double ymax = std::max(hDataSPlot->GetMaximum(), hMCsig->GetMaximum());
        hDataSPlot->SetMinimum(0.0);
        hDataSPlot->SetMaximum(1.35 * ymax);

        TCanvas* c = new TCanvas(Form("c_sp_%s", tag.Data()), "", 760, 650);
        c->SetLeftMargin(0.14);
        hDataSPlot->Draw("E");
        hMCsigBand->Draw("E2 SAME");
        hMCsig->Draw("HIST SAME");
        hDataSPlot->Draw("E SAME");

        TLegend* leg = new TLegend(0.65, 0.70, 0.90, 0.89);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetHeader(particleLabel(treeName), "C");
        leg->AddEntry(hDataSPlot, "sPlot (sig)", "lep");
        leg->AddEntry(hMCsig, "MC", "l");
        leg->Draw();

        c->SaveAs(Form("./SPLOT/%s/%s_splot.pdf", treeName.Data(), tag.Data()));

        fout->cd();
        hDataSPlot->Write();
        hDataBkgSPlot->Write();
        hMCsig->Write();

        delete leg;
        delete c;
        delete hDataSPlot;
        delete hDataBkgSPlot;
        delete hMCsigBand;
        delete hMCsig;

    }

    fout->Close();
    delete frame;
    delete cMass;
    if (fitRes) delete fitRes;
    fData->Close();
    fMC->Close();
    fModel->Close();
}

void run_compare_methods(
    TString sidebandFile,
    TString splotFile,
    TString treeName,
    TString suffix = "")
{
    gSystem->mkdir("./COMPARE/", true);
    gSystem->mkdir(Form("./COMPARE/%s", treeName.Data()), true);

    TFile* fSB = TFile::Open(sidebandFile, "READ");
    TFile* fSP = TFile::Open(splotFile, "READ");
    if (!fSB || fSB->IsZombie() || !fSP || fSP->IsZombie()) {
        Error("run_compare_methods", "Could not open sideband/splot output files");
        return;
    }

    auto vars = getSignalVars(treeName);
    for (const auto& v : vars) {
        TString baseVar = baseVarFromExpr(v.expr);
        if (baseVar == "Bnorm_trk1Dxy" || baseVar == "Balpha") continue;
        TString exprPlot = v.absVal ? Form("abs(%s)", v.expr.Data()) : v.expr;
        TString tag = makeTag(exprPlot);
        TString label = v.expr;

        TH1D* hSBraw = (TH1D*)fSB->Get(Form("hDataSub_%s", tag.Data()));
        TH1D* hMCraw = (TH1D*)fSB->Get(Form("hMCsig_%s",   tag.Data()));
        TH1D* hSBbkgRaw = (TH1D*)fSB->Get(Form("hDataSB_%s", tag.Data()));
        TH1D* hSPraw = (TH1D*)fSP->Get(Form("hDataSPlot_%s", tag.Data()));
        TH1D* hSPbkgRaw = (TH1D*)fSP->Get(Form("hDataBkgSPlot_%s", tag.Data()));
        if (!hSBraw || !hMCraw || !hSPraw || !hSBbkgRaw || !hSPbkgRaw) continue;

        TH1D* hSB = (TH1D*)hSBraw->Clone(Form("hSB_cmp_%s", tag.Data()));
        TH1D* hMC = (TH1D*)hMCraw->Clone(Form("hMC_cmp_%s", tag.Data()));
        TH1D* hSP = (TH1D*)hSPraw->Clone(Form("hSP_cmp_%s", tag.Data()));
        TH1D* hSBbkg = (TH1D*)hSBbkgRaw->Clone(Form("hSBbkg_cmp_%s", tag.Data()));
        TH1D* hSPbkg = (TH1D*)hSPbkgRaw->Clone(Form("hSPbkg_cmp_%s", tag.Data()));

        if (hSB->Integral() > 0) hSB->Scale(1.0 / hSB->Integral());
        if (hMC->Integral() > 0) hMC->Scale(1.0 / hMC->Integral());
        if (hSP->Integral() > 0) hSP->Scale(1.0 / hSP->Integral());
        if (hSBbkg->Integral() > 0) hSBbkg->Scale(1.0 / hSBbkg->Integral());
        if (hSPbkg->Integral() > 0) hSPbkg->Scale(1.0 / hSPbkg->Integral());

        hMC->SetLineColor(kOrange + 7);
        hMC->SetLineWidth(2);
        hSB->SetLineColor(kBlue + 1);
        hSB->SetMarkerColor(kBlue + 1);
        hSB->SetMarkerStyle(20);
        hSB->SetLineWidth(2);
        hSP->SetLineColor(kRed + 1);
        hSP->SetMarkerColor(kRed + 1);
        hSP->SetMarkerStyle(24);
        hSP->SetLineWidth(2);
        hSBbkg->SetLineColorAlpha(kBlack, 0.35);
        hSBbkg->SetMarkerColorAlpha(kBlack, 0.35);
        hSBbkg->SetMarkerStyle(21);
        hSBbkg->SetLineWidth(1);
        hSPbkg->SetLineColorAlpha(kGray + 1, 0.45);
        hSPbkg->SetMarkerColorAlpha(kGray + 1, 0.45);
        hSPbkg->SetMarkerStyle(25);
        hSPbkg->SetLineWidth(1);
        TH1D* hMCBand = makeMCUncBand(hMC, Form("hMCBand_cmp_%s", tag.Data()));

        TH1D* hDiffSB = (TH1D*)hSB->Clone(Form("hDiffSB_%s", tag.Data()));
        TH1D* hDiffSP = (TH1D*)hSP->Clone(Form("hDiffSP_%s", tag.Data()));
        hDiffSB->Add(hMC, -1.0);
        hDiffSP->Add(hMC, -1.0);

        hDiffSB->SetLineColor(kBlue + 1);
        hDiffSB->SetMarkerColor(kBlue + 1);
        hDiffSB->SetMarkerStyle(20);
        hDiffSB->SetLineWidth(2);
        hDiffSP->SetLineColor(kRed + 1);
        hDiffSP->SetMarkerColor(kRed + 1);
        hDiffSP->SetMarkerStyle(24);
        hDiffSP->SetLineWidth(2);

        double ymax = std::max(hMC->GetMaximum(), std::max(hSB->GetMaximum(), hSP->GetMaximum()));
        hSB->SetMaximum(1.35 * ymax);
        hSB->SetMinimum(0.0);

        double dmax = std::max(std::abs(hDiffSB->GetMaximum()), std::abs(hDiffSB->GetMinimum()));
        dmax = std::max(dmax, std::abs(hDiffSP->GetMaximum()));
        dmax = std::max(dmax, std::abs(hDiffSP->GetMinimum()));
        if (dmax <= 0) dmax = 0.05;

        TCanvas* c = new TCanvas(Form("c_cmp_%s", tag.Data()), "", 760, 650);
        TPad* pTop = new TPad(Form("pTop_%s", tag.Data()), "", 0.0, 0.30, 1.0, 1.0);
        TPad* pBot = new TPad(Form("pBot_%s", tag.Data()), "", 0.0, 0.00, 1.0, 0.30);
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
        hSB->GetXaxis()->SetLabelSize(0.0);
        hSB->GetXaxis()->SetTitleSize(0.0);
        hSB->GetXaxis()->SetTickLength(0.0);
        hSB->Draw("E");
        hSP->Draw("E SAME");
        hSBbkg->Draw("E SAME");
        hSPbkg->Draw("E SAME");
        hMCBand->Draw("E2 SAME");
        hMC->Draw("HIST SAME");
        hSB->Draw("E SAME");
        hSP->Draw("E SAME");

        TLegend* leg = new TLegend(0.7, 0.68, 0.9, 0.9);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetHeader(particleLabel(treeName), "C");
        leg->AddEntry(hMC, "MC", "l");
        leg->AddEntry(hSB, "Sideband Subtraction (sig)", "lep");
        leg->AddEntry(hSP, "sPlot (sig)", "lep");
        leg->AddEntry(hSBbkg, "Sideband bkg", "lep");
        leg->AddEntry(hSPbkg, "sPlot bkg", "lep");
        leg->Draw();

        pBot->cd();
        hDiffSB->SetTitle("");
        hDiffSB->GetYaxis()->SetTitle("Data - MC");
        hDiffSB->GetYaxis()->SetTitleSize(0.09);
        hDiffSB->GetYaxis()->SetLabelSize(0.08);
        hDiffSB->GetYaxis()->SetTitleOffset(0.7);
        hDiffSB->GetYaxis()->SetNdivisions(304);
        hDiffSB->GetXaxis()->SetTitleSize(0.11);
        hDiffSB->GetXaxis()->SetLabelSize(0.10);
        hDiffSB->GetXaxis()->SetTitleOffset(1.1);
        hDiffSB->GetYaxis()->SetRangeUser(-1.2 * dmax, 1.2 * dmax);
        hDiffSB->Draw("E");
        hDiffSP->Draw("E SAME");

        TLine* l0 = new TLine(hDiffSB->GetXaxis()->GetXmin(), 0.0, hDiffSB->GetXaxis()->GetXmax(), 0.0);
        l0->SetLineStyle(2);
        l0->SetLineWidth(2);
        l0->SetLineColor(kOrange + 7);
        l0->Draw("SAME");

        c->SaveAs(Form("./COMPARE/%s/%s%s_compare.pdf", treeName.Data(), tag.Data(), suffix.Data()));

        delete l0;
        delete leg;
        delete pTop;
        delete pBot;
        delete c;
        delete hDiffSB;
        delete hDiffSP;
        delete hSB;
        delete hSBbkg;
        delete hMCBand;
        delete hMC;
        delete hSPbkg;
        delete hSP;
    }

    fSB->Close();
    fSP->Close();
}

void DataSIGNAL_VS_MC(
    TString dataPath  = "/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_DATA_wScore.root",
    TString mcPath    = "/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_MC_wScore_X3872.root",
    TString modelPath = "../fitER/ROOTfiles/nominalFitModel_ntmix_ppRef.root",
    TString baseCut   = "xgb_score > 0.61 && BQvalue<0.15",
    TString treeName  = "ntmix")
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gSystem->mkdir("./SPLOT/", true);
    gSystem->mkdir("./BKGsub/", true);

    std::cout << "Running Sideband and sPlot from one macro..." << std::endl;
    run_sideband_method(dataPath, mcPath, modelPath, baseCut, treeName);
    run_splot_method(dataPath, mcPath, modelPath, baseCut, treeName);
    run_compare_methods(Form("./bkgSub_%s.root", treeName.Data()), Form("./splot_%s.root", treeName.Data()), treeName, "");
    std::cout << "Done. Outputs: ./bkgSub_" << treeName << ".root, ./splot_" << treeName << ".root and ./COMPARE/" << treeName << "/*.pdf" << std::endl;
}
