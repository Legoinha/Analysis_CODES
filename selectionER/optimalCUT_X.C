void optimalCUT_X(TString system = "ppRef") {
    gStyle->SetOptStat(0);

    TFile* fileData = TFile::Open(Form("/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_%s_DATA_wScore.root", system.Data()));
    TFile* fileX =   TFile::Open(Form("/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_%s_MC_wScore_X3872.root", system.Data()));

    TTree *data = nullptr, *mcX = nullptr;
    fileData->GetObject("ntmix", data);
    fileX->GetObject("ntmix", mcX);

    TString sideband = "((Bmass > 3.95 && Bmass < 4.0))";
    double Fs_X = 1800 / 67709.00;
    double bkgScale = (23500) / 32000.00;
    double bestThr = 0., bestFom = -1., ymax = 0.;

    TH1F htmp("htmp", "", 1, 0, 1);
    TGraph gX;
    int ip = 0;

    for (double thr = 0.; thr <= 1.0001; thr += 0.01) {
        TString sel = Form("xgb_score > %.3f", thr);
        TString selData = Form("%s && %s", sideband.Data(), sel.Data());

        htmp.Reset();
        double sx = mcX->Project("htmp", "xgb_score", sel);
        htmp.Reset();
        double b = data->Project("htmp", "xgb_score", selData) ;

        double denX = Fs_X * sx + bkgScale * b;
        double fomX = (denX > 0.) ? Fs_X * sx / TMath::Sqrt(denX) : 0.;
        if (!std::isfinite(fomX)) fomX = 0.;
        fomX = std::round(fomX * 10.) / 10.;

        gX.SetPoint(ip, thr, fomX);
        if (fomX > bestFom) {
            bestFom = fomX;
            bestThr = thr;
        }
        if (fomX > ymax) ymax = fomX;
        // Print each iteration with Bkg*fb and Signal*fs
        std::cout << Form("thr = %.3f, FOM_X = %.2f, Bkg_DATA = %.2f, Signal_MC = %.2f", thr, fomX, b, sx) << std::endl;
        ip++;
    }

    TCanvas c("c", "", 800, 600);
    TH1F* frame = c.DrawFrame(0., 0., 1., 1.2 * ymax);
    frame->SetTitle(";xgb_score;FOM");

    gX.SetLineColor(kOrange - 3);
    gX.SetMarkerColor(kOrange - 3);
    gX.SetMarkerStyle(20);
    gX.SetLineWidth(2);
    gX.Draw("LP SAME");

    TLine lThr(bestThr, 0., bestThr, bestFom);
    TLine lFom(0., bestFom, bestThr, bestFom);
    lThr.SetLineColor(kOrange - 3);
    lFom.SetLineColor(kOrange - 3);
    lThr.SetLineStyle(2);
    lFom.SetLineStyle(2);
    lThr.Draw();
    lFom.Draw();

    TLegend leg(0.6, 0.75, 0.88, 0.86);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(&gX, "X(3872) FOM", "lp");
    leg.Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.15, 0.85, system);
    latex.DrawLatex(0.15, 0.80, Form("Best X(3872) thr. = %.2f", bestThr));
    c.SaveAs("bdt_optimization.pdf");

    std::cout << Form("%s: best X(3872) thr. = %.2f, FOM = %.1f", system.Data(), bestThr, bestFom) << std::endl;
}

int main(int argc, char** argv) {
    TString system = "ppRef";
    if (argc > 1) system = argv[1];
    optimalCUT_X(system);
    return 0;
}
