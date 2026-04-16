void optimalCUT() {

    TString TREE = "ntmix"; // "ntmix", "ntKp", "ntKstar", "ntphi"
    // Input files 
    TString path_to_file = "/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_DATA_wScore.root";
    TString path_to_MC   = "/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_MC_wScore.root";

    TFile *file_Data = TFile::Open(path_to_file.Data());
    TFile *file_MC_X3872 = TFile::Open(path_to_MC.Data());

    // Get trees (expect name "TREE")
    TTree *tree_DATA  = nullptr;
    TTree *tree_X3872 = nullptr;
    file_Data->GetObject(TREE, tree_DATA);
    file_MC_X3872->GetObject(TREE, tree_X3872);

    TString cut = "BQvalue < 0.15"; // 
    //TString cut = "1"; // 
    //TString cut = "BQvalue < 0.13 && (Bpt > 15 && Bpt < 50) && abs(By) < 1.6 && CentBin < 90 "; // Example final cuts from run2
    
    // Sideband definition for data (user requested): Bmass > 3.95 || (Bmass > 3.75 && Bmass < 3.8)
    TString sb_def = "((Bmass > 3.95 && Bmass < 4.0) || (Bmass > 3.75 && Bmass < 3.80))";

    // Scan BDTScore thresholds
    double best_thr = 0.03;
    double best_fom = -1.0;

    // Range and step for threshold scan
    double thr_min = 0.0;
    double thr_max = 1;
    double step = 0.005;

    std::cout<<"Starting BDT cut scan from "<<thr_min<<" to "<<thr_max<<" step "<<step<<std::endl;

    // Temporary histogram for projecting
    TH1F *htmp = new TH1F("htmp","htmp",1,-10,10);

    // Prepare arrays for graph
    int nsteps = int((thr_max - thr_min)/step + 1.5);
    double *xs = new double[nsteps];
    double *ys = new double[nsteps];
    int ipoint = 0;

    for (double thr = thr_min; thr <= thr_max + 1e-12; thr += step) {
        // Build selection strings
        TString sel_mc = Form("isX3872 == 1 && %s && (xgb_score > %g)", cut.Data(), thr);
        TString sel_data = Form("%s && %s && (xgb_score > %g)", sb_def.Data(), cut.Data(), thr);

        // Count MC signal events (S_mc)
        htmp->Reset();
        Int_t nS = tree_X3872->Project("htmp", "xgb_score", sel_mc.Data());

        // Count data sideband events (B_sb)
        htmp->Reset();
        Int_t nB = tree_DATA->Project("htmp", "xgb_score", sel_data.Data());

        double S_mc = static_cast<double>(nS);
        double B_sb = static_cast<double>(nB)/2;
        double fom = 0.0;
        double num =  2250.0/62847.0 * S_mc;
        double den = (2250.0/62847.0 * S_mc) + (11500.0/10132.0 * B_sb);
        if (den > 0.0) fom = num / TMath::Sqrt(den);
        if (!std::isfinite(fom)) fom = 0.0;

        std::cout<<Form("thr=%5.3f  S_mc=%.1f  B_sb=%.1f  FOM=%.2f", thr, S_mc, B_sb, fom)<<std::endl;

        // store for graph
        if (ipoint < nsteps) {
            xs[ipoint] = thr;
            ys[ipoint] = fom;
            ipoint++;
        }

        if (fom > best_fom) {
            best_fom = fom;
            best_thr = thr;
        }
    }

    delete htmp;

    std::cout<<"\nBest threshold: "<<best_thr<<"  Best FOM: "<<best_fom<<std::endl;

    // Create and save a graph of FOM vs threshold
    TCanvas *c = new TCanvas("c_optim","BDT optimization",800,600);
    TGraph *g = new TGraph(ipoint, xs, ys);

    // Build an explicit frame so PDF is never blank
    double ymax = 0.0;
    for (int i = 0; i < ipoint; ++i) {
        if (std::isfinite(ys[i]) && ys[i] > ymax) ymax = ys[i];
    }
    if (ymax <= 0.0) ymax = 1.0;

    TH1F *frame = c->DrawFrame(thr_min, 0.0, thr_max, 1.2 * ymax);
    frame->SetTitle(";BDTScore ;FOM");

    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetLineWidth(2);
    g->Draw("LP SAME");
    // draw vertical line at best thr
    TLine *l = new TLine(best_thr, 0, best_thr, 1.2 * ymax);
    l->SetLineColor(kRed);
    l->SetLineStyle(2);
    l->Draw();

    // Print best values directly on PDF
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.6, 0.85, Form("Best threshold = %.3f", best_thr));
    latex.DrawLatex(0.6, 0.82, Form("Best FOM = %.5f", best_fom));
    if (ipoint == 0) latex.DrawLatex(0.15, 0.77, "No scan points stored");

    c->SaveAs("bdt_optimization.pdf");

    delete[] xs; delete[] ys;
    delete c; delete g; delete l;

}
