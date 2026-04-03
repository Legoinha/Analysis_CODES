void optimalCUT() {

    TString TREE = "tree"; // "ntmix", "ntKp", "ntKstar", "ntphi"
    // Input files 
    TString path_to_file = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/selected_events/bchi2prob_optuna5/DATA_with_score.root";
    TString path_to_MC   = "/eos/home-l/leyao/pbpb_work/X_analysis/XGBoost/selected_events/bchi2prob_optuna5/MC_with_score.root";


    TFile *file_Data = TFile::Open(path_to_file.Data());
    TFile *file_MC_X3872 = TFile::Open(path_to_MC.Data());
    if (!file_Data || file_Data->IsZombie()) { std::cout<<"Cannot open data file"<<std::endl; return; }
    if (!file_MC_X3872 || file_MC_X3872->IsZombie()) { std::cout<<"Cannot open MC file"<<std::endl; return; }

    // Get trees (expect name "TREE")
    TTree *tree_DATA  = nullptr;
    TTree *tree_X3872 = nullptr;
    file_Data->GetObject(TREE, tree_DATA);
    file_MC_X3872->GetObject(TREE, tree_X3872);
    if (!tree_DATA) {std::cout << "Data tree " << TREE << " not found"<<std::endl; return;}
    if (!tree_X3872){std::cout << "MC tree "   << TREE << " not found"<<std::endl; return;}

    TString cut = "BQvalue < 0.15"; // 
    //TString cut = "1"; // 
    //TString cut = "BQvalue < 0.13 && (Bpt > 15 && Bpt < 50) && abs(By) < 1.6 && CentBin < 90 "; // Example final cuts from run2
    
    // Sideband definition for data (user requested): Bmass > 3.95 || (Bmass > 3.75 && Bmass < 3.8)
    TString sb_def = "((Bmass > 3.95 && Bmass < 4.0) || (Bmass > 3.75 && Bmass < 3.80))";

    // Scan BDTScore thresholds
    double best_thr = 0.03;
    double best_fom = -1.0;

    // Range and step for threshold scan
    double thr_min = 0.4;
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
        TString sel_mc = Form("%s && (xgb_score > %g)"  , cut.Data(), thr);
        TString sel_data = Form("%s && %s && (xgb_score > %g)", sb_def.Data(), cut.Data(), thr);

        // Count MC signal events (S_mc)
        htmp->Reset();
        Int_t nS = tree_X3872->Project("htmp", "xgb_score", sel_mc.Data());

        // Count data sideband events (B_sb)
        htmp->Reset();
        Int_t nB = tree_DATA->Project("htmp", "xgb_score", sel_data.Data());

        double S_mc = static_cast<double>(nS);
        double B_sb = static_cast<double>(nB)/2;

        float fom = 0.0;
        float num = 3600.0/90000.0 * S_mc;
        float den = (3600.0/90000.0 * S_mc) + B_sb;
        fom = num / TMath::Sqrt(den);

        std::cout<<Form("thr=%5.3f  S_mc=%.1f  B_sb=%.1f  FOM=%.5f", thr, S_mc, B_sb, fom)<<std::endl;

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

    // Save result to a small text file
    TString outname = "bdt_optimization_result.txt";
    std::ofstream ofs(outname.Data());
    if (ofs.is_open()) {
        ofs<<"best_thr "<<best_thr<<"\n";
        ofs<<"best_fom "<<best_fom<<"\n";
        ofs.close();
        std::cout<<"Results written to "<<outname<<std::endl;
    }

    // Create and save a graph of FOM vs threshold
    TCanvas *c = new TCanvas("c_optim","BDT optimization",800,600);
    TGraph *g = new TGraph(ipoint, xs, ys);
    g->SetTitle("FOM vs BDTScore threshold;BDTScore threshold;FOM");
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetLineWidth(2);
    g->Draw("ALP");
    // draw vertical line at best thr
    TLine *l = new TLine(best_thr, 0, best_thr, g->GetHistogram() ? g->GetHistogram()->GetMaximum() : 1.0);
    l->SetLineColor(kRed);
    l->SetLineStyle(2);
    l->Draw();
    c->SaveAs("bdt_optimization.pdf");

    // Save graph to a ROOT file
    TFile fout("bdt_optimization.root","RECREATE");
    g->Write("g_fom_vs_thr");
    fout.Close();

    delete[] xs; delete[] ys;
    delete c; delete g; delete l;

}
