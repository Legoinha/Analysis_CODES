void plot_xgb_score() {
    gStyle->SetOptStat(0);

    // Open files
    TFile *file_data = TFile::Open("DATA_with_score.root");
    TFile *file_mc = TFile::Open("MC_with_score.root");
    if (!file_data || file_data->IsZombie() || !file_mc || file_mc->IsZombie()) {
        std::cout << "Error opening input ROOT files" << std::endl;
        return;
    }
    
    // Get trees (adjust tree name if needed)
    TTree *tree_data = (TTree*)file_data->Get("tree");
    TTree *tree_mc = (TTree*)file_mc->Get("tree");
    if (!tree_data || !tree_mc) {
        std::cout << "Tree 'tree' not found in one of the files" << std::endl;
        return;
    }

    // Create a single canvas and overlay both distributions
    TCanvas *c = new TCanvas("c", "XGB Score Comparison", 900, 700);
    c->SetMargin(0.12, 0.04, 0.12, 0.06);
    
    // Fill histograms
    tree_data->Draw("xgb_score >> h_data(100, 0, 1)", "", "hist");
    TH1F *h_data = (TH1F*)gDirectory->Get("h_data");
    tree_mc->Draw("xgb_score >> h_mc(100, 0, 1)", "", "hist");
    TH1F *h_mc = (TH1F*)gDirectory->Get("h_mc");

    // Normalize to unit area
    double int_data = h_data->Integral();
    double int_mc = h_mc->Integral();
    if (int_data > 0) h_data->Scale(1.0 / int_data);
    if (int_mc > 0) h_mc->Scale(1.0 / int_mc);

    // Style
    h_data->SetTitle("xgb_score comparison (normalized); xgb_score; ");
    h_data->SetLineColor(kBlack);
    h_data->SetMarkerColor(kBlack);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.7);
    h_data->SetLineWidth(2);
    h_data->SetStats(0);
    
    h_mc->SetLineColor(kRed);
    h_mc->SetLineWidth(2);
    h_mc->SetStats(0);

    // Draw on same pad
    double ymax = TMath::Max(h_data->GetMaximum(), h_mc->GetMaximum());
    h_data->SetMaximum(1.2 * ymax);
    h_data->Draw("E");
    h_mc->Draw("HIST SAME");

    // Legend
    TLegend *leg = new TLegend(0.62, 0.74, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_data, "DATA", "lep");
    leg->AddEntry(h_mc, "MC", "l");
    leg->Draw();
    
    c->SaveAs("xgb_score_comparison.pdf");
    
    file_data->Close();
    file_mc->Close();
}