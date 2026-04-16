void plot_xgb_score() {
    gStyle->SetOptStat(0);

    // Open files
    TFile *file_data = TFile::Open("/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_DATA_wScore.root");
    TFile *file_mc = TFile::Open("/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_MC_wScore.root");    
    // Get trees (adjust tree name if needed)
    TTree *tree_data = (TTree*)file_data->Get("ntmix");
    TTree *tree_mc = (TTree*)file_mc->Get("ntmix");
    // Create a single canvas and overlay both distributions
    TCanvas *c = new TCanvas("c", "XGB Score Comparison", 900, 700);
    c->SetMargin(0.12, 0.04, 0.12, 0.06);
    
    // Fill histograms
    tree_data->Draw("xgb_score >> h_data(100, 0, 1)", "", "hist");
    TH1F *h_data = (TH1F*)gDirectory->Get("h_data");
    tree_mc->Draw("xgb_score >> h_mc_x(100, 0, 1)", "isX3872==1", "hist");
    TH1F *h_mc_x = (TH1F*)gDirectory->Get("h_mc_x");
    tree_mc->Draw("xgb_score >> h_mc_psi(100, 0, 1)", "isX3872==0", "hist");
    TH1F *h_mc_psi = (TH1F*)gDirectory->Get("h_mc_psi");

    // Normalize to unit area
    double int_data = h_data->Integral();
    double int_mc_x = h_mc_x->Integral();
    double int_mc_psi = h_mc_psi->Integral();
    if (int_data > 0) h_data->Scale(1.0 / int_data);
    if (int_mc_x > 0) h_mc_x->Scale(1.0 / int_mc_x);
    if (int_mc_psi > 0) h_mc_psi->Scale(1.0 / int_mc_psi);

    // Style
    h_data->SetTitle("");
    h_data->SetLineColor(kBlack);
    h_data->SetMarkerColor(kBlack);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.7);
    h_data->SetLineWidth(2);
    h_data->SetStats(0);
    
    h_mc_x->SetLineColor(kOrange-3);
    h_mc_x->SetLineWidth(2);
    h_mc_x->SetStats(0);

    h_mc_psi->SetLineColor(kAzure+2);
    h_mc_psi->SetLineStyle(2);
    h_mc_psi->SetLineWidth(2);
    h_mc_psi->SetStats(0);

    // Draw on same pad
    double ymax_mc = TMath::Max(h_mc_x->GetMaximum(), h_mc_psi->GetMaximum());
    double ymax = TMath::Max(h_data->GetMaximum(), ymax_mc);
    h_data->SetMaximum(1.2 * ymax);
    h_data->Draw("E");
    h_mc_x->Draw("HIST SAME");
    h_mc_psi->Draw("HIST SAME");

    // Legend
    TLegend *leg = new TLegend(0.58, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h_data, "DATA", "lep");
    leg->AddEntry(h_mc_x, "MC X(3872)", "l");
    leg->AddEntry(h_mc_psi, "MC #psi(2S)", "l");
    leg->Draw();
    
    c->SaveAs("xgb_score_comparison.pdf");
    
    file_data->Close();
    file_mc->Close();
}