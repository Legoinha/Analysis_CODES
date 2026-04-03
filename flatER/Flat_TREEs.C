#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <vector>
#include <iostream>
#include <cmath>
#include "../plotER/aux/masses.h"

void Flat_TREEs(TString treename="ntmix", TString P_vs_NP="prompt", TString SYSTEM="PbPb23", bool isMC = true)
{

    bool Fid_region  = true;  // apply fiducial region selection
    bool PreSel_cuts = true;  // apply quality cuts (+ pre-selection)
    bool KeepGenLv = true;     // whether to copy gen-level tree (ntGen) to output file (only for MC)
    TString specCASES = "";    // for extra naming in output file (e.g. to distinguish different selection cuts)
    if (treename == "ntmix" && isMC && P_vs_NP == "NP") { specCASES = "_nonPrompt";}
    if (KeepGenLv == false) { specCASES = "_RECOonly";}

    std::cout << "Flattening tree: " << treename << " , System: " << SYSTEM << " , isMC: " << isMC << std::endl;
    bool isPP     = (SYSTEM == "ppRef");
    bool isPbPb23 = (SYSTEM == "PbPb23");
    bool isPbPb24 = (SYSTEM == "PbPb24");
    bool isPbPb25 = (SYSTEM == "PbPb25");
    bool isPbPb   = (isPbPb23 || isPbPb24 || isPbPb25);

    // --------------------------------------------------
    // Input files using TChain
    // --------------------------------------------------
    TChain *tin = new TChain("Bfinder/" + treename);
    
    if(isPP){ // MC ppRef
        if (isMC){ 
            if (treename == "ntmix"){
                if (P_vs_NP == "NP"){
                    tin->Add("/eos/user/h/hmarques/MC_ppRef_X3872/nonprompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder/260226_100453/0000/*.root");   //PSI2S (non-prompt, small)
                    tin->Add("/eos/user/h/hmarques/MC_ppRef_X3872/nonprompt_X3872_to_Jpsi_Rho_phat5_Bfinder/260225_232446/0000/*.root");    //X3872 (non-prompt, small)
                } else {
                    tin->Add("/eos/user/h/hmarques/MC_ppRef_X3872/prompt_X3872_to_Jpsi_Rho_phat5_Bfinder/260228_194148/0000/*.root");    //X3872 (small)
                    tin->Add("/eos/user/h/hmarques/MC_ppRef_X3872/prompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder/260228_194229/0000/*.root");   //PSI2S (small)
                    tin->Add("/eos/user/h/hmarques/MC_ppRef_X3872/prompt_X3872_to_Jpsi_Rho_phat5_Bfinder/260228_181540/0000/*.root");    //X3872 (large)
                    tin->Add("/eos/user/h/hmarques/MC_ppRef_X3872/prompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder/260228_181719/0000/*.root");   //PSI2S (large)
                }
            }
            else if (treename == "ntphi"){  tin->Add("/eos/user/h/hmarques/MC_ppRef_Bmesons/MC_ppRef_X/Bsubs_to_JPsiPHI_pthat10/260224_023406/0000/*.root");}
            else if (treename == "ntKstar"){tin->Add("/eos/user/h/hmarques/MC_ppRef_Bmesons/MC_ppRef_X/Bzero_to_JpsiKstar_pthat10/260224_023043/0000/*.root");}
            else if (treename == "ntKp"){   tin->Add("/eos/user/h/hmarques/MC_ppRef_Bmesons/MC_ppRef_X/Bplus_to_JPsiK_pthat10/260224_125939/0000/*.root");}
        }
        else{ // DATA ppRef
            tin->Add("/eos/user/h/hmarques/PPRef2024/PPRefDoubleMuon0/PPRefDoubleMuon0/260225_125246/0000/*"); 
            tin->Add("/eos/user/h/hmarques/PPRef2024/PPRefDoubleMuon1/PPRefDoubleMuon1/260225_125254/0000/*"); 
            tin->Add("/eos/user/h/hmarques/PPRef2024/PPRefDoubleMuon2/PPRefDoubleMuon2/260225_125302/0000/*");  
            tin->Add("/eos/user/h/hmarques/PPRef2024/PPRefDoubleMuon3/PPRefDoubleMuon3/260225_125311/0000/*"); 
        }
    }
    else if(isPbPb23){
         if(isMC){
            //tin->Add("/eos/user/h/hmarques/MC_PbPb23_X3872/MC_PbPb_X3872/prompt_X3872_to_Jpsi_Rho_phat5_Bfinder/260322_030836/0000/*.root");      //X3872 (small)
            //tin->Add("/eos/user/h/hmarques/MC_PbPb23_X3872/MC_PbPb_X3872/prompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder/260322_030701/0000/*.root");     //PSI2S (small)
            tin->Add("/gstore/t3cms/store/user/hmarques/MC_PbPb23_X3872/MC_PbPb_X3872/prompt_X3872_to_Jpsi_Rho_phat5_Bfinder_large/260331_100051/0000/*.root");      //X3872 (lip, large)
            tin->Add("/gstore/t3cms/store/user/hmarques/MC_PbPb23_X3872/MC_PbPb_X3872/prompt_PSI2S_to_Jpsi_pipi_phat5_Bfinder_large/260331_095521/0000/*.root");     //PSI2S (lip, large)
        
        }    
        else{}
    }

    std::cout << "Total files added: "      << tin->GetNtrees() << std::endl;
    std::cout << "Total entries in chain: " << tin->GetEntries() << std::endl;

    // --------------------------------------------------
    // Output file
    // --------------------------------------------------
    TString datatype  = "";
    if(isMC){datatype = "MC";}
    else    {datatype = "DATA";}

    TString outputFile = "flat_" + treename + "_" + SYSTEM + "_" + datatype + specCASES + ".root";  //
    TFile *fout = new TFile(outputFile, "RECREATE");
    TTree *tout = new TTree(treename, "Flattened tree");

    // --------------------------------------------------
    // Prepare gen-level input chain and output tree
    // --------------------------------------------------
    TChain *tgen = new TChain("Bfinder/ntGen");
    // Reuse exactly the same file list used by the reco chain
    TObjArray *fileList = tin->GetListOfFiles();
    if (fileList) {
        TIter next(fileList);
        TObject *obj = nullptr;
        while ((obj = next())) {
            const char *fname = obj->GetTitle();
            tgen->Add(fname);
        }
    }

    TTree *tgenOut = nullptr;

    // --------------------------------------------------
    // Input reconstructed branches 
    // --------------------------------------------------
    const Int_t MAXCAND = 25000;
    Int_t Bsize, nSelectedChargedTracks, CentBin;
    Float_t Bmass[MAXCAND], Bpt[MAXCAND], By[MAXCAND], Bgen[MAXCAND];
    Bool_t Bmu1isTriggered[MAXCAND];
    Bool_t Bmu2isTriggered[MAXCAND];
    Bool_t Bmu1SoftMuID[MAXCAND];
    Bool_t Bmu2SoftMuID[MAXCAND];
    Bool_t Bmu1HybridSoftMuID[MAXCAND];
    Bool_t Bmu2HybridSoftMuID[MAXCAND];
    Bool_t Bmu1isAcc[MAXCAND];
    Bool_t Bmu2isAcc[MAXCAND];
    Bool_t  Btrk1highPurity[MAXCAND];
    Bool_t  Btrk2highPurity[MAXCAND];
    Float_t Btrk1Pt[MAXCAND];
    Float_t Btrk1Eta[MAXCAND];
    Float_t Btrk1PtErr[MAXCAND];
    Float_t Btrk1Chi2ndf[MAXCAND];
    Float_t Btrk2Pt[MAXCAND];
    Float_t Btrk2Eta[MAXCAND];
    Float_t Btrk2PtErr[MAXCAND];
    Float_t Btrk2Chi2ndf[MAXCAND];
    Float_t Btktkmass[MAXCAND];
    Float_t Bujmass[MAXCAND];
    Float_t BujvProb[MAXCAND];
    Float_t Bchi2Prob[MAXCAND];
    Float_t BtrkPtimb[MAXCAND];
    Float_t Btrk1dR[MAXCAND];
    Float_t Btrk2dR[MAXCAND];
    Float_t Btktkpt[MAXCAND];
    Float_t BtktkvProb[MAXCAND];
    Float_t BQvalue[MAXCAND];
    Float_t Bnorm_svpvDistance_2D[MAXCAND];
    Float_t BsvpvDistance_2D[MAXCAND];
    Float_t Balpha[MAXCAND];
    Float_t Bdtheta[MAXCAND];
    Float_t Bcos_dtheta[MAXCAND];
    Float_t Bnorm_trk1Dxy[MAXCAND];
    Float_t Bnorm_trk2Dxy[MAXCAND];
    Float_t Btrk1nPixelLayer[MAXCAND];
    Float_t Btrk1nStripLayer[MAXCAND];
    Float_t Btrk2nPixelLayer[MAXCAND];  
    Float_t Btrk2nStripLayer[MAXCAND];
    Float_t BLxy[MAXCAND];

    tin->SetBranchAddress("Bsize", &Bsize);
    tin->SetBranchAddress("nSelectedChargedTracks", &nSelectedChargedTracks);
    tin->SetBranchAddress("CentBin", &CentBin);
    tin->SetBranchAddress("Bmass", Bmass);
    tin->SetBranchAddress("Bpt", Bpt);
    tin->SetBranchAddress("By", By);
    tin->SetBranchAddress("Bgen", Bgen);
    tin->SetBranchAddress("Bmu1isTriggered", Bmu1isTriggered);
    tin->SetBranchAddress("Bmu2isTriggered", Bmu2isTriggered);
    tin->SetBranchAddress("Bmu1SoftMuID", Bmu1SoftMuID);
    tin->SetBranchAddress("Bmu2SoftMuID", Bmu2SoftMuID);
    tin->SetBranchAddress("Bmu1HybridSoftMuID", Bmu1HybridSoftMuID);
    tin->SetBranchAddress("Bmu2HybridSoftMuID", Bmu2HybridSoftMuID);
    tin->SetBranchAddress("Bmu1isAcc", Bmu1isAcc);
    tin->SetBranchAddress("Bmu2isAcc", Bmu2isAcc);  
    tin->SetBranchAddress("Btrk1highPurity", Btrk1highPurity);
    tin->SetBranchAddress("Btrk2highPurity", Btrk2highPurity);
    tin->SetBranchAddress("Btrk1Pt", Btrk1Pt);
    tin->SetBranchAddress("Btrk1Eta", Btrk1Eta);
    tin->SetBranchAddress("Btrk1PtErr", Btrk1PtErr);
    tin->SetBranchAddress("Btrk1Chi2ndf", Btrk1Chi2ndf);
    tin->SetBranchAddress("Btrk2Pt", Btrk2Pt);
    tin->SetBranchAddress("Btrk2Eta", Btrk2Eta);
    tin->SetBranchAddress("Btrk2PtErr", Btrk2PtErr);
    tin->SetBranchAddress("Btrk2Chi2ndf", Btrk2Chi2ndf);
    tin->SetBranchAddress("Btktkmass", Btktkmass);
    tin->SetBranchAddress("Bujmass", Bujmass);
    tin->SetBranchAddress("BujvProb", BujvProb);
    tin->SetBranchAddress("Bchi2Prob", Bchi2Prob);
    tin->SetBranchAddress("BtrkPtimb", BtrkPtimb);
    tin->SetBranchAddress("Btrk1dR", Btrk1dR);
    tin->SetBranchAddress("Btrk2dR", Btrk2dR);
    tin->SetBranchAddress("Btktkpt", Btktkpt);
    tin->SetBranchAddress("BtktkvProb", BtktkvProb);
    tin->SetBranchAddress("BQvalue", BQvalue);
    tin->SetBranchAddress("Bnorm_svpvDistance_2D", Bnorm_svpvDistance_2D);
    tin->SetBranchAddress("BsvpvDistance_2D", BsvpvDistance_2D);
    tin->SetBranchAddress("Balpha", Balpha);
    tin->SetBranchAddress("Bdtheta", Bdtheta);
    tin->SetBranchAddress("Bcos_dtheta", Bcos_dtheta);
    tin->SetBranchAddress("Bnorm_trk1Dxy", Bnorm_trk1Dxy);
    tin->SetBranchAddress("Bnorm_trk2Dxy", Bnorm_trk2Dxy);
    tin->SetBranchAddress("Btrk1nPixelLayer", Btrk1nPixelLayer);
    tin->SetBranchAddress("Btrk1nStripLayer", Btrk1nStripLayer);
    tin->SetBranchAddress("Btrk2nPixelLayer", Btrk2nPixelLayer);  
    tin->SetBranchAddress("Btrk2nStripLayer", Btrk2nStripLayer);
    tin->SetBranchAddress("BLxy", BLxy);
    // --------------------------------------------------
    // Output branches (for Selection // Analysis)
    // --------------------------------------------------
    float Bmass_out, Bpt_out, By_out, Bchi2Prob_out, Btrk1dR_out, Btrk2dR_out, BtrkPtimb_out, Btktkpt_out, Bujmass_out, Bnorm_svpvDistance_2D_out;
    float BQvalue_out, Bnorm_trk1Dxy_out, Bnorm_trk2Dxy_out, Balpha_out, Bdtheta_out, Bcos_dtheta_out, Btktkmass_out, Btrk1Pt_out, Btrk2Pt_out;
    float Bgen_out, BsvpvDistance_2D_out, BtktkvProb_out, BLxy_out;
    int CentBin_out, nSelectedChargedTracks_out;
    Bool_t isX3872_out;

    tout->Branch("nSelectedChargedTracks", &nSelectedChargedTracks_out, "nSelectedChargedTracks/I");
    tout->Branch("CentBin", &CentBin_out, "CentBin/I");
    tout->Branch("Bmass", &Bmass_out, "Bmass/F");
    tout->Branch("Bpt", &Bpt_out, "Bpt/F");
    tout->Branch("By", &By_out, "By/F");
    tout->Branch("Bchi2Prob", &Bchi2Prob_out, "Bchi2Prob/F");
    tout->Branch("Btrk1dR", &Btrk1dR_out, "Btrk1dR/F");
    tout->Branch("Btrk2dR", &Btrk2dR_out, "Btrk2dR/F");
    tout->Branch("BtrkPtimb", &BtrkPtimb_out, "BtrkPtimb/F");
    tout->Branch("Btktkpt", &Btktkpt_out, "Btktkpt/F");
    tout->Branch("Bujmass", &Bujmass_out, "Bujmass/F");
    tout->Branch("Bnorm_svpvDistance_2D", &Bnorm_svpvDistance_2D_out, "Bnorm_svpvDistance_2D/F");
    tout->Branch("BsvpvDistance_2D", &BsvpvDistance_2D_out, "BsvpvDistance_2D/F");
    tout->Branch("BQvalue", &BQvalue_out, "BQvalue/F");
    tout->Branch("Bnorm_trk1Dxy", &Bnorm_trk1Dxy_out, "Bnorm_trk1Dxy/F");
    tout->Branch("Bnorm_trk2Dxy", &Bnorm_trk2Dxy_out, "Bnorm_trk2Dxy/F");
    tout->Branch("Balpha", &Balpha_out, "Balpha/F");
    tout->Branch("Bdtheta", &Bdtheta_out, "Bdtheta/F");
    tout->Branch("Bcos_dtheta", &Bcos_dtheta_out, "Bcos_dtheta/F");
    tout->Branch("Btktkmass", &Btktkmass_out, "Btktkmass/F");
    tout->Branch("Btrk1Pt", &Btrk1Pt_out, "Btrk1Pt/F");
    tout->Branch("Btrk2Pt", &Btrk2Pt_out, "Btrk2Pt/F");
    tout->Branch("BtktkvProb", &BtktkvProb_out, "BtktkvProb/F");
    tout->Branch("BLxy", &BLxy_out, "BLxy/F");
    if (isMC){
        tout->Branch("Bgen", &Bgen_out, "Bgen/F");
        if (treename == "ntmix") tout->Branch("isX3872", &isX3872_out, "isX3872/O");
    }

    // --------------------------------------------------
    // Event loop
    // --------------------------------------------------
    Long64_t nentries = tin->GetEntries();
    std::cout << "Processing " << nentries << " entries..." << std::endl;
    for(Long64_t ev=0; ev<nentries; ++ev)
    {
        // keep track of progress
        if (ev % 100000 == 0) {
            std::cout << "Processing event " << ev << " / " << nentries 
                      << " (" << (100.0*ev/nentries) << "%)" << std::endl;
        }
        tin->GetEntry(ev);

        // identify source by filename
        TString fname = tin->GetCurrentFile()->GetName();
        isX3872_out = fname.Contains("X3872_") ? 1 : 0;
    
        // Event-level variables
        CentBin_out = CentBin;
        nSelectedChargedTracks_out = nSelectedChargedTracks;

        // --------------------------------------------------
        // Candidate loop
        // --------------------------------------------------
        for(int i=0; i<Bsize; ++i)
        {
            // Pre-Selection (quality) Cuts
            // DATA is already filtered (except for the cuts marked with <-------- ** ); 
            // MC is not.
            if (PreSel_cuts){
                // Muons
                // Muon Trigger matching   <-------- **
                if(isPP && !(Bmu1isTriggered[i] && Bmu2isTriggered[i] )) continue;
                // Soft vs HybridSoft Muon <-------- **
                if((isPP   && !(Bmu1SoftMuID[i] && Bmu2SoftMuID[i] )) ||
                   (isPbPb && !(Bmu1HybridSoftMuID[i] && Bmu2HybridSoftMuID[i]))) continue;
                // Muon Acceptance
                if( !(Bmu1isAcc[i] && Bmu2isAcc[i])) continue;
                // diMuon system 
                if( !(abs(Bujmass[i] - JPSI_MASS) < 0.15)) continue;
                if( !(BujvProb[i] > 0.01)) continue;

                // Tracks
                // single Track channel
                if(abs(Btrk1Eta[i]) > 2.4) continue;
                if((isPP && !(Btrk1Pt[i] > 0.5)) || (isPbPb && !(Btrk1Pt[i] > 0.9)) ) continue;
                if(((Btrk1PtErr[i] / Btrk1Pt[i]) >= 0.1)) continue;
                if(((Btrk1nPixelLayer[i] + Btrk1nStripLayer[i]) <= 10)) continue;
                if(((Btrk1Chi2ndf[i]/(Btrk1nPixelLayer[i] + Btrk1nStripLayer[i])) >= 0.18)) continue;
                if(!Btrk1highPurity[i]) continue;
                // diTrack channel 
                if (treename !=  "ntKp"){ 
                    if(abs(Btrk2Eta[i]) > 2.4) continue;
                    if((isPP && !(Btrk2Pt[i] > 0.5)) || (isPbPb && !(Btrk2Pt[i] > 0.9)) ) continue;
                    if(((Btrk2PtErr[i] / Btrk2Pt[i]) >= 0.1)) continue;
                    if(((Btrk2nPixelLayer[i] + Btrk2nStripLayer[i]) <= 10)) continue;
                    if(((Btrk2Chi2ndf[i]/(Btrk2nPixelLayer[i] + Btrk2nStripLayer[i])) >= 0.18)) continue;
                    if(!Btrk2highPurity[i]) continue;
                    //diTrack system
                    //if (treename == "ntphi"   && !(abs(Btktkmass[i] - PHI_MASS)  <0.015)) continue;
                    //if (treename == "ntKstar" && !(abs(Btktkmass[i] - KSTAR_MASS)<0.150)) continue;
                }

                // Candidate
                if(Bchi2Prob[i] < 0.005) continue;
                if(treename == "ntmix" && Bpt[i]<4) continue; 
                else if (Bpt[i]<1) continue;
            }

            // Fiducial Region
            if (Fid_region){
                if ( ( (Bpt[i] < 5) && abs(By[i]) > 2.4) ) continue; 
                if ( (treename != "ntmix") && !( (Bpt[i] < 10  && abs(By[i]) > 1.5) || (Bpt[i] > 10  && Bpt[i] < 50) ) ) continue;
                if ( (treename == "ntmix") && !( (Bpt[i] < 7.5 && abs(By[i]) > 1.4) || (Bpt[i] > 7.5 && Bpt[i] < 50) ) ) continue;
            }

            // keep MC signal only
            Bgen_out = Bgen[i];
            if( isMC && !(Bgen_out==23333 || Bgen_out==24333 || Bgen_out==23433 || Bgen_out==24433 || Bgen_out==41000 )) continue;

            Bmass_out   = Bmass[i];
            Bpt_out     = Bpt[i];
            By_out      = By[i];
            Bchi2Prob_out = Bchi2Prob[i];
            Btrk1dR_out = Btrk1dR[i];
            Btrk2dR_out = Btrk2dR[i];
            BtrkPtimb_out = BtrkPtimb[i];
            Btktkpt_out   = Btktkpt[i];
            Bujmass_out   = Bujmass[i];
            Bnorm_svpvDistance_2D_out = Bnorm_svpvDistance_2D[i];
            BsvpvDistance_2D_out = BsvpvDistance_2D[i];
            BQvalue_out = BQvalue[i];
            Bnorm_trk1Dxy_out = Bnorm_trk1Dxy[i];
            Bnorm_trk2Dxy_out = Bnorm_trk2Dxy[i];
            Balpha_out  = Balpha[i];
            BtktkvProb_out = BtktkvProb[i];
            Bdtheta_out = Bdtheta[i];
            Bcos_dtheta_out = Bcos_dtheta[i];
            Btktkmass_out = Btktkmass[i];
            Btrk1Pt_out   = Btrk1Pt[i];
            Btrk2Pt_out   = Btrk2Pt[i];
            BLxy_out      = BLxy[i];

            if(!std::isfinite(Bmass_out) || !std::isfinite(Bpt_out) || !std::isfinite(By_out) || !std::isfinite(Bnorm_trk1Dxy_out) ||
            !std::isfinite(CentBin_out) || !std::isfinite(Bchi2Prob_out) || !std::isfinite(Btrk1dR_out) || !std::isfinite(Bnorm_svpvDistance_2D_out)) continue;
            tout->Fill();
        }
    }

    // --------------------------------------------------
    // Copy ntGen tree (flat, one row per candidate in mass window)
    // --------------------------------------------------
    if (isMC && KeepGenLv) {
        fout->cd();

        // Gen-level input branches
        Int_t Gsize;
        Int_t GpdgId[MAXCAND];
        Float_t Gmu1eta[MAXCAND];
        Float_t Gmu1pt[MAXCAND];
        Float_t Gmu2eta[MAXCAND];
        Float_t Gmu2pt[MAXCAND];
        Float_t Gtk1pt[MAXCAND];
        Float_t Gtk1eta[MAXCAND];
        Float_t Gtk2pt[MAXCAND];
        Float_t Gtk2eta[MAXCAND];
        Float_t Gpt[MAXCAND];
        Float_t Gy[MAXCAND];
        Int_t GisSignal[MAXCAND];

        tgen->SetBranchAddress("Gsize", &Gsize);
        tgen->SetBranchAddress("GisSignal", GisSignal);
        tgen->SetBranchAddress("Gmu1eta", Gmu1eta);
        tgen->SetBranchAddress("Gmu1pt", Gmu1pt);
        tgen->SetBranchAddress("Gmu2eta", Gmu2eta);
        tgen->SetBranchAddress("Gmu2pt", Gmu2pt);
        tgen->SetBranchAddress("Gtk1pt", Gtk1pt);
        tgen->SetBranchAddress("Gtk1eta", Gtk1eta);
        tgen->SetBranchAddress("Gtk2pt", Gtk2pt);
        tgen->SetBranchAddress("Gtk2eta", Gtk2eta);
        tgen->SetBranchAddress("GpdgId", GpdgId);
        tgen->SetBranchAddress("Gpt", Gpt);
        tgen->SetBranchAddress("Gy", Gy);

        // Output tree (flat)
        tgenOut = new TTree("ntGen", "Gen-level (flat, filtered)");
        Float_t Gmu1eta_out, Gmu1pt_out, Gmu2eta_out, Gmu2pt_out;
        Float_t Gtk1pt_out, Gtk1eta_out, Gtk2pt_out, Gtk2eta_out, Gpt_out, Gy_out;
        Int_t GpdgId_out;
        tgenOut->Branch("Gmu1eta", &Gmu1eta_out, "Gmu1eta/F");
        tgenOut->Branch("Gmu1pt", &Gmu1pt_out, "Gmu1pt/F");
        tgenOut->Branch("Gmu2eta", &Gmu2eta_out, "Gmu2eta/F");
        tgenOut->Branch("Gmu2pt", &Gmu2pt_out, "Gmu2pt/F");
        tgenOut->Branch("Gtk1pt", &Gtk1pt_out, "Gtk1pt/F");
        tgenOut->Branch("Gtk1eta", &Gtk1eta_out, "Gtk1eta/F");
        tgenOut->Branch("Gtk2pt", &Gtk2pt_out, "Gtk2pt/F");
        tgenOut->Branch("Gtk2eta", &Gtk2eta_out, "Gtk2eta/F");
        tgenOut->Branch("GpdgId", &GpdgId_out, "GpdgId/I");
        tgenOut->Branch("Gpt", &Gpt_out, "Gpt/F");
        tgenOut->Branch("Gy", &Gy_out, "Gy/F");

        const Long64_t ngen = tgen->GetEntries();
        for (Long64_t ev=0; ev<ngen; ++ev) {
            tgen->GetEntry(ev);
            for (int i=0; i<Gsize; ++i) {

                if (treename == "ntmix" && GisSignal[i] != 7) continue;              // keep only X(3872), both prompt and non-prompt
                else if (treename == "ntphi"   && abs(GpdgId[i]) != 531) continue;  
                else if (treename == "ntKstar" && abs(GpdgId[i]) != 511) continue;  
                else if (treename == "ntKp"    && abs(GpdgId[i]) != 521) continue;

                // Fiducial Region
                if (Fid_region){
                    if ( !( (Gpt[i] > 5) && abs(Gy[i]) < 2.4) ) continue; 
                    if ( (treename != "ntmix") && !( (Gpt[i] < 10  && abs(Gy[i]) > 1.5) || (Gpt[i] > 10  && Gpt[i] < 50) ) ) continue;
                    if ( (treename == "ntmix") && !( (Gpt[i] < 7.5 && abs(Gy[i]) > 1.4) || (Gpt[i] > 7.5 && Gpt[i] < 50) ) ) continue;
                }

                GpdgId_out  = GpdgId[i];
                Gmu1eta_out = Gmu1eta[i];
                Gmu1pt_out  = Gmu1pt[i];
                Gmu2eta_out = Gmu2eta[i];
                Gmu2pt_out  = Gmu2pt[i];
                Gtk1pt_out  = Gtk1pt[i];
                Gtk1eta_out = Gtk1eta[i];
                Gtk2pt_out  = Gtk2pt[i];
                Gtk2eta_out = Gtk2eta[i];
                Gpt_out = Gpt[i];
                Gy_out  = Gy[i];
                tgenOut->Fill();

            }
        }
    }

    std::cout << "Output tree has " << tout->GetEntries() << " entries" << std::endl;
    // Write trees (only two top-level trees)
    tout->Write();
    if (isMC && tgenOut) tgenOut->Write();
    fout->Close();
    delete tin;
    std::cout << "Done & Saved -> " << outputFile << "\n";
}