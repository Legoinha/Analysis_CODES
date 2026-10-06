#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TVector2.h>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "../plotER/aux/masses.h"

void Flat_TREEs( TString FILEIN="", TString NUN="", TString TREENAME="ntmix", TString SYSTEM="ppRef", TString KIND="MC", TString PARTICLE="", TString PVSNP="", TString MCNORMFILE="config/mc_normalization.csv")
{
    bool isMC = (KIND == "MC");
    bool Fid_region  = true;       // apply fiducial region selection
    TString specCASES = PARTICLE + PVSNP;   // explicit case tag passed by the runner after filelist preselection

    std::cout << "Flattening tree: " << TREENAME << " , System: " << SYSTEM << " , kind: " << KIND << std::endl;
    bool isPP     = (SYSTEM == "ppRef") ;
    bool isPbPb   = (SYSTEM == "PbPb23" || SYSTEM == "PbPb24" || SYSTEM == "PbPb25");

    // --------------------------------------------------
    // Input files: one ROOT file (or wildcard) per line of FILEIN; lines starting with # are skipped
    // --------------------------------------------------
    TChain *tin = new TChain("Bfinder/" + TREENAME);
    std::ifstream fin(FILEIN.Data());
    std::string line;
    while (std::getline(fin, line))
        if (line[0] != '#') tin->Add(line.c_str());
    std::cout << "Total files added: "      << tin->GetNtrees() << std::endl;
    std::cout << "Total entries in chain: " << tin->GetEntries() << std::endl;

    // Gen-level input chain, on exactly the same files as the reco chain
    TChain *tgen = new TChain("Bfinder/ntGen");
    for (TObject *file : *tin->GetListOfFiles()) tgen->Add(file->GetTitle());

    // --------------------------------------------------
    // MC pThat weights (config/README.md). The pThat-X campaigns of this system, tree, particle and
    // promptness are the rows of MCNORMFILE, sorted by threshold X. An event with generator pThat p
    // can come from every campaign j with X_j < p, so
    //     pThatreweight(p) = 1 / sum_{X_j < p} L_j,     L_j = n_gen_j / (xsec_j * filter_eff_j),
    // with n_gen_j the Bfinder/ntGen events of the input files of campaign j (matched by
    // path_pattern). A campaign without input files has n_gen = 0 and adds no luminosity.
    // pThats: the threshold X of the campaign of each input file, the same for reco and gen chains.
    // --------------------------------------------------
    struct Campaign { int pthat; TString pattern; double xsec, filterEff; Long64_t nGen = 0; int files = 0; };
    std::vector<Campaign> campaigns;
    std::vector<double> cumulativeLuminosity;   // sum of L_j over the campaigns up to each one
    std::vector<Int_t> campaignPthat;           // per input file
    if (isMC) {
        const TString particle = TREENAME == "ntKp" ? "bplus" : TREENAME == "ntphi" ? "bs" : TREENAME == "ntKstar" ? "bzero"
                               : PARTICLE.Contains("psi2s", TString::kIgnoreCase) ? "psi2s" : "x3872";
        const TString promptness = TREENAME != "ntmix" ? "inclusive"
                                 : PVSNP.Contains("nonprompt", TString::kIgnoreCase) ? "nonprompt" : "prompt";
        // Columns: system,tree,particle,promptness,pthat,path_pattern,xsec_pb,filter_eff. The header and
        // the # comments never match the system.
        std::ifstream table(MCNORMFILE.Data());
        while (std::getline(table, line)) {
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream row(line);
            TString rowSystem, rowTree, rowParticle, rowPromptness;
            Campaign campaign;
            row >> rowSystem >> rowTree >> rowParticle >> rowPromptness >> campaign.pthat >> campaign.pattern
                >> campaign.xsec >> campaign.filterEff;
            if (rowSystem.EqualTo(SYSTEM, TString::kIgnoreCase) && rowTree.EqualTo(TREENAME, TString::kIgnoreCase) &&
                rowParticle == particle && rowPromptness == promptness)
                campaigns.push_back(campaign);
        }
        std::sort(campaigns.begin(), campaigns.end(), [](const Campaign &a, const Campaign &b) { return a.pthat < b.pthat; });

        tgen->GetEntries();   // opens every file and fills the per-file entry offsets
        for (Int_t i = 0; i < tgen->GetNtrees(); ++i) {
            const TString path = tgen->GetListOfFiles()->At(i)->GetTitle();
            Campaign &campaign = campaigns.at(std::find_if(campaigns.begin(), campaigns.end(), [&path](const Campaign &c) {
                return path.Contains(c.pattern, TString::kIgnoreCase); }) - campaigns.begin());
            campaign.nGen += tgen->GetTreeOffset()[i + 1] - tgen->GetTreeOffset()[i];
            ++campaign.files;
            campaignPthat.push_back(campaign.pthat);
        }
        double luminosity = 0.;
        for (const auto &campaign : campaigns) {
            luminosity += campaign.nGen / (campaign.xsec * campaign.filterEff);
            cumulativeLuminosity.push_back(luminosity);
        }

        std::cout << "MC normalization table: " << MCNORMFILE << std::endl;
        std::cout << "pThat weights, w = 1 / sum of L_j over samples with pThat_j below the event pThat:" << std::endl;
        for (std::size_t k = 0; k < campaigns.size(); ++k) {
            std::cout << "  pThat > " << std::setw(2) << campaigns[k].pthat
                      << "  files=" << std::setw(4) << campaigns[k].files
                      << "  xsec(pb)=" << campaigns[k].xsec
                      << "  filterEff=" << campaigns[k].filterEff
                      << "  nGen=" << campaigns[k].nGen
                      << "  L=" << campaigns[k].nGen / (campaigns[k].xsec * campaigns[k].filterEff)
                      << "  weight for events in this pThat range=" << std::setprecision(12)
                      << 1. / cumulativeLuminosity[k] << std::setprecision(6) << std::endl;
        }
    }

    // --------------------------------------------------
    // Output file
    // --------------------------------------------------
    TString outputFile = "flat_" + TREENAME + "_" + SYSTEM + "_" + KIND + specCASES + NUN + ".root";  //
    TFile *fout = new TFile(outputFile, "RECREATE");
    TTree *tout = new TTree(TREENAME + PARTICLE, "Flattened tree");
    // Do not leave intermediate autosaved cycles in chunk files. Besides making
    // the file larger, such cycles can retain references to a chunk's old path.
    tout->SetAutoSave(0);
    TTree *tevent = new TTree("ntEvent", "Event-level tree");
    tevent->SetAutoSave(0);

    TTree *tgenOut = nullptr;


    // --------------------------------------------------
    // Input reconstructed branches 
    // --------------------------------------------------
    const Int_t MAXCAND = 40000;
    Int_t RunNo, LumiNo, Bsize, CentBin;
    ULong64_t EvtNo;
    Float_t PVx, PVy, PVz, PVxE, PVyE, PVzE, PVnchi2;
    Float_t nChargedTracks, nChargedTracks_LOOSE, nChargedTracks_TIGHT;
    static Int_t Bindex[MAXCAND];
    static Float_t Bmass[MAXCAND], Bpt[MAXCAND], By[MAXCAND], Bgen[MAXCAND];
    static Float_t BvtxX[MAXCAND], BvtxY[MAXCAND];
    static Bool_t Bmu1isTriggered[MAXCAND];
    static Bool_t Bmu2isTriggered[MAXCAND];
    static Bool_t Bmu1SoftMuID[MAXCAND];
    static Bool_t Bmu2SoftMuID[MAXCAND];
    static Bool_t Bmu1HybridSoftMuID[MAXCAND];
    static Bool_t Bmu2HybridSoftMuID[MAXCAND];
    static Bool_t Bmu1isAcc[MAXCAND];
    static Bool_t Bmu2isAcc[MAXCAND];
    static Bool_t  Btrk1highPurity[MAXCAND];
    static Bool_t  Btrk2highPurity[MAXCAND];
    static Float_t Btrk1Pt[MAXCAND];
    static Float_t Btrk1Eta[MAXCAND];
    static Float_t Btrk1Phi[MAXCAND];
    static Float_t Btrk1PtErr[MAXCAND];
    static Float_t Btrk1Chi2ndf[MAXCAND];
    static Float_t Btrk2Pt[MAXCAND];
    static Float_t Btrk2Eta[MAXCAND];
    static Float_t Btrk2Phi[MAXCAND];
    static Float_t Btrk2PtErr[MAXCAND];
    static Float_t Btrk2Chi2ndf[MAXCAND];
    static Float_t Btktkmass[MAXCAND];
    static Float_t Bujmass[MAXCAND];
    static Float_t BujvProb[MAXCAND];
    static Float_t Bchi2Prob[MAXCAND];
    static Float_t BtrkPtimb[MAXCAND];
    static Float_t Btrk1dR[MAXCAND];
    static Float_t Btrk2dR[MAXCAND];
    static Float_t Btktkpt[MAXCAND];
    static Float_t BtktkvProb[MAXCAND];
    static Float_t BQvalue[MAXCAND];
    static Float_t Bnorm_svpvDistance_2D[MAXCAND];
    static Float_t BsvpvDistance_2D[MAXCAND];
    static Float_t BsvpvDisErr_2D[MAXCAND];
    static Float_t Balpha[MAXCAND];
    static Float_t Bdtheta[MAXCAND];
    static Float_t Bcos_dtheta[MAXCAND];
    static Float_t Bnorm_trk1Dxy[MAXCAND];
    static Float_t Bnorm_trk2Dxy[MAXCAND];
    static Float_t Btrk1nPixelLayer[MAXCAND];
    static Float_t Btrk1nStripLayer[MAXCAND];
    static Float_t Btrk2nPixelLayer[MAXCAND];  
    static Float_t Btrk2nStripLayer[MAXCAND];
    static Float_t BLxy[MAXCAND];
    static Float_t Bmu1pt[MAXCAND];
    static Float_t Bmu2pt[MAXCAND];
    static Float_t Bmu1eta[MAXCAND];
    static Float_t Bmu2eta[MAXCAND];
    static Float_t Bmu1phi[MAXCAND];
    static Float_t Bmu2phi[MAXCAND];
    static Float_t Bujpt[MAXCAND];
    static Float_t Bujeta[MAXCAND];
    static Float_t Bujphi[MAXCAND];
    static Float_t Bujlxy[MAXCAND];
    static Bool_t BdiTrackFitValid[MAXCAND];
    static Float_t Btrk1Dz1[MAXCAND];
    static Float_t Btrk2Dz1[MAXCAND];
    static Float_t Btrk1DzError1[MAXCAND];
    static Float_t Btrk2DzError1[MAXCAND];
    static Float_t Btrk1Dxy1[MAXCAND];
    static Float_t Btrk2Dxy1[MAXCAND];
    static Float_t Btrk1DxyError1[MAXCAND];
    static Float_t Btrk2DxyError1[MAXCAND];
    static Float_t Btktketa[MAXCAND];
    static Float_t Btktkphi[MAXCAND];
    static Float_t Btktky[MAXCAND];
    static Float_t Bdoubletpt[MAXCAND];
    static Float_t Bdoubleteta[MAXCAND];
    static Float_t Bdoubletphi[MAXCAND];
    static Float_t Bdoublety[MAXCAND];
    static Float_t Bnorm_trk1Dz[MAXCAND];
    static Float_t Bnorm_trk2Dz[MAXCAND];

    tin->SetBranchAddress("Bsize", &Bsize);
    tin->SetBranchAddress("RunNo", &RunNo);
    tin->SetBranchAddress("EvtNo", &EvtNo);
    tin->SetBranchAddress("LumiNo", &LumiNo);
    tin->SetBranchAddress("PVx", &PVx);
    tin->SetBranchAddress("PVy", &PVy);
    tin->SetBranchAddress("PVz", &PVz);
    tin->SetBranchAddress("PVxE", &PVxE);
    tin->SetBranchAddress("PVyE", &PVyE);
    tin->SetBranchAddress("PVzE", &PVzE);
    tin->SetBranchAddress("PVnchi2", &PVnchi2);
    tin->SetBranchAddress("nChargedTracks", &nChargedTracks);
    tin->SetBranchAddress("nChargedTracks_LOOSE", &nChargedTracks_LOOSE);
    tin->SetBranchAddress("nChargedTracks_TIGHT", &nChargedTracks_TIGHT);
    tin->SetBranchAddress("CentBin", &CentBin);
    Float_t pthat = 0.;
    if (isMC) tin->SetBranchAddress("pthat", &pthat);
    tin->SetBranchAddress("Bindex", Bindex);
    tin->SetBranchAddress("Bmass", Bmass);
    tin->SetBranchAddress("BvtxX", BvtxX);
    tin->SetBranchAddress("BvtxY", BvtxY);
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
    tin->SetBranchAddress("Btrk1Phi", Btrk1Phi);
    tin->SetBranchAddress("Btrk1PtErr", Btrk1PtErr);
    tin->SetBranchAddress("Btrk1Chi2ndf", Btrk1Chi2ndf);
    tin->SetBranchAddress("Btrk2Pt", Btrk2Pt);
    tin->SetBranchAddress("Btrk2Eta", Btrk2Eta);
    tin->SetBranchAddress("Btrk2Phi", Btrk2Phi);
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
    tin->SetBranchAddress("BsvpvDisErr_2D", BsvpvDisErr_2D);
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
    tin->SetBranchAddress("Bmu1pt", Bmu1pt);
    tin->SetBranchAddress("Bmu2pt", Bmu2pt);
    tin->SetBranchAddress("Bmu1eta", Bmu1eta);
    tin->SetBranchAddress("Bmu2eta", Bmu2eta);
    tin->SetBranchAddress("Bmu1phi", Bmu1phi);
    tin->SetBranchAddress("Bmu2phi", Bmu2phi);
    tin->SetBranchAddress("Bujpt", Bujpt);
    tin->SetBranchAddress("Bujeta", Bujeta);
    tin->SetBranchAddress("Bujphi", Bujphi);
    tin->SetBranchAddress("Bujlxy", Bujlxy);
    tin->SetBranchAddress("BdiTrackFitValid", BdiTrackFitValid);
    tin->SetBranchAddress("Btrk1Dz1", Btrk1Dz1);
    tin->SetBranchAddress("Btrk2Dz1", Btrk2Dz1);
    tin->SetBranchAddress("Btrk1DzError1", Btrk1DzError1);
    tin->SetBranchAddress("Btrk2DzError1", Btrk2DzError1);
    tin->SetBranchAddress("Btrk1Dxy1", Btrk1Dxy1);
    tin->SetBranchAddress("Btrk2Dxy1", Btrk2Dxy1);
    tin->SetBranchAddress("Btrk1DxyError1", Btrk1DxyError1);
    tin->SetBranchAddress("Btrk2DxyError1", Btrk2DxyError1);
    tin->SetBranchAddress("Btktketa", Btktketa);
    tin->SetBranchAddress("Btktkphi", Btktkphi);
    tin->SetBranchAddress("Btktky", Btktky);
    tin->SetBranchAddress("Bdoubletpt", Bdoubletpt);
    tin->SetBranchAddress("Bdoubleteta", Bdoubleteta);
    tin->SetBranchAddress("Bdoubletphi", Bdoubletphi);
    tin->SetBranchAddress("Bdoublety", Bdoublety);
    tin->SetBranchAddress("Bnorm_trk1Dz", Bnorm_trk1Dz);
    tin->SetBranchAddress("Bnorm_trk2Dz", Bnorm_trk2Dz);
    // --------------------------------------------------
    // Output branches (for Selection // Analysis)
    // --------------------------------------------------
    int RunNo_out, LumiNo_out, Bsize_out, Bindex_out, candidateIndex_out;
    ULong64_t EvtNo_out;
    Long64_t eventIndex_out;
    int nCandPassPreselection_out, nCandPassFiducial_out, nCandWritten_out;
    float PVx_out, PVy_out, PVz_out, PVxE_out, PVyE_out, PVzE_out, PVnchi2_out;
    float nChargedTracks_out, nChargedTracks_LOOSE_out, nChargedTracks_TIGHT_out;
    float Bmass_out, Bpt_out, By_out, Bchi2Prob_out, Btrk1dR_out, Btrk2dR_out, BtrkPtimb_out, Btktkpt_out, Bujmass_out, Bnorm_svpvDistance_2D_out;
    float BQvalue_out, Bnorm_trk1Dxy_out, Bnorm_trk2Dxy_out, Balpha_out, Bdtheta_out, Bcos_dtheta_out, Btktkmass_out, Btrk1Pt_out, Btrk2Pt_out;
    float Bgen_out, BsvpvDistance_2D_out, BtktkvProb_out, BLxy_out;
    float BvtxX_out, BvtxY_out, BsvpvDisErr_2D_out;
    float Btrk1Eta_out, Btrk2Eta_out, Btrk1Phi_out, Btrk2Phi_out, Btrk1PtErr_out, Btrk2PtErr_out;
    float BujvProb_out, Bmu1pt_out, Bmu2pt_out, Bnorm_trk1Dz_out, Bnorm_trk2Dz_out;
    float Bmu1eta_out, Bmu2eta_out, Bmu1phi_out, Bmu2phi_out;
    float Bujpt_out, Bujeta_out, Bujphi_out, Bujlxy_out;
    float Btrk1Dz1_out, Btrk2Dz1_out, Btrk1DzError1_out, Btrk2DzError1_out;
    float Btrk1Dxy1_out, Btrk2Dxy1_out, Btrk1DxyError1_out, Btrk2DxyError1_out;
    float Btktketa_out, Btktkphi_out, Btktky_out;
    float Bdoubletpt_out, Bdoubleteta_out, Bdoubletphi_out, Bdoublety_out;
    float Btrk1Chi2ndf_out, Btrk2Chi2ndf_out;
    float Btrk1nPixelLayer_out, Btrk2nPixelLayer_out, Btrk1nStripLayer_out, Btrk2nStripLayer_out;
    bool BdiTrackFitValid_out, Bmu1isTriggered_out, Bmu2isTriggered_out;
    bool Bmu1SoftMuID_out, Bmu2SoftMuID_out, Bmu1HybridSoftMuID_out, Bmu2HybridSoftMuID_out;
    bool Bmu1isAcc_out, Bmu2isAcc_out, Btrk1highPurity_out, Btrk2highPurity_out;
    int CentBin_out;
    // The trk1/trk2 and mu1/mu2 labels are not ordered, so these symmetric versions are kept
    // too: pT-ordered values, pT fractions of the B, and the largest and smallest Delta R and
    // |eta|. Muon Delta R is to the J/psi, like Btrk1dR and Btrk2dR. ntKp has one track only.
    const bool twoTrack = (TREENAME != "ntKp");
    float BtrkLeadPt_out, BtrkSubPt_out, BtrkLeadPtFrac_out, BtrkSubPtFrac_out, BtrkMaxdR_out, BtrkMindR_out, BtrkMaxAbsEta_out;
    float BmuLeadPt_out, BmuSubPt_out, BmuMaxdR_out, BmuMindR_out, BmuMaxAbsEta_out;
    Double_t pThatreweight_out = 1.;
    Int_t pThats_out = 0;

    tout->Branch("RunNo", &RunNo_out, "RunNo/I");
    tout->Branch("EvtNo", &EvtNo_out, "EvtNo/l");
    tout->Branch("LumiNo", &LumiNo_out, "LumiNo/I");
    tout->Branch("eventIndex", &eventIndex_out, "eventIndex/L");
    tout->Branch("candidateIndex", &candidateIndex_out, "candidateIndex/I");
    tout->Branch("Bindex", &Bindex_out, "Bindex/I");
    tout->Branch("Bsize", &Bsize_out, "Bsize/I");
    tout->Branch("PVx", &PVx_out, "PVx/F");
    tout->Branch("PVy", &PVy_out, "PVy/F");
    tout->Branch("PVz", &PVz_out, "PVz/F");
    tout->Branch("PVxE", &PVxE_out, "PVxE/F");
    tout->Branch("PVyE", &PVyE_out, "PVyE/F");
    tout->Branch("PVzE", &PVzE_out, "PVzE/F");
    tout->Branch("PVnchi2", &PVnchi2_out, "PVnchi2/F");
    tout->Branch("nChargedTracks", &nChargedTracks_out, "nChargedTracks/F");
    tout->Branch("nChargedTracks_LOOSE", &nChargedTracks_LOOSE_out, "nChargedTracks_LOOSE/F");
    tout->Branch("nChargedTracks_TIGHT", &nChargedTracks_TIGHT_out, "nChargedTracks_TIGHT/F");
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
    tout->Branch("BujvProb", &BujvProb_out, "BujvProb/F");
    tout->Branch("Bnorm_svpvDistance_2D", &Bnorm_svpvDistance_2D_out, "Bnorm_svpvDistance_2D/F");
    tout->Branch("BsvpvDistance_2D", &BsvpvDistance_2D_out, "BsvpvDistance_2D/F");
    tout->Branch("BsvpvDisErr_2D", &BsvpvDisErr_2D_out, "BsvpvDisErr_2D/F");
    tout->Branch("BQvalue", &BQvalue_out, "BQvalue/F");
    tout->Branch("Bnorm_trk1Dxy", &Bnorm_trk1Dxy_out, "Bnorm_trk1Dxy/F");
    tout->Branch("Bnorm_trk2Dxy", &Bnorm_trk2Dxy_out, "Bnorm_trk2Dxy/F");
    tout->Branch("Balpha", &Balpha_out, "Balpha/F");
    tout->Branch("Bdtheta", &Bdtheta_out, "Bdtheta/F");
    tout->Branch("Bcos_dtheta", &Bcos_dtheta_out, "Bcos_dtheta/F");
    tout->Branch("Btktkmass", &Btktkmass_out, "Btktkmass/F");
    tout->Branch("Btrk1Pt", &Btrk1Pt_out, "Btrk1Pt/F");
    tout->Branch("Btrk2Pt", &Btrk2Pt_out, "Btrk2Pt/F");
    tout->Branch("Btrk1Eta", &Btrk1Eta_out, "Btrk1Eta/F");
    tout->Branch("Btrk2Eta", &Btrk2Eta_out, "Btrk2Eta/F");
    tout->Branch("Btrk1Phi", &Btrk1Phi_out, "Btrk1Phi/F");
    tout->Branch("Btrk2Phi", &Btrk2Phi_out, "Btrk2Phi/F");
    tout->Branch("Btrk1PtErr", &Btrk1PtErr_out, "Btrk1PtErr/F");
    tout->Branch("Btrk2PtErr", &Btrk2PtErr_out, "Btrk2PtErr/F");
    tout->Branch("Btrk1Chi2ndf", &Btrk1Chi2ndf_out, "Btrk1Chi2ndf/F");
    tout->Branch("Btrk2Chi2ndf", &Btrk2Chi2ndf_out, "Btrk2Chi2ndf/F");
    tout->Branch("Btrk1nPixelLayer", &Btrk1nPixelLayer_out, "Btrk1nPixelLayer/F");
    tout->Branch("Btrk2nPixelLayer", &Btrk2nPixelLayer_out, "Btrk2nPixelLayer/F");
    tout->Branch("Btrk1nStripLayer", &Btrk1nStripLayer_out, "Btrk1nStripLayer/F");
    tout->Branch("Btrk2nStripLayer", &Btrk2nStripLayer_out, "Btrk2nStripLayer/F");
    tout->Branch("Btrk1highPurity", &Btrk1highPurity_out, "Btrk1highPurity/O");
    tout->Branch("Btrk2highPurity", &Btrk2highPurity_out, "Btrk2highPurity/O");
    tout->Branch("BtktkvProb", &BtktkvProb_out, "BtktkvProb/F");
    tout->Branch("BLxy", &BLxy_out, "BLxy/F");
    tout->Branch("BvtxX", &BvtxX_out, "BvtxX/F");
    tout->Branch("BvtxY", &BvtxY_out, "BvtxY/F");
    tout->Branch("Bmu1pt", &Bmu1pt_out, "Bmu1pt/F");
    tout->Branch("Bmu2pt", &Bmu2pt_out, "Bmu2pt/F");
    tout->Branch("Bmu1eta", &Bmu1eta_out, "Bmu1eta/F");
    tout->Branch("Bmu2eta", &Bmu2eta_out, "Bmu2eta/F");
    tout->Branch("Bmu1phi", &Bmu1phi_out, "Bmu1phi/F");
    tout->Branch("Bmu2phi", &Bmu2phi_out, "Bmu2phi/F");
    tout->Branch("Bmu1isTriggered", &Bmu1isTriggered_out, "Bmu1isTriggered/O");
    tout->Branch("Bmu2isTriggered", &Bmu2isTriggered_out, "Bmu2isTriggered/O");
    tout->Branch("Bmu1SoftMuID", &Bmu1SoftMuID_out, "Bmu1SoftMuID/O");
    tout->Branch("Bmu2SoftMuID", &Bmu2SoftMuID_out, "Bmu2SoftMuID/O");
    tout->Branch("Bmu1HybridSoftMuID", &Bmu1HybridSoftMuID_out, "Bmu1HybridSoftMuID/O");
    tout->Branch("Bmu2HybridSoftMuID", &Bmu2HybridSoftMuID_out, "Bmu2HybridSoftMuID/O");
    tout->Branch("Bmu1isAcc", &Bmu1isAcc_out, "Bmu1isAcc/O");
    tout->Branch("Bmu2isAcc", &Bmu2isAcc_out, "Bmu2isAcc/O");
    tout->Branch("Bujpt", &Bujpt_out, "Bujpt/F");
    tout->Branch("Bujeta", &Bujeta_out, "Bujeta/F");
    tout->Branch("Bujphi", &Bujphi_out, "Bujphi/F");
    tout->Branch("Bujlxy", &Bujlxy_out, "Bujlxy/F");
    tout->Branch("BdiTrackFitValid", &BdiTrackFitValid_out, "BdiTrackFitValid/O");
    tout->Branch("Btrk1Dz1", &Btrk1Dz1_out, "Btrk1Dz1/F");
    tout->Branch("Btrk2Dz1", &Btrk2Dz1_out, "Btrk2Dz1/F");
    tout->Branch("Btrk1DzError1", &Btrk1DzError1_out, "Btrk1DzError1/F");
    tout->Branch("Btrk2DzError1", &Btrk2DzError1_out, "Btrk2DzError1/F");
    tout->Branch("Btrk1Dxy1", &Btrk1Dxy1_out, "Btrk1Dxy1/F");
    tout->Branch("Btrk2Dxy1", &Btrk2Dxy1_out, "Btrk2Dxy1/F");
    tout->Branch("Btrk1DxyError1", &Btrk1DxyError1_out, "Btrk1DxyError1/F");
    tout->Branch("Btrk2DxyError1", &Btrk2DxyError1_out, "Btrk2DxyError1/F");
    tout->Branch("Btktketa", &Btktketa_out, "Btktketa/F");
    tout->Branch("Btktkphi", &Btktkphi_out, "Btktkphi/F");
    tout->Branch("Btktky", &Btktky_out, "Btktky/F");
    tout->Branch("Bdoubletpt", &Bdoubletpt_out, "Bdoubletpt/F");
    tout->Branch("Bdoubleteta", &Bdoubleteta_out, "Bdoubleteta/F");
    tout->Branch("Bdoubletphi", &Bdoubletphi_out, "Bdoubletphi/F");
    tout->Branch("Bdoublety", &Bdoublety_out, "Bdoublety/F");
    tout->Branch("Bnorm_trk1Dz", &Bnorm_trk1Dz_out, "Bnorm_trk1Dz/F");
    tout->Branch("Bnorm_trk2Dz", &Bnorm_trk2Dz_out, "Bnorm_trk2Dz/F");
    if (twoTrack) {
        tout->Branch("BtrkLeadPt", &BtrkLeadPt_out, "BtrkLeadPt/F");
        tout->Branch("BtrkSubPt", &BtrkSubPt_out, "BtrkSubPt/F");
        tout->Branch("BtrkLeadPtFrac", &BtrkLeadPtFrac_out, "BtrkLeadPtFrac/F");
        tout->Branch("BtrkSubPtFrac", &BtrkSubPtFrac_out, "BtrkSubPtFrac/F");
        tout->Branch("BtrkMaxdR", &BtrkMaxdR_out, "BtrkMaxdR/F");
        tout->Branch("BtrkMindR", &BtrkMindR_out, "BtrkMindR/F");
        tout->Branch("BtrkMaxAbsEta", &BtrkMaxAbsEta_out, "BtrkMaxAbsEta/F");
    }
    tout->Branch("BmuLeadPt", &BmuLeadPt_out, "BmuLeadPt/F");
    tout->Branch("BmuSubPt", &BmuSubPt_out, "BmuSubPt/F");
    tout->Branch("BmuMaxdR", &BmuMaxdR_out, "BmuMaxdR/F");
    tout->Branch("BmuMindR", &BmuMindR_out, "BmuMindR/F");
    tout->Branch("BmuMaxAbsEta", &BmuMaxAbsEta_out, "BmuMaxAbsEta/F");
    if (isMC){
        tout->Branch("Bgen", &Bgen_out, "Bgen/F");
        tout->Branch("pThatreweight", &pThatreweight_out, "pThatreweight/D");
        tout->Branch("pthat", &pthat, "pthat/F");
        tout->Branch("pThats", &pThats_out, "pThats/I");
    }

    tevent->Branch("RunNo", &RunNo_out, "RunNo/I");
    tevent->Branch("EvtNo", &EvtNo_out, "EvtNo/l");
    tevent->Branch("LumiNo", &LumiNo_out, "LumiNo/I");
    tevent->Branch("eventIndex", &eventIndex_out, "eventIndex/L");
    tevent->Branch("Bsize", &Bsize_out, "Bsize/I");
    tevent->Branch("PVx", &PVx_out, "PVx/F");
    tevent->Branch("PVy", &PVy_out, "PVy/F");
    tevent->Branch("PVz", &PVz_out, "PVz/F");
    tevent->Branch("PVxE", &PVxE_out, "PVxE/F");
    tevent->Branch("PVyE", &PVyE_out, "PVyE/F");
    tevent->Branch("PVzE", &PVzE_out, "PVzE/F");
    tevent->Branch("PVnchi2", &PVnchi2_out, "PVnchi2/F");
    tevent->Branch("nChargedTracks", &nChargedTracks_out, "nChargedTracks/F");
    tevent->Branch("nChargedTracks_LOOSE", &nChargedTracks_LOOSE_out, "nChargedTracks_LOOSE/F");
    tevent->Branch("nChargedTracks_TIGHT", &nChargedTracks_TIGHT_out, "nChargedTracks_TIGHT/F");
    tevent->Branch("CentBin", &CentBin_out, "CentBin/I");
    tevent->Branch("nCandPassPreselection", &nCandPassPreselection_out, "nCandPassPreselection/I");
    tevent->Branch("nCandPassFiducial", &nCandPassFiducial_out, "nCandPassFiducial/I");
    tevent->Branch("nCandWritten", &nCandWritten_out, "nCandWritten/I");
    if (isMC) {
        tevent->Branch("pThatreweight", &pThatreweight_out, "pThatreweight/D");
        tevent->Branch("pthat", &pthat, "pthat/F");
        tevent->Branch("pThats", &pThats_out, "pThats/I");
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
        if (isMC) {
            pThats_out = campaignPthat[tin->GetTreeNumber()];
            std::size_t k = 0;
            while (k + 1 < campaigns.size() && pthat > campaigns[k + 1].pthat) ++k;
            pThatreweight_out = 1. / cumulativeLuminosity[k];
        }

        // Event-level variables
        RunNo_out = RunNo;
        EvtNo_out = EvtNo;
        LumiNo_out = LumiNo;
        eventIndex_out = ev;
        Bsize_out = Bsize;
        PVx_out = PVx;
        PVy_out = PVy;
        PVz_out = PVz;
        PVxE_out = PVxE;
        PVyE_out = PVyE;
        PVzE_out = PVzE;
        PVnchi2_out = PVnchi2;
        CentBin_out = CentBin;
        nChargedTracks_out = nChargedTracks;
        nChargedTracks_LOOSE_out = nChargedTracks_LOOSE;
        nChargedTracks_TIGHT_out = nChargedTracks_TIGHT;
        nCandPassPreselection_out = 0;
        nCandPassFiducial_out = 0;
        nCandWritten_out = 0;

        // --------------------------------------------------
        // Candidate loop
        // --------------------------------------------------
        for(int i=0; i<Bsize; ++i)
        {
            // Pre-Selection (quality) Cuts
            // DATA is already filtered (except for the cuts marked with <-------- ** ); 
            // MC is not.
            if (true){
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
                if (TREENAME !=  "ntKp"){ 
                    if(abs(Btrk2Eta[i]) > 2.4) continue;
                    if((isPP && !(Btrk2Pt[i] > 0.5)) || (isPbPb && !(Btrk2Pt[i] > 0.9)) ) continue;
                    if(((Btrk2PtErr[i] / Btrk2Pt[i]) >= 0.1)) continue;
                    if(((Btrk2nPixelLayer[i] + Btrk2nStripLayer[i]) <= 10)) continue;
                    if(((Btrk2Chi2ndf[i]/(Btrk2nPixelLayer[i] + Btrk2nStripLayer[i])) >= 0.18)) continue;
                    if(!Btrk2highPurity[i]) continue;
                    //diTrack system
                    if (TREENAME == "ntphi"   && !(abs(Btktkmass[i] - PHI_MASS)  <0.015)) continue;
                    if (TREENAME == "ntKstar" && !(abs(Btktkmass[i] - KSTAR_MASS)<0.150)) continue;
                }

                // Candidate
                if (TREENAME == "ntmix" && (Bmass[i] > 4.0 || Bmass[i] < 3.6)) continue;
                if(Bchi2Prob[i] < 0.005) continue;
                if(TREENAME == "ntmix" && Bpt[i]<4) continue; 
                else if (Bpt[i]<1) continue;
            }
            ++nCandPassPreselection_out;

            // Fiducial Region
            if (Fid_region){
                if (!((Bpt[i] >= 5) && abs(By[i]) <= 2.4)) continue;
            }
            ++nCandPassFiducial_out;

            // keep MC signal only
            Bgen_out = Bgen[i];
            if( isMC && !(Bgen_out==23333 || Bgen_out==24333 || Bgen_out==23433 || Bgen_out==24433 || Bgen_out==41000 )) continue;

            candidateIndex_out = i;
            Bindex_out = Bindex[i];
            Bmass_out   = Bmass[i];
            Bpt_out     = Bpt[i];
            By_out      = By[i];
            Bchi2Prob_out = Bchi2Prob[i];
            Btrk1dR_out = Btrk1dR[i];
            Btrk2dR_out = Btrk2dR[i];
            BtrkPtimb_out = BtrkPtimb[i];
            Btktkpt_out   = Btktkpt[i];
            Bujmass_out   = Bujmass[i];
            BujvProb_out  = BujvProb[i];
            Bnorm_svpvDistance_2D_out = Bnorm_svpvDistance_2D[i];
            BsvpvDistance_2D_out = BsvpvDistance_2D[i];
            BsvpvDisErr_2D_out = BsvpvDisErr_2D[i];
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
            Btrk1Eta_out = Btrk1Eta[i];
            Btrk2Eta_out = Btrk2Eta[i];
            Btrk1Phi_out = Btrk1Phi[i];
            Btrk2Phi_out = Btrk2Phi[i];
            Btrk1PtErr_out = Btrk1PtErr[i];
            Btrk2PtErr_out = Btrk2PtErr[i];
            Btrk1Chi2ndf_out = Btrk1Chi2ndf[i];
            Btrk2Chi2ndf_out = Btrk2Chi2ndf[i];
            Btrk1nPixelLayer_out = Btrk1nPixelLayer[i];
            Btrk2nPixelLayer_out = Btrk2nPixelLayer[i];
            Btrk1nStripLayer_out = Btrk1nStripLayer[i];
            Btrk2nStripLayer_out = Btrk2nStripLayer[i];
            Btrk1highPurity_out = Btrk1highPurity[i];
            Btrk2highPurity_out = Btrk2highPurity[i];
            BvtxX_out = BvtxX[i];
            BvtxY_out = BvtxY[i];
            BLxy_out  = BLxy[i];
            Bmu1pt_out = Bmu1pt[i];
            Bmu2pt_out = Bmu2pt[i];
            Bmu1eta_out = Bmu1eta[i];
            Bmu2eta_out = Bmu2eta[i];
            Bmu1phi_out = Bmu1phi[i];
            Bmu2phi_out = Bmu2phi[i];
            Bmu1isTriggered_out = Bmu1isTriggered[i];
            Bmu2isTriggered_out = Bmu2isTriggered[i];
            Bmu1SoftMuID_out = Bmu1SoftMuID[i];
            Bmu2SoftMuID_out = Bmu2SoftMuID[i];
            Bmu1HybridSoftMuID_out = Bmu1HybridSoftMuID[i];
            Bmu2HybridSoftMuID_out = Bmu2HybridSoftMuID[i];
            Bmu1isAcc_out = Bmu1isAcc[i];
            Bmu2isAcc_out = Bmu2isAcc[i];
            Bujpt_out = Bujpt[i];
            Bujeta_out = Bujeta[i];
            Bujphi_out = Bujphi[i];
            Bujlxy_out = Bujlxy[i];
            BdiTrackFitValid_out = BdiTrackFitValid[i];
            Btrk1Dz1_out = Btrk1Dz1[i];
            Btrk2Dz1_out = Btrk2Dz1[i];
            Btrk1DzError1_out = Btrk1DzError1[i];
            Btrk2DzError1_out = Btrk2DzError1[i];
            Btrk1Dxy1_out = Btrk1Dxy1[i];
            Btrk2Dxy1_out = Btrk2Dxy1[i];
            Btrk1DxyError1_out = Btrk1DxyError1[i];
            Btrk2DxyError1_out = Btrk2DxyError1[i];
            Btktketa_out = Btktketa[i];
            Btktkphi_out = Btktkphi[i];
            Btktky_out = Btktky[i];
            Bdoubletpt_out = Bdoubletpt[i];
            Bdoubleteta_out = Bdoubleteta[i];
            Bdoubletphi_out = Bdoubletphi[i];
            Bdoublety_out = Bdoublety[i];
            Bnorm_trk1Dz_out = Bnorm_trk1Dz[i];
            Bnorm_trk2Dz_out = Bnorm_trk2Dz[i];

            if (twoTrack) {
                BtrkLeadPt_out = std::max(Btrk1Pt_out, Btrk2Pt_out);
                BtrkSubPt_out = std::min(Btrk1Pt_out, Btrk2Pt_out);
                BtrkLeadPtFrac_out = BtrkLeadPt_out / Bpt_out;
                BtrkSubPtFrac_out = BtrkSubPt_out / Bpt_out;
                BtrkMaxdR_out = std::max(Btrk1dR_out, Btrk2dR_out);
                BtrkMindR_out = std::min(Btrk1dR_out, Btrk2dR_out);
                BtrkMaxAbsEta_out = std::max(std::abs(Btrk1Eta_out), std::abs(Btrk2Eta_out));
            }
            const float mu1dR = std::hypot(Bmu1eta_out - Bujeta_out, TVector2::Phi_mpi_pi(Bmu1phi_out - Bujphi_out));
            const float mu2dR = std::hypot(Bmu2eta_out - Bujeta_out, TVector2::Phi_mpi_pi(Bmu2phi_out - Bujphi_out));
            BmuLeadPt_out = std::max(Bmu1pt_out, Bmu2pt_out);
            BmuSubPt_out = std::min(Bmu1pt_out, Bmu2pt_out);
            BmuMaxdR_out = std::max(mu1dR, mu2dR);
            BmuMindR_out = std::min(mu1dR, mu2dR);
            BmuMaxAbsEta_out = std::max(std::abs(Bmu1eta_out), std::abs(Bmu2eta_out));

            if(!std::isfinite(Bmass_out) || !std::isfinite(Bpt_out) || !std::isfinite(By_out) || !std::isfinite(Bnorm_trk1Dxy_out) ||
            !std::isfinite(CentBin_out) || !std::isfinite(Bchi2Prob_out) || !std::isfinite(Btrk1dR_out) || !std::isfinite(Bnorm_svpvDistance_2D_out)) continue;
            ++nCandWritten_out;
            tout->Fill();
        }
        tevent->Fill();
    }

    // --------------------------------------------------
    // Copy ntGen tree (flat, one row per candidate in mass window)
    // --------------------------------------------------
    if (isMC) {
        fout->cd();

        // Gen-level input branches
        Int_t Gsize;
        static Int_t GpdgId[MAXCAND];
        static Float_t Gmu1eta[MAXCAND];
        static Float_t Gmu1pt[MAXCAND];
        static Float_t Gmu2eta[MAXCAND];
        static Float_t Gmu2pt[MAXCAND];
        static Float_t Gtk1pt[MAXCAND];
        static Float_t Gtk1eta[MAXCAND];
        static Float_t Gtk2pt[MAXCAND];
        static Float_t Gtk2eta[MAXCAND];
        static Float_t Gpt[MAXCAND];
        static Float_t Gy[MAXCAND];
        static Int_t GisSignal[MAXCAND];

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
        Float_t pthat_gen = 0.;
        tgen->SetBranchAddress("pthat", &pthat_gen);

        // Output tree (flat)
        tgenOut = new TTree("ntGen", "Gen-level (flat, filtered)");
        tgenOut->SetAutoSave(0);
        Float_t Gmu1eta_out, Gmu1pt_out, Gmu2eta_out, Gmu2pt_out;
        Float_t Gtk1pt_out, Gtk1eta_out, Gtk2pt_out, Gtk2eta_out, Gpt_out, Gy_out;
        Int_t GpdgId_out;
        Double_t pThatreweight_gen_out = 1.;
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
        tgenOut->Branch("pThatreweight", &pThatreweight_gen_out, "pThatreweight/D");
        tgenOut->Branch("pthat", &pthat_gen, "pthat/F");
        Int_t pThats_gen_out = 0;
        tgenOut->Branch("pThats", &pThats_gen_out, "pThats/I");

        const Long64_t ngen = tgen->GetEntries();
        for (Long64_t ev=0; ev<ngen; ++ev) {
            tgen->GetEntry(ev);
            pThats_gen_out = campaignPthat[tgen->GetTreeNumber()];
            std::size_t k = 0;
            while (k + 1 < campaigns.size() && pthat_gen > campaigns[k + 1].pthat) ++k;
            pThatreweight_gen_out = 1. / cumulativeLuminosity[k];
            for (int i=0; i<Gsize; ++i) {

                if (TREENAME == "ntmix" && GisSignal[i] != 7) continue;               // keep only signal ntmix candidates for the selected MC sample
                else if (TREENAME == "ntphi"   && abs(GpdgId[i]) != 531) continue;  
                else if (TREENAME == "ntKstar" && abs(GpdgId[i]) != 511) continue;  
                else if (TREENAME == "ntKp"    && abs(GpdgId[i]) != 521) continue;

                // Fiducial Region
                if (Fid_region){
                    if ( !( (Gpt[i] >= 5) && abs(Gy[i]) <= 2.4) ) continue; 
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
    // Write top-level trees
    tout->Write("", TObject::kOverwrite);
    tevent->Write("", TObject::kOverwrite);
    if (isMC) tgenOut->Write("", TObject::kOverwrite);
    fout->Purge();
    fout->Close();
    delete tin;
    delete tgen;
    delete fout;
    std::cout << "Done & Saved -> " << outputFile << "\n";
}
