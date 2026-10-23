# Step 1: read the flat trees once, apply the study selections, cache the needed columns as .npz.
import ROOT, numpy as np, sys
ROOT.gROOT.SetBatch(True); ROOT.EnableImplicitMT(8)
X = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/"
BM = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/Bmesons/"
PROMPT = "BLxy*(Bmass/Bpt) < 0.005"
SAMPLES = {
  "ppRef_psi2S": dict(data=(X+"ppRef24/flat_ntmix_ppRef_DATA.root","ntmix"), mc=(X+"ppRef24/flat_ntmix_ppRef_MC_PSI2S.root","ntmix_PSI2S"),
       sel="BQvalue < 0.15 && Btrk1dR < 0.5 && Btrk2dR < 0.5 && Bpt >= 7.5 && Bpt < 50 && "+PROMPT, win=(3.62,3.76)),
  "ppRef_X3872": dict(data=(X+"ppRef24/flat_ntmix_ppRef_DATA.root","ntmix"), mc=(X+"ppRef24/flat_ntmix_ppRef_MC_X3872.root","ntmix_X3872"),
       sel="BQvalue < 0.15 && Btrk1dR < 0.5 && Btrk2dR < 0.5 && Bpt >= 7.5 && Bpt < 50 && "+PROMPT, win=(3.80,3.96)),
  "PbPb23_psi2S": dict(data=(X+"PbPb23/flat_ntmix_PbPb23_DATA.root","ntmix"), mc=(X+"PbPb23/flat_ntmix_PbPb23_MC_PSI2S.root","ntmix_PSI2S"),
       sel="Bpt > 15 && Bpt < 50 && BQvalue < 0.15 && abs(By) < 1.6 && Btrk1dR < 0.35 && Btrk2dR < 0.35 && "+PROMPT, win=(3.62,3.76)),
  "PbPb23_Bp": dict(data=(BM+"PbPb23/flat_ntKp_PbPb23_DATA.root","ntKp"), mc=(BM+"PbPb23/flat_ntKp_PbPb23_MC.root","ntKp"),
       sel="Bpt > 10 && Bpt < 50 && abs(By) < 2.4 && Balpha < 0.05 && BsvpvDistance_2D > 0.02 && Btrk1Pt > 1.0", win=(5.05,5.60)),
  # 2024 data compared with the 2023 MC (no 2024 B+ MC available): relative year comparison only
  "PbPb24_Bp": dict(data=(BM+"PbPb24/flat_ntKp_PbPb24_DATA.root","ntKp"), mc=(BM+"PbPb23/flat_ntKp_PbPb23_MC.root","ntKp"),
       sel="Bpt > 10 && Bpt < 50 && abs(By) < 2.4 && Balpha < 0.05 && BsvpvDistance_2D > 0.02 && Btrk1Pt > 1.0", win=(5.05,5.60)),
  "ppRef_Bp": dict(data=(BM+"ppRef/flat_ntKp_ppRef_DATA.root","ntKp"), mc=(BM+"ppRef/flat_ntKp_ppRef_MC.root","ntKp"),
       sel="Bpt > 10 && Bpt < 50 && abs(By) < 2.4 && Balpha < 0.05 && BsvpvDistance_2D > 0.02 && Btrk1Pt > 1.0", win=(5.05,5.60)),
}
COMMON = ["Bmass","Bchi2Prob","BujvProb","PVnchi2","EvtNo","Bpt","nChargedTracks","Bujmass","Btrk1Pt","Btrk1Eta",
          "Btrk1Dxy1","Btrk1DxyError1","Btrk1Dz1","Btrk1DzError1"]
TWO_TRK = ["BtktkvProb","Btrk2Pt","Btrk2Eta","Btrk2Dxy1","Btrk2DxyError1","Btrk2Dz1","Btrk2DzError1"]
name = sys.argv[1]; s = SAMPLES[name]
cols = COMMON + ([] if name.endswith("_Bp") else TWO_TRK)
win = f"Bmass > {s['win'][0]} && Bmass < {s['win'][1]}"
def frame(f, t):
    df = ROOT.RDataFrame(t, f)
    # EvtNo is only used for the even/odd closure split; fall back to the entry number when absent
    return df if "EvtNo" in [str(c) for c in df.GetColumnNames()] else df.Define("EvtNo", "(double)rdfentry_")
d = frame(*s["data"]).Filter(s["sel"]).Filter(win)
m = frame(*s["mc"]).Filter(s["sel"]).Filter(win).Filter("Bgen==23333")
D = d.AsNumpy(cols); M = m.AsNumpy(cols + ["pThatreweight"])
np.savez_compressed(f"cache/{name}.npz", **{"d_"+k: np.asarray(v, dtype=np.float64) for k,v in D.items()},
                    **{"m_"+k: np.asarray(v, dtype=np.float64) for k,v in M.items()})
print(name, "data", len(D["Bmass"]), "MC", len(M["Bmass"]))
