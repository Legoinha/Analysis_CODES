#!/usr/bin/env python3
import os
import sys
import ROOT
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt




def get_particle_label(tree_name):
    PARTICLE_LABELS = {
    "ntKp": "B+",
    "ntKstar": "B0",
    "ntphi": "Bs",}
    return PARTICLE_LABELS.get(tree_name, tree_name)

def plot_correlation_dataMC(TREE="ntmix", systemNAME="ppRef"):
    gROOT = ROOT.gROOT
    gROOT.SetBatch(True)

    # VARIABLES to correlate
    variables = [ "Btrk1Pt", "Btrk2Pt", "Btrk1dR", "Btrk2dR", "BtrkPtimb", "BtktkvProb", "Bchi2Prob", "BLxy", "Bmass", "BQvalue", "Bpt"]
    # VARIABLES to correlate

    # PATHS to data/MC files
    path_to_data = ""
    path_to_mc = ""

    if systemNAME == "PbPb23":
        if TREE == "ntmix":
            path_to_data = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_DATA.root"
            path_to_mc = "/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23/flat_ntmix_PbPb23_MC.root"
    else:
        if TREE == "ntmix": #X3872 ppRef
            path_to_data = "/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_DATA_wScore.root"
            path_to_mc = "/eos/user/h/hmarques/Analysis_CODES/flatER/X3872/flat_ntmix_ppRef_MC_wScore.root"
        elif TREE == "ntphi": #B0s ppRef
            path_to_data = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntphi_ppRef_DATA.root"
            path_to_mc = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntphi_ppRef_MC.root"
        elif TREE == "ntKp": #B+ ppRef
            path_to_data = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKp_ppRef_DATA.root"
            path_to_mc = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKp_ppRef_MC.root"
        elif TREE == "ntKstar": #B0 ppRef
            path_to_data = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_DATA.root"
            path_to_mc = "/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_ntKstar_ppRef_MC.root"

    # HARDCODE your selection here
    selection = "1"
    #if TREE == "ntmix":
    #    selection = "xgb_score > 0.5 && BQvalue < 0.15"
    #else:
    #    selection = "Bnorm_svpvDistance_2D > 4"

    # DATA sideband definition
    sideband = "1"
    if TREE == "ntmix":
        sideband = "(Bmass > 3.95 || (Bmass > 3.75 && Bmass < 3.8))"
    else: 
        sideband = "(Bmass > 5.5)"

    data_cut = f"({selection}) && ({sideband})"
    data_file = ROOT.TFile.Open(path_to_data)
    mc_file = ROOT.TFile.Open(path_to_mc)
    data_tree = data_file.Get(TREE)
    mc_tree = mc_file.Get(TREE)

    common_vars = []
    for var in variables:
        has_data = data_tree.GetListOfBranches().FindObject(var)
        has_mc = mc_tree.GetListOfBranches().FindObject(var)
        if has_data and has_mc:
            common_vars.append(var)

    print("TREE:", TREE)
    print("systemNAME:", systemNAME)
    print("Selection:", selection)
    print("Sideband:", sideband)

    rdf_data = ROOT.RDataFrame(TREE, path_to_data).Filter(data_cut)
    arr_data = rdf_data.AsNumpy(common_vars)
    n_data = len(arr_data[common_vars[0]])
    mat_data = np.column_stack([np.asarray(arr_data[var], dtype=float) for var in common_vars])
    corr_data = np.corrcoef(mat_data, rowvar=False)
    out_dir = f"./correlations/{TREE}"
    os.makedirs(out_dir, exist_ok=True)

    draw_matrix(
        corr_data,
        common_vars,
        f"data sideband {systemNAME}",
        f"{out_dir}/data_sideband_corr_{systemNAME}_{TREE}.pdf"
    )

    mc_cases = [("1", get_particle_label(TREE), TREE)]
    if TREE == "ntmix":
        mc_cases = [
            ("isX3872==1", "X(3872)", TREE),
            ("isX3872==0", "psi(2S)", f"{TREE}_psi2s"),
        ]

    for mc_sel, particle_label, output_tag in mc_cases:
        mc_cut = f"({selection}) && ({mc_sel})"
        print("MC selection:", mc_sel)
        rdf_mc = ROOT.RDataFrame(TREE, path_to_mc).Filter(mc_cut)
        arr_mc = rdf_mc.AsNumpy(common_vars)
        n_mc = len(arr_mc[common_vars[0]])
        print("n_data:", n_data)
        print("n_mc:", n_mc)
        mat_mc = np.column_stack([np.asarray(arr_mc[var], dtype=float) for var in common_vars])
        corr_mc = np.corrcoef(mat_mc, rowvar=False)
        draw_matrix(
            corr_mc,
            common_vars,
            f"{particle_label} signal MC {systemNAME}",
            f"{out_dir}/mc_corr_{systemNAME}_{output_tag}.pdf"
        )

def draw_matrix(matrix, labels, title, out_name, vmax=1.0):
    fig_size = max(8, int(0.7 * len(labels)))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)
    ax.set_title(title)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=8)

    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(out_name)
    plt.close(fig)
    print(f"Saved correlation matrix to {out_name}")

if __name__ == "__main__":
    tree = "ntmix"
    system = "ppRef"
    if len(sys.argv) > 1:
        tree = sys.argv[1]
    if len(sys.argv) > 2:
        system = sys.argv[2]
    plot_correlation_dataMC(tree, system)

## python3 correlation_dataMC.py ntphi ppRef
