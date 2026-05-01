#ifndef EFFER_AUX_UTI_H
#define EFFER_AUX_UTI_H

#include "TString.h"

inline TString GetMCEffPath(TString treename, TString system)
{
    if (system.Contains("PbPb23")) {
        if (treename == "ntmix")       return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb23_scored_MC_X3872.root";
        if (treename == "ntmix_psi2s") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb23_scored_MC_PSI2S.root";
    } else if (system.Contains("PbPb24")) {
        if (treename == "ntmix")       return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb24_scored_MC_X3872.root";
        if (treename == "ntmix_psi2s") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb24_scored_MC_PSI2S.root";
    } else if (system.Contains("PbPb")) {
        if (treename == "ntmix")       return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb_scored_MC_X3872.root";
        if (treename == "ntmix_psi2s") return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb_scored_MC_PSI2S.root";
    } else {
        if (treename == "ntmix")        return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_ppRef_scored_MC_X3872.root";
        if (treename == "ntmix_psi2s")  return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_ppRef_scored_MC_PSI2S.root";
        if (treename == "ntphi")        return Form("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_%s_%s_MC.root", treename.Data(), system.Data());
        if (treename == "ntKp")         return Form("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_%s_%s_MC.root", treename.Data(), system.Data());
        if (treename == "ntKstar")      return Form("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_%s_%s_MC.root", treename.Data(), system.Data());
    }

    return "";
}

inline TString GetGenEffPath(TString treename, TString system)
{
    TString sharedSystem = system;
    if (system == "ppRef") sharedSystem = "ppRef24";

    if (treename == "ntmix") {
        if (system == "ppRef") return Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/%s/flat_ntmix_ppRef_MC_X3872.root", sharedSystem.Data());
        return Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/%s/flat_ntmix_%s_MC.root", sharedSystem.Data(), system.Data());
    }

    if (treename == "ntmix_psi2s") {
        if (system == "ppRef") return Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/%s/flat_ntmix_ppRef_MC_PSI2S.root", sharedSystem.Data());
        return Form("/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/%s/flat_ntmix_%s_MC.root", sharedSystem.Data(), system.Data());
    }

    return GetMCEffPath(treename, system);
}






inline TString GetDataEffPath(TString treename, TString system)
{
    if (treename == "ntmix" || treename == "ntmix_psi2s") {
        if (system.Contains("PbPb23")) return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb23_scored_DATA.root";
        if (system.Contains("PbPb24")) return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb24_scored_DATA.root";
        if (system.Contains("PbPb"))   return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_PbPb_scored_DATA.root";
        return "/eos/user/h/hmarques/Analysis_CODES/selectionER/scored_samples/flat_ntmix_ppRef_scored_DATA.root";
    }

    if (treename == "ntphi" || treename == "ntKp" || treename == "ntKstar") {
        return Form("/eos/user/h/hmarques/Analysis_CODES/flatER/Bmeson/flat_%s_%s_DATA.root", treename.Data(), system.Data());
    }

    return "";
}

#endif
