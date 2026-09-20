#pragma once

#include "TH1D.h"
#include "TString.h"

struct EffCase {
    TString suffix;
    TString label;
};

struct EffResult {
    EffCase method;
    TH1D* hAvg;
    TH1D* hYield;
};
