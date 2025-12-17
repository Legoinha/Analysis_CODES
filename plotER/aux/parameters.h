using namespace RooFit;
using namespace std;



int nbinsmasshisto = 100;

double minhisto = 0;
double maxhisto = 999;
double minhisto_B=5.;
double maxhisto_B=6.;
double minhisto_X=3.6;
double maxhisto_X=4.0;


const int N_pt_Bins_X = 6;
std::vector<double> ptbinsvec_X = {5, 7.5, 10, 15, 25, 50, 100};

///FOR TESTING TESTING
//const int N_pt_Bins_X = 1;
//std::vector<double> ptbinsvec_X = {15, 16};
///FOR TESTING TESTING

const int N_pt_Bins_B = 8;
std::vector<double> ptbinsvec_B = {2, 4, 6, 8, 10, 15, 20, 30, 60};

const int N_y_Bins_X = 4;
std::vector<double> ybinsvec = {0.0, 0.8, 1.5, 2.0, 2.4};

const int N_mult_Bins_X = 5;
std::vector<double> nmbinsvec = {0, 20, 40, 60, 80, 100};

const int N_cent_Bins_X = 5;
std::vector<double> centbinsvec = {0, 20, 40, 60, 80, 100};

