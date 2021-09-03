
#include <TROOT.h>

const int dim    = 45;  // number of fit parameters (here const later depends on number 
                        // of data samples)

const int CommonPar = 5; // Number of Parameters common to all samples
const int SSpecPar  = 8; // Number of Sample specific parameters

// 'pointers' to the parameters
const int mpeak    = 0;  // position of mass peak  
const int mscale   = 1;  // scale factor for calculated sigma of mass
const int ctscale  = 2;  // scale factor for calculated sigma of ctau
const int ctau     = 3;  // lifetime of B-meson
const int tailfrac = 4; 
const int mconst   = 5;  // const of linear background in the mass distribution 
                         // (slope follows from normalisation condition)
const int bgrfrac  = 6;  // background fraction  
const int lbdapos  = 7;  // lambda+ of right side exponential of bgr ctau distribution
const int fracpos  = 8;  // f+ of right side exponential of bgr tau distribution
const int lbdaneg  = 9;  // lambda- of left side exponential 
const int fracneg  = 10;  // f- of left side exponential
const int lbdapos2 = 11; // lambda++ of 2nd right side exponential of bgr ctau distribution
const int fracpos2 = 12; // f++ of 2nd right side exponential og bgr ctau distribution

const Double_t sqrt2    = sqrt(2.0);
const Double_t sqrt2pi  = 2.5066282746;
// 
// The following are some predefined cuts and constants
//
const Double_t PDGBmass    = 5.279;           // PDG value of B-meson mass in GeV
//const Double_t PDGBmass  = 5.624; // lamda_b
const Double_t masscut     = 8.0;             // cut for normalized mass distribution (# of sigmas) 
const Double_t ctaucut     = 0.01;            // ctau cut for mass histograms 
const Double_t masswin     = 0.150;           // width of mass window in GeV
const Double_t sbmasswin   = 0.100;           // width of side band mass windows in GeV


