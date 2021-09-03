#ifndef BEXCLFIT_FCTS_H
#define BEXCLFIT_FCTS_H

#include <TROOT.h>


Double_t signalmassdis(Double_t x, Double_t sig, Double_t *parms);
Double_t bgrmassdis(Double_t x, Double_t *parms, double mmin, double mmax);
Double_t backgr(Double_t *x,Double_t *par);
Double_t bgrctdis(Double_t x, Double_t sig, Double_t *parms);
Double_t signalctdis(Double_t x, Double_t sig, Double_t *parms);
Double_t sigdis(Double_t x, Double_t sig, Double_t *parms);

void testfun(void);
#endif   // BEXCLFIT_FCTS
