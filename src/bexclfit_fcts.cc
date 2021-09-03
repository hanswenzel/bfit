#include <iostream>
#include "bexclfit_fcts.h"
#include "TMath.h"
#include "global.h"

extern int offset,combine_flg;
//extern Double_t sbmasswin,masswin;

void testfun()
{
  printf("offset = %d\n",offset);
  printf("combine_flg = %d\n",combine_flg);
}

Double_t signalmassdis(Double_t x, Double_t sig, Double_t *parms)
{ 
  //---------------------------------------------------------------------------
  //  This is a simple function to parameterize the pure B-signal in the 
  //  mass distribution. This is parametrized by a simple normalized gaussian 
  //  distribution. 
  //---------------------------------------------------------------------------
  Double_t result,m=x,sigma;
  sigma = parms[mscale]*sig;
  result = (exp(-((m-parms[mpeak])*(m-parms[mpeak]))/(2.0*sigma*sigma)))/(sigma*sqrt2pi);
  return result;
} 


Double_t bgrmassdis(Double_t x, Double_t *parms, double mmin, double mmax)
{ 
  //---------------------------------------------------------------------------
  //  This is a simple function to parameterize the background in the 
  //  mass distribution. This is parametrized by a simple linear function. 
  //---------------------------------------------------------------------------
  Double_t result,m=x;
  Double_t a;
 
  a = 2.0/(mmax*mmax-mmin*mmin) - (2.0 *parms[mconst+offset])/(mmax+mmin);  
  result = a * m + parms[mconst+offset];
  return result;
} 


Double_t backgr(Double_t *x,Double_t *par)
{ 
  //---------------------------------------------------------------------------
  //  This is a simple function to parameterize the shape of the bgr in the 
  //  exclusive lifetime distribution. This function is used in the binned fit. 
  //
  //       par[0]: NORMALIZATION FACTOR
  //       par[1]: SIGMA OF GAUSSIAN
  //       par[2]: LAMBDA1 OF RIGHT SIDE EXPONENTIAL
  //       par[3]: F1 OF RIGHT SIDE EXPONENTIAL
  //       par[4]: LAMBDA2 OF LEFT SIDE EXPONENTIAL
  //       par[5]: F2 OF LEFT SIDE EXPONENTIAL
  //       par[6]: LAMBDA3 OF RIGHT SIDE EXPONENTIAL 2
  //       par[7]: F3 OF RIGHT SIDE EXPONENTIAL 2
  //---------------------------------------------------------------------------
  Double_t result,l=*x;

  result = exp((-0.5*l*l)/(par[1]*par[1]))*(par[0]/(par[1]*sqrt2pi));
  //if (l>=0) {
  // result=(1-par[3]-par[5])*result+par[0]*(par[3]/par[2])*exp(-l/par[2]);
  //} else {
  // result=(1-par[3]-par[5])*result+par[0]*(par[5]/par[4])*exp(l/par[4]);
  //}
    if (l>=0) {
    result=(1-par[3]-par[5]-par[7])*result+par[0]*(par[3]/par[2])*exp(-l/par[2])+par[0]*(par[7]/par[6])*exp(-l/par[6]);
 } else {
    result=(1-par[3]-par[5])*result+par[0]*(par[5]/par[4])*exp(l/par[4]);
    }


  return result;
} 




Double_t bgrctdis(Double_t x, Double_t sig, Double_t *parms)
{ 
  //----------------------------------------------------------------------------
  //  This is a simple function to parameterize the shape of the bgr in the 
  //  exclusive lifetime distribution. 
  //----------------------------------------------------------------------------
  Double_t result,l=x, sigma=sig;
  
  sigma = sigma*parms[ctscale];
  
  Double_t norm=1.0/(sigma*sqrt2pi);
  result =norm*exp((-0.5*l*l)/(sigma*sigma));
  if (l>=0) {
    result=(1-parms[fracpos+offset]-parms[fracneg+offset]-parms[fracpos2+offset])*result+(parms[fracpos+offset]/parms[lbdapos+offset])*exp(-l/parms[lbdapos+offset])+(parms[fracpos2+offset]/parms[lbdapos2+offset])*exp(-l/parms[lbdapos2+offset]);
   
  } else {

    result=(1-parms[fracpos+offset]-parms[fracneg+offset]-parms[fracpos2+offset])*result+(parms[fracneg+offset]/parms[lbdaneg+offset])*exp(l/parms[lbdaneg+offset]); 
  }
  
  return result;
} 


/*
Double_t signalctdis(Double_t x, Double_t sig, Double_t *parms)
{
  //---------------------------------------------------------------------------
  //  This is a simple function to parameterize the shape of the signal 
  //  in the exclusive lifetime  distribution.
  //---------------------------------------------------------------------------
  Double_t result,l=x;
  Double_t sigma = sig;

  sigma = sigma*parms[ctscale];

  result = ((0.5/parms[ctau])*exp((0.5*sigma*sigma)/(parms[ctau]*parms[ctau])-l/parms[ctau])*(1+TMath::Erf(l/(sqrt2*sigma)-sigma/(parms[ctau]*sqrt2))));

  return result;
}
*/

Double_t signalctdis(Double_t x, Double_t sig, Double_t *parms)
{
  //---------------------------------------------------------------------------
  //  This is a simple function to parameterize the shape of the signal 
  //  in the exclusive lifetime  distribution.
  //---------------------------------------------------------------------------
  Double_t result,l=x;
  Double_t sigma = sig;
  Double_t s2 = 0.100;
  sigma = sigma*parms[ctscale];


  //  result = ((0.5/parms[ctau])*exp((0.5*sigma*sigma)/(parms[ctau]*parms[ctau])-l/parms[ctau])*(1+TMath::Erf(l/(sqrt2*sigma)-sigma/(parms[ctau]*sqrt2))));

  result = (0.5/parms[ctau])*(
    (1-parms[tailfrac])*(exp((0.5*sigma*sigma)/(parms[ctau]*parms[ctau])-l/parms[ctau])*(1+TMath::Erf(l/(sqrt2*sigma)-sigma/(parms[ctau]*sqrt2))))
    +parms[tailfrac]*(exp((0.5*s2*s2)/(parms[ctau]*parms[ctau])-l/parms[ctau])*(1+TMath::Erf(l/(sqrt2*s2)-s2/(parms[ctau]*sqrt2))))
    );
  return result;
}




