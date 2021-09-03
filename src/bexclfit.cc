//--------------------------------------------------------------------------
// File and Version Information:  
// bexclfit.cc
//
// Description:
//
//      This Program performs binned and unbinned fits to the lifetime 
//      and mass distribution of exclusively reconstructed B-mesons in 
//      the  B-->J/Psi K mode. 
//      Finally, a combined fit of mass and lifetime distributions 
//      and the statistical combination of different decay modes of charged 
//      and neutral B-mesons is done.
// 
//      The meaning of the Parameters obtained by the unbinned fitting procedure is: 
//      described in bexclfit_fcts.h
//
// Environment:
//      Software developed for the CDF Collaboration
//
// Author List:
//         Hans Wenzel    Fermilab        email::wenzel@fnal.gov
//         Andreas Heiss  IEKP Karlsruhe  email::heiss@fnal.gov
//
//------------------------------------------------------------------------
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <TROOT.h>
#include <TFile.h>
#include <TRint.h>
#include <TH1.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TSystem.h>
#include <TPaveStats.h>
#include <TStyle.h>  
#include <TKey.h>
#include <TGraph.h>
#include "BDecay.h"
#include "bexclfit_fcts.h"
#include "global.h"
extern "C" {
#include <libconf.h>
}   


const char par_name[dim][20]={
                       "mpeak  ",
                       "mscale ",
                       "ctscale",
                       "ctau   ",
		       "tailfrac",
                       "mconst0","bgrfrac0","lbdapos0","fracpos0","lbdapneg0","fracneg0",
		       "lbda2pos0","frac2pos0",
                       "mconst1","bgrfrac1","lbdapos1","fracpos1","lbdapneg1","fracneg1",
		       "lbda2pos1","frac2pos1",
                       "mconst2","bgrfrac2","lbdapos2","fracpos2","lbdapneg2","fracneg2",
                       "lbda2pos2","frac2pos2",
		       "mconst3","bgrfrac3","lbdapos3","fracpos3","lbdapneg3","fracneg3",
                       "lbda2pos3","frac2pos3",
		       "mconst4","bgrfrac4","lbdapos4","fracpos4","lbdapneg4","fracneg4",
                       "lbda2pos4","frac2pos4"};


Double_t params[dim],errparams[dim];
static Double_t step[dim] = {0.001,0.001,0.001,0.0005,0.001,
                             0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
			     0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
			     0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
			     0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,
			     0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};

Double_t sspars[10][13],sserrpars[10][13],sfpar[13],errsfpar[13];   
Double_t m_bgrfrac[10];
Double_t sigma;
Double_t binwidth  = 0.0; 
Double_t mbinwidth = 0.0;


int count_bdecay=0,samplesize[10];
BDecay  *BtoJPSIK_array[10],*Bsample,*BtoJPSIK;
std::vector<Measurement> *junker;                  // temporary hans 
TF1 *backufcn, *sigufcn, *lifeufcn, *massufcn, *backfit;
TF1 *gaussfcn;
TCanvas *c[10];
TCanvas *cc,*cs;
int cnum;

char mfitstat[10][15],bfitstat[10][15],fitstat[10][15];
char combfitstat[15],mcombfitstat[15];

int ierflg = 0;
using namespace std;

// Declare global variables also known within bexclfit_fcts.cc 
int combine_flg; // 0 = single sample fit; 1 = combined fit
int offset;      // This is always 0 in the single sample fit loop, will be used
                 // in the combined fit



void printsummary(int l)
{
  int ii;

  for (ii=0; ii<count_bdecay; ii++) {
    cout << char(27) << "[34m";
    printf("\n---------- Sample no. %2d : %s ----------\n",ii,BtoJPSIK_array[ii]->GetName());
    cout << char(27) << "[0m";
    printf("Status of massfit, bgrfit, final fit : ");
    if (!strncmp(mfitstat[ii],"CONVERGED",9)) {cout << char(27) << "[32m";}
    else {cout << char(27) << "[31m";}
    printf("%s ",mfitstat[ii]);
    cout << char(27) << "[0m" << ","; 
    if (!strncmp(bfitstat[ii],"CONVERGED",9)) {cout << char(27) << "[32m";}
    else {cout << char(27) << "[31m";}
    printf("%s ",bfitstat[ii]);
    cout << char(27) << "[0m" << ","; 
    if (!strncmp(fitstat[ii],"CONVERGED",9)) {cout << char(27) << "[32m";}
    else {cout << char(27) << "[31m";}
    printf("%s \n",fitstat[ii]);
    cout << char(27) << "[0m";
    printf("mpeak    %10.5f  +/- %10.7f\n",sspars[ii][0],sserrpars[ii][0]);
    printf("mscale   %10.5f  +/- %10.7f\n",sspars[ii][1],sserrpars[ii][1]);
    printf("ctscale  %10.5f  +/- %10.7f\n",sspars[ii][2],sserrpars[ii][2]);
    printf("ctau     %10.5f  +/- %10.5f\n",10000*sspars[ii][3],10000*sserrpars[ii][3]);
    printf("tailfrac %10.5f  +-  %10.7f\n",sspars[ii][4],sserrpars[ii][4]);
    printf("mconst   %10.5f  +/- %10.7f\n",sspars[ii][5],sserrpars[ii][5]);
    printf("bgrfrac  %10.5f  +/- %10.7f\n",sspars[ii][6],sserrpars[ii][6]);
    printf("lbdapos  %10.5f  +/- %10.5f\n",10000*sspars[ii][7],10000*sserrpars[ii][7]);
    printf("fracpos  %10.5f  +/- %10.7f\n",sspars[ii][8],sserrpars[ii][8]);
    printf("lbdaneg  %10.5f  +/- %10.5f\n",10000*sspars[ii][9],10000*sserrpars[ii][9]);
    printf("fracneg  %10.5f  +/- %10.7f\n",sspars[ii][10],sserrpars[ii][10]);
    printf("lbda2pos %10.5f  +/- %10.5f\n",10000*sspars[ii][11],10000*sserrpars[ii][11]);
    printf("frac2pos %10.5f  +/- %10.7f\n",sspars[ii][12],sserrpars[ii][12]);
    cout << char(27) << "[31m";
   printf("Number of candidates : %f +/- %f\n\n",(1.-sspars[ii][bgrfrac])*(float)samplesize[ii],(float)samplesize[ii]*sserrpars[ii][bgrfrac]);
   cout << char(27) << "[0m";
  }
  if (l>0) {
    printf("Combined lifetime : %10.5f  +/- %10.5f\n",10000*params[ctau],10000*errparams[ctau]);
    printf("Combined mass     : %10.8f  +/- %10.8f\n",params[mpeak],errparams[mpeak]);
    printf("ctscale           : %10.8f  +/- %10.8f\n",params[ctscale],errparams[ctscale]);
    printf("mscale            : %10.8f  +/- %10.8f\n",params[mscale],errparams[mscale]);
    printf("tailfrac          : %10.8f  +/- %10.8f\n",params[tailfrac],errparams[tailfrac]);
    printf("Status of combined fit (mass, final) : ");
    if (!strncmp(mcombfitstat,"CONVERGED",9)) {cout << char(27) << "[32m";}
    else {cout << char(27) << "[31m";}
    printf("%s  ",mcombfitstat);
    cout << char(27) << "[0m";
    if (!strncmp(combfitstat,"CONVERGED",9)) {cout << char(27) << "[32m";}
    else {cout << char(27) << "[31m";}
    printf("%s ",combfitstat);
    cout << char(27) << "[0m";
    printf("\n"); 
  }
  printf("\n");  
}

void massfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *params, Int_t iflag)
{
  //---------------------------------------------------------------------------
  // this is the function used by minuit to do the unbinned fit to the mass 
  // distribution
  //----------------------------------------------------------------------------
  double mass;   
  int nn;          // this stuff is not implemented yet !

  offset=0; // Offset of sample specific parameters used in combined fit only
  f = 0.0;

  if (combine_flg) {
    nn = count_bdecay;
  } else {
    nn = 1; // Only 'loop' once
    Bsample = BtoJPSIK;
    offset=0;
  }
  for (int ii=0; ii<nn; ii++) {
    if (combine_flg==1) {
      Bsample = BtoJPSIK_array[ii];  // set pointer to current sample
      offset = ii*SSpecPar;          // calculate offset of sample specific parameters
    }
  std::vector<Measurement>::iterator iterator = Bsample->beginiter();   
    while (iterator != Bsample->enditer()) {
      if (iterator->Getctau() > ctaucut) {
	mass =  params[bgrfrac+offset]*
	  bgrmassdis(iterator->Getmass(),params,Bsample->Getmin_mass(),Bsample->Getmax_mass())+
	  (1.0-params[bgrfrac+offset])*signalmassdis(iterator->Getmass(),iterator->Getsigmamass(),params);
	f    = log(mass)+f;
      }
      ++iterator;
    }
  } // next BDecay object

  f= -2.0*f;
  return ;
}


Double_t draw_massfcn(Double_t *x,Double_t *params)
{
  //------------------------------------------------------------
  // this is the function used to draw the unbinned fit results 
  // for the meaning of the fit parameters see function massdis
  //------------------------------------------------------------
  Double_t mass;   
  Double_t sum = 0.0;
  Double_t m=*x;
  int nn;

  if (combine_flg) {
    nn = count_bdecay;
  } else {
    nn = 1; // Only 'loop' once
    Bsample = BtoJPSIK; 
    offset = 0; 
  }
  for (int ii=0;ii<nn;ii++) {
    if (combine_flg) {
      Bsample = BtoJPSIK_array[ii];
      offset=ii*SSpecPar;
    }
   std::vector<Measurement>::iterator iterator = Bsample->beginiter();   
    double mmin = Bsample->Getmin_mass();
    double mmax = Bsample->Getmax_mass();
    while (iterator != Bsample->enditer()) {
      if (iterator->Getctau() > ctaucut) {
	mass =  params[bgrfrac+offset]*bgrmassdis(m,params,mmin,mmax)+
	  (1.0-params[bgrfrac+offset])*signalmassdis(m,iterator->Getsigmamass(),params);
	sum = sum + mass;
      }
      ++iterator;
    }
  }    // next BDecay object
  sum = sum * mbinwidth;
  return sum;
}



void bgrctfcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *params, Int_t iflag)
{
  //----------------------------------------------------------------------
  // This is the function used by minuit to do the unbinned fit 
  // to the background ctau distribution.
  //----------------------------------------------------------------------
  // THIS FUNCTION CANNOT BE USED FOR A COMBINED FIT !
  //----------------------------------------------------------------------
  f=0.0;
  offset=0; // just to make sure
  Double_t bgr=0.0;
   std::vector<Measurement>::iterator bgr_iterator = junker->begin();
   while(bgr_iterator!= junker->end())
     {      
       bgr = bgrctdis(bgr_iterator->Getctau(),bgr_iterator->Getsigmactau(),params);
       f   = log(bgr)+f;
       ++bgr_iterator;
     }
   f=   -2.0 * f;
}

Double_t draw_bgrctfcn(Double_t *x,Double_t *params)
{
  //----------------------------------------------------------------------
  // This is the function used to display the result of the unbinned fit 
  // to the background ctau distribution.
  //----------------------------------------------------------------------
  // THIS FUNCTION CANNOT BE USED FOR A COMBINED FIT !
  //----------------------------------------------------------------------

  Double_t bgr;   
  Double_t sum = 0.0;
  Double_t l=*x;   
  std::vector<Measurement>::iterator bgr_iterator = junker->begin();
  while(bgr_iterator!= junker->end())
    {
      bgr = bgrctdis(l,bgr_iterator->Getsigmactau(),params);
      sum = sum + bgr;  
      ++bgr_iterator;
    }
  sum = sum * binwidth;
  return sum;
}


void mass_ct_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *params, Int_t iflag)
{
  //----------------------------------------------------------------------
  // This is the function used by minuit to do the unbinned fit 
  // to the signal ctau distribution.
  //----------------------------------------------------------------------
  int nn;
  double result;

  f = 0.0;
  
  if (combine_flg) {
    nn = count_bdecay;
  } else {
    nn = 1;  // Only 'loop' once
    Bsample = BtoJPSIK;
    offset = 0;
  }
  for (int ii=0; ii<nn; ii++) {
    if (combine_flg) {
      Bsample = BtoJPSIK_array[ii];  // set pointer to current sample
      offset = ii*SSpecPar;          // calculate offset of sample specific parameters
    }
    std::vector<Measurement>::iterator iterator = Bsample->beginiter();    
    while (iterator != Bsample->enditer()) {
      result = 
	params[bgrfrac+offset]*
	bgrctdis(iterator->Getctau(),iterator->Getsigmactau(),params)*
	bgrmassdis(iterator->Getmass(),params,Bsample->Getmin_mass(),Bsample->Getmax_mass())
	+(1.0-params[bgrfrac+offset])*signalctdis(iterator->Getctau(),iterator->Getsigmactau(),params)*signalmassdis(iterator->Getmass(),iterator->Getsigmamass(),params);
      f = log(result)+f;
      iterator++;
    }
  }
  f = -2.0*f;
  return;
}



Double_t draw_mass_ct_fcn(Double_t *x,Double_t *params)
{
  //----------------------------------------------------------------------
  // This is the function used to display the result of the unbinned fit 
  // to the signal ctau distribution.
  //----------------------------------------------------------------------
  int nn;
  Double_t result=0.0;
  Double_t sum=0.0;
  Double_t sum2=0.;   
  Double_t l=*x,sig;

  if (combine_flg) {
    nn = count_bdecay;
  } else {
    nn = 1;  // Only 'loop' once
    Bsample = BtoJPSIK;
    offset = 0;
  }
  for (int ii=0; ii<nn; ii++) {
    sum2=0.;
    if (combine_flg) {
      Bsample = BtoJPSIK_array[ii];  // set pointer to current sample
      offset = ii*SSpecPar;          // calculate offset of sample specific parameters
    }
    std::vector<Measurement>::iterator iterator = Bsample->beginiter();    
    while (iterator != Bsample->enditer()) {
      sig = iterator->Getsigmactau()*params[ctscale];
      result = params[bgrfrac+offset] * bgrctdis(l,sig,params) + 
	(1-params[bgrfrac+offset])*signalctdis(l,sig,params);
      sum    = sum + result; sum2=sum2+result;
      ++iterator;   
    }
  }
  sum = sum * binwidth;
  return sum;
}



Double_t draw_sigctdis(Double_t *x,Double_t *params)
{
  //----------------------------------------------------------------------
  // This is the function used to display the result of the unbinned fit 
  // to the signal ctau distribution.
  //----------------------------------------------------------------------
  Double_t result,sum;   
  Double_t l=*x,sig;
  int nn;
  
  result=0.;
  sum=0.;
  
  if (combine_flg) {
    nn = count_bdecay;
  } else {
    nn = 1;  // Only 'loop' once
    Bsample = BtoJPSIK;
    offset = 0;
  }
  for (int ii=0; ii<nn; ii++) {
    if (combine_flg) {
      Bsample = BtoJPSIK_array[ii];  // set pointer to current sample
      offset = ii*SSpecPar;          // calculate offset of sample specific parameters
    }
    std::vector<Measurement>::iterator iterator = Bsample->beginiter();    
    while (iterator != Bsample->enditer()) {
      sig = iterator->Getsigmactau();
      result = (1-params[bgrfrac+offset])*signalctdis(l,sig,params);
      //result = result * binwidth;
      sum    = sum + result;
      ++iterator;      
    }
  }
  sum = sum * binwidth;
  return sum;
}



int main(int argc, char **argv)
{
 static TROOT exclusivefit("exclusivefit","B lifetime fitting");
 static TRint app("app",&argc,argv,NULL,0);   
 int i,ii;
 int draw1=1,draw2=1,draw3=1,draw4=1,draw5=1,draw6=1,cdiv1=3,cdiv2=2;
 int cdraw1=1,cdraw2=1,cdraw3=1,cdraw4=1,ccdiv1=2,ccdiv2=2,doprint=0,doprintps=0;
 int csizx=1000,csizy=800,ccsizx=1000,ccsizy=800,binfac=1,binfit=1;
 float Ctau,Mass,Emass;
 double mmin,mmax,mmin2,mmax2;
 char datafile[100],confname[30],canvasname[10],canvastitle[30],psfn[20];
 Double_t stp[5],stp_n[5];  // startparameter for binned mass fits
 Double_t stp_l[8];  // startparameter for binned lifetime fit (background)


 // Read config file
 ierflg = conf_read("bexclfit.conf");
 if (ierflg != 0) {
   printf("ERROR reading config file !\n");
   exit(0);
 }          
 
 conf_readkeyword("DRAWHIS","%d %d %d %d %d %d",&draw1,&draw2,&draw3,&draw4,&draw5,&draw6);
 conf_readkeyword("CDIV","%d %d",&cdiv1,&cdiv2);
 conf_readkeyword("C_DRAWHIS","%d %d %d %d",&cdraw1,&cdraw2,&cdraw3,&cdraw4);
 conf_readkeyword("C_CDIV","%d %d",&ccdiv1,&ccdiv2);
 conf_readkeyword("CSIZE","%d %d",&csizx,&csizy);
 conf_readkeyword("C_CSIZE","%d %d",&ccsizx,&ccsizy);
 conf_readkeyword("PRINT","%d %d",&doprint,&doprintps);
 conf_readkeyword("BINFAC","%d",&binfac);
 conf_readkeyword("BINFIT","%d",&binfit);

 //conf_readkeyword("MASSWIN","%F",&masswin);
 //conf_readkeyword("SBMASSWIN","%F",&sbmasswin);

 // Read data file
 conf_readkeyword("DATAFILE","%s",datafile); // get data filename from configfile 
 TFile *f = new TFile(datafile);
 f->Print();
 f ->ls();

 // cout << " the system id is "<< gSystem->GetPid()<<endl;
 // gROOT->pwd();   // show current directory
 gROOT->cd();  // back to program scope 
 TIter next(f->GetListOfKeys());
 TKey *key;
 TObject *obj;

 // Check, how many different BDecay objects we have in the file
 while ((key = (TKey*)next())) {	 
   obj = (TObject*)f->Get(key->GetName());
   if(obj->InheritsFrom(BDecay::Class())) {    
     BtoJPSIK_array[count_bdecay] = (BDecay*)obj;
     BtoJPSIK_array[count_bdecay]->Print();       
     count_bdecay++;
   }
 }
 cout << char(27) << "[1m" << char(27) << "[34m";
 printf("\nNumber of subsamples in this file   %2d \n\n",count_bdecay);
 cout << char(27) << "[0m";

 //
 // now book some histograms:
 //
 TH1F *masshis     = new TH1F("masshis"     ,"mass distribution (ctau > 100 #mum)",70,5.15,5.4);
 //TH1F *masshis     = new TH1F("masshis"     ,"mass distribution",70,5.491,5.791);
 TH1F *normasshis  = new TH1F("normasshis"  ,"normalized mass distribution (ctau>100 #mum)",70,-masscut,masscut);
 //TH1F *ubmasshis   = new TH1F("ubmasshis"   ,"mass distribution (c#tau>100 #mum)",60,5.129,5.429);
 TH1F *ubmasshis   = new TH1F("ubmasshis"   ,"",60,5.129,5.429);
 //TH1F *ubmasshis   = new TH1F("ubmasshis"   ,"mass distribution (ctau>100 micron)",70,5.491,5.791);
 TH1F *fctubmasshis= new TH1F("fctubmasshis","fitted mass distribution (unbinned)",480,5.129,5.429);
 
 TH1F *ctauhis     = new TH1F("ctauhis"     ,"signal c#tau  distribution",140,-0.3,0.4);
 //TH1F *ubctauhis   = new TH1F("ubctauhis"   ,"signal c#tau  distribution",140,-0.3,0.4);
 TH1F *ubctauhis   = new TH1F("ubctauhis"   ,"",140,-0.3,0.4);
 TH1F *lifeufcnhis = new TH1F("lifeufcnhis ","b contribution to fitted ctau  distribution  (unbinned)",70*binfac,-0.3,0.4);
 TH1F *l2fcnhis = new TH1F("l2fcnhis ","b contribution to fitted ctau  distribution  (unbinned)",140,-0.3,0.4);
 TH1F *sigufcnhis  = new TH1F("sigufcnhis"  ,"fitted c#tau  distribution (unbinned)",70*binfac,-0.3,0.4);
 TH1F *ctaubgr     = new TH1F("ctaubgr"     ,"bgr c#tau  distribution",140,-0.3,0.4);
 //TH1F *ubctaubgr   = new TH1F("ubctaubgr"   ,"bgr c#tau  distribution",140,-0.3,0.4);
 TH1F *ubctaubgr   = new TH1F("ubctaubgr"   ,"",140,-0.3,0.4);
 //TH1F *ubctaubgr   = new TH1F("ubctaubgr"   ,"bgr c#tau  distribution",140,-0.3,0.4);

 TH1F *backufcnhis = new TH1F("backufcnhis" ,"fitted bgr c#tau  distribution (unbinned)",70*binfac,-0.3,0.4);

 // For syst. errors
 TH1F *bgrsubctauhis = new TH1F("bgrsubctauhis","",140,-0.3,0.4);
 TH1F *nullhis = new TH1F("nullhis","",140,-0.3,0.4);
 TH1F *leftctaubgr = new TH1F("leftbgr","left bgr",140,-0.3,0.4);
 TH1F *rightctaubgr = new TH1F("rightbgr","right bgr",140,-0.3,0.4);
 TH1F *gausshis = new TH1F("gausshis","b contribution to fitted ctau  distribution  (unbinned)",1400,-0.3,0.4);

 gStyle->SetCanvasColor(0);
 gROOT->SetStyle("Plain"); // Just don't draw any Pads and fancy stuff

 //
 // ********** Loop over samples and fit seperately *********
 //
 combine_flg = 0;  // set flag for fit of single samples

 for (int ii=0; ii<count_bdecay; ii++) {
   cout << char(27) << "[1m" << char(27) << "[31m";
   printf("\n\n**********************************************************\n");
   printf("* Processing sample No.  %2d ...                          *\n",ii);
   printf("**********************************************************\n");
   cout << char(27) << "[0m";

  
   BtoJPSIK= BtoJPSIK_array[ii];   // set pointer to current BDecay object
   mmin = BtoJPSIK->Getmin_mass();
   mmax = BtoJPSIK->Getmax_mass();
   gROOT->pwd();   // show current directory
  
   cout << " Capacity of sample:  " << BtoJPSIK->Capacity()<< endl;
   cout << " size of sample    :  " << BtoJPSIK ->Size() << endl;
   cout << " Mass min, max     :  " << mmin << ", " << mmax << endl; 
   cout << " Masswindow        :  " << masswin << endl;
   cout << " Sideband masswin  :  " << sbmasswin << endl;

   mmin2 = PDGBmass-masswin;
   mmax2 = PDGBmass+masswin;
   cout << " mmin2 mmax2 = " << mmin2 << " " << mmax2 << endl;
   if (mmin2>mmin && mmax2<mmax) {
     BtoJPSIK ->Prune(mmin2,mmax2);
     mmin = BtoJPSIK->Getmin_mass();
     mmax = BtoJPSIK->Getmax_mass();
     cout << " Changed mass min, max     :  " << mmin << ", " << mmax << endl; 
   } else if (mmin2<mmin || mmax2>mmax) {
     cout << "Masswindow too big !!!" << endl;
     exit(1);
   }
   //cout << " Capacity of sample " << ii << " :  " << BtoJPSIK->Capacity()<< endl;
   //cout << " size of sample " << ii << " :      " << BtoJPSIK->Size() << endl;
   
   samplesize[ii] = BtoJPSIK ->Size();
   cnum = 0; // canvas pad number   
   // Reset histograms
   masshis->Reset();
   normasshis->Reset();
   ubmasshis->Reset();
   fctubmasshis->Reset();
   ctauhis->Reset();
   ubctauhis->Reset();
   lifeufcnhis->Reset();
   sigufcnhis->Reset();
   ctaubgr->Reset();
   ubctaubgr->Reset();
   backufcnhis->Reset();
   l2fcnhis->Reset();


   std::vector<Measurement>::iterator iterator = BtoJPSIK->beginiter();    
   while (iterator != BtoJPSIK->enditer()) {
     if (iterator->Getmass()>mmin && iterator->Getmass()<mmax) {
       Ctau  = iterator->Getctau();
       Mass  = iterator->Getmass();
       Emass = iterator->Getsigmamass();
       ctauhis  ->Fill(Ctau);
       ubctauhis->Fill(Ctau);
       //       if (Mass<(mmin+sbmasswin)||Mass>(mmax-sbmasswin)) { // fill histos for ctau sideband fit
       //  ctaubgr  ->Fill(Ctau);
       //  ubctaubgr->Fill(Ctau);
       //}
       if (Ctau>ctaucut) {                           // fill histos for mass fit
	 masshis  ->Fill(Mass);
	 ubmasshis->Fill(Mass);
	 if (TMath::Abs((Mass-PDGBmass)/Emass)<masscut){  // fill normalized mass distribution
	   normasshis->Fill((Mass-PDGBmass)/Emass);
	 }
       }
       if (Mass>5.2 && Mass<5.23) {
	 leftctaubgr->Fill(Ctau);
       }
       if (Mass<5.38 && Mass>5.35) {
	 rightctaubgr->Fill(Ctau);
       }
     }
       ++iterator;
   }

   junker =BtoJPSIK->Sidebands(mmin,mmax,sbmasswin);
   cout << " Capacity of sample:  " << junker->capacity()<< endl;
   cout << " size of sample    :  " << junker->size()<< endl;
   cout << " Mass min, max     :  " << mmin << ", " << mmax << endl; 
   std::vector<Measurement>::iterator bgr_iterator = junker->begin();
   while(bgr_iterator!= junker->end())
     {
	 ctaubgr  ->Fill(bgr_iterator->Getctau());
	 ubctaubgr->Fill(bgr_iterator->Getctau());
	 ++bgr_iterator;
     }
   //   delete junker;


   //       if (Mass<(mmin+sbmasswin)||Mass>(mmax-sbmasswin)) { // fill histos for ctau sideband fit
   //	 ctaubgr  ->Fill(Ctau);
   //	 ubctaubgr->Fill(Ctau);
   //    }
   // For systematic error
   ctauhis->Sumw2();
   ctaubgr->Sumw2();
   bgrsubctauhis->Add(ctauhis,ctaubgr,1.,-1.5);

   //bgrsubctauhis->Add(ctauhis,leftctaubgr,1.,-10.);

   // Create and subdivide canvas
   sprintf(canvastitle,"Fit results subsample %d",ii);
   sprintf(canvasname,"c%d",ii);
   c[ii] = new TCanvas(canvasname,canvastitle,200,10,csizx,csizy);
   c[ii]->Divide(cdiv1,cdiv2);
   if (draw1) c[ii]->cd(++cnum);

   //
   // First do a binned likelihood fit to the mass:
   //  

   // Read startparameters from config file
   sprintf(confname,"SP%d",ii+1);
   conf_readkeyword(confname,"%F %F %F %F %F",&stp[0],&stp[1],&stp[2],&stp[3],&stp[4]); 
   
   TF1 *massfit = new TF1("massfit","[0] + [1]*x+[2]*exp(-0.5*((x-[3])/[4])^2)"); 
   gStyle->SetOptFit(1011);
   gStyle->SetOptStat(111111);
   massfit->SetParNames("p0","p1","const","mean","sigma"); 
   massfit->SetParameters(stp[0],stp[1],stp[2],stp[3],stp[4]);
   if (binfit) 
     if (draw1) {masshis->Fit("massfit");} else {masshis->Fit("massfit","N");}
   c[ii]->Update();

   //
   // perform a binned fit to the normalized mass distribution:
   // the sigma of this fit serves as start value of the error scale factor
   // of the unbinned mass fit
   // 
   if (draw2) c[ii]->cd(++cnum);
   // Read startparameters from config file
   sprintf(confname,"NSP%d",ii+1);
   conf_readkeyword(confname,"%F %F %F %F %F",&stp_n[0],&stp_n[1],&stp_n[2],&stp_n[3],&stp_n[4]); 
  
   gStyle->SetOptFit(1011);
   gStyle->SetOptStat(11);
   TF1 *nmassfit = new TF1("nmassfit","[0] + [1]*x+[2]*exp(-0.5*((x-[3])/[4])^2)");
   nmassfit->SetParNames("p0","p1","const","mean","sigma");
   nmassfit->SetParameters(stp_n[0],stp_n[1],stp_n[2],stp_n[3],stp_n[4]);
   if (draw2) {normasshis->Fit("nmassfit");} else {normasshis->Fit("nmassfit","N");}
   c[ii]->Update();

   //
   // do an unbinned likelihood fit to the mass
   //
   gStyle->SetOptStat(0); // No statistic box
   if (draw3) c[ii]->cd(++cnum);
   Float_t nsignal = (massfit->GetParameter(2)*massfit->GetParameter(4)*sqrt2pi)/masshis->GetBinWidth(1);
   Float_t entries = (Float_t)(masshis->GetEntries());
   printf("NSIGNAL, NENTRIES = %f, %f\n",nsignal,entries);
   // Set start values . Use results from bined fit.
   sfpar[mpeak] = massfit->GetParameter(3);
   //sfpar[mscale] = nmassfit->GetParameter(4);
   sfpar[mscale] = 1.2;
   sfpar[bgrfrac] = (entries-nsignal)/entries;
   //sfpar[bgrfrac] = 0.8;
   sfpar[mconst] = 15.;

   TMinuit *gmMinuit = new TMinuit(13);  // CAREFUL !  Hardcoded 13 !  
   gmMinuit->SetFCN(massfcn);

   for (i = 0; i<13; i++) {   // CAREFUL !  Hardcoded 13 !  
     gmMinuit->mnparm(i,par_name[i],sfpar[i], step[i], 0,0,ierflg);
   }
   
   gmMinuit->FixParameter(lbdapos);
   gmMinuit->FixParameter(fracpos);
   gmMinuit->FixParameter(lbdaneg);
   gmMinuit->FixParameter(fracneg);
   gmMinuit->FixParameter(lbdapos2);
   gmMinuit->FixParameter(fracpos2);
   gmMinuit->FixParameter(ctscale);
   gmMinuit->FixParameter(ctau);
   gmMinuit->FixParameter(tailfrac);
   
   gmMinuit->Migrad(); 

   // Store fit status of mass fit   
   strcpy(mfitstat[ii],gMinuit->fCstatu); strcat(mfitstat[ii],""); 

                       
   for (i=0; i<13; i++) {
     gmMinuit->GetParameter(i,sfpar[i],errsfpar[i]);
   }
   // Store the background fraction for later use in combined mass fit
   m_bgrfrac[ii] = sfpar[bgrfrac];
   
   gmMinuit->GetParameter(mconst,sspars[ii][mconst],sserrpars[ii][mconst]);
   gmMinuit->GetParameter(mpeak,sspars[ii][mpeak],sserrpars[ii][mpeak]);
   gmMinuit->GetParameter(mscale,sspars[ii][mscale],sserrpars[ii][mscale]);

   //delete gmMinuit;

   // Now draw the result
   mbinwidth = ubmasshis->GetBinWidth(1);
   massufcn  = new TF1("massufcn",draw_massfcn,mmin,mmax,13);
   //massufcn->SetParNames("peak","sigma","f1","f2");
   massufcn->SetParameters(sfpar);
   massufcn->SetFillColor(38);
   massufcn->SetFillStyle(1001);
   massufcn->SetLineWidth(2);
   ubmasshis->SetXTitle("m_{B} [GeV/c^{2}]");
   ubmasshis->SetYTitle("Entries per 5 MeV/c^{ 2}");
   ubmasshis->SetLabelSize(0.06,"X"); ubmasshis->SetLabelSize(0.06,"Y");
   ubmasshis->SetTitleSize(0.06,"X");  ubmasshis->SetTitleSize(0.06,"Y");
   ubmasshis->SetTitleOffset(1.1,"X");
   ubmasshis->SetTitleOffset(0.7,"Y");
   ubmasshis->SetMaximum(ubmasshis->GetMaximum()+.1*ubmasshis->GetMaximum());
   
   if (draw3) ubmasshis->Draw();
   //fctubmasshis->SetLineColor(4);
   fctubmasshis->SetLineColor(1);
   fctubmasshis->Eval(massufcn, "A");
   if (draw3) fctubmasshis->Draw("LSAME");
   c[ii]->Update();
   
   //
   // perform a binned likelihood fit to the ct bgr distribution:
   //  
   if (draw4) {c[ii]->cd(++cnum);
   gPad->SetLogy();}

   sprintf(confname,"BCTSP%d",ii+1);
   conf_readkeyword(confname,"%F %F %F %F %F %F %F %F",&stp_l[0],&stp_l[1],&stp_l[2],&stp_l[3],&stp_l[4],&stp_l[5],&stp_l[6],&stp_l[7]);

   /* backfit = new TF1("backfit",backgr,-.3,.4,6);
   backfit->SetParNames("Const","sigma","lambdaplus","fplus","lambdaminus","fminus");
   backfit->SetParameters(stp_l[0],stp_l[1],stp_l[2],stp_l[3],stp_l[4],stp_l[5]);
   ctaubgr->Fit("backfit","LE");
   */
   delete backfit; // If don't delete, the fit doesn't work the second time !
   backfit = new TF1("backfit",backgr,-.3,.4,8);
   backfit->SetParNames("Const","sigma","lambdaplus","fplus","lambdaminus","fminus","lamdaplus2","fplus2");
   backfit->SetParameters(stp_l[0],stp_l[1],stp_l[2],stp_l[3],stp_l[4],stp_l[5],stp_l[6],stp_l[7]);
   if (binfit)
     if (draw4) {ctaubgr->Fit("backfit","LE");} else {ctaubgr->Fit("backfit","LEN");}
   // Get start parameters for unbinned fit
   sfpar[lbdaneg]  = backfit->GetParameter(4);
   sfpar[fracneg]  = backfit->GetParameter(5);
   sfpar[lbdapos]  = backfit->GetParameter(2);
   sfpar[fracpos]  = backfit->GetParameter(3);
   sfpar[lbdapos2] = backfit->GetParameter(6);
   sfpar[fracpos2] = backfit->GetParameter(7);
   sfpar[ctscale]  = 1.1;
   
   c[ii]->Update();
   if (draw5) {c[ii]->cd(++cnum);
   gPad->SetLogy();}
   //
   // now estimate the background fraction and normalize the sideband 
   // distributions to the signal mass window (scale).
   //
   
   float scale = masswin/sbmasswin;  
   sfpar[bgrfrac] = (scale*ctaubgr->GetEntries())/(ubctauhis->GetEntries()); 
   //sfpar[bgrfrac] = 0.8;
   


   gmMinuit = new TMinuit(13);   // CAREFUL : Hardcoded 13    
   gmMinuit->SetFCN(bgrctfcn);
   //  
   // Set starting values and step sizes for parameters
   //
   for (int i = 0;i<13;i++) {  // CAREFUL : Hardcoded 13 
       gmMinuit->mnparm(i,par_name[i],sfpar[i],step[i],0,0,ierflg);
   }

   gmMinuit->FixParameter(bgrfrac);
   gmMinuit->FixParameter(mpeak);
   gmMinuit->FixParameter(mscale);
   gmMinuit->FixParameter(mconst);
   gmMinuit->FixParameter(ctau);   
   gmMinuit->FixParameter(tailfrac);
   gmMinuit->Migrad();

   // Store fit status of bgr fit   
   strcpy(bfitstat[ii],gMinuit->fCstatu); strcat(bfitstat[ii],""); 

   for (int i = 0;i<13;i++) {
     gmMinuit->GetParameter(i,sfpar[i],errsfpar[i]);
   } 

   binwidth = ubctaubgr->GetBinWidth(1);
   backufcn = new TF1("backufcn",draw_bgrctfcn,-.3,.4,dim);
   backufcn->SetParNames("lambdaplus","fplus","lambdaminus","fminus","scale");
   backufcn->SetParameters(sfpar);
   ubctaubgr->SetXTitle("c#tau [cm]");
   ubctaubgr->SetYTitle("Entries per 0.1 mm");
   ubctaubgr->SetLabelSize(0.06,"X"); ubctaubgr->SetLabelSize(0.06,"Y");
   ubctaubgr->SetTitleSize(0.06,"X"); ubctaubgr->SetTitleSize(0.06,"Y");
   ubctaubgr->SetTitleOffset(1.1,"X"); ubctaubgr->SetTitleOffset(0.7,"Y"); 
   ubctaubgr->SetOption("PE1");
   ubctaubgr->SetMarkerColor(2);
   ubctaubgr->SetMarkerSize(0.4);
   ubctaubgr->SetMarkerStyle(22);
    if (draw5) ubctaubgr->Draw();
   backufcnhis->Eval(backufcn, "A");
   backufcnhis->SetLineWidth(1);
   backufcnhis->SetLineColor(1);
   if (draw5) backufcnhis->Draw("HISTSSAME");
   c[ii]->Update();
   delete junker;
   
   gmMinuit = new TMinuit(13); // CAREFUL: Hardcoded 13 !! 
   gmMinuit->mncler();
   gmMinuit->SetFCN(mass_ct_fcn);
   //
   // Set starting values and step sizes for parameters
   //
   sfpar[ctau] = 0.045;
   //sfpar[ctau] = 0.037;
   sfpar[tailfrac] = 0.;

   
   // Stuff to test uncertainty due to wrong background parametrization
   //   sfpar[fracpos] = sfpar[fracneg]; // Set left and right side bgr exp. equal
   //sfpar[lbdapos] = sfpar[lbdaneg]; //  only works for fixed parameters
   //sfpar[lbdapos2] = 0.03; // wrong !!!
   //sfpar[fracneg]=0.05;

   //sfpar[ctscale] = 1.1;

   
   

   for (int i = 0;i<13;i++) {
       gmMinuit->mnparm(i,par_name[i],sfpar[i],step[i],0,0,ierflg);
   } 
   //
   // First we vary only ctau and keep everything else fixed
   //
   for (int i = 0;i<13;i++)  { 
       gmMinuit->FixParameter(i);
   }
   gmMinuit->Release(ctau);
   gmMinuit->SetMaxIterations(100);
   gmMinuit->Migrad();
   gmMinuit->Release(ctscale);
   gmMinuit->Migrad();
   gmMinuit->Release(bgrfrac);
   gmMinuit->Migrad();
   //
   // Then we vary everything :=)  ... or not.
   //
   gmMinuit->Release(mpeak);
   gmMinuit->Release(mscale);
   gmMinuit->Release(mconst);
   gmMinuit->Release(lbdapos);
   gmMinuit->Release(fracpos);
   gmMinuit->Release(lbdaneg);
   gmMinuit->Release(fracneg);
   gmMinuit->Release(lbdapos2);
   gmMinuit->Release(fracpos2);
   //  gmMinuit->Release(tailfrac);
   gmMinuit->Migrad();
   
   // Store fit status of main fit   
   strcpy(fitstat[ii],gMinuit->fCstatu); strcat(fitstat[ii],""); 
   
   for (int i = 0;i<13;i++) {
     gmMinuit->GetParameter(i,sfpar[i],errsfpar[i]);
   } 
   
   gmMinuit->GetParameter(ctau,sspars[ii][ctau],sserrpars[ii][ctau]);
   gmMinuit->GetParameter(ctscale,sspars[ii][ctscale],sserrpars[ii][ctscale]);
   gmMinuit->GetParameter(lbdapos,sspars[ii][lbdapos],sserrpars[ii][lbdapos]);
   gmMinuit->GetParameter(fracpos,sspars[ii][fracpos],sserrpars[ii][fracpos]);
   gmMinuit->GetParameter(lbdaneg,sspars[ii][lbdaneg],sserrpars[ii][lbdaneg]);
   gmMinuit->GetParameter(fracneg,sspars[ii][fracneg],sserrpars[ii][fracneg]);
   gmMinuit->GetParameter(lbdapos2,sspars[ii][lbdapos2],sserrpars[ii][lbdapos2]);
   gmMinuit->GetParameter(fracpos2,sspars[ii][fracpos2],sserrpars[ii][fracpos2]);
   gmMinuit->GetParameter(bgrfrac,sspars[ii][bgrfrac],sserrpars[ii][bgrfrac]);
   gmMinuit->GetParameter(tailfrac,sspars[ii][tailfrac],sserrpars[ii][tailfrac]);

   //double emat [dim*dim];
   //
   // just for later if we want to store the covariance matrix 
   //
   // gmMinuit->mnemat(emat,dim);
   // for (int i = 0;i<10;i++)
   //   {
   //     for (int j = 0;j<10;j++) 
   //	 { 
   //	   cout <<"cov["<<i<<","<<j<<"]=   " << emat [i+j*dim] << endl;
   //	   cout << " " <<emat [i+j*dim]; 
   // }
   //   cout <<endl;
   // }
   //cout << "emat"<<emat<<endl;

   //FitRes *result = new FitRes("name","name","bfcn",par_name,params,errparams,ii);
   //result->PrintFitRes();
   
   //  gmMinuit->Command("scan 10");
   

   // Now draw everything into one pad
   if (draw6) {c[ii]->cd(++cnum);
   gPad->SetLogy();}
   sigufcn = new TF1("sigufcn",draw_mass_ct_fcn,-.3,.4,dim);
   sigufcn->SetParNames("ctau","scale");
   sigufcn->SetParameters(sfpar);
   sigufcn->SetLineWidth(1);

   ubctauhis->SetOption("PE1");
   ubctauhis->SetMarkerColor(2);
   ubctauhis->SetMarkerSize(0.4);
   ubctauhis->SetMarkerStyle(22);
   ubctauhis->SetFillStyle(1001);
   if (draw6) ubctauhis->Draw("PE1");
   
   lifeufcn = new TF1("lifeufcn",draw_sigctdis,-.3,.4,dim);
   lifeufcn->SetParNames("ctau","scale");
   lifeufcn->SetParameters(sfpar);
   lifeufcn->SetLineWidth(1);
   ubctauhis->SetXTitle("c#tau [cm]");
   ubctauhis->SetYTitle("Entries per 0.1 mm");
   ubctauhis->SetLabelSize(0.06,"X"); ubctauhis->SetLabelSize(0.06,"Y");
   ubctauhis->SetTitleSize(0.06,"X"); ubctauhis->SetTitleSize(0.06,"Y");
   ubctauhis->SetTitleOffset(1.1,"X"); ubctauhis->SetTitleOffset(0.7,"Y"); 
   ubctauhis->Draw("PE1HISTSAME");
   lifeufcnhis->Eval(lifeufcn, "A");
   l2fcnhis->Eval(lifeufcn, "A"); // same as lifeufcnhis but with other binning
   lifeufcnhis->SetFillColor(38);
   lifeufcnhis->SetFillStyle(3005);
   lifeufcnhis->SetLineWidth(1);
   lifeufcnhis->SetLineColor(4);
   if (draw6) lifeufcnhis->Draw("HISTSSAME");
   sigufcnhis->Eval(sigufcn, "A");
   sigufcnhis->SetLineWidth(1);
   sigufcnhis->SetLineColor(4);
   if (draw6) sigufcnhis->Draw("HISTSSAME");
   c[ii]->Update();

   // **** inserted macro for PhD-Thesis
   if (1) {
//=========Macro generated from canvas: c0/Fit results subsample 0
//=========  (Fri Apr 20 14:49:01 2001) by ROOT version3.00/00
   TCanvas *c0 = new TCanvas("c0", "Fit results subsample 0",200,10,800,700);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   c0->SetHighLightColor(2);
   c0->Range(0,0,1,1);
   c0->SetFillColor(0);
   c0->SetBorderMode(0);
   c0->SetBorderSize(2);
  
// ------------>Primitives in pad: c0_1
//   TPad *c0_1 = new TPad("c0_1", "c0_1",0.01,0.676667,0.99,0.99);
   TPad *c0_1 = new TPad("c0_1", "c0_1",0.01,0.51,0.44,0.99);
   c0_1->Draw();
   c0_1->cd();
   c0_1->Range(5.09148,-8.69472,5.46652,36.4253);
   c0_1->SetFillColor(0);
   c0_1->SetBorderMode(0);
   c0_1->SetBorderSize(2);
   c0_1->SetTopMargin(0.05);
   c0_1->SetLeftMargin(0.176);
   c0_1->SetRightMargin(0.085);
   c0_1->SetBottomMargin(0.13);
   c0_1->SetFrameBorderMode(0);
   
   ubmasshis->GetXaxis()->SetLabelSize(0.05);
   ubmasshis->GetXaxis()->SetTitleSize(0.05);
   ubmasshis->GetXaxis()->SetTitleOffset(1.2);
   ubmasshis->GetYaxis()->SetLabelSize(0.05);
   ubmasshis->GetYaxis()->SetTitleSize(0.05);
   ubmasshis->GetYaxis()->SetTitleOffset(2.);
   ubmasshis->GetYaxis()->SetLabelOffset(0.01);
   ubmasshis->SetLineWidth(2);
   ubmasshis->Draw("");
   
   fctubmasshis->SetLineWidth(2);
   fctubmasshis->Draw("lsame");
   c0_1->Modified();
   c0->cd();
  
// ------------>Primitives in pad: c0_2
   //TPad *c0_2 = new TPad("c0_2", "c0_2",0.01,0.343333,0.99,0.656667);
   TPad *c0_2 = new TPad("c0_2", "c0_2",0.45,0.51,0.99,0.99);
   c0_2->Draw();
   c0_2->cd();
   c0_2->Range(-0.387545,-1.19386,0.487545,3.43195);
   c0_2->SetFillColor(0);
   c0_2->SetBorderMode(0);
   c0_2->SetBorderSize(2);
   c0_2->SetLogy();
   c0_2->SetTopMargin(0.05);
   c0_2->SetBottomMargin(0.13);
   c0_2->SetRightMargin(0.09);
   c0_2->SetLeftMargin(0.);
   c0_2->SetFrameBorderMode(0);
   
   ubctaubgr->SetMarkerColor(2);
   ubctaubgr->SetMarkerStyle(20);
   ubctaubgr->SetMarkerSize(0.4);
   ubctaubgr->GetXaxis()->SetLabelSize(0.05);
   ubctaubgr->GetXaxis()->SetTitleSize(0.05);
   ubctaubgr->GetXaxis()->SetTitleOffset(1.2);
   ubctaubgr->GetYaxis()->SetLabelSize(0.05);
   ubctaubgr->GetYaxis()->SetTitleSize(0.05);
   ubctaubgr->GetYaxis()->SetTitleOffset(1.4);
   ubctaubgr->Draw("");
   
   backufcnhis->SetLineWidth(2);
   backufcnhis->Draw("histssame");
   c0_2->Modified();
   c0->cd();
  
// ------------>Primitives in pad: c0_3
   TPad *c0_3 = new TPad("c0_3", "c0_3",0.01,0.01,0.99,0.49);
   c0_3->Draw();
   c0_3->cd();
   c0_3->Range(-0.387545,-1.22542,0.487545,3.59513);
   c0_3->SetFillColor(0);
   c0_3->SetBorderMode(0);
   c0_3->SetBorderSize(2);
   c0_3->SetLogy();
   c0_3->SetTopMargin(0.05);
   c0_3->SetBottomMargin(0.13);
   c0_3->SetLeftMargin(0.078);
   c0_3->SetRightMargin(0.05);
   c0_3->SetFrameBorderMode(0);
   
   ubctauhis->SetMarkerColor(2);
   ubctauhis->SetMarkerStyle(20);
   ubctauhis->SetMarkerSize(.4);
   ubctauhis->GetXaxis()->SetLabelSize(0.05);
   ubctauhis->GetXaxis()->SetTitleSize(0.05);
   ubctauhis->GetXaxis()->SetTitleOffset(1.2);
   ubctauhis->GetYaxis()->SetLabelSize(0.05);
   ubctauhis->GetYaxis()->SetTitleSize(0.05);
   ubctauhis->GetYaxis()->SetTitleOffset(0.85);
   ubctauhis->Draw("pe1");
   
   
   
   lifeufcnhis ->SetFillColor(1);
   lifeufcnhis ->SetFillStyle(3005);
   lifeufcnhis ->SetLineColor(1);
   lifeufcnhis->SetLineWidth(2);
   lifeufcnhis ->Draw("histssame");
   
   sigufcnhis->SetLineColor(1);
   sigufcnhis->SetLineWidth(2);
   sigufcnhis->Draw("histssame");
   c0_3->Modified();
   c0->cd();
   c0->Modified();
   c0->cd();
}

   // *******************************************************************

   // For systematic error 
   /*
   nullhis->Add(bgrsubctauhis,l2fcnhis,1.,-1.);
   c[ii]->cd(5);
   gPad->SetLogy(0);
   bgrsubctauhis->Draw();
   lifeufcnhis->Draw("HISTSAME");
   c[ii]->cd(6);
   gPad->SetLogy(0);
   nullhis->Draw();
   */

   // don't forget to put all the necessary delete statements here 
   // just to be a good citizen.
   //delete massfit;
   //delete nmassfit;
   //delete backfit;
   //delete sigufcn;
   //delete lifeufcn;

   if (doprint) {
     sprintf(psfn,"sample%d.eps",ii);
     c[ii]->Print(psfn);
   }
   if (doprintps) {
     sprintf(psfn,"sample%d.ps",ii);
     c[ii]->Print(psfn);
   }

   /*
   TGraph *contour  = (TGraph*)gmMinuit->Contour(5, 0, 3) ;
   TCanvas *c5 = new TCanvas("c5","Contours",10,10,800,600);
   gStyle->SetPalette(1);
   //TH2D *conthis = new TH2D("","",100,5.27,5.282,100,.35,.48);
   contour->SetMarkerStyle(21);
   contour->SetTitle("");
   //contour->GetHistogram()->GetXaxis()->SetTitle("blah");
   //contour->Draw("alp");
   //conthis->SetXTitle("m(B) [GeV/c^{2}]");
   //conthis->SetYTitle("c#tau [cm]");
   //conthis->GetXaxis()->SetLabelSize(0.04);
   //conthis->GetXaxis()->SetTitleSize(0.04);
   //conthis->GetXaxis()->SetTitleOffset(1.2);
   //conthis->GetYaxis()->SetLabelSize(0.04);
   //conthis->GetYaxis()->SetTitleSize(0.04);
   //conthis->GetYaxis()->SetTitleOffset(0.85);
   //conthis->Draw();
   contour->Draw("al");
   //contour->Draw("alp");
   */
         
 } // Next BDecay object  ******* end of loop over samples *******
 
   
 if (count_bdecay < 2) {  // If there's only one sample, we're done.
   printsummary(0);
    
   app.Run();   
   return 0;
 }



 //
 // ********** Combined fit starts  here ! **********
//

 cout << char(27) << "[1m" << char(27) << "[31m";
 printf("\n\n**********************************************************\n");
 printf("***                 COMBINED FIT                       ***\n");
 printf("**********************************************************\n");
 cout << char(27) << "[0m";

 // Fill parameter array
 for (i=0; i<4; i++) {   // Common Parameters: Use results of 1st sample
   params[i] = sspars[0][i];
   printf("Parameter %d : %10.7f\n",i,params[i]);  // test printout
 }
 for (ii=0; ii<count_bdecay; ii++) {
   for (i=0; i<SSpecPar; i++) {
     params[CommonPar+ii*SSpecPar+i] = sspars[ii][CommonPar+i];
     printf("Parameter %d : %10.7f\n",CommonPar+ii*SSpecPar+i,params[CommonPar+ii*SSpecPar+i]);  // test printout
   }
 }

 BtoJPSIK= BtoJPSIK_array[0];
 mmin = BtoJPSIK->Getmin_mass();
 mmax = BtoJPSIK->Getmax_mass(); 
 
 // Reset histograms
 masshis->Reset();
 normasshis->Reset();
 ubmasshis->Reset();
 fctubmasshis->Reset();
 ctauhis->Reset();
 ubctauhis->Reset();
 lifeufcnhis->Reset();
 sigufcnhis->Reset();
 ctaubgr->Reset();
 ubctaubgr->Reset();
 backufcnhis->Reset();
 bgrsubctauhis->Reset();

 // Fill histos 
 for (ii=0; ii<count_bdecay; ii++) {
   BtoJPSIK= BtoJPSIK_array[ii]; 
   cout << " Capacity of sample " << ii << " :  " << BtoJPSIK->Capacity()<< endl;
   cout << " size of sample " << ii << " :      " << BtoJPSIK->Size() << endl;
   
   std::vector<Measurement>::iterator iterator = BtoJPSIK->beginiter();    
   while (iterator != BtoJPSIK->enditer()) {
     if (iterator->Getmass()>mmin&&iterator->Getmass()<mmax) {
       float Ctau  = iterator->Getctau();
       float Mass  = iterator->Getmass();
       float Emass = iterator->Getsigmamass();
       ctauhis  ->Fill(Ctau);
       ubctauhis->Fill(Ctau);
       if (Mass<(mmin+sbmasswin)||Mass>(mmax-sbmasswin)) {  // fill histos for ctau sideband fit
	 ctaubgr  ->Fill(Ctau);
	 ubctaubgr->Fill(Ctau);
       }
       if (Ctau>ctaucut) {                          // fill histos for mass fit
	 masshis  ->Fill(Mass);
	 ubmasshis->Fill(Mass);
	 if (TMath::Abs((Mass-PDGBmass)/Emass)<masscut) {   // fill normalized mass distribution
	   normasshis->Fill((Mass-PDGBmass)/Emass);
	 }
       }
     }
     ++iterator;
   }
 } // next subsample

 // For systematic error 
 bgrsubctauhis->Add(ctauhis,ctaubgr,1.,-1.5);

 gStyle->SetOptFit(1011);
 gStyle->SetOptStat(11);
 
 cc = new TCanvas("cc","Results of combined fits",200,10,ccsizx,ccsizy);
 cc->Divide(ccdiv1,ccdiv2);
 cnum=0;
 if (cdraw1) cc->cd(++cnum);
 //
 // First do a binned likelihood fit to the mass:
 //  
 TF1 *massfit = new TF1("massfit","[0] + [1]*x+[2]*exp(-0.5*((x-[3])/[4])^2)"); 
 gStyle->SetOptFit(1011);
 gStyle->SetOptStat(111111);
 massfit->SetParNames("p0","p1","const","mean","sigma");
 massfit->SetParameters(1.21839e+02, -1.93711e+01, 100., PDGBmass, 0.010);
 if (cdraw1) {masshis->Fit("massfit");} else {masshis->Fit("massfit","N");}

 //
 // perform a binned fit to the normalized mass distribution:
 // the sigma of this fit serves as start value of the error scale factor
 // of the unbinned mass fit
 // 
 cc->Update();
 if (cdraw2) cc->cd(++cnum);  
 gStyle->SetOptFit(1011);
 gStyle->SetOptStat(11);
 TF1 *nmassfit = new TF1("nmassfit","[0] + [1]*x+[2]*exp(-0.5*((x-[3])/[4])^2)");
 nmassfit->SetParNames("p0","p1","const","mean","sigma");
 nmassfit->SetParameters(6.0, 0.0, 37., 0.0,1.0);
 if (cdraw2) {normasshis->Fit("nmassfit");} else {normasshis->Fit("nmassfit","N");}
 if (cdraw2) normasshis->Draw();
 cc->Update();

 //
 // do an unbinned likelihood fit to the mass
 //
 gStyle->SetOptStat(0); // No statistic box
 combine_flg = 1; // Switch likelihood functions to combine mode

 if (cdraw3) cc->cd(++cnum);
 Float_t nsignal = (massfit->GetParameter(2)*massfit->GetParameter(4)*sqrt2pi)/masshis->GetBinWidth(1);
 Float_t entries = (Float_t)(masshis->GetEntries());

 params[bgrfrac] = (entries-nsignal)/entries;
 params[mpeak] = massfit->GetParameter(3);
 params[mscale] = nmassfit->GetParameter(4);
 //params[mconst] = 9.0;

 for (i=0; i<count_bdecay; i++) {
   params[mconst+i*SSpecPar] = sspars[i][mconst];
   // Use the stored backgound fractions from the massfits (ct>100um)
   params[bgrfrac+i*SSpecPar] = m_bgrfrac[i]; 
 }

 TMinuit *gmMinuit = new TMinuit(CommonPar+count_bdecay*SSpecPar); 
 gmMinuit->SetFCN(massfcn);
 
 for (i=0; i<(CommonPar+count_bdecay*SSpecPar); i++) {
   gmMinuit->mnparm(i,par_name[i],params[i], step[i], 0,0,ierflg);
 }
 for (i=0; i<count_bdecay; i++) {
   gmMinuit->FixParameter(lbdapos+i*SSpecPar);
   gmMinuit->FixParameter(fracpos+i*SSpecPar);
   gmMinuit->FixParameter(lbdaneg+i*SSpecPar);
   gmMinuit->FixParameter(fracneg+i*SSpecPar);
   gmMinuit->FixParameter(lbdapos2+i*SSpecPar);
   gmMinuit->FixParameter(fracpos2+i*SSpecPar);
 }
 gmMinuit->FixParameter(ctscale);
 gmMinuit->FixParameter(ctau);
 gmMinuit->Migrad(); 
 // store fit status
 strcpy(mcombfitstat,gMinuit->fCstatu); strcat(mcombfitstat,"");

 gmMinuit->GetParameter(bgrfrac,params[bgrfrac],errparams[bgrfrac]);
 gmMinuit->GetParameter(mpeak,params[mpeak],errparams[mpeak]);
 gmMinuit->GetParameter(mscale,params[mscale],errparams[mscale]);
 gmMinuit->GetParameter(mconst,params[mconst],errparams[mconst]);
 for (int i=0; i<(CommonPar+count_bdecay*SSpecPar); i++) {
   gmMinuit->GetParameter(i,params[i],errparams[i]);
 } 

 // Draw everything
 mbinwidth = ubmasshis->GetBinWidth(1);
 massufcn  = new TF1("massufcn",draw_massfcn,mmin,mmax,dim);
 //massufcn->SetParNames("peak","sigma","f1","f2");
 massufcn->SetParameters(params);
 massufcn->SetFillColor(38);
 massufcn->SetFillStyle(1001);
 massufcn->SetLineWidth(1);
 ubmasshis->SetXTitle("m(B) [GeV/c^{2}]");
 if (cdraw3) ubmasshis->Draw();
 fctubmasshis->Eval(massufcn, "A");
 if (cdraw3) fctubmasshis->Draw("LSAME");
 cc->Update();

     


 // bgrfrac was overwritten in the massfit. We get the result of the
 // single sample fits as start parameters
 for (i=0; i<count_bdecay; i++) {
   params[bgrfrac+i*SSpecPar] = sspars[i][bgrfrac];
 }
 params[ctau] = 0.046;  

 params[tailfrac] = 0.0;

 gmMinuit = new TMinuit(CommonPar+count_bdecay*SSpecPar); 
 gmMinuit->SetFCN(mass_ct_fcn);
 // Set starting values and step sizes for parameters
 for (int i = 0; i<(CommonPar+count_bdecay*SSpecPar); i++) {
   gmMinuit->mnparm(i,par_name[i],params[i],step[i],0,0,ierflg);
 } 
 
 // First we vary only ctau and keep everything else fixed
 for (int i=0; i<(CommonPar+count_bdecay*SSpecPar); i++)  { 
   gmMinuit->FixParameter(i);
 }
 gmMinuit->Release(ctau);

 gmMinuit->SetMaxIterations(100);
 gmMinuit->Migrad();

 gmMinuit->Release(ctscale);
 gmMinuit->Migrad();

 for (int i=0; i<count_bdecay; i++)   // release bgrfracs
   gmMinuit->Release(bgrfrac+i*SSpecPar);  
 gmMinuit->Migrad();

 //
 // Then we vary everything :=)  ... or not.
 //
 gmMinuit->Release(mpeak);
 gmMinuit->Release(mscale);
 // gmMinuit->Release(tailfrac);  // Fraction of second gaussian
 for (int i=0; i<count_bdecay; i++) {
   gmMinuit->Release(mconst+i*SSpecPar);
   gmMinuit->Release(lbdapos+i*SSpecPar);
   gmMinuit->Release(fracpos+i*SSpecPar);
   gmMinuit->Release(lbdaneg+i*SSpecPar);
   gmMinuit->Release(fracneg+i*SSpecPar);
   gmMinuit->Release(lbdapos2+i*SSpecPar);
   gmMinuit->Release(fracpos2+i*SSpecPar);
 }
 gmMinuit->Migrad();

 strcpy(combfitstat,gMinuit->fCstatu); strcat(combfitstat,"");

 for (int i=0; i<(CommonPar+count_bdecay*SSpecPar); i++) {
   gmMinuit->GetParameter(i,params[i],errparams[i]);
 } 

 // Now draw everything into one pad
 if (cdraw4) {cc->cd(++cnum);
 gPad->SetLogy();}
 sigufcn = new TF1("sigufcn",draw_mass_ct_fcn,-.3,.4,dim);
 //sigufcn->SetParNames("ctau","scale");
 sigufcn->SetParameters(params);
 sigufcn->SetLineWidth(1);
 

 ubctauhis->SetMarkerColor(2);
 ubctauhis->SetMarkerStyle(20);
 ubctauhis->SetMarkerSize(.4);
 ubctauhis->GetXaxis()->SetLabelSize(0.038);
 ubctauhis->GetXaxis()->SetTitleSize(0.038);
 ubctauhis->GetXaxis()->SetTitleOffset(1.2);
 ubctauhis->GetYaxis()->SetLabelSize(0.038);
 ubctauhis->GetYaxis()->SetTitleSize(0.038);
 ubctauhis->GetYaxis()->SetTitleOffset(1.1);
 ubctauhis->SetFillStyle(1001);
 ubctauhis->SetXTitle("c#tau [cm]");
 if (cdraw4) ubctauhis->Draw("PE1");
 
 lifeufcn = new TF1("lifeufcn",draw_sigctdis,-.3,.4,dim);
 lifeufcn->SetParNames("ctau","scale");
 lifeufcn->SetParameters(params);
 lifeufcn->SetLineWidth(1);
 if (draw4) ubctauhis->Draw("PE1HISTSAME");
 lifeufcnhis->Eval(lifeufcn, "A");
 lifeufcnhis ->SetFillColor(1);
 lifeufcnhis ->SetFillStyle(3005);
 lifeufcnhis ->SetLineColor(1);
 lifeufcnhis->SetLineWidth(1);
 if (cdraw4) lifeufcnhis->Draw("HISTSSAME");
 sigufcnhis->Eval(sigufcn, "A");
 sigufcnhis->SetLineWidth(1);
 sigufcnhis->SetLineColor(1);
 if (cdraw4) sigufcnhis->Draw("HISTSSAME");
 cc->Update();

 if (doprint) {
   cc->Print("combined.eps");
 }
 if (doprintps) {
   cc->Print("combined.ps");
 }


 // For systematic error
 /* 
 gaussfcn = new TF1("gauss","7/(6.2831854*0.10)*exp(-x*x/0.016)",-.3,.4);
 gausshis->Eval(gaussfcn, "A");
 cs = new TCanvas("cs","",200,10,600,250);
 cs->Divide(2,1);
 //bgrsubctauhis->Sumw2();
 l2fcnhis->Sumw2();
 //nullhis->Sumw2();
 nullhis->Add(bgrsubctauhis,l2fcnhis,1.,-1.);
 nullhis->SetMarkerColor(2);
 nullhis->SetMarkerStyle(20);
 nullhis->SetMarkerSize(.4);
 cs->cd(1);
 gPad->SetLogy(0);
 bgrsubctauhis->SetMarkerColor(2);
 bgrsubctauhis->SetMarkerStyle(20);
 bgrsubctauhis->SetMarkerSize(.4);
 bgrsubctauhis->GetXaxis()->SetLabelSize(.048);
 bgrsubctauhis->GetYaxis()->SetLabelSize(.048);
 bgrsubctauhis->GetXaxis()->SetTitleSize(.048);
 bgrsubctauhis->GetYaxis()->SetTitleSize(.048);
 bgrsubctauhis->GetYaxis()->SetTitle("entries per 0.1 mm");
 bgrsubctauhis->GetXaxis()->SetTitle("c#tau [cm]");
 bgrsubctauhis->Draw("PE1");
 lifeufcnhis->Draw("HISTSAME");
 cs->cd(2);
 gPad->SetLogy(0);
 nullhis->GetXaxis()->SetLabelSize(.048);
 nullhis->GetYaxis()->SetLabelSize(.048);
 nullhis->GetXaxis()->SetTitleSize(.048);
 nullhis->GetYaxis()->SetTitleSize(.048);
 nullhis->GetYaxis()->SetTitle("entries per 0.1 mm");
 nullhis->GetXaxis()->SetTitle("c#tau [cm]");
 nullhis->Draw("PE1");
 gausshis->Draw("HISTSAME");
 */
 
 /*
 TGraph *contour  = (TGraph*)gmMinuit->Contour(50, 0, 3) ;
 TCanvas *c5 = new TCanvas("c5","",10,10,800,600);
 gStyle->SetPalette(1);
 contour->SetMarkerStyle(20);
 contour->SetMarkerSize(0.1);
 contour->SetTitle(""); 
 contour->Draw("al");  
 */
 

 printsummary(1); // print final result
 
                   
 

 app.Run();   
 return 0;
}
