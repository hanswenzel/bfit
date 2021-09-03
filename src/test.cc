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
#include <memory>
#include <TROOT.h>
#include <TFile.h>
#include <TRint.h>
#include <TH1.h>
#include <TF1.h>
#include <TF2.h>
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


int count_bdecay=0,samplesize[10];
BDecay  *BtoJPSIK_array[10],*Bsample,*BtoJPSIK;
std::vector<Measurement> *junker;

// Declare global variables also known within bexclfit_fcts.cc 
int combine_flg; // 0 = single sample fit; 1 = combined fit
int offset;      // This is always 0 in the single sample fit loop, will be used
                 // in the combined fit

int main(int argc, char **argv)
{
 static TROOT exclusivefit("exclusivefit","B lifetime fitting");
 static TRint app("app",&argc,argv,NULL,0);   
 TH1F *ubctaubgr   = new TH1F("ubctaubgr"   ,"",140,-0.3,0.4);
 TH1F *masshis     = new TH1F("masshis"     ,"mass distribution (ctau > 100 #mum)",70,5.15,5.4);
 double mmin,mmax;
 TFile *f = new TFile("data/BtoP2SJPsiK_3_4.5_0.75.root");
 f->Print();
 f ->ls();
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
 std::cout << char(27) << "[1m" << char(27) << "[34m";
 printf("\nNumber of subsamples in this file   %2d \n\n",count_bdecay);
 std::cout << char(27) << "[0m";
 for (int ii=0; ii<count_bdecay; ii++) {
   std::cout << char(27) << "[1m" << char(27) << "[31m";
   printf("\n\n**********************************************************\n");
   printf("* Processing sample No.  %2d ...                          *\n",ii);
   printf("**********************************************************\n");
   std::cout << char(27) << "[0m";
   BtoJPSIK= BtoJPSIK_array[ii];   // set pointer to current BDecay object
   mmin = BtoJPSIK->Getmin_mass();
   mmax = BtoJPSIK->Getmax_mass();
   gROOT->pwd();   // show current directory
   std::cout << " Capacity of sample:  " << BtoJPSIK->Capacity()<< std::endl;
   std::cout << " size of sample    :  " << BtoJPSIK ->Size() << std::endl;
   std::cout << " Mass min, max     :  " << mmin << ", " << mmax << std::endl; 
   //   BtoJPSIK ->Begin();
   //   std::vector<Measurement>::iterator iter3 = BtoJPSIK->enditer();
   //std::vector<Measurement>::iterator iterator = BtoJPSIK->beginiter();
   //   std::cout<<"test:  "<< iter3->Getctau()<< std::endl;
   //   std::cout<<"test2:  "<< iter4->Getctau()<< std::endl;
   //std::cout<<"sorting:  " << std::endl;
   //BtoJPSIK->Sort();
   //  while(iterator!= BtoJPSIK->enditer()) {
   //    std::cout << iterator->Getmass()<<std::endl;
   //    iterator++;
   //	   }
   //BtoJPSIK->Prune(5.135,5.135);   
   mmin = BtoJPSIK->Getmin_mass();
   mmax = BtoJPSIK->Getmax_mass();
   //   auto_ptr<std::vector<Measurement> > junker(BtoJPSIK->Sidebands(mmin,mmax,sbmasswin));
   junker =BtoJPSIK->Sidebands(mmin,mmax,sbmasswin);
   std::vector<Measurement>::iterator iterator = junker->begin();
   while(iterator!= junker->end())
     {
	 ubctaubgr->Fill(iterator->Getctau());
         masshis->Fill(iterator->Getmass());
	 iterator++;
     }

   gROOT->pwd();   // show current directory
   std::cout << " Capacity of sample:  " << junker->capacity()<< std::endl;
   std::cout << " size of sample    :  " << junker->size()<< std::endl;
   std::cout << " Mass min, max     :  " << mmin << ", " << mmax << std::endl; 
   //   cc = new TCanvas("cc","Results of combined fits",200,10,ccsizx,ccsizy);   
    masshis->Draw();
//   std::cout<<"test2: "<< *iter3 << std::endl;
   
   //   while (BtoJPSIK ->iter != BtoJPSIK->enditer()) {
   //  std::cout <<  BtoJPSIK->iter->Getctau()<<std::endl;
   //  BtoJPSIK->iter++;
   //}
 }
 app.Run();   
 return 0;
}

