//--------------------------------------------------------------------------
// File and Version Information:  
// mkbdecay.cc V2.0  June, 14th 2000
//
// Description:
// Convert root-Ntuples into BDecay-Object files and apply cuts
//
// Environment:
//      Software developed for the CDF Collaboration
//
// Author List:
//         Hans Wenzel    IEKP Karlsruhe and Fermilab 
//                        email::wenzel@fnal.gov
//         Andreas Heiss  IEKP Karlsruhe
//
//------------------------------------------------------------------------
#define VERSION "3.1"

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRint.h>
#include <TH1.h>
#include <TF1.h>
#include <TF2.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TMinuit.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <vector>
#include "h1000.h"
#include "Measurement.h"
#include "BDecay.h"
extern "C" {
#include <libconf.h>
}    
using namespace std;

// 
// The following are some predefined cuts and constants
//
const Double_t PDGBmass    = 5.279;           // PDG value of B-meson mass in GeV
//const Double_t PDGBmass = 5.624;
// Set some defaults
Double_t masscut     = 8.0;  // cut for normalized mass distribution (# of sigmas) 
Double_t masswin     = 0.250;       // width of mass window in GeV
Float_t Chiprob_cut  = 0.01;        // chi^2 propability cut
Float_t Chi2d_cut    = 100.;        // chi^2 (2D) cut
Float_t Pt_cut       = 4.5;         // Cut on B pt
Float_t Kpt_cut      = 1.0;         // Cut on Kaon pt
Float_t Ectau_cut    = 0.0250;      // Cut on ectau  
Float_t lmmin        = 5.229;       // signal region for listfile
Float_t lmmax        = 5.329;       //   "      "     "     "
Float_t lctau        = 0.0;         // minimum ctau for list
int nsvxhits = 2;
int checkhits = 0;
int checksvxacc = 0;
Double_t mmin,mmax;
FILE *lf;

// offset and combine_flg are only needed for compatibility with BDecay.
// This sucks and should be fixed soon !
int offset,combine_flg;
// same for sbmasswin !
Double_t sbmasswin;

TTree *newtree,*tree ;
char cutname[200],filename[250],conffile[250],oname[20];
char outfile[200],listfile[200],buffer[100]; 
//----------------------------------------------
Float_t Run, Event;
Float_t Chiprob, Chi2d;
Float_t Mass;
Float_t Emass;
Float_t Ctau;
Float_t Ectau;
Float_t Isolz;
Float_t Bitcode;
Float_t Pt;
Float_t Kpt;
Float_t Svx1,Svx2;

//Float_t Kchrg;
//----------------------------------------------

int main(int argc, char **argv)
{
  Int_t bits,i,nbad;
  int ierr,bcode,clist[100],cbad[100],lindex,j,flag,lflag;
  int l0,l1,l2,l3,svxacc,svxh1,svxh2;

  static TROOT mkBDecay("mkBDecay","Convert ntuples into BDecay objects file");
  static TRint app("app",&argc,argv,NULL,0);


  printf("\n%c[31mmkBDecay Version %s %c[0m\n",27,VERSION,27);
  strcpy(filename,""); strcpy(conffile,"");
  strcpy(oname,"Bstuff"); strcpy(outfile,"BDecay.root");
  for (i=0; i<100; i++) {cbad[i]=0; clist[i]=0;}
  i=0; 
  // Command line arguments
  while(i<(argc-1)) {
    i++;
    if (!strncmp(argv[i],"-c",2) && i<(argc-1)) {
      strcpy(conffile,argv[++i]);
      continue;
    }   
    if (strncmp(argv[i],"-",1)) strcpy(filename,argv[i]);
  }

  if (strlen(filename)<2) {
    printf("\nFilename of root ntuple file: ");
    scanf("%s",filename);
  }
  if (strlen(conffile)<2) strcpy(conffile,"mkbdecay.conf");

  printf("\nNtuple   : %c[34m%s%c[0m  \nConffile : %c[34m%s%c[0m\n\n",27,filename,27,27,conffile,27);

  // Read config file
  ierr = conf_read(conffile);
  if (ierr != 0) {
    printf("ERROR reading config file !\n");
    exit(0);
  }          
  // Now read some stuff. Ignore error code.
  conf_readkeyword("BITCODE","%d",&bcode);
  conf_readkeyword("MASSWIN","%F",&masswin);
  conf_readkeyword("OBJECTNAME","%s",oname);
  conf_readkeyword("OUTFILE","%s",outfile);
  conf_readkeyword("CHISQR","%f %f",&Chiprob_cut,&Chi2d_cut);
  conf_readkeyword("PTB","%f",&Pt_cut);
  conf_readkeyword("PTK","%f",&Kpt_cut);
  conf_readkeyword("SIGMASS","%f %f",&lmmin,&lmmax);
  conf_readkeyword("SIGCTAU","%f",&lctau);
  conf_readkeyword("ECTAU","%f",&Ectau_cut);
  conf_readkeyword("NSVXHITS","%d",&nsvxhits);
  conf_readkeyword("CHECKHITS","%d",&checkhits);
  conf_readkeyword("CHECKSVXACC","%d",&checksvxacc);

  if (!conf_readkeyword("LISTFILE","%s",listfile)) { 
    lf = fopen(listfile,"w");       // User wants an event list
    lflag=1;                        // so let's open the file and set a flag
    printf("An event list will be written to file %s\n",listfile);
    printf("The signal region is defined as %8.5f - %8.5f GeV\n",lmmin,lmmax);
    printf("Minimum ctau for signal events: %f\n",lctau);
  } else {
    lflag=0;
  }

  if (bcode==0) printf("%c[31mWARNING ! Accepting all decays ! %c[0m\n",27,27); 
  printf("\nBitcode=%d  masswin=%f\n",bcode,masswin);
  printf("Objectname = %s\n",oname);
  printf("Chicuts       : %f  %f\n",Chiprob_cut,Chi2d_cut);
  printf("pt cuts       : B %f     K %f\n",Pt_cut, Kpt_cut);
  printf("ectau cut     : %f\n",Ectau_cut);
  printf("min. svx hits : %d\n",nsvxhits);
  printf("check svx hit pattern for < 3 hits : %d\n",checkhits);
  printf("check svxacc  : %d\n",checksvxacc);

  mmin=PDGBmass-masswin; mmax=PDGBmass+masswin;
  TFile *ff = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  if (!ff) {
    ff = new TFile(filename);
  }
  tree = (TTree*)gDirectory->Get("h1000");
  
  // shuffle the TTree into memory just using the branches we are interested in  
  // and just accepting events passing a few pre selection cuts. 
  // (would be nice if we could actually prune branches instead of just
  // making them inactive.
  h1000 *t= new h1000(tree);
  
  t->fChain->SetBranchStatus("*",0);
  t->fChain->SetBranchStatus("Mass" ,1);
  t->fChain->SetBranchStatus("Emass",1);
  t->fChain->SetBranchStatus("Pt",1);
  t->fChain->SetBranchStatus("Kpt",1);
  t->fChain->SetBranchStatus("Chiprob",1);
  t->fChain->SetBranchStatus("Chi2d",1);
  t->fChain->SetBranchStatus("Ctau",1);
  t->fChain->SetBranchStatus("Ectau",1);
  t->fChain->SetBranchStatus("Isolz",1);
  t->fChain->SetBranchStatus("Bitcode",1);
  t->fChain->SetBranchStatus("Run",1);
  t->fChain->SetBranchStatus("Event",1);
  t->fChain->SetBranchStatus("Svx1",1);
  t->fChain->SetBranchStatus("Svx2",1);

  t->fChain->SetBranchStatus("Kchrg",1);

  gROOT->cd();    // back to program scope 
  gROOT->pwd();   // show current directory
  //  cout << gROOT->GetPath() << endl;
  
  // now construct the character string defining the preselection cut:
  // this way we don't have hardcoded cuts in the code :-)  
  sprintf(cutname,"Chiprob>%f&&Chi2d<%f&&Mass>%f&&Mass<%f&&Pt>%f&&Kpt>%f&&abs(Ectau)<%f",Chiprob_cut,Chi2d_cut,mmin,mmax,Pt_cut,Kpt_cut,Ectau_cut);
  newtree =t->fChain->CopyTree(cutname);
  
  newtree->SetBranchAddress("Run",&Run);
  newtree->SetBranchAddress("Event",&Event);
  newtree->SetBranchAddress("Mass",&Mass); 
  newtree->SetBranchAddress("Emass",&Emass); 
  newtree->SetBranchAddress("Pt",&Pt);
  newtree->SetBranchAddress("Kpt",&Kpt);
  newtree->SetBranchAddress("Chiprob",&Chiprob); 
  newtree->SetBranchAddress("Chi2d",&Chi2d); 
  newtree->SetBranchAddress("Ctau",&Ctau);
  newtree->SetBranchAddress("Ectau",&Ectau); 
  newtree->SetBranchAddress("Isolz",&Isolz); 
  newtree->SetBranchAddress("Bitcode",&Bitcode); 
  newtree->SetBranchAddress("Svx1",&Svx1);
  newtree->SetBranchAddress("Svx2",&Svx2);
  
  // Now loop over newtree and fill the vector:
  nbad = 0; lindex=0; 
  std::vector<Measurement>  m;
  m.reserve(newtree->GetEntries());
  cout << " Capacity of vector:  " << m.capacity()<< endl;
  for (i=0; i<newtree->GetEntries();i++) {
    newtree->GetEntry(i);   
    bits = (Int_t)Bitcode;
    if (bits != bcode && bcode!=0) {
      nbad++; flag=0;
      for (j=1; j<=lindex; j++) {
	if (bits==clist[j]) {
	  cbad[j]++;
	  flag=1;
	}
      }
      if (!flag) {
	lindex++;
	clist[lindex]=bits;
	cbad[lindex]=1;
      }
      continue; // Skip event with wrong bitcode
    }
    
    // Check no. of svx hits and hit pattern
    // 1st lepton
    svxacc = ((int)Svx1 & 128)/128;
    l0 = (int)Svx1 & 1;
    l1 = ((int)Svx1 & 2)/2;
    l2 = ((int)Svx1 & 4)/4;
    l3 = ((int)Svx1 & 8)/8;
    svxh1 = l0+l1+l2+l3;   // No. of svx hits

    //printf("%d: %d %d %d %d\n",svxh1,l0,l1,l2,l3);
    if (svxh1 < nsvxhits) continue; // Not enough hits
    if ((svxh1 < 3) && (checkhits)) {
      if (! ((l0 && l1) || (l0 && l2) || (l1 && l2))) continue;
    }
    if (checksvxacc && (! svxacc)) continue;

    // 2nd lepton
    svxacc = ((int)Svx2 & 128)/128;
    l0 = (int)Svx2 & 1;
    l1 = ((int)Svx2 & 2)/2;
    l2 = ((int)Svx2 & 4)/4;
    l3 = ((int)Svx2 & 8)/8;
    svxh2 = l0+l1+l2+l3;   // No. of svx hits

    if (svxh2 < nsvxhits) continue; // Not enough hits
    if ((svxh2 < 3) && (checkhits)) {
      if (! ((l0 && l1) || (l0 && l2) || (l1 && l2))) continue;
    }
    if (checksvxacc && (! svxacc)) continue;

    m.push_back(Measurement(Ctau,Ectau,Mass,Emass,Isolz,Pt,Kpt));
    if (lflag) {  // if lflag is set, we write event to the list file
      // First, check if candidate is in the defined signal region
      if (Mass>=lmmin && Mass<=lmmax && Ctau>=lctau) { 
	sprintf(buffer,"%9d %9d   %9.6f\n",(int)Run,(int)Event,Mass);
	fputs(buffer,lf);
      }
    }
  }
 
  delete newtree;
  BDecay *B;
  B = new BDecay(oname,"- ",m); 

  B->SetCuts(mmin,mmax,Chiprob_cut,Chi2d_cut,Pt_cut,Kpt_cut);
  //B->Setbitcode(bcode);

  TFile *f = new TFile(outfile, "RECREATE"); 
  B->Write(oname);
  f ->ls();
  delete f;
  printf("\nNo. of events with wrong bitcode : %d\n",nbad);
  for (j=1; j<=lindex; j++) {
    printf("Bitcode %6d :   No. of times : %7d\n",clist[j],cbad[j]);
  }
  if (lflag) fclose(lf);  // close event list file
  //app.Run();   
  return 0;
}
