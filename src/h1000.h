//////////////////////////////////////////////////////////
//   This class has been automatically generated 
//     (Tue Jan 30 15:04:10 2001 by ROOT version3.00/00)
//   from TTree h1000/Bd_new
//   found on file: 1000_b0_diele_corr_nocon_3001.root
//////////////////////////////////////////////////////////


#ifndef h1000_h
#define h1000_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class h1000 {
   public :
   TTree          *fChain;   //pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //current Tree number in a TChain
//Declaration of leaves types
   Float_t         Run;
   Float_t         Event;
   Float_t         Lum;
   Float_t         Bitcode;
   Float_t         Beamx;
   Float_t         Beamy;
   Float_t         Vertx;
   Float_t         Verty;
   Float_t         Vertz;
   Float_t         Chiprob;
   Float_t         Chi2d;
   Float_t         Mass;
   Float_t         Emass;
   Float_t         Normmass;
   Float_t         Pt;
   Float_t         Kmass;
   Float_t         Kpt;
   Float_t         Kchrg;
   Float_t         Kvertx;
   Float_t         Kverty;
   Float_t         Kvertz;
   Float_t         Lxy;
   Float_t         Elxy;
   Float_t         Ctau;
   Float_t         Ectau;
   Float_t         Isol;
   Float_t         Isolz;
   Float_t         Svx1;
   Float_t         Svx2;

//List of branches
   TBranch        *b_Run;
   TBranch        *b_Event;
   TBranch        *b_Lum;
   TBranch        *b_Bitcode;
   TBranch        *b_Beamx;
   TBranch        *b_Beamy;
   TBranch        *b_Vertx;
   TBranch        *b_Verty;
   TBranch        *b_Vertz;
   TBranch        *b_Chiprob;
   TBranch        *b_Chi2d;
   TBranch        *b_Mass;
   TBranch        *b_Emass;
   TBranch        *b_Normmass;
   TBranch        *b_Pt;
   TBranch        *b_Kmass;
   TBranch        *b_Kpt;
   TBranch        *b_Kchrg;
   TBranch        *b_Kvertx;
   TBranch        *b_Kverty;
   TBranch        *b_Kvertz;
   TBranch        *b_Lxy;
   TBranch        *b_Elxy;
   TBranch        *b_Ctau;
   TBranch        *b_Ectau;
   TBranch        *b_Isol;
   TBranch        *b_Isolz;
   TBranch        *b_Svx1;
   TBranch        *b_Svx2;

   h1000(TTree *tree=0);
   ~h1000();
   Int_t  Cut(Int_t entry);
   Int_t  GetEntry(Int_t entry);
   Int_t  LoadTree(Int_t entry);
   void   Init(TTree *tree);
   void   Loop();
   Bool_t Notify();
   void   Show(Int_t entry = -1);
};

#endif

#ifdef h1000_cxx
h1000::h1000(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("default");
      if (!f) {
         f = new TFile("default");
      }
      tree = (TTree*)gDirectory->Get("h1000");

   }
   Init(tree);
}

h1000::~h1000()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t h1000::GetEntry(Int_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t h1000::LoadTree(Int_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Int_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void h1000::Init(TTree *tree)
{
//   Set branch addresses
   if (tree == 0) return;
   fChain    = tree;
   fCurrent = -1;

   fChain->SetBranchAddress("Run",&Run);
   fChain->SetBranchAddress("Event",&Event);
   fChain->SetBranchAddress("Lum",&Lum);
   fChain->SetBranchAddress("Bitcode",&Bitcode);
   fChain->SetBranchAddress("Beamx",&Beamx);
   fChain->SetBranchAddress("Beamy",&Beamy);
   fChain->SetBranchAddress("Vertx",&Vertx);
   fChain->SetBranchAddress("Verty",&Verty);
   fChain->SetBranchAddress("Vertz",&Vertz);
   fChain->SetBranchAddress("Chiprob",&Chiprob);
   fChain->SetBranchAddress("Chi2d",&Chi2d);
   fChain->SetBranchAddress("Mass",&Mass);
   fChain->SetBranchAddress("Emass",&Emass);
   fChain->SetBranchAddress("Normmass",&Normmass);
   fChain->SetBranchAddress("Pt",&Pt);
   fChain->SetBranchAddress("Kmass",&Kmass);
   fChain->SetBranchAddress("Kpt",&Kpt);
   fChain->SetBranchAddress("Kchrg",&Kchrg);
   fChain->SetBranchAddress("Kvertx",&Kvertx);
   fChain->SetBranchAddress("Kverty",&Kverty);
   fChain->SetBranchAddress("Kvertz",&Kvertz);
   fChain->SetBranchAddress("Lxy",&Lxy);
   fChain->SetBranchAddress("Elxy",&Elxy);
   fChain->SetBranchAddress("Ctau",&Ctau);
   fChain->SetBranchAddress("Ectau",&Ectau);
   fChain->SetBranchAddress("Isol",&Isol);
   fChain->SetBranchAddress("Isolz",&Isolz);
   fChain->SetBranchAddress("Svx1",&Svx1);
   fChain->SetBranchAddress("Svx2",&Svx2);
   Notify();
}

Bool_t h1000::Notify()
{
//   called when loading a new file
//   get branch pointers
   b_Run = fChain->GetBranch("Run");
   b_Event = fChain->GetBranch("Event");
   b_Lum = fChain->GetBranch("Lum");
   b_Bitcode = fChain->GetBranch("Bitcode");
   b_Beamx = fChain->GetBranch("Beamx");
   b_Beamy = fChain->GetBranch("Beamy");
   b_Vertx = fChain->GetBranch("Vertx");
   b_Verty = fChain->GetBranch("Verty");
   b_Vertz = fChain->GetBranch("Vertz");
   b_Chiprob = fChain->GetBranch("Chiprob");
   b_Chi2d = fChain->GetBranch("Chi2d");
   b_Mass = fChain->GetBranch("Mass");
   b_Emass = fChain->GetBranch("Emass");
   b_Normmass = fChain->GetBranch("Normmass");
   b_Pt = fChain->GetBranch("Pt");
   b_Kmass = fChain->GetBranch("Kmass");
   b_Kpt = fChain->GetBranch("Kpt");
   b_Kchrg = fChain->GetBranch("Kchrg");
   b_Kvertx = fChain->GetBranch("Kvertx");
   b_Kverty = fChain->GetBranch("Kverty");
   b_Kvertz = fChain->GetBranch("Kvertz");
   b_Lxy = fChain->GetBranch("Lxy");
   b_Elxy = fChain->GetBranch("Elxy");
   b_Ctau = fChain->GetBranch("Ctau");
   b_Ectau = fChain->GetBranch("Ectau");
   b_Isol = fChain->GetBranch("Isol");
   b_Isolz = fChain->GetBranch("Isolz");
   b_Svx1 = fChain->GetBranch("Svx1");
   b_Svx2 = fChain->GetBranch("Svx2");
   return kTRUE;
}

void h1000::Show(Int_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t h1000::Cut(Int_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef h1000_cxx

