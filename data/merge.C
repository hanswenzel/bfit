{
gSystem.Load("BDecay.so");

TFile *f = new TFile("../data/b0k0s/B0toJPsiK0S_3_4.5_1.00.root");
S1 = (BDecay *)f->Get("B0toJPsiK0S");
S1 -> Size();
delete f;

TFile *f2 = new TFile("../data/b0kstar/B0toJPsiKstar_3_5.0_1.25.root");
S2 = (BDecay *)f2->Get("B0toJPsiKstar");
S2 -> Size();
delete f2;

TFile *out = new TFile("B0_2404.root","RECREATE");
S1->Write("B0toJPsiK0S");
S2->Write("B0toJPsiKstar");
out->ls();
delete out;
}
