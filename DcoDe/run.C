{
//TFile *f=new TFile("40.root");
TChain *c=new TChain("hiporoot");
c->Add("40.root");
c->Process("DCoDe.C");

}
