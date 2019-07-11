{
TFile *f=new TFile("/home/justind/DATA/40.root");
TChain *c=new TChain("hipotree");
c->Add("/home/justind/DATA/40.root");
c->Process("DcoDe.C");

}
