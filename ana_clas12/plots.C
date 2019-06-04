{
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSelector.h"
#include "TCanvas.h"
using namespace std;

TFile *f=new TFile("ana_test.root");
TTree *e=(TTree *)f->Get("el");
TTree *p=(TTree *)f->Get("pDVCS");
TTree *n=(TTree *)f->Get("nDVCS");
TCanvas *c=new TCanvas("plots","plots",600,600);

c->Divide(3,2);
c->cd(1);
e->Draw("El_Theta:El_Phi","","colz");
c->cd(2);
p->Draw("Pr_Theta:Pr_Phi","","colz");
c->cd(3);
n->Draw("N_Theta:N_Phi","","colz");
c->cd(4);
e->Draw("Q2:Xbj","Xbj>0 && Xbj<1.5 && Q2>1","colz");
c->cd(5);
p->Draw("Ph_Theta:Ph_Phi","","colz");

c->Update();
c->Draw();

}
