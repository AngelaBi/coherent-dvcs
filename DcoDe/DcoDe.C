#define DcoDe_cxx
// The class definition in DcoDe.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("DcoDe.C")
// root> T->Process("DcoDe.C","some options")
// root> T->Process("DcoDe.C+")
//


#include "DcoDe.h"
//#include <TH2.h>
//#include <TStyle.h>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "Rtypes.h"
#include "/home/justind/Clas12Tool/Clas12Root/ParticleTree.h"
#include "/home/justind/Clas12Tool/Clas12Root/HipoRootAction.h"
#include "/home/justind/Clas12Tool/Clas12Banks3/clas12reader.h"
using namespace std;

void DcoDe::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   outfile = new TFile("test.root","RECREATE");

   AddBranches();
   firstCall= kTRUE;
   Pmass = 0.938272;
   Nmass = 0.93957;
   Dmass = 1.8756;
   LightSpeed = 29.9792458;
   Ebeam = 10.2;

   tstart = 124.25;

   ElectronBeam.SetXYZT(0,0,Ebeam,Ebeam);
   Target_Vec.SetXYZT(0,0,0,Dmass);
   PTarget_Vec.SetXYZT(0,0,0,Pmass);
   NTarget_Vec.SetXYZT(0,0,0,Nmass);

   firstcall = kFALSE;
}

void DcoDe::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t DcoDe::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);
   nelec=0;
   i_el = -1;
   nprot=0;
   i_pr = -1;
   nneut=0;
   i_n=0;
   nphotEC=0;
   i_phEC = -1;
   nphotFT=0;
   i_phFT = -1;
   nphot=0;
   i_ph=-1;
   npion=0;
   FD_proton=0;
   CD_proton=0;
   FD_neutron=0;
   CD_neutron=0;

   int size = PBANK.Pid.GetSize();

for(int j=0;j<size;j++)
    {
      Float_t pmom[size];
      pmom[j]=TMath::Sqrt((PBANK.Px)[j]*(PBANK.Px)[j]+(PBANK.Py)[j]*(PBANK.Py)[j]+(PBANK.Pz)[j]*(PBANK.Pz)[j]);

 //electron ID
      if((PBANK.Charge)[j]<0 && (PBANK.Pid)[j]==11 && (PBANK.Status)[j]>1000)
        {
          nelec++;
          i_el=j;
          El_Vec.SetPxPyPzE((PBANK.Px)[j],(PBANK.Py)[j],(PBANK.Pz)[j],pmom[j]);
          Q2 = 4*Ebeam*El_Vec.P()*TMath::Power(TMath::Sin(El_Vec.Theta()/2),2.);;
          Float_t nu = Ebeam-El_P;
          W = TMath::Sqrt(Pmass*Pmass + 2*Pmass*nu - Q2);
          El_P = pmom[j];
          Xbj = Q2/(2*Pmass*(Ebeam-El_P));
          El_Theta = El_Vec.Theta()*180./TMath::Pi();
          El_Phi = El_Vec.Phi()*180./TMath::Pi();
          El_vz = (PBANK.Vz)[j];
          El_vx = (PBANK.Vx)[j];
          El_vy = (PBANK.Vy)[j];
        }
 if(nelec>=1)
        {
          //photon ID 
          if(j!=i_el && (PBANK.Charge)[j]==0 && (PBANK.Pid)[j]==22)
            {
              nphot++;
              i_ph=j;
              Ph_Vec.SetPxPyPzE((PBANK.Px)[j],(PBANK.Py)[j],(PBANK.Pz)[j],pmom[j]);
              Ph_P = Ph_Vec.E();
              Ph_Theta = Ph_Vec.Theta()*180./TMath::Pi();
              Ph_Phi = Ph_Vec.Phi()*180./TMath::Pi();
            }
if(j!=i_el && j!=i_ph && (PBANK.Pid)[j]==2212)
            {
              //proton ID
              //nprot++;
              i_pr=j;
              Pr_Vec.SetPxPyPzE((PBANK.Px)[j],(PBANK.Py)[j],(PBANK.Pz)[j],TMath::Sqrt(pmom[j]*pmom[j]+Pmass*Pmass));
              Pr_P = pmom[j];
              Pr_Theta = Pr_Vec.Theta()*180./TMath::Pi();
              Pr_Phi = Pr_Vec.Phi()*180./TMath::Pi();
              Pr_vz = (PBANK.vz)[j];
              if((PBANK.Status)[j]>=2000 && (PBANK.Status)[j]<4000)FD_proton=1;
              else if((PBANK.Status)[j]>=4000)CD_proton=1;
              if(FD_proton==1 || CD_proton==1)nprot++;
            }
          if(j!=i_el && j!=i_ph && TMath::Abs((PBANK.Pid)[j])==211)
            {
              npion++;
            }
 if(j!=i_el && j!=i_ph && j!=i_pr && (PBANK.Charge)[j]==0 && (PBANK.Pid)[j]==2112)
            {
        //nneut++;
              i_n=j;
              N_Vec.SetPxPyPzE((PBANK.Px)[j],(PBANK.Py)[j],(PBANK.Pz)[j],TMath::Sqrt(pmom[j]*pmom[j]+Nmass*Nmass));
              N_P = pmom[j];
              N_Theta = N_Vec.Theta()*180./TMath::Pi();
              N_Phi = N_Vec.Phi()*180./TMath::Pi();
              if((PBANK.Status)[j]>=2000 && (PBANK.Status)[j]<4000)FD_neutron=1;
              else if((PBANK.Status)[j]>=4000)CD_neutron=1;
              if(FD_neutron==1 || CD_neutron==1)nneut++;
            }
        }
    }
if(nelec>=1)el_tree->Fill();
  if(nelec>=1 && nprot==1 && nphot>=1)
    {
      Calc_kine_P();
      pDVCS_tree->Fill();
    }
  if(nelec>=1 && nneut==1 && nphot>=1 && nprot==0 && npion==0)
    {
      Calc_kine_N();
      nDVCS_tree->Fill();
    }
  if (nelec==1 && nprot==1){
      cout<<"hello?"<<endl;
      Calc_kine_P();
      ep_tree->Fill();
}
   return kTRUE;
}

void DcoDe::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void DcoDe::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  outfile->cd();
  el_tree->Write();
  pDVCS_tree->Write();
  nDVCS_tree->Write();
  ep_tree->Write();
  //  MC_tree->Write();
  outfile->Write();
  outfile->Close();

}

void DcoDe::Calc_kine_P()
{
    TVector3 VelectronIn,VelectronOut,VprotonOut,VphotonOut,Vlepto,Vhadro,VhadroPP,Vvirtualphoton;
    TVector3 VX;

    VelectronIn    = ElectronBeam.Vect();
    VelectronOut   = El_Vec.Vect();
    VprotonOut     = Pr_Vec.Vect();
    VphotonOut     = Ph_Vec.Vect();
    Vvirtualphoton = (ElectronBeam-El_Vec).Vect();

    Vlepto         = VelectronIn.Cross(VelectronOut);
    Vhadro         = VprotonOut.Cross(Vvirtualphoton);
    VhadroPP       = VprotonOut.Cross(VphotonOut);

    Phi_Pr = 180./TMath::Pi()*Vlepto.Angle(Vhadro);
    Phi_Ph = 180./TMath::Pi()*Vlepto.Angle(VhadroPP);

    if(Vlepto.Dot(VprotonOut)>0.)  Phi_Pr = 360.-Phi_Pr;
    if(Vlepto.Dot(VphotonOut)<0.)  Phi_Ph = 360.-Phi_Ph;

    t_Pr  = (Pr_Vec-PTarget_Vec).M2();
    t_Ph  = (ElectronBeam-El_Vec-Ph_Vec).M2();

    TLorentzVector BalV = ElectronBeam+Target_Vec-Ph_Vec-El_Vec-Pr_Vec;
    mm2_epg = (ElectronBeam+Target_Vec-Pr_Vec-El_Vec-Ph_Vec).M2();
    Xbal  = BalV.X();
    Ybal  = BalV.Y();
    Zbal  = BalV.Z();
    Ebal  = BalV.E();
}

void DcoDe::Calc_kine_N()
{

    TVector3 VelectronIn,VelectronOut,VneutronOut,VphotonOut,Vlepto,Vhadro,VhadroPP,Vvirtualphoton;
    TVector3 VX;

    VelectronIn    = ElectronBeam.Vect();
    VelectronOut   = El_Vec.Vect();
    VneutronOut     = N_Vec.Vect();
    VphotonOut     = Ph_Vec.Vect();
    Vvirtualphoton = (ElectronBeam-El_Vec).Vect();

    Vlepto         = VelectronIn.Cross(VelectronOut);
    Vhadro         = VneutronOut.Cross(Vvirtualphoton);
    VhadroPP       = VneutronOut.Cross(VphotonOut);

    Phi_N = 180./TMath::Pi()*Vlepto.Angle(Vhadro);
    Phi_Ph = 180./TMath::Pi()*Vlepto.Angle(VhadroPP);

    if(Vlepto.Dot(VneutronOut)>0.)  Phi_N = 360.-Phi_N;
    if(Vlepto.Dot(VphotonOut)<0.)  Phi_Ph = 360.-Phi_Ph;

    t_N  = (N_Vec-NTarget_Vec).M2();
    t_Ph  = (ElectronBeam-El_Vec-Ph_Vec).M2();

    TLorentzVector BalV = ElectronBeam+Target_Vec-Ph_Vec-El_Vec-N_Vec;

    mm2_eng = (ElectronBeam+Target_Vec-N_Vec-El_Vec-Ph_Vec).M2();
    Xbal  = BalV.X();
    Ybal  = BalV.Y();
    Zbal  = BalV.Z();
    Ebal  = BalV.E();
}

void DcoDe::AddBranches()
{
   el_tree = (TTree*) new TTree("el","el");
   el_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
   el_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
   //   el_tree->Branch("bpola",&bpola,"bpola/F");
   el_tree->Branch("helicity",&helicity,"helicity/I");
   el_tree->Branch("TarPol",&TarPol,"TarPol/F");
   el_tree->Branch("Q2",&Q2,"Q2/F");
   el_tree->Branch("Xbj", &Xbj, "Xbj/F");
   el_tree->Branch("W", &W, "W/F");
   el_tree->Branch("El_P",&El_P,"El_P/F");
   el_tree->Branch("El_Theta",&El_Theta,"El_Theta/F");
   el_tree->Branch("El_Phi",&El_Phi,"El_Phi/F");
   el_tree->Branch("El_vx",&El_vx,"El_vx/F");
   el_tree->Branch("El_vy",&El_vy,"El_vy/F");
   el_tree->Branch("El_vz",&El_vz,"El_vz/F");
   el_tree->Branch("El_vx_MC",&El_vx_MC,"El_vx_MC/F");
   el_tree->Branch("El_vy_MC",&El_vy_MC,"El_vy_MC/F");
   el_tree->Branch("El_vz_MC",&El_vz_MC,"El_vz_MC/F");
   el_tree->Branch("El_P_MC",&El_P_MC,"El_P_MC/F");
   el_tree->Branch("El_Theta_MC",&El_Theta_MC,"El_Theta_MC/F");
   el_tree->Branch("El_Phi_MC",&El_Phi_MC,"El_Phi_MC/F");
   el_tree->Branch("Q2_MC",&Q2_MC,"Q2_MC/F");
   el_tree->Branch("Xbj_MC", &Xbj_MC, "Xbj_MC/F");
   el_tree->Branch("W_MC", &W_MC, "W_MC/F");

   pDVCS_tree = (TTree*) new TTree("pDVCS","pDVCS");
   pDVCS_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
   pDVCS_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
   pDVCS_tree->Branch("helicity",&helicity,"helicity/I");
   pDVCS_tree->Branch("TarPol",&TarPol,"TarPol/F");
   pDVCS_tree->Branch("t_Pr",&t_Pr,"t_Pr/F");
   pDVCS_tree->Branch("t_Ph",&t_Ph,"t_Ph/F");
   pDVCS_tree->Branch("Phi_Pr",&Phi_Pr,"Phi_Pr/F");
   pDVCS_tree->Branch("Phi_Ph",&Phi_Ph,"Phi_Ph/F");
   pDVCS_tree->Branch("Q2",&Q2,"Q2/F");
   pDVCS_tree->Branch("Xbj", &Xbj, "Xbj/F");
   pDVCS_tree->Branch("W", &W, "W/F");
   pDVCS_tree->Branch("Xbal",&Xbal,"Xbal/F");
   pDVCS_tree->Branch("Ybal",&Ybal,"Ybal/F");
   pDVCS_tree->Branch("Zbal",&Zbal,"Zbal/F");
   pDVCS_tree->Branch("Ebal",&Ebal,"Ebal/F");
   pDVCS_tree->Branch("nelec",&nelec,"nelec/I");
   pDVCS_tree->Branch("El_P",&El_P,"El_P/F");
   pDVCS_tree->Branch("El_Theta",&El_Theta,"El_Theta/F");
   pDVCS_tree->Branch("El_Phi",&El_Phi,"El_Phi/F");
   pDVCS_tree->Branch("El_vz",&El_vz,"El_vz/F");
   pDVCS_tree->Branch("El_vx",&El_vx,"El_vx/F");
   pDVCS_tree->Branch("El_vy",&El_vy,"El_vy/F");
   pDVCS_tree->Branch("nprot",&nprot,"nprot/I");
   pDVCS_tree->Branch("Pr_P",&Pr_P,"Pr_P/F");
   pDVCS_tree->Branch("Pr_Theta",&Pr_Theta,"Pr_Theta/F");
   pDVCS_tree->Branch("Pr_Phi",&Pr_Phi,"Pr_Phi/F");
   pDVCS_tree->Branch("Pr_vz",&Pr_vz,"Pr_vz/F");
   pDVCS_tree->Branch("FD_proton",&FD_proton,"FD_proton/I");
   pDVCS_tree->Branch("CD_proton",&CD_proton,"CD_proton/I");
   pDVCS_tree->Branch("nphotEC",&nphotEC,"nphotEC/I");
   pDVCS_tree->Branch("Ph_EC_P",Ph_EC_P,"Ph_EC_P[nphotEC]/F");
   pDVCS_tree->Branch("Ph_EC_Theta",Ph_EC_Theta,"Ph_EC_Theta[nphotEC]/F");
   pDVCS_tree->Branch("Ph_EC_Phi",Ph_EC_Phi,"Ph_EC_Phi[nphotEC]/F");
   pDVCS_tree->Branch("nphotFT",&nphotFT,"nphotFT/I");
   pDVCS_tree->Branch("Ph_FT_P",Ph_FT_P,"Ph_FT_P[nphotFT]/F");
   pDVCS_tree->Branch("Ph_FT_Theta",Ph_FT_Theta,"Ph_FT_Theta[nphotFT]/F");
   pDVCS_tree->Branch("Ph_FT_Phi",Ph_FT_Phi,"Ph_FT_Phi[nphotFT]/F");
   pDVCS_tree->Branch("Ph_det",&Ph_det,"Ph_det/I");
   pDVCS_tree->Branch("nphot",&nphot,"nphot/I");
   pDVCS_tree->Branch("Ph_P",&Ph_P,"Ph_P/F");
   pDVCS_tree->Branch("Ph_Theta",&Ph_Theta,"Ph_Theta/F");
   pDVCS_tree->Branch("Ph_Phi",&Ph_Phi,"Ph_Phi/F");
   pDVCS_tree->Branch("mm2_epg",&mm2_epg,"mm2_epg/F");
   pDVCS_tree->Branch("El_P_MC",&El_P_MC,"El_P_MC/F");
   pDVCS_tree->Branch("El_Theta_MC",&El_Theta_MC,"El_Theta_MC/F");
   pDVCS_tree->Branch("El_Phi_MC",&El_Phi_MC,"El_Phi_MC/F");
   pDVCS_tree->Branch("Pr_P_MC",&Pr_P_MC,"Pr_P_MC/F");
   pDVCS_tree->Branch("Pr_Theta_MC",&Pr_Theta_MC,"Pr_Theta_MC/F");
   pDVCS_tree->Branch("Pr_Phi_MC",&Pr_Phi_MC,"Pr_Phi_MC/F");
   pDVCS_tree->Branch("Ph_P_MC",&Ph_P_MC,"Ph_P_MC/F");
   pDVCS_tree->Branch("Ph_Theta_MC",&Ph_Theta_MC,"Ph_Theta_MC/F");
   pDVCS_tree->Branch("Ph_Phi_MC",&Ph_Phi_MC,"Ph_Phi_MC/F");
   pDVCS_tree->Branch("Q2_MC",&Q2_MC,"Q2_MC/F");
   pDVCS_tree->Branch("Xbj_MC", &Xbj_MC, "Xbj_MC/F");
   pDVCS_tree->Branch("W_MC", &W_MC, "W_MC/F");
   pDVCS_tree->Branch("t_Pr_MC",&t_Pr_MC,"t_Pr_MC/F");
   pDVCS_tree->Branch("t_N_MC",&t_N_MC,"t_N_MC/F");
   pDVCS_tree->Branch("t_Ph_MC",&t_Ph_MC,"t_Ph_MC/F");
   pDVCS_tree->Branch("Phi_MC",&Phi_MC,"Phi_MC/F");
   pDVCS_tree->Branch("n_spec",&n_spec,"n_spec/I");
   pDVCS_tree->Branch("p_spec",&p_spec,"p_spec/I");

   nDVCS_tree = (TTree*) new TTree("nDVCS","nDVCS");
   nDVCS_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
   nDVCS_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
   nDVCS_tree->Branch("helicity",&helicity,"helicity/I");
   nDVCS_tree->Branch("TarPol",&TarPol,"TarPol/F");
   nDVCS_tree->Branch("Q2",&Q2,"Q2/F");
   nDVCS_tree->Branch("Xbj", &Xbj, "Xbj/F");
   nDVCS_tree->Branch("W", &W, "W/F");
   nDVCS_tree->Branch("t_Ph",&t_Ph,"t_Ph/F");
   nDVCS_tree->Branch("t_N",&t_N,"t_N/F");
   nDVCS_tree->Branch("Phi_Ph",&Phi_Ph,"Phi_Ph/F");
   nDVCS_tree->Branch("Phi_N",&Phi_N,"Phi_N/F");
   nDVCS_tree->Branch("Xbal",&Xbal,"Xbal/F");
   nDVCS_tree->Branch("Ybal",&Ybal,"Ybal/F");
   nDVCS_tree->Branch("Zbal",&Zbal,"Zbal/F");
   nDVCS_tree->Branch("Ebal",&Ebal,"Ebal/F");
   nDVCS_tree->Branch("nelec",&nelec,"nelec/I");
   nDVCS_tree->Branch("El_P",&El_P,"El_P/F");
   nDVCS_tree->Branch("El_Theta",&El_Theta,"El_Theta/F");
   nDVCS_tree->Branch("El_Phi",&El_Phi,"El_Phi/F");
   nDVCS_tree->Branch("El_vz",&El_vz,"El_vz/F");
   nDVCS_tree->Branch("El_vx",&El_vx,"El_vx/F");
   nDVCS_tree->Branch("El_vy",&El_vy,"El_vy/F");
   nDVCS_tree->Branch("nprot",&nprot,"nprot/I");
   nDVCS_tree->Branch("FD_neutron",&FD_neutron,"FD_neutron/I");
   nDVCS_tree->Branch("CD_neutron",&CD_neutron,"CD_neutron/I");
   nDVCS_tree->Branch("N_P",&N_P,"N_P/F");
   nDVCS_tree->Branch("N_Theta",&N_Theta,"N_Theta/F");
   nDVCS_tree->Branch("N_Phi",&N_Phi,"N_Phi/F");
   nDVCS_tree->Branch("nphotEC",&nphotEC,"nphotEC/I");
   nDVCS_tree->Branch("Ph_EC_P",Ph_EC_P,"Ph_EC_P[nphotEC]/F");
   nDVCS_tree->Branch("Ph_EC_Theta",Ph_EC_Theta,"Ph_EC_Theta[nphotEC]/F");
   nDVCS_tree->Branch("Ph_EC_Phi",Ph_EC_Phi,"Ph_EC_Phi[nphotEC]/F");
   nDVCS_tree->Branch("nphotFT",&nphotFT,"nphotFT/I");
   nDVCS_tree->Branch("Ph_FT_P",Ph_FT_P,"Ph_FT_P[nphotFT]/F");
   nDVCS_tree->Branch("Ph_FT_Theta",Ph_FT_Theta,"Ph_FT_Theta[nphotFT]/F");
   nDVCS_tree->Branch("Ph_FT_Phi",Ph_FT_Phi,"Ph_FT_Phi[nphotFT]/F");
   nDVCS_tree->Branch("Ph_det",&Ph_det,"Ph_det/I");
   nDVCS_tree->Branch("nphot",&nphot,"nphot/I");
   nDVCS_tree->Branch("Ph_P",&Ph_P,"Ph_P/F");
   nDVCS_tree->Branch("Ph_Theta",&Ph_Theta,"Ph_Theta/F");
   nDVCS_tree->Branch("Ph_Phi",&Ph_Phi,"Ph_Phi/F");
   nDVCS_tree->Branch("mm2_eng",&mm2_eng,"mm2_eng/F");
   nDVCS_tree->Branch("El_P_MC",&El_P_MC,"El_P_MC/F");
   nDVCS_tree->Branch("El_Theta_MC",&El_Theta_MC,"El_Theta_MC/F");
   nDVCS_tree->Branch("El_Phi_MC",&El_Phi_MC,"El_Phi_MC/F");
   nDVCS_tree->Branch("Pr_P_MC",&Pr_P_MC,"Pr_P_MC/F");
   nDVCS_tree->Branch("Pr_Theta_MC",&Pr_Theta_MC,"Pr_Theta_MC/F");
   nDVCS_tree->Branch("Pr_Phi_MC",&Pr_Phi_MC,"Pr_Phi_MC/F");
   nDVCS_tree->Branch("N_P_MC",&N_P_MC,"N_P_MC/F");
   nDVCS_tree->Branch("N_Theta_MC",&N_Theta_MC,"N_Theta_MC/F");
   nDVCS_tree->Branch("N_Phi_MC",&N_Phi_MC,"N_Phi_MC/F");
   nDVCS_tree->Branch("Ph_P_MC",&Ph_P_MC,"Ph_P_MC/F");
   nDVCS_tree->Branch("Ph_Theta_MC",&Ph_Theta_MC,"Ph_Theta_MC/F");
   nDVCS_tree->Branch("Ph_Phi_MC",&Ph_Phi_MC,"Ph_Phi_MC/F");
   nDVCS_tree->Branch("Q2_MC",&Q2_MC,"Q2_MC/F");
   nDVCS_tree->Branch("Xbj_MC", &Xbj_MC, "Xbj_MC/F");
   nDVCS_tree->Branch("W_MC", &W_MC, "W_MC/F");
   nDVCS_tree->Branch("t_Pr_MC",&t_Pr_MC,"t_Pr_MC/F");
   nDVCS_tree->Branch("t_N_MC",&t_N_MC,"t_N_MC/F");
   nDVCS_tree->Branch("t_Ph_MC",&t_Ph_MC,"t_Ph_MC/F");
   nDVCS_tree->Branch("Phi_MC",&Phi_MC,"Phi_MC/F");
   nDVCS_tree->Branch("n_spec",&n_spec,"n_spec/I");
   nDVCS_tree->Branch("p_spec",&p_spec,"p_spec/I");

   ep_tree = (TTree*) new TTree("ep","ep");
   ep_tree->Branch("El_P",&El_P,"El_P/F");
   ep_tree->Branch("El_Theta",&El_Theta,"El_Theta/F");
   ep_tree->Branch("El_vz",&El_vz,"El_vz/F");
   ep_tree->Branch("El_vx",&El_vx,"El_vx/F");
   ep_tree->Branch("El_vy",&El_vy,"El_vy/F");
   ep_tree->Branch("El_Phi",&El_Phi,"El_Phi/F");
   ep_tree->Branch("Pr_P",&Pr_P,"Pr_P/F");
   ep_tree->Branch("Pr_Theta",&Pr_Theta,"Pr_Theta/F");
   ep_tree->Branch("Pr_Phi",&Pr_Phi,"Pr_Phi/F");
   ep_tree->Branch("Q2",&Q2,"Q2/F");
   ep_tree->Branch("Xbj", &Xbj, "Xbj/F");
   ep_tree->Branch("W", &W, "W/F");
}

}
