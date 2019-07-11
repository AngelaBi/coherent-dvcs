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
#include <TH2.h>
#include <TStyle.h>
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TSelector.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "Rtypes.h"

using namespace std;

void DcoDe::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   outfile = new TFile("DcoDe_test.root","RECREATE");
   AddBranches();

   ElectronBeam.SetXYZT(0,0,EBEAM,EBEAM);
   Target_Vec.SetXYZT(0,0,0,MDEUT);
   PTarget_Vec.SetXYZT(0,0,0,MPROT);
   NTarget_Vec.SetXYZT(0,0,0,MNEUT);

}

void DcoDe::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}
//defining functions for pDVCS and nDVCS
void DcoDe::Calc_kine_P()
{
    TVector3 VelectronIn,VelectronOut,VprotonOut,VphotonOut,Vlepto,Vhadro,VhadroPP,Vvirtualphoton;
    TVector3 VX;

    VelectronIn    = ElectronBeam.Vect();
    VelectronOut   = Elec_Vec.Vect();
    VprotonOut     = Prot_Vec.Vect();
    VphotonOut     = Ph_Vec.Vect();
    Vvirtualphoton = (ElectronBeam-Elec_Vec).Vect();

    Vlepto         = VelectronIn.Cross(VelectronOut);
    Vhadro         = VprotonOut.Cross(Vvirtualphoton);
    VhadroPP       = VprotonOut.Cross(VphotonOut);

    Phi_Prot = 180./TMath::Pi()*Vlepto.Angle(Vhadro);
    Phi_Ph = 180./TMath::Pi()*Vlepto.Angle(VhadroPP);

    if(Vlepto.Dot(VprotonOut)>0.)  Phi_Prot = 360.-Phi_Prot;
    if(Vlepto.Dot(VphotonOut)<0.)  Phi_Ph = 360.-Phi_Ph;

    t_Prot  = (Prot_Vec-PTarget_Vec).M2();
    t_Ph  = (ElectronBeam-Elec_Vec-Ph_Vec).M2();

    TLorentzVector BalV = ElectronBeam+Target_Vec-Ph_Vec-Elec_Vec-Prot_Vec;
    mm2_epg = (ElectronBeam+Target_Vec-Prot_Vec-Elec_Vec-Ph_Vec).M2();
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
    VelectronOut   = Elec_Vec.Vect();
    VneutronOut     = Neut_Vec.Vect();
    VphotonOut     = Ph_Vec.Vect();
    Vvirtualphoton = (ElectronBeam-Elec_Vec).Vect();

    Vlepto         = VelectronIn.Cross(VelectronOut);
    Vhadro         = VneutronOut.Cross(Vvirtualphoton);
    VhadroPP       = VneutronOut.Cross(VphotonOut);

    Phi_Neut = 180./TMath::Pi()*Vlepto.Angle(Vhadro);
    Phi_Ph = 180./TMath::Pi()*Vlepto.Angle(VhadroPP);

    if(Vlepto.Dot(VneutronOut)>0.)  Phi_Neut = 360.-Phi_Neut;
    if(Vlepto.Dot(VphotonOut)<0.)  Phi_Ph = 360.-Phi_Ph;

    t_Neut  = (Neut_Vec-NTarget_Vec).M2();
    t_Ph  = (ElectronBeam-Elec_Vec-Ph_Vec).M2();

    TLorentzVector BalV = ElectronBeam+Target_Vec-Ph_Vec-Elec_Vec-Neut_Vec;

    mm2_eng = (ElectronBeam+Target_Vec-Neut_Vec-Elec_Vec-Ph_Vec).M2();
    Xbal  = BalV.X();
    Ybal  = BalV.Y();
    Zbal  = BalV.Z();
    Ebal  = BalV.E();
}

void DcoDe::AddBranches(){
 el_tree = (TTree*) new TTree("el","el");
 el_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
 //el_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
 //el_tree->Branch("helicity",&helicity,"helicity/I");
 //el_tree->Branch("TarPol",&TarPol,"TarPol/F");
 el_tree->Branch("Q2",&Q2,"Q2/F");
 el_tree->Branch("Xbj", &Xbj, "Xbj/F");
 el_tree->Branch("W", &W, "W/F");
 el_tree->Branch("Elec_P",&Elec_P,"Elec_P/F");
 el_tree->Branch("Elec_Theta",&Elec_Theta,"Elec_Theta/F");
 el_tree->Branch("Elec_Phi",&Elec_Phi,"Elec_Phi/F");
 el_tree->Branch("Elec_Vx",&Elec_Vx,"Elec_Vx/F");
 el_tree->Branch("Elec_Vy",&Elec_Vy,"Elec_Vy/F");
 el_tree->Branch("Elec_Vz",&Elec_Vz,"Elec_Vz/F");

 pDVCS_tree = (TTree*) new TTree("pDVCS","pDVCS");
 pDVCS_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
 // pDVCS_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
 // pDVCS_tree->Branch("helicity",&helicity,"helicity/I");
 // pDVCS_tree->Branch("TarPol",&TarPol,"TarPol/F");
 pDVCS_tree->Branch("t_Prot",&t_Prot,"t_Prot/F");
 pDVCS_tree->Branch("t_Ph",&t_Ph,"t_Ph/F");
 pDVCS_tree->Branch("Phi_Prot",&Phi_Prot,"Phi_Prot/F");
 pDVCS_tree->Branch("Phi_Ph",&Phi_Ph,"Phi_Ph/F");
 pDVCS_tree->Branch("Q2",&Q2,"Q2/F");
 pDVCS_tree->Branch("Xbj", &Xbj, "Xbj/F");
 pDVCS_tree->Branch("W", &W, "W/F");
 pDVCS_tree->Branch("Xbal",&Xbal,"Xbal/F");
 pDVCS_tree->Branch("Ybal",&Ybal,"Ybal/F");
 pDVCS_tree->Branch("Zbal",&Zbal,"Zbal/F");
 pDVCS_tree->Branch("Ebal",&Ebal,"Ebal/F");
 pDVCS_tree->Branch("nElec",&nElec,"nElec/I");
 pDVCS_tree->Branch("Elec_P",&Elec_P,"Elec_P/F");
 pDVCS_tree->Branch("Elec_Theta",&Elec_Theta,"Elec_Theta/F");
 pDVCS_tree->Branch("Elec_Phi",&Elec_Phi,"Elec_Phi/F");
 pDVCS_tree->Branch("Elec_Vz",&Elec_Vz,"Elec_Vz/F");
 pDVCS_tree->Branch("Elec_Vx",&Elec_Vx,"Elec_Vx/F");
 pDVCS_tree->Branch("Elec_Vy",&Elec_Vy,"Elec_Vy/F");
 pDVCS_tree->Branch("nProt",&nProt,"nProt/I");
 pDVCS_tree->Branch("Prot_P",&Prot_P,"Prot_P/F");
 pDVCS_tree->Branch("Prot_Theta",&Prot_Theta,"Prot_Theta/F");
 pDVCS_tree->Branch("Prot_Phi",&Prot_Phi,"Prot_Phi/F");
 pDVCS_tree->Branch("Prot_Vz",&Prot_Vz,"Prot_Vz/F");
 pDVCS_tree->Branch("FD_proton",&FD_Prot,"FD_proton/I");
 pDVCS_tree->Branch("CD_proton",&CD_Prot,"CD_proton/I");
 pDVCS_tree->Branch("nphotElec",&nphotElec,"nphotElec/I");
 //pDVCS_tree->Branch("Ph_EC_P",Ph_EC_P,"Ph_EC_P[nphotEC]/F");
 //pDVCS_tree->Branch("Ph_EC_Theta",Ph_EC_Theta,"Ph_EC_Theta[nphotEC]/F");
 //pDVCS_tree->Branch("Ph_EC_Phi",Ph_EC_Phi,"Ph_EC_Phi[nphotEC]/F");
 //pDVCS_tree->Branch("nphotFT",&nphotFT,"nphotFT/I");
 //pDVCS_tree->Branch("Ph_FT_P",Ph_FT_P,"Ph_FT_P[nphotFT]/F");
 //pDVCS_tree->Branch("Ph_FT_Theta",Ph_FT_Theta,"Ph_FT_Theta[nphotFT]/F");
 //pDVCS_tree->Branch("Ph_FT_Phi",Ph_FT_Phi,"Ph_FT_Phi[nphotFT]/F");
 //pDVCS_tree->Branch("Ph_det",&Ph_det,"Ph_det/I");
 pDVCS_tree->Branch("nPh",&nPh,"nPh/I");
 pDVCS_tree->Branch("Ph_P",&Ph_P,"Ph_P/F");
 pDVCS_tree->Branch("Ph_Theta",&Ph_Theta,"Ph_Theta/F");
 pDVCS_tree->Branch("Ph_Phi",&Ph_Phi,"Ph_Phi/F");
 pDVCS_tree->Branch("mm2_epg",&mm2_epg,"mm2_epg/F");

 nDVCS_tree = (TTree*) new TTree("nDVCS","nDVCS");
 nDVCS_tree->Branch("RunNumber",&RunNumber,"RunNumber/I");
 //nDVCS_tree->Branch("pol_sign",&pol_sign,"pol_sign/I");
 //nDVCS_tree->Branch("helicity",&helicity,"helicity/I");
 //nDVCS_tree->Branch("TarPol",&TarPol,"TarPol/F");
 nDVCS_tree->Branch("Q2",&Q2,"Q2/F");
 nDVCS_tree->Branch("Xbj", &Xbj, "Xbj/F");
 nDVCS_tree->Branch("W", &W, "W/F");
 nDVCS_tree->Branch("t_Ph",&t_Ph,"t_Ph/F");
 nDVCS_tree->Branch("t_Neut",&t_Neut,"t_Neut/F");
 nDVCS_tree->Branch("Phi_Ph",&Phi_Ph,"Phi_Ph/F");
 nDVCS_tree->Branch("Phi_Neut",&Phi_Neut,"Phi_Neut/F");
 nDVCS_tree->Branch("Xbal",&Xbal,"Xbal/F");
 nDVCS_tree->Branch("Ybal",&Ybal,"Ybal/F");
 nDVCS_tree->Branch("Zbal",&Zbal,"Zbal/F");
 nDVCS_tree->Branch("Ebal",&Ebal,"Ebal/F");
 nDVCS_tree->Branch("nElec",&nElec,"nElec/I");
 nDVCS_tree->Branch("Elec_P",&Elec_P,"Elec_P/F");
 nDVCS_tree->Branch("Elec_Theta",&Elec_Theta,"Elec_Theta/F");
 nDVCS_tree->Branch("Elec_Phi",&Elec_Phi,"Elec_Phi/F");
 nDVCS_tree->Branch("Elec_Vz",&Elec_Vz,"Elec_Vz/F");
 nDVCS_tree->Branch("Elec_Vx",&Elec_Vx,"Elec_Vx/F");
 nDVCS_tree->Branch("Elec_Vy",&Elec_Vy,"Elec_Vy/F");
 nDVCS_tree->Branch("nProt",&nProt,"nProt/I");
 nDVCS_tree->Branch("FD_neutron",&FD_Neut,"FD_neutron/I");
 nDVCS_tree->Branch("CD_neutron",&CD_Neut,"CD_neutron/I");
 nDVCS_tree->Branch("Neut_P",&Neut_P,"Neut_P/F");
 nDVCS_tree->Branch("Neut_Theta",&Neut_Theta,"Neut_Theta/F");
 nDVCS_tree->Branch("Neut_Phi",&Neut_Phi,"Neut_Phi/F");
 nDVCS_tree->Branch("nphotElec",&nphotElec,"nphotElec/I");
 //nDVCS_tree->Branch("Ph_EC_P",Ph_EC_P,"Ph_EC_P[nphotEC]/F");
 //nDVCS_tree->Branch("Ph_EC_Theta",Ph_EC_Theta,"Ph_EC_Theta[nphotEC]/F");
 //nDVCS_tree->Branch("Ph_EC_Phi",Ph_EC_Phi,"Ph_EC_Phi[nphotEC]/F");
 //nDVCS_tree->Branch("nphotFT",&nphotFT,"nphotFT/I");
 //nDVCS_tree->Branch("Ph_FT_P",Ph_FT_P,"Ph_FT_P[nphotFT]/F");
 //nDVCS_tree->Branch("Ph_FT_Theta",Ph_FT_Theta,"Ph_FT_Theta[nphotFT]/F");
 //nDVCS_tree->Branch("Ph_FT_Phi",Ph_FT_Phi,"Ph_FT_Phi[nphotFT]/F");
 //nDVCS_tree->Branch("Ph_det",&Ph_det,"Ph_det/I");
 nDVCS_tree->Branch("nPh",&nPh,"nPh/I");
 nDVCS_tree->Branch("Ph_P",&Ph_P,"Ph_P/F");
 nDVCS_tree->Branch("Ph_Theta",&Ph_Theta,"Ph_Theta/F");
 nDVCS_tree->Branch("Ph_Phi",&Ph_Phi,"Ph_Phi/F");
 nDVCS_tree->Branch("mm2_eng",&mm2_eng,"mm2_eng/F");

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
   int size = CLAS12Data_PBANK_Pid.GetSize();



   for(int i=0;i<size;i++){

Float_t pmom[size];
pmom[i]=TMath::Sqrt((CLAS12Data_PBANK_Px)[i]*(CLAS12Data_PBANK_Px)[i]+(CLAS12Data_PBANK_Py)[i]*(CLAS12Data_PBANK_Py)[i]+(CLAS12Data_PBANK_Pz)[i]*(CLAS12Data_PBANK_Pz)[i]);

	//electron info
	if( (CLAS12Data_PBANK_Pid)[i]==11 && (CLAS12Data_PBANK_Charge)[i]<0){
		nElec++;
		i_Elec=i;
		TVector3 Elec_ptemp((CLAS12Data_PBANK_Px)[i],(CLAS12Data_PBANK_Py)[i],(CLAS12Data_PBANK_Pz)[i]);
		Elec_Vec.SetVectM(Elec_ptemp,MELE);
		Float_t nu = EBEAM - pmom[i];
		Q2 = 4*EBEAM*Elec_Vec.P()*TMath::Power(TMath::Sin(Elec_Vec.Theta()/2),2.);
		W = TMath::Sqrt(MPROT*MPROT + 2*MPROT*nu -Q2);
		Xbj = Q2/(2*MPROT*(EBEAM-pmom[i]));
		Elec_Vx = (CLAS12Data_PBANK_Vx)[i];
		Elec_Vy = (CLAS12Data_PBANK_Vy)[i];
		Elec_Vy = (CLAS12Data_PBANK_Vy)[i];
		Elec_P = pmom[i];
		Elec_Theta = Elec_Vec.Theta()*180./TMath::Pi();
		Elec_Phi = Elec_Vec.Phi()*180./TMath::Pi();
	}
	if (nElec>=1){
		//photon info
		if(i!=i_Elec && (CLAS12Data_PBANK_Charge)[i]==0 && (CLAS12Data_PBANK_Pid)[i]==22){
			nPh++;
			i_Ph=i;
			TVector3 Ph_ptemp((CLAS12Data_PBANK_Px)[i],(CLAS12Data_PBANK_Py)[i],(CLAS12Data_PBANK_Pz)[i]);
			Ph_Vec.SetVectM(Ph_ptemp,pmom[i]);
			Ph_P = Ph_Vec.E();
			Ph_Theta = Ph_Vec.Theta()*180./TMath::Pi();
			Ph_Phi = Ph_Vec.Phi()*180./TMath::Pi();
		}
		//proton info
		if(i!=i_Elec && i!=i_Ph && (CLAS12Data_PBANK_Pid)[i]==2212){
			i_Prot=i;
			TVector3 Prot_ptemp((CLAS12Data_PBANK_Px)[i],(CLAS12Data_PBANK_Py)[i],(CLAS12Data_PBANK_Pz)[i]);
			Prot_Vec.SetVectM(Prot_ptemp,MPROT);
			Prot_P = pmom[i];
			Prot_Theta = Prot_Vec.Theta()*180./TMath::Pi();
			Prot_Phi = Prot_Vec.Phi()*180./TMath::Pi();
			Prot_Pz = (CLAS12Data_PBANK_Pz)[i];
			Prot_Vz = (CLAS12Data_PBANK_Vz)[i];
		}
		//pion info
		if(i!=i_Elec && i!=i_Ph && i!=i_Prot && TMath::Abs((CLAS12Data_PBANK_Pid)[i]==211)){
				nPion++;		
		}
		//neutron info
		if(i!=i_Elec && i!=i_Ph && i!=i_Prot && (CLAS12Data_PBANK_Charge)[i]==0 && (CLAS12Data_PBANK_Pid)[i]==2112){
			i_Neut=i;
			TVector3 Neut_ptemp((CLAS12Data_PBANK_Px)[i],(CLAS12Data_PBANK_Py)[i],(CLAS12Data_PBANK_Pz)[i]);
			Neut_Vec.SetVectM(Neut_ptemp,MNEUT);
			Neut_P = pmom[i];
			Neut_Pz = (CLAS12Data_PBANK_Pz)[i];
			Neut_Theta = Neut_Vec.Theta()*180./TMath::Pi();
			Neut_Phi = Neut_Vec.Phi()*180./TMath::Pi();
			if((CLAS12Data_PBANK_Status)[i]>=2000 && (CLAS12Data_PBANK_Status)[i]<4000)FD_Neut=1;
			else if ((CLAS12Data_PBANK_Status)[i]>=4000)CD_Neut=1;
			if (FD_Neut==1 || CD_Neut==1) nNeut++;
		}
		el_tree->Fill();

		if(nProt==1 && nPh>=1){
			Calc_kine_P();
			pDVCS_tree->Fill();
		}
		if(nNeut==1 && nPh>=1 && nProt==0 && nPion==0){
			Calc_kine_N();
			nDVCS_tree->Fill();
		}
	}
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
	outfile->Write();
	outfile->Close();

}
