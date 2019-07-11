//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 14 11:06:19 2019 by ROOT version 6.16/00
// from TTree hipotree/tree converted from hipo data
// found on file: 40.root
//////////////////////////////////////////////////////////

#ifndef MyDcoDe_h
#define MyDcoDe_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>



class MyDcoDe : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Pid = {fReader, "CLAS12Data.PBANK_Pid"};
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Charge = {fReader, "CLAS12Data.PBANK_Charge"};
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Status = {fReader, "CLAS12Data.PBANK_Status"};
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Px = {fReader, "CLAS12Data.PBANK_Px"};
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Py = {fReader, "CLAS12Data.PBANK_Py"};
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Pz = {fReader, "CLAS12Data.PBANK_Pz"};
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Vx = {fReader, "CLAS12Data.PBANK_Vx"};
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Vy = {fReader, "CLAS12Data.PBANK_Vy"};
   TTreeReaderArray<Float_t> CLAS12Data_PBANK_Vz = {fReader, "CLAS12Data.PBANK_Vz"};
   TTreeReaderArray<Float_t> CLAS12Data_P_Time = {fReader, "CLAS12Data.P_Time"};
   TTreeReaderArray<Float_t> CLAS12Data_P_Theta = {fReader, "CLAS12Data.P_Theta"};
   TTreeReaderArray<Float_t> CLAS12Data_P_Phi = {fReader, "CLAS12Data.P_Phi"};
   TTreeReaderArray<Float_t> CLAS12Data_CTOF_Time = {fReader, "CLAS12Data.CTOF_Time"};
   TTreeReaderArray<Float_t> CLAS12Data_P_Region = {fReader, "CLAS12Data.P_Region"};
   TTreeReaderArray<Float_t> CLAS12Data_P_P = {fReader, "CLAS12Data.P_P"};
   TTreeReaderArray<Int_t> CLAS12Data_EVNT_RunNumber = {fReader, "CLAS12Data.EVNT_RunNumber"};
   TTreeReaderArray<Int_t> CLAS12Data_EVNT_EventNumber = {fReader, "CLAS12Data.EVNT_EventNumber"};
   TTreeReaderArray<Int_t> CLAS12Data_EVNT_Helicity = {fReader, "CLAS12Data.EVNT_Helicity"};
   TTreeReaderArray<Int_t> CLAS12Data_EVNT_Type = {fReader, "CLAS12Data.EVNT_Type"};
   TTreeReaderArray<Float_t> CLAS12Data_EVNT_StartTime = {fReader, "CLAS12Data.EVNT_StartTime"};
   TTreeReaderArray<Float_t> CLAS12Data_EVNT_RFTime = {fReader, "CLAS12Data.EVNT_RFTime"};
   TTreeReaderArray<Long64_t> CLAS12Data_EVNT_Trigger = {fReader, "CLAS12Data.EVNT_Trigger"};
   TTreeReaderArray<Int_t> CLAS12Data_CND1_Index = {fReader, "CLAS12Data.CND1_Index"};
   TTreeReaderArray<Int_t> CLAS12Data_CND1_Pindex = {fReader, "CLAS12Data.CND1_Pindex"};
   TTreeReaderArray<Float_t> CLAS12Data_CND1_Status = {fReader, "CLAS12Data.CND1_Status"};
   TTreeReaderArray<Int_t> CLAS12Data_CND1_Sector = {fReader, "CLAS12Data.CND1_Sector"};
   TTreeReaderArray<Int_t> CLAS12Data_CND1_Layer = {fReader, "CLAS12Data.CND1_Layer"};
   TTreeReaderArray<Float_t> CLAS12Data_CND1_Energy = {fReader, "CLAS12Data.CND1_Energy"};
   TTreeReaderArray<Float_t> CLAS12Data_CND1_Time = {fReader, "CLAS12Data.CND1_Time"};
   TTreeReaderArray<Float_t> CLAS12Data_CND1_X = {fReader, "CLAS12Data.CND1_X"};
   TTreeReaderArray<Float_t> CLAS12Data_CND1_Y = {fReader, "CLAS12Data.CND1_Y"};
   TTreeReaderArray<Float_t> CLAS12Data_CND1_Z = {fReader, "CLAS12Data.CND1_Z"};
   TTreeReaderArray<Int_t> CLAS12Data_ECIN_Index = {fReader, "CLAS12Data.ECIN_Index"};
   TTreeReaderArray<Int_t> CLAS12Data_ECIN_Pindex = {fReader, "CLAS12Data.ECIN_Pindex"};
   TTreeReaderArray<Int_t> CLAS12Data_ECIN_Sector = {fReader, "CLAS12Data.ECIN_Sector"};
   TTreeReaderArray<Int_t> CLAS12Data_ECIN_Layer = {fReader, "CLAS12Data.ECIN_Layer"};
   TTreeReaderArray<Int_t> CLAS12Data_ECIN_Status = {fReader, "CLAS12Data.ECIN_Status"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Energy = {fReader, "CLAS12Data.ECIN_Energy"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Time = {fReader, "CLAS12Data.ECIN_Time"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Path = {fReader, "CLAS12Data.ECIN_Path"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_X = {fReader, "CLAS12Data.ECIN_X"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Y = {fReader, "CLAS12Data.ECIN_Y"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Z = {fReader, "CLAS12Data.ECIN_Z"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Hx = {fReader, "CLAS12Data.ECIN_Hx"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Hy = {fReader, "CLAS12Data.ECIN_Hy"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Hz = {fReader, "CLAS12Data.ECIN_Hz"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Du = {fReader, "CLAS12Data.ECIN_Du"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Dv = {fReader, "CLAS12Data.ECIN_Dv"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_Dw = {fReader, "CLAS12Data.ECIN_Dw"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_M2u = {fReader, "CLAS12Data.ECIN_M2u"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_M2v = {fReader, "CLAS12Data.ECIN_M2v"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_M2w = {fReader, "CLAS12Data.ECIN_M2w"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_M3u = {fReader, "CLAS12Data.ECIN_M3u"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_M3v = {fReader, "CLAS12Data.ECIN_M3v"};
   TTreeReaderArray<Float_t> CLAS12Data_ECIN_M3w = {fReader, "CLAS12Data.ECIN_M3w"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Index = {fReader, "CLAS12Data.ECOUT_Index"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Pindex = {fReader, "CLAS12Data.ECOUT_Pindex"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Sector = {fReader, "CLAS12Data.ECOUT_Sector"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Layer = {fReader, "CLAS12Data.ECOUT_Layer"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Status = {fReader, "CLAS12Data.ECOUT_Status"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Energy = {fReader, "CLAS12Data.ECOUT_Energy"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Time = {fReader, "CLAS12Data.ECOUT_Time"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Path = {fReader, "CLAS12Data.ECOUT_Path"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_X = {fReader, "CLAS12Data.ECOUT_X"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Y = {fReader, "CLAS12Data.ECOUT_Y"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Z = {fReader, "CLAS12Data.ECOUT_Z"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Hx = {fReader, "CLAS12Data.ECOUT_Hx"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Hy = {fReader, "CLAS12Data.ECOUT_Hy"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Hz = {fReader, "CLAS12Data.ECOUT_Hz"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Du = {fReader, "CLAS12Data.ECOUT_Du"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Dv = {fReader, "CLAS12Data.ECOUT_Dv"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_Dw = {fReader, "CLAS12Data.ECOUT_Dw"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_M2u = {fReader, "CLAS12Data.ECOUT_M2u"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_M2v = {fReader, "CLAS12Data.ECOUT_M2v"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_M2w = {fReader, "CLAS12Data.ECOUT_M2w"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_M3u = {fReader, "CLAS12Data.ECOUT_M3u"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_M3v = {fReader, "CLAS12Data.ECOUT_M3v"};
   TTreeReaderArray<Float_t> CLAS12Data_ECOUT_M3w = {fReader, "CLAS12Data.ECOUT_M3w"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Index = {fReader, "CLAS12Data.PCAL_Index"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Pindex = {fReader, "CLAS12Data.PCAL_Pindex"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Sector = {fReader, "CLAS12Data.PCAL_Sector"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Layer = {fReader, "CLAS12Data.PCAL_Layer"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Status = {fReader, "CLAS12Data.PCAL_Status"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Energy = {fReader, "CLAS12Data.PCAL_Energy"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Time = {fReader, "CLAS12Data.PCAL_Time"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Path = {fReader, "CLAS12Data.PCAL_Path"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_X = {fReader, "CLAS12Data.PCAL_X"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Y = {fReader, "CLAS12Data.PCAL_Y"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Z = {fReader, "CLAS12Data.PCAL_Z"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Hx = {fReader, "CLAS12Data.PCAL_Hx"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Hy = {fReader, "CLAS12Data.PCAL_Hy"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Hz = {fReader, "CLAS12Data.PCAL_Hz"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Du = {fReader, "CLAS12Data.PCAL_Du"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Dv = {fReader, "CLAS12Data.PCAL_Dv"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_Dw = {fReader, "CLAS12Data.PCAL_Dw"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_M2u = {fReader, "CLAS12Data.PCAL_M2u"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_M2v = {fReader, "CLAS12Data.PCAL_M2v"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_M2w = {fReader, "CLAS12Data.PCAL_M2w"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_M3u = {fReader, "CLAS12Data.PCAL_M3u"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_M3v = {fReader, "CLAS12Data.PCAL_M3v"};
   TTreeReaderArray<Float_t> CLAS12Data_PCAL_M3w = {fReader, "CLAS12Data.PCAL_M3w"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Index = {fReader, "CLAS12Data.HTCC_Index"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Pindex = {fReader, "CLAS12Data.HTCC_Pindex"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Sector = {fReader, "CLAS12Data.HTCC_Sector"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Status = {fReader, "CLAS12Data.HTCC_Status"};
   TTreeReaderArray<Int_t> CLAS12Data_HTCC_Nphe = {fReader, "CLAS12Data.HTCC_Nphe"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Time = {fReader, "CLAS12Data.HTCC_Time"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Path = {fReader, "CLAS12Data.HTCC_Path"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_X = {fReader, "CLAS12Data.HTCC_X"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Y = {fReader, "CLAS12Data.HTCC_Y"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Z = {fReader, "CLAS12Data.HTCC_Z"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Theta = {fReader, "CLAS12Data.HTCC_Theta"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Phi = {fReader, "CLAS12Data.HTCC_Phi"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_Dtheta = {fReader, "CLAS12Data.HTCC_Dtheta"};
   TTreeReaderArray<Float_t> CLAS12Data_HTCC_DPhi = {fReader, "CLAS12Data.HTCC_DPhi"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Index = {fReader, "CLAS12Data.FTCAL_Index"};
   TTreeReaderArray<Int_t> CLAS12Data_FTCAL_Pindex = {fReader, "CLAS12Data.FTCAL_Pindex"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Detector = {fReader, "CLAS12Data.FTCAL_Detector"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Energy = {fReader, "CLAS12Data.FTCAL_Energy"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Time = {fReader, "CLAS12Data.FTCAL_Time"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Path = {fReader, "CLAS12Data.FTCAL_Path"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Chi2 = {fReader, "CLAS12Data.FTCAL_Chi2"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_X = {fReader, "CLAS12Data.FTCAL_X"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Y = {fReader, "CLAS12Data.FTCAL_Y"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Z = {fReader, "CLAS12Data.FTCAL_Z"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Dx = {fReader, "CLAS12Data.FTCAL_Dx"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Dy = {fReader, "CLAS12Data.FTCAL_Dy"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Radius = {fReader, "CLAS12Data.FTCAL_Radius"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Size = {fReader, "CLAS12Data.FTCAL_Size"};
   TTreeReaderArray<Float_t> CLAS12Data_FTCAL_Status = {fReader, "CLAS12Data.FTCAL_Status"};


   MyDcoDe(TTree * /*tree*/ =0) { }
   virtual ~MyDcoDe() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(MyDcoDe,0);

};

#endif

#ifdef MyDcoDe_cxx
void MyDcoDe::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t MyDcoDe::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef MyDcoDe_cxx
