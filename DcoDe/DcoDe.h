//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun  7 13:53:20 2019 by ROOT version 6.14/06
// from TTree hipotree/tree converted from hipo data
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef DcoDe_h
#define DcoDe_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TVirtualIndex.h>
#include <iostream>             // std::cout, std::endl
#include <fstream>              // std::ifstream
#include <sstream>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TSystem.h>
#include <cstdlib>
#include <map>
// Headers needed by this particular selector
#include <vector>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using namespace std;

class DcoDe : public TSelector {
public :
 //   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   TTree *el_tree,*pDVCS_tree,*nDVCS_tree,*ep_tree;
   TTree *MC_tree;
   TFile *outfile;
   Bool_t   firstCall;
   void    AddBranches();
   void    AddMCBranches();
   void    Calc_kine_P();
   void    Calc_kine_N();
   void    Calc_kine_MC();


// Fixed size dimensions of array or collections stored in the TTree if any.
  Int_t nelec,nprot,nneut,nphotFT,nphotEC,nphot,npion,CD_proton,FD_proton,CD_neutron,FD_neutron;
  Int_t i_el,i_pr,i_n,i_n_cnd,i_phEC,i_phFT,i_ph;
  Int_t isortEC[40],isortIC[40];
  Int_t GoodECPhotIndex,GoodICPhotIndex;
  Int_t RunNumber,FileNumber;
  Int_t pol_sign,helicity;
  Int_t Ph_det;
  Float_t TarPol;
  Float_t Ebeam, Pmass, Nmass, Dmass,LightSpeed;
  Float_t Q2, Xbj, W,t_Pr,t_Ph,t_N,El_P,Phi_Pr,Phi_Ph, Phi_N,Angl_X_g,Angl_hg_hp,Xbal,Ybal,Zbal,Ebal,El_Theta, El_Phi,El_vx,El_vy,El_vz,Pr_P,Pr_Theta,Pr_Phi,Pr_vz,N_P,N_Theta,N_Phi,Ph_EC_P[40],Ph_EC_Theta[40],Ph_EC_Phi[40],Ph_FT_P[40],Ph_FT_Theta[40],Ph_FT_Phi[40],Ph_P,Ph_Theta,Ph_Phi;
  Float_t mm2_epg,mm2_eng,mm2_ep,mm2_eg;
  Float_t El_vx_MC,El_vy_MC,El_vz_MC;
  Float_t El_P_MC,El_Theta_MC,El_Phi_MC;
  Float_t Pr_P_MC,Pr_Theta_MC,Pr_Phi_MC;
  Float_t N_P_MC,N_Theta_MC,N_Phi_MC;
  Float_t Ph_P_MC,Ph_Theta_MC,Ph_Phi_MC;
  Float_t Q2_MC, Xbj_MC, W_MC,t_Pr_MC,t_Ph_MC,t_N_MC,Phi_MC;
  Float_t tstart;
  Int_t n_spec,p_spec;

  TLorentzVector ElectronBeam, Target_Vec,PTarget_Vec,NTarget_Vec;
  TLorentzVector El_Vec,Pr_Vec,N_Vec,Ph_IC_Vec[40],Ph_EC_Vec[40],Ph_Vec;
  TLorentzVector El_Vec_MC,Pr_Vec_MC,N_Vec_MC,Ph_Vec_MC;

   TVector3 pp;
   TVector3 PP[4];


   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<float> CLAS12Data_P_Time = {fReader, "CLAS12Data.P_Time"};
   TTreeReaderArray<float> CLAS12Data_P_Region = {fReader, "CLAS12Data.P_Region"};
   TTreeReaderArray<float> CLAS12Data_P_Theta = {fReader, "CLAS12Data.P_Theta"};
   TTreeReaderArray<float> CLAS12Data_P_Phi = {fReader, "CLAS12Data.P_Phi"};
   TTreeReaderArray<float> CLAS12Data_ECIN_Energy = {fReader, "CLAS12Data.ECIN_Energy"};
   TTreeReaderArray<float> CLAS12Data_ECOUT_Energy = {fReader, "CLAS12Data.ECOUT_Energy"};
   TTreeReaderArray<float> CLAS12Data_CTOF_Time = {fReader, "CLAS12Data.CTOF_Time"};
   TTreeReaderArray<float> CLAS12Data_HTCC_Nphe = {fReader, "CLAS12Data.HTCC_Nphe"};
   TTreeReaderArray<float> CLAS12Data_P_P = {fReader, "CLAS12Data.P_P"};


   TTreeReaderArray<long> CLAS12Data_CND1_hits_id = {fReader, "CLAS12Data.CND1_hits_id"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_status = {fReader, "CLAS12Data.CND1_hits_status"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_trkID = {fReader, "CLAS12Data.CND1_hits_trkID"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_sector = {fReader, "CLAS12Data.CND1_hits_sector"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_layer = {fReader, "CLAS12Data.CND1_hits_layer"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_component = {fReader, "CLAS12Data.CND1_hits_component"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_indexLadc = {fReader, "CLAS12Data.CND1_hits_indexLadc"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_indexRadc = {fReader, "CLAS12Data.CND1_hits_indexRadc"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_indexLtdc = {fReader, "CLAS12Data.CND1_hits_indexLtdc"};
   TTreeReaderArray<long> CLAS12Data_CND1_hits_indexRtdc = {fReader, "CLAS12Data.CND1_hits_indexRtdc"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_energy = {fReader, "CLAS12Data.CND1_hits_energy"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_time = {fReader, "CLAS12Data.CND1_hits_time"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_energy_unc = {fReader, "CLAS12Data.CND1_hits_energy_unc"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_time_unc = {fReader, "CLAS12Data.CND1_hits_time_unc"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_x = {fReader, "CLAS12Data.CND1_hits_x"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_y = {fReader, "CLAS12Data.CND1_hits_y"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_z = {fReader, "CLAS12Data.CND1_hits_z"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_x_unc = {fReader, "CLAS12Data.CND1_hits_x_unc"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_y_unc = {fReader, "CLAS12Data.CND1_hits_y_unc"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_z_unc = {fReader, "CLAS12Data.CND1_hits_z_unc"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_tx = {fReader, "CLAS12Data.CND1_hits_tx"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_ty = {fReader, "CLAS12Data.CND1_hits_ty"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_tz = {fReader, "CLAS12Data.CND1_hits_tz"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_tlength = {fReader, "CLAS12Data.CND1_hits_tlength"};
   TTreeReaderArray<float> CLAS12Data_CND1_hits_pathlength = {fReader, "CLAS12Data.CND1_hits_pathlength"};

//   TTreeReaderArray<long> RICH_tdc_sector = {fReader, "RICH_tdc_sector"};
//   TTreeReaderArray<long> RICH_tdc_layer = {fReader, "RICH_tdc_layer"};
//   TTreeReaderArray<long> RICH_tdc_component = {fReader, "RICH_tdc_component"};
//   TTreeReaderArray<long> RICH_tdc_order = {fReader, "RICH_tdc_order"};
//   TTreeReaderArray<long> RICH_tdc_TDC = {fReader, "RICH_tdc_TDC"};
//
//   TTreeReaderArray<long> BAND_hits_id = {fReader, "BAND_hits_id"};
//   TTreeReaderArray<long> BAND_hits_sector = {fReader, "BAND_hits_sector"};
//   TTreeReaderArray<long> BAND_hits_layer = {fReader, "BAND_hits_layer"};
//   TTreeReaderArray<long> BAND_hits_component = {fReader, "BAND_hits_component"};
//   TTreeReaderArray<float> BAND_hits_meantimeTdc = {fReader, "BAND_hits_meantimeTdc"};
//   TTreeReaderArray<float> BAND_hits_meantimeFadc = {fReader, "BAND_hits_meantimeFadc"};
//   TTreeReaderArray<float> BAND_hits_difftimeTdc = {fReader, "BAND_hits_difftimeTdc"};
//   TTreeReaderArray<float> BAND_hits_difftimeFadc = {fReader, "BAND_hits_difftimeFadc"};
//   TTreeReaderArray<float> BAND_hits_adcLcorr = {fReader, "BAND_hits_adcLcorr"};
//   TTreeReaderArray<float> BAND_hits_adcRcorr = {fReader, "BAND_hits_adcRcorr"};
//   TTreeReaderArray<float> BAND_hits_tFadcLcorr = {fReader, "BAND_hits_tFadcLcorr"};
//   TTreeReaderArray<float> BAND_hits_tFadcRcorr = {fReader, "BAND_hits_tFadcRcorr"};
//   TTreeReaderArray<float> BAND_hits_tTdcLcorr = {fReader, "BAND_hits_tTdcLcorr"};
//   TTreeReaderArray<float> BAND_hits_tTdcRcorr = {fReader, "BAND_hits_tTdcRcorr"};
//   TTreeReaderArray<float> BAND_hits_x = {fReader, "BAND_hits_x"};
//   TTreeReaderArray<float> BAND_hits_y = {fReader, "BAND_hits_y"};
//   TTreeReaderArray<float> BAND_hits_z = {fReader, "BAND_hits_z"};
//   TTreeReaderArray<float> BAND_hits_ux = {fReader, "BAND_hits_ux"};
//   TTreeReaderArray<float> BAND_hits_uy = {fReader, "BAND_hits_uy"};
//   TTreeReaderArray<float> BAND_hits_uz = {fReader, "BAND_hits_uz"};
//
//   TTreeReaderArray<long> RUN_config_run = {fReader, "RUN_config_run"};
//   TTreeReaderArray<long> RUN_config_event = {fReader, "RUN_config_event"};
//   TTreeReaderArray<long> RUN_config_unixtime = {fReader, "RUN_config_unixtime"};
//   TTreeReaderArray<long> RUN_config_trigger = {fReader, "RUN_config_trigger"};
//   TTreeReaderArray<long> RUN_config_timestamp = {fReader, "RUN_config_timestamp"};
//   TTreeReaderArray<long> RUN_config_type = {fReader, "RUN_config_type"};
//   TTreeReaderArray<long> RUN_config_mode = {fReader, "RUN_config_mode"};
//   TTreeReaderArray<float> RUN_config_torus = {fReader, "RUN_config_torus"};
//   TTreeReaderArray<float> RUN_config_solenoid = {fReader, "RUN_config_solenoid"};

   TTreeReaderArray<long> CLAS12Data_header_rn = {fReader, "CLAS12Data.header_rn"};
   TTreeReaderArray<long> CLAS12Data_header_en = {fReader, "CLAS12Data.header_en"};
   TTreeReaderArray<long> CLAS12Data_header_ty = {fReader, "CLAS12Data.header_ty"};
   TTreeReaderArray<long> CLAS12Data_header_ec = {fReader, "CLAS12Data.header_ec"};
   TTreeReaderArray<long> CLAS12Data_header_np = {fReader, "CLAS12Data.header_np"};
   TTreeReaderArray<long> CLAS12Data_header_trg = {fReader, "CLAS12Data.header_trg"};
   TTreeReaderArray<long> CLAS12Data_header_hel = {fReader, "CLAS12Data.header_hel"};
   TTreeReaderArray<float> CLAS12Data_header_et = {fReader, "CLAS12Data.header_et"};
   TTreeReaderArray<float> CLAS12Data_header_bcg = {fReader, "CLAS12Data.header_bcg"};
   TTreeReaderArray<float> CLAS12Data_header_lt = {fReader, "CLAS12Data.header_lt"};
   TTreeReaderArray<float> CLAS12Data_header_st = {fReader, "CLAS12Data.header_st"};
   TTreeReaderArray<float> CLAS12Data_header_rf = {fReader, "CLAS12Data.header_rf"};
   TTreeReaderArray<float> CLAS12Data_header_pt = {fReader, "CLAS12Data.header_pt"};

   TTreeReaderArray<long> CLAS12Data_PBANK_Pid = {fReader, "CLAS12Data.PBANK_Pid"};
   TTreeReaderArray<long> CLAS12Data_PBANK_Charge = {fReader, "CLAS12Data.PBANK_Charge"};
   TTreeReaderArray<long> CLAS12Data_PBANK_Status = {fReader, "CLAS12Data_PBANK_Status"};
   TTreeReaderArray<float> CLAS12Data_PBANK_Px = {fReader, "CLAS12Data.PBANK_Px"};
   TTreeReaderArray<float> CLAS12Data_PBANK_Py = {fReader, "CLAS12Data.PBANK_Py"};
   TTreeReaderArray<float> CLAS12Data_PBANK_Pz = {fReader, "CLAS12Data.PBANK_Pz"};
   TTreeReaderArray<float> CLAS12Data_PBANK_Vx = {fReader, "CLAS12Data.PBANK_Vx"};
   TTreeReaderArray<float> CLAS12Data_PBANK_Vy = {fReader, "CLAS12Data.PBANK_Vy"};
   TTreeReaderArray<float> CLAS12Data_PBANK_Vz = {fReader, "CLAS12Data.PBANK_Vz"};
   TTreeReaderArray<float> CLAS12_PBANK_Beta = {fReader, "CLAS12Data.PBANK_Beta"};
   TTreeReaderArray<float> CLAS12Data_PBANK_Chi2pid = {fReader, "CLAS12Data.PBANK_Chi2pid"};

   TTreeReaderArray<long> CLAS12Data_calorimeter_index = {fReader, "CLAS12Data.calorimeter_index"};
   TTreeReaderArray<long> CLAS12Data_calorimeter_pindex = {fReader, "CLAS12Data.calorimeter_pindex"};
   TTreeReaderArray<long> CLAS12Data_calorimeter_detector = {fReader, "CLAS12Data.calorimeter_detector"};
   TTreeReaderArray<long> CLAS12Data_calorimeter_sector = {fReader, "CLAS12Data.calorimeter_sector"};
   TTreeReaderArray<long> CLAS12Data_calorimeter_layer = {fReader, "CLAS12Data.calorimeter_layer"};
   TTreeReaderArray<long> CLAS12Data_calorimeter_status = {fReader, "CLAS12Data.calorimeter_status"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_energy = {fReader, "CLAS12Data.calorimeter_energy"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_time = {fReader, "CLAS12Data.calorimeter_time"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_path = {fReader, "CLAS12Data.calorimeter_path"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_chi2 = {fReader, "CLAS12Data.calorimeter_chi2"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_x = {fReader, "CLAS12Data.calorimeter_x"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_y = {fReader, "CLAS12Data.calorimeter_y"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_z = {fReader, "CLAS12Data.calorimeter_z"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_hx = {fReader, "CLAS12Data.calorimeter_hx"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_hy = {fReader, "CLAS12Data.calorimeter_hy"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_hz = {fReader, "CLAS12Data.calorimeter_hz"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_lu = {fReader, "CLAS12Data.calorimeter_lu"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_lv = {fReader, "CLAS12Data.calorimeter_lv"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_lw = {fReader, "CLAS12Data.calorimeter_lw"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_du = {fReader, "CLAS12Data.calorimeter_du"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_dv = {fReader, "CLAS12Data.calorimeter_dv"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_dw = {fReader, "CLAS12Data.calorimeter_dw"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_m2u = {fReader, "CLAS12Data.calorimeter_m2u"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_m2v = {fReader, "CLAS12Data.calorimeter_m2v"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_m2w = {fReader, "CLAS12Data.calorimeter_m2w"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_m3u = {fReader, "CLAS12Data.calorimeter_m3u"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_m3v = {fReader, "CLAS12Data.calorimeter_m3v"};
   TTreeReaderArray<float> CLAS12Data_calorimeter_m3w = {fReader, "CLAS12Data.calorimeter_m3w"};

   TTreeReaderArray<long> CLAS12Data_HTCC_index = {fReader, "CLAS12Data.HTCC_index"};
   TTreeReaderArray<long> CLAS12Data_HTCC_pindex = {fReader, "CLAS12Data.HTCC_pindex"};
   TTreeReaderArray<long> CLAS12Data_HTCC_detector = {fReader, "CLAS12Data.HTCC_detector"};
   TTreeReaderArray<long> CLAS12Data_HTCC_sector = {fReader, "CLAS12Data.HTCC_sector"};
   TTreeReaderArray<long> CLAS12Data_HTCC_status = {fReader, "CLAS12Data.HTCC_status"};
   TTreeReaderArray<float> CLAS12Data_HTCC_nphe = {fReader, "CLAS12Data.HTCC_nphe"};
   TTreeReaderArray<float> CLAS12Data_HTCC_time = {fReader, "CLAS12Data.HTCC_time"};
   TTreeReaderArray<float> CLAS12Data_HTCC_path = {fReader, "CLAS12Data.HTCC_path"};
   TTreeReaderArray<float> CLAS12Data_HTCC_chi2 = {fReader, "CLAS12Data.HTCC_chi2"};
   TTreeReaderArray<float> CLAS12Data_HTCC_x = {fReader, "CLAS12Data.HTCC_x"};
   TTreeReaderArray<float> CLAS12Data_HTCC_y = {fReader, "CLAS12Data.HTCC_y"};
   TTreeReaderArray<float> CLAS12Data_HTCC_z = {fReader, "CLAS12Data.HTCC_z"};
   TTreeReaderArray<float> CLAS12Data_HTCC_theta = {fReader, "CLAS12Data.HTCC_theta"};
   TTreeReaderArray<float> CLAS12Data_HTCC_phi = {fReader, "CLAS12Data.HTCC_phi"};
   TTreeReaderArray<float> CLAS12Data_HTCC_dtheta = {fReader, "CLAS12Data.HTCC_dtheta"};
   TTreeReaderArray<float> CLAS12Data_HTCC_dphi = {fReader, "CLAS12Data.HTCC_dphi"};

   TTreeReaderArray<long> CLAS12Data_forwardtagger_index = {fReader, "CLAS12Data.forwardtagger_index"};
   TTreeReaderArray<long> CLAS12Data_forwardtagger_pindex = {fReader, "CLAS12Data.forwardtagger_pindex"};
   TTreeReaderArray<long> CLAS12Data_forwardtagger_detector = {fReader, "CLAS12Data.forwardtagger_detector"};
   TTreeReaderArray<long> CLAS12Data_forwardtagger_size = {fReader, "CLAS12Data.forwardtagger_size"};
   TTreeReaderArray<long> CLAS12Data_forwardtagger_status = {fReader, "CLAS12Data.forwardtagger_status"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_energy = {fReader, "CLAS12Data.forwardtagger_energy"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_time = {fReader, "CLAS12Data.forwardtagger_time"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_path = {fReader, "CLAS12Data.forwardtagger_path"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_chi2 = {fReader, "CLAS12Data.forwardtagger_chi2"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_x = {fReader, "CLAS12Data.forwardtagger_x"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_y = {fReader, "CLAS12Data.forwardtagger_y"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_z = {fReader, "CLAS12Data.forwardtagger_z"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_dx = {fReader, "CLAS12Data.forwardtagger_dx"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_dy = {fReader, "CLAS12Data.forwardtagger_dy"};
   TTreeReaderArray<float> CLAS12Data_forwardtagger_radius = {fReader, "CLAS12Data.forwardtagger_radius"};

   TTreeReaderArray<long> CLAS12Data_Scintillator_index = {fReader, "CLAS12Data.Scintillator_index"};
   TTreeReaderArray<long> CLAS12Data_Scintillator_pindex = {fReader, "CLAS12Data.Scintillator_pindex"};
   TTreeReaderArray<long> CLAS12Data_Scintillator_detector = {fReader, "CLAS12Data.Scintillator_detector"};
   TTreeReaderArray<long> CLAS12Data_Scintillator_sector = {fReader, "CLAS12Data.Scintillator_sector"};
   TTreeReaderArray<long> CLAS12Data_Scintillator_layer = {fReader, "CLAS12Data.Scintillator_layer"};
   TTreeReaderArray<long> CLAS12Data_Scintillator_component = {fReader, "CLAS12Data.Scintillator_component"};
   TTreeReaderArray<long> CLAS12Data_Scintillator_status = {fReader, "CLAS12Data.Scintillator_status"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_energy = {fReader, "CLAS12Data.Scintillator_energy"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_time = {fReader, "CLAS12Data.Scintillator_time"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_path = {fReader, "CLAS12Data.Scintillator_path"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_chi2 = {fReader, "CLAS12Data.Scintillator_chi2"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_x = {fReader, "CLAS12Data.Scintillator_x"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_y = {fReader, "CLAS12Data.Scintillator_y"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_z = {fReader, "CLAS12Data.Scintillator_z"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_hx = {fReader, "CLAS12Data.Scintillator_hx"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_hy = {fReader, "CLAS12Data.Scintillator_hy"};
   TTreeReaderArray<float> CLAS12Data_Scintillator_hz = {fReader, "CLAS12Data.Scintillator_hz"};

   TTreeReaderArray<long> CLAS12Data_Track_index = {fReader, "CLAS12Data.Track_index"};
   TTreeReaderArray<long> CLAS12Data_Track_pindex = {fReader, "CLAS12Data.Track_pindex"};
   TTreeReaderArray<long> CLAS12Data_Track_detector = {fReader, "CLAS12Data.Track_detector"};
   TTreeReaderArray<long> CLAS12Data_Track_sector = {fReader, "CLAS12Data.Track_sector"};
   TTreeReaderArray<long> CLAS12Data_Track_status = {fReader, "CLAS12Data.Track_status"};
   TTreeReaderArray<long> CLAS12Data_Track_q = {fReader, "CLAS12Data.Track_q"};
   TTreeReaderArray<long> CLAS12Data_Track_NDF = {fReader, "CLAS12Data.Track_NDF"};
   TTreeReaderArray<long> CLAS12Data_Track_NDFnomm = {fReader, "CLAS12Data.Track_NDFnomm"};
   TTreeReaderArray<float> CLAS12Data_Track_chi2 = {fReader, "CLAS12Data.Track_chi2"};
   TTreeReaderArray<float> CLAS12Data_Track_pxnomm = {fReader, "CLAS12Data.Track_pxnomm"};
   TTreeReaderArray<float> CLAS12Data_Track_pynomm = {fReader, "CLAS12Data.Track_pynomm"};
   TTreeReaderArray<float> CLAS12Data_Track_pznomm = {fReader, "CLAS12Data.Track_pznomm"};
   TTreeReaderArray<float> CLAS12Data_Track_vxnomm = {fReader, "CLAS12Data.Track_vxnomm"};
   TTreeReaderArray<float> CLAS12Data_Track_vynomm = {fReader, "CLAS12Data.Track_vynomm"};
   TTreeReaderArray<float> CLAS12Data_Track_vznomm = {fReader, "CLAS12Data.Track_vznomm"};
   TTreeReaderArray<float> CLAS12Data_Track_chi2nomm = {fReader, "CLAS12Data.Track_chi2nomm"};

   TTreeReaderArray<long> CLAS12Data_covmatrix_index = {fReader, "CLAS12Data.covmatrix_index"};
   TTreeReaderArray<long> CLAS12Data_covmatrix_pindex = {fReader, "CLAS12Data.covmatrix_pindex"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C11 = {fReader, "CLAS12Data.covmatrix_C11"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C12 = {fReader, "CLAS12Data.covmatrix_C12"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C13 = {fReader, "CLAS12Data.covmatrix_C13"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C14 = {fReader, "CLAS12Data.covmatrix_C14"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C15 = {fReader, "CLAS12Data.covmatrix_C15"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C22 = {fReader, "CLAS12Data.covmatrix_C22"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C23 = {fReader, "CLAS12Data.covmatrix_C23"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C24 = {fReader, "CLAS12Data.covmatrix_C24"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C25 = {fReader, "CLAS12Data.covmatrix_C25"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C33 = {fReader, "CLAS12Data.covmatrix_C33"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C34 = {fReader, "CLAS12Data.covmatrix_C34"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C35 = {fReader, "CLAS12Data.covmatrix_C35"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C44 = {fReader, "CLAS12Data.covmatrix_C44"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C45 = {fReader, "CLAS12Data.covmatrix_C45"};
   TTreeReaderArray<float> CLAS12Data_covmatrix_C55 = {fReader, "CLAS12Data.covmatrix_C55"};

   TTreeReaderArray<long> CLAS12Data_traj_pindex = {fReader, "CLAS12Data.traj_pindex"};
   TTreeReaderArray<long> CLAS12Data_traj_index = {fReader, "CLAS12Data.traj_index"};
   TTreeReaderArray<long> CLAS12Data_traj_detId = {fReader, "CLAS12Data.traj_detId"};
   TTreeReaderArray<long> CLAS12Data_traj_q = {fReader, "CLAS12Data.traj_q"};
   TTreeReaderArray<float> CLAS12Data_traj_x = {fReader, "CLAS12Data.traj_x"};
   TTreeReaderArray<float> CLAS12Data_traj_y = {fReader, "CLAS12Data.traj_y"};
   TTreeReaderArray<float> CLAS12Data_traj_z = {fReader, "CLAS12Data.traj_z"};
   TTreeReaderArray<float> CLAS12Data_traj_cx = {fReader, "CLAS12Data.traj_cx"};
   TTreeReaderArray<float> CLAS12Data_traj_cy = {fReader, "CLAS12Data.traj_cy"};
   TTreeReaderArray<float> CLAS12Data_traj_cz = {fReader, "CLAS12Data.traj_cz"};
   TTreeReaderArray<float> CLAS12Data_traj_pathlength = {fReader, "CLAS12Data.traj_pathlength"};


   DcoDe(TTree * /*tree*/ =0) { }
   virtual ~DcoDe() { }
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

   ClassDef(DcoDe,0);

};

#endif

#ifdef DcoDe_cxx
void DcoDe::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t DcoDe::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef DcoDe_cxx
