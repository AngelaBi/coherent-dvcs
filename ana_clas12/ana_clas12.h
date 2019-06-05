//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 26 11:09:57 2018 by ROOT version 6.02/02
// from TTree hipo2root/CLAS12 banks in ROOT
// found on file: out__171218_095035_events_deut_run31127632.dat.root
//////////////////////////////////////////////////////////

#ifndef ana_clas12_h
#define ana_clas12_h

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
// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

using namespace std;

class ana_clas12 : public TSelector {
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

   TTreeReaderArray<long> CND_hits_id = {fReader, "CND_hits_id"};
   TTreeReaderArray<long> CND_hits_status = {fReader, "CND_hits_status"};
   TTreeReaderArray<long> CND_hits_trkID = {fReader, "CND_hits_trkID"};
   TTreeReaderArray<long> CND_hits_sector = {fReader, "CND_hits_sector"};
   TTreeReaderArray<long> CND_hits_layer = {fReader, "CND_hits_layer"};
   TTreeReaderArray<long> CND_hits_component = {fReader, "CND_hits_component"};
   TTreeReaderArray<long> CND_hits_indexLadc = {fReader, "CND_hits_indexLadc"};
   TTreeReaderArray<long> CND_hits_indexRadc = {fReader, "CND_hits_indexRadc"};
   TTreeReaderArray<long> CND_hits_indexLtdc = {fReader, "CND_hits_indexLtdc"};
   TTreeReaderArray<long> CND_hits_indexRtdc = {fReader, "CND_hits_indexRtdc"};
   TTreeReaderArray<float> CND_hits_energy = {fReader, "CND_hits_energy"};
   TTreeReaderArray<float> CND_hits_time = {fReader, "CND_hits_time"};
   TTreeReaderArray<float> CND_hits_energy_unc = {fReader, "CND_hits_energy_unc"};
   TTreeReaderArray<float> CND_hits_time_unc = {fReader, "CND_hits_time_unc"};
   TTreeReaderArray<float> CND_hits_x = {fReader, "CND_hits_x"};
   TTreeReaderArray<float> CND_hits_y = {fReader, "CND_hits_y"};
   TTreeReaderArray<float> CND_hits_z = {fReader, "CND_hits_z"};
   TTreeReaderArray<float> CND_hits_x_unc = {fReader, "CND_hits_x_unc"};
   TTreeReaderArray<float> CND_hits_y_unc = {fReader, "CND_hits_y_unc"};
   TTreeReaderArray<float> CND_hits_z_unc = {fReader, "CND_hits_z_unc"};
   TTreeReaderArray<float> CND_hits_tx = {fReader, "CND_hits_tx"};
   TTreeReaderArray<float> CND_hits_ty = {fReader, "CND_hits_ty"};
   TTreeReaderArray<float> CND_hits_tz = {fReader, "CND_hits_tz"};
   TTreeReaderArray<float> CND_hits_tlength = {fReader, "CND_hits_tlength"};
   TTreeReaderArray<float> CND_hits_pathlength = {fReader, "CND_hits_pathlength"};
   TTreeReaderArray<long> RICH_tdc_sector = {fReader, "RICH_tdc_sector"};
   TTreeReaderArray<long> RICH_tdc_layer = {fReader, "RICH_tdc_layer"};
   TTreeReaderArray<long> RICH_tdc_component = {fReader, "RICH_tdc_component"};
   TTreeReaderArray<long> RICH_tdc_order = {fReader, "RICH_tdc_order"};
   TTreeReaderArray<long> RICH_tdc_TDC = {fReader, "RICH_tdc_TDC"};
   TTreeReaderArray<long> BAND_hits_id = {fReader, "BAND_hits_id"};
   TTreeReaderArray<long> BAND_hits_sector = {fReader, "BAND_hits_sector"};
   TTreeReaderArray<long> BAND_hits_layer = {fReader, "BAND_hits_layer"};
   TTreeReaderArray<long> BAND_hits_component = {fReader, "BAND_hits_component"};
   TTreeReaderArray<float> BAND_hits_meantimeTdc = {fReader, "BAND_hits_meantimeTdc"};
   TTreeReaderArray<float> BAND_hits_meantimeFadc = {fReader, "BAND_hits_meantimeFadc"};
   TTreeReaderArray<float> BAND_hits_difftimeTdc = {fReader, "BAND_hits_difftimeTdc"};
   TTreeReaderArray<float> BAND_hits_difftimeFadc = {fReader, "BAND_hits_difftimeFadc"};
   TTreeReaderArray<float> BAND_hits_adcLcorr = {fReader, "BAND_hits_adcLcorr"};
   TTreeReaderArray<float> BAND_hits_adcRcorr = {fReader, "BAND_hits_adcRcorr"};
   TTreeReaderArray<float> BAND_hits_tFadcLcorr = {fReader, "BAND_hits_tFadcLcorr"};
   TTreeReaderArray<float> BAND_hits_tFadcRcorr = {fReader, "BAND_hits_tFadcRcorr"};
   TTreeReaderArray<float> BAND_hits_tTdcLcorr = {fReader, "BAND_hits_tTdcLcorr"};
   TTreeReaderArray<float> BAND_hits_tTdcRcorr = {fReader, "BAND_hits_tTdcRcorr"};
   TTreeReaderArray<float> BAND_hits_x = {fReader, "BAND_hits_x"};
   TTreeReaderArray<float> BAND_hits_y = {fReader, "BAND_hits_y"};
   TTreeReaderArray<float> BAND_hits_z = {fReader, "BAND_hits_z"};
   TTreeReaderArray<float> BAND_hits_ux = {fReader, "BAND_hits_ux"};
   TTreeReaderArray<float> BAND_hits_uy = {fReader, "BAND_hits_uy"};
   TTreeReaderArray<float> BAND_hits_uz = {fReader, "BAND_hits_uz"};
   TTreeReaderArray<long> RUN_config_run = {fReader, "RUN_config_run"};
   TTreeReaderArray<long> RUN_config_event = {fReader, "RUN_config_event"};
   TTreeReaderArray<long> RUN_config_unixtime = {fReader, "RUN_config_unixtime"};
   TTreeReaderArray<long> RUN_config_trigger = {fReader, "RUN_config_trigger"};
   TTreeReaderArray<long> RUN_config_timestamp = {fReader, "RUN_config_timestamp"};
   TTreeReaderArray<long> RUN_config_type = {fReader, "RUN_config_type"};
   TTreeReaderArray<long> RUN_config_mode = {fReader, "RUN_config_mode"};
   TTreeReaderArray<float> RUN_config_torus = {fReader, "RUN_config_torus"};
   TTreeReaderArray<float> RUN_config_solenoid = {fReader, "RUN_config_solenoid"};
   TTreeReaderArray<long> REC_Event_NRUN = {fReader, "REC_Event_NRUN"};
   TTreeReaderArray<long> REC_Event_NEVENT = {fReader, "REC_Event_NEVENT"};
   TTreeReaderArray<long> REC_Event_TYPE = {fReader, "REC_Event_TYPE"};
   TTreeReaderArray<long> REC_Event_EvCAT = {fReader, "REC_Event_EvCAT"};
   TTreeReaderArray<long> REC_Event_NPGP = {fReader, "REC_Event_NPGP"};
   TTreeReaderArray<long> REC_Event_TRG = {fReader, "REC_Event_TRG"};
   TTreeReaderArray<long> REC_Event_Helic = {fReader, "REC_Event_Helic"};
   TTreeReaderArray<float> REC_Event_EVNTime = {fReader, "REC_Event_EVNTime"};
   TTreeReaderArray<float> REC_Event_BCG = {fReader, "REC_Event_BCG"};
   TTreeReaderArray<float> REC_Event_LT = {fReader, "REC_Event_LT"};
   TTreeReaderArray<float> REC_Event_STTime = {fReader, "REC_Event_STTime"};
   TTreeReaderArray<float> REC_Event_RFTime = {fReader, "REC_Event_RFTime"};
   TTreeReaderArray<float> REC_Event_PTIME = {fReader, "REC_Event_PTIME"};
   TTreeReaderArray<long> REC_Particle_pid = {fReader, "REC_Particle_pid"};
   TTreeReaderArray<long> REC_Particle_charge = {fReader, "REC_Particle_charge"};
   TTreeReaderArray<long> REC_Particle_status = {fReader, "REC_Particle_status"};
   TTreeReaderArray<float> REC_Particle_px = {fReader, "REC_Particle_px"};
   TTreeReaderArray<float> REC_Particle_py = {fReader, "REC_Particle_py"};
   TTreeReaderArray<float> REC_Particle_pz = {fReader, "REC_Particle_pz"};
   TTreeReaderArray<float> REC_Particle_vx = {fReader, "REC_Particle_vx"};
   TTreeReaderArray<float> REC_Particle_vy = {fReader, "REC_Particle_vy"};
   TTreeReaderArray<float> REC_Particle_vz = {fReader, "REC_Particle_vz"};
   TTreeReaderArray<float> REC_Particle_beta = {fReader, "REC_Particle_beta"};
   TTreeReaderArray<float> REC_Particle_chi2pid = {fReader, "REC_Particle_chi2pid"};
   TTreeReaderArray<long> REC_Calorimeter_index = {fReader, "REC_Calorimeter_index"};
   TTreeReaderArray<long> REC_Calorimeter_pindex = {fReader, "REC_Calorimeter_pindex"};
   TTreeReaderArray<long> REC_Calorimeter_detector = {fReader, "REC_Calorimeter_detector"};
   TTreeReaderArray<long> REC_Calorimeter_sector = {fReader, "REC_Calorimeter_sector"};
   TTreeReaderArray<long> REC_Calorimeter_layer = {fReader, "REC_Calorimeter_layer"};
   TTreeReaderArray<long> REC_Calorimeter_status = {fReader, "REC_Calorimeter_status"};
   TTreeReaderArray<float> REC_Calorimeter_energy = {fReader, "REC_Calorimeter_energy"};
   TTreeReaderArray<float> REC_Calorimeter_time = {fReader, "REC_Calorimeter_time"};
   TTreeReaderArray<float> REC_Calorimeter_path = {fReader, "REC_Calorimeter_path"};
   TTreeReaderArray<float> REC_Calorimeter_chi2 = {fReader, "REC_Calorimeter_chi2"};
   TTreeReaderArray<float> REC_Calorimeter_x = {fReader, "REC_Calorimeter_x"};
   TTreeReaderArray<float> REC_Calorimeter_y = {fReader, "REC_Calorimeter_y"};
   TTreeReaderArray<float> REC_Calorimeter_z = {fReader, "REC_Calorimeter_z"};
   TTreeReaderArray<float> REC_Calorimeter_hx = {fReader, "REC_Calorimeter_hx"};
   TTreeReaderArray<float> REC_Calorimeter_hy = {fReader, "REC_Calorimeter_hy"};
   TTreeReaderArray<float> REC_Calorimeter_hz = {fReader, "REC_Calorimeter_hz"};
   TTreeReaderArray<float> REC_Calorimeter_lu = {fReader, "REC_Calorimeter_lu"};
   TTreeReaderArray<float> REC_Calorimeter_lv = {fReader, "REC_Calorimeter_lv"};
   TTreeReaderArray<float> REC_Calorimeter_lw = {fReader, "REC_Calorimeter_lw"};
   TTreeReaderArray<float> REC_Calorimeter_du = {fReader, "REC_Calorimeter_du"};
   TTreeReaderArray<float> REC_Calorimeter_dv = {fReader, "REC_Calorimeter_dv"};
   TTreeReaderArray<float> REC_Calorimeter_dw = {fReader, "REC_Calorimeter_dw"};
   TTreeReaderArray<float> REC_Calorimeter_m2u = {fReader, "REC_Calorimeter_m2u"};
   TTreeReaderArray<float> REC_Calorimeter_m2v = {fReader, "REC_Calorimeter_m2v"};
   TTreeReaderArray<float> REC_Calorimeter_m2w = {fReader, "REC_Calorimeter_m2w"};
   TTreeReaderArray<float> REC_Calorimeter_m3u = {fReader, "REC_Calorimeter_m3u"};
   TTreeReaderArray<float> REC_Calorimeter_m3v = {fReader, "REC_Calorimeter_m3v"};
   TTreeReaderArray<float> REC_Calorimeter_m3w = {fReader, "REC_Calorimeter_m3w"};
   TTreeReaderArray<long> REC_Cherenkov_index = {fReader, "REC_Cherenkov_index"};
   TTreeReaderArray<long> REC_Cherenkov_pindex = {fReader, "REC_Cherenkov_pindex"};
   TTreeReaderArray<long> REC_Cherenkov_detector = {fReader, "REC_Cherenkov_detector"};
   TTreeReaderArray<long> REC_Cherenkov_sector = {fReader, "REC_Cherenkov_sector"};
   TTreeReaderArray<long> REC_Cherenkov_status = {fReader, "REC_Cherenkov_status"};
   TTreeReaderArray<float> REC_Cherenkov_nphe = {fReader, "REC_Cherenkov_nphe"};
   TTreeReaderArray<float> REC_Cherenkov_time = {fReader, "REC_Cherenkov_time"};
   TTreeReaderArray<float> REC_Cherenkov_path = {fReader, "REC_Cherenkov_path"};
   TTreeReaderArray<float> REC_Cherenkov_chi2 = {fReader, "REC_Cherenkov_chi2"};
   TTreeReaderArray<float> REC_Cherenkov_x = {fReader, "REC_Cherenkov_x"};
   TTreeReaderArray<float> REC_Cherenkov_y = {fReader, "REC_Cherenkov_y"};
   TTreeReaderArray<float> REC_Cherenkov_z = {fReader, "REC_Cherenkov_z"};
   TTreeReaderArray<float> REC_Cherenkov_theta = {fReader, "REC_Cherenkov_theta"};
   TTreeReaderArray<float> REC_Cherenkov_phi = {fReader, "REC_Cherenkov_phi"};
   TTreeReaderArray<float> REC_Cherenkov_dtheta = {fReader, "REC_Cherenkov_dtheta"};
   TTreeReaderArray<float> REC_Cherenkov_dphi = {fReader, "REC_Cherenkov_dphi"};
   TTreeReaderArray<long> REC_ForwardTagger_index = {fReader, "REC_ForwardTagger_index"};
   TTreeReaderArray<long> REC_ForwardTagger_pindex = {fReader, "REC_ForwardTagger_pindex"};
   TTreeReaderArray<long> REC_ForwardTagger_detector = {fReader, "REC_ForwardTagger_detector"};
   TTreeReaderArray<long> REC_ForwardTagger_size = {fReader, "REC_ForwardTagger_size"};
   TTreeReaderArray<long> REC_ForwardTagger_status = {fReader, "REC_ForwardTagger_status"};
   TTreeReaderArray<float> REC_ForwardTagger_energy = {fReader, "REC_ForwardTagger_energy"};
   TTreeReaderArray<float> REC_ForwardTagger_time = {fReader, "REC_ForwardTagger_time"};
   TTreeReaderArray<float> REC_ForwardTagger_path = {fReader, "REC_ForwardTagger_path"};
   TTreeReaderArray<float> REC_ForwardTagger_chi2 = {fReader, "REC_ForwardTagger_chi2"};
   TTreeReaderArray<float> REC_ForwardTagger_x = {fReader, "REC_ForwardTagger_x"};
   TTreeReaderArray<float> REC_ForwardTagger_y = {fReader, "REC_ForwardTagger_y"};
   TTreeReaderArray<float> REC_ForwardTagger_z = {fReader, "REC_ForwardTagger_z"};
   TTreeReaderArray<float> REC_ForwardTagger_dx = {fReader, "REC_ForwardTagger_dx"};
   TTreeReaderArray<float> REC_ForwardTagger_dy = {fReader, "REC_ForwardTagger_dy"};
   TTreeReaderArray<float> REC_ForwardTagger_radius = {fReader, "REC_ForwardTagger_radius"};
   TTreeReaderArray<long> REC_Scintillator_index = {fReader, "REC_Scintillator_index"};
   TTreeReaderArray<long> REC_Scintillator_pindex = {fReader, "REC_Scintillator_pindex"};
   TTreeReaderArray<long> REC_Scintillator_detector = {fReader, "REC_Scintillator_detector"};
   TTreeReaderArray<long> REC_Scintillator_sector = {fReader, "REC_Scintillator_sector"};
   TTreeReaderArray<long> REC_Scintillator_layer = {fReader, "REC_Scintillator_layer"};
   TTreeReaderArray<long> REC_Scintillator_component = {fReader, "REC_Scintillator_component"};
   TTreeReaderArray<long> REC_Scintillator_status = {fReader, "REC_Scintillator_status"};
   TTreeReaderArray<float> REC_Scintillator_energy = {fReader, "REC_Scintillator_energy"};
   TTreeReaderArray<float> REC_Scintillator_time = {fReader, "REC_Scintillator_time"};
   TTreeReaderArray<float> REC_Scintillator_path = {fReader, "REC_Scintillator_path"};
   TTreeReaderArray<float> REC_Scintillator_chi2 = {fReader, "REC_Scintillator_chi2"};
   TTreeReaderArray<float> REC_Scintillator_x = {fReader, "REC_Scintillator_x"};
   TTreeReaderArray<float> REC_Scintillator_y = {fReader, "REC_Scintillator_y"};
   TTreeReaderArray<float> REC_Scintillator_z = {fReader, "REC_Scintillator_z"};
   TTreeReaderArray<float> REC_Scintillator_hx = {fReader, "REC_Scintillator_hx"};
   TTreeReaderArray<float> REC_Scintillator_hy = {fReader, "REC_Scintillator_hy"};
   TTreeReaderArray<float> REC_Scintillator_hz = {fReader, "REC_Scintillator_hz"};
   TTreeReaderArray<long> REC_Track_index = {fReader, "REC_Track_index"};
   TTreeReaderArray<long> REC_Track_pindex = {fReader, "REC_Track_pindex"};
   TTreeReaderArray<long> REC_Track_detector = {fReader, "REC_Track_detector"};
   TTreeReaderArray<long> REC_Track_sector = {fReader, "REC_Track_sector"};
   TTreeReaderArray<long> REC_Track_status = {fReader, "REC_Track_status"};
   TTreeReaderArray<long> REC_Track_q = {fReader, "REC_Track_q"};
   TTreeReaderArray<long> REC_Track_NDF = {fReader, "REC_Track_NDF"};
   TTreeReaderArray<long> REC_Track_NDF_nomm = {fReader, "REC_Track_NDF_nomm"};
   TTreeReaderArray<float> REC_Track_chi2 = {fReader, "REC_Track_chi2"};
   TTreeReaderArray<float> REC_Track_px_nomm = {fReader, "REC_Track_px_nomm"};
   TTreeReaderArray<float> REC_Track_py_nomm = {fReader, "REC_Track_py_nomm"};
   TTreeReaderArray<float> REC_Track_pz_nomm = {fReader, "REC_Track_pz_nomm"};
   TTreeReaderArray<float> REC_Track_vx_nomm = {fReader, "REC_Track_vx_nomm"};
   TTreeReaderArray<float> REC_Track_vy_nomm = {fReader, "REC_Track_vy_nomm"};
   TTreeReaderArray<float> REC_Track_vz_nomm = {fReader, "REC_Track_vz_nomm"};
   TTreeReaderArray<float> REC_Track_chi2_nomm = {fReader, "REC_Track_chi2_nomm"};
   TTreeReaderArray<long> REC_CovMat_index = {fReader, "REC_CovMat_index"};
   TTreeReaderArray<long> REC_CovMat_pindex = {fReader, "REC_CovMat_pindex"};
   TTreeReaderArray<float> REC_CovMat_C11 = {fReader, "REC_CovMat_C11"};
   TTreeReaderArray<float> REC_CovMat_C12 = {fReader, "REC_CovMat_C12"};
   TTreeReaderArray<float> REC_CovMat_C13 = {fReader, "REC_CovMat_C13"};
   TTreeReaderArray<float> REC_CovMat_C14 = {fReader, "REC_CovMat_C14"};
   TTreeReaderArray<float> REC_CovMat_C15 = {fReader, "REC_CovMat_C15"};
   TTreeReaderArray<float> REC_CovMat_C22 = {fReader, "REC_CovMat_C22"};
   TTreeReaderArray<float> REC_CovMat_C23 = {fReader, "REC_CovMat_C23"};
   TTreeReaderArray<float> REC_CovMat_C24 = {fReader, "REC_CovMat_C24"};
   TTreeReaderArray<float> REC_CovMat_C25 = {fReader, "REC_CovMat_C25"};
   TTreeReaderArray<float> REC_CovMat_C33 = {fReader, "REC_CovMat_C33"};
   TTreeReaderArray<float> REC_CovMat_C34 = {fReader, "REC_CovMat_C34"};
   TTreeReaderArray<float> REC_CovMat_C35 = {fReader, "REC_CovMat_C35"};
   TTreeReaderArray<float> REC_CovMat_C44 = {fReader, "REC_CovMat_C44"};
   TTreeReaderArray<float> REC_CovMat_C45 = {fReader, "REC_CovMat_C45"};
   TTreeReaderArray<float> REC_CovMat_C55 = {fReader, "REC_CovMat_C55"};
   TTreeReaderArray<long> REC_Traj_pindex = {fReader, "REC_Traj_pindex"};
   TTreeReaderArray<long> REC_Traj_index = {fReader, "REC_Traj_index"};
   TTreeReaderArray<long> REC_Traj_detId = {fReader, "REC_Traj_detId"};
   TTreeReaderArray<long> REC_Traj_q = {fReader, "REC_Traj_q"};
   TTreeReaderArray<float> REC_Traj_x = {fReader, "REC_Traj_x"};
   TTreeReaderArray<float> REC_Traj_y = {fReader, "REC_Traj_y"};
   TTreeReaderArray<float> REC_Traj_z = {fReader, "REC_Traj_z"};
   TTreeReaderArray<float> REC_Traj_cx = {fReader, "REC_Traj_cx"};
   TTreeReaderArray<float> REC_Traj_cy = {fReader, "REC_Traj_cy"};
   TTreeReaderArray<float> REC_Traj_cz = {fReader, "REC_Traj_cz"};
   TTreeReaderArray<float> REC_Traj_pathlength = {fReader, "REC_Traj_pathlength"};

   ana_clas12(TTree * /*tree*/ =0) { }
   virtual ~ana_clas12() { }
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

   ClassDef(ana_clas12,0);
};

#endif

#ifdef ana_clas12_cxx
void ana_clas12::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

 
   fReader.SetTree(tree);

}


Bool_t ana_clas12::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef ana_clas12_cxx
