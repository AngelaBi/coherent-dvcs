//particleTree4 /dir/file.hipo out.root Ex4_TreeMaker.C
{


//   treemaker.Branch("BAND");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");


//   treemaker.Branch("Config");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");


//   treemaker.Branch("Event");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");



//   treemaker.Branch("Scintillator");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");


//   treemaker.Branch("Track");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");

   
//   treemaker.Branch("CovMat");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");


//   treemaker.Branch("Traj");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");
//   treemaker.Branch("");



  treemaker.Branch("PBANK.Pid/I"); 
  treemaker.Branch("PBANK.Charge/I");
  treemaker.Branch("PBANK.Status/I");
  treemaker.Branch("PBANK.Px/F");
  treemaker.Branch("PBANK.Py/F");
  treemaker.Branch("PBANK.Pz/F");
  treemaker.Branch("PBANK.Vx/F");
  treemaker.Branch("PBANK.Vy/F");
  treemaker.Branch("PBANK.Vz/F");
////   treemaker.Branch("PBANK.Beta/F");
////   treemaker.Branch("PBANK.Chi2pid/F");
  treemaker.Branch("P.Time/F");
  treemaker.Branch("P.Theta/F");
  treemaker.Branch("P.Phi/F");
  treemaker.Branch("CTOF.Time/F");
  treemaker.Branch("P.Region/F");
  treemaker.Branch("P.P/F");

  treemaker.Branch("EVNT.RunNumber/I");
  treemaker.Branch("EVNT.EventNumber/I");
 // treemaker.Branch("EVNT.Helicity/I");
 // treemaker.Branch("HEL.helicity/I");
  treemaker.Branch("EVNT.Type/I");
  treemaker.Branch("EVNT.StartTime/F");
  treemaker.Branch("EVNT.RFTime/F");
//  treemaker.Branch("EVNT.EVNime/F");
//  treemaker.Branch("EVNT.BCG/F");
//  treemaker.Branch("EVNT.NPGP/I");
  //treemaker.Branch("EVNT.LT/F"); //not working wrong type
  treemaker.Branch("EVNT.Trigger/L");

  treemaker.Branch("CND1.Index/I");
  treemaker.Branch("CND1.Pindex/I");
  treemaker.Branch("CND1.Status/F");
  treemaker.Branch("CND1.Sector/I");
  treemaker.Branch("CND1.Layer/I");
  treemaker.Branch("CND1.Energy/F");
  treemaker.Branch("CND1.Time/F");
  treemaker.Branch("CND1.X/F");
  treemaker.Branch("CND1.Y/F");
  treemaker.Branch("CND1.Z/F");

  treemaker.Branch("ECIN.Index/I");
  treemaker.Branch("ECIN.Pindex/I");
//  treemaker.Branch("ECIN.Detector/F");
  treemaker.Branch("ECIN.Sector/I");
  treemaker.Branch("ECIN.Layer/I");
  treemaker.Branch("ECIN.Status/I");
  treemaker.Branch("ECIN.Energy/F");
  treemaker.Branch("ECIN.Time/F");
  treemaker.Branch("ECIN.Path/F");
  //   treemaker.Branch("ECIN.Chi2/F");
  treemaker.Branch("ECIN.X/F");
  treemaker.Branch("ECIN.Y/F");
  treemaker.Branch("ECIN.Z/F");
  treemaker.Branch("ECIN.Hx/F");
  treemaker.Branch("ECIN.Hy/F");
  treemaker.Branch("ECIN.Hz/F");
//  treemaker.Branch("ECIN.Lu/F");
//  treemaker.Branch("ECIN.Lv/F");
//  treemaker.Branch("ECIN.Lw/F");
  treemaker.Branch("ECIN.Du/F");
  treemaker.Branch("ECIN.Dv/F");
  treemaker.Branch("ECIN.Dw/F");
  treemaker.Branch("ECIN.M2u/F");
  treemaker.Branch("ECIN.M2v/F");
  treemaker.Branch("ECIN.M2w/F");
  treemaker.Branch("ECIN.M3u/F");
  treemaker.Branch("ECIN.M3v/F");
  treemaker.Branch("ECIN.M3w/F");


  treemaker.Branch("ECOUT.Index/F");
  treemaker.Branch("ECOUT.Pindex/F");
//  treemaker.Branch("ECOUT.Detector/F");
  treemaker.Branch("ECOUT.Sector/F");
  treemaker.Branch("ECOUT.Layer/F");
  treemaker.Branch("ECOUT.Status/F");
  treemaker.Branch("ECOUT.Energy/F");
  treemaker.Branch("ECOUT.Time/F");
  treemaker.Branch("ECOUT.Path/F");
  //   treemaker.Branch("ECOUT.Chi2/F");
  treemaker.Branch("ECOUT.X/F");
  treemaker.Branch("ECOUT.Y/F");
  treemaker.Branch("ECOUT.Z/F");
  treemaker.Branch("ECOUT.Hx/F");
  treemaker.Branch("ECOUT.Hy/F");
  treemaker.Branch("ECOUT.Hz/F");
//  treemaker.Branch("ECOUT.Lu/F");
//  treemaker.Branch("ECOUT.Lv/F");
//  treemaker.Branch("ECOUT.Lw/F");
  treemaker.Branch("ECOUT.Du/F");
  treemaker.Branch("ECOUT.Dv/F");
  treemaker.Branch("ECOUT.Dw/F");
  treemaker.Branch("ECOUT.M2u/F");
  treemaker.Branch("ECOUT.M2v/F");
  treemaker.Branch("ECOUT.M2w/F");
  treemaker.Branch("ECOUT.M3u/F");
  treemaker.Branch("ECOUT.M3v/F");
  treemaker.Branch("ECOUT.M3w/F");


  treemaker.Branch("PCAL.Index/F");
  treemaker.Branch("PCAL.Pindex/F");
//  treemaker.Branch("PCAL.Detector/F");
  treemaker.Branch("PCAL.Sector/F");
  treemaker.Branch("PCAL.Layer/F");
  treemaker.Branch("PCAL.Status/F");
  treemaker.Branch("PCAL.Energy/F");
  treemaker.Branch("PCAL.Time/F");
  treemaker.Branch("PCAL.Path/F");
  //   treemaker.Branch("PCAL.Chi2/F");
  treemaker.Branch("PCAL.X/F");
  treemaker.Branch("PCAL.Y/F");
  treemaker.Branch("PCAL.Z/F");
  treemaker.Branch("PCAL.Hx/F");
  treemaker.Branch("PCAL.Hy/F");
  treemaker.Branch("PCAL.Hz/F");
//  treemaker.Branch("PCAL.Lu/F");
//  treemaker.Branch("PCAL.Lv/F");
//  treemaker.Branch("PCAL.Lw/F");
  treemaker.Branch("PCAL.Du/F");
  treemaker.Branch("PCAL.Dv/F");
  treemaker.Branch("PCAL.Dw/F");
  treemaker.Branch("PCAL.M2u/F");
  treemaker.Branch("PCAL.M2v/F");
  treemaker.Branch("PCAL.M2w/F");
  treemaker.Branch("PCAL.M3u/F");
  treemaker.Branch("PCAL.M3v/F");
  treemaker.Branch("PCAL.M3w/F");


  treemaker.Branch("HTCC.Index/F");
  treemaker.Branch("HTCC.Pindex/F");
//  treemaker.Branch("HTCC.Detector/F");
  treemaker.Branch("HTCC.Sector/F");
  treemaker.Branch("HTCC.Status/F");
  treemaker.Branch("HTCC.Nphe/I");
  treemaker.Branch("HTCC.Time/F");
  treemaker.Branch("HTCC.Path/F");
  //   treemaker.Branch("HTCC.Chi2/F");
  treemaker.Branch("HTCC.X/F");
  treemaker.Branch("HTCC.Y/F");
  treemaker.Branch("HTCC.Z/F");
  treemaker.Branch("HTCC.Theta/F");
  treemaker.Branch("HTCC.Phi/F");
  treemaker.Branch("HTCC.Dtheta/F");
  treemaker.Branch("HTCC.DPhi/F");

  treemaker.Branch("FTCAL.Index/F");
  treemaker.Branch("FTCAL.Pindex/I");
  treemaker.Branch("FTCAL.Detector/F");
  treemaker.Branch("FTCAL.Energy/F");
  treemaker.Branch("FTCAL.Time/F");
  treemaker.Branch("FTCAL.Path/F");
  treemaker.Branch("FTCAL.Chi2/F");
  treemaker.Branch("FTCAL.X/F");
  treemaker.Branch("FTCAL.Y/F");
  treemaker.Branch("FTCAL.Z/F");
  treemaker.Branch("FTCAL.Dx/F");
  treemaker.Branch("FTCAL.Dy/F");
  treemaker.Branch("FTCAL.Radius/F");
  treemaker.Branch("FTCAL.Size/F");
  treemaker.Branch("FTCAL.Status/F");



//  //e.g. Only save electron information
  // treemaker.AddParticleCut("PBANK.Pid==11");

  //Event topology cuts
  // treemaker.AddAtLeastPid(211,1); //at least 1 pi+
  // treemaker.AddExactPid(11,1);    //exactly 1 electron
  // treemaker.AddZeroOfRestPid();  //nothing else, if not this line any of anything else
 
  treemaker.Fill();
}
