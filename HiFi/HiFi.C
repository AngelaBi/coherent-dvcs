//particleTree4 /dir/file.hipo out.root Ex4_TreeMaker.C
{

  treemaker.Branch("P.Time/F");
  treemaker.Branch("PBANK.Px/F");
  treemaker.Branch("PBANK.Py/F");
  treemaker.Branch("PBANK.Pz/F");
  treemaker.Branch("PBANK.Vx/F");
  treemaker.Branch("PBANK.Vy/F");
  treemaker.Branch("PBANK.Vz/F");
  treemaker.Branch("PBANK.Pid/I");
  treemaker.Branch("P.Region/I");
  treemaker.Branch("P.Theta/F");
  //treemaker.Branch("P.Phi/F");
  //treemaker.Branch("ECIN.Energy/F");
 // treemaker.Branch("ECOUT.Energy/F");
  //treemaker.Branch("CTOF.Time/F");
  treemaker.Branch("HTCC.Nphe/I");
 // treemaker.Branch("P.P/F");


  //e.g. Only save electron information
  // treemaker.AddParticleCut("PBANK.Pid==11");

  //Event topology cuts
  // treemaker.AddAtLeastPid(211,1); //at least 1 pi+
  // treemaker.AddExactPid(11,1);    //exactly 1 electron
  // treemaker.AddZeroOfRestPid();  //nothing else, if not this line any of anything else
 
  treemaker.Fill();
}
