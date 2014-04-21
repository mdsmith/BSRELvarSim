fileToExe = "/home/martin/Software/Simulation/BranchSiteRELnew.bf";
_inDirectoryPaths = { "0": "/home/martin/Software/Simulation/Simplex1Seq_reduced.nex"};
_fileLine = 0;
inputRedirect = {};
inputRedirect["01"]="Universal";
inputRedirect["02"]="Yes";
inputRedirect [ "03" ]   = _inDirectoryPaths[_fileLine];
inputRedirect["04"]="Y";
inputRedirect["05"]= _inDirectoryPaths[_fileLine] + ".output";  // Add this so you have the output somewhere

 ExecuteAFile ( fileToExe, inputRedirect );
