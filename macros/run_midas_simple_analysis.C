// JK 18.6.2024

// example:
// root -l 'macros/run_midas_simple_analysis.C++("/scratch/WCTE/2023/root_files/root_run_000407.root")'



#include "midas_data.C"

// midas_data_D300
// midas_data_D301



void run_midas_simple_analysis(TString runno = "409", TString treename = "midas_data_D302")
{

  TString infilename = "/scratch/WCTE/2023/root_files/root_run_000" + runno + ".root";
  auto *infile = TFile::Open(infilename);
  auto tree = (TTree*) infile -> Get(treename);
  cout << infile << " " << tree << endl;
  cout << infile->GetName() << " " << tree->GetName() << endl;
  midas_data *mdana = new midas_data(tree);
  mdana -> Loop(runno);
  cout << "DONE!" << endl;
}
