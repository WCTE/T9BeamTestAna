#define midas_data_cxx
#include "midas_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>


midas_data::midas_data(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root_run_000409.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root_run_000409.root");
      }
      f->GetObject("midas_data_D302",tree);

   }
   Init(tree);
}

midas_data::~midas_data()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t midas_data::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t midas_data::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void midas_data::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Channel0 = 0;
   Channel1 = 0;
   Channel2 = 0;
   Channel3 = 0;
   Channel4 = 0;
   Channel5 = 0;
   Channel6 = 0;
   Channel7 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("spillNumber", &spillNumber, &b_spillNumber);
   fChain->SetBranchAddress("timeStamp", &timeStamp, &b_timeStamp);
   fChain->SetBranchAddress("triggerTime", &triggerTime, &b_triggerTime);
   fChain->SetBranchAddress("Channel0", &Channel0, &b_Channel0);
   fChain->SetBranchAddress("Channel1", &Channel1, &b_Channel1);
   fChain->SetBranchAddress("Channel2", &Channel2, &b_Channel2);
   fChain->SetBranchAddress("Channel3", &Channel3, &b_Channel3);
   fChain->SetBranchAddress("Channel4", &Channel4, &b_Channel4);
   fChain->SetBranchAddress("Channel5", &Channel5, &b_Channel5);
   fChain->SetBranchAddress("Channel6", &Channel6, &b_Channel6);
   fChain->SetBranchAddress("Channel7", &Channel7, &b_Channel7);
   Notify();
}

Bool_t midas_data::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void midas_data::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t midas_data::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}



void midas_data::Loop(TString runno)
{
//   In a ROOT session, you can do:
//      root> .L midas_data.C
//      root> midas_data t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TString hname = "h2_ch2";
   TH2D *h2_ch2 = new TH2D(hname, ";Time bin;Voltage [ADC]", 100, 0, 300, 100, -100, 17000);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      double max = -999;
      double min = 99999;
      for (unsigned int i = 0; i < Channel2 -> size(); ++i) {
	double val = Channel2 -> at(i);
	//cout << i << " " << val << endl;
	if (val > max)
	  max = val;
	if (val < min)
	  min = val;
	h2_ch2 -> Fill(i,  val);
      }
      //cout << "min=" << min << " max=" << max << endl;
      

      
   }

   gStyle -> SetPalette(1);
   h2_ch2 -> SetStats(0);
   // h2_ch2 -> Scale(1./h2_ch2 -> GetEntries());
   
   TString cn = "heatMap_" + runno;
   auto can = new TCanvas(cn, cn);
   can -> cd(1);
   gPad -> SetLogz(1);
   h2_ch2 -> Draw("colz");
   auto tex = new TLatex(0.15, 0.915, TString("Run ") + runno + ", " + Form(" %1.2f M events", h2_ch2 -> GetEntries() / 1e6));
   tex -> SetNDC();
   tex -> Draw();
   
   can -> Print(TString(can -> GetName()) + TString(".png"));
   can -> Print(TString(can -> GetName()) + TString(".pdf"));
}
