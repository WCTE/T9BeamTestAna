//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 18 15:46:17 2024 by ROOT version 6.30/04
// from TTree midas_data_D302/Digitizer D302
// found on file: root_run_000409.root
//////////////////////////////////////////////////////////

#ifndef midas_data_h
#define midas_data_h

#include <iostream>

using namespace std;

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class midas_data {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          eventNumber;
   UInt_t          spillNumber;
   UInt_t          timeStamp;
   UInt_t          triggerTime;
   vector<double>  *Channel0;
   vector<double>  *Channel1;
   vector<double>  *Channel2;
   vector<double>  *Channel3;
   vector<double>  *Channel4;
   vector<double>  *Channel5;
   vector<double>  *Channel6;
   vector<double>  *Channel7;

   // List of branches
   TBranch        *b_eventNumber;   //!
   TBranch        *b_spillNumber;   //!
   TBranch        *b_timeStamp;   //!
   TBranch        *b_triggerTime;   //!
   TBranch        *b_Channel0;   //!
   TBranch        *b_Channel1;   //!
   TBranch        *b_Channel2;   //!
   TBranch        *b_Channel3;   //!
   TBranch        *b_Channel4;   //!
   TBranch        *b_Channel5;   //!
   TBranch        *b_Channel6;   //!
   TBranch        *b_Channel7;   //!

   midas_data(TTree *tree=0);
   virtual ~midas_data();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString runnno);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
