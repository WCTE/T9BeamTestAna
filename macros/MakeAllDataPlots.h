#ifndef MakeAllDataPlots_h
#define MakeAllDataPlots_h

// JK 2023--2024

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TString.h"

#include <string>
#include <vector>
#include <iostream>


#include "EventInfo.h"
#include "channelReadClass.h"

#include "tofUtil.C"

#include <string>

using namespace std;

const int nMaxChannels = 32;


class MakeAllDataPlots
{


 private:

  bool _isHodoscopeRun;
  bool _noAct1Cuts; // for ToF fits
  bool _useWindowIntCharge;
  
  Double_t peakVoltage[nMaxChannels][1];
  Double_t peakTime[nMaxChannels][1];
  Double_t signalTime[nMaxChannels][1];
  Double_t intCharge[nMaxChannels][1];
  Double_t pedestal[nMaxChannels];
  Double_t pedestalSigma[nMaxChannels];
  Double_t nPeaks[nMaxChannels];
  

  Int_t _lastSpillNumber;
  Int_t _nSpills;

  map<TString,double> _NPeaksA; 
  map<TString,double> _NPeaksC; 
  map<TString,int> _PeakIDA; // amplitudes, p.e. non-corrected
  map<TString,int> _PeakIDC; // charged base, usually p.e. corrected and window integrated
  map<TString,double> _Amplitudes; // amplitude
  map<TString,double> _Charges; // charge
  map<TString,double> _SignalTimes; // time

  // ranges
  double _tofmin;
  double _tofmax;
  double _tofmaxhigh;
  int    _ntofbins;

  double _tofminlow;
  double _tofmaxlow;
  int    _ntofbinslow;

  int _ntofbins2d;

  double _actChargeMin;
  double _actChargeMax;

  double _actAmplitudeMin;
  double _actAmplitudeMax;

  double _PbGChargeMin;
  double _PbGChargeMax;

  double _PbGAmplitudeMin;
  double _PbGAmplitudeMax;

  double _trigScintChargeMin;
  double _trigScintChargeMax;
  double _trigScintAmplitudeMin;
  double _trigScintAmplitudeMax;

  double _trigScint0_CweightedX;
  double _trigScint0_CweightedY;
  double _trigScint1_CweightedX;
  double _trigScint1_CweightedY;

  double _trigScint0_TweightedX;
  double _trigScint0_TweightedY;
  double _trigScint1_TweightedX;
  double _trigScint1_TweightedY;

  double _xL;
  double _xR;
  double _yD;
  double _yU;

  double _ACTwidth;
  double _act0x;
  double _act1x;
  double _act2x;
  double _act3x;
  
  // cuts
  map<int, map<TString,double > > _cutsMap;

  // IO
  
  string _fileName;
  int _momentum;
  TString _peakMode;
  TFile *_infile;
  TFile *_outFile;
  EventInfo *_eventInfo;

  int _nChannels;
  vector<TString> _treeNames;
  vector<int> _channelToHodoscope;
  
  // standard per channel
  vector<TH1D> _hCharge;
  vector<TH1D> _hVoltage;
  vector<TH1D> _hPedestalSigma;
  vector<TH1D> _hTime;
  vector<TH1D> _hnPeaksA;
  vector<TH1D> _hnPeaksC;
  
  map<TString,TH1D*> _histos1d;
  map<TString,TH2D*> _histos2d;

  int _Nmin;
  int _ent[nMaxChannels];
  channelReadClass *_reader[nMaxChannels];
  TTree *_trees[nMaxChannels];
  map<TString, channelReadClass*> _readerMap;

  int _debug;
  double _t0;
  double _t1;
  double _tof;

  // for peak multiplicity:
  bool _onePeakInAllACTs;
  bool _onePeakInAllToFs;
  bool _onePeakInAll;
  bool _onePeakInPbGlass;
  
  bool _moreThanOnePeakInAllACTs;
  bool _moreThanOnePeakInAllToFs;
  bool _moreThanOnePeakInAll;
  
  bool _PbGlassAboveElectronLevel;
  bool _ACT23AboveElectronLevel;

  // for tof fit:
  bool _isdACT23pb;
  bool _ispACT23pb;
  bool _isMuACT23pb;
  bool _isElACT23pb;

  tofUtil* _tofutil;


  // ACTs
  double _act0c;
  double _act1c;
  double _act2c;
  double _act3c;
  
  double _act0a;
  double _act1a;
  double _act2a;
  double _act3a;
  
  double _act23aAver;
  double _act23cAver;
  
  // hole counters and lead glass
  
  double _hc0c;
  double _hc0a;
  
  double _hc1c;
  double _hc1a;
  
  double _pbc;
  double _pba;
  
  double _trigScintA;
  double _trigScintC;
  
  double _trigScint0A;
  double _trigScint0C;
  double _trigScint0LA;
  double _trigScint0LC;
  double _trigScint0RA;
  double _trigScint0RC;
  
  double _trigScint1A;
  double _trigScint1C;
  double _trigScint1LA;
  double _trigScint1LC;
  double _trigScint1RA;
  double _trigScint1RC;
  
 public:
  
  MakeAllDataPlots(string fileName, int momentum, bool isHodoscopeRun, TString peakMode = "", bool useWindowIntCharge = false);
  ~MakeAllDataPlots();

  int getHighestPeakIndex(channelReadClass *reader, bool useCharges);
  void Init(bool noAct1Cuts);
  void InitReaders();
  void InitTofHistos();
  void InitGeneralHistos();
  void InitTrigScintHistos(TString dirname, TString selTag, TString selTit);
  void InitChargedHistos();
  void InitHodoscopeHistos();

  void ReadChannels();
  bool PassedPeakCuts();
  void FillChannels();
  void ComputeChargesAndAmplitudes();
  void FillTrigScintHistos(TString selTag);
  void FillTofHistos();
  void FillChargedHistos();
  void FillHodoscopeHistos();

  
  void Loop( int verbose = 10000, int debug = 0);
  void Terminate();

};

#endif
