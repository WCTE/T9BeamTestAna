
#include "MakeAllDataPlots.h"
#include "Tools.C"

// jiri kvita 2023--2024

// new
using namespace std;


// ______________________________________________________________


MakeAllDataPlots::MakeAllDataPlots(string fileName, int momentum, bool isHodoscopeRun, TString peakMode, bool useWindowIntCharge ) {

  _fileName = fileName;
  _momentum = momentum;
  _peakMode = peakMode;
  _isHodoscopeRun = isHodoscopeRun;
  _useWindowIntCharge = useWindowIntCharge;
  _debug = 0;

  cout << "MakeAllDataPlots::MakeAllDataPlot: Configured as:" << endl;
  cout << " fileName: " << _fileName.c_str() << endl 
       << " momentum: " << _momentum << endl 
       << " peakMode: " << _peakMode << endl 
       << " isHodoscopeRun: " << _isHodoscopeRun << endl
       << " useWindowIntCharge: " << useWindowIntCharge << endl
       << endl;


  
}
// ______________________________________________________________


MakeAllDataPlots::~MakeAllDataPlots() {

}
// ______________________________________________________________



void MakeAllDataPlots::Init(bool noAct1Cuts)
{

  _noAct1Cuts = noAct1Cuts;

  tofmin = 10.;
  tofmax = 40.;
  ntofbins = 200;

  tofminlow = 10.;
  tofmaxlow = 20.;
  ntofbinslow = 150;

  ntofbins2d = 400;

  actChargeMin = 0.0;
  actChargeMax = 1.1; // 2.
  actAmplitudeMax =  0.8; // 2.

  PbGAmplitudeMin =  0.; // 2.
  PbGAmplitudeMax =  2.; // 2.
  
  PbGChargeMin =  0.; // 2.
  PbGChargeMax =  2.; // 2. // 1.4

  trigScintChargeMin = 0.;
  trigScintChargeMax = 18.;
  trigScintAmplitudeMin = 0.;
  trigScintAmplitudeMax = 10.;
  
  gSystem->Exec("mkdir -p histos/");

  _infile = new TFile(_fileName.c_str(), "READ");
  
  
  TString peakModeTag = "";
  if (_peakMode != "")
    peakModeTag = "_" + _peakMode;
  TString outFileName = TString(_fileName.substr(0, _fileName.size()-5).c_str()) + "_plots" + peakModeTag + ".root";
  outFileName = outFileName.ReplaceAll("output/", "histos/").ReplaceAll("ntuple_files/","histos/").ReplaceAll("windowpe_analyzed/","histos/windowpe_analyzed/");
  outFileName = outFileName.ReplaceAll("data", "histos");

  _outFile = new TFile(outFileName.Data(), "RECREATE");
  _outFile -> cd();


  _cutsMap[900] = {   { "tof_t0_cut", 4.9}, 
		      { "tof_t1_cut", 6}, 
		      { "act23_pi_minA", 0.4}, 
		      { "act23_pi_maxA", 2.4}, 
		      { "pb_min", 0.1}, 
		      { "actThresh", 0.5},
  };

  cout << "Initialized with: " << endl;
  cout << " _noAct1Cuts: " << _noAct1Cuts << endl;
    
}


// ______________________________________________________________


void MakeAllDataPlots::InitReaders()
{
  cout << "In InitReaders" << endl;
  _eventInfo = new EventInfo(_infile, "EventInfo");
  cout << "Initializing the tree readers..." << endl;
  for (int ich = 0; ich < _nChannels; ++ich) {
    //cout << "Initializing " <<  _treeNames[ich] << endl;
    _reader[ich] = new channelReadClass(_infile, _treeNames[ich], _useWindowIntCharge);
    _trees[ich] = _reader[ich] -> fChain;
    _ent[ich] = _trees[ich] -> GetEntries();
    _readerMap[_treeNames[ich]] = _reader[ich];
    cout << "  ...reader for " << _treeNames[ich].Data() << ": " << _ent[ich] << " entries." << endl;
  }
  cout << "done Init" << endl;

  // jk 16.11.2023
  _Nmin = 999999999;
  for (int ich = 0; ich < _nChannels; ++ich) {
    if (_ent[ich] < _Nmin)
      _Nmin = _ent[ich];
  }
  cout << "Minimal entries over trees: " << _Nmin << endl;

  
}

// ______________________________________________________________

void MakeAllDataPlots::InitGeneralHistos() {
  
  _nChannels = 32;
  cout << "InitGeneralHistos" << endl;
  
  _outFile -> mkdir("General");
  _outFile -> cd("General");
  
 for(int i = 0; i < _nChannels; i++) {
    string name1 = "hRef_Charge" + to_string(i);
    string name2 = "hRef_Voltage" + to_string(i);
    string name3 = "hRef_Hits" + to_string(i);
    string name4 = "hRef_PedestalSigma" + to_string(i);
    string name5 = "hRef_Time" + to_string(i);
    string name6 = "hRef_nPeaks" + to_string(i);
    string name7 = "hRef_Pedestal" + to_string(i);
    string name8 = "hRef_PedestalNbPeaks" + to_string(i);

    string title1 = "Channel " + to_string(i) + "; Charge [nC]; Triggers";
    string title2 = "Channel " + to_string(i) + "; Total Amplitude [V]; Triggers";
    string title3 = "Channel " + to_string(i) + "; Hits per trigger; Triggers";
    string title4 = "Channel " + to_string(i) + "; #sigma_{ped} [V]; Triggers";
    string title5 = "Channel " + to_string(i) + "; Time [ns]; Triggers";
    string title6 = "Channel " + to_string(i) + "; Number of peaks; Triggers";
    string title7 = "Channel " + to_string(i) + "; Pedestal Amplitude; Triggers/1mV";
    string title8 = "Channel " + to_string(i) + "; Pedestal Amplitude; Number of peaks";

    TH1D temp1(name1.c_str(), title1.c_str(), 400, actChargeMin, actChargeMax);
    TH1D temp2(name2.c_str(), title2.c_str(), 400, 0., 2);
    TH1D temp3(name3.c_str(), title3.c_str(), 5, -0.5, 4.5);
    TH1D temp4(name4.c_str(), title4.c_str(), 200, 0., 0.01);
    TH1D temp5(name5.c_str(), title5.c_str(), 270, 0., 540.);
    TH1D temp6(name6.c_str(), title6.c_str(), 20, 0., 20.);
    TH1D temp7(name7.c_str(), title7.c_str(), 1000., 1.65, 1.65+1000*0.0012207);
    TH2D temp8(name8.c_str(), title8.c_str(), 200, 0., 3*0.8,  5, 0., 5.);

    hCharge.push_back(temp1);
    hVoltage.push_back(temp2);
    hPedestalSigma.push_back(temp4);
    hTime.push_back(temp5);
    hnPeaks.push_back(temp6);
  }
 _outFile -> cd("../");
 cout << "done" << endl;
 
}

// ______________________________________________________________

void MakeAllDataPlots::InitHodoscopeHistos() {
  cout << "In initHodoscopeHistos" << endl;
  _channelToHodoscope = {8, 9, 10, 11, 12, 13, 14, 0, 1, 2, 3, 4, 5, 6, 7};
  _nChannels = 31;
  // the order here is not directly translatable to channel number
  // due to the missing channel 7 in digi0
  _treeNames = {
    "ACT0L",    "ACT0R",
    "ACT1L",    "ACT1R",
    "ACT3L", 	"ACT3R",
    "TriggerScint",
    "TOF00", 	"TOF01", 	"TOF02", 	"TOF03",
    "TOF10", 	"TOF11", 	"TOF12", 	"TOF13",
    "PbGlass",
    "HD8",   "HD9",   "HD10",   "HD11",   "HD12", "HD13",   "HD14",
    "HD0",   "HD1",   "HD2",  "HD3",  "HD4",  "HD5",  "HD6",  "HD7"
  };
  _outFile -> mkdir("Hodoscope");
  _outFile -> cd("Hodoscope");
  
  // lead glass A vs hodoscope occupancy with some amplitude cuts
  // subject to mV callibration!!
  // jiri on shift 28.7.2023
  TString name = "LeadGlassPhotonAVsPositronHodoOcc";
  TString title = name + ";HD Channel ID;A^{#gamma}_{Pb}";
  _histos2d[name] = new TH2D(name, title, 15, 0, 15, 125, 0, 0.12);

  name = "LeadGlassPhotonAVsPositronMaxHodoOcc";
  title = name + ";HD Max. Channel ID;A^{#gamma}_{Pb}";
  _histos2d[name] = new TH2D(name, title, 15, 0, 15, 125, 0, 0.12);

  name = "HodoOccScatter";
  title = name + ";HD channel;HD channel;entries";
  _histos2d[name] = new TH2D(name, title, 15, 0, 15, 15, 0, 15);
   
  name = "HodoOccScatterFrac";
  title = name + ";HD channel;HD channel;fractions";
  _histos2d[name] = new TH2D(name, title, 15, 0, 15, 15, 0, 15);
  _histos1d["hnHitsHodoscope"] = new TH1D("hnHitsHodoscope", ";Hodoscope channel; Number of hits", 15, 0, 15);
  
  _outFile -> cd("../");
  cout << "done" << endl;

}


// ______________________________________________________________

void MakeAllDataPlots::InitTofHistos()
{
  cout << "In InitTofHistos" << endl;

  _outFile -> mkdir("TOF");
  _outFile -> cd("TOF");
  // TOF 1D
  _histos1d["hTOFAll"] = new TH1D("hTOFAll", ";t_{TOF}^{All} [ns]", 120, tofmin, tofmax);
  _histos1d["hTOFAllWide"] = new TH1D("hTOFAllWide", ";t_{TOF}^{All} [ns]", 2*ntofbins, tofmin, 2*tofmax);
  _histos1d["hTOFEl"] = new TH1D("hTOFEl", ";t_{TOF}^{e} [ns]", ntofbins, tofmin, tofmax);
  _histos1d["hTOFOther"] = new TH1D("hTOFOther", ";t_{TOF}^{non-e} [ns]", ntofbins, tofmin, tofmax);
  _histos1d["hTOFOther_act1cuts"] = new TH1D("hTOFOther_act1cuts", ";t_{TOF}^{non-e} [ns]", ntofbins, tofmin, tofmax);
  
  
  _histos1d["hTOFAllLow"] = new TH1D("hTOFAllLow", ";t_{TOF}^{All} [ns]", ntofbinslow, tofminlow, tofmaxlow);
  _histos1d["hTOFElLow"] = new TH1D("hTOFElLow", ";t_{TOF}^{e} [ns]", ntofbinslow, tofminlow, tofmaxlow);
  _histos1d["hTOFOtherLow"] = new TH1D("hTOFOtherLow", ";t_{TOF}^{non-e} [ns]", ntofbinslow, tofminlow, tofmaxlow);
  _histos1d["hTOFOtherLow_act1cuts"] = new TH1D("hTOFOtherLow_act1cuts", ";t_{TOF}^{non-e} [ns]", ntofbinslow, tofminlow, tofmaxlow);

  _histos1d["hT0"] = new TH1D("hRef_T0", "", 270, 50, 320);
  _histos1d["hT1"] = new TH1D("hRef_T1", "", 270, 50, 320);

  // jiri
  _histos1d["hTimeReso0"] = new TH1D("hTimeReso0", "", 200, -100, 100);
  _histos1d["hTimeReso1"] = new TH1D("hTimeReso1", "", 200, -100, 100);
  _histos1d["hTimeReso0_zoom"] = new TH1D("hTimeReso0_zoom", "", 160, 20, 30);
  _histos1d["hTimeReso1_zoom"] = new TH1D("hTimeReso1_zoom", "", 160, -5, 5);

  // 2023 time offset analysis
  _histos1d["hTimeDiffTOF01"] = new TH1D("hTimeDiffTOF01", "hTimeDiffTOF01", 100, -12.,12.);
  _histos1d["hTimeDiffTOF02"] = new TH1D("hTimeDiffTOF02", "hTimeDiffTOF02", 100, -12.,12.);
  _histos1d["hTimeDiffTOF03"] = new TH1D("hTimeDiffTOF03", "hTimeDiffTOF03", 100, -12.,12.);

  _histos1d["hTimeDiffTOF11"] = new TH1D("hTimeDiffTOF11", "hTimeDiffTOF11", 100, -12.,12.);
  _histos1d["hTimeDiffTOF12"] = new TH1D("hTimeDiffTOF12", "hTimeDiffTOF12", 100, -12.,12.);
  _histos1d["hTimeDiffTOF13"] = new TH1D("hTimeDiffTOF13", "hTimeDiffTOF13", 100, -12.,12.);


  //acraplet TOF analysis
  _histos1d["hTimeTOF0"] = new TH1D("hTimeTOF0", "; hTimeTOF0", 100, 0.,50.);
  _histos1d["hTimeTOF1"] = new TH1D("hTimeTOF1", "; hTimeTOF1", 100, 0.,50.);
  _histos1d["hTimeTOF2"] = new TH1D("hTimeTOF2", "; hTimeTOF2", 100, 0.,50.);
  _histos1d["hTimeTOF3"] = new TH1D("hTimeTOF3", "; hTimeTOF3", 100, 0.,50.);

  _outFile -> cd("../");

}

// ______________________________________________________________

void MakeAllDataPlots::InitChargedHistos()
{

  cout << "InitChargedHistos" << endl;
  _nChannels = 19; 

  // https://docs.google.com/spreadsheets/d/1QBHKEbpC_roTHyY5QJSExFnMmDGn6aLyWtHHhrKORZA/edit?usp=sharing

  
  _treeNames = {
    "ACT0L",    "ACT0R",
    "ACT1L",    "ACT1R",
    "ACT2L", 	"ACT2R",
    "ACT3L", 	"ACT3R",
    "TOF00", 	"TOF01", 	"TOF02", 	"TOF03",
    "TOF10", 	"TOF11", 	"TOF12", 	"TOF13",
    "Hole0", 	"Hole1", 	"PbGlass"
  };

  double ACT0Gain = 5.;

  _outFile -> mkdir("Charged");
  _outFile -> cd("Charged");
  
  //lead glass vs act 2 and 3 - identify particles
  _histos2d["hRef_pbA_act23A"] = new TH2D("hRef_pbA_act23A", "; Pb-glass Amplitude ; (ACT2+ACT3)/2 Amplitude", 200, 0., PbGAmplitudeMax, 400, 0., 2*actAmplitudeMax);
  _histos2d["hRef_pbC_act23C"] = new TH2D("hRef_pbC_act23C", "; Pb-glass Charge ; (ACT2+ACT3)/2 Charge", 200, actChargeMin, actChargeMax, 400, 0., 2*actAmplitudeMax);

  _histos2d["hRef_pbC_act23A"] = new TH2D("hRef_pbC_act23A", "; Pb-glass Charge ; (ACT2+ACT3)/2 Amplitude", 200, 0., actAmplitudeMax, 400, 0., 2*actAmplitudeMax);
  _histos2d["hRef_pbA_act23C"] = new TH2D("hRef_pbA_act23C", "; Pb-glass Charge ; (ACT2+ACT3)/2 Amplitude", 200, PbGChargeMin, PbGChargeMax, 400, 0., 2*actAmplitudeMax);

  _histos2d["hRef_pbA_act0A"] = new TH2D("hRef_pbA_act0A", "; Pb-glass Amplitude ; ACT0 Amplitude", 200, 0., PbGAmplitudeMax, 400, 0., ACT0Gain*actAmplitudeMax);
  _histos2d["hRef_pbC_act0C"] = new TH2D("hRef_pbC_act0C", "; Pb-glass Charge ; ACT1 Charge", 200, PbGChargeMin, PbGChargeMax, 400, 0., ACT0Gain*actChargeMax);
  _histos2d["hRef_pbA_act1A"] = new TH2D("hRef_pbA_act1A", "; Pb-glass Amplitude ; ACT1 Amplitude", 200, 0., PbGAmplitudeMax, 400, 0., actAmplitudeMax);
  _histos2d["hRef_pbC_PbG1C"] = new TH2D("hRef_pbC_PbG1C", "; Pb-glass Charge ; PBG1 Charge", 200, actChargeMin, PbGChargeMax, 400, 0., PbGChargeMax);
  _histos2d["hRef_pbA_act1C"] = new TH2D("hRef_pbA_act1C", "; Pb-glass Amplitude ; ACT1 Charge", 200, 0., PbGAmplitudeMax, 400,actChargeMin, actChargeMax);

  _histos2d["hRef_pbC_act1C"] = new TH2D("hRef_pbC_act1C", "; Pb-glass Charge ; ACT1 Charge", 200, actChargeMin, actChargeMax, 400, 0., actAmplitudeMax);
  // act1cuts
  _histos2d["hRef_pbA_act1A_act1cuts"] = new TH2D("hRef_pbA_act1A_act1cuts", "; Pb-glass Amplitude ; ACT1 Amplitude", 200, 0., PbGAmplitudeMax, 400, 0., actAmplitudeMax);


  // 5.3.2024
  // Trigger scintillators (unfortunatelly labelled as TOF all through out the code;-)
  // amplitude and charge vs tof:

  _outFile -> cd("../");

  _outFile -> mkdir("TrigScint");
  _outFile -> cd("TrigScint");
  // both trig scintil.:
  _histos2d["hRef_pbC_TrigScintC"] = new TH2D("hRef_pbC_TrigScintC", "; Pb-glass Charge ; Trig. scint. Charge", 200, PbGChargeMin, PbGChargeMax, 400, trigScintChargeMin, trigScintChargeMax);
  _histos2d["hRef_pbA_TrigScintC"] = new TH2D("hRef_pbA_TrigScintC", "; Pb-glass Amplitude ; Trig. scint. Charge", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintChargeMax);
  _histos2d["hRef_pbC_TrigScintA"] = new TH2D("hRef_pbC_TrigScintA", "; Pb-glass Charge ; Trig. scint. Amplitude", 200, PbGChargeMin, PbGChargeMax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax);
  _histos2d["hRef_pbA_TrigScintA"] = new TH2D("hRef_pbA_TrigScintA", "; Pb-glass Amplitude ; Trig. scint. Amplitude", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintAmplitudeMax);

  _histos2d["hRef_TOF_TrigScintC"] = new TH2D("hRef_TOF_TrigScintC", "; t_{1}-t_{0} [ns] ; Trig. scint. Charge", ntofbins2d, tofmin, tofmax, 400, trigScintChargeMin, trigScintChargeMax);
  _histos2d["hRef_TOF_TrigScintA"] = new TH2D("hRef_TOF_TrigScintA", "; t_{1}-t_{0} [ns] ; Trig. scint. Amplitude", ntofbins2d, tofmin, tofmax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax);

  // TOF0X
  _histos2d["hRef_pbC_TrigScint0C"] = new TH2D("hRef_pbC_TrigScint0C", "; Pb-glass Charge ; Trig. scint. 0 Charge", 200, PbGChargeMin, PbGChargeMax, 400, trigScintChargeMin, trigScintChargeMax/2.);
  _histos2d["hRef_pbA_TrigScint0C"] = new TH2D("hRef_pbA_TrigScint0C", "; Pb-glass Amplitude ; Trig. scint. 0 Charge", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintChargeMax/2.);
  _histos2d["hRef_pbC_TrigScint0A"] = new TH2D("hRef_pbC_TrigScint0A", "; Pb-glass Charge ; Trig. scint. 0 Amplitude", 200, PbGChargeMin, PbGChargeMax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/2.);
  _histos2d["hRef_pbA_TrigScint0A"] = new TH2D("hRef_pbA_TrigScint0A", "; Pb-glass Amplitude ; Trig. scint. 0 Amplitude", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintAmplitudeMax/2.);

  _histos2d["hRef_TOF_TrigScint0C"] = new TH2D("hRef_TOF_TrigScint0C", "; t_{1}-t_{0} [ns] ; Trig. scint. 0 Charge", ntofbins2d, tofmin, tofmax, 400, trigScintChargeMin, trigScintChargeMax/2.);
  _histos2d["hRef_TOF_TrigScint0A"] = new TH2D("hRef_TOF_TrigScint0A", "; t_{1}-t_{0} [ns] ; Trig. scint. 0 Amplitude", ntofbins2d, tofmin, tofmax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/2.);

  // TOF0 0+1
  _histos2d["hRef_pbC_TrigScint001C"] = new TH2D("hRef_pbC_TrigScint001C", "; Pb-glass Charge ; Trig. scint. 0 0+1 Charge", 200, PbGChargeMin, PbGChargeMax, 400, trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_pbA_TrigScint001C"] = new TH2D("hRef_pbA_TrigScint001C", "; Pb-glass Amplitude ; Trig. scint. 0 0+1 Charge", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_pbC_TrigScint001A"] = new TH2D("hRef_pbC_TrigScint001A", "; Pb-glass Charge ; Trig. scint. 0 0+1 Amplitude", 200, PbGChargeMin, PbGChargeMax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/4.);
  _histos2d["hRef_pbA_TrigScint001A"] = new TH2D("hRef_pbA_TrigScint001A", "; Pb-glass Amplitude ; Trig. scint. 0 0+1 Amplitude", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintAmplitudeMax/4.);

  _histos2d["hRef_TOF_TrigScint001C"] = new TH2D("hRef_TOF_TrigScint001C", "; t_{1}-t_{0} [ns] ; Trig. scint. 0 0+1 Charge", ntofbins2d, tofmin, tofmax, 400, trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_TOF_TrigScint001A"] = new TH2D("hRef_TOF_TrigScint001A", "; t_{1}-t_{0} [ns] ; Trig. scint. 0 0+1 Amplitude", ntofbins2d, tofmin, tofmax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/4.);

   // TOF0 2+3
  _histos2d["hRef_pbC_TrigScint023C"] = new TH2D("hRef_pbC_TrigScint023C", "; Pb-glass Charge ; Trig. scint. 0 2+3 Charge", 200, PbGChargeMin, PbGChargeMax, 400, trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_pbA_TrigScint023C"] = new TH2D("hRef_pbA_TrigScint023C", "; Pb-glass Amplitude ; Trig. scint. 0 2+3 Charge", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_pbC_TrigScint023A"] = new TH2D("hRef_pbC_TrigScint023A", "; Pb-glass Charge ; Trig. scint. 0 2+3 Amplitude", 200, PbGChargeMin, PbGChargeMax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/4.);
  _histos2d["hRef_pbA_TrigScint023A"] = new TH2D("hRef_pbA_TrigScint023A", "; Pb-glass Amplitude ; Trig. scint. 0 2+3 Amplitude", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintAmplitudeMax/4.);

  _histos2d["hRef_TOF_TrigScint023C"] = new TH2D("hRef_TOF_TrigScint023C", "; t_{1}-t_{0} [ns] ; Trig. scint. 0 2+3 Charge", ntofbins2d, tofmin, tofmax, 400, trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_TOF_TrigScint023A"] = new TH2D("hRef_TOF_TrigScint023A", "; t_{1}-t_{0} [ns] ; Trig. scint. 0 2+3 Amplitude", ntofbins2d, tofmin, tofmax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/4.);

  // TOF1X
  _histos2d["hRef_pbC_TrigScint1C"] = new TH2D("hRef_pbC_TrigScint1C", "; Pb-glass Charge ; Trig. scint. 1 Charge", 200, PbGChargeMin, PbGChargeMax, 400, trigScintChargeMin, trigScintChargeMax/2.);
  _histos2d["hRef_pbA_TrigScint1C"] = new TH2D("hRef_pbA_TrigScint1C", "; Pb-glass Amplitude ; Trig. scint. 1 Charge", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintChargeMax/2.);
  _histos2d["hRef_pbC_TrigScint1A"] = new TH2D("hRef_pbC_TrigScint1A", "; Pb-glass Charge ; Trig. scint. 1 Amplitude", 200, PbGChargeMin, PbGChargeMax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/2.);
  _histos2d["hRef_pbA_TrigScint1A"] = new TH2D("hRef_pbA_TrigScint1A", "; Pb-glass Amplitude ; Trig. scint. 1 Amplitude", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintAmplitudeMax/2.);

  _histos2d["hRef_TOF_TrigScint1C"] = new TH2D("hRef_TOF_TrigScint1C", "; t_{1}-t_{0} [ns] ; Trig. scint. 1 Charge", ntofbins2d, tofmin, tofmax, 400, trigScintChargeMin, trigScintChargeMax/2.);
  _histos2d["hRef_TOF_TrigScint1A"] = new TH2D("hRef_TOF_TrigScint1A", "; t_{1}-t_{0} [ns] ; Trig. scint. 1 Amplitude", ntofbins2d, tofmin, tofmax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/2.);
  
  // TOF1 0+1
  _histos2d["hRef_pbC_TrigScint101C"] = new TH2D("hRef_pbC_TrigScint101C", "; Pb-glass Charge ; Trig. scint.1 0+1 Charge", 200, PbGChargeMin, PbGChargeMax, 400, trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_pbA_TrigScint101C"] = new TH2D("hRef_pbA_TrigScint101C", "; Pb-glass Amplitude ; Trig. scint.1 0+1 Charge", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_pbC_TrigScint101A"] = new TH2D("hRef_pbC_TrigScint101A", "; Pb-glass Charge ; Trig. scint.1 0+1 Amplitude", 200, PbGChargeMin, PbGChargeMax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/4.);
  _histos2d["hRef_pbA_TrigScint101A"] = new TH2D("hRef_pbA_TrigScint101A", "; Pb-glass Amplitude ; Trig. scint.1 0+1 Amplitude", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintAmplitudeMax/4.);

  _histos2d["hRef_TOF_TrigScint101C"] = new TH2D("hRef_TOF_TrigScint101C", "; t_{1}-t_{0} [ns] ; Trig. scint.1 0+1 Charge", ntofbins2d, tofmin, tofmax, 400, trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_TOF_TrigScint101A"] = new TH2D("hRef_TOF_TrigScint101A", "; t_{1}-t_{0} [ns] ; Trig. scint.1 0+1 Amplitude", ntofbins2d, tofmin, tofmax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/4.);

  // TOF1 2+3
  _histos2d["hRef_pbC_TrigScint123C"] = new TH2D("hRef_pbC_TrigScint123C", "; Pb-glass Charge ; Trig. scint.1 2+3 Charge", 200, PbGChargeMin, PbGChargeMax, 400, trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_pbA_TrigScint123C"] = new TH2D("hRef_pbA_TrigScint123C", "; Pb-glass Amplitude ; Trig. scint.1 2+3 Charge", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_pbC_TrigScint123A"] = new TH2D("hRef_pbC_TrigScint123A", "; Pb-glass Charge ; Trig. scint.1 2+3 Amplitude", 200, PbGChargeMin, PbGChargeMax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/4.);
  _histos2d["hRef_pbA_TrigScint123A"] = new TH2D("hRef_pbA_TrigScint123A", "; Pb-glass Amplitude ; Trig. scint.1 2+3 Amplitude", 200, 0., PbGAmplitudeMax, 400,trigScintChargeMin, trigScintAmplitudeMax/4.);

  _histos2d["hRef_TOF_TrigScint123C"] = new TH2D("hRef_TOF_TrigScint123C", "; t_{1}-t_{0} [ns] ; Trig. scint.1 2+3 Charge", ntofbins2d, tofmin, tofmax, 400, trigScintChargeMin, trigScintChargeMax/4.);
  _histos2d["hRef_TOF_TrigScint123A"] = new TH2D("hRef_TOF_TrigScint123A", "; t_{1}-t_{0} [ns] ; Trig. scint.1 2+3 Amplitude", ntofbins2d, tofmin, tofmax, 400, trigScintAmplitudeMin, trigScintAmplitudeMax/4.);
  
  _outFile -> cd("../");


  _outFile -> cd("Charged");

  // L-R studies:
  // 17.11.2023
  _histos2d["hRef_act0LA_act0RA_noneZero"] = new TH2D("hRef_act0LA_act0RA_nonZero", "; ACT0L Amplitude;ACT0R Amplitude", 200, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);
  _histos2d["hRef_act1LA_act1RA_noneZero"] = new TH2D("hRef_act1LA_act1RA_nonZero", "; ACT1L Amplitude;ACT1R Amplitude", 200, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);
  _histos2d["hRef_act2LA_act2RA_noneZero"] = new TH2D("hRef_act2LA_act2RA_nonZero", "; ACT2L Amplitude;ACT2R Amplitude", 200, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);
  _histos2d["hRef_act3LA_act3RA_noneZero"] = new TH2D("hRef_act3LA_act3RA_nonZero", "; ACT3L Amplitude;ACT3R Amplitude", 200, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);

  _histos2d["hRef_act0LA_act1LA_noneZero"] = new TH2D("hRef_act0LA_act1LA_nonZero", "; ACT0L Amplitude;ACT1L Amplitude", 200, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);
  _histos2d["hRef_act0RA_act1RA_noneZero"] = new TH2D("hRef_act0RA_act1RA_nonZero", "; ACT0R Amplitude;ACT1R Amplitude", 200, 0., actAmplitudeMax, 400, 0., actAmplitudeMax);


  // (ACT2+ACT3)/2 vs TOF plots
  _histos2d["hRef_TOFACT23A"] = new TH2D("hRef_TOFACT23A", "; t_{1}-t_{0} [ns]; (ACT2+ACT3)/2 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., 2*actAmplitudeMax);
  _histos2d["hRef_TOFACT23C"] = new TH2D("hRef_TOFACT23C", "; t_{1}-t_{0} [ns]; (ACT2+ACT3)/2 Charge", ntofbins2d, tofmin, tofmax, 200, actChargeMin, 2*actChargeMax);

  // also ACT 0 and 1, separately:
  /* seems they were already defined below...
  _histos2d["hRef_TOFACT0A"] = new TH2D("hRef_TOFACT0A", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  _histos2d["hRef_TOFACT1A"] = new TH2D("hRef_TOFACT1A", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  _histos2d["hRef_TOFACT0C"] = new TH2D("hRef_TOFACT0C", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, actChargeMin, actChargeMax);
  _histos2d["hRef_TOFACT1C"] = new TH2D("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, actChargeMin, actChargeMax);
  */

  //TOF vs Pb-glass plots
  _histos2d["hRef_PbATOF"] = new TH2D("hRef_PbATOF", "; Pb-glass Amplitude; t_{1}-t_{0} [ns]", 200, 0., PbGAmplitudeMax, ntofbins2d, tofmin, tofmax);
  _histos2d["hRef_PbCTOF"] = new TH2D("hRef_PbCTOF", "; Pb-glass Charge; t_{1}-t_{0} [ns]", 200, PbGChargeMin, PbGChargeMax, ntofbins2d, tofmin, tofmax);
  _histos2d["hRef_TOFPbA"] = new TH2D("hRef_TOFPbA", "; t_{1}-t_{0} [ns]; Pb-glass Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., PbGAmplitudeMax);
  _histos2d["hRef_TOFPbC"] = new TH2D("hRef_TOFPbC", "; t_{1}-t_{0} [ns]; Pb-glass Charge", ntofbins2d, tofmin, tofmax, 200, 0., PbGChargeMax);
  
  //acraplet - investigate "weird electrons"
  _histos2d["hHC0AHC1A"] = new TH2D("hweirdE_HC0AHC1A", "; Hole Counter 0 Amplitude; Hole Counter 1 Amplitude", 200, 0., 1000, 200, 0., 1000.);
  _histos2d["hHC0CHC1C"] = new TH2D("hweirdE_HC0CHC1C", "; Hole Counter 0 Charge; Hole Counter 1 Charge", 200, 0., 1., 200, 0., 1.);

  // no cuts
  _histos2d["hRef_TOFACT0A"] = new TH2D("hRef_TOFACT0A", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  _histos2d["hRef_TOFACT1A"] = new TH2D("hRef_TOFACT1A", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  _histos2d["hRef_TOFACT2A"] = new TH2D("hRef_TOFACT2A", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);
  _histos2d["hRef_TOFACT3A"] = new TH2D("hRef_TOFACT3A", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", ntofbins2d, tofmin, tofmax, 200, 0., actAmplitudeMax);

  _histos2d["hRef_TOFACT0C"] = new TH2D("hRef_TOFACT0C", "; t_{1}-t_{0} [ns]; ACT0 Charge", ntofbins2d, tofmin, tofmax, 200, actChargeMin, actChargeMax);
  _histos2d["hRef_TOFACT1C"] = new TH2D("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", ntofbins2d, tofmin, tofmax, 200, actChargeMin, actChargeMax);
  _histos2d["hRef_TOFACT2C"] = new TH2D("hRef_TOFACT2C", "; t_{1}-t_{0} [ns]; ACT2 Charge", ntofbins2d, tofmin, tofmax, 200, actChargeMin, actChargeMax);
  _histos2d["hRef_TOFACT3C"] = new TH2D("hRef_TOFACT3C", "; t_{1}-t_{0} [ns]; ACT3 Charge", ntofbins2d, tofmin, tofmax, 200, actChargeMin, actChargeMax);


  // ACT2+ACT3 cut
  _histos1d["hTOF_act2act3cut"] = new TH1D("hTOF_act2act3cut", "; t_{1}-t_{0} [ns];", 120, tofmin, tofmax);

  // 2D ACT charges
  _histos2d["hACT2CACT1C"] = new TH2D("hRef_ACT2CACT1C", "; ACT2 Charge; ACT1 Charge", 200, actChargeMin, actChargeMax, 200, actChargeMin, actChargeMax);
  _histos2d["hACT3CACT2C"] = new TH2D("hRef_ACT3CACT2C", "; ACT3 Charge; ACT2 Charge", 200, actChargeMin, actChargeMax, 200, actChargeMin, actChargeMax);
  _histos2d["hACT1CACT3C"] = new TH2D("hRef_ACT1CACT3C", "; ACT1 Charge; ACT3 Charge", 200, actChargeMin, actChargeMax, 200, actChargeMin, actChargeMax);

  _outFile -> cd("../");

  // nPeak 2D plots;)
  int nbn = 16.;
  double n1 = 0.;
  double n2 = 4.;

  _outFile -> mkdir("nPeaks");
  _outFile -> cd("nPeaks");

  _histos2d["hnPeaksACT23vsnPeaksToF"] = new TH2D("hnPeaksACT23vsnPeaksToF", "hnPeaksACT23vsnPeaksToF;<n_{Peaks}^{ToF}>;<n_{Peaks}^{ACT23}>", nbn, n1, n2, nbn, n1, n2);
  _histos2d["hnPeaksToF1vsnPeaksToF0"] = new TH2D("hnPeaksToF1vsnPeaksToF0", "hnPeaksToF1vsnPeaksToF0;<n_{Peaks}^{ToF0}>;<n_{Peaks}^{ToF1}>", nbn, n1, n2, nbn, n1, n2);
  _histos2d["hnPeaksACT3vsnPeaksACT2"] = new TH2D("hnPeaksACT3vsnPeaksACT2", "hnPeaksACT3vsnPeaksACT2;<n_{Peaks}^{ACT2}>;<n_{Peaks}^{ACT3}>", nbn/2, n1, n2, nbn/2, n1, n2);
  
  _histos2d["hnPeaksACT23vsToF"] = new TH2D("hnPeaksACT23vsToF", "hnPeaksACT23vsToF;t_{TOF};<n_{Peaks}^{ACT23}>", ntofbins2d/4, tofmin, tofmax, nbn, n1, n2);
  _histos2d["hnPeaksACT23vsToFlow"] = new TH2D("hnPeaksACT23vsToFlow", "hnPeaksACT23vsToF;t_{TOF};<n_{Peaks}^{ACT23}>", ntofbins2d/4, tofminlow, tofmaxlow, nbn, n1, n2);
  _histos2d["hnPeaksToFvsToF"] = new TH2D("hnPeaksToFvsToF", "hnPeaksToFvsToF;t_{TOF};<n_{Peaks}^{ToF}>", ntofbins2d/4, tofmin, tofmax, nbn, n1, n2);
  _histos2d["hnPeaksToFvsToFlow"] = new TH2D("hnPeaksToFvsToFlow", "hnPeaksToFvsToF;t_{TOF};<n_{Peaks}^{ToF}>", ntofbins2d/4, tofminlow, tofmaxlow, nbn, n1, n2);
  _histos2d["hnPeaksACT23vsLeadGlassA"] = new TH2D("hnPeaksACT23vsLeadGlassA", "hnPeaksACT23vsLeadGlassA;lead glass A;<n_{Peaks}^{ACT23}>", 100,  0., actAmplitudeMax/2., nbn, n1, n2);
  _histos2d["hnPeaksToFvsLeadGlassA"] = new TH2D("hnPeaksToFvsLeadGlassA", "hnPeaksToFvsLeadGlassA;lead glass A;<n_{Peaks}^{ToF}>", 100,  0., actAmplitudeMax/2., nbn, n1, n2);
  n1 = 0.;
  n2 = 10.;
  _histos2d["hnPeaksLeadGlassvsLeadGlassA"] = new TH2D("hnPeaksLeadGlassvsLeadGlassA", "hnPeaksLeadGlassvsLeadGlassA;lead glass A;n_{Peaks}^{Pb}", 100,  0., actAmplitudeMax/2., int(n2-n1), n1, n2);
  _outFile -> cd("../");
  
  cout << "done" << endl;
  
    
} // InitHistos



// ______________________________________________________________

int MakeAllDataPlots::getHighestPeakIndex(channelReadClass *reader)
 {
   int imax = -1;
   double maxA = -999;
   double a;
   for (int ipeak = 0; ipeak < reader -> nPeaks; ++ipeak) {
     a = reader -> PeakVoltage[ipeak];
     if (a > maxA) {
       maxA = a;
       imax = ipeak;
     }
   }
   return imax;
   //   return 0; 
 }

// ______________________________________________________________
// peakMode: "", a, b, c, d, e, f, g, h, i

void MakeAllDataPlots::Loop(int verbose, int debug) {

  _debug = debug;
  // +-------------------------------+
  // |         event loop            |
  // +-------------------------------+

  cout << "Event loop!" << endl;
  
  // TODO:
  // check also the number of entries in the trees?

  _nSpills = 0;
  _lastSpillNumber = -999;
  
  // for(int ientry = 0; ientry < _ent[0]; ientry++) {
  for(int ientry = 0; ientry < _Nmin; ientry++) {

    if (ientry % verbose == 0) {
      cout << "processing " << ientry << " / " << _ent[0] << endl;
    }
    _eventInfo -> LoadTree(ientry);
    _eventInfo -> GetEntry(ientry);
    Long64_t  RunNumber = _eventInfo->RunNumber;
    Int_t EventNumber = _eventInfo -> EventNumber;
    Int_t SpillNumber = _eventInfo -> SpillNumber;
    if (_lastSpillNumber != SpillNumber) {
      _nSpills++;
      _lastSpillNumber = SpillNumber;
    }
    
    //    cout << " RunNumber=" << RunNumber << " EventNumber=" << EventNumber << " SpillNumber=" << SpillNumber << endl;
    
    for (int ich = 0; ich < _nChannels; ++ich) {
      if (_debug)	cout << "getting entry for " <<  _treeNames[ich] << endl;
      _reader[ich] -> LoadTree(ientry);
      _reader[ich] -> GetEntry(ientry);
    }
    if (_debug)      cout << "done" << endl;


    this -> ReadChannels();
    
    // peak cuts on demand
    if (!_isHodoscopeRun) {
      if (! this -> PassedPeakCuts())
	continue;
    } // not hodoscope run

    this -> FillChannels();
    

    this -> FillTofHistos();
    if (_debug) cout << "filled tof" << endl;
    if (!_isHodoscopeRun) {
      if (_debug) cout << "filling charged" << endl;
      this -> FillChargedHistos();
    } else {
      if (_debug) cout << "filling hodoscope" << endl;
      this -> FillHodoscopeHistos();
    }
    
  } // entries


} // Run


// ______________________________________________________________
bool MakeAllDataPlots::PassedPeakCuts()
{

  //    vector<int> indices(_nChannels, 0);
      
  _onePeakInAllACTs = true;
  _onePeakInAllToFs = true;
  _onePeakInAll = true;
      
  _moreThanOnePeakInAllACTs = true;
  _moreThanOnePeakInAllToFs = true;
  _moreThanOnePeakInAll = true;
      
  _PbGlassAboveElectronLevel = true;
  _ACT23AboveElectronLevel = true;
      
  double PbGlassElectronThreshA = 5;
  double PbGlassElectronUpperThreshA = 6.5;
      
  double ACTC23ElectronThreshA = 1.5;
  double ACTC23ElectronUpperThreshA = 3.5;
      
  if (_debug)      cout << "point a" << endl;
  _onePeakInPbGlass = (NPeaks["PbGlass"] == 1);
  if (_debug)      cout << "point b" << endl;
      
  // this is WCTE TB 2023 Run1 (4xACTSm trigger tofs lead glas; no hodoscope) specific!
  // to be updated for the highest peak
  // so all [0] need to be changed to appropriate PeakID
  for(int j = 0; j < _nChannels; j++) {
    if (j < 16) {
      _onePeakInAll = _onePeakInAll && (_reader[j] -> nPeaks == 1);
      // by All we mean trigger tofs and ACTs, not leadglas nor hole counters
      _moreThanOnePeakInAll = _moreThanOnePeakInAll && (_reader[j] -> nPeaks > 1);
    }
    if (j < 8) {
      _onePeakInAllACTs = _onePeakInAllACTs && (_reader[j] -> nPeaks == 1);
      _moreThanOnePeakInAllACTs = _moreThanOnePeakInAllACTs && (_reader[j] -> nPeaks > 1);
    }
    if (j >= 8 && j < 16) {
      _onePeakInAllToFs = _onePeakInAllToFs && (_reader[j] -> nPeaks == 1);
      _moreThanOnePeakInAllToFs = _moreThanOnePeakInAllToFs && (_reader[j] -> nPeaks > 1);
    }
    
    // dirty add-on to select electrons
    if (j == 18) {
      _PbGlassAboveElectronLevel = _PbGlassAboveElectronLevel && (_reader[j] -> PeakVoltage[0] > PbGlassElectronThreshA);
      _PbGlassAboveElectronLevel = _PbGlassAboveElectronLevel && (_reader[j] -> PeakVoltage[0] < PbGlassElectronUpperThreshA);
    }
    if (j == 4) {
      _ACT23AboveElectronLevel = _ACT23AboveElectronLevel && ((_reader[j] -> PeakVoltage[0] + _reader[j+1] -> PeakVoltage[0] + _reader[j+2] -> PeakVoltage[0] + _reader[j+3] -> PeakVoltage[0])/2. > ACTC23ElectronThreshA);
	  
      _ACT23AboveElectronLevel = _ACT23AboveElectronLevel && ((_reader[j] -> PeakVoltage[0] + _reader[j+1] -> PeakVoltage[0] + _reader[j+2] -> PeakVoltage[0] + _reader[j+3] -> PeakVoltage[0])/2. < ACTC23ElectronUpperThreshA);
    }
  } // channels
      
      
  if (_peakMode == "a" && ! (_onePeakInAll) )
    return false;
  if (_peakMode == "b" && (! _moreThanOnePeakInAllACTs) )
    return false;
  if (_peakMode == "c" && ! (_onePeakInAllACTs && _moreThanOnePeakInAllToFs) )
    return false;
  if (_peakMode == "d" && ! (_moreThanOnePeakInAllACTs && _onePeakInAllToFs) )
    return false;
  if (_peakMode == "e" && ! (_onePeakInAllACTs) )
    return false;
  if (_peakMode == "f" && ! (_onePeakInAllToFs) )
    return false;
  if (_peakMode == "g" && ( ! (_onePeakInAllToFs) || ! (_ACT23AboveElectronLevel) || ! (_PbGlassAboveElectronLevel)))
    return false;
  if (_peakMode == "h" && ( ! (_onePeakInAllToFs && _onePeakInPbGlass)) )
    return false;

  if (_peakMode == "i" && (  (Charges["TOF10"] < 0.04 && Charges["TOF10"] > 0.03) || (Amplitudes["TOF00"] > 1.43 && Amplitudes["TOF00"] < 1.435) ) )
    return false;


  
  return true;

}

// ______________________________________________________________
void MakeAllDataPlots::ReadChannels()
{

  // +-------------------------------------------------------+
  // |   read all channels information for all waveforms!    |
  // +-------------------------------------------------------+

    for (int ich = 0; ich < _nChannels; ++ich) {
      TString chname = _treeNames[ich];
      if (_debug)      cout << "point c, " << chname.Data() << endl;
      PeakID[chname] = getHighestPeakIndex(_readerMap[chname]);
      int ipeak = PeakID[chname];
      int iwindowpeak = 0; // !!!
      NPeaks[chname] = _readerMap[chname] -> nPeaks;
      if ( ipeak >= 0 && ipeak < _readerMap[chname] -> nPeaks) {
	Amplitudes[chname]  = _readerMap[chname] -> PeakVoltage[ipeak];
	if (_useWindowIntCharge) // && !chname.Contains("TOF"))
	  Charges[chname]     = _readerMap[chname] -> WindowIntCharge[iwindowpeak];
	else
	  Charges[chname]     = _readerMap[chname] -> IntCharge[ipeak];
	SignalTimes[chname] = _readerMap[chname] -> SignalTime[ipeak];
      } else {
	Amplitudes[chname]  = 0.;
	Charges[chname]     = 0.;
	SignalTimes[chname] = 0.;
      }
      if (_debug)      cout << "point d" << endl;

    } // channels
}



// ______________________________________________________________
void MakeAllDataPlots::FillChannels()
{

  for (int ich = 0; ich < _nChannels; ++ich) {
      TString chname = _treeNames[ich];
      if (_debug)      cout << "point c, " << chname.Data() << endl;
      int ipeak = PeakID[chname];
      if ( ipeak >= 0 && ipeak < _readerMap[chname] -> nPeaks) {

	// can be simplified using the above maps
	
	hCharge.at(ich).Fill(Charges[chname]);
	hVoltage.at(ich).Fill(Amplitudes[chname]);
	hTime.at(ich).Fill(SignalTimes[chname]);
	hPedestalSigma.at(ich).Fill(_reader[ich] -> PedestalSigma);
	hnPeaks.at(ich).Fill(_reader[ich] -> nPeaks);

	/*
	// histograms over all channels
	hCharge.at(ich).Fill(_reader[ich] -> IntCharge[ipeak]);
	hVoltage.at(ich).Fill(_reader[ich] -> PeakVoltage[ipeak]);
	hTime.at(ich).Fill(_reader[ich] -> SignalTime[ipeak]);
	//hNbPeaks.at(ich).Fill(_reader[ich] -> nPeaks);
	hPedestalSigma.at(ich).Fill(_reader[ich] -> PedestalSigma);
	hnPeaks.at(ich).Fill(_reader[ich] -> nPeaks);
	*/
	
      } 
    } // channels
}
// ______________________________________________________________

void MakeAllDataPlots::FillTofHistos()
{

  // TOF trigger scintilators

    double t00 = SignalTimes["TOF00"];
    double t01 = SignalTimes["TOF01"];
    double t02 = SignalTimes["TOF02"];
    double t03 = SignalTimes["TOF03"];

    double t10 = SignalTimes["TOF10"];
    double t11 = SignalTimes["TOF11"];
    double t12 = SignalTimes["TOF12"];
    double t13 = SignalTimes["TOF13"];

    if (_debug)      cout << "point e" << endl;

    // JK's time resolution of 2022, diagonal combinations
    double t0a = (t00 + t03) / 2.;
    double t0b = (t01 + t02) / 2.;
    double t1a = (t10 + t13) / 2.;
    double t1b = (t11 + t12) / 2.;

    // time diffs for time resolution histogramme 2022:
    double t0diff = t0a - t0b;
    double t1diff = t1a - t1b;

    // jiri
    // Fill resolution histograms of 2022:
    _histos1d["hTimeReso0"]->Fill(t0diff);
    _histos1d["hTimeReso1"]->Fill(t1diff);
    _histos1d["hTimeReso0_zoom"]->Fill(t0diff);
    _histos1d["hTimeReso1_zoom"]->Fill(t1diff);

    _histos1d["hTimeDiffTOF01"]->Fill(t00 - t01); // carefull, this particular histogram is a difference of hit times,
    _histos1d["hTimeDiffTOF02"]->Fill(t00 - t02); // not a TOF per se but instead the cable length difference
    _histos1d["hTimeDiffTOF03"]->Fill(t00 - t03); // plus the photon travel time though the panel
 
    _histos1d["hTimeDiffTOF11"]->Fill(t11 - t10);
    _histos1d["hTimeDiffTOF12"]->Fill(t11 - t12);
    _histos1d["hTimeDiffTOF13"]->Fill(t11 - t13);
    
    // acraplet
    // compare the hit times for the same event recorded by trigger PMTs on the same side (up/down, left/right)
    // assumption: the light travel time trough the trigger counter to the PMT should be about the same (if the beam is well aligned)
    // then we can check if the TOF is constant for a run
    // idea: using tof dofference we can triangulate the position of a given pulse on the trigger! could be a fun thing to check
    // this is after the calibration
    _histos1d["hTimeTOF0"]->Fill(t11 - t00); // positive tof
    _histos1d["hTimeTOF1"]->Fill(t13 - t02);
    _histos1d["hTimeTOF2"]->Fill(t10 - t01);
    _histos1d["hTimeTOF3"]->Fill(t12 - t03);

    // Time of flight!
    _t0 = ( t00 + t01 + t02 + t03 ) / 4.;
    _t1 = ( t10 + t11 + t12 + t13 ) / 4.;
    _tof = _t1 - _t0;

    if (_debug)      cout << "done tof" << endl;
}

// ______________________________________________________________

void MakeAllDataPlots::FillChargedHistos()
{
    if (_debug)      cout << "charged a" << endl;

    // ACTs
    double act0c = Charges["ACT0L"] + Charges["ACT0R"];
    double act1c = Charges["ACT1L"] + Charges["ACT1R"];
    double act2c = Charges["ACT2L"] + Charges["ACT2R"];
    double act3c = Charges["ACT3L"] + Charges["ACT3R"];

    double act0a = Amplitudes["ACT0L"] + Amplitudes["ACT0R"];
    double act1a = Amplitudes["ACT1L"] + Amplitudes["ACT1R"];
    double act2a = Amplitudes["ACT2L"] + Amplitudes["ACT2R"];
    double act3a = Amplitudes["ACT3L"] + Amplitudes["ACT3R"];

    double act23aAver = (act2a + act3a) / 2.;
    double act23cAver = (act2c + act3c) / 2.;

    // hole counters and lead glass
    
    double hc0c = Charges["Hole0"];
    double hc0a = Amplitudes["Hole0"];

    double hc1c = Charges["Hole1"];
    double hc1a = Amplitudes["Hole1"];

    double pbc = Charges["PbGlass"];
    double pba = Amplitudes["PbGlass"];

    double trigScintA = Charges["TOF00"] + Charges["TOF01"] + Charges["TOF02"] + Charges["TOF03"] + Charges["TOF10"] + Charges["TOF11"] + Charges["TOF12"] + Charges["TOF13"];
    double trigScintC = Amplitudes["TOF00"] + Amplitudes["TOF01"] + Amplitudes["TOF02"] + Amplitudes["TOF03"] + Amplitudes["TOF10"] + Amplitudes["TOF11"] + Amplitudes["TOF12"] + Amplitudes["TOF13"];

    double trigScint0A = Charges["TOF00"] + Charges["TOF01"] + Charges["TOF02"] + Charges["TOF03"];
    double trigScint0C = Amplitudes["TOF00"] + Amplitudes["TOF01"] + Amplitudes["TOF02"] + Amplitudes["TOF03"];
    double trigScint001A = Charges["TOF00"] + Charges["TOF01"];
    double trigScint001C = Amplitudes["TOF00"] + Amplitudes["TOF01"];
    double trigScint023A = Charges["TOF02"] + Charges["TOF03"];
    double trigScint023C = Amplitudes["TOF02"] + Amplitudes["TOF03"];

    double trigScint1A = Charges["TOF10"] + Charges["TOF11"] + Charges["TOF12"] + Charges["TOF13"];
    double trigScint1C = Amplitudes["TOF10"] + Amplitudes["TOF11"] + Amplitudes["TOF12"] + Amplitudes["TOF13"];
    double trigScint101A = Charges["TOF10"] + Charges["TOF11"];
    double trigScint101C = Amplitudes["TOF10"] + Amplitudes["TOF11"];
    double trigScint123A = Charges["TOF12"] + Charges["TOF13"];
    double trigScint123C = Amplitudes["TOF12"] + Amplitudes["TOF13"];


    // cout << " trigScintA=" << trigScintA << " trigScintC=" << trigScintC << endl;

    // HACK!
    //    pba = _readerMap["PbGlass"] -> PeakVoltage[0];

    // 17.11.2023:
    // A's are zero if no peaks;)
    if (Amplitudes["ACT0L"] > 0.05 &&  Amplitudes["ACT0R"] > 0.05) {
      _histos2d["hRef_act0LA_act0RA_noneZero"]->Fill(Amplitudes["ACT0L"], Amplitudes["ACT0R"] );
    }
    if (Amplitudes["ACT1L"] > 0.01 &&  Amplitudes["ACT1R"] > 0.01) {
      _histos2d["hRef_act1LA_act1RA_noneZero"]->Fill(Amplitudes["ACT1L"], Amplitudes["ACT1R"] );
    }
    if (Amplitudes["ACT2L"] > 0.01 &&  Amplitudes["ACT2R"] > 0.01) {
      _histos2d["hRef_act2LA_act2RA_noneZero"]->Fill(Amplitudes["ACT2L"], Amplitudes["ACT2R"] );
    }
    if (Amplitudes["ACT3L"] > 0.01 &&  Amplitudes["ACT3R"] > 0.01) {
      _histos2d["hRef_act3LA_act3RA_noneZero"]->Fill(Amplitudes["ACT3L"], Amplitudes["ACT3R"] );
    }

    if (Amplitudes["ACT0L"] > 0.05 &&  Amplitudes["ACT1L"] > 0.01) {
      _histos2d["hRef_act0LA_act1LA_noneZero"]->Fill(Amplitudes["ACT0L"], Amplitudes["ACT1L"] );
    }
    if (Amplitudes["ACT0R"] > 0.05 &&  Amplitudes["ACT1R"] > 0.01) {
      _histos2d["hRef_act0RA_act1RA_noneZero"]->Fill(Amplitudes["ACT0R"], Amplitudes["ACT1R"] );
    }

    
    if (_debug)      cout << "charged b" << endl;
    
    // amplitudes vs tof
    _histos2d["hRef_TOFACT0A"]->Fill(_tof, act0a);
    _histos2d["hRef_TOFACT1A"]->Fill(_tof, act1a);
    _histos2d["hRef_TOFACT2A"]->Fill(_tof, act2a);
    _histos2d["hRef_TOFACT3A"]->Fill(_tof, act3a);

    // Trigger scintillators (unfortunatelly labelled as TOF all through out the code;-)
    // amplitude and charge vs tof:
    // TODO
    
    
    // lead glass vs acts and tof
    _histos2d["hRef_pbA_act23A"]->Fill(pba, act23aAver);
    _histos2d["hRef_pbC_act23A"]->Fill(pbc, act23aAver);
    _histos2d["hRef_TOFACT23A"]->Fill(_tof, act23aAver);

    _histos2d["hRef_pbA_act0A"]->Fill(pba, act0a);
    _histos2d["hRef_pbA_act1A"]->Fill(pba, act1a);
    _histos2d["hRef_pbA_act1C"]->Fill(pba, act1c);
    
    _histos2d["hRef_PbATOF"]->Fill(pba, _tof);
    _histos2d["hRef_TOFPbA"]->Fill(_tof, pba);
    _histos2d["hRef_TOFPbC"]->Fill(_tof, pbc);

    // charges vs tof
    _histos2d["hRef_TOFACT0C"]->Fill(_tof, act0c);
    _histos2d["hRef_TOFACT1C"]->Fill(_tof, act1c);
    _histos2d["hRef_TOFACT2C"]->Fill(_tof, act2c);
    _histos2d["hRef_TOFACT3C"]->Fill(_tof, act3c);

    // lead glass vs acts and tof
    _histos2d["hRef_pbC_act23C"]->Fill(pbc, act23cAver);
    _histos2d["hRef_pbA_act23C"]->Fill(pba, act23cAver);
    _histos2d["hRef_TOFACT23C"]->Fill(_tof, act23cAver);
    
    _histos2d["hRef_pbC_act0C"]->Fill(pbc, act0c);
    _histos2d["hRef_TOFACT0C"]->Fill(_tof, act0c);
    _histos2d["hRef_pbC_act1C"]->Fill(pbc, act1c);
    _histos2d["hRef_TOFACT1C"]->Fill(_tof, act1c);
    
    _histos2d["hRef_PbCTOF"]->Fill(pbc, _tof);

    // act 2d plots
    _histos2d["hACT1CACT3C"]->Fill(act1c, act3c);
    _histos2d["hACT3CACT2C"]->Fill(act3c, act2c);
    _histos2d["hACT2CACT1C"]->Fill(act2c, act1c);

    // acraplet - weird electrons which do not see anything in the ACT
    if (act23aAver != 1.5 && _tof >= 13.5 && _tof <= 16.5) {
      _histos2d["hHC0AHC1A"]->Fill(hc0a, hc1a);
      _histos2d["hHC0CHC1C"]->Fill(hc0c, hc1c);
    }


    // JK 5.3.2024
    // both trigger scintillators:
    _histos2d["hRef_pbC_TrigScintC"]->Fill(pbc, trigScintC);
    _histos2d["hRef_pbA_TrigScintC"]->Fill(pba, trigScintC);
    _histos2d["hRef_pbC_TrigScintA"]->Fill(pbc, trigScintA);
    _histos2d["hRef_pbA_TrigScintA"]->Fill(pba, trigScintA);
    
    _histos2d["hRef_TOF_TrigScintC"]->Fill(_tof, trigScintC);
    _histos2d["hRef_TOF_TrigScintA"]->Fill(_tof, trigScintA);

    // TOF0X trigger scintillators:
    _histos2d["hRef_pbC_TrigScint0C"]->Fill(pbc, trigScint0C);
    _histos2d["hRef_pbA_TrigScint0C"]->Fill(pba, trigScint0C);
    _histos2d["hRef_pbC_TrigScint0A"]->Fill(pbc, trigScint0A);
    _histos2d["hRef_pbA_TrigScint0A"]->Fill(pba, trigScint0A);
    
    _histos2d["hRef_TOF_TrigScint0C"]->Fill(_tof, trigScint0C);
    _histos2d["hRef_TOF_TrigScint0A"]->Fill(_tof, trigScint0A);

    // TOF0 0+1 trigger scintillators:
    _histos2d["hRef_pbC_TrigScint001C"]->Fill(pbc, trigScint001C);
    _histos2d["hRef_pbA_TrigScint001C"]->Fill(pba, trigScint001C);
    _histos2d["hRef_pbC_TrigScint001A"]->Fill(pbc, trigScint001A);
    _histos2d["hRef_pbA_TrigScint001A"]->Fill(pba, trigScint001A);
    
    _histos2d["hRef_TOF_TrigScint001C"]->Fill(_tof, trigScint001C);
    _histos2d["hRef_TOF_TrigScint001A"]->Fill(_tof, trigScint001A);

    // TOF0 2+3 trigger scintillators:
    _histos2d["hRef_pbC_TrigScint023C"]->Fill(pbc, trigScint023C);
    _histos2d["hRef_pbA_TrigScint023C"]->Fill(pba, trigScint023C);
    _histos2d["hRef_pbC_TrigScint023A"]->Fill(pbc, trigScint023A);
    _histos2d["hRef_pbA_TrigScint023A"]->Fill(pba, trigScint023A);
    
    _histos2d["hRef_TOF_TrigScint023C"]->Fill(_tof, trigScint023C);
    _histos2d["hRef_TOF_TrigScint023A"]->Fill(_tof, trigScint023A);

    
    // TOF1X trigger scintillators:
    _histos2d["hRef_pbC_TrigScint1C"]->Fill(pbc, trigScint1C);
    _histos2d["hRef_pbA_TrigScint1C"]->Fill(pba, trigScint1C);
    _histos2d["hRef_pbC_TrigScint1A"]->Fill(pbc, trigScint1A);
    _histos2d["hRef_pbA_TrigScint1A"]->Fill(pba, trigScint1A);
    
    _histos2d["hRef_TOF_TrigScint1C"]->Fill(_tof, trigScint1C);
    _histos2d["hRef_TOF_TrigScint1A"]->Fill(_tof, trigScint1A);

    // TOF1 0+1 trigger scintillators:
    _histos2d["hRef_pbC_TrigScint101C"]->Fill(pbc, trigScint101C);
    _histos2d["hRef_pbA_TrigScint101C"]->Fill(pba, trigScint101C);
    _histos2d["hRef_pbC_TrigScint101A"]->Fill(pbc, trigScint101A);
    _histos2d["hRef_pbA_TrigScint101A"]->Fill(pba, trigScint101A);
    
    _histos2d["hRef_TOF_TrigScint101C"]->Fill(_tof, trigScint101C);
    _histos2d["hRef_TOF_TrigScint101A"]->Fill(_tof, trigScint101A);

    // TOF1 2+3 trigger scintillators:
    _histos2d["hRef_pbC_TrigScint123C"]->Fill(pbc, trigScint123C);
    _histos2d["hRef_pbA_TrigScint123C"]->Fill(pba, trigScint123C);
    _histos2d["hRef_pbC_TrigScint123A"]->Fill(pbc, trigScint123A);
    _histos2d["hRef_pbA_TrigScint123A"]->Fill(pba, trigScint123A);
    
    _histos2d["hRef_TOF_TrigScint123C"]->Fill(_tof, trigScint123C);
    _histos2d["hRef_TOF_TrigScint123A"]->Fill(_tof, trigScint123A);


    // times
    _histos1d["hT0"]->Fill(_t0);
    _histos1d["hT1"]->Fill(_t1);

    if (_debug)      cout << "charged c" << endl;

    // 2D <nPeaks> studies
    int nPeaksToFAver = (_readerMap["TOF00"]->nPeaks + _readerMap["TOF01"]->nPeaks + _readerMap["TOF02"]->nPeaks + _readerMap["TOF03"]->nPeaks + _readerMap["TOF10"]->nPeaks + _readerMap["TOF11"]->nPeaks + _readerMap["TOF12"]->nPeaks + _readerMap["TOF13"]->nPeaks) / 8;
    int nPeaksACT23Aver = (_readerMap["ACT2L"]->nPeaks + _readerMap["ACT2R"]->nPeaks + _readerMap["ACT3L"]->nPeaks + _readerMap["ACT3R"]->nPeaks) / 4;

    int nPeaksToF0Aver = (_readerMap["TOF00"]->nPeaks + _readerMap["TOF01"]->nPeaks + _readerMap["TOF02"]->nPeaks + _readerMap["TOF03"]->nPeaks ) / 4;
    int nPeaksToF1Aver = (_readerMap["TOF10"]->nPeaks + _readerMap["TOF11"]->nPeaks + _readerMap["TOF12"]->nPeaks + _readerMap["TOF13"]->nPeaks) / 4;
    int nPeaksACT2Aver = (_readerMap["ACT2L"]->nPeaks + _readerMap["ACT2R"]->nPeaks ) / 2;
    int nPeaksACT3Aver = (_readerMap["ACT3L"]->nPeaks + _readerMap["ACT3R"]->nPeaks) / 2;

    int nPeaksLeadGlass = _readerMap["PbGlass"]->nPeaks;
   
    _histos2d["hnPeaksACT23vsnPeaksToF"]->Fill(nPeaksToFAver, nPeaksACT23Aver);

    _histos2d["hnPeaksToF1vsnPeaksToF0"]->Fill(nPeaksToF1Aver, nPeaksToF0Aver);
    _histos2d["hnPeaksACT3vsnPeaksACT2"]->Fill(nPeaksACT2Aver, nPeaksACT3Aver);

    _histos2d["hnPeaksACT23vsToF"]->Fill(_tof,  nPeaksACT23Aver);
    _histos2d["hnPeaksACT23vsToFlow"]->Fill(_tof,  nPeaksACT23Aver);
    _histos2d["hnPeaksToFvsToF"]->Fill(_tof, nPeaksToFAver);
    _histos2d["hnPeaksToFvsToFlow"]->Fill(_tof, nPeaksToFAver);
    _histos2d["hnPeaksLeadGlassvsLeadGlassA"]->Fill(pba, nPeaksLeadGlass);

    _histos2d["hnPeaksACT23vsLeadGlassA"]->Fill(pba, nPeaksACT23Aver);
    _histos2d["hnPeaksToFvsLeadGlassA"]->Fill(pba, nPeaksToFAver);

    // selections:

    bool isEl = false;

    bool passed_act1a_cuts = false;
    bool passed_act2a_cuts = false;

    //    switch(_momentum)	  {

    if (_noAct1Cuts) {
      // originally designed on 900 MeV/c
      // Alie:
      /*
      if (_tof >= _cutsMap[900]["tof_t1_cut"] && _onePeakInAllToFs) {
	_isdACT23pb = true;
      } else if (_tof >= _cutsMap[900]["tof_t0_cut"] && _onePeakInAllToFs){
	_ispACT23pb = true;
      } else if ( (act2a+act3a)/2 < _cutsMap[900]["act23_pi_minA"] && _tof <= _cutsMap[900]["tof_t0_cut"] && _onePeakInAllToFs) {
	_isMuACT23pb = true;
      } else if ((act2a+act3a)/2 > _cutsMap[900]["act23_pi_minA"] && _tof <= _cutsMap[900]["tof_t0_cut"] && _onePeakInAllToFs) {
	_isElACT23pb = true;
      }
      */
      

   
    } else {

      // TODO
      // to move to a map, too
      /*
      // also ACT1 cuts
      // originally designed on 420 MeV/c
      if ( pbc > 0.9 || pbc < 0.1 || act1c > 0.25) {  // custom electron removal cut (act2a + act3a)/2.
	  isEl = true;
      }
      if ((act2a+act3a)/2 > y0_cut - y0_cut/x0_cut * pba && pba > pb_min && _onePeakInAllToFs) {
	std::cout << _onePeakInAllToFs << std::endl;
	isElACT23pb = true;
      }
      else if ( (act2a+act3a)/2 > act23_pi_maxA && pba > pb_min && _onePeakInAllToFs) {
	isMuACT23pb = true;
      }
      else if ( (act2a+act3a)/2 > act23_pi_minA && pba > pb_min && _onePeakInAllToFs) {
	isPiACT23pb = true;
      }
      else if ( pba < pb_min && !_onePeakInAllToFs) {
	ispACT23pb = true;
      }
      */

      // JK looing at histos/windowpe_analyzed/windowPE_-16ns_45ns_run000352_plots_f.root
      //if (act1a < 0.15 && pba < 0.25) {

      // cuts III
      //if (act1a < 0.15 && pba < 0.27 && act23aAver < 0.55) {
      // cuts II
      //if (act1a < 0.15 && pba < 0.27) {
      // cuts I
      if (act1a < 0.15) {
	// this cut is too destructive:  act0a < 0.15
	//also too destructive: pba < 0.15
	passed_act1a_cuts = true;
	_histos2d["hRef_pbA_act1A_act1cuts"]->Fill(pba, act1a);
	_histos1d["hTOFOther_act1cuts"]->Fill(_tof);
	_histos1d["hTOFOtherLow_act1cuts"]->Fill(_tof);
      
	
      }
      
      

    }

	
	 //default: {
	//        if (ientry < 10)
        //  cout << "WARNING: Using default settings for the " << _momentum << " MeV/c beam" << endl;
	// add ACT1 cuts?
        // JK if ( (act2a + act3a)/2. > 3.) { // custom electron removal cut
        //  isEl = true;
        //}

	 //  } // default
	 //    } // case

    if (isEl) {
      // electrons
      _histos1d["hTOFEl"]->Fill(_tof);
      _histos1d["hTOFElLow"]->Fill(_tof);
    }
    else {
      // non-electrons
      _histos1d["hTOFOther"]->Fill(_tof);
      _histos1d["hTOFOtherLow"]->Fill(_tof);
    } // non-electrons

    _histos1d["hTOFAll"]->Fill(_tof);
    _histos1d["hTOFAllWide"]->Fill(_tof);
    _histos1d["hTOFAllLow"]->Fill(_tof);

    if (_debug)      cout << "done charged" << endl;
    
}

// ______________________________________________________________
void MakeAllDataPlots::FillHodoscopeHistos() {


  // Jiri
  double pbanew = Amplitudes["LeadGlass"];
  int maxhdchid = -1;
  double maxA = -1;
  bool hits[15];
  for(int j = 0; j < 15; ++j)
    hits[j] = false;
  for(int j = 0; j < _nChannels; j++) {
    TString chname = _treeNames[j];
    if (! chname.Contains("HD"))
      continue;
    int hdchid = -1;
    double Aj = Amplitudes[chname];
    if (j >= 17 && j < 24)
      hdchid = j - 9;
    else
      hdchid = j - 24;
    if (Aj > 0.12) {
      hits[hdchid] = true;
      _histos2d["LeadGlassPhotonAVsPositronHodoOcc"] -> Fill(hdchid, pbanew);
      if (Aj > maxA) {
	maxA = Aj;
	maxhdchid = hdchid;
      }
    }  // A cut
  } // channels
  
  if (maxhdchid >= 0)
    _histos2d["LeadGlassPhotonAVsPositronMaxHodoOcc"] -> Fill(maxhdchid, pbanew);
  for(int j = 0; j < 15; ++j) {
    for(int k = 0; k <= j; ++k) {
      if (hits[j] && hits[k])
	_histos2d["HodoOccScatter"] -> Fill(j,k); 
    }
  }
  
  
  // Alie:
  // JK: the cut to check!
  double threshHodoscopeHit = 0.12; //Threshold estimated by hand using the online hist as reference 
  //each channel has a different 1pe peak (in particular ch11 has a low one) 
  // the online analysis threshold of 400mV corresponds to 1.9 in these units but we are cutting a lot of hits especially in channel 11, using 1.5 is better there  
  //careful ! position of the detectors on the digitiser have moved!!!
  for (int i = 0; i < 15; i++){
    TString hdname = Form("HD%i", i);
    if (Amplitudes[hdname] >= threshHodoscopeHit){
      _histos1d["hnHitsHodoscope"] -> Fill(i);
      //_histos1d["hnHitsHodoscope"] -> Fill(_channelToHodoscope[i]);
    }
  }
  


  
}

// ______________________________________________________________


void MakeAllDataPlots::Terminate()
{

  cout << "End of event loop!" << endl;
  // normalization of the fraction 2D occupancy map
  if (_isHodoscopeRun) {
    for(int j = 0; j < 15; ++j) {
      double diag = _histos2d["HodoOccScatter"] -> GetBinContent(j,j);
      //cout << "diag=" << diag << endl;
      for(int k = 0; k <= j; ++k) {
	if (diag > 0.)
	  _histos2d["HodoOccScatterFrac"] -> SetBinContent(j, k, _histos2d["HodoOccScatter"] -> GetBinContent(j,k) / diag);
      } // k
    } // j
    _histos2d["HodoOccScatterFrac"] -> Scale(1.);
    //    _histos2d["HodoOccScatter"] -> Scale(1.);
  }

  _outFile -> cd();
  TH1D *spills_h = new TH1D("nspills", "nspills", 1, 0, 1);
  spills_h -> Fill(0.5, _nSpills);
  
  _outFile -> Write();

}


// ______________________________________________________________
