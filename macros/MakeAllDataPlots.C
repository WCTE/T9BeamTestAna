#include "MakeAllDataPlots.h"
#include "Tools.C"

// jiri kvita 2023--2024

using namespace std;

// ______________________________________________________________


MakeAllDataPlots::MakeAllDataPlots(string fileName, int momentum, bool isHodoscopeRun, TString peakMode, bool useWindowIntCharge ) {

  _fileName = fileName;
  _momentum = momentum;
  _peakMode = peakMode;
  _isHodoscopeRun = isHodoscopeRun;
  _useWindowIntCharge = useWindowIntCharge;
  _debug = 0;

  // Luan, 7.3.2024
  // nominal diameter of the PMTs is 2.85 cm
  // thickness of the scintillator tile is 1 cm

  double dpmt = 2*16.73; //31.85; //28.5; // mm
  double delta = 11.55; // mm
  double sf = 1.;
  double rx = 86./2; // mm
  double ry = 86./2; // mm
  _xL = -rx*sf; // guess, mm
  _xR = +rx*sf; // guess, mm
  _yD = (-ry+(dpmt + delta)/2.)*sf; // guess, mm
  _yU = (+ry-(dpmt + delta)/2.)*sf; // guess, mm
  /*
    _xL = -rx*sf; // guess, mm
    _xR = +rx*sf; // guess, mm
    _yD = (-ry+dpmt/2.)*sf; // guess, mm
    _yU = (+ry-dpmt/2.)*sf; // guess, mm
  */


  // ACT
  _ACTwidth = 100; // mm
  
  cout << "MakeAllDataPlots::MakeAllDataPlot: Configured as:" << endl;
  cout << " fileName           : " << _fileName.c_str() << endl 
       << " momentum           : " << _momentum << endl 
       << " peakMode           : " << _peakMode << endl 
       << " isHodoscopeRun     : " << _isHodoscopeRun << endl
       << " useWindowIntCharge : " << useWindowIntCharge << endl
       << endl;

  cout << "Estimated Trigger Scintillators PMTs be apart by " << _yU - _yD << " mm" << endl;
  cout << "Estimated Trigger Scintillators dimensions: " << 2*rx << endl;
  
  _tofutil = new tofUtil();
  
}
// ______________________________________________________________


MakeAllDataPlots::~MakeAllDataPlots() {

}
// ______________________________________________________________

void MakeAllDataPlots::Init(bool noAct1Cuts)
{

  _noAct1Cuts = noAct1Cuts;

  _tofmin = 10.;
  _tofmax = 40.;
  _tofmaxhigh = 85.;
  _ntofbins = 200;

  _tofminlow = 10.;
  _tofmaxlow = 20.;
  _ntofbinslow = 150;

  _ntofbins2d = 400;

  // scale factors for p.e. charges limits
  double PEsfTOF = 150;
  double PEsfPbG = 300.;
  double PEsfACT = 10.;

  if (!_useWindowIntCharge) {
    PEsfTOF = 1.;
    PEsfPbG = 1.;
    PEsfACT = 1.;
  }
  
  _actChargeMin = 0.0;
  _actChargeMax = 2.*PEsfACT; // /2. 1.1* 2.*
  _actAmplitudeMax =  2.; // 2.

  _PbGAmplitudeMin =  0.; // 2.
  _PbGAmplitudeMax =  2.; // 2.

  _PbGChargeMin =  0.; // 2.
  _PbGChargeMax =  2*PEsfPbG; // 2. // 1.4

  _trigScintChargeMin = 0.;
  _trigScintChargeMax = 2.*PEsfTOF;
  _trigScintAmplitudeMin = 0.;
  _trigScintAmplitudeMax = 2.;//10.*PEsfTOF;
  
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
  cout << "  noAct1Cuts: " << _noAct1Cuts << endl;
    
}


// ______________________________________________________________


void MakeAllDataPlots::InitReaders()
{
  cout << "In InitReaders" << endl;
  _eventInfo = new EventInfo(_infile, "EventInfo");
  cout << "Initializing the tree readers" << endl;
  for (int ich = 0; ich < _nChannels; ++ich) {
    cout << "  initializing " <<  _treeNames[ich] << endl;
    _reader[ich] = new channelReadClass(_infile, _treeNames[ich], _useWindowIntCharge);
    _trees[ich] = _reader[ich] -> fChain;
    _ent[ich] = _trees[ich] -> GetEntries();
    _readerMap[_treeNames[ich]] = _reader[ich];
    cout << "    ...reader for " << _treeNames[ich].Data() << ": " << _ent[ich] << " entries" << endl;
  }
  cout << "Done Init" << endl;

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
  
  //_nChannels = 32;
  cout << "InitGeneralHistos" << endl;
  
  _outFile -> mkdir("General");
  _outFile -> cd("General");

  for (int i = 0; i < _nChannels; ++i) {
    TString chname = _treeNames[i];
    // cout << "Preparing histos for channel " << chname.Data() << endl;
    TString name1 = Form("hRef_Charge%i", i);
    TString name2 = Form("hRef_Voltage%i", i);
    TString name3 = Form("hRef_Hits%i", i);
    TString name4 = Form("hRef_PedestalSigma%i", i);
    TString name5 = Form("hRef_Time%i", i);
    TString name6a = Form("hRef_nPeaksA%i", i);
    TString name6c = Form("hRef_nPeaksC%i", i);
    TString name7 = Form("hRef_Pedestal%i", i);
    TString name8 = Form("hRef_PedestalNbPeaks%i", i);

    TString title1 = Form("Channel %i", i) + TString("; Charge [nC]; Triggers"); 
    TString title2 = Form("Channel %i", i) + TString("; Total Amplitude [V]; Triggers");
    TString title3 = Form("Channel %i", i) + TString("; Hits per trigger; Triggers");
    TString title4 = Form("Channel %i", i) + TString("; #sigma_{ped} [V]; Triggers");
    TString title5 = Form("Channel %i", i) + TString("; Time [ns]; Triggers");
    TString title6a = Form("Channel %i", i) + TString("; Number of ampl. peaks; Triggers");
    TString title6c = Form("Channel %i", i) + TString("; Number of charg. peaks; Triggers");
    TString title7 = Form("Channel %i", i) + TString("; Pedestal Amplitude; Triggers/1mV");
    TString title8 = Form("Channel %i", i) + TString("; Pedestal Amplitude; Number of peaks");

    double cmin = _actChargeMin;
    double cmax = _actChargeMax;
    double amin = _actAmplitudeMin;
    double amax = _actAmplitudeMax;
    if (chname.Contains("TOF")) {
      //      cout << "Setting histo limits to TOF" << endl;
      cmin = _trigScintChargeMin;
      cmax = _trigScintChargeMax;
      amin = _trigScintAmplitudeMin;
      amax = _trigScintAmplitudeMax;
    } else if (chname.Contains("PbG")) {
      //      cout << "Setting histo limits to PbG" << endl;
      cmin = _PbGChargeMin;
      cmax = _PbGChargeMax;
      amin = _PbGAmplitudeMin;
      amax = _PbGAmplitudeMax;
    }
    
    TH1D temp1(name1, title1, 400, cmin, cmax);
    TH1D temp2(name2, title2, 400, amin, amax);
    TH1D temp3(name3, title3, 5, -0.5, 4.5);
    TH1D temp4(name4, title4, 200, 0., 0.01);
    TH1D temp5(name5, title5, 270, 0., 540.);
    TH1D temp6a(name6a, title6a, 20, 0., 20.);
    TH1D temp6c(name6c, title6c, 20, 0., 20.);
    TH1D temp7(name7, title7, 1000., 1.65, 1.65+1000*0.0012207);
    TH2D temp8(name8, title8, 200, 0., 3*0.8,  5, 0., 5.);

    _hCharge.push_back(temp1);
    _hVoltage.push_back(temp2);
    _hPedestalSigma.push_back(temp4);
    _hTime.push_back(temp5);
    _hnPeaksA.push_back(temp6a);
    _hnPeaksC.push_back(temp6c);
  }
 
 _outFile -> cd("../");
 // cout << "Done" << endl;
 
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
  //cout << "Done" << endl;

}


// ______________________________________________________________

void MakeAllDataPlots::InitTofHistos()
{
  cout << "In InitTofHistos" << endl;

  _outFile -> mkdir("TOF");
  _outFile -> cd("TOF");
  // TOF 1D
  _histos1d["hTOFAll"] = new TH1D("hTOFAll", ";t_{TOF}^{All} [ns]", 120, _tofmin, _tofmax);
  _histos1d["hTOFAllWide"] = new TH1D("hTOFAllWide", ";t_{TOF}^{All} [ns]", 2*_ntofbins, _tofmin, 2*_tofmax);
  _histos1d["hTOFEl"] = new TH1D("hTOFEl", ";t_{TOF}^{e} [ns]", _ntofbins, _tofmin, _tofmax);
  _histos1d["hTOFOther"] = new TH1D("hTOFOther", ";t_{TOF}^{non-e} [ns]", _ntofbins, _tofmin, _tofmax);
  _histos1d["hTOFOther_act1cuts"] = new TH1D("hTOFOther_act1cuts", ";t_{TOF}^{non-e} [ns]", _ntofbins, _tofmin, _tofmax);
  
  
  _histos1d["hTOFAllLow"] = new TH1D("hTOFAllLow", ";t_{TOF}^{All} [ns]", _ntofbinslow, _tofminlow, _tofmaxlow);
  _histos1d["hTOFElLow"] = new TH1D("hTOFElLow", ";t_{TOF}^{e} [ns]", _ntofbinslow, _tofminlow, _tofmaxlow);
  _histos1d["hTOFOtherLow"] = new TH1D("hTOFOtherLow", ";t_{TOF}^{non-e} [ns]", _ntofbinslow, _tofminlow, _tofmaxlow);
  _histos1d["hTOFOtherLow_act1cuts"] = new TH1D("hTOFOtherLow_act1cuts", ";t_{TOF}^{non-e} [ns]", _ntofbinslow, _tofminlow, _tofmaxlow);

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

// TDirectory name for ourput file, and selection tak for histo name and selection title

void MakeAllDataPlots::InitTrigScintHistos(TString dirname, TString selTag, TString selTit)
{

  _outFile -> mkdir(dirname);
  _outFile -> cd(dirname);

			    
  int nbs = 400;
  double tofmax = _tofmax;
  if (selTag.Contains("D-") || selTag.Contains("T-") ) {
    tofmax = _tofmaxhigh;
  }
  
 // both trig scintil.:
  _histos2d["hRef_pbC_TrigScintC" + selTag] = new TH2D("hRef_pbC_TrigScintC" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. Charge",
						       nbs, _PbGChargeMin, _PbGChargeMax, nbs, 8*_trigScintChargeMin, 8*_trigScintChargeMax);
  _histos2d["hRef_pbA_TrigScintC" + selTag] = new TH2D("hRef_pbA_TrigScintC" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. Charge",
						       nbs, 0., _PbGAmplitudeMax, nbs,8*_trigScintChargeMin, 8*_trigScintChargeMax);
  _histos2d["hRef_pbC_TrigScintA" + selTag] = new TH2D("hRef_pbC_TrigScintA" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. Amplitude",
						       nbs, _PbGChargeMin, _PbGChargeMax, nbs, 8*_trigScintAmplitudeMin, 8*_trigScintAmplitudeMax);
  _histos2d["hRef_pbA_TrigScintA" + selTag] = new TH2D("hRef_pbA_TrigScintA" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. Amplitude",
						       nbs, 0., _PbGAmplitudeMax, nbs, 8*_trigScintChargeMin, 8*_trigScintAmplitudeMax);

  _histos2d["hRef_TOF_TrigScintC" + selTag] = new TH2D("hRef_TOF_TrigScintC" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. Charge",
						       _ntofbins2d, _tofmin, tofmax, nbs, 8*_trigScintChargeMin, 8*_trigScintChargeMax);
  _histos2d["hRef_TOF_TrigScintA" + selTag] = new TH2D("hRef_TOF_TrigScintA" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. Amplitude",
						       _ntofbins2d, _tofmin, tofmax, nbs, 8*_trigScintAmplitudeMin, 8*_trigScintAmplitudeMax);

  // TOF0X
  _histos2d["hRef_pbC_TrigScint0C" + selTag] = new TH2D("hRef_pbC_TrigScint0C" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 0 Charge",
							nbs, _PbGChargeMin, _PbGChargeMax, nbs, 4*_trigScintChargeMin, 4*_trigScintChargeMax);
  _histos2d["hRef_pbA_TrigScint0C" + selTag] = new TH2D("hRef_pbA_TrigScint0C" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 0 Charge",
							nbs, 0., _PbGAmplitudeMax, nbs, 4*_trigScintChargeMin, 4*_trigScintChargeMax);
  _histos2d["hRef_pbC_TrigScint0A" + selTag] = new TH2D("hRef_pbC_TrigScint0A" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 0 Amplitude",
							nbs, _PbGChargeMin, _PbGChargeMax, nbs, 4*_trigScintAmplitudeMin, 4*_trigScintAmplitudeMax);
  _histos2d["hRef_pbA_TrigScint0A" + selTag] = new TH2D("hRef_pbA_TrigScint0A" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 0 Amplitude",
							nbs, 0., _PbGAmplitudeMax, nbs, 4*_trigScintChargeMin, 4*_trigScintAmplitudeMax);

  _histos2d["hRef_TOF_TrigScint0C" + selTag] = new TH2D("hRef_TOF_TrigScint0C" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 0 Charge",
							_ntofbins2d, _tofmin, tofmax, nbs, 4*_trigScintChargeMin, 4*_trigScintChargeMax);
  _histos2d["hRef_TOF_TrigScint0A" + selTag] = new TH2D("hRef_TOF_TrigScint0A" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 0 Amplitude",
							_ntofbins2d, _tofmin, tofmax, nbs, 4*_trigScintAmplitudeMin, 4*_trigScintAmplitudeMax);

  // TOF0 L
  _histos2d["hRef_pbC_TrigScint0LC" + selTag] = new TH2D("hRef_pbC_TrigScint0LC" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 0L Charge",
							 nbs, _PbGChargeMin, _PbGChargeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_pbA_TrigScint0LC" + selTag] = new TH2D("hRef_pbA_TrigScint0LC" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 0L Charge",
							 nbs, 0., _PbGAmplitudeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_pbC_TrigScint0LA" + selTag] = new TH2D("hRef_pbC_TrigScint0LA" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 0L Amplitude",
							 nbs, _PbGChargeMin, _PbGChargeMax, nbs, 2*_trigScintAmplitudeMin, 2*_trigScintAmplitudeMax);
  _histos2d["hRef_pbA_TrigScint0LA" + selTag] = new TH2D("hRef_pbA_TrigScint0LA" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 0L Amplitude",
							 nbs, 0., _PbGAmplitudeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintAmplitudeMax);

  _histos2d["hRef_TOF_TrigScint0LC" + selTag] = new TH2D("hRef_TOF_TrigScint0LC" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 0L Charge",
							 _ntofbins2d, _tofmin, tofmax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_TOF_TrigScint0LA" + selTag] = new TH2D("hRef_TOF_TrigScint0LA" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 0L Amplitude",
							 _ntofbins2d, _tofmin, tofmax, nbs, 2*_trigScintAmplitudeMin, 2*_trigScintAmplitudeMax);

   // TOF0 R
  _histos2d["hRef_pbC_TrigScint0RC" + selTag] = new TH2D("hRef_pbC_TrigScint0RC" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 0R Charge",
							 nbs, _PbGChargeMin, _PbGChargeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_pbA_TrigScint0RC" + selTag] = new TH2D("hRef_pbA_TrigScint0RC" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 0R Charge",
							 nbs, 0., _PbGAmplitudeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_pbC_TrigScint0RA" + selTag] = new TH2D("hRef_pbC_TrigScint0RA" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 0R Amplitude",
							 nbs, _PbGChargeMin, _PbGChargeMax, nbs, 2*_trigScintAmplitudeMin, 2*_trigScintAmplitudeMax);
  _histos2d["hRef_pbA_TrigScint0RA" + selTag] = new TH2D("hRef_pbA_TrigScint0RA" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 0R Amplitude",
							 nbs, 0., _PbGAmplitudeMax, nbs,_trigScintChargeMin, _trigScintAmplitudeMax/4.);

  _histos2d["hRef_TOF_TrigScint0RC" + selTag] = new TH2D("hRef_TOF_TrigScint0RC" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 0R Charge",
							 _ntofbins2d, _tofmin, tofmax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_TOF_TrigScint0RA" + selTag] = new TH2D("hRef_TOF_TrigScint0RA" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 0R Amplitude",
							 _ntofbins2d, _tofmin, tofmax, nbs, 2*_trigScintAmplitudeMin, 2*_trigScintAmplitudeMax);

  // L-R
  _histos2d["hRef_TrigScint0RC_TrigScint0LC" + selTag] = new TH2D("hRef_TrigScint0RC_TrigScint0LC" + selTag, "; " + selTit + " Trig. scint. 0 R Charge; " + selTit + " Trig. scint. 0 L Charge", nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);


  // TOF1X
  _histos2d["hRef_pbC_TrigScint1C" + selTag] = new TH2D("hRef_pbC_TrigScint1C" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 1 Charge",
							nbs, _PbGChargeMin, _PbGChargeMax, nbs, 4*_trigScintChargeMin, 4*_trigScintChargeMax);
  _histos2d["hRef_pbA_TrigScint1C" + selTag] = new TH2D("hRef_pbA_TrigScint1C" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 1 Charge",
							nbs, 0., _PbGAmplitudeMax, nbs,4*_trigScintChargeMin, 4*_trigScintChargeMax);
  _histos2d["hRef_pbC_TrigScint1A" + selTag] = new TH2D("hRef_pbC_TrigScint1A" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 1 Amplitude",
							nbs, _PbGChargeMin, _PbGChargeMax, nbs, _trigScintAmplitudeMin, _trigScintAmplitudeMax/2.);
  _histos2d["hRef_pbA_TrigScint1A" + selTag] = new TH2D("hRef_pbA_TrigScint1A" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 1 Amplitude",
							nbs, 0., _PbGAmplitudeMax, nbs,_trigScintChargeMin, _trigScintAmplitudeMax/2.);

  _histos2d["hRef_TOF_TrigScint1C" + selTag] = new TH2D("hRef_TOF_TrigScint1C" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 1 Charge",
							_ntofbins2d, _tofmin, tofmax, nbs, 4*_trigScintChargeMin, 4*_trigScintChargeMax);
  _histos2d["hRef_TOF_TrigScint1A" + selTag] = new TH2D("hRef_TOF_TrigScint1A" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 1 Amplitude",
							_ntofbins2d, _tofmin, tofmax, nbs, _trigScintAmplitudeMin, _trigScintAmplitudeMax/2.);
  
  // TOF1 L
  _histos2d["hRef_pbC_TrigScint1LC" + selTag] = new TH2D("hRef_pbC_TrigScint1LC" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 1L Charge",
							 nbs, _PbGChargeMin, _PbGChargeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_pbA_TrigScint1LC" + selTag] = new TH2D("hRef_pbA_TrigScint1LC" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 1L Charge",
							 nbs, 0., _PbGAmplitudeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_pbC_TrigScint1LA" + selTag] = new TH2D("hRef_pbC_TrigScint1LA" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 1L Amplitude",
							 nbs, _PbGChargeMin, _PbGChargeMax, nbs, 2*_trigScintAmplitudeMin, 2*_trigScintAmplitudeMax);
  _histos2d["hRef_pbA_TrigScint1LA" + selTag] = new TH2D("hRef_pbA_TrigScint1LA" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 1L Amplitude",
							 nbs, 0., _PbGAmplitudeMax, nbs,_trigScintChargeMin, _trigScintAmplitudeMax/4.);

  _histos2d["hRef_TOF_TrigScint1LC" + selTag] = new TH2D("hRef_TOF_TrigScint1LC" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 1L Charge",
							 _ntofbins2d, _tofmin, tofmax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_TOF_TrigScint1LA" + selTag] = new TH2D("hRef_TOF_TrigScint1LA" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 1L Amplitude",
							 _ntofbins2d, _tofmin, tofmax, nbs, 2*_trigScintAmplitudeMin, 2*_trigScintAmplitudeMax);

  // TOF1 R
  _histos2d["hRef_pbC_TrigScint1RC" + selTag] = new TH2D("hRef_pbC_TrigScint1RC" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 1R Charge",
							 nbs, _PbGChargeMin, _PbGChargeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_pbA_TrigScint1RC" + selTag] = new TH2D("hRef_pbA_TrigScint1RC" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 1R Charge",
							 nbs, 0., _PbGAmplitudeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_pbC_TrigScint1RA" + selTag] = new TH2D("hRef_pbC_TrigScint1RA" + selTag, "; " + selTit + "  Pb-glass Charge; " + selTit + " Trig. scint. 1R Amplitude",
							 nbs, _PbGChargeMin, _PbGChargeMax, nbs, 2*_trigScintAmplitudeMin, 2*_trigScintAmplitudeMax);
  _histos2d["hRef_pbA_TrigScint1RA" + selTag] = new TH2D("hRef_pbA_TrigScint1RA" + selTag, "; " + selTit + "  Pb-glass Amplitude; " + selTit + " Trig. scint. 1R Amplitude",
							 nbs, 0., _PbGAmplitudeMax, nbs,_trigScintChargeMin, _trigScintAmplitudeMax/4.);

  _histos2d["hRef_TOF_TrigScint1RC" + selTag] = new TH2D("hRef_TOF_TrigScint1RC" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 1R Charge",
							 _ntofbins2d, _tofmin, tofmax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  _histos2d["hRef_TOF_TrigScint1RA" + selTag] = new TH2D("hRef_TOF_TrigScint1RA" + selTag, "; " + selTit + "  t_{1}-t_{0} [ns]; " + selTit + " Trig. scint. 1R Amplitude",
							 _ntofbins2d, _tofmin, tofmax, nbs, 2*_trigScintAmplitudeMin, 2*_trigScintAmplitudeMax);
  
  // L-R
  _histos2d["hRef_TrigScint1RC_TrigScint1LC" + selTag] = new TH2D("hRef_TrigScint1RC_TrigScint1LC" + selTag, "; " + selTit + " Trig. scint. 1 R Charge ; " + selTit + "  Trig. scint. 1L Charge",
								  nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  // L-L between trig scinti 0 and 1
  _histos2d["hRef_TrigScint0LC_TrigScint1LC" + selTag] = new TH2D("hRef_TrigScint0LC_TrigScint1LC" + selTag, "; " + selTit + " Trig. scint. 0 L Charge ; " + selTit + "  Trig. scint. 1L Charge",
								  nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);
  
  // R-R between trig scinti 0 and 1
  _histos2d["hRef_TrigScint0RC_TrigScint1RC" + selTag] = new TH2D("hRef_TrigScint0RC_TrigScint1RC" + selTag, "; " + selTit + " Trig. scint. 0 R Charge ; " + selTit + "  Trig. scint. 1R Charge",
								  nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax, nbs, 2*_trigScintChargeMin, 2*_trigScintChargeMax);

  // extrapolated map, charge-weighted
  double sf = 1.2;// 0.24;
  _histos2d["hRef_TrigScint0_Cweighted_xymap" + selTag] = new TH2D("hRef_TrigScint0_Cweighted_xymap" + selTag, "; " + selTit + " Trig. scint. 0 Ch.-weighted X [quasi mm];" + selTit + " Trig. scint. 0 Ch.-weighted Y [quasi mm];", 100, sf*_xL, sf*_xR, 100, sf*_yD, sf*_yU);
  _histos2d["hRef_TrigScint1_Cweighted_xymap" + selTag] = new TH2D("hRef_TrigScint1_Cweighted_xymap" + selTag, "; " + selTit + " Trig. scint. 1 Ch.-weighted X [quasi mm];" + selTit + " Trig. scint. 1 Ch.-weighted Y [quasi mm];", 100, sf*_xL, sf*_xR, 100, sf*_yD, sf*_yU);
  
  _histos2d["hRef_TrigScint0_Tweighted_xymap" + selTag] = new TH2D("hRef_TrigScint0_Tweighted_xymap" + selTag, "; " + selTit + " Trig. scint. 0 time-weighted X [quasi mm];" + selTit + " Trig. scint. 0 time-weighted Y [quasi mm];", 100, sf*_xL, sf*_xR, 100, sf*_yD, sf*_yU);
  _histos2d["hRef_TrigScint1_Tweighted_xymap" + selTag] = new TH2D("hRef_TrigScint1_Tweighted_xymap" + selTag, "; " + selTit + " Trig. scint. 1 time-weighted X [quasi mm];" + selTit + " Trig. scint. 1 time-weighted Y [quasi mm];", 100, sf*_xL, sf*_xR, 100, sf*_yD, sf*_yU);
  
  
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
  _histos2d["hRef_pbA_act23A"] = new TH2D("hRef_pbA_act23A", "; Pb-glass Amplitude ; (ACT2+ACT3)/2 Amplitude", 200, 0., _PbGAmplitudeMax, 400, 0., 2*_actAmplitudeMax);
  _histos2d["hRef_pbC_act23C"] = new TH2D("hRef_pbC_act23C", "; Pb-glass Charge ; (ACT2+ACT3)/2 Charge", 200, _actChargeMin, _actChargeMax, 400, 0., 2*_actAmplitudeMax);

  _histos2d["hRef_pbC_act23A"] = new TH2D("hRef_pbC_act23A", "; Pb-glass Charge ; (ACT2+ACT3)/2 Amplitude", 200, 0., _actAmplitudeMax, 400, 0., 2*_actAmplitudeMax);
  _histos2d["hRef_pbA_act23C"] = new TH2D("hRef_pbA_act23C", "; Pb-glass Charge ; (ACT2+ACT3)/2 Amplitude", 200, _PbGChargeMin, _PbGChargeMax, 400, 0., 2*_actAmplitudeMax);

  _histos2d["hRef_pbA_act0A"] = new TH2D("hRef_pbA_act0A", "; Pb-glass Amplitude ; ACT0 Amplitude", 200, 0., _PbGAmplitudeMax, 400, 0., ACT0Gain*_actAmplitudeMax);
  _histos2d["hRef_pbC_act0C"] = new TH2D("hRef_pbC_act0C", "; Pb-glass Charge ; ACT1 Charge", 200, _PbGChargeMin, _PbGChargeMax, 400, 0., ACT0Gain*_actChargeMax);
  _histos2d["hRef_pbA_act1A"] = new TH2D("hRef_pbA_act1A", "; Pb-glass Amplitude ; ACT1 Amplitude", 200, 0., _PbGAmplitudeMax, 400, 0., _actAmplitudeMax);
  _histos2d["hRef_pbC_PbG1C"] = new TH2D("hRef_pbC_PbG1C", "; Pb-glass Charge ; PBG1 Charge", 200, _actChargeMin, _PbGChargeMax, 400, 0., _PbGChargeMax);
  _histos2d["hRef_pbA_act1C"] = new TH2D("hRef_pbA_act1C", "; Pb-glass Amplitude ; ACT1 Charge", 200, 0., _PbGAmplitudeMax, 400,_actChargeMin, _actChargeMax);

  _histos2d["hRef_pbC_act1C"] = new TH2D("hRef_pbC_act1C", "; Pb-glass Charge ; ACT1 Charge", 200, _actChargeMin, _actChargeMax, 400, 0., _actAmplitudeMax);
  // act1cuts
  _histos2d["hRef_pbA_act1A_act1cuts"] = new TH2D("hRef_pbA_act1A_act1cuts", "; Pb-glass Amplitude ; ACT1 Amplitude", 200, 0., _PbGAmplitudeMax, 400, 0., _actAmplitudeMax);

  _outFile -> cd("../");

  // 5.3.2024
  // Trigger scintillators (unfortunatelly labelled as TOF all through out the code;-)
  // amplitude and charge vs tof:
  // all PID:
  this -> InitTrigScintHistos("TrigScint", "", "");
  // also for proton-like selection
  this -> InitTrigScintHistos("TrigScint_p", "_p-like", " p-like");
  this -> InitTrigScintHistos("TrigScint_D", "_D-like", " D-like");
  this -> InitTrigScintHistos("TrigScint_T", "_T-like", " T-like");
  this -> InitTrigScintHistos("TrigScint_e", "_e-like", " e-like");
  this -> InitTrigScintHistos("TrigScint_mu", "_mu-like", " mu-like");
  this -> InitTrigScintHistos("TrigScint_pi", "_pi-like", " pi-like");


  _outFile -> cd("Charged");

  // L-R studies:
  // 17.11.2023
  double act0sf = 2.;
  _histos2d["hRef_act0LA_act0RA_nonZero"] = new TH2D("hRef_act0LA_act0RA_nonZero", "; ACT0L Amplitude;ACT0R Amplitude", 200, 0., _actAmplitudeMax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_act1LA_act1RA_nonZero"] = new TH2D("hRef_act1LA_act1RA_nonZero", "; ACT1L Amplitude;ACT1R Amplitude", 200, 0., _actAmplitudeMax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_act2LA_act2RA_nonZero"] = new TH2D("hRef_act2LA_act2RA_nonZero", "; ACT2L Amplitude;ACT2R Amplitude", 200, 0., _actAmplitudeMax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_act3LA_act3RA_nonZero"] = new TH2D("hRef_act3LA_act3RA_nonZero", "; ACT3L Amplitude;ACT3R Amplitude", 200, 0., _actAmplitudeMax, 200, 0., _actAmplitudeMax);

  _histos2d["hRef_act0LA_act1LA_nonZero"] = new TH2D("hRef_act0LA_act1LA_nonZero", "; ACT0L Amplitude;ACT1L Amplitude", 200, 0., act0sf*_actAmplitudeMax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_act0RA_act1RA_nonZero"] = new TH2D("hRef_act0RA_act1RA_nonZero", "; ACT0R Amplitude;ACT1R Amplitude", 200, 0., act0sf*_actAmplitudeMax, 200, 0., _actAmplitudeMax);

  _histos2d["hRef_act0LC_act0RC_nonZero"] = new TH2D("hRef_act0LC_act0RC_nonZero", "; ACT0L Charge;ACT0R Charge", 200, 0., act0sf*_actChargeMax, 200, 0., act0sf*_actChargeMax);
  _histos2d["hRef_act1LC_act1RC_nonZero"] = new TH2D("hRef_act1LC_act1RC_nonZero", "; ACT1L Charge;ACT1R Charge", 200, 0., _actChargeMax, 200, 0., _actChargeMax);
  _histos2d["hRef_act2LC_act2RC_nonZero"] = new TH2D("hRef_act2LC_act2RC_nonZero", "; ACT2L Charge;ACT2R Charge", 200, 0., _actChargeMax, 200, 0., _actChargeMax);
  _histos2d["hRef_act3LC_act3RC_nonZero"] = new TH2D("hRef_act3LC_act3RC_nonZero", "; ACT3L Charge;ACT3R Charge", 200, 0., _actChargeMax, 200, 0., _actChargeMax);

  // not good given the different gain...
  _histos2d["hRef_act0LC_act1LC_nonZero"] = new TH2D("hRef_act0LC_act1LC_nonZero", "; ACT0L Charge;ACT1L Charge", 200, 0., act0sf*_actChargeMax, 200, 0., _actChargeMax);
  _histos2d["hRef_act0RC_act1RC_nonZero"] = new TH2D("hRef_act0RC_act1RC_nonZero", "; ACT0R Charge;ACT1R Charge", 200, 0., act0sf*_actChargeMax, 200, 0., _actChargeMax);

  // diff, charges only
  double dsf = 0.75;
  _histos1d["hRef_act0LC_minus_act0RC_nonZero"] = new TH1D("hRef_act0LC_minus_act0RC_nonZero", "; ACT0L Charge - ACT0R Charge", 200, -act0sf*_actChargeMax*dsf, act0sf*_actChargeMax*dsf);
  _histos1d["hRef_act1LC_minus_act1RC_nonZero"] = new TH1D("hRef_act1LC_minus_act1RC_nonZero", "; ACT1L Charge - ACT1R Charge", 200, -_actChargeMax*dsf, _actChargeMax*dsf);
  _histos1d["hRef_act2LC_minus_act2RC_nonZero"] = new TH1D("hRef_act2LC_minus_act2RC_nonZero", "; ACT2L Charge - ACT2R Charge", 200, -_actChargeMax*dsf, _actChargeMax*dsf);
  _histos1d["hRef_act3LC_minus_act3RC_nonZero"] = new TH1D("hRef_act3LC_minus_act3RC_nonZero", "; ACT3L Charge - ACT3R Charge", 200, -_actChargeMax*dsf, _actChargeMax*dsf);

  // 1D x map, in mm:
  double xrange = 70; // mm
  _histos1d["hRef_act0Xmap_nonZero"] = new TH1D("hRef_act0Xmap_nonZero", "; ACT0 weighted X [quasi mm]", 100, -xrange, xrange );
  _histos1d["hRef_act1Xmap_nonZero"] = new TH1D("hRef_act1Xmap_nonZero", "; ACT1 weighted X [quasi mm]", 100, -xrange, xrange );
  _histos1d["hRef_act2Xmap_nonZero"] = new TH1D("hRef_act2Xmap_nonZero", "; ACT2 weighted X [quasi mm]", 100, -xrange, xrange );
  _histos1d["hRef_act3Xmap_nonZero"] = new TH1D("hRef_act3Xmap_nonZero", "; ACT3 weighted X [quasi mm]", 100, -xrange, xrange );
  
  // not good given the different gain...
  _histos1d["hRef_act0LC_minus_act1LC_nonZero"] = new TH1D("hRef_act0LC_minus_act1LC_nonZero", "; ACT0L Charge - ACT1L Charge", 200, -_actChargeMax*dsf, _actChargeMax*dsf);
  _histos1d["hRef_act0RC_minus_act1RC_nonZero"] = new TH1D("hRef_act0RC_minus_act1RC_nonZero", "; ACT0R Charge - ACT1R Charge", 200, -_actChargeMax*dsf, _actChargeMax*dsf);

  // (ACT2+ACT3)/2 vs TOF plots
  _histos2d["hRef_TOFACT23A"] = new TH2D("hRef_TOFACT23A", "; t_{1}-t_{0} [ns]; (ACT2+ACT3)/2 Amplitude", _ntofbins2d, _tofmin, _tofmax, 200, 0., 2*_actAmplitudeMax);
  _histos2d["hRef_TOFACT23C"] = new TH2D("hRef_TOFACT23C", "; t_{1}-t_{0} [ns]; (ACT2+ACT3)/2 Charge", _ntofbins2d, _tofmin, _tofmax, 200, _actChargeMin, 2*_actChargeMax);

  // also ACT 0 and 1, separately:
  /* seems they were already defined below...
  _histos2d["hRef_TOFACT0A"] = new TH2D("hRef_TOFACT0A", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", _ntofbins2d, _tofmin, _tofmax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_TOFACT1A"] = new TH2D("hRef_TOFACT1A", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", _ntofbins2d, _tofmin, _tofmax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_TOFACT0C"] = new TH2D("hRef_TOFACT0C", "; t_{1}-t_{0} [ns]; ACT0 Charge", _ntofbins2d, _tofmin, _tofmax, 200, _actChargeMin, _actChargeMax);
  _histos2d["hRef_TOFACT1C"] = new TH2D("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", _ntofbins2d, _tofmin, _tofmax, 200, _actChargeMin, _actChargeMax);
  */

  //TOF vs Pb-glass plots
  _histos2d["hRef_PbATOF"] = new TH2D("hRef_PbATOF", "; Pb-glass Amplitude; t_{1}-t_{0} [ns]", 200, 0., _PbGAmplitudeMax, _ntofbins2d, _tofmin, _tofmax);
  _histos2d["hRef_PbCTOF"] = new TH2D("hRef_PbCTOF", "; Pb-glass Charge; t_{1}-t_{0} [ns]", 200, _PbGChargeMin, _PbGChargeMax, _ntofbins2d, _tofmin, _tofmax);
  _histos2d["hRef_TOFPbA"] = new TH2D("hRef_TOFPbA", "; t_{1}-t_{0} [ns]; Pb-glass Amplitude", _ntofbins2d, _tofmin, _tofmax, 200, 0., _PbGAmplitudeMax);
  _histos2d["hRef_TOFPbC"] = new TH2D("hRef_TOFPbC", "; t_{1}-t_{0} [ns]; Pb-glass Charge", _ntofbins2d, _tofmin, _tofmax, 200, 0., _PbGChargeMax);
  
  //acraplet - investigate "weird electrons"
  _histos2d["hHC0AHC1A"] = new TH2D("hweirdE_HC0AHC1A", "; Hole Counter 0 Amplitude; Hole Counter 1 Amplitude", 200, 0., 1000, 200, 0., 1000.);
  _histos2d["hHC0CHC1C"] = new TH2D("hweirdE_HC0CHC1C", "; Hole Counter 0 Charge; Hole Counter 1 Charge", 200, 0., 1., 200, 0., 1.);

  // no cuts
  _histos2d["hRef_TOFACT0A"] = new TH2D("hRef_TOFACT0A", "; t_{1}-t_{0} [ns]; ACT0 Amplitude", _ntofbins2d, _tofmin, _tofmax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_TOFACT1A"] = new TH2D("hRef_TOFACT1A", "; t_{1}-t_{0} [ns]; ACT1 Amplitude", _ntofbins2d, _tofmin, _tofmax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_TOFACT2A"] = new TH2D("hRef_TOFACT2A", "; t_{1}-t_{0} [ns]; ACT2 Amplitude", _ntofbins2d, _tofmin, _tofmax, 200, 0., _actAmplitudeMax);
  _histos2d["hRef_TOFACT3A"] = new TH2D("hRef_TOFACT3A", "; t_{1}-t_{0} [ns]; ACT3 Amplitude", _ntofbins2d, _tofmin, _tofmax, 200, 0., _actAmplitudeMax);

  _histos2d["hRef_TOFACT0C"] = new TH2D("hRef_TOFACT0C", "; t_{1}-t_{0} [ns]; ACT0 Charge", _ntofbins2d, _tofmin, _tofmax, 200, _actChargeMin, _actChargeMax);
  _histos2d["hRef_TOFACT1C"] = new TH2D("hRef_TOFACT1C", "; t_{1}-t_{0} [ns]; ACT1 Charge", _ntofbins2d, _tofmin, _tofmax, 200, _actChargeMin, _actChargeMax);
  _histos2d["hRef_TOFACT2C"] = new TH2D("hRef_TOFACT2C", "; t_{1}-t_{0} [ns]; ACT2 Charge", _ntofbins2d, _tofmin, _tofmax, 200, _actChargeMin, _actChargeMax);
  _histos2d["hRef_TOFACT3C"] = new TH2D("hRef_TOFACT3C", "; t_{1}-t_{0} [ns]; ACT3 Charge", _ntofbins2d, _tofmin, _tofmax, 200, _actChargeMin, _actChargeMax);


  // ACT2+ACT3 cut
  _histos1d["hTOF_act2act3cut"] = new TH1D("hTOF_act2act3cut", "; t_{1}-t_{0} [ns];", 120, _tofmin, _tofmax);

  // 2D ACT charges
  _histos2d["hACT2CACT1C"] = new TH2D("hRef_ACT2CACT1C", "; ACT2 Charge; ACT1 Charge", 200, _actChargeMin, _actChargeMax, 200, _actChargeMin, _actChargeMax);
  _histos2d["hACT3CACT2C"] = new TH2D("hRef_ACT3CACT2C", "; ACT3 Charge; ACT2 Charge", 200, _actChargeMin, _actChargeMax, 200, _actChargeMin, _actChargeMax);
  _histos2d["hACT1CACT3C"] = new TH2D("hRef_ACT1CACT3C", "; ACT1 Charge; ACT3 Charge", 200, _actChargeMin, _actChargeMax, 200, _actChargeMin, _actChargeMax);

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
  
  _histos2d["hnPeaksACT23vsToF"] = new TH2D("hnPeaksACT23vsToF", "hnPeaksACT23vsToF;t_{TOF};<n_{Peaks}^{ACT23}>", _ntofbins2d/4, _tofmin, _tofmax, nbn, n1, n2);
  _histos2d["hnPeaksACT23vsToFlow"] = new TH2D("hnPeaksACT23vsToFlow", "hnPeaksACT23vsToF;t_{TOF};<n_{Peaks}^{ACT23}>", _ntofbins2d/4, _tofminlow, _tofmaxlow, nbn, n1, n2);
  _histos2d["hnPeaksToFvsToF"] = new TH2D("hnPeaksToFvsToF", "hnPeaksToFvsToF;t_{TOF};<n_{Peaks}^{ToF}>", _ntofbins2d/4, _tofmin, _tofmax, nbn, n1, n2);
  _histos2d["hnPeaksToFvsToFlow"] = new TH2D("hnPeaksToFvsToFlow", "hnPeaksToFvsToF;t_{TOF};<n_{Peaks}^{ToF}>", _ntofbins2d/4, _tofminlow, _tofmaxlow, nbn, n1, n2);
  _histos2d["hnPeaksACT23vsLeadGlassA"] = new TH2D("hnPeaksACT23vsLeadGlassA", "hnPeaksACT23vsLeadGlassA;lead glass A;<n_{Peaks}^{ACT23}>", 100,  0., _actAmplitudeMax/2., nbn, n1, n2);
  _histos2d["hnPeaksToFvsLeadGlassA"] = new TH2D("hnPeaksToFvsLeadGlassA", "hnPeaksToFvsLeadGlassA;lead glass A;<n_{Peaks}^{ToF}>", 100,  0., _actAmplitudeMax/2., nbn, n1, n2);
  n1 = 0.;
  n2 = 10.;
  _histos2d["hnPeaksLeadGlassvsLeadGlassA"] = new TH2D("hnPeaksLeadGlassvsLeadGlassA", "hnPeaksLeadGlassvsLeadGlassA;lead glass A;n_{Peaks}^{Pb}", 100,  0., _actAmplitudeMax/2., int(n2-n1), n1, n2);
  _outFile -> cd("../");
  
  cout << "done" << endl;
  
    
} // InitHistos



// ______________________________________________________________

int MakeAllDataPlots::getHighestPeakIndex(channelReadClass *reader, bool useCharges)
 {
   int imax = -1;
   double maxA = -999;
   double a = -999;
   int nn = _useWindowIntCharge ? reader -> nWindowPeaks : reader -> nPeaks;
   for (int ipeak = 0; ipeak < nn; ++ipeak) {
     if (useCharges) {
       if (_useWindowIntCharge) {
	 // a = reader -> WindowIntCharge[ipeak];
	 a = reader -> WindowIntPE[ipeak];
       } else
	 a = reader -> IntCharge[ipeak];
     } else {
       a = reader -> PeakVoltage[ipeak];       
     }
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

    // READ!
    this -> ReadChannels();
    
    // peak cuts on demand
    if (!_isHodoscopeRun) {
      if (! this -> PassedPeakCuts())
	continue;
    } // not hodoscope run

    // FILL!
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

  cout << "End of event loop!" << endl;
} // Loop


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
  _onePeakInPbGlass = (_NPeaksC["PbGlass"] == 1);
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

  if (_peakMode == "i" && (  (_Charges["TOF10"] < 0.04 && _Charges["TOF10"] > 0.03) || (_Amplitudes["TOF00"] > 1.43 && _Amplitudes["TOF00"] < 1.435) ) )
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

      // read amplitudes:
      _PeakIDA[chname] = getHighestPeakIndex(_readerMap[chname], false);
      int ipeak = _PeakIDA[chname];
      _NPeaksA[chname] = _readerMap[chname] -> nPeaks;
      if ( ipeak >= 0 && ipeak < _readerMap[chname] -> nPeaks) {
	_Amplitudes[chname]  = _readerMap[chname] -> PeakVoltage[ipeak];
      } else {
	_Amplitudes[chname]  = 0.;
      }

      // read window charges:
      _PeakIDC[chname] = getHighestPeakIndex(_readerMap[chname], true);
      ipeak = _PeakIDC[chname];
      _NPeaksC[chname] = _readerMap[chname] -> nPeaks;
      if ( ipeak >= 0 && ipeak < _readerMap[chname] -> nPeaks) {
	if (_useWindowIntCharge) { // && !chname.Contains("TOF"))
	  // preferred, to compare PMTs
	  if (chname != "PbGlass")
	    _Charges[chname]     = _readerMap[chname] -> WindowIntPE[ipeak];// ?!?!?!
	  else
	    _Charges[chname]     = _readerMap[chname] -> WholeWaveformInt;// ?!?!?!
	  // NOT preferred!:
	  // _Charges[chname]     = _readerMap[chname] -> WindowIntCharge[ipeak]; 
	  _SignalTimes[chname] = _readerMap[chname] -> SignalTimeCorrected[ipeak];
	} else {
	  _Charges[chname]     = _readerMap[chname] -> IntCharge[ipeak];
	  _SignalTimes[chname] = _readerMap[chname] -> SignalTime[ipeak];
	}
      } else {
	_Charges[chname]     = 0.;
	_SignalTimes[chname] = 0.;
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
      
      int ipeak = _PeakIDA[chname];
      _hnPeaksA.at(ich).Fill(_reader[ich] -> nPeaks);
      if ( ipeak >= 0 && ipeak < _readerMap[chname] -> nPeaks) {
	// can be simplified using the above maps
	_hVoltage.at(ich).Fill(_Amplitudes[chname]);
	_hPedestalSigma.at(ich).Fill(_reader[ich] -> PedestalSigma);
      }

      ipeak = _PeakIDC[chname];
      int nn =  _useWindowIntCharge ? _readerMap[chname] -> nWindowPeaks : _readerMap[chname] -> nPeaks;
      _hnPeaksC.at(ich).Fill(nn);
      if ( ipeak >= 0 && ipeak < nn) {
	// can be simplified using the above maps
	_hCharge.at(ich).Fill(_Charges[chname]);
	_hTime.at(ich).Fill(_SignalTimes[chname]);
      }
      
      // depricated:
	/*
	// histograms over all channels
	_hCharge.at(ich).Fill(_reader[ich] -> IntCharge[ipeak]);
	_hVoltage.at(ich).Fill(_reader[ich] -> PeakVoltage[ipeak]);
	_hTime.at(ich).Fill(_reader[ich] -> SignalTime[ipeak]);
	//hNbPeaks.at(ich).Fill(_reader[ich] -> nPeaks);
	_hPedestalSigma.at(ich).Fill(_reader[ich] -> PedestalSigma);
	_hnPeaks.at(ich).Fill(_reader[ich] -> nPeaks);
	*/

    } // channels
}

// ______________________________________________________________

void MakeAllDataPlots::FillTofHistos()
{

  // TOF trigger scintilators

    double t00 = _SignalTimes["TOF00"];
    double t01 = _SignalTimes["TOF01"];
    double t02 = _SignalTimes["TOF02"];
    double t03 = _SignalTimes["TOF03"];

    double t10 = _SignalTimes["TOF10"];
    double t11 = _SignalTimes["TOF11"];
    double t12 = _SignalTimes["TOF12"];
    double t13 = _SignalTimes["TOF13"];

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


void MakeAllDataPlots::ComputeChargesAndAmplitudes() {

  if (_debug)      cout << "charged a" << endl;

    // ACTs
    _act0c = _Charges["ACT0L"] + _Charges["ACT0R"];
    _act1c = _Charges["ACT1L"] + _Charges["ACT1R"];
    _act2c = _Charges["ACT2L"] + _Charges["ACT2R"];
    _act3c = _Charges["ACT3L"] + _Charges["ACT3R"];

    _act0a = _Amplitudes["ACT0L"] + _Amplitudes["ACT0R"];
    _act1a = _Amplitudes["ACT1L"] + _Amplitudes["ACT1R"];
    _act2a = _Amplitudes["ACT2L"] + _Amplitudes["ACT2R"];
    _act3a = _Amplitudes["ACT3L"] + _Amplitudes["ACT3R"];

    _act23aAver = (_act2a + _act3a) / 2.;
    _act23cAver = (_act2c + _act3c) / 2.;

    // hole counters and lead glass
    
    _hc0c = _Charges["Hole0"];
    _hc0a = _Amplitudes["Hole0"];

    _hc1c = _Charges["Hole1"];
    _hc1a = _Amplitudes["Hole1"];

    _pbc = _Charges["PbGlass"];
    _pba = _Amplitudes["PbGlass"];

    _trigScintC = _Charges["TOF00"] + _Charges["TOF01"] + _Charges["TOF02"] + _Charges["TOF03"] + _Charges["TOF10"] + _Charges["TOF11"] + _Charges["TOF12"] + _Charges["TOF13"];
    _trigScintA = _Amplitudes["TOF00"] + _Amplitudes["TOF01"] + _Amplitudes["TOF02"] + _Amplitudes["TOF03"] + _Amplitudes["TOF10"] + _Amplitudes["TOF11"] + _Amplitudes["TOF12"] + _Amplitudes["TOF13"];


    // Bruno, 7.3.2024:
    // 0 and 1 are up
    // 2 and 3 are down
    // 'TOF' trigger scintillator 0:
    // Channels 1 and 3 are on the right side looking downstream and 0 & 2 on the left side.
    // 'TOF' trigger scintillator 1:
    // Channels 0 and 2 are on the right side looking downstream and 1 & 3 on the left side.
    
    _trigScint0C = _Charges["TOF00"] + _Charges["TOF01"] + _Charges["TOF02"] + _Charges["TOF03"];
    _trigScint0A = _Amplitudes["TOF00"] + _Amplitudes["TOF01"] + _Amplitudes["TOF02"] + _Amplitudes["TOF03"];
    _trigScint0LC = _Charges["TOF00"] + _Charges["TOF02"];
    _trigScint0LA = _Amplitudes["TOF00"] + _Amplitudes["TOF02"];
    _trigScint0RC = _Charges["TOF01"] + _Charges["TOF03"];
    _trigScint0RA = _Amplitudes["TOF01"] + _Amplitudes["TOF03"];

    _trigScint1C = _Charges["TOF10"] + _Charges["TOF11"] + _Charges["TOF12"] + _Charges["TOF13"];
    _trigScint1A = _Amplitudes["TOF10"] + _Amplitudes["TOF11"] + _Amplitudes["TOF12"] + _Amplitudes["TOF13"];
    _trigScint1LC = _Charges["TOF11"] + _Charges["TOF13"];
    _trigScint1LA = _Amplitudes["TOF11"] + _Amplitudes["TOF13"];
    _trigScint1RC = _Charges["TOF10"] + _Charges["TOF12"];
    _trigScint1RA = _Amplitudes["TOF10"] + _Amplitudes["TOF12"];


    // Charge-weighted hit positions in T0 and T1
    _trigScint0_CweightedX =  ( _xL*_Charges["TOF00"] + _xL*_Charges["TOF02"] + _xR*_Charges["TOF01"] + _xR*_Charges["TOF03"] ) / _trigScint0C;
    _trigScint1_CweightedX =  ( _xL*_Charges["TOF11"] + _xL*_Charges["TOF13"] + _xR*_Charges["TOF10"] + _xR*_Charges["TOF12"] ) / _trigScint1C; 
    _trigScint0_CweightedY =  ( _yU*_Charges["TOF00"] + _yU*_Charges["TOF01"] +  _yD*_Charges["TOF02"] + _yD*_Charges["TOF03"] ) / _trigScint0C;
    _trigScint1_CweightedY =  ( _yU*_Charges["TOF10"] + _yU*_Charges["TOF11"] +  _yD*_Charges["TOF12"] + _yD*_Charges["TOF13"] ) / _trigScint1C;
    
    // Time-weighted hit positions in T0 and T1
    
    double trigScint0SumTime = _SignalTimes["TOF00"] + _SignalTimes["TOF01"] + _SignalTimes["TOF02"] + _SignalTimes["TOF03"];
    double trigScint1SumTime = _SignalTimes["TOF10"] + _SignalTimes["TOF11"] + _SignalTimes["TOF12"] + _SignalTimes["TOF13"];
    
    _trigScint0_TweightedX =  ( _xL*_SignalTimes["TOF00"] + _xL*_SignalTimes["TOF02"] + _xR*_SignalTimes["TOF01"] + _xR*_SignalTimes["TOF03"] ) / trigScint0SumTime;
    _trigScint1_TweightedX =  ( _xL*_SignalTimes["TOF11"] + _xL*_SignalTimes["TOF13"] + _xR*_SignalTimes["TOF10"] + _xR*_SignalTimes["TOF12"] ) / trigScint1SumTime;
    _trigScint0_TweightedY =  ( _yU*_SignalTimes["TOF00"] + _yU*_SignalTimes["TOF01"] + _yD*_SignalTimes["TOF02"] + _yD*_SignalTimes["TOF03"] ) / trigScint0SumTime;
    _trigScint1_TweightedY =  ( _yU*_SignalTimes["TOF10"] + _yU*_SignalTimes["TOF11"] + _yD*_SignalTimes["TOF12"] + _yD*_SignalTimes["TOF13"] ) / trigScint1SumTime;
    

      
}

// ______________________________________________________________

void MakeAllDataPlots::FillTrigScintHistos(TString selTag)
{
    _histos2d["hRef_pbC_TrigScintC" + selTag]->Fill(_pbc, _trigScintC);
    _histos2d["hRef_pbA_TrigScintC" + selTag]->Fill(_pba, _trigScintC);
    _histos2d["hRef_pbC_TrigScintA" + selTag]->Fill(_pbc, _trigScintA);
    _histos2d["hRef_pbA_TrigScintA" + selTag]->Fill(_pba, _trigScintA);
    
    _histos2d["hRef_TOF_TrigScintC" + selTag]->Fill(_tof, _trigScintC);
    _histos2d["hRef_TOF_TrigScintA" + selTag]->Fill(_tof, _trigScintA);

    // TOF0X trigger scintillators:
    _histos2d["hRef_pbC_TrigScint0C" + selTag]->Fill(_pbc, _trigScint0C);
    _histos2d["hRef_pbA_TrigScint0C" + selTag]->Fill(_pba, _trigScint0C);
    _histos2d["hRef_pbC_TrigScint0A" + selTag]->Fill(_pbc, _trigScint0A);
    _histos2d["hRef_pbA_TrigScint0A" + selTag]->Fill(_pba, _trigScint0A);
    
    _histos2d["hRef_TOF_TrigScint0C" + selTag]->Fill(_tof, _trigScint0C);
    _histos2d["hRef_TOF_TrigScint0A" + selTag]->Fill(_tof, _trigScint0A);

    // TOF0 L trigger scintillators:
    _histos2d["hRef_pbC_TrigScint0LC" + selTag]->Fill(_pbc, _trigScint0LC);
    _histos2d["hRef_pbA_TrigScint0LC" + selTag]->Fill(_pba, _trigScint0LC);
    _histos2d["hRef_pbC_TrigScint0LA" + selTag]->Fill(_pbc, _trigScint0LA);
    _histos2d["hRef_pbA_TrigScint0LA" + selTag]->Fill(_pba, _trigScint0LA);
    
    _histos2d["hRef_TOF_TrigScint0LC" + selTag]->Fill(_tof, _trigScint0LC);
    _histos2d["hRef_TOF_TrigScint0LA" + selTag]->Fill(_tof, _trigScint0LA);
    
    // TOF0 R trigger scintillators:
    _histos2d["hRef_pbC_TrigScint0RC" + selTag]->Fill(_pbc, _trigScint0RC);
    _histos2d["hRef_pbA_TrigScint0RC" + selTag]->Fill(_pba, _trigScint0RC);
    _histos2d["hRef_pbC_TrigScint0RA" + selTag]->Fill(_pbc, _trigScint0RA);
    _histos2d["hRef_pbA_TrigScint0RA" + selTag]->Fill(_pba, _trigScint0RA);
    
    _histos2d["hRef_TOF_TrigScint0RC" + selTag]->Fill(_tof, _trigScint0RC);
    _histos2d["hRef_TOF_TrigScint0RA" + selTag]->Fill(_tof, _trigScint0RA);

    // L-R:
    _histos2d["hRef_TrigScint0RC_TrigScint0LC" + selTag]->Fill(_trigScint0RC, _trigScint0LC);

    // TOF1X trigger scintillators:
    _histos2d["hRef_pbC_TrigScint1C" + selTag]->Fill(_pbc, _trigScint1C);
    _histos2d["hRef_pbA_TrigScint1C" + selTag]->Fill(_pba, _trigScint1C);
    _histos2d["hRef_pbC_TrigScint1A" + selTag]->Fill(_pbc, _trigScint1A);
    _histos2d["hRef_pbA_TrigScint1A" + selTag]->Fill(_pba, _trigScint1A);
    
    _histos2d["hRef_TOF_TrigScint1C" + selTag]->Fill(_tof, _trigScint1C);
    _histos2d["hRef_TOF_TrigScint1A" + selTag]->Fill(_tof, _trigScint1A);

    // TOF1 L trigger scintillators:
    _histos2d["hRef_pbC_TrigScint1LC" + selTag]->Fill(_pbc, _trigScint1LC);
    _histos2d["hRef_pbA_TrigScint1LC" + selTag]->Fill(_pba, _trigScint1LC);
    _histos2d["hRef_pbC_TrigScint1LA" + selTag]->Fill(_pbc, _trigScint1LA);
    _histos2d["hRef_pbA_TrigScint1LA" + selTag]->Fill(_pba, _trigScint1LA);
    
    _histos2d["hRef_TOF_TrigScint1LC" + selTag]->Fill(_tof, _trigScint1LC);
    _histos2d["hRef_TOF_TrigScint1LA" + selTag]->Fill(_tof, _trigScint1LA);

    // TOF1 R trigger scintillators:
    _histos2d["hRef_pbC_TrigScint1RC" + selTag]->Fill(_pbc, _trigScint1RC);
    _histos2d["hRef_pbA_TrigScint1RC" + selTag]->Fill(_pba, _trigScint1RC);
    _histos2d["hRef_pbC_TrigScint1RA" + selTag]->Fill(_pbc, _trigScint1RA);
    _histos2d["hRef_pbA_TrigScint1RA" + selTag]->Fill(_pba, _trigScint1RA);
    
    _histos2d["hRef_TOF_TrigScint1RC" + selTag]->Fill(_tof, _trigScint1RC);
    _histos2d["hRef_TOF_TrigScint1RA" + selTag]->Fill(_tof, _trigScint1RA);
    
    // L-R:
    _histos2d["hRef_TrigScint1RC_TrigScint1LC" + selTag]->Fill(_trigScint1RC, _trigScint1LC);

    // L-L and R-R
    _histos2d["hRef_TrigScint0LC_TrigScint1LC" + selTag]->Fill(_trigScint0LC, _trigScint1LC);
    _histos2d["hRef_TrigScint0RC_TrigScint1RC" + selTag]->Fill(_trigScint0RC, _trigScint1RC);

    // the interpolated weoghted hits positions
    // weighted by the charges
     _histos2d["hRef_TrigScint0_Cweighted_xymap" + selTag]-> Fill(_trigScint0_CweightedX,  _trigScint0_CweightedY);
     _histos2d["hRef_TrigScint1_Cweighted_xymap" + selTag]-> Fill(_trigScint1_CweightedX,  _trigScint1_CweightedY);
     // weighted by the time
     _histos2d["hRef_TrigScint0_Tweighted_xymap" + selTag]-> Fill(_trigScint0_TweightedX,  _trigScint0_TweightedY);
     _histos2d["hRef_TrigScint1_Tweighted_xymap" + selTag]-> Fill(_trigScint1_TweightedX,  _trigScint1_TweightedY);
    
}

// ______________________________________________________________


void MakeAllDataPlots::FillChargedHistos()
{


  this -> ComputeChargesAndAmplitudes();

    // cout << " trigScintA=" << trigScintA << " trigScintC=" << trigScintC << endl;

    // HACK!
    //    pba = _readerMap["PbGlass"] -> PeakVoltage[0];

    // 17.11.2023:
    // A's are zero if no peaks;)
    if (_Amplitudes["ACT0L"] > 0.05 &&  _Amplitudes["ACT0R"] > 0.05) {
      _histos2d["hRef_act0LA_act0RA_nonZero"]->Fill(_Amplitudes["ACT0L"], _Amplitudes["ACT0R"] );
    }
    if (_Amplitudes["ACT1L"] > 0.01 &&  _Amplitudes["ACT1R"] > 0.01) {
      _histos2d["hRef_act1LA_act1RA_nonZero"]->Fill(_Amplitudes["ACT1L"], _Amplitudes["ACT1R"] );
    }
    if (_Amplitudes["ACT2L"] > 0.01 &&  _Amplitudes["ACT2R"] > 0.01) {
      _histos2d["hRef_act2LA_act2RA_nonZero"]->Fill(_Amplitudes["ACT2L"], _Amplitudes["ACT2R"] );
    }
    if (_Amplitudes["ACT3L"] > 0.01 &&  _Amplitudes["ACT3R"] > 0.01) {
      _histos2d["hRef_act3LA_act3RA_nonZero"]->Fill(_Amplitudes["ACT3L"], _Amplitudes["ACT3R"] );
    }

    if (_Amplitudes["ACT0L"] > 0.05 &&  _Amplitudes["ACT1L"] > 0.01) {
      _histos2d["hRef_act0LA_act1LA_nonZero"]->Fill(_Amplitudes["ACT0L"], _Amplitudes["ACT1L"] );
    }
    if (_Amplitudes["ACT0R"] > 0.05 &&  _Amplitudes["ACT1R"] > 0.01) {
      _histos2d["hRef_act0RA_act1RA_nonZero"]->Fill(_Amplitudes["ACT0R"], _Amplitudes["ACT1R"] );
    }


    
    _act0x = 0;
    _act1x = 0;
    _act2x = 0;
    _act3x = 0;
    double minCharge = 0.025; // 0.02
    
    // charges and their diffs
    if (_Charges["ACT0L"] > minCharge &&  _Charges["ACT0R"] > minCharge) {
      _histos2d["hRef_act0LC_act0RC_nonZero"]->Fill(_Charges["ACT0L"], _Charges["ACT0R"] );
      _histos1d["hRef_act0LC_minus_act0RC_nonZero"]->Fill(_Charges["ACT0L"] - _Charges["ACT0R"] );
      _act0x = (-_ACTwidth/2.*_Charges["ACT0L"] + _ACTwidth/2.*_Charges["ACT0R"]) / (_Charges["ACT0L"] + _Charges["ACT0R"]);
      _histos1d["hRef_act0Xmap_nonZero"]->Fill(_act0x);
    }
    if (_Charges["ACT1L"] > minCharge &&  _Charges["ACT1R"] > minCharge) {
      _histos2d["hRef_act1LC_act1RC_nonZero"]->Fill(_Charges["ACT1L"], _Charges["ACT1R"] );
      _histos1d["hRef_act1LC_minus_act1RC_nonZero"]->Fill(_Charges["ACT1L"] - _Charges["ACT1R"] );
      _act1x = (-_ACTwidth/2.*_Charges["ACT1L"] + _ACTwidth/2.*_Charges["ACT1R"]) / (_Charges["ACT1L"] + _Charges["ACT1R"]);
      _histos1d["hRef_act1Xmap_nonZero"]->Fill(_act1x);
	    
    }
    if (_Charges["ACT2L"] > minCharge &&  _Charges["ACT2R"] > minCharge) {
      _histos2d["hRef_act2LC_act2RC_nonZero"]->Fill(_Charges["ACT2L"], _Charges["ACT2R"] );
      _histos1d["hRef_act2LC_minus_act2RC_nonZero"]->Fill(_Charges["ACT2L"] - _Charges["ACT2R"] );
      _act2x = (-_ACTwidth/2.*_Charges["ACT2L"] + _ACTwidth/2.*_Charges["ACT2R"]) / (_Charges["ACT2L"] + _Charges["ACT2R"]);
      _histos1d["hRef_act2Xmap_nonZero"]->Fill(_act2x);
    }
    if (_Charges["ACT3L"] > minCharge &&  _Charges["ACT3R"] > minCharge) {
      _histos2d["hRef_act3LC_act3RC_nonZero"]->Fill(_Charges["ACT3L"], _Charges["ACT3R"] );
      _histos1d["hRef_act3LC_minus_act3RC_nonZero"]->Fill(_Charges["ACT3L"] - _Charges["ACT3R"] );
      _act3x = (-_ACTwidth/2.*_Charges["ACT3L"] + _ACTwidth/2.*_Charges["ACT3R"]) / (_Charges["ACT3L"] + _Charges["ACT3R"]);
      _histos1d["hRef_act3Xmap_nonZero"]->Fill(_act3x);
    }

    
    if (_Charges["ACT0L"] > minCharge &&  _Charges["ACT1L"] > minCharge) {
    _histos2d["hRef_act0LC_act1LC_nonZero"]->Fill(_Charges["ACT0L"], _Charges["ACT1L"] );
    _histos1d["hRef_act0LC_minus_act1LC_nonZero"]->Fill(_Charges["ACT0L"] - _Charges["ACT1L"] );
    }
    if (_Charges["ACT0R"] > minCharge &&  _Charges["ACT1R"] > minCharge) {
      _histos2d["hRef_act0RC_act1RC_nonZero"]->Fill(_Charges["ACT0R"], _Charges["ACT1R"] );
      _histos1d["hRef_act0RC_minus_act1RC_nonZero"]->Fill(_Charges["ACT0R"] - _Charges["ACT1R"] );
    }

    if (_debug)      cout << "charged b" << endl;
    
    // amplitudes vs tof
    _histos2d["hRef_TOFACT0A"]->Fill(_tof, _act0a);
    _histos2d["hRef_TOFACT1A"]->Fill(_tof, _act1a);
    _histos2d["hRef_TOFACT2A"]->Fill(_tof, _act2a);
    _histos2d["hRef_TOFACT3A"]->Fill(_tof, _act3a);

    // Trigger scintillators (unfortunatelly labelled as TOF all through out the code;-)
    // amplitude and charge vs tof:
    // TODO
    
    
    // lead glass vs acts and tof
    _histos2d["hRef_pbA_act23A"]->Fill(_pba, _act23aAver);
    _histos2d["hRef_pbC_act23A"]->Fill(_pbc, _act23aAver);
    _histos2d["hRef_TOFACT23A"]->Fill(_tof, _act23aAver);

    _histos2d["hRef_pbA_act0A"]->Fill(_pba, _act0a);
    _histos2d["hRef_pbA_act1A"]->Fill(_pba, _act1a);
    _histos2d["hRef_pbA_act1C"]->Fill(_pba, _act1c);
    
    _histos2d["hRef_PbATOF"]->Fill(_pba, _tof);
    _histos2d["hRef_TOFPbA"]->Fill(_tof, _pba);
    _histos2d["hRef_TOFPbC"]->Fill(_tof, _pbc);

    // ACT charges vs tof
    _histos2d["hRef_TOFACT0C"]->Fill(_tof, _act0c);
    _histos2d["hRef_TOFACT1C"]->Fill(_tof, _act1c);
    _histos2d["hRef_TOFACT2C"]->Fill(_tof, _act2c);
    _histos2d["hRef_TOFACT3C"]->Fill(_tof, _act3c);

    // lead glass vs acts and tof
    _histos2d["hRef_pbC_act23C"]->Fill(_pbc, _act23cAver);
    _histos2d["hRef_pbA_act23C"]->Fill(_pba, _act23cAver);
    _histos2d["hRef_TOFACT23C"]->Fill(_tof, _act23cAver);
    
    _histos2d["hRef_pbC_act0C"]->Fill(_pbc, _act0c);
    _histos2d["hRef_TOFACT0C"]->Fill(_tof, _act0c);
    _histos2d["hRef_pbC_act1C"]->Fill(_pbc, _act1c);
    _histos2d["hRef_TOFACT1C"]->Fill(_tof, _act1c);
    
    _histos2d["hRef_PbCTOF"]->Fill(_pbc, _tof);

    // act 2d plots
    _histos2d["hACT1CACT3C"]->Fill(_act1c, _act3c);
    _histos2d["hACT3CACT2C"]->Fill(_act3c, _act2c);
    _histos2d["hACT2CACT1C"]->Fill(_act2c, _act1c);

    // acraplet - weird electrons which do not see anything in the ACT
    if (_act23aAver != 1.5 && _tof >= 13.5 && _tof <= 16.5) {
      _histos2d["hHC0AHC1A"]->Fill(_hc0a, _hc1a);
      _histos2d["hHC0CHC1C"]->Fill(_hc0c, _hc1c);
    }


    // JK 5.3.2024
    // trigger scintillators (also used for tof:)
    // changes and amplitudes study

    this -> FillTrigScintHistos("");
    // and now also for p-like selection:
    // expected tof of protons:
    double ptofExp = _tofutil -> getTof("p", _momentum);
    double DtofExp = _tofutil -> getTof("D", _momentum);
    double TtofExp = _tofutil -> getTof("T", _momentum);
    double etofExp = _tofutil -> getTof("e", _momentum);
    double mutofExp = _tofutil -> getTof("mu", _momentum);
    double pitofExp = _tofutil -> getTof("pi", _momentum);
    double tsigma = 1.2; // ns
    // cout << " momentum: " << _momentum << " tof=" << _tof << " ptofExp=" << ptofExp <<  " " << fabs(_tof - ptofExp) << " " <<  tsigma << endl;
    if ( fabs(_tof - ptofExp) < tsigma) {
      this -> FillTrigScintHistos("_p-like");
    } else 
    if ( fabs(_tof - DtofExp) < tsigma) {
      this -> FillTrigScintHistos("_D-like");
    } else 
    if ( fabs(_tof - TtofExp) < tsigma) {
      this -> FillTrigScintHistos("_T-like");
    } else
      if ( fabs(_tof - etofExp) < tsigma) {
      this -> FillTrigScintHistos("_e-like");
    } else
      if ( fabs(_tof - mutofExp) < tsigma) {
      this -> FillTrigScintHistos("_mu-like");
    } else
      if ( fabs(_tof - pitofExp) < tsigma) {
      this -> FillTrigScintHistos("_pi-like");
    } 

    

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
    _histos2d["hnPeaksLeadGlassvsLeadGlassA"]->Fill(_pba, nPeaksLeadGlass);

    _histos2d["hnPeaksACT23vsLeadGlassA"]->Fill(_pba, nPeaksACT23Aver);
    _histos2d["hnPeaksToFvsLeadGlassA"]->Fill(_pba, nPeaksToFAver);

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
      // to move to a map, too?
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
      //if (_act1a < 0.15 && _pba < 0.25) {

      // cuts III
      //if (_act1a < 0.15 && _pba < 0.27 && _act23aAver < 0.55) {
      // cuts II
      //if (_act1a < 0.15 && _pba < 0.27) {
      // cuts I
      if (_act1a < 0.15) {
	// this cut is too destructive:  _act0a < 0.15
	//also too destructive: _pba < 0.15
	passed_act1a_cuts = true;
	_histos2d["hRef_pbA_act1A_act1cuts"]->Fill(_pba, _act1a);
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
  double pbanew = _Amplitudes["LeadGlass"];
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
    double Aj = _Amplitudes[chname];
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
    if (_Amplitudes[hdname] >= threshHodoscopeHit){
      _histos1d["hnHitsHodoscope"] -> Fill(i);
      //_histos1d["hnHitsHodoscope"] -> Fill(_channelToHodoscope[i]);
    }
  }
  


  
}

// ______________________________________________________________


void MakeAllDataPlots::Terminate()
{

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
