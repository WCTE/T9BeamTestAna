#include "tofUtil.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

tofUtil::tofUtil() {

    conv = 1.e9; // to ns
    //TB2022 l = 2.90 // m
    //TB2023:
    Ltof = 3.49; // m
    c = 299792458; // m/s

    this -> MakeDicts();
}


tofUtil::~tofUtil() {}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void tofUtil::MakeDicts()
{

  ms["e"] = 0.511;
  ms["mu"] = 105.66;
  ms["pi"] = 139.57;
  ms["K"] = 493.7;
  ms["p"] = 938.3;
  ms["D"] = 1876.;
  ms["T"] = 3.01604928*931.494; // isotope mass in Da

  pcols["e"] = kRed;
  pcols["mu"] = kBlue;
  pcols["pi"] = kGreen+2;
  pcols["K"] = kGray+1;
  pcols["p"] = kBlack;
  pcols["D"] = kViolet;
  pcols["T"] = kCyan+2;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tofUtil::getTof(TString part, int momentum)
{
  //cout << "mass: " << ms[part] << endl;
  return this -> getTof(ms[part], momentum);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tofUtil::getTof(double m, int momentum)
{
  return Ltof / c*sqrt(1.+pow(m/momentum,2))*conv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
double tofUtil::TofToMomentum(double tof, double m)
{
    //the tof needs to be the absolute flying time
  double p = 0.;
  double val = pow((tof) * c / (conv * Ltof), 2) - 1;
  if (val > 0.) 
    p = m/sqrt(val);
  return p;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// input: tofdiff, i.e. time of flight difference of a particle and tof of electrons!
// then also mass of the particle and uncertainties in the emasured times of the particle and of electrons (not the tof resolution or fitted leaks widths!)
std::pair<double,double> tofUtil::TofDiffToMomentum(double tofdiff, double m, double sigmate, double sigmatParticle)
{
  //the tof needs to be the tof subtracted by the electrons TOF
  double coL = c / (conv * Ltof); // c over L ;) [s^{-1}]
  double val = pow((tofdiff) * coL + 1, 2) - 1;
  double p = 0.;
  double perr = 0.;
  if (val > 0.) {
    p = m/sqrt(val);
    if (sigmate > 0. || sigmatParticle > 0.) {
      double sigmatSq = pow(sigmate,2) + pow(sigmatParticle,2);
      if (sigmatSq > 0.) {
	double sigmat = sqrt(sigmatSq);
	perr = pow(p,3) / pow(m,2) * ( tofdiff*coL + 1) * coL * sigmat;
      }
    }
  }
  return std::pair<double,double>(std::make_pair(p, perr));
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double tofUtil::getTofDiff(TString particle1, TString particle2, int momentum)
{
  double m1 = ms[particle1];
  double m2 = ms[particle2];
  double t1 = getTof(m1, momentum);
  double t2 = getTof(m2, momentum);
  return t2 - t1;
}
		  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
