// JK 6.3.2024 based on tofUtil.py

#include <iostream>
#include <map>
#include <cmath>

using namespace std;

class tofUtil
{


 private:
  double Ltof;
  double c;
  double conv;
  
  map<TString,double> ms;
  map<TString,double> pcols;
  
 public:

  tofUtil();
  ~tofUtil();
  
  void MakeDicts();
  double getTof(double m, int momentum);
  double getTof(TString part, int momentum);
  double TofToMomentum(double tof, double m);
  std::pair<double,double> TofDiffToMomentum(double tofdiff, double m, double sigmate = 0., double sigmatParticle = 0.);
  double getTofDiff(TString particle1, TString particle2, int momentum);

};




