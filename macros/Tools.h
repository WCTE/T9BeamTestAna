#ifndef tools_h
# define tools_h



#include <iostream>


double GetBeta(double mass, double momentum);

TObjArray * tokenizeTString(const TString& inputString, const TString& delimiter);

#endif
