#include "Tools.h"


double GetBeta(double mass, double momentum) {
  double bg = momentum/mass;
  double beta = sqrt(bg*bg/(1+bg*bg));
  return beta;
}


TObjArray * tokenizeTString(const TString& inputString, const TString& delimiter) {
    TObjArray *tokens = inputString.Tokenize(delimiter);
    //cout << "gonna print tokens..." << endl;
    /*
      int it = 0;
    if (tokens) {
      //cout << "printing tokens:" << endl;
      TIter next(tokens);
      while (TObject *obj = next()) {
	TObjString *token = dynamic_cast<TObjString*>(obj);
	cout << "token: " << token << endl;
	if (token) {
	  //cout << "token" << it << " \"" << token->GetString().Data() << "\"" << endl;
	  it++;
	}
      }
    } else {
      cerr << "Error tokenizing the output file name!" << endl;
    }
    */
    return tokens;
}
