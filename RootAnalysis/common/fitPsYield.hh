#ifndef FITPSYIELD_HH
#define FITPSYIELD_HH

#include <utility>
#include <vector>
#include <string>
#include <map>

#include "TH1.h"
#include "TH2.h"

// ----------------------------------------------------------------------
// -- fitPsYield
// -------
//
// class to fit n TH2D filled with w8=1 for different prescales
//
// ----------------------------------------------------------------------

class fitPsYield {
public:
  fitPsYield(std::string hname, int verbose = 1);
  ~fitPsYield();

private:
  int         fVerbose;
  std::string fBaseName;

  std::vector<TH2D*> fHists;
};


#endif
