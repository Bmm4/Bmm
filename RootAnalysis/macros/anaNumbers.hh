#ifndef ANANUMBERS
#define ANANUMBERS

#include <string>
#include <utility>
#include <vector>
#include <map>

// the number of mass region windows: low sideband, high sideband, B0 window, Bs window
#define NWIN 4


class number {
public:
  number() : val(0), estat(0), esyst(0), etot(0), name("") {};
  ~number() {};
  void clear(std::string sname) {val = estat = esyst = etot = 0.; name = sname;}
  double val;
  double estat, esyst, etot;
  std::string name;
};

class anaNumbers {
public:
  anaNumbers(std::string name = "", int ichan = 0);
  ~anaNumbers() {};
  void clear();
  std::string fName;
  int fChan;
  // -- vectors for mass window yields
  std::vector<number> fMcYield, fObsYield, fFitYield;
  // -- overall results of fits
  number fSignalFit;
  // -- yield numbers for efficiency and acceptance
  number fGenAccFileYield, fGenAccYield, fGenFileYield, fGenYield, fCombGenYield, fProdGenYield;
  number fChanYield, fRecoYield, fCandYield, fMuidYield, fTrigYield;
  number fEffRecoPt;
  // -- efficiency and acceptance
  number fAcc, fEffCand, fEffAna, fEffTrigMC, fEffMuidMC, fEffTot, fEffProdMC;
  number fEffFilter; // this is the MC filter efficiency, used for eq lumi calculation!
  // -- cross feed: pxx signal in signal, pyx other in signal, etc.
  number fPxx, fPyx, fPlx, fPhx;
};

#endif
