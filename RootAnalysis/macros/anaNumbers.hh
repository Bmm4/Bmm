#ifndef ANANUMBERS
#define ANANUMBERS

#include <string>
#include <utility>
#include <vector>
#include <map>

// the number of mass region windows: low sideband, high sideband, B0 window, Bs window
#define NWIN 5


class number {
public:
  number() : val(0), estat(0), esyst(0), etot(0), name("") {};
  ~number() {};
  void clear(std::string sname) {val = estat = esyst = etot = 0.; name = sname;}
  void setErrors(double sta, double sys);
  void setErrors(double sta, double sys, double tot);
  void add2Errors(number &);
  void calcEtot();
  double val;
  double estat, esyst, etot;
  std::string name;
};

class anaNumbers {
public:
  anaNumbers(std::string name = "", int ichan = 0);
  ~anaNumbers() {};
  void clear();
  std::string fName, fNameMc, fNameDa, fSel;
  int fChan;
  // -- NWIN-element vectors for mass window yields: lo, Bd, Bs, hi, all
  std::vector<number> fMcYield, fObsYield, fFitYield;
  // -- NWIN-element vector for signal migration
  std::vector<number> fFrac;
  // -- overall results of fits and scaled yields
  number fSignalFit, fScaledYield;
  number fW8SignalFit; // with (correction) weights applied
  double fScaleFactor;
  // -- yield numbers for dataset
  number fGenFileYield, // number of (B) events in file/dataset
    fGenYield,          // number of (B) events generated (correcting for gen-level filter, but NOT for effFilter!)
    fTotGenYield, fProdGenYield;
  // -- yield numbers for corresponding ACC dataset
  number fAccGenFileYield, // number of (B) events in ACC file/dataset
    fAccGenYield;          // number of (B) events in ACC generated (correcting for gen-level filter, but NOT for effFilter!)
  number fAccRecoYield, fAccCandYield;
  number fCandYield, fAnaYield, fMuidYield, fTrigYield;
  // -- efficiency and acceptance
  number fAcc, fEffCand, fEffAna, fEffAccAna, fEffTrigMC, fEffMuidMC, fEffTot, fEffProdMC;
  number fEffFilter; // this is the MC filter efficiency, used for eq lumi calculation!
  number fEffGenSel; // this is the gen-level filter efficiency, needed to go from non-Acc samples to acc samples! Obtained from ratio of effFilter's
};

#endif
