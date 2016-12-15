#ifndef FITPSYIELD_HH
#define FITPSYIELD_HH

#include <utility>
#include <vector>
#include <string>
#include <map>

#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"

#include "initFunc.hh"

// ----------------------------------------------------------------------
// -- fitPsYield for a SINGLE channel!
// -------------
//
// ----------------------------------------------------------------------

class Results {
public:
  Results() {};
  ~Results() {};
  void clear() {
    fSgSigma = fSgPeak = 0.;
    fSg  = fSgE = 0.;
    fBg  = fBgE = 0.;
    fBg0 = fBg0E = 0.;
  }
  double fSgSigma, fSgPeak; // whatever sigma and peak mean
  double fSg,      fSgE;    // integrated over signal region (whatever that means)
  double fBg,      fBgE;    // integrated over complete histogram
  double fBg0,     fBg0E;   // in signal region
};

class psd {
 public:
  psd(): fEntries(0), fPs(0), fChi2(0.), fProb(0.), fNdof(0), fH1(0) {fResults.clear();};
  ~psd() {delete fH1;}
  int      fEntries;
  int      fPs;
  double   fChi2, fProb;
  int      fNdof;
  TH1D*    fH1;
  Results  fResults;
};

class fitPsYield {
public:
  fitPsYield(std::string hname = "hNo_cnc_bupsikData_chan0", TDirectory *pD = 0, int verbose = 1);
  ~fitPsYield();

  double getSignalYield() {return fSummary.fSg; }
  double getSignalError() {return fSummary.fSgE;}

  void printSummary();

  TH1D* getUnweightedCombination() {return fCombined;}
  TH1D* getWeightedCombination()   {return fCombinedW8;}

  void fitBu2JpsiKp(int limitpars, std::string pdfprefix);

  void  fit0_Bu2JpsiKp(psd *res, int limitpars = 0, std::string pdfprefix = ".");
  void  fit0_Bs2JpsiPhi(int limitpars = 0);

private:
  int         fVerbose;
  std::string fBaseName;
  TDirectory  *fpDir;
  initFunc    *fpIF;

  TH2D        *fH2;
  TH1D        *fCombined, *fCombinedW8;
  std::vector<psd*> fData;                 // vectors for different prescales with data
  psd         *fW8Combined, *fUnW8Combined;
  Results     fSummary;

  std::vector<double> fPar, fParE;
  double         fCombMax, fCombS2All, fCombS2AllE;
};


#endif
