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
  fitPsYield(TH2D *h2, int verbose);
  ~fitPsYield();
  void initFromHist(TH2D *h);

  double getSignalRMS() {return fSummary.fSgSigma; }
  double getSignalPeak() {return fSummary.fSgPeak; }

  double getSignalYield() {return fSummary.fSg; }
  double getSignalError() {return fSummary.fSgE;}

  double getBackgroundYield() {return fSummary.fBg; }
  double getBackgroundError() {return fSummary.fBgE;}

  double getSignalW8Yield() {return fW8Combined->fResults.fSg; }
  double getSignalW8Error() {return fW8Combined->fResults.fSgE;}

  double getSignalUnW8Yield() {return fUnW8Combined->fResults.fSg; }
  double getSignalUnW8Error() {return fUnW8Combined->fResults.fSgE;}

  double getSignalUnW8Chi2() {return fUnW8Combined->fChi2; }
  double getSignalUnW8Ndof() {return fUnW8Combined->fNdof; }


  void printSummary();

  TH1D* getUnweightedCombination() {return fCombined;}
  TH1D* getWeightedCombination()   {return fCombinedW8;}


  void fitBu2JpsiKp(int limitpars, std::string pdfprefix, int whichfit = 0, double lo = 5.0, double hi = 5.8, double sigma = 0.020);
  void fit0_Bu2JpsiKp(psd *res, int limitpars = 0, std::string pdfprefix = ".", double lo = 5.0, double hi = 5.8, double sigma = 0.020, bool keepFunctions = false);
  void fit0_Bu2JpsiKp(TH1D *h, int limitpars = 0, std::string pdfprefix = ".", double lo = 5.0, double hi = 5.8, double sigma = 0.020);
  void fit1_Bu2JpsiKp(psd *res, int limitpars = 0, std::string pdfprefix = ".", double lo = 5.0, double hi = 5.8, double sigma = 0.020, bool keepFunctions = false);
  void fit1_Bu2JpsiKp(TH1D *h, int limitpars = 0, std::string pdfprefix = ".", double lo = 5.0, double hi = 5.8, double sigma = 0.020);

  void fitBs2JpsiPhi(int limitpars, std::string pdfprefix, int whichfit = 0, double lo = 5.0, double hi = 5.8, double sigma = 0.020);
  void fit0_Bs2JpsiPhi(psd *res, int limitpars = 0, std::string pdfprefix = ".", double lo = 5.0, double hi = 5.8, double sigma = 0.020, bool keepFunctions = false);
  void fit0_Bs2JpsiPhi(TH1D *h, int limitpars = 0, std::string pdfprefix = ".", double lo = 5.0, double hi = 5.8, double sigma = 0.020);
  void fit1_Bs2JpsiPhi(psd *res, int limitpars = 0, std::string pdfprefix = ".", double lo = 5.0, double hi = 5.8, double sigma = 0.020, bool keepFunctions = false);
  void fit1_Bs2JpsiPhi(TH1D *h, int limitpars = 0, std::string pdfprefix = ".", double lo = 5.0, double hi = 5.8, double sigma = 0.020);

  void fitBd2JpsiKstar(int limitpars, std::string pdfprefix, int whichfit = 0);
  void fit0_Bd2JpsiKstar(psd *res, int limitpars = 0, std::string pdfprefix = ".", bool keepFunctions = false);
  void fit0_Bd2JpsiKstar(TH1D *h, int limitpars = 0, std::string pdfprefix = ".");


  TF1* getFunction(std::string name);
  TF1* listFunctions();

private:
  int         fVerbose;
  std::string fBaseName;
  initFunc    *fpIF;

  TH2D        *fH2;
  TH1D        *fCombined, *fCombinedW8;
  std::vector<psd*> fData;                 // vectors for different prescales with data
  psd         *fW8Combined, *fUnW8Combined;
  Results     fSummary;

  std::vector<double> fPar, fParE;
  double         fCombMax, fCombS2All, fCombS2AllE;
  std::map<std::string, TF1*> fFunctions;
};


#endif
