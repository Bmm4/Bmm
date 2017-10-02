#ifndef PLOTRESULTS_h
#define PLOTRESULTS_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"
#include "anaNumbers.hh"

#include <fstream>
#include <iostream>


// ----------------------------------------------------------------------
class plotResults: public plotClass {

public :
  plotResults(std::string dir = "results",
	      std::string files = "plotResults.2016.files",
	      std::string cuts = "baseCuts.cuts",
	      std::string setup = "");
  virtual ~plotResults();

  // -- Main analysis methods
  virtual void makeAll(std::string what = "all");
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void bookHist(std::string dsname);

  void dumpDatasets();
  void genSummary(std::string dsname, std::string dir);

  // -- main results
  void fillAndSaveHistograms(int nevents = -1, int start = 0);
  void saveHistograms(std::string smode);
  void resetHistograms(bool deleteThem = false);
  void initNumbers(anaNumbers &a);
  void otherNumbers(string ds);
  void rareBgHists(string ds, int nevents = -1);
  void scaleYield(anaNumbers &aSig, anaNumbers &aNorm, double pRatio = 1.0, bool useW8 = true);
  void getAccAndEffFromEffTree(string ds, anaNumbers &a, cuts &b, int proc);
  void numbersFromHist(anaNumbers &a, std::string syst);
  std::string rareAccName(std::string sname);
  bool skipThisBg(std::string name);
  void calculateNumbers(std::string mode);
  void calculateB2JpsiNumbers(anaNumbers &a);
  void calculateSgNumbers(anaNumbers &a);
  void calculateCombBgNumbers(anaNumbers &a, int mode = 0, double lo = 5.45, double hi = 5.9);
  void calculateRareBgNumbers(int chan);
  void calculateBs2Bu(int chan);
  void calculatePerformance(int chan);

  enum INTMODE {LO=0, BD, BS, HI, ALL};
  double massIntegral(TH1* h, INTMODE i, int chan);
  void printNumbers(anaNumbers &a, ostream &OUT);
  void dumpTex(number &a, const std::string &s, int ndigits = 3, int sgn = 0);


  // -- code for loops
  void   loopFunction1();
  void   loopFunction2();
  void   loopFunction3();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

private:

  static const int NBINS = 100;
  double fMassLo, fMassHi
    , fBgLo, fBgHi
    , fSgLo, fSgHi
    , fNoLo, fNoHi
    , fCsLo, fCsHi
    ;

  std::string fHistWithAllCuts;

  std::vector<std::string> fHistStrings;
  std::map<std::string, std::vector<TH2D*> > fhAccAll, fhAccPass;
  std::map<std::string, std::vector<TH1D*> > fhAccPtAll, fhAccPtPass, fhAccEtaAll, fhAccEtaPass;

  std::map<std::string, std::vector<TH1D*> > fhGenAndAccNumbers;
  std::map<std::string, std::vector<TH1D*> > fhMassAbsNoCuts;
  std::map<std::string, std::vector<TH1D*> > fhMassNoCuts;
  std::map<std::string, std::vector<TH1D*> > fhMassWithAnaCuts;
  std::map<std::string, std::vector<TH1D*> > fhMassWithMuonCuts;
  std::map<std::string, std::vector<TH1D*> > fhMassWithTriggerCuts;
  std::map<std::string, std::vector<TH1D*> > fhMassWithAllCuts;
  std::map<std::string, std::vector<TH1D*> > fhMassWithAllCutsBlind;
  std::map<std::string, std::vector<TH1D*> > fhMassWithAllCutsSeagull;
  std::map<std::string, std::vector<TH1D*> > fhMassWithAllCutsSeagullBlind;
  std::map<std::string, std::vector<TH1D*> > fhMassWithAllCutsCowboy;
  std::map<std::string, std::vector<TH1D*> > fhMassWithAllCutsCowboyBlind;

  std::map<std::string, std::vector<TH1D*> > fhW8MassWithAllCuts;
  std::map<std::string, std::vector<TH1D*> > fhW8MassWithAllCutsSeagull;
  std::map<std::string, std::vector<TH1D*> > fhW8MassWithAllCutsCowboy;

  std::map<std::string, std::vector<TH1D*> > fhMassWithMassCuts;
  std::map<std::string, std::vector<TH2D*> > fhNorm, fhNormC;
  std::map<std::string, std::vector<TH2D*> > fhW8Norm, fhW8NormC;

  bool fSaveSmallTree;

  // -- vector for each channel
  std::vector<anaNumbers>
  fNoNumbers,        // B+ normalization
    fCsNumbers,      // Bs normalization
    fB0Numbers,      // B0 -> J/psi K*0
    fBsmmNumbers,    // Bs -> mu mu signal
    fBdmmNumbers,    // B0 -> mu mu signal
    fHhNumbers,      // peaking background
    fSlNumbers,      // semileptonic background
    fCombNumbers,    // combinatorial background (from data)
    fNpNumbers,      // non-peaking background: SL (from MC) and combinatorial background (from data)
    fBgNumbers,      // background: SL and HH (from MC) and  combinatorial background (from data)
    fSgAndBgNumbers; // background and signal
  std::map<std::string, std::vector<anaNumbers*> > fRareNumbers;
  // -- map (name) -> (vector with per-channel errors) for all relative systematic uncertainties
  std::map<std::string, std::vector<double> > fSystematics;
  // ----------------------------------------------------------------------
  ClassDef(plotResults,1)

};


#endif
