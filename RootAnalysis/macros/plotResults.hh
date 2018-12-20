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
	      std::string files = "plotResults.2016BF.files",
	      std::string cuts = "baseCuts.2016.cuts",
	      std::string setup = "",
	      int year = 0);
  virtual ~plotResults();

  // -- Main analysis methods
  virtual void makeAll(std::string what = "all");
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void bookHist(std::string dsname);

  void dumpDatasets();
  void genSummary(std::string dsname, std::string dir);

  // -- main results
  void fillAndSaveHistograms(int start = 0, int nevents = -1);
  void saveHistograms(std::string smode, TFile *f = 0);
  void resetHistograms(bool deleteThem = false);
  void initNumbers(anaNumbers &a);
  void otherNumbers(string ds);
  void rareBgHists(string ds, int nevents = -1, int nbg = -1);
  void scaleYield(anaNumbers &aSig, anaNumbers &aNorm, double pRatio = 1.0, bool useW8 = true);
  void getAccAndEffFromEffTree(string ds, anaNumbers &a, cuts &b, int proc);
  void numbersFromHist(anaNumbers &a, std::string syst);
  std::string rareAccName(std::string sname);
  bool skipThisBg(std::string name);
  void calculateNumbers(std::string mode, int chan);
  void calculateB2JpsiNumbers(anaNumbers &a);
  void calculateSgNumbers(anaNumbers &a);
  void calculateCombBgNumbers(anaNumbers &a, int mode = 0, double lo = 5.45, double hi = 5.9);
  void calculateRareBgNumbers(int chan);
  void calculateBs2Bu(int chan);
  void calculatePerformance(int chan);
  void triggerEfficiency(int chan, double &eff, double &effE, std::string mode, std::string sample = "bg");
  void scanBDT(std::string fname, bool createTexFile = true);
  void scanBDTEffRatio(std::string fname, std::string string1, std::string string2);
  void sysAna(std::string sample1, std::string sample2, std::string massc = "C", int ichan = 0, int buversion = 0, int bsversion = 1);
  void sysNorm(std::string sample1, int ichan = 0, int version = 0);
  void sysBsVsBu();
  void sysEffLifetime(std::string cutname = "WithAnaCuts");
  void sysDoubleRatio(std::string sample1, std::string sample2, std::string chansel, std::string var, std::string cutlevel, double cut);
  void displayScanBDT(std::string what = "CSBF", int mode = 0, int chan = 0);
  void showScanBDT(std::string what = "all");
  double findVarValue(std::string varName, std::vector<std::string> &lines);
  double findVarEtot(std::string varName, std::vector<std::string> &lines);
  double findVarEstat(std::string varName, std::vector<std::string> &lines);

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
  std::map<std::string, std::vector<TH1D*> > fhBdtCrossCheck;
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

  // -- anti-muon sample
  std::map<std::string, std::vector<TH1D*> > fhAmMassWithAllCuts;
  std::map<std::string, std::vector<TH1D*> > fhAmMassWithMuonCuts;
  std::map<std::string, std::vector<TH1D*> > fhAmW8MassWithAllCuts;
  // -- real (anti-)muon sample
  std::map<std::string, std::vector<TH1D*> > fhRMMassWithMuonCuts,  fhRMMassWithAllCuts;
  std::map<std::string, std::vector<TH1D*> > fhRAmMassWithMuonCuts, fhRAmMassWithAllCuts;

  std::map<std::string, std::vector<TH1D*> > fhMassWithMassCuts;
  std::map<std::string, std::vector<TH2D*> > fhNorm, fhNormC;
  std::map<std::string, std::vector<TH2D*> > fhW8Norm, fhW8NormC;
  std::map<std::string, std::vector<TH2D*> > fhMassSysStep0;
  std::map<std::string, std::vector<TH2D*> > fhMassSysStep0C;
  std::map<std::string, std::vector<TH2D*> > fhMassSysStep1;
  std::map<std::string, std::vector<TH2D*> > fhMassSysStep1C;

  std::vector<TH1D*> fhBdt;

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
