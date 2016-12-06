#ifndef PLOTRESULTS_h
#define PLOTRESULTS_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"

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
  void fillAndSaveHistograms(int nevents = -1);
  void saveHistograms(std::string smode);
  void resetHistograms();


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

  std::vector<TH2D*> fhAccAll, fhAccPass;
  std::vector<TH1D*> fhAccPtAll, fhAccPtPass, fhAccEtaAll, fhAccEtaPass;

  std::vector<TH1D*> fhGenAndAccNumbers;
  std::vector<TH1D*> fhMassAbsNoCuts;
  std::vector<TH1D*> fhMassNoCuts;
  std::vector<TH1D*> fhMassWithAnaCuts;
  std::vector<TH1D*> fhMassWithMuonCuts;
  std::vector<TH1D*> fhMassWithTriggerCuts;
  std::vector<TH1D*> fhMassWithAllCuts;
  std::vector<TH1D*> fhMassWithAllCutsBlind;

  std::vector<TH1D*> fhW8MassWithAllCuts;
  std::vector<TH1D*> fhW8MassWithAllCutsBlind;

  std::vector<TH1D*> fhMassWithMassCuts;
  std::vector<TH1D*> fhNo, fhNoC;
  std::vector<TH1D*> fhCs, fhCsC;
  std::vector<TH1D*> fhB0, fhB0C;

  // ----------------------------------------------------------------------
  ClassDef(plotResults,1)

};


#endif
