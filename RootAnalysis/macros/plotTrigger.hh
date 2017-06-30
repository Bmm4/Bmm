#ifndef PLOTTRIGGER_h
#define PLOTTRIGGER_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotTrigger: public plotClass {

public :
  plotTrigger(std::string dir = "results",
	   std::string files = "plotResults.2017.files",
	   std::string cuts = "baseCuts.cuts",
	   std::string setup = "");
  virtual        ~plotTrigger();

  // -- Main analysis methods
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void makeAll(std::string what = "all");

  // -- trigger efficiency studies
  void plotTisEfficiency(std::string dsname);
  void runTisEfficiency(std::string dsname);
  void plotL1Seeds(std::string dsname);

  void refTrgEfficiency(std::string selection, std::string dsname = "bupsikMc");

  // -- code for loops
  void   loopFunction1();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

private:


  // ----------------------------------------------------------------------
  TH1D *fpHnorm, *fpHpass;
  TH1D *fpHL1s0, *fpHL1s1, *fpHL1s2, *fpHL1s3, *fpHL1s4, *fpHL1s5, *fpHL1s6, *fpHL1All;
  std::vector<TH1D *> fHma, fHmc;
  double fS, fSE, fN, fNE, fW, fWE, fB, fBE, fChi2Dof;
  double fSsigma, fSRMS;
  double fEntries;

  std::vector<TH1D*> fHmass0, fHmass1;
  std::vector<TH1D*> fHBd0, fHBd1;
  std::vector<TH1D*> fHBs0, fHBs1;
  std::vector<TH1D*> fHBg0, fHBg1;
  TProfile *fpPmass0, *fpPmass1;

  int fRefTrigger;

  ClassDef(plotTrigger,1)

};


#endif
