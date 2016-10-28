#ifndef PLOTSTUFF_h
#define PLOTSTUFF_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotStuff: public plotClass {

public :
  plotStuff(std::string dir = "results",
	   std::string files = "plotResults.2016.files",
	   std::string cuts = "baseCuts.cuts",
	   std::string setup = "");
  virtual        ~plotStuff();

  // -- Main analysis methods
  virtual void bookHist(std::string dsname);
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void makeAll(std::string what = "all");

  // -- PV studies
  void pvStudy(std::string dsname, std::string selection = "&& fl1 > 0.01", std::string fmod = "");

  // -- yield stability
  void yieldStability(std::string dsname, std::string trg = "HLT");


  // -- code for loops
  void   loopFunction1();
  void   loopFunction2();
  void   loopFunction3();
  void   loopFunction4();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);



private:


  // ----------------------------------------------------------------------
  TH1D *fpHnorm, *fpHpass;
  TH1D *fpHL1s0, *fpHL1s1, *fpHL1s2, *fpHL1s3, *fpHL1s4, *fpHL1s5, *fpHL1s6, *fpHL1All;
  std::vector<TH1D *> fHma, fHmc;
  double fS, fSE, fN, fNE, fW, fWE, fB, fBE, fChi2Dof;
  double fSsigma, fSRMS;
  double fEntries;

  std::map<std::string, TH2D*> fYieldRTR, fYieldHLT;

  int fRefTrigger;

  ClassDef(plotStuff,1)

};


#endif
