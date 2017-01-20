#ifndef PLOTSTUFF_h
#define PLOTSTUFF_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"
#include "pvTreeData.hh"

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

  // -- PV studies: our PV choice or best pointing angle
  void pvStudy(std::string dsname, std::string selection = "&& fl1 > 0.01", std::string fmod = "");
  // -- Pileup studies: what happens when another PV is close
  void puStudy(std::string dsname);

  // -- yield stability
  void yieldStability(std::string dsname, std::string trg = "HLT");
  void yieldStabilityOld(std::string dsname, std::string trg = "HLT");
  void yieldStabilityRatios(std::string trg = "HLT");

  // -- code for loops
  void setupPvTree(TTree *t);
  void loopFunction1();
  void loopFunction2();
  void loopFunction3();
  void loopFunction4();

  void loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);



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

  struct pvTreeData fpv;
  static const int NCHAN = 5;
  TH1D *fpHmultFar[NCHAN], *fpHmultClose05[NCHAN]
    , *fpHdzmin[NCHAN]
    , *fpHmultFar2[NCHAN], *fpHmultClose2[NCHAN]
    , *fpHlz1[NCHAN], *fpHtmlz1[NCHAN]
    , *fpHlz2[NCHAN]
    ;

  TProfile *fpP1Mult[NCHAN], *fpP1flsxy[NCHAN], *fpP1fls3d[NCHAN], *fpP1fl3d[NCHAN], *fpP1dfl3d[NCHAN], *fpP1tau[NCHAN], *fpP1dtau[NCHAN];
  TProfile *fpP2Mult[NCHAN], *fpP2flsxy[NCHAN], *fpP2fls3d[NCHAN], *fpP2fl3d[NCHAN], *fpP2dfl3d[NCHAN], *fpP2tau[NCHAN], *fpP2dtau[NCHAN];
  TProfile *fpP3tau[NCHAN],  *fpP3dtau[NCHAN];

  ClassDef(plotStuff,1)

};


#endif
