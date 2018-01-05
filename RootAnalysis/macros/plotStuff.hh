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
	    std::string setup = "",
	    int year = 0);
  virtual        ~plotStuff();

  // -- Main analysis methods
  virtual void bookHist(std::string dsname);
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void makeAll(std::string what = "all");

  // -- PV studies: our PV choice or best pointing angle
  void pvStudy(std::string dsname, std::string selection = "&& fl1 > 0.01", std::string fmod = "");
  // -- Pileup studies: what happens when another PV is close
  void puStudy(std::string dsname, std::string dsname2 = "");
  void o2Profile(TProfile *p, std::string dsname, std::string dsname2, int i);

  // -- yield stability
  void yieldStability(std::string dsname, std::string trg = "HLT");
  void yieldStabilityRatios(std::string trg = "HLT");
  void yieldStudy(int run = 278273, string ds = "bmmData");
  void runStudy(std::string dsname, std::string mode = "ana");

  // -- misreconstructed background
  void wrongReco(std::string ds1, std::string mode, std::string selection = "hlt");
  void plotWrongReco(std::string var, int nbin, double min, double max, std::string selection,
		     std::string wds, std::string wdir,
		     std::string cds, std::string cdir);

  // -- variety of methods
  void massResolution(std::string file1, std::string file2);
  void tauEfficiency(string varname, string cut, string otherSelection, string dsname);


  // -- code for loops
  void setupPvTree(TTree *t);
  void loopFunction1();
  void loopFunction2();
  void loopFunction3();
  void loopFunction4();
  void loopFunction5();

  void loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);


  void plot2016Eras(double ymax);

private:


  // ----------------------------------------------------------------------
  TH1D *fpHnorm, *fpHpass;
  TH1D *fpHL1s0, *fpHL1s1, *fpHL1s2, *fpHL1s3, *fpHL1s4, *fpHL1s5, *fpHL1s6, *fpHL1All;
  std::vector<TH1D *> fHma, fHmc;
  double fS, fSE, fN, fNE, fW, fWE, fB, fBE, fChi2Dof;
  double fSsigma, fSRMS;
  double fEntries;

  std::map<std::string, TH2D*> fYieldHists;
  std::map<std::string, TH1D*> fPlots;
  int fSplitRun;

  int fRefTrigger;

  struct pvTreeData fpv;
  static const int NCHAN = 5;
  TH1D *fpHmultFar[NCHAN], *fpHmultClose05[NCHAN]
    , *fpHdzmin[NCHAN]    , *fpHpvn[NCHAN]
    , *fpHmultFar2[NCHAN], *fpHmultClose2[NCHAN]
    , *fpHlz1[NCHAN],      *fpHlz2[NCHAN]
    , *fpHlzs1[NCHAN],     *fpHlzs2[NCHAN]
    ;

  // -- vs dzmin
  TProfile *fpP0Mult[NCHAN], *fpP0fls3d[NCHAN], *fpP0fl3d[NCHAN], *fpP0dfl3d[NCHAN], *fpP0tau[NCHAN], *fpP0dtau[NCHAN];
  TProfile *fpP0l1[NCHAN], *fpP0l2[NCHAN], *fpP0npv[NCHAN], *fpP0Iso[NCHAN], *fpP0Ntrk[NCHAN], *fpP0PvIps[NCHAN], *fpP0Bdt[NCHAN];

  // -- vs lz(1)
  TProfile *fpP1Mult[NCHAN], *fpP1fls3d[NCHAN], *fpP1fl3d[NCHAN], *fpP1dfl3d[NCHAN], *fpP1tau[NCHAN], *fpP1dtau[NCHAN];
  TProfile *fpP1l2[NCHAN], *fpP1dzmin[NCHAN], *fpP1npv[NCHAN], *fpP1Iso[NCHAN], *fpP1Ntrk[NCHAN], *fpP1PvIps[NCHAN], *fpP1Bdt[NCHAN];

  // -- vs lz(2)
  TProfile *fpP2Mult[NCHAN], *fpP2fls3d[NCHAN], *fpP2fl3d[NCHAN], *fpP2dfl3d[NCHAN], *fpP2tau[NCHAN], *fpP2dtau[NCHAN];
  TProfile *fpP2l1[NCHAN], *fpP2dzmin[NCHAN], *fpP2npv[NCHAN], *fpP2Iso[NCHAN], *fpP2Ntrk[NCHAN], *fpP2PvIps[NCHAN], *fpP2Bdt[NCHAN];

  // -- vs NPV
  TProfile *fpP3Mult[NCHAN], *fpP3fls3d[NCHAN], *fpP3fl3d[NCHAN], *fpP3dfl3d[NCHAN], *fpP3tau[NCHAN], *fpP3dtau[NCHAN];
  TProfile *fpP3l1[NCHAN], *fpP3l2[NCHAN], *fpP3dzmin[NCHAN], *fpP3Iso[NCHAN], *fpP3Ntrk[NCHAN], *fpP3PvIps[NCHAN], *fpP3Bdt[NCHAN];

  // -- mass resolution
  std::vector<TH1D*> fHmass0, fHmass1;

  std::vector<int> fLargeRuns;

  std::map<std::string, TH2D*> fHLs0, fHLs1;
  std::map<std::string, TProfile*> fProf;
  std::map<std::string, TH1D *> fvHists;

  ClassDef(plotStuff,1)

};


#endif
