#ifndef PLOTWORK_h
#define PLOTWORK_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotWork: public plotClass {

public :
  plotWork(std::string dir = "results",
	   std::string files = "plotWork.2016.files",
	   std::string cuts = "plotClass.2016.cuts",
	   std::string setup = "");
  virtual        ~plotWork();

  // -- Main analysis methods
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void makeAll(std::string what = "all");

  // -- PV studies
  void pvStudy(std::string filename, std::string selection = "");

  // -- yield stability
  void yieldStability(std::string dsname, std::string trg = "HLT");

  // -- genSummary plots
  void genSummary(std::string dsname, std::string dir);

  // -- fitting studies
  void fitStudies(std::string ds1, std::string tag, int nevt = -1, int nstart = -1);
  void  fitStudiesFit0(TH1D *h1, int chan);
  double MKKLO, MKKHI, DR, PTK1, PTK2, PTPSI;

  // -- misreconstructed background
  void wrongReco(std::string ds1, std::string mode, std::string selection = "hlt");
  void plotWrongReco(std::string var, int nbin, double min, double max, std::string selection,
		     std::string wds, std::string wdir,
		     std::string cds, std::string cdir);

  // -- trigger efficiency studies
  std::string selectionString(int imode, int itrig);
  void plotTisEfficiency(std::string dsname);
  void runTisEfficiency(std::string dsname);
  void plotL1Seeds(std::string dsname);

  void refTrgEfficiency(std::string selection, std::string dsname = "bupsikMc");
  void efficiencyVariable(std::string var, std::string effvar = "hlt", int iselection = 10,
			  int nbin = 20, double xmin = 0., double xmax = 20., std::string dsname = "bupsikMc");




  std::string removeVarFromSelection(std::string var, std::string selection);
  void bookHist(std::string dsname);

  /// NOTE: This works with the output of genAnalysis!!!!
  void prodSummary(string ds1, int year);


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

  std::map<int, TH2D*> fYieldRTR, fYieldHLT;

  int fRefTrigger;

  ClassDef(plotWork,1)

};


#endif
