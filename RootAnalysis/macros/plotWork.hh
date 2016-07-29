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
  virtual void makeAll(int bitmask = 0);

  // -- yield stability
  void yieldStability(std::string dsname, std::string trg);

  // -- fitting studies
  void fitJpsiPhi(string ds1);


  // -- trigger efficiency studies
  std::string selectionString(int imode, int itrig);
  void plotTisEfficiency(std::string dsname);
  void runTisEfficiency(std::string dsname);

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

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

private:


  // ----------------------------------------------------------------------
  ClassDef(plotWork,1)

  TH1D *fpHnorm, *fpHpass;

  std::map<int, TH1D*> fYieldRTR, fYieldHLT;

  int fRefTrigger;

};


#endif
