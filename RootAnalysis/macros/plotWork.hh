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

  void plotTisEfficiency(std::string dsname);
  void runTisEfficiency(std::string dsname);


  void bookHist(std::string dsname);
  TH1D *fpHnorm, *fpHpass;


  /// NOTE: This works with the output of genAnalysis!!!!
  void prodSummary(string ds1, int year);


  // -- code for loops
  void   loopFunction1();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

private:


  // ----------------------------------------------------------------------
  ClassDef(plotWork,1)

};


#endif
