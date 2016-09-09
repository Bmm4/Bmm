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

  // -- code for loops
  void   loopFunction1();
  void   loopFunction2();
  void   loopFunction3();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

private:


  // ----------------------------------------------------------------------
  ClassDef(plotResults,1)

};


#endif
