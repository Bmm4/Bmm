#ifndef PLOTWORK_h
#define PLOTWORK_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotWork: public plotClass {

public :
  plotWork(std::string dir = "results", std::string files = "plotWork.files", std::string setup = "default");
  virtual        ~plotWork();

  void   setCuts(std::string cuts);

  // -- Main analysis methods
  void   loadFiles(std::string afiles);
  void   makeAll(int bitmask = 0);
  void   bookHist(int mode);

  // -- validate gen production
  void   prodSummary(string ds1, int year = 2014);

  // -- code for loops
  void   loopFunction1();


  void   setupTree(TTree *t);

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);
  void   candAnalysis();

private:

  TTree* fTree;

  struct redTreeData fb;


  bool fGoodCand;
  double PTLO;


  // ----------------------------------------------------------------------
  ClassDef(plotWork,1)

};


#endif
