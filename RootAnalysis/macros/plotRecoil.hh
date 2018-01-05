#ifndef PLOTRECOIL_h
#define PLOTRECOIL_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotRecoil: public plotClass {

public :
  plotRecoil(std::string dir = "results",
	     std::string files = "plotRecoil.files",
	     std::string cuts = "plotClass.2016.cuts",
	     std::string setup = "default",
	     int year = 0);
  virtual        ~plotRecoil();

  void   setCuts(std::string cuts);

  // -- Main analysis methods
  void   loadFiles(std::string afiles);
  void   makeAll(int bitmask = 0);
  void   bookHist(int mode);

  void   recoil0();
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
  ClassDef(plotRecoil,1)

};


#endif
