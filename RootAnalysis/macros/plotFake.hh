#ifndef PLOTFAKE_h
#define PLOTFAKE_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotFake: public plotClass {

public :
  plotFake(std::string dir = "results",
	   std::string files = "plotFake.files",
	   std::string cuts = "plotClass.2016.cuts",
	   std::string setup = "default");
  virtual        ~plotFake();

  void   setCuts(std::string cuts);

  // -- Main analysis methods
  void   loadFiles(std::string afiles);
  void   makeAll(int bitmask = 0);
  void   bookHist(int mode);

  // -- Fake rate from light hadrons
  void   fakeRate(std::string var = "pt", std::string dataset = "data_charmonium", std::string particle = "pion");
  void   fitKs(TH1D*);
  void   fitPhi(TH1D*);
  void   fitLambda(TH1D*);

  void   loopFunction1();


  void   setupTree(TTree *t);

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);
  void   candAnalysis();

private:

  TTree* fTree;

  struct redTreeData fb;

  bool fGoodCand;
  double PTLO;

  double fYield, fYieldE;


  // ----------------------------------------------------------------------
  ClassDef(plotFake,1)

};


#endif
