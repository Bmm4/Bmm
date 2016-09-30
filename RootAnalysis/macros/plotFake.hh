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

  // -- Main analysis methods
  virtual void bookHist(int mode);
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void makeAll(int bitmask = 0);

  // -- Fake rate from light hadrons
  void   playKs(std::string cuts = "flxy<4&&flsxy>15&&pvips<3&&prob>0.01&&abs(m1ips)>5&&abs(m2ips)>5&&m1q*m2q<0");
  void   playPhi(std::string cuts = "m1q*m2q<0&&prob>0.01&&maxdoca<0.002&&m2pt>4&&m1pv==m2pv&&m1bpixl1&&m2bpixl1");
  void   playLambda(std::string cuts = "flxy<4&&flsxy>15&&fls3d>15&&chi2<3&&pvips<3");
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
