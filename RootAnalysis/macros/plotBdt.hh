#ifndef PLOTBDT_h
#define PLOTBDT_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"


struct bdtSetup {
  int NTrees, nEventsMin, MaxDepth, NNodesMax;
  int nCuts;
  float AdaBoostBeta;
};


// ----------------------------------------------------------------------
class plotBdt: public plotClass {

public :
  plotBdt(std::string dir = "results",
	   std::string files = "plotResults.2016.files",
	   std::string cuts = "baseCuts.cuts",
	   std::string setup = "");
  virtual        ~plotBdt();

  // -- Main analysis methods
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void makeAll(std::string what = "all");
  virtual void bookHist(std::string dsname);

  void  make(int offset, string filename, int evt, int clean);
  void  createInputFiles(std::string filename = "tmva-trees.root", int randomseed = -1);
  void  setupTree(TTree *t, redTreeData &b);
  void  train(std::string oname = "TMVA-0", std::string filename = "/scratch/ursl/tmva-trees.root", int nsg = -1, int nbg = -1);
  void  apply(const char *fname);

  TH1D* getRanking(std::string fname, std::string prefix, std::string type);
  void  mvas(std::string fname);
  void  analyze(std::string fname);
  void  cleanup(std::string fname);

  void  setVars(string vars) {fVariables = vars; }
  void  setBDTParameters(string pars) {fBDTParameters = pars; }
  void  setApply0() {fApplyOn0 = true;  fApplyOn1 = false; fApplyOn2 = false;};
  void  setApply1() {fApplyOn0 = false; fApplyOn1 = true;  fApplyOn2 = false;};
  void  setApply2() {fApplyOn0 = false; fApplyOn1 = false; fApplyOn2 = true; };
  void  setTrainAntiMuon(bool yes) {fTrainAntiMuon = yes;};
  void  setChannel(int channel) {fChannel = channel;};

  void  writeOut(TFile*, TH1*);
  void  redrawStats(double x, double y, const char *newname, int color);

  // -- code for loops
  void  loopFunction1();
  void  loopFunction2();
  void  loopFunction3();
  void  loopFunction4();
  void  loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

private:

  double fLumiScale;

  int  fChannel;
  std::vector<TMVA::Reader*> fReader;
  double fBDT, fBDT0, fBDT1, fBDT2;
  bool   fApplyOn0, fApplyOn1, fApplyOn2;

  std::string fVariables, fBDTParameters;
  bdtSetup    fBdtSetup;
  bool        fTrainAntiMuon;


  ClassDef(plotBdt,1)

};


#endif
