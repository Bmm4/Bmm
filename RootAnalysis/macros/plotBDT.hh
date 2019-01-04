#ifndef PLOTBDT_h
#define PLOTBDT_h

#include "plotClass.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotBDT: public plotClass {

public :
  plotBDT(std::string dir = "results",
	  std::string files = "plotResults.2016.files",
	  std::string cuts = "nada",
	  std::string setup = "",
	  int year = 0);
  virtual ~plotBDT();

  // -- Main analysis methods
  virtual void init();
  //  virtual void loadFiles(std::string afiles);
  virtual void makeAll(std::string what = "all");

  void tmvaControlPlots(int ichan);
  void dumpParameters(std::string type);
  void tmvaPlots(std::string type);

  void apply(std::string what = "ana");
  void opspt(std::string what = "hi");
  TH1D* ssb(TH2D *hS, TH2D *hB, std::string what = "hi");


  std::string replaceLabelWithTex(string label);
  void SetFrameStyle( TH1* frame, Float_t scale = 1.0);
  void NormalizeHists( TH1* sig, TH1* bkg = 0);
  void SetSignalAndBackgroundStyle( TH1* sig, TH1* bkg, TH1* all = 0);
  void GetMethodTitle( TString & name, TKey * ikey);
  void GetMethodName( TString & name, TKey * mkey);
  void GetMethodTitle( TString & name, TDirectory * idir);
  int  GetNumberOfTargets( TDirectory *dir);
  int  GetNumberOfInputVariables( TDirectory *dir);
  std::string parseXmlOption(std::string line);
  int  bdtString2Channel(std::string s);
  void setBdtStrings(int ichan);
  void readLogFile(std::string sfile);

  void xmlParsingVariables(std::string weightfile);
  void xmlParsingReadTree(std::string xmlfile);
  double getMaximum(TH1 *h1, TH1 *h2);


  // -- migrated from tmva1 (TLF = TMVA logfile)
  void getTLFEventNumbers();
  void getTLFRanking(std::string prefix = "BDT", std::string type = "events0");
  void getTLFParameters(std::string prefix = "BDT", std::string type = "events0");

  // -- analysis of BDT optimization
  void bdtOptMakeTree(std::string logfile);
  void bdtOptAnaTree(std::string rootfile);

  // -- code for loops
  void   loopFunction1();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);

private:
  TFile *fRootFile;

  std::string fBdtString, fBdtLogFile;
  std::vector<std::string> fLogFileLines;
  std::map<std::string, std::vector<std::string> > fReaderVariables;

  // -- XML parsing
  std::map<int, std::string> fVariableMap;
  TH1D *fhBdtVariables, *fhBdtVariablesW8, *fhBdtNodes, *fhBdtNodesW8;
  std::vector<TH1D*> fhBdtVariableCuts;
  std::vector<TH1D*> fhBdtVariableCutsW8;

  double fBDT0, fBDT1, fBDT2;
  ClassDef(plotBDT,1)

};


#endif
