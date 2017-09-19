#ifndef PLOTCLASS_h
#define PLOTCLASS_h


#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include "common/dataset.hh"
#include "common/selpoint.hh"
#include "common/AnalysisCuts.hh"
#include "common/PidTable.hh"
#include "common/initFunc.hh"

#include "redTreeData.hh"
#include "cuts.hh"
#include "anaNumbers.hh"
#include "ReaderData.hh"

// -- TMVA related
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"


// struct readerData {
//   float pt, eta, m1eta, m2eta, m1pt, m2pt;
//   float fls3d, alpha, maxdoca, pvip, pvips, iso, docatrk, chi2dof, closetrk;
//   float m;
//   float closetrks1, closetrks2, closetrks3;
//   float m1iso, m2iso;
//   float pvdchi2, othervtx;
//   float pv2lip, pv2lips;
// };


// ----------------------------------------------------------------------
class plotClass: public TObject {

public :
  plotClass(std::string dir = "results",
	    std::string files = "plotClass.files",
	    std::string cuts = "nada",
	    std::string setup = "");
  virtual        ~plotClass();

  enum MODE {UNSET, BMM, BDMM, BSMM, RARE, BU2JPSIKP, BD2JPSIKSTAR, BS2JPSIPHI, BS2JPSIF, FAKEKS, FAKEPHI, FAKELAMBDA, FAKEPSI};

  // -- stuff to run over the tree from any derived class
  virtual void   setup(std::string ds);
  TTree*         getTree(std::string ds, std::string dir = "", std::string tree = "events");
  virtual void   setupTree(TTree *t, std::string mode = "");
  virtual void   candAnalysis();
  virtual void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);
  virtual void   loopFunction1();
  virtual void   loopFunction2();

  // -- physics utilities
  virtual double getValueByLabel(TH1D *h, std::string label);
  void           closeHistFile();
  virtual void   readCuts(std::string filename);
  void           readFile(std::string filename, std::vector<std::string> &lines);
  virtual void   printCuts(ostream &OUT);

  virtual void   loadFiles(std::string afiles);
  TFile*         loadFile(std::string afiles);

  void           changeSetup(std::string dir, std::string name, std::string setup);
  void           insertDataset(std::string dsname, dataset *);
  void           muonBdtSetup(TH1D *h1, std::string &prefixB, double &cutB, std::string &prefixE, double &cutE);

  std::string    era(int run);
  int            iera(int run);

  // -- Main analysis methods
  virtual void   makeAll(int bitmask = 0);
  virtual void   init();
  virtual void   treeAnalysis();

  // -- overlays and normalizing histograms
  void           normHist(TH1 *, std::string ds="", int method = NONORM);

  // -- overlay 3
  void           overlay(TH1* h1, std::string f1, TH1 *h2, std::string f2, TH1* h3 = 0, std::string f3 = "",
			 int method = NONORM, bool log = false, bool legend = true, double xleg = 0.4, double yleg = 0.6);
  void           overlay(std::string h1name, std::string f1, std::string h2name, std::string f2, std::string h3name = "", std::string f3 = "",
			 int method = NONORM, bool log = false, bool legend = true, double xleg = 0.4, double yleg = 0.6);

  virtual void   bookHist(std::string name);


  // -- display utilities
  TStyle *       setTdrStyle();
  void           cd(std::string dataset, std::string dir = "");
  void           replaceAll(std::string &sInput, const std::string &oldString, const std::string &newString);
  void           newLegend(double x1, double y1, double x2, double y2, std::string title = "");
  void           makeCanvas(int i = 3);
  void           setItalic();
  void           setRoman();
  void           savePad(std::string name, TCanvas *c = 0);

  // -- histograms
  std::map<std::string, TH1*> fHists;


  // -- cuts
  int fChan, fNchan, fYear;
  std::string fCutsFileName;
  std::vector<cuts*> fCuts;

  double fAccPt, fAccEtaGen, fAccEtaRec;

  // -- BDT reader setup
  //  virtual TMVA::Reader* setupReader(std::string xmlFile, readerData &rd);
  virtual void calcBDT();

  struct redTreeData fb;
  std::vector<TMVA::Reader*> fReaderEvents0;
  std::vector<TMVA::Reader*> fReaderEvents1;
  std::vector<TMVA::Reader*> fReaderEvents2;
  TString fMvaMethod;
  bool fIsMC, fIsSignal;
  double fBDT;
  // readerData frd;
  ReaderData frd;


  // -- PidTables
  PidTable *fptT1, *fptT2, *fptM;
  PidTable *fptT1MC, *fptT2MC, *fptMMC;

  // -- split into seagull and cowboys
  PidTable *fptSgT1, *fptSgT2, *fptSgM;
  PidTable *fptSgT1MC, *fptSgT2MC, *fptSgMMC;

  PidTable *fptCbT1, *fptCbT2, *fptCbM;
  PidTable *fptCbT1MC, *fptCbT2MC, *fptCbMMC;

  PidTable *fptFakePosKaons, *fptFakePosPions, *fptFakePosProtons, *fptPosMuons;
  PidTable *fptFakeNegKaons, *fptFakeNegPions, *fptFakeNegProtons, *fptNegMuons;


  // -- setup and cuts
  double MASSMIN, MASSMAX, SIGBOXMIN, SIGBOXMAX, BGLBOXMIN, BGLBOXMAX, BGHBOXMIN, BGHBOXMAX;

  bool fGoodAcceptance, fPreselection, fPreselectionBDT, fWideMass, fGoodHLT, fGoodMuonsID, fGoodGlobalMuons,
    fGoodBdtPt, fGoodMuonsPt, fGoodMuonsEta, fGoodTracks, fGoodTracksPt, fGoodTracksEta;
  bool fGoodQ, fGoodPvAveW8, fGoodLip, fGoodLipS, fGoodIp, fGoodIpS, fGoodMaxDoca,
    fGoodPt, fGoodEta, fGoodAlpha, fGoodFLS, fGoodChi2, fGoodIso, fGoodM1Iso, fGoodM2Iso;
  bool fGoodCloseTrack, fGoodCloseTrackS1, fGoodCloseTrackS2, fGoodCloseTrackS3,
    fGoodDocaTrk, fGoodJpsiCuts, fGoodCNC, fGoodBDT, fGoodDcand;

  bool fIsCowboy;

  double fW8, fW8MisId, fW8MmuID, fW8Mtrig, fW8DmuID, fW8Dtrig;

  int fRunMin, fRunMax; // if you want to look at a specific run range

  AnalysisCuts fCncCuts;
  AnalysisCuts fBdtCuts;



  std::string fTexFileName;
  std::ofstream fTEX;

  std::string fHistFileName;
  TFile      *fHistFile;

  enum HistNorm {NONORM,     // do not touch the normalization
		 SOMETHING,  // normalize to what is given in fNorm
		 UNITY,      // normalize all to 1
		 XSECTION,   // the resulting histograms will be cross sections [ pb]!!
		 LUMI        // according to the number provided in fLumi       [/fb]!!
  };
  double fNorm, fLumi; // [fLumi] = 1/fb!!!
  number fFsfu;
  double fCrossSection;
  double fBfPsiMuMu, fBfPsiMuMuE,
    fBfPhiKpKm, fBfPhiKpKmE,
    fBfKstarKpPim, fBfKstarKpPimE;

  bool   fDBX, fDoUseBDT;
  int    fVerbose;
  double fEpsilon;

  std::string fDirectory, fSetup, fSuffix, fSuffixSel, fSample, fNumbersFileName, fTreeDir;
  int   fSetupInt;
  // -- datasets (files and associated information)
  std::map<std::string, dataset*> fDS;
  // -- current dataset for analysis
  dataset* fCds;
  initFunc *fIF;

  // -- Display utilities
  std::map<std::string, std::string> fVarToTex, fVarToTexSymbol;
  std::string fStampString, fStampCms, fStampLumi;
  int fFont;
  double fSize;
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5;
  TLatex *tl;
  TBox *box;
  TArrow *pa;
  TLine *pl;
  TLegend *legg;
  TLegendEntry *legge;

  enum MODE fMode;


  static const int MAXPS = 20;

  // ----------------------------------------------------------------------
  ClassDef(plotClass,1)

};

#endif
