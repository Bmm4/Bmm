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

#include "dataset.hh"
#include "selpoint.hh"

// ----------------------------------------------------------------------
class plotClass: public TObject {

public :
                 plotClass(std::string dir = "hpt0", std::string files = "plotLq.files", std::string setup = "m");
  virtual        ~plotClass();
  void           closeHistFile();

  virtual void   loadFiles(std::string afiles);
  TFile*         loadFile(std::string afiles);

  // -- Main analysis methods 
  virtual void   makeAll(int bitmask = 0);
  virtual void   treeAnalysis(); 

  // -- overlays and normalizing histograms
  void           normHist(TH1 *, std::string ds="", int method = NONORM); 
  virtual void   overlayAll();
  // -- overlay 2
  void           overlay(TH1* h1, std::string f1, TH1 *h2, std::string f2, 
			 int method = NONORM, bool log = false, bool legend = true, double xleg = 0.4, double yleg = 0.6);
  void           overlay(std::string h1name, std::string f1, std::string h2name, std::string f2, 
			 int method = NONORM, bool log = false, bool legend = true, double xleg = 0.4, double yleg = 0.6);

  // -- overlay 3
  void           overlay(TH1* h1, std::string f1, TH1 *h2, std::string f2, TH1* h3, std::string f3, 
			 int method = NONORM, bool log = false, bool legend = true, double xleg = 0.4, double yleg = 0.6);
  void           overlay(std::string h1name, std::string f1, std::string h2name, std::string f2, std::string h3name, std::string f3, 
			 int method = NONORM, bool log = false, bool legend = true, double xleg = 0.4, double yleg = 0.6);

  virtual void   bookHist(std::string name); 

  TTree*         getTree(std::string ds, std::string dir = ""); 
  virtual void   setupTree(TTree *t); 
  virtual void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  virtual void   loopFunction1(); 
  virtual void   loopFunction2(); 

  void           cd(std::string dataset, std::string dir = "");
  void           replaceAll(std::string &sInput, const std::string &oldString, const std::string &newString);
  void           newLegend(double x1, double y1, double x2, double y2, std::string title = "");
  void           makeCanvas(int i = 3);


  // -- histograms
  std::map<std::string, TH1*> fHists;

protected:

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

  bool   fDBX;
  int    fVerbose; 
  double fEpsilon; 
  bool   fIsMC;

  std::string fDirectory, fSetup, fSuffix;   

  // -- datasets (files and associated information)
  std::map<std::string, dataset*> fDS; 
  // -- current dataset for analysis
  std::string fCds; 

  // -- Display utilities
  int fFont; 
  double fSize; 
  TCanvas *c0, *c1, *c2, *c3, *c4, *c5;
  TLatex *tl; 
  TBox *box;
  TArrow *pa; 
  TLine *pl; 
  TLegend *legg;
  TLegendEntry *legge;


  // ----------------------------------------------------------------------
  ClassDef(plotClass,1) 

};

#endif
