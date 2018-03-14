#ifndef TMVA1
#define TMVA1

#include <iostream>
#include <fstream>
#include <vector>
#include <map>


#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"

#include "redTreeData.hh"
#include "preselection.hh"
#include "common/util.hh"

#include "ReaderData.hh"

struct files {
  std::string dname;
  std::string sname;
  std::string rname;
};


class tmva1: public TObject {

  public:

  tmva1(int year = 2016, std::string vars = "pt:eta:alpha:fls3d:maxdoca:pvip:pvips:iso:m1iso:m2iso",
	std::string pars = "NTrees=800");
  ~tmva1();

  TCanvas* getC0();
  void train(std::string oname = "TMVA-0", std::string filename = "/scratch/ursl/bdt/tmva-trees-0-2016.root", int nsg = -1, int nbg = -1,
	     std::string cut = "");
  void makeAll(int offset, std::string filename = "", int chan = -1);
  void make(int offset, std::string filename, int evt);

  void createInputFile(std::string fname, std::string sfile, std::string dfile, int randomSeed = -1);

  void setBDTParameters(std::string pars) {fBDTParameters = pars;}
  void newLegend(double x1, double y1, double x2, double y2, std::string title = "");
  void writeOut(TFile*, TH1*);

  void setApply0() {fApplyOn0 = true;  fApplyOn1 = false; fApplyOn2 = false;};  // apply on even events, train on odd
  void setApply1() {fApplyOn0 = false; fApplyOn1 = true;  fApplyOn2 = false;};  // apply on odd events, train on even
  void setApply2() {fApplyOn0 = false; fApplyOn1 = false; fApplyOn2 = true; };  // apply on odd events, train on even
  void setTrainAntiMuon(bool yes) {fTrainAntiMuon = yes;};
  void setChannel(int channel) {fChannel = channel;};

  files fInputFiles;

  std::vector<double> fKS, fMaxSSB, fMaxBdt;
  TH1D *fH1s, *fH1b, *fH1r;

  bool fApplyOn0, fApplyOn1, fApplyOn2;
  bool fTrainAntiMuon;
  int fChannel, fYear;
  double fRsigma, fLumiScale;
  std::string fVariables, fBDTParameters;

  redTreeData ftd;
  ReaderData frd;
  double fBDT, fBDT0, fBDT1, fBDT2;
  std::vector<TMVA::Reader*> fReader;
  presel fPresel;

  TLegend *legg;
  TLegendEntry *legge;
  TLatex *tl;

  ClassDef(tmva1,1) //Testing tmva1
};

// void setTitles(TH1 *h, const char *sx, const char *sy,
// 	       float size = 0.05, float xoff = 1.1, float yoff = 1.1, float lsize = 0.05, int font = 42);
// void setHist(TH1 *h, int color = kBlack, int symbol = 20, double size = 1., double width = 2.);
// void shrinkPad(double b = 0.1, double l = 0.1, double r = 0.1, double t = 0.1);

#endif
