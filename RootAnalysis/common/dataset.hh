#ifndef DATASET_H
#define DATASET_H

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


class dataset {
public:
  dataset();
  void   cd(std::string dir) {fF->cd(dir.c_str()); }
  TFile* getFile() {return fF;}
  TH1D*  getHist(std::string name, bool clone = true);
  TH2D*  getHist2(std::string name, bool clone = true);
  void   setHistStyle(TH1 *);

  TFile *fF;
  std::string fName, fFullName;
  // -- decay/process specifics
  double fXsec, fBf, fMass, fLambda;
  // -- generation information
  double fFilterEff, fLumi;
  // -- display
  int fColor, fLcolor, fFcolor, fSymbol, fFillStyle;
  double fSize, fWidth;

  int fNclone;
};

#endif
