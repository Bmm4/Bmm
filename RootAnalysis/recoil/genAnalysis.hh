#ifndef GENANALYSIS_H
#define GENANALYSIS_H

#include <iostream>
#include <vector>
#include <utility>

#include <TString.h>
#include <TChain.h>

#include "t1Reader.hh"

class genAnalysis : public t1Reader {

public:
  genAnalysis(TChain *tree, TString evtClassName);
  ~genAnalysis();

  void         bookHist();
  void         startAnalysis();
  void         eventProcessing();
  void         endAnalysis();
  void         fillHist();
  void         readCuts(TString filename, int dump = 1);
  void         initVariables();

  int          muonType(TGenCand *pCand);
  void         printBdecays();
  void         printB2JpsiXdecays();
  void         recoilValidation();
  bool         decayModeValidation(TGenCand *pGen, int mode);

  int CANDTYPE, CANDTRUTH;
  int BRECOTYPE, BRECOTRUTH;

  int fRecId, fRecMode, fCanId, fCanMode;
  double fRecPt, fRecEta, fRecPhi, fCanPt, fCanEta, fCanPhi;

  std::map<int, std::vector<int> > fDecayModes;

};

#endif
