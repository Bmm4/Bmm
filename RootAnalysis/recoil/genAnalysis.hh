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
  void         plotBdecays(int btype = 531);
  void         printB2JpsiXdecays();
  void         recoilValidation();
  bool         decayModeValidation(TGenCand *pGen, int mode);
  bool         isAncestor(TGenCand *pAnc, TGenCand *pCand);

  int CANDTYPE, CANDTRUTH;
  int BRECOTYPE, BRECOTRUTH;

  int fRecId, fRecMode, fCanId, fCanMode;
  double fRecPt, fRecEta, fRecPhi, fCanPt, fCanEta, fCanPhi;

  double fMu1Pt, fMu1Eta, fMu1Phi, fMu2Pt, fMu2Eta, fMu2Phi;
  double fKa1Pt, fKa1Eta, fKa1Phi, fKa2Pt, fKa2Eta, fKa2Phi;

  std::map<int, std::vector<int> > fDecayModes;

};

#endif
