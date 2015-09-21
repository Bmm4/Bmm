#ifndef GENANALYSIS_H
#define GENANALYSIS_H

#include <iostream>
#include <vector>
#include <utility>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include "rootio/TAna01Event.hh"
#include "rootio/TGenCand.hh"
#include "rootio/TAnaCand.hh"
#include "rootio/TAnaTrack.hh"
#include "rootio/TAnaJet.hh"
#include "rootio/TAnaVertex.hh"

#include "treeReader01.hh"

class genAnalysis : public treeReader01 {

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
  void         bbbarCrossSection();

  void         printBdecays(); 
  void         genB(); 
  void         compare2PDG(TGenCand *pCand, int year = 2014, bool finalize = false); 
  double       isolation(TGenCand *pCand); 

  int          NTOTAL; 
  double       XSECTION;
  
  double       f511Mass, f521Mass, f531Mass, f5122Mass, f541Mass; 

};

#endif
