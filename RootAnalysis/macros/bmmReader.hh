#ifndef BMMREADER_H
#define BMMREADER_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>

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

#include "common/JSON.hh"

#include "treeReader01.hh"

#define DR      57.29577951

class PidTable; 
class candAna;

class bmmReader : public treeReader01 {

public:
  bmmReader(TChain *tree, TString evtClassName);
  ~bmmReader();

  virtual void   startAnalysis();
  virtual void   endAnalysis();
  virtual void   eventProcessing();
  virtual void   readCuts(TString filename, int dump = 1);
  virtual void   bookHist();

  virtual void   processType();
  virtual void   setYear(int year) {fYear = year;}

  std::vector<candAna*> lCandAnalysis;


  int fProcessType; 
  int fYear; 
};

#endif
