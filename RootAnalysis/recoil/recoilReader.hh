#ifndef RECOILREADER_H
#define RECOILREADER_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <set>
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

#include "t1Reader.hh"

#define DR      57.29577951


class candAna;

// -- Class derived from t1Reader
//    setup various candAna* classes to analyse an event stored in a T1 tree
class recoilReader : public t1Reader {

public:
  recoilReader(TChain *tree, TString evtClassName);
  ~recoilReader();

  virtual void   startAnalysis();
  virtual void   endAnalysis();
  virtual void   eventProcessing();
  virtual void   readCuts(TString filename, int dump = 1);
  virtual void   bookHist();

  std::vector<candAna*> lCandAnalysis;

  std::set<int> fCandTypes;
};

#endif
