#ifndef BMMREADER_H
#define BMMREADER_H

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

#include "treeReader01.hh"

#define DR      57.29577951


struct hltPathInfo {
  bool result, wasRun, error, v;
  int  prescale;
};


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

  virtual void   processTypePythia6();
  virtual void   processTypePythia8();

  virtual std::string splitTrigRange(std::string tl, int &r1, int &r2);

  std::vector<candAna*> lCandAnalysis;

  int fProcessType;
  std::set<int> fCandTypes;
  std::map<std::string, hltPathInfo> fHltPathInfo;
};

#endif
