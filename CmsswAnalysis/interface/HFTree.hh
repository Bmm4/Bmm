#ifndef _HFTREE_h_
#define _HFTREE_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TFile;
class TTree;
class TH1D;
class TAna01Event;

// ----------------------------------------------------------------------
class HFTree : public edm::EDAnalyzer {
 public:
  explicit HFTree(const edm::ParameterSet&);
  ~HFTree();

 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  TFile        *fFile;
  TTree        *fTree;
  TAna01Event  *fEvent;
  TH1D         *fH1;

  bool         fRequireCand, fFullGenBlock;

  std::string  fFileName;
  std::string  fTreeName;

  int          fVerbose;
  int          fPrintFrequency;
  int          fEventCounter;

};

#endif
