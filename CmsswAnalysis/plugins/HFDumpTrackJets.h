#ifndef  HFDumpTrackJetsH
#define  HFDumpTrackJetsH

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


class TFile;
class TTree;
class TAna01Event;

class TrackAssociatorBase;

// ----------------------------------------------------------------------
class HFDumpTrackJets : public edm::EDAnalyzer {
 public:
  explicit HFDumpTrackJets(const edm::ParameterSet&);
  ~HFDumpTrackJets();
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int           fVerbose;
  int           fDoFlavorTagging;
  std::string   fJetsLabel;
  std::string   fJetsTagLabel;
  std::string   fTracksLabel;
  std::string   fGenCandidatesLabel;
  edm::InputTag fsourceByRefer;
 

 
  int nevt;
 

  

};

#endif
