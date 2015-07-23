#ifndef  HFInclBMuonTrackJetsH
#define  HFInclBMuonTrackJetsH

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class TFile;
class TTree;
class TAna01Event;

const double mmuon = 0.10565837;

// ----------------------------------------------------------------------
class HFInclBMuonTrackJets : public edm::EDAnalyzer {
 public:
  explicit HFInclBMuonTrackJets(const edm::ParameterSet&);
  ~HFInclBMuonTrackJets();
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int         fVerbose;
  double      fJetMatch;
  double      fJetEtMin; 
  std::string fMuonLabel; 
  std::string fJetsLabel;
  std::string fTracksLabel;
  std::string fVertexLabel;
  std::string fSimVertexLabel;
  double      fVertexMinNdof;
  bool        fAllowGlobalOnly; 
  
  int         fNevt;

};

#endif
