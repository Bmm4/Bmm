#ifndef GUARD_HFDEBUG_H
#define GUARD_HFDEBUG_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h" 
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class TFile;
class TTree;
class TAna01Event;


// ----------------------------------------------------------------------
class HFDebug : public edm::EDAnalyzer {
 public:
  explicit HFDebug(const edm::ParameterSet&);
  ~HFDebug();
  
 private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run &, const edm::EventSetup &);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void endRun(edm::Run const&, edm::EventSetup const&);

  int           fNevt;
  int           fVerbose;
  
  edm::InputTag fTracksLabel;
  std::string   fHLTProcessName; 
  edm::InputTag fL1GTReadoutRecordLabel; 
  edm::InputTag fL1GTmapLabel;
  edm::InputTag fL1MuonsLabel;
  edm::InputTag fTriggerEventLabel;
  edm::InputTag fHLTResultsLabel;
  
  HLTConfigProvider fHltConfig;
  bool              fValidHLTConfig;


  edm::EDGetTokenT<edm::View<reco::Track> > fTokenTrack;

  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord>  fTokenL1GlobalTriggerReadoutRecord;
  edm::EDGetTokenT<trigger::TriggerEvent> fTokenTriggerEvent;
  edm::EDGetTokenT<edm::TriggerResults> fTokenTriggerResults;
  edm::EDGetTokenT<L1GlobalTriggerObjectMap>  fTokenL1GlobalTriggerObjectMap;
  edm::EDGetTokenT<L1GlobalTriggerObjectMapRecord> fTokenL1GlobalTriggerObjectMapRecord;

  edm::Handle<edm::View<reco::Track> > fTracksHandle;

  // new setup
  edm::InputTag                                          fL1MuonsNewLabel,
                                                         fL2MuonsNewLabel,
                                                         fL3MuonsNewLabel;
  edm::EDGetTokenT<l1extra::L1MuonParticleCollection>    l1extramuToken_;
  edm::EDGetTokenT<reco::RecoChargedCandidateCollection> MuCandTag2Token_, MuCandTag3Token_;


};

#endif
