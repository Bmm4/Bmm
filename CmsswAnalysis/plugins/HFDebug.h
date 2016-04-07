#ifndef GUARD_HFDEBUG_H
#define GUARD_HFDEBUG_H

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

// ----------------------------------------------------------------------
class HFDebug : public HFVirtualDecay {
 public:
  explicit HFDebug(const edm::ParameterSet&);
  ~HFDebug();
  
 protected:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run &, const edm::EventSetup &);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void endRun(edm::Run const&, edm::EventSetup const&);

  virtual void dumpConfiguration();
  double matchTrigger(std::vector<int> &trigIndices, const trigger::TriggerObjectCollection &trigObjs, 
		      edm::Handle<trigger::TriggerEvent> &triggerEvent, const reco::Muon &mu);
  TrajectoryStateOnSurface cylExtrapTrkSam(reco::TrackRef track, double rho);
  TrajectoryStateOnSurface surfExtrapTrkSam(reco::TrackRef track, double z);
  FreeTrajectoryState freeTrajStateMuon(reco::TrackRef track);

  
  edm::InputTag fTriggerEventLabel;
  edm::EDGetTokenT<trigger::TriggerEvent> fTokenTriggerEvent;

  edm::InputTag fHLTResultsLabel;
  edm::EDGetTokenT<edm::TriggerResults> fTokenTriggerResults;

  std::string   fHLTProcessName; 

  std::vector<std::string> fTriggerNames;
  std::vector<int>         fTriggerIndices;

  edm::ESHandle<Propagator> fPropagatorAlong;
  edm::ESHandle<Propagator> fPropagatorOpposite;
  
  HLTConfigProvider fHltConfig;
  bool              fValidHLTConfig;

};

#endif
