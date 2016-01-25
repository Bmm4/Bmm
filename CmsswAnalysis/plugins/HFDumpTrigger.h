#ifndef GUARD_HFDUMPTRIGGER_H
#define GUARD_HFDUMPTRIGGER_H

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
//76 #include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

class TFile;
class TTree;
class TAna01Event;


// ----------------------------------------------------------------------
class HFDumpTrigger : public edm::EDAnalyzer {
 public:
  explicit HFDumpTrigger(const edm::ParameterSet&);
  ~HFDumpTrigger();
  
 private:
  virtual void beginJob() ;
  virtual void beginRun(const edm::Run &, const edm::EventSetup &);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void endRun(edm::Run const&, edm::EventSetup const&);

  int           fVerbose;
  int           fNevt;

  std::string   fHLTProcessName; 
  edm::InputTag fL1GTReadoutRecordLabel; 
  edm::InputTag fL1GTmapLabel;
  edm::InputTag fL1MuonsLabel;
  edm::InputTag fTriggerEventLabel;
  edm::InputTag fHLTResultsLabel;
  
  //  HLTPrescaleProvider  fHltPrescaleProvider

  HLTConfigProvider fHltConfig;
  bool              fValidHLTConfig;

  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord>  fTokenL1GlobalTriggerReadoutRecord;
  edm::EDGetTokenT<trigger::TriggerEvent> fTokenTriggerEvent;
  edm::EDGetTokenT<edm::TriggerResults> fTokenTriggerResults;
  edm::EDGetTokenT<L1GlobalTriggerObjectMap>  fTokenL1GlobalTriggerObjectMap;
  edm::EDGetTokenT<L1GlobalTriggerObjectMapRecord> fTokenL1GlobalTriggerObjectMapRecord;
};

#endif
