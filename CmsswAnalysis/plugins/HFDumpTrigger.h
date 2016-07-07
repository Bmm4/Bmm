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

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

// -- legacy L1 information
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

// -- 2016 L1 information
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

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
  virtual void resetRun(bool changed);

  int           fVerbose;
  int           fNevt;

  std::string   fHLTProcessName;
  edm::InputTag fL1MuonsLabel;

  edm::InputTag fTriggerEventLabel;
  edm::EDGetTokenT<trigger::TriggerEvent> fTokenTriggerEvent;

  edm::InputTag fHLTResultsLabel;
  edm::EDGetTokenT<edm::TriggerResults> fTokenTriggerResults;

  // -- legacy L1 information
  edm::InputTag fL1TriggerReadoutLabel;
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> fTokenL1TriggerReadoutRecord;


  // -- 2016 L1 information
  edm::InputTag   fAlgInputTag;
  edm::InputTag   fExtInputTag;
  edm::EDGetTokenT<BXVector<GlobalAlgBlk> > fAlgToken;
  edm::EDGetTokenT<BXVector<GlobalExtBlk> > fExtToken;
  l1t::L1TGlobalUtil *fGtUtil;


  HLTPrescaleProvider  fHltPrescaleProvider;

  HLTConfigProvider fHltConfig;
  bool              fValidHLTConfig;


};

#endif
