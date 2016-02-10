// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFDebug
// ------------
//
// 2016/01/27  Urs Langenegger      first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDebug.h"

#include <sstream>
#include <iostream>
#include <bitset>

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

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TH1.h"
#include "TFile.h"
#include "TDirectory.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"
#include "Bmm/RootAnalysis/rootio/TAnaMuon.hh"
#include "Bmm/RootAnalysis/rootio/TTrgObj.hh"
#include "Bmm/RootAnalysis/rootio/TTrgObjv2.hh"



// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;


// ----------------------------------------------------------------------
HFDebug::HFDebug(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fNevt(0),  

  fTracksLabel(iConfig.getUntrackedParameter<edm::InputTag>("tracksLabel", edm::InputTag("generalTracks"))),
  fTokenTrack(consumes<View<Track> >(fTracksLabel)), 

  fL1MuonsLabel(iConfig.getUntrackedParameter<InputTag>("L1MuonsLabel")),

  fHLTProcessName(iConfig.getUntrackedParameter<string>("HLTProcessName")),

  fTriggerEventLabel(iConfig.getUntrackedParameter<InputTag>("TriggerEventLabel")),
  fTokenTriggerEvent(consumes<TriggerEvent>(fTriggerEventLabel)), 

  fHLTResultsLabel(iConfig.getUntrackedParameter<InputTag>("HLTResultsLabel")),
  fTokenTriggerResults(consumes<TriggerResults>(fHLTResultsLabel)),

  fHltPrescaleProvider(iConfig, consumesCollector(), *this),

  fL1MuonsNewLabel(iConfig.getUntrackedParameter<edm::InputTag>("l1muonsLabel", edm::InputTag("hltL1extraParticles"))),
  fL2MuonsNewLabel(iConfig.getUntrackedParameter<edm::InputTag>("l2muonsLabel", edm::InputTag("hltL2MuonCandidates"))),
  fL3MuonsNewLabel(iConfig.getUntrackedParameter<edm::InputTag>("l3muonsLabel", edm::InputTag("hltL3MuonCandidates"))) {  

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDebug constructor" << endl;
  cout << "--- Verbose                     : " << fVerbose << endl;
  cout << "--- HLT process name            : " << fHLTProcessName << endl;
  cout << "--- L1 Muons Label              : " << fL1MuonsLabel << endl;
  cout << "--- HLTResultsLabel             : " << fHLTResultsLabel << endl;
  cout << "--- Trigger Event Label         : " << fTriggerEventLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;


  l1extramuToken_  = consumes<l1extra::L1MuonParticleCollection>(fL1MuonsNewLabel); 
  MuCandTag2Token_ = consumes<reco::RecoChargedCandidateCollection>(fL2MuonsNewLabel);
  MuCandTag3Token_ = consumes<reco::RecoChargedCandidateCollection>(fL3MuonsNewLabel);
}


// ----------------------------------------------------------------------
HFDebug::~HFDebug() {
  
}


// ----------------------------------------------------------------------
void HFDebug::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  fNevt++;

  cout << "======================================================================" << endl;
  stringstream evtstring;
  evtstring << "run: " << iEvent.id().run()
	    << " event: " << iEvent.id().event()
	    << " LS: " <<  iEvent.luminosityBlock();
  cout << evtstring.str() << endl;
  cout << "======================================================================" << endl;

  
  // ----------------------------------------------------------------------
  // -- tracks
  // ----------------------------------------------------------------------
  try {
    iEvent.getByToken(fTokenTrack, fTracksHandle);
    fTracksHandle.isValid();
  } catch(cms::Exception&){
    cout << "could not get tracks, giving up" << endl;
    return;
  }

  int i(0), hiPurity(0);
  for (View<reco::Track>::const_iterator itTrack = fTracksHandle->begin(); itTrack != fTracksHandle->end(); ++itTrack) {

    reco::TrackBase::TrackQuality trackQualityhighPur = reco::TrackBase::qualityByName("highPurity");
    if (itTrack->quality(trackQualityhighPur)) {
      hiPurity = 1;
    } else {
      hiPurity = 0;
    }

    cout << Form("%3d", i) << " track p/t = "
	 << itTrack->pt() << ", " << itTrack->eta() << ", " << itTrack->phi()
	 << " algo = " << itTrack->algo()
	 << " high purity = " << hiPurity
	 << endl;
    ++i; 
  }

    
  
  // ----------------------------------------------------------------------
  // -- HLT results
  // ----------------------------------------------------------------------
  vector<string> validTriggerNames;
  if (fValidHLTConfig) validTriggerNames = fHltConfig.triggerNames();
  else cerr << "==> HFDebug: No valid Trigger configuration!!!" << endl;
  //can assert?!  fHltConfig.dump("PrescaleTable");

  if (validTriggerNames.size() < 1) {
    cout << "==>HFDebug: NO valid trigger names returned by HLT config provided!!??" << endl;
    return;
  }

  Handle<TriggerResults> hHLTresults;
  bool hltF = true;

  try {
    iEvent.getByToken(fTokenTriggerResults, hHLTresults);
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==>HFDebug> Triggerresults  " << fHLTResultsLabel.encode() << " not found " << endl;
    hltF = false;
    return;
  }
  
  Handle<trigger::TriggerEvent> triggerSummary;
  hltF = true;
  try {
    iEvent.getByToken(fTokenTriggerEvent, triggerSummary);
    hltF = true;
  } catch (const cms::Exception& e) {
    hltF = false;
    cout<<"Error!! No TriggerEvent with label " << fTriggerEventLabel << endl;
    return;
  }

  if (!hltF) return;
  
  const TriggerNames &trigName = iEvent.triggerNames(*hHLTresults);
  cout << "hHLTresults->size() = " << hHLTresults->size() << " and HLT accept = " << gHFEvent->fHLTDecision << endl;
  
  unsigned int index(999); 
  bool wasrun(false), result(false);
  int prescale(1); 
  int psSet = -1;
  psSet = fHltPrescaleProvider.prescaleSet(iEvent, iSetup);
  
  // Loop over all HLT-paths
  for (unsigned int it = 0; it < validTriggerNames.size(); ++it) {
    index    = trigName.triggerIndex(validTriggerNames[it]); 
    const unsigned int triggerIndex = fHltConfig.triggerIndex(validTriggerNames[it]); //dk
    if(index!=triggerIndex) cout<<" something wrong 1 "<<index<<" "<<triggerIndex<<endl;
    
    if (index >= hHLTresults->size())
      continue;
    
    TString validTriggerNamesT = TString(validTriggerNames[it]);
    result   = hHLTresults->accept(index);
    wasrun   = hHLTresults->wasrun(index);
    if (psSet > -1) {
      prescale = fHltConfig.prescaleValue(psSet, validTriggerNames[it]);
    } else {
      //	cout << "==>HFDebug> error in prescale set!?" << endl;
      prescale = 0;
    }
    
    const vector<string>& moduleLabels(fHltConfig.moduleLabels(index));
    const unsigned int moduleIndex(hHLTresults->index(index));
    
    cout << "--- new HLT path " << evtstring.str() << " ------------" << endl;
    cout << validTriggerNames[it] << " result = " << result << " "
	 << " prescale = " << prescale
	 << " wasrun = " << wasrun
	 << " moduleLabels.size() = "  << moduleLabels.size()
	 << " moduleIndex = "  << moduleIndex
	 << " triggerIndex = " << triggerIndex
	 << " sizeFilters = " << triggerSummary->sizeFilters()
	 << endl;
    
    // loop over all modules 
    for (unsigned int j=0;j<=moduleIndex;j++) {
      const string& moduleLabel = moduleLabels[j]; 
      const string& type = fHltConfig.moduleType(moduleLabel);
      const unsigned int filterIndex = triggerSummary->filterIndex(InputTag(moduleLabel,"","HLT"));
      cout << "-> " << j << " label " << moduleLabel << " type " << type << " index "
	   << filterIndex
	   << endl;
      
      if (filterIndex < triggerSummary->sizeFilters()) {
	
	const trigger::Vids& VIDS (triggerSummary->filterIds(filterIndex));
	const trigger::Keys& KEYS(triggerSummary->filterKeys(filterIndex));
	const size_type nI(VIDS.size());
	const size_type nK(KEYS.size());
	assert(nI==nK);
	const size_type n(max(nI,nK));
	
	cout << "   " << j << " label " << moduleLabel << " type " << type << " index "
	     << filterIndex << " " << n << " accepted TRIGGER objects found: "
	     << endl;
	
	const trigger::TriggerObjectCollection& TOC = triggerSummary->getObjects();
	for (size_type i=0; i!=n; ++i) {
	  const trigger::TriggerObject& TO=TOC[KEYS[i]];
	  cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
	       << TO.id() << " " << TO.pt() << " " << TO.eta() << " " 
	       << TO.phi() << " " << TO.mass()<< endl;                                                                 
	} 
      } 
    } // for j, modul loop
  } 

  // -- HLT L1 muon candidates
  string L1NameCollection("hltL1extraParticles");
  trigger::size_type Index(0);
  Index = triggerSummary->collectionIndex(edm::InputTag(L1NameCollection, "", fTriggerEventLabel.process()));
  
  if (Index < triggerSummary->sizeCollections()) {
    TString label = TString(L1NameCollection.c_str());
    const trigger::Keys& Keys(triggerSummary->collectionKeys());
    const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
    const trigger::size_type n1 (Keys.at(Index));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj( triggerSummary->getObjects().at(i) );
      cout << "===> L1: " << i << " id: " << obj.id() 
	   << " m = " << obj.mass() << " pT,eta,phi = " << obj.pt() << "," 
	   <<  obj.eta() << "," << obj.phi()
	   << " label = " << label
	   << endl;
      
      TTrgObj *pTO = gHFEvent->addTrgObj();
      pTO->fP.SetPtEtaPhiE(obj.pt(), 
			   obj.eta(), 
			   obj.phi(), 
			   obj.energy()
			   ); 
      pTO->fID     = obj.id(); 
      pTO->fLabel  = label;
      
    } 
    cout << "===> Found L1 trigger collection -> " << L1NameCollection << " " << (n1-n0)
	 << endl;
  }

  // -- HLT L2 muon candidates 
  string L2NameCollection("hltL2MuonCandidates");
  Index = triggerSummary->collectionIndex(edm::InputTag(L2NameCollection, "", fTriggerEventLabel.process()));
  
  if (Index < triggerSummary->sizeCollections()) {
    TString label = TString(L2NameCollection.c_str());
    const trigger::Keys& Keys(triggerSummary->collectionKeys());
    const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
    const trigger::size_type n1 (Keys.at(Index));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj( triggerSummary->getObjects().at(i) );
      if(fVerbose>10) 
	cout << "===> L2: " << i << " id: " << obj.id() << " m = " 
	     << obj.mass() << " pT,eta,phi = " << obj.pt() << "," 
	     <<  obj.eta() << "," << obj.phi()
	     << " label = " << label
	     << endl;
      
      TTrgObj *pTO = gHFEvent->addTrgObj();
      pTO->fP.SetPtEtaPhiE(obj.pt(), 
			   obj.eta(), 
			   obj.phi(), 
			   obj.energy()
			   ); 
      pTO->fID     = obj.id(); 
      pTO->fLabel  = label;
    }
    cout << "===> Found L2 trigger collection -> " << L2NameCollection << " " 
	 << (n1-n0)
	 << endl;
  }
  
  // -- HLT L3 muon candidates
  string L3NameCollection("hltL3MuonCandidates");
  Index = triggerSummary->collectionIndex(edm::InputTag(L3NameCollection, "", fTriggerEventLabel.process()));
  if (Index < triggerSummary->sizeCollections()) {
    TString label = TString(L3NameCollection.c_str());
    const trigger::Keys& Keys(triggerSummary->collectionKeys());
    const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
    const trigger::size_type n1 (Keys.at(Index));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj( triggerSummary->getObjects().at(i) );
      cout << "===> L3: " << i << " id: " << obj.id() << " m = " 
	   << obj.mass() << " pT,eta,phi = " << obj.pt() << "," 
	   <<  obj.eta() << "," << obj.phi()
	   << " label = " << label
	   << endl;
      
      TTrgObj *pTO = gHFEvent->addTrgObj();
      pTO->fP.SetPtEtaPhiE(obj.pt(), 
			   obj.eta(), 
			   obj.phi(), 
			   obj.energy()
			   ); 
      pTO->fID     = obj.id(); 
      pTO->fLabel  = label;
      
    }
    cout << "===> Found L3 trigger collection -> " << L3NameCollection << " " 
	 << (n1-n0)
	 << endl;
  }


  // -- new setup, based on
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/HLTrigger/HLTanalyzers/interface/HLTAnalyzer.h
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/HLTrigger/HLTanalyzers/src/HLTAnalyzer.cc
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/HLTrigger/HLTanalyzers/src/HLTMuon.cc
  // https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/HLTrigger/HLTanalyzers/python/HLTAnalyser_cfi.py
  /*  
  Handle<l1extra::L1MuonParticleCollection> l1extmu; // hltL1extraParticles
  Handle<RecoChargedCandidateCollection> mucands2, mucands3; // hltL2MuonCandidates hltL3MuonCandidates

  try {
    iEvent.getByToken(l1extramuToken_, l1extmu);
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==>HFDebug> l1extra::L1MuonParticleCollection  " << fHLTResultsLabel.encode() << " not found " << endl;
    return;
  }

  try {
    iEvent.getByToken(MuCandTag2Token_, mucands2);
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==>HFDebug> RecoChargedCandidateCollection  " << fL3MuonsNewLabel.encode() << " not found " << endl;
    return;
  }

  typedef reco::RecoChargedCandidateCollection::const_iterator cand;
  int idx(0); 
  //  reco::RecoChargedCandidateCollection myMucands2 = ;
  for (cand i = (*mucands2).begin(); i != (*mucands2).end(); i++) {
    reco::TrackRef tk = i->get<reco::TrackRef>();
    cout << "L2 " << idx << " pt = " << tk->pt() << " eta = " << tk->eta() << " phi = " <<  tk->phi() << endl;
  }

  
  try {
    iEvent.getByToken(MuCandTag2Token_, mucands3);
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==>HFDebug> RecoChargedCandidateCollection  " << fL3MuonsNewLabel.encode() << " not found " << endl;
    return;
  }
  */

  
  
} 


// ----------------------------------------------------------------------
void  HFDebug::beginRun(const Run &run, const EventSetup &iSetup) {
  bool hasChanged;
  fValidHLTConfig = fHltConfig.init(run,iSetup,fHLTProcessName,hasChanged);
  cout << fHltConfig.tableName() << endl;
  vector<string> pds = fHltConfig.datasetNames();
  
  TDirectory *pDir = gDirectory; 
  gHFFile->cd();
  TH1D *h1 = new TH1D(Form("pd_run%d", run.run()), fHltConfig.tableName().c_str(), pds.size(), 0., pds.size()); 
  h1->SetDirectory(gHFFile); 

  for (unsigned int i = 0; i < pds.size(); ++i) {
    cout << "   " << pds[i] << endl;
    h1->GetXaxis()->SetBinLabel(i+1, pds[i].c_str()); 
  }

  h1->Write();
  delete h1; 


  for (unsigned int ipd = 0; ipd < pds.size(); ++ipd) {
    vector<string> pdTriggers = fHltConfig.datasetContent(pds[ipd]);
    cout << "  --> pd: " << pds[ipd] << endl;
    h1 = new TH1D(Form("triggers_%s_run%d", pds[ipd].c_str(), run.run()), 
		  Form("triggers_%s_run%d (%s)", pds[ipd].c_str(), run.run(), fHltConfig.tableName().c_str()), 
		  pdTriggers.size(), 0., pdTriggers.size()); 
    h1->SetDirectory(gHFFile); 
    for (unsigned int it = 0; it < pdTriggers.size(); ++it) {
      cout << "           " << pdTriggers[it] << endl;
      h1->GetXaxis()->SetBinLabel(it+1, pdTriggers[it].c_str()); 
    }
    h1->Write();
    delete h1; 
  }
  
  pDir->cd(); 


}

// ----------------------------------------------------------------------
void HFDebug::endRun(Run const&, EventSetup const&) {
  fValidHLTConfig = false;
} // HFDebug::endRun()

// ----------------------------------------------------------------------
void  HFDebug::beginJob() {

}

// ----------------------------------------------------------------------
void  HFDebug::endJob() {

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDebug);
