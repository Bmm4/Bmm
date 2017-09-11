// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFDumpTrigger
// ------------
// 2016/09/28 Urs Langenegger      removed TTrgObj
// 2016/07/07 Urs Langenegger      add L1 (both for stage-2 and stage-1)
// 2016/03/22 Danek Kotlinski      modify TTrgObjv2 to include all modules in the path
// 2016/01/22 Urs Langenegger      migrate to "consumes"
// 2016/01/13 Danek Kotlinski      L1 work
// 2015/10/01 Urs Langenegger      all primary datasets and triggers
// 2015/08/25 Danek Kotlinski      TTrgObjv2 and more protection for *Mu* names
// stone age  Urs Langenegger      first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpTrigger.h"

#include <iostream>
#include <bitset>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

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
#include "Bmm/RootAnalysis/rootio/TTrgObjv2.hh"



// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;
using namespace l1t;


// ----------------------------------------------------------------------
HFDumpTrigger::HFDumpTrigger(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fNevt(0),

  fHLTProcessName(iConfig.getUntrackedParameter<string>("HLTProcessName")),
  fL1MuonsLabel(iConfig.getUntrackedParameter<InputTag>("L1MuonsLabel")),

  fTriggerEventLabel(iConfig.getUntrackedParameter<InputTag>("TriggerEventLabel")),
  fTokenTriggerEvent(consumes<TriggerEvent>(fTriggerEventLabel)),

  fHLTResultsLabel(iConfig.getUntrackedParameter<InputTag>("HLTResultsLabel")),
  fTokenTriggerResults(consumes<TriggerResults>(fHLTResultsLabel)),

  fL1TriggerReadoutLabel(iConfig.getUntrackedParameter<InputTag>("gtCollection", InputTag("gtDigis"))),
  fTokenL1TriggerReadoutRecord(consumes<L1GlobalTriggerReadoutRecord>(fL1TriggerReadoutLabel)),

  fAlgInputTag(iConfig.getUntrackedParameter<InputTag>("AlgInputTag", InputTag("gtStage2Digis"))),
  fExtInputTag(iConfig.getUntrackedParameter<InputTag>("ExtInputTag", InputTag("gtStage2Digis"))),

  fAlgToken(consumes<BXVector<GlobalAlgBlk> >(fAlgInputTag)),
  fExtToken(consumes<BXVector<GlobalExtBlk> >(fExtInputTag)),

  fHltPrescaleProvider(iConfig, consumesCollector(), *this) {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTrigger constructor" << endl;
  cout << "--- Verbose                     : " << fVerbose << endl;
  cout << "--- HLT process name            : " << fHLTProcessName << endl;
  cout << "--- L1 Muons Label              : " << fL1MuonsLabel << endl;
  cout << "--- HLTResultsLabel             : " << fHLTResultsLabel << endl;
  cout << "--- Trigger Event Label         : " << fTriggerEventLabel << endl;
  cout << "--- Trigger Digi Label          : " << fL1TriggerReadoutLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fGtUtil      = new L1TGlobalUtil(iConfig, consumesCollector(), *this, fAlgInputTag, fExtInputTag);
}


// ----------------------------------------------------------------------
HFDumpTrigger::~HFDumpTrigger() {

}


// ----------------------------------------------------------------------
void HFDumpTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  fNevt++;

  // -- 2016 L1
  Handle<BXVector<GlobalAlgBlk>> alg;
  iEvent.getByToken(fAlgToken, alg);

  if (alg.isValid()) {
    Handle<BXVector<GlobalExtBlk>> ext;
    iEvent.getByToken(fExtToken, ext);

    fGtUtil->retrieveL1(iEvent, iSetup, fAlgToken);

    // grab the map for the final decisions
    const std::vector<std::pair<std::string, bool> > initialDecisions = fGtUtil->decisionsInitial();
    const std::vector<std::pair<std::string, bool> > intermDecisions = fGtUtil->decisionsInterm();
    const std::vector<std::pair<std::string, bool> > finalDecisions = fGtUtil->decisionsFinal();
    const std::vector<std::pair<std::string, int> >  prescales = fGtUtil->prescales();
    // const std::vector<std::pair<std::string, bool> > masks; // = fGtUtil->masks();
    // const std::vector<std::pair<std::string, bool> > vetoMasks; // = fGtUtil->vetoMasks();

    if (0) {
      cout << "    Bit                  Algorithm Name                  Init    aBXM  Final   PS Factor     Masked    Veto " << endl;
      cout << "============================================================================================================" << endl;
    }
    for (unsigned int i = 0; i < initialDecisions.size(); ++i) {
      if (i >= NL1T) break;
      // get the name and trigger result
      std::string name = (initialDecisions.at(i)).first;
      if (name == "NULL") continue;

      bool resultInit = (initialDecisions.at(i)).second;

      // get prescaled and final results (need some error checking here)
      bool resultInterm = (intermDecisions.at(i)).second;
      bool resultFin = (finalDecisions.at(i)).second;

      // get the prescale and mask (needs some error checking here)
      int prescale = (prescales.at(i)).second;
      // bool mask    = (masks.at(i)).second;
      // bool veto    = (vetoMasks.at(i)).second;

      if (0) cout << std::dec << setfill(' ') << "   " << setw(5) << i << "   "
		  << setw(40) << name.c_str() << "   " << setw(7) << resultInit
		  << setw(7) << resultInterm << setw(7) << resultFin << setw(10) << prescale
		  // << setw(11) << mask << setw(9) << veto
		  << endl;

      gHFEvent->fL1TNames[i]    = name.c_str();
      gHFEvent->fL1TPrescale[i] = prescale;
      gHFEvent->fL1TResult[i]   = resultFin;
      // gHFEvent->fL1TMask[i]     = mask;
      // gHFEvent->fL1TError[i]    = veto;

    }
    bool finOR = fGtUtil->getFinalOR();
    if (0) {
      cout << "--> FinalOR = " << finOR << endl;
      cout << "===========================================================================================================" << endl;
    }
    gHFEvent->fL1TDecision =  finOR;
  } else {
    // -- legacy L1
    edm::ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
    const L1GtTriggerMenu* menu = menuRcd.product();

    edm::Handle<L1GlobalTriggerReadoutRecord> hL1Readout;
    iEvent.getByToken(fTokenL1TriggerReadoutRecord, hL1Readout);

    if (hL1Readout.isValid()) {
      DecisionWord const& dWord(hL1Readout->decisionWord());
      if (dWord.size() > 0) {
	bool finOR(false);
	if (0) cout << "----------------------------------------------------------------------" << endl;
	int cnt(0);
	for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
	  if (cnt >= NL1T) break;
	  bool result = menu->gtAlgorithmResult((algo->second).algoName(), dWord);
	  if (result) finOR = true;
	  if (0) cout << cnt << " L1Menu Name: " << (algo->second).algoName() << " Alias: " << (algo->second).algoAlias()
		      << " result: " << result
		      << endl;

	  gHFEvent->fL1TNames[cnt]    = (algo->second).algoName().c_str();
	  gHFEvent->fL1TPrescale[cnt] = 1;
	  gHFEvent->fL1TResult[cnt]   = result;
	  gHFEvent->fL1TMask[cnt]     = 1;
	  gHFEvent->fL1TError[cnt]    = 0;
	  ++cnt;
	}
	if (0) cout << "--> FinalOR = " << finOR << endl;
      }
    }
  }

  // ----------------------------------------------------------------------
  // -- HLT results
  // ----------------------------------------------------------------------

  if (fVerbose > 99) cout << "Resetting all trigger arrays" << endl;
  for (int i = 0; i < NHLT; ++i) {
    gHFEvent->fHLTPrescale[i] = gHFEvent->fHLTResult[i] = gHFEvent->fHLTWasRun[i] = gHFEvent->fHLTError[i] = 0;
  }

  vector<string> validTriggerNames;
  if (fValidHLTConfig) validTriggerNames = fHltConfig.triggerNames();
  else cerr << "==> HFDumpTrigger: No valid Trigger configuration!!!" << endl;
  //can assert?!  fHltConfig.dump("PrescaleTable");

  if (validTriggerNames.size() < 1) {
    cout << "==>HFDumpTrigger: NO valid trigger names returned by HLT config provided!!??" << endl;
    return;
  }

  Handle<TriggerResults> hHLTresults;
  bool hltF = true;
  bool muonTrigger = false, muonObjects=false, l3muon=false;
  int countSelectedMuonObjects=0, countLabels=0;
  int lastModuleIndex=-1;
  string lastModuleLabel="";
  string lastModuleType="";
  const int TEMP_SIZE = 100;
  int lastMuonID[TEMP_SIZE], lastMuonIndex[TEMP_SIZE];
  float lastMuonPt[TEMP_SIZE], lastMuonEta[TEMP_SIZE], lastMuonPhi[TEMP_SIZE],lastMuonE[TEMP_SIZE];
  for(int i=0; i<TEMP_SIZE;++i) lastMuonIndex[i]=-1;

  try {
    iEvent.getByToken(fTokenTriggerResults, hHLTresults);
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==>HFDumpTrigger> Triggerresults  " << fHLTResultsLabel.encode() << " not found " << endl;
    hltF = false;
    return;
  }

  Handle<trigger::TriggerEvent> triggerSummary;
  hltF = true;
  try {
    iEvent.getByToken(fTokenTriggerEvent, triggerSummary);
  } catch (const cms::Exception& e) {
    hltF = false;
    cout<<"Error!! No TriggerEvent with label " << fTriggerEventLabel << endl;
    return;
  }

  if (hltF) {
    const TriggerNames &trigName = iEvent.triggerNames(*hHLTresults);
    gHFEvent->fHLTDecision = hHLTresults->accept();
    if (fVerbose > 1)
      cout << "hHLTresults->size() = " << hHLTresults->size() << " and HLT accept = " << gHFEvent->fHLTDecision << endl;

    unsigned int index(999);
    bool wasrun(false), result(false), error(false);
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
      error    = hHLTresults->error(index);
      if (psSet > -1) {
	prescale = fHltConfig.prescaleValue(psSet, validTriggerNames[it]);
      } else {
	//	cout << "==>HFDumpTrigger> error in prescale set!?" << endl;
	prescale = 0;
      }

      gHFEvent->fHLTNames[index]    = validTriggerNamesT;
      gHFEvent->fHLTResult[index]   = result;
      gHFEvent->fHLTWasRun[index]   = wasrun;
      gHFEvent->fHLTError[index]    = error;
      gHFEvent->fHLTPrescale[index] = prescale;

      const vector<string>& moduleLabels(fHltConfig.moduleLabels(index));
      const unsigned int moduleIndex(hHLTresults->index(index));

      if ( (fVerbose > 99) ) {
	cout << " HLTName = " << gHFEvent->fHLTNames[index]
	     << " Index = " << index<<"/"<<triggerIndex
	     << " Result = " <<  gHFEvent->fHLTResult[index]
	     << " WasRun = " <<  gHFEvent->fHLTWasRun[index]
	     << " Error = " <<  gHFEvent->fHLTError[index]
	     << " Prescale = " <<  gHFEvent->fHLTPrescale[index]
	     << endl;
      }

      if(result) {  // Do only for passed HLT
	if (fVerbose > 2)
	  cout<<" passed "<<validTriggerNames[it]<<" "
	      <<moduleLabels.size()<<" "
	      <<moduleIndex<<" Index "<<index<<" "<<it<<" "<<triggerIndex<<endl;

	if(1) {  // store all, does not cost much space

	  bool isMuon =
	    (validTriggerNamesT.Contains("mu")) ||
	    (validTriggerNamesT.Contains("Mu")) ||
	    (validTriggerNamesT.Contains("MU"));  // select only muon triggers
	  muonTrigger = muonTrigger
	    || (isMuon && // accumulate if several passed HLTs
		(!(validTriggerNamesT.Contains("Multiplicity")) && // ignore
		 !(validTriggerNamesT.Contains("L1simulation_step")) && //ignore
		 !(validTriggerNamesT.Contains("noMu")) && // ignore those
		 !(validTriggerNamesT.Contains("AlCa")) ));

	  lastModuleIndex=-1;
	  lastModuleLabel="";
	  lastModuleType="";
	  int numMuons = 0;

	  // loop over all modules
	  for(unsigned int j=0;j<=moduleIndex;j++) {
	    const string& moduleLabel = moduleLabels[j];
	    const string& type = fHltConfig.moduleType(moduleLabel);
	    const unsigned int filterIndex = triggerSummary->filterIndex(InputTag(moduleLabel,"","HLT"));
	    if(fVerbose>98) {
	      cout<<j<<" label "<<moduleLabel<<" type "<<type<<" index "
		  <<filterIndex<<endl;
	    }

	    if(filterIndex < triggerSummary->sizeFilters()) {

	      const trigger::Vids& VIDS (triggerSummary->filterIds(filterIndex));
	      const trigger::Keys& KEYS(triggerSummary->filterKeys(filterIndex));
	      const size_type nI(VIDS.size());
	      const size_type nK(KEYS.size());
	      assert(nI==nK);
	      const size_type n(max(nI,nK));

	      TString moduleLabelT = TString(moduleLabel);
	      //bool muLabel = (moduleLabelT.Contains("mu") ||
	      //	      moduleLabelT.Contains("Mu")) ||
	      //!moduleLabelT.Contains("MuonNo");
	      //bool muLabel = true;  // override to save all modules, is it a big increae?
	      //if( (n>0) && muLabel) {
	      if(n>0) { // save all modules, not a big increae.
		countLabels++;
		muonObjects = true;
		lastModuleIndex = filterIndex;
		lastModuleLabel = moduleLabel;
		lastModuleType  = type;

		// loop over objects
		numMuons=0;
		const trigger::TriggerObjectCollection& TOC = triggerSummary->getObjects();
		for (size_type i=0; i!=n; ++i) {
		  const trigger::TriggerObject& TO=TOC[KEYS[i]];

		  if(fVerbose>11)
		    cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
			 << TO.id() << " " << TO.pt() << " " << TO.eta() << " "
			 << TO.phi() << " " << TO.mass()<< endl;

		  if(i<TEMP_SIZE) {
		    numMuons++;
		    lastMuonIndex[i] = KEYS[i];
		    lastMuonID[i] = TO.id();
		    lastMuonPt[i] = TO.pt();
		    lastMuonEta[i] = TO.eta();
		    lastMuonPhi[i] = TO.phi();
		    lastMuonE[i] = TO.energy();
		  } else {
		    // cout<<" more than "<<TEMP_SIZE<<" muons in trigger object? size too small "
		    // 	<<i<<" "<<n<<endl;
		    // cout<<validTriggerNamesT<<" "<<index<<" "<<moduleLabelT<<" "<<type<<" "
		    // 	<<filterIndex<<" "<<n<<" "<<KEYS[i]<<" "<<TO.pt()<<" "<<TO.eta()<<endl;
		  }

		} // end for(i)

		// save every module
		if(numMuons>0) { // do only for modules with at least 1 kinemtics
		  TTrgObjv2 *pTO = gHFEvent->addTrgObjv2();
		  pTO->fHltPath  = validTriggerNames[it];
		  pTO->fHltIndex  = index;
		  pTO->fLabel  = lastModuleLabel;
		  pTO->fType  = lastModuleType;
		  pTO->fNumber  = lastModuleIndex;
		  TLorentzVector v;
		  for(int m=0;m<numMuons;++m) {
		    v.SetPtEtaPhiE(lastMuonPt[m],lastMuonEta[m],lastMuonPhi[m],lastMuonE[m]);
		    pTO->fP.push_back(v);
		    pTO->fID.push_back(lastMuonID[m]);
		    pTO->fIndex.push_back(lastMuonIndex[m]);
		    lastMuonIndex[m]=-1;
		  }  //  for(m)
		} // numMuons>0
	      } // end if n>0
	    }  // end if
	  } // for j, modul loop
	} // if mu
      } // if hlt passed

      // Do we need this?
      if ( moduleIndex < moduleLabels.size() &&
	   fHltConfig.moduleType(moduleLabels[moduleIndex]) == "HLTPrescaler" ){
	if (fVerbose > 99) cout << " HLTPrescaler  " << endl;
	int tmp = gHFEvent->fHLTError[index];
	gHFEvent->fHLTError[index] = (tmp<<2);
	gHFEvent->fHLTError[index] |= 1;
      }

      if ( (gHFEvent->fHLTError[index] & 1)  && (fVerbose > 99) )
	cout << " Last active module type =  " << fHltConfig.moduleType(moduleLabels[moduleIndex]) << endl;


    } // for it
  } //  if hltf


  // ----------------------------------------------------------------------
  // -- Get trigger muon  objects
  // ----------------------------------------------------------------------

  if (hltF) {
    // -- HLT L1 muon candidates FIXME: Fills NOTHING
    //    string L1NameCollection("hltL1extraParticles");
    string L1NameCollection("l1extraParticles");
    trigger::size_type Index(0);
    //    Index = triggerSummary->collectionIndex(edm::InputTag(L1NameCollection, "", fTriggerEventLabel.process()));
    Index = triggerSummary->collectionIndex(edm::InputTag(L1NameCollection, "", "RECO"));

    if (Index < triggerSummary->sizeCollections()) {
      TString label = TString(L1NameCollection.c_str());
      const trigger::Keys& Keys(triggerSummary->collectionKeys());
      const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
      const trigger::size_type n1 (Keys.at(Index));
      for (trigger::size_type i = n0; i != n1; ++i) {
	const trigger::TriggerObject& obj( triggerSummary->getObjects().at(i) );
	if(fVerbose>10)
	  cout << "===> For L1: " << i << " id: " << obj.id()
	       << " m = " << obj.mass() << " pT,eta,phi = " << obj.pt() << ","
	       <<  obj.eta() << "," << obj.phi() << endl;

	TTrgObjv2 *pTO = gHFEvent->addTrgObjv2();
	TLorentzVector v;
	v.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
	pTO->fP.push_back(v);
	pTO->fID.push_back(obj.id());
	pTO->fP.push_back(v);
	pTO->fLabel = label;
	pTO->fType = "l1muon";
      }
      if(fVerbose>9)
	cout << "===> Found L1 trigger collection -> " << L1NameCollection << " " << (n1-n0)<<endl;
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
	  cout << "===> For L2: " << i << " id: " << obj.id() << " m = "
	       << obj.mass() << " pT,eta,phi = " << obj.pt() << ","
	       <<  obj.eta() << "," << obj.phi() << endl;

	TTrgObjv2 *pTO = gHFEvent->addTrgObjv2();
	TLorentzVector v;
	v.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
	pTO->fP.push_back(v);
	pTO->fID.push_back(obj.id());
	pTO->fLabel  = label;
	pTO->fType = "l2muon";
      }
      if(fVerbose>9) cout << "===> Found L2 trigger collection -> " << L2NameCollection << " "
			  << (n1-n0)<< endl;
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
	l3muon=true;
	if(fVerbose>10)
	  cout << "===> For L3: " << i << " id: " << obj.id() << " m = "
	       << obj.mass() << " pT,eta,phi = " << obj.pt() << ","
	       <<  obj.eta() << "," << obj.phi() << endl;

	TTrgObjv2 *pTO = gHFEvent->addTrgObjv2();
	TLorentzVector v;
	v.SetPtEtaPhiE(obj.pt(), obj.eta(), obj.phi(), obj.energy());
	pTO->fP.push_back(v);
	pTO->fID.push_back(obj.id());
	pTO->fLabel  = label;
	pTO->fType = "l3muon";
      }
      if(fVerbose>9) cout << "===> Found L3 trigger collection -> " << L3NameCollection << " "
			  << (n1-n0)<< endl;
    }
  }


  if (fVerbose > 0)  {
    cout << "for event " << fNevt
	 << " muon trigger " << muonTrigger
	 << " muon objects " << muonObjects << " labels" << countLabels << " selected "
	 << countSelectedMuonObjects << " l3muon " << l3muon
	 << endl;

    if(fVerbose >1) {
      cout << " stored objects " << gHFEvent->nTrgObjv2() << endl;
      for(int i=0; i<gHFEvent->nTrgObjv2();++i) {
	TTrgObjv2 *pTO = gHFEvent->getTrgObjv2(i);
	pTO->dump();

	cout<<i<<" "<<pTO->fHltPath<<" "<<pTO->fHltIndex<<" "
	    <<pTO->fLabel<<" "<<pTO->fType<<" "<<pTO->fNumber<<endl;

	vector<int> muonIndex = pTO->fIndex;
	vector<int> muonID = pTO->fID;
	vector<TLorentzVector> muonP = pTO->fP;
	int num = muonIndex.size();
	for(int n=0;n<num;++n) {
	  int index = muonIndex[n];
	  int id = muonID[n];
	  TLorentzVector p = muonP[n];
	  cout<<n<<" "<<index<<" "<<id<<" "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;
	}
      }
    }
  }


}

// ----------------------------------------------------------------------
void  HFDumpTrigger::beginRun(const Run &run, const EventSetup &iSetup) {
  bool hasChanged;
  fValidHLTConfig = fHltConfig.init(run,iSetup,fHLTProcessName,hasChanged);


  // -- Initialize fHltConfig_
  bool changed = true;
  if (fHltPrescaleProvider.init(run, iSetup, fHLTResultsLabel.process(), changed)){
    fHltConfig =  fHltPrescaleProvider.hltConfigProvider();
    resetRun(changed);
  } else {
    cout << "HFDumpTrigger::beginRun> " << "HLTPrescaleProvider initialization failed!" << endl;
    return;
  }

  cout << "HFDumpMuons::beginRun: " << run.run() << endl;
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
void HFDumpTrigger::endRun(Run const&, EventSetup const&) {
  fValidHLTConfig = false;
} // HFDumpTrigger::endRun()

// ----------------------------------------------------------------------
void  HFDumpTrigger::beginJob() {

}

// ----------------------------------------------------------------------
void  HFDumpTrigger::endJob() {

}

// ----------------------------------------------------------------------
void HFDumpTrigger::resetRun(bool changed) {

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTrigger);
