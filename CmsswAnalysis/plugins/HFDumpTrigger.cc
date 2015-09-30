#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpTrigger.h"

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
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
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

#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"


// -- Yikes!
extern TAna01Event *gHFEvent;
extern TFile       *gHFFile;

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;


// ----------------------------------------------------------------------
HFDumpTrigger::HFDumpTrigger(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fHLTProcessName(iConfig.getUntrackedParameter<string>("HLTProcessName")),
  fL1GTReadoutRecordLabel(iConfig.getUntrackedParameter<InputTag>("L1GTReadoutRecordLabel", edm::InputTag("gtDigis"))),
  fL1GTmapLabel(iConfig.getUntrackedParameter<InputTag>("hltL1GtObjectMap")),
  fL1MuonsLabel(iConfig.getUntrackedParameter<InputTag>("L1MuonsLabel")),
  fTriggerEventLabel(iConfig.getUntrackedParameter<InputTag>("TriggerEventLabel")),
  fHLTResultsLabel(iConfig.getUntrackedParameter<InputTag>("HLTResultsLabel"))
{

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTrigger constructor" << endl;
  cout << "--- Verbose                     : " << fVerbose << endl;
  cout << "--- HLT process name            : " << fHLTProcessName << endl;
  cout << "--- L1 GT Readout Record Label  : " << fL1GTReadoutRecordLabel << endl;
  cout << "--- L1 GT Object Map Label      : " << fL1GTmapLabel << endl;
  cout << "--- L1 Muons Label              : " << fL1MuonsLabel << endl;
  cout << "--- HLTResultsLabel             : " << fHLTResultsLabel << endl;
  cout << "--- Trigger Event Label         : " << fTriggerEventLabel << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fNevt = 0; 
  
}


// ----------------------------------------------------------------------
HFDumpTrigger::~HFDumpTrigger() {
  
}


// ----------------------------------------------------------------------
void HFDumpTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  fNevt++;

  // TESTING
  //if(fNevt<100) return;
  //if( fNevt!=169 && fNevt!=511 && fNevt!=940  ) return;
  // This will not be needed once we switch to TrgObjv2
  // class muonTrigObject {
  // public:
  //   int hltIndex;
  //   string hltPath;
  //   int lastModuleIndex;
  //   int lastModuleLevel;
  //   string lastModuleLabel;
  //   string lastModuleType;
  //   vector<int> muonIndex;
  //   vector<int> muonID;
  //   vector<float> muonPt;
  //   vector<float> muonEta;
  //   vector<float> muonPhi;
  // };
  // vector<muonTrigObject> muonTrigObjects;
  // vector<muonTrigObject>::iterator  imto;
  // muonTrigObject oneMuonTrigObject;


  // ----------------------------------------------------------------------
  // -- L1 results: physics and technical triggers
  // ----------------------------------------------------------------------
  if (fVerbose > 99) cout << "Resetting all trigger arrays" << endl;
  for (int i = 0; i < NL1T; ++i) {
    gHFEvent->fL1TPrescale[i] = gHFEvent->fL1TResult[i] = gHFEvent->fL1TMask[i] = gHFEvent->fL1TError[i] = 0; 
  }

  for (int i = 0; i < NLTT; ++i) {
    gHFEvent->fLTTPrescale[i] = gHFEvent->fLTTResult[i] = gHFEvent->fLTTMask[i] = gHFEvent->fLTTError[i] = 0; 
  }

  for (int i = 0; i < NHLT; ++i) {
    gHFEvent->fHLTPrescale[i] = gHFEvent->fHLTResult[i] = gHFEvent->fHLTWasRun[i] = gHFEvent->fHLTError[i] = 0; 
  }
    
  if (fVerbose > 5) cout << "Retrieving trigger records" << endl;
  Handle<L1GlobalTriggerReadoutRecord> L1GTRR;
  iEvent.getByLabel(fL1GTReadoutRecordLabel,L1GTRR);
  Handle<L1GlobalTriggerObjectMapRecord> hL1GTmap; 
  iEvent.getByLabel("hltL1GtObjectMap", hL1GTmap);

  if (fVerbose > 5) cout << "Retrieving L1GtUtils" << endl;
  L1GtUtils l1GtUtils;
  l1GtUtils.retrieveL1EventSetup(iSetup);
  // cout << "L1 trigger menu: ";
  // cout << l1GtUtils.l1TriggerMenu() << endl;

  if (fVerbose > 5) cout << "Get L1GtTriggerMenu" << endl;
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();

  string algoname; 
  int    algobit(-1); 
  bool   result(false); 
  bool   resultBeforeMask(false); // not really used, needed by interface which is by ref
  int    prescale(0); 
  int    mask(0); 
  int    iErrorCode(0); 

  for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
    algoname = (algo->second).algoName();
    algobit  = (algo->second).algoBitNumber();
    //result   = l1GtUtils.decisionAfterMask(iEvent, algoname, iErrorCode);
    //mask     = l1GtUtils.triggerMask(iEvent, algoname, iErrorCode);
    //prescale = l1GtUtils.prescaleFactor(iEvent, algoname, iErrorCode);
    // this does the same in one go, moreover the three calls above use this - three times of course...
    iErrorCode = l1GtUtils.l1Results(iEvent, algoname, resultBeforeMask, result, prescale, mask);

    gHFEvent->fL1TNames[algobit]    = TString(algoname);
    gHFEvent->fL1TResult[algobit]   = result;
    gHFEvent->fL1TMask[algobit]     = mask;
    gHFEvent->fL1TError[algobit]    = iErrorCode;
    gHFEvent->fL1TPrescale[algobit] = prescale;
  }

  for (CItAlgo algo = menu->gtTechnicalTriggerMap().begin(); algo != menu->gtTechnicalTriggerMap().end(); ++algo) {
    algoname = (algo->second).algoName();
    algobit  = (algo->second).algoBitNumber();
    result   = l1GtUtils.decisionAfterMask(iEvent, algoname, iErrorCode);
    mask     = l1GtUtils.triggerMask(iEvent, algoname, iErrorCode);
    prescale = l1GtUtils.prescaleFactor(iEvent, algoname, iErrorCode);
    
    gHFEvent->fLTTNames[algobit]    = TString(algoname);
    gHFEvent->fLTTResult[algobit]   = result;
    gHFEvent->fLTTMask[algobit]     = mask;
    gHFEvent->fLTTError[algobit]    = iErrorCode;
    gHFEvent->fLTTPrescale[algobit] = prescale;
  }

  // ----------------------------------------------------------------------
  // -- HLT results
  // ----------------------------------------------------------------------

  vector<string> validTriggerNames;
  if (validHLTConfig) validTriggerNames = hltConfig.triggerNames();
  else cerr << "==> HFDumpTrigger: No valid Trigger configuration!!!" << endl;
  //can assert?!  hltConfig.dump("PrescaleTable");

  if (validTriggerNames.size() < 1) {
    cout << "==>HFDumpTrigger: NO valid trigger names returned by HLT config provided!!??" << endl;
    return;
  }

  Handle<TriggerResults> hHLTresults;
  bool hltF = true;
  int selected=0;  // to store the right trigger objects
  string selectedObj[10];
  int selectedObjIndex[10];
  bool muonTrigger = false, muonLabels = false, muonObjects=false, l3muon=false;  
  int countSelectedMuonObjects=0, countMuonLabels=0;
  //int countL1Modules=0, countL2Modules=0, countL3Modules=0, countLNoModules=0;
  int lastModuleIndex=-1;
  //int lastModuleLevel=-1;
  string lastModuleLabel="";
  string lastModuleType="";
  const int TEMP_SIZE = 100;
  int lastMuonID[TEMP_SIZE], lastMuonIndex[TEMP_SIZE];
  float lastMuonPt[TEMP_SIZE], lastMuonEta[TEMP_SIZE], lastMuonPhi[TEMP_SIZE],lastMuonE[TEMP_SIZE];
  for(int i=0; i<TEMP_SIZE;++i) lastMuonIndex[i]=-1;

  try {
    iEvent.getByLabel(fHLTResultsLabel, hHLTresults);
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==>HFDumpTrigger> Triggerresults  " << fHLTResultsLabel.encode() << " not found " << endl;
    hltF = false;
    return;
  }
  
  Handle<trigger::TriggerEvent> triggerSummary;
  hltF = true;
  try {
    iEvent.getByLabel(fTriggerEventLabel, triggerSummary);
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
    psSet = hltConfig.prescaleSet(iEvent, iSetup);
    //    cout << "validTriggerNames.size() = " << validTriggerNames.size() << endl;

    // Loop over all HLT-paths
    for (unsigned int it = 0; it < validTriggerNames.size(); ++it) {
      index    = trigName.triggerIndex(validTriggerNames[it]); 
      const unsigned int triggerIndex = hltConfig.triggerIndex(validTriggerNames[it]); //dk
      if(index!=triggerIndex) cout<<" something wrong 1 "<<index<<" "<<triggerIndex<<endl;

      if (index >= hHLTresults->size())
	continue;
		
      TString validTriggerNamesT = TString(validTriggerNames[it]);
      result   = hHLTresults->accept(index);
      wasrun   = hHLTresults->wasrun(index);
      error    = hHLTresults->error(index);
      if (psSet > -1) {
	prescale = hltConfig.prescaleValue(psSet, validTriggerNames[it]);
      } else {
	//	cout << "==>HFDumpTrigger> error in prescale set!?" << endl;
	prescale = 0;
      }
      
      gHFEvent->fHLTNames[index]    = validTriggerNamesT;
      gHFEvent->fHLTResult[index]   = result;
      gHFEvent->fHLTWasRun[index]   = wasrun;
      gHFEvent->fHLTError[index]    = error;
      gHFEvent->fHLTPrescale[index] = prescale;

      const vector<string>& moduleLabels(hltConfig.moduleLabels(index));
      const unsigned int moduleIndex(hHLTresults->index(index));

      if ( (fVerbose > 99) ) {
	cout << " HLTName = " << gHFEvent->fHLTNames[index]
	     << " Index = " << index<<"/"<<triggerIndex
	     << " Result = " <<  gHFEvent->fHLTResult[index]
	     << " WasRun = " <<  gHFEvent->fHLTWasRun[index]
	     << " Error = " <<  gHFEvent->fHLTError[index]
	     << " Presacle = " <<  gHFEvent->fHLTPrescale[index] 
	     << endl;
      }

      // Find L3 Muon modules 
      // Do only for passed HLT
      if(result) {

	if (fVerbose > 2) 
	  cout<<" passed "<<validTriggerNames[it]<<" "
	      <<moduleLabels.size()<<" "
	      <<moduleIndex<<" Index "<<index<<" "<<it<<" "<<triggerIndex<<endl;

	bool lookAt = 
	  (validTriggerNamesT.Contains("mu")) ||
	  (validTriggerNamesT.Contains("Mu")) ||
	  (validTriggerNamesT.Contains("MU"));  // select only muon triggers

	if(lookAt || fVerbose>98) { 

	  muonTrigger = muonTrigger 
	    || (lookAt && // accumulate if several passed HLTs
		(!(validTriggerNamesT.Contains("Multiplicity")) && // ignore
		 !(validTriggerNamesT.Contains("L1simulation_step")) && //ignore
		 !(validTriggerNamesT.Contains("noMu")) && // ignore those
		 !(validTriggerNamesT.Contains("AlCa")) ));

	  
	  lastModuleIndex=-1;
	  //lastModuleLevel=-1;
	  lastModuleLabel="";
	  lastModuleType="";
	  int numMuons = 0;

	  // loop over all modules 
	  for(unsigned int j=0;j<=moduleIndex;j++) {
	    const string& moduleLabel = moduleLabels[j]; 
	    const string& type = hltConfig.moduleType(moduleLabel);
	    const unsigned int filterIndex = triggerSummary->
	      filterIndex(InputTag(moduleLabel,"","HLT"));
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
	      bool muLabel = (moduleLabelT.Contains("mu") ||
			      moduleLabelT.Contains("Mu")) ||
		              !moduleLabelT.Contains("MuonNo");

	      if( (n>0) && muLabel) {

		countMuonLabels++;
		muonLabels = muonLabels || muLabel; // accumulate if several hlts
		if(fVerbose>11) 
		  cout<<j<<" label "<<moduleLabel<<" type "<<type<<" index "
		      <<filterIndex<<" "<< n  << " accepted TRIGGER objects found: " << endl;

		muonObjects = true;

		lastModuleIndex = filterIndex;
		lastModuleLabel = moduleLabel;
		lastModuleType  = type;

		// TESTING
		// if(moduleLabelT.Contains("L3") || moduleLabelT.Contains("l3")) {
		//   countL3Modules++;
		//   lastModuleLevel  = 3;
		//   //cout<<" L3 module "<<endl;
		// } else if(moduleLabelT.Contains("L2") || moduleLabelT.Contains("l2")) {
		//   countL2Modules++;
		//   lastModuleLevel  = 2;
		//   //cout<<" L2 module "<<endl;
		// } else if(moduleLabelT.Contains("L1") || moduleLabelT.Contains("l1")) {
		//   countL1Modules++; 
		//   lastModuleLevel  = 1;
		//   //cout<<" L1 module "<<endl;
	        // } else {
		//   countLNoModules++;
		//   lastModuleLevel  = 3;
		//   //cout<<" non L module "<<endl;
		// }

		// loop over muons
		numMuons=0;
		const trigger::TriggerObjectCollection& TOC = triggerSummary->getObjects();
		for (size_type i=0; i!=n; ++i) {
		  const trigger::TriggerObject& TO=TOC[KEYS[i]];

		  if(i<TEMP_SIZE) {
		    numMuons++;
		    lastMuonIndex[i] = KEYS[i];
		    lastMuonID[i] = TO.id();
		    lastMuonPt[i] = TO.pt();
		    lastMuonEta[i] = TO.eta();
		    lastMuonPhi[i] = TO.phi();
		    lastMuonE[i] = TO.energy();
		  } else {
		    cout<<" more than "<<TEMP_SIZE<<" muons in trigger object? size too small "
			<<i<<" "<<n<<endl;
		    cout<<validTriggerNamesT<<" "<<index<<" "<<moduleLabelT<<" "<<type<<" "
			<<filterIndex<<" "<<n<<" "<<KEYS[i]<<" "<<TO.pt()<<" "<<TO.eta()<<endl;
		  }

		  if(fVerbose>11)  
		    cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
			 << TO.id() << " " << TO.pt() << " " << TO.eta() << " " 
			 << TO.phi() << " " << TO.mass()<< endl;                                                                 
		} // end for
		
	      } // end if n>0
	    }  // end if

	    // do we still need this? maybe keep it for the trigger confirmation 
	    if(type == "HLTMuonDimuonL3Filter") {
	      if(selected<10) {
		selectedObj[selected] = moduleLabels[j];
		selectedObjIndex[selected] = j+1;
		selected++;
	      } else {
		cout<<"ERROR: array too small "<<selected<<endl;
	      }
	    } // if type


	  } // for j, modul loop

	  // testing
	  // oneMuonTrigObject.hltIndex = index;
	  // oneMuonTrigObject.hltPath = validTriggerNames[it];
	  // oneMuonTrigObject.muonIndex.clear();
	  // oneMuonTrigObject.muonID.clear();
	  // oneMuonTrigObject.muonPt.clear();
	  // oneMuonTrigObject.muonEta.clear();
	  // oneMuonTrigObject.muonPhi.clear();
	  // //oneMuonTrigObject.muonE.clear();
	  // oneMuonTrigObject.lastModuleIndex = lastModuleIndex;
	  // oneMuonTrigObject.lastModuleLabel = lastModuleLabel;
	  // oneMuonTrigObject.lastModuleType  = lastModuleType;
	  // oneMuonTrigObject.lastModuleLevel = lastModuleLevel;
	  // for(int m=0;m<numMuons;++m) {	    
	  //   oneMuonTrigObject.muonID.push_back(lastMuonID[m]);
	  //   oneMuonTrigObject.muonIndex.push_back(lastMuonIndex[m]);
	  //   oneMuonTrigObject.muonPt.push_back(lastMuonPt[m]);
	  //   oneMuonTrigObject.muonEta.push_back(lastMuonEta[m]);
	  //   oneMuonTrigObject.muonPhi.push_back(lastMuonPhi[m]);
	  //   //oneMuonTrigObject.muonE.push_back(lastMuonE[m]);
	  // }
	  // muonTrigObjects.push_back(oneMuonTrigObject);
	  //cout<<" save "<<oneMuonTrigObject.hltPath<<" "<<oneMuonTrigObject.hltIndex<<" "
	  //  <<lastModuleIndex<<" "<<lastModuleLabel<<" "
	  //  <<lastModuleType<<" "<<lastModuleLevel<<" "<<numMuons<<endl;

	  // load the trigger information to TrgObjv2
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
	  }


	} // if mu
      } // if hlt passed

      //  cout << " Last active module - label/type: "
      //   << moduleLabels[moduleIndex] << "/" << hltConfig.moduleType(moduleLabels[moduleIndex])
      //   << endl;

      if ( moduleIndex < moduleLabels.size() && hltConfig.moduleType(moduleLabels[moduleIndex]) == "HLTPrescaler" ){
	if (fVerbose > 99) cout << " HLTPrescaler  " << endl;
	int tmp = gHFEvent->fHLTError[index];
	gHFEvent->fHLTError[index] = (tmp<<2);
	gHFEvent->fHLTError[index] |= 1;
      }

      // cout << "gHFEvent->fHLTError[index] = " << gHFEvent->fHLTError[index] << endl;
      if ( (gHFEvent->fHLTError[index] & 1)  && (fVerbose > 99) )
	cout << " Last active module type =  " << hltConfig.moduleType(moduleLabels[moduleIndex]) << endl;


    } // for it 
  } //  if hltf
  

  // Testing only 
  if (fVerbose > 99)  {
    cout<<"Selected HLT modules "<<selected<<endl;
    for(int i=0; i<selected; ++i) {
      cout<<i<<" "<<selectedObj[i]<<" "<<selectedObjIndex[i]<<endl;
    }
    if(selected>1) cout<<" MORE THAN ONE OBJECT SELECTED "<<endl;
  } // if verbose 

  // ----------------------------------------------------------------------
  // -- Get trigger muon  objects
  // ----------------------------------------------------------------------
  
  if (hltF) {
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
	if(fVerbose>10) 
	   cout << "===> For L1: " << i << " id: " << obj.id() 
	   << " m = " << obj.mass() << " pT,eta,phi = " << obj.pt() << "," 
	   <<  obj.eta() << "," << obj.phi() << endl;
			  
	TTrgObj *pTO = gHFEvent->addTrgObj();
	pTO->fP.SetPtEtaPhiE(obj.pt(), 
			     obj.eta(), 
			     obj.phi(), 
			     obj.energy()
			     ); 
	pTO->fID     = obj.id(); 
	pTO->fLabel  = label;
			  
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
			  
	TTrgObj *pTO = gHFEvent->addTrgObj();
	pTO->fP.SetPtEtaPhiE(obj.pt(), 
			     obj.eta(), 
			     obj.phi(), 
			     obj.energy()
			     ); 
	pTO->fID     = obj.id(); 
	pTO->fLabel  = label;
			  
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
			
	TTrgObj *pTO = gHFEvent->addTrgObj();
	pTO->fP.SetPtEtaPhiE(obj.pt(), 
			     obj.eta(), 
			     obj.phi(), 
			     obj.energy()
			     ); 
	pTO->fID     = obj.id(); 
	pTO->fLabel  = label;
			
      }
      if(fVerbose>9) cout << "===> Found L3 trigger collection -> " << L3NameCollection << " " 
			  << (n1-n0)<< endl;
    }

    // -- muon filter objects
    TriggerObjectCollection allObjects = triggerSummary->getObjects();
    if(fVerbose>9) cout << "===> Found muon filter objects -> " 
			<< triggerSummary->sizeFilters() << endl;
    for (int i=0; i < triggerSummary->sizeFilters(); i++){         
      //cout<<i<<endl;
      Keys keys = triggerSummary->filterKeys(i);
      //cout<<keys.size()<<endl;
      if (keys.size() > 0) {

	//cout<<triggerSummary->filterLabel(i)  // crash 
	//<<" "<<triggerSummary->filterTagEncoded(i)
	//<<" "<<triggerSummary->collectionTagEncoded(i)
	//<<" "<<triggerSummary->filterTag(i)  
	//<<" "<<triggerSummary->filterIndex(triggerSummary->filterTag(i))
	//  <<" "<<triggerSummary->collectionTag(i)
	//<<" "<<triggerSummary->collectionIndex(triggerSummary->collectionTag(i))
	//  <<endl;

	TString label = TString(triggerSummary->filterTag(i).label()+":"+triggerSummary->filterTag(i).process()+":"+triggerSummary->filterTag(i).instance()+":");
	
	if (fVerbose > 99) 
	  cout <<i<<" object "<<triggerSummary->filterTag(i).label()<<endl; 
	
	// -- the following removes cross trigger filter objects even when they have Mu in the name!
	//if (label.Contains("Jet")) continue; // diable dk.,there are mu+jet cross triggers 
	//if (label.Contains("HT")) continue;  // same
	//if (label.Contains("Tau")) continue; // same
	//if (label.Contains("EG")) continue; // unused
	//if (label.Contains("Pi0")) continue;
	//if (label.Contains("AlCa")) continue;
	//if (label.Contains("Multiplicity")) continue;
	// if (label.Contains("Mu") || label.Contains("mu")) {
	//   // fill
	// } else {
	//   continue;
	// }

	if (fVerbose > 12) 
	  cout <<i<<" object "<<triggerSummary->filterTag(i).label()<<endl; 

	// Check if this path was in the fired HLT path 
	bool matched=false;
	for(int n=0;n<selected; ++n) {
	  if ( triggerSummary->filterTag(i).label() == selectedObj[n] ) {
	    selectedObjIndex[n]=0;
	    matched=true;
	    if (fVerbose > 12) 
	      cout<<" This is the triggered object --> "
		  <<i<<" "<<triggerSummary->filterTag(i).label()<<endl;
	    break;
	  }
	}

	if(matched) {  // save only matched objects 
	  countSelectedMuonObjects++;
	  // loop over particles (muons?) in this object and save them
	  for (unsigned int j=0; j<keys.size(); j++){
	    TTrgObj *pTO = gHFEvent->addTrgObj();
	    pTO->fP.SetPtEtaPhiE(allObjects[keys[j]].pt(), 
				 allObjects[keys[j]].eta(), 
				 allObjects[keys[j]].phi(), 
				 allObjects[keys[j]].energy()
				 ); 
	    pTO->fID     = allObjects[keys[j]].id(); 
	    pTO->fLabel  = label;
	    pTO->fNumber = i;  // marked selected objects for later analysis  
	    if (fVerbose > 12) 
	      cout << " pt = " <<  allObjects[keys[j]].pt() 
		   << " eta = " << allObjects[keys[j]].eta() 
		   << " phi = " << allObjects[keys[j]].phi() 
		   << " e = " << allObjects[keys[j]].energy() 
		   << " id = " << allObjects[keys[j]].id() 
		//<< " label: " << pTO->fLabel
		   << " number:" << pTO->fNumber
		   << endl;
	  } // end for j
	} // if matched 
      } // if size>0
    } // for i
    
  }  // if hltf
  
  if (fVerbose > 0)  {
    cout<<" for event "<< fNevt
	<< " muon trigger "<<muonTrigger<<" modules "<<muonLabels
	<<" muon objects "<<muonObjects<<" labels"<<countMuonLabels<<" selected "
	<<countSelectedMuonObjects<<" l3muon "<<l3muon<<endl;
      //<<" L1 "<<countL1Modules<<" L2 "<<countL2Modules<<" L3 "<<countL3Modules<<" LNo "
      //<<countLNoModules<<" : "<<" objects "<<countMuonObjects<<" "

    if(fVerbose >1) {

      // TESTING
      // cout<<" my trigger summary "<<muonTrigObjects.size()<<endl;
      // for(imto=muonTrigObjects.begin(); imto!=muonTrigObjects.end();++imto) {
      // 	int hltIndex = imto->hltIndex;
      // 	string hltPath = imto->hltPath;
      // 	int lastModuleIndex = imto->lastModuleIndex;
      // 	int lastModuleLevel = imto->lastModuleLevel;
      // 	string lastModuleLabel = imto->lastModuleLabel;
      // 	string lastModuleType = imto->lastModuleType;
      // 	vector<int> muonIndex = imto->muonIndex;
      // 	vector<int> muonID= imto->muonID;
      // 	vector<float> muonPt = imto->muonPt;
      // 	vector<float> muonEta = imto->muonEta;
      // 	vector<float> muonPhi = imto->muonPhi;
      // 	int numOfMuons = muonIndex.size();
      // 	//if(fVerbose>9 || numOfMuons>2) {
      // 	  cout<<" hlt "<<hltPath<<" index "<<hltIndex
      // 	      <<" last module index "<<lastModuleIndex<<" label "
      // 	      <<lastModuleLabel<<" type "<<lastModuleType<<" level "
      // 	      <<lastModuleLevel<<" num of muons "<<numOfMuons
      // 	      <<endl;
      // 	  for(int i=0;i<numOfMuons;++i) {
      // 	    cout<<"muon "<<i<<" index,id,pt,eta,phi "<<muonIndex[i]<<","
      // 		<<muonID[i]<<","<<muonPt[i]<<","
      // 		<<muonEta[i]<<","<<muonPhi[i]<<endl;
      // 	  } // for 
      // 	  //} // if 
      // } // for
      //}
      
      
      
      cout<<" stored objects "<<gHFEvent->nTrgObjv2()<<" "<<gHFEvent->nTrgObj()<<endl;
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
    } // fVerbose

    if (fVerbose > 99)  {
      cout<<"Selected HLT modules "<<selected<<endl;
      for(int i=0; i<selected; ++i) {
	cout<<selectedObj[i]<<" "<<selectedObjIndex[i]<<endl;
	if(selectedObjIndex[i]!=0) cout<<" found non matched selected object "<<selectedObj[i]<<endl;
      }
    }   

  } // verbose?


  // Check trigger information consistenct 
  if(muonTrigger && 
     (muonTrigger&&muonLabels&&muonObjects) != 
     (muonTrigger||muonLabels||muonObjects) ) {
    cout<<" Inconsistent muon trigger information evt= "<<fNevt
	<< " muon trigger "<<muonTrigger<<" modules "<<muonLabels
	<<" muon objects "<<muonObjects<<" labels"<<countMuonLabels<<" selected "
	<<countSelectedMuonObjects<<" l3muon "<<l3muon<<endl;
    //<<" L1 "<<countL1Modules<<" L2 "<<countL2Modules<<" L3 "<<countL3Modules<<" LNo "
    //<<countLNoModules<<" : "<<" l3muon "<<l3muon
    //<<" objects "<<muonFilterObjects
    //<<" num of mu filter objects selected "<<countSelectedMuonObjects
    //	<<endl;

    // cout<<" my trigger summary "<<muonTrigObjects.size()<<endl; 
    // for(imto=muonTrigObjects.begin(); imto!=muonTrigObjects.end();++imto) {
    //   int hltIndex = imto->hltIndex;
    //   string hltPath = imto->hltPath;
    //   int lastModuleIndex = imto->lastModuleIndex;
    //   int lastModuleLevel = imto->lastModuleLevel;
    //   string lastModuleLabel = imto->lastModuleLabel;
    //   string lastModuleType = imto->lastModuleType;
    //   vector<int> muonIndex = imto->muonIndex;
    //   vector<int> muonID= imto->muonID;
    //   vector<float> muonPt = imto->muonPt;
    //   vector<float> muonEta = imto->muonEta;
    //   vector<float> muonPhi = imto->muonPhi;
    //   int numOfMuons = muonIndex.size();

    //   cout<<" hlt "<<hltPath<<" index "<<hltIndex
    // 	  <<" last module index "<<lastModuleIndex<<" label "
    // 	  <<lastModuleLabel<<" type "<<lastModuleType<<" level "
    // 	  <<lastModuleLevel<<" num of muons "<<numOfMuons
    // 	  <<endl;
    //   for(int i=0;i<numOfMuons;++i) {
    // 	cout<<"muon "<<i<<" index,id,pt,eta,phi "<<muonIndex[i]<<","
    // 	    <<muonID[i]<<","<<muonPt[i]<<","
    // 	    <<muonEta[i]<<","<<muonPhi[i]<<endl;
    //   } // for 
    // } // for

  } // end if
  
} // the end

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpTrigger::beginRun(const Run &run, const EventSetup &iSetup)
{
  bool hasChanged;
  validHLTConfig = hltConfig.init(run,iSetup,fHLTProcessName,hasChanged);
  //  cout << "HFDumpTrigger::beginRun> hltConfig.tableName() = " << hltConfig.tableName() << endl;
  vector<string> v = hltConfig.datasetNames();
  
  TDirectory *pDir = gDirectory; 
  gHFFile->cd();
  TH1D *h1 = new TH1D(Form("pd%d", run.run()), hltConfig.tableName().c_str(), v.size()+1, 0., v.size()+1); 
  h1->SetDirectory(gHFFile); 

  for (unsigned int i = 0; i < v.size(); ++i) {
    //    cout << "                           " << v[i] << endl;
    h1->GetXaxis()->SetBinLabel(i+1, v[i].c_str()); 
  }

  h1->Write();
  delete h1; 

  pDir->cd(); 
}

void HFDumpTrigger::endRun(Run const&, EventSetup const&)
{
  validHLTConfig = false;
} // HFDumpTrigger::endRun()

// ------------ method called once each job just before starting event loop  ------------
void  HFDumpTrigger::beginJob() {

}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpTrigger::endJob() {

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTrigger);
