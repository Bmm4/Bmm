#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <Math/VectorUtil.h>
#include <TMath.h>

#include "HFDumpJets.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/JetFloatAssociation.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "PhysicsTools/JetMCUtils/interface/CandMCTag.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaJet.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/getRef.h"


// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace edm;
using namespace reco;
using namespace JetMCTagUtils;
using namespace CandMCTagUtils;




// ----------------------------------------------------------------------
HFDumpJets::HFDumpJets(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fDoFlavorTagging(iConfig.getUntrackedParameter<int>("doflavortagging", 0)),
  fJetsLabel(iConfig.getUntrackedParameter<string>("jetsLabel", string("MCJetCorJetIcone5"))),
 fJetsTagLabel(iConfig.getUntrackedParameter<string>("jetsTagLabel", string("simpleSecondaryVertexBJetTags"))),
  fGenCandidatesLabel(iConfig.getUntrackedParameter<string>("genparticlesLabel", string("genParticles"))), 
  fsourceByRefer(iConfig.getParameter<edm::InputTag>("sourceByRefer"))
 
{
  
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpJets constructor" << endl;
  cout << "--- " << fJetsLabel.c_str() << endl;
  cout << "----------------------------------------------------------------------" << endl;
  
}


// ----------------------------------------------------------------------
HFDumpJets::~HFDumpJets() {
  
}


// ----------------------------------------------------------------------
void HFDumpJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  nevt++;
  
  // -- get the collection of CaloJets 
  edm::Handle<reco::CaloJetCollection> caloJets;
  iEvent.getByLabel( fJetsLabel.c_str(), caloJets ); 
 if( !caloJets.isValid()) { if(fVerbose > 0) cout<<"****** no "<<fJetsLabel<<endl; return; }
  // btag   
edm::Handle<reco::JetTagCollection> bTagH;
iEvent.getByLabel(fJetsTagLabel.c_str(), bTagH);
 bool btaginfo=false;
 if( bTagH.isValid()) {
//const reco::JetTagCollection & bTags = *(bTagHandle.product());
if (fVerbose > 0) cout << "==>HFDumpTrackJets> bjetstag =" <<  bTagH->size() << endl;
 btaginfo=true;
 }else {
   if(fVerbose > 0) cout << " ***** not found " <<fJetsTagLabel<< endl;
 }


 if(fVerbose > 0) cout<<"  HFDumpJets====> "<< caloJets->size()<<endl;
  //Loop over CaloJets (with Jet Flavor Tagging)
  int jetIndex=0; 
  bool first = true;
  TAnaJet *pCaloJet;
  JetMatchedPartonsCollection::const_iterator j;
  for ( reco::CaloJetCollection::const_iterator cal = caloJets->begin(); cal != caloJets->end(); cal ++ ) {
    
    pCaloJet = gHFEvent->addCaloJet();
    pCaloJet->fIndex            = jetIndex;  

    pCaloJet->fPlab.SetPtEtaPhi(cal->pt(),
				cal->eta(),
				cal->phi()
				);
 
    pCaloJet->fQ                = cal->charge();
    pCaloJet->fE                = cal->energy();
    pCaloJet->fEt               = cal->et();
    pCaloJet->fM                = cal->mass();
    pCaloJet->fMt               = cal->mt();

    pCaloJet->fEMEnergy         = cal->emEnergyFraction();
    pCaloJet->fHADEnergy        = cal->energyFractionHadronic();
    pCaloJet->finvisibleEnergy  = -1;

    pCaloJet->fn60              = cal->n60();
    pCaloJet->fn90              = cal->n90();
     
    pCaloJet->fJetFlavorAlgo    = -9999;
    pCaloJet->fJetFlavorPhys    = -9999;
    pCaloJet->fJetFlavorEne     = -9999;

    pCaloJet->fD1               = -9999; //not used
    pCaloJet->fD2               = -9999; //fJetFlavorNear2
    pCaloJet->fD3               = -9999; //fJetFlavorNear3
    pCaloJet->fD4               = -9999; //fJetFlavorHeavy
    pCaloJet->fD5               = -9999; //not used
    pCaloJet->fD6               = -9999; //not used
    pCaloJet->fD7               = -9999; //not used
    pCaloJet->ftrjpvindx       = 0; //not used
    pCaloJet->fbtag         = -9999; // btag SSV output


   
    // btag info
    if(btaginfo) {
 const reco::JetTagCollection & tagColl = *(bTagH.product());
 double rmin = 0.5; // jets match cone
 for (JetTagCollection::const_iterator ijt = tagColl.begin();
        ijt != tagColl.end(); ++ijt) {
   // match with actual jet
   TVector3 jetbcand;
   jetbcand.SetPtEtaPhi(ijt->first->pt(), ijt->first->eta(), ijt->first->phi());  
  double r = (pCaloJet->fPlab).DeltaR(jetbcand);
  if (fVerbose > 2)  cout<<" jettag "<<r<<"  "<< ijt->second<<" "<<ijt->first->pt()<<" "<<cal->pt()<<endl;
  if(r<rmin) {
pCaloJet->fbtag= ijt->second;
 if (fVerbose > 1) cout<<r<<" btag "<<pCaloJet->fbtag<<endl; 
  }
 }
    }//btaginfo

    if (fDoFlavorTagging == 1) { 

      // -- get the collection for flavor matching
      edm::Handle<reco::JetMatchedPartonsCollection> theTagByRef;
      iEvent.getByLabel (fsourceByRefer , theTagByRef);
      
      // -- get the collection of GenParticles 
      edm::Handle<GenParticleCollection> genParticlesH;
      genParticlesH.clear();
      iEvent.getByLabel (fGenCandidatesLabel.c_str(), genParticlesH);
      
      //GenParticles for Jet Flavor Tagging
      std::vector<const GenParticle *> cands;
      cands.clear();
      vector<const GenParticle *>::const_iterator found = cands.begin();
      for (GenParticleCollection::const_iterator p = genParticlesH->begin(); p != genParticlesH->end(); ++p) {
	cands.push_back( & * p );
      }
      
      if (fVerbose > 0) cout << "==>HFDumpJets> nJets =" << caloJets->size() << endl;
      if (caloJets->size() != theTagByRef->size()) {
	if (fVerbose > 0) cout << "==>HFDumpJets>ERROR: Different Size of JetCollections "  << endl; 
      }

      if (first) {
	j  = theTagByRef->begin();
	first = false;
      }
	

      const Jet *aJet             = (*j).first.get();
      const MatchedPartons aMatch = (*j).second; 
      
      if (cal->eta() != aJet->eta()) {
	if (fVerbose > 0) cout << "==>HFDumpJets>ERROR: Different jets in JetCollections" << endl; 
      }
    
      const GenParticleRef theHeaviest = aMatch.heaviest() ;
      if(theHeaviest.isNonnull()) {
	int index = -1;
	found = find(cands.begin(), cands.end(), theHeaviest.get());
	if (found != cands.end()) {
	  index = found - cands.begin();
	  pCaloJet->fD4 = index;
	} 
	if (fVerbose > 0) {
	  cout << "theHeaviest flav idx (p,eta,phi)= " 
	       << theHeaviest.get()->pdgId() << " " 
	       << index
	       << " (" << theHeaviest.get()->p() 
	       << ","  << theHeaviest.get()->eta() 
	       << ","  << theHeaviest.get()->phi() << ") " << endl; 
	}
      }
      const GenParticleRef theNearest2 = aMatch.nearest_status2() ;
      if(theNearest2.isNonnull()) {
	int index = -1;
	found = find(cands.begin(), cands.end(), theNearest2.get());
	if (found != cands.end()) {
	  index = found - cands.begin();
	  pCaloJet->fD2 = index;
	}
	if (fVerbose > 0) {
	  cout << "theNearest Stat2  flav idx (p,eta,phi)= " 
	       << theNearest2.get()->pdgId() << " " 
	       << index
	       << " (" << theNearest2.get()->p() 
	       << ","  << theNearest2.get()->eta() 
	       << ","  << theNearest2.get()->phi() << ") " << endl; 
	}
      }
      const GenParticleRef theNearest3 = aMatch.nearest_status3() ;
      if(theNearest3.isNonnull()) {
	int index = -1;
	found = find(cands.begin(), cands.end(), theNearest3.get());
	if (found != cands.end()) {
	  index = found - cands.begin();
	  pCaloJet->fD3 = index;
	}
	if (fVerbose > 0) {
	  cout << "theNearest Stat3  flav idx (p,eta,phi)= " 
	       << theNearest3.get()->pdgId() << " " 
	       << index
	       << " (" << theNearest3.get()->p() 
	       << ","  << theNearest3.get()->eta() 
	       << ","  << theNearest3.get()->phi() << ") " << endl; 
	}
      }
      const GenParticleRef thePhyDef = aMatch.physicsDefinitionParton() ;
      if(thePhyDef.isNonnull()) {
	int index = -1;
	found = find(cands.begin(), cands.end(), thePhyDef.get());
	if (found != cands.end()) {
	  index = found - cands.begin();
	  pCaloJet->fJetFlavorPhys = index;
	}
	if (fVerbose > 0) {
	  cout << "thePhysDefinition flav idx (p,eta,phi)= " 
	       << thePhyDef.get()->pdgId() << " " 
	       << index 
	       << " (" << thePhyDef.get()->p() 
	       << ","  << thePhyDef.get()->eta() 
	       << ","  << thePhyDef.get()->phi() << ") " << endl; 
	}
	
      }
      const GenParticleRef theAlgDef = aMatch.algoDefinitionParton() ;
      if(theAlgDef.isNonnull()) {
	int index = -1;
	found = find(cands.begin(), cands.end(), theAlgDef.get());
	if (found != cands.end()) {
	  index = found - cands.begin();
	  pCaloJet->fJetFlavorAlgo = index;
	}
	if (fVerbose > 0) {
	  cout << "theAlgoDefinition flav idx (p,eta,phi)= " 
	       << theAlgDef.get()->pdgId() << " " 
	       << index 
	       << " (" << theAlgDef.get()->p() 
	       << ","  << theAlgDef.get()->eta() 
	       << ","  << theAlgDef.get()->phi() << ") " << endl; 
	} 
	
      }
    }
    
    jetIndex++;
    j++;
 

  if (fVerbose > 0) {cout<<" calojet"<<endl; pCaloJet->dump(); }  
  } 

  return;

}

// ------------ method called once each job just before starting event loop  ------------
//void  HFDumpJets::beginJob(const edm::EventSetup& setup) {
void  HFDumpJets::beginJob() {
  nevt=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpJets::endJob() { 

  cout << "HFDumpJet>        Summary: Events processed: " << nevt << endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpJets);
