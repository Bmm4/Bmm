#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <Math/VectorUtil.h>
#include <TMath.h>

#include "HFDumpTrackJets.h"

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/BasicJet.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"

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





// ----------------------------------------------------------------------
HFDumpTrackJets::HFDumpTrackJets(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fDoFlavorTagging(iConfig.getUntrackedParameter<int>("doflavortagging", 0)),
  fJetsLabel(iConfig.getUntrackedParameter<string>("jetsLabel", string("ic5TrackJets"))),
  fJetsTagLabel(iConfig.getUntrackedParameter<string>("jetsTagLabel", string("simpleSecondaryVertexBJetTags"))),
  fTracksLabel(iConfig.getUntrackedParameter<string>("tracksLabel", string("trackCandidates"))),
  fGenCandidatesLabel(iConfig.getUntrackedParameter<string>("genparticlesLabel", string("genParticles"))),
  fsourceByRefer(iConfig.getParameter<edm::InputTag>("sourceByRefer"))


{

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTrackJets constructor" << endl;
  cout << "--- " << fJetsLabel.c_str() << endl;
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDumpTrackJets::~HFDumpTrackJets() {

}


// ----------------------------------------------------------------------
void HFDumpTrackJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if(fVerbose > 0) cout<<" ===HFDumpTrackJets "<<endl;
  nevt++;

  // handle to 0.5 cone ctf track jets
  edm::Handle<BasicJetCollection> jetsH;
  iEvent.getByLabel(fJetsLabel.c_str(),jetsH);
  const BasicJetCollection *jets   = jetsH.product();
  if( !jetsH.isValid()) {
    if (fVerbose > 0) cout << "**** HFDumpTrackJets>  no " << fJetsLabel << endl;
    return;
  }

  // get btag
  // Get b tag information
  edm::Handle<reco::JetTagCollection> bTagH;
  iEvent.getByLabel(fJetsTagLabel.c_str(), bTagH);
  bool btaginfo=false;
  if( bTagH.isValid()) {
    if (fVerbose > 0) cout << "==>HFDumpTrackJets> bjetstag =" <<  bTagH->size() << endl;
    btaginfo=true;
  }else{ if(fVerbose > 0) cout << "**** not found " <<fJetsTagLabel.c_str()<< endl;
  }

  if (fVerbose > 0) cout << "==>HFDumpTrackJets> nTrackJets =" << jetsH->size() << endl;

  //tracks (jet constituents)
  Handle<reco::CandidateView> candidates1Handle;
  iEvent.getByLabel(fTracksLabel.c_str(), candidates1Handle);

  int jetIndex=0;
  bool first = true;
  TAnaJet *pTrackJet;

  JetMatchedPartonsCollection::const_iterator j;
  for ( BasicJetCollection::const_iterator it = jets->begin(); it != jets->end(); it ++ ) {

    //for(pat::JetCollection::const_iterator it = jetsH->begin(); it != jetsH->end(); ++it) {



    pTrackJet = gHFEvent->addTrackJet();
    pTrackJet->fIndex            = jetIndex;

    pTrackJet->fPlab.SetPtEtaPhi(it->pt(),
				 it->eta(),
				 it->phi()
				 );

    pTrackJet->fQ                = it->charge();
    pTrackJet->fE                = it->energy();
    pTrackJet->fEt               = it->et();
    pTrackJet->fM                = it->mass();
    pTrackJet->fMt               = it->mt();

    pTrackJet->fEMEnergy         = -1;//not used
    pTrackJet->fHADEnergy        = -1;//not used
    pTrackJet->finvisibleEnergy  = -1;//not used

    pTrackJet->fn60              = -9999;//not used
    pTrackJet->fn90              = -9999;//not used

    pTrackJet->fJetFlavorAlgo    = -9999;
    pTrackJet->fJetFlavorPhys    = -9999;
    pTrackJet->fJetFlavorEne     = -9999;

    pTrackJet->fD1               = -9999; //
    pTrackJet->fD2               = -9999; //fJetFlavorNear2
    pTrackJet->fD3               = -9999; //fJetFlavorNear3
    pTrackJet->fD4               = -9999; //fJetFlavorHeavy
    pTrackJet->fD5               = -9999; //not used
    pTrackJet->fD6               = -9999; //not used
    pTrackJet->fD7               = -9999; //not used

    // not working...
    //    pTrackJet->ftrjpvindx            = it->primaryVertex().index(); //TrackJet primary vertex index (O-first from HI)
    pTrackJet->fbtag         = -9999; // btag SSV output


    //    pTrackJet->fVtx.fPoint.SetX(it->primaryVertex().x());
    //    pTrackJet->fVtx.fPoint.SetY(it->primaryVertex().y());
    //    pTrackJet->fVtx.fPoint.SetZ(it->primaryVertex().z());

    if (fVerbose > 0) pTrackJet->dump();

    // btag info
    if(btaginfo) {
      const reco::JetTagCollection & tagColl = *(bTagH.product());
      double rmin = 0.1; // jets match cone
      for (JetTagCollection::const_iterator ijt = tagColl.begin();
	   ijt != tagColl.end(); ++ijt) {
	// match with actual jet
	TVector3 jetbcand;
	jetbcand.SetPtEtaPhi(ijt->first->pt(), ijt->first->eta(), ijt->first->phi());
	double r = (pTrackJet->fPlab).DeltaR(jetbcand);
	if (fVerbose > 1) cout<<" btag match "<<r<<endl;
	if(r<rmin) pTrackJet->fbtag= ijt->second;
      }
    }//btaginfo

    //jet constituents
    std::vector< const reco::Candidate * > Constituent = it->getJetConstituentsQuick();
    if (fVerbose > 0) cout<<jetIndex<<" TrackJets constituents "<<Constituent.size()<<endl;

    for (unsigned int i=0; i< Constituent.size(); i++) {

      int index  = -1;
      const reco::Candidate * consti = Constituent[i];
      if (consti) {
	for (unsigned int j = 0; j < candidates1Handle->size(); ++ j ) {
	  const Candidate &p2 = (*candidates1Handle)[j];
	  const Candidate *  p1 = &p2;
	  if (consti->pt() == p1->pt() && consti->phi() == p1->phi() && consti->eta() == p1->eta() )
	    index = j;

	}
	pTrackJet->addTrack(index);
	if (fVerbose) cout << index << " pt " << consti->pt() << " phi " << consti->phi() << " eta " << consti->eta() << endl;
      }

    }

    if (fDoFlavorTagging == 1) {


      //////////////////////////////////////////
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

      if (jetsH->size() != theTagByRef->size()) {
	if (fVerbose > 0) cout << "==>HFDumpTrackJets>ERROR: Different Size of JetCollections " << endl;
      }

      ///////////////////////////////////////////

      if (first) {
	j  = theTagByRef->begin();
	first = false;
      }

      const Jet *aJet             = (*j).first.get();
      const MatchedPartons aMatch = (*j).second;

      if (it->eta() != aJet->eta()) {
	if (fVerbose > 0) cout << "==>HFDumpTrackJets>ERROR: Different jets in JetCollections " << endl;
      }

      //jet tag
      const GenParticleRef theHeaviest = aMatch.heaviest() ;
      if(theHeaviest.isNonnull()) {
	int index = -1;
	found = find(cands.begin(), cands.end(), theHeaviest.get());
	if (found != cands.end()) {
	  index = found - cands.begin();
	  pTrackJet->fD4 = index;
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
	  pTrackJet->fD2 = index;
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
	  pTrackJet->fD3 = index;
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
	  pTrackJet->fJetFlavorPhys = index;
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
	  pTrackJet->fJetFlavorAlgo = index;
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
    if (fVerbose > 0) pTrackJet->dump();
    jetIndex++;
    j++;
  }
}

// ------------ method called once each job just before starting event loop  ------------
//void  HFDumpTrackJets::beginJob(const edm::EventSetup& setup) {
void  HFDumpTrackJets::beginJob() {
  nevt=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpTrackJets::endJob() {
  cout << "HFDumpTrackJet>   Summary: Events processed: " << nevt << endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTrackJets);
