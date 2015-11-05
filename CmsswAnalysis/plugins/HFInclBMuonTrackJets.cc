#include <iostream>

#include "HFInclBMuonTrackJets.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaMuon.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"


#include "DataFormats/JetReco/interface/BasicJetCollection.h"
#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include <SimDataFormats/Vertex/interface/SimVertex.h>
#include <SimDataFormats/Vertex/interface/SimVertexContainer.h>
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"



#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace reco;
using namespace edm;


// ----------------------------------------------------------------------
HFInclBMuonTrackJets::HFInclBMuonTrackJets(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fJetMatch(iConfig.getUntrackedParameter<double>("jetmatch", 0.5)),
  fJetEtMin(iConfig.getUntrackedParameter<double>("jetetmin", 1)),
  fMuonLabel(iConfig.getUntrackedParameter<string>("muonLabel", string("globalMuons"))),
  fJetsLabel(iConfig.getUntrackedParameter<string>("jetsLabel", string("myak5TrackJets"))),
  fTracksLabel(iConfig.getUntrackedParameter<string>("tracksLabel", string("trackCandidates"))),
  fVertexLabel(iConfig.getUntrackedParameter<string>("vertexLabel", string("offlinePrimaryVerticesWithBS"))),
  fSimVertexLabel(iConfig.getUntrackedParameter<string>("simvertexLabel", string("g4SimHits"))),
  fVertexMinNdof(iConfig.getUntrackedParameter<double>("vertexMinNdof", 4.))
{
  using namespace std;
  fAllowGlobalOnly=iConfig.getUntrackedParameter<bool>("allowGlobalOnly", false);
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFInclBMuonTrackJets constructor" << endl;
  cout << "--- muonLabel         " << fMuonLabel.c_str() << endl;
  cout << "--- jetsLabel         " << fJetsLabel.c_str() << endl;
  cout << "--- tracksLabel       " << fTracksLabel.c_str() << endl;
  cout << "--- vertexLabel       " << fVertexLabel.c_str() << endl;
  cout << "--- simvertexLabel    " << fSimVertexLabel.c_str() << endl;
  if (fAllowGlobalOnly){
    cout << "!!! global-only muons allowed,  not required to be tracker muons !!!!!" << endl;
  }else{
    cout << "--- all muons will be required to be tracker muons " << endl;
  }
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFInclBMuonTrackJets::~HFInclBMuonTrackJets() {
  
}


// ----------------------------------------------------------------------
void HFInclBMuonTrackJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  fNevt++;

  int printout(0); 
  
  if (printout) cout << "----------------------------------------------------------------------" << endl;

  //muon collection
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel( fMuonLabel.c_str(), muons);
  if (!muons.isValid()) { cout<<"****** no "<<fMuonLabel<<endl; return; }

  //primary vertex
  Handle<reco::VertexCollection> primaryVertex;
  iEvent.getByLabel(fVertexLabel.c_str(),primaryVertex); 
  if (!primaryVertex.isValid()) { cout<<"****** no "<<fVertexLabel<<endl; return; }

  // use beamspot, or 0 if not available
  math::XYZPoint bs = math::XYZPoint(0.,0.,0.);
  edm::Handle<reco::BeamSpot> beamSpotCollection;
  iEvent.getByLabel("offlineBeamSpot", beamSpotCollection);
  if (beamSpotCollection.isValid()){ bs = beamSpotCollection->position();}
  else { cout<<"****no beamspot "<<endl;}
 
  // use simvertexes
  math::XYZPoint simv0 = math::XYZPoint(0.,0.,0.);
  bool simvFound =false;
  Handle<SimVertexContainer> simVtxs;
  if (iEvent.getByLabel( fSimVertexLabel.c_str(), simVtxs) && simVtxs.isValid()){
    simvFound = (simVtxs->size() != 0); 
    if(simvFound)  simv0=(*simVtxs)[0].position(); 
  }

  // TrackJets 
//   Handle<TrackJetCollection> jetsH;
//   iEvent.getByLabel(fJetsLabel.c_str(),jetsH);
//   if (!jetsH.isValid()) { cout<<"****** HFInclBMuonTrackJets> no "<<fJetsLabel<<endl;  }

  Handle<BasicJetCollection> jetsH;
  iEvent.getByLabel(fJetsLabel.c_str(),jetsH);
  if (!jetsH.isValid()) { cout<<"****** HFInclBMuonTrackJets> no "<<fJetsLabel<<endl;  }

  //tracks (jet constituents)
  Handle<reco::CandidateView> candidates1Handle;
  iEvent.getByLabel(fTracksLabel.c_str(), candidates1Handle); 
  if (!candidates1Handle.isValid()) { cout<<"****** HFInclBMuonTrackJets> no "<<fTracksLabel<<endl;  }

  //transient track builder
  edm::ESHandle<TransientTrackBuilder> builder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",builder);  

  TAnaTrack *pTrack; 
  int index (0), matchedJetIndex(-1);
  if (fVerbose > 5) cout << "==>HFInclBMuonTrackJets> nMuons =" << muons->size() << endl;

  for (reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++ muon ) { 
    
    if (!muon::isGoodMuon(*muon, muon::AllGlobalMuons)) continue;
    if (muon->track()->pt() < 4.0) continue;
    
    if (printout) cout << "==>HFInclBMuonTrackJets> found global muon " << endl;
    if (printout) cout << "MUON pt/eta/phi = " << muon->track()->pt() << "/" <<  muon->track()->eta()  << "/" <<  muon->track()->phi() << endl;
    
    if (jetsH.isValid()) {
      const BasicJetCollection *jets = jetsH.product();  
      if (!candidates1Handle.isValid()) continue;
      if ((muon->track()).index() > candidates1Handle->size() ) continue;
      double rmin = fJetMatch; 
      BasicJet* matchedjet=0;
      bool found = false;
      int indj=0;
      if (printout) cout << " number of trackjets " << jets->size() << endl;
      
      for (BasicJetCollection::const_iterator it = jets->begin(); it != jets->end(); it ++ ) {
	TVector3 jetcand;
	jetcand.SetPtEtaPhi(it->pt(), it->eta(), it->phi());
	TVector3 a; 
	a.SetPtEtaPhi(muon->track()->pt(), muon->track()->eta(), muon->track()->phi());
	double r = a.DeltaR(jetcand);
	if (printout) cout << "   " << indj << " pt/eta/phi: " << it->pt() << "/" << it->eta() << "/" << it->phi() 
			   << " ntrk = " << it->nConstituents()
			   << " r = " << r; 
	bool muonOnlyJet = (r < 1.e-6) && (it->nConstituents() == 1); 
	if (!muonOnlyJet && r<rmin && it->et() > fJetEtMin) {
	  rmin    = r;
	  matchedjet   = (*it).clone();
	  matchedJetIndex = indj;
	  found = true;
	  if (printout) cout << " matched!";
	}
	if (muonOnlyJet) {
	  if (printout) cout << " this is a muon-only jet, do not match"; 
	}
	if (printout) cout << endl;
	indj++;
      } 
	  
      // -- found a nearby track jet: remove the muon
      if (found) {
	TLorentzVector vect;
	vect.SetPtEtaPhiE(matchedjet->pt(), matchedjet->eta(), matchedjet->phi(), matchedjet->energy());
	bool foundmuon = false;
	std::vector< const reco::Candidate * > Constituent = matchedjet->getJetConstituentsQuick();
	if (printout) cout << "  found a trackjet with nconsti = " << Constituent.size() << " close to muon, jet pt/eta/phi = " 
			   << matchedjet->pt() << "/" << matchedjet->eta() << "/" << matchedjet->phi()
			   << endl;
	
	for (unsigned int i=0; i< Constituent.size(); i++) {
	  unsigned int idx  = 99999;
	  const reco::Candidate * consti = Constituent[i];
	  if (consti) {
    	    for (unsigned int j = 0; j < candidates1Handle->size(); ++ j ) {
	      const Candidate &p2 = (*candidates1Handle)[j];
	      const Candidate *  p1 = &p2;
	      if ((TMath::Abs(1. - consti->pt()/p1->pt()) < 0.001) && (TMath::Abs(1. - consti->eta()/p1->eta()) < 0.001)) {
		idx = j;
		break;
	      }
	    }

	    if (printout) cout << "            consti " << i 
			       << " with trk index " << idx
			       << " and trk pt/eta/phi: " << consti->pt() << "/" << consti->eta() << "/" << consti->phi() << " " 
			       << endl;

		
	    if (idx == (muon->track()).index()) {
	      foundmuon = true;
	    }
	  } 
	} // jet consti
	
	    
	if (foundmuon) {
	  pTrack           = gHFEvent->addSigTrack();
	  pTrack->fInt1    = 100100; // type
	  pTrack->fInt2    = matchedJetIndex;
	  
	  pTrack->fDouble1 = -99.;   // ptrel
	  TVector3 a; 
	  a.SetPtEtaPhi(muon->track()->pt(), muon->track()->eta(), muon->track()->phi());
	  pTrack->fPlab    = a;
	  
	  TLorentzVector muvect;
	  muvect.SetPtEtaPhiM(muon->track()->pt(), muon->track()->eta(), muon->track()->phi(), mmuon);
	  vect = vect - muvect;
	  pTrack->fDouble1 = muvect.Perp(vect.Vect()); //Ptrel in respct to the corrected jets direction
	  if (printout) cout << "--> matched muon with ptrel = " << muvect.Perp(vect.Vect()) << endl;
	  
	  //define direction
	  GlobalVector direction(vect.X(),vect.Y(),vect.Z());
	  const TransientTrack & transientTrack = builder->build(&(*(muon->track()))); 
	  
	  pTrack->fIndex = muon->track().index();
	  pTrack->fQ     = transientTrack.charge(); 
	  TAnaMuon *pM = gHFEvent->getSimpleTrackMuon(muon->track().index()); 
	  if (pM) {
	    pTrack->fMuID = pM->fMuID;
	  } else {
	    pTrack->fMuID = 0;
	  }
	} // found muon	in jet
	if (printout) cout << "--> done with muon " << muon->track().index() << endl;
      } // found jet
    } // jetsH valid
    index++;  
  } // muons  
}


// ------------ method called once each job just before starting event loop  ------------
//void  HFInclBMuonTrackJets::beginJob(const edm::EventSetup& setup) {
void  HFInclBMuonTrackJets::beginJob() {
  fNevt=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFInclBMuonTrackJets::endJob() { 
  cout << "HFInclBMuonTrackJets>     Summary: Events processed: " << fNevt << endl;
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFInclBMuonTrackJets);
