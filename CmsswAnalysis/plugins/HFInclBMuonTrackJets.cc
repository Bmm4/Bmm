#include <iostream>

#include "HFInclBMuonTrackJets.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
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
  
  if (fVerbose > -1) cout << "----------------------------------------------------------------------" << endl;

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
  int index = 0;
  if (fVerbose > 0) cout << "==>HFInclBMuonTrackJets> nMuons =" << muons->size() << endl;

  for (  reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++ muon ) { 
    //if (muon->isGlobalMuon()&&muon->isTrackerMuon()) { // add only if global and tracker muons
    if (muon->isGlobalMuon()&& (muon->isTrackerMuon() || fAllowGlobalOnly)) { // add all global muons for tests
      
      //      reco::TrackRef gt=muon->innerTrack();

      if (fVerbose > 0) cout << "==>HFInclBMuonTrackJets> found global muon " << endl;
      pTrack           = gHFEvent->addSigTrack();
      pTrack->fInt1    = 100100; // type
      pTrack->fInt2    = -1;     // jet index

      pTrack->fDouble1 = -99.;   // ptrel
      TVector3 a; 
      a.SetPtEtaPhi(muon->track()->pt(), muon->track()->eta(), muon->track()->phi());
      pTrack->fPlab    = a;

      

      cout << "MUON pt/eta = " << muon->track()->pt() << "/" <<  muon->track()->eta() << endl;


      if(jetsH.isValid()) {
	//	const TrackJetCollection *jets   = jetsH.product();  
	const BasicJetCollection *jets   = jetsH.product();  
	if(candidates1Handle.isValid() ){
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	  if ( (muon->track()).index() <= candidates1Handle->size() ) {
	  
	    double rmin = fJetMatch; 
	    //	    TrackJet* matchedjet=0;
	    BasicJet* matchedjet=0;
	    bool found = false;
	    int indj=0;
	    if (fVerbose > 0) cout<<" trackjets "<<jets->size()<<endl;
	    //	    for ( TrackJetCollection::const_iterator it = jets->begin(); it != jets->end(); it ++ ) {
	    for ( BasicJetCollection::const_iterator it = jets->begin(); it != jets->end(); it ++ ) {
            
	      TVector3 jetcand;
	      jetcand.SetPtEtaPhi(it->pt(), it->eta(), it->phi());
	      double r = (pTrack->fPlab).DeltaR(jetcand);
	    
	      if (r<rmin && it->et() > fJetEtMin) {
		rmin    = r;
		matchedjet   = (*it).clone();
		found = true;
		pTrack->fInt2 = indj;
	      }
	      indj++;
	    } 
	    if (found) {//found a track jet
	    
	      if (fVerbose > 0) cout << "==>HFInclBMuonTrackJets> found a trackjet close to muon, jet pt/eta = " 
				     << matchedjet->pt() << "/" << matchedjet->eta()
				     << endl;
	      //subtract muon 
	      TLorentzVector vect;
	      vect.SetPtEtaPhiE(matchedjet->pt(), matchedjet->eta(), matchedjet->phi(), matchedjet->energy());
	    
	      bool foundmuon = false;
	      std::vector< const reco::Candidate * > Constituent = matchedjet->getJetConstituentsQuick();
	      if (fVerbose > 0) cout<<" matchedjet tracks "<<Constituent.size()<<endl;
	      //jet constituents
	      for (unsigned int i=0; i< Constituent.size(); i++) {
	      
		unsigned int idx  = 99999;
		const reco::Candidate * consti = Constituent[i];
		if (consti) {
		
		  //get index of consti in all tracks
		  for (unsigned int j = 0; j < candidates1Handle->size(); ++ j ) {
		    const Candidate &p2 = (*candidates1Handle)[j];
		    const Candidate *  p1 = &p2;
		    if ((TMath::Abs(1. - consti->pt()/p1->pt()) < 0.001) && (TMath::Abs(1. - consti->eta()/p1->eta()) < 0.001)) {
		      idx = j;
		    }


		  }
		  //check whether the index of consti in tracks is the same as the index of the muon track
		  if (idx == (muon->track()).index()) {
		    foundmuon = true;
		  }
		} 
	      } // jet consti
	  
	   	    
	      if (foundmuon) {
		TLorentzVector muvect;
		muvect.SetPtEtaPhiM(muon->track()->pt(), muon->track()->eta(), muon->track()->phi(), mmuon);
		vect = vect - muvect;
		pTrack->fDouble1 = muvect.Perp(vect.Vect()); //Ptrel in respct to the corrected jets direction
		if (fVerbose) cout << "matched muon " 
				   << " pt " << muon->pt()  << " eta " << muon->eta() << " phi " << muon->phi()
				   << " ptrel = " << muvect.Perp(vect.Vect())
				   << endl;
	      }
	      //define direction
	      GlobalVector direction(vect.X(),vect.Y(),vect.Z());
	      const TransientTrack & transientTrack = builder->build(&(*(muon->track()))); 
	    
	      /*FIXME
		pTrack->fTip3d    = IPTools::signedImpactParameter3D(transientTrack,direction,*pVertex).second.value();
		pTrack->fTip3dE   = IPTools::signedImpactParameter3D(transientTrack, direction, *pVertex).second.error();
		pTrack->fTip      = IPTools::signedTransverseImpactParameter(transientTrack, direction, *pVertex).second.value();
		pTrack->fTipE     = IPTools::signedTransverseImpactParameter(transientTrack, direction, *pVertex).second.error();	     
		FIXME*/
	    }  //found trackjet in dR	  
	  } // tracks=muon check
	} // tracks valid
      } // trackjets valid
      
      
      // find a vertex to which muon belongs and remove it from vertex track list, get myVertex
      // match by dR<0.1
      double drmin=0.1;
      bool findmuon=false;
      double dzmax=10000.;
      vector<TransientTrack> mytracks;
      const  Vertex* pvr0=0;

      if( primaryVertex.isValid()){
	int ipv=0;
	for (reco::VertexCollection::const_iterator pv =  primaryVertex->begin(); pv!= primaryVertex->end(); ++pv){  
	  pvr0=0;
	  double ddz= muon->innerTrack()->dz(pv->position());
	  if(abs(ddz)<dzmax) {
	    dzmax=ddz;
	    pvr0 =&(*pv);
	  }  

	  for(std::vector<TrackBaseRef>::const_iterator pvt = pv->tracks_begin(); pvt != pv->tracks_end(); pvt++) {
	    TVector3 trkcand;
	    TrackRef pvtr=pvt->castTo<TrackRef>();
	    trkcand.SetPtEtaPhi(pvtr->pt(), pvtr->eta(), pvtr->phi());
	      
	    double ddr = (pTrack->fPlab).DeltaR(trkcand);
	    //ul if (fVerbose > 0) cout<<"ipv "<<ipv<<" "<<ddr<<"= dr  pvtr pt"<< pvtr->pt()<<" pTrack->fPlabpb  "<<  pTrack->fPlab.Perp()<<endl;  
	      
	    if(ddr<drmin) {
	      // found muon
	      findmuon=true;
	    } else {
	      // fill mytracks
	      TransientTrack  transientTrack =  builder->build(pvtr);
	      if (beamSpotCollection.isValid()) {
		reco::BeamSpot vertexBeamSpot_= *beamSpotCollection;
		transientTrack.setBeamSpot(vertexBeamSpot_);
	      }
	      mytracks.push_back(transientTrack);
	    } //ddr
	  } //pvt
	    

	  // build vertex wo muon
	  if(findmuon) {
	    if (fVerbose > 0) cout<<" muon is found in PV"<<endl;
	    AdaptiveVertexFitter* theFitter=new AdaptiveVertexFitter();
	    // this can be used for rererecos when bs is well defined
	    // if (beamSpotCollection.isValid()) {
	    //reco::BeamSpot vertexBeamSpot_= *beamSpotCollection;
	    // TransientVertex myVertex = theFitter->vertex(mytracks, vertexBeamSpot_); //fit with beam-constraint
	    //     }
	    TransientVertex myVertex = theFitter->vertex(mytracks); // fit without beam-constraint
	    if(myVertex.isValid()&&myVertex.degreesOfFreedom()>1){
	      // require at least 2 degree of freedom
	      //FIXME math::XYZPoint vpnm = math::XYZPoint(myVertex.position().x(),myVertex.position().y(),myVertex.position().z());
	      //FIXME pTrack->fDxypvnm = muon->innerTrack()->dxy(vpnm); // transverse IP without muon
		
	      // using z beam  direction
	      GlobalVector direction(0.,0.,1.);
	      const TransientTrack & transientTrack = builder->build(&(*(muon->track()))); 
	      const Vertex vv(myVertex); //PV without muon
	      /* FIXME
		 pTrack->fTip3d    = IPTools::signedImpactParameter3D(transientTrack,direction,vv).second.value();
		 pTrack->fTip3dE   = IPTools::signedImpactParameter3D(transientTrack, direction, vv).second.error();

		 pTrack->fLip      = IPTools::signedTransverseImpactParameter(transientTrack, direction, vv).second.value();
		 pTrack->fLipE     = IPTools::signedTransverseImpactParameter(transientTrack, direction, vv).second.error();	 
		 FIXME */
	    } // valid myvertex
	    delete theFitter;
	    break; // break after the first found vertex with muon
	  }//foundmuon
	  ipv++;
	} //pv
      } //pv valid
	
	// use vertex closest in z if muon not found in pv collection
	// in this case the pTrack->fLip stays indefined in order to distinguish
      if(!findmuon&&pvr0) {
	//FIXME math::XYZPoint vpnm = math::XYZPoint(pvr0->position().x(),pvr0->position().y(),pvr0->position().z());
	//FIXME pTrack->fDxypvnm = muon->innerTrack()->dxy(vpnm); // transverse IP for closest Iin z
	//FIXME pTrack->fLip =-99999.; // different from default -9999
      }

      if (fVerbose > 0) {
	/* FIXME
	   pTrack->dump(); 
	   cout << " dump muon" << endl;
	   cout << "pTrack->fTip3d=" << pTrack->fTip3d << "+-" << pTrack->fTip3dE << endl;
	   cout << "pTrack->fTip=" << pTrack->fTip << "+-" << pTrack->fTipE << endl;
	   cout << "pTrack->fLip=" << pTrack->fLip << "+-" << pTrack->fLipE << endl;
	   cout << "pTrack->fDxypv= " <<pTrack->fDxypv << endl;
	   cout << "pTrack->fDxypvnm= " <<pTrack->fDxypvnm << endl;
	   cout << "pTrack->fDxysimv= " <<pTrack->fDxysimv << endl;
	   FIXME */
      }
	
    }//is global muon
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    index++;  
    
    
  } //muon
  
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
