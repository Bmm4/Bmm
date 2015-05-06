#include <iostream>

#include "HFDumpSignal.h"

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
HFDumpSignal::HFDumpSignal(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fJetMatch(iConfig.getUntrackedParameter<double>("jetmatch", 0.5)),
  fJetEtMin(iConfig.getUntrackedParameter<double>("jetetmin", 1)),
  fMuonLabel(iConfig.getUntrackedParameter<string>("muonLabel", string("globalMuons"))),
  fJetsLabel(iConfig.getUntrackedParameter<string>("jetsLabel", string("ak5TrackJets"))),
  fTracksLabel(iConfig.getUntrackedParameter<string>("tracksLabel", string("trackCandidates"))),
  fVertexLabel(iConfig.getUntrackedParameter<string>("vertexLabel", string("offlinePrimaryVerticesWithBS"))),
  fSimVertexLabel(iConfig.getUntrackedParameter<string>("simvertexLabel", string("g4SimHits"))),
  fVertexMinNdof(iConfig.getUntrackedParameter<double>("vertexMinNdof", 4.))
   {
  using namespace std;
  fAllowGlobalOnly=iConfig.getUntrackedParameter<bool>("allowGlobalOnly", false);
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpSignal constructor" << endl;
  cout << "--- " << fMuonLabel.c_str() << endl;
  if (fAllowGlobalOnly){
    cout << "!!! global-only muons allowed,  not required to be tracker muons !!!!!" << endl;
  }else{
    cout << "all muons will be required to be tracker muons " << endl;
  }
  cout << "----------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
HFDumpSignal::~HFDumpSignal() {
  
}


// ----------------------------------------------------------------------
void HFDumpSignal::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  int RunNumber=iEvent.id().run();
  int nn=iEvent.id().event();
  bool inspect=(RunNumber==133874)&&((nn==32751514)||(nn==20949710)||(nn==11880401));
  nevt++;
  

  //muon collection
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel( fMuonLabel.c_str(), muons);
  if( !muons.isValid()) { cout<<"****** no "<<fMuonLabel<<endl; return; }

  //primary vertex
  Handle<reco::VertexCollection> primaryVertex;
  iEvent.getByLabel(fVertexLabel.c_str(),primaryVertex); 
 if( !primaryVertex.isValid()) { cout<<"****** no "<<fVertexLabel<<endl; return; }


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
 // iEvent.getByLabel( fSimVertexLabel.c_str(), simVtxs);
 if( iEvent.getByLabel( fSimVertexLabel.c_str(), simVtxs) && simVtxs.isValid()){
   simvFound = (simVtxs->size() != 0); 
   if(simvFound)  simv0=(*simVtxs)[0].position(); 
 }


 // TrackJets 
 Handle<BasicJetCollection> jetsH;
 iEvent.getByLabel(fJetsLabel.c_str(),jetsH);
 // bool jetvalid=true;
 if( !jetsH.isValid()) { cout<<"****** no "<<fJetsLabel<<endl;  }

 //tracks (jet constituents)
 Handle<reco::CandidateView> candidates1Handle;
 iEvent.getByLabel(fTracksLabel.c_str(), candidates1Handle); 
 if( !candidates1Handle.isValid()) { cout<<"****** no "<<fTracksLabel<<endl;  }
 

 //transient track builder
 edm::ESHandle<TransientTrackBuilder> builder;
 iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",builder);  
 

 /////////////////////   refitted track stuff   ///////////////////////
 Handle<reco::TrackCollection> generalTracksH;
 Handle<reco::TrackCollection> refittedTracksH;
 std::map<double ,const reco::Track *> g2r;
 //FIXME bool useRefittedTracks=false;
 bool hasRefittedTracks=iEvent.getByLabel("TrackRefitter",refittedTracksH) && refittedTracksH.isValid();
 bool hasTracks=iEvent.getByLabel("generalTracks",generalTracksH) && generalTracksH.isValid();


 if(
    hasRefittedTracks && hasTracks
    &&(generalTracksH->size() == refittedTracksH->size())
    ){
   //std::cout << "HFDumpSignal: found refitted tracks " << generalTracksH->size() << " = " << refittedTracksH->size() << std::endl;
   for(unsigned int i=0; i<generalTracksH->size(); i++){
     g2r[generalTracksH->at(i).pt()] = &(refittedTracksH->at(i));
   }
   //FIXME useRefittedTracks=true;
 }else if(hasTracks){
   for(unsigned int i=0; i<generalTracksH->size(); i++){
     g2r[generalTracksH->at(i).pt()] = &(generalTracksH->at(i));  // map to itself if no refitted tracks available
   }
 }else{
   cout << "HFDumpSignal: no tracks, what are we doing here? " << endl;
 }
 /////////////////////   refitted track stuff   ///////////////////////

  TAnaTrack *pTrack; 
  int index = 0;
  if (fVerbose > 0) cout << "==>HFDumpSignal> nMuons =" << muons->size() << endl;




  if(inspect){
  int muindex=0;
  cout << "==>HFDumpSignal> nMuons =" << muons->size() << endl;
  for (  reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++ muon ) { 
    cout << muindex << ")  g=" << muon->isGlobalMuon() << " t=" << muon->isTrackerMuon() << ":";
    muindex++;
    cout << Form(" ( %8.1f %4.1f %4.1f )   ", muon->pt(), muon->phi(), muon->eta());
    if(muon->isGlobalMuon()&&muon->isTrackerMuon()){
      reco::TrackRef InnerTrack=muon->innerTrack();
      cout << Form("%3d    %2d    %5.1f   %5.1f   %6.4f   %6.4f   %2d  ",
                   int((muon->globalTrack())->hitPattern().numberOfValidStripHits()+(muon->globalTrack())->hitPattern().numberOfValidPixelHits()),
                   int((muon->globalTrack())->hitPattern().pixelLayersWithMeasurement()),
                   float((muon->globalTrack())->chi2()/(muon->globalTrack())->ndof()),
                   float(muon->innerTrack()->normalizedChi2()),
                   float(0),
                   float(0),
                   int((muon->globalTrack())->hitPattern().numberOfValidMuonHits()));
    }
    cout  << endl;
  }
  }



  for (  reco::MuonCollection::const_iterator muon = muons->begin(); muon != muons->end(); ++ muon ) { 
    //if (muon->isGlobalMuon()&&muon->isTrackerMuon()) { // add only if global and tracker muons
    if (muon->isGlobalMuon()&& (muon->isTrackerMuon() || fAllowGlobalOnly)) { // add all global muons for tests
      
      //      reco::TrackRef gt=muon->innerTrack();

      if (fVerbose > 0) cout << "==>HFDumpSignal> found global muon " << endl;
      pTrack            = gHFEvent->addSigTrack();

      /* this goes/went into HFDUmpMuons!

      if (muon->isTrackerMuon()){  pTrack->fKaID   = 1; pTrack->fMuType   = 0;}else{ pTrack->fKaID   = 0; pTrack->fMuType   = 5;}
      if (muon->isGlobalMuon()){ pTrack->fMCID   = 1;}else{pTrack->fMCID=0;}

      pTrack->fMuType   = 0;                            //0=RECO, 1=L1, 2=HLTL2, 3=HLTL3 
      pTrack->fMuID     = (muon->track()).index();      //index of muon track in RECO track block
      pTrack->fIndex    = index;                        //index in muon block
      pTrack->fGenIndex = -1;                           //not used here
      pTrack->fQ        = muon->charge();               //charge
      pTrack->fPlab.SetPtEtaPhi(muon->pt(),
				muon->eta(),
				muon->phi()
				);

      
      pTrack->fChi2        = (muon->globalTrack())->chi2();
      pTrack->fDof         = int((muon->globalTrack())->ndof());
      pTrack->fHits        = (muon->globalTrack())->numberOfValidHits();  

      pTrack->fMuonSelector = 0;
      bool laststation = muon::isGoodMuon(*muon, muon::TMLastStationAngTight);
      bool tightmuon = muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight);
      bool matched = muon->numberOfMatches()>1;
      if (laststation) pTrack->fMuonSelector +=1;
      if (tightmuon)  pTrack->fMuonSelector  +=2;
      if (matched)  pTrack->fMuonSelector  +=4;
      if (muon->isTrackerMuon()) pTrack->fMuonSelector +=8;
      
      pTrack->fMuonCSCHits = (muon->globalTrack())->hitPattern().numberOfValidMuonCSCHits();
      pTrack->fMuonDTHits  = (muon->globalTrack())->hitPattern().numberOfValidMuonDTHits();
      pTrack->fMuonRPCHits = (muon->globalTrack())->hitPattern().numberOfValidMuonRPCHits();
      pTrack->fMuonHits    = (muon->globalTrack())->hitPattern().numberOfValidMuonHits();
      
      pTrack->fBPIXHits    = (muon->globalTrack())->hitPattern().numberOfValidPixelBarrelHits(); 
      pTrack->fFPIXHits    = (muon->globalTrack())->hitPattern().numberOfValidPixelEndcapHits();
      pTrack->fPixelHits   = (muon->globalTrack())->hitPattern().numberOfValidPixelHits();
      
      pTrack->fStripHits   = (muon->globalTrack())->hitPattern().numberOfValidStripHits(); 
      pTrack->fTECHits     = (muon->globalTrack())->hitPattern().numberOfValidStripTECHits(); 
      pTrack->fTIBHits     = (muon->globalTrack())->hitPattern().numberOfValidStripTIBHits();
      pTrack->fTIDHits     = (muon->globalTrack())->hitPattern().numberOfValidStripTIDHits();
      pTrack->fTOBHits     = (muon->globalTrack())->hitPattern().numberOfValidStripTOBHits(); 
      
      pTrack->fBPIXLayers  = (muon->globalTrack())->hitPattern().pixelBarrelLayersWithMeasurement();  
      pTrack->fFPIXLayers  = (muon->globalTrack())->hitPattern().pixelEndcapLayersWithMeasurement();  
      pTrack->fPixelLayers = (muon->globalTrack())->hitPattern().pixelLayersWithMeasurement();  
      pTrack->fStripLayers = (muon->globalTrack())->hitPattern().stripLayersWithMeasurement(); 
      pTrack->fTECLayers   = (muon->globalTrack())->hitPattern().stripTECLayersWithMeasurement();
      pTrack->fTIBLayers   = (muon->globalTrack())->hitPattern().stripTIBLayersWithMeasurement(); 
      pTrack->fTIDLayers   = (muon->globalTrack())->hitPattern().stripTIDLayersWithMeasurement();
      pTrack->fTOBLayers   = (muon->globalTrack())->hitPattern().stripTOBLayersWithMeasurement();
     
      pTrack->fTrChi2norm  = muon->innerTrack()->normalizedChi2(); // chi2 of the tracker track
      pTrack->fExpectedHitsInner = (muon->innerTrack())->trackerExpectedHitsInner().numberOfHits();
      pTrack->fExpectedHitsOuter = (muon->innerTrack())->trackerExpectedHitsOuter().numberOfHits();
      pTrack->fLostHits    = (muon->innerTrack())->hitPattern().numberOfLostTrackerHits();
      pTrack->fValidHitInFirstPixelBarrel=(muon->innerTrack())->hitPattern().hasValidHitInFirstPixelBarrel();
 
      pTrack->fDxybs       = muon->innerTrack()->dxy(bs);          // Dxy relative to the beam spot
      pTrack->fDzbs          = muon->innerTrack()->dz(bs);        // dz relative to bs or o if not available


      // find the nearest primary vertex
      double dmin=9999.;
      double ztrack=muon->vertex().Z();
      const reco::Vertex* pVertex =&(*primaryVertex->begin()); // if selection below fails, take the first one
      for (reco::VertexCollection::const_iterator iv = primaryVertex->begin(); iv != primaryVertex->end(); ++iv) {
	if (iv->ndof()>fVertexMinNdof){
	  if (fabs(ztrack-iv->z()) < dmin){
	    dmin=fabs(ztrack-iv->z());
	    pVertex=&(*iv);
	  }
	}
      }

      // impact parameter of the tracker track relative  to the PV  
      pTrack->fDxypv        = muon->innerTrack()->dxy(pVertex->position()); 
      pTrack->fDzpv         = muon->innerTrack()->dz(pVertex->position());       
 
      pTrack->fDxyE        = muon->innerTrack()->dxyError();  // error on Dxy
      pTrack->fDzE          = muon->innerTrack()->dzError();

      pTrack->fMuonIsosumPT    = muon->isolationR03().sumPt;  // muon isolations
      pTrack->fMuonIsoEcal     = muon->isolationR03().emEt;
      pTrack->fMuonIsoHcal     = muon->isolationR03().hadEt;
      pTrack->fMuonCalocomp     = muon->caloCompatibility();
      pTrack->fMuonCalocomp     = muon->caloCompatibility();
      pTrack->fVertex.SetXYZ(muon->vertex().X(),
				muon->vertex().Y(),
				muon->vertex().Z()
				);
      pTrack->fMuonVertexChi2=muon->vertexChi2();
 
      pTrack->fDxysimv = -9999;  // first  sim vertex 
      pTrack->fDzsimv = -9999;  
      if(simvFound) {
	pTrack->fDxysimv =   muon->innerTrack()->dxy(simv0); 
 	pTrack->fDzsimv =  muon->innerTrack()->dz(simv0); 
      }

      if(hasRefittedTracks){
	if(useRefittedTracks){
	  //const reco::Track* refitted = g2r[muon->innerTrack()->pt()];
	  const reco::Track refitted = generalTracksH->at((muon->innerTrack()).index());
	  pTrack->fDxypv        = refitted.dxy(bs);
	  pTrack->fDxysimv      = refitted.dxy(simv0);
	}else{
	  pTrack->fDxypv        = -99999.;
	  pTrack->fDxysimv      = -99999.;
	}
      }

     //
      pTrack->fLip      = -9999;
      pTrack->fLipE     = -9999;
      pTrack->fTip      = -9999;
      pTrack->fTipE     = -9999;   // 3d and transverse impact parameters
      pTrack->fTip3d    = -9999;
      pTrack->fTip3dE   = -9999;   // 3d and transverse impact parameters 
      pTrack->fIndxj    = -9999;   //   index of the closest jets after selection
      pTrack->fMuonPtrel    = -9999;   //  ptrel to this jet, can be racalculated latter
      pTrack->fDxypvnm = -9999;  // dx pTrack->fDxypvnm = muon->innerTrack()->dxy(vpnm);y relative to vertex without the muon
 
  */

 
      //calculate signed transverse ip
      //get closest jet, use it to calculate 3dIP 

      if(jetsH.isValid()) {
	const BasicJetCollection *jets   = jetsH.product();  
	if(candidates1Handle.isValid() ){
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if ( (muon->track()).index() <= candidates1Handle->size() ) {
	  
	  double rmin = fJetMatch; 
	  BasicJet* matchedjet=0;
	  bool found = false;
          int indj=0;
	  if (fVerbose > 0) cout<<" trackjets "<<jets->size()<<endl;
	  for ( BasicJetCollection::const_iterator it = jets->begin(); it != jets->end(); it ++ ) {
            
	    TVector3 jetcand;
	    jetcand.SetPtEtaPhi(it->pt(), it->eta(), it->phi());
	    double r = (pTrack->fPlab).DeltaR(jetcand);
	    
	    if (r<rmin && it->et() > fJetEtMin) {
	      rmin    = r;
	      matchedjet   = (*it).clone();
	      found = true;
	      //FIXME              pTrack->fIndxj=indj;
	    }
	    indj++;
	  } 
	  if (found) {//found a track jet
	    
	    if (fVerbose > 0) cout << "==>HFDumpSignal> found a trackjet close to muon " << endl;
	    //subtract muon 
	    TLorentzVector vect;
	    vect.SetPtEtaPhiE(matchedjet->pt(), matchedjet->eta(), matchedjet->phi(), matchedjet->energy());
	    
	    bool foundmuon = false;
	    std::vector< const reco::Candidate * > Constituent = matchedjet->getJetConstituentsQuick();
	    if (fVerbose > 0) cout<<" matchedjet tracks "<<Constituent.size()<<endl;
	    //5555555555555555555555555555555555555555555
	    //jet constituents
	    for (unsigned int i=0; i< Constituent.size(); i++) {
	      
	      unsigned int idx  = 99999;
	      const reco::Candidate * consti = Constituent[i];
	      if (consti) {
		
		//get index of consti in all tracks
		for (unsigned int j = 0; j < candidates1Handle->size(); ++ j ) {
		  const Candidate &p2 = (*candidates1Handle)[j];
		  const Candidate *  p1 = &p2;
		  if (consti->pt() == p1->pt() && consti->phi() == p1->phi() && consti->eta() == p1->eta() ) 
		    idx = j;
		  
		}
		//check whether the index of consti in tracks is the same as the index of the muon track
		if (idx == (muon->track()).index()) {
		  foundmuon = true;
		  if (fVerbose) cout << "muon " << idx << " pt " << consti->pt() << " phi " << consti->phi() << " eta " << consti->eta() << endl;
		}
	      } 
	    } // jet const
	  
	   	    
	    if (foundmuon) {
	      TLorentzVector muvect;
	      muvect.SetPtEtaPhiM(muon->track()->pt(), muon->track()->eta(), muon->track()->phi(), mmuon);
	      vect = vect - muvect;
	      //FIXME pTrack->fMuonPtrel=muvect.Perp(vect.Vect()); //Ptrel in respct to the corrected jets direction
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


	//// find a vertex to which muon belongs and remove it from vertex track list, get myVertex
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
	      if (fVerbose > 0) cout<<"ipv "<<ipv<<" "<<ddr<<"= dr  pvtr pt"<< pvtr->pt()<<" pTrack->fPlabpb  "<<  pTrack->fPlab.Perp()<<endl;  
	      
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
	pTrack->dump(); 
	/* FIXME
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
  
  if (fVerbose > 0) cout<<" --->HFDumpSignal "<< index << endl;
}

// ------------ method called once each job just before starting event loop  ------------
//void  HFDumpSignal::beginJob(const edm::EventSetup& setup) {
void  HFDumpSignal::beginJob() {
  nevt=0;
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFDumpSignal::endJob() { 
  cout << "HFDumpSignal>     Summary: Events processed: " << nevt << endl;
 
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpSignal);
