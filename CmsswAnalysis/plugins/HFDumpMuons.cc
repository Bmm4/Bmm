// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFDumpMuons
// -----------
//
// 2016/01/21 Urs Langenegger      derive from HFVirtualDecay,
//                                 migrate to "consumes",
//                                 remove calomuons and 
// stone age  Urs Langenegger      first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpMuons.h"

#include <iostream>
#include <cmath>

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"
#include "Bmm/RootAnalysis/rootio/TTrgObj.hh"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "CommonTools/Statistics/interface/ChiSquared.h"

#include <TFile.h>
#include <TH1.h>

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace edm;
using namespace reco;

// -- sort the vector with xpTracks
static bool dist_less(const xpTrack &x, const xpTrack &y) {
  return (x.dist < y.dist); 
}


// ----------------------------------------------------------------------
HFDumpMuons::HFDumpMuons(const edm::ParameterSet& iConfig):
  HFVirtualDecay(iConfig),
  fCaloMuonsLabel(iConfig.getUntrackedParameter<InputTag>("calomuonsLabel")),
  fMaxTrackDistToStore(iConfig.getUntrackedParameter<double>("maxTrackDist",0.1)),
  fDocaVertex(iConfig.getUntrackedParameter<double>("docaVertex",0.05)),
  fKeepBest(iConfig.getUntrackedParameter<int>("keepBest",3)),
  fMaxCandTracks(iConfig.getUntrackedParameter<int>("maxCandTracks",3)),
  fpropM1(iConfig.getParameter<edm::ParameterSet>("propM1")),
  fpropM2(iConfig.getParameter<edm::ParameterSet>("propM2")) {
  dumpConfiguration();
  
  fTokenCaloMuon      = consumes<CaloMuonCollection>(fCaloMuonsLabel);
  
}


// ----------------------------------------------------------------------
void HFDumpMuons::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpMuons configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  caloMuonsLabel           " << fCaloMuonsLabel << endl;
  cout << "---  fMaxTrackDistToStore:    " << fMaxTrackDistToStore << endl;
  cout << "---  docaVertex:              " << fDocaVertex << endl;
  cout << "---  keepBest:                " << fKeepBest << endl;
  cout << "---  maxCandTracks:           " << fMaxCandTracks << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
HFDumpMuons::~HFDumpMuons() {

}

// ----------------------------------------------------------------------
void HFDumpMuons::beginRun(const Run& iRun, const EventSetup& iSetup) {
  cout << "HFDumpMuons::beginRun" << endl;
  fpropM1.init(iSetup);
  fpropM2.init(iSetup);
}

// ----------------------------------------------------------------------
void HFDumpMuons::analyze(const Event& iEvent, const EventSetup& iSetup) {
  if (1) cout << "--- HFDumpMuons -------------------------------------------------------------------" << endl;

  if (fVerbose > 0) cout << "==>HFDumpMuons> new event " << endl;

  try {
    HFVirtualDecay::analyze(iEvent, iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFfDumpMuons> " << e.fMsg << endl;
    return;
  }

  fListBuilder->setMinPt(-1.);

  extrapolateTracks();
  
  int im(0);
  for (MuonCollection::const_iterator iMuon = fMuonCollection->begin(); iMuon != fMuonCollection->end(); ++ iMuon ) {
    fillMuon(*iMuon, im); 
    ++im;
  }

  if (fVerbose > 0) {
    for (int im = 0; im < gHFEvent->nMuons(); ++im) {
      gHFEvent->getMuon(im)->dump();
      cout << "px = " << gHFEvent->getMuon(im)->fPlab.Px() 
	   << " py = " << gHFEvent->getMuon(im)->fPlab.Py()
	   << " pz = " << gHFEvent->getMuon(im)->fPlab.Pz()
	   << endl;
    }
  }
  
}


// ----------------------------------------------------------------------
void HFDumpMuons::fillMuon(const reco::Muon& rm, int im) {

  TrackRef gTrack = rm.globalTrack();
  TrackRef iTrack = rm.innerTrack();
  TrackRef oTrack = rm.outerTrack();

  TAnaMuon *pM = gHFEvent->addMuon();    

  
  if (rm.innerTrack().isNonnull()) {
    Track trk(*iTrack);
    fillAnaTrack(pM, trk, rm.innerTrack().index(), -2, &fVertexCollection, fMuonCollection, &fBeamSpot); 
  } else {
    pM->fIndex = -23;
  }
  pM->fMuIndex = im;
  pM->fMuID    = muonID(rm);
  pM->fQ       = rm.charge();

  pM->fMuonChi2         = rm.combinedQuality().trkKink;
  pM->fTimeInOut        = rm.time().timeAtIpInOut; 
  pM->fTimeInOutE       = rm.time().timeAtIpInOutErr; 
  pM->fTimeOutIn        = rm.time().timeAtIpOutIn; 
  pM->fTimeOutInE       = rm.time().timeAtIpOutInErr; 
  pM->fTimeNdof         = rm.time().nDof;
  pM->fNmatchedStations = rm.numberOfMatchedStations();

  bool isGlobalMuon = muon::isGoodMuon(rm, muon::AllGlobalMuons);
  if (isGlobalMuon != rm.isGlobalMuon()) {
    cout << "?????????? error ????  isGlobalMuon() != muon::isGoodMuon(rm, muon::AllGlobalMuons) ???" << endl;
  }
  
  // -- variables for MVA muon ID
  if (gTrack.isNonnull() && iTrack.isNonnull()) {
    const HitPattern track_hp  = iTrack->hitPattern();
    //  changes in 72X
    //const HitPattern exp_track_out_hp = iTrack->trackerExpectedHitsOuter();
    reco::MuonQuality muQuality = rm.combinedQuality();
    
    pM->fItrkValidFraction       = iTrack->validFraction(); //1
    pM->fGtrkNormChi2            = gTrack->normalizedChi2(); //2
    pM->fChi2LocalPosition       = muQuality.chi2LocalPosition; //3
    //pM->fNumberOfLostTrkHits     = exp_track_out_hp.numberOfLostTrackerHits(); //4 change for 72X
    // should I use this 
    //pM->fNumberOfLostTrkHits     = track_hp.numberOfLostTrackerHits(HitPattern::TRACK_HITS); //4
    // OR this 
    pM->fNumberOfLostTrkHits     = track_hp.numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS); //4
    pM->fSegmentComp             = muon::segmentCompatibility(rm); //5
    pM->fGtrkProb                = muQuality.glbTrackProbability; //6
    pM->fChi2LocalMomentum       = muQuality.chi2LocalMomentum; //6
    pM->fNumberOfValidTrkHits    = track_hp.numberOfValidTrackerHits(); //7
  }


  if (gTrack.isNonnull()) {
    Track trk(*gTrack);
    pM->fGlobalPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
    pM->fNvalidMuonHits = gTrack->hitPattern().numberOfValidMuonHits();


    pM->fNhitsDT  = gTrack->hitPattern().numberOfValidMuonDTHits();
    pM->fNhitsCSC = gTrack->hitPattern().numberOfValidMuonCSCHits();
    pM->fNhitsRPC = gTrack->hitPattern().numberOfValidMuonRPCHits();
    
  }

  if (iTrack.isNonnull()) {
    Track trk(*iTrack);
    pM->fInnerPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
    //    cout << " ipt = " << trk.pt() << " " <<  trk.eta() << " "  <<  trk.phi() << endl;
  }

  if (oTrack.isNonnull()) {
    Track trk(*oTrack);
    pM->fOuterPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
  }

  // -- propagate muons to muon system to get their impact point
  TVector3 muPosM1; 
  bool validM1(false); 

  if (isGlobalMuon && doExtrapolate(rm.pt(), rm.eta())) {
    TrajectoryStateOnSurface prop_M1 = fpropM1.extrapolate(rm);
    TrajectoryStateOnSurface prop_M2;
    if (oTrack.isNonnull()) {
      //      Track trk(*oTrack);
      prop_M2 = fpropM2.extrapolate(*oTrack);
    }
    
    if (prop_M1.isValid()) {
      pM->fPositionAtM1.SetXYZ(prop_M1.globalPosition().x(), prop_M1.globalPosition().y(), prop_M1.globalPosition().z());
    }
    if (prop_M2.isValid()) {
      pM->fPositionAtM2.SetXYZ(prop_M2.globalPosition().x(), prop_M2.globalPosition().y(), prop_M2.globalPosition().z());
      cout << "===> M2:   "
	   << Form(" %2d id: %+2d m = %4.1f pT/eta/phi = %6.3f/%+5.4f/%+5.4f, expol: %6.3f/%+5.4f/%+5.4f",
		   1, 3, 0.105, rm.pt(), rm.eta(), rm.phi(), rm.pt(), prop_M2.globalPosition().eta(), static_cast<float>(prop_M2.globalPosition().phi()))
	   << endl;

    }

    if (oTrack.isNonnull()) {
      TrajectoryStateOnSurface propOuter = fpropM1.extrapolate(*oTrack);
      if (propOuter.isValid()) {
	// -- FIXME: not sure whether the next two lines would not be better in prop_M1.isValid() block?!
	validM1 = true; 
	muPosM1.SetXYZ(propOuter.globalPosition().x(),propOuter.globalPosition().y(),propOuter.globalPosition().z());
	pM->fMuonTrackPosAtM1.SetXYZ(propOuter.globalPosition().x(),propOuter.globalPosition().y(),propOuter.globalPosition().z());
	pM->fMuonTrackPlabAtM1.SetXYZ(propOuter.globalMomentum().x(),propOuter.globalMomentum().y(),propOuter.globalMomentum().z());
      }
    }
  }

  // -- doca of close tracks to muon
  if (isGlobalMuon && iTrack.isNonnull() && fTTB.isValid()) {
    TwoTrackMinimumDistance md(TwoTrackMinimumDistance::SlowMode);
    Track trkMuon(*iTrack);
    TransientTrack transTrkMuon = fTTB->build(trkMuon);
	  
    for (size_t k = 0; k < fTracksHandle->size(); k++) {
      if (k == iTrack.index()) continue; // own track
		  
      TrackBaseRef bRefTrk(fTracksHandle, k);
      Track trk(*bRefTrk);
      TransientTrack transTrk = fTTB->build(trk);
		  
      md.calculate(transTrkMuon.initialFreeState(), transTrk.initialFreeState());
      if (md.distance() < fMaxTrackDistToStore) {
	pM->fNstTracks.insert(std::make_pair(k,md.distance()));
      }
    }
  }


  if (isGlobalMuon && validM1) {
    vector<xpTrack> xvec; 
    for (unsigned int i = 0; i < fXpTracks.size(); ++i) {
      xpTrack x = fXpTracks[i]; 
      if (x.idx == static_cast<int>(iTrack.index())) continue;
      x.dist = (muPosM1 - x.r).Mag();  
      xvec.push_back(x); 
    }
    
    // -- sort the vector & keep only the first TAnaMuon::NXPTRACKS
    sort(xvec.begin(), xvec.end(), dist_less);
    if (TAnaMuon::NXPTRACKS < xvec.size()) {
      for (unsigned int ii = 0; ii < TAnaMuon::NXPTRACKS; ++ii) pM->fXpTracks[ii] = xvec[ii]; 
    } else {
      for (unsigned int ii = 0; ii < xvec.size(); ++ii) pM->fXpTracks[ii] = xvec[ii]; 
    }
  }

  // do the muon vertex analysis
  if (isGlobalMuon && iTrack.isNonnull()) {
    std::set<unsigned> trks;
    double prob = -1.0;
    trks.insert(iTrack.index());
    findVertex(pM,&trks,&prob);
    pM->fVtxProb = prob;
    pM->fVtxTracks = trks;
  }
}

// ----------------------------------------------------------------------
void HFDumpMuons::fillCaloMuon(const reco::CaloMuon& rm, int im) {

  TAnaMuon *pM = gHFEvent->addMuon();    
  pM->fMuID    = 0;  // this assumes that fillCaloMuon is independent from the main muons d.k.

  if (rm.innerTrack().isNonnull()) {
    pM->fIndex = rm.innerTrack().index();
    pM->fQ        = rm.charge();
  } else {
    pM->fIndex = -23;
    pM->fQ     = 0;
  }
  pM->fMuIndex = im; 
  pM->fMuID   |= 0x1<<15;

  pM->fNhitsDT  = 0; 
  pM->fNhitsCSC = 0; 
  pM->fNhitsRPC = 0; 
  

  TrackRef iTrack = rm.innerTrack();

  if (iTrack.isNonnull()) {
    Track trk(*iTrack);
    pM->fInnerPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
  } 

}

// ----------------------------------------------------------------------
void HFDumpMuons::findVertex(TAnaMuon *anaMu, std::set<unsigned> *trkIcs, double *prob) {
  std::vector<TransientTrack> transTracks;
  std::vector<std::pair<double,unsigned> > bestTracks;
  std::map<int,float>::const_iterator mapIt;
  std::set<unsigned>::const_iterator it;
  std::set<unsigned> resultIcs;
  KalmanVertexFitter kvf;
  double best;
  unsigned ix;
	
  // build the transient tracks with 'trkIcs'
  for (it = trkIcs->begin(); it != trkIcs->end(); ++it) {
    TrackBaseRef baseRef(fTracksHandle, *it);
    Track trk(*baseRef);
    TransientTrack ttrack = fTTB->build(trk);
    transTracks.push_back(ttrack);
  }
	
  for (mapIt = anaMu->fNstTracks.begin(); mapIt != anaMu->fNstTracks.end(); ++mapIt) {
		
    if (mapIt->second >= fDocaVertex)
      continue;
		
    if (trkIcs->count(mapIt->first) > 0)
      continue; // already included
		
    trkIcs->insert(mapIt->first);
		
    TrackBaseRef baseRef(fTracksHandle, mapIt->first);
    Track trk(*baseRef);
    TransientTrack ttrack = fTTB->build(trk);
    transTracks.push_back(ttrack);
		
    TransientVertex vtx = kvf.vertex(transTracks);
    ChiSquared chi(vtx.totalChiSquared(), vtx.degreesOfFreedom());
    best = chi.probability();
    if (!TMath::IsNaN(best))
      bestTracks.push_back(make_pair(chi.probability(),mapIt->first));
		
    trkIcs->erase(mapIt->first);
    transTracks.pop_back();
  }
	
  // only iterate the most promosing 'keep'
  std::sort(bestTracks.begin(),bestTracks.end());
  if (bestTracks.size() > fKeepBest) bestTracks.erase(bestTracks.begin(),bestTracks.end()-fKeepBest);
	
  best = *prob;
  resultIcs = *trkIcs;
  for (ix = 0; ix < bestTracks.size(); ix++) {
		
    std::set<unsigned> curTracks = *trkIcs;
    double result = bestTracks[ix].first;
		
    curTracks.insert(bestTracks[ix].second);
		
    if (curTracks.size() < fMaxCandTracks)
      findVertex(anaMu,&curTracks,&result);
		
    if (best < result) {
      best = result;
      resultIcs = curTracks;
    }
  }
	
  // save
  *prob = best;
  *trkIcs = resultIcs;
} // findVertex()


// ----------------------------------------------------------------------
void HFDumpMuons::extrapolateTracks() {

  fXpTracks.clear(); 

  for (unsigned int k = 0; k < fTracksHandle->size(); ++k) {
    TrackBaseRef bRefTrk(fTracksHandle, k);
    Track trk(*bRefTrk);

    double eta = TMath::Abs(trk.eta()); 
    double pt  = trk.pt();
    
    // -- skip regions where extrapolation is futile
    if (!doExtrapolate(pt, eta)) continue;

    TrajectoryStateOnSurface tsos = fpropM1.extrapolate(trk);
    if (tsos.isValid()) {
      xpTrack x; 
      x.idx = k; 
      x.r = TVector3(tsos.globalPosition().x(), tsos.globalPosition().y(), tsos.globalPosition().z()); 
      x.p = TVector3(tsos.globalMomentum().x(), tsos.globalMomentum().y(), tsos.globalMomentum().z());
      x.dist = 9999.;
      fXpTracks.push_back(x); 
    }
  }

}


// ----------------------------------------------------------------------
bool HFDumpMuons::doExtrapolate(double pt, double eta) { 
  if (pt < 0.8) return false;
  if (eta > 2.4) return false;
  if (eta < 0.8 && pt < 3.4) return false;
  if (eta > 0.8 && eta < 1.4 && (pt < -eta + 4.2)) return false;
  if (eta > 1.4 && eta < 2.4 && (pt < -eta + 2.8)) return false;
  return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpMuons);
