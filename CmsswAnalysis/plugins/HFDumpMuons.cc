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
  fpropM2(iConfig.getParameter<edm::ParameterSet>("propM2")),
  fOutwardPropM1(iConfig.getParameter<edm::ParameterSet>("OutwardPropM1")),
  fInwardPropM1(iConfig.getParameter<edm::ParameterSet>("InwardPropM1")),
  fweightFileBarrel(iConfig.getUntrackedParameter<edm::FileInPath>("weightFileBarrel")),
  fweightFileEndcap(iConfig.getUntrackedParameter<edm::FileInPath>("weightFileEndcap")) {
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
  fOutwardPropM1.init(iSetup);
  fInwardPropM1.init(iSetup);
}


// ----------------------------------------------------------------------
void HFDumpMuons::beginJob() {
  cout << "HFDumpMuons::beginJob" << endl;
  barrelBDT.setWeightFile(fweightFileBarrel);
  barrelBDT.setupReader();
  endcapBDT.setWeightFile(fweightFileEndcap);
  endcapBDT.setupReader();
}

// ----------------------------------------------------------------------
void HFDumpMuons::analyze(const Event& iEvent, const EventSetup& iSetup) {
  if (0) cout << "--- HFDumpMuons -------------------------------------------------------------------" << endl;

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
    pM->fNumberOfValidPixHits    = track_hp.numberOfValidPixelHits();
    pM->fInnerChi2               = iTrack->normalizedChi2();
  }


  if (gTrack.isNonnull()) {
    Track trk(*gTrack);
    pM->fGlobalPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
    pM->fNvalidMuonHits = gTrack->hitPattern().numberOfValidMuonHits();


    pM->fNhitsDT  = gTrack->hitPattern().numberOfValidMuonDTHits();
    pM->fNhitsCSC = gTrack->hitPattern().numberOfValidMuonCSCHits();
    pM->fNhitsRPC = gTrack->hitPattern().numberOfValidMuonRPCHits();

    fillMuonDetHits(pM,gTrack);
    pM->fGlbKinkFinder  = rm.combinedQuality().glbKink;
    pM->fTrkRelChi2     = rm.combinedQuality().trkRelChi2;
    pM->fStaRelChi2     = rm.combinedQuality().staRelChi2;
    pM->fGlbDeltaEtaPhi = rm.combinedQuality().globalDeltaEtaPhi;
    pM->fStaTrkMult     = getTrackMultiplicity(rm,true,false);
    pM->fTmTrkMult      = getTrackMultiplicity(rm,false,false);

    barrelBDT.getMuon()->fillBDTmuon(rm, &fVertexCollection, &fBeamSpot, pM->fStaTrkMult, pM->fTmTrkMult);
    endcapBDT.getMuon()->fillBDTmuon(rm, &fVertexCollection, &fBeamSpot, pM->fStaTrkMult, pM->fTmTrkMult);
    pM->fBarrelBDTresponse = barrelBDT.evaluate();
    pM->fEndcapBDTresponse = endcapBDT.evaluate();
  }

  if (iTrack.isNonnull()) {
    Track trk(*iTrack);
    pM->fInnerPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
    //    cout << " ipt = " << trk.pt() << " " <<  trk.eta() << " "  <<  trk.phi() << endl;
  }

  if (oTrack.isNonnull()) {
    Track trk(*oTrack);
    pM->fOuterPlab.SetPtEtaPhi(trk.pt(), trk.eta(), trk.phi());
    pM->fOuterChi2 = oTrack->normalizedChi2();
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
      if (0) cout << "===> M2:   "
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

// ----------------------------------------------------------------------
int HFDumpMuons::getTrackMultiplicity(const reco::Muon& Mu, bool MuType, bool verbose) {
  //MuType=true -->STA
  //MuType=false -->TrackerMuons
  int TrkMult(0);

  bool muonFlag = muon::isGoodMuon(Mu, muon::AllGlobalMuons);
  bool HPflag = false;
  if (muonFlag)
    {HPflag = (Mu.innerTrack())->quality(Track::highPurity);}
  if (!muonFlag || !HPflag)
    {return -1;}

  if (verbose)
    {
      cout << "HM: starting process looking for: " << MuType << " ("
	   << (fMuonCollection->end() - fMuonCollection->begin())
	   << " candidates)" << endl;
    }

  for (MuonCollection::const_iterator MuIt = fMuonCollection->begin(); MuIt != fMuonCollection->end();++MuIt)
    {
      if (verbose) {
	cout << "HM: Muon " << (MuIt - fMuonCollection->begin()) << endl;
	cout << "HM: muonType: " << (*MuIt).type() << endl;
      }
      if (MuType && !muon::isGoodMuon(*MuIt,muon::AllStandAloneMuons))
	{continue;}
      else if (!MuType && !muon::isGoodMuon(*MuIt,muon::AllTrackerMuons))
	{continue;}

      if (verbose) {cout << "HM: Searching in cone for " << MuType << endl;}

      if (MuType && (*MuIt).outerTrack().isNull())
	{
	  if (verbose) {cout << "HM: Invalid STA." << endl;}
	  continue;
	}

      if (MuType)
	{if ( tracksAreEqual(Mu.outerTrack(),(*MuIt).outerTrack()) ) {continue;}}
      else
	{
	  //to avoid double counting with TM is tricky, as GM and TM are produced by different algorithms.
	}

      double max_distance(0);
      //optimal max_distance determined using a BDT
      if (MuType) {max_distance = 150;}
      else {max_distance = 100;}

      double distance = getDistanceM1(Mu,(*MuIt), MuType,verbose);
      if ( (distance < max_distance) && (distance>=0) )
	{TrkMult++;}
      if (verbose) {cout << "HM: Moving on to next muon." << endl;}
    }

  return TrkMult;
}

double HFDumpMuons::getDistanceM1(const reco::Muon& muon, const reco::Muon& sec, bool MuType, bool verbose) {

  if (verbose) {cout << "HM: Filling the HMhitInfo." << endl;}
  if (muon.innerTrack().isNull())
    {if (verbose) {cout << "HM: GM track is not valid." << endl;}return -1;}

  TrackRef track;
  if (MuType) {track = sec.outerTrack();}
  else {track = sec.innerTrack();}
  if (track.isNull())
    {
      if (verbose) {cout << "HM: Invalid (secondary) track." << endl;}
      return -1;
    }

  TrajectoryStateOnSurface main1TSOS = fOutwardPropM1.extrapolate(*muon.globalTrack());
  TrajectoryStateOnSurface sec1TSOS;
  if (MuType)
    {
      sec1TSOS = fInwardPropM1.extrapolate(*track);
      if (verbose) {cout << "HM: Created TSOS using outer track." << endl;}
    }
  else
    {
      sec1TSOS = fpropM1.extrapolate(*track);
      if (verbose) {cout << "HM: Created TSOS using inner track." << endl;}
    }

  if ( !main1TSOS.isValid() || !sec1TSOS.isValid() )
    {
      if (verbose) {cout << "HM: At least one TSOS is not valid." << endl;}
      return -1;
    }
  TVector3 Vmain1;
  TVector3 Vsec1;
  Vmain1.SetXYZ(main1TSOS.globalPosition().x(),main1TSOS.globalPosition().y(),main1TSOS.globalPosition().z());
  Vsec1.SetXYZ(sec1TSOS.globalPosition().x(),sec1TSOS.globalPosition().y(),sec1TSOS.globalPosition().z());

  if (verbose) {cout << "HM: 3D dist: " << (Vmain1-Vsec1).Mag() << endl;}

  return (Vmain1-Vsec1).Mag();
}

bool HFDumpMuons::tracksAreEqual(const TrackRef& main,const TrackRef& STA) {
  double margin = 0.01;

  if (main.isNull())
    {return false;}
  if (STA.isNull())
    {return false;}

  double dpt = TMath::Abs(main->pt() - STA->pt());
  double deta = TMath::Abs(main->eta() - STA->eta());
  double dphi = TMath::Abs(main->phi() - STA->phi());
  if (dpt < margin && deta < margin && dphi < margin)
    {return true;}
  return false;
}

void HFDumpMuons::fillMuonDetHits(TAnaMuon* pM, TrackRef gTrack) {

  //Fill valid hits per muon detector and per station from HitPattern
  //The invalid hits are not filled. It is not very productive using the global track
  //as they are filtered out during the global muon reconstruction
  //Nevertheless, the hits are checked for validity.
  const reco::HitPattern &GMpattern = gTrack->hitPattern();
  pM->fvDThits = {0,0,0,0};
  pM->fvRPChits = {0,0,0,0};
  pM->fvCSChits = {0,0,0,0};

  for (int i=0;i<gTrack->hitPattern().numberOfHits(HitPattern::TRACK_HITS);i++)
    {
      uint32_t hit = GMpattern.getHitPattern(HitPattern::TRACK_HITS, i);
      if ( !GMpattern.validHitFilter(hit) )
	{continue;}

      if (GMpattern.getMuonStation(hit) == 1)
	{
	  if (GMpattern.muonDTHitFilter(hit)) pM->fvDThits[0]++;
	  if (GMpattern.muonRPCHitFilter(hit)) pM->fvRPChits[0]++;
	  if (GMpattern.muonCSCHitFilter(hit)) pM->fvCSChits[0]++;
	}
      else if (GMpattern.getMuonStation(hit) == 2)
	{
	  if (GMpattern.muonDTHitFilter(hit)) pM->fvDThits[1]++;
	  if (GMpattern.muonRPCHitFilter(hit)) pM->fvRPChits[1]++;
	  if (GMpattern.muonCSCHitFilter(hit)) pM->fvCSChits[1]++;
	}
      else if (GMpattern.getMuonStation(hit) == 3)
	{
	  if (GMpattern.muonDTHitFilter(hit)) pM->fvDThits[2]++;
	  if (GMpattern.muonRPCHitFilter(hit)) pM->fvRPChits[2]++;
	  if (GMpattern.muonCSCHitFilter(hit)) pM->fvCSChits[2]++;
	}
      else if (GMpattern.getMuonStation(hit) == 4)
	{
	  if (GMpattern.muonDTHitFilter(hit)) pM->fvDThits[3]++;
	  if (GMpattern.muonRPCHitFilter(hit)) pM->fvRPChits[3]++;
	  if (GMpattern.muonCSCHitFilter(hit)) pM->fvCSChits[3]++;
	}
    }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpMuons);
