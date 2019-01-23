#include "Bmm/CmsswAnalysis/plugins/HfrJpsiXRecoilProducer.h"
#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"

#include <memory>
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"

#include <TLorentzVector.h>


#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"

using namespace std;
using namespace edm;
using namespace reco;

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
// Utility routine to sort the mindoca array of the candidate...
static bool docaLess(pair<int, double> x, pair<int,double> y) {
  return x.second < y.second;
}


// ----------------------------------------------------------------------
HfrJpsiXRecoilProducer::HfrJpsiXRecoilProducer(const ParameterSet& iConfig) :
  HfrBaseProducer(iConfig),
  fVtxProb(iConfig.getUntrackedParameter<double>("vtxProb", 0.01)) {
  dumpConfiguration();
  produces<vector<int> >(); // vector of track indices
  produces<int>();          // possibly a PV index
  produces<VertexCollection>();
}


// ----------------------------------------------------------------------
void HfrJpsiXRecoilProducer::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HfrJpsiXRecoilProducer::dumpConfiguration()" << endl;
  HfrBaseProducer::dumpConfiguration();
  cout << "----------------------------------------------------------------------" << endl;
}

// ----------------------------------------------------------------------
void HfrJpsiXRecoilProducer::put(edm::Event &iEvent,
			       std::unique_ptr<std::vector<int> > &trkIdxColl,
			       std::unique_ptr<int> &pvIdx,
			       std::unique_ptr<VertexCollection> &vertexColl) {
  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
  iEvent.put(std::move(vertexColl));
}

// ----------------------------------------------------------------------
void HfrJpsiXRecoilProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  auto trkIdxColl = std::make_unique<vector<int> >();
  auto pvIdx = std::make_unique<int>(1234);
  std::unique_ptr<VertexCollection> vertexColl(new VertexCollection());

  if (fVerbose > 0)  cout << "======================================================================"
			  << endl
			  << "=== HfrJpsiXRecoilProducer run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event() << endl
			  << "----------------------------------------------------------------------"
			  << endl;

  HfrBaseProducer::analyze(iEvent, iSetup);

  // -- find all muons
  vector<int> midx;
  muon::SelectionType muonType = muon::selectionTypeFromString(fMuonQualityString);
  reco::MuonCollection::const_iterator muonIt;
  for (muonIt = fMuonCollection->begin(); muonIt != fMuonCollection->end(); ++muonIt) {
    if (!muon::isGoodMuon(*muonIt, muonType)) {
      //      cout << "==>HfrJpsiXRecoilProducer> not a global muon " << endl;
      continue;
    }
    int ixMu = muonIt->track().index();
    if (ixMu < 0) continue;
    reco::TrackBaseRef rTrackView(fTracksHandle, ixMu);
    const reco::Track track(*rTrackView);
    if (track.pt() < fMuonMinPt) continue;
    midx.push_back(ixMu);
  }

  if (midx.size() < 2) {
    if (0) cout << "==>HfrJpsiXRecoilProducer> less than 2 global muons in event, not putting any track list into event " << endl
		<< "   midx.size() = " << midx.size() << endl
		<< "   fMuonCollection->size() = " << fMuonCollection->size()
		<< endl;
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    return;
  }

  // -- simply use the best J/psi based on the mass
  TLorentzVector p4m1, p4m2, p4psi, m1best, m2best;
  double mbest(-999.);
  int im1(-1), im2(-1);
  for (unsigned int im = 0; im < midx.size(); ++im) {
    for (unsigned int in = im+1; in < midx.size(); ++in) {
      reco::TrackBaseRef m1(fTracksHandle, midx[im]);
      reco::Track t1(*m1);
      reco::TrackBaseRef m2(fTracksHandle, midx[in]);
      reco::Track t2(*m2);
      if (t1.charge()*t2.charge() > 0) {
	if (0) cout << "==>HfrJpsiXRecoilProducer> charge correlation failed: "  << t1.charge() << " x " << t2.charge()
		    << " for muon idxs = " << midx[im] << " and " << midx[in]
		    << endl;
	continue;
      }
      p4m1.SetXYZM(t1.px(), t1.py(), t1.pz(), MMUON);
      p4m2.SetXYZM(t2.px(), t2.py(), t2.pz(), MMUON);
      p4psi = p4m1 + p4m2;
      if (TMath::Abs(p4psi.M() - MJPSI) < TMath::Abs(mbest - MJPSI)) {
	// cout << "HfrJpsiXRecoilProducer> replacing  m = " << mbest << " with " << p4psi.M()
	//      << " from muon idxs = " << midx[im] << " and " << midx[in] << endl;
	mbest = p4psi.M();
	m1best = p4m1;
	m2best = p4m2;
	im1 = midx[im];
	im2 = midx[in];
      }
    }
  }

  if ((im1 < 0) || (im2 < 0)) {
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    if (0) cout << "==>HfrJpsiXRecoilProducer> charge correlation failed "  << endl;
    return;
  }

  if ((mbest < 2.5) || (mbest > 4.5)) {
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    if (0) cout << "==>HfrJpsiXRecoilProducer> psi mass outside of range: " << mbest << endl;
    return;
  }

  // -- setup vertex fit
  vector<TransientTrack>vtt;
  TrackBaseRef rm1(fTracksHandle, im1);
  vtt.push_back(fTTB->build(*rm1));
  TrackBaseRef rm2(fTracksHandle, im2);
  vtt.push_back(fTTB->build(*rm2));

  // -- find best PV with J/psi only
  RefCountedKinematicTree psiTree;
  pair<int, double> psiPV = findBestPV(vtt, psiTree);
  RefCountedKinematicVertex   psiVertex = psiTree->currentDecayVertex();
  Vertex                       myVertex = mkVertex(psiTree);


  // -- Vtx probability and displacement from beamspot(!) copied from HLTrigger/btau/src/HLTDisplacedmumuFilter.cc
  double vtxProb = 0.0;
  if ((psiVertex->chiSquared() >= 0.0) && (psiVertex->degreesOfFreedom() > 0)) {
    vtxProb = TMath::Prob(psiVertex->chiSquared(), psiVertex->degreesOfFreedom());
  }
  if (0) cout << "==>HfrJpsiXRecoilProducer> J/psi mass = " << mbest << " from muons with trk idx = " << im1 << " and " << im2
	      << " and prob = " << vtxProb
	      << " PV(idx, |lip|) = " << psiPV.first << ", " << psiPV.second
	      << endl;
  if (vtxProb < fVtxProb) {
    // cout << "skipping vertex with probability " << vtxProb << " and mass = " << mbest << endl;
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    return;
  }

  // -- create vector<int, double> with hadronic tracks, sorted by doca to dimuon vertex
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TrajectoryStateOnSurface tsos;
  VertexDistance3D a3d;
  vector<pair<int, double> > trkIdxDoca;
  for (unsigned int ix = 0; ix < fTracksHandle->size(); ix++) {
    // -- skip muons, already added above to be the first in the list
    if (ix == static_cast<unsigned int>(im1)) continue;
    if (ix == static_cast<unsigned int>(im2)) continue;
    reco::TrackBaseRef rTrackView(fTracksHandle, ix);
    const reco::Track track(*rTrackView);
    TransientTrack transTrack = fTTB->build(*rTrackView);
    if (!track.quality(reco::TrackBase::qualityByName(fTrackQualityString))) continue;
    if (track.pt() < fTrackMinPt) continue;

    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), psiVertex->vertexState().position());
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), psiVertex->vertexState());
    if (doca.value() > 0.1) continue;
    trkIdxDoca.push_back(make_pair(ix, doca.value()));
  }
  sort(trkIdxDoca.begin(), trkIdxDoca.end(), docaLess);

  // -- create track index list with J/psi muon indices and all qualifying hadronic tracks.
  RefCountedKinematicTree psitTree;
  double psitProb(vtxProb);
  double psitlip(psiPV.second);
  vector<int> trkList;
  trkList.push_back(im1);
  trkList.push_back(im2);
  for (unsigned int ix = 0; ix < trkIdxDoca.size(); ++ix) {
    if (trkIdxDoca[ix].second > 0.020) {
      if (0) cout << "reached doca = " << trkIdxDoca[ix].second << ", breaking at track " << trkIdxDoca[ix].first << " at position " << ix << endl;
      break;
    }
    reco::TrackBaseRef baseRef(fTracksHandle, trkIdxDoca[ix].first);
    vtt.push_back(fTTB->build(*baseRef));
    pair<int, double> psitPV = findBestPV(vtt, psitTree);
    RefCountedKinematicVertex psitVertex = psitTree->currentDecayVertex();

    double prob = 0.0;
    if ((psitVertex->chiSquared() >= 0.0) && (psitVertex->degreesOfFreedom() > 0)) {
      prob = TMath::Prob(psitVertex->chiSquared(), psitVertex->degreesOfFreedom());
    }
    if (0 && prob < psitProb) {
      cout << "skipping track " << trkIdxDoca[ix].first << " vertex with probability " << prob << " because J/psi vtx had prob = " << psitProb << endl;
      vtt.pop_back();
      continue;
    }
    if (0 && prob < fVtxProb) {
      cout << "skipping track " << trkIdxDoca[ix].first << " vertex with probability " << prob << " because fVtxProb = " << fVtxProb << endl;
      vtt.pop_back();
      continue;
    }

    if (psitPV.first == psiPV.first && psitPV.second > psitlip) {
      if (0) cout << "skipping track " << trkIdxDoca[ix].first << " |lip| increased to " << psitPV.second << " compared to = " << psitlip << endl;
      vtt.pop_back();
      continue;
    }

    if (psitPV.first != psiPV.first) {
      if (0) cout << "skipping track " << trkIdxDoca[ix].first << " PV changed to " << psitPV.first << " compared to = " << psiPV.first << endl;
      vtt.pop_back();
      continue;
    }


    psitlip  = psitPV.second;
    psitProb = prob;

    trkList.push_back(trkIdxDoca[ix].first);

    if (0) cout << "keeping track " << trkIdxDoca[ix].first << " in vertex fit, with doca = " << trkIdxDoca[ix].second
		<< " pt = " << baseRef->pt()
		<< " prob = " << prob << " idx(pv) = " << psitPV.first << " |lip| = " << psitPV.second
		<< endl;
  }

  Vertex mypVertex = fVertexCollection[psiPV.first];
  GlobalPoint pVertex(mypVertex.position().x(), mypVertex.position().y(), mypVertex.position().z());

  // ----------------------------------------------------------------------
  // -- store JpsiX cand in gHFEvent (cf. HFListCand)
  // ----------------------------------------------------------------------
  TAnaCand *pCand = fillCand(trkList, psiPV.first);
  pCand->fPvIdx = psiPV.first;
  pCand->fPvLip = psiPV.second;

  // ----------------------------------------------------------------------
  // -- create list of recoil tracks (cf. HfrSimpleRecoilProducer)
  // ----------------------------------------------------------------------
  vector<int> recoilList, missList;
  bool skip(false);
  for (unsigned int ix = 0; ix < fTracksHandle->size(); ix++) {
    skip = false;
    // -- skip tracks from J/psi X
    if (trkList.end() != find(trkList.begin(), trkList.end(), static_cast<int>(ix))) {
      continue;
    }
    TrackBaseRef rTrackView(fTracksHandle, ix);
    const Track track(*rTrackView);
    TransientTrack transTrack = fTTB->build(*rTrackView);
    if (!track.quality(reco::TrackBase::qualityByName(fTrackQualityString))) {
      skip = true;
    }
    if (track.pt() < fTrackMinPt) {
      skip = true;
    }

    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), pVertex);
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), mypVertex);
    if (doca.value() > fDocaMaxPv) {
      skip = true;
    }

    if (false == skip) {
      recoilList.push_back(ix);
      trkIdxColl->push_back(ix);
    } else {
      int trkPvIdx = getPv(ix, &fVertexCollection);
      // miss list contains high-purity tracks from nearby
      if ((doca.value() < 1.0) && (track.quality(reco::TrackBase::qualityByName(fTrackQualityString)))
	  ) {
	missList.push_back(static_cast<int>(ix));
      }
    }
  }

  pCand = fillCand(recoilList, psiPV.first);
  pCand->fType  = fType+1; // override the value set
  pCand->fPvIdx = psiPV.first;
  pCand->fPvLip = psiPV.second;

  pCand = fillCand(missList, psiPV.first);
  pCand->fType  = fType+2; // override the value set
  pCand->fPvIdx = psiPV.first;
  pCand->fPvLip = psiPV.second;

  pvIdx = std::make_unique<int>(psiPV.first);

  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
  iEvent.put(std::move(vertexColl));
}

// ----------------------------------------------------------------------
void HfrJpsiXRecoilProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HfrJpsiXRecoilProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}

// ----------------------------------------------------------------------
RefCountedKinematicTree HfrJpsiXRecoilProducer::fitTree(std::vector<reco::TransientTrack> &vtt) {
  // -- fit vertex
  KinematicParticleFactoryFromTransientTrack factory;
  ParticleMass muon_mass = 0.1056583;
  float muon_sigma = 0.0000000001;
  float chi = 0.;
  float ndf = 0.;
  vector<RefCountedKinematicParticle> particles;
  for(vector<TransientTrack>::const_iterator i = vtt.begin();  i != vtt.end(); ++i) {
    particles.push_back(factory.particle(*i, muon_mass, chi, ndf, muon_sigma));
  }
  KinematicParticleVertexFitter fitter;
  RefCountedKinematicTree kinTree = fitter.fit(particles);
  kinTree->movePointerToTheTop();
  return kinTree;
}

// ----------------------------------------------------------------------
pair<int, double> HfrJpsiXRecoilProducer::findBestPV(vector<TransientTrack> &vtt, RefCountedKinematicTree &kinTree) {
  // -- fit vertex
  KinematicParticleFactoryFromTransientTrack factory;
  ParticleMass muon_mass = 0.1056583;
  float muon_sigma = 0.0000000001;
  float chi = 0.;
  float ndf = 0.;
  vector<RefCountedKinematicParticle> particles;
  for(vector<TransientTrack>::const_iterator i = vtt.begin();  i != vtt.end(); ++i) {
    particles.push_back(factory.particle(*i, muon_mass, chi, ndf, muon_sigma));
  }
  KinematicParticleVertexFitter fitter;
  kinTree = fitter.fit(particles);
  kinTree->movePointerToTheTop();

  // -- setup extrapolators
  RefCountedKinematicParticle kinParticle = kinTree->currentParticle();
  RefCountedKinematicVertex   kinVertex = kinTree->currentDecayVertex();
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  int bestPV(-1);
  double lip(999.), bestLip(9999.);
  // -- loop over PVs and determine best one
  for (unsigned int iv = 0; iv < fVertexCollection.size(); ++iv) {
    Vertex currentPV = fVertexCollection[iv];

    TrajectoryStateOnSurface tsos1 = extrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),
							      RecoVertex::convertPos(currentPV.position()));
    std::pair<bool,Measurement1D> currentIp = IPTools::signedDecayLength3D(tsos1, GlobalVector(0,0,1), currentPV);
    lip = currentIp.second.value();
    //    cout << "pv " << iv << " lip: " << lip << endl;
    if (TMath::Abs(lip) < bestLip) {
      bestLip = TMath::Abs(lip);
      bestPV = iv;
    }
  }
  return make_pair(bestPV, bestLip);
}


// ----------------------------------------------------------------------
Vertex HfrJpsiXRecoilProducer::mkVertex(RefCountedKinematicTree &kinTree) {
  RefCountedKinematicVertex kinVertex = kinTree->currentDecayVertex();
  return Vertex(RecoVertex::convertPos(kinVertex->vertexState().position()),
		kinVertex->vertexState().error().matrix(),
		kinVertex->chiSquared(),
		kinVertex->degreesOfFreedom(),
		kinTree->daughterParticles().size());
}


// ----------------------------------------------------------------------
TAnaCand* HfrJpsiXRecoilProducer::fillCand(vector<int> trkList, int pvIdx) {
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fType = fType;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + trkList.size() - 1;
  TAnaTrack *pTrack(0);

  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TrajectoryStateOnSurface tsos;

  Vertex mypVertex = fVertexCollection[pvIdx];
  GlobalPoint pVertex(mypVertex.position().x(), mypVertex.position().y(), mypVertex.position().z());
  VertexDistance3D a3d;

  TLorentzVector p4t, p4tot;
  p4tot.SetXYZT(0., 0., 0., 0.);
  for (vector<int>::const_iterator trkIt = trkList.begin(); trkIt != trkList.end(); ++trkIt) {
    int idx = *trkIt;
    //    cout << "idx = " << ix << " -> track idx = " << idx << " pv idx = " << pvidx << endl;
    TSimpleTrack *sTrack = gHFEvent->getSimpleTrack(idx);
    int gidx = sTrack->getGenIndex();
    pTrack = gHFEvent->addSigTrack();

    TrackBaseRef baseRef(fTracksHandle, idx);
    const Track trackView(*baseRef);
    const BeamSpot *pBeamSpot = &fBeamSpot;

    TransientTrack transTrack = fTTB->build(*baseRef);
    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), pVertex);
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), mypVertex);

    fillAnaTrack(pTrack, trackView, idx, gidx, &fVertexCollection, fMuonCollection, pBeamSpot);
    pTrack->fDouble1 = doca.value();

    p4t.SetXYZM(trackView.px(), trackView.py(), trackView.pz(), MPION);
    p4tot += p4t;
  }
  pCand->fMass  = p4tot.M();
  pCand->fPlab  = p4tot.Vect();
  return pCand;
}



DEFINE_FWK_MODULE(HfrJpsiXRecoilProducer);
