#include "Bmm/CmsswAnalysis/plugins/HFJpsiXListProducer.h"
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

using namespace std;
using namespace edm;
using namespace reco;

// ----------------------------------------------------------------------
HFJpsiXListProducer::HFJpsiXListProducer(const ParameterSet& iConfig) :
  HfrBaseProducer(iConfig),
  fVtxProb(iConfig.getUntrackedParameter<double>("vtxProb", 0.01)) {
  dumpConfiguration();
  produces<vector<int> >(); // vector of track indices
  produces<int>();          // possibly a PV index
  produces<VertexCollection>();
}


// ----------------------------------------------------------------------
void HFJpsiXListProducer::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFJpsiXListProducer::dumpConfiguration()" << endl;
  HfrBaseProducer::dumpConfiguration();
  cout << "----------------------------------------------------------------------" << endl;
}

// ----------------------------------------------------------------------
void HFJpsiXListProducer::put(edm::Event &iEvent,
			       std::unique_ptr<std::vector<int> > &trkIdxColl,
			       std::unique_ptr<int> &pvIdx,
			       std::unique_ptr<VertexCollection> &vertexColl) {
  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
  iEvent.put(std::move(vertexColl));
}

// ----------------------------------------------------------------------
void HFJpsiXListProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  auto trkIdxColl = std::make_unique<vector<int> >();
  auto pvIdx = std::make_unique<int>(1234);
  std::unique_ptr<VertexCollection> vertexColl(new VertexCollection());

  if (fVerbose > 0)  cout << "======================================================================"
			  << endl
			  << "=== HFJpsiXListProducer run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event() << endl
			  << "----------------------------------------------------------------------"
			  << endl;

  HfrBaseProducer::analyze(iEvent, iSetup);

  // -- find J/psi
  vector<int> midx;
  muon::SelectionType muonType = muon::selectionTypeFromString(fMuonQualityString);
  reco::MuonCollection::const_iterator muonIt;
  for (muonIt = fMuonCollection->begin(); muonIt != fMuonCollection->end(); ++muonIt) {
    if (!muon::isGoodMuon(*muonIt, muonType)) {
      //      cout << "==>HFJpsiXListProducer> not a global muon " << endl;
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
    cout << "==>HFJpsiXListProducer> less than 2 global muons in event, not putting any track list into event " << endl
	 << "   midx.size() = " << midx.size() << endl
	 << "   fMuonCollection->size() = " << fMuonCollection->size()
	 << endl;
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    return;
  }

  // -- at this point, simply use the best J/psi based on the mass.
  TLorentzVector p4m1, p4m2, p4psi, m1best, m2best;
  double mbest(-999.);
  int im1(-1), im2(-1);
  for (unsigned int im = 0; im < midx.size(); ++im) {
    for (unsigned int in = im+1; in < midx.size(); ++in) {
      reco::TrackBaseRef m1(fTracksHandle, midx[im]);
      reco::Track t1(*m1);
      reco::TrackBaseRef m2(fTracksHandle, midx[in]);
      reco::Track t2(*m2);
      if (t1.charge()*t2.charge() > 0) continue;
      p4m1.SetXYZM(t1.px(), t1.py(), t1.pz(), MMUON);
      p4m2.SetXYZM(t2.px(), t2.py(), t2.pz(), MMUON);
      p4psi = p4m1 + p4m2;
      if (TMath::Abs(p4psi.M() - MJPSI) < TMath::Abs(mbest - MJPSI)) {
	cout << "HFJpsiXListProducer> replacing  m = " << mbest << " with " << p4psi.M()
	     << " from muon idxs = " << midx[im] << " and " << midx[in] << endl;
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
    return;
  }

  if ((mbest < 2.5) || (mbest > 4.5)) {
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    return;
  }
  cout << "HFJpsiXListProducer> J/psi mass = " << mbest << " from muons with trk idx = " << im1 << " and " << im2 << endl;

  // -- setup vertex fit
  vector<TransientTrack>vtt;
  TrackBaseRef rm1(fTracksHandle, im1);
  vtt.push_back(fTTB->build(*rm1));
  TrackBaseRef rm2(fTracksHandle, im2);
  vtt.push_back(fTTB->build(*rm2));

  KalmanVertexFitter fitter;
  TransientVertex myVertex = fitter.vertex(vtt);
  Vertex secVertex = myVertex;

  // -- find best PV
  RefCountedKinematicTree kinTree = fitTree(vtt);
  kinTree->movePointerToTheTop();
  RefCountedKinematicParticle kinParticle = kinTree->currentParticle();

  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  int bestPV(-1);
  double lip(999.), bestLip(9999.);
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

  // -- Vtx probability and displacement from beamspot(!) copied from HLTrigger/btau/src/HLTDisplacedmumuFilter.cc
  double vtxProb = 0.0;
  if ((secVertex.chi2() >= 0.0) && (secVertex.ndof() > 0)) vtxProb = TMath::Prob(secVertex.chi2(), secVertex.ndof());
  if (vtxProb < fVtxProb) {
    // cout << "skipping vertex with probability " << vtxProb << " and mass = " << mbest << endl;
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    return;
  }

  //  const reco::Vertex::Point& vpoint = secVertex.position();
  reco::Vertex::Error verr = secVertex.error();
  GlobalError err(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2));
  GlobalPoint displacementFromBeamspot( -1.*((fBeamSpot.x0() -  secVertex.x()) +  (secVertex.z() - fBeamSpot.z0()) * fBeamSpot.dxdz()),
					-1.*((fBeamSpot.y0() - secVertex.y())+  (secVertex.z() - fBeamSpot.z0()) * fBeamSpot.dydz()),
					0.);

  // float flxy  = displacementFromBeamspot.perp();
  // float flxye = sqrt(err.rerr(displacementFromBeamspot));
  // double flsxy = flxy/flxye;
  // cout << "J/psi vtx = " << myVertex.position() << " with flxy/flxye = " << flxy << "/" << flxye << " = " << flsxy << endl;

  vertexColl->push_back(secVertex);

  pvIdx = std::make_unique<int>(bestPV);
  trkIdxColl->push_back(im1);
  trkIdxColl->push_back(im2);

  // -- create and fill the track collections
  TrajectoryStateOnSurface tsos;
  VertexDistance3D a3d;
  int ntracks(0);
  for (unsigned int ix = 0; ix < fTracksHandle->size(); ix++) {
    // -- skip muons, already added above to be the first in the list
    if (ix == static_cast<unsigned int>(im1)) continue;
    if (ix == static_cast<unsigned int>(im2)) continue;
    reco::TrackBaseRef rTrackView(fTracksHandle, ix);
    const reco::Track track(*rTrackView);
    TransientTrack transTrack = fTTB->build(*rTrackView);
    ++ntracks;
    if (!track.quality(reco::TrackBase::qualityByName(fTrackQualityString))) continue;
    if (track.pt() < fTrackMinPt) continue;

    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), myVertex.position());
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), myVertex);
    if (doca.value() > 0.1) continue;

    trkIdxColl->push_back(ix);
    //    cout << " adding track idx = " << ix << " with pT = " << track.pt() << " and doca = " << doca.value() << endl;
  }
  cout << "==>HFJpsiXListProducer> put into event tracklist with " << trkIdxColl->size()
       << " tracks, removed "  << ntracks - trkIdxColl->size()
       << " first two indices: " << trkIdxColl->at(0) << " " << trkIdxColl->at(1)
       << " best PV idx: " << bestPV << " at z = " << fVertexCollection[bestPV].z()
       << endl;


  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
  iEvent.put(std::move(vertexColl));
}

// ----------------------------------------------------------------------
void HFJpsiXListProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HFJpsiXListProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}

// ----------------------------------------------------------------------
RefCountedKinematicTree HFJpsiXListProducer::fitTree(vector<TransientTrack> &vtt) {
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


DEFINE_FWK_MODULE(HFJpsiXListProducer);
