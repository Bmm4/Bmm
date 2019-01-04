#include "Bmm/CmsswAnalysis/plugins/HfrJpsiXRecoilProducer.h"
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
    cout << "==>HfrJpsiXRecoilProducer> less than 2 global muons in event, not putting any track list into event " << endl
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
	cout << "==>HfrJpsiXRecoilProducer> charge correlation failed: "  << t1.charge() << " x " << t2.charge()
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
    cout << "==>HfrJpsiXRecoilProducer> charge correlation failed "  << endl;
    return;
  }

  if ((mbest < 2.5) || (mbest > 4.5)) {
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    cout << "==>HfrJpsiXRecoilProducer> psi mass outside of range: " << mbest << endl;
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
  cout << "==>HfrJpsiXRecoilProducer> J/psi mass = " << mbest << " from muons with trk idx = " << im1 << " and " << im2 << " and prob = " << vtxProb << endl;
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

    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), psiVertex->vertexState().position());
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), psiVertex->vertexState());
    if (doca.value() > 0.1) continue;
    trkIdxDoca.push_back(make_pair(ix, doca.value()));
  }
  sort(trkIdxDoca.begin(), trkIdxDoca.end(), docaLess);

  RefCountedKinematicTree psitTree;
  double psitProb(vtxProb);
  for (unsigned int ix = 0; ix < trkIdxDoca.size(); ++ix) {
    if (trkIdxDoca[ix].second > 0.025) {
      cout << "reached doca = " << trkIdxDoca[ix].second << ", breaking" << endl;
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
    if (prob < psitProb) {
      cout << "skipping vertex with probability " << prob << " because J/psi vtx had prob = " << psitProb << endl;
      vtt.pop_back();
      continue;
    }
    if (prob < fVtxProb) {
      cout << "skipping vertex with probability " << prob << " because fVtxProb = " << fVtxProb << endl;
      vtt.pop_back();
      continue;
    }

    psitProb = prob;
    cout << "keeping track " << trkIdxDoca[ix].first << " in vertex fit, with doca = " << trkIdxDoca[ix].second
	 << " pt = " << baseRef->pt()
	 << " prob = " << prob << " idx(pv) = " << psitPV.first
	 << endl;
  }


  vertexColl->push_back(mkVertex(psiTree));

  pvIdx = std::make_unique<int>(psiPV.first);

  trkIdxColl->push_back(im1);
  trkIdxColl->push_back(im2);


  // -- now try to add the closest tracks and check that
  //    - the vertex fit is good
  //    - whether the PV changed



  cout << "==>HfrJpsiXRecoilProducer> put into event tracklist with " << trkIdxColl->size()
       << " tracks, removed "  << ntracks - trkIdxColl->size()
       << " first two indices: " << trkIdxColl->at(0) << " " << trkIdxColl->at(1)
       << " J/psi PV idx: " << psiPV.first << " with |lip| = " << psiPV.second
       << endl;



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


DEFINE_FWK_MODULE(HfrJpsiXRecoilProducer);


  // KalmanVertexFitter fitter;
  // TransientVertex myVertex = fitter.vertex(vtt);
  // Vertex secVertex = myVertex;


  // AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  // int psiPV(-1);
  // double lip(999.), psiLip(9999.);
  // for (unsigned int iv = 0; iv < fVertexCollection.size(); ++iv) {
  //   Vertex currentPV = fVertexCollection[iv];

  //   TrajectoryStateOnSurface tsos1 = extrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),
  // 							      RecoVertex::convertPos(currentPV.position()));
  //   std::pair<bool,Measurement1D> currentIp = IPTools::signedDecayLength3D(tsos1, GlobalVector(0,0,1), currentPV);
  //   lip = currentIp.second.value();
  //   //    cout << "pv " << iv << " lip: " << lip << endl;
  //   if (TMath::Abs(lip) < psiLip) {
  //     psiLip = TMath::Abs(lip);
  //     psiPV = iv;
  //   }
  // }

  //  const reco::Vertex::Point& vpoint = secVertex.position();
  // reco::Vertex::Error verr = secVertex.error();
  // GlobalError err(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2));
  // GlobalPoint displacementFromBeamspot( -1.*((fBeamSpot.x0() -  secVertex.x()) +  (secVertex.z() - fBeamSpot.z0()) * fBeamSpot.dxdz()),
  // 					-1.*((fBeamSpot.y0() - secVertex.y())+  (secVertex.z() - fBeamSpot.z0()) * fBeamSpot.dydz()),
  // 					0.);
  // float flxy  = displacementFromBeamspot.perp();
  // float flxye = sqrt(err.rerr(displacementFromBeamspot));
  // double flsxy = flxy/flxye;
  // cout << "J/psi vtx = " << myVertex.position() << " with flxy/flxye = " << flxy << "/" << flxye << " = " << flsxy << endl;
