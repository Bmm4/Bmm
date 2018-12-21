#include "Bmm/CmsswAnalysis/plugins/HfrJpsiXListProducer.h"
#include "Bmm/RootAnalysis/common/HFMasses.hh"

#include <memory>
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include <TLorentzVector.h>

using namespace std;
using namespace edm;
using namespace reco;

// ----------------------------------------------------------------------
HfrJpsiXListProducer::HfrJpsiXListProducer(const ParameterSet& iConfig) :
  HfrBaseProducer(iConfig),
  fVtxProb(iConfig.getUntrackedParameter<double>("vtxProb", 0.01)) {
  produces<vector<int> >(); // vector of track indices
  produces<int>();          // possibly a PV index
  produces<VertexCollection>();
}


// ----------------------------------------------------------------------
void HfrJpsiXListProducer::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HfrBaseProducer::dumpConfiguration()" << endl;
  HfrBaseProducer::dumpConfiguration();
}

// ----------------------------------------------------------------------
void HfrJpsiXListProducer::put(edm::Event &iEvent,
			       std::unique_ptr<std::vector<int> > &trkIdxColl,
			       std::unique_ptr<int> &pvIdx,
			       std::unique_ptr<VertexCollection> &vertexColl) {
  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
  iEvent.put(std::move(vertexColl));
}

// ----------------------------------------------------------------------
void HfrJpsiXListProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  auto trkIdxColl = std::make_unique<vector<int> >();
  auto pvIdx = std::make_unique<int>(1234);
  std::unique_ptr<VertexCollection> vertexColl(new VertexCollection());

  if (fVerbose > 0)  cout << " ======================================================================"
			  << endl
			  << "=== HfrJpsiXListProducer run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event() << endl
			  << " ----------------------------------------------------------------------"
			  << endl;

  HfrBaseProducer::analyze(iEvent, iSetup);

  // -- find J/psi
  vector<int> midx;
  reco::MuonCollection::const_iterator muonIt;
  muon::SelectionType muonType = muon::selectionTypeFromString(fMuonQualityString);
  for (muonIt = fMuonCollection->begin(); muonIt != fMuonCollection->end(); ++muonIt) {
    if (!muon::isGoodMuon(*muonIt, muonType)) continue;
    int ixMu = muonIt->track().index();
    if (ixMu < 0) continue;
    reco::TrackBaseRef rTrackView(fTracksHandle, ixMu);
    const reco::Track track(*rTrackView);
    if (track.pt() < fMuonMinPt) continue;
    midx.push_back(ixMu);
  }

  if (midx.size() < 2) {
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    return;
  }
  // -- at this point, simply use the best J/psi based on the mass.
  TLorentzVector p4m1, p4m2, p4psi, m1best, m2best;
  double mbest(-999.);
  int im1(-1), im2(-1);
  for (unsigned int im = 0; im < midx.size(); ++im) {
    for (unsigned int in = im; in < midx.size(); ++in) {
      reco::TrackBaseRef m1(fTracksHandle, midx[im]);
      reco::Track t1(*m1);
      reco::TrackBaseRef m2(fTracksHandle, midx[in]);
      reco::Track t2(*m2);
      if (t1.charge()*t2.charge() > 0) continue;
      p4m1.SetXYZM(t1.px(), t1.py(), t1.pz(), MMUON);
      p4m2.SetXYZM(t2.px(), t2.py(), t2.pz(), MMUON);
      p4psi = p4m1 + p4m2;
      if (TMath::Abs(p4psi.M() - MJPSI) < TMath::Abs(mbest - MJPSI)) {
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

  // -- setup vertex fit
  vector<TransientTrack>vtt;
  TrackBaseRef rm1(fTracksHandle, im1);
  vtt.push_back(fTTB->build(*rm1));
  TrackBaseRef rm2(fTracksHandle, im2);
  vtt.push_back(fTTB->build(*rm2));

  KalmanVertexFitter fitter;
  TransientVertex myVertex = fitter.vertex(vtt);
  Vertex secVertex = myVertex;

  // -- Vtx probability and displacement from beamspot(!) copied from HLTrigger/btau/src/HLTDisplacedmumuFilter.cc
  double vtxProb = 0.0;
  if ((secVertex.chi2() >= 0.0) && (secVertex.ndof() > 0)) vtxProb = TMath::Prob(secVertex.chi2(), secVertex.ndof());
  if (vtxProb < fVtxProb) {
    // cout << "skipping vertex with probability " << vtxProb << " and mass = " << mbest << endl;
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    return;
  }

  if ((mbest < 2.5) || (mbest > 4.5)) {
    put(iEvent, trkIdxColl, pvIdx, vertexColl);
    return;
  }
  //  cout << "J/psi mass = " << mbest << " from muons with trk idx = " << im1 << " and " << im2 << endl;


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
  trkIdxColl->push_back(im1);
  trkIdxColl->push_back(im2);

  // -- create and fill the track collections
  TrajectoryStateOnSurface tsos;
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
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
  cout << "==>HfrJpsiXListProducer> put into event tracklist with " << trkIdxColl->size()
       << " tracks, removed "  << ntracks - trkIdxColl->size()
       << endl;


  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
  iEvent.put(std::move(vertexColl));
}

// ----------------------------------------------------------------------
void HfrJpsiXListProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HfrJpsiXListProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}

DEFINE_FWK_MODULE(HfrJpsiXListProducer);
