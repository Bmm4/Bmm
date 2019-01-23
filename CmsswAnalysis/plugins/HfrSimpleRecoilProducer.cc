#include "Bmm/CmsswAnalysis/plugins/HfrSimpleRecoilProducer.h"
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

#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include <TLorentzVector.h>

using namespace std;
using namespace edm;
using namespace reco;


// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HfrSimpleRecoilProducer::HfrSimpleRecoilProducer(const ParameterSet& iConfig) :
  HfrBaseProducer(iConfig),
  fVtxProb(iConfig.getUntrackedParameter<double>("vtxProb", 0.01)),
  fCandType(iConfig.getUntrackedParameter<int>("candtype", 0)) {
  dumpConfiguration();
  produces<vector<int> >(); // vector of track indices
  produces<int>();          // possibly a PV index
}


// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HfrSimpleRecoilProducer::dumpConfiguration()" << endl;
  HfrBaseProducer::dumpConfiguration();
  cout << "----------------------------------------------------------------------" << endl;
}

// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::put(edm::Event &iEvent,
				  std::unique_ptr<std::vector<int> > &trkIdxColl,
				  std::unique_ptr<int> &pvIdx) {
  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
}

// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  auto trkIdxColl = std::make_unique<vector<int> >();
  auto pvIdx = std::make_unique<int>(1234);

  if (fVerbose > 0)  cout << " ======================================================================"
			  << endl
			  << "=== HfrSimpleRecoilProducer run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event() << endl
			  << " ----------------------------------------------------------------------"
			  << endl;

  HfrBaseProducer::analyze(iEvent, iSetup);

  // -- get candidate and PV from gHFEvent
  TAnaCand *pCand(0);
  TAnaTrack *pS(0);
  int ipv(-1);
  vector<int> brecoIdx;
  for (int ic = 0; ic < gHFEvent->nCands(); ++ic) {
    pCand = gHFEvent->getCand(ic);
    if (fCandType == pCand->fType) {
      ipv = pCand->fPvIdx;
      cout << "found one " << pCand->fType << " with PV idx = " << ipv << endl;
      for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
	pS = gHFEvent->getSigTrack(i);
	pS->dump();
	brecoIdx.push_back(pS->fIndex);
      }
    }
  }
  if (ipv < 0) {
      put(iEvent, trkIdxColl, pvIdx);
      return;
  }
  pvIdx = std::make_unique<int>(ipv);

  Vertex myVertex = fVertexCollection[ipv];
  GlobalPoint pVertex(myVertex.position().x(), myVertex.position().y(), myVertex.position().z());

  // -- create and fill the track collections
  TrajectoryStateOnSurface tsos;
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  VertexDistance3D a3d;
  int ntracks(0), cnt(0);
  for (unsigned int ix = 0; ix < fTracksHandle->size(); ix++) {
    // -- skip tracks from BRECO
    for (unsigned int ib = 0; ib < brecoIdx.size(); ++ib) {
      if (static_cast<int>(ix) == brecoIdx[ib]) continue;
    }
    reco::TrackBaseRef rTrackView(fTracksHandle, ix);
    const reco::Track track(*rTrackView);
    TransientTrack transTrack = fTTB->build(*rTrackView);
    ++ntracks;
    if (!track.quality(reco::TrackBase::qualityByName(fTrackQualityString))) continue;
    if (track.pt() < fTrackMinPt) continue;

    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), pVertex);
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), myVertex);
    if (doca.value() > 0.1) continue;

    trkIdxColl->push_back(ix);
    ++cnt;
  }
  cout << "==>HfrSimpleRecoilProducer> put into event tracklist with " << trkIdxColl->size()
       << " tracks, removed "  << ntracks - trkIdxColl->size()
       << endl;

  // -- put into event
  put(iEvent, trkIdxColl, pvIdx);
}

// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}

DEFINE_FWK_MODULE(HfrSimpleRecoilProducer);
