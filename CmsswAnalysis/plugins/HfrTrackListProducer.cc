#include "Bmm/CmsswAnalysis/plugins/HfrTrackListProducer.h"

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

using namespace std;
using namespace edm;
using namespace reco;

// ----------------------------------------------------------------------
HfrTrackListProducer::HfrTrackListProducer(const ParameterSet& iConfig) :
  HfrBaseProducer(iConfig) {
  produces<vector<int> >(); // vector of track indices
  produces<int>();          // possibly a PV index
}


// ----------------------------------------------------------------------
void HfrTrackListProducer::dumpConfiguration() {
  using namespace std;
  cout << "---  verbose                     " << fVerbose << endl;
  cout << "---  tracksLabel                 " << fTracksLabel << endl;
  cout << "---  trackQualityString          " << fTrackQualityString << endl;
  cout << "---  PrimaryVertexLabel          " << fPrimaryVertexLabel << endl;
  cout << "---  BeamSpotLabel               " << fBeamSpotLabel << endl;
  cout << "---  muonsLabel                  " << fMuonsLabel << endl;
  cout << "---  muonQuality                 " << fMuonQualityString << endl;
  cout << "---  trackMinPt                  " << fTrackMinPt << endl;
}



// ----------------------------------------------------------------------
void HfrTrackListProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  if (fVerbose > 0)  cout << " ======================================================================"
			  << "=== HfrTrackListProducer run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event() << endl
			  << " ----------------------------------------------------------------------"
			  << endl;

  HfrBaseProducer::analyze(iEvent, iSetup);

  // -- create and fill the track collections
  auto trkIdxColl = std::make_unique<vector<int> >();
  int ntracks(0);
  for (unsigned int ix = 0; ix < fTracksHandle->size(); ix++) {
    reco::TrackBaseRef rTrackView(fTracksHandle, ix);
    const reco::Track track(*rTrackView);
    ++ntracks;
    if (!track.quality(reco::TrackBase::qualityByName(fTrackQualityString))) continue;
    if (track.pt() < fTrackMinPt) continue;
    trkIdxColl->push_back(ix);
    cout << " adding track idx = " << ix << " with pT = " << track.pt() << " to hfrtracks list" << endl;
  }
  cout << "==>HfrTrackListProducer> put into event tracklist with " << trkIdxColl->size()
       << " tracks, removed "  << ntracks - trkIdxColl->size()
       << endl;

  iEvent.put(std::move(trkIdxColl));

  auto pvIdx = std::make_unique<int>(1234);
  iEvent.put(std::move(pvIdx));
}

// ----------------------------------------------------------------------
void HfrTrackListProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HfrTrackListProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}

DEFINE_FWK_MODULE(HfrTrackListProducer);
