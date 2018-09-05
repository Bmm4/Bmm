#include "Bmm/CmsswAnalysis/interface/HfrTrackListProducer.hh"

#include <memory>
#include "RecoParticleFlow/PFTracking/interface/PFTrackProducer.h"
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrackFwd.h"
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
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("generalTracks"))),
  fTrackQualityString(iConfig.getUntrackedParameter<string>("trackQualityString", string("highPurity"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fBeamSpotLabel(iConfig.getUntrackedParameter<InputTag>("BeamSpotLabel", InputTag("offlineBeamSpot"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel", InputTag("muons"))),
  fMuonQualityString(iConfig.getUntrackedParameter<string>("muonQuality", string("AllGlobalMuons"))),
  fTrackMinPt(iConfig.getUntrackedParameter<double>("trackMinPt", 3.0)) {
  produces<TrackCollection>();
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
  // -- magnetic field
  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  fMagneticField = fieldHandle.product();

  // -- primary vertices
  Handle<VertexCollection> hVertexCollection;
  try {
    iEvent.getByToken(fTokenVertex, hVertexCollection);
    hVertexCollection.isValid();
    fVertexCollection = *hVertexCollection;
  } catch(cms::Exception&){
    //    throw HFSetupException("No primary vertex collection found, skipping");
  }

  // -- beam spot
  Handle<BeamSpot> hBeamSpot;
  try {
    iEvent.getByToken(fTokenBeamSpot, hBeamSpot);
    hBeamSpot.isValid();
    fBeamSpot = *hBeamSpot;
  } catch(cms::Exception&){
    //    throw HFSetupException("No beam spot collection found, skipping");
  }

  // -- tracks
  try {
    iEvent.getByToken(fTokenTrack, fTracksHandle);
    fTracksHandle.isValid();
  } catch(cms::Exception&){
    //    throw HFSetupException(Form("No valid TrackCollection with label '%s' found, skipping", fTracksLabel.encode().c_str()));
  }

  // -- muons
  Handle<MuonCollection> hMuonCollection;
  try {
    iEvent.getByToken(fTokenMuon, hMuonCollection);
    hMuonCollection.isValid();
    fMuonCollection = hMuonCollection.product();
  } catch(cms::Exception&){
    //    throw HFSetupException(Form("No valid MuonCollection with label '%s' found, skipping",fMuonsLabel.encode().c_str()));
  }


  // -- create and fill the track collections
  auto trackColl = std::make_unique<reco::TrackCollection>();
  for (unsigned int ix = 0; ix < fTracksHandle->size(); ix++) {
    reco::TrackBaseRef rTrackView(fTracksHandle, ix);
    const reco::Track track(*rTrackView);
    if (!track.quality(reco::TrackBase::qualityByName(fTrackQualityString))) continue;
    if (track.pt() < fTrackMinPt) continue;
    trackColl->push_back(track);
  }



  iEvent.put(std::move(trackColl));
}

// ----------------------------------------------------------------------
void HfrTrackListProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HfrTrackListProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}
