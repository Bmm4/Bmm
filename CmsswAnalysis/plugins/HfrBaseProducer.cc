#include "Bmm/CmsswAnalysis/plugins/HfrBaseProducer.h"

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
HfrBaseProducer::HfrBaseProducer(const ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fType(iConfig.getUntrackedParameter<int>("type", 0)),
  fTracksLabel(iConfig.getUntrackedParameter<InputTag>("tracksLabel", InputTag("generalTracks"))),
  fTrackQualityString(iConfig.getUntrackedParameter<string>("trackQualityString", string("highPurity"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<InputTag>("PrimaryVertexLabel", InputTag("offlinePrimaryVertices"))),
  fBeamSpotLabel(iConfig.getUntrackedParameter<InputTag>("BeamSpotLabel", InputTag("offlineBeamSpot"))),
  fMuonsLabel(iConfig.getUntrackedParameter<InputTag>("muonsLabel", InputTag("muons"))),
  fMuonQualityString(iConfig.getUntrackedParameter<string>("muonQuality", string("AllGlobalMuons"))),
  fTrackMinPt(iConfig.getUntrackedParameter<double>("trackMinPt", 0.4)),
  fMuonMinPt(iConfig.getUntrackedParameter<double>("muonMinPt", 4.0)),
  fDocaMaxPv(iConfig.getUntrackedParameter<double>("docaMaxPv", 0.1)),
  fDocaMaxTk(iConfig.getUntrackedParameter<double>("docaMaxTk", 0.1)) {

  // produces<vector<int> >(); // vector of track indices
  // produces<int>();          // possibly a PV index

  fTokenBeamSpot    = consumes<BeamSpot>(fBeamSpotLabel);
  fTokenTrack       = consumes<edm::View<reco::Track> >(fTracksLabel) ;
  fTokenMuon       = consumes<MuonCollection>(fMuonsLabel) ;
  fTokenVertex      = consumes<VertexCollection>(fPrimaryVertexLabel);

}


// ----------------------------------------------------------------------
void HfrBaseProducer::dumpConfiguration() {
  cout << "---  verbose                     " << fVerbose << endl;
  cout << "---  type                        " << fType << endl;
  cout << "---  tracksLabel                 " << fTracksLabel << endl;
  cout << "---  trackQualityString          " << fTrackQualityString << endl;
  cout << "---  PrimaryVertexLabel          " << fPrimaryVertexLabel << endl;
  cout << "---  BeamSpotLabel               " << fBeamSpotLabel << endl;
  cout << "---  muonsLabel                  " << fMuonsLabel << endl;
  cout << "---  muonQuality                 " << fMuonQualityString << endl;
  cout << "---  trackMinPt                  " << fTrackMinPt << endl;
  cout << "---  muonMinPt                   " << fMuonMinPt << endl;
  cout << "---  docaMaxPv                   " << fDocaMaxPv << endl;
  cout << "---  docaMaxTk                   " << fDocaMaxTk << endl;
}



// ----------------------------------------------------------------------
// -- this is called from all derived classes
void HfrBaseProducer::analyze(Event& iEvent, const EventSetup& iSetup) {
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
    cout << "ERROR: did not get track collection from event!" << endl;
    //    throw HFSetupException(Form("No valid TrackCollection with label '%s' found, skipping", fTracksLabel.encode().c_str()));
  }

  // -- load the transient track builder
  edm::ESHandle<TransientTrackBuilder> ttbHandle;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttbHandle);
  if (!ttbHandle.isValid()) {
    cout << "Error: no TransientTrackBuilder handle found" << endl;
    return;
  }
  fTTB = ttbHandle.product();

  // -- muons
  Handle<MuonCollection> hMuonCollection;
  try {
    iEvent.getByToken(fTokenMuon, hMuonCollection);
    hMuonCollection.isValid();
    fMuonCollection = hMuonCollection.product();
  } catch(cms::Exception&){
    //    throw HFSetupException(Form("No valid MuonCollection with label '%s' found, skipping",fMuonsLabel.encode().c_str()));
  }

}


// ----------------------------------------------------------------------
// -- this is NOT called from any derived class!
void HfrBaseProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  if (fVerbose > 0)  cout << " ======================================================================"
			  << "=== HfrBaseProducer run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event() << endl
			  << " ----------------------------------------------------------------------"
			  << endl;
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
    cout << "ERROR: did not get track collection from event!" << endl;
    //    throw HFSetupException(Form("No valid TrackCollection with label '%s' found, skipping", fTracksLabel.encode().c_str()));
  }

  // -- load the transient track builder
  edm::ESHandle<TransientTrackBuilder> ttbHandle;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttbHandle);
  if (!ttbHandle.isValid()) {
    cout << "Error: no TransientTrackBuilder handle found" << endl;
    return;
  }
  fTTB = ttbHandle.product();

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
  cout << "==>HfrBaseProducer> put into event tracklist with " << trkIdxColl->size()
       << " tracks, removed "  << ntracks - trkIdxColl->size()
       << endl;

  iEvent.put(std::move(trkIdxColl));

  auto pvIdx = std::make_unique<int>(1234);
  iEvent.put(std::move(pvIdx));
}

// ----------------------------------------------------------------------
void HfrBaseProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HfrBaseProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}

DEFINE_FWK_MODULE(HfrBaseProducer);
