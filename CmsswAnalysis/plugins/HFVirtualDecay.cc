// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFVirtualDecay
// --------------
//
// 2016/01/20 Urs Langenegger      migrate to "consumes"
// 2013/01/31 Christoph Naegeli    first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "HFVirtualDecay.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/Common/interface/Handle.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <TString.h>

using namespace std;
using namespace edm;
using namespace reco;

// ----------------------------------------------------------------------
HFVirtualDecay::HFVirtualDecay(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  fFilterLabel(iConfig.getUntrackedParameter<edm::InputTag>("filterLabel", edm::InputTag(""))),
  fTracksLabel(iConfig.getUntrackedParameter<edm::InputTag>("tracksLabel", edm::InputTag("generalTracks"))),
  fTrackQualityString(iConfig.getUntrackedParameter<std::string>("trackQualityString",std::string("highPurity"))),
  fPrimaryVertexLabel(iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertexLabel", edm::InputTag("offlinePrimaryVertices"))),
  fBeamSpotLabel(iConfig.getUntrackedParameter<edm::InputTag>("BeamSpotLabel", edm::InputTag("offlineBeamSpot"))),
  fMuonsLabel(iConfig.getUntrackedParameter<edm::InputTag>("muonsLabel", edm::InputTag("muons"))),
  fMuonQualityString(iConfig.getUntrackedParameter<std::string>("muonQualityString",std::string("AllGlobalMuons"))),
  fTrackPt(iConfig.getUntrackedParameter<double>("trackPt", 3.0)),
  fMuonPt(iConfig.getUntrackedParameter<double>("muonPt", 1.0)),
  fChi2(iConfig.getUntrackedParameter<double>("chi2", 100.)),
  fCandPt(iConfig.getUntrackedParameter<double>("candpt", -1.)),
  fCandLo(iConfig.getUntrackedParameter<double>("candlo", -1.)),
  fCandHi(iConfig.getUntrackedParameter<double>("candhi", -1.)),
  fDimuonPt(iConfig.getUntrackedParameter<double>("dimuonpt", -1.)),
  fPvIpS(iConfig.getUntrackedParameter<double>("pvips", 99999.0)),
  fFlsxy(iConfig.getUntrackedParameter<double>("flsxy", 3.0)),
  fFlxy(iConfig.getUntrackedParameter<double>("flxy", 9999.0)),
  fFls3d(iConfig.getUntrackedParameter<double>("fls3d", -1.0)),
  fMassE(iConfig.getUntrackedParameter<double>("masserr", 9999.0)),
  fDeltaR(iConfig.getUntrackedParameter<double>("deltaR", 99.0)),
  fMaxDoca(iConfig.getUntrackedParameter<double>("maxDoca", 0.2)),
  fMaxD0(iConfig.getUntrackedParameter<double>("maxD0", 999.)),
  fMaxDz(iConfig.getUntrackedParameter<double>("maxDz", 999.)),
  fPvWeight(iConfig.getUntrackedParameter<double>("pvWeight", 0.6)),
  fType(iConfig.getUntrackedParameter<int>("type")),
  fUseBeamspotConstraint(iConfig.getUntrackedParameter<bool>("useBeamspotConstraint", true)),
  fDoFilter(false) {

  fTokenBeamSpot    = consumes<BeamSpot>(fBeamSpotLabel);
  fTokenTrack       = consumes<edm::View<reco::Track> >(fTracksLabel) ;
  fTokenMuon       = consumes<MuonCollection>(fMuonsLabel) ;
  fTokenVertex      = consumes<VertexCollection>(fPrimaryVertexLabel);
  if (fFilterLabel.encode() != "") {
    fDoFilter = true;
    fTokenTrkIdx      = consumes<vector<int> >(fFilterLabel);
    fTokenPvIdx       = consumes<int>(fFilterLabel) ;
  }
}


// ----------------------------------------------------------------------
void HFVirtualDecay::dumpConfiguration() {
  using namespace std;
  cout << "---  verbose                     " << fVerbose << endl;
  cout << "---  filterLabel                 " << fFilterLabel << endl;
  cout << "---  tracksLabel                 " << fTracksLabel << endl;
  cout << "---  trackQualityString          " << fTrackQualityString << endl;
  cout << "---  PrimaryVertexLabel          " << fPrimaryVertexLabel << endl;
  cout << "---  BeamSpotLabel               " << fBeamSpotLabel << endl;
  cout << "---  muonsLabel                  " << fMuonsLabel << endl;
  cout << "---  muonQuality                 " << fMuonQualityString << endl;

  cout << "---  trackPt                     " << fTrackPt << endl;
  cout << "---  muonPt                      " << fMuonPt << endl;
  cout << "---  chi2                        " << fChi2 << endl;
  cout << "---  candpt                      " << fCandPt << endl;
  cout << "---  candlo                      " << fCandLo << endl;
  cout << "---  candhi                      " << fCandHi << endl;
  cout << "---  dimuonpt                    " << fDimuonPt << endl;
  cout << "---  pvips                       " << fPvIpS << endl;
  cout << "---  flsxy                       " << fFlsxy << endl;
  cout << "---  flxy                        " << fFlxy << endl;
  cout << "---  fls3d                       " << fFls3d << endl;
  cout << "---  masserr                     " << fMassE << endl;

  cout << "---  deltaR                      " << fDeltaR << endl;
  cout << "---  maxDoca                     " << fMaxDoca << endl;
  cout << "---  maxD0                       " << fMaxD0 << endl;
  cout << "---  maxDz                       " << fMaxDz << endl;
  cout << "---  pvWeight                    " << fPvWeight << endl;

  cout << "---  type                        " <<  fType << endl;
  cout << "---  useBeamSpotConstraint       " <<  fUseBeamspotConstraint << endl;
}

// ----------------------------------------------------------------------
void HFVirtualDecay::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using std::cout; using std::endl;

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
    if (fVertexCollection.size() == 0) throw HFSetupException("Primary vertex collection is empty, skipping");
  } catch(cms::Exception&){
    throw HFSetupException("No primary vertex collection found, skipping");
  }

  // -- beam spot
  Handle<BeamSpot> hBeamSpot;
  try {
    iEvent.getByToken(fTokenBeamSpot, hBeamSpot);
    hBeamSpot.isValid();
    fBeamSpot = *hBeamSpot;
  } catch(cms::Exception&){
    throw HFSetupException("No beam spot collection found, skipping");
  }

  // -- tracks
  try {
    iEvent.getByToken(fTokenTrack, fTracksHandle);
    fTracksHandle.isValid();
  } catch(cms::Exception&){
    throw HFSetupException(Form("No valid TrackCollection with label '%s' found, skipping", fTracksLabel.encode().c_str()));
  }

  // -- muons
  Handle<MuonCollection> hMuonCollection;
  try {
    iEvent.getByToken(fTokenMuon, hMuonCollection);
    hMuonCollection.isValid();
    fMuonCollection = hMuonCollection.product();
  } catch(cms::Exception&){
    throw HFSetupException(Form("No valid MuonCollection with label '%s' found, skipping",fMuonsLabel.encode().c_str()));
  }

  // -- setup filtering for recoil physics
  int pvidx(0);
  vector<int> vidx;
  if (fDoFilter) {
    // -- recoil track indices
    try {
      iEvent.getByToken(fTokenTrkIdx, fTrkIdxHandle);
      fTrkIdxHandle.isValid();
    } catch(cms::Exception&){
      throw HFSetupException(Form("No valid TrkIdxCollection with label '%s' found, skipping", fFilterLabel.encode().c_str()));
    }

    // -- recoil pv index
    try {
      iEvent.getByToken(fTokenPvIdx, fPvIdxHandle);
      fPvIdxHandle.isValid();
    } catch(cms::Exception&){
      throw HFSetupException(Form("No valid PvIndex with label '%s' found, skipping", fFilterLabel.encode().c_str()));
    }
    pvidx = (*fPvIdxHandle);
    vidx  = (*fTrkIdxHandle);
  }

  // -- load the transient track builder
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", fTTB);
  if (!fTTB.isValid()) throw HFSetupException("Error: no TransientTrackBuilder found");

  // -- construct a list builder
  fListBuilder.reset(new  HFTrackListBuilder(fTracksHandle, fMuonCollection, fTTB.product(), fVerbose));
  fListBuilder->setMaxD0(fMaxD0);
  fListBuilder->setMaxDz(fMaxDz);
  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setTrackQuality(fTrackQualityString);

  if (fDoFilter) {
    fListBuilder->setRecoilTrkIdx(vidx);
    fListBuilder->setRecoilPvIdx(pvidx);
    if (fVerbose > 0) {
      cout << "HFVirtualDecay: recoil trk list size = " << vidx.size() << " pv idx = " << pvidx << endl;
    }
  }

  muon::SelectionType muonType = muon::selectionTypeFromString(fMuonQualityString);
  fListBuilder->setMuonQuality(muonType);

  // -- construct the sequential vertex fitter
  fSequentialFitter.reset(new HFSequentialVertexFit(fTracksHandle, fMuonCollection, fTTB.product(), hVertexCollection, fMagneticField, fBeamSpot, fVerbose));
  fSequentialFitter->setPvW8(fPvWeight);
  fSequentialFitter->useBeamspotConstraint(fUseBeamspotConstraint);
}
