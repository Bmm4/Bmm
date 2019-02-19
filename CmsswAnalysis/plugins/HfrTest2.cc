#include "FWCore/Framework/interface/MakerMacros.h"
#include "HfrTest2.h"

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


#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"

using namespace std;
using namespace edm;
using namespace reco;

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HfrTest2::HfrTest2(const edm::ParameterSet& iConfig) :
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 1)),
  fTracksLabel(iConfig.getUntrackedParameter<edm::InputTag>("tracksLabel", edm::InputTag("generalTracks"))),
  fPVLabel(iConfig.getUntrackedParameter<edm::InputTag>("pvLabel", edm::InputTag("offlinePrimaryVertices"))),
  fBeamSpotLabel(iConfig.getUntrackedParameter<edm::InputTag>("BeamSpotLabel", edm::InputTag("offlineBeamSpot"))),
  fMuonsLabel(iConfig.getUntrackedParameter<edm::InputTag>("muonsLabel", edm::InputTag("muons"))),
  fMuonQualityString(iConfig.getUntrackedParameter<std::string>("muonQualityString",std::string("AllGlobalMuons"))),
  fTrkIdxLabel(iConfig.getUntrackedParameter<edm::InputTag>("trkIdxLabel", edm::InputTag("hfrtracks"))),
  fPvIdxLabel(iConfig.getUntrackedParameter<edm::InputTag>("pvIdxLabel", edm::InputTag("hfrtracks"))),
  fVertexLabel(iConfig.getUntrackedParameter<edm::InputTag>("vertexLabel", edm::InputTag("hfrtracks"))),
  fLeadingTrackPt(iConfig.getUntrackedParameter<double>("leadingTrackPt", -1.0)) {
  dumpConfiguration();
  fTokenTrack       = consumes<edm::View<reco::Track> >(fTracksLabel) ;
  fTokenPV          = consumes<VertexCollection>(fPVLabel);
  fTokenBeamSpot    = consumes<BeamSpot>(fBeamSpotLabel);
  fTokenMuon        = consumes<MuonCollection>(fMuonsLabel) ;
  fTokenTrkIdx      = consumes<vector<int> >(fTrkIdxLabel);
  fTokenPvIdx       = consumes<int>(fPvIdxLabel) ;
}


// ----------------------------------------------------------------------
void HfrTest2::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HfrTest2 constructor" << endl;
  cout << "---  muonsLabel                  " << fMuonsLabel << endl;
  cout << "---  muonQuality                 " << fMuonQualityString << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HfrTest2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (fVerbose > 0)  cout << "=== HfrTest2 run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event()
			  << " ==================================================================="
			  << endl;


  // -- tracks
  try {
    iEvent.getByToken(fTokenTrack, fTracksHandle);
    fTracksHandle.isValid();
  } catch(cms::Exception&){
    throw HFSetupException(Form("No valid TrackCollection with label '%s' found, skipping", fTracksLabel.encode().c_str()));
  }

  // -- load the transient track builder
  edm::ESHandle<TransientTrackBuilder> ttbHandle;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttbHandle);
  if (!ttbHandle.isValid()) {
    cout << "Error: no TransientTrackBuilder handle found" << endl;
    return;
  }
  fTTB = ttbHandle.product();

  // -- magnetic field
  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  fMagneticField = fieldHandle.product();

  // -- primary vertices
  Handle<VertexCollection> hPVCollection;
  try {
    iEvent.getByToken(fTokenPV, hPVCollection);
    hPVCollection.isValid();
    fPVCollection = hPVCollection.product();
    if (fPVCollection->size() == 0) throw HFSetupException("Primary vertex collection is empty, skipping");
  } catch(cms::Exception&){
    throw HFSetupException("No primary vertex collection found, skipping");
  }

  // -- beam spot
  Handle<BeamSpot> hBeamSpot;
  try {
    iEvent.getByToken(fTokenBeamSpot, hBeamSpot);
    hBeamSpot.isValid();
    fBeamSpot = hBeamSpot.product();
  } catch(cms::Exception&){
    throw HFSetupException("No beam spot collection found, skipping");
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


  // -- recoil track indices
  try {
    iEvent.getByToken(fTokenTrkIdx, fTrkIdxHandle);
    fTrkIdxHandle.isValid();
  } catch(cms::Exception&){
    throw HFSetupException(Form("No valid TrkIdxCollection with label '%s' found, skipping", fTrkIdxLabel.encode().c_str()));
  }

  // -- recoil pv index
  try {
    iEvent.getByToken(fTokenPvIdx, fPvIdxHandle);
    fPvIdxHandle.isValid();
  } catch(cms::Exception&){
    throw HFSetupException(Form("No valid PvIndex with label '%s' found, skipping", fTrkIdxLabel.encode().c_str()));
  }

  //  int pvidx = (*fPvIdxHandle);
  vector<int> vidx = (*fTrkIdxHandle);

  // -- must have at least two muons (filled in first two places)
  if (vidx.size() < 2) return;

  // -- vertex the two muons once again
  vector<TransientTrack>vtt;
  TrackBaseRef rm1(fTracksHandle, vidx[0]); vtt.push_back(fTTB->build(*rm1));
  TrackBaseRef rm2(fTracksHandle, vidx[1]); vtt.push_back(fTTB->build(*rm2));

  KalmanVertexFitter fitter;
  TransientVertex myVertex = fitter.vertex(vtt);
  Vertex secVertex = myVertex;

  // -- Vtx probability and displacement from beamspot(!) copied from HLTrigger/btau/src/HLTDisplacedmumuFilter.cc
  double vtxProb = 0.0;
  if ((secVertex.chi2() >= 0.0) && (secVertex.ndof() > 0)) vtxProb = TMath::Prob(secVertex.chi2(), secVertex.ndof());

  reco::Vertex::Error verr = secVertex.error();
  GlobalError err(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2));
  GlobalPoint displacementFromBeamspot( -1.*((fBeamSpot->x0() -  secVertex.x()) +  (secVertex.z() - fBeamSpot->z0()) * fBeamSpot->dxdz()),
					-1.*((fBeamSpot->y0() - secVertex.y())+  (secVertex.z() - fBeamSpot->z0()) * fBeamSpot->dydz()),
					0.);

  float flxy  = displacementFromBeamspot.perp();
  float flxye = sqrt(err.rerr(displacementFromBeamspot));

  // double flsxy = flxy/flxye;
  //  cout << "J/psi vtx = " << myVertex.position() << " with flxy/flxye = " << flxy << "/" << flxye << " = " << flsxy << endl;

  TrajectoryStateOnSurface tsos;
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  VertexDistance3D a3d;
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fType = 1;
  pCand->fVtx.fDxy = flxy;
  pCand->fVtx.fDxyE = flxye;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + vidx.size() - 1;
  TAnaTrack *pTrack(0);

  // -- calculate dimuon mass and momentum
  TLorentzVector p4m1, p4m2, p4psi;
  TrackBaseRef m1Ref(fTracksHandle, vidx[0]);
  Track t1(*m1Ref);
  TrackBaseRef m2Ref(fTracksHandle, vidx[1]);
  Track t2(*m2Ref);
  p4m1.SetXYZM(t1.px(), t1.py(), t1.pz(), MMUON);
  p4m2.SetXYZM(t2.px(), t2.py(), t2.pz(), MMUON);
  p4psi = p4m1 + p4m2;

  pCand->fMass = p4psi.M();
  pCand->fPlab = p4psi.Vect();

  // -- loop over track from produced list
  TLorentzVector p4t, p4tot;
  for (unsigned int ix = 0; ix < vidx.size(); ix++) {
    int idx = vidx[ix];
    //    cout << "idx = " << ix << " -> track idx = " << idx << " pv idx = " << pvidx << endl;
    TrackBaseRef baseRef(fTracksHandle, idx);
    TransientTrack transTrack = fTTB->build(*baseRef);
    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), myVertex.position());
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), myVertex);

    TSimpleTrack *sTrack = gHFEvent->getSimpleTrack(idx);
    pTrack = gHFEvent->addSigTrack();
    int gidx = sTrack->getGenIndex();
    Track trackView(*baseRef);
    fillAnaTrack(pTrack, trackView, idx, gidx, fPVCollection, fMuonCollection, fBeamSpot);
    pTrack->fDouble1 = doca.value();

    p4t.SetXYZM(trackView.px(), trackView.py(), trackView.pz(), MKAON);
    p4tot = p4m1 + p4m2 + p4t;
    pTrack->fBsTip = p4tot.M();

    // -- check compatibility in another vertex fit
    if (ix  > 1) {
      vtt.push_back(fTTB->build(*baseRef));
      TransientVertex tVertex = fitter.vertex(vtt);
      Vertex tmpVertex = tVertex;

      double prob(0.);
      if ((tmpVertex.chi2() >= 0.0) && (tmpVertex.ndof() > 0)) prob = TMath::Prob(tmpVertex.chi2(), tmpVertex.ndof());
      pTrack->fDouble2 = vtxProb;
      pTrack->fDouble3 = prob;
      pTrack->fRefChi2 = tmpVertex.chi2();
      pTrack->fRefDof  = tmpVertex.ndof();

      // -- remove temporary track from vtt
      vtt.pop_back();
    } else {
      pTrack->fDouble2 = vtxProb;
      pTrack->fDouble3 = vtxProb;
      pTrack->fRefChi2 = -1.;
      pTrack->fRefDof  = -1;
    }
  }
  cout << "HfrTest2 cand dump: " << endl;
  pCand->dump();
  // for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
  //   gHFEvent->getSigTrack(is)->dump();
  // }
}


// ----------------------------------------------------------------------
int HfrTest2::idFromMass(double mass) {
  if (TMath::Abs(mass - MELECTRON) < 0.0001) return 11;
  if (TMath::Abs(mass - MMUON) < 0.0001) return 13;
  if (TMath::Abs(mass - MPION) < 0.0001) return 211;
  if (TMath::Abs(mass - MKAON) < 0.0001) return 321;
  if (TMath::Abs(mass - MPROTON) < 0.0001) return 2212;
  cout << "%%%> HfrTest2: mass " << mass << " not associated to any stable particle" << endl;
  return 0;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HfrTest2);
