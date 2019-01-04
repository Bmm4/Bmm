#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFJpsiXCand.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"

using namespace std;
using namespace edm;
using namespace reco;

// -- Yikes!
extern TAna01Event *gHFEvent;

enum { dimension = 5 };
typedef math::Error<dimension>::type CovarianceMatrix;

// ----------------------------------------------------------------------
// Utility routine to sort the mindoca array of the candidate...
static bool docaLess(pair<int, double> x, pair<int,double> y) {
  return x.second < y.second;
}

// ----------------------------------------------------------------------
HFJpsiXCand::HFJpsiXCand(const edm::ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fRefType(iConfig.getUntrackedParameter<int>("refType", 3000068)) {
  dumpConfiguration();
}


// ----------------------------------------------------------------------
void HFJpsiXCand::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFJpsiXCand constructor" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  refType                     " << fRefType << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFJpsiXCand::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (fVerbose > 0)  cout << "=== HFJpsiXCand run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event()
			  << " ==================================================================="
			  << endl;
  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFJpsiXCand::analyze> " << e.fMsg << endl;
    return;
  }


  // -- get candidate and PV from gHFEvent
  TAnaCand *tCand(0);
  TAnaTrack *pS(0);
  int ipv(-1);
  vector<int> brecoIdx;
  cout << "----------------------------------------------------------------------" << endl;
  for (int ic = 0; ic < gHFEvent->nCands(); ++ic) {
    tCand = gHFEvent->getCand(ic);
    if (fRefType == tCand->fType) {
      ipv = tCand->fPvIdx;
      cout << "found one with PV idx = " << ipv << " at z = " << gHFEvent->getPV(ipv)->fPoint.Z()
	   << (((*fPvIdxHandle) == ipv && ((*fPvIdxHandle) != 1234))? " ": " XXX ")
	   << endl;

      for (int i = tCand->fSig1; i <= tCand->fSig2; ++i) {
	pS = gHFEvent->getSigTrack(i);
	pS->dump();
      }
    }
  }
  cout << "---------------------------" << endl;


  fListBuilder->setCallerName("HFJpsiXCand");
  fListBuilder->setMinPt(fMuonPt);
  fListBuilder->setMaxDocaToTracks(fMaxDoca);

  cout << "==>HFJpsiXCand> fListBuilder->getMuonList()" << endl;
  vector<int> midx = fListBuilder->getMuonList();
  cout << "==>HFJpsiXCand> #muons = " << midx.size() << endl;
  if (midx.size() < 2) return;

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
	cout << "HFJpsiXCand> replacing  m = " << mbest << " with " << p4psi.M()
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
    return;
  }

  cout << "==>HFJpsiXCand> J/psi muons at " << im1 << " and " << im2 << endl;

  vector<TransientTrack>vtt;
  TrackBaseRef rm1(fTracksHandle, im1); vtt.push_back(fTTB->build(*rm1));
  TrackBaseRef rm2(fTracksHandle, im2); vtt.push_back(fTTB->build(*rm2));

  KalmanVertexFitter fitter;

  // -- TEST
  if (0) {
    TwoTrackMinimumDistance md;
    reco::Track t1(*rm1);
    reco::Track t2(*rm2);
    CovarianceMatrix cm1 = t1.covariance();
    CovarianceMatrix cm2 = t2.covariance();
    //    cout << cm << endl;
    //    cm*= 5.;
    //    cout << cm << endl;
    Track tmpt1(t1.chi2(), t1.ndof(), t1.referencePoint(), t1.momentum(), t1.charge(), cm1);
    Track tmpt2(t2.chi2(), t2.ndof(), t2.referencePoint(), t2.momentum(), t2.charge(), cm2);
    // CovarianceMatrix cmt = tmpt.covariance();
    // cout << cmt << endl;
    TransientTrack tt1  = fTTB->build(*rm1);
    TransientTrack tt1a = fTTB->build(tmpt1);
    TransientTrack tt2  = fTTB->build(*rm2);
    TransientTrack tt2a = fTTB->build(tmpt2);
    // -- calculate with 'original' tracks
    vector<TransientTrack> vtt;
    vtt.push_back(tt1);
    vtt.push_back(tt2);
    TransientVertex myVertex = fitter.vertex(vtt);
    AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
    TrajectoryStateOnSurface tsos1 = extrapolator.extrapolate(tt1.initialFreeState(), myVertex.position());
    TrajectoryStateOnSurface tsos2 = extrapolator.extrapolate(tt2.initialFreeState(), myVertex.position());
    md.calculate(tsos1, tsos2);
    cout << "XXXXXXXXX doca: " << md.distance() << endl;
    // -- calculate with 'your own' tracks
    vector<TransientTrack> vttt;
    vttt.push_back(tt1a);
    vttt.push_back(tt2a);
    TransientVertex myVertext = fitter.vertex(vttt);
    TrajectoryStateOnSurface tsos1t = extrapolator.extrapolate(tt1a.initialFreeState(), myVertext.position());
    TrajectoryStateOnSurface tsos2t = extrapolator.extrapolate(tt2a.initialFreeState(), myVertext.position());
    md.calculate(tsos1t, tsos2t);
    cout << "XXXXXXXXX doca: " << md.distance() << endl;
  }

  TransientVertex myVertex = fitter.vertex(vtt);
  Vertex secVertex = myVertex;

  // -- Vtx probability and displacement from beamspot(!) copied from HLTrigger/btau/src/HLTDisplacedmumuFilter.cc
  // double vtxProb = 0.0;
  // if ((secVertex.chi2() >= 0.0) && (secVertex.ndof() > 0)) vtxProb = TMath::Prob(secVertex.chi2(), secVertex.ndof());

  reco::Vertex::Error verr = secVertex.error();
  GlobalError err(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2));
  GlobalPoint displacementFromBeamspot(-1.*((fBeamSpot.x0() - secVertex.x()) + (secVertex.z() - fBeamSpot.z0())*fBeamSpot.dxdz()),
				       -1.*((fBeamSpot.y0() - secVertex.y()) + (secVertex.z() - fBeamSpot.z0())*fBeamSpot.dydz()),
				       0.);

  float flxy  = displacementFromBeamspot.perp();
  float flxye = sqrt(err.rerr(displacementFromBeamspot));

  // double flsxy = flxy/flxye;
  //  cout << "J/psi vtx = " << myVertex.position() << " with flxy/flxye = " << flxy << "/" << flxye << " = " << flsxy << endl;


  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setMaxDocaToTracks(fMaxDoca);

  cout << "==>HFJpsiXCand> fListBuilder->getTrackList()" << endl;
  vector<int> vidx = fListBuilder->getTrackList();
  cout << "==>HFJpsiXCand> #tracks = " << vidx.size() << endl;
  cout << "==>HFJpsiXCand> PV idx  = " << (*fPvIdxHandle) << endl;

  TrajectoryStateOnSurface tsos;
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  VertexDistance3D a3d;
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fType = fType;
  pCand->fVtx.fDxy = flxy;
  pCand->fVtx.fDxyE = flxye;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + vidx.size() - 1;

  // -- calculate dimuon mass and momentum
  Track t1(*rm1);
  Track t2(*rm2);
  p4m1.SetXYZM(t1.px(), t1.py(), t1.pz(), MMUON);
  p4m2.SetXYZM(t2.px(), t2.py(), t2.pz(), MMUON);
  p4psi = p4m1 + p4m2;
  pCand->fMass = p4psi.M();
  pCand->fPlab = p4psi.Vect();

  // -- create vector<int, double> with hadronic tracks, sorted by doca to dimuon vertex
  vector<pair<int, double> > trkIdxDoca;
  for (unsigned int ix = 0; ix < vidx.size(); ix++) {
    int idx = vidx[ix];
    // -- skip muons
    if (idx == im1 || idx == im2) continue;
    //    cout << "idx = " << ix << " -> track idx = " << idx << " pv idx = " << pvidx << endl;
    TrackBaseRef baseRef(fTracksHandle, idx);
    TransientTrack transTrack = fTTB->build(*baseRef);
    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), myVertex.position());
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), myVertex);
    trkIdxDoca.push_back(make_pair(idx, doca.value()));
  }
  sort(trkIdxDoca.begin(), trkIdxDoca.end(), docaLess);
  cout << "sorted: " << endl;
  for (vector<pair<int, double> >::iterator it = trkIdxDoca.begin(); it != trkIdxDoca.end(); ++it) {
    cout << Form("idx = %3d, doca = %5.3g", it->first, it->second) << endl;
  }


}




//define this as a plug-in
DEFINE_FWK_MODULE(HFJpsiXCand);
