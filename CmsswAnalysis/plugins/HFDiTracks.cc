#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDiTracks.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"

using std::cout;
using std::endl;

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFDiTracks::HFDiTracks(const edm::ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fLeadingTrackPt(iConfig.getUntrackedParameter<double>("leadingTrackPt", -1.0)),
  fTrack1Mass(iConfig.getUntrackedParameter<double>("track1Mass", MMUON)),
  fTrack2Mass(iConfig.getUntrackedParameter<double>("track2Mass", MMUON)),
  fExtra(iConfig.getUntrackedParameter<double>("extra", 0.3)),
  fUnlikeCharge(iConfig.getUntrackedParameter<bool>("unlikeCharge", true)),
  fNbrMuons(iConfig.getUntrackedParameter<int>("nbrMuons",2)),
  fCloseToMuons(iConfig.getUntrackedParameter<bool>("closeToMuons",false)) {
  dumpConfiguration();
}


// ----------------------------------------------------------------------
void HFDiTracks::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDiTracks constructor" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  leadingTrackPt:           " << fLeadingTrackPt << endl;
  cout << "---  track1Mass:               " << fTrack1Mass << endl;
  cout << "---  track2Mass:               " << fTrack2Mass << endl;
  cout << "---  extra:                    " << fExtra << endl;
  cout << "---  nbrMuons:                 " << fNbrMuons << endl;
  cout << "---  closeToMuons:             " << fCloseToMuons << endl;
  cout << "---  unlikeCharge:             " << fUnlikeCharge << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFDiTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (fVerbose > 0)  cout << "=== HFDiTracks run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event()
			  << " ==================================================================="
			  << endl;
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;

  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  }
  catch (HFSetupException e) {
    cout << "==>HFDiTracks> " << e.fMsg << endl;
    return;
  }

  std::vector<int> trkList, ltrkList;
  std::vector<int> muonList;

  if (fNbrMuons > 0 || fCloseToMuons) {
    fListBuilder->setMinPt(fMuonPt);
    muonList = fListBuilder->getMuonList();
  }
  if (fNbrMuons < 2) {
    fListBuilder->setMinPt(fTrackPt);
    if (fCloseToMuons) {
      fListBuilder->setMaxDocaToTracks(fMaxDoca);
      fListBuilder->setCloseTracks(&muonList);
    }
    trkList = fListBuilder->getTrackList();
    // -- second list with leading track pt cut
    if (fLeadingTrackPt > 0.) {
      fListBuilder->setMinPt(fLeadingTrackPt);
      if (fCloseToMuons) {
	fListBuilder->setMaxDocaToTracks(fMaxDoca);
	fListBuilder->setCloseTracks(&muonList);
      }
      ltrkList = fListBuilder->getTrackList();
    }
  }

  HFTwoParticleCombinatoricsNew a(fTracksHandle, fVerbose);
  HFTwoParticleCombinatoricsSet candSet;
  if (fNbrMuons < 1) {
    candSet = a.combine((fLeadingTrackPt > 0.? ltrkList: trkList), fTrack1Mass,
			trkList, fTrack2Mass,
			fCandLo - fExtra, fCandHi + fExtra,
			TMath::Abs(fTrack1Mass - fTrack2Mass) < 0.0001);
  } else if (fNbrMuons < 2) {
    candSet = a.combine(muonList, fTrack1Mass,
			trkList, fTrack2Mass,
			fCandLo - fExtra, fCandHi + fExtra,
			TMath::Abs(fTrack1Mass - fTrack2Mass) < 0.0001);
  } else {
    candSet = a.combine(muonList, fTrack1Mass,
			muonList, fTrack2Mass,
			fCandLo - fExtra, fCandHi + fExtra,
			TMath::Abs(fTrack1Mass - fTrack2Mass) < 0.0001);
  }
  if (fVerbose > 0) cout << "==>HFDiTracks> candidate list size: " << candSet.size() << endl;

  TAnaCand *pCand(0);
  for (HFTwoParticleCombinatoricsNew::iterator trkIt = candSet.begin(); trkIt != candSet.end(); ++trkIt) {
    reco::TrackBaseRef rTrackView1(fTracksHandle, trkIt->first);
    reco::Track t1(*rTrackView1);
    reco::TrackBaseRef rTrackView2(fTracksHandle, trkIt->second);
    reco::Track t2(*rTrackView2);
    if (fUnlikeCharge && (t1.charge() * t2.charge() > 0)) continue;
    if (fDeltaR < 90.) {
      TLorentzVector p4t1, p4t2;
      p4t1.SetXYZM(t1.px(), t1.py(), t1.pz(), fTrack1Mass);
      p4t2.SetXYZM(t2.px(), t2.py(), t2.pz(), fTrack2Mass);
      if (p4t1.DeltaR(p4t2) > fDeltaR) {
	if (fVerbose) cout << fType << " skipping because p4t1.DeltaR(p4t2) = " << p4t1.DeltaR(p4t2) << " > " << fDeltaR << endl;
	continue;
      }
    }

    HFDecayTree theTree(fType, true, 0, false);
    theTree.addTrack(trkIt->first, idFromMass(fTrack1Mass));
    theTree.addTrack(trkIt->second, idFromMass(fTrack2Mass));
    theTree.addNodeCut(&HFDecayTree::passMaxDoca,   -1., fMaxDoca, "maxdoca");
    theTree.addNodeCut(&HFDecayTree::passMass, fCandLo,   fCandHi, "mass");
    theTree.addNodeCut(&HFDecayTree::passPt,   fCandPt,     1.e9, "pt");
    if (fFlxy > 0)  theTree.addNodeCut(&HFDecayTree::passFlxy,      -1.,    fFlxy, "flxy");
    if (fFlsxy > 0) theTree.addNodeCut(&HFDecayTree::passFlsxy,  fFlsxy,     1.e9, "flsxy");
    if (fPvIpS > 0) theTree.addNodeCut(&HFDecayTree::passPvips,      -1,   fPvIpS, "pvips");
    if (fChi2 > 0)  theTree.addNodeCut(&HFDecayTree::passChi2Dof,   -1.,    fChi2, "chi2dof");

    fSequentialFitter->doFit(&theTree);
    if (fVerbose) {
      theTree.dump();
      pCand = theTree.getAnaCand();
      if (0 != pCand) {
	cout << "==> filled this candidate!" << endl;
      }
    }
  }
}
// ----------------------------------------------------------------------
int HFDiTracks::idFromMass(double mass) {
  if (TMath::Abs(mass - MELECTRON) < 0.0001) return 11;
  if (TMath::Abs(mass - MMUON) < 0.0001) return 13;
  if (TMath::Abs(mass - MPION) < 0.0001) return 211;
  if (TMath::Abs(mass - MKAON) < 0.0001) return 321;
  if (TMath::Abs(mass - MPROTON) < 0.0001) return 2212;
  cout << "%%%> HFDiTracks: mass " << mass << " not associated to any stable particle" << endl;
  return 0;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDiTracks);
