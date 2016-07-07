#include "HFBs2JpsiPhi.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"

#include <iostream>

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace reco;
using namespace edm;

// ----------------------------------------------------------------------
HFBs2JpsiPhi::HFBs2JpsiPhi(const edm::ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)),
  fPsiLo(iConfig.getUntrackedParameter<double>("psiLo", 2.8)),
  fPsiHi(iConfig.getUntrackedParameter<double>("psiHi", 3.3)),
  fPhiLo(iConfig.getUntrackedParameter<double>("phiLo", 0.98)),
  fPhiHi(iConfig.getUntrackedParameter<double>("phiHi", 1.06)) {
  dumpConfiguration();
}

// ----------------------------------------------------------------------
void HFBs2JpsiPhi::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBs2JpsiPhi configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiLo:                    " << fPsiLo << endl;
  cout << "---  psiHi:                    " << fPsiHi << endl;
  cout << "---  phiLo:                    " << fPhiLo << endl;
  cout << "---  phiHi:                    " << fPhiHi << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFBs2JpsiPhi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;

  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  }
  catch (HFSetupException e) {
    cout << "==>HFBs2JpsiPhi> " << e.fMsg << endl;
    return;
  }

  fListBuilder->setMinPt(fMuonPt);
  vector<int> muonList = fListBuilder->getMuonList();
  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setMaxDocaToTracks(fMaxDoca);
  fListBuilder->setCloseTracks(&muonList);
  vector<int> trkList = fListBuilder->getTrackList();

  if (muonList.size() < static_cast<unsigned int>(fPsiMuons)) return;
  HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
  HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList),MMUON,
						     (fPsiMuons < 2 ? trkList : muonList),MMUON,
						     fPsiLo-0.2, fPsiHi+0.2, 1);
  HFTwoParticleCombinatoricsSet phiList = a.combine(trkList, MKAON, trkList, MKAON, 0.97, 1.10, 1);

  if (fVerbose > 0) cout << "==>HFBs2JpsiPhi> J/psi list size: " << psiList.size() << endl;
  if (fVerbose > 0) cout << "==>HFBs2JpsiPhi> phi list size: " << phiList.size() << endl;

  // -- Build J/psi + phi
  TLorentzVector psi, phi, m1, m2, ka1, ka2, bs;
  double mass(0.);
  TAnaCand *pCand(0);
  HFDecayTreeIterator iterator;
  HFDecayTree theTree(300531, true, MBS, false, -1.0, true);
  for (HFTwoParticleCombinatoricsNew::iterator psiIt = psiList.begin(); psiIt != psiList.end(); ++psiIt) {
    unsigned int iMuon1 = psiIt->first;
    unsigned int iMuon2 = psiIt->second;

    reco::TrackBaseRef mu1TrackView(fTracksHandle, iMuon1);
    reco::Track tMuon1(*mu1TrackView);
    if (tMuon1.pt() < fMuonPt)  continue;
    m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON);

    reco::TrackBaseRef mu2TrackView(fTracksHandle, iMuon2);
    reco::Track tMuon2(*mu2TrackView);
    if (tMuon2.pt() < fMuonPt)  continue;
    m2.SetPtEtaPhiM(tMuon2.pt(), tMuon2.eta(), tMuon2.phi(), MMUON);

    if ((fMaxD0 < 90.) && (tMuon1.d0() > fMaxD0)) continue;
    if ((fMaxDz < 90.) && (tMuon1.dz() > fMaxDz)) continue;
    if ((fMaxD0 < 90.) && (tMuon2.d0() > fMaxD0)) continue;
    if ((fMaxDz < 90.) && (tMuon2.dz() > fMaxDz)) continue;

    psi = m1 + m2;

    for (HFTwoParticleCombinatoricsNew::iterator phiIt = phiList.begin(); phiIt != phiList.end(); ++phiIt) {
      unsigned int iKaon1 = phiIt->first;
      unsigned int iKaon2 = phiIt->second;

      if (iKaon1 == iMuon1 || iKaon1 == iMuon2) continue;
      if (iKaon2 == iMuon1 || iKaon2 == iMuon2) continue;
      reco::TrackBaseRef rTrackView1(fTracksHandle, iKaon1);
      reco::Track tKaon1(*rTrackView1);
      if (tKaon1.pt() < fTrackPt) continue;
      if ((fMaxD0 < 90.) && (tKaon1.d0() > fMaxD0)) continue;
      if ((fMaxDz < 90.) && (tKaon1.dz() > fMaxDz)) continue;
      ka1.SetXYZM(tKaon1.px(), tKaon1.py(), tKaon1.pz(), MKAON);
      if ((fDeltaR < 90.) && (psi.DeltaR(ka1) > fDeltaR)) continue;
      reco::TrackBaseRef rTrackView2(fTracksHandle, iKaon2);
      reco::Track tKaon2(*rTrackView2);
      if (tKaon2.pt() < fTrackPt) continue;
      if ((fMaxD0 < 90.) && (tKaon2.d0() > fMaxD0)) continue;
      if ((fMaxDz < 90.) && (tKaon2.dz() > fMaxDz)) continue;
      ka2.SetXYZM(tKaon2.px(), tKaon2.py(), tKaon2.pz(), MKAON);
      if ((fDeltaR < 90.) && (psi.DeltaR(ka2) > fDeltaR)) continue;

      phi = ka1 + ka2;

      bs = psi + phi;
      mass = bs.M();
      if (mass < (fCandLo-0.3)) continue;
      if (mass > (fCandHi+0.3)) continue;

      // -- sequential fit: J/Psi kaons
      theTree.clear(300531, true, MBS, false, -1.0, true);
      iterator = theTree.addDecayTree(300443, false, MJPSI, false);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca,     -1., fMaxDoca, "maxdoca");
      iterator->addNodeCut(&HFDecayTree::passMass,     fPsiLo,   fPsiHi,    "mass");
      iterator->addNodeCut(&HFDecayTree::passPt,    fDimuonPt,     1.e9,      "pt");

      iterator = theTree.addDecayTree(300333, false, MPHI, false);
      iterator->addTrack(iKaon1,321);
      iterator->addTrack(iKaon2,321);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca, -1., fMaxDoca, "maxdoca");
      iterator->addNodeCut(&HFDecayTree::passMass, fPhiLo,   fPhiHi,    "mass");

      theTree.addNodeCut(&HFDecayTree::passMaxDoca,  -1., fMaxDoca, "maxdoca");
      theTree.addNodeCut(&HFDecayTree::passMass, fCandLo,  fCandHi, "mass");
      theTree.addNodeCut(&HFDecayTree::passFlsxy, fFlsxy,     1.e9, "flsxy");
      theTree.addNodeCut(&HFDecayTree::passPvips,     -1,   fPvIpS, "pvips");

      fSequentialFitter->doFit(&theTree);
      pCand = theTree.getAnaCand();
      if (0 == pCand) continue;

      // -- sequential fit: J/Psi (constrained) phi (unconstrained) WILL NOT BE FILLED INTO T1!
      theTree.clear(400531, true, MBS, false, -1.0, true);
      iterator = theTree.addDecayTree(400443, true, MJPSI, true);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca,     -1., fMaxDoca, "maxdoca");
      iterator->addNodeCut(&HFDecayTree::passMass,     fPsiLo,   fPsiHi,    "mass");
      iterator->addNodeCut(&HFDecayTree::passPt,    fDimuonPt,     1.e9,      "pt");

      // FIXME: try to simply call set_massConstraint() and change it instead of setting up a new decay tree!
      iterator = theTree.addDecayTree(400333, false, MPHI, false);
      iterator->addTrack(iKaon1,321);
      iterator->addTrack(iKaon2,321);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca, -1., fMaxDoca, "maxdoca");
      iterator->addNodeCut(&HFDecayTree::passMass, fPhiLo,   fPhiHi,    "mass");

      theTree.addNodeCut(&HFDecayTree::passMaxDoca,  -1., fMaxDoca, "maxdoca");
      theTree.addNodeCut(&HFDecayTree::passMass, fCandLo,  fCandHi, "mass");
      theTree.addNodeCut(&HFDecayTree::passFlsxy, fFlsxy,     1.e9, "flsxy");
      theTree.addNodeCut(&HFDecayTree::passPvips,     -1,   fPvIpS, "pvips");
      theTree.addNodeCut(&HFDecayTree::passNever,     1.,       1., "never");

      fSequentialFitter->doFit(&theTree);
      // -- but we store its relevant information into the unconstrained candidate, saved above
      pCand->fDouble1 = theTree.fTV.mass;
      pCand->fDouble2 = theTree.fTV.masserr;
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFBs2JpsiPhi);
