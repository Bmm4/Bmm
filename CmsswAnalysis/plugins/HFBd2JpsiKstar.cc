#include "HFBd2JpsiKstar.h"

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
HFBd2JpsiKstar::HFBd2JpsiKstar(const edm::ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)),
  fPsiLo(iConfig.getUntrackedParameter<double>("psiLo", 2.8)),
  fPsiHi(iConfig.getUntrackedParameter<double>("psiHi", 3.3)),
  fBdWindow(iConfig.getUntrackedParameter<double>("BdWindow", 0.8)),
  fKstarLo(iConfig.getUntrackedParameter<double>("kstarLo", 0.7)),
  fKstarHi(iConfig.getUntrackedParameter<double>("kstarHi", 1.1)){
  dumpConfiguration();
}


// ----------------------------------------------------------------------
void HFBd2JpsiKstar::dumpConfiguration() {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBd2JpsiKstar configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiLo:                    " << fPsiLo << endl;
  cout << "---  psiHi:                    " << fPsiHi << endl;
  cout << "---  kstarLo:                  " << fKstarLo << endl;
  cout << "---  kstarHi:                  " << fKstarHi << endl;
  cout << "---  BdWindow:                 " << fBdWindow << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFBd2JpsiKstar::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;

  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch (HFSetupException e) {
    cout << "==>HFBd2JpsiKstar> " << e.fMsg << endl;
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
  HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList), MMUON,
						     (fPsiMuons < 2 ? trkList : muonList), MMUON,
						     fPsiLo-0.2, fPsiHi+0.2, 1);

  HFTwoParticleCombinatoricsSet kstarList = a.combine(trkList, MKAON, trkList, MPION, fKstarLo, fKstarHi, 0);

  if (fVerbose > 0) cout << "==>HFBd2JpsiKstar> J/psi list size: " << psiList.size() << endl;
  if (fVerbose > 0) cout << "==>HFBd2JpsiKstar> kstar list size: " << kstarList.size() << endl;

  // -- Build J/psi + kstar
  TLorentzVector psi, m1, m2, kstar, ka, pi, b0;
  double mass(0.);
  TAnaCand *pCand(0);
  HFDecayTreeIterator iterator;
  HFDecayTree theTree(300511, true, MB_0, false, -1.0, true);
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

    for (HFTwoParticleCombinatoricsNew::iterator kstarIt = kstarList.begin(); kstarIt != kstarList.end(); ++kstarIt) {
      unsigned int iKaon = kstarIt->first;
      unsigned int iPion = kstarIt->second;

      if (iKaon == iMuon1 || iKaon  == iMuon2) continue;
      if (iPion  == iMuon1 || iPion == iMuon2) continue;
      reco::TrackBaseRef rTrackView1(fTracksHandle, iKaon);
      reco::Track tKaon(*rTrackView1);
      if (tKaon.pt() < fTrackPt) continue;
      if ((fMaxD0 < 90.) && (tKaon.d0() > fMaxD0)) continue;
      if ((fMaxDz < 90.) && (tKaon.dz() > fMaxDz)) continue;

      reco::TrackBaseRef rTrackView2(fTracksHandle, iPion);
      reco::Track tPion(*rTrackView2);
      if (tPion.pt() < fTrackPt) continue;
      if ((fMaxD0 < 90.) && (tPion.d0() > fMaxD0)) continue;
      if ((fMaxDz < 90.) && (tPion.dz() > fMaxDz)) continue;
      if (tKaon.charge()*tPion.charge() > 0) continue;

      ka.SetXYZM(tKaon.px(), tKaon.py(), tKaon.pz(), MKAON);
      if ((fDeltaR < 90.) && (psi.DeltaR(ka) > fDeltaR)) continue;
      pi.SetXYZM(tPion.px(), tPion.py(), tPion.pz(), MPION);
      if ((fDeltaR < 90.) && (psi.DeltaR(pi) > fDeltaR)) continue;

      kstar = ka + pi;

      b0 = psi + kstar;
      mass = b0.M();
      if (mass < (fCandLo-0.3)) continue;
      if (mass > (fCandHi+0.3)) continue;

      // -- sequential fit: J/Psi kaons
      theTree.clear(300511, true, MB_0, false, -1.0, true);

      iterator = theTree.addDecayTree(300443, false, MJPSI, false);
      iterator->addTrack(iMuon1, 13);
      iterator->addTrack(iMuon2, 13);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca,     -1., fMaxDoca, "maxdoca");
      iterator->addNodeCut(&HFDecayTree::passMass,     fPsiLo,   fPsiHi,    "mass");
      iterator->addNodeCut(&HFDecayTree::passPt,    fDimuonPt,     1.e9,      "pt");

      iterator = theTree.addDecayTree(300313, false, MKSTAR, false);
      iterator->addTrack(iKaon, 321);
      iterator->addTrack(iPion, 211);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca,  -1., fMaxDoca, "maxdoca");
      iterator->addNodeCut(&HFDecayTree::passMass, fKstarLo,  fKstarHi,    "mass");

      theTree.addNodeCut(&HFDecayTree::passMaxDoca,  -1., fMaxDoca, "maxdoca");
      theTree.addNodeCut(&HFDecayTree::passMass, fCandLo,  fCandHi, "mass");
      theTree.addNodeCut(&HFDecayTree::passFlsxy, fFlsxy,     1.e9, "flsxy");
      theTree.addNodeCut(&HFDecayTree::passPvips,     -1,   fPvIpS, "pvips");

      fSequentialFitter->doFit(&theTree);
      pCand = theTree.getAnaCand();
      if (0 == pCand) continue;

      // -- sequential fit: J/Psi (constrained) phi (unconstrained) WILL NOT BE FILLED INTO T1!
      theTree.clear(400511, true, MB_0, false, -1.0, true);
      iterator = theTree.addDecayTree(400443, true, MJPSI, true);
      iterator->addTrack(iMuon1, 13);
      iterator->addTrack(iMuon2, 13);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca,     -1., fMaxDoca, "maxdoca");
      iterator->addNodeCut(&HFDecayTree::passMass,     fPsiLo,   fPsiHi,    "mass");
      iterator->addNodeCut(&HFDecayTree::passPt,    fDimuonPt,     1.e9,      "pt");

      iterator = theTree.addDecayTree(400313, false, MKSTAR, false);
      iterator->addTrack(iKaon, 321);
      iterator->addTrack(iPion, 211);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca,  -1., fMaxDoca, "maxdoca");
      iterator->addNodeCut(&HFDecayTree::passMass, fKstarLo,  fKstarHi,    "mass");

      theTree.addNodeCut(&HFDecayTree::passMaxDoca,  -1., fMaxDoca, "maxdoca");
      theTree.addNodeCut(&HFDecayTree::passMass, fCandLo,  fCandHi, "mass");
      theTree.addNodeCut(&HFDecayTree::passFlsxy, fFlsxy,     1.e9, "flsxy");
      theTree.addNodeCut(&HFDecayTree::passPvips,     -1,   fPvIpS, "pvips");
      theTree.addNodeCut(&HFDecayTree::passNever,     1.,       1., "never");

      fSequentialFitter->doFit(&theTree);

      // -- but we store its relevant information into the unconstrained candidate, saved above
      pCand->fMassC  = theTree.fTV.mass;
      pCand->fMassCE = theTree.fTV.masserr;

    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFBd2JpsiKstar);
