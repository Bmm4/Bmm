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
  fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)),
  fPhiWindow(iConfig.getUntrackedParameter<double>("phiWindow", 0.2)),
  fBsWindow(iConfig.getUntrackedParameter<double>("BsWindow", 0.8)) {
  dumpConfiguration();
} // HFBs2JpsiPhi()

// ----------------------------------------------------------------------
void HFBs2JpsiPhi::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBs2JpsiPhi configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiWindow:                " << fPsiWindow << endl;
  cout << "---  phiWindow:                " << fPhiWindow << endl;
  cout << "---  BsWindow:                 " << fBsWindow << endl;
  cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()


// ----------------------------------------------------------------------
void HFBs2JpsiPhi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const double MB0(4.8), MB1(6.0), MJPSI0(3.0), MJPSI1(3.2), MPHI0(0.98), MPHI1(1.06);
  cout << "=== HFBs2JpsiPhi ===================================================================" << endl;
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
						     MJPSI-fPsiWindow, MJPSI+fPsiWindow, 1);
  HFTwoParticleCombinatoricsSet phiList = a.combine(trkList, MKAON, trkList, MKAON, MPHI-fPhiWindow, MPHI+fPhiWindow, 1);
	
  if (fVerbose > 0) cout << "==>HFBs2JpsiPhi> J/psi list size: " << psiList.size() << endl;
  if (fVerbose > 0) cout << "==>HFBs2JpsiPhi> phi list size: " << phiList.size() << endl;
	
  // -- Build J/psi + phi
  TLorentzVector psi, phi, m1, m2, ka1, ka2, bs;
  TAnaCand *pCand(0); 
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
		
    psi = m1 + m2;
    if ((TMath::Abs(psi.M() - MJPSI) > fPsiWindow)) continue;
		
    for (HFTwoParticleCombinatoricsNew::iterator phiIt = phiList.begin(); phiIt != phiList.end(); ++phiIt) {
      unsigned int iKaon1 = phiIt->first;
      unsigned int iKaon2 = phiIt->second;
			
      if (iKaon1 == iMuon1 || iKaon1 == iMuon2) continue;
      if (iKaon2 == iMuon1 || iKaon2 == iMuon2) continue;
      reco::TrackBaseRef rTrackView1(fTracksHandle, iKaon1);
      reco::Track tKaon1(*rTrackView1);
      if (tKaon1.pt() < fTrackPt) continue;
      ka1.SetXYZM(tKaon1.px(), tKaon1.py(), tKaon1.pz(), MKAON); 
      if (psi.DeltaR(ka1) > fDeltaR) continue;
			
      reco::TrackBaseRef rTrackView2(fTracksHandle, iKaon2);
      reco::Track tKaon2(*rTrackView2);
      if (tKaon2.pt() < fTrackPt) continue;
      ka2.SetXYZM(tKaon2.px(), tKaon2.py(), tKaon2.pz(), MKAON); 
      if (psi.DeltaR(ka2) > fDeltaR) continue;
			
      phi = ka1 + ka2;
      if (TMath::Abs(phi.M() - MPHI) > fPhiWindow) continue;
      
      bs = psi + phi; 
      if (TMath::Abs(bs.M() - MBS) > fBsWindow) continue;
			
      // HFDecayTree:
      // addDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false);
      //        clear(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false); 

      // -- sequential fit: J/Psi kaons
      cout << "fitting 300531 (m = " << bs.M() << ") with tracks " << iMuon1 << " " << iMuon2 << " " << iKaon1 << " " << iKaon2 << endl;
      theTree.clear(300531, true, MBS, false, -1.0, true);
      HFDecayTreeIterator iterator = theTree.addDecayTree(300443, false, MJPSI, false);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.maxDoca), &(iterator->fTV.maxDocaV),  -1., fMaxDoca, "300443 maxdoca"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.mass),    &(iterator->fTV.massV),  MJPSI0,   MJPSI1, "300443 J/psi mass"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.pt),      &(iterator->fTV.ptV), fDimuonPt,     1.e9, "300443 pt"));
			
      iterator = theTree.addDecayTree(300333, false, MPHI, false);
      iterator->addTrack(iKaon1,321);
      iterator->addTrack(iKaon2,321);
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.maxDoca), &(iterator->fTV.maxDocaV),  -1., fMaxDoca, "300333 maxdoca"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.mass),    &(iterator->fTV.massV),   MPHI0,    MPHI1, "300333 phi mass"));
      
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.maxDoca),  &(theTree.fTV.maxDocaV),    -1., fMaxDoca, "300531 maxdoca"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.mass),     &(theTree.fTV.massV),       MB0,      MB1, "300531 Bs mass"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.flsxy),    &(theTree.fTV.flsxyV),   fFlsxy,     1.e9, "300531 flsxy"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.pvips),    &(theTree.fTV.pvipsV),       -1,   fPvIpS, "300531 pvips"));
      
      fSequentialFitter->doFit(&theTree);
      pCand = theTree.getAnaCand();
      if (0 == pCand) continue;

      // -- sequential fit: J/Psi (constrained) phi (unconstrained) WILL NOT BE FILLED INTO T1!
      theTree.clear(400531, true, MBS, false, -1.0, true);
      iterator = theTree.addDecayTree(400443, true, MJPSI, true);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.maxDoca), &(iterator->fTV.maxDocaV),  -1., fMaxDoca, "400443 maxdoca"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.mass),    &(iterator->fTV.massV),  MJPSI0,   MJPSI1, "400443 J/psi mass"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.pt),      &(iterator->fTV.ptV), fDimuonPt,     1.e9, "400443 pt"));
			
      iterator = theTree.addDecayTree(400333, false, MPHI, false);
      iterator->addTrack(iKaon1,321);
      iterator->addTrack(iKaon2,321);
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.maxDoca), &(iterator->fTV.maxDocaV), -1., fMaxDoca, "400333 maxdoca"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.mass),    &(iterator->fTV.massV),  MPHI0,    MPHI1, "400333 phi mass"));
			
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.maxDoca), &(theTree.fTV.maxDocaV),  -1., fMaxDoca, "400531 maxdoca"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.mass),    &(theTree.fTV.massV),     MB0,      MB1, "400531 Bs mass"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.flsxy),   &(theTree.fTV.flsxyV), fFlsxy,     1.e9, "400531 flsxy"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.pvips),    &(theTree.fTV.pvipsV),    -1,   fPvIpS, "400531 pvips"));
      // -- the following is never true, therefore the 400531 candidate will not be filled into the T1 tree
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.zero),    &(theTree.fTV.zeroV),      1.,       2., "400531 zero"));

      fSequentialFitter->doFit(&theTree);
      // -- but we store its relevant information into the unconstrained candidate, saved above
      pCand->fDouble1 = theTree.fTV.mass;
      pCand->fDouble2 = theTree.fTV.masserr;
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFBs2JpsiPhi);
