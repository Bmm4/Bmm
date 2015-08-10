#include "HFBd2JpsiKstar.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"

#include <iostream>

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFBd2JpsiKstar::HFBd2JpsiKstar(const edm::ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)),
  fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)),
  fK0Window(iConfig.getUntrackedParameter<double>("k0Window", 0.2)),
  fB0Window(iConfig.getUntrackedParameter<double>("B0Window", 0.8)) {
  dumpConfiguration();
} // HFBd2JpsiKstar()

// ----------------------------------------------------------------------
void HFBd2JpsiKstar::dumpConfiguration() {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBd2JpsiKstar configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiWindow:                " << fPsiWindow << endl;
  cout << "---  k0Window:                 " << fK0Window << endl;
  cout << "---  BdWindow:                 " << fB0Window << endl;
  cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()


// ----------------------------------------------------------------------
void HFBd2JpsiKstar::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
	
  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  }
  catch (HFSetupException e) {
    cout << "==>HFBd2JpsiKstar> " << e.fMsg << endl;
    return;
  }
	
  fListBuilder->setMinPt(fMuonPt); // work with muon pt and not with track pt
  vector<int> muonList = fListBuilder->getMuonList();
  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setMaxDocaToTracks(fMaxDoca);
  fListBuilder->setCloseTracks(&muonList);
  vector<int> trkList = fListBuilder->getTrackList();
	
  if (muonList.size() < static_cast<unsigned int>(fPsiMuons)) return; // not enough muons
  HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
  HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList),MMUON,
						     (fPsiMuons < 2 ? trkList : muonList),MMUON,
						     MJPSI-fPsiWindow,MJPSI+fPsiWindow, 1);
  
  //////adapted  by jmonroy 8-7-2015 

  HFTwoParticleCombinatoricsSet kstarList = a.combine(trkList,MKAON,trkList,MPION,MKSTAR-fK0Window,MKSTAR+fK0Window, 0);

  if (fVerbose > 0) cout << "==>HFBs2JpsiKstar> J/psi list size: " << psiList.size() << endl;
  if (fVerbose > 0) cout << "==>HFBd2JpsiKstar> kstar list size: " << kstarlist.size() << endl;
	
  // -- Build J/psi + kstar
  TLorentzVector psi, m1, m2, kstar, ka, pi, b0;
 

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
		

    for (HFTwoParticleCombinatoricsNew::iterator kstarIt = kstarList.begin(); kstarIt != kstarList.end(); ++kstarIt) {
      unsigned int iKaon = kstarIt->first;
      unsigned int iPion = kstarIt->second;
      
      if (iKaon == iMuon1 || iKaon  == iMuon2) continue;

      reco::TrackBaseRef rTrackView1(fTracksHandle, iKaon);
      reco::Track tKaon(*rTrackView1);
      if (tKaon.pt() < fTrackPt) continue;
      ka.SetXYZM(tKaon.px(), tKaon.py(), tKaon.pz(), MKAON); 
      if (psi.DeltaR(ka) > fDeltaR) continue;


      if (iPion  == iMuon1 || iPion == iMuon2) continue;

      reco::TrackBaseRef rTrackView2(fTracksHandle, iPion);
      reco::Track tPion(*rTrackView2);
      if (tPion.pt() < fTrackPt) continue;
      pi.SetXYZM(tPion.px(), tPion.py(), tPion.pz(), MPION); 
      if (psi.DeltaR(pi) > fDeltaR) continue;

      kstar = ka + pi;
      if ((TMath::Abs(kstar.M() - MKSTAR) > fK0Window)) continue;

      b0 = psi + kstar; 
      if (TMath::Abs(b0.M() - MB_0) > fB0Window) continue;

      // -- sequential fit: J/Psi kaon + pion
      HFDecayTree theTree(310511, true, MB_0, false, -1.0, true);
			
      HFDecayTreeIterator iterator = theTree.addDecayTree(300443, false, MJPSI, false); // Don't use kinematic particle for the Psi
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      iterator = theTree.addDecayTree(300892, false, MKSTAR, false);
      iterator->addTrack(iKaon,321);
      iterator->addTrack(iPion,211);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      fSequentialFitter->doFit(&theTree);
			
      // -- sequential fit: J/Psi (constraint) kstar (unconstraint)
      theTree.clear(410511, true, MB_0, false, -1.0, true);
      iterator = theTree.addDecayTree(400443, true, MJPSI, true);
      iterator->addTrack(iMuon1,13);
      iterator->addTrack(iMuon2,13);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      iterator = theTree.addDecayTree(400892, false, MKSTAR, false);
      iterator->addTrack(iKaon,321);
      iterator->addTrack(iPion,211);
      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
			
      fSequentialFitter->doFit(&theTree);

    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFBd2JpsiKstar);
