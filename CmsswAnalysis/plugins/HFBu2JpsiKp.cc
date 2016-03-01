#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBu2JpsiKp.h"

#include <iostream>
#include <utility>

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"

#include "Bmm/CmsswAnalysis/interface/HFSequentialVertexFit.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFDecayTree.hh"
#include "Bmm/CmsswAnalysis/interface/HFTrackListBuilder.hh"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace reco;
using namespace edm;


// ----------------------------------------------------------------------
HFBu2JpsiKp::HFBu2JpsiKp(const ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fPsiMuons(iConfig.getUntrackedParameter<int>("psiMuons", 2)),
  fPsiWindow(iConfig.getUntrackedParameter<double>("psiWindow", 0.3)),
  fBuWindow(iConfig.getUntrackedParameter<double>("BuWindow", 0.8)) {
  dumpConfiguration();
}

// ----------------------------------------------------------------------
void HFBu2JpsiKp::dumpConfiguration() {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBu2JpsiKp configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiWindow:                " << fPsiWindow << endl;
  cout << "---  BuWindow:                " << fBuWindow << endl;
  cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()

// ----------------------------------------------------------------------
void HFBu2JpsiKp::analyze(const Event& iEvent, const EventSetup& iSetup) {
  const double MB0(4.8), MB1(6.0), MJPSI0(3.0), MJPSI1(3.2);
  cout << "=== HFBu2JpsiKp ===================================================================" << endl;
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;

  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFBu2JpsiKp> " << e.fMsg << endl;
    return;
  }

  fListBuilder->setMinPt(fMuonPt); 
  vector<int> muonList = fListBuilder->getMuonList();
  if (muonList.size() < static_cast<unsigned int>(fPsiMuons)) return; // not enough muons
	
  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setMaxDocaToTracks(fMaxDoca);
  fListBuilder->setCloseTracks(&muonList);
  vector<int> trkList = fListBuilder->getTrackList();
	
  HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
  HFTwoParticleCombinatoricsSet psiList = a.combine( (fPsiMuons < 1 ? trkList : muonList), MMUON,
						     (fPsiMuons < 2 ? trkList : muonList), MMUON,
						     MJPSI - fPsiWindow, MJPSI + fPsiWindow, 1 );

  if (fVerbose > 0) cout << "==>HFBu2JpsiKp> J/psi list size: " << psiList.size() << endl;

  // -- Build J/psi + track
  TAnaCand *pCand(0); 
  TLorentzVector psi, cpsi, m1, m2, ka, bu;
  HFDecayTree theTree(300521, true, MBPLUS, false, -1.0, true);
  for (HFTwoParticleCombinatoricsNew::iterator psiIt = psiList.begin(); psiIt != psiList.end(); ++psiIt) {
    int iMuon1 = psiIt->first; 
    int iMuon2 = psiIt->second; 

    reco::TrackBaseRef mu1TrackView(fTracksHandle, iMuon1);
    reco::Track tMuon1(*mu1TrackView);
    if (tMuon1.pt() < fMuonPt)  continue;
    m1.SetPtEtaPhiM(tMuon1.pt(), tMuon1.eta(), tMuon1.phi(), MMUON); 

    TrackBaseRef mu2TrackView(fTracksHandle, iMuon2);
    Track tMuon2(*mu2TrackView);
    if (tMuon2.pt() < fMuonPt)  continue;
    m2.SetPtEtaPhiM(tMuon2.pt(), tMuon2.eta(), tMuon2.phi(), MMUON); 

    psi = m1 + m2; 
    if ((TMath::Abs(psi.M() - MJPSI) > fPsiWindow)) continue;
		
    for (vector<int>::const_iterator trkIt = trkList.begin(); trkIt != trkList.end(); ++trkIt) {
      if (*trkIt == iMuon1 || *trkIt == iMuon2) continue; 
      TrackBaseRef rTrackView(fTracksHandle, *trkIt);
      Track tKaon(*rTrackView);
      if (tKaon.d0() > fMaxD0) continue;
      if (tKaon.dz() > fMaxDz) continue;
      if (tKaon.pt() < fTrackPt) continue;
      ka.SetXYZM(tKaon.px(), tKaon.py(), tKaon.pz(), MKAON); 
      if (psi.DeltaR(ka) > fDeltaR) continue; 

      bu = ka + psi; 
      if (TMath::Abs(bu.M() - MBPLUS) > fBuWindow) continue;

      // HFDecayTree:
      // addDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false);
      //        clear(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false); 

      // -- sequential fit: J/Psi kaon
      cout << "fitting 300521 (m = " << bu.M() << ") with tracks " << iMuon1 << " " << iMuon2 << " " << *trkIt << endl;
      theTree.clear(300521, true, MBPLUS, false, -1.0, true);
      HFDecayTreeIterator iterator = theTree.addDecayTree(300443, false, MJPSI, false);
      iterator->addTrack(iMuon1, 13);
      iterator->addTrack(iMuon2, 13);
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.maxDoca), &(iterator->fTV.maxDocaV),  -1., fMaxDoca, "300443 maxdoca"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.mass),    &(iterator->fTV.massV),  MJPSI0,   MJPSI1, "300443 J/psi mass"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.pt),      &(iterator->fTV.ptV), fDimuonPt,     1.e9, "300443 pt"));

      theTree.addTrack(*trkIt, 321);
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.maxDoca), &(theTree.fTV.maxDocaV),  -1., fMaxDoca, "300521 maxdoca"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.mass),    &(theTree.fTV.massV),     MB0,      MB1, "300521 B+ mass"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.flsxy),   &(theTree.fTV.flsxyV), fFlsxy,     1.e9, "300521 flsxy"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.pvips),   &(theTree.fTV.pvipsV), -1,       fPvIpS, "300521 pvips"));

      fSequentialFitter->doFit(&theTree);
      pCand = theTree.getAnaCand();
      if (0 == pCand) continue;
      
      // -- sequential fit: (mass-constrained) J/Psi kaon WILL NOT BE FILLED INTO T1!
      theTree.clear(400521, true, MBPLUS, false, -1.0, true);
      iterator = theTree.addDecayTree(400443, true, MJPSI, true);
      iterator->addTrack(iMuon1, 13);
      iterator->addTrack(iMuon2, 13);
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.maxDoca), &(iterator->fTV.maxDocaV), -1., fMaxDoca, "400443 maxdoca"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.mass),    &(iterator->fTV.massV), MJPSI0,   MJPSI1, "400443 J/psi mass"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.pt),      &(iterator->fTV.ptV), fDimuonPt,     1.e9, "300443 pt"));
      if (fVerbose > 5) cout << "==>HFBu2JpsiKp> sequential fit with mass constraint" << endl;

      theTree.addTrack(*trkIt, 321);
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.maxDoca), &(theTree.fTV.maxDocaV),  -1., fMaxDoca, "400521 maxdoca"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.mass),    &(theTree.fTV.massV),     MB0,      MB1, "400521 B+ mass"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.flsxy),   &(theTree.fTV.flsxyV), fFlsxy,     1.e9, "400521 flsxy"));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.pvips),   &(theTree.fTV.pvipsV), -1,       fPvIpS, "400521 pvips"));
      // -- the following is never true, therefore the 400521 candidate will not be filled into the T1 tree
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.zero),    &(theTree.fTV.zeroV),      1.,       2., "400521 zero"));

      fSequentialFitter->doFit(&theTree);
      // -- but we store its relevant information into the unconstrained candidate, saved above
      pCand->fDouble1 = theTree.fTV.mass;
      pCand->fDouble2 = theTree.fTV.masserr;
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFBu2JpsiKp);
