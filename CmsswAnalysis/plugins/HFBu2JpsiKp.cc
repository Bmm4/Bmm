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
  fPsiLo(iConfig.getUntrackedParameter<double>("psiLo", 2.8)),
  fPsiHi(iConfig.getUntrackedParameter<double>("psiHi", 3.3)) {
  dumpConfiguration();
}

// ----------------------------------------------------------------------
void HFBu2JpsiKp::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBu2JpsiKp configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  psiMuons:                 " << fPsiMuons << endl;
  cout << "---  psiLo:                    " << fPsiLo << endl;
  cout << "---  psiHi:                    " << fPsiHi << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFBu2JpsiKp::analyze(const Event& iEvent, const EventSetup& iSetup) {
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
						     fPsiLo-0.2, fPsiHi+0.2, 1 );

  if (fVerbose > 0) cout << "==>HFBu2JpsiKp> J/psi list size: " << psiList.size() << endl;

  // -- Build J/psi + track
  TAnaCand *pCand(0); 
  double mass(0.1); 
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

    if ((fMaxD0 < 90.) && (tMuon1.d0() > fMaxD0)) continue;
    if ((fMaxDz < 90.) && (tMuon1.dz() > fMaxDz)) continue;
    if ((fMaxD0 < 90.) && (tMuon2.d0() > fMaxD0)) continue;
    if ((fMaxDz < 90.) && (tMuon2.dz() > fMaxDz)) continue;
    
    psi = m1 + m2; 
		
    for (vector<int>::const_iterator trkIt = trkList.begin(); trkIt != trkList.end(); ++trkIt) {
      if (*trkIt == iMuon1 || *trkIt == iMuon2) continue; 
      TrackBaseRef rTrackView(fTracksHandle, *trkIt);
      Track tKaon(*rTrackView);
      if ((fMaxD0 < 90.) && (tKaon.d0() > fMaxD0)) continue;
      if ((fMaxDz < 90.) && (tKaon.dz() > fMaxDz)) continue;
      if (tKaon.pt() < fTrackPt) continue;
      ka.SetXYZM(tKaon.px(), tKaon.py(), tKaon.pz(), MKAON); 
      if ((fDeltaR < 90.) && (psi.DeltaR(ka) > fDeltaR)) continue; 
    
      bu = ka + psi; 
      mass = bu.M();
      
      if (mass < (fCandLo-0.3)) continue;
      if (mass > (fCandHi+0.3)) continue;

      // HFDecayTree:
      // addDecayTree(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false);
      //        clear(int pID, bool doVertexing, double mass, bool massConstraint, double massSigma = -1.0, bool daughtersToPV = false); 

      // -- sequential fit: J/Psi kaon
      theTree.clear(300521, true, MBPLUS, false, -1.0, true);
      HFDecayTreeIterator iterator = theTree.addDecayTree(300443, false, MJPSI, false);
      iterator->addTrack(iMuon1, 13);
      iterator->addTrack(iMuon2, 13);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca,     -1., fMaxDoca, "maxdoca"); 
      iterator->addNodeCut(&HFDecayTree::passMass,     fPsiLo,   fPsiHi,    "mass");
      iterator->addNodeCut(&HFDecayTree::passPt,    fDimuonPt,     1.e9,      "pt"); 

      theTree.addTrack(*trkIt, 321);
      theTree.addNodeCut(&HFDecayTree::passMaxDoca,  -1., fMaxDoca, "maxdoca");
      theTree.addNodeCut(&HFDecayTree::passMass, fCandLo,  fCandHi, "mass");
      theTree.addNodeCut(&HFDecayTree::passFlsxy, fFlsxy,     1.e9, "flsxy");
      theTree.addNodeCut(&HFDecayTree::passPvips,     -1,   fPvIpS, "pvips");

      fSequentialFitter->doFit(&theTree);
      pCand = theTree.getAnaCand();
      if (0 == pCand) continue;
      
      // -- sequential fit: (mass-constrained) J/Psi kaon WILL NOT BE FILLED INTO T1!
      theTree.clear(400521, true, MBPLUS, false, -1.0, true);
      iterator = theTree.addDecayTree(400443, true, MJPSI, true);
      iterator->addTrack(iMuon1, 13);
      iterator->addTrack(iMuon2, 13);
      iterator->addNodeCut(&HFDecayTree::passMaxDoca,     -1., fMaxDoca, "maxdoca"); 
      iterator->addNodeCut(&HFDecayTree::passMass,     fPsiLo,   fPsiHi,    "mass");
      iterator->addNodeCut(&HFDecayTree::passPt,    fDimuonPt,     1.e9,      "pt"); 

      theTree.addTrack(*trkIt, 321);
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
DEFINE_FWK_MODULE(HFBu2JpsiKp);
