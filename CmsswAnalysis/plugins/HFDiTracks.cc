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
  fMassLow(iConfig.getUntrackedParameter<double>("massLow", 0.0)),
  fMassHigh(iConfig.getUntrackedParameter<double>("massHigh", 12.0)),
  fNbrMuons(iConfig.getUntrackedParameter<int>("nbrMuons",2)),
  fCloseToMuons(iConfig.getUntrackedParameter<bool>("closeToMuons",false)) {
  dumpConfiguration();
} // HFDiTracks()


// ----------------------------------------------------------------------
void HFDiTracks::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDiTracks constructor" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  leadingTrackPt:           " << fLeadingTrackPt << endl;
  cout << "---  track1Mass:               " << fTrack1Mass << endl;
  cout << "---  track2Mass:               " << fTrack2Mass << endl;
  cout << "---  massLow:                  " << fMassLow << endl;
  cout << "---  massHigh:                 " << fMassHigh << endl;
  cout << "---  nbrMuons:                 " << fNbrMuons << endl;
  cout << "---  closeToMuons:             " << fCloseToMuons << endl;
  cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()


// ----------------------------------------------------------------------
void HFDiTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  cout << "=== HFDiTracks ===================================================================" << endl;
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
			fMassLow, fMassHigh,
			TMath::Abs(fTrack1Mass - fTrack2Mass) < 0.0001);
  } else if (fNbrMuons < 2) {
    candSet = a.combine(muonList, fTrack1Mass,
			trkList, fTrack2Mass,
			fMassLow, fMassHigh,
			TMath::Abs(fTrack1Mass - fTrack2Mass) < 0.0001);
  } else {
    candSet = a.combine(muonList, fTrack1Mass,
			muonList, fTrack2Mass,
			fMassLow, fMassHigh,
			TMath::Abs(fTrack1Mass - fTrack2Mass) < 0.0001);
  } 
  if (fVerbose > 0) cout << "==>HFDiTracks> candidate list size: " << candSet.size() << endl;
	
  for (HFTwoParticleCombinatoricsNew::iterator trkIt = candSet.begin(); trkIt != candSet.end(); ++trkIt) {
    cout << "fitting tracks " << trkIt->first << " " << trkIt->second << endl;
    HFDecayTree theTree(fType, true, 0, false);
    theTree.addTrack(trkIt->first, idFromMass(fTrack1Mass));
    theTree.addTrack(trkIt->second, idFromMass(fTrack2Mass));
    theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.maxDoca), &(theTree.fTV.maxDocaV),   -1.,  fMaxDoca, Form("ditracks maxdoca %d", fType))); 
    theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.mass),    &(theTree.fTV.massV), fMassLow, fMassHigh, Form("ditracks mass %d", fType)));
    theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.pt),      &(theTree.fTV.ptV),    fCandPt,      999., Form("ditracks pt %d", fType))); 
    theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.flxy),    &(theTree.fTV.flxyV),      -1.,     fFlxy, Form("ditracks flxy %d", fType)));
    theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.flsxy),   &(theTree.fTV.flsxyV),  fFlsxy,      1.e9, Form("ditracks flsxy %d", fType)));
    theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.pvips),   &(theTree.fTV.pvipsV),     -1.,    fPvIpS, Form("ditracks pvips %d", fType)));
		
    fSequentialFitter->doFit(&theTree);
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
