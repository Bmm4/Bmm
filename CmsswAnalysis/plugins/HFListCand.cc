#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFListCand.h"

#include <iostream>
#include <utility>

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"

#include "Bmm/CmsswAnalysis/interface/HFSequentialVertexFit.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"
#include "Bmm/CmsswAnalysis/interface/HFDecayTree.hh"
#include "Bmm/CmsswAnalysis/interface/HFTrackListBuilder.hh"
#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

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
HFListCand::HFListCand(const ParameterSet& iConfig) :
  HFVirtualDecay(iConfig) {
  dumpConfiguration();
}

// ----------------------------------------------------------------------
void HFListCand::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFListCand configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFListCand::analyze(const Event& iEvent, const EventSetup& iSetup) {
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;

  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFListCand> " << e.fMsg << endl;
    return;
  }

  fListBuilder->setCallerName("HFListCand");
  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setMaxDocaToTracks(fMaxDoca);
  vector<int> trkList = fListBuilder->getTrackList();

  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fType = fType;
  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + trkList.size() - 1;
  TAnaTrack *pTrack(0);

  TLorentzVector p4t, p4tot;
  p4tot.SetXYZT(0., 0., 0., 0.);
  for (vector<int>::const_iterator trkIt = trkList.begin(); trkIt != trkList.end(); ++trkIt) {
    int idx = *trkIt;
    //    cout << "idx = " << ix << " -> track idx = " << idx << " pv idx = " << pvidx << endl;
    TSimpleTrack *sTrack = gHFEvent->getSimpleTrack(idx);
    pTrack = gHFEvent->addSigTrack();
    TrackBaseRef baseRef(fTracksHandle, idx);
    int gidx = sTrack->getGenIndex();
    Track trackView(*baseRef);
    const reco::BeamSpot *pBeamSpot = &fBeamSpot;
    fillAnaTrack(pTrack, trackView, idx, gidx, &fVertexCollection, fMuonCollection, pBeamSpot);
    p4t.SetXYZM(trackView.px(), trackView.py(), trackView.pz(), MPION);
    p4tot += p4t;
  }
  pCand->fMass = p4tot.M();
  pCand->fPlab = p4tot.Vect();

  cout << "HFListCand dump: " << endl;
  pCand->dump();
  for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
    gHFEvent->getSigTrack(is)->dump();
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFListCand);
