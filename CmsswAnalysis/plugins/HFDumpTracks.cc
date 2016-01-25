// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFDumpTracks
// ------------
//
// 2016/01/20 Urs Langenegger      derive from HFVirtualDecay, migrate to "consumes"
// stone age  Urs Langenegger      first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFDumpTracks.h"
#include "HFDumpMuons.h"

#include <iostream>

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h" 
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/TrackerRecHit2D/interface/SiTrackerGSRecHit2D.h" 

#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TSimpleTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"
#include "Bmm/RootAnalysis/rootio/TAnaMuon.hh"
#include "Bmm/RootAnalysis/rootio/TTrgObj.hh"

#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

#include <TFile.h>
#include <TH1.h>

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace edm;
using namespace reco;


// ----------------------------------------------------------------------
HFDumpTracks::HFDumpTracks(const ParameterSet &iConfig) : 
  HFVirtualDecay(iConfig),
  fDoTruthMatching(iConfig.getUntrackedParameter<int>("doTruthMatching", 1)),
  fDumpSimpleTracks(iConfig.getUntrackedParameter<bool>("dumpSimpleTracks", true)),
  fDumpRecTracks(iConfig.getUntrackedParameter<bool>("dumpRecTracks", false))  {
  dumpConfiguration();
}


// ----------------------------------------------------------------------
void HFDumpTracks::dumpConfiguration() {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDumpTracks configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  dumpSimpleTracks:        " << fDumpSimpleTracks << endl;
  cout << "---  dumpRecTracks:           " << fDumpRecTracks << endl;
  cout << "---  doTruthMatching:         " << fDoTruthMatching << endl;  // 0 = nothing, 1 = TrackingParticles, 2 = FAMOS, 3 = TAna01Event
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
HFDumpTracks::~HFDumpTracks() {
  
}


// ----------------------------------------------------------------------
void HFDumpTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (fVerbose > 0) cout << "==>HFDumpTracks> new event " << endl;

  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFfDumpTracks> " << e.fMsg << endl;
    return;
  }

  fListBuilder->setMinPt(-1.);
  
  
  int genIdx(-1); 
  
  vector<int> muonList = fListBuilder->getMuonList();
  vector<int> trkList  = fListBuilder->getTrackList();
  for (unsigned int i = 0; i < fTracksHandle->size(); ++i) {
    TrackBaseRef rTrackView(fTracksHandle, i);
    Track track(*rTrackView);
		    
    // -- Muon?
    int mid = 0; 
    for (unsigned int im = 0; im < muonList.size(); ++im) {
      if (i == static_cast<unsigned int>(muonList[im])) {
	mid = 1; 
	break;
      }
    }
    
    // -- truth matching with TAna01Event::getGenIndexWithDeltaR(...)
    genIdx = -1; 
    if (3 == fDoTruthMatching) {    
      genIdx  = gHFEvent->getGenIndexWithDeltaR(track.pt(), track.eta(), track.phi(), track.charge()); 
    }
    
    // -- fill the tracks
    if (fDumpSimpleTracks) {
      TSimpleTrack *st = gHFEvent->addSimpleTrack();
      fillSimpleTrack(st, track, i, mid, genIdx, &fVertexCollection); 
    } 

    if (fDumpRecTracks) {
      TAnaTrack *at = gHFEvent->addRecTrack();
      fillAnaTrack(at, track, i, genIdx, &fVertexCollection, fMuonCollection, &fBeamSpot); 
    }
  }

  cleanupTruthMatching(fTracksHandle, fMagneticField);

  if (fVerbose > 0) {
    if (fDumpSimpleTracks) {
      for (int it = 0; it < gHFEvent->nSimpleTracks(); ++it) {
	gHFEvent->getSimpleTrack(it)->dump();
      }
    }

    if (fDumpRecTracks) {
      for (int it = 0; it < gHFEvent->nRecTracks(); ++it) {
	gHFEvent->getRecTrack(it)->dump();
      }
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDumpTracks);
