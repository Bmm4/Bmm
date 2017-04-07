#include <memory>

#include "DataFormats/TrackReco/interface/Track.h"
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


#include "Bmm/RootAnalysis/common/JSON.hh"

#include <iostream>
#include <vector>

using namespace std;
using namespace edm;
using namespace reco;

class HFSkipEvents : public EDFilter {
public:
  explicit HFSkipEvents(const ParameterSet&);
  ~HFSkipEvents();

private:
  virtual void beginJob();
  virtual void endJob() ;
  virtual bool filter(Event&, const EventSetup&);

  int          fVerbose;

  int          filterOnPrimaryVertex;
  InputTag     fPrimaryVertexCollectionLabel;
  edm::EDGetTokenT<reco::VertexCollection> fTokenVertex;

  int          filterOnTracks;
  InputTag     fTrackCollectionLabel;

  int          filterOnMuons;
  InputTag     fMuonCollectionLabel;

  int          filterOnJSON;
  edm::FileInPath fJSONFile;
  JSON           *fpJSON;


  int fNpv, fNtk, fNmu, fNjson;
  int fNfailed, fNpassed;
  int fEvent;

};

// ----------------------------------------------------------------------
HFSkipEvents::HFSkipEvents(const edm::ParameterSet& iConfig):
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 0)),
  filterOnPrimaryVertex(iConfig.getUntrackedParameter<int>("filterOnPrimaryVertex", 0)),
  fPrimaryVertexCollectionLabel(iConfig.getUntrackedParameter<InputTag>("primaryVertexCollectionLabel", edm::InputTag("offlinePrimaryVertices"))),
  filterOnTracks(iConfig.getUntrackedParameter<int>("filterOnTrackMaximum", 0)),
  fTrackCollectionLabel(iConfig.getUntrackedParameter<InputTag>("TrackCollectionLabel", edm::InputTag("generalTracks"))),
  filterOnMuons(iConfig.getUntrackedParameter<int>("filterOnMuonMinimum", 0)),
  fMuonCollectionLabel(iConfig.getUntrackedParameter<InputTag>("MuonCollectionLabel", edm::InputTag("muons"))),
  filterOnJSON(iConfig.getUntrackedParameter<int>("filterOnJson", 0)),
  fJSONFile(iConfig.getUntrackedParameter<edm::FileInPath>("JSONFile"))
{
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFSkipEvents constructor" << endl;
  cout << "---  verbose:                         " << fVerbose << endl;
  cout << "---  filterOnPrimaryVertex:           " << filterOnPrimaryVertex << endl;
  cout << "---  primaryVertexCollectionLabel:    " << fPrimaryVertexCollectionLabel << endl;
  cout << "---  filterOnTrackMaximum:            " << filterOnTracks << endl;
  cout << "---  filterOnMuonMinimum:             " << filterOnMuons << endl;
  cout << "---  TrackCollectionLabel:            " << fTrackCollectionLabel << endl;
  cout << "---  filterOnJSON:                    " << filterOnJSON << endl;
  cout << "---  JSONFile:                        " << fJSONFile.fullPath() << endl;
  cout << "----------------------------------------------------------------------" << endl;

  fpJSON = new JSON(fJSONFile.fullPath().c_str(), 1);

  fNpv = fNtk = fNmu = fNjson = 0;
  fEvent = fNfailed = fNpassed = 0;

  fTokenVertex      = consumes<VertexCollection>(fPrimaryVertexCollectionLabel);

}


// ----------------------------------------------------------------------
HFSkipEvents::~HFSkipEvents() {

}



// ----------------------------------------------------------------------
bool HFSkipEvents::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  bool result(true);
  ++fEvent;

  int run = iEvent.id().run();
  int ls  = iEvent.luminosityBlock();

  bool goodJson(false);
  if (filterOnJSON > 0) {
    goodJson = fpJSON->good(run, ls);
    if (goodJson) {
      ++fNjson;
    } else {
      cout << "==> HFSkipEvent> rejecting run/ls = " << run << "/" << ls << endl;
      result = false;
    }
  } else {
    goodJson = true;
    ++fNjson;
  }

  edm::Handle<reco::VertexCollection> hVertices;
  int goodVertices(-1);
  try{
    iEvent.getByLabel(fPrimaryVertexCollectionLabel, hVertices);
    const reco::VertexCollection vertices = *(hVertices.product());
    goodVertices = 0;
    for (unsigned int i = 0; i < vertices.size(); ++i) {
      if (vertices[i].isFake()) {
      } else {
	++goodVertices;
      }
    }
  } catch (cms::Exception &ex) {
    if (fVerbose > 1) cout << "No Vertex collection with label " << fPrimaryVertexCollectionLabel << endl;
  }

  edm::Handle<MuonCollection> hMuons;
  int goodMuons(-1);
  try{
    iEvent.getByLabel(fMuonCollectionLabel, hMuons);
    goodMuons = (*(hMuons.product())).size();
  } catch (cms::Exception &ex)  {
    if (fVerbose > 1) cout << "No Muon collection with label " << fMuonCollectionLabel << endl;
  }
  bool goodMu(false);
  if (filterOnMuons > 0) {
    goodMu = (goodMuons >= filterOnMuons);
    if (goodMu) {
      //      result = true;
      ++fNmu;
    } else {
      result = false;
    }
  } else {
    goodMu = true;
    ++fNmu;
  }


  edm::Handle<reco::TrackCollection> hTracks;
  int goodTracks(-1);
  try{
    iEvent.getByLabel(fTrackCollectionLabel, hTracks);
    goodTracks = (*(hTracks.product())).size();
  } catch (cms::Exception &ex) {
    if (fVerbose > 1) cout << "No Track collection with label " << fTrackCollectionLabel << endl;
  }


  bool goodPv(false);
  if (filterOnPrimaryVertex > 0)  {
    goodPv = (goodVertices >= filterOnPrimaryVertex);
    if (goodPv) {
      ++fNpv;
    } else {
      result = false;
    }
  } else {
    goodPv = true;
    ++fNpv;
  }

  bool goodTk(false);
  if (filterOnTracks > 0) {
    goodTk = (goodTracks < filterOnTracks);
    if (goodTk) {
      ++fNtk;
    } else {
      result = false;
    }
  } else {
    ++fNtk;
    goodTk = true;
  }



  if (fVerbose > 0) {
    char line[20];
    sprintf(line, "%7d", fEvent);
    cout << "HFSkipEvents: " << line
	 << " " << iEvent.id().event()
	 << " result: " << (result? "true ":"false")
	 << " PV: " << (goodPv?1:0) << "(" << goodVertices << ")"
	 << " Tk: " << (goodTk?1:0) << "(" << goodTracks << ")"
	 << " Mu: " << (goodMu?1:0) << "(" << goodMuons << ")"
	 << " JSON: " << (goodJson?1:0) << "(run/ls = " << run << "/" << ls << ")"
	 << endl;
  }

  if (result) {
    ++fNpassed;
  } else {
    ++fNfailed;
  }

  return result;
}


// --------------------------------------------------------------------------------
void  HFSkipEvents::beginJob() {
}


// --------------------------------------------------------------------------------
void  HFSkipEvents::endJob() {

    cout << "HFSkipEvents: "
	 << " passed events: " << fNpassed
	 << " failed events: " << fNfailed
	 << " good fNpv: " << fNpv
	 << " fNtk: " << fNtk
	 << " fNmu: " << fNmu
	 << " fNjson: " << fNjson
	 << endl;

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFSkipEvents);
