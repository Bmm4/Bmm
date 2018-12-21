#ifndef _HFRTEST_h_
#define _HFRTEST_h_

#include <string>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "Bmm/CmsswAnalysis/interface/HFSequentialVertexFit.hh"

struct HFSetupException {
HFSetupException(const char* msg) : fMsg(msg) {}
  std::string fMsg;
};


class HfrTest : public edm::EDAnalyzer {
 public:
  explicit HfrTest(const edm::ParameterSet&);

 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();

  virtual int  idFromMass(double mass);

  int     fVerbose;
  edm::InputTag fTracksLabel;
  edm::InputTag fPVLabel;
  edm::InputTag fBeamSpotLabel;
  edm::InputTag fMuonsLabel;
  std::string   fMuonQualityString;

  edm::InputTag fTrkIdxLabel;
  edm::InputTag fPvIdxLabel;
  edm::InputTag fVertexLabel;
  double fLeadingTrackPt;

  edm::Handle<edm::View<reco::Track> > fTracksHandle;
  edm::EDGetTokenT<edm::View<reco::Track> > fTokenTrack;
  edm::EDGetTokenT<reco::VertexCollection> fTokenPV;
  edm::EDGetTokenT<reco::BeamSpot> fTokenBeamSpot;
  edm::EDGetTokenT<reco::MuonCollection>  fTokenMuon;

  edm::Handle<std::vector<int> > fTrkIdxHandle;
  edm::EDGetTokenT<std::vector<int> > fTokenTrkIdx;

  edm::Handle<int> fPvIdxHandle;
  edm::EDGetTokenT<int> fTokenPvIdx;

  const reco::VertexCollection        *fPVCollection;
  const reco::MuonCollection          *fMuonCollection;
  const reco::BeamSpot                *fBeamSpot;

  const MagneticField                 *fMagneticField;
  const TransientTrackBuilder         *fTTB;

};

#endif
