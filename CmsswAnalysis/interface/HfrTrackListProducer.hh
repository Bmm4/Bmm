#ifndef HFRTRACKLISTPRODUCER_H
#define HFRTRACKLISTPRODUCER_H

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


#include <memory>
#include <vector>

class Trajectory;

class HfrTrackListProducer : public edm::stream::EDProducer<> {
public:

  void dumpConfiguration();
  explicit HfrTrackListProducer(const edm::ParameterSet&);

private:
  void beginRun(const edm::Run&,const edm::EventSetup&) override;
  void endRun(const edm::Run&,const edm::EventSetup&) override;

  ///Produce the PFRecTrack collection
  void produce(edm::Event&, const edm::EventSetup&) override;

  int fVerbose;

  edm::InputTag fTracksLabel;
  std::string   fTrackQualityString;
  edm::InputTag fPrimaryVertexLabel;
  edm::InputTag fBeamSpotLabel;
  edm::InputTag fMuonsLabel;
  std::string   fMuonQualityString;

  double fTrackMinPt;

  const MagneticField                 *fMagneticField;
  const reco::MuonCollection          *fMuonCollection;
  reco::VertexCollection               fVertexCollection;
  reco::BeamSpot                       fBeamSpot;

  edm::Handle<edm::View<reco::Track> > fTracksHandle;
  edm::ESHandle<TransientTrackBuilder> fTTB;

  edm::EDGetTokenT<reco::BeamSpot> fTokenBeamSpot;
  edm::EDGetTokenT<edm::View<reco::Track> > fTokenTrack;
  edm::EDGetTokenT<reco::MuonCollection>  fTokenMuon;
  edm::EDGetTokenT<reco::VertexCollection> fTokenVertex;


};
#endif
