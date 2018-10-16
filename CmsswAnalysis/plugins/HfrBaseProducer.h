#ifndef HFRBASEPRODUCER_H
#define HFRBASEPRODUCER_H

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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"


#include <memory>
#include <vector>

class HfrBaseProducer : public edm::stream::EDProducer<> {
public:

  void dumpConfiguration();
  explicit HfrBaseProducer(const edm::ParameterSet&);

protected:
  void beginRun(const edm::Run&,const edm::EventSetup&) ;
  void endRun(const edm::Run&,const edm::EventSetup&) ;

  void analyze(edm::Event&, const edm::EventSetup&) ;
  void produce(edm::Event&, const edm::EventSetup&) ;

  int fVerbose;

  edm::InputTag fTracksLabel;
  std::string   fTrackQualityString;
  edm::InputTag fPrimaryVertexLabel;
  edm::InputTag fBeamSpotLabel;
  edm::InputTag fMuonsLabel;
  std::string   fMuonQualityString;

  double fTrackMinPt;
  double fMuonMinPt;

  const MagneticField                 *fMagneticField;
  const reco::MuonCollection          *fMuonCollection;
  reco::VertexCollection               fVertexCollection;
  reco::BeamSpot                       fBeamSpot;

  edm::Handle<edm::View<reco::Track> > fTracksHandle;

  edm::EDGetTokenT<reco::BeamSpot> fTokenBeamSpot;
  edm::EDGetTokenT<edm::View<reco::Track> > fTokenTrack;
  edm::EDGetTokenT<reco::MuonCollection>  fTokenMuon;
  edm::EDGetTokenT<reco::VertexCollection> fTokenVertex;

  const TransientTrackBuilder         *fTTB;

};
#endif
