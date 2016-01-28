#ifndef HFVIRTUALDECAY_H
#define HFVIRTUALDECAY_H

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

#include "Bmm/CmsswAnalysis/interface/HFTrackListBuilder.hh"
#include "Bmm/CmsswAnalysis/interface/HFSequentialVertexFit.hh"

struct HFSetupException {
HFSetupException(const char* msg) : fMsg(msg) {}
  std::string fMsg;
};


// ----------------------------------------------------------------------
class HFVirtualDecay : public edm::EDAnalyzer {
	
 public:
  explicit HFVirtualDecay(const edm::ParameterSet& iConfig);
	
 protected:
  virtual void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup);
		
  virtual void dumpConfiguration();
 protected:
  int fVerbose;
		
  edm::InputTag fTracksLabel;
  std::string   fTrackQualityString;
  edm::InputTag fPrimaryVertexLabel;
  edm::InputTag fBeamSpotLabel;
  edm::InputTag fMuonsLabel;
  std::string   fMuonQualityString;
		
  double fTrackPt;
  double fMuonPt;
  double fDeltaR;
  double fMaxDoca;
  double fMaxD0;
  double fMaxDz;
  double fPvWeight;
		
  int fType;
	
 protected:
  const MagneticField                 *fMagneticField;
  const reco::MuonCollection          *fMuonCollection;
  reco::VertexCollection               fVertexCollection;
  reco::BeamSpot                       fBeamSpot;
		
  edm::Handle<edm::View<reco::Track> > fTracksHandle;
  edm::ESHandle<TransientTrackBuilder> fTTB;
		
  std::auto_ptr<HFTrackListBuilder>    fListBuilder;
  std::auto_ptr<HFSequentialVertexFit> fSequentialFitter;

  edm::EDGetTokenT<reco::BeamSpot> fTokenBeamSpot;
  edm::EDGetTokenT<edm::View<reco::Track> > fTokenTrack;
  edm::EDGetTokenT<reco::MuonCollection>  fTokenMuon;
  edm::EDGetTokenT<reco::VertexCollection> fTokenVertex;

};

#endif
