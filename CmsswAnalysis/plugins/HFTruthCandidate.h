#ifndef _HFTRUTHCANDIDATE_h_
#define _HFTRUTHCANDIDATE_h_

#include <memory>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class HFTruthCandidate : public edm::EDAnalyzer {

public:

  explicit HFTruthCandidate(const edm::ParameterSet&);
  ~HFTruthCandidate();
  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

private:

  edm::InputTag                             fTracksLabel;
  edm::EDGetTokenT<edm::View<reco::Track> > fTokenTrack;

  edm::InputTag                             fPrimaryVertexLabel;
  edm::EDGetTokenT<reco::VertexCollection>  fTokenVertex;

  edm::InputTag                             fBeamSpotLabel;
  edm::EDGetTokenT<reco::BeamSpot>          fTokenBeamSpot;

  edm::InputTag                             fMuonsLabel;
  edm::EDGetTokenT<reco::MuonCollection>    fTokenMuon;

  int                fType, fGenType, fMotherID;
  bool               fPartialDecayMatching;

  double             fMaxDoca;
  int                fVerbose; 
  
  int                fStableDaughters; 
  std::multiset<int> fDaughtersSet;
  std::multiset<int> fDaughtersGammaSet;
  std::multiset<int> fDaughtersGamma2Set;

  std::vector<int>   fDaughtersID;

  edm::ESHandle<TransientTrackBuilder> fTTB;

};

#endif
