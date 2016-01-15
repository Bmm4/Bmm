#ifndef _HFDUMPSTUFF_h_
#define _HFDUMPSTUFF_h_

#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
#include "DataFormats/Common/interface/ConditionsInEdm.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

class TFile;
class TTree;
class TAna01Event;

class TrackAssociatorBase;

// ----------------------------------------------------------------------
class HFDumpStuff : public edm::EDAnalyzer {
 public:
  explicit HFDumpStuff(const edm::ParameterSet&);
  ~HFDumpStuff();
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  int                  fVerbose; 
  std::string	       fCandidates1Label, fCandidates2Label, fCandidates3Label;
  edm::InputTag        fLumiSummaryLabel, fBeamSpotLabel, fPrimaryVertexLabel, fPrimaryVertexTracksLabel; 

  edm::EDGetTokenT<LumiSummary> fTokenLumiSummary;
  edm::EDGetTokenT<reco::BeamSpot> fTokenBeamSpot;
  edm::EDGetTokenT<std::vector<reco::Track> > fTokenTrack;
  edm::EDGetTokenT<reco::VertexCollection> fTokenVertex;
  
};

#endif
