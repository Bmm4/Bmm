#ifndef _HFBU2MUTAUK_h_
#define _HFBU2MUTAUK_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "Bmm/RootAnalysis/rootio/TGenCand.hh"

class HFBu2MuTauK : public HFVirtualDecay {

 public:
  explicit HFBu2MuTauK(const edm::ParameterSet&);

 protected:
  void analyze(const edm::Event&, const edm::EventSetup&);
  void dumpConfiguration();
  // -- find from gHFEvent the pointers to the charged stable final state particles
  void genMatch();
  // -- produce a vector of int that matches the list of truth-identified tracks
  std::vector<int> truthIndices();
  // -- starting from a vector of int (mu, ka, had1, had2, had3) build a candidate of type
  TAnaCand* buildCandidate(std::vector<int>, int type, bool applyCuts = true);
  // -- write a candidate to the tree (if all cuts are passing)
  TAnaCand* fillCand(int type, RefCountedKinematicTree &tree, std::vector<int> trkIdx, std::vector<int> trkId, reco::Vertex &vtx);
  // -- given something, calculate IP(3D) wrt a vertex
  std::pair<bool, Measurement1D> doca2Vtx(reco::TransientTrack &, reco::Vertex &);
  std::pair<bool, Measurement1D> doca2Vtx(RefCountedKinematicTree &, reco::Vertex &);
  std::pair<bool, Measurement1D> doca2Vtx(TrajectoryStateOnSurface &, reco::Vertex &);
  int    fMcType;
  double fMuKaVtxProb, fMuKaMaxDoca;
  double fTauVtxProb, fTauMaxDoca, fTauMaxIP;

  TGenCand *fpGenB, *fpGenMu, *fpGenKa, *fpGenTau, *fpGenHad1, *fpGenHad2, *fpGenHad3;

};

#endif
