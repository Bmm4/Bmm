#ifndef HFJPSIXLISTPRODUCER_H
#define HFJPSIXLISTPRODUCER_H

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "Bmm/CmsswAnalysis/plugins/HfrBaseProducer.h"

#include <memory>
#include <vector>

class HFJpsiXListProducer : public HfrBaseProducer {
public:

  void dumpConfiguration();
  explicit HFJpsiXListProducer(const edm::ParameterSet&);

private:
  void beginRun(const edm::Run&,const edm::EventSetup&) override;
  void endRun(const edm::Run&,const edm::EventSetup&) override;
  RefCountedKinematicTree fitTree(std::vector<reco::TransientTrack> &);

  ///Produce the PFRecTrack collection
  void produce(edm::Event&, const edm::EventSetup&) override;
  void put(edm::Event&, std::unique_ptr<std::vector<int> >&, std::unique_ptr<int>&, std::unique_ptr<reco::VertexCollection>&);
  double fVtxProb;

};
#endif
