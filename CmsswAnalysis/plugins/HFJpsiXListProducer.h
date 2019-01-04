#ifndef HFRJPSIXLISTPRODUCER_H
#define HFRJPSIXLISTPRODUCER_H

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Bmm/CmsswAnalysis/plugins/HfrBaseProducer.h"

#include <memory>
#include <vector>

class HfrJpsiXListProducer : public HfrBaseProducer {
public:

  void dumpConfiguration();
  explicit HfrJpsiXListProducer(const edm::ParameterSet&);

private:
  void beginRun(const edm::Run&,const edm::EventSetup&) override;
  void endRun(const edm::Run&,const edm::EventSetup&) override;

  ///Produce the PFRecTrack collection
  void produce(edm::Event&, const edm::EventSetup&) override;
  void put(edm::Event&, std::unique_ptr<std::vector<int> >&, std::unique_ptr<int>&, std::unique_ptr<reco::VertexCollection>&);
  double fVtxProb;

};
#endif
