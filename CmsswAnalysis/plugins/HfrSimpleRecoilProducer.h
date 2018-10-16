#ifndef HFRSIMPLERECOILPRODUCER_H
#define HFRSIMPLERECOILPRODUCER_H

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Bmm/CmsswAnalysis/plugins/HfrBaseProducer.h"

#include <memory>
#include <vector>

class HfrSimpleRecoilProducer : public HfrBaseProducer {
public:

  void dumpConfiguration();
  explicit HfrSimpleRecoilProducer(const edm::ParameterSet&);

private:
  void beginRun(const edm::Run&,const edm::EventSetup&) override;
  void endRun(const edm::Run&,const edm::EventSetup&) override;

  ///Produce the PFRecTrack collection
  void produce(edm::Event&, const edm::EventSetup&) override;
  void put(edm::Event&, std::unique_ptr<std::vector<int> >&, std::unique_ptr<int>&);
  double fVtxProb;
  int    fCandType;
};
#endif
