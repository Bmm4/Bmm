#ifndef HFSIMPLELISTPRODUCER_H
#define HFSIMPLELISTPRODUCER_H

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Bmm/CmsswAnalysis/plugins/HfrBaseProducer.h"

#include <memory>
#include <vector>

// ----------------------------------------------------------------------
// Create a track list for a candidate
// - Takes from gHFEvent a candidate 'type'
// - produces a track list of the candidate's tracks
// - put the list into the event
// ----------------------------------------------------------------------


class HFSimpleListProducer : public HfrBaseProducer {
public:

  void dumpConfiguration();
  explicit HFSimpleListProducer(const edm::ParameterSet&);

 private:
  void beginRun(const edm::Run&,const edm::EventSetup&) override;
  void endRun(const edm::Run&,const edm::EventSetup&) override;

  void produce(edm::Event&, const edm::EventSetup&) override;
  void put(edm::Event&, std::unique_ptr<std::vector<int> >&, std::unique_ptr<int>&);
};
#endif
