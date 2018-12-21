#ifndef HFRTRACKLISTPRODUCER_H
#define HFRTRACKLISTPRODUCER_H

#include "Bmm/CmsswAnalysis/plugins/HfrBaseProducer.h"

#include <memory>
#include <vector>

class HfrTrackListProducer : public HfrBaseProducer {
public:

  void dumpConfiguration();
  explicit HfrTrackListProducer(const edm::ParameterSet&);

private:
  void beginRun(const edm::Run&,const edm::EventSetup&) override;
  void endRun(const edm::Run&,const edm::EventSetup&) override;

  ///Produce the PFRecTrack collection
  void produce(edm::Event&, const edm::EventSetup&) override;

  int fVerbose;

};
#endif
