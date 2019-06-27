#ifndef _HFRTEST_h_
#define _HFRTEST_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"


class HfrTest : public HFVirtualDecay {
 public:
  explicit HfrTest(const edm::ParameterSet&);

 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();

};

#endif
