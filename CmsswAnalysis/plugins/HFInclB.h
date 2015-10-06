#ifndef _HFINCLB_h_
#define _HFINCLB_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

class HFInclB : public HFVirtualDecay {
  
 public:
  explicit HFInclB(const edm::ParameterSet&);
  
 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();
  
  int           fNmuons;
};

#endif
