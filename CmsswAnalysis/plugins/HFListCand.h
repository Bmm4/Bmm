#ifndef _HFLISTCAND_h_
#define _HFLISTCAND_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

class HFListCand : public HFVirtualDecay {

 public:
  explicit HFListCand(const edm::ParameterSet&);

 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();

};

#endif
