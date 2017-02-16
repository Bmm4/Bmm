#ifndef _HFBS2JPSIF0_h_
#define _HFBS2JPSIF0_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

// ----------------------------------------------------------------------
class HFBs2Jpsif0 : public HFVirtualDecay {
 public:
  explicit HFBs2Jpsif0(const edm::ParameterSet&);

 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();

  int           fPsiMuons;
  double        fPsiLo, fPsiHi, ff0Lo, ff0Hi;
};

#endif
