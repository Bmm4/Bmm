#ifndef _HFBD2JPSIKSTAR_h_
#define _HFBD2JPSIKSTAR_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

// ----------------------------------------------------------------------
class HFBd2JpsiKstar : public HFVirtualDecay {
 public:
  explicit HFBd2JpsiKstar(const edm::ParameterSet&);

 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();

  int           fPsiMuons;
  double        fPsiLo, fPsiHi, fPhiWindow, fBdWindow, fKstarLo, fKstarHi;

};

#endif
