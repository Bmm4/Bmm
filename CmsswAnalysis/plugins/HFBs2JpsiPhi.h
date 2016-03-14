#ifndef _HFBS2JPSIPHI_h_
#define _HFBS2JPSIPHI_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

// ----------------------------------------------------------------------
class HFBs2JpsiPhi : public HFVirtualDecay {
 public:
  explicit HFBs2JpsiPhi(const edm::ParameterSet&);
  
 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();
  
  int           fPsiMuons;
  double        fPsiLo, fPsiHi, fPhiLo, fPhiHi;
};

#endif
