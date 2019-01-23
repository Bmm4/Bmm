#ifndef _HFBU2MUTAUK_h_
#define _HFBU2MUTAUK_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

class HFBu2MuTauK : public HFVirtualDecay {

 public:
  explicit HFBu2MuTauK(const edm::ParameterSet&);

 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();

  double fMuKaVtxProb;
  double fTauVtxProb;

};

#endif
