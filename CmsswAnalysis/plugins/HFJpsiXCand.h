#ifndef _HFJPSIXCAND_h_
#define _HFJPSIXCAND_h_

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"


class HFJpsiXCand : public HFVirtualDecay {
 public:
  explicit HFJpsiXCand(const edm::ParameterSet&);

 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  //  virtual void addTrack(int idx);
  virtual void dumpConfiguration();

  int    fRefType;

};

#endif
