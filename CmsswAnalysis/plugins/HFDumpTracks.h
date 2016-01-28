#ifndef _HFDUMPTRACKS_h_
#define _HFDUMPTRACKS_h_

#include <memory>

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

#define HFMAXTRACK 10000

class HFDumpTracks : public HFVirtualDecay {
  
 public:
  explicit HFDumpTracks(const edm::ParameterSet&);
  ~HFDumpTracks();

 private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();
  
  int                  fDoTruthMatching;
  bool                 fDumpSimpleTracks, fDumpRecTracks;

};

#endif
