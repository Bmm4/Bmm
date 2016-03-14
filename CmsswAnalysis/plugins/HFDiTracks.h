#ifndef _HFDITRACKS_h_
#define _HFDITRACKS_h_

#include <string>

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

class HFDiTracks : public HFVirtualDecay {
 public:
  explicit HFDiTracks(const edm::ParameterSet&);
  
 protected:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void dumpConfiguration();
  
  virtual int  idFromMass(double mass);
  
  double fLeadingTrackPt; 
  double fTrack1Mass;
  double fTrack2Mass;
  double fExtra;
  
  int    fNbrMuons;
  bool   fCloseToMuons;
};

#endif
