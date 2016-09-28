#ifndef _HFDUMPMUONS_h_
#define _HFDUMPMUONS_h_

#include <memory>

#include "Bmm/CmsswAnalysis/plugins/HFVirtualDecay.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"


#include "Bmm/RootAnalysis/rootio/TAnaMuon.hh"
#include "Bmm/CmsswAnalysis/interface/HFBDT.hh" 


// ----------------------------------------------------------------------
class HFDumpMuons : public HFVirtualDecay {
 public:
  explicit                  HFDumpMuons(const edm::ParameterSet&);
  ~HFDumpMuons();
  
 private:
  virtual void              analyze(const edm::Event&, const edm::EventSetup&);
  virtual void              dumpConfiguration();
  virtual void              beginRun(const edm::Run&, const edm::EventSetup&);

  void                      fillMuon(const reco::Muon& tr, int type);
  void                      fillCaloMuon(const reco::CaloMuon& tr, int type);
  void                      extrapolateTracks(); 
  bool                      doExtrapolate(double pt, double eta);
  void                      findVertex(TAnaMuon *anaMu, std::set<unsigned> *trkIcs, double *prob);

  int                       getTrackMultiplicity(const reco::Muon& muon, bool STA, bool verbose);
  double                    getDistanceM1(const reco::Muon& muon, const reco::Muon& sec, bool MuType, bool verbose);
  bool                      tracksAreEqual(const reco::TrackRef& main,const reco::TrackRef& STA);
  void                      fillMuonDetHits(TAnaMuon* pM, reco::TrackRef gTrack);


  edm::InputTag             fCaloMuonsLabel;
  
  double                    fMaxTrackDistToStore;
  double                    fDocaVertex;           // try vertexing only with tracks closer than this
  unsigned                  fKeepBest;             // number of candidates to keep for iterative vertex search
  unsigned                  fMaxCandTracks;        // max number of tracks for the muon candidate vertex

  PropagateToMuon           fpropM1, fpropM2, fOutwardPropM1, fInwardPropM1;
  std::vector<xpTrack>      fXpTracks;

  HFBDT barrelBDT,endcapBDT;
  edm::FileInPath           fweightFileBarrel;
  edm::FileInPath           fweightFileEndcap;

  edm::EDGetTokenT<reco::CaloMuonCollection> fTokenCaloMuon;

};

#endif
