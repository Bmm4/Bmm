#ifndef TANAMUON
#define TANAMUON


#include "TObject.h"
#include "TArrayC.h"
#include "TArrayF.h"
#include "TVector3.h"

#include "TAnaTrack.hh"

#include <set>
#include <map>
#include <vector>

class xpTrack {
public: 
  xpTrack(int i = -1) {idx = i;}
  virtual ~xpTrack() {}
  void clear() {idx = -1;}
  int idx; 
  float dist; 
  TVector3 r, p;
  ClassDef(xpTrack,1)
};


class TAnaMuon : public TAnaTrack {

public:

  TAnaMuon(int index = -1);
  ~TAnaMuon();
  void clear(); 
  void dump();

  // ----------------------------------------------------------------------
  TVector3  fInnerPlab, fOuterPlab, fGlobalPlab;
  TVector3  fPositionAtM1, fPositionAtM2;
  TVector3  fMuonTrackPosAtM1, fMuonTrackPlabAtM1;

  // see AN2008_098 "Muon identification in CMS"
  float     fMuonChi2, fMuonZ, fMuonR;
  float     fCaloComp, fSegmComp, fIsolation; 
  int       fNvalidMuonHits, fNhitsCSC, fNhitsDT,fNhitsRPC; 
  int       fNmatchedStations;
  float     fTimeInOut, fTimeInOutE, fTimeOutIn, fTimeOutInE; 
  int       fTimeNdof;

  // -- new variables for MVA muon id
  float  fItrkValidFraction; 
  float  fGtrkNormChi2;    
  float  fChi2LocalPosition, fChi2LocalMomentum; 
  int    fNumberOfLostTrkHits, fNumberOfValidTrkHits; 
  float  fSegmentComp;    
  float  fGtrkProb;           
 
  // information on tracks near this candidate
  std::map<int,float> fNstTracks; // usage: (trackIx, doca)
  // information on tracks extrapolated to muon station 1 (position and momentum): tkIdx, dist to muon, tk pos, tk mom
  static const unsigned int NXPTRACKS = 10; 
  xpTrack fXpTracks[NXPTRACKS];
  // information on best vertex of muon
  double fVtxProb;
  std::set<unsigned> fVtxTracks;

  // -- reserve
  int            fInt1,    fInt2,    fInt3; 
  double         fDouble1, fDouble2, fDouble3;

  int STATrkMult_150, TMTrkMult_100;
  int fNumberOfValidPixHits;
  double finnerChi2;
  double fouterChi2;
  std::vector<int> fvDThits;
  std::vector<int> fvRPChits;
  std::vector<int> fvCSChits;
  float fglbKinkFinder;
  float ftrkRelChi2;
  float fstaRelChi2;
  float fglbDeltaEtaPhi;
  //muonID BDT
  double barrelBDTresponse, endcapBDTresponse;

private:

  ClassDef(TAnaMuon,1)

};

#endif
