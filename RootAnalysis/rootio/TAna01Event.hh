#ifndef TANA01EVENT
#define TANA01EVENT


#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "TGenCand.hh"
#include "TSimpleTrack.hh"
#include "TAnaTrack.hh"
#include "TAnaMuon.hh"
#include "TTrgObj.hh"
#include "TAnaCand.hh"
#include "TAnaVertex.hh"
#include "TAnaJet.hh"

class TAna01Event : public TObject {

public:

  TAna01Event() { };
  TAna01Event(Int_t Option);
  virtual ~TAna01Event() { };
  virtual void  Clear(Option_t *option ="");
  void dump();

  // ----------------------------------------------------------------------

  // -- Generator block
  int                 nGenCands() {return fnGenCands;}
  TGenCand*           getGenCand(int n);
  virtual TGenCand*   addGenCand();
  void                dumpGenBlock();
  void                clearGenBlock();
  // For a given 'SimTrack', find the index in the generator block by doing a (p,id) matching
  int                 getGenIndex(double px, double py, double pz, int charge, double precision = 0.05);
  int                 getGenIndexWithDeltaR(double pt, double eta, double phi, double charge);
  int                 getGenIndexWithDeltaR(const TLorentzVector &tlv, double charge, double dRthrsh = 99999.9);
  int                 getGenIndexWithDeltaR(const TLorentzVector &tlv, const TVector3 &vtx, 
					    double charge, double dRthrsh = 99999.9, double dVtxRatioThrsh = 1.3);
  bool                isAncestor(TGenCand *pMother, TGenCand *pParticle); 

  // -- reduced/truncated generator block, filled by plugins/HFTruthCandidate.cc
  int                 nGenT() {return fnGenT;}
  TGenCand*           getGenTWithIndex(int n); // return TGenCand with fNumber = n
  TGenCand*           getGenT(int n);          // return TGenCand at position n
  virtual TGenCand*   addGenT();
  void                fillGenT(int n); // fill all daugthers and mothers of this genIdx into truncated gen blocj

  // -- RecTracks
  int                 nRecTracks() {return fnRecTracks;}
  TAnaTrack*          getRecTrack(int n);
  virtual TAnaTrack*  addRecTrack();

  // -- SimpleTracks
  int                    nSimpleTracks() {return fnSimpleTracks;}
  TSimpleTrack*          getSimpleTrack(int n);
  virtual TSimpleTrack*  addSimpleTrack();
  // -- these helper functions need to be here, and not in TSimpleTrack ...
  int                    getSimpleTrackMCID(int n);
  int                    getSimpleTrackMuonIdx(int n);
  TAnaMuon*              getSimpleTrackMuon(int n); 

  // -- Signal Tracks (e.g. refitted tracks forming cands in the next block)
  int                 nSigTracks() {return fnSigTracks;}
  TAnaTrack*          getSigTrack(int n);
  virtual TAnaTrack*  addSigTrack();

  // -- Signal Candidates (e.g. B_s)
  int                 nCands() {return fnCandidates;}
  TAnaCand*           getCand(int n);
  virtual TAnaCand*   addCand();

  // -- Muons
  int                 nMuons() {return fnMuons;}
  TAnaMuon*           getMuon(int n);
  virtual TAnaMuon*   addMuon();

  // -- Trigger Objects
  int                 nTrgObj() {return fnTrgObj;}
  TTrgObj*            getTrgObj(int n);
  virtual TTrgObj*    addTrgObj();

  // -- Primary vertices
  int                 nPV()    {return fnPV;}
  TAnaVertex*         bestPV() {return getPV(fBestPV); }
  TAnaVertex*         getPV(int n);
  virtual TAnaVertex* addPV();

  // -- Beamspot
  TAnaVertex          fBeamSpot; 

  // -- GenJets
  int                 nGenJets() {return fnGenJets;}
  TAnaJet*            getGenJet(int n);
  virtual TAnaJet*    addGenJet();

  // -- CaloJets
  int                 nCaloJets() {return fnCaloJets;}
  TAnaJet*            getCaloJet(int n);
  virtual TAnaJet*    addCaloJet();

  // -- TrackJets
  int                 nTrackJets() {return fnTrackJets;}
  TAnaJet*            getTrackJet(int n);
  virtual TAnaJet*    addTrackJet();


  // ----------------------------------------------------------------------
  // -- Basic event and detector information
  Long64_t          fRunNumber, fEventNumber;
  int               fEventBits;
  int               fDetectorStatus; 
  int               fEventTag;
  int               fBestPV;   // be careful with this...

  // -- MC event/generation information
  int               fProcessID;
  double            fXsec;
  double            fFilterEff;
  double            fEventWeight;
  double            fPtHat;

  // -- Lumi
  double            fLumi, fLumiInt; 
  int               fLumiSection; 
  int               fOrbit; 
  int               fBx; 
  
  unsigned int      fTimeLo, fTimeHi; 

  // -- Trigger words
  bool              fL1TDecision, fHLTDecision;
  #define NL1T 128
  #define NLTT 64
  #define NHLT 1024
  // -- L1 trigger
  TString           fL1TNames[NL1T];
  int               fL1TPrescale[NL1T];
  bool              fL1TResult[NL1T];
  bool              fL1TMask[NL1T];
  bool              fL1TError[NL1T];

  // -- L1 technical trigger
  TString           fLTTNames[NLTT];
  int               fLTTPrescale[NLTT];
  bool              fLTTResult[NLTT];
  bool              fLTTMask[NLTT];
  bool              fLTTError[NLTT];

  // -- HLT
  TString           fHLTNames[NHLT];
  int               fHLTPrescale[NHLT];
  bool              fHLTResult[NHLT];
  bool              fHLTWasRun[NHLT];
  bool              fHLTError[NHLT];

  // -- MET
  TVector3          fGenMET, fMET0, fMET1;  // only x and y component are relevant. z could contain type information. 

  // -- Reserve variables
  int               fIntRes[10]; 
  double            fDoubleRes[10]; 

private:

  int               fnGenCands;
  TClonesArray      *fGenCands;

  int               fnGenT;
  TClonesArray      *fGenT;

  int               fnSimpleTracks;
  TClonesArray      *fSimpleTracks;

  int               fnRecTracks;
  TClonesArray      *fRecTracks;

  int               fnMuons;
  TClonesArray      *fMuons;

  int               fnTrgObj;
  TClonesArray      *fTrgObj;

  int               fnSigTracks;
  TClonesArray      *fSigTracks;

  int               fnCaloJets;
  TClonesArray      *fCaloJets;

  int               fnGenJets;
  TClonesArray      *fGenJets;

  int               fnTrackJets;
  TClonesArray      *fTrackJets;

  int               fnCandidates;
  TClonesArray      *fCandidates;

  int               fnPV;
  TClonesArray      *fPV;

  ClassDef(TAna01Event,1)

};

#endif
