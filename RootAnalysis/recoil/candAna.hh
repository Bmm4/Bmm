#ifndef CANDANA_H
#define CANDANA_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

#include <TH1.h>
#include <TH2.h>

// -- TMVA related
#include "TMVA/Reader.h"

#include "recoilReader.hh"
#include "common/AnalysisCuts.hh"

#include "cuts.hh"


#define NTRKMAX 1000

struct mvaMuonIDData {
  // -- Luca's original setup
  float trkValidFract, glbNChi2;
  float pt, eta;
  float segComp, chi2LocMom, chi2LocPos, glbTrackProb;
  float NTrkVHits, NTrkEHitsOut;
  // -- Stephan's other variables
  float glbTrackTailProb, glbDeltaEtaPhi, iValFrac;
  float LWH; //int
  float dxyRef, dzRef;
  float kinkFinder, glbKinkFinder, glbKinkFinderLOG, timeAtIpInOutErr, outerChi2;
  float valPixHits, TMTrkMult100; //int
  float innerChi2, trkRelChi2;
  float vMuonHitComb;
  float Qprod; //int
  float spectatorDummy;
};


class TTree;
class TDirectory;

// ----------------------------------------------------------------------
class candAna {

public:
  candAna(recoilReader *pReader, std::string name, std::string cutsFile);
  virtual ~candAna();

  virtual void          evtAnalysis(TAna01Event *evt);
  virtual void          readCuts(std::string fileName, int dump = 1);
  virtual void          printCuts(ostream &OUT);
  virtual void          readFile(std::string fileName, std::vector<std::string> &lines);
  virtual TMVA::Reader* setupMuonMvaReader(std::string xmlFile, mvaMuonIDData &rd);
  virtual void          setupReducedTree(TTree *);

  virtual void          genMatch();
  virtual void          recoMatch();
  virtual void          candMatch();

  virtual void          resetAllData();
  virtual void          genAnalysis();
  virtual void          brecoAnalysis();
  virtual void          candAnalysis();
  virtual void          candEvaluation();
  virtual void          endAnalysis();

  virtual void          bookHist();

  virtual void          play();

  TH1D*                 getHist(std::string);

  std::pair<TVector3, TVector3> parallelAndPerp(TVector3 direction, TVector3 momVis);
  std::pair<TVector3, TVector3> parallelAndPerp2(TVector3 direction, TVector3 momVis);
  std::pair<double, double>  nuRecoMom0(double compVisPar, double compVisPerp, double eVis, double mVis, double mTot = 1.777);
  double      minMassPair(std::vector<TLorentzVector> a);
  double      maxMassPair(std::vector<TLorentzVector> a);
  // -- calculate opening angle between three vectors (in fact the sine of half opening angle)
  double      oaTriplet(std::vector<TLorentzVector> a);


  // -- general setup
  std::string fName;
  std::string fCutFile;
  TDirectory *fHistDir;
  recoilReader *fpReader;
  TTree *fTree, *fEffTree;
  TAna01Event *fpEvt;
  TAnaCand *fpCand, *fpCandTruth, *fpBrecoCand, *fpRecoilCand, *fpRecoilMissCand, *fpOsCand;

  // -- MVA muon ID
  TMVA::Reader *fMvaMuonID;
  mvaMuonIDData mrd;
  double  fMuBDT;

  int fVerbose, fDbx;

  AnalysisCuts fAnaCuts;
  std::vector<cuts*> fCuts;

  // -- "Cuts"
  int CANDTYPE, CANDTRUTH, DATACAND, BLIND, NOPRESELECTION;

  double       MASSMIN,   MASSMAX;
  double       SIGBOXMIN, SIGBOXMAX;
  double       BGLBOXMIN, BGLBOXMAX;
  double       BGHBOXMIN, BGHBOXMAX;

  double TRACKPTLO, TRACKPTHI, TRACKETALO, TRACKETAHI
    , TRACKTIP, TRACKLIP
    , MUPTLO, MUPTHI
    , MUETALO, MUETAHI, MUIP, MUBDT
    ;


  // -- BRECO analysis
  int         BRECOTYPE;

  double      fBrecoMass, fBrecoPt, fBrecoEta, fBrecoPhi, fBrecoPvZ;
  int         fBrecoPvIdx;
  double      fMu1BrecoPt, fMu1BrecoEta, fMu1BrecoPhi, fMu2BrecoPt, fMu2BrecoEta, fMu2BrecoPhi;

  // -- Signal cand analysis
  int         fCandPvI, fBrecoPvI;
  int         fMissedTracks;

  // -- TM
  int     fGenBTmi;
  int     fCandTmi;  // index of truth-matched candidate
  bool    fRecoTm;   // flag whether the final state particles have reconstructed (simple) tracks
  bool    fTm;       // flag whether the current candidate is truth-matched
  bool    fTruthRm;  // flag whether the current candidate is reco-matched (for truth-identified cands)
  int     fNGenPhotons;
  int     fGenBpartial;
  int     fProcessType;

  TGenCand    *fpGenB;
  TVector3    fVtxBProd;   // from truth-matching
  TVector3    fPV;         // from fPvIdx
  TVector3    fVtxBDecay;

  // -- Variables
  bool   fPreselection;
  int    fNchan;
  int    fPvN;
  double fPvX, fPvY, fPvZ;

  // -- inherited from t1Reader/recoilReader
  int fIsMC;
  Long64_t fRun, fEvt;
  int fLS, fChainEvent;
  int fYear;
  std::string fEra;
  double fLumi;
  bool fJSON;

};

#endif //  CANDANA_H
