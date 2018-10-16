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
  virtual void          candAnalysis();
  virtual void          candEvaluation();
  virtual void          endAnalysis();

  virtual void          genMatch();
  virtual void          recoMatch();
  virtual void          candMatch();

  virtual void          bookHist();

  virtual void        play();

  // -- general setup
  std::string fName;
  std::string fCutFile;
  TDirectory *fHistDir;
  recoilReader *fpReader;
  TTree *fTree, *fEffTree;
  TAna01Event *fpEvt;
  TAnaCand *fpCand, *fpOsCand;

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

  // -- TM
  int     fGenBTmi;
  int     fNGenPhotons;
  int     fCandTmi;
  int     fGenBpartial;
  int     fProcessType;

  // -- Variables
  int fPreselection, fNchan;
  int fCandIdx;

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
