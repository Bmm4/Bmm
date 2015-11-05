#ifndef CANDANA_H
#define CANDANA_H

#include <string>
#include <vector>
#include <map>

#include <TH1.h>
#include <TH2.h>

#include "inclbReader.hh"
#include "RedTreeData.hh"

class TTree; 
class TDirectory; 

// ----------------------------------------------------------------------
class candAna {
  
public:
  candAna(inclbReader *pReader, std::string name, std::string cutsFile);
  virtual ~candAna();

  virtual void        evtAnalysis(TAna01Event *evt);
  virtual void        candAnalysis();
  virtual void        endAnalysis();
  virtual void        setupRedTree(TTree *);
  virtual void        fillRedTreeData();
  void                triggerSelection();
  void                genAnalysis();
  
  virtual void        bookHist();

  virtual void        readCuts(std::string fileName, int dump = 1);
  virtual void        readFile(std::string fileName, std::vector<std::string> &lines);

  void                getSigTracks(std::vector<int> &v, TAnaCand *pC);
  std::string         splitTrigRange(std::string tl, int &r1, int &r2);
  bool                doTriggerMatching(TAnaTrack *fp1);

  std::string fName; 
  std::string fCutFile; 
  TDirectory *fHistDir; 
  inclbReader *fpReader; 
  TTree *fTree, *fAmsTree; 
  TAna01Event *fpEvt;
  TAnaTrack *fpSigTrack;
  TSimpleTrack *fpTrack;
  TAnaMuon *fpMuon;
  TGenCand *fpGenMuon;
  TAnaJet *fpSigJet;

  int fVerbose;
  int fIsMC;

  Long64_t fRun, fEvt;
  int fLS;
  int fEvent; 
  int fRunRange;
  int fYear;
  int fProcessType, fMuonProcessType;

  std::map<std::string, std::pair<int, int> > HLTRANGE;

  // -- cuts and setup
  int TYPE, NOPRESELECTION, IGNORETRIGGER, SELMODE; 
  int TRUTHCAND;
  double MUPTLO, MUPTHI; 
  double MUETALO, MUETAHI; 
  

  // -- variables for reduced tree, they are from fpCand
  bool    fDoFill;
  bool    fJSON;
  bool    fMuId, fGoodHLT;
  int     fHltType, fHltPs; 
  std::string fHLTPath;

  int     fType; 
  double  fMuPt, fMuEta, fMuPhi, fMuPtRel; 
  double  fJetPt, fJetEta, fJetPhi; 
  double  fGenMuPt, fGenMuEta, fGenMuPhi;
  struct RedTreeData fRTD;

};

#endif //  CANDANA_H
