#ifndef CANDANADSTAR_H
#define CANDANADSTAR_H

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include "candAna.hh"

class candAnaDstar : public candAna {
  
public:
  candAnaDstar(bmmReader *pReader, std::string name, std::string cutsFile);
  ~candAnaDstar();

  //void        evtAnalysis(TAna01Event *evt);
  bool        anaMC(bool);
  void        candAnalysis();
  void        moreBasicCuts();

  int         truthMatch(TAnaCand *pC, int verbose = 0); 
  void        dumpHFTruthCand(TAnaCand *pC); 
  void        dumpHFDstarCand(TAnaCand *pC); 
  void        readCuts(string filename, int dump);
  void        bookHist();
  void        genMatch();
  void        recoMatch();
  void        candMatch();
  void        efficiencyCalculation();
  int         doTest(TAnaCand *pC, int mode =-1); 
  int         getJpsi(int &idx1, int &idx2);
  bool        matchToTriggeredMuon(TAnaTrack *pt, int &idx);
  int         doMuonTriggerMatching(void); // match all offline muons to hlt muons 
  double      doTriggerMatchingTest(int &idx, int muon, int trig = 0); // match a single track to HLT
  double      findBestMatchingTrack(TLorentzVector,int,int &); // match trigger track to all reco tracks
  int         getHLTId(TString); // get the id (integer) from the HLT name 
  void        dumpAll(); // print everything
  bool        doTriggerMatchingForDs(double &dr1, double &dr2); // match hlt muons to offline 
  bool        doTriggerVeto(TAnaTrack *pt, bool singleMatch,
			    bool matchPt, bool anyModule, float drCut, int histoOffset=0); 
  // match the 2 muons from the dimuon to HLT
  //virtual bool        doTriggerMatching(TAnaTrack *pt1, TAnaTrack *pt2);
  // match a single track to HLT
  //virtual bool        doTriggerMatching(TAnaTrack *pt, bool anyTrig = false,
  //					bool muonsOnly=true, bool anyModule=false);
  // To return the full deltaR not just a bool
  virtual double      doTriggerMatchingR(TAnaTrack *pt, bool anyTrig = false,
					 bool muonsOnly=true, bool anyModule=false);
  // Veto track too close to a trigger object
  //virtual bool        doTriggerVeto(TAnaTrack *pt, bool singleMatch,
  //				    bool matchPt, bool anyModule, float drCut, int histoOffset=0);



  bool        analyzeHltInfo(bool allTrigClean=false);
  void        candEvaluation();
  void        triggerHLT();

private:
  TTree * tree;

  // reduced tree variables 
  int ftm, fnclose;
  bool fmuid1, fmuid2, fmumat1, fmumat2;
  float fmds, fmdz;
  float ffls3d,fchi2,falpha,fqpis,fdr;
  float fpt,fptdz,fptpis,fptpi,fptk;
  float fpvd, fiso;
  float feta, fetapi, fetak;
  float fchipi, fchik;
  float fmudr1, fmudr2;  
  float fmatch1dr, fmatch2dr;
  bool fmuidmva1, fmuidmva2; // MVA muon id
  bool fveto;
  double fmva1, fmva2; // MVA
  // test
  float fmatch1dr1, fmatch2dr1,fmatch1dr2, fmatch2dr2,fmatch1dr3, fmatch2dr3,
    fmatch1dr4, fmatch2dr4, fmatch1dr5, fmatch2dr5;
  bool fmatchTrigs, fb2, fb3,fb4,fb5,fb6,fb7,fb8,fb9;
  bool ftmp1,ftmp2,ftmp3,ftmp4,ftmp5,ftmp6,ftmp7,ftmp8,ftmp9; 
  int fitmp1, fitmp2, fitmp3, fitmp4, fitmp5;

 int     fhltType; // to hold the HLT information d.k
  // mc
  float fmcmds, fmcmdz;
  float fmcpt,fmcptdz,fmcptpis,fmcptpi,fmcptk;

  map<unsigned int, unsigned int, less<unsigned int> > hltObjMap;
  vector<int> hltMatchedMuons;
  vector<int> hltInfo;

};

#endif
