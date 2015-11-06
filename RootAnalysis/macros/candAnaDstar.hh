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
  bool        anaMC();
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
  int         doMuonTriggerMatching(TAnaCand *pC, int mode=0); // match all offline muons 
  double      doTriggerMatchingTest(TAnaTrack *pt, int trig = 0); // match a single track to HLT

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
  float fmatch1dr1, fmatch2dr1,fmatch1dr2, fmatch2dr2,fmatch1dr3, fmatch2dr3,fmatch1dr4, fmatch2dr4;
  bool fb1, fb2, fb3;

  vector<int> hltMatchedMuons;

};

#endif
