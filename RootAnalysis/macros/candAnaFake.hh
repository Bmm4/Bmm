#ifndef CANDANAFAKE_H
#define CANDANAFAKE_H

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

class candAnaFake : public candAna {

public:
  candAnaFake(bmmReader *pReader, std::string name, std::string cutsFile);
  ~candAnaFake();

  void        candAnalysis();
  // void        hhAnalysis();
  void        efficiencyCalculation();

  void        processType();
  void        genMatch();
  void        genMatchOld();
  void        recoMatch();
  void        candMatch();

  void        bookHist();
  void        readCuts(string filename, int dump);

  void        evtAnalysis(TAna01Event *evt);
  bool        anaMC(TAna01Event *evt);


  int         truthMatch(TAnaCand *pC, int verbose = 0);
  void        dumpHFTruthCand(TAnaCand *pC);
  // void        dumpHFHhCand(TAnaCand *pC);


private:

  static const int NTRKMAX = 10;
  TTree *fFakeTree;
  double  fCandMKS, fCandMconv, fCandMLambda, fCandMEta;
  int    fFakeNtrk, fFakeId[NTRKMAX], fFakeQ[NTRKMAX], fFakeGm[NTRKMAX];
  float  fFakePt[NTRKMAX], fFakeEta[NTRKMAX], fFakePhi[NTRKMAX], fFakeBdt[NTRKMAX], fFakeDistTrig[NTRKMAX], fFakeDistMuon[NTRKMAX];
  bool   fFakeHP[NTRKMAX], fFakeTis[NTRKMAX];

  float  fFakeTip[NTRKMAX], fFakeLip[NTRKMAX]
    , fFakeInnerChi2[NTRKMAX]
    , fFakeOuterChi2[NTRKMAX]
    , fFakeChi2LocalPosition[NTRKMAX]
    , fFakeChi2LocalMomentum[NTRKMAX]
    , fFakeStaTrkMult[NTRKMAX]
    , fFakeTmTrkMult[NTRKMAX]
    , fFakeDeltaR[NTRKMAX]
    , fFakeDxyRef[NTRKMAX]
    , fFakeDzRef[NTRKMAX]
    , fFakeItrkValidFraction[NTRKMAX]
    , fFakeSegmentComp[NTRKMAX]
    , fFakeGtrkNormChi2[NTRKMAX]
    , fFakeDz[NTRKMAX]
    , fFakeGtrkProb[NTRKMAX]
    , fFakeMuonChi2[NTRKMAX]
    , fFakeGlbKinkFinder[NTRKMAX]
    , fFakeStaRelChi2[NTRKMAX]
    , fFakeTrkRelChi2[NTRKMAX]
    , fFakeGlbDeltaEtaPhi[NTRKMAX]
    , fFakeTimeInOut[NTRKMAX]
    , fFakeTimeInOutE[NTRKMAX]
    , fFakeTimeInOutS[NTRKMAX]
    ;

  int fFakeNvalidMuonHits[NTRKMAX]
    , fFakeNmatchedStations[NTRKMAX]
    , fFakeLayersWithHits[NTRKMAX]
    , fFakeNumberOfValidTrkHits[NTRKMAX]
    , fFakeNumberOfLostTrkHits[NTRKMAX]
    , fFakeNumberOfValidPixHits[NTRKMAX]
    , fFakeRPChits[NTRKMAX]
    , fFakeRPChits1[NTRKMAX]
    , fFakeRPChits2[NTRKMAX]
    , fFakeRPChits3[NTRKMAX]
    , fFakeRPChits4[NTRKMAX]
    , fFakeMuonDetectorHitsCombination[NTRKMAX]

    ;


};

#endif
