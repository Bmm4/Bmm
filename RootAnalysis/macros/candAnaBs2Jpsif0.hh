#ifndef CANDANABS2JPSIF0_H
#define CANDANABS2JPSIF0_H

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

class candAnaBs2Jpsif0 : public candAna {

public:
  candAnaBs2Jpsif0(bmmReader *pReader, std::string name, std::string cutsFile);
  virtual ~candAnaBs2Jpsif0();

  void        candAnalysis();
  void        efficiencyCalculation();
  void        moreBasicCuts();
  void        moreReducedTree(TTree *);

  void        genMatch();
  void        genMatchOld();
  void        recoMatch();
  void        candMatch();

  void        readCuts(string filename, int dump);
  void        bookHist();
  void        fillCandidateHistograms(int offset);

  int          JPSITYPE;
  double       JPSIMASSLO, JPSIMASSHI;

  double       fJpsiMass, fJpsiPt, fJpsiEta, fJpsiPhi, fJpsiCosA, fJpsiMaxDoca, fJpsiFLSxy, fJpsiVtxProb;

  bool         fGoodJpsiMass;

  double       fPi1Pt, fPi1Eta, fPi1Phi;
  double       fPi2Pt, fPi2Eta, fPi2Phi;
  double       ff0Pt, ff0Eta, ff0Phi;

  double       fPi1PtNrf, fPi1EtaNrf;
  double       fPi2PtNrf, fPi2EtaNrf;

  double       fPi1PtGen, fPi1EtaGen, fPi2PtGen, fPi2EtaGen;
  int          fPi1GenID, fPi2GenID;
  int          fPi1TkQuality, fPi2TkQuality;

  // -- TM
  int                     fGenPi1Tmi, fGenPi2Tmi;
  int                     fRecPi1Tmi, fRecPi2Tmi;

  // -- effTree
  float fETpi1pt, fETpi1eta, fETg3pt, fETg3eta;
  float fETpi2pt, fETpi2eta, fETg4pt, fETg4eta;
  int   fETpi1q,  fETpi2q;
  bool  fETpi1gt, fETpi2gt;

  // -- Additional variables and cuts for Bs -> J/psi f0
  int               F0TYPE;
  double            MPIPILO, MPIPIHI, DELTAR;
  double            ff0DeltaR, fMPiPi;
  bool              fGoodDeltaR, fGoodMPIPI;

  double            fBsJpsiPhiMass;

};

#endif
