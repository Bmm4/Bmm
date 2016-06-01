#ifndef CANDANARECOIL_H
#define CANDANARECOIL_H

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

class candAnaRecoil : public candAna {

public:
  candAnaRecoil(bmmReader *pReader, std::string name, std::string cutsFile);
  virtual ~candAnaRecoil();

  void        evtAnalysis(TAna01Event *);
  void        candAnalysis();

  void        processType();
  void        genMatch();
  void        recoMatch();
  void        candMatch();

  void        bookHist();
  void        readCuts(string filename, int dump);
  void        moreReducedTree(TTree *);
  void        efficiencyCalculation();

  int         BRECOTYPE, BRECOTRUTH;

  double      fBrecoMass, fBrecoPt, fBrecoEta, fBrecoPhi, fBrecoPvZ;
  int         fBrecoPvIdx;
  double      fMu1BrecoPt, fMu1BrecoEta, fMu1BrecoPhi, fMu2BrecoPt, fMu2BrecoEta, fMu2BrecoPhi;

  TAnaCand    *fpBreco;
  int         fBrecoIdx;

  int         fGenBrecoBTmi;
  int         fGenBrecoM1Tmi, fGenBrecoM2Tmi, fNGenBrecoPhotons, fGenBrecoK1Tmi;
  int         fMu1GenBrecoID, fMu2GenBrecoID, fK1GenBrecoID;
  double      fGenBrecoLifeTime;

  int         fRecBrecoM1Tmi, fRecBrecoM2Tmi, fRecBrecoK1Tmi;
  int         fCandBrecoTmi;

  TVector3    fGenPV, fGenBrecoPV;

};

#endif
