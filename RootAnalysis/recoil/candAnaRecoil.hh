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

#define NTRKMAX 1000

class candAnaRecoil : public candAna {

public:
  candAnaRecoil(recoilReader *pReader, std::string name, std::string cutsFile);
  virtual ~candAnaRecoil();

  void        candAnalysis();
  void        brecoAnalysis();

  void        bookHist();
  void        moreReducedTree(TTree *);
  void        readCuts(string filename, int dump);
  void        dump();
  void        dumpBkmt();
  void        printGenBDecays();
  bool        decayModeValidation(TGenCand *pCand, int mode);

  std::map<int, std::vector<int> > fDecayModes;

  int         BRECOTYPE, BRECOTRUTH;

  double      fBrecoMass, fBrecoPt, fBrecoEta, fBrecoPhi, fBrecoPvZ;
  int         fBrecoPvIdx;
  double      fMu1BrecoPt, fMu1BrecoEta, fMu1BrecoPhi, fMu2BrecoPt, fMu2BrecoEta, fMu2BrecoPhi;
  double      fJpsiFlsxy;
  int         fNTrk;
  double      fDoca[NTRKMAX], fProb1[NTRKMAX], fProb2[NTRKMAX], fChi2[NTRKMAX], fDof[NTRKMAX], fTrkMass[NTRKMAX];
  bool        fCorrect[NTRKMAX];

  TAnaCand    *fpBreco;
  std::vector<int> fGenIndices;
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
