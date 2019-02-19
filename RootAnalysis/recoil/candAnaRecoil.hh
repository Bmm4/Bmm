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

  void        genAnalysis();
  void        candAnalysis();
  void        brecoAnalysis();

  void        bookHist();
  void        moreReducedTree(TTree *);
  void        readCuts(string filename, int dump);

  std::pair<TVector3, TVector3> parallelAndPerp(TVector3 direction, TVector3 momVis);
  std::pair<TVector3, TVector3> parallelAndPerp2(TVector3 direction, TVector3 momVis);
  std::pair<double, double>  nuRecoMom0(double compVisPar, double compVisPerp, double eVis, double mVis, double mTot = 1.777);
  double      minMassPair(std::vector<TLorentzVector> a);
  double      maxMassPair(std::vector<TLorentzVector> a);

  void        dump();
  TH1D*       getHist(std::string);
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

  // -------------------
  // -- recoil = signal!
  // -------------------

  // -- pointers to HepMC cands and vertices
  TGenCand    *fpGenB,  *fpGenMu,  *fpGenKa, *fpGenTau, *fpGenHad1, *fpGenHad2, *fpGenHad3, *fpGenNu, *fpGenGa1, *fpGenGa2;
  TVector3    fGenVtxBProd, fGenVtxBDecay, fGenVtxTauDecay, fGenDirectionTau;
  // -- calculated/derived quantities
  TLorentzVector f4GenTauAtDecay, f4GenTauHad, f4GenNur0, f4GenNur0Pos, f4GenNur0Neg, f4GenBr0, f4GenBr0Pos, f4GenBr0Neg;
  // -- reco'ed/derived quantities
  TVector3    fVtxBProd, fVtxBDecay, fVtxTauDecay, fDirectionTau;
  TLorentzVector f4Muon, f4Kaon, f4TauHad, f4Nur0Pos, f4Nur0Neg, f4Br0Pos, f4Br0Neg;
  int            fNgamma;

  // -- gen quantities
  double         fGenBPt,
    fGenBr0Mass, fGenBr0PosMass, fGenBr0NegMass,
    fGenTauPt, fGenTauDecTime, fGenTauPara, fGenTauPerp,
    fGenTaur0Pt, fGenTaur0PosPt, fGenTaur0NegPt,
    fGenTauHadPt, fGenTauHadMass, fGenTauHadAngle, fGenTauHadPerp, fGenTauHadPara,
    fGen2HadMinMass, fGen2HadMaxMass,
    fGen2HadSSMass, fGen2HadOSminMass, fGen2HadOSmaxMass,
    fGenHad1Pt, fGenHad2Pt, fGenHad3Pt,
    fGenMuPt, fGenKaPt, fGenMuKaMass,
    fGenNur0PosPara, fGenNur0NegPara,
    fGenNuPara, fGenNuPerp
    ;


  // -- reco quantities
  double fPvX, fPvY, fPvZ,
    fBPt,
    fBrMass, fBr0PosMass, fBr0NegMass,
    fTauPt, fTauDecTime,
    fTaur0Pt, fTaur0PosPt, fTaur0NegPt,
    fTauHadPt, fTauHadMass, fTauHadAngle, fTauHadPerp, fTauHadPara, fTauFl, fTauFls,
    fBTauAngle,
    f2HadMinMass,f2HadMaxMass,
    f2HadSSMass, f2HadOSminMass, f2HadOSmaxMass,
    fHad1Pt, fHad2Pt, fHad3Pt, fHad1Eta, fHad2Eta, fHad3Eta,
    fMuPt, fKaPt, fMuEta, fKaEta, fMuKaMass,
    fNur0PosPara, fNur0NegPara
    ;



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
