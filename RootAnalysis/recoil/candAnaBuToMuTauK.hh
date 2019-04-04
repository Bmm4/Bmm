#ifndef CANDANABUTOMUTAUK_H
#define CANDANABUTOMUTAUK_H

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

class candAnaBuToMuTauK : public candAna {

public:
  candAnaBuToMuTauK(recoilReader *pReader, std::string name, std::string cutsFile);
  virtual ~candAnaBuToMuTauK();

  void        genAnalysis();
  void        candAnalysis();
  void        candEvaluation();

  void        genMatch();
  void        recoMatch();
  void        candMatch();


  void        resetAllData();
  void        bookHist();
  void        moreReducedTree(TTree *);
  void        readCuts(string filename, int dump);

  void        dump();
  void        dumpBkmt();
  void        printGenBDecays();

  // -------------------
  // -- recoil = signal!
  // -------------------

  // -- pointers to SIGNAL HepMC cands and vertices
  TGenCand    *fpGenMu,  *fpGenKa, *fpGenTau, *fpGenHad1, *fpGenHad2, *fpGenHad3, *fpGenNu, *fpGenGa1, *fpGenGa2;
  TVector3    fGenVtxBProd, fGenVtxBDecay, fGenVtxTauDecay, fGenDirectionTau;
  // -- calculated/derived quantities
  TLorentzVector f4GenTauAtDecay, f4GenTauHad, f4GenNur0, f4GenNur0Pos, f4GenNur0Neg, f4GenBr0, f4GenBr0Pos, f4GenBr0Neg;
  // -- reco'ed/derived quantities
  TVector3    fVtxTauDecay, fDirectionB, fDirectionTau;
  TLorentzVector f4Muon, f4Kaon, f4MuKa, f4Had, f4Nur0Pos, f4Nur0Neg, f4Br0Pos, f4Br0Neg;

  // -- pointers to SIGNAL simpleTracks
  TSimpleTrack *fpTmMu,  *fpTmKa, *fpTmHad1, *fpTmHad2, *fpTmHad3;
  // -- vector of all signal tracks: mu, ka, had1, had2, had3
  std::vector<TAnaTrack*> fSignalTracks;

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
  double fBPt,
    fBrMass, fBr0PosMass, fBr0NegMass,
    fMuKaDoca3D, fMuKaLip, fMuKaDocaMax, fMuKaFl, fMuKaFls,
    fHadPt, fHadDecTime, fHadFl, fHadFls, fHadMass, fHadDoca3D, fHadLip, fHadDocaMax,
    fTaur0Pt, fTaur0PosPt, fTaur0NegPt,
    fTauHadAngle, fTauHadPerp, fTauHadPara,
    fBTauAngle, fBMuKaAngle,
    f2HadMinMass,f2HadMaxMass,
    f2HadSSMass, f2HadOSminMass, f2HadOSmaxMass,
    fHad1Pt, fHad2Pt, fHad3Pt, fHad1Eta, fHad2Eta, fHad3Eta,
    fMuPt, fKaPt, fMuEta, fKaEta, fMuKaMass,
    fNur0PosPara, fNur0NegPara
    ;


  // -- tracks of signal cand
  int         fNTrk;
  double      fDoca[NTRKMAX], fPt[NTRKMAX];
  bool        fCorrect[NTRKMAX], fInRecoil[NTRKMAX];



};

#endif
