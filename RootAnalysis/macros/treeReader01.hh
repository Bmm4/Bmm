#ifndef TREEREADER01_H
#define TREEREADER01_H

#include <iostream>
#include <string>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TTimeStamp.h>

#include "rootio/TAna01Event.hh"
#include "rootio/TGenCand.hh"
#include "rootio/TAnaCand.hh"
#include "rootio/TAnaTrack.hh"
#include "rootio/TAnaVertex.hh"

#include "common/JSON.hh"
#include "common/Lumi.hh"
#include "common/PdTrigger.hh"

#define DR      57.29577951

class treeReader01 {
public:
  treeReader01(TChain *tree, TString evtClassName);
  virtual      ~treeReader01();
  virtual void       init(TString evtClassName);

  virtual void       openHistFile(TString filename);
  virtual void       closeHistFile();
  virtual void       bookHist();
  virtual void       readCuts(TString filename, int dump = 1);

  virtual void       startAnalysis();
  virtual void       endAnalysis();
  virtual int        loop(int nevents = 1, int start = -1);
  virtual TFile*     getFile() {return fpChain->GetCurrentFile();}
  virtual void       eventProcessing();
  virtual void       initVariables();
  virtual void       fillHist();
  virtual bool       goodRun();
  virtual void       setVerbosity(int f) {std::cout << Form("setVerbosity(%d)", f) << std::endl;  fVerbose = f;}
  virtual void       setYear(int f) {std::cout << Form("setYear(%d)", f) << std::endl;  fYear = f;}
  virtual void       setEra(std::string s) {std::cout << Form("setEra(%s)", s.c_str()) << std::endl;  fEra = s;}
  virtual void       setMC(int f) {std::cout << Form("setMC(%d)", f) << std::endl; fIsMC = f;}
  virtual void       setCheckCandTypes(bool f) {std::cout << "checking cand types: " << f << std::endl; fCheckCandTypes = f;}
  virtual void       runBlind() {std::cout << "running blinded" << std::endl; BLIND = 1;}
  virtual int        numberOfBPixLayers(TAnaTrack *t);
  virtual int        numberOfPixLayers(TAnaTrack *t);
  virtual int        numberOfPixelHits(TAnaTrack *pTrack);
  virtual int        numberOfBPixLayer1Hits(TAnaTrack *t);
  virtual int        numberOfTrackerLayers(TAnaTrack *t);
  virtual void       setJSONFile(const char *name) {JSONFILE = name;}
  virtual void       setLumiFile(const char *name) {LUMIFILE = name;}
  virtual void       ignoreJSON()  {std::cout << "ignoreJSON() called" << std::endl; fIgnoreJson = true;}
  virtual PdTrigger* pdTrigger() {return fpPdTrigger;}

  int fVerbose;
  int fYear;
  std::string fEra;

protected:

  TChain      *fpChain;        // pointer to the analyzed TTree or TChain
  TFile       *fpHistFile;     // for output histograms and reduced trees
  TString      fChainFileName; // the name of the chain file
  TString      fCutFile;       // contains file with the cut definitions

  TAna01Event *fpEvt;

  // -- Pre-filled variables
  int          fNentries;      // number of events in chain; filled in treeReader01::treeReader01()
  int          fEvent;         // current sequential event number in chain; filled in treeReader01::loop()
  long int     fEvt;           // current event number; filled in treeReader01::loop()
  long int     fRun;           // current run number; filled in treeReader01::loop()
  int          fLS;            // current lumi section; filled in treeReader01::loop()


  // -- Histogram pointers
  TTree       *fTree;

  JSON *fpJSON;
  PdTrigger *fpPdTrigger;
  Lumi *fpLumi;

  bool fCheckCandTypes;

  // -- Cut values
  double
      PTLO
    , PTHI
    , ETALO
    , ETAHI
    ;
  int TYPE;
  int BLIND;
  std::string JSONFILE, DSNAME, LUMIFILE;
  bool fIgnoreJson;

  int fIsMC;
};

// ----------------------------------------------------------------------
inline void mk4Vector(TLorentzVector &p4, const Double_t p, const Double_t t, const Double_t f, const Double_t m) {
  p4.SetXYZM(p*TMath::Sin(t)*TMath::Cos(f), p*TMath::Sin(t)*TMath::Sin(f), p*TMath::Cos(t), m);
}

#endif
