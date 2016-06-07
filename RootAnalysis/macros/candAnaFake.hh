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


};

#endif
