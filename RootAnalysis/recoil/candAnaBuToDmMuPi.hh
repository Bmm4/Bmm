#ifndef CANDANABUTODMMUPI_H
#define CANDANABUTODMMUPI_H

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

class candAnaBuToDmMuPi : public candAna {

public:
  candAnaBuToDmMuPi(recoilReader *pReader, std::string name, std::string cutsFile);
  virtual ~candAnaBuToDmMuPi();

  void        candAnalysis();
  void        moreReducedTree(TTree *);

  void        genMatch();
  void        recoMatch();
  void        candMatch();

  void        readCuts(string filename, int dump);
  void        bookHist();
  void        fillCandidateHistograms(int offset);



  // -- effTree

};

#endif
