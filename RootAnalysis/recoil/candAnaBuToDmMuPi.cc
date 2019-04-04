#include "candAnaBuToDmMuPi.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"
#include "common/ana.hh"

#include "common/AnalysisDistribution.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBuToDmMuPi::candAnaBuToDmMuPi(recoilReader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  fAnaCuts.setAcName("candAnaBuToDmMuPi");
  BLIND = 0;
  cout << "==> candAnaBuToDmMuPi: name = " << name << ", reading cutsfile " << cutsFile << endl;
  readCuts(cutsFile, 1);
}


// ----------------------------------------------------------------------
candAnaBuToDmMuPi::~candAnaBuToDmMuPi() {
  cout << "==> candAnaBuToDmMuPi: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBuToDmMuPi::candAnalysis() {

  candAna::candAnalysis();

  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(10);
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(3);

}

// ----------------------------------------------------------------------
void candAnaBuToDmMuPi::candEvaluation() {

}


// ----------------------------------------------------------------------
// -- version for PYTHIA8 and new-style EvtGen daughter insertion
void candAnaBuToDmMuPi::genMatch() {

}

// ----------------------------------------------------------------------
void candAnaBuToDmMuPi::recoMatch() {

}


// ----------------------------------------------------------------------
void candAnaBuToDmMuPi::candMatch() {

}


// ----------------------------------------------------------------------
void candAnaBuToDmMuPi::bookHist() {
  cout << "==>candAnaBuToDmMuPi: bookHist" << endl;
  candAna::bookHist();

  moreReducedTree(fTree);

  // -- Additional effTree variables

}



// ----------------------------------------------------------------------
void candAnaBuToDmMuPi::moreReducedTree(TTree *t) {
  // -- Additional reduced tree variables
}


// ----------------------------------------------------------------------
void candAnaBuToDmMuPi::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump);

  fCutFile = filename;

  if (dump) cout << "==> candAnaBuToDmMuPi: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines;
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  string cstring = "B cand";

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str());

    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);
  }
}
