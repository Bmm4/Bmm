#include "t1Reader.hh"

#include "TRandom.h"

using namespace std;

#include "t1Reader.icc"

// ----------------------------------------------------------------------
void t1Reader::startAnalysis() {
  cout << "t1Reader: startAnalysis: ..." << endl;
  fpJSON = new JSON(JSONFILE.c_str());
  fpLumi = new Lumi(LUMIFILE.c_str());
}

// ----------------------------------------------------------------------
void t1Reader::endAnalysis() {
  cout << "t1Reader: endAnalysis: ..." << endl;
}


// ----------------------------------------------------------------------
bool t1Reader::goodRun() {
  return true;
}

// ----------------------------------------------------------------------
// this function is normally not called as it is overridden in recoilReader
void t1Reader::eventProcessing() {
  TAnaTrack *pTrack;
  TAnaCand *pCand;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "Found " << fpEvt->nGenCands() << " gen cands in event" << endl;
  cout << "Found " << fpEvt->nSigTracks() << " sig tracks in event" << endl;
  cout << "Found " << fpEvt->nRecTracks() << " rec tracks in event" << endl;
  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nRecTracks());

  cout << "------------------------------" << endl;
  ((TH1D*)fpHistFile->Get("h20"))->Fill(fpEvt->nRecTracks());
  for (int it = 0; it < fpEvt->nRecTracks(); ++it) {
    pTrack = fpEvt->getRecTrack(it);
    ((TH1D*)fpHistFile->Get("h10"))->Fill(pTrack->fPlab.Perp());
    cout << "R: "; pTrack->dump();
  }

  cout << "------------------------------" << endl;
  for (int it = 0; it < fpEvt->nCands(); ++it) {
    pCand = fpEvt->getCand(it);
    cout << "C: " << pCand->fType << " "; pCand->dump();
    ((TH1D*)fpHistFile->Get("h100"))->Fill(pCand->fMass);
  }

  fpHistFile->cd();
  fillHist();
  fTree->Fill();

}


// ----------------------------------------------------------------------
void t1Reader::fillHist() {


}

// ----------------------------------------------------------------------
void t1Reader::bookHist() {
  cout << "==> t1Reader: bookHist> " << endl;

  new TH1D("h1", "nTrk", 40, 0., 40.);
  new TH1D("h10", "pT", 40, 0., 20.);
  new TH1D("h20", "ntrk", 20, 0, 20.);

  new TH1D("h100", "m", 40, 2.8, 3.4);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun ,"run/I");

}


// ----------------------------------------------------------------------
void t1Reader::initVariables() {
  cout << "t1Reader: initVariables: ..." << endl;

  fRun = -1;
  BLIND = 0;
  JSONFILE = "";
  LUMIFILE = "";
}
