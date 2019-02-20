#include "recoilReader.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"
#include "common/util.hh"

#include "candAna.hh"
#include "candAnaBuToMuTauK.hh"

using namespace std;

// ----------------------------------------------------------------------
recoilReader::recoilReader(TChain *tree, TString evtClassName): t1Reader(tree, evtClassName) {
  cout << "==> recoilReader: constructor..." << endl;
  fVerbose = 0;
}


// ----------------------------------------------------------------------
recoilReader::~recoilReader() {
  cout << "==> recoilReader: destructor..." << endl;
}


// ----------------------------------------------------------------------
void recoilReader::startAnalysis() {
  cout << "==> recoilReader: fVerbose = " << fVerbose << endl;
  if (JSONFILE == "") {
    cout << "No JSONFILE provided, ignoring JSON" << endl;
    fpJSON = 0;
    fIgnoreJson = true;
  } else {
    fpJSON = new JSON(JSONFILE.c_str(), fVerbose);
  }
  if (LUMIFILE == "") {
    cout << "No LUMIFILE provided, deriving LUMIFILE from " << JSONFILE << endl;
    LUMIFILE = JSONFILE;
    replaceAll(LUMIFILE, "txt", "lumi");
    if (LUMIFILE == "") {
      cout << "Cannot derive LUMIFILE, no lumi information" << endl;
      fpLumi = 0;
    } else {
      fpLumi = new Lumi(LUMIFILE.c_str(), fVerbose);
    }
  }
}

// ----------------------------------------------------------------------
void recoilReader::endAnalysis() {
  cout << "==> recoilReader: endAnalysis() destroying cand analysis modules" << endl;
  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    lCandAnalysis[i]->endAnalysis();
  }

}


// ----------------------------------------------------------------------
void recoilReader::eventProcessing() {
  ((TH1D*)fpHistFile->Get("monEvents"))->Fill(0);
  static bool json = false;
  static int oldRun(-1);
  static int oldLS(-1);
  static double rlumi(-1.);
  if (fIsMC) {
    json = true;
  } else if (fIgnoreJson) {
    json = true;
  } else {
    if ((fLS != oldLS) || (fRun != oldRun)) {
      oldLS = fLS;
      // do not reset oldRun, this is done below!
      json = fpJSON->good(fRun, fLS);
    }
  }

  if (fRun != oldRun) {
    oldRun = fRun;
    if (!fIsMC && json && !fIgnoreJson) {
      if (0 != fpLumi->contains(fRun)) {
	rlumi = fpLumi->lumi(fRun);
      } else {
	if (fVerbose > 100) {
	  cout << "Run " << fRun << " has no lumi information" << endl;
	}
      }
    }
  }

  if (fCheckCandTypes) {
    fCandTypes.clear();
    for (int i = 0; i < fpEvt->nCands(); ++i) {
      fCandTypes.insert(fpEvt->getCand(i)->fType);
    }
  }

  TH1D *hbs = (TH1D*)fpHistFile->Get("bsz");
  hbs->Fill(fpEvt->fBeamSpot.fPoint.Z());
  TH1D *h1 = (TH1D*)fpHistFile->Get("npv");
  h1->Fill(fpEvt->nPV());
  TH1D *h2 = (TH1D*)fpHistFile->Get("pvz");
  TH1D *h0 = (TH1D*)fpHistFile->Get("pv0z");
  for (int i = 0; i < fpEvt->nPV(); ++i) {
    double z = fpEvt->getPV(i)->fPoint.Z();
    if (0 == i) h0->Fill(fpEvt->getPV(i)->fPoint.Z());
    h2->Fill(z);

  }

  unsigned int ilm = lCandAnalysis.size();
  if (ilm < 1) return;

  // -- call candidate analyses
  // cout << "recoilReader: " << fChainEvent << endl;
  // cout << "lCandAnalysis.size() = " << lCandAnalysis.size() << endl;
  for (unsigned int i = 0; i < ilm; ++i) {
    if (fCheckCandTypes) {
      if (fCandTypes.find(lCandAnalysis[i]->CANDTYPE) == fCandTypes.end()) {
	continue;
      }
    }

    lCandAnalysis[i]->fIsMC        = fIsMC;
    lCandAnalysis[i]->fJSON        = json;
    lCandAnalysis[i]->fRun         = fRun;
    lCandAnalysis[i]->fEvt         = fEvt;
    lCandAnalysis[i]->fLS          = fLS;
    lCandAnalysis[i]->fLumi        = rlumi;
    lCandAnalysis[i]->fChainEvent  = fChainEvent;

    lCandAnalysis[i]->evtAnalysis(fpEvt);
  }

}


// ----------------------------------------------------------------------
void recoilReader::bookHist() {
  fpHistFile->cd();
  TH1D *h(0);
  h = new TH1D("monEvents", "monEvents", 20, 0., 20.);
  h = new TH1D("ntracks", "ntracks", 100, 0., 1000.);
  h = new TH1D("pt0", "pt(tracks)", 100, 0., 10.);
  h = new TH1D("pt1", "pt(tracks)", 100, 0., 50.);
  h = new TH1D("eta", "eta(tracks)", 60, -3.0, 3.0);
  h = new TH1D("phi", "phi(tracks)", 50, -3.15, 3.15);
  h = new TH1D("npv", "npv", 100, 0, 100.);
  h = new TH1D("pvz", "pvz", 100, -25., 25.);
  h = new TH1D("pv0z","pv0z", 100, -25., 25.);
  h = new TH1D("bsz", "bsz", 100, -25., 25.);

  (void)h; // make compiler warning go away

  for (unsigned int i = 0; i < lCandAnalysis.size(); ++i) {
    lCandAnalysis[i]->bookHist();
  }

}


// ----------------------------------------------------------------------
// Changed setup compared to previous times:
// directly pass all candAna* setups after -C, separated by :
// The first line in the candAna* file should contain CLASSNAME for proper instantiation
// This avoids duplicating all cut files for the reader and candAna*
// ----------------------------------------------------------------------
void recoilReader::readCuts(TString filenames, int dump) {
  if (dump) cout << "==> recoilReader: Reading " << filenames << " for classes setup" << endl;

  vector<string> cutfiles = split(string(filenames), ':');
  ifstream INS;
  for (unsigned int i = 0; i < cutfiles.size(); ++i) {
    cout << "cutfile i = " << i << ": " << cutfiles[i] << endl;
    string sline("nada");
    INS.open(cutfiles[i]);
    while (getline(INS, sline)) {
      if (string::npos != sline.find("CLASS")) {
	replaceAll(sline, "CLASSNAME", "");
	cleanupString(sline);
	break;
      }
    }
    INS.close();

    if (string::npos != sline.find("candAnaBuToMuTauK")) {
      candAna *a = new candAnaBuToMuTauK(this, sline, cutfiles[i]);
      a->BLIND = BLIND;
      lCandAnalysis.push_back(a);
    }
  }

  cout << "Added " << lCandAnalysis.size() << " candidate analysis modules to lCandAnalysis" << endl;
}
