#include "genAnalysis.hh"
#include "common/HFMasses.hh"
#include "common/util.hh"

#include "TRandom.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

#include <algorithm>
#include <cmath>
#include <map>

using namespace std;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: bin/runBmm -r genAnalysis -f test.root
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// ----------------------------------------------------------------------
genAnalysis::genAnalysis(TChain *tree, TString evtClassName): t1Reader(tree, evtClassName) {
  cout << "==> genAnalysis: constructor..." << endl;
  // -- signal modes
  vector<int> v80; v80.push_back(531); v80.push_back(13); v80.push_back(13);
  fDecayModes.insert(make_pair(80, v80));
  vector<int> v200; v200.push_back(521); v200.push_back(321); v200.push_back(15); v200.push_back(211); v200.push_back(211); v200.push_back(211); v200.push_back(13);
  fDecayModes.insert(make_pair(200, v200));
  vector<int> v300; v300.push_back(521); v200.push_back(321); v200.push_back(15); v200.push_back(211); v200.push_back(211); v200.push_back(211); v200.push_back(13);
  fDecayModes.insert(make_pair(200, v200));

  // -- BRECO modes
  vector<int> v68; v68.push_back(521); v68.push_back(443); v68.push_back(13); v68.push_back(13); v68.push_back(321);
  fDecayModes.insert(make_pair(68, v68));
  vector<int> v69; v69.push_back(511); v69.push_back(443); v69.push_back(13); v69.push_back(13); v69.push_back(321); v69.push_back(321); v69.push_back(211);
  fDecayModes.insert(make_pair(69, v69));
}

// ----------------------------------------------------------------------
genAnalysis::~genAnalysis() {
  cout << "==> genAnalysis: destructor..." << endl;
}

// ----------------------------------------------------------------------
void genAnalysis::startAnalysis() {
  cout << "==> genAnalysis: Starting analysis..." << endl;
}


// ----------------------------------------------------------------------
void genAnalysis::eventProcessing() {
  if (0) fpEvt->dumpGenBlock();

  if (1) printBdecays();
  if (0) printB2JpsiXdecays();
  if (0) recoilValidation();

  fTree->Fill();
}


// ----------------------------------------------------------------------
void genAnalysis::endAnalysis() {

}


// ----------------------------------------------------------------------
void genAnalysis::recoilValidation() {

  TGenCand *pCand, *pD;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "New event: " << fChainEvent << " gen block with " << fpEvt->nGenCands() << " gen cands" << endl;
  int nreco(0), ncand(0);
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if ((CANDTRUTH == TMath::Abs(pCand->fID)) && ((pCand->fDau2 - pCand->fDau1) > 0)) {
      if (!decayModeValidation(pCand, 80)) {
	//	cout << "Failed decayModeValidation(pCand, 80)" << endl;
	continue;
      }
      ++ncand;
      fCanId = pCand->fID;
      fCanMode = 68;
      fCanPt = pCand->fP.Perp();
      fCanEta = pCand->fP.Eta();
      fCanPhi = pCand->fP.Phi();

      if (0) {
	cout << "signal decay: " << endl;
	pCand->dump(1);
	for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	  pD = fpEvt->getGenCand(iD);
	  pD->dump(1);
	}
      }
    }
    if (BRECOTRUTH == TMath::Abs(pCand->fID)) {
      if (!decayModeValidation(pCand, 68)) {
	//	cout << "Failed decayModeValidation(pCand, 68)" << endl;
	continue;
      }
      ++nreco;
      fRecId = pCand->fID;
      fRecMode = 68;
      fRecPt = pCand->fP.Perp();
      fRecEta = pCand->fP.Eta();
      fRecPhi = pCand->fP.Phi();
      if (0) {
	cout << "RECO decay: " << endl;
	pCand->dump(1);
	for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	  pD = fpEvt->getGenCand(iD);
	  pD->dump(1);
	}
      }
    }
  }

  if (0 == nreco || 0 == ncand) {
    cout << "XXXXXXXXXXXXXXXXXX no recoil gen event!!!!!!" << endl;
  }
}

// ----------------------------------------------------------------------
bool genAnalysis::decayModeValidation(TGenCand *pCand, int mode) {
  vector<int> tDaughters = fDecayModes[mode];
  if (tDaughters[0] != TMath::Abs(pCand->fID)) {
    return false;
  }
  TGenCand *pD(0);
  vector<int> cDaughters;
  cDaughters.push_back(TMath::Abs(pCand->fID));
  for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
    pD = fpEvt->getGenCand(iD);
    if (22 == pD->fID) continue;
    cDaughters.push_back(TMath::Abs(pD->fID));
  }

  if (cDaughters.size() != tDaughters.size()) {
    return false;
  }
  for (unsigned int i = 0; i < tDaughters.size(); ++i) {
    if (tDaughters[i] != cDaughters[i]) return false;
  }

  return true;
}

// ----------------------------------------------------------------------
void genAnalysis::printBdecays() {

  TGenCand *pCand, *pD, *pDD;
  cout << "--- Event " << fChainEvent << " ------------------------------------------------------------------" << endl;
  cout << "gen block with " << fpEvt->nGenCands() << " gen cands" << endl;
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if ((TMath::Abs(pCand->fID) == 531)
	|| (TMath::Abs(pCand->fID) == 511)
	|| (TMath::Abs(pCand->fID) == 521)
	|| (TMath::Abs(pCand->fID) == 541)
	|| (TMath::Abs(pCand->fID) == 5122)
	|| (TMath::Abs(pCand->fID) == 5112)
	|| (TMath::Abs(pCand->fID) == 5212)
	|| (TMath::Abs(pCand->fID) == 5222)
	) {
      cout << "B meson" << endl;
      pCand->dump(1);
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = fpEvt->getGenCand(iD);
	pD->dump(1);
	// -- tau decays are not part of the B candidate daughters!? ONLY if the tau is the last particle in the decay chain!
	// if (15 == TMath::Abs(pD->fID)) {
	//   for (int idd = pD->fDau1; idd <= pD->fDau2; ++idd) {
	//     pDD = fpEvt->getGenCand(idd);
	//     pDD->dump(1);
	//   }
	// }
      }
    }
  }

  return;
}

// ----------------------------------------------------------------------
void genAnalysis::printB2JpsiXdecays() {
  TGenCand *pCand, *pD, *pDD;
  cout << "--- Event " << fChainEvent << " ------------------------------------------------------------------" << endl;
  cout << "gen block with " << fpEvt->nGenCands() << " gen cands" << endl;
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if ((TMath::Abs(pCand->fID) == 531)
	|| (TMath::Abs(pCand->fID) == 511)
	|| (TMath::Abs(pCand->fID) == 521)
	|| (TMath::Abs(pCand->fID) == 541)
	|| (TMath::Abs(pCand->fID) == 5122)
	|| (TMath::Abs(pCand->fID) == 5112)
	|| (TMath::Abs(pCand->fID) == 5212)
	|| (TMath::Abs(pCand->fID) == 5222)
	) {
      int OK(0);
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = fpEvt->getGenCand(iD);
	if (TMath::Abs(pD->fID) == 443) {
	  OK = 1;
	}
      }
      if (0 == OK) {
	cout << "B meson" << endl;
	pCand->dump(1);
	continue;
      }
      cout << "B -> J/psiX meson" << endl;
      pCand->dump(1);
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = fpEvt->getGenCand(iD);
	pD->dump(1);
	// -- tau decays are not part of the B candidate daughters!? ONLY if the tau is the last particle in the decay chain!
	// if (15 == TMath::Abs(pD->fID)) {
	//   for (int idd = pD->fDau1; idd <= pD->fDau2; ++idd) {
	//     pDD = fpEvt->getGenCand(idd);
	//     pDD->dump(1);
	//   }
	// }
      }
    }
  }

  return;
}

// ----------------------------------------------------------------------
void genAnalysis::initVariables() {
}


// ----------------------------------------------------------------------
void genAnalysis::fillHist() {
}

// ----------------------------------------------------------------------
void genAnalysis::bookHist() {
  cout << "==> genAnalysis: bookHist " << endl;
  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",     &fRun,    "run/I");
  fTree->Branch("recid",   &fRecId,  "recid/I");
  fTree->Branch("recpt",   &fRecPt,  "recpt/D");
  fTree->Branch("receta",  &fRecEta, "receta/D");
  fTree->Branch("recphi",  &fRecPhi, "recphi/D");
  fTree->Branch("canid",   &fCanId,  "canid/I");
  fTree->Branch("canpt",   &fCanPt,  "canpt/D");
  fTree->Branch("caneta",  &fCanEta, "caneta/D");
  fTree->Branch("canphi",  &fCanPhi, "canphi/D");

}

// ----------------------------------------------------------------------
void genAnalysis::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> genAnalysis: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> genAnalysis: Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  int ibin(-1);
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "CANDTYPE")) {
      CANDTYPE = int(CutValue);
      if (dump) cout << "CANDTYPE:       " << CANDTYPE << endl;
      ibin = 2;
      hcuts->SetBinContent(ibin, CANDTYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: CANDTYPE       :: %i", CutName, CANDTYPE));
      ok = 1;
    }

    if (!strcmp(CutName, "CANDTRUTH")) {
      CANDTRUTH = int(CutValue);
      if (dump) cout << "CANDTRUTH:      " << CANDTRUTH << endl;
      ibin = 3;
      hcuts->SetBinContent(ibin, CANDTRUTH);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: CANDTRUTH      :: %i", CutName, CANDTRUTH));
      ok = 1;
    }

    if (!strcmp(CutName, "BRECOTYPE")) {
      BRECOTYPE = int(CutValue);
      if (dump) cout << "BRECOTYPE:      " << BRECOTYPE << endl;
      ibin = 4;
      hcuts->SetBinContent(ibin, BRECOTYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BRECOTYPE      :: %i", CutName, BRECOTYPE));
      ok = 1;
    }

    if (!strcmp(CutName, "BRECOTRUTH")) {
      BRECOTRUTH = int(CutValue);
      if (dump) cout << "BRECOTRUTH:     " << BRECOTRUTH << endl;
      ibin = 5;
      hcuts->SetBinContent(ibin, BRECOTRUTH);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BRECOTRUTH     :: %i", CutName, BRECOTRUTH));
      ok = 1;
    }


    if (!ok) cout << "==> genAnalysis: ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}
