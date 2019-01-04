#include "candAnaRecoil.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"
#include "common/util.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaRecoil::candAnaRecoil(recoilReader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaRecoil: constructor..." << endl;
  readCuts(cutsFile, 1);

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
candAnaRecoil::~candAnaRecoil() {
  cout << "==> candAnaRecoil: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaRecoil::dumpBkmt() {
  cout << "event " << fEvt << endl;
  fGenIndices.clear();
  TGenCand *pCand, *pD;
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if (521 == TMath::Abs(pCand->fID)) {
      //      if (decayModeValidation(pCand, 200)) {
	for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	  pD = fpEvt->getGenCand(iD);
	  pD->dump();
	  // if (isStableCharged(pD->fID)) {
	  //   if (fVerbose > -1) cout << "gen idx = " << pD->fNumber << " with pT = " << pD->fP.Perp() << " eta = " << pD->fP.Eta() << endl;
	  //   fGenIndices.push_back(pD->fNumber);
	  // }
	}
	//      }
    }
  }

}


// ----------------------------------------------------------------------
void candAnaRecoil::dump() {
  fGenIndices.clear();

  TGenCand *pCand, *pD;
  if (fVerbose > 1)  {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "New event: " << fChainEvent << " gen block with " << fpEvt->nGenCands() << " gen cands" << endl;
  }
  int nreco(0), ncand(0);
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if (521 == TMath::Abs(pCand->fID)) {
      if (decayModeValidation(pCand, 68)) {
	for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	  pD = fpEvt->getGenCand(iD);
	  if (isStableCharged(pD->fID)) {
	    if (fVerbose > 1) cout << "gen idx = " << pD->fNumber << " with pT = " << pD->fP.Perp() << " eta = " << pD->fP.Eta() << endl;
	    fGenIndices.push_back(pD->fNumber);
	  }
	}
      }
    }
  }
  if (fVerbose > 1) {
    for (unsigned int i = 0; i < fGenIndices.size(); ++i) {
      cout << "gen index: " << fGenIndices[i] << endl;
    }
  }

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    if (fVerbose > 2) pCand->dump();

    if (BRECOTYPE == pCand->fType) {
      for (int i = 0; i < NTRKMAX; ++i) {
	fDoca[i] = fDoca[i] = -1.;
      }
      int itrk(-1), nn(0);
      for (int it = pCand->fSig1; it <= pCand->fSig2; ++it) {
	TAnaTrack *p0 = fpEvt->getSigTrack(it);
	int found(0);
	for (unsigned int ig = 0; ig < fGenIndices.size(); ++ig) {
	  if (p0->fGenIndex == fGenIndices[ig]) {
	    found = 1;
	    break;
	  }
	}
	if (fVerbose > 1) cout << p0->fIndex << " pT = " << p0->fPlab.Perp() << " eta = " << p0->fPlab.Eta() << " fDouble1 = " << p0->fDouble1
	     << " mcidx = " << p0->fGenIndex
	     << " found = " << found
	     << endl;
	if (found) {
	  ++itrk;
	  fDoca[itrk] = p0->fDouble1;
	  fCorrect[itrk] = true;
	  if (TMath::Abs(p0->fMCID) != 13) {
	    ((TH1D*)fHistDir->Get("doca0"))->Fill(p0->fDouble1);
	  }
	} else {
	  ++itrk;
	  fDoca[itrk] = p0->fDouble1;
	  fCorrect[itrk] = false;
	  ((TH1D*)fHistDir->Get("doca2"))->Fill(p0->fDouble1);
	}
	((TH1D*)fHistDir->Get("doca1"))->Fill(p0->fDouble1);
	++nn;
      }
      fJpsiFlsxy = pCand->fVtx.fDxy/pCand->fVtx.fDxyE;
      fNTrk = itrk;
    }
  }
}


// ----------------------------------------------------------------------
void candAnaRecoil::brecoAnalysis() {
  TAnaTrack *pS(0);
  vector<int> trueTrkIdx;
  fNTrk = 0;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    // -- this is the "true"/cheating breco cand
    if (BRECOTYPE == pCand->fType) {
      cout << "===> BRECOTYPE" << endl;
      pCand->dump();
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	pS = fpEvt->getSigTrack(is);
	trueTrkIdx.push_back(pS->fIndex);
	pS->dump();
      }
    }


    // -- study the JpsiX BRECO candidate
    if (1) break;

    int itrk(0);
    for (int i = 0; i < NTRKMAX; ++i) {
      fDoca[i] = fProb1[i] = fProb2[i] = fChi2[i] = fDof[i] = fTrkMass[i] = -1.;
      fCorrect[i] = false;
    }

    if (1 == pCand->fType) {
      cout << "===> TYPE = 1" << endl;
      pCand->dump();
      int icnt(0);
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	pS = fpEvt->getSigTrack(is);
	cout << pS->fIndex;
	if (trueTrkIdx.end() != find(trueTrkIdx.begin(), trueTrkIdx.end(), pS->fIndex)) {
	  cout << "*" <<  "(" << pS->fDouble1 << ")";
	  if (icnt > 1) {
	    fDoca[itrk]    = pS->fDouble1;
	    fProb1[itrk]   = pS->fDouble2;
	    fProb2[itrk]   = pS->fDouble3;
	    fChi2[itrk]    = pS->fRefChi2;
	    fDof[itrk]     = pS->fRefDof;
	    fTrkMass[itrk] = pS->fBsTip;
	    fCorrect[itrk] = true;
	    ++itrk;
	  }
	} else {
	  cout << "(" << pS->fDouble1 << ")";
	  if (icnt > 1) {
	    fDoca[itrk]    = pS->fDouble1;
	    fProb1[itrk]   = pS->fDouble2;
	    fProb2[itrk]   = pS->fDouble3;
	    fChi2[itrk]    = pS->fRefChi2;
	    fDof[itrk]     = pS->fRefDof;
	    fTrkMass[itrk] = pS->fBsTip;
	    fCorrect[itrk] = false;
	    ++itrk;
	    if (0 && pS->fDouble1 < 0.02 && pS->fDouble3 > 0.7) {
	      cout << " doca = " << pS->fDouble1 << endl;
	      fpEvt->dump();
	    }
	  }
	}
	++icnt;
      }
      cout << endl;
      fNTrk = itrk;
      fJpsiFlsxy = pCand->fVtx.fDxy/pCand->fVtx.fDxyE;
    }


  }

}

// ----------------------------------------------------------------------
void candAnaRecoil::candAnalysis() {
  //  dumpBkmt();
  cout << "----------------------------------------------------------------------" << endl;
  cout << "candAnaRecoil::candAnalysis() event =  " << fChainEvent << endl;
  brecoAnalysis();

  TAnaTrack *pS(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    if (CANDTYPE == pCand->fType) {
      cout << "===> CANDTYPE" << endl;
      pCand->dump();
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	pS = fpEvt->getSigTrack(is);
	cout << pS->fIndex << " ";
      }
      cout << endl;
    }
  }

  candAna::candAnalysis();
}


// ----------------------------------------------------------------------
void candAnaRecoil::moreReducedTree(TTree *t) {
  t->Branch("flsxy",      &fJpsiFlsxy,  "jpsiflsxy/D");
  t->Branch("ntrk",   &fNTrk,   "ntrk/I");
  t->Branch("correct",&fCorrect,"correct[ntrk]/O");
  t->Branch("doca",   &fDoca,   "doca[ntrk]/D");
  t->Branch("prob1",  &fProb1,  "prob1[ntrk]/D");
  t->Branch("prob2",  &fProb2,  "prob2[ntrk]/D");
  t->Branch("chi2",   &fChi2,   "chi2[ntrk]/D");
  t->Branch("dof",    &fDof,    "dof[ntrk]/D");
  t->Branch("trkmass",&fTrkMass,"trkmass[ntrk]/D");
}

// ----------------------------------------------------------------------
void candAnaRecoil::bookHist() {
  cout << "==>candAnaRecoil: bookHist" << endl;
  fHistDir->cd();
  new TH1D("mrecoil", "mrecoil", 50, 5., 6.);
  new TH1D("m", "m (breco)", 50, 5., 6.);
  new TH1D("doca0", "doca(cand trk, J/psi vertex)", 100, 0., 0.1);
  new TH1D("doca1", "doca(all trk, J/psi vertex)", 100, 0., 0.1);
  new TH1D("doca2", "doca(other trk, J/psi vertex)", 100, 0., 0.1);

  new TH2D("candPv", "candPv", 10, -1., 9., 10, -1., 9.);

  new TH1D("brecoPv", "brecoPv", 10, -1., 9.);
  new TH1D("recoilPv", "recoilPv", 10, -1., 9.);

  new TH1D("recoilTRPVz", "recoil: truth - reco PV z ", 100, -1., 1.);
  new TH1D("brecoTRPVz", "breco: truth - reco PV z ", 100, -1., 1.);

  new TH2D("recoilTRPVzVsY", "recoil: truth - reco PV z vs y(B)", 40, -2.0, 2.0, 100, -1., 1.);
  new TH2D("brecoTRPVzVsY", "breco: truth - reco PV z vs y(B)", 40, -2.0, 2.0, 100, -1., 1.);

  candAna::bookHist();
  moreReducedTree(fTree);
}


// ----------------------------------------------------------------------
void candAnaRecoil::readCuts(string filename, int dump) {

  candAna::readCuts(filename, dump);


  fCutFile = filename;

  if (dump) cout << "==> candAnaRecoil: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines;
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  int ibin;
  string cstring = "B cand";

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str());

    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "BRECOTYPE")) {
      BRECOTYPE = static_cast<int>(CutValue);
      if (dump) cout << "BRECOTYPE:      " << BRECOTYPE << endl;
      ibin = 210;
      hcuts->SetBinContent(ibin, BRECOTYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BRECOTYPE :: %d", CutName, BRECOTYPE));
    }

    if (!strcmp(CutName, "BRECOTRUTH")) {
      BRECOTRUTH = static_cast<int>(CutValue);
      if (dump) cout << "BRECOTRUTH:      " << BRECOTRUTH << endl;
      ibin = 211;
      hcuts->SetBinContent(ibin, BRECOTRUTH);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BRECOTRUTH :: %d", CutName, BRECOTRUTH));
    }

    if (!strcmp(CutName, "CANDTYPE")) {
      CANDTYPE = static_cast<int>(CutValue);
      if (dump) cout << "CANDTYPE:      " << CANDTYPE << endl;
      ibin = 212;
      hcuts->SetBinContent(ibin, CANDTYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BRECOTYPE :: %d", CutName, CANDTYPE));
    }

    if (!strcmp(CutName, "CANDTRUTH")) {
      CANDTRUTH = static_cast<int>(CutValue);
      if (dump) cout << "CANDTRUTH:      " << CANDTRUTH << endl;
      ibin = 213;
      hcuts->SetBinContent(ibin, CANDTRUTH);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: CANDTRUTH :: %d", CutName, CANDTRUTH));
    }
  }

}


// ----------------------------------------------------------------------
void candAnaRecoil::printGenBDecays() {

  TGenCand *pCand, *pD;
  cout << "--- Event " << fChainEvent << " ------------------------------------------------------------------" << endl;
  cout << "gen block with " << fpEvt->nGenCands() << " gen cands" << endl;
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if ((TMath::Abs(pCand->fID) == 531) || (TMath::Abs(pCand->fID) == 511) || (TMath::Abs(pCand->fID) == 521)) {
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
bool candAnaRecoil::decayModeValidation(TGenCand *pCand, int mode) {
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
