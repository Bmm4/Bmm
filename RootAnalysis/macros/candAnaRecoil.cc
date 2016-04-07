#include "candAnaRecoil.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaRecoil::candAnaRecoil(bmmReader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaRecoil: constructor..." << endl;
  readCuts(cutsFile, 1); 
}


// ----------------------------------------------------------------------
candAnaRecoil::~candAnaRecoil() {
  cout << "==> candAnaRecoil: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaRecoil::efficiencyCalculation() {
}

// ----------------------------------------------------------------------
// BRECO selection/reconstruction
void candAnaRecoil::evtAnalysis(TAna01Event *evt) {
  
  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(1);
  fpEvt = evt; 
  fBadEvent = false;

  if (fIsMC) {
    genMatch(); 
    recoMatch(); 
    candMatch(); 
    if (fBadEvent) {
      cout << "XXXXXXX BAD EVENT XXXXXX SKIPPING XXXXX" << endl;
      return;
    }
    efficiencyCalculation();
  } 

  triggerSelection();
  runRange(); 

  TAnaCand *pCand(0);
  bool fillNoCand(true); 
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    if (BRECOTYPE != pCand->fType) {
      if (fVerbose > 39) cout << "  skipping candidate at " << iC << " which is of type " << pCand->fType 
			      << " looking for type " << BRECOTYPE << endl;
      continue;
    }

    if (fVerbose > 2) cout << "  taking candidate at " << iC << " which is of type " << pCand->fType << endl;


    fillNoCand = false; 
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(3);
    fpBreco = pCand;
    fBrecoIdx = iC; 

    candAnalysis();
    if (0 == fpCand) {
      continue;
    }
    
    if (0) {
      if (fBrecoPt < 30 && fBrecoPhi > -2.9 && fBrecoPhi < -2.7 && fCandPhi > 0.35 && fCandPhi < 0.55) {
	cout << "breco pt = " << fBrecoPt << " phi = " << fBrecoPhi << " eta = " << fBrecoEta << endl;
      }
    }
    
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(5);
    if (fIsMC) {
      fTree->Fill(); 
    } else {  // DATA
      if (NOPRESELECTION) {
	fPreselection = true; 
      }
      
      if (fPreselection) { 
	fTree->Fill(); 
      }
    } 
  }  

  // -- fill events with no passing BRECO candidate (one entry per event)
  if (fillNoCand) ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(2);

}



// ----------------------------------------------------------------------
void candAnaRecoil::candAnalysis() {
  // -- fpCand is to be set. fpBreco is set in candAnaRecoil::evtAnalysis!!
  fpCand = 0; 
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    if (TYPE == pCand->fType) {
      fpCand = pCand;
      fCandIdx = iC; 
    } 
    if (fpCand) break;
  }

  if (!fpCand) {
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(4);
    return;
  }
  
  candAna::candAnalysis();
  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(5);


  
  ((TH1D*)fHistDir->Get("m"))->Fill(fpCand->fMass);
  ((TH1D*)fHistDir->Get("mrecoil"))->Fill(fpBreco->fMass);
  
  fBrecoMass  = fpBreco->fMass;
  fBrecoPt    = fpBreco->fPlab.Perp();
  fBrecoEta   = fpBreco->fPlab.Eta();
  fBrecoPhi   = fpBreco->fPlab.Phi();
  if (fpBreco->fPvIdx > -1 && fpBreco->fPvIdx < fpEvt->nPV()) {
    fBrecoPvIdx = fpBreco->fPvIdx;
    fBrecoPvZ   = fpEvt->getPV(fBrecoPvIdx)->fPoint.Z();
  } else {
    fBrecoPvIdx = -99; 
    fBrecoPvZ   = -99.;
  }
  // -- histogram the PV to which they are associated
  ((TH2D*)fHistDir->Get("candPv"))->Fill(fpBreco->fPvIdx, fpCand->fPvIdx);
  
  TAnaTrack *p0, *p1(0), *p2(0); 

  // -- daughters of fpCand
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     
    ((TH1D*)fHistDir->Get("recoilPv"))->Fill(p0->fPvIdx);
  }
  
  
  // -- get daughters of breco
  for (int it = fpBreco->fSig1; it <= fpBreco->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);     

    ((TH1D*)fHistDir->Get("brecoPv"))->Fill(p0->fPvIdx);

    if (TMath::Abs(p0->fMCID) != 13) continue;
    if (0 == p1) {
      p1 = p0; 
    } else {
      p2 = p0; 
    }
  }

  // -- order muons with pT
  if (p1->fPlab.Perp() < p2->fPlab.Perp()) {
    p0 = p1;
    p1 = p2;
    p2 = p0;
  }
  
  fMu1BrecoPt  = p1->fPlab.Perp(); 
  fMu1BrecoEta = p1->fPlab.Eta(); 
  fMu1BrecoPhi = p1->fPlab.Phi(); 
  fMu2BrecoPt  = p2->fPlab.Perp(); 
  fMu2BrecoEta = p2->fPlab.Eta(); 
  fMu2BrecoPhi = p2->fPlab.Phi(); 
  
  if (0) cout << "candidates found; pCand (recoil) = " << fpCand << " with m = " << fpCand->fMass
	      << " and fpBreco = " << fpBreco << " with m = " << fpBreco->fMass
	      << " fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi
	      << " fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi
	      << " fCandTmi = " << fCandTmi
	      << endl;
  
}

// ----------------------------------------------------------------------
void candAnaRecoil::processType() {

}


// ----------------------------------------------------------------------
void candAnaRecoil::genMatch() {
  fGenM1Tmi = fGenM2Tmi = -1; 
  fNGenPhotons = 0; 

  int id1(13), id2(13); 

  // -- modifications for rare backgrounds
  if (1000082 == BRECOTYPE) {id1 = 321; id2 = 321;}
  if (1000083 == BRECOTYPE) {id1 = 321; id2 = 211;}
  if (1000084 == BRECOTYPE) {id1 = 211; id2 = 211;}
  if (1000085 == BRECOTYPE) {id1 = 211; id2 = 13;}
  if (1000086 == BRECOTYPE) {id1 = 321; id2 = 13;}

  if (1000091 == BRECOTYPE) {id1 = 211; id2 = 211;}
  if (1000092 == BRECOTYPE) {id1 = 321; id2 = 211;}
  if (1000093 == BRECOTYPE) {id1 = 321; id2 = 321;}
  if (1000095 == BRECOTYPE) {id1 = 211; id2 = 13;}
  if (1000098 == BRECOTYPE) {id1 = 211; id2 = 211;}

  if (1000060 == BRECOTYPE) {id1 = 2212; id2 = 211;}
  if (1000061 == BRECOTYPE) {id1 = 2212; id2 = 321;}
  if (1000062 == BRECOTYPE) {id1 = 2212; id2 = 13;}

  TGenCand *pC(0), *pM1(0), *pM2(0), *pB(0); 
  bool goodMatch(false); 
  for (int i = 0; i < fpEvt->nGenT(); ++i) {
    pC = fpEvt->getGenT(i); 
    if (BRECOTRUTH == TMath::Abs(pC->fID)) {
      pM1 = pM2 = 0; 
      pB = pC;
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenTWithIndex(id); 
	if (id1 == TMath::Abs(pC->fID) || id2 == TMath::Abs(pC->fID)) {
	  if (0 == pM1) {
	    pM1 = fpEvt->getGenTWithIndex(id); 
	  } else {
	    pM2 = fpEvt->getGenTWithIndex(id); 
	  }
	}
      }
      if (0 != pM1 && 0 != pM2) {
	goodMatch = true; 
	fNGenPhotons = pB->fDau2 - pB->fDau1 - 1; 
	if (fVerbose > 10) {
	  cout << "found gen match for B gen idx = " << pB->fNumber << endl;
	}
	break;
      }
    }
  }

  
  fGenBTmi = -1; 
  if (goodMatch) {
    fGenBTmi = pB->fNumber; 
    double m = pB->fP.Mag();
    double p = pB->fP.P();
    // Meson pointer
    TGenCand *pM = fpEvt->getGenTWithIndex(pB->fMom1); 
    // the meson is the original except if it oscillated
    if (BRECOTRUTH != TMath::Abs(pM->fID)) pM = pB;
    double x = (pM1->fV - pM->fV).Mag(); 
    fGenLifeTime = x*m/p/TMath::Ccgs();

    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fMu1GenID = pM1->fID;
      fMu2GenID = pM2->fID;
      fGenM1Tmi = pM1->fNumber; 
      fGenM2Tmi = pM2->fNumber; 
    } else {
      fMu1GenID = pM2->fID;
      fMu2GenID = pM1->fID;
      fGenM1Tmi = pM2->fNumber; 
      fGenM2Tmi = pM1->fNumber; 
    }
  } else {
    fGenM1Tmi = -1; 
    fGenM2Tmi = -1; 
  }
  
  if (fVerbose > 10) {
    cout << "fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fGenM2Tmi = " << fGenM2Tmi << endl;

  }

}


// ----------------------------------------------------------------------
void candAnaRecoil::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = -1; 
  TSimpleTrack *pT(0);
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i); 
    if (pT->getGenIndex() < 0) continue;

    // -- muon 1
    if (pT->getGenIndex() == fGenM1Tmi) {
      fRecM1Tmi = i; 
    }

    // -- muon 2
    if (pT->getGenIndex() == fGenM2Tmi) {
      fRecM2Tmi = i; 
    }

    // -- skip rest if both matches found
    if (fRecM1Tmi > -1 && fRecM2Tmi > -1) break;
  }


  if (fVerbose > 10) {
    cout << "fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaRecoil::candMatch() {
  fCandTmi = -1;   
  int idx(-1); 
  int d1Matched(0), d2Matched(0); 
  TAnaCand *pCand(0);
  if (fRecM1Tmi < 0 ||  fRecM2Tmi < 0) return; 
  
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC); 
    if (BRECOTYPE != pCand->fType)  continue;
    
    d1Matched = d2Matched = 0; 
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex; 
      if (fVerbose > 10) {
	cout << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << endl;
      }
      if (idx == fRecM1Tmi) {
	d1Matched = 1; 
      }
      if (idx == fRecM2Tmi) {
	d2Matched = 1; 
      }
    }
    
    if (d1Matched && d2Matched) {
      fCandTmi = iC;
      break;
    }
  }
  if (fVerbose > 10) {
    cout << "fCandTmi = " << fCandTmi << endl;
  }
}



// ----------------------------------------------------------------------
void candAnaRecoil::bookHist() {
  cout << "==>candAnaRecoil: bookHist" << endl;
  fHistDir->cd();
  new TH1D("mrecoil", "mrecoil", 50, 5., 6.);
  new TH1D("m", "m (breco)", 50, 5., 6.);

  new TH2D("candPv", "candPv", 10, -1., 9., 10, -1., 9.);

  new TH1D("brecoPv", "brecoTrkPv", 10, -1., 9.);
  new TH1D("recoilPv", "recoilTrkPv", 10, -1., 9.);

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
      ibin = 210;
      hcuts->SetBinContent(ibin, BRECOTRUTH);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BRECOTRUTH :: %d", CutName, BRECOTRUTH));
    }

  }
  
}


// ----------------------------------------------------------------------
void candAnaRecoil::moreReducedTree(TTree *t) {
  // -- Additional reduced tree variables
  t->Branch("mbreco",     &fBrecoMass,    "mbreco/D");
  t->Branch("ptbreco",    &fBrecoPt,      "ptbreco/D");
  t->Branch("etabreco",   &fBrecoEta,     "etabreco/D");
  t->Branch("phibreco",   &fBrecoPhi,     "phibreco/D");
  t->Branch("pvidxbreco", &fBrecoPvIdx,   "pvidxbreco/I");
  t->Branch("pvzbreco",   &fBrecoPvZ,     "pvzbreco/D");
  t->Branch("m1ptbreco",  &fMu1BrecoPt,   "m1ptbreco/D");
  t->Branch("m1etabreco", &fMu1BrecoEta,  "m1etabreco/D");
  t->Branch("m1phibreco", &fMu1BrecoPhi,  "m1phibreco/D");
  t->Branch("m2ptbreco",  &fMu2BrecoPt,   "m2ptbreco/D");
  t->Branch("m2etabreco", &fMu2BrecoEta,  "m2etabreco/D");
  t->Branch("m2phibreco", &fMu2BrecoPhi,  "m2phibreco/D");

}
