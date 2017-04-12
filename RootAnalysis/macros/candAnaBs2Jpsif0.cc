#include "candAnaBs2Jpsif0.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"

#include "common/AnalysisDistribution.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBs2Jpsif0::candAnaBs2Jpsif0(bmmReader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  fAnaCuts.setAcName("candAnaBs2Jpsif0");
  fGenPi1Tmi = fGenPi2Tmi = fRecPi1Tmi = fRecPi2Tmi = -1;
  BLIND = 0;
  cout << "==> candAnaBs2Jpsif0: name = " << name << ", reading cutsfile " << cutsFile << endl;
  readCuts(cutsFile, 1);
}


// ----------------------------------------------------------------------
candAnaBs2Jpsif0::~candAnaBs2Jpsif0() {
  cout << "==> candAnaBs2Jpsif0: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::candAnalysis() {

  fGoodTracks     = false;
  fGoodTracksPt   = false;
  fGoodTracksEta  = false;

  fGoodAcceptance = false;
  fGoodJpsiCuts   = false;


  if (0 == fpCand) return;

  // -- check for overlap with a Bs -> J/psi phi candidate
  vector<int> idx0, idx1;
  getSigTracks(idx0, fpCand);
  int overlap(0);
  fBsJpsiPhiMass = fBdJpsiKstarMass = -99.;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pC = fpEvt->getCand(iC);
    if ((300531 == pC->fType) || (300511 == pC->fType)) {
      // stay in the loop
    } else {
      continue;
    }
    idx1.clear();
    getSigTracks(idx1, pC);
    // -- check for the same tracks
    overlap = 0;
    for (unsigned int i0 = 0; i0 < idx0.size(); ++i0) {
      for (unsigned int i1 = 0; i1 < idx1.size(); ++i1) {
	if (idx0[i0] == idx1[i1]) {
	  ++overlap;
	}
      }
    }
    // -- if all 4 track overlap, get kstar mass of the other cand
    if (4 == overlap) {
      if (300531 == pC->fType) fBsJpsiPhiMass   = pC->fMass;
      if (300511 == pC->fType) fBdJpsiKstarMass = pC->fMass;
    }
  }

  // -- Check for J/psi mass
  TAnaCand *pD = 0;
  fGoodJpsiMass = false;
  double chi2(0.), ndof(0.);
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i);
    //    cout << "i = " << i << " pD = " << pD << endl;
    if (pD->fType == JPSITYPE) {
      if ((JPSIMASSLO < pD->fMass) && (pD->fMass < JPSIMASSHI)) fGoodJpsiMass = true;
      fJpsiMass = pD->fMass;
      fJpsiPt   = pD->fPlab.Perp();
      fJpsiEta  = pD->fPlab.Eta();
      fJpsiPhi  = pD->fPlab.Phi();

      fJpsiCosA = TMath::Cos(pD->fAlpha);
      fJpsiMaxDoca = pD->fMaxDoca;
      fJpsiFLSxy   = pD->fVtx.fDxy/pD->fVtx.fDxyE;
      fJpsiVtxProb = pD->fVtx.fProb;

      chi2 = pD->fVtx.fChi2;
      ndof = pD->fVtx.fNdof;
    }
    if (pD->fType == F0TYPE) {
      if ((MPIPILO < pD->fMass) && (pD->fMass < MPIPIHI)) fGoodMPIPI = true;
      fMPiPi  = pD->fMass;
      ff0Pt   = pD->fPlab.Perp();
      ff0Eta  = pD->fPlab.Eta();
      ff0Phi  = pD->fPlab.Phi();
    }
  }

  // -- Get Kaons
  TAnaTrack *p0;
  TAnaTrack *p1(0);
  TAnaTrack *p2(0);

  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);
    if (TMath::Abs(p0->fMCID) != 211) continue;
    if (0 == p1) {
      p1 = p0;
    } else {
      p2 = p0;
    }
  }

  if (0 == p1) {
    cout << "candAnaBs2Jpsif0::candAnalysis  no pion 1 found " << endl;
    return;
  }
  if (0 == p2) {
    cout << "candAnaBs2Jpsif0::candAnalysis  no pion 2 found " << endl;
    return;
  }

  // -- order the pions according to (refitted) track pT
  if (p2->fRefPlab.Perp() > p1->fRefPlab.Perp()) {
    p0 = p2;
    p2 = p1;
    p1 = p0;
  }

  fPi1Pt        = p1->fRefPlab.Perp();
  fPi1Eta       = p1->fRefPlab.Eta();
  fPi1Phi       = p1->fRefPlab.Phi();
  fPi1TkQuality = highPurity(p1);
  fPi1PtNrf     = p1->fPlab.Perp();
  fPi1EtaNrf    = p1->fPlab.Eta();

  fPi2Pt        = p2->fRefPlab.Perp();
  fPi2Eta       = p2->fRefPlab.Eta();
  fPi2Phi       = p2->fRefPlab.Phi();
  fPi2TkQuality = highPurity(p2);
  fPi2PtNrf     = p2->fPlab.Perp();
  fPi2EtaNrf    = p2->fPlab.Eta();

  if (fCandTmi > -1 && fCandTmi == fpCand->fIndex) {
    TGenCand *pg1 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(p1->fIndex)->getGenIndex());
    fPi1PtGen     = pg1->fP.Perp();
    fPi1EtaGen    = pg1->fP.Eta();
    TGenCand *pg2 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(p2->fIndex)->getGenIndex());
    fPi2PtGen     = pg2->fP.Perp();
    fPi2EtaGen    = pg2->fP.Eta();
  } else {
    fPi1PtGen     = -99.;
    fPi1EtaGen    = -99.;
    fPi2PtGen     = -99.;
    fPi2EtaGen    = -99.;
  }

  ff0DeltaR  = p1->fPlab.DeltaR(p2->fPlab);

  TLorentzVector ka1, ka2, pi1, pi2;
  ka1.SetPtEtaPhiM(fPi1Pt, fPi1Eta, fPi1Phi, MKAON);
  ka2.SetPtEtaPhiM(fPi2Pt, fPi2Eta, fPi2Phi, MKAON);
  pi1.SetPtEtaPhiM(fPi1Pt, fPi1Eta, fPi1Phi, MPION);
  pi2.SetPtEtaPhiM(fPi2Pt, fPi2Eta, fPi2Phi, MPION);

  TLorentzVector f0Cand = ka1 + ka2;
  ff0Pt   = f0Cand.Pt();
  ff0Eta  = f0Cand.Eta();
  ff0Phi  = f0Cand.Phi();

  fGoodDeltaR   = (ff0DeltaR < DELTAR);
  fGoodMPIPI    = ((MPIPILO < fMPiPi ) && (fMPiPi < MPIPIHI));

  candAna::candAnalysis();


  // -- overwrite specific variables
  fCandChi2    = chi2;
  fCandDof     = ndof;
  fCandChi2Dof = chi2/ndof;

  fGoodTracks    = fGoodTracks    && fPi1TkQuality && fPi2TkQuality;
  fGoodTracksPt  = fGoodTracksPt  && ((TRACKPTLO < fPi1Pt)   && (fPi1Pt < TRACKPTHI)   && (TRACKPTLO < fPi2Pt)   && (fPi2Pt < TRACKPTHI));
  fGoodTracksEta = fGoodTracksEta && ((TRACKETALO < fPi1Eta) && (fPi1Eta < TRACKETAHI) && (TRACKETALO < fPi2Eta) && (fPi2Eta < TRACKETAHI));

  fGoodAcceptance = fGoodAcceptance && fGoodTracks    && fGoodTracksPt && fGoodTracksEta;
  fGoodJpsiCuts   = fGoodJpsiMass   && fGoodMPIPI     && fGoodDeltaR;
  //  fGoodJpsiCuts   = fGoodJpsiCuts   && (fJpsiPt > 7.) && (fJpsiCosA > 0.9) &&  (fJpsiFLSxy > 3.) && (fJpsiVtxProb > 0.1);
  fGoodJpsiCuts   = fGoodJpsiCuts   && (fJpsiPt > 7.) && (fJpsiVtxProb > 0.1);

  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(10);
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(4);
}

// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::moreBasicCuts() {
  cout << "   candAnaBs2Jpsif0: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodJpsiMass", "m(J/psi)", fGoodJpsiMass);
  fAnaCuts.addCut("fGoodDeltaR", "Delta R(KK)", fGoodDeltaR);
  fAnaCuts.addCut("fGoodMPIPI", "m(PiPi) [GeV]", fGoodMPIPI);
}

// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::genMatch() {

  fGenM1Tmi = fGenM2Tmi = fGenPi1Tmi = fGenPi2Tmi = -1;
  fNGenPhotons = 0;

  TGenCand *pC(0), *pB(0), *pPsi(0), *pf0(0), *pM1(0), *pM2(0), *pK1(0), *pK2(0), *pTmp(0);
  int ngamma(0);
  bool goodMatch(false);
  for (int i = 0; i < fpEvt->nGenT(); ++i) {
    pC = fpEvt->getGenT(i);
    if (531 == TMath::Abs(pC->fID)) {
      pB = pC;
      ngamma = 0;
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenTWithIndex(id);
	if (22 == TMath::Abs(pC->fID)) ++ngamma;
	if (443 == TMath::Abs(pC->fID)) {
	  pPsi = pC;
	  pM1 = pM2 = 0;
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd);
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenTWithIndex(idd);
	      } else {
		pM2 = fpEvt->getGenTWithIndex(idd);
	      }
	    }
	  }
	} else if (332 == TMath::Abs(pC->fID)) {
	  pf0 = fpEvt->getGenTWithIndex(id);
	  pK1 = pK2 = 0;
	  for (int idd = pf0->fDau1; idd <= pf0->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd);
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (211 == TMath::Abs(pC->fID)) {
	      if (0 == pK1) {
		pK1 = fpEvt->getGenTWithIndex(idd);
	      } else {
		pK2 = fpEvt->getGenTWithIndex(idd);
	      }
	    }
	  }
	}
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK1 && 0 != pK2 && (pPsi->fMom1 == pf0->fMom1)) {
	// -- check that there are no other direct daughters than J/psi K (plus possibly photons)
	int nDaughters(0);
	for (int ij = 0; ij < fpEvt->nGenT(); ++ij) {
	  pTmp = fpEvt->getGenT(ij);
	  if (pTmp->fMom1 == pB->fNumber) {
	    if (pTmp->fID != 22) ++nDaughters;
	  }
	}
	if (2 == nDaughters) {
	  goodMatch = true;
	  fNGenPhotons = ngamma;
	  break;
	}
      }
    }
  }

  if (!goodMatch) {
    if (fVerbose > 2) cout << "No matched signal decay found" << endl;
    return;
  }

  fGenBTmi = -1;
  fPi1GenID = -99999;
  fPi2GenID = -99999;
  if (goodMatch) {
    fMu1GenID = pM1->fID;
    fMu2GenID = pM2->fID;
    fPi1GenID = pK1->fID;
    fPi2GenID = pK2->fID;
    fGenBTmi = pB->fNumber;
    double m = pB->fP.Mag();
    double p = pB->fP.P();
    // mother pointer
    TGenCand *pM = fpEvt->getGenTWithIndex(pB->fMom1);
    // use the mother if it has the same PDGID (it oscillated)
    if (TMath::Abs(pB->fID) != TMath::Abs(pM->fID)) pM = pB;
    double x = (pM1->fV - pM->fV).Mag();
    fGenFl3d = x;
    fGenLifeTime = x*m/p/TMath::Ccgs();
    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fGenM1Tmi = pM1->fNumber;
      fGenM2Tmi = pM2->fNumber;
    } else {
      fGenM1Tmi = pM2->fNumber;
      fGenM2Tmi = pM1->fNumber;
    }
    if (pK1->fP.Perp() > pK2->fP.Perp()) {
      fGenPi1Tmi = pK1->fNumber;
      fGenPi2Tmi = pK2->fNumber;
    } else {
      fGenPi1Tmi = pK2->fNumber;
      fGenPi2Tmi = pK1->fNumber;
    }
  } else {
    fGenM1Tmi = -1;
    fGenM2Tmi = -1;
    fGenPi1Tmi = -1;
    fGenPi2Tmi = -1;
  }

}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::genMatchOld() {

  fGenM1Tmi = fGenM2Tmi = fGenPi1Tmi = fGenPi2Tmi = -1;
  fNGenPhotons = 0;

  TGenCand *pC(0), *pB(0), *pPsi(0), *pf0(0), *pM1(0), *pM2(0), *pPi1(0), *pPi2(0);
  int nb(0), ngamma(0);
  bool goodMatch(false);
  for (int i = 0; i < fpEvt->nGenT(); ++i) {
    pC = fpEvt->getGenT(i);
    if (531 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = pB->fDau2 - pB->fDau1 + 1;
      if (nb > 2) continue; // skip B decays where more than J/psi and f0 came from Bs
      ngamma = 0;
      for (int id = pB->fDau1; id <= pB->fDau2; ++id) {
	pC = fpEvt->getGenTWithIndex(id);
	if (22 == TMath::Abs(pC->fID)) ++ngamma;
	if (443 == TMath::Abs(pC->fID)) {
	  pPsi = pC;
	  pM1 = pM2 = 0;
	  for (int idd = pPsi->fDau1; idd <= pPsi->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd);
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (13 == TMath::Abs(pC->fID)) {
	      if (0 == pM1) {
		pM1 = fpEvt->getGenTWithIndex(idd);
	      } else {
		pM2 = fpEvt->getGenTWithIndex(idd);
	      }
	    }
	  }
	} else if (332 == TMath::Abs(pC->fID)) {
	  pf0 = fpEvt->getGenTWithIndex(id);
	  pPi1 = pPi2 = 0;
	  for (int idd = pf0->fDau1; idd <= pf0->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd);
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (211 == TMath::Abs(pC->fID)) {
	      if (0 == pPi1) {
		pPi1 = fpEvt->getGenTWithIndex(idd);
	      } else {
		pPi2 = fpEvt->getGenTWithIndex(idd);
	      }
	    }
	  }
	}
      }
      if (0 != pM1 && 0 != pM2 && 0 != pPi1 && 0 != pPi2
	  && (pPsi->fMom1 == pf0->fMom1)
	  ) {
	goodMatch = true;
	fNGenPhotons = ngamma;
	break;
      }
    }
  }

  if (!goodMatch) {
    if (fVerbose > 2) cout << "No matched signal decay found" << endl;
    return;
  }

  fGenBTmi = -1;
  fPi1GenID = -99999;
  fPi2GenID = -99999;
  if (goodMatch) {
    fMu1GenID = pM1->fID;
    fMu2GenID = pM2->fID;
    fPi1GenID = pPi1->fID;
    fPi2GenID = pPi2->fID;
    fGenBTmi = pB->fNumber;
    double m = pB->fP.Mag();
    double p = pB->fP.P();
    // Meson pointer
    TGenCand *pM = fpEvt->getGenTWithIndex(pB->fMom1);
    // the meson is the original except if it oscillated
    if (531 != TMath::Abs(pM->fID)) pM = pB;
    double x = (pM1->fV - pM->fV).Mag();
    fGenFl3d = x;
    fGenLifeTime = x*m/p/TMath::Ccgs();
    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fGenM1Tmi = pM1->fNumber;
      fGenM2Tmi = pM2->fNumber;
    } else {
      fGenM1Tmi = pM2->fNumber;
      fGenM2Tmi = pM1->fNumber;
    }
    if (pPi1->fP.Perp() > pPi2->fP.Perp()) {
      fGenPi1Tmi = pPi1->fNumber;
      fGenPi2Tmi = pPi2->fNumber;
    } else {
      fGenPi1Tmi = pPi2->fNumber;
      fGenPi2Tmi = pPi1->fNumber;
    }
  } else {
    fGenM1Tmi = -1;
    fGenM2Tmi = -1;
    fGenPi1Tmi = -1;
    fGenPi2Tmi = -1;
  }

}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = fRecPi1Tmi = fRecPi2Tmi =-1;
  TSimpleTrack *pT(0);
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i);
    if (pT->getGenIndex() < 0) continue;

    // -- muon 1
    if (fGenM1Tmi > -1 && pT->getGenIndex() == fGenM1Tmi) {
      fRecM1Tmi = i;
    }

    // -- muon 2
    if (fGenM2Tmi > -1 && pT->getGenIndex() == fGenM2Tmi) {
      fRecM2Tmi = i;
    }

    // -- kaon 1
    if (fGenPi1Tmi > -1 && pT->getGenIndex() == fGenPi1Tmi) {
      fRecPi1Tmi = i;
    }

    // -- kaon 2
    if (fGenPi2Tmi > -1 && pT->getGenIndex() == fGenPi2Tmi) {
      fRecPi2Tmi = i;
    }

    // -- skip rest if all matches found
    if (fRecM1Tmi > -1 && fRecM2Tmi > -1 && fRecPi1Tmi > -1 && fRecPi2Tmi > -1) break;
  }


  if (fVerbose > 10) {
    cout << "fRecM1Tmi =  " << fRecM1Tmi  << " matched to fGenM1Tmi =  " << fGenM1Tmi << endl;
    cout << "fRecM2Tmi =  " << fRecM2Tmi  << " matched to fGenM2Tmi =  " << fGenM2Tmi << endl;
    cout << "fRecPi1Tmi = " << fRecPi1Tmi << " matched to fGenPi1Tmi = " << fGenPi1Tmi << endl;
    cout << "fRecPi2Tmi = " << fRecPi2Tmi << " matched to fGenPi2Tmi = " << fGenPi2Tmi << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::candMatch() {

  fCandTmi = -1;
  int idx(-1), type(-1);
  int d1Matched(0), d2Matched(0), d3Matched(0), d4Matched(0);
  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);
    if (TYPE != pCand->fType) continue;

    d1Matched = d2Matched = d3Matched = d4Matched = 0;
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex;
      type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
      //       if (fGenM1Tmi > -1) cout << "  --> " << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecK1Tmi << endl;
      if (fVerbose > 10) {
	cout << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecPi1Tmi << " " << fRecPi2Tmi << endl;
      }
      if (fRecM1Tmi > -1 && type == 13 && idx == fRecM1Tmi) {
	d1Matched = 1;
      }
      if (fRecM2Tmi > -1 && type == 13 && idx == fRecM2Tmi) {
	d2Matched = 1;
      }
      if (fRecPi1Tmi > -1 && type == 211 && idx == fRecPi1Tmi) {
	d3Matched = 1;
      }
      if (fRecPi2Tmi > -1 && type == 211 && idx == fRecPi2Tmi) {
	d4Matched = 1;
      }
    }

    if (d1Matched && d2Matched && d3Matched && d4Matched) {
      fCandTmi = iC;
      break;
    }
  }
  if (fVerbose > 10) {
    cout << "fCandTmi = " << fCandTmi << " matched to rec tracks " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecPi1Tmi  << " " << fRecPi2Tmi
	 << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::bookHist() {
  cout << "==>candAnaBs2Jpsif0: bookHist" << endl;
  candAna::bookHist();

  moreReducedTree(fTree);

  // -- Additional effTree variables
  fEffTree->Branch("pi1pt",  &fETpi1pt,           "pi1pt/F");
  fEffTree->Branch("g3pt",   &fETg3pt,            "g3pt/F");
  fEffTree->Branch("pi1eta", &fETpi1eta,          "pi1eta/F");
  fEffTree->Branch("g3eta",  &fETg3eta,           "g3eta/F");
  fEffTree->Branch("pi1q",   &fETpi1q,            "pi1q/I");
  fEffTree->Branch("pi1gt",  &fETpi1gt,           "pi1gt/O");

  fEffTree->Branch("pi2pt",  &fETpi2pt,           "pi2pt/F");
  fEffTree->Branch("g4pt",   &fETg4pt,            "g4pt/F");
  fEffTree->Branch("pi2eta", &fETpi2eta,          "pi2eta/F");
  fEffTree->Branch("g4eta",  &fETg4eta,           "g4eta/F");
  fEffTree->Branch("pi2q",   &fETpi2q,            "pi2q/I");
  fEffTree->Branch("pi2gt",  &fETpi2gt,           "pi2gt/O");

}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::moreReducedTree(TTree *t) {

  // -- Additional reduced tree variables
  t->Branch("mpsi",        &fJpsiMass,    "mpsi/D");
  t->Branch("psipt",       &fJpsiPt,      "psipt/D");
  t->Branch("psieta",      &fJpsiEta,     "psieta/D");
  t->Branch("psiphi",      &fJpsiPhi,     "psiphi/D");
  t->Branch("psicosa",     &fJpsiCosA,    "psicosa/D");
  t->Branch("psimaxdoca",  &fJpsiMaxDoca, "psimaxdoca/D");
  t->Branch("psiflsxy",    &fJpsiFLSxy,   "psiflsxy/D");
  t->Branch("psiprob",     &fJpsiVtxProb, "psiprob/D");

  t->Branch("mpipi", &fMPiPi,     "mpipi/D");
  t->Branch("f0pt",  &ff0Pt,      "f0pt/D");
  t->Branch("f0eta", &ff0Eta,     "f0eta/D");
  t->Branch("f0phi", &ff0Phi,     "f0phi/D");
  t->Branch("f0dr",  &ff0DeltaR,  "f0dr/D");
  t->Branch("mbspsiphi", &fBsJpsiPhiMass,  "mbspsiphi/D");
  t->Branch("mbdpsikstar", &fBdJpsiKstarMass,  "mbdpsikstar/D");



  t->Branch("pi1pt",  &fPi1Pt,       "pi1pt/D");
  t->Branch("pi1eta", &fPi1Eta,      "pi1eta/D");
  t->Branch("pi1phi", &fPi1Phi,      "pi1phi/D");
  t->Branch("pi1gt",  &fPi1TkQuality,"pi1gt/I");
  t->Branch("pi2pt",  &fPi2Pt,       "pi2pt/D");
  t->Branch("pi2eta", &fPi2Eta,      "pi2eta/D");
  t->Branch("pi2phi", &fPi2Phi,      "pi2phi/D");
  t->Branch("pi2gt",  &fPi2TkQuality,"pi2gt/I");



  t->Branch("t3pt",  &fPi1PtNrf, "t3pt/D");
  t->Branch("t3eta", &fPi1EtaNrf,"t3eta/D");

  t->Branch("t4pt",  &fPi2PtNrf, "t4pt/D");
  t->Branch("t4eta", &fPi2EtaNrf,"t4eta/D");

  t->Branch("g3pt", &fPi1PtGen,  "g3pt/D");
  t->Branch("g3eta",&fPi1EtaGen, "g3eta/D");
  t->Branch("g3id", &fPi1GenID,  "g3id/I");

  t->Branch("g4pt", &fPi2PtGen,  "g4pt/D");
  t->Branch("g4eta",&fPi2EtaGen, "g4eta/D");
  t->Branch("g4id", &fPi2GenID,  "g4id/I");
}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::fillCandidateHistograms(int offset) {
  candAna::fillCandidateHistograms(offset);
}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::efficiencyCalculation() {
  // -- gen level
  TGenCand *pB(0), *pM1(0), *pM2(0), *pPi1(0), *pPi2(0);
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi || -1 == fGenPi1Tmi || -1 == fGenPi2Tmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }
  pB  = fpEvt->getGenTWithIndex(fGenBTmi);
  pM1 = fpEvt->getGenTWithIndex(fGenM1Tmi);
  pM2 = fpEvt->getGenTWithIndex(fGenM2Tmi);
  pPi1 = fpEvt->getGenTWithIndex(fGenPi1Tmi);
  pPi2 = fpEvt->getGenTWithIndex(fGenPi2Tmi);

  // -- reco level
  TSimpleTrack *prM1(0), *prM2(0), *prPi1(0), *prPi2(0);
  double bla(0);
  int m1Matched(0), m2Matched(0), pi1Matched(0), pi2Matched(0), m1ID(0), m1mvaID(0), m2ID(0), m2mvaID(0),
    m1GT(0), m2GT(0), pi1GT(0), pi2GT(0);
  if (fRecM1Tmi > -1) {
    m1Matched = 1;
    prM1 = fpEvt->getSimpleTrack(fRecM1Tmi);
    if (mvaMuon(prM1, bla)) m1mvaID = 1;
    if (prM1->getHighPurity()) {
      m1GT = 1;
    } else {
      m1GT = 0;
    }
  }

  if (fRecM2Tmi > -1) {
    m2Matched = 1;
    prM2 = fpEvt->getSimpleTrack(fRecM2Tmi);
    if (mvaMuon(prM2, bla)) m2mvaID = 1;
    if (prM2->getHighPurity()) {
      m2GT = 1;
    } else {
      m2GT = 0;
    }
  }

  if (fRecPi1Tmi > -1) {
    pi1Matched = 1;
    prPi1 = fpEvt->getSimpleTrack(fRecPi1Tmi);
    if (prPi1->getHighPurity()) {
      pi1GT = 1;
    } else {
      pi1GT = 0;
    }
  }

  if (fRecPi2Tmi > -1) {
    pi2Matched = 1;
    prPi2 = fpEvt->getSimpleTrack(fRecPi2Tmi);
    if (prPi2->getHighPurity()) {
      pi2GT = 1;
    } else {
      pi2GT = 0;
    }
  }

  // -- cand level
  TAnaCand *pCand(0);
  if (fCandTmi > -1) {
    pCand = fpEvt->getCand(fCandTmi);
  }

  m1ID = m1mvaID;
  m2ID = m2mvaID;

  // -- EffTree filling for all events with a signal decay
  fETgm    = pB->fP.M();
  fETgpt   = pB->fP.Perp();
  fETgtau  = fGenLifeTime;
  fETgeta  = pB->fP.Eta();
  fETg1pt  = pM1->fP.Perp();
  fETg1eta = pM1->fP.Eta();
  fETg2pt  = pM2->fP.Perp();
  fETg2eta = pM2->fP.Eta();
  fETg3pt  = pPi1->fP.Perp();
  fETg3eta = pPi1->fP.Eta();
  fETg4pt  = pPi2->fP.Perp();
  fETg4eta = pPi2->fP.Eta();
  if (m1Matched) {
    fETm1pt  = prM1->getP().Perp();
    fETm1eta = prM1->getP().Eta();
    fETm1q   = prM1->getCharge();
    fETm1gt  = (m1GT>0?true:false);
    fETm1id  = (m1ID>0?true:false);
    fETm1mvaid = (m1mvaID>0?true:false);
  } else {
    fETm1pt  = -99.;
    fETm1eta = -99.;
    fETm1q   = -99;
    fETm1gt  = false;
    fETm1id  = false;
    fETm1mvaid = false;
  }
  if (m2Matched) {
    fETm2pt  = prM2->getP().Perp();
    fETm2eta = prM2->getP().Eta();
    fETm2q   = prM2->getCharge();
    fETm2gt  = (m2GT>0?true:false);
    fETm2id  = (m2ID>0?true:false);
    fETm2mvaid = (m2mvaID>0?true:false);
  } else {
    fETm2pt  = -99.;
    fETm2eta = -99.;
    fETm2q   = -99;
    fETm2gt  = false;
    fETm2id  = false;
    fETm2mvaid = false;
  }
  if (m1Matched && m2Matched) {
    fETchan = detChan(fETm1eta, fETm2eta);
  } else {
    fETchan = -1;
  }
  if (pi1Matched) {
    fETpi1pt  = prPi1->getP().Perp();
    fETpi1eta = prPi1->getP().Eta();
    fETpi1q   = prPi1->getCharge();
    fETpi1gt  = (pi1GT>0?true:false);
  } else {
    fETpi1pt  = -99.;
    fETpi1eta = -99.;
    fETpi1q   = -99;
    fETpi1gt  = false;
  }
  if (pi2Matched) {
    fETpi2pt  = prPi2->getP().Perp();
    fETpi2eta = prPi2->getP().Eta();
    fETpi2q   = prPi2->getCharge();
    fETpi2gt  = (pi2GT>0?true:false);
  } else {
    fETpi2pt  = -99.;
    fETpi2eta = -99.;
    fETpi2q   = -99;
    fETpi2gt  = false;
  }
  if (pCand) {
    fETcandMass = pCand->fMass;
    fETtau      = pCand->fTau3d;
  } else {
    fETcandMass = -99.;
    fETtau      = -99.;
  }

  fEffTree->Fill();

}


// ----------------------------------------------------------------------
void candAnaBs2Jpsif0::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump);

  fCutFile = filename;

  if (dump) cout << "==> candAnaBs2Jpsif0: Reading " << fCutFile << " for cut settings" << endl;
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

    if (!strcmp(CutName, "JPSITYPE")) {
      JPSITYPE = static_cast<int>(CutValue);
      if (dump) cout << "JPSITYPE:      " << JPSITYPE << endl;
      ibin = 210;
      hcuts->SetBinContent(ibin, JPSITYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: J/#psi ID :: %d", CutName, JPSITYPE));
    }

    if (!strcmp(CutName, "JPSIMASSLO")) {
      JPSIMASSLO = CutValue;
      if (dump) cout << "JPSIMASSLO:      " << JPSIMASSLO << endl;
      ibin = 211;
      hcuts->SetBinContent(ibin, JPSIMASSLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m(J/#psi) :: %3.1f", CutName, JPSIMASSLO));
    }

    if (!strcmp(CutName, "JPSIMASSHI")) {
      JPSIMASSHI = CutValue;
      if (dump) cout << "JPSIMASSLO:      " << JPSIMASSHI << endl;
      ibin = 212;
      hcuts->SetBinContent(ibin, JPSIMASSHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m(J/#psi) :: %3.1f", CutName, JPSIMASSHI));
    }

    if (!strcmp(CutName, "F0TYPE")) {
      F0TYPE = static_cast<int>(CutValue);
      if (dump) cout << "F0TYPE:      " << F0TYPE << endl;
      ibin = 213;
      hcuts->SetBinContent(ibin, F0TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: f0 ID :: %d", CutName, F0TYPE));
    }

    if (!strcmp(CutName, "MPIPILO")) {
      MPIPILO = CutValue;
      if (dump) cout << "MPIPILO:           " << MPIPILO << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MPIPILO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{min}(#pi#pi) :: %3.1f", CutName, MPIPILO));
    }

    if (!strcmp(CutName, "MPIPIHI")) {
      MPIPIHI = CutValue;
      if (dump) cout << "MPIPIHI:           " << MPIPIHI << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MPIPIHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{max}(#pi#pi) :: %3.1f", CutName, MPIPIHI));
    }

    if (!strcmp(CutName, "DELTAR")) {
      DELTAR = CutValue;
      if (dump) cout << "DELTAR:           " << DELTAR << endl;
      ibin = 301;
      hcuts->SetBinContent(ibin, DELTAR);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #Delta R(KK) :: %3.1f", CutName, DELTAR));
    }



  }

}
