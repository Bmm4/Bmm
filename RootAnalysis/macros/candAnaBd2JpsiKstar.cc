#include "candAnaBd2JpsiKstar.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"

#include "common/AnalysisDistribution.hh"

using namespace std;

// ----------------------------------------------------------------------
candAnaBd2JpsiKstar::candAnaBd2JpsiKstar(bmmReader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  fGenKTmi = fGenPiTmi = fRecKTmi = fRecPiTmi = -1;
  BLIND = 0;
  cout << "==> candAnaBd2JpsiKstar: name = " << name << ", reading cutsfile " << cutsFile << endl;
  readCuts(cutsFile, 1);
}


// ----------------------------------------------------------------------
candAnaBd2JpsiKstar::~candAnaBd2JpsiKstar() {
  cout << "==> candAnaBd2JpsiKstar: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::candAnalysis() {

  if (0 == fpCand) return;

  TAnaCand *pC(0), *pD(0);
  // -- Check for another candidate with the same tracks that is closer to the PDG K*0 mass
  double mkstar(0.), mkstarOther(0.);
  for (int i = fpCand->fDau1; i <= fpCand->fDau2; ++i) {
    if (i < 0) break;
    pD = fpEvt->getCand(i);
    if (pD->fType == KSTARTYPE) {
      mkstar     = pD->fMass;
      break;
    }
  }
  vector<int> idx0, idx1;
  getSigTracks(idx0, fpCand);
  int overlap(0);
  fKstarFail = false;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pC = fpEvt->getCand(iC);
    if (pC == fpCand) continue;
    if (3000070 != pC->fType) {
      if (pC->fType != fpCand->fType) continue;
    }
    idx1.clear();
    getSigTracks(idx1, pC);
    // -- check for the same tracks
    overlap = 0;
    for (unsigned int i0 = 0; i0 < idx0.size(); ++i0) {
      if (idx0[i0] == idx1[i0]) {
	++overlap;
      }
    }
    // -- if all 4 track overlap, get kstar mass of the other cand
    if (4 == overlap) {
      for (int i = pC->fDau1; i <= pC->fDau2; ++i) {
	if (i < 0) break;
	pD = fpEvt->getCand(i);
	if (pD->fType == KSTARTYPE) {
	  mkstarOther = pD->fMass;
	  break;
	}
      }
      if (TMath::Abs(mkstarOther - MKSTAR) < TMath::Abs(mkstar - MKSTAR)) {
	if (0) cout << "other cand " << pC->fIndex << " track indices has better kstar mass: " << mkstarOther << " (" << TMath::Abs(mkstarOther - MKSTAR)
		    << ") then this cand: " << mkstar  << " (" << TMath::Abs(mkstar - MKSTAR) << ")"
		    << endl;
	if (3000070 == pC->fType) {
	  fKstarFail = true;
	} else {
	  fKstarFail = true;
	  return;
	}
      }
    }
  }

  // -- Check for J/psi mass
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

    ///////added by jmonroy //////////////////////////////

    if (pD->fType == KSTARTYPE ) {                                            //is there a KSTARTYPE ???? I defined it in the .hh
      if ((MKPILO < pD->fMass) && (pD->fMass < MKPIHI)) fGoodMKPI = true;
      fMKPI     = pD->fMass;
      fKstarPt   = pD->fPlab.Perp();
      fKstarEta  = pD->fPlab.Eta();
      fKstarPhi  = pD->fPlab.Phi();
      }

  }


  // -- Get Kaons
  TAnaTrack *p0;
  TAnaTrack *p1(0);
  TAnaTrack *p2(0);

  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);
    if (0 == p0) {
      cout << "candAnaBd2JpsiKstar::candAnalysis problem with sigtrack?? " << endl;
      return;
    }
    if (321 == TMath::Abs(p0->fMCID)) {
      p1 = p0;
    }
    if (211 == TMath::Abs(p0->fMCID)) {
      p2 = p0;
    }
  }
  if (0 == p1) {
    cout << "candAnaBd2JpsiKstar::candAnalysis  no kaon found " << endl;
    return;
  }
  if (0 == p2) {
    cout << "candAnaBd2JpsiKstar::candAnalysis  no pion found " << endl;
    return;
  }

  fKaPt        = p1->fRefPlab.Perp();
  fKaEta       = p1->fRefPlab.Eta();
  fKaPhi       = p1->fRefPlab.Phi();
  fKaTkQuality = highPurity(p1);
  fKaPtNrf     = p1->fPlab.Perp();
  fKaEtaNrf    = p1->fPlab.Eta();

  fPiPt        = p2->fRefPlab.Perp();
  fPiEta       = p2->fRefPlab.Eta();
  fPiPhi       = p2->fRefPlab.Phi();
  fPiTkQuality = highPurity(p2);
  fPiPtNrf     = p2->fPlab.Perp();
  fPiEtaNrf    = p2->fPlab.Eta();

  fKaMissid = tightMuon(p1);  // true for tight  muons
  fPiMissid = tightMuon(p2);
  fKaMuMatch = doTriggerMatching(p1); // see if it matches HLT muon
  fPiMuMatch = doTriggerMatching(p2); // see if it matches HLT muon

  if(0) { // special tests d.k.
    double mva=0;

    fKaMissid2 = mvaMuon(p1,mva);  // true for tight  muons
    fKaMuMatch = doTriggerMatching(p1,false); // see if it matches HLT muon
    fKaMuMatch2 = doTriggerMatching(p1,true); // see if it matches HLT muon
    fKaMuMatchR = doTriggerMatchingR(p1,false); // see if it matches HLT muon
    fKaMuMatchR2 = doTriggerMatchingR(p1,true); // see if it matches HLT muon
    fKaMuMatchR3 = matchToMuon(p1,true); // see if it matches HLT muon

    // Pion
    fPiMissid2 = mvaMuon(p2,mva);  // true for tight  muons
    fPiMuMatch = doTriggerMatching(p2,false); // see if it matches HLT muon
    fPiMuMatch2 = doTriggerMatching(p2,true); // see if it matches HLT muon
    fPiMuMatchR = doTriggerMatchingR(p2,false); // see if it matches HLT muon
    fPiMuMatchR2 = doTriggerMatchingR(p2,true); // see if it matches HLT muon
    fPiMuMatchR3 = matchToMuon(p2,true); // see if it matches HLT muon

    //if(fKa1Missid || fKa2Missid) cout<<"missid "<<fKa1Missid<<" "<<fKa2Missid<<" "<<fKa1MuMatch<<" "<<fKa2MuMatch<<endl;
  }

  if (fCandTmi > -1 && fCandTmi == fpCand->fIndex) {
    if (fpEvt->getSimpleTrack(p1->fIndex)->getGenIndex() < fpEvt->nGenT()) {
      TGenCand *pg1 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(p1->fIndex)->getGenIndex());
      fKaPtGen     = pg1->fP.Perp();
      fKaEtaGen    = pg1->fP.Eta();
    }
    if (fpEvt->getSimpleTrack(p2->fIndex)->getGenIndex() < fpEvt->nGenT()) {
      TGenCand *pg2 = fpEvt->getGenTWithIndex(fpEvt->getSimpleTrack(p2->fIndex)->getGenIndex());
      pg2->dump();
      fPiPtGen     = pg2->fP.Perp();
      fPiEtaGen    = pg2->fP.Eta();
    }
  } else {
    fKaPtGen     = -99.;
    fKaEtaGen    = -99.;
    fPiPtGen     = -99.;
    fPiEtaGen    = -99.;
  }

  fDeltaR  = p1->fPlab.DeltaR(p2->fPlab);

  TLorentzVector ka, pi;
  ka.SetPtEtaPhiM(fKaPt, fKaEta, fKaPhi, MKAON);
  pi.SetPtEtaPhiM(fPiPt, fPiEta, fPiPhi, MPION);

  TLorentzVector kstarCand = ka + pi;
  fKstarPt   = kstarCand.Pt();
  fKstarEta  = kstarCand.Eta();
  fKstarPhi  = kstarCand.Phi();

  fGoodDeltaR = (fDeltaR < DELTAR);
  fGoodMKPI    = ((MKPILO < fMKPI ) && (fMKPI < MKPIHI));

  candAna::candAnalysis();
  fGoodTracksPt = fGoodTracksPt && ((fKaPt > TRACKPTLO) && (fPiPt > TRACKPTLO));

  fPreselection = fPreselection && fGoodJpsiMass && fGoodMKPI && fGoodDeltaR &&fGoodTracksPt;
  fPreselection = fPreselection && fWideMass;

  // -- overwrite specific variables
  fCandChi2    = chi2;
  fCandDof     = ndof;
  fCandChi2Dof = chi2/ndof;


  if(0) { // misid test d.k.
    if( (p1->fIndex == fpMuon1->fIndex) || (p1->fIndex ==fpMuon2->fIndex) )
      cout<<" Kaon is a MUON "<<fEvt<<" "<<fpCand<<" "<<p1->fIndex<<" "<<fpMuon1->fIndex<<" "<<fpMuon2->fIndex<<" "<<fEvt<<endl;
    if( (p2->fIndex == fpMuon1->fIndex) || (p2->fIndex ==fpMuon2->fIndex) )
      cout<<" Pion is a MUON "<<fEvt<<" "<<fpCand<<" "<<p2->fIndex<<" "<<fpMuon1->fIndex<<" "<<fpMuon2->fIndex<<" "<<fEvt<<endl;

    TVector3 muonMom1 = fpMuon1->fPlab;
    TVector3 muonMom2 = fpMuon2->fPlab;

    TVector3 trackMom = p1->fPlab;  // test track momentum
    double dR1 = muonMom1.DeltaR(trackMom);
    double dR2 = muonMom2.DeltaR(trackMom);

    if(dR1<dR2) { fKaMuMatchR4 = dR1; fKaMuMatchR5 = dR2;}
    else        { fKaMuMatchR4 = dR2; fKaMuMatchR5 = dR1;}

    trackMom = p2->fPlab;  // test track momentum
    dR1 = muonMom1.DeltaR(trackMom);
    dR2 = muonMom2.DeltaR(trackMom);

    if(dR1<dR2) { fPiMuMatchR4 = dR1; fPiMuMatchR5 = dR2;}
    else        { fPiMuMatchR4 = dR2; fPiMuMatchR5 = dR1;}
  } // end special test

  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(10);
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(4);
}

// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::moreBasicCuts() {
  cout << "   candAnaBd2JpsiKstar: more basic cuts" << endl;
  fAnaCuts.addCut("fGoodJpsiMass", "m(J/psi)", fGoodJpsiMass);
  fAnaCuts.addCut("fGoodDeltaR", "Delta R(KPI)", fGoodDeltaR);
  fAnaCuts.addCut("fGoodMKPI", "m(KPI) [GeV]", fGoodMKPI);
}

// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::genMatch() {

  fGenM1Tmi = fGenM2Tmi = fGenKTmi = -1;
  fNGenPhotons = 0;

  TGenCand *pC(0), *pB(0), *pPsi(0), *pKstar(0), *pM1(0), *pM2(0), *pK(0), *pPi(0), *pTmp(0);
  int ngamma(0);
  bool goodMatch(false);
  for (int i = 0; i < fpEvt->nGenT(); ++i) {
    pC = fpEvt->getGenT(i);
    if (511 == TMath::Abs(pC->fID)) {
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
	} else if (313 == TMath::Abs(pC->fID)) {
	  pKstar = fpEvt->getGenTWithIndex(id);
	  pK = pPi = 0;
	  for (int idd = pKstar->fDau1; idd <= pKstar->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd);
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (321 == TMath::Abs(pC->fID)) {
	      pK = fpEvt->getGenTWithIndex(idd);
	    }
	    if (211 == TMath::Abs(pC->fID)) {
	      pPi = fpEvt->getGenTWithIndex(idd);
	    }
	  }
	}
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK && 0 != pPi && pPsi != 0 && pKstar != 0 && (pPsi->fMom1 == pKstar->fMom1)) {
	// -- check that there are no other direct daughters than J/psi Kstar (plus possibly photons)
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
  fKaGenID = -99999;
  fPiGenID = -99999;
  if (goodMatch) {
    fMu1GenID = pM1->fID;
    fMu2GenID = pM2->fID;
    fKaGenID = pK->fID;
    fPiGenID = pPi->fID;
    fGenBTmi = pB->fNumber;
    double m = pB->fP.Mag();
    double p = pB->fP.P();
    // mother pointer
    TGenCand *pM = fpEvt->getGenTWithIndex(pB->fMom1);
    // use the mother if it has the same PDGID (it oscillated)
    if (TMath::Abs(pB->fID) != TMath::Abs(pM->fID)) pM = pB;
    double x = (pM1->fV - pM->fV).Mag();
    fGenLifeTime = x*m/p/TMath::Ccgs();
    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fGenM1Tmi = pM1->fNumber;
      fGenM2Tmi = pM2->fNumber;
    } else {
      fGenM1Tmi = pM2->fNumber;
      fGenM2Tmi = pM1->fNumber;
    }
    fGenKTmi = pK->fNumber;
    fGenPiTmi = pPi->fNumber;
  } else {
    fGenM1Tmi = -1;
    fGenM2Tmi = -1;
    fGenKTmi = -1;
    fGenPiTmi = -1;
  }
}


// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::genMatchOld() {

  fGenM1Tmi = fGenM2Tmi = fGenKTmi = -1;
  fNGenPhotons = 0;

  TGenCand *pC(0), *pB(0), *pPsi(0), *pKstar(0), *pM1(0), *pM2(0), *pK(0), *pPi(0);
  int nb(0), ngamma(0);
  bool goodMatch(false);
  for (int i = 0; i < fpEvt->nGenT(); ++i) {
    pC = fpEvt->getGenT(i);
    if (511 == TMath::Abs(pC->fID)) {
      pB = pC;
      nb = pB->fDau2 - pB->fDau1 + 1;
      if (nb > 2) continue; // skip B decays where more than J/psi and kstar came from B
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
	} else if (313 == TMath::Abs(pC->fID)) {
	  pKstar = fpEvt->getGenTWithIndex(id);
	  pK = pPi = 0;
	  for (int idd = pKstar->fDau1; idd <= pKstar->fDau2; ++idd) {
	    pC = fpEvt->getGenTWithIndex(idd);
	    if (22 == TMath::Abs(pC->fID)) ++ngamma;
	    if (321 == TMath::Abs(pC->fID)) {
	      if (0 == pK) {
		pK = fpEvt->getGenTWithIndex(idd);
	      } else {
		pPi = fpEvt->getGenTWithIndex(idd);
	      }
	    }
	  }
	}
      }
      if (0 != pM1 && 0 != pM2 && 0 != pK && 0 != pPi
	  && (pPsi->fMom1 == pKstar->fMom1)
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
  fKaGenID = -99999;
  fPiGenID = -99999;
  if (goodMatch) {
    fMu1GenID = pM1->fID;
    fMu2GenID = pM2->fID;
    fKaGenID = pK->fID;
    fPiGenID = pPi->fID;
    fGenBTmi = pB->fNumber;
    double m = pB->fP.Mag();
    double p = pB->fP.P();
    // Meson pointer
    TGenCand *pM = fpEvt->getGenTWithIndex(pB->fMom1);
    // the meson is the original except if it oscillated
    if (511 != TMath::Abs(pM->fID)) pM = pB;
    double x = (pM1->fV - pM->fV).Mag();
    fGenLifeTime = x*m/p/TMath::Ccgs();
    if (pM1->fP.Perp() > pM2->fP.Perp()) {
      fGenM1Tmi = pM1->fNumber;
      fGenM2Tmi = pM2->fNumber;
    } else {
      fGenM1Tmi = pM2->fNumber;
      fGenM2Tmi = pM1->fNumber;
    }
    // the code below is wrong!
    if (pK->fP.Perp() > pPi->fP.Perp()) {
      fGenKTmi = pK->fNumber;
      fGenPiTmi = pPi->fNumber;
    } else {
      fGenKTmi = pPi->fNumber;
      fGenPiTmi = pK->fNumber;
    }
  } else {
    fGenM1Tmi = -1;
    fGenM2Tmi = -1;
    fGenKTmi = -1;
    fGenPiTmi = -1;
  }

}


// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::recoMatch() {

  fRecM1Tmi = fRecM2Tmi = fRecKTmi = fRecPiTmi =-1;
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

    // -- kaon
    if (fGenKTmi > -1 && pT->getGenIndex() == fGenKTmi) {
      fRecKTmi = i;
    }

    // -- pion
    if (fGenPiTmi > -1 && pT->getGenIndex() == fGenPiTmi) {
      fRecPiTmi = i;
    }

    // -- skip rest if all matches found
    if (fRecM1Tmi > -1 && fRecM2Tmi > -1 && fRecKTmi > -1 && fRecPiTmi > -1) break;
  }


  if (fVerbose > 10) {
    cout << "fRecM1Tmi = " << fRecM1Tmi << " matched to fGenM1Tmi = " << fGenM1Tmi << endl;
    cout << "fRecM2Tmi = " << fRecM2Tmi << " matched to fGenM2Tmi = " << fGenM2Tmi << endl;
    cout << "fRecKTmi = " << fRecKTmi << " matched to fGenKTmi = " << fGenKTmi << endl;
    cout << "fRecPiTmi = " << fRecPiTmi << " matched to fGenPiTmi = " << fGenPiTmi << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::candMatch() {

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
	cout << idx << " " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecKTmi << endl;
      }
      if (fRecM1Tmi > -1 && type == 13 && idx == fRecM1Tmi) {
	d1Matched = 1;
      }
      if (fRecM2Tmi > -1 && type == 13 && idx == fRecM2Tmi) {
	d2Matched = 1;
      }
      if (fRecKTmi > -1 && type == 321 && idx == fRecKTmi) {
	d3Matched = 1;
      }
      if (fRecPiTmi > -1 && type == 211 && idx == fRecPiTmi) {
	d4Matched = 1;
      }
    }

    if (d1Matched && d2Matched && d3Matched && d4Matched) {
      fCandTmi = iC;
      break;
    }
  }
  if (fVerbose > 10) {
    cout << "fCandTmi = " << fCandTmi << " matched to rec tracks " << fRecM1Tmi << " " << fRecM2Tmi << " " << fRecKTmi  << " " << fRecPiTmi
	 << endl;
  }

}


// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::bookHist() {
  cout << "==>candAnaBd2JpsiKstar: bookHist" << endl;
  candAna::bookHist();

  moreReducedTree(fTree);

  // -- Additional effTree variables
  fEffTree->Branch("kpt",   &fETkpt,            "kpt/F");
  fEffTree->Branch("g3pt",  &fETg3pt,           "g3pt/F");
  fEffTree->Branch("keta",  &fETketa,           "keta/F");
  fEffTree->Branch("g3eta", &fETg3eta,          "g3eta/F");
  fEffTree->Branch("kq",    &fETkq,             "kq/I");
  fEffTree->Branch("kgt",   &fETkgt,            "kgt/O");

  fEffTree->Branch("pipt",   &fETpipt,          "pipt/F");
  fEffTree->Branch("g4pt",   &fETg4pt,          "g4pt/F");
  fEffTree->Branch("pieta",  &fETpieta,         "pieta/F");
  fEffTree->Branch("g4eta",  &fETg4eta,         "g4eta/F");
  fEffTree->Branch("piq",    &fETpiq,           "piq/I");
  fEffTree->Branch("pigt",   &fETpigt,          "pigt/O");

}


// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::moreReducedTree(TTree *t) {

  // -- Additional reduced tree variables
  t->Branch("mpsi",        &fJpsiMass,    "mpsi/D");
  t->Branch("psipt",       &fJpsiPt,      "psipt/D");
  t->Branch("psieta",      &fJpsiEta,     "psieta/D");
  t->Branch("psiphi",      &fJpsiPhi,     "psiphi/D");
  t->Branch("psicosa",     &fJpsiCosA,    "psicosa/D");
  t->Branch("psimaxdoca",  &fJpsiMaxDoca, "psimaxdoca/D");
  t->Branch("psiflsxy",    &fJpsiFLSxy,   "psiflsxy/D");
  t->Branch("psiprob",     &fJpsiVtxProb, "psiprob/D");

  t->Branch("mkpi",      &fMKPI,      "mkpi/D");
  t->Branch("kstarpt",   &fKstarPt,    "kstarpt/D");
  t->Branch("kstareta",  &fKstarEta,   "kstareta/D");
  t->Branch("kstarphi",  &fKstarPhi,   "kstarphi/D");
  t->Branch("dr",        &fDeltaR,   "dr/D");
  t->Branch("kstarfail", &fKstarFail,   "kstarfail/O");

  t->Branch("kpt",  &fKaPt,    "kpt/D");
  t->Branch("keta", &fKaEta,   "keta/D");
  t->Branch("kphi", &fKaPhi,   "kphi/D");
  t->Branch("kgt",  &fKaTkQuality,"kgt/I");
  t->Branch("pipt",  &fPiPt,    "pipt/D");
  t->Branch("pieta", &fPiEta,   "pieta/D");
  t->Branch("piphi", &fPiPhi,   "piphi/D");
  t->Branch("pigt",  &fPiTkQuality,"pigt/I");



  t->Branch("t3pt",  &fKaPtNrf, "t3pt/D");
  t->Branch("t3eta", &fKaEtaNrf,"t3eta/D");

  t->Branch("t4pt",  &fPiPtNrf, "t4pt/D");
  t->Branch("t4eta", &fPiEtaNrf,"t4eta/D");

  t->Branch("g3pt", &fKaPtGen,  "g3pt/D");
  t->Branch("g3eta",&fKaEtaGen, "g3eta/D");
  t->Branch("g3id", &fKaGenID,  "g3id/I");

  t->Branch("g4pt", &fPiPtGen,  "g4pt/D");
  t->Branch("g4eta",&fPiEtaGen, "g4eta/D");
  t->Branch("g4id", &fPiGenID,  "g4id/I");

  t->Branch("kmissid",  &fKaMissid,    "kmissid/O");
  t->Branch("pimissid",  &fPiMissid,    "pimissid/O");
  t->Branch("kmumatch", &fKaMuMatch,  "kmumatch/O");
  t->Branch("pimumatch", &fPiMuMatch,  "pimumatch/O");

  if(0) { // for testing d.k.

    t->Branch("kmissid2",  &fKaMissid2,    "kmissid2/O");
    t->Branch("kmumatch2", &fKaMuMatch2,    "kmumatch2/O");

    t->Branch("kmumatchr", &fKaMuMatchR,    "kmumatchr/F");
    t->Branch("kmumatchr2", &fKaMuMatchR2,    "kmumatchr2/F");
    t->Branch("kmumatchr3", &fKaMuMatchR3,    "kmumatchr3/F");
    t->Branch("kmumatchr4", &fKaMuMatchR4,    "kmumatchr4/F");
    t->Branch("kmumatchr5", &fKaMuMatchR5,    "kmumatchr5/F");

    t->Branch("pimissid2",  &fPiMissid2,    "pimissid2/O");
    t->Branch("pimumatch2", &fPiMuMatch2,   "pimumatch2/O");

    t->Branch("pimumatchr",  &fPiMuMatchR,   "pimumatchr/F");
    t->Branch("pimumatchr2", &fPiMuMatchR2,  "pimumatchr2/F");
    t->Branch("pimumatchr3", &fPiMuMatchR3,  "pimumatchr3/F");
    t->Branch("pimumatchr4", &fPiMuMatchR4,  "pimumatchr4/F");
    t->Branch("pimumatchr5", &fPiMuMatchR5,  "pimumatchr5/F");
  } // end testing

}


// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::fillCandidateHistograms(int offset) {
  candAna::fillCandidateHistograms(offset);
}


// ----------------------------------------------------------------------
void candAnaBd2JpsiKstar::efficiencyCalculation() {

  fGoodEffCand = false;

  // -- gen level
  TGenCand *pB(0), *pM1(0), *pM2(0), *pK(0), *pPi(0);
  if (-1 == fGenM1Tmi || -1 == fGenM2Tmi || -1 == fGenKTmi || -1 == fGenPiTmi) {
    if (fVerbose > 2 ) cout << "--------------------> No matched signal decay found" << endl;
    return;
  }
  pB  = fpEvt->getGenTWithIndex(fGenBTmi);
  pM1 = fpEvt->getGenTWithIndex(fGenM1Tmi);
  pM2 = fpEvt->getGenTWithIndex(fGenM2Tmi);
  pK = fpEvt->getGenTWithIndex(fGenKTmi);
  pPi = fpEvt->getGenTWithIndex(fGenPiTmi);

  // -- reco level
  TSimpleTrack *prM1(0), *prM2(0), *prK(0), *prPi(0);
  double bla(0);
  int m1Matched(0), m2Matched(0), kMatched(0), piMatched(0), m1ID(0), m1tmID(0), m1mvaID(0), m2ID(0), m2tmID(0), m2mvaID(0),
    m1GT(0), m2GT(0), kGT(0), piGT(0);
  if (fRecM1Tmi > -1) {
    m1Matched = 1;
    prM1 = fpEvt->getSimpleTrack(fRecM1Tmi);
    if (tightMuon(prM1)) m1tmID = 1;
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
    if (tightMuon(prM2)) m2tmID = 1;
    if (mvaMuon(prM2, bla)) m2mvaID = 1;
    if (prM2->getHighPurity()) {
      m2GT = 1;
    } else {
      m2GT = 0;
    }
  }

  if (fRecKTmi > -1) {
    kMatched = 1;
    prK = fpEvt->getSimpleTrack(fRecKTmi);
    if (prK->getHighPurity()) {
      kGT = 1;
    } else {
      kGT = 0;
    }
  }

  if (fRecPiTmi > -1) {
    piMatched = 1;
    prPi = fpEvt->getSimpleTrack(fRecPiTmi);
    if (prPi->getHighPurity()) {
      piGT = 1;
    } else {
      piGT = 0;
    }
  }

  // -- cand level
  TAnaCand *pCand(0);
  if (fCandTmi > -1) {
    pCand = fpEvt->getCand(fCandTmi);
  }

  m1ID = m1tmID;
  m2ID = m2tmID;

  // -- EffTree filling for all events with a signal decay
  fETgm    = pB->fP.M();
  fETgpt   = pB->fP.Perp();
  fETgtau  = fGenLifeTime;
  fETgeta  = pB->fP.Eta();
  fETg1pt  = pM1->fP.Perp();
  fETg1eta = pM1->fP.Eta();
  fETg2pt  = pM2->fP.Perp();
  fETg2eta = pM2->fP.Eta();
  fETg3pt  = pK->fP.Perp();
  fETg3eta = pK->fP.Eta();
  fETg4pt  = pPi->fP.Perp();
  fETg4eta = pPi->fP.Eta();
  if (m1Matched) {
    fETm1pt  = prM1->getP().Perp();
    fETm1eta = prM1->getP().Eta();
    fETm1q   = prM1->getCharge();
    fETm1gt  = (m1GT>0?true:false);
    fETm1id  = (m1ID>0?true:false);
    fETm1tmid  = (m1tmID>0?true:false);
    fETm1mvaid = (m1mvaID>0?true:false);
  } else {
    fETm1pt  = -99.;
    fETm1eta = -99.;
    fETm1q   = -99;
    fETm1gt  = false;
    fETm1id  = false;
    fETm1tmid  = false;
    fETm1mvaid = false;
  }
  if (m2Matched) {
    fETm2pt  = prM2->getP().Perp();
    fETm2eta = prM2->getP().Eta();
    fETm2q   = prM2->getCharge();
    fETm2gt  = (m2GT>0?true:false);
    fETm2id  = (m2ID>0?true:false);
    fETm2tmid  = (m2tmID>0?true:false);
    fETm2mvaid = (m2mvaID>0?true:false);
  } else {
    fETm2pt  = -99.;
    fETm2eta = -99.;
    fETm2q   = -99;
    fETm2gt  = false;
    fETm2id  = false;
    fETm2tmid  = false;
    fETm2mvaid = false;
  }
  if (m1Matched && m2Matched) {
    fETchan = detChan(fETm1eta, fETm2eta);
  } else {
    fETchan = -1;
  }
  if (kMatched) {
    fETkpt  = prK->getP().Perp();
    fETketa = prK->getP().Eta();
    fETkq   = prK->getCharge();
    fETkgt  = (kGT>0?true:false);
  } else {
    fETkpt  = -99.;
    fETketa = -99.;
    fETkq   = -99;
    fETkgt  = false;
  }
  if (piMatched) {
    fETpipt  = prPi->getP().Perp();
    fETpieta = prPi->getP().Eta();
    fETpiq   = prPi->getCharge();
    fETpigt  = (piGT>0?true:false);
  } else {
    fETpipt  = -99.;
    fETpieta = -99.;
    fETpiq   = -99;
    fETpigt  = false;
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
void candAnaBd2JpsiKstar::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump);

  fCutFile = filename;

  if (dump) cout << "==> candAnaBd2JpsiKstar: Reading " << fCutFile << " for cut settings" << endl;
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

    if (!strcmp(CutName, "KSTARTYPE")) {
      KSTARTYPE = static_cast<int>(CutValue);
      if (dump) cout << "KSTARTYPE:      " << KSTARTYPE << endl;
      ibin = 213;
      hcuts->SetBinContent(ibin, KSTARTYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #ktar ID :: %d", CutName, KSTARTYPE));
    }

    if (!strcmp(CutName, "MKPILO")) {
      MKPILO = CutValue;
      if (dump) cout << "MKPILO:           " << MKPILO << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKPILO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{min}(KPI) :: %3.1f", CutName, MKPILO));
    }

    if (!strcmp(CutName, "MKPIHI")) {
      MKPIHI = CutValue;
      if (dump) cout << "MKPIHI:           " << MKPIHI << endl;
      ibin = 300;
      hcuts->SetBinContent(ibin, MKPIHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: m^{max}(KPI) :: %3.1f", CutName, MKPIHI));
    }

    if (!strcmp(CutName, "DELTAR")) {
      DELTAR = CutValue;
      if (dump) cout << "DELTAR:           " << DELTAR << endl;
      ibin = 301;
      hcuts->SetBinContent(ibin, DELTAR);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #Delta R(KPI) :: %3.1f", CutName, DELTAR));
    }



  }

}
