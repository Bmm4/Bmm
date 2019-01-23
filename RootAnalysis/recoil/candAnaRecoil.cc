#include "candAnaRecoil.hh"
#include <cmath>
#include <string>
#include <stdarg.h>

#include "common/HFMasses.hh"
#include "common/util.hh"
#include "common/ana.hh"

using namespace std;

string streamTVector3(TVector3 v) {
  return string(Form("(%+4.3f, %+4.3f, %+4.3f)", v.X(),  v.Y(),  v.Z()));
}

// ----------------------------------------------------------------------
candAnaRecoil::candAnaRecoil(recoilReader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaRecoil: constructor..." << endl;
  readCuts(cutsFile, 1);

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
void candAnaRecoil::genAnalysis() {
  TGenCand *pCand, *pD, *pT;
  // -- Force initialization!
  fpGenB = fpGenMu = fpGenKa = fpGenTau = fpGenHad1 = fpGenHad2 = fpGenHad3 = fpGenNu = fpGenGa1 = fpGenGa2 = 0;
  fGenVtxBProd = fGenVtxBDecay = fGenVtxTauDecay = fGenDirectionTau = TVector3(0., 0., 0.);
  f4GenTauAtDecay = f4GenNur0Pos = f4GenNur0Neg = f4GenBr0Pos = f4GenBr0Neg = f4GenBr0 = f4GenNur0 = TLorentzVector(0., 0., 0., 0.);

  fGenBPt =
    fGenBr0Mass = fGenBr0PosMass = fGenBr0NegMass =
    fGenTauPt = fGenTauDecTime =
    fGenTauHadPt = fGenTauHadMass = fGenTauHadAngle = fGenTauHadPerp = fGenTauHadPara =
    fGen2HadMinMass = fGen2HadMaxMass =
    fGenHad1Pt = fGenHad2Pt = fGenHad3Pt =
    fGenMuPt = fGenKaPt = fGenMuKaMass =
    fGenNur0PosPara =
    fGenNur0NegPara =
    fGenNuPara = fGenNuPerp
    = -99.;

  fNgamma = 0;

  // -- now try to find signal decay
  bool goodMatch(false);
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if (521 == TMath::Abs(pCand->fID)) {
      fpGenB = fpGenMu = fpGenKa = fpGenTau = fpGenHad1 = fpGenHad2 = fpGenHad3 = fpGenNu = fpGenGa1 = fpGenGa2 = 0;
      fGenVtxBProd = fGenVtxBDecay = fGenVtxTauDecay = fGenDirectionTau = TVector3(0., 0., 0.);
      fpGenB = pCand;
      fGenVtxBProd = pCand->fV;
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = fpEvt->getGenCand(iD);
	if (321 == TMath::Abs(pD->fID)) fpGenKa = pD;
	if (13 == TMath::Abs(pD->fID)) fpGenMu = pD;
	if (15 == TMath::Abs(pD->fID)) {
	  fpGenTau = pD;
	  fGenVtxBDecay = pD->fV; // pD->fV is the production vertex for pD, i.e. the decay vertex of its mother
	  for (int iT = pD->fDau1; iT <= pD->fDau2; ++iT) {
	    pT = fpEvt->getGenCand(iT);
	    if (16 == TMath::Abs(pT->fID)) fpGenNu = pT;
	    if (0 != pT->fQ) {
	      if (0 == fpGenHad1) {
		fpGenHad1 = pT;
		fGenVtxTauDecay = pT->fV;
	      } else if (0 == fpGenHad2) {
		fpGenHad2 = pT;
	      } else if (0 == fpGenHad3) {
		fpGenHad3 = pT;
	      } else {
		cout << "XXXXXXXXXXXXXXXX what should I do here?????" << endl;
	      }
	    }
	    if (22 == TMath::Abs(pT->fID)) {
	      if (0 == fpGenGa1) {
		fpGenGa1 = pT;
	      } else if (0 == fpGenGa2) {
		fpGenGa2 = pT;
	      } else {
		cout << "XXXXXXXXXXXXXXXXXXXXXXXX photon overflow!!!!!!!!" << endl;
	      }
	    }
	  }
	}

      }
      if (fpGenB && fpGenMu && fpGenKa && fpGenTau && fpGenHad1 && fpGenHad2 && fpGenHad3) {
	goodMatch = true;
	break;
      }
    }
  }

  if (goodMatch) {
    // -- sort according to pT of hadronic tracks
    double pt1 = fpGenHad1->fP.Perp();
    double pt2 = fpGenHad2->fP.Perp();
    double pt3 = fpGenHad3->fP.Perp();
    TGenCand *tgc1(0), *tgc2(0), *tgc3(0);
    if (pt1 > pt2 && pt1 > pt3) {
      tgc1 = fpGenHad1;
      if (pt2 > pt3) {
	tgc2 = fpGenHad2;
	tgc3 = fpGenHad3;
      } else {
	tgc2 = fpGenHad3;
	tgc3 = fpGenHad2;
      }
    } else if (pt2 > pt1 && pt2 > pt3) {
      tgc1 = fpGenHad2;
      if (pt1 > pt3) {
	tgc2 = fpGenHad1;
	tgc3 = fpGenHad3;
      } else {
	tgc2 = fpGenHad3;
	tgc3 = fpGenHad1;
      }
    } else if (pt3 > pt1 && pt3 > pt2) {
      tgc1 = fpGenHad3;
      if (pt1 > pt2) {
	tgc2 = fpGenHad1;
	tgc3 = fpGenHad2;
      } else {
	tgc2 = fpGenHad2;
	tgc3 = fpGenHad1;
      }
    } else {
      cout << "THIS SHOULD NOT HAPPEN!!!" << endl;
    }
    fpGenHad1 = tgc1;
    fpGenHad2 = tgc2;
    fpGenHad3 = tgc3;

    f4GenTauHad   = fpGenHad1->fP + fpGenHad2->fP + fpGenHad3->fP;
    fGenDirectionTau = fGenVtxTauDecay - fGenVtxBDecay;
    fGenTauDecTime = fGenDirectionTau.Mag()*MTAU/fpGenTau->fP.Vect().Mag()/TMath::Ccgs();
    fGenTauHadAngle = fGenDirectionTau.Angle(f4GenTauHad.Vect());

    f4GenTauAtDecay = fpGenHad1->fP + fpGenHad2->fP + fpGenHad3->fP  + fpGenNu->fP;
    if (fpGenGa1) {
      f4GenTauAtDecay += fpGenGa1->fP;
      ++fNgamma;
    }
    if (fpGenGa2) {
      f4GenTauAtDecay += fpGenGa2->fP;
      ++fNgamma;
    }

    pair<TVector3, TVector3> nucomp = parallelAndPerp(fGenDirectionTau, fpGenNu->fP.Vect());
    fGenNuPara = nucomp.first.Mag();
    fGenNuPerp = nucomp.second.Mag();

    pair<TVector3, TVector3> decomp = parallelAndPerp(fGenDirectionTau, f4GenTauHad.Vect());
    double evis  = f4GenTauHad.E();
    double mvis  = f4GenTauHad.M();
    double ppar  = decomp.first.Mag();
    double pperp = decomp.second.Mag();
    pair<double, double> twosol = nuRecoMom0(ppar, pperp, evis, mvis, MTAU);
    fGenTauHadPerp = ppar;
    fGenTauHadPara = pperp;

    // -- calculate reconstructed neutrino 4-momentum
    TVector3 vDirUnit = fGenDirectionTau.Unit();
    TVector3 vRecoPos = twosol.first*vDirUnit - decomp.second;
    TVector3 vRecoNeg = twosol.second*vDirUnit - decomp.second;
    f4GenNur0Pos.SetVectM(vRecoPos, 0.);
    f4GenNur0Neg.SetVectM(vRecoNeg, 0.);

    fGenNur0PosPara = twosol.first;
    fGenNur0NegPara = twosol.second;


    f4GenBr0Pos = f4GenNur0Pos + f4GenTauHad + fpGenMu->fP + fpGenKa->fP;
    f4GenBr0Neg = f4GenNur0Neg + f4GenTauHad + fpGenMu->fP + fpGenKa->fP;

    fGenBPt        = fpGenB->fP.Perp();
    fGenTauHadMass = f4GenTauHad.M();

    fGenTaur0PosPt = (f4GenNur0Pos + f4GenTauHad).Perp();
    fGenTaur0NegPt = (f4GenNur0Neg + f4GenTauHad).Perp();

    fGenBr0PosMass = f4GenBr0Pos.M();
    fGenBr0NegMass = f4GenBr0Neg.M();
    TVector3 vTrue = fpGenNu->fP.Vect();
    double diffPos = (vRecoPos - vTrue).Mag();
    double diffNeg = (vRecoNeg - vTrue).Mag();
    if (diffPos < diffNeg) {
      fGenTaur0Pt  = fGenTaur0PosPt;
      fGenBr0Mass  = fGenBr0PosMass;
    } else {
      fGenTaur0Pt  = fGenTaur0NegPt;
      fGenBr0Mass  = fGenBr0NegMass;
    }

    fGenTauPt      = fpGenTau->fP.Perp();
    fGenTauHadPt   = f4GenTauHad.Perp();

    fGenHad1Pt     = fpGenHad1->fP.Perp();
    fGenHad2Pt     = fpGenHad2->fP.Perp();
    fGenHad3Pt     = fpGenHad3->fP.Perp();

    vector<TLorentzVector> vhad;
    vhad.push_back(fpGenHad1->fP);
    vhad.push_back(fpGenHad2->fP);
    vhad.push_back(fpGenHad3->fP);
    fGen2HadMinMass = minMassPair(vhad);
    fGen2HadMaxMass = maxMassPair(vhad);

    fGenMuPt       = fpGenMu->fP.Perp();
    fGenKaPt       = fpGenKa->fP.Perp();
    fGenMuKaMass   = (fpGenMu->fP + fpGenKa->fP).M();



  }

  ((TH1D*)fHistDir->Get("ngamma"))->Fill(fNgamma);

}


// ----------------------------------------------------------------------
pair<TVector3,TVector3> candAnaRecoil::parallelAndPerp(TVector3 direction, TVector3 momVis) {
  TVector3 par, perp;

  TVector3 unitDir = direction.Unit();
  double compPar = momVis.Dot(unitDir);
  par = compPar*unitDir;

  perp = momVis - par;
  return make_pair(par, perp);
}


// ----------------------------------------------------------------------
pair<double, double> candAnaRecoil::nuRecoMom0(double compVisPar, double compVisPerp, double eVis, double mVis, double mTot) {
  double mTot2  = mTot*mTot;
  double eVis2  = eVis*eVis;
  double mVis2  = mVis*mVis;
  double pperp2 = compVisPerp * compVisPerp;
  double ppar2  = compVisPar * compVisPar;

  double paren  = (mTot2 - mVis2 - 2.*pperp2);

  double a = (paren*compVisPar)/(2.*(ppar2 - eVis2));

  double r = (paren*paren*eVis2)/(4.*(ppar2 - eVis2)*(ppar2 - eVis2)) + (eVis2*pperp2)/(ppar2 - eVis2);

  double solPlus(0.), solMinus(0.);
  if (r > 0) {
    solPlus  = -a + TMath::Sqrt(r);
    solMinus = -a - TMath::Sqrt(r);
  }

  return make_pair(solPlus, solMinus);
}


// ----------------------------------------------------------------------
void candAnaRecoil::candAnalysis() {
  //  dumpBkmt();
  cout << "----------------------------------------------------------------------" << endl;
  cout << "candAnaRecoil::candAnalysis() event =  " << fChainEvent << endl;
  brecoAnalysis();

  TAnaTrack *pS(0);
  // -- TRUTH J/psi K*0 candidate
  vector<int> brecoIdx;
  int nmatch(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    if (3000069 == pCand->fType) {
      cout << "===> J/psi Kstar0: Truth" << endl;
      pCand->dump();
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	pS = fpEvt->getSigTrack(is);
	cout << pS->fIndex << " ";
	brecoIdx.push_back(pS->fIndex);
      }
      cout << endl;
    }
  }
  // -- J/psi X candidate
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    if (10000 == pCand->fType) {
      int found(0);
      cout << "===> J/psi X: Cand" << endl;
      pCand->dump();
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	pS = fpEvt->getSigTrack(is);
	if (brecoIdx.end() == find(brecoIdx.begin(), brecoIdx.end(), pS->fIndex)) {
	  cout << pS->fIndex << " ";
	} else {
	  cout << pS->fIndex << "* ";
	  ++nmatch;
	}
      }
      cout << endl;
      if (nmatch == brecoIdx.size()) cout << "** completely J/psi X matched" << endl;
    }
  }

  // -- TRUTH candidate B+ -> mu tau K
  vector<int> signalIdx;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    if (CANDTYPE == pCand->fType) {
      cout << "===> CANDTYPE: " << CANDTYPE << endl;
      pCand->dump();
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	pS = fpEvt->getSigTrack(is);
	cout << pS->fIndex <<  Form("(pv:%d/pt:%3.2f) ", pS->fPvIdx, pS->fPlab.Perp());
	signalIdx.push_back(pS->fIndex);
      }
      cout << endl;
    }
  }

  // -- recoil of J/psi X candidate
  nmatch = 0;
  vector<int> recoilIdx;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    if (10001 == pCand->fType) {
      cout << "===> J/psi X: RECOIL" << endl;
      pCand->dump();
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	pS = fpEvt->getSigTrack(is);
	if (signalIdx.end() == find(signalIdx.begin(), signalIdx.end(), pS->fIndex)) {
	  cout << pS->fIndex << " ";
	} else {
	  cout << pS->fIndex << "* ";
	  ++nmatch;
	}
	recoilIdx.push_back(pS->fIndex);
      }
      cout << endl;
      if (nmatch == signalIdx.size()) cout << "** completely signal matched" << endl;
    }
  }
  // -- missed recoil of J/psi X candidate
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    TAnaCand *pCand = fpEvt->getCand(iC);
    if (10002 == pCand->fType) {
      cout << "===> J/psi X: MISSED RECOIL" << endl;
      pCand->dump();
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	pS = fpEvt->getSigTrack(is);
	cout << pS->fIndex << Form("(pv:%d/pt:%3.2f/doca:%4.3f) ", pS->fPvIdx, pS->fPlab.Perp(), pS->fDouble1);
	if ((pS->fPlab.Perp() > 0.5)
	    && (recoilIdx.end() == find(recoilIdx.begin(), recoilIdx.end(), pS->fIndex))
	    && highPurity(pS)
	    ) {
	  ((TH1D*)fHistDir->Get("doca"))->Fill(pS->fDouble1);
	}
      }
      cout << endl;
    }
  }

  candAna::candAnalysis();
}


// ----------------------------------------------------------------------
void candAnaRecoil::moreReducedTree(TTree *t) {
  // -- gen level
  t->Branch("gngamma",     &fNgamma,         "gngamma/I");
  t->Branch("gptb",        &fGenBPt,         "gptb/D");
  t->Branch("gmtauhad",    &fGenTauHadMass,  "gmtauhad/D");
  t->Branch("gm2hadmax",   &fGen2HadMaxMass, "gm2hadmax/D");
  t->Branch("gm2hadmin",   &fGen2HadMinMass, "gm2hadmin/D");
  t->Branch("gmbr0pos",    &fGenBr0PosMass,  "gmbr0pos/D");
  t->Branch("gmbr0neg",    &fGenBr0NegMass,  "gmbr0neg/D");
  t->Branch("gmbr0",       &fGenBr0Mass,     "gmbr0/D");
  t->Branch("gt0tau",      &fGenTauDecTime,  "gt0tau/D");
  t->Branch("gpttau",      &fGenTauPt,       "gpttau/D");
  t->Branch("gpttauhad",   &fGenTauHadPt,    "gpttauhad/D");
  t->Branch("gatauhad",    &fGenTauHadAngle, "gatauhad/D");
  t->Branch("gparanur0pos",&fGenNur0PosPara, "gparanur0pos/D");
  t->Branch("gparanur0neg",&fGenNur0NegPara, "gparanur0neg/D");
  t->Branch("gparanu",     &fGenNuPara,      "gparanu/D");
  t->Branch("gperpnu",     &fGenNuPerp,      "gperpnu/D");
  t->Branch("gperptauhad", &fGenTauHadPerp,  "gperptauhad/D");
  t->Branch("gparatauhad", &fGenTauHadPara,  "gparatauhad/D");
  t->Branch("gpttaur0pos", &fGenTaur0PosPt,  "gpttaur0pos/D");
  t->Branch("gpttaur0neg", &fGenTaur0NegPt,  "gpttaur0neg/D");
  t->Branch("gpttaur0",    &fGenTaur0Pt,     "gpttaur0/D");
  t->Branch("gptmu",       &fGenMuPt,        "gptmu/D");
  t->Branch("gptka",       &fGenKaPt,        "gptka/D");
  t->Branch("gmmuka",      &fGenMuKaMass,    "gmmuka/D");
  t->Branch("gpthad1",     &fGenHad1Pt,      "gpthad1/D");
  t->Branch("gpthad2",     &fGenHad2Pt,      "gpthad2/D");
  t->Branch("gpthad3",     &fGenHad3Pt,      "gpthad3/D");

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
  new TH1D("ngamma", "ngamma", 5, 0., 5.);
  new TH1D("doca", "doca", 100, 0., 2.0);
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


// ----------------------------------------------------------------------
double candAnaRecoil::minMassPair(vector<TLorentzVector> v ) {
  double minMass(99999.);
  for (unsigned int i = 0; i < v.size(); ++i) {
    for (unsigned int j = i+1; j < v.size(); ++j) {
      double mass = (v[i] + v[j]).M();
      if (mass < minMass) minMass = mass;
    }
  }
  return minMass;
}

// ----------------------------------------------------------------------
double candAnaRecoil::maxMassPair(vector<TLorentzVector> v ) {
  double maxMass(-99999.);
  for (unsigned int i = 0; i < v.size(); ++i) {
    for (unsigned int j = i+1; j < v.size(); ++j) {
      double mass = (v[i] + v[j]).M();
      if (mass > maxMass) maxMass = mass;
    }
  }
  return maxMass;
}
