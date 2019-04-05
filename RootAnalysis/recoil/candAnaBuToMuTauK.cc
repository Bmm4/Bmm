#include "candAnaBuToMuTauK.hh"
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
candAnaBuToMuTauK::candAnaBuToMuTauK(recoilReader *pReader, std::string name, std::string cutsFile) : candAna(pReader, name, cutsFile) {
  cout << "==> candAnaBuToMuTauK: constructor..." << endl;
  readCuts(cutsFile, 1);

}


// ----------------------------------------------------------------------
candAnaBuToMuTauK::~candAnaBuToMuTauK() {
  cout << "==> candAnaBuToMuTauK: destructor..." << endl;
}


// ----------------------------------------------------------------------
void candAnaBuToMuTauK::resetAllData() {
  fSignalTracks.clear();

  // -- pointers to SIGNAL HepMC cands and vertices
  fpGenMu=  fpGenKa= fpGenTau= fpGenHad1= fpGenHad2= fpGenHad3= fpGenNu= fpGenGa1= fpGenGa2 = 0;
  fGenVtxBProd= fGenVtxBDecay= fGenVtxTauDecay= fGenDirectionTau= TVector3(0., 0., 0.);
  // -- calculated/derived quantities
  f4GenTauAtDecay= f4GenTauHad= f4GenNur0= f4GenNur0Pos= f4GenNur0Neg= f4GenBr0= f4GenBr0Pos= f4GenBr0Neg = TLorentzVector(0., 0., 0., 0.);
  // -- reco'ed/derived quantities
  fVtxTauDecay= fDirectionTau = fDirectionB = TVector3(0., 0., 0.);
  f4Muon= f4Kaon= f4MuKa = f4Had= f4Nur0Pos= f4Nur0Neg= f4Br0Pos= f4Br0Neg = TLorentzVector(0., 0., 0., 0.);

  // -- pointers to SIGNAL simpleTracks
  fpTmMu=  fpTmKa= fpTmHad1= fpTmHad2= fpTmHad3 = 0;

  fCandPvI = fBrecoPvI = -99;;
  fNTrk = 0;
  for (int i = 0; i < NTRKMAX; ++i) {
    fCorrect[i] = fInRecoil[i] = false;
    fDoca[i] = -99.;
    fPt[i] = -99.;
  }
  // -- gen quantities
  fGenBPt=
    fGenBr0Mass= fGenBr0PosMass= fGenBr0NegMass=
    fGenTauPt= fGenTauDecTime= fGenTauPara= fGenTauPerp=
    fGenTaur0Pt= fGenTaur0PosPt= fGenTaur0NegPt=
    fGenTauHadPt= fGenTauHadMass= fGenTauHadAngle= fGenTauHadPerp= fGenTauHadPara=
    fGen2HadMinMass= fGen2HadMaxMass=
    fGen2HadSSMass= fGen2HadOSminMass= fGen2HadOSmaxMass=
    fGenHad1Pt= fGenHad2Pt= fGenHad3Pt=
    fGenMuPt= fGenKaPt= fGenMuKaMass=
    fGenNur0PosPara= fGenNur0NegPara=
    fGenNuPara= fGenNuPerp= -99.
    ;


  // -- reco quantities
  fBrMass= fBr0PosMass= fBr0NegMass=
    fMuKaDoca3D= fMuKaLip= fMuKaDocaMax=fMuKaFl= fMuKaFls=
    fHadPt= fHadMass= fHadDecTime= fHadFl= fHadFls= fHadDoca3D= fHadLip= fHadDocaMax= fHadOa=
    fTaur0Pt= fTaur0PosPt= fTaur0NegPt=
    fTauHadAngle= fTauHadPerp= fTauHadPara=
    fBTauAngle= fBMuKaAngle=
    f2HadMinMass=f2HadMaxMass=
    f2HadSSMass= f2HadOSminMass= f2HadOSmaxMass=
    fHad1Pt= fHad2Pt= fHad3Pt= fHad1Eta= fHad2Eta= fHad3Eta=
    fMuPt= fKaPt= fMuEta= fKaEta=
    fMuKaDecTime= fMuKaMass= fMuKaPt= fMuKaHadMass=
    fNur0PosPara= fNur0NegPara=-99.
    ;


  candAna::resetAllData();

}


// ----------------------------------------------------------------------
void candAnaBuToMuTauK::genMatch() {
  TGenCand *pCand, *pD, *pT;
  // -- Force initialization!
  fpGenB = fpGenMu = fpGenKa = fpGenTau = fpGenHad1 = fpGenHad2 = fpGenHad3 = fpGenNu = fpGenGa1 = fpGenGa2 = 0;
  fGenVtxBProd = fGenVtxBDecay = fGenVtxTauDecay = fGenDirectionTau = TVector3(0., 0., 0.);
  f4GenTauAtDecay = f4GenNur0Pos = f4GenNur0Neg = f4GenBr0Pos = f4GenBr0Neg = f4GenBr0 = f4GenNur0 = TLorentzVector(0., 0., 0., 0.);

  fGenBPt =
    fGenBr0Mass = fGenBr0PosMass = fGenBr0NegMass =
    fGenTauPt = fGenTauDecTime = fGenTauPara = fGenTauPerp =
    fGenTauHadPt = fGenTauHadMass = fGenTauHadAngle = fGenTauHadPerp = fGenTauHadPara =
    fGen2HadMinMass = fGen2HadMaxMass =
    fGen2HadSSMass = fGen2HadOSminMass = fGen2HadOSmaxMass =
    fGenHad1Pt = fGenHad2Pt = fGenHad3Pt =
    fGenMuPt = fGenKaPt = fGenMuKaMass =
    fGenNur0PosPara =
    fGenNur0NegPara =
    fGenNuPara = fGenNuPerp
    = -99.;

  fNGenPhotons = 0;

  fGenBTmi = -1;

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
	fGenBTmi = fpGenB->fNumber;
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


	break;
      }
    }
  }

  if (0 && fGenBTmi > -1) {
    cout << "genMatch: pT = "
      << fpGenMu->fP.Perp() << " "
      << fpGenKa->fP.Perp() << " "
      << fpGenHad1->fP.Perp() << " "
      << fpGenHad2->fP.Perp() << " "
      << fpGenHad3->fP.Perp() << " "
      << endl;
  } else {
    //    cout << "candAnaBuToMuTauK::genMatch()" << endl;
  }
}

// ----------------------------------------------------------------------
void candAnaBuToMuTauK::recoMatch() {
  fRecoTm = false;
  // -- search for true tracks
  fpTmMu = fpTmKa = fpTmHad1 = fpTmHad2 = fpTmHad3 = 0;
  TSimpleTrack *pT(0);
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i);
    if (pT->getGenIndex() < 0) continue;
    // -- muon
    if (0 != fpGenMu && pT->getGenIndex() == fpGenMu->fNumber) {
      fpTmMu = pT;
    }

    // -- kaon
    if (0 != fpGenKa && pT->getGenIndex() == fpGenKa->fNumber) {
      fpTmKa = pT;
    }

    // -- had1
    if (0 != fpGenHad1 && pT->getGenIndex() == fpGenHad1->fNumber) {
      fpTmHad1 = pT;
    }

    // -- had2
    if (0 != fpGenHad2 && pT->getGenIndex() == fpGenHad2->fNumber) {
      fpTmHad2 = pT;
    }

    // -- had3
    if (0 != fpGenHad3 && pT->getGenIndex() == fpGenHad3->fNumber) {
      fpTmHad3 = pT;
    }


    // -- skip rest if all matches found
    if (0!= fpTmMu && 0 != fpTmKa && 0 != fpTmHad1 && 0 != fpTmHad2 && 0 != fpTmHad3) {
      fRecoTm = true;
      break;
    }
  }
  if (fRecoTm) {
    cout << "recoMatch:pT = "
	 << fpTmMu->getP().Perp() << " "
	 << fpTmKa->getP().Perp() << " "
	 << fpTmHad1->getP().Perp() << " "
	 << fpTmHad2->getP().Perp() << " "
	 << fpTmHad3->getP().Perp() << " "
	 << " indices: "
	 << fpTmMu->getIndex() << " "
	 << fpTmKa->getIndex() << " "
	 << fpTmHad1->getIndex() << " "
	 << fpTmHad2->getIndex() << " "
	 << fpTmHad3->getIndex() << " "
      	 << endl;
  } else {
      //    cout << "candAnaBuToMuTauK::recoMatch()" << endl;
  }

}



// ----------------------------------------------------------------------
// -- search for true candidate
void candAnaBuToMuTauK::candMatch() {
  fpCandTruth = 0;
  fCandTmi = -1;
  int idx(-1), type(-1);
  int muMatched(0), kaMatched(0), had1Matched(0), had2Matched(0), had3Matched(0);
  TAnaCand *pCand(0), *pTau(0);

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    if (CANDTYPE != pCand->fType) continue;

    muMatched = kaMatched = 0;
    for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
      idx = fpEvt->getSigTrack(i)->fIndex;
      type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
      if (0 != fpTmMu && type == 13 && idx == fpTmMu->getIndex()) {
        muMatched = 1;
      }
      if (0 != fpTmKa && type == 321 && idx == fpTmKa->getIndex()) {
        kaMatched = 1;
      }
    }

    for (int iCC = iC+1; iCC < fpEvt->nCands(); ++iCC) {
      pTau = fpEvt->getCand(iCC);
      if (pTau->fMom != iC) continue;

      had1Matched = had2Matched = had3Matched = 0;
      for (int ii = pTau->fSig1; ii <= pTau->fSig2; ++ii) {
	idx = fpEvt->getSigTrack(ii)->fIndex;
	type = TMath::Abs(fpEvt->getSigTrack(ii)->fMCID);
	if (0 != fpTmHad1 && type == 211 && idx == fpTmHad1->getIndex()) {
	  had1Matched = 1;
	}
	if (0 != fpTmHad2 && type == 211 && idx == fpTmHad2->getIndex()) {
	  had2Matched = 1;
	}
	if (0 != fpTmHad3 && type == 211 && idx == fpTmHad3->getIndex()) {
	  had3Matched = 1;
	}
      }
      if (had1Matched && had2Matched  && had3Matched) break;
    }

    if (muMatched && kaMatched  && had1Matched && had2Matched && had3Matched) {
      fCandTmi = iC;
      fpCandTruth = fpEvt->getCand(iC);
      break;
    }
  }

  if (muMatched && kaMatched  && had1Matched && had2Matched && had3Matched) {
    cout << "Matched cand idx = " << fCandTmi << " with all tracks to correct TGC, fCandTmi = " << fCandTmi << endl;
  } else {
    if (fRecoTm) {
      cout << "candAnaBuToMuTauK::candMatch() FAILED" << endl;
    }
  }
}

// ----------------------------------------------------------------------
void candAnaBuToMuTauK::dump() {
}



// ----------------------------------------------------------------------
void candAnaBuToMuTauK::genAnalysis() {
  int verbose(0);
  if (fGenBTmi > -1) {
    f4GenTauHad   = fpGenHad1->fP + fpGenHad2->fP + fpGenHad3->fP;
    fGenDirectionTau = fGenVtxTauDecay - fGenVtxBDecay;
    fGenTauDecTime = fGenDirectionTau.Mag()*MTAU/fpGenTau->fP.Vect().Mag()/TMath::Ccgs();
    fGenTauHadAngle = fGenDirectionTau.Angle(f4GenTauHad.Vect());

    f4GenTauAtDecay = fpGenHad1->fP + fpGenHad2->fP + fpGenHad3->fP  + fpGenNu->fP;
    if (fpGenGa1) {
      f4GenTauAtDecay += fpGenGa1->fP;
      ++fNGenPhotons;
    }
    if (fpGenGa2) {
      f4GenTauAtDecay += fpGenGa2->fP;
      ++fNGenPhotons;
    }

    // if (fNGenPhotons > 0) {
    //   // cout << "photons encountered, skipping event!" << endl;
    //   // return;
    // }

    if (verbose > 0) cout << "gen test" << endl;
    pair<TVector3, TVector3> nucomp = parallelAndPerp(fGenDirectionTau, fpGenNu->fP.Vect());

    fGenNuPara = nucomp.first.Mag();
    fGenNuPerp = nucomp.second.Mag();

    pair<TVector3, TVector3> taucomp = parallelAndPerp(fGenDirectionTau, fpGenTau->fP.Vect());
    fGenTauPara = taucomp.first.Mag();
    fGenTauPerp = taucomp.second.Mag();

    pair<TVector3, TVector3> hcomp = parallelAndPerp(fGenDirectionTau, f4GenTauHad.Vect());
    double evis  = f4GenTauHad.E();
    double mvis  = f4GenTauHad.M();
    double ppar  = hcomp.first.Mag();
    double pperp = hcomp.second.Mag();


    // -- tests and checks ...
    TVector3 v3tau = fpGenTau->fP.Vect();
    TVector3 v3nu  = fpGenNu->fP.Vect();
    TVector3 v3had = f4GenTauHad.Vect();
    TVector3 manu  = v3had + v3nu;

    TLorentzVector v4tau = fpGenTau->fP;
    TLorentzVector v4nu  = fpGenNu->fP;
    TLorentzVector v4had = f4GenTauHad;

    if (verbose > 0) {
      cout << "along tau direction: tau = " << taucomp.first.Mag()
	   << " sum(had + nu) = " << hcomp.first.Mag() + nucomp.first.Mag()
	   << " had = " << hcomp.first.Mag() << " neutrino = " << nucomp.first.Mag()
	   << endl;
      cout << "orthogonal to tau direction: tau = " << taucomp.second.Mag()
	   << " sum(had + nu) = " << hcomp.second.Mag() - nucomp.second.Mag()
	   << " had = " << hcomp.second.Mag() << " neutrino = " << nucomp.second.Mag()
	   << endl;
    }

    // -- tests and checks ...

    pair<double, double> twosol = nuRecoMom0(ppar, pperp, evis, mvis, MTAU);
    fGenTauHadPerp = ppar;
    fGenTauHadPara = pperp;

    // -- calculate reconstructed neutrino 4-momentum
    TVector3 vDirUnit = fGenDirectionTau.Unit();
    TVector3 vRecoPos = twosol.first*vDirUnit - hcomp.second;
    TVector3 vRecoNeg = twosol.second*vDirUnit - hcomp.second;
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

    vector<TLorentzVector> vHad;
    vHad.push_back(fpGenHad1->fP);
    vHad.push_back(fpGenHad2->fP);
    vHad.push_back(fpGenHad3->fP);
    fGen2HadMinMass = minMassPair(vHad);
    fGen2HadMaxMass = maxMassPair(vHad);

    fGenMuPt       = fpGenMu->fP.Perp();
    fGenKaPt       = fpGenKa->fP.Perp();
    fGenMuKaMass   = (fpGenMu->fP + fpGenKa->fP).M();


    double ma =  (vHad[0] + vHad[1]).M();
    double mb =  (vHad[0] + vHad[2]).M();
    double mc =  (vHad[1] + vHad[2]).M();
    if (fpGenHad1->fQ == fpGenHad2->fQ) {
      fGen2HadSSMass = ma;
      if (mb < mc) {
	fGen2HadOSminMass = mb;
	fGen2HadOSmaxMass = mc;
      } else {
	fGen2HadOSminMass = mc;
	fGen2HadOSmaxMass = mb;
      }
    } else if (fpGenHad1->fQ == fpGenHad3->fQ) {
      fGen2HadSSMass = mb;
      if (ma < mc) {
	fGen2HadOSminMass = ma;
	fGen2HadOSmaxMass = mc;
      } else {
	fGen2HadOSminMass = mc;
	fGen2HadOSmaxMass = ma;
      }
    } else if (fpGenHad2->fQ == fpGenHad3->fQ) {
      fGen2HadSSMass = mc;
      if (ma < mb) {
	fGen2HadOSminMass = ma;
	fGen2HadOSmaxMass = mb;
      } else {
	fGen2HadOSminMass = mb;
	fGen2HadOSmaxMass = ma;
      }
    } else {
      cout << "GEN confused?!       " << fpGenHad1->fQ << " " << fpGenHad2->fQ << " " << fpGenHad3->fQ << endl;
    }
  }

  ((TH1D*)fHistDir->Get("ngamma"))->Fill(fNGenPhotons);

}


// ----------------------------------------------------------------------
void candAnaBuToMuTauK::candAnalysis() {
  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(10);
  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(2);

  fTm = (fpCand == fpCandTruth);
  fPvX = fPvY = fPvZ =
    fBr0PosMass = fBr0NegMass =
    fMuKaDoca3D= fMuKaLip= fMuKaDocaMax=fMuKaFl= fMuKaFls=
    fHadPt = fHadDecTime = fHadMass = fHadFl = fHadFls = fHadDoca3D = fHadLip = fHadDocaMax= fHadOa=
    fTauHadAngle = fTauHadPerp = fTauHadPara =
    f2HadMinMass = f2HadMaxMass =
    f2HadSSMass = f2HadOSminMass = f2HadOSmaxMass =
    fHad1Pt = fHad2Pt = fHad3Pt = fHad1Eta = fHad2Eta = fHad3Eta =
    fMuPt = fKaPt = fMuEta = fKaEta =
    fMuKaDecTime = fMuKaMass = fMuKaPt = fMuKaHadMass =
    fNur0PosPara = fNur0NegPara =
    -99.;
  fPvN = -99;
  fMissedTracks = 0;

  fVtxBProd = fPV = fVtxBDecay = fVtxTauDecay = fDirectionTau = fDirectionB = TVector3(0., 0., 0.);
  f4Had = f4Nur0Pos = f4Nur0Neg = f4Br0Pos = f4Br0Neg = TLorentzVector(0., 0., 0., 0.);

  // -- SIGNAL candidate B+ -> mu tau K
  TLorentzVector tlv, tlw;
  if (0) {
    cout << "dump signal candidate " << fpCand->fType << endl;
    dumpCand(fpCand, fpEvt);
  }
  TAnaTrack *pMuon = fpEvt->getSigTrack(fpCand->fSig1);
  TAnaTrack *pKaon = fpEvt->getSigTrack(fpCand->fSig1+1);
  TAnaCand  *pTau  = fpEvt->getCand(fpCand->fDau1);
  TAnaTrack *pHad1 = fpEvt->getSigTrack(pTau->fSig1);
  TAnaTrack *pHad2 = fpEvt->getSigTrack(pTau->fSig1+1);
  TAnaTrack *pHad3 = fpEvt->getSigTrack(pTau->fSig1+2);

  if (pHad1->fPlab.Perp() < pHad2->fPlab.Perp()) {
    TAnaTrack *pHad = pHad1;
    pHad1 = pHad2;
    pHad2 = pHad;
  }
  if (pHad1->fPlab.Perp() < pHad3->fPlab.Perp()) {
    TAnaTrack *pHad = pHad1;
    pHad1 = pHad3;
    pHad3 = pHad;
  }
  if (pHad2->fPlab.Perp() < pHad3->fPlab.Perp()) {
    TAnaTrack *pHad = pHad2;
    pHad2 = pHad3;
    pHad3 = pHad;
  }
  fSignalTracks.clear();
  fSignalTracks.push_back(pMuon);
  fSignalTracks.push_back(pKaon);
  fSignalTracks.push_back(pHad1);
  fSignalTracks.push_back(pHad2);
  fSignalTracks.push_back(pHad3);

  fCandPvI  = fpCand->fPvIdx;
  if (fpRecoilCand) {
    bool foundTrack(false);
    fNTrk = fSignalTracks.size();
    fBrecoPvI = fpRecoilCand->fPvIdx;
    for (unsigned int is = 0; is < fSignalTracks.size(); ++is) {
      TAnaTrack *ps = fSignalTracks[is];
      foundTrack = false;
      double docar(0.);
      for (int ir = fpRecoilCand->fSig1; ir <= fpRecoilCand->fSig2; ++ir) {
	TAnaTrack *pr = fpEvt->getSigTrack(ir);
	if (pr->fIndex == ps->fIndex) {
	  docar = pr->fDouble1;
	  foundTrack = true;
	  break;
	}
      }
      fDoca[is] = ps->fDouble1;
      fPt[is]   = ps->fPlab.Perp();
      if (false == foundTrack) {
	double ntrk = TMath::Power(2, static_cast<int>(is));
	fMissedTracks += ntrk;
	fInRecoil[is] = false;
      } else {
	fInRecoil[is] = true;
      }
    }
  }

  f4Muon.SetXYZM(pMuon->fPlab.X(), pMuon->fPlab.Y(), pMuon->fPlab.Z(), MMUON);
  f4Kaon.SetXYZM(pKaon->fPlab.X(), pKaon->fPlab.Y(), pKaon->fPlab.Z(), MKAON);
  f4MuKa = f4Muon + f4Kaon;

  vector<TLorentzVector> vHad;
  tlv.SetXYZM(pHad1->fPlab.X(), pHad1->fPlab.Y(), pHad1->fPlab.Z(), MPION); vHad.push_back(tlv);
  tlv.SetXYZM(pHad2->fPlab.X(), pHad2->fPlab.Y(), pHad2->fPlab.Z(), MPION); vHad.push_back(tlv);
  tlv.SetXYZM(pHad3->fPlab.X(), pHad3->fPlab.Y(), pHad3->fPlab.Z(), MPION); vHad.push_back(tlv);
  f4Had = vHad[0] + vHad[1] + vHad[2];

  // big guess: take the closest to the gen PV
  int    minI(-1);
  double minDist(99.);
  for (int ip = 0; ip < fpEvt->nPV(); ++ip) {
    TAnaVertex *pv = fpEvt->getPV(ip);
    //    cout << "PV " << ip << ": " << formatTVector3(pv->fPoint, 1) << endl;
    double pvdist = (fGenVtxBProd - pv->fPoint).Mag();
    if (pvdist < minDist) {
      minDist = pvdist;
      minI = ip;
    }
  }
  fVtxBProd     = fpEvt->getPV(minI)->fPoint;
  if (fpCand->fPvIdx > -1 && fpCand->fPvIdx < fpEvt->nPV()) {
    fPV           = fpEvt->getPV(fpCand->fPvIdx)->fPoint;
  }
  fPvX = fPV.X();
  fPvY = fPV.Y();
  fPvZ = fPV.Z();

  fMuKaDoca3D   = fpCand->fPvIP3d;
  fMuKaLip      = fpCand->fPvLip;
  fMuKaDocaMax  = fpCand->fMaxDoca;
  fMuKaFl       = fpCand->fVtx.fD3d;
  fMuKaFls      = fpCand->fVtx.fD3d/fpCand->fVtx.fD3dE;

  fVtxTauDecay  = pTau->fVtx.fPoint;
  fVtxBDecay    = fpCand->fVtx.fPoint;
  fDirectionB   = fVtxBDecay - fPV;
  fDirectionTau = fVtxTauDecay - fVtxBDecay;
  fHadFl        = pTau->fVtx.fD3d;
  fHadFls       = pTau->fVtx.fD3d / pTau->fVtx.fD3dE;
  fHadDoca3D    = pTau->fPvIP3d;
  fHadDocaMax   = pTau->fMaxDoca;
  fHadLip       = pTau->fPvLip;

  //??  fHadDecTime   = fDirectionTau.Mag()*MTAU/f4Had.Rho()/TMath::Ccgs();
  fBTauAngle    = fDirectionB.Angle(fDirectionTau);
  fBMuKaAngle   = fDirectionB.Angle(f4MuKa.Vect());

  if (fTauHadAngle > 2.9) {
    cout << "fTauHadAngle = " << fTauHadAngle
	 << " d(tv, sv) = " << formatTVector3(fVtxBDecay-fVtxTauDecay) << Form(" r = %5.4f", (fVtxBDecay-fVtxTauDecay).Mag())
	 << " in tree = " << pTau->fVtx.fD3d << "/" << pTau->fVtx.fD3dE
	 << endl;
    cout << "PV: reco" << formatTVector3(fVtxBProd, 1) << " gen" << formatTVector3(fGenVtxBProd, 1) << endl;
    cout << "SV: reco" << formatTVector3(fVtxBDecay, 1) << " gen" << formatTVector3(fGenVtxBDecay, 1) << endl;
    cout << "TV: reco" << formatTVector3(fVtxTauDecay, 1) << " gen" << formatTVector3(fGenVtxTauDecay, 1) << endl;
    cout << "muon: reco" << formatTVector3(f4Muon.Vect(), 1);
    if (fpGenMu) {
      cout << " gen" << formatTVector3(fpGenMu->fP.Vect(), 1) << endl;
    } else {
      cout << endl;
    }
    cout << "kaon: reco" << formatTVector3(f4Kaon.Vect(), 1);
    if (fpGenKa) {
      cout << " gen" << formatTVector3(fpGenKa->fP.Vect(), 1) << endl;
    } else {
      cout << endl;
    }
    cout << "had1: reco" << formatTVector3(vHad[0].Vect(), 1);
    if (fpGenHad1) {
      cout << " gen" << formatTVector3(fpGenHad1->fP.Vect(), 1) << endl;
    } else {
      cout << endl;
    }
    cout << "had2: reco" << formatTVector3(vHad[1].Vect(), 1);
    if (fpGenHad2) {
      cout << " gen" << formatTVector3(fpGenHad2->fP.Vect(), 1) << endl;
    } else {
      cout << endl;
    }
    cout << "had3: reco" << formatTVector3(vHad[2].Vect(), 1);
    if (fpGenHad3) {
      cout << " gen" << formatTVector3(fpGenHad3->fP.Vect(), 1) << endl;
    } else {
      cout << endl;
    }
  }

  double ma =  (vHad[0] + vHad[1]).M();
  double mb =  (vHad[0] + vHad[2]).M();
  double mc =  (vHad[1] + vHad[2]).M();
  if (pHad1->fQ == pHad2->fQ) {
    f2HadSSMass = ma;
    if (mb < mc) {
      f2HadOSminMass = mb;
      f2HadOSmaxMass = mc;
    } else {
      f2HadOSminMass = mc;
      f2HadOSmaxMass = mb;
    }
  } else if (pHad1->fQ == pHad3->fQ) {
    f2HadSSMass = mb;
    if (ma < mc) {
      f2HadOSminMass = ma;
      f2HadOSmaxMass = mc;
    } else {
      f2HadOSminMass = mc;
      f2HadOSmaxMass = ma;
    }
  } else if (pHad2->fQ == pHad3->fQ) {
    f2HadSSMass = mc;
    if (ma < mb) {
      f2HadOSminMass = ma;
      f2HadOSmaxMass = mb;
    } else {
      f2HadOSminMass = mb;
      f2HadOSmaxMass = ma;
    }
  } else {
    cout << "RECO confused?!       " << pHad1->fQ << " " << pHad2->fQ << " " << pHad3->fQ << endl;
  }

  // cout << "dump signal truth-identified candidate " << fpCand->fType << endl;
  // dumpCand(fpCand, fpEvt);


  pair<TVector3, TVector3> hcomp = parallelAndPerp(fDirectionTau, f4Had.Vect());
  double evis  = f4Had.E();
  double mvis  = f4Had.M();
  double ppar  = hcomp.first.Mag();
  double pperp = hcomp.second.Mag();
  pair<double, double> twosol = nuRecoMom0(ppar, pperp, evis, mvis, MTAU);
  fTauHadPerp = ppar;
  fTauHadPara = pperp;

  // -- calculate reconstructed neutrino 4-momentum
  TVector3 vDirUnit = fDirectionTau.Unit();
  TVector3 vRecoPos = twosol.first*vDirUnit - hcomp.second;
  TVector3 vRecoNeg = twosol.second*vDirUnit - hcomp.second;
  f4Nur0Pos.SetVectM(vRecoPos, 0.);
  f4Nur0Neg.SetVectM(vRecoNeg, 0.);

  fNur0PosPara = twosol.first;
  fNur0NegPara = twosol.second;

  f4Br0Pos = f4Nur0Pos + f4Had + f4Muon + f4Kaon;
  f4Br0Neg = f4Nur0Neg + f4Had + f4Muon + f4Kaon;

  fBr0PosMass = f4Br0Pos.M();
  fBr0NegMass = f4Br0Neg.M();

  fHadDecTime = pTau->fTau3d;
  fTauHadAngle  = fDirectionTau.Angle(f4Had.Vect());
  fHadOa   = oaTriplet(vHad);
  fHadPt = pTau->fPlab.Perp();
  fHadMass = pTau->fMass;
  f2HadMinMass = minMassPair(vHad);
  f2HadMaxMass = maxMassPair(vHad);
  fHad1Pt  = pHad1->fPlab.Perp();
  fHad2Pt  = pHad2->fPlab.Perp();
  fHad3Pt  = pHad3->fPlab.Perp();
  fHad1Eta = pHad1->fPlab.Eta();
  fHad2Eta = pHad2->fPlab.Eta();
  fHad3Eta = pHad3->fPlab.Eta();
  fMuPt    = pMuon->fPlab.Perp();
  fMuEta   = pMuon->fPlab.Eta();
  fKaPt    = pKaon->fPlab.Perp();
  fKaEta   = pKaon->fPlab.Eta();

  fMuKaDecTime = fpCand->fTau3d;
  fMuKaMass    = fpCand->fMass;
  fMuKaPt      = fpCand->fPlab.Perp();

  fMuKaHadMass = (f4MuKa + f4Had).M();

  //  -- check whether this (TRUTHCAND) candidate is also a reco candidate
  if (CANDTYPE == CANDTRUTH) {
    fTruthRm = false;
    int idx(-1), type(-1);
    int muMatched(0), kaMatched(0), had1Matched(0), had2Matched(0), had3Matched(0);
    TAnaCand *pCand(0), *pTau(0);
    for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
      pCand = fpEvt->getCand(iC);
      if (10032 == pCand->fType) {

	muMatched = kaMatched = 0;
	for (int i = pCand->fSig1; i <= pCand->fSig2; ++i) {
	  idx = fpEvt->getSigTrack(i)->fIndex;
	  type = TMath::Abs(fpEvt->getSigTrack(i)->fMCID);
	  if (0 != fpTmMu && type == 13 && idx == fpTmMu->getIndex()) {
	    muMatched = 1;
	  }
	  if (0 != fpTmKa && type == 321 && idx == fpTmKa->getIndex()) {
	    kaMatched = 1;
	  }
	}

	for (int iCC = iC+1; iCC < fpEvt->nCands(); ++iCC) {
	  pTau = fpEvt->getCand(iCC);
	  if (pTau->fMom != iC) continue;

	  had1Matched = had2Matched = had3Matched = 0;
	  for (int ii = pTau->fSig1; ii <= pTau->fSig2; ++ii) {
	    idx = fpEvt->getSigTrack(ii)->fIndex;
	    type = TMath::Abs(fpEvt->getSigTrack(ii)->fMCID);
	    if (0 != fpTmHad1 && type == 211 && idx == fpTmHad1->getIndex()) {
	      had1Matched = 1;
	    }
	    if (0 != fpTmHad2 && type == 211 && idx == fpTmHad2->getIndex()) {
	      had2Matched = 1;
	    }
	    if (0 != fpTmHad3 && type == 211 && idx == fpTmHad3->getIndex()) {
	      had3Matched = 1;
	    }
	  }
	  if (had1Matched && had2Matched  && had3Matched) break;
	}

	if (muMatched && kaMatched  && had1Matched && had2Matched && had3Matched) {
	  fTruthRm = true;
	  break;
	}
      }

    }
  }
  candAna::candAnalysis();
}


// ----------------------------------------------------------------------
void candAnaBuToMuTauK::candEvaluation() {

}

// ----------------------------------------------------------------------
void candAnaBuToMuTauK::moreReducedTree(TTree *t) {
  cout << "candAnaBuToMuTauK::moreReducedTree for fTree = " << t << " in directory " << fHistDir->GetName() << endl;
  // -- gen level
  t->Branch("gptb",        &fGenBPt,           "gptb/D");
  t->Branch("gmtauhad",    &fGenTauHadMass,    "gmtauhad/D");
  t->Branch("gm2hadmax",   &fGen2HadMaxMass,   "gm2hadmax/D");
  t->Branch("gm2hadmin",   &fGen2HadMinMass,   "gm2hadmin/D");
  t->Branch("gm2hadss",    &fGen2HadSSMass,    "gm2hadss/D");
  t->Branch("gm2hadosmin", &fGen2HadOSminMass, "gm2hadosmin/D");
  t->Branch("gm2hadosmax", &fGen2HadOSmaxMass, "gm2hadosmax/D");
  t->Branch("gmbr0pos",    &fGenBr0PosMass,    "gmbr0pos/D");
  t->Branch("gmbr0neg",    &fGenBr0NegMass,    "gmbr0neg/D");
  t->Branch("gmbr0",       &fGenBr0Mass,       "gmbr0/D");
  t->Branch("gt0tau",      &fGenTauDecTime,    "gt0tau/D");
  t->Branch("gpttau",      &fGenTauPt,         "gpttau/D");
  t->Branch("gpttauhad",   &fGenTauHadPt,      "gpttauhad/D");
  t->Branch("gatauhad",    &fGenTauHadAngle,   "gatauhad/D");
  t->Branch("gparanur0pos",&fGenNur0PosPara,   "gparanur0pos/D");
  t->Branch("gparanur0neg",&fGenNur0NegPara,   "gparanur0neg/D");
  t->Branch("gparanu",     &fGenNuPara,        "gparanu/D");
  t->Branch("gperpnu",     &fGenNuPerp,        "gperpnu/D");
  t->Branch("gparatau",    &fGenTauPara,       "gparatau/D");
  t->Branch("gperptau",    &fGenTauPerp,       "gperptau/D");
  t->Branch("gperptauhad", &fGenTauHadPerp,    "gperptauhad/D");
  t->Branch("gparatauhad", &fGenTauHadPara,    "gparatauhad/D");
  t->Branch("gpttaur0pos", &fGenTaur0PosPt,    "gpttaur0pos/D");
  t->Branch("gpttaur0neg", &fGenTaur0NegPt,    "gpttaur0neg/D");
  t->Branch("gpttaur0",    &fGenTaur0Pt,       "gpttaur0/D");
  t->Branch("gptmu",       &fGenMuPt,          "gptmu/D");
  t->Branch("gptka",       &fGenKaPt,          "gptka/D");
  t->Branch("gmmuka",      &fGenMuKaMass,      "gmmuka/D");
  t->Branch("gpthad1",     &fGenHad1Pt,        "gpthad1/D");
  t->Branch("gpthad2",     &fGenHad2Pt,        "gpthad2/D");
  t->Branch("gpthad3",     &fGenHad3Pt,        "gpthad3/D");

  t->Branch("tm",          &fTm,             "tm/O");
  t->Branch("rm",          &fTruthRm,        "rm/O");

  t->Branch("t0muka",      &fMuKaDecTime,    "t0muka/D");
  t->Branch("mmukahad",    &fMuKaHadMass,    "mmukahad/D");

  t->Branch("doca3dmuka",  &fMuKaDoca3D,     "doca3dmuka/D");
  t->Branch("docamaxmuka", &fMuKaDocaMax,    "docamaxmuka/D");
  t->Branch("lipmuka",     &fMuKaLip,        "lipmuka/D");

  t->Branch("ptmuka",      &fMuKaPt,         "ptmuka/D");
  t->Branch("mmuka",       &fMuKaMass,       "mmuka/D");
  t->Branch("flmuka",      &fMuKaFl,         "flmuka/D");
  t->Branch("flsmuka",     &fMuKaFls,        "flsmuka/D");

  t->Branch("mhad",        &fHadMass,        "mhad/D");
  t->Branch("m2hadmax",    &f2HadMaxMass,    "m2hadmax/D");
  t->Branch("m2hadmin",    &f2HadMinMass,    "m2hadmin/D");
  t->Branch("m2hadss",     &f2HadSSMass,     "m2hadss/D");
  t->Branch("m2hadosmin",  &f2HadOSminMass,  "m2hadosmin/D");
  t->Branch("m2hadosmax",  &f2HadOSmaxMass,  "m2hadosmax/D");
  t->Branch("mbr0pos",     &fBr0PosMass,     "mbr0pos/D");
  t->Branch("mbr0neg",     &fBr0NegMass,     "mbr0neg/D");
  t->Branch("abtau",       &fBTauAngle,      "abtau/D");
  t->Branch("t0had",       &fHadDecTime,     "t0had/D");
  t->Branch("flhad",       &fHadFl,          "flhad/D");
  t->Branch("flshad",      &fHadFls,         "flshad/D");
  t->Branch("pthad",       &fHadPt,          "pthad/D");
  t->Branch("doca3dhad",   &fHadDoca3D,      "doca3dhad/D");
  t->Branch("docamaxhad",  &fHadDocaMax,     "docamaxhad/D");
  t->Branch("oahad",       &fHadOa,          "oahad/D");
  t->Branch("liphad",      &fHadLip,         "liphad/D");
  t->Branch("pthad",       &fHadPt,          "pthad/D");
  t->Branch("atauhad",     &fTauHadAngle,    "atauhad/D");
  t->Branch("abmuka",      &fBMuKaAngle,     "abmuka/D");
  t->Branch("paranur0pos", &fNur0PosPara,    "paranur0pos/D");
  t->Branch("paranur0neg", &fNur0NegPara,    "paranur0neg/D");
  t->Branch("perptauhad",  &fTauHadPerp,     "perptauhad/D");
  t->Branch("paratauhad",  &fTauHadPara,     "paratauhad/D");
  t->Branch("pttaur0pos",  &fTaur0PosPt,     "pttaur0pos/D");
  t->Branch("pttaur0neg",  &fTaur0NegPt,     "pttaur0neg/D");
  t->Branch("pttaur0",     &fTaur0Pt,        "pttaur0/D");
  t->Branch("ptmu",        &fMuPt,           "ptmu/D");
  t->Branch("etamu",       &fMuEta,          "etamu/D");
  t->Branch("ptka",        &fKaPt,           "ptka/D");
  t->Branch("etaka",       &fKaEta,          "etaka/D");
  t->Branch("pthad1",      &fHad1Pt,         "pthad1/D");
  t->Branch("pthad2",      &fHad2Pt,         "pthad2/D");
  t->Branch("pthad3",      &fHad3Pt,         "pthad3/D");
  t->Branch("etahad1",     &fHad1Eta,        "etahad1/D");
  t->Branch("etahad2",     &fHad2Eta,        "etahad2/D");
  t->Branch("etahad3",     &fHad3Eta,        "etahad3/D");
  t->Branch("missedtracks",&fMissedTracks,   "missedtracks/I");

  t->Branch("ntrk",    &fNTrk,   "ntrk/I");
  t->Branch("tcorrect",&fCorrect,"tcorrect[ntrk]/O");
  t->Branch("inrecoil",&fInRecoil,"inrecoil[ntrk]/O");
  t->Branch("doca",    &fDoca,    "doca[ntrk]/D");
  t->Branch("pt",      &fPt,      "pt[ntrk]/D");
}

// ----------------------------------------------------------------------
void candAnaBuToMuTauK::bookHist() {
  cout << "==>candAnaBuToMuTauK: bookHist in fHistDir = " << fHistDir->GetName() << endl;
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
void candAnaBuToMuTauK::readCuts(string filename, int dump) {

  candAna::readCuts(filename, dump);


  fCutFile = filename;

  if (dump) cout << "==> candAnaBuToMuTauK: Reading " << fCutFile << " for cut settings" << endl;
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

  }

}


// ----------------------------------------------------------------------
void candAnaBuToMuTauK::printGenBDecays() {

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
