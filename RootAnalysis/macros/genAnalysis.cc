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
genAnalysis::genAnalysis(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> genAnalysis: constructor..." << endl;

  f511Mass = f521Mass = f531Mass = f5122Mass = f541Mass = 0;
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

  if (0) yVsEta();
  if (0) genB();
  if (0) printBdecays();
  if (0) bbbarCrossSection();
  if (1) validateLb2PMuNu();



}


// ----------------------------------------------------------------------
void genAnalysis::yVsEta() {

  static int first(1);
  if (1 == first) {
    first = 0;
    new TH2D("yeta1", "yeta1", 50, 0., 2.5, 50, 0., 2.5);
    new TH2D("yeta2", "yeta2", 50, 0., 2.5, 50, 0., 2.5);
    new TH2D("yetaMin", "yetaMin", 50, 0., 2.5, 50, 0., 2.5);
    new TH2D("yetaMax", "yetaMax", 50, 0., 2.5, 50, 0., 2.5);
    new TH2D("yetaMin3", "yetaMin3", 50, 0., 2.5, 50, 0., 2.5);
    new TH2D("yetaMax3", "yetaMax3", 50, 0., 2.5, 50, 0., 2.5);
  }


  TGenCand *pCand(0);
  int aid(0);
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    aid = TMath::Abs(pCand->fID);
    if (521 == aid) {
      TGenCand *pD(0), *pJpsi(0), *pK(0), *pMu1(0), *pMu2(0);
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = fpEvt->getGenCand(iD);
	if (443 == TMath::Abs(pD->fID)) {
	  pJpsi = pD;
	}
	if (321 == TMath::Abs(pD->fID)) {
	  pK = pD;
	}
      }
      if (pJpsi) {
	for (int iD = pJpsi->fDau1; iD <= pJpsi->fDau2; ++iD) {
	  pD = fpEvt->getGenCand(iD);
	  if (13 == pD->fID) {
	    pMu1 = pD;
	  }
	  if (-13 == pD->fID) {
	    pMu2 = pD;
	  }
	}

	if (pCand && pMu1 && pMu2 && pK) {
	  double yB = TMath::Abs(pCand->fP.Rapidity());
	  double eta1 = TMath::Abs(pMu1->fP.PseudoRapidity());
	  double eta2 = TMath::Abs(pMu2->fP.PseudoRapidity());
	  double etaK = TMath::Abs(pK->fP.PseudoRapidity());
	  double etaMin = (eta1 < eta2? eta1: eta2);
	  double etaMax = (eta1 > eta2? eta1: eta2);
	  double etaMin3 = (etaK < etaMin? etaK: etaMin);
	  double etaMax3 = (etaK > etaMax? etaK: etaMin);

	  ((TH1D*)fpHistFile->Get("yeta1"))->Fill(yB, eta1);
	  ((TH1D*)fpHistFile->Get("yeta2"))->Fill(yB, eta2);
	  ((TH1D*)fpHistFile->Get("yetaMin"))->Fill(yB, etaMin);
	  ((TH1D*)fpHistFile->Get("yetaMax"))->Fill(yB, etaMax);
	  ((TH1D*)fpHistFile->Get("yetaMin3"))->Fill(yB, etaMin3);
	  ((TH1D*)fpHistFile->Get("yetaMax3"))->Fill(yB, etaMax3);
	  return;
	}

      }



    }
  }

}




// ----------------------------------------------------------------------
void genAnalysis::endAnalysis() {
  compare2PDG(0, 2014, true);

  for (unsigned int i = 0; i < fRunEvents.size(); ++i) {
    //    cout << fRunEvents[i].first << " " << fRunEvents[i].second << endl;
  }
}

// ----------------------------------------------------------------------
void genAnalysis::genB() {

  static int first(1);
  if (1 == first) {
    first = 0;
    static const double aparticles[] = {511, 521, 531, 5122};
    vector<int> particles(aparticles, aparticles + sizeof(aparticles)/sizeof(aparticles[0]));
    for (unsigned int i = 0; i < particles.size(); ++i) {
      new TH1D(Form("pt%d", particles[i]), Form("pt%d", particles[i]), 100, 0., 50.);
      new TH1D(Form("cpt%d", particles[i]), Form("pt%d", particles[i]), 30, 0., 30.);
      new TH1D(Form("eta%d", particles[i]), Form("eta%d", particles[i]), 50, -10., 10.);
      new TH1D(Form("ceta%d", particles[i]), Form("ceta%d", particles[i]), 60, -3., 3.);
      new TH1D(Form("mom%d", particles[i]), Form("mom%d", particles[i]), 500, 500., 1000.);
      new TH1D(Form("iso%d", particles[i]), Form("iso%d", particles[i]), 51, 0., 1.02);
    }
  }

  TGenCand *pCand, *pD;
  int mom(0), aid(0);
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    aid = TMath::Abs(pCand->fID);

    // -- ignore the oscillated ones (in this loop look at produced B,
    //    in compare2PDG the mixed ones will be used for the lifetime)
    if (531 == aid || 511 == aid) {
      if (aid == TMath::Abs(fpEvt->getGenCand(pCand->fMom1)->fID)) {
	continue;
      }
    }

    if (531 == aid || 521 == aid || 511 == aid || 5122 == aid) {
      compare2PDG(pCand, 2014, false);
      double iso = isolation(pCand);
      pD = fpEvt->getGenCand(pCand->fMom1);
      mom = TMath::Abs(pD->fID);
      ((TH1D*)fpHistFile->Get(Form("pt%d", aid)))->Fill(pCand->fP.Perp());
      ((TH1D*)fpHistFile->Get(Form("cpt%d", aid)))->Fill(pCand->fP.Perp());
      ((TH1D*)fpHistFile->Get(Form("eta%d", aid)))->Fill(pCand->fP.Eta());
      ((TH1D*)fpHistFile->Get(Form("ceta%d", aid)))->Fill(pCand->fP.Eta());
      ((TH1D*)fpHistFile->Get(Form("mom%d", aid)))->Fill(mom);
      // -- ignore underflow!
      if (iso > -0.5) ((TH1D*)fpHistFile->Get(Form("iso%d", aid)))->Fill(iso);

      if (511 == aid) f511Mass = pCand->fP.M();
      if (521 == aid) f521Mass = pCand->fP.M();
      if (531 == aid) f531Mass = pCand->fP.M();
      if (541 == aid) f541Mass = pCand->fP.M();
      if (5122 == aid) f5122Mass = pCand->fP.M();

    }
  }

}


// ----------------------------------------------------------------------
void genAnalysis::printBdecays() {

  TGenCand *pCand, *pD;
  cout << "gen block with " << fpEvt->nGenCands() << " gen cands" << endl;
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);

    if (TMath::Abs(pCand->fID) == 531 || TMath::Abs(pCand->fID) == 511) {
      cout << "--- Event " << fEvent << " ------------------------------------------------------------------" << endl;
      pCand->dump(2);
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = fpEvt->getGenCand(iD);
	pD->dump(2);
	//	cout << "pT: " << pD->fP.Perp() << " eta: " << pD->fP.Eta() << endl;
      }
    }
  }


  return;

  TAnaCand *pC;
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pC = fpEvt->getCand(iC);
    cout << "typ " << pC->fType << endl;
  }

  if (fpEvt->nRecTracks() > 0) {
    TAnaTrack *pT;
    cout << "rec tracks: " << fpEvt->nRecTracks() << endl;
    for (int it = 0; it < fpEvt->nRecTracks(); ++it) {
      pT = fpEvt->getRecTrack(it);
      pT->dump();
    }
  }

  TString a;
  int ps(0);
  bool result(false), wasRun(false), error(false);

  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i];
    ps = fpEvt->fHLTPrescale[i];
    wasRun = fpEvt->fHLTWasRun[i];
    result = fpEvt->fHLTResult[i];
    error  = fpEvt->fHLTError[i];

    if (wasRun && result) {
      cout << a << " " << ps << endl;
    }
  }


}


// ----------------------------------------------------------------------
void genAnalysis::validateLb2PMuNu() {

  static bool first(true);
  if (first) {
    first = false;
    new TH1D("prpt", "pupt", 100, 0., 50.);
    new TH1D("mupt", "mupt", 100, 0., 50.);
    new TH1D("preta", "pueta", 52, -2.6, 2.6);
    new TH1D("mueta", "mueta", 52, -2.6, 2.6);
    new TH1D("mprmu", "mprmu", 120, 0., 6.);
  }

  int lbx, lby;

  TGenCand *pCand, *pD, *pMom;
  //  cout << "gen block with " << fpEvt->nGenCands() << " gen cands" << endl;
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);

    if (TMath::Abs(pCand->fID) == 5122) {
      lbx = static_cast<int>(1000*pCand->fP.X());
      lby = static_cast<int>(1000*pCand->fP.Y());
      pair<int, int> runEvent(lbx, lby);
      //      cout << lbx << " " << lby << endl;
      bool duplicate(false);
      for (unsigned int i = 0; i < fRunEvents.size(); ++i) {
	if (fRunEvents[i].first == lbx && fRunEvents[i].second == lby) {
	  duplicate = true;
	  cout << "duplicate proton" << endl;
	}
      }
      fRunEvents.push_back(runEvent);
      // cout << "--- Event " << fEvent << " ------------------------------------------------------------------" << endl;
      // pCand->dump(2);
      TGenCand *pPr(0), *pMu(0), *pNu(0);
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = fpEvt->getGenCand(iD);
	pMom = fpEvt->getGenCand(pD->fMom1);
	// skip daughters that are not direct daughters of Lb (this removes Lb -> Lc X)
	if (TMath::Abs(pMom->fID) != 5122) continue;
	if (TMath::Abs(pD->fID) == 13) {
	  pMu = pD;
	}
	if (TMath::Abs(pD->fID) == 2212) {
	  pPr = pD;
	}
	if (TMath::Abs(pD->fID) == 14) {
	  pNu = pD;
	}
      }
      if (pPr && pMu && pNu) {
	if (pPr->fP.Perp() < 3.5 || pMu->fP.Perp() < 3.5) {
	  cout << "--- Event " << fEvent << " ------------------------------------------------------------------" << endl;
	  pCand->dump(2);
	  pPr->dump(2);
	  pMu->dump(2);
	  pNu->dump(2);
	}
	((TH1D*)gDirectory->Get("mupt"))->Fill(pMu->fP.Perp());
	((TH1D*)gDirectory->Get("prpt"))->Fill(pPr->fP.Perp());
	((TH1D*)gDirectory->Get("mueta"))->Fill(pMu->fP.Eta());
	((TH1D*)gDirectory->Get("preta"))->Fill(pPr->fP.Eta());
	TLorentzVector m = pPr->fP + pMu->fP;
	((TH1D*)gDirectory->Get("mprmu"))->Fill(m.M());
      }
    }
  }


  return;

}


// ----------------------------------------------------------------------
void genAnalysis::bbbarCrossSection() {

  TGenCand *pCand;
  int muType(0), evtType(0), nevt(0), bevt(0), bacc(0);
  bool acc(false);
  double pt(0.), eta(0.);
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);

    if (TMath::Abs(pCand->fID) == 13) {
      //      cout << "muon at " << iC << endl;
      muType = muonType(pCand);
      pt = pCand->fP.Perp();
      eta= pCand->fP.Eta();

      evtType |= muType;

      if (pt > 5 && eta < 2.1 && eta > -2.1) {
	acc = true;
      } else {
	acc = false;
      }

      ((TH1D*)fpHistFile->Get("h1"))->Fill(muType);
      if (muType > 0) ((TH1D*)fpHistFile->Get("h100"))->Fill(pt);
      if (muType ==1) ((TH1D*)fpHistFile->Get("h110"))->Fill(pt);
      if (muType ==1 && acc) ((TH1D*)fpHistFile->Get("h112"))->Fill(pt);

      if (muType > 0)         ++nevt;
      if (muType == 1)        ++bevt;
      if (muType == 1 && acc) ++bacc;

    }

  }

  ((TH1D*)fpHistFile->Get("e1"))->Fill(evtType);
  if (evtType > 0) ((TH1D*)fpHistFile->Get("e100"))->Fill(pt);
  if (evtType ==1) ((TH1D*)fpHistFile->Get("e110"))->Fill(pt);
  if (bevt    > 0) ((TH1D*)fpHistFile->Get("e111"))->Fill(pt);
  if (bacc    > 0) ((TH1D*)fpHistFile->Get("e112"))->Fill(pt);

}


// ----------------------------------------------------------------------
int genAnalysis::muonType(TGenCand *pCand) {
  int id(999);
  int rest(0), ganz(0);
  int momI = pCand->fMom1;
  TGenCand *pMom;
  bool foundT(false), foundC(false), foundB(false), light(false);

  int result(1024);

  string str("");
  str += Form(" %i", pCand->fID);

  while (momI > 1 && momI < fpEvt->nGenCands()) {
    pMom = fpEvt->getGenCand(momI);
    momI = pMom->fMom1;

    id  = TMath::Abs(pMom->fID);
    str += Form(" %i", pMom->fID);

    rest =  id%1000;
    ganz =  id/1000;

    // -- hit a string
    if (rest == 92) {
      break;
    }

    // -- hit a quark
    if (rest < 6) {
      break;
    }

    if (ganz == 5) {
      foundB = 1;
      break;
    }

    if (ganz == 4) {
      foundC = 1;
      break;
    }

    if (rest == 15) {
      foundT = 1;
    }

    if (15 < rest && rest < 300) {
      light = 1;
    }

    if (rest > 399 && rest < 499) {
      foundC = 1;
    }

    if (rest > 499 && rest < 599) {
      foundB = 1;
      break;
    }
  }

  result = 0;

  if (foundB) result += 1;
  if (foundC) result += 2;
  if (foundT) result += 4;
  if (light)  result += 8;

//   if (result > -1) {
//     cout << fEvent << " " << str << " --> " << result << endl;
//   }

  return result;
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

  new TH1D("h1",   "mu type", 20, 0., 20.);
  new TH1D("h100", "mu pt 0", 100, 0., 20.);
  new TH1D("h110", "mu pt ==1", 100, 0., 20.);
  new TH1D("h112", "mu pt ==1", 100, 0., 20.);

  new TH1D("e1",  "evt type", 20, 0., 20.);
  new TH1D("e100", "evt pt 0", 100, 0., 20.);
  new TH1D("e110", "evt pt ==1", 100, 0., 20.);
  new TH1D("e111", "evt pt ==1", 100, 0., 20.);
  new TH1D("e112", "evt acc ==1", 100, 0., 20.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",     &fRun,     "run/I");

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
  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "NTOTAL")) {
      NTOTAL = int(CutValue); ok = 1;
      if (dump) cout << "NTOTAL:           " << NTOTAL << endl;
    }

    if (!strcmp(CutName, "XSECTION")) {
      XSECTION = double(CutValue); ok = 1;
      if (dump) cout << "XSECTION:           " << XSECTION << endl;
    }

    if (!ok) cout << "==> genAnalysis: ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;

}


// ----------------------------------------------------------------------
double genAnalysis::isolation(TGenCand *pCand) {

  // -- not interested outside of acceptance
  if (TMath::Abs(pCand->fP.Eta()) > 2.4) return -1.;

  TGenCand *pC;
  double iso(0.);
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pC = fpEvt->getGenCand(iC);
    if (isStableCharged(pC->fID)) {
      if (fpEvt->isAncestor(pCand, pC)) continue;
      if (pCand->fP.DeltaR(pC->fP) < 0.7 && pC->fP.Perp() > 0.9) {
	iso += pC->fP.Perp();
      }
    }
  }

  return (pCand->fP.Perp()/(pCand->fP.Perp() + iso));
}


// ----------------------------------------------------------------------
void genAnalysis::compare2PDG(TGenCand *pCand, int year, bool finalize) {

  static double M511(0.), M521(0.), M531(0.), M541(0.), M5122(0.);
  static double L511(0.), L521(0.), L531(0.), L541(0.), L5122(0.);

  if (2014 == year) {
    M511 = 5.27958;
    M521 = 5.27926;
    M531 = 5.36677;
    M541 = 6.2756;
    M5122= 5.6195;

    L511 = 455.4;
    L521 = 491.1;
    L531 = 453.3;
    L541 = 135.5;
    L5122= 435;
  }

  static int first(1);
  if (1 == first) {
    first = 0;
    vector<int> particles = defVector(5, 511, 521, 531, 541, 5122);
    for (unsigned int i = 0; i < particles.size(); ++i) {
      new TH1D(Form("t%d", particles[i]), Form("t%d", particles[i]), 50, 0., 5000.);
      new TH1D(Form("tpos%d", particles[i]), Form("tpos%d", particles[i]), 50, 0., 5000.);
      new TH1D(Form("tneg%d", particles[i]), Form("tneg%d", particles[i]), 50, 0., 5000.);
    }

    new TH1D("m511",  "m511",  100, 5.2794, 5.2796);
    new TH1D("m521",  "m521",  100, 5.279, 5.280);
    new TH1D("m531",  "m531",  100, 5.360, 5.370);
    new TH1D("m541",  "m541",  100, 6.270, 6.280);
    new TH1D("m5122", "m5122", 100, 5.610, 5.630);
  }


  if (finalize) {
    // -- Masses
    cout << "Masses" << endl;
    cout << Form("B0: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
		 ((TH1D*)fpHistFile->Get("m511"))->GetMean(),
		 M511,
		 ((TH1D*)fpHistFile->Get("m511"))->GetMean() - M511
		 )
	 << endl;

    cout << Form("B+: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
		 ((TH1D*)fpHistFile->Get("m521"))->GetMean(),
		 M521,
		 ((TH1D*)fpHistFile->Get("m521"))->GetMean() - M521
		 )
	 << endl;

    cout << Form("Bs: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
		 ((TH1D*)fpHistFile->Get("m531"))->GetMean(),
		 M531,
		 ((TH1D*)fpHistFile->Get("m531"))->GetMean() - M531
		 )
	 << endl;

    cout << Form("Lb: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
		 ((TH1D*)fpHistFile->Get("m5122"))->GetMean(),
		 M5122,
		 ((TH1D*)fpHistFile->Get("m5122"))->GetMean() - M5122
		 )
	 << endl;

    if (((TH1D*)fpHistFile->Get("m541"))->GetEntries() > 1)
      cout << Form("Bc: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
		   ((TH1D*)fpHistFile->Get("m541"))->GetMean(),
		   M541,
		   ((TH1D*)fpHistFile->Get("m541"))->GetMean() - M541
		   )
	   << endl;

    // -- Lifetimes
    cout << "Lifetime" << endl;
    double t(1.), tE(1.), chi2(0);
    TH1D *h = (TH1D*)fpHistFile->Get("t521");
    TF1 *f(0);
    gStyle->SetOptFit(1);
    if (h->GetEntries() > 100) {
      h->Fit("expo", "lq");
      //      gPad->SaveAs("t521.pdf");
      f = (TF1*)h->GetFunction("expo");
      chi2 = f->GetChisquare()/f->GetNDF();
      t    = -1./f->GetParameter(1);
      tE   = -t*f->GetParError(1)/f->GetParameter(1);
      cout << Form("B+: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L521, (t-L521)/tE, chi2) << endl;
    }
    h = (TH1D*)fpHistFile->Get("t511");
    if (h->GetEntries() > 100) {
      h->Fit("expo", "ql");
      //      gPad->SaveAs("t511.pdf");
      f = (TF1*)h->GetFunction("expo");
      chi2 = f->GetChisquare()/f->GetNDF();
      t    = -1./f->GetParameter(1);
      tE   = -t*f->GetParError(1)/f->GetParameter(1);
      cout << Form("B0: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L511, (t-L511)/tE, chi2) << endl;
    }

    h = (TH1D*)fpHistFile->Get("t531");
    if (h->GetEntries() > 100) {
      h->Fit("expo", "ql");
      //      gPad->SaveAs("t531.pdf");
      f = (TF1*)h->GetFunction("expo");
      chi2 = f->GetChisquare()/f->GetNDF();
      t    = -1./f->GetParameter(1);
      tE   = -t*f->GetParError(1)/f->GetParameter(1);
      cout << Form("Bs: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L531, (t-L531)/tE, chi2) << endl;
    }

    h = (TH1D*)fpHistFile->Get("t5122");
    if (h->GetEntries() > 100) {
      h->Fit("expo", "ql");
      //      gPad->SaveAs("t5122.pdf");
      f = (TF1*)h->GetFunction("expo");
      chi2 = f->GetChisquare()/f->GetNDF();
      t    = -1./f->GetParameter(1);
      tE   = -t*f->GetParError(1)/f->GetParameter(1);
      cout << Form("Lb: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L5122, (t-L5122)/tE, chi2) << endl;
    }

    h = (TH1D*)fpHistFile->Get("t541");
    if (h->GetEntries() > 100) {
      h->Fit("expo", "ql");
      //      gPad->SaveAs("t541.pdf");
      f = (TF1*)h->GetFunction("expo");
      chi2 = f->GetChisquare()/f->GetNDF();
      t    = -1./f->GetParameter(1);
      tE   = -t*f->GetParError(1)/f->GetParameter(1);
      cout << Form("Bc: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L541, (t-L541)/tE, chi2) << endl;
    }

    //  -- save at the end to remove the intermittent root printout
    ((TH1D*)fpHistFile->Get("t521"))->Draw();
    gPad->SaveAs("t521.pdf");

    ((TH1D*)fpHistFile->Get("t511"))->Draw();
    gPad->SaveAs("t511.pdf");

    ((TH1D*)fpHistFile->Get("t531"))->Draw();
    gPad->SaveAs("t531.pdf");

    ((TH1D*)fpHistFile->Get("t5122"))->Draw();
    gPad->SaveAs("t5122.pdf");

    ((TH1D*)fpHistFile->Get("t541"))->Draw();
    gPad->SaveAs("t541.pdf");

    return;
  }

  int aid = TMath::Abs(pCand->fID);

  if (511 == aid || 521 == aid || 531 == aid || 541 == aid || 5122 == aid) {
  } else {
    cout << "wrong aid for compare2PDG, returning ... " << endl;
    return;
  }


  ((TH1D*)fpHistFile->Get(Form("m%d", aid)))->Fill(pCand->fP.M());

  // -- get daughter, and check for oscillated cand
  TGenCand *pDau  = fpEvt->getGenCand(pCand->fDau1);
  if (aid == TMath::Abs(pDau->fID)) pDau  = fpEvt->getGenCand(pDau->fDau1);

  double x        = (pCand->fV - pDau->fV).Mag();
  double lifetime = x * pCand->fP.M() / pCand->fP.Rho() / TMath::Ccgs();
  lifetime *= 1.e6*299792458;
  if (0 && 531 == aid && lifetime < 20) {
    cout << "----------------------------------------------------------------------" << endl;
    fpEvt->dumpGenBlock();
    cout << "pCand at " << pCand->fNumber << " pDau at " << pDau->fNumber << endl;
    cout << "----------------------------------------------------------------------" << endl;
  }

  ((TH1D*)fpHistFile->Get(Form("t%d", aid)))->Fill(lifetime);
  if (pDau->fID > 0) ((TH1D*)fpHistFile->Get(Form("tpos%d", aid)))->Fill(lifetime);
  if (pDau->fID < 0) ((TH1D*)fpHistFile->Get(Form("tneg%d", aid)))->Fill(lifetime);

}
