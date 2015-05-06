#include "genAnalysis.hh"
#include "common/HFMasses.hh"
#include "common/util.hh"

#include "TRandom.h"
#include <cmath>
#include <map>

using namespace std;

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: bin/runInclB -r genAnalysis -f test.root
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// ----------------------------------------------------------------------
genAnalysis::genAnalysis(TChain *tree, TString evtClassName): treeReader01(tree, evtClassName) {
  cout << "==> genAnalysis: constructor..." << endl;
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

  printMuons();
  //  bbbarCrossSection();
  // -- initialize all variables
  //  initVariables(); 

}



// ----------------------------------------------------------------------
void genAnalysis::endAnalysis() {

  double ftilde(0.11); 
  double f = 2.*ftilde*(1-ftilde) + ftilde*ftilde; 
  double n0 = ((TH1D*)fpHistFile->Get("e111"))->GetSumOfWeights(); 
  double n1 = ((TH1D*)fpHistFile->Get("e112"))->GetSumOfWeights();
  double c0 = ((TH1D*)fpHistFile->Get("h110"))->GetSumOfWeights(); 
  double c1 = ((TH1D*)fpHistFile->Get("h112"))->GetSumOfWeights();

  if (NTOTAL > 0) {
    cout << "XSECTION:            " << XSECTION << endl;
    cout << "NTOTAL:              " << NTOTAL << endl;
    cout << "b->mu events:        " << n0 << endl;
    cout << "b->mu decays:        " << c0 << endl;
    cout << "b->mu events in acc: " << n1 << endl;
    cout << "b->mu decays in acc: " << c1 << endl;
    
    cout << "bbbar xsection:       " << n0/(f*NTOTAL)*XSECTION << endl;
    cout << "bbbar xsection:       " << c0/(2.*ftilde*NTOTAL)*XSECTION << endl;
    cout << "b-mu  xsection:       " << n0/(NTOTAL)*XSECTION << endl;
    cout << "b-mu in acc xsection: " << n1/(NTOTAL)*XSECTION << endl;
  }
}

// ----------------------------------------------------------------------
void genAnalysis::printMuons() {

  static int first(1); 
  if (1 == first) {
    first = 0; 
    TH1D *h0 = new TH1D("muMom", "muon direct mother", 30, -1., 1.); 
    h0->SetLabelSize(0.03, "X");
    TH1D *h1 = new TH1D("muAnc", "muon ancestor", 64, 0., 64.); 
    vector<int> mid; 
    static const int arr[] = {511, 521, 531, 5122, 411, 421, 431, 4122, 130};
    vector<int> mId (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    int center = h0->FindBin(0.); 
    for (unsigned int i = 0; i < mId.size(); ++i) {
      h0->GetXaxis()->SetBinLabel(center+1+i, Form("%d", mId[i]));
      h0->GetXaxis()->SetBinLabel(center-1-i, Form("%d", -mId[i]));
    }
    h0->SetLabelSize(0.03, "X");
    if (0) cout << h1 << endl;

    for (int i = 1; i < h0->GetNbinsX(); ++i) {
      cout << Form("%3d  ", i) 
	   << h0->GetXaxis()->GetBinLabel(i) << " -> " << atoi(h0->GetXaxis()->GetBinLabel(i))
	   << endl;
    }
  }    
  
  TGenCand *pCand; 
  cout << "======================================================================" << endl;
  cout << "gen block with " << fpEvt->nGenCands() << " gen cands" << endl;
  fpEvt->dumpGenBlock(); 
  cout << "======================================================================" << endl;
  cout << "Check for b and bBar" << endl;
  TGenCand *bMeson(0), *bBarMeson(0); 
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if (TMath::Abs(pCand->fID) == 5 && (pCand->fDau1 != pCand->fDau2)) {
      pCand->dump();
      TGenCand *pC; 
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pC = fpEvt->getGenCand(iD);
	if (isBeautyMeson(pC->fID) || isBeautyBaryon(pC->fID)) {
	  if (pC->fID > 0) {
	    bMeson = pC; 
	    cout << "bMeson    "; 
	  } else {
	    bBarMeson = pC; 
	    cout << "bBarMeson "; 
	  }
	} else {
	    cout << "          "; 
	}
	pC->dump();
      }
    }
    if (bMeson && bBarMeson) break;
  }

  TH1D *h0 = (TH1D*)fpHistFile->Get("muMom"); 
  TH1D *h1 = (TH1D*)fpHistFile->Get("muAnc"); 
  for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
    pCand = fpEvt->getGenCand(iC);
    if (13 == TMath::Abs(pCand->fID)) {
      int muMom(-1), muType(0); 
      muType = muonType(pCand); 
      h1->Fill(muType);
      muMom = (pCand->fMom1 > -1 && pCand->fMom1 < fpEvt->nGenCands() ? fpEvt->getGenCand(pCand->fMom1)->fID: -1);
      muonBinHist(muMom, h0);
      pCand->dump();
      cout << "MUON from bMeson: " << fpEvt->isAncestor(bMeson, pCand) 
	   << " or from bBarMeson: " << fpEvt->isAncestor(bBarMeson, pCand) 
	   << " muon type = " << muType
	   << endl;
    }  
  }    

  if (bMeson) {
    cout << "bMeson daughters" << endl;
    bMeson->dump(); 
    for (int iC = bMeson->fDau1; iC <= bMeson->fDau2; ++iC) {
      pCand = fpEvt->getGenCand(iC);
      pCand->dump(); 
//       if (isBeautyMesonWeak(pCand->fID)) {
// 	TGenCand *pC; 
// 	for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
// 	  pC = fpEvt->getGenCand(iD);
// 	  cout << " "; pC->dump(); 
// 	}
//       }
    }
  }

  if (bBarMeson) {
    cout << "bBarMeson daughters" << endl;
    bBarMeson->dump(); 
    for (int iC = bBarMeson->fDau1; iC <= bBarMeson->fDau2; ++iC) {
      pCand = fpEvt->getGenCand(iC);
      pCand->dump(); 
//       if (isBeautyMesonWeak(pCand->fID)) {
// 	TGenCand *pC; 
// 	for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
// 	  pC = fpEvt->getGenCand(iD);
// 	  cout << " "; pC->dump(); 
// 	}
//       }
    }
  }

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
  bool foundT(false), foundC(false), foundB(false), light(false), foundX(false); 

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

    // -- hit a b baryon
    if (ganz == 5) {
      foundB = 1; 
      break;
    }
 
    // -- hit a c baryon
    if (ganz == 4) {
      foundC = 1; 
    }
    
    // -- muon from tau decay
    if (rest == 15) {
      foundT = 1; 
    }

    if (15 < rest && rest < 400) {
      light = 1; 
    }

    if (rest > 399 && rest < 499) {
      foundC = 1; 
    }

    if (rest > 499 && rest < 599) {
      foundB = 1; 
      break;
    }

    if (rest > 1999 && rest < 5999) {
      foundX = 1; 
      //      break;
    }
  }

  result = 0; 

  if (foundB) result +=  1; 
  if (foundC) result +=  2; 
  if (foundT) result +=  4; 
  if (light)  result +=  8; 
  if (foundX) result += 16; 

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
  
  TH1D *h; 
  h = new TH1D("h1",   "mu type", 20, 0., 20.);
  h = new TH1D("h100", "mu pt 0", 100, 0., 20.);
  h = new TH1D("h110", "mu pt ==1", 100, 0., 20.);
  h = new TH1D("h112", "mu pt ==1", 100, 0., 20.);

  h = new TH1D("e1",  "evt type", 20, 0., 20.);
  h = new TH1D("e100", "evt pt 0", 100, 0., 20.);
  h = new TH1D("e110", "evt pt ==1", 100, 0., 20.);
  h = new TH1D("e111", "evt pt ==1", 100, 0., 20.);
  h = new TH1D("e112", "evt acc ==1", 100, 0., 20.);

  if (0) cout << h << endl;

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


// ---------------------------------------------------------------------------------------
void genAnalysis::muonBinHist(int id, TH1D *h) {
  //  double bin(0); 
  //   h->Fill(Form("%d", id), 1.); 
  //   if (id == -130) {
  //     cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX -130" << endl;
  //   }
  //   return;
  for (int i = 1; i < h->GetNbinsX(); ++i) {
    if (id == atoi(h->GetXaxis()->GetBinLabel(i))) {
      h->AddBinContent(i); 
      break;
    }
  }
}
