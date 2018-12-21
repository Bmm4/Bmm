#include "candAna.hh"

#include <TProfile.h>

#include "common/HFMasses.hh"
#include "common/AnalysisDistribution.hh"
#include "common/util.hh"

using namespace std;

struct near_track_t {
  int ix;
  float doca;
  float p;
  float pt;
  float pt_rel;
  float deltaR;
};

// ----------------------------------------------------------------------
candAna::candAna(recoilReader *pReader, string name, string cutsFile) {
  fAnaCuts.setAcName("candAna");
  fpReader = pReader;
  fVerbose = fpReader->fVerbose;
  fDbx     = -2;
  fYear    = fpReader->fYear;
  fEra     = fpReader->fEra;
  fName    = name;
  cout << "======================================================================" << endl;
  cout << "==> candAna: name = " << name << ", " << cutsFile << ", year " << fYear << ", era ->" << fEra << "<-" << endl;

  fGenBTmi = fNGenPhotons = fCandTmi = -1;

  fHistDir = gFile->mkdir(fName.c_str());

  cout << "======================================================================" << endl;

}


// ----------------------------------------------------------------------
candAna::~candAna() {
  cout << "==> candAna: destructor..." << endl;
}

// ----------------------------------------------------------------------
void candAna::endAnalysis() {
  cout << "This was for year " << fYear << " and era ->" << fEra << "<-" << endl;
  TH1D *h1 = ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())));
  if (h1) {
    cout << Form("==> mon%s: events seen    = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(h1->FindBin(1.)))) << endl;
    cout << Form("==> mon%s: cands analysed = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(h1->FindBin(10.)))) << endl;
    cout << Form("==> mon%s: cands passed   = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(h1->FindBin(11.)))) << endl;
    cout << Form("==> mon%s: cands failed   = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(h1->FindBin(20.)))) << endl;
    if (h1->GetBinContent(2) < 1) {
      cout << Form("==> mon%s: error, no events seen!", fName.c_str()) << endl;
    }
  } else {
    cout << Form("==> mon%s: error, histogram not found!", fName.c_str()) << endl;
  }

}


// ----------------------------------------------------------------------
void candAna::evtAnalysis(TAna01Event *evt) {
  cout << "  candAna: " << fChainEvent << endl;

  fpEvt = evt;

  if (0) {
    TAnaCand *pCand(0);
    TAnaTrack *pT(0);
    static int ngen(0), nrec(0);
    for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
      pCand = fpEvt->getCand(iC);
      if (-200 == pCand->fType) ++ngen;
      if (2000200 == pCand->fType) {
	++nrec;
	for (int iD = pCand->fSig1; iD <= pCand->fSig2; ++iD) {
	  pT = fpEvt->getSigTrack(iD);
	  pT->dump();
	}
	cout << Form("%4d", iC) << " cand -> " << pCand->fType << ", nrec/ngen = " << nrec << "/" << ngen << " = " << static_cast<double>(nrec)/ngen << endl;
      }
    }
    return;
  }


  bool fillNoCand(true);
  TAnaCand *pCand(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    if (fVerbose == -66) {
      cout << Form("%4d", iC) << " cand -> " << pCand->fType << endl;
      continue;
    }

    if (CANDTYPE != pCand->fType) {
      continue;
    } else {
      //      cout << "candAna: found cand with type " << CANDTYPE << endl;
    }
    // -- call derived functions (will jump back into candAna::candAnalysis for the common stuff!)
    candAnalysis();
    candEvaluation();

    if (fIsMC) {
      fTree->Fill();
    } else {  // DATA
      if (BLIND && (fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX)) {
	// -- do nothing
      } else {
	if (0 == NOPRESELECTION) {
	  // -- use fPreselection as filled in candEvaluation (below)
	} else if (1 == NOPRESELECTION) {
	  // -- original NOPRESELECTION,  this likely leads to very big reduced trees (too big for Chandi)
	  fPreselection = true;
	} else {
	  //  what here?
	}
	((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(32);

	if (fPreselection) {
	  if (fJSON) {
	    fTree->Fill();
	    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(11);
	  } else {
	    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(20);
	  }

	  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(12);
	  if (fJSON) ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(33);
	} else {
	  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(20);
	  if ( fVerbose > 9 ) cout << " failed preselection" << endl;

	} // if preselection
      } // if blind
    } // if MC
  }  // loop over cands

  // -- fill events with no passing candidate (one entry per event)
  if (fillNoCand) ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(2);

}

// ----------------------------------------------------------------------
void candAna::candAnalysis() {

  if (0 == fpCand) return;

  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(1);

}

// ----------------------------------------------------------------------
void candAna::candEvaluation() {
  fAnaCuts.update();

}


// ----------------------------------------------------------------------
void candAna::genMatch() {
  cout << "candAna::genMatch()  wrong function" << endl;
}


// ----------------------------------------------------------------------
void candAna::recoMatch() {
  cout << "candAna::recoMatch()  wrong function" << endl;
}


// ----------------------------------------------------------------------
void candAna::candMatch() {
  cout << "candAna::candMatch()  wrong function" << endl;
}


// ----------------------------------------------------------------------
void candAna::bookHist() {

  fHistDir->cd();

  TH1D *h11(0);
  (void)h11;
  h11 = new TH1D(Form("mon%s", fName.c_str()), Form("mon%s", fName.c_str()), 50, 0., 50.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupReducedTree(fTree);

  // -- Efficiency/Acceptance Tree
  fEffTree = new TTree("effTree", "effTree");
  fEffTree->Branch("run",    &fRun,               "run/L");
  fEffTree->Branch("evt",    &fEvt,               "evt/L");

  // -- Analysis distributions
  TH1D *h = new TH1D("analysisDistributions", "analysisDistributions", 10000, 0., 10000.);
  (void)h;

}


// ----------------------------------------------------------------------
void candAna::setupReducedTree(TTree *t) {

  t->Branch("run",     &fRun,               "run/L");
  t->Branch("evt",     &fEvt,               "evt/L");
  t->Branch("ls",      &fLS,                "ls/I");

}

// ----------------------------------------------------------------------
void candAna::readCuts(string fileName, int dump) {

  dump = 1;
  // -- set up cut sequence for analysis
  // basicCuts();
  // moreBasicCuts();
  // candidateCuts();
  // moreCandidateCuts();

  fCutFile = fileName;

  if (dump) cout << "==> candAna: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines;
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[1000]; //  XmlName[1000];
  fHistDir->cd();
  if (dump) {
    cout << "gDirectory: ";
    fHistDir->pwd();
  }
  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fCutFile.c_str());
  int ibin;
  string cstring = "B cand";

  int ok(0);
  float cutvalue;
  string cutname("nada");

  // -- determine fNchan
  cuts *a = 0;
  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    if (string::npos != cutLines[i].find("nchan")) {
      cleanupString(cutLines[i]);
      vector<string> lineItems = split(cutLines[i], ' ');
      fNchan = atoi(lineItems[1].c_str());
    }
  }

  if (fNchan < 1) {
    cout << "no analysis channels found?!" << endl;
  } else {
    cout << "creating " << fNchan << " analysis channels" << endl;
  }

  for (int i = 0; i < fNchan; ++i) {
    a = new cuts;
    a->index = i;
    fCuts.push_back(a);
  }


  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    // -- read the baseCuts file to get the channel definition and cuts used for (possible) preselection
    cleanupString(cutLines[i]);
    vector<string> lineItems = split(cutLines[i], ' ');
    if (lineItems.size() == 0) {
      continue;
    }
    cutname  = lineItems[0];
    ok = 0;
    if (static_cast<unsigned int>(fNchan + 1) == lineItems.size()) {
      for (unsigned int j = 1; j < lineItems.size(); ++j) {
	if (cutname == "index") {
	  // -- do nothing, indices already assigned above
	  ok = 1;
	}

	if (cutname == "mBdLo") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->mBdLo = cutvalue; ok = 1;
	}

	if (cutname == "mBdHi") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->mBdHi = cutvalue; ok = 1;
	}

	if (cutname == "mBuLo") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->mBuLo = cutvalue; ok = 1;
	}

	if (cutname == "mBuHi") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->mBuHi = cutvalue; ok = 1;
	}

	if (cutname == "mBsLo") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->mBsLo = cutvalue; ok = 1;
	}

	if (cutname == "mBsHi") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->mBsHi = cutvalue; ok = 1;
	}

	if (cutname == "metaMin") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->metaMin = cutvalue; ok = 1;
	}

	if (cutname == "metaMax") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->metaMax = cutvalue; ok = 1;
	}

	if (cutname == "muonbdt") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->muonbdt = cutvalue; ok = 1;
	}

	if (cutname == "m1pt") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->m1pt = cutvalue; ok = 1;
	}

	if (cutname == "m2pt") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->m2pt = cutvalue; ok = 1;
	}

	if (cutname == "etaMin") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->etaMin = cutvalue; ok = 1;
	}

	if (cutname == "etaMax") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->etaMax = cutvalue; ok = 1;
	}

	if (cutname == "phiMin") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->phiMin = cutvalue; ok = 1;
	}

	if (cutname == "phiMax") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->phiMax = cutvalue; ok = 1;
	}

	if (cutname == "pt") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pt = cutvalue; ok = 1;
	}

	if (cutname == "alpha") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->alpha = cutvalue; ok = 1;
	}

      	if (cutname == "fls3d") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->fls3d = cutvalue; ok = 1;
	}

	if (cutname == "flsxy") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->flsxy = cutvalue; ok = 1;
	}

	if (cutname == "flxyLo") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->flxyLo = cutvalue; ok = 1;
	}

	if (cutname == "flxyHi") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->flxyHi = cutvalue; ok = 1;
	}

	if (cutname == "chi2dof") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->chi2dof = cutvalue; ok = 1;
	}

	if (cutname == "iso") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->iso = cutvalue; ok = 1;
	}
	if (cutname == "m1iso") {
	  cutvalue = cutvalue; ok = 1;
	  fCuts[j-1]->m1iso = cutvalue; ok = 1;
	}
	if (cutname == "m2iso") {
	  cutvalue = cutvalue; ok = 1;
	  fCuts[j-1]->m2iso = cutvalue; ok = 1;
	}

	if (cutname == "docatrk") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->docatrk = cutvalue; ok = 1;
	}

	if (cutname == "closetrk") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->closetrk = cutvalue; ok = 1;
	}
	if (cutname == "closetrks1") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->closetrks1 = cutvalue; ok = 1;
	}
	if (cutname == "closetrks2") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->closetrks2 = cutvalue; ok = 1;
	}
	if (cutname == "closetrks3") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->closetrks3 = cutvalue; ok = 1;
	}

	if (cutname == "iso") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->iso = cutvalue; ok = 1;
	}

	if (cutname == "maxdoca") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->maxdoca = cutvalue; ok = 1;
	}

	if (cutname == "pvlip") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pvlip = cutvalue; ok = 1;
	}

	if (cutname == "pvlips") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pvlips = cutvalue; ok = 1;
	}

	if (cutname == "pvip") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pvip = cutvalue; ok = 1;
	}

	if (cutname == "pvips") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pvips = cutvalue; ok = 1;
	}

	if (cutname == "pv2lip") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pv2lip = cutvalue; ok = 1;
	}

	if (cutname == "pv2lips") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pv2lips = cutvalue; ok = 1;
	}

	if (0 == ok) {
	  cout << "XXXX unknown cut or cannot parse ->" << cutLines[i] << "<-" << endl;
	}

      }

    }

    // -- now back to the original cut reading
    sprintf(buffer, "%s", cutLines[i].c_str());

    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "CANDTYPE")) {
      CANDTYPE = int(CutValue);
      if (dump) cout << "TYPE:           " << CANDTYPE << endl;
      if (1313 == CANDTYPE) cstring = "#mu^{+}#mu^{-}";
      if (301313 == CANDTYPE) cstring = "#mu^{+}#mu^{-}";
      if (200521 == CANDTYPE) cstring = "#mu^{+}#mu^{-}K^{+}";
      if (300521 == CANDTYPE) cstring = "#mu^{+}#mu^{-}K^{+}";
      ibin = 1;
      hcuts->SetBinContent(ibin, CANDTYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
    }

    if (!strcmp(CutName, "CANDTRUTH")) {
      CANDTRUTH = int(CutValue);
      if (dump) cout << "CANDTRUTH:           " << CANDTRUTH << endl;
      ibin = 4;
      hcuts->SetBinContent(ibin, CANDTRUTH);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: CAND TRUTH", CutName));
    }

    if (!strcmp(CutName, "DATACAND")) {
      DATACAND = int(CutValue);
      if (dump) cout << "DATACAND:           " << DATACAND << endl;
      ibin = 4;
      hcuts->SetBinContent(ibin, DATACAND);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: DATA CAND", CutName));
    }


    if (!strcmp(CutName, "NOPRESELECTION")) {
      NOPRESELECTION = int(CutValue);
      if (dump) cout << "NOPRESELECTION     " << NOPRESELECTION << endl;
      ibin = 6;
      hcuts->SetBinContent(ibin, NOPRESELECTION);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Ignore preselection :: %i", CutName, NOPRESELECTION));
    }

    if (!strcmp(CutName, "SIGBOXMIN")) {
      SIGBOXMIN = CutValue;
      if (dump) cout << "SIGBOXMIN:           " << SIGBOXMIN << endl;
      ibin = 90;
      hcuts->SetBinContent(ibin, SIGBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: SIGBOXMIN :: %6.3f", CutName, SIGBOXMIN));
    }

    if (!strcmp(CutName, "SIGBOXMAX")) {
      SIGBOXMAX = CutValue;
      if (dump) cout << "SIGBOXMAX:           " << SIGBOXMAX << endl;
      ibin = 91;
      hcuts->SetBinContent(ibin, SIGBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: SIGBOXMAX :: %6.3f", CutName, SIGBOXMAX));
    }

    if (!strcmp(CutName, "BGLBOXMIN")) {
      BGLBOXMIN = CutValue;
      if (dump) cout << "BGLBOXMIN:           " << BGLBOXMIN << endl;
      ibin = 92;
      hcuts->SetBinContent(ibin, BGLBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGLBOXMIN :: %6.3f", CutName, BGLBOXMIN));
    }

    if (!strcmp(CutName, "BGLBOXMAX")) {
      BGLBOXMAX = CutValue;
      if (dump) cout << "BGLBOXMAX:           " << BGLBOXMAX << endl;
      ibin = 93;
      hcuts->SetBinContent(ibin, BGLBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGLBOXMAX :: %6.3f", CutName, BGLBOXMAX));
    }

    if (!strcmp(CutName, "BGHBOXMIN")) {
      BGHBOXMIN = CutValue;
      if (dump) cout << "BGHBOXMIN:           " << BGHBOXMIN << endl;
      ibin = 94;
      hcuts->SetBinContent(ibin, BGHBOXMIN);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGHBOXMIN :: %6.3f", CutName, BGHBOXMIN));
    }

    if (!strcmp(CutName, "BGHBOXMAX")) {
      BGHBOXMAX = CutValue;
      if (dump) cout << "BGHBOXMAX:           " << BGHBOXMAX << endl;
      ibin = 95;
      hcuts->SetBinContent(ibin, BGHBOXMAX);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BGHBOXMAX :: %6.3f", CutName, BGHBOXMAX));
    }

    // -- Tracks
    if (!strcmp(CutName, "TRACKPTLO")) {
      TRACKPTLO = CutValue;
      if (dump) cout << "TRACKPTLO:           " << TRACKPTLO << " GeV" << endl;
      ibin = 101;
      hcuts->SetBinContent(ibin, TRACKPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(track) :: %3.1f", CutName, TRACKPTLO));
    }

    if (!strcmp(CutName, "TRACKPTHI")) {
      TRACKPTHI = CutValue;
      if (dump) cout << "TRACKPTHI:           " << TRACKPTHI << " GeV" << endl;
      ibin = 102;
      hcuts->SetBinContent(ibin, TRACKPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{max}(track) :: %3.1f", CutName, TRACKPTHI));
    }

    if (!strcmp(CutName, "TRACKTIP")) {
      TRACKTIP = CutValue;
      if (dump) cout << "TRACKTIP:           " << TRACKTIP << " cm" << endl;
      ibin = 103;
      hcuts->SetBinContent(ibin, TRACKTIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: doca_{xy}(track) :: %3.1f", CutName, TRACKTIP));
    }

    if (!strcmp(CutName, "TRACKLIP")) {
      TRACKLIP = CutValue;
      if (dump) cout << "TRACKLIP:           " << TRACKLIP << " cm" << endl;
      ibin = 104;
      hcuts->SetBinContent(ibin, TRACKLIP);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: doca_{z}(track) :: %3.1f", CutName, TRACKLIP));
    }

    if (!strcmp(CutName, "TRACKETALO")) {
      TRACKETALO = CutValue;
      if (dump) cout << "TRACKETALO:           " << TRACKETALO << " " << endl;
      ibin = 105;
      hcuts->SetBinContent(ibin, TRACKETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta_{min}(track) :: %3.1f", CutName, TRACKETALO));
    }

    if (!strcmp(CutName, "TRACKETAHI")) {
      TRACKETAHI = CutValue;
      if (dump) cout << "TRACKETAHI:           " << TRACKETAHI << " " << endl;
      ibin = 106;
      hcuts->SetBinContent(ibin, TRACKETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta_{max}(track) :: %3.1f", CutName, TRACKETAHI));
    }

    if (!strcmp(CutName, "MUBDTXML")) {
      char xml[1000];
      sscanf(buffer, "%s %s", CutName, xml);
      string tl(xml);
      if (dump) {
	cout << "MUBDTXML:       " << xml << endl;
      }
      ibin = 207;
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, xml));
      fMvaMuonID = setupMuonMvaReader(string(xml), mrd);
    }

    if (!strcmp(CutName, "MUBDT")) {
      MUBDT = CutValue;
      if (dump) cout << "MUBDT:           " << MUBDT << endl;
      ibin = 208;
      hcuts->SetBinContent(ibin, MUBDT);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: BDT(#mu) :: %3.1f", CutName, MUBDT));
    }

  }

  printCuts(cout);
}

// ----------------------------------------------------------------------
void candAna::readFile(string filename, vector<string> &lines) {
  cout << "    readFile " << filename << endl;
  char  buffer[200];
  ifstream is(filename.c_str());
  if (!is) {
    exit(1);
  }
  string sbuffer;
  while (is.getline(buffer, 200, '\n')) {
    sbuffer = string(buffer);
    if (buffer[0] != '+') {
      lines.push_back(sbuffer);
    } else {
      if (string::npos != sbuffer.find("YEARERA")) {
	replaceAll(sbuffer, "YEARERA", Form("%d%s", fYear, fEra.c_str()));
      }
      if (string::npos != sbuffer.find("YEAR")) {
	replaceAll(sbuffer, "YEAR", Form("%d", fYear));
      }
      replaceAll(sbuffer, "+input ", "");
      readFile(sbuffer, lines);
    }
  }

}



// ----------------------------------------------------------------------
void candAna::play() {

}


// ----------------------------------------------------------------------
void candAna::printCuts(ostream &OUT) {

  OUT << "----------------------------------------------------------------------" << endl;
  OUT << "channel    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  OUT << Form("%10d", fCuts[i]->index);
  OUT << endl;

  OUT << "metaMin    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->metaMin);
  }
  OUT << endl;

  OUT << "metaMax    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->metaMax);
  }
  OUT << endl;

  OUT << "muonBDT    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->muonbdt);
  }
  OUT << " (ignored) " << endl;

  OUT << "l1seeds     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    for (int is = fCuts[i]->l1seeds.size()*2; is < 10; ++is) OUT << " ";
    for (unsigned is = 0; is < fCuts[i]->l1seeds.size(); ++is) OUT << Form("%d ", fCuts[i]->l1seeds[is]);
  }
  OUT << endl;

  OUT << "mBdLo      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    OUT << Form("%10.3f", fCuts[i]->mBdLo);
  }
  OUT << endl;

  OUT << "mBdHi      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    OUT << Form("%10.3f", fCuts[i]->mBdHi);
  }
  OUT << endl;

  OUT << "mBsLo      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->mBsLo);
  }
  OUT << endl;

  OUT << "mBsHi      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->mBsHi);
  }
  OUT << endl;

  OUT << "mBuLo      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->mBuLo);
  }
  OUT << endl;

  OUT << "mBuHi      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->mBuHi);
  }
  OUT << endl;


  OUT << "m1pt       ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->m1pt);
  }
  OUT << endl;

  OUT << "m2pt       ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->m2pt);
  }
  OUT << endl;


  OUT << "etaMin     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->etaMin);
  }
  OUT << endl;

  OUT << "etaMax     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->etaMax);
  }
  OUT << endl;

  OUT << "pt         ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pt);
  }
  OUT << endl;

  OUT << "iso        ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->iso);
  }
  OUT << endl;

  OUT << "m1iso      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->m1iso);
  }
  OUT << endl;

  OUT << "m2iso      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->m2iso);
  }
  OUT << endl;

  OUT << "chi2dof    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->chi2dof);
  }
  OUT << endl;

  OUT << "alpha      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->alpha);
  }
  OUT << endl;

  OUT << "fls3d      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->fls3d);
  }
  OUT << endl;

  OUT << "flsxy      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->flsxy);
  }
  OUT << endl;

  OUT << "docatrk    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->docatrk);
  }
  OUT << endl;

  OUT << "closetrk   ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->closetrk);
  }
  OUT << endl;

  OUT << "closetrks1 ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->closetrks1);
  }
  OUT << endl;

  OUT << "closetrks2 ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->closetrks2);
  }
  OUT << endl;

  OUT << "closetrks3 ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->closetrks3);
  }
  OUT << endl;

  OUT << "maxdoca    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->maxdoca);
  }
  OUT << endl;

  OUT << "pvip       ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pvip);
  }
  OUT << endl;

  OUT << "pvips      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pvips);
  }
  OUT << endl;

  OUT << "pvlip      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pvlip);
  }
  OUT << endl;

  OUT << "pvlips     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pvlips);
  }
  OUT << endl;

  OUT << "pv2lip     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pv2lip);
  }
  OUT << endl;

  OUT << "pv2lips    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pv2lips);
  }
  OUT << endl;

  OUT << "bdtXml     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10s", fCuts[i]->bdtXml.c_str());
  }
  OUT << endl;

  OUT << "bdtCut     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->bdtCut);
  }
  OUT << endl;

  OUT << "bdtMuPt    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->bdtMuPt);
  }
  OUT << endl;


  OUT << "----------------------------------------------------------------------" << endl;

  OUT.flush();

  return;
}


// ----------------------------------------------------------------------
TMVA::Reader* candAna::setupMuonMvaReader(string xmlFile, mvaMuonIDData &d) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  TString dir    = "weights/";
  TString methodNameprefix = "BDT";
  // -- read in variables from weight file
  vector<string> allLines;
  char  buffer[2000];
  string weightFile = "weights/TMVA-" + xmlFile + ".weights.xml";
  cout << "setupMuonMvaReader, open file " << weightFile << endl;
  ifstream is(weightFile.c_str());
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  int nvars(-1);
  string::size_type m1, m2;
  string stype;
  cout << "  read " << allLines.size() << " lines " << endl;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\"");
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2);
      cout << "  " << stype << " variables" << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
        m1 = allLines[j].find("Expression=\"")+10;
        m2 = allLines[j].find("\" Label=\"");
        stype = allLines[j].substr(m1+2, m2-m1-2);
        if (stype == "segComp") {
          //      cout << "  adding segComp" << endl;
          reader->AddVariable("segComp", &d.segComp);
          continue;
        }
        if (stype == "chi2LocMom") {
          //      cout << "  adding chi2LocMom" << endl;
          reader->AddVariable("chi2LocMom", &d.chi2LocMom);
          continue;
        }
        if (stype == "chi2LocPos") {
          //      cout << "  adding chi2LocPos" << endl;
          reader->AddVariable("chi2LocPos", &d.chi2LocPos);
          continue;
        }
        if (stype == "NTrkVHits") {
          //      cout << "  adding NTrkVHits" << endl;
          reader->AddVariable("NTrkVHits", &d.NTrkVHits);
          continue;
        }

        if (stype == "trkValidFract") {
          //      cout << "  adding trkValidFract" << endl;
          reader->AddVariable( "trkValidFract", &d.trkValidFract);
          continue;
        }
        if (stype == "glbNChi2") {
          reader->AddVariable("glbNChi2", &d.glbNChi2);
          //      cout << "  adding glbNChi2" << endl;
          continue;
        }
        if (stype == "eta") {
          //      cout << "  adding eta" << endl;
          reader->AddVariable("eta", &d.eta);
          continue;
        }
        if (stype == "pt") {
          //      cout << "  adding pt" << endl;
          reader->AddVariable("pt", &d.pt);
          continue;
        }
        if (stype == "glbTrackProb") {
          //      cout << "  adding glbTrackProb" << endl;
          reader->AddVariable("glbTrackProb", &d.glbTrackProb);
          continue;
        }
        if (stype == "NTrkEHitsOut") {
          //      cout << "  adding NTrkEHitsOut" << endl;
          reader->AddVariable("NTrkEHitsOut", &d.NTrkEHitsOut);
          continue;
        }

        // -- variables for SW's BDT
        if (stype == "glbTrackTailProb") {
          //      cout << "  adding glbTrackTailProb" << endl;
          reader->AddVariable("glbTrackTailProb", &d.glbTrackTailProb);
          continue;
        }
        if (stype == "glbDeltaEtaPhi") {
          //      cout << "  adding glbDeltaEtaPhi" << endl;
          reader->AddVariable("glbDeltaEtaPhi", &d.glbDeltaEtaPhi);
          continue;
        }
        if (stype == "iValFrac") {
          //      cout << "  adding iValFrac (aka trkValidFract)" << endl;
          reader->AddVariable("iValFrac", &d.trkValidFract);
          continue;
        }
        if (stype == "LWH") {
          //      cout << "  adding LWH" << endl;
          reader->AddVariable("LWH", &d.LWH);
          continue;
        }
        if (stype == "dxyRef") {
          //      cout << "  adding dxyRef" << endl;
          reader->AddVariable("dxyRef", &d.dxyRef);
          continue;
        }
        if (stype == "kinkFinder") {
          //      cout << "  adding kinkFinder" << endl;
          reader->AddVariable("kinkFinder", &d.kinkFinder);
          continue;
        }
        if (stype == "dzRef") {
          //      cout << "  adding dzRef" << endl;
          reader->AddVariable("dzRef", &d.dzRef);
          continue;
        }
        if (stype == "glbKinkFinder") {
          //      cout << "  adding glbKinkFinder" << endl;
          reader->AddVariable("glbKinkFinder", &d.glbKinkFinder);
          continue;
        }
        if (stype == "TMath::Log(2+glbKinkFinder)") {
          //      cout << "  adding TMath::Log(2+glbKinkFinder)" << endl;
          reader->AddVariable("TMath::Log(2+glbKinkFinder)", &d.glbKinkFinderLOG);
          continue;
        }
        if (stype == "timeAtIpInOutErr") {
          //      cout << "  adding timeAtIpInOutErr" << endl;
          reader->AddVariable("timeAtIpInOutErr", &d.timeAtIpInOutErr);
          continue;
        }
        if (stype == "outerChi2") {
          //      cout << "  adding outerChi2" << endl;
          reader->AddVariable("outerChi2", &d.outerChi2);
          continue;
        }
        if (stype == "valPixHits") {
          //      cout << "  adding valPixHits" << endl;
          reader->AddVariable("valPixHits", &d.valPixHits);
          continue;
        }
        if (stype == "TMTrkMult100") {
          //      cout << "  adding TMTrkMult100" << endl;
          reader->AddVariable("TMTrkMult100", &d.TMTrkMult100);
          continue;
        }
        if (stype == "innerChi2") {
          //      cout << "  adding innerChi2" << endl;
          reader->AddVariable("innerChi2", &d.innerChi2);
          continue;
        }
        if (stype == "trkRelChi2") {
          //      cout << "  adding trkRelChi2" << endl;
          reader->AddVariable("trkRelChi2", &d.trkRelChi2);
          continue;
        }
        if (stype == "vMuonHitComb") {
          //      cout << "  adding vMuonHitComb" << endl;
          reader->AddVariable("vMuonHitComb", &d.vMuonHitComb);
          continue;
        }
        if (stype == "Qprod") {
          //      cout << "  adding Qprod" << endl;
          reader->AddVariable("Qprod", &d.Qprod);
          continue;
        }
      }
      break;
    }
  }

  nvars = -1;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add spectators
    if (string::npos != allLines[i].find("Spectators NSpec")) {
      m1 = allLines[i].find("=\"");
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2);
      //      cout << "==> " << stype << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
        m1 = allLines[j].find("Expression=\"")+10;
        m2 = allLines[j].find("\" Label=\"");
        stype = allLines[j].substr(m1+2, m2-m1-2);
        //      cout << "ivar " << j-i << " spectator string: ->" << stype << "<-" << endl;
        //      cout << " adding " << stype << " as a spectator dummy" << endl;
        reader->AddSpectator( stype.c_str(), &d.spectatorDummy);
      }
      break;
    }
  }

  reader->BookMVA("BDT", TString(weightFile.c_str()));
  return reader;
}
