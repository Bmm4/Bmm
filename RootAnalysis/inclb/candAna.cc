#include "candAna.hh"

#include "common/HFMasses.hh"
//#include "fastjet/ClusterSequence.hh"

#include "common/ana.hh"
#include "common/util.hh"

using namespace std;
//using namespace fastjet;

// ----------------------------------------------------------------------
candAna::candAna(inclbReader *pReader, string name, string cutsFile) {
  fpReader = pReader; 
  fVerbose = fpReader->fVerbose;	 
  fYear    = fpReader->fYear; 
  fName    = name; 
  cout << "======================================================================" << endl;
  cout << "==> candAna: name = " << name << ", reading cutsfile " << cutsFile << " setup for year " << fYear << " verbose = " << fVerbose << endl;
  fHistDir = gFile->mkdir(fName.c_str());
  readCuts(cutsFile, 1); 
  cout << "======================================================================" << endl;	 

  // // -- test FASTJET
  //   vector<PseudoJet> particles;
  //   // an event with three particles:   px    py  pz      E
  //   particles.push_back( PseudoJet(   99.0,  0.1,  0, 100.0) ); 
  //   particles.push_back( PseudoJet(    4.0, -0.1,  0,   5.0) ); 
  //   particles.push_back( PseudoJet(  -99.0,    0,  0,  99.0) );
  
  //   // choose a jet definition
  //   double R = 0.7;
  //   JetDefinition jet_def(antikt_algorithm, R);
  
  //   // run the clustering, extract the jets
  //   ClusterSequence cs(particles, jet_def);
  //   vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
  
  //   // print out some infos
  //   cout << "Clustering with " << jet_def.description() << endl;
  
  //   // print the jets
  //   cout <<   "        pt y phi" << endl;
  //   for (unsigned i = 0; i < jets.size(); i++) {
  //     cout << "jet " << i << ": "<< jets[i].pt() << " " 
  // 	 << jets[i].rap() << " " << jets[i].phi() << endl;
  //     vector<PseudoJet> constituents = jets[i].constituents();
  //     for (unsigned j = 0; j < constituents.size(); j++) {
  //       cout << "    constituent " << j << "'s pt: " << constituents[j].pt()
  //            << endl;
  //     }
  //   }
  
}


// ----------------------------------------------------------------------
candAna::~candAna() {
  cout << "==> candAna: destructor..." << endl;
}

// ----------------------------------------------------------------------
void candAna::endAnalysis() {
  TH1D *h1 = ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())));
  if (h1) {
    cout << Form("==> mon%s: events seen    = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(2))) << endl;
    cout << Form("==> mon%s: cands analysed = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(11))) << endl;
    cout << Form("==> mon%s: cands passed   = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(12))) << endl;
    cout << Form("==> mon%s: cands failed   = %d", fName.c_str(), static_cast<int>(h1->GetBinContent(21))) << endl;
    if (h1->GetBinContent(2) < 1) {
      cout << Form("==> mon%s: error, no events seen!", fName.c_str()) << endl; 
    }
  } else {
    cout << Form("==> mon%s: error, histogram not found!", fName.c_str()) << endl; 
  }    
}



// ----------------------------------------------------------------------
void candAna::evtAnalysis(TAna01Event *evt) {
  
  fpEvt = evt; 

  TAnaTrack *pSigTrack(0);
  if (fVerbose > 39) { 
    cout << "---- Evt: " << fEvt << " n(sig tracks) = " << fpEvt->nSigTracks() << endl;
  }
  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(1);
  fpSigTrack = 0; 
  fpSigJet = 0; 
  for (int iC = 0; iC < fpEvt->nSigTracks(); ++iC) {
    pSigTrack = fpEvt->getSigTrack(iC);
    
    if (TYPE != pSigTrack->fInt1) {
      if (fVerbose > 39) cout << "  skipping sig track at " << iC << " which is of type " << pSigTrack->fInt1 << endl;
      continue;
    } else {
      if (fVerbose > 39) cout << "  analyzing sig track at " << iC << " which is of type " << pSigTrack->fInt1 
			      << " with ptrel = " << pSigTrack->fDouble1
			      << " jet index = " << pSigTrack->fInt2
			      << endl;
    }
    
    if (pSigTrack->fInt2 < 0) continue;

    fpSigTrack = pSigTrack; 
    fpSigJet   = fpEvt->getTrackJet(pSigTrack->fInt2);

    fpMuon = fpEvt->getSimpleTrackMuon(fpSigTrack->fIndex); 
    fpTrack = fpEvt->getSimpleTrack(fpSigTrack->fIndex); 

    if (fVerbose > 99) {
      TVector3 vj = fpSigJet->fPlab - fpSigTrack->fPlab; 
      vj.SetMag(1.); 
      TVector3 vm = fpSigTrack->fPlab; 
      double ptrel = vj.Cross(vm).Mag();
      if (1) cout << "LS " << fLS << " sigTrack " << iC << " type " << fpSigTrack->fInt1 << " pt " << fpSigTrack->fPlab.Perp()
		  << " jet idx: " << pSigTrack->fInt2
		  << " jet pt/eta/phi = " << fpSigJet->fPlab.Perp() << "/" << fpSigJet->fPlab.Eta() << "/" << fpSigJet->fPlab.Phi() 
		  << " ntrk = " << fpSigJet->getNtracks()
		  << " ptrel = " << fpSigTrack->fDouble1 << " == " << ptrel
		  << endl;
    }

    genAnalysis();
    candAnalysis();
    triggerSelection();
    fillRedTreeData();

    bool doFill(true);
    if ((3 == fIsMC) && (0 < fProcessType)) {
      doFill = false;
    }
    if ((4 == fIsMC) && ((499 < fProcessType && fProcessType < 600) || (0 == fProcessType))) {
      doFill = false;
    }
    if ((5 == fIsMC) && ((399 < fProcessType && fProcessType < 500) || (0 == fProcessType))) {
      doFill = false;
    }
    if (doFill) {
      fTree->Fill(); 
    }
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(11);
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(31);
  } 

}


// ----------------------------------------------------------------------
void candAna::genAnalysis() {

  fProcessType = 0; 

  if (0 == fIsMC) {
    //    cout << "not MC" << endl;
    return;
  }

  if (fpTrack->getGenIndex() < 0) {
    //    cout << "no matching gen cand found" << endl;
    return;
  }

  //  fpEvt->dumpGenBlock();

  fpGenMuon = fpEvt->getGenCand(fpTrack->getGenIndex());
  TGenCand *pM(fpGenMuon); 
  int iMom(2), aid(0); 
  bool beauty(false), 
    charm(false), 
    light(false), 
    ccbar(false), 
    tau(false), 
    fake(false); 
  
  int cnt(-1); 
  while (iMom > 1 && iMom < fpEvt->nGenCands()) {
    ++cnt;
    //    cout << "cnt = " << cnt << " iMom = " << iMom << " pM = " << pM << endl;
    iMom = pM->fMom1;
    pM   = fpEvt->getGenCand(iMom);

    aid = TMath::Abs(pM->fID); 
    if ((211 == aid) 
	|| (321 == aid) 
	|| (2212 == aid) 
	) fake = true;

    if ((111 == aid)        // pi0
	|| (221 == aid)     // eta
	|| (113 == aid)     // rho0
	|| (213 == aid)     // rho+
	|| (223 == aid)     // omega
	|| (333 == aid)     // phi
	) light = true; 

    if ((441 == aid )       // eta_c
	|| (443 == aid)     // J/psi
	|| (100443 == aid)  // J/psi
	) ccbar = true; 

    if (15 == aid)          // tau
      tau = true;

    if (isCharmMesonWeak(aid)) charm = true; 

    if (isBeautyMesonWeak(aid)) {
      beauty = true; 
      break;
    }

    if (21 == aid) break;
    if (5 == aid) break;
    if (4 == aid) break;
    if (3 == aid) break;
    if (2 == aid) break;
    if (1 == aid) break;
  }

  if (beauty) {
    if (charm) {
      if (tau) {
	fProcessType = 550; 
      } else if (light) {
	fProcessType = 560; 
      } else {
	fProcessType = 510; 
      }
    } else if (ccbar) {
      fProcessType = 530; 
    } else if (tau) {
      fProcessType = 540; 
    } else if (light) {
      fProcessType = 520; 
    } else if (fake) {
      fProcessType = 590; 
    } else {
      fProcessType = 500; 
    }
  } else if (charm) {
    if (tau) {
      fProcessType = 440; 
    } else if (light) {
      fProcessType = 420; 
    } else if (fake) {
      fProcessType = 490; 
    } else {
      fProcessType = 400; 
    }
  }


  //  cout << "==> fProcessType = " << fProcessType 
  //       << endl;
}


// ----------------------------------------------------------------------
void candAna::candAnalysis() {

  if (!fpSigTrack) return;
  if (!fpTrack) return;
  if (!fpMuon) return;
  if (!fpSigJet) return;


  fMuId    = tightMuon(fpMuon); 

  fType = fpSigTrack->fInt1; 

  fMuPtRel = fpSigTrack->fDouble1; 

  fMuPt    = fpSigTrack->fPlab.Perp(); 
  fMuEta   = fpSigTrack->fPlab.Eta(); 
  fMuPhi   = fpSigTrack->fPlab.Phi(); 

  fJetPt   = fpSigJet->fPlab.Perp(); 
  fJetEta  = fpSigJet->fPlab.Eta(); 
  fJetPhi  = fpSigJet->fPlab.Phi(); 

}


// ----------------------------------------------------------------------
void candAna::bookHist() {

  fHistDir->cd();

  TH1D *h11(0); 
  (void)h11; 
  h11 = new TH1D(Form("mon%s", fName.c_str()), Form("mon%s", fName.c_str()), 50, 0., 50.); 


  TH2D *h22(0); 
  (void)h22; 

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupRedTree(fTree); 
  
  // -- Analysis distributions
  TH1D *h = new TH1D("analysisDistributions", "analysisDistributions", 10000, 0., 10000.); 
  (void)h;

}


// ----------------------------------------------------------------------
void candAna::setupRedTree(TTree *t) {

  t->Branch("run",     &fRTD.run,          "run/L");
  t->Branch("evt",     &fRTD.evt,          "evt/L");
  t->Branch("ls",      &fRTD.ls,           "ls/I");

  t->Branch("muid",    &fRTD.muid,         "muid/O");
  t->Branch("hlt",     &fRTD.hlt,          "hlt/O");
  t->Branch("hltmatch",&fRTD.hltmatch,     "hltmatch/O");
  t->Branch("json",    &fRTD.json,         "json/O");

  t->Branch("type",    &fRTD.type,         "type/I");

  t->Branch("pt",      &fRTD.pt,           "pt/F");
  t->Branch("eta",     &fRTD.eta,          "eta/F");
  t->Branch("phi",     &fRTD.phi,          "phi/F");
  t->Branch("ptrel",   &fRTD.ptrel,        "ptrel/F");

  t->Branch("jpt",     &fRTD.jpt,          "jpt/F");
  t->Branch("jeta",    &fRTD.jeta,         "jeta/F");
  t->Branch("jphi",    &fRTD.jphi,         "jphi/F");

}


// ----------------------------------------------------------------------
void candAna::readCuts(string fileName, int dump) {

  // -- define default values for some cuts
  NOPRESELECTION = 0; 
  IGNORETRIGGER  = 0; 

  // -- set up cut sequence for analysis
  fCutFile = fileName;

  if (dump) cout << "==> candAna: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines; 
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[1000];
  fHistDir->cd();
  if (dump) cout << "gDirectory: "; fHistDir->pwd();
  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fCutFile.c_str());
  int ibin; 
  string cstring = "B cand"; 

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); 
      if (dump) cout << "TYPE:           " << TYPE << endl;
      if (100100 == TYPE) cstring = "TrackJets";
      ibin = 1;
      hcuts->SetBinContent(ibin, TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
    }

    if (!strcmp(CutName, "TRIGRANGE")) {
      char triggerlist[1000]; 
      sscanf(buffer, "%s %s", CutName, triggerlist);
      string tl(triggerlist); 
      int r1(0), r2(0); 
      string hlt = splitTrigRange(tl, r1, r2); 
      HLTRANGE.insert(make_pair(hlt, make_pair(r1, r2))); 
      if (dump) {
	cout << "HLTRANGE:       " << hlt << " from " << r1 << " to " << r2 << endl; 
      }
      ibin = 3; 
      hcuts->SetBinContent(ibin, 1);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: %s", CutName, triggerlist));
    }


    if (!strcmp(CutName, "IGNORETRIGGER")) {
      IGNORETRIGGER = int(CutValue); 
      if (dump) cout << "IGNORETRIGGER      " << IGNORETRIGGER << endl;
      ibin = 5;
      hcuts->SetBinContent(ibin, IGNORETRIGGER);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Ignore trigger :: %i", CutName, IGNORETRIGGER));
    }

    if (!strcmp(CutName, "NOPRESELECTION")) {
      NOPRESELECTION = int(CutValue); 
      if (dump) cout << "NOPRESELECTION     " << NOPRESELECTION << endl;
      ibin = 6;
      hcuts->SetBinContent(ibin, NOPRESELECTION);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Ignore preselection :: %i", CutName, NOPRESELECTION));
    }


    if (!strcmp(CutName, "MUPTLO")) {
      MUPTLO = CutValue; 
      if (dump) cout << "MUPTLO:           " << MUPTLO << " GeV" << endl;
      ibin = 202;
      hcuts->SetBinContent(ibin, MUPTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{min}(#mu) :: %3.1f", CutName, MUPTLO));
    }

    if (!strcmp(CutName, "MUPTHI")) {
      MUPTHI = CutValue; 
      if (dump) cout << "MUPTHI:           " << MUPTHI << " GeV" << endl;
      ibin = 203;
      hcuts->SetBinContent(ibin, MUPTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: p_{T}^{max}(#mu) :: %3.1f", CutName, MUPTHI));
    }

    if (!strcmp(CutName, "MUETALO")) {
      MUETALO = CutValue; 
      if (dump) cout << "MUETALO:           " << MUETALO << endl;
      ibin = 204;
      hcuts->SetBinContent(ibin, MUETALO);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{min}(#mu) :: %3.1f", CutName, MUETALO));
    }

    if (!strcmp(CutName, "MUETAHI")) {
      MUETAHI = CutValue; 
      if (dump) cout << "MUETAHI:           " << MUETAHI << endl;
      ibin = 205;
      hcuts->SetBinContent(ibin, MUETAHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: #eta^{max}(#mu) :: %3.1f", CutName, MUETAHI));
    }

  }

  if (dump)  cout << "------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void candAna::readFile(string filename, vector<string> &lines) {
  cout << "    readFile " << filename << endl;
  char  buffer[200];
  ifstream is(filename.c_str());
  if (!is) {
    exit(1);
  }
  char input[1000]; 
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] != '+') {
      lines.push_back(string(buffer));
    } else {
      sscanf(buffer, "+input %s", input);
      readFile(input, lines); 
    }
  }

}


// ----------------------------------------------------------------------
void candAna::getSigTracks(vector<int> &v, TAnaCand *pC) {
  TAnaCand *pD; 
  TAnaTrack *pT; 
  vector<int> bla; 

  // -- loop over daughters
  if (pC->fDau1 > -1) {
    for (int j = pC->fDau1; j <= pC->fDau2; ++j) {
      pD = fpEvt->getCand(j); 
      getSigTracks(bla, pD); 
    }

    for (unsigned j = 0; j < bla.size(); ++j) v.push_back(bla[j]);
  }

  // -- add direct sigtracks
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getSigTrack(i); 
    if (v.end() == find(v.begin(), v.end(), pT->fIndex)) {
      v.push_back(pT->fIndex); 
    }
  }

}


// ----------------------------------------------------------------------
void candAna::fillRedTreeData() {
  fRTD.run       = fRun;
  fRTD.evt       = fEvt;
  fRTD.ls        = fLS;
  fRTD.type      = fProcessType;

  fRTD.muid      = fMuId;
  fRTD.hlt       = fGoodHLT;
  fRTD.hltmatch  = fHLTmatch;
  fRTD.json      = fJSON;

  fRTD.pt        = fMuPt; 
  fRTD.eta       = fMuEta; 
  fRTD.phi       = fMuPhi; 
  fRTD.ptrel     = fMuPtRel; 

  fRTD.jpt        = fJetPt; 
  fRTD.jeta       = fJetEta; 
  fRTD.jphi       = fJetPhi; 
}


// ----------------------------------------------------------------------
string candAna::splitTrigRange(string tl, int &r1, int &r2) {

  string::size_type id1 = tl.find_first_of("("); 
  string::size_type id2 = tl.find_first_of(":"); 
  string::size_type id3 = tl.find_first_of(")"); 

  //cout << "tl: " << tl << endl;
  string hlt = tl.substr(0, id1);
  //cout << "hlt: " << hlt << endl;
  string a   = tl.substr(id1+1, id2-id1-1);
  r1 = atoi(a.c_str());
  //cout << "1st a: " << a << " -> r1 = " << r1 << endl;
  a  = tl.substr(id2+1, id3-id2-1); 
  r2 = atoi(a.c_str());
  //cout << "2nd a: " << a << " -> r2 = " << r2 << endl;

  return hlt; 

}


// ----------------------------------------------------------------------
void candAna::triggerSelection() {

  fGoodHLT = false; 

  TString a; 
  int ps(0); 
  bool result(false), wasRun(false), error(false); 

  if (0 && HLTRANGE.begin()->first == "NOTRIGGER" ) { 
    cout << "NOTRIGGER requested... " << endl;
    fGoodHLT = true; 
    return;
  }

  if (fVerbose == -31) {
    cout << "--------------------  L1" << endl;
    for (int i = 0; i < NL1T; ++i) {
      result = wasRun = error = false;
      a = fpEvt->fL1TNames[i]; 
      ps = fpEvt->fL1TPrescale[i]; 
      result = fpEvt->fL1TResult[i]; 
      error  = fpEvt->fL1TMask[i]; 
      //if (a.Contains("Mu")) {
      if (result ) {
	cout << a <<  " mask: " << error << " result: " << result << " ps: " << ps << endl;
      }
    }
  }
  
  if (fVerbose == -32) cout<<" event " << fEvt << endl;

  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i]; 
    ps = fpEvt->fHLTPrescale[i]; 
    wasRun = fpEvt->fHLTWasRun[i]; 
    result = fpEvt->fHLTResult[i]; 
    error  = fpEvt->fHLTError[i]; 

    if (-31 == fVerbose && wasRun > 0) cout << a 
					    << " was run: " << wasRun 
					    << " result: " << result 
					    << " ps: " << ps 
					    << endl;

    string spath; 
    for (map<string, pair<int, int> >::iterator imap = HLTRANGE.begin(); imap != HLTRANGE.end(); ++imap) {  
      spath = imap->first; 
      if (!a.CompareTo(imap->first.c_str())) {
	fGoodHLT = result; 
	if (fVerbose > 1 || -32 == fVerbose  ) cout << "exact match: " << imap->first.c_str() << " HLT: " << a << " result: " << result << endl;
      }
    }
  }      

  if (false == fGoodHLT) {
    if (fVerbose > 1 || fVerbose == -32) cout << "------->  event NOT triggered (pt = " << fMuPt << ", ptrel = " << fMuPtRel << ")" << endl;
  } else {
    if (fVerbose > 1 || fVerbose == -32) cout << "------->  event     triggered (pt = " << fMuPt << ", ptrel = " << fMuPtRel << ")" << endl;
  }
}
