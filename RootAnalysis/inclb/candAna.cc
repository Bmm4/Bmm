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

    // -- MC filling?
    fDoFill = true;
    if ((3 == fIsMC) && (0 < fMuonProcessType)) {
      fDoFill = false;
    }
    if ((4 == fIsMC) && ((499 < fMuonProcessType && fMuonProcessType < 600) || (0 == fMuonProcessType))) {
      fDoFill = false;
    }
    if ((5 == fIsMC) && ((399 < fMuonProcessType && fMuonProcessType < 500) || (0 == fMuonProcessType))) {
      fDoFill = false;
    }


    // -- for data only fill HLT-triggered (and matched!) muons from good runs
    if (0 == fIsMC) {
      if (!fGoodHLT) fDoFill = false; 
      if (!fJSON) fDoFill = false; 
      if (0) cout << "goodHLT: " << fGoodHLT << " JSON: " << fJSON 
		  << " run: " << fRun << " evt: " << fEvt << " ls: " << fLS
		  << " fDoFill: " << fDoFill
		  << endl;
    }


    fillRedTreeData();
    
    if (fDoFill) {
      fTree->Fill(); 
    }
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(11);
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(31);
  } 

}


// ----------------------------------------------------------------------
void candAna::genAnalysis() {

  
  fMuonProcessType = 0; 

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
  fGenMuPt  = fpGenMuon->fP.Perp();
  fGenMuEta = fpGenMuon->fP.Eta();
  fGenMuPhi = fpGenMuon->fP.Phi();

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
	fMuonProcessType = 550; 
      } else if (light) {
	fMuonProcessType = 560; 
      } else {
	fMuonProcessType = 510; 
      }
    } else if (ccbar) {
      fMuonProcessType = 530; 
    } else if (tau) {
      fMuonProcessType = 540; 
    } else if (light) {
      fMuonProcessType = 520; 
    } else if (fake) {
      fMuonProcessType = 590; 
    } else {
      fMuonProcessType = 500; 
    }
  } else if (charm) {
    if (tau) {
      fMuonProcessType = 440; 
    } else if (light) {
      fMuonProcessType = 420; 
    } else if (fake) {
      fMuonProcessType = 490; 
    } else {
      fMuonProcessType = 400; 
    }
  }


  //  cout << "==> fMuonProcessType = " << fMuonProcessType 
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

  new TH1D(Form("mon%s", fName.c_str()), Form("mon%s", fName.c_str()), 50, 0., 50.); 

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupRedTree(fTree); 
  
  // -- Analysis distributions
  new TH1D("analysisDistributions", "analysisDistributions", 10000, 0., 10000.); 


  cout << "booking hists for " << fName << endl;
  TH1::SetDefaultSumw2(kTRUE);
  double MAXPTREL(10.);
  int NBINS(50); 
  new TH1D("ptrel_4_6", "ptrel_4_6", NBINS, 0., MAXPTREL); 
  new TH1D("ptrel_6_8", "ptrel_6_8", NBINS, 0., MAXPTREL); 
  new TH1D("ptrel_8_10", "ptrel_8_10", NBINS, 0., MAXPTREL); 
  new TH1D("ptrel_10_15", "ptrel_10_15", NBINS, 0., MAXPTREL); 
  new TH1D("ptrel_15_20", "ptrel_15_20", NBINS, 0., MAXPTREL); 
  new TH1D("ptrel_20_30", "ptrel_20_30", NBINS, 0., MAXPTREL); 
  new TH1D("ptrel_30_40", "ptrel_30_40", NBINS, 0., MAXPTREL); 
  new TH1D("ptrel_40_50", "ptrel_40_50", NBINS, 0., MAXPTREL); 
  
  TH1D *h(0);
  TH2D *h2(0); 
  char hname[1000], htitle[1000];
  // -- level: 
  //    5 trigger-matched muon + jet
  //    (6 trigger-matched muon with impact parameter)
  int i(5); 
  // -- tag
  //    0 data???
  //    1 b 
  //    2 c
  //    3 uds
  //    ... more to come ...
  vector<int> vTag; 
  vTag.push_back(0); 
  vTag.push_back(1); 
  vTag.push_back(2); 
  vTag.push_back(3); 
  vTag.push_back(8);  // FCR
  vTag.push_back(9);  // FEC
  vTag.push_back(10); // GS
  for (unsigned int jv = 0; jv < vTag.size(); ++jv) {
    int j = vTag[jv];
    sprintf(hname, "RECO_%d_%d_ptrelvsmuoneta", i, j); 
    sprintf(htitle, "RECO ptrel vs muon eta: level %d tag %d", i, j);
    h2 = new TH2D(hname, htitle, 80, -4., 4., 100, 0., 10.); 
    setTitles(h2, "#eta^{muon}", "p_{T}^{rel} [GeV]"); 
    h2->Sumw2(); 
    
    sprintf(hname, "RECO_%d_%d_ptrelvsmuonpt", i, j); 
    sprintf(htitle, "RECO ptrel vs muon pt: level %d tag %d", i, j);
    h2 = new TH2D(hname, htitle, 100, 0., 100., 100, 0., 10.); 
    setTitles(h2, "pt^{muon}[GeV]", "p_{T}^{rel} [GeV]"); 
    h2->Sumw2();
    
    for (int m = 1; m < 11; ++m) {
      sprintf(hname, "RECO_%d_%d_ptrelvsmuonpt%c", i, j, 96+m); 
      sprintf(htitle, "RECO ptrel vs muon pt: level %d tag %d",i,j);
      h2 = new TH2D(hname,htitle, 100, 0., 100., 100, 0., 10.); 
      setTitles(h2, "pt^{muon}[GeV]", "p_{T}^{rel} [GeV]"); 
      h2->Sumw2();
      
      sprintf(hname, "RECO_%d_%d_ptrelvsmuoneta%c", i, j, 96+m); 
      sprintf(htitle, "RECO ptrel vs muon eta: level %d tag %d", i, j);
      h2 = new TH2D(hname,htitle, 80, -4., 4., 100, 0., 10.); 
      setTitles(h2, "#eta^{muon}", "p_{T}^{rel} [GeV]"); 
      h2->Sumw2();
    }
    
    
    sprintf(hname, "GEN_%d_%d_muon_pt", i, j); 
    sprintf(htitle, "GEN muon pt: level %d tag %d", i, j);
    h = new TH1D(hname, htitle, 100, 0., 100.); 
    setTitles(h, "p_{T} [GeV]", "events/bin"); 
    h->Sumw2();
    
    sprintf(hname, "GEN_%d_%d_muon_eta", i, j); 
    sprintf(htitle, "GEN muon #eta: level %d tag %d", i, j);
    h = new TH1D(hname,htitle, 80, -4., 4.); 
    setTitles(h, "#eta", "events/bin"); 
    h->Sumw2(); 
    
    sprintf(hname, "RECO_%d_%d_muon_pt", i, j); 
    sprintf(htitle,"RECO muon pt: level %d tag %d", i, j);
    h = new TH1D(hname,htitle, 100, 0., 100.); 
    setTitles(h, "p_{T} [GeV]", "events/bin"); 
    h->Sumw2();
    
    sprintf(hname, "RECO_%d_%d_muon_eta", i, j); 
    sprintf(htitle, "RECO muon #eta: level %d tag %d", i, j);
    h = new TH1D(hname, htitle, 80, -4., 4.); 
    setTitles(h, "#eta", "events/bin"); 
    h->Sumw2();
  }

}


// ----------------------------------------------------------------------
void candAna::setupRedTree(TTree *t) {

  t->Branch("run",     &fRTD.run,          "run/L");
  t->Branch("evt",     &fRTD.evt,          "evt/L");
  t->Branch("ls",      &fRTD.ls,           "ls/I");

  t->Branch("muid",    &fRTD.muid,         "muid/O");
  t->Branch("hlt",     &fRTD.hlt,          "hlt/O");
  t->Branch("hltt",    &fHltType,          "hltt/I");
  t->Branch("ps",      &fHltPs,            "ps/I");
  t->Branch("json",    &fRTD.json,         "json/O");

  t->Branch("type",    &fRTD.type,         "type/I");
  t->Branch("procid",  &fRTD.procid,       "procid/I");

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
  fRTD.procid    = fProcessType;
  fRTD.type      = fMuonProcessType;

  fRTD.muid      = fMuId;
  fRTD.hlt       = fGoodHLT;
  fRTD.json      = fJSON;

  fRTD.pt        = fMuPt; 
  fRTD.eta       = fMuEta; 
  fRTD.phi       = fMuPhi; 
  fRTD.ptrel     = fMuPtRel; 

  fRTD.jpt        = fJetPt; 
  fRTD.jeta       = fJetEta; 
  fRTD.jphi       = fJetPhi; 

  fHistDir->cd();
  
  if (fDoFill && fGoodHLT) {
    if ( 4 <= fMuPt && fMuPt <  6) { ((TH1D*)fHistDir->Get("ptrel_4_6"))->Fill(fMuPtRel); }
    if ( 6 <= fMuPt && fMuPt <  8) { ((TH1D*)fHistDir->Get("ptrel_6_8"))->Fill(fMuPtRel); }
    if ( 8 <= fMuPt && fMuPt < 10) { ((TH1D*)fHistDir->Get("ptrel_8_10"))->Fill(fMuPtRel); }
    if (10 <= fMuPt && fMuPt < 15) { ((TH1D*)fHistDir->Get("ptrel_10_15"))->Fill(fMuPtRel); }
    if (15 <= fMuPt && fMuPt < 20) { ((TH1D*)fHistDir->Get("ptrel_15_20"))->Fill(fMuPtRel); }
    if (20 <= fMuPt && fMuPt < 30) { ((TH1D*)fHistDir->Get("ptrel_20_30"))->Fill(fMuPtRel); }
    if (30 <= fMuPt && fMuPt < 40) { ((TH1D*)fHistDir->Get("ptrel_30_40"))->Fill(fMuPtRel); }
    if (40 <= fMuPt && fMuPt < 50) { ((TH1D*)fHistDir->Get("ptrel_40_50"))->Fill(fMuPtRel); }


    if (0 == fIsMC) {
      ((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuoneta", 0)))->Fill(fMuEta, fMuPtRel); 
      ((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuonpt", 0)))->Fill(fMuPt, fMuPtRel); 
      for (int m = 1; m < 11; ++m) {
	double ptrelm = fMuPtRel*(1.02 - m*0.01);
	((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuoneta%c", 0, 96+m)))->Fill(fMuEta, ptrelm); 
	((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuonpt%c", 0, 96+m)))->Fill(fMuPt, ptrelm); 
      }
      ((TH1D*)fHistDir->Get(Form("RECO_5_%d_muon_pt", 0)))->Fill(fMuPt); 
      ((TH1D*)fHistDir->Get(Form("RECO_5_%d_muon_eta", 0)))->Fill(fMuEta); 
      ((TH1D*)fHistDir->Get(Form("GEN_5_%d_muon_pt", 0)))->Fill(fGenMuPt); 
      ((TH1D*)fHistDir->Get(Form("GEN_5_%d_muon_eta", 0)))->Fill(fGenMuEta); 
    }

    if (5 == fIsMC) {
      ((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuoneta", 1)))->Fill(fMuEta, fMuPtRel); 
      ((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuonpt", 1)))->Fill(fMuPt, fMuPtRel); 
      for (int m = 1; m < 11; ++m) {
	double ptrelm = fMuPtRel*(1.02 - m*0.01);
	((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuoneta%c", 1, 96+m)))->Fill(fMuEta, ptrelm); 
	((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuonpt%c", 1, 96+m)))->Fill(fMuPt, ptrelm); 
      }
      ((TH1D*)fHistDir->Get(Form("RECO_5_%d_muon_pt", 1)))->Fill(fMuPt); 
      ((TH1D*)fHistDir->Get(Form("RECO_5_%d_muon_eta", 1)))->Fill(fMuEta); 
      ((TH1D*)fHistDir->Get(Form("GEN_5_%d_muon_pt", 1)))->Fill(fGenMuPt); 
      ((TH1D*)fHistDir->Get(Form("GEN_5_%d_muon_eta", 1)))->Fill(fGenMuEta); 
      // -- production processes
      if (40 == fProcessType) {
	((TH2D*)fHistDir->Get(Form("RECO_8_%d_ptrelvsmuoneta", 1)))->Fill(fMuEta, fMuPtRel); 
	((TH2D*)fHistDir->Get(Form("RECO_8_%d_ptrelvsmuonpt", 1)))->Fill(fMuPt, fMuPtRel);
      }

      if (41 == fProcessType) {
	((TH2D*)fHistDir->Get(Form("RECO_9_%d_ptrelvsmuoneta", 1)))->Fill(fMuEta, fMuPtRel); 
	((TH2D*)fHistDir->Get(Form("RECO_9_%d_ptrelvsmuonpt", 1)))->Fill(fMuPt, fMuPtRel);
      }

      if (42 == fProcessType) {
	((TH2D*)fHistDir->Get(Form("RECO_10_%d_ptrelvsmuoneta", 1)))->Fill(fMuEta, fMuPtRel); 
	((TH2D*)fHistDir->Get(Form("RECO_10_%d_ptrelvsmuonpt", 1)))->Fill(fMuPt, fMuPtRel);
      }

    }

    if (4 == fIsMC) {
      ((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuoneta", 2)))->Fill(fMuEta, fMuPtRel); 
      ((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuonpt", 2)))->Fill(fMuPt, fMuPtRel); 
      for (int m = 1; m < 11; ++m) {
	double ptrelm = fMuPtRel*(1.02 - m*0.01);
	((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuoneta%c", 2, 96+m)))->Fill(fMuEta, ptrelm); 
	((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuonpt%c", 2, 96+m)))->Fill(fMuPt, ptrelm); 
      }
      ((TH1D*)fHistDir->Get(Form("RECO_5_%d_muon_pt", 2)))->Fill(fMuPt); 
      ((TH1D*)fHistDir->Get(Form("RECO_5_%d_muon_eta", 2)))->Fill(fMuEta); 
      ((TH1D*)fHistDir->Get(Form("GEN_5_%d_muon_pt", 2)))->Fill(fGenMuPt); 
      ((TH1D*)fHistDir->Get(Form("GEN_5_%d_muon_eta", 2)))->Fill(fGenMuEta); 
    }

    if (3 == fIsMC) {
      ((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuoneta", 3)))->Fill(fMuEta, fMuPtRel); 
      ((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuonpt", 3)))->Fill(fMuPt, fMuPtRel); 
      for (int m = 1; m < 11; ++m) {
	double ptrelm = fMuPtRel*(1.02 - m*0.01);
	((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuoneta%c", 3, 96+m)))->Fill(fMuEta, ptrelm); 
	((TH2D*)fHistDir->Get(Form("RECO_5_%d_ptrelvsmuonpt%c", 3, 96+m)))->Fill(fMuPt, ptrelm); 
      }
      ((TH1D*)fHistDir->Get(Form("RECO_5_%d_muon_pt", 3)))->Fill(fMuPt); 
      ((TH1D*)fHistDir->Get(Form("RECO_5_%d_muon_eta", 3)))->Fill(fMuEta); 
      ((TH1D*)fHistDir->Get(Form("GEN_5_%d_muon_pt", 3)))->Fill(fGenMuPt); 
      ((TH1D*)fHistDir->Get(Form("GEN_5_%d_muon_eta", 3)))->Fill(fGenMuEta); 
    }

  }
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
  fHltType = 0; 
  fHltPs = 0; 

  TString a; 
  int ps(0); 
  bool result(false), wasRun(false), error(false); 

  // NOTRIGGER, just accept the event
  if ( HLTRANGE.begin()->first == "NOTRIGGER" ) { 
    if(fVerbose>2) cout << "NOTRIGGER requested... " << endl;
    fGoodHLT = true; 
    return;
  }

  // Check HLT
  // For every passed HLT look for a matching tigger from our list.
  // If if it confirmed by our list than match it with an object in the TrgObjv2 list
  // Mark the TrigObjv2 object my add a large number to the index. 
  // Like this it can be recogised in the track match search. 
  if(fVerbose==-32) cout<<" event "<<fEvt<<endl;
  fHltType=0; // reset this variable 
  int foundNumHltObjects=0;
  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i]; 
    ps = fpEvt->fHLTPrescale[i]; 
    wasRun = fpEvt->fHLTWasRun[i]; 
    result = fpEvt->fHLTResult[i]; 
    error  = fpEvt->fHLTError[i]; 
    if (wasRun && result) { // passed
      if (fVerbose>2  || (-32 == fVerbose) ) cout << "passed: " << a << endl;
      if ((a == "digitisation_step") 
	  || (a == "L1simulation_step") 
	  || (a == "digi2raw_step") 
	  || (a == "HLTriggerFinalPath") 
	  || (a == "raw2digi_step") 
	  || (a == "reconstruction_step") 
	  ) { 
	//	cout<<" does this ever happen? " <<a<<endl;
	continue; // skip, go to the next triggr
      }
      
      // loop over our list of HLTs and look for matching 
      // We assume that an event can be matched to only one trigger from our list,
      // that is they have to be exclusive.
      bool good = false;

      string spath; 
      int rmin, rmax; 
      for (map<string, pair<int, int> >::iterator imap = HLTRANGE.begin(); 
	   imap != HLTRANGE.end(); ++imap) {  
	spath = imap->first; 
	rmin = imap->second.first; 
	rmax = imap->second.second; 
	if (!a.CompareTo(imap->first.c_str())) {
	  good=true;
	  if (fVerbose > 2 || -32 == fVerbose  ) 
	    cout << "exact match: " << imap->first.c_str() << " HLT: " << a 
		 << " result: " << result << endl;
	  break; // can we skip the rest?
	}	  
	if (a.Contains(spath.c_str()) && (rmin <= fRun) && (fRun <= rmax)) {
	  good=true;
	  if (fVerbose > 2 || -32 == fVerbose) 
	    cout << "close match: " << imap->first.c_str() << " HLT: " << a 
		 << " result: " << result << " in run " << fRun << endl;
	  break; // can we skip the rest?
	} // end if
	
      } // end for loop 
	
      if(good) {  // for matched hlt paths select the trigger object 	  
	// this trigger matched one in our list
	fGoodHLT = true;
	fHltPs = ps; 
	fHltType = i;  
	fHLTPath = a.Data();
	
	bool foundHltObject = false;
	TTrgObjv2 *pTO;     
	if( (fVerbose>9) || (fVerbose==-32)) 
	  cout<<" TTrgObjv2 objects, size= "<<fpEvt->nTrgObjv2()<<endl;
	for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {  // loop over all saved hlt objects 
	  pTO = fpEvt->getTrgObjv2(i); 
	  
	  if(a == pTO->fHltPath) { // found the right one, matched the hlt-name  
	    // this trig object matches a passed and selected trigger 
	    foundHltObject=true;
	    foundNumHltObjects++;
	    
	    int hltIndex = pTO->fHltIndex;
	    // mark it as selected by adding a large number, so >1000.
	    if(hltIndex<1000) {pTO->fHltIndex = hltIndex + (foundNumHltObjects*1000);}
	    else cout<<" hltIndex>1000 "<<hltIndex<<" problem marking it"<<endl;
	    
	    vector<int> muonIndex = pTO->fIndex;
	    vector<int> muonID = pTO->fID;
	    vector<TLorentzVector> muonP = pTO->fP;
	    int num = muonIndex.size();
	    
	    if( (fVerbose>9) || (fVerbose==-32)) {
	      cout<<" matched: "<<pTO->fHltPath<<" hlt-index: "<<pTO->fHltIndex<<" module label: "
		  <<pTO->fLabel<<" type: "<<pTO->fType<<" num of particles: "<<num<<endl;
	      pTO->dump();
	    }
	    
	    break; // found already, skip the rest
	  } // if matched 
	} // end for loop 
	
	if(!foundHltObject) 
	  cout<<"Warning: canAns::triggerSelection: matching trigger module not found! "
	      <<a<<endl;
	} // end if fGoodHLT
    } // if passed      
  } // end for loop hlt
  if( (fVerbose>9) || (fVerbose==-32)) 
    cout<<" number of found mathcing hlt objects: "<<foundNumHltObjects<<endl;

  bool goodMatch(false); 
  if (fGoodHLT) goodMatch = doTriggerMatching(fpMuon); 
  
  if (fGoodHLT) {
    if (fVerbose > 1 || fVerbose == -32) {
      cout << "------->  event     triggered! mu pt/eta = " << fMuPt << "/" << fMuEta;
      if (goodMatch) {
	cout << " and matched!";
      } else {
	cout << " but NOT matched!"; 
      }
      cout << endl;
    }
  } else {
    if (fVerbose > 1 || fVerbose == -32) cout << "------->  event NOT triggered! mu pt/eta = " << fMuPt << "/" << fMuEta << endl;
  }
  fGoodHLT = fGoodHLT && goodMatch; 
}


// ----------------------------------------------------------------------
// A simple trigger matcher based on deltaR (from Frank) + pt matching.
// check 2 muons, use only the selected hlt objects which correspond to triggers
// which passed and were on out trigger list. 
// Only consider trig objects which match our trigger list. 
// The main cuts are: deltaRthr for DR and deltaPtMatch for pt 
// uses TTrgObjv2
bool candAna::doTriggerMatching(TAnaTrack *fp1) { // call the normal version with (true)
  int indx1=-1;
  const double deltaRthr(0.02); // final cut, Frank had 0.5, change 0.020
  const double deltaPtMatch(0.15); // the pt matching cut 
  const int verboseThr = 30;
  //const bool localPrint = false;
  bool localPrint = (fVerbose==-32) || (fVerbose > verboseThr);
  int mu1match(-1);
  string hlt1;
  double deltaRmin1(100);
  double trigMatchDeltaPt1 = 99.;
  bool match=false;
  TTrgObjv2 *pTO;
  TLorentzVector tlvMu1;
   
  if (localPrint) {
    cout << "mu1: pt,eta,phi: " << fp1->fPlab.Perp() << " " << fp1->fPlab.Eta() << " " << fp1->fPlab.Phi()<< endl;
  }
  
  tlvMu1.SetPtEtaPhiM(fp1->fPlab.Perp(),fp1->fPlab.Eta(),fp1->fPlab.Phi(),MMUON); // assume a muon

  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    //pTO->dump();
    int hltIndex = pTO->fHltIndex;
    if(hltIndex>1000) { // this object was selected, matches our trigger list
      if(localPrint) cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
                         <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<endl;
     
      bool match1=false;
      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
        int index = muonIndex[n];  
        int id = muonID[n];  
        TLorentzVector p = muonP[n];  

        if( abs(id) != 13 ) { // if not muon trigger skip 
          if(localPrint) cout<<" a none hlt-muon found in a trigger object, skip it, id= "
                             <<id<<" "<<pTO->fHltPath<<" "<<pTO->fLabel<<" "<<pTO->fType<<endl;
          continue;  // skip checking non-muon objects 
        }

        // check direction matching
        double deltaR1 = p.DeltaR(tlvMu1);

        if(localPrint) {
          cout<<" particle"<<n<<" index "<<index<<" id "<<id
              <<" pt/eta/phi "<<p.Pt()<<"/"<<p.Eta()<<"/"<<p.Phi()<<" i/n "<<i<<"/"<<n 
              <<" dr "<<deltaR1<<endl;
        }

        // muon 1
        if(deltaR1<deltaRmin1) {
          deltaRmin1=deltaR1;  // best match until now
          if (fVerbose > verboseThr || localPrint) {cout << " mu1 selected "<< deltaR1 <<endl;}
            // check now the pt matching 
          double trigMatchDeltaPt=999.;
          if (fp1->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp1->fPlab.Mag())/fp1->fPlab.Mag(); 
          if( trigMatchDeltaPt < deltaPtMatch ) {  // check if it is good enough
            if (deltaR1<deltaRthr) {
              trigMatchDeltaPt1=trigMatchDeltaPt;
              mu1match = n;
              hlt1 = pTO->fLabel;
              indx1=i;
              match1=true;
            } // if delta 
          } // if pt match 
        } // if direction match 

      } // end for loop n

      // check that at least one module matched both
      match = match || match1; 

    } // end if valid module 

  } // loop over all modules


  if (localPrint) 
    cout << " best match "
         <<indx1<<" "<< deltaRmin1 << " "<<mu1match<<" "<<hlt1<<" "<<trigMatchDeltaPt1<<" "<<endl;
  
  bool HLTmatch = false;
  if(match && mu1match>-1) {
    if(localPrint) cout<<" matching OK "<<indx1<<endl;
    HLTmatch=true;
  }

  return HLTmatch;
}
