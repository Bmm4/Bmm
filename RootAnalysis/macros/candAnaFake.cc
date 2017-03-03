#include "candAnaFake.hh"
#include <cmath>
#include <string>

#include "common/HFMasses.hh"


using namespace std;


// ----------------------------------------------------------------------
candAnaFake::candAnaFake(bmmReader *pReader, std::string name, std::string cutsFile) :
  candAna(pReader, name, cutsFile) {
  cout << "==> candAnaFake: name = " << name << ", reading cutsfile " << cutsFile << endl;
  BLIND = 0;
  readCuts(cutsFile, 1);

}


// ----------------------------------------------------------------------
candAnaFake::~candAnaFake() {
  cout << "==> candAnaFake: destructor..." << endl;
}


// ----------------------------------------------------------------------
bool candAnaFake::anaMC(TAna01Event *evt) {
  return true;
}


//------------------------------------------------------------------------------------------
void candAnaFake::evtAnalysis(TAna01Event *evt) {
  candAna::evtAnalysis(evt);
  return;
}


// ----------------------------------------------------------------------
void candAnaFake::candAnalysis() {

  TAnaMuon  *pm[NTRKMAX];
  TSimpleTrack *pt[NTRKMAX];

  bool skip(false);
  if (!fIsMC) {
    // -- fill arrays based on cands
    if (0 == fpCand) {
      cout << "No candidate with type = " << TYPE << " found, skipping" << endl;
      return;
    }

    candAna::candAnalysis();

    if (!fJSON) {
      //    cout << "failed JSON" << endl;
      skip = true;
    }
    if (fCandChi2 > 3) {
      //    cout << "failed chi2/dof" << endl;
      skip = true;
    }
    if (fCandDoca > 0.01) {
      //    cout << "failed candDoca" << endl;
      skip = true;
    }
    if (!fGoodTracks) {
      //    cout << "failed goodTracks" << endl;
      skip = true;
    }
    if (!fGoodQ) {
      //    cout << "failed goodQ" << endl;
      skip = true;
    }
    // FIXME: These should be replaced with the proper signed 3D impact parameter (to be filled!)
    //  if (fpMuon1->fLip/fpMuon1->fLipE < 1.) skip = true;
    //  if (fpMuon2->fLip/fpMuon2->fLipE < 1.) skip = true;

    if (skip) return;

    if (310 == TRUTHCAND) {
      // see bmm4.org::160217a
      //    if (fCandA > 0.01)        skip = true;
      if (fCandFLxy > 4.0)   skip = true;
      if (fCandFLSxy < 15)   skip = true;
      if (fCandFLS3d < 15)   skip = true;
      if (fCandPvIpS > 5.0)  skip = true;
      if (fCandDoca > 0.004) skip = true;

      // -- calculate Lambda mass
      TLorentzVector pr, pi, la;
      pr.SetPtEtaPhiM(fpMuon1->fPlab.Perp(), fpMuon1->fPlab.Eta(), fpMuon1->fPlab.Phi(), MPROTON);
      pi.SetPtEtaPhiM(fpMuon2->fPlab.Perp(), fpMuon2->fPlab.Eta(), fpMuon2->fPlab.Phi(), MPION);
      la = pr + pi;
      fCandMLambda = la.M();
      //      if (TMath::Abs(mla - MLAMBDA_0) < 0.006) skip = true;

      // -- calculate eta -> mu mu  mass
      TLorentzVector m1, m2, eta;
      m1.SetPtEtaPhiM(fpMuon1->fPlab.Perp(), fpMuon1->fPlab.Eta(), fpMuon1->fPlab.Phi(), MMUON);
      m2.SetPtEtaPhiM(fpMuon2->fPlab.Perp(), fpMuon2->fPlab.Eta(), fpMuon2->fPlab.Phi(), MMUON);
      eta = m1 + m2;
      fCandMEta = eta.M();

      fFakeId[0] = 211;
      fFakeId[1] = 211;
    } else if (333 == TRUTHCAND) {
      // see bmm4.org::160217a
      if (fCandDoca > 0.004) skip = true;
      fFakeId[0] = 321;
      fFakeId[1] = 321;
    } else if (3122 == TRUTHCAND) {
      // root [227] events->Draw("m", "1.0<m&&m<1.2&&fls3d>15&&flxy<4&&maxdoca<0.005&&pvips<5", "")
      //      if (fCandFLxy > 4.0)   skip = true;
      if (fCandFLSxy < 15)   skip = true;
      if (fCandFLS3d < 15)   skip = true;
      if (fCandDoca > 0.004) skip = true;
      //      if (fCandPvIpS > 5.0)  skip = true;

      // -- calculate KS masses
      TLorentzVector pi1, pi2, ks;
      pi1.SetPtEtaPhiM(fpMuon1->fPlab.Perp(), fpMuon1->fPlab.Eta(), fpMuon1->fPlab.Phi(), MPION);
      pi2.SetPtEtaPhiM(fpMuon2->fPlab.Perp(), fpMuon2->fPlab.Eta(), fpMuon2->fPlab.Phi(), MPION);
      ks = pi1 + pi2;
      fCandMKS = ks.M();
      if (TMath::Abs(fCandMKS - MKSHORT) < 0.025) {
	//      cout << " mks = " << mks << endl;
	//	skip = true;
      }
      // -- calculate gamma conversion masses
      TLorentzVector el1, el2, conv;
      el1.SetPtEtaPhiM(fpMuon1->fPlab.Perp(), fpMuon1->fPlab.Eta(), fpMuon1->fPlab.Phi(), MELECTRON);
      el2.SetPtEtaPhiM(fpMuon2->fPlab.Perp(), fpMuon2->fPlab.Eta(), fpMuon2->fPlab.Phi(), MELECTRON);
      conv = el1 + el2;
      fCandMconv = conv.M();
      // -- proton has higher pT than pion
      if (211 == TMath::Abs(fpMuon1->fMCID)) {
	skip = true;
      }
      fFakeId[0] = 2212;
      fFakeId[1] = 211;

    }
    if (skip) return;

    ((TH1D*)fHistDir->Get("mall"))->Fill(fpCand->fMass);

    pt[0] = fpEvt->getSimpleTrack(fpMuon1->fIndex);
    pt[1] = fpEvt->getSimpleTrack(fpMuon2->fIndex);

    pm[0] = (fpMuon1->fMuIndex>-1?fpEvt->getMuon(fpMuon1->fMuIndex):0);
    pm[1] = (fpMuon2->fMuIndex>-1?fpEvt->getMuon(fpMuon2->fMuIndex):0);

    fFakeNtrk = 2;
  } else {
    // -- fill arrays based on MC truth
    TSimpleTrack *pT(0);
    int cnt(0);
    for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
      pT = fpEvt->getSimpleTrack(i);
      if (pT->getP().Perp() < 4.0) continue;
      if (TMath::Abs(pT->getP().Eta()) > 2.4) continue;
      // -- check if simple track is matched to gen level
      if (pT->getGenIndex() < 0) continue;
      // -- check if simple track has entry in muon list
      TAnaMuon *pM(0);
      if (pT->getMuonID() > 0) {
	bool isGM(false);
	int idx = pT->getIndex();
	for (int im = 0; im < fpEvt->nMuons(); ++im) {
	  pM = fpEvt->getMuon(im);
	  if (idx == pM->fIndex) {
	    // -- found the corresponding muon, now check if it is a global muon
	    isGM = ((pM->fMuID & 2) == 2);
	    // cout << "check if simple track idx = " << i << " is a global muon: " << (isGM?"yes":"no")
	    //  	 << ": genidx = " << pT->getGenIndex()
	    // 	 << " pT(ST) = " << pT->getP().Perp() << " pT(TAnaMuon) = " << pM->fPlab.Perp()
	    //    	 << endl;
	    if (isGM) break;
	  } else {
	    pM = 0;
	  }
	}
      }
      pt[cnt]      = pT;
      pm[cnt]      = pM;
      fFakeId[cnt] = TMath::Abs(fpEvt->getGenCand(pT->getGenIndex())->fID);
      ++cnt;
      if (cnt == NTRKMAX-1) break;
    }
    fFakeNtrk = cnt;
  }

  // -- now fill into tree
  double result(0.);
  double d2trig, dm1, dm2;
  for (int im = 0; im < fFakeNtrk; ++im) {
    result = 0.;

    fFakePt[im]                = pt[im]->getP().Perp();
    fFakeEta[im]               = pt[im]->getP().Eta();
    fFakePhi[im]               = pt[im]->getP().Phi();
    fFakeQ[im]                 = pt[im]->getCharge();
    fFakeHP[im]                = pt[im]->getHighPurity();
    dist2PdTrigger(pt[im], d2trig, dm1, dm2);
    fFakeDistTrig[im]          = d2trig;
    fFakeDistMuon[im]          = distToMuon(pt[im]);
    fFakeTis[im]               = tis(fpCand);

    if (0 == pm[im]) {
      fFakeGm[im]                     = -99;

      fFakeBdt[im]                    = -99.;
      fFakeTip[im]                    = -99.;
      fFakeLip[im]                    = -99.;
      fFakeInnerChi2[im]              = -99.;
      fFakeOuterChi2[im]              = -99.;
      fFakeChi2LocalPosition[im]      = -99.;
      fFakeChi2LocalMomentum[im]      = -99.;
      fFakeStaTrkMult[im]             = -99.;
      fFakeTmTrkMult[im]              = -99.;
      fFakeDeltaR[im]                 = -99.;
      fFakeItrkValidFraction[im]      = -99.;
      fFakeSegmentComp[im]            = -99.;
      fFakeGtrkNormChi2[im]           = -99.;
      fFakeDxyRef[im]                 = -99.;
      fFakeDzRef[im]                 = -99.;
      fFakeDz[im]                     = -99.;
      fFakeGtrkProb[im]               = -99.;
      fFakeNumberOfValidTrkHits[im]   = -99.;
      fFakeNumberOfLostTrkHits[im]    = -99.;
      fFakeMuonChi2[im]               = -99.;
      fFakeGlbKinkFinder[im]          = -99.;
      fFakeStaRelChi2[im]             = -99.;
      fFakeTrkRelChi2[im]             = -99.;
      fFakeGlbDeltaEtaPhi[im]         = -99.;
      fFakeTimeInOut[im]              = -99.;
      fFakeTimeInOutE[im]             = -99.;
      fFakeTimeInOutS[im]             = -99.;

      fFakeNvalidMuonHits[im]              = -99;
      fFakeNmatchedStations[im]            = -99;
      fFakeLayersWithHits[im]              = -99;
      fFakeNumberOfValidPixHits[im]        = -99;
      fFakeRPChits[im]                     = -99;
      fFakeRPChits1[im]                    = -99;
      fFakeRPChits2[im]                    = -99;
      fFakeRPChits3[im]                    = -99;
      fFakeRPChits4[im]                    = -99;
      fFakeMuonDetectorHitsCombination[im] = -99;
      continue;
    } else {
      fFakeGm[im]                     = ((pm[im]->fMuID & 2) == 2);

      mvaMuon(pm[im], result, false);
      fFakeBdt[im]                    = result;

      fFakeTip[im]                    = pm[im]->fTip;
      fFakeLip[im]                    = pm[im]->fLip;
      fFakeInnerChi2[im]              = pm[im]->fInnerChi2;
      fFakeOuterChi2[im]              = (pm[im]->fOuterChi2 < 100.?pm[im]->fOuterChi2:101.);
      fFakeChi2LocalPosition[im]      = (pm[im]->fChi2LocalPosition < 100?pm[im]->fChi2LocalPosition:101.);
      fFakeChi2LocalMomentum[im]      = (pm[im]->fChi2LocalMomentum < 100.?pm[im]->fChi2LocalMomentum:101.);
      fFakeStaTrkMult[im]             = pm[im]->fStaTrkMult;
      fFakeTmTrkMult[im]              = pm[im]->fTmTrkMult;
      fFakeDeltaR[im]                 = pm[im]->fInnerPlab.DeltaR(pm[im]->fOuterPlab);
      fFakeItrkValidFraction[im]      = pm[im]->fItrkValidFraction;
      fFakeSegmentComp[im]            = pm[im]->fSegmentComp;
      fFakeGtrkNormChi2[im]           = (pm[im]->fGtrkNormChi2 < 10?pm[im]->fGtrkNormChi2:11.);
      fFakeDxyRef[im]                 = pm[im]->fTip;
      fFakeDzRef[im]                  = pm[im]->fLip;
      fFakeDz[im]                     = pm[im]->fdz;
      fFakeLip[im]                    = pm[im]->fLip;
      fFakeGtrkProb[im]               = (pm[im]->fGtrkProb < 100?pm[im]->fGtrkProb:101);
      fFakeNumberOfValidTrkHits[im]   = pm[im]->fNumberOfValidTrkHits;
      fFakeNumberOfLostTrkHits[im]    = pm[im]->fNumberOfLostTrkHits;
      fFakeMuonChi2[im]               = pm[im]->fMuonChi2;
      fFakeGlbKinkFinder[im]          = pm[im]->fGlbKinkFinder;
      fFakeStaRelChi2[im]             = (pm[im]->fStaRelChi2< 100.?pm[im]->fStaRelChi2:101);
      fFakeTrkRelChi2[im]             = pm[im]->fTrkRelChi2;
      fFakeGlbDeltaEtaPhi[im]         = pm[im]->fGlbDeltaEtaPhi;
      fFakeTimeInOut[im]              = pm[im]->fTimeInOut;
      fFakeTimeInOutE[im]             = pm[im]->fTimeInOutE;
      fFakeTimeInOutS[im]             = (pm[im]->fTimeInOutE > 1.e-15?pm[im]->fTimeInOut/pm[im]->fTimeInOutE :-98.);

      fFakeNvalidMuonHits[im]         = pm[im]->fNvalidMuonHits;
      fFakeNmatchedStations[im]       = pm[im]->fNmatchedStations;
      fFakeLayersWithHits[im]         = pm[im]->fLayersWithHits;
      fFakeNumberOfValidPixHits[im]   = pm[im]->fNumberOfValidPixHits;
      unsigned int size = pm[im]->fvRPChits.size();
      if (size == 4) {
	fFakeRPChits1[im] = pm[im]->fvRPChits[0];
	fFakeRPChits2[im] = pm[im]->fvRPChits[1];
	fFakeRPChits3[im] = pm[im]->fvRPChits[2];
	fFakeRPChits4[im] = pm[im]->fvRPChits[3];

	fFakeRPChits[im]  = fFakeRPChits1[im];
	fFakeRPChits[im]  += fFakeRPChits2[im];
	fFakeRPChits[im]  += fFakeRPChits3[im];
	fFakeRPChits[im]  += fFakeRPChits4[im];

	fFakeMuonDetectorHitsCombination[im] = getDetVarComb(pm[im]);
      }
    }
  }

  fFakeTree->Fill();

  bool muid(false);
  // -- leg 1
  if (fMu1Pt >  3.5) {
    muid = ((fpMuon1->fMuID & 2) == 2); // global muon!
    ((TH1D*)fHistDir->Get("pt1all"))->Fill(fMu1Pt, fpCand->fMass);
    ((TH1D*)fHistDir->Get("eta1all"))->Fill(TMath::Abs(fMu1Eta), fpCand->fMass);
    if (muid) {
      ((TH1D*)fHistDir->Get("pt1muo"))->Fill(fMu1Pt, fpCand->fMass);
      ((TH1D*)fHistDir->Get("eta1muo"))->Fill(TMath::Abs(fMu1Eta), fpCand->fMass);
    } else {
      ((TH1D*)fHistDir->Get("pt1nmu"))->Fill(fMu1Pt, fpCand->fMass);
      ((TH1D*)fHistDir->Get("eta1nmu"))->Fill(TMath::Abs(fMu1Eta), fpCand->fMass);
    }
  }

  // -- leg 2
  if (fMu2Pt >  3.5) {
    muid = ((fpMuon2->fMuID & 2) == 2); // global muon!
    ((TH1D*)fHistDir->Get("pt2all"))->Fill(fMu2Pt, fpCand->fMass);
    ((TH1D*)fHistDir->Get("eta2all"))->Fill(TMath::Abs(fMu2Eta), fpCand->fMass);
    if (muid) {
      ((TH1D*)fHistDir->Get("pt2muo"))->Fill(fMu2Pt, fpCand->fMass);
      ((TH1D*)fHistDir->Get("eta2muo"))->Fill(TMath::Abs(fMu2Eta), fpCand->fMass);
    } else {
      ((TH1D*)fHistDir->Get("pt2nmu"))->Fill(fMu2Pt, fpCand->fMass);
      ((TH1D*)fHistDir->Get("eta2nmu"))->Fill(TMath::Abs(fMu2Eta), fpCand->fMass);
    }
  }

  return;
}

// ----------------------------------------------------------------------
void candAnaFake::dumpHFTruthCand(TAnaCand *pC) {
  TAnaTrack *pT(0);
  if( pC->fSig1 == -1 && pC->fSig2==-1 ) return;
  for (int i = pC->fSig1; i <= pC->fSig2; ++i) {
    pT = fpEvt->getRecTrack(fpEvt->getSigTrack(i)->fIndex);
    pT->dump();
  }
}

// ----------------------------------------------------------------------
int candAnaFake::truthMatch(TAnaCand *pCand, int verbose) {
  return 0;
}


// ----------------------------------------------------------------------
void candAnaFake::bookHist() {
  cout << "==>candAnaFake: bookHist" << endl;

  candAna::bookHist();

  fHistDir->cd();

  fFakeTree = new TTree("fakeTree", "fakeTree");
  fFakeTree->Branch("run",     &fRun,               "run/L");
  fFakeTree->Branch("evt",     &fEvt,               "evt/L");

  fFakeTree->Branch("m",       &fCandM,             "m/D");
  fFakeTree->Branch("mks",     &fCandMKS,           "mks/D");
  fFakeTree->Branch("mconv",   &fCandMconv,         "mconv/D");
  fFakeTree->Branch("mlambda", &fCandMLambda,       "mlambda/D");
  fFakeTree->Branch("meta",    &fCandMEta,          "meta/D");

  fFakeTree->Branch("pvip",    &fCandPvIp,          "pvip/D");
  fFakeTree->Branch("pvips",   &fCandPvIpS,         "pvips/D");

  fFakeTree->Branch("chi2dof", &fCandChi2Dof,       "chi2dof/D");
  fFakeTree->Branch("fls3d",   &fCandFLS3d,         "fls3d/D");
  fFakeTree->Branch("fl3d",    &fCandFL3d,          "fl3d/D");
  fFakeTree->Branch("flxy",    &fCandFLxy,          "flxy/D");
  fFakeTree->Branch("fl3dE",   &fCandFL3dE,         "fl3dE/D");
  fFakeTree->Branch("flsxy",   &fCandFLSxy,         "flsxy/D");
  fFakeTree->Branch("maxdoca", &fCandDoca,          "maxdoca/D");
  fFakeTree->Branch("tis",     &fTIS,               "tis/O");
  fFakeTree->Branch("cowboy",  &fCowboy,            "cowboy/O");

  fFakeTree->Branch("ntrk",    &fFakeNtrk,          "ntrk/I");
  fFakeTree->Branch("id",      fFakeId,             "id[ntrk]/I");
  fFakeTree->Branch("q",       fFakeQ,              "q[ntrk]/I");
  fFakeTree->Branch("gm",      fFakeGm,             "gm[ntrk]/I");
  fFakeTree->Branch("pt",      fFakePt,             "pt[ntrk]/F");
  fFakeTree->Branch("eta",     fFakeEta,            "eta[ntrk]/F");
  fFakeTree->Branch("phi",     fFakePhi,            "phi[ntrk]/F");
  fFakeTree->Branch("hp",      fFakeHP,             "hp[ntrk]/O");
  fFakeTree->Branch("tis",     fFakeTis,            "tis[ntrk]/F");
  fFakeTree->Branch("dtrig",   fFakeDistTrig,       "dtrig[ntrk]/F");
  fFakeTree->Branch("dmuon",   fFakeDistMuon,       "dmuon[ntrk]/F");

  fFakeTree->Branch("bdt",         fFakeBdt,               "bdt[ntrk]/F");

  fFakeTree->Branch("deltar",               fFakeDeltaR,               "deltar[ntrk]/F");                    // tree->deltaR = swMu->fInnerPlab.DeltaR(swMu->fOuterPlab);
  fFakeTree->Branch("gtrknormchi2",         fFakeGtrkNormChi2,         "gtrknormchi2[ntrk]/F");              // tree->gNchi2 = swMu->fGtrkNormChi2;
  fFakeTree->Branch("nvalidmuonhits",       fFakeNvalidMuonHits,       "nvalidmuonhits[ntrk]/I");            // tree->vMuHits = swMu->fNvalidMuonHits;
  fFakeTree->Branch("nmatchedstations",     fFakeNmatchedStations,     "nmatchedstations[ntrk]/I");          // tree->mMuStations = swMu->fNmatchedStations;
  fFakeTree->Branch("dxyref",               fFakeDxyRef,               "dxyref[ntrk]/F");                    // tree->dxyRef = swMu->fTip;
  fFakeTree->Branch("dzref",                fFakeDzRef,                "dzref[ntrk]/F");                     // tree->dzRef = swMu->fLip;
  fFakeTree->Branch("layerswithhits",       fFakeLayersWithHits,       "layerswithhits[ntrk]/I");            // tree->LWH = swMu->fLayersWithHits;
  fFakeTree->Branch("numberofvalidpixhits", fFakeNumberOfValidPixHits, "numberofvalidpixhits[ntrk]/I");      // tree->valPixHits = swMu->fNumberOfValidPixHits;
  fFakeTree->Branch("tip",                  fFakeTip,                  "tip[ntrk]/F");                       // tree->innerChi2 = swMu->fTip;
  fFakeTree->Branch("lip",                  fFakeLip,                  "lip[ntrk]/F");                       // tree->innerChi2 = swMu->fLip;
  fFakeTree->Branch("innerchi2",            fFakeInnerChi2,            "innerchi2[ntrk]/F");                 // tree->innerChi2 = swMu->fInnerChi2;
  fFakeTree->Branch("outerchi2",            fFakeOuterChi2,            "outerchi2[ntrk]/F");                 // tree->outerChi2 = swMu->fOuterChi2;
  fFakeTree->Branch("itrkvalidfraction",    fFakeItrkValidFraction,    "itrkValidFraction[ntrk]/F");         // tree->iValFrac = swMu->fItrkValidFraction;
  fFakeTree->Branch("segmentcomp",          fFakeSegmentComp,          "segmentComp[ntrk]/F");               // tree->segComp = swMu->fSegmentComp;
  fFakeTree->Branch("chi2localposition",    fFakeChi2LocalPosition,    "chi2localposition[ntrk]/F");         // tree->chi2LocPos = swMu->fChi2LocalPosition;
  fFakeTree->Branch("chi2localmomentum",    fFakeChi2LocalMomentum,    "chi2localmomentum[ntrk]/F");         // tree->chi2LocMom = swMu->fChi2LocalMomentum;
  fFakeTree->Branch("gtrktailprob",         fFakeGtrkProb,             "gtrktailprob[ntrk]/F");              // tree->glbTrackTailProb = swMu->fGtrkProb;
  fFakeTree->Branch("numberofvalidtrkhits", fFakeNumberOfValidTrkHits, "numberofvalidtrkhits[ntrk]/I");      // tree->NTrkVHits = swMu->fNumberOfValidTrkHits;
  fFakeTree->Branch("muonchi2",             fFakeMuonChi2,             "muonchi2[ntrk]/F");                  // tree->kinkFinder = swMu->fMuonChi2;
  fFakeTree->Branch("rpchits",              fFakeRPChits,              "rpchits[ntrk]/I");                   //  tree->vRPC = swMu->fvRPChits;
  // no entry for DT hits, needs to be changed in BDT setup
  fFakeTree->Branch("glbkinkfinder",        fFakeGlbKinkFinder,        "glbkinkfinder[ntrk]/F");             // tree->glbKinkFinder = swMu->fGlbKinkFinder;
  fFakeTree->Branch("starelchi2",           fFakeStaRelChi2,           "starelchi2[ntrk]/F");                // tree->staRelChi2 = swMu->fStaRelChi2;
  fFakeTree->Branch("trkrelchi2",           fFakeTrkRelChi2,           "trkrelchi2[ntrk]/F");                // tree->trkRelChi2 = swMu->fTrkRelChi2;
  fFakeTree->Branch("glbdeltaetaphi",       fFakeGlbDeltaEtaPhi,       "glbdeltaetaphi[ntrk]/F");            // tree->glbDeltaEtaPhi = swMu->fGlbDeltaEtaPhi;
  fFakeTree->Branch("timeinout",            fFakeTimeInOut,            "timeinout[ntrk]/F");                 // tree->timeAtIpInOut = swMu->fTimeInOut;
  fFakeTree->Branch("timeinoute",           fFakeTimeInOutE,           "timeinoute[ntrk]/F");                // tree->timeAtIpInOutErr = swMu->fTimeInOutE;
  fFakeTree->Branch("timeinouts",           fFakeTimeInOutS,           "timeinouts[ntrk]/F");
  // no entry for DT hits per layer, needs to be changed in BDT setup
  fFakeTree->Branch("rpchits1",             fFakeRPChits1,             "rpchits1[ntrk]/I");                  // tree->vRPC_1 = (swMu->fvRPChits)[0];
  fFakeTree->Branch("rpchits2",             fFakeRPChits2,             "rpchits2[ntrk]/I");                  // tree->vRPC_1 = (swMu->fvRPChits)[1];
  fFakeTree->Branch("rpchits3",             fFakeRPChits3,             "rpchits3[ntrk]/I");                  // tree->vRPC_1 = (swMu->fvRPChits)[2];
  fFakeTree->Branch("rpchits4",             fFakeRPChits4,             "rpchits4[ntrk]/I");                  // tree->vRPC_1 = (swMu->fvRPChits)[3];
  fFakeTree->Branch("statrkmult",           fFakeStaTrkMult,           "statrkmult[ntrk]/F");                // tree->STATrkMult150 = swMu->fStaTrkMult;
  fFakeTree->Branch("tmtrkmult",            fFakeTmTrkMult,            "tmtrkmult[ntrk]/F");                 // tree->TMTrkMult100 = swMu->fTmTrkMult;
  fFakeTree->Branch("mudethitscomb",        fFakeMuonDetectorHitsCombination, "mudethitscomb[ntrk]/I");      // the combination

  fFakeTree->Branch("numberoflosttrkhits",  fFakeNumberOfLostTrkHits,  "numberoflosttrkhits[ntrk]/I");













  const double xbins[] = {0., 4., 6., 10., 15., 25., 50.};
  const int nbins(6);

  string name("K_{S}");
  double mmin(0.45), mmax(0.55);
  double mbin(50);
  TH1::SetDefaultSumw2(kTRUE);
  vector<string> mode;
  mode.push_back("all");
  mode.push_back("nmu");
  mode.push_back("muo");
  if (310 == TRUTHCAND) {
    new TH1D("mall", Form("%s", name.c_str()), mbin, mmin, mmax);
    for (unsigned i = 0; i < mode.size(); ++i) {
      new TH2D(Form("pt1%s", mode[i].c_str()), Form("pt1 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("pt2%s", mode[i].c_str()), Form("pt2 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("eta1%s", mode[i].c_str()), Form("eta1 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
      new TH2D(Form("eta2%s", mode[i].c_str()), Form("eta2 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
    }
  } else if (333 == TRUTHCAND) {
    name = "#phi";
    mmin = 0.98;
    mmax = 1.06;
    new TH1D("mall", Form("%s", name.c_str()), mbin, mmin, mmax);
    for (unsigned i = 0; i < mode.size(); ++i) {
      new TH2D(Form("pt1%s", mode[i].c_str()), Form("pt1 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("pt2%s", mode[i].c_str()), Form("pt2 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("eta1%s", mode[i].c_str()), Form("eta1 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
      new TH2D(Form("eta2%s", mode[i].c_str()), Form("eta2 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
    }
  } else if (443 == TRUTHCAND) {
    name = "#psi";
    mmin = 2.7;
    mmax = 3.4;
    new TH1D("mall", Form("%s", name.c_str()), mbin, mmin, mmax);
    for (unsigned i = 0; i < mode.size(); ++i) {
      new TH2D(Form("pt1%s", mode[i].c_str()), Form("pt1 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("pt2%s", mode[i].c_str()), Form("pt2 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("eta1%s", mode[i].c_str()), Form("eta1 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
      new TH2D(Form("eta2%s", mode[i].c_str()), Form("eta2 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
    }
  } else if (3122 == TRUTHCAND) {
    name = "#Lambda";
    mmin = 1.05;
    mmax = 1.20;
    new TH1D("mall", Form("%s", name.c_str()), mbin, mmin, mmax);
    for (unsigned i = 0; i < mode.size(); ++i) {
      new TH2D(Form("pt1%s", mode[i].c_str()), Form("pt1 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("pt2%s", mode[i].c_str()), Form("pt2 pion (%s)", mode[i].c_str()), nbins, xbins, mbin, mmin, mmax);
      new TH2D(Form("eta1%s", mode[i].c_str()), Form("eta1 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
      new TH2D(Form("eta2%s", mode[i].c_str()), Form("eta2 pion (%s)", mode[i].c_str()), 6, 0., 2.4, mbin, mmin, mmax);
    }
  }

}

// ----------------------------------------------------------------------
void candAnaFake::readCuts(string filename, int dump) {
  candAna::readCuts(filename, dump);

  fCutFile = filename;

  if (dump) cout << "==> candAnaFake: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines;
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[200];
  fHistDir->cd();
  TH1D *hcuts = (TH1D*)fHistDir->Get("hcuts");
  hcuts->GetXaxis()->SetBinLabel(200, fCutFile.c_str());
  string cstring = "B cand";
  int ibin;

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str());

    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

  } // end for

}


// ----------------------------------------------------------------------
void candAnaFake::genMatch() {
  //  cout << "genMatch()  function" << endl;
}

// ----------------------------------------------------------------------
void candAnaFake::genMatchOld() {
  //  cout << "genMatchOld()  function" << endl;
}


// ----------------------------------------------------------------------
void candAnaFake::recoMatch() {
  //  cout << "recoMatch() function" << endl;
}


// ----------------------------------------------------------------------
void candAnaFake::candMatch() {
  //  cout << "candMatch()  function" << endl;
}


// ----------------------------------------------------------------------
void candAnaFake::efficiencyCalculation() {
  //  cout << "efficiencyCalculation() function" << endl;
}

// ----------------------------------------------------------------------
void candAnaFake::processType() {
  //  cout << "processType() function" << endl;
}
