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
candAna::candAna(bmmReader *pReader, string name, string cutsFile) {
  fpReader = pReader;
  fVerbose = fpReader->fVerbose;
  fYear    = fpReader->fYear;
  fNchan   = -1;
  fName    = name;
  cout << "======================================================================" << endl;
  cout << "==> candAna: name = " << name << ", reading cutsfile " << cutsFile << " setup for year " << fYear << endl;

  MASSMIN = 4.9;
  MASSMAX = 5.9;
  BLIND = 0;

  fL1Seeds = 0;
  fGenBTmi = fGenM1Tmi = fGenM2Tmi = fNGenPhotons = fRecM1Tmi = fRecM2Tmi = fCandTmi = -1;

  fHistDir = gFile->mkdir(fName.c_str());
  pvStudy(true);

  cout << "======================================================================" << endl;

}


// ----------------------------------------------------------------------
candAna::~candAna() {
  cout << "==> candAna: destructor..." << endl;
}

// ----------------------------------------------------------------------
void candAna::endAnalysis() {
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

  h1 = ((TH1D*)fHistDir->Get("hmc"));
  if (h1) h1->Write();
  h1 = ((TH1D*)fHistDir->Get("hda"));
  if (h1) h1->Write();

}



// ----------------------------------------------------------------------
void candAna::evtAnalysis(TAna01Event *evt) {

  ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(1);
  fpEvt = evt;
  fBadEvent = false;

  if (1233 == fVerbose) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "new event: " << fEvt << endl;
    cout << "----------------------------------------------------------------------" << endl;
    for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
      fpEvt->getCand(iC)->dump();
    }
    return;
  }


  if (1234 == fVerbose) {
    fpEvt->dump();
    return;
  }

  if (1235 == fVerbose) {
    play();
    return;
  }

  if (1236 == fVerbose) {
    triggerEff("HLT_Dimuon16_Jpsi_v2", "HLT_DoubleMu4_3_Jpsi_Displaced_v2", 1);
    triggerEff("HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_v2", "HLT_Dimuon0er16_Jpsi_NoVertexing_v2", 2);
    return;
  }

  if (1237 == fVerbose) {
    cout << "--- event " << fEvent << " -------------------------------------------------------------------" << endl;
    fpEvt->dumpGenBlock();
    return;
  }

  if (1238 == fVerbose) {
    fHistDir->cd();
    TAnaCand *pCand(0);
    TAnaTrack *p1(0), *p2(0);
    static int first(1);
    static TH1D *hmc(0), *hda(0);
    if (first == 1) {
      first = 0;
      // hmc = new TH1D("hmc", "", 100, 1.0, 1.05);
      // hda = new TH1D("hda", "", 100, 1.0, 1.05);
      // hmc = new TH1D("hmc", "", 100, 0.45, 0.55);
      // hda = new TH1D("hda", "", 100, 0.45, 0.55);
      hmc = new TH1D("hmc", "", 100, 1.0, 1.2);
      hda = new TH1D("hda", "", 100, 1.0, 1.2);
    }
    int GETYPE(3122);
    int DATYPE(113122);
    //    fpEvt->dumpGenBlock();
    TGenCand *pGC(0), *pG1(0), *pG2(0);
    int nmatch(0);
    for (int iC = 0; iC < fpEvt->nGenCands(); ++iC) {
      pGC = fpEvt->getGenCand(iC);
      if (pGC->fID == GETYPE) {
	pG1 = fpEvt->getGenCand(pGC->fDau1);
	pG2 = fpEvt->getGenCand(pGC->fDau2);
	if (pG1->fP.Perp() > 3.5 && pG2->fP.Perp() > 1.0) {
	  cout << "======================================================================" << endl;
	  cout << "--> found " << GETYPE << " in evt = " << fEvt
	       << " with flight length = " << pG1->fV.Mag()
	       << endl;
	  pGC->dump();
	  pG1->dump();
	  pG2->dump();
	  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
	    int idx = fpEvt->getSimpleTrack(i)->getGenIndex();
	    if (idx == pGC->fDau1) {
	      cout << "reco track " << i << " matches " << pGC->fDau1 << endl;
	      fpEvt->getSimpleTrack(i)->dump();
	      ++nmatch;
	    }
	    if (idx == pGC->fDau2) {
	      cout << "reco track " << i << " matches " << pGC->fDau2 << endl;
	      fpEvt->getSimpleTrack(i)->dump();
	      ++nmatch;
	    }
	  }
	  if (nmatch == 2) break;
	}
      }
    }

    TAnaTrack *q1(0), *q2(0);
    for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
      pCand = fpEvt->getCand(iC);
      if (pCand->fType == DATYPE) {
	q1 = fpEvt->getSigTrack(pCand->fSig1);
	q2 = fpEvt->getSigTrack(pCand->fSig2);
	if (q1->fPlab.Perp() > 3.5 && q2->fPlab.Perp() > 1.5) {
	  hda->Fill(pCand->fMass);
	  cout << "    ==> found " << DATYPE << "!!!" << endl;
	  cout << "pCand->fSig1 = " << pCand->fSig1 << " pCand->fSig2 = " << pCand->fSig2 << endl;
	  q1->dump();
	  q2->dump();
	  cout << "----------------------------------------------------------------------" << endl;
	  break;
	}
      }
    }
    return;
  }

  if (fVerbose>0) {
    cout << "======================================================================" << endl;
    cout << " event: " << fEvt << " run: " << fRun << " LS: " << fLS << " JSON: " << fJSON << " cands: " << fpEvt->nCands() << " verbose: "
	 << fVerbose << " MC: " << fIsMC << endl;
  }

  // -- trigger selection (trigger matching is done AFTER [or in] candAnalysis)
  fTIS = fTOS = fGoodHLT1 = fRefTrigger = false;
  triggerHLT();
  triggerL1T();

  triggerSelection();
  runRange();

  if (fIsMC) {
    genMatch();
    recoMatch();
    candMatch();
    pvStudy();
    if (fBadEvent) {
      cout << "XXXXXXX BAD EVENT XXXXXX SKIPPING XXXXX" << endl;
      return;
    }
    efficiencyCalculation();

  }

  // -- print HLT path(s) that fired
  if (0) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "fired HLT paths" << endl;
    for (int i = 0; i < NHLT; ++i) {
      if (fpEvt->fHLTResult[i]) {
	cout << "   " << fpEvt->fHLTNames[i] << endl;
      }
    }
    cout << "fired L1 bits" << endl;
    for (int i = 0; i < NL1T; ++i) {
      if (fpEvt->fL1TResult[i]) {
	cout << "   " << fpEvt->fL1TNames[i] << endl;
      }
    }
    return;
  }


  // -- special call for non-candidate analysis (e.g. candAnaFake for MC)
  if (TYPE == -1) {
    candAnalysis();
    return;
  }

  TAnaCand *pCand(0);
  if (fVerbose == -66) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << " event: " << fEvt << " run: " << fRun << " nSimpleTracks(): " << fpEvt->nSimpleTracks() << endl;
  }
  bool fillNoCand(true);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    if (fVerbose == -66) {
      cout << Form("%4d", iC) << " cand -> " << pCand->fType << endl;
      continue;
    }

    if (TYPE != pCand->fType) {
      if (fVerbose > 39) cout << "  skipping candidate at " << iC << " which is of type " << pCand->fType
			      << " looking for type "<<TYPE<<endl;
      continue;
    }

    fillNoCand = false;
    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(3);
    fpCand = pCand;
    fCandIdx = iC;

    // -- call derived functions
    candAnalysis();

    // Trigger matching
    fHLTmatch=false;
    if (fGoodHLT && fpMuon1 != NULL && fpMuon2 != NULL){ // do only when 2 muons exist
      // check matching for both muons in parallel
      fHLTmatch = doTriggerMatching(fpMuon1, fpMuon2);
    }

    if (0)
      cout << " cand " << fpCand->fType
	   << " run " << fRun << " ls " << fLS << " event " << fEvt << " json " << fJSON
	   << " chan = " << fChan
	   << " hlt1 " <<fGoodHLT1 << " TOS " << fTOS
	   << " presel = " << fPreselection
	   << endl;

    if (fIsMC) {
      fTree->Fill();

      ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(11);
      ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(31);
      if (fJSON) ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(33);
      if (fJSON&&fGoodHLT) ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(34);
      if (fJSON&&fGoodHLT&&fHLTmatch) ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(35);
    } else {  // DATA
      if (NOPRESELECTION) {
	fPreselection = true;
	((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(32);
      }
      if (BLIND
	  && (fMu1MvaId && fMu2MvaId)
	  && (fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX)
	  ) {
      } else {
	if (fPreselection) {
	  if (fJSON) {
	    fTree->Fill();
	    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(11);
	  } else {
	    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(20);
	  }

	  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(12);
	  //((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(31);
	  if (fJSON) ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(33);
	  if (fJSON&&fGoodHLT) ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(34);
	  if (fJSON&&fGoodHLT&&fHLTmatch) {
	    ((TH1D*)fHistDir->Get(Form("mon%s", fName.c_str())))->Fill(35);
	    ((TH1D*)fHistDir->Get("run1"))->Fill(fRun);
	  }
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

  fpMuon1 = fpMuon2 = 0;

  ((TH1D*)fHistDir->Get("../monEvents"))->Fill(1);

  TAnaVertex *pVtx;
  fPvN = 0;
  for (int i = 0; i < fpEvt->nPV(); ++i) {
    pVtx = fpEvt->getPV(i);
    if (0 == pVtx->fStatus) {
      ++fPvN;
    } else {
      //      cout << "skipping  fake vertex" << endl;
    }
  }

  if (fpCand->fPvIdx > -1 && fpCand->fPvIdx < fpEvt->nPV()) {
    TAnaVertex *pv = fpEvt->getPV(fpCand->fPvIdx);
    fPvIdx = fpCand->fPvIdx;
    fPvX = pv->fPoint.X();
    fPvY = pv->fPoint.Y();
    fPvZ = pv->fPoint.Z();
    fPvNtrk = pv->getNtracks();
    fPvNdof = pv->fNdof;
    fPvAveW8 = ((fPvNdof+2.)/2.)/fPvNtrk;
  } else {
    fPvIdx = -99;
    fPvX = -99.;
    fPvY = -99.;
    fPvZ = -99.;
    fPvNtrk = -99;
  }

  if (fIsMC) {
    if (fCandTmi == fCandIdx) {
      if (fNGenPhotons) {
	fCandTM = 2;
      } else {
	fCandTM = 1;
      }
    } else {
      fCandTM = 0;
    }
  } else {
    fCandTM = 0;
  }

  fCandType     = fpCand->fType;
  fCandPt       = fpCand->fPlab.Perp();
  fCandP        = fpCand->fPlab.Mag();
  fCandEta      = fpCand->fPlab.Eta();
  fCandPhi      = fpCand->fPlab.Phi();
  fCandM        = fpCand->fMass;
  fCandME       = fpCand->fMassE;
  fCandDoca     = fpCand->fMaxDoca;

  // -- values of cand wrt PV
  fCandPvTip    = fpCand->fPvTip;
  fCandPvTipE   = fpCand->fPvTipE;
  fCandPvTipS   = fCandPvTip/fCandPvTipE;
  fCandPvLip    = fpCand->fPvLip;
  fCandPvLipE   = fpCand->fPvLipE;
  fCandPvLipS   = fCandPvLip/fCandPvLipE;
  fCandPv2Lip   = fpCand->fPv2Lip;
  fCandPv2LipS  = fpCand->fPv2Lip/fpCand->fPv2LipE;
  if (fpCand->fPv2Lip > 999) fCandPv2LipS = 999.; // in this case no 2nd best PV was found, reset to 'good' state
  fCandPv12Lip  = fCandPvLip/fpCand->fPv2Lip;
  fCandPv12LipE = fCandPvLipE/fpCand->fPv2LipE;
  fCandPv12LipS = fCandPvLipS/(fpCand->fPv2Lip/fpCand->fPv2LipE);

  // -- 3d impact parameter wrt PV
  if (0) {
    fCandPvIp     = TMath::Sqrt(fCandPvLip*fCandPvLip + fCandPvTip*fCandPvTip);
    fCandPvIpE    = (fCandPvLip*fCandPvLip*fCandPvLipE*fCandPvLipE + fCandPvTip*fCandPvTip*fCandPvTipE*fCandPvTipE)/(fCandPvIp*fCandPvIp);
    fCandPvIpE    = TMath::Sqrt(fCandPvIpE);
    fCandPvIpS    = fCandPvIp/fCandPvIpE;
    if (TMath::IsNaN(fCandPvIpS)) fCandPvIpS = -1.;
  }
  // -- new version directly from CMSSSW and no longer patched...
  fCandPvIp     = fpCand->fPvIP3d;
  fCandPvIpE    = fpCand->fPvIP3dE;
  fCandPvIpS    = fpCand->fPvIP3d/fpCand->fPvIP3dE;
  if (TMath::IsNaN(fCandPvIpS)) fCandPvIpS = -1.;

  fCandPvIp3D   = fpCand->fPvIP3d;
  fCandPvIpE3D  = fpCand->fPvIP3dE;
  fCandPvIpS3D  = fpCand->fPvIP3d/fpCand->fPvIP3dE;
  if (TMath::IsNaN(fCandPvIpS3D)) fCandPvIpS3D = -1.;

  //old  fCandM2 = constrainedMass();
  //  fCandM2 = fpCand->fDouble1;
  fCandM2 = fpCand->fMassC;

  // -- new variables
  fCandPvDeltaChi2 = fpCand->fDeltaChi2;

  TAnaTrack *p0;
  TAnaTrack *p1(0);
  TAnaTrack *p2(0);

  fCandQ    = 0;

  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    p0 = fpEvt->getSigTrack(it);
    fCandQ += p0->fQ;
    if (TMath::Abs(p0->fMCID) != 13) continue;
    if (0 == p1) {
      p1 = p0;
    } else {
      p2 = p0;
    }
  }

  // -- for rare backgrounds there are no "true" muons
  if (fpCand->fType > 1000000 && fpCand->fType < 2000000) {
    p1 = fpEvt->getSigTrack(fpCand->fSig1);
    p2 = fpEvt->getSigTrack(fpCand->fSig2);
  }

  // -- candidates for fake rate determination from light hadrons
  if (fpCand->fType == 11310 || fpCand->fType == 11333 || fpCand->fType == 113122) {
    p1 = fpEvt->getSigTrack(fpCand->fSig1);
    p2 = fpEvt->getSigTrack(fpCand->fSig2);
  }
  if (fpCand->fType == 2000100 || fpCand->fType == 2000101 || fpCand->fType == 2000102) {
    p1 = fpEvt->getSigTrack(fpCand->fSig1);
    p2 = fpEvt->getSigTrack(fpCand->fSig2);
  }

  //  Dstar (prompt and from Bd)
  if (fpCand->fType == 54 || fpCand->fType == 300054 || fpCand->fType == 300031) {
    return;
  }

  if (fVerbose>98) {
    cout<<" mu1 "<<p1->fRefPlab.Perp()<<" "<<p1->fRefPlab.Perp()<<" "<<p1->fRefPlab.Perp()<<endl;
    cout<<" mu2 "<<p2->fRefPlab.Perp()<<" "<<p2->fRefPlab.Perp()<<" "<<p2->fRefPlab.Perp()<<endl;
  }

  // -- Bd2DstarPi
  if (fpCand->fType == 600030 || fpCand->fType == 3000030) {
    if (fpCand->fDau1 < 0 || fpCand->fDau1 > fpEvt->nCands()) return;
    TAnaCand *pD = fpEvt->getCand(fpCand->fDau1);
    pD = fpEvt->getCand(pD->fDau1);
    // -- the pion and kaon of the D0 assume the "muon" roles (??)
    p1 = fpEvt->getSigTrack(pD->fSig1);
    p2 = fpEvt->getSigTrack(pD->fSig2);
  }


  if (p1->fPlab.Perp() < p2->fPlab.Perp()) {
    p0 = p1;
    p1 = p2;
    p2 = p0;
  }


  fpMuon1 = p1;
  fpMuon2 = p2;

  fDeltaR = fpMuon1->fRefPlab.DeltaR(fpMuon2->fRefPlab);


  muScaleCorrectedMasses();

  fMu1TrkLayer  = fpReader->numberOfTrackerLayers(p1);
  fMu1GmId      = ((p1->fMuID & 2) == 2);
  fMu1TmId      = tightMuon(p1);

  fMu1BDT       = -1.;
  fMu1MvaId     = mvaMuon(p1, fMu1BDT);
  fMu1rTmId     = tightMuon(p1, false);
  fMu1rBDT      = -1.;
  fMu1rMvaId    = mvaMuon(p1, fMu1rBDT, false);

  fMu1BDTLM     = -1.;
  fMu1MvaIdLM   = mvaMuonLM(p1, fMu1BDTLM);
  fMu1rBDTLM    = -1.;
  fMu1rMvaIdLM  = mvaMuonLM(p1, fMu1rBDTLM, false);

  fTrigMatchDeltaPt = 99.;
  //cout<<" do trigger matching for muon 1 "<<endl;
  // false - consider only selecetd triggers, true - match to muon trigger objects
  fMu1TrigM     = doTriggerMatchingR(p1, false, true, false);
  if (fTrigMatchDeltaPt > 0.1) fMu1TrigM *= -1.;

  //  fMu1Id        = fMu1MvaId && (fMu1TrigM < 0.1) && (fMu1TrigM > 0);
  fMu1Id        = fMu1MvaId;
  if (HLTRANGE.begin()->first == "NOTRIGGER") fMu1Id = true;

  fMu1Pt        = p1->fRefPlab.Perp();
  fMu1Eta       = p1->fRefPlab.Eta();
  fMu1Phi       = p1->fRefPlab.Phi();
  fMu1PtNrf     = p1->fPlab.Perp();
  fMu1EtaNrf    = p1->fPlab.Eta();
  fMu1PhiNrf    = p1->fPlab.Phi();
  fMu1TkQuality = highPurity(p1);
  fMu1Q         = p1->fQ;
  fMu1Pix       = fpReader->numberOfPixLayers(p1);
  fMu1BPix      = fpReader->numberOfBPixLayers(p1);
  fMu1BPixL1    = fpReader->numberOfBPixLayer1Hits(p1);
  fMu1PV        = p1->fPvIdx;
  fMu1IP        = p1->fBsTip;
  fMu1IPE       = p1->fBsTipE;
  fMu1IPS       = p1->fTip/p1->fTipE;

  if (p1->fMuIndex > -1) {
    TAnaMuon *pm = fpEvt->getMuon(p1->fMuIndex);
    fillMuonData(fMu1Data, pm);
    fMu1Data.mbdt = fMu1rBDT;
    fMu1Chi2     = pm->fMuonChi2;
    fMu1Iso      = isoMuon(fpCand, pm);
    fMu1VtxProb  = pm->fVtxProb;
  } else {
    fillMuonData(fMu1Data, 0);
    fMu1Data.mbdt = -98.;
    fMu1Chi2 = -98.;
    fMu1Iso  = -98.;
    fMu1VtxProb = 99.;
  }

  if (fCandTM && fGenM1Tmi < 0) fpEvt->dump();

  TGenCand *pg1(0), *pg2(0);
  if (fCandTmi > -1) {
    pg1           = fpEvt->getGenTWithIndex(p1->fGenIndex);
    fMu1PtGen     = pg1->fP.Perp();
    fMu1EtaGen    = pg1->fP.Eta();
    fMu1PhiGen    = pg1->fP.Phi();
  } else {
    fMu1PtGen     = -99.;
    fMu1EtaGen    = -99.;
    fMu1PhiGen    = -99.;
  }

  //cout<<" check muon 2"<<endl;
  //  fMu2Id        = goodMuon(p2);
  fMu2TrkLayer  = fpReader->numberOfTrackerLayers(p2);
  fMu2GmId      = ((p2->fMuID & 2) == 2);
  fMu2TmId      = tightMuon(p2);

  fMu2BDT       = -1.;
  fMu2MvaId     = mvaMuon(p2, fMu2BDT);
  fMu2rTmId     = tightMuon(p2, false);
  fMu2rBDT      = -1.;
  fMu2rMvaId    = mvaMuon(p2, fMu2rBDT, false);

  fMu2BDTLM     = -1.;
  fMu2MvaIdLM   = mvaMuonLM(p2, fMu2BDTLM);
  fMu2rBDTLM    = -1.;
  fMu2rMvaIdLM  = mvaMuonLM(p2, fMu2rBDTLM, false);

  fTrigMatchDeltaPt = 99.;
  //cout<<" do trigger matching for muon 2 "<<endl;
  fMu2TrigM     = doTriggerMatchingR(p2, false, true, false);

  if (fTrigMatchDeltaPt > 0.1) fMu2TrigM *= -1.;
  //  fMu2Id        = fMu2MvaId && (fMu2TrigM < 0.1) && (fMu2TrigM > 0);
  fMu2Id        = fMu2MvaId;
  if (HLTRANGE.begin()->first == "NOTRIGGER") fMu2Id = true;


  fMu2Pt        = p2->fRefPlab.Perp();
  fMu2Eta       = p2->fRefPlab.Eta();
  fMu2Phi       = p2->fRefPlab.Phi();
  fMu2PtNrf     = p2->fPlab.Perp();
  fMu2EtaNrf    = p2->fPlab.Eta();
  fMu2PhiNrf    = p2->fPlab.Phi();
  fMu2TkQuality = highPurity(p2);
  fMu2Q         = p2->fQ;
  fMu2Pix       = fpReader->numberOfPixLayers(p2);
  fMu2BPix      = fpReader->numberOfBPixLayers(p2);
  fMu2BPixL1    = fpReader->numberOfBPixLayer1Hits(p2);
  fMu2PV        = p2->fPvIdx;
  fMu2IP        = p2->fBsTip;
  fMu2IPE       = p2->fBsTipE;
  fMu2IPS       = p2->fTip/p2->fTipE;

  // -- fill tree for muon id MVA studies
  if (0 && (fMu1rTmId || fMu2rTmId)) {
    TAnaMuon *pt(0);
    for (int i = 0; i < 2; ++i) {
      pt = 0;
      if (0 == i) {
	if (!fMu1rTmId) continue;
	int idx = p1->fMuIndex;
	if (idx > -1 && idx < fpEvt->nMuons()) {
	  pt = fpEvt->getMuon(idx);
	}
      }
      if (1 == i) {
	if (!fMu2rTmId) continue;
	int idx = p2->fMuIndex;
	if (idx > -1 && idx < fpEvt->nMuons()) {
	  pt = fpEvt->getMuon(idx);
	}
      }

      if (0 == pt) {
	cout << "no TAnaMuon found despite fMu1TmId != 0!!" << endl;
	continue;
      }

      //      fillMuonData(fMuonData, pt);

      fMuonData.pt            = pt->fPlab.Perp();
      fMuonData.eta           = pt->fPlab.Eta();
      fMuonData.validMuonHits    = 0;
      fMuonData.glbNChi2         = pt->fGtrkNormChi2;
      fMuonData.nMatchedStations = pt->fNmatchedStations;
      fMuonData.validPixelHits   = fpReader->numberOfPixLayers(pt); // FIXME, kind of correct
      fMuonData.trkLayerWithHits = fpReader->numberOfTrackerLayers(pt);

      fMuonData.trkValidFract = pt->fItrkValidFraction;
      fMuonData.segComp       = pt->fSegmentComp;
      fMuonData.chi2LocMom    = pt->fChi2LocalMomentum;
      fMuonData.chi2LocPos    = pt->fChi2LocalPosition;
      fMuonData.glbTrackProb  = pt->fGtrkProb;
      fMuonData.NTrkVHits     = static_cast<float>(pt->fNumberOfValidTrkHits);
      fMuonData.NTrkEHitsOut  = static_cast<float>(pt->fNumberOfLostTrkHits);

      fMuonData.kink          = pt->fMuonChi2;

      fMuonData.dpt           = pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag();
      fMuonData.dptrel        = TMath::Abs(pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag())/pt->fInnerPlab.Mag();
      if (pt->fOuterPlab.Mag() > 3.) {
	fMuonData.deta          = pt->fInnerPlab.Eta() - pt->fOuterPlab.Eta();
	fMuonData.dphi          = pt->fInnerPlab.DeltaPhi(pt->fOuterPlab);
	fMuonData.dr            = pt->fInnerPlab.DeltaR(pt->fOuterPlab);
      } else {
	fMuonData.deta          = -99.;
	fMuonData.dphi          = -99.;
	fMuonData.dr            = -99.;
      }

      fMuonIdTree->Fill();
    }
  }

  // -- cut on fMuIndex so that fake muons (from rare backgrounds) can be treated above as real muons
  if (p1->fMuIndex > -1 && p2->fMuIndex > -1) {
    TVector3 rm1  = fpEvt->getMuon(p1->fMuIndex)->fPositionAtM2;
    TVector3 rm2  = fpEvt->getMuon(p2->fMuIndex)->fPositionAtM2;

    if (rm1.Mag() > 0.1 && rm2.Mag() > 0.1) {
      TVector3 rD   = rm2-rm1;
      fMuDist   = rD.Mag();
      fMuDeltaR =  rm1.DeltaR(rm2);
    } else {
      fMuDist   = -99.;
      fMuDeltaR = -99.;
    }

    if (fVerbose > 10) cout << "dist: " << fMuDist << " dr = " << fMuDeltaR << endl;
  } else {
    fMuDist   = -99.;
    fMuDeltaR = -99.;
  }

  if (p2->fMuIndex > -1) {
    TAnaMuon *pm = fpEvt->getMuon(p2->fMuIndex);
    fillMuonData(fMu2Data, pm);
    fMu2Data.mbdt = fMu2rBDT;
    fMu2Chi2      = pm->fMuonChi2;
    fMu2Iso       = isoMuon(fpCand, pm);
    fMu2VtxProb   = pm->fVtxProb;
  } else {
    fillMuonData(fMu2Data, 0);
    fMu2Data.mbdt = -98.;
    fMu2Chi2 = -98.;
    fMu2Iso  = -98.;
    fMu2VtxProb = 99.;
  }


  xpDistMuons();

  double dphi = p1->fPlab.DeltaPhi(p2->fPlab);
  fCowboy = (p1->fQ*dphi > 0);

  fMu1W8Mu = fMu2W8Mu = fMu1W8Tr = fMu2W8Tr = -2.;

  if (fCandTmi > -1) {
    pg2           = fpEvt->getGenTWithIndex(p2->fGenIndex);
    fMu2PtGen     = pg2->fP.Perp();
    fMu2EtaGen    = pg2->fP.Eta();
    fMu2PhiGen    = pg2->fP.Phi();
  } else {
    fMu2PtGen     = -99.;
    fMu2EtaGen    = -99.;
    fMu2PhiGen    = -99.;
  }

  fGenMass = -99.;
  if (0 != pg1 && 0 != pg2) {
    TLorentzVector gendimuon = pg1->fP + pg2->fP;
    fGenMass = gendimuon.M();
  }
  //  cout << "m(mu,mu) = " << fGenMass << " n(photons) = " << fNGenPhotons << endl;

  fCandW8Mu     = fMu1W8Mu*fMu2W8Mu;
  if (TMath::Abs(fCandW8Mu) > 1.) fCandW8Mu = 0.2; // FIXME correction for missing entries at low pT
  fCandW8Tr     = fMu1W8Tr*fMu2W8Tr;
  if (TMath::Abs(fCandW8Tr) > 1.) fCandW8Tr = 0.2; // FIXME correction for missing entries at low pT

  // -- FIXME ????
  int pvidx = (fpCand->fPvIdx > -1? fpCand->fPvIdx : 0);
  // -- this is from the full candidate
  TAnaVertex sv = fpCand->fVtx;
  // // -- this is from the dimuon vertex
  // TAnaCand *pD;
  // TAnaVertex sv2m;
  // bool good2m(false);
  // for (int id = fpCand->fDau1; id <= fpCand->fDau2; ++id) {
  //   if (id < 0) break;
  //   pD = fpEvt->getCand(id);
  //   if (300443 == pD->fType) {
  //     sv2m = pD->fVtx;
  //     good2m = true;
  //     break;
  //   }
  // }

  // -- go back to original!
  sv = fpCand->fVtx;

  fCandA      = fpCand->fAlpha;
  fCandCosA   = TMath::Cos(fCandA);

  //  virtual double      isoClassicWithDOCA(TAnaCand*, double dca, double r = 1.0, double ptmin = 0.9);
  double iso = isoClassicWithDOCA(fpCand, 0.05, 0.7, 0.9); // 500um DOCA cut
  fCandIso      = iso;
  fCandPvTrk    = fCandI0trk;
  fCandIsoTrk   = fCandI2trk;
  pair<int, int> pclose;
  pclose = nCloseTracks(fpCand, 0.03, 1, 0.5);
  fCandCloseTrk = pclose.first;
  fCandCloseTrkS1 = pclose.second;
  pclose = nCloseTracks(fpCand, 0.03, 2, 0.5);
  fCandCloseTrkS2 = pclose.second;
  pclose = nCloseTracks(fpCand, 0.03, 3, 0.5);
  fCandCloseTrkS3 = pclose.second;
  pclose = nCloseTracks(fpCand, 0.03, 4, 0.5);
  fCandCloseTrkS4 = pclose.second;
  pclose = nCloseTracks(fpCand, 0.03, 5, 0.5);
  fCandCloseTrkS5 = pclose.second;

  fCandChi2    = sv.fChi2;
  fCandDof     = sv.fNdof;
  fCandProb    = sv.fProb;
  fCandChi2Dof = fCandChi2/fCandDof;

  fCandOtherVtx = TMath::Max(fMu1VtxProb, fMu2VtxProb) - sv.fProb;

  fCandFL3d  = sv.fD3d;
  fCandFL3dE = sv.fD3dE;
  fCandFLS3d = sv.fD3d/sv.fD3dE;
  if (TMath::IsNaN(fCandFLS3d)) fCandFLS3d = -1.;
  fCandFLxy  = sv.fDxy;
  fCandFLSxy = sv.fDxy/sv.fDxyE;
  if (TMath::IsNaN(fCandFLSxy)) fCandFLSxy = -1.;

  fCandTau   = fpCand->fTau3d;
  fCandTauE  = fpCand->fTau3dE;

  fCandTauxy   = fpCand->fTauxy;
  fCandTauxyE  = fpCand->fTauxyE;

  // -- variables for production mechanism studies
  //  fpOsCand      = osCand(fpCand);
  fOsIso        = osIsolation(fpCand, 1.0, 0.9);
  fOsRelIso     = fOsIso/fCandPt;
  int osm       = osMuon(fpCand, 0.);
  fOsMuonPt     = (osm > -1? fpEvt->getSimpleTrack(osm)->getP().Perp():-1.);
  fOsMuonPtRel  = (osm > -1? fpEvt->getSimpleTrack(osm)->getP().Perp(fpCand->fPlab):-1);
  fOsMuonDeltaR =  (osm > -1? fpCand->fPlab.DeltaR(fpEvt->getSimpleTrack(osm)->getP()):-1.);

  // // -- dimuon vertex version
  // if (good2m) {
  //   fmmChi2    = sv2m.fChi2;
  //   fmmDof     = sv2m.fNdof;
  //   fmmProb    = sv2m.fProb;
  //   fmmFL3d    = sv2m.fD3d;
  //   fmmFL3dE   = sv2m.fD3dE;
  //   fmmFLS3d   = sv2m.fD3d/sv2m.fD3dE;
  //   if (TMath::IsNaN(fmmFLS3d)) fmmFLS3d = -1.;
  //   fmmFLSxy   = sv2m.fDxy/sv2m.fDxyE;
  //   if (TMath::IsNaN(fmmFLSxy)) fmmFLSxy = -1.;
  //   // -- pD should be/is still pointing to the J/psi daughter IF good2m == true
  //   fmmCosA    = TMath::Cos(pD->fAlpha);
  //   fmmMaxDoca = pD->fMaxDoca;
  //   fmmPt      = pD->fP.Perp();
  // } else {
  //   fmmChi2  = -1.;
  //   fmmDof   = -1;
  //   fmmProb  = -1.;
  //   fmmFL3d  = -1.;
  //   fmmFL3dE = -1.;
  //   fmmFLS3d = -1.;
  //   fmmFLSxy = -1.;
  //   fmmCosA    = -1.;
  //   fmmMaxDoca = -1.;
  //   fmmPt      = -1.;
  // }

  if (fpCand->fNstTracks.size() == 0) {
    //    cout << "HHHHEEEELLLLPPPP" << endl;
    fCandDocaTrk = 99.;
    fCandDocaTrkBdt = 99.;
  } else {
    fCandDocaTrk    = fpCand->fNstTracks[0].second.first;
    fCandDocaTrkBdt = fpCand->fNstTracks[0].second.first;

    int nsize(fpCand->fNstTracks.size());
    TSimpleTrack *ps;
    for (int i = 0; i<nsize; ++i) {
      int trkId = fpCand->fNstTracks[i].first;
      ps = fpEvt->getSimpleTrack(trkId);
      // -- check that any track associated with a definitive vertex is from the same or the closest (compatible) other PV
      if ((ps->getPvIndex() > -1) && (ps->getPvIndex() != pvidx)) continue;

      fCandDocaTrk = fpCand->fNstTracks[i].second.first;
      break;
    }

  }

  // if (BLIND && fpCand->fMass > SIGBOXMIN && fpCand->fMass < SIGBOXMAX  && fCandIso < 0.7) {
  //   calcBDT();
  // } else {
  //   calcBDT();
  // }

  fChan = detChan(fMu1Eta, fMu2Eta);

  if (0 && fMu2Pt < 4.0) {
    cout << "======== chan = " << fChan << " hlt1 = " << fGoodHLT1 << " m2pt = " << fMu2Pt << " " << fL1SeedString << endl;
  }
  fTOS = tos(fpCand);
  fTIS = tis(fpCand);
  fRefTrigger = refTrigger(fpCand, "HLT_Mu7p5_Track3p5_Jpsi_v2");

  // -- check for good channel, else the cuts array cannot be used!
  if (fChan < 0) {
    //    cout << "Failed chan, l1 seeds = " << fL1SeedString << endl;
    return;
  }


  fWideMass       = ((fpCand->fMass > MASSMIN) && (fpCand->fMass < MASSMAX));

  fGoodMuonsID    = (fMu1Id && fMu2Id);
  fGoodMuonsTmID  = (fMu1TmId && fMu2TmId);
  fGoodMuonsMvaID = (fMu1MvaId && fMu2MvaId);
  fGoodMuonsPt    = ((fMu1Pt > fCuts[fChan]->m1pt) && (fMu1Pt < 14000.) && (fMu2Pt > fCuts[fChan]->m2pt) && (fMu2Pt < 14000.));
  double etaLead(TMath::Abs(fMu1Eta));
  if (TMath::Abs(fMu2Eta) > etaLead) etaLead = TMath::Abs(fMu2Eta);
  fGoodMuonsEta   = ((fCuts[fChan]->metaMin < etaLead) && (etaLead < fCuts[fChan]->metaMax));
  fGoodTracks     = (highPurity(p1) && highPurity(p2));
  fGoodTracksPt   = ((fMu1Pt > TRACKPTLO) && (fMu1Pt < TRACKPTHI) && (fMu2Pt > TRACKPTLO) && (fMu2Pt < TRACKPTHI));
  fGoodTracksEta  = ((fMu1Eta > TRACKETALO) && (fMu1Eta < TRACKETAHI) && (fMu2Eta > TRACKETALO) && (fMu2Eta < TRACKETAHI));

  fGoodQ          = (fMu1Q*fMu2Q < 0);
  fGoodPvAveW8    = (fPvAveW8 > PVAVEW8);
  fGoodPvLip      = (TMath::Abs(fCandPvLip) < fCuts[fChan]->pvlip);
  fGoodPvLipS     = (TMath::Abs(fCandPvLipS) < fCuts[fChan]->pvlips);
  fGoodPv2Lip     = (TMath::Abs(fCandPv2Lip) > fCuts[fChan]->pv2lip);
  fGoodPv2LipS    = (TMath::Abs(fCandPv2LipS) > fCuts[fChan]->pv2lips);
  fGoodMaxDoca    = (TMath::Abs(fCandDoca) < fCuts[fChan]->maxdoca);
  fGoodIp         = (TMath::Abs(fCandPvIp) < fCuts[fChan]->pvip);
  fGoodIpS        = (TMath::Abs(fCandPvIpS) < fCuts[fChan]->pvips);

  fGoodPt         = (fCandPt > fCuts[fChan]->pt);
  fGoodEta        = ((fCandEta > fCuts[fChan]->etaMin) && (fCandEta < fCuts[fChan]->etaMax));
  fGoodAlpha      = (fCandA < fCuts[fChan]->alpha);
  fPreselAlpha    = (fCandA < 0.2);
  fGoodChi2       = (fCandChi2/fCandDof < fCuts[fChan]->chi2dof);
  fGoodFLS        =  ((fCandFLS3d > fCuts[fChan]->fls3d) && (fCandFLSxy > fCuts[fChan]->flsxy));
  if (TMath::IsNaN(fCandFLS3d)) fGoodFLS = false;

  fGoodCloseTrack   = (fCandCloseTrk < fCuts[fChan]->closetrk);
  fGoodCloseTrackS1 = (fCandCloseTrkS1 < fCuts[fChan]->closetrks1);
  fGoodCloseTrackS2 = (fCandCloseTrkS2 < fCuts[fChan]->closetrks2);
  fGoodCloseTrackS3 = (fCandCloseTrkS3 < fCuts[fChan]->closetrks3);
  fGoodIso          = (fCandIso > fCuts[fChan]->iso);
  fGoodM1Iso        = (fMu1Iso > fCuts[fChan]->m1iso);
  fGoodM2Iso        = (fMu2Iso > fCuts[fChan]->m2iso);
  fGoodDocaTrk      = (fCandDocaTrk > fCuts[fChan]->docatrk);
  fGoodLastCut      = true;

  fAnaCuts.update();

  fillRedTreeData();

  // -- to be consistent with the BDT traning
  ((TH1D*)fHistDir->Get("test3"))->Fill(1.);

  fPreselection = fGoodQ && fGoodMuonsEta && fWideMass && fPreselAlpha;
  if (0) cout << "fGoodQ = " << fGoodQ
	      << " fGoodMuonsEta = " << fGoodMuonsEta
	      << Form(" (%3.1f, %3.1f)", fMu1Eta, fMu2Eta)
	      << " fGoodMuonsPt = " << fGoodMuonsPt
	      << Form(" (%3.1f, %3.1f)", fMu1Pt, fMu2Pt)
	      << " fWideMass = " << fWideMass
	      << " fGoodAlpha = " <<  fGoodAlpha
	      << endl;
  if (fPreselection) ((TH1D*)fHistDir->Get("test3"))->Fill(2.);

  fPreselection = fPreselection && fGoodHLT1;
  if (fPreselection) ((TH1D*)fHistDir->Get("test3"))->Fill(3.);

}

// ----------------------------------------------------------------------
void candAna::fillCandidateHistograms(int offset) {
}


// ----------------------------------------------------------------------
int candAna::detChan(double m1eta, double m2eta) {
  int mode(1);
  // mode 0: channels 0 .. n are simply increasingly more forward regions for the most-forward muon
  // mode 1: channels 0 .. n are arbitrary eta regions for the most-forward muon combined with L1SEED requirements


  double m1 = TMath::Abs(m1eta);
  double m2 = TMath::Abs(m2eta);

  if (0 == mode) {
    int im1(-1), im2(-1);
    for (int ichan = 0; ichan < fNchan; ++ichan) {
      if ((m1 > fCuts[ichan]->metaMin) && (m1 < fCuts[ichan]->metaMax)) {
	im1 = ichan;
	break;
      }
    }
    for (int ichan = 0; ichan < fNchan; ++ichan) {
      if ((m2 > fCuts[ichan]->metaMin) && (m2 < fCuts[ichan]->metaMax)) {
	im2 = ichan;
	break;
      }
    }
    if ((im1 < 0) || (im2 < 0)) return -1;
    return (im1>im2?im1:im2);
  } else if (1 == mode) {
    if (m2 > m1) m1 = m2;
    bool found(false);
    for (int ichan = 0; ichan < fNchan; ++ichan) {
      if ((m1 > fCuts[ichan]->metaMin) && (m1 < fCuts[ichan]->metaMax)) {
	for (unsigned int is = 0; is < fCuts[ichan]->l1seeds.size(); ++is) {
	  found = false;
	  if (fL1Seeds & (0x1<<fCuts[ichan]->l1seeds[is])) {
	    found = true;
	    break;
	  }
	}
	return ichan;
      }
    }
    return -1;
  }
  return -1;
}

// ----------------------------------------------------------------------
void candAna::basicCuts() {
  cout << "    candAna basic cuts" << endl;
  fAnaCuts.addCut("fWideMass", "m(B candidate) [GeV]", fWideMass);
  fAnaCuts.addCut("fGoodHLT", "HLT", fGoodHLT1);
  fAnaCuts.addCut("fGoodMuonsID", "lepton ID", fGoodMuonsID);
  fAnaCuts.addCut("fGoodMuonsPt", "p_{T,#mu} [GeV]", fGoodMuonsPt);
  fAnaCuts.addCut("fGoodMuonsEta", "#eta_{#mu}", fGoodMuonsEta);
  fAnaCuts.addCut("fGoodTracks", "good tracks", fGoodTracks);
  fAnaCuts.addCut("fGoodTracksPt", "p_{T,trk} [GeV]", fGoodTracksPt);
  fAnaCuts.addCut("fGoodTracksEta", "#eta_{trk} ", fGoodTracksEta);
}


// ----------------------------------------------------------------------
void candAna::moreBasicCuts() {
  cout << "    candAna more basic cuts?" << endl;

}


// ----------------------------------------------------------------------
void candAna::candidateCuts() {
  cout << "    candAna candidate cuts" << endl;
  fAnaCuts.addCut("fGoodQ", "q_{1} 1_{2}", fGoodQ);
  fAnaCuts.addCut("fGoodPvAveW8", "<w8>", fGoodPvAveW8);
  fAnaCuts.addCut("fGoodPvLip", "LIP(PV)", fGoodPvLip);
  fAnaCuts.addCut("fGoodPvLipS", "LIPS(PV)", fGoodPvLipS);
  fAnaCuts.addCut("fGoodPv2Lip", "LIP(PV2)", fGoodPv2Lip);
  fAnaCuts.addCut("fGoodPv2LipS", "LIPS(PV2)", fGoodPv2LipS);
  fAnaCuts.addCut("fGoodIp", "IP", fGoodIp);
  fAnaCuts.addCut("fGoodIpS", "IPS", fGoodIpS);
  fAnaCuts.addCut("fGoodMaxDoca", "MAXDOCA", fGoodMaxDoca);
  fAnaCuts.addCut("fGoodPt", "p_{T,B}", fGoodPt);
  fAnaCuts.addCut("fGoodEta", "#eta_{B}", fGoodEta);
  fAnaCuts.addCut("fGoodAlpha", "#alpha", fGoodAlpha);
  fAnaCuts.addCut("fGoodFLS", "l/#sigma(l)", fGoodFLS);
  fAnaCuts.addCut("fGoodChi2", "#chi^{2}", fGoodChi2);
  fAnaCuts.addCut("fGoodIso", "I_{trk}", fGoodIso);
  fAnaCuts.addCut("fGoodCloseTrack", "close track veto", fGoodCloseTrack);
  fAnaCuts.addCut("fGoodDocaTrk", "d_{ca}(trk)", fGoodDocaTrk);
  fAnaCuts.addCut("fGoodLastCut", "lastCut", fGoodLastCut);
}


// ----------------------------------------------------------------------
void candAna::moreCandidateCuts() {
  cout << "    candAna more candidate cuts?" << endl;

}



// ----------------------------------------------------------------------
bool candAna::highPurity(TAnaTrack *pt) {

  if (pt->fTrackQuality & 4) return true;
  return false;


//   if (TRACKQUALITY > 0 && (0 == (pt->fTrackQuality & TRACKQUALITY))) {
//     if (fVerbose > 5) cout << "track " << pt->fIndex << " failed track quality: " << pt->fTrackQuality << endl;
//     return false;
//   }

//   if (TMath::Abs(pt->fTip) > TRACKTIP) {
//     if (fVerbose > 5) cout << "track " << pt->fIndex << " failed tip: " << pt->fTip
// 			   << " pointing to PV = "  << pt->fPvIdx  << endl;
//     return false;
//   }

//   if (TMath::Abs(pt->fLip) > TRACKLIP) {
//     if (fVerbose > 5) cout << "track " << pt->fIndex << " failed lip: " << pt->fLip
// 			   << " pointing to PV = "  << pt->fPvIdx  << endl;
//     return false;
//   }

  return true;
}


// ----------------------------------------------------------------------
void candAna::genMatch() {
  cout << "candAna::genMatch()  wrong function" << endl;
}

// ----------------------------------------------------------------------
void candAna::genMatchOld() {
  cout << "candAna::genMatchOld()  wrong function" << endl;
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
void candAna::efficiencyCalculation() {
  cout << "candAna::efficiencyCalculation()  wrong function" << endl;
}

// ----------------------------------------------------------------------
void candAna::triggerHLT() {
  fGoodHLT1 = false;
  fhltType = -1;
  fHLT1Path = "nada";

  TString a;
  int ps(0);
  bool result(false), wasRun(false), error(false);

  // -- NOTRIGGER, just accept the event
  if (HLTRANGE.begin()->first == "NOTRIGGER") {
    if (fVerbose>2) cout << "NOTRIGGER requested... " << endl;
    fGoodHLT1 = true;
    fHltPrescale = 1;
    return;
  }

  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i];
    ps = fpEvt->fHLTPrescale[i];
    wasRun = fpEvt->fHLTWasRun[i];
    result = fpEvt->fHLTResult[i];
    error  = fpEvt->fHLTError[i];

    if (wasRun && result) { // passed
      if (fVerbose>1  || (-32 == fVerbose) ) cout << "passed: " << a << endl;

      bool good = false;
      string spath;
      int rmin, rmax;
      for (map<string, pair<int, int> >::iterator imap = HLTRANGE.begin(); imap != HLTRANGE.end(); ++imap) {
	spath = imap->first;
	rmin = imap->second.first;
	rmax = imap->second.second;
	if (!a.CompareTo(imap->first.c_str()) && (rmin <= fRun) && (fRun <= rmax) ) {
	  // cout << "----------------------------------------------------------------------" << endl;
	  // cout << a << " result = " << result << " prescale: " << ps << " wasRun = " << wasRun
	  //      << " event: " << fEvt << " run: " << fRun
	  //      << endl;
	  good = true;
	  if (fVerbose > 1 || -32 == fVerbose  )
	    cout << "exact match: " << imap->first.c_str() << " HLT: " << a
		 << " result: " << result << endl;
	  break;
	}
	if (a.Contains(spath.c_str()) && (rmin <= fRun) && (fRun <= rmax)) {
	  good = true;
	  if (fVerbose > 1 || -32 == fVerbose)
	    cout << "close match: " << imap->first.c_str() << " HLT: " << a
		 << " result: " << result << " in run " << fRun << endl;
	  break;
	}
      }
      if (good) {
	fHltPrescale = ps;
	fHLT1Path    = a;
	fGoodHLT1    = true;
	return;
      }
    }
  }
}

// ----------------------------------------------------------------------
void candAna::triggerL1T() {
  fL1Seeds = 0;
  fL1SeedString = "";
  for (int i = 0; i < NL1T; ++i) {
    if (!fpEvt->fL1TResult[i]) continue;
    if ("L1_DoubleMu0er1p6_dEtaMax1p8" == fpEvt->fL1TNames[i]) {
      fL1Seeds |= 0x1;
      fL1SeedString += fpEvt->fL1TNames[i];
      fL1SeedString += " ";
      continue;
    } else if ("L1_DoubleMu0er1p6_dEta_Max1p8_OS" == fpEvt->fL1TNames[i]) {
      fL1Seeds |= 0x1<<1;
      fL1SeedString += fpEvt->fL1TNames[i];
      fL1SeedString += " ";
      continue;
    } else if ("L1_DoubleMu0er1p4_dEta_Max1p8_OS" == fpEvt->fL1TNames[i]) {
      fL1Seeds |= 0x1<<2;
      fL1SeedString += fpEvt->fL1TNames[i];
      fL1SeedString += " ";
      continue;
    } else  if ("L1_DoubleMu_10_0_dEta_Max1p8" == fpEvt->fL1TNames[i]) {
      fL1Seeds |= 0x1<<3;
      fL1SeedString += fpEvt->fL1TNames[i];
      fL1SeedString += " ";
      continue;
    } else if ("L1_DoubleMu_11_4" == fpEvt->fL1TNames[i]) {
      fL1Seeds |= 0x1<<4;
      fL1SeedString += fpEvt->fL1TNames[i];
      fL1SeedString += " ";
      continue;
    } else if ("L1_DoubleMu_12_5" == fpEvt->fL1TNames[i]) {
      fL1Seeds |= 0x1<<5;
      fL1SeedString += fpEvt->fL1TNames[i];
      fL1SeedString += " ";
      continue;
    } else if ("L1_DoubleMu_13_6" == fpEvt->fL1TNames[i]) {
      fL1Seeds |= 0x1<<6;
      fL1SeedString += fpEvt->fL1TNames[i];
      fL1SeedString += " ";
      continue;
    }
  }

}


// ----------------------------------------------------------------------
void candAna::triggerSelection() {
  const bool skipL12 = true;
  fGoodHLT = false;
  fhltType = -1;
  fHLTPath = "";
  hltObjMap.clear(); // reset the trig object map

  TString a;
  int ps(0);
  bool result(false), wasRun(false), error(false);

  // NOTRIGGER, just accept the event
  if ( HLTRANGE.begin()->first == "NOTRIGGER" ) {
    if (fVerbose>2) cout << "NOTRIGGER requested... " << endl;
    fGoodHLT = true;
    return;
  }

  // NOT USED (but maybe could be usefull in the future)
  // Any trigger fired (we should check it versus the valid hlt list for this data set
  //if ( HLTRANGE.begin()->first == "ALLTRIGGER") { // extend to ALLTRIGGER 16/1/13 d.k.
  //if (fVerbose > 2) cout << "ALLTRIGGER requested... " << endl;
    // For data make always true, the event has been triggered anyway
    //if (!fIsMC) {fGoodHLT = true;}
    //else { // for MC check if there was any HLT firing
    // for (int i = 0; i < NHLT; ++i) {
    //   result = wasRun = false;
    //   //a = fpEvt->fHLTNames[i];
    //   //ps = fpEvt->fHLTPrescale[i];
    //   wasRun = fpEvt->fHLTWasRun[i];
    //   result = fpEvt->fHLTResult[i];
    //   //error  = fpEvt->fHLTError[i];
    //   if (wasRun && result) {fGoodHLT=true; break;}
    // } // end HLT for loop
    //} // if-else
  //} // ALLTRIGGER

  bool pdTrigger = false;
  // Any trigger fired (we should check it versus the valid hlt list for this data set
  if ( HLTRANGE.begin()->first == "PDTRIGGER") { //
    if (fVerbose > 2) cout << "PDTRIGGER requested... " << endl;
    fpReader->pdTrigger()->setHLTKey(fRun, fpReader->getFile());  // needed for the new way
    pdTrigger=true;
  } // PDTRIGGER

  // Check L1, not really used
  if ( (fVerbose>9) || (fVerbose == -31)) {
    cout << "--------------------  L1: " << endl;
    for (int i = 0; i < NL1T; ++i) {
      result = wasRun = error = false;
      a = fpEvt->fL1TNames[i];
      ps = fpEvt->fL1TPrescale[i];
      result = fpEvt->fL1TResult[i];
      error  = fpEvt->fL1TMask[i];
      //if (a.Contains("Mu")) {
      //      if (result) {
	cout << a <<  " mask: " << error << " result: " << result << " ps: " << ps << endl;
	//      }
    }
  }

  // Check HLT
  // For every passed HLT look for a matching tigger from our list.
  // If if it confirmed by our list than match it with an object in the TrgObjv2 list
  // Mark the TrigObjv2 object my add a large number to the index.
  // Like this it can be recogised in the track match search.
  if ( (fVerbose>9) || (fVerbose==-32) ) cout<<" event "<<fEvt<<endl;
  bool isMuonTrigger=false; // just for diagnostics
  int foundNumHltObjects=0;
  int foundNumHlts=0;
  for (int i = 0; i < NHLT; ++i) {
    result = wasRun = error = false;
    a = fpEvt->fHLTNames[i];
    ps = fpEvt->fHLTPrescale[i];
    wasRun = fpEvt->fHLTWasRun[i];
    result = fpEvt->fHLTResult[i];
    error  = fpEvt->fHLTError[i];

    if (wasRun && result) { // passed

      if (fVerbose>1  || (-32 == fVerbose) ) cout << "passed: " << a << endl;
      if ((a == "digitisation_step")
	  || (a == "L1simulation_step")
	  || (a == "digi2raw_step")
	  || (a == "HLTriggerFinalPath")
	  || (a == "raw2digi_step")
	  || (a == "reconstruction_step")
	  ) {
	//cout<<" does this ever happen? " <<a<<endl;
	continue; // skip, go to the next triggr
      }

      // loop over our list of HLTs and look for matching
      // We assume that an event can be matched to only one trigger from our list,
      // that is they have to be exclusive.
      bool good = false;
      if (pdTrigger) { // PDTRIGGER mode - accept all triggers which fire and are on the DS list

	// check that this trigger belongs to our DS
	bool rightDS = fpReader->pdTrigger()->triggerInPd(DSNAME, a.Data());
	if (fVerbose>9) cout<<" check hlt-path "<<a.Data()<<" DS name "<<DSNAME<<" included? "<<rightDS<<endl;
	if (rightDS) { // hlt_path in this DS
	  // check if this is a L1/L2 type trigger or L3, assume that L1/L2 is always in the definition
	  bool isL1L2 = a.Contains("L1") || a.Contains("L2");
	  if (skipL12 && isL1L2) {
	    if (fVerbose>-1) cout<<" HIT-path os L1/L2 type, skip "<<a.Data()<<endl;
	    continue;
	  } else good=true; // accept, assume it is L3

	} else  { // hlt_path not in this DS
	  if (fVerbose>1)
	    cout<<" HIT-path not in this DS, skip it: "<<a.Data()<<" "<<rightDS<<" DS name: "<<DSNAME<<endl;
	  continue;
	} // end if

      } else { // TRIGRANGE mode - select only the trigger from our list

	string spath;
	int rmin, rmax;
	for (map<string, pair<int, int> >::iterator imap = HLTRANGE.begin();
	     imap != HLTRANGE.end(); ++imap) {
	  spath = imap->first;
	  rmin = imap->second.first;
	  rmax = imap->second.second;
	  if (!a.CompareTo(imap->first.c_str())) {
	    good=true;
	    if (fVerbose > 1 || -32 == fVerbose  )
	      cout << "exact match: " << imap->first.c_str() << " HLT: " << a
		   << " result: " << result << endl;
	    break; // can we skip the rest?
	  }
	  if (a.Contains(spath.c_str()) && (rmin <= fRun) && (fRun <= rmax)) {
	    good=true;
	    if (fVerbose > 1 || -32 == fVerbose)
	      cout << "close match: " << imap->first.c_str() << " HLT: " << a
		   << " result: " << result << " in run " << fRun << endl;
	    break; // can we skip the rest?
	  } // end if

	} // end for loop
      } // if pdTrigger

      if (good) {  // for matched hlt paths select the trigger object
	// this trigger matched one in our list
	foundNumHlts++;
	fGoodHLT = true;
	fhltType = i;
	fHLTPath = a.Data();
	isMuonTrigger = a.Contains("Mu") || a.Contains("mu") || a.Contains("MU");
	//if (ps!=1) cout<<"prescale not one "<<a.Data()<<" "<<ps<<endl;

	continue;

	bool foundHltObject = false;
	int countModules=0, lastIndex=-1;
	TTrgObjv2 *pTO;
	if ( (fVerbose>9) || (fVerbose==-32))
	  cout<<" TTrgObjv2 objects, size= "<<fpEvt->nTrgObjv2()<<endl;
	for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {  // loop over all saved hlt objects
	  pTO = fpEvt->getTrgObjv2(i);

	  if (a == pTO->fHltPath) { // found the right one, matched the hlt-name
	    // this trig object matches a passed and selected trigger
	    foundHltObject=true;
	    foundNumHltObjects++;
	    countModules++;
            lastIndex=i;
	    int hltIndex = pTO->fHltIndex; // HLT path index

            // mark all module as selected by adding a large number, so >1000.
            // do it either here (ALL MODULES)  or below (LAST MODULE)
	    hltObjMap[i]= (hltIndex & 0x7FFFFFFF); // make sure the uper bit is free

	    vector<int> muonIndex = pTO->fIndex;
	    vector<int> muonID = pTO->fID;
	    vector<TLorentzVector> muonP = pTO->fP;
	    int num = muonIndex.size();

	    if ( (fVerbose>9) || (fVerbose==-32)) {
	      cout<<" matched: "<<pTO->fHltPath<<" hlt-index: "<<pTO->fHltIndex<<" module label: "
		  <<pTO->fLabel<<" type: "<<pTO->fType<<" num of particles: "<<num<<endl;
	      pTO->dump();
	    }

	  } // if matched
	} // end for loop

	// only mark the last module for matching, do it by setting the highest bit
        if (lastIndex>-1) {
          pTO = fpEvt->getTrgObjv2(lastIndex);
	  hltObjMap[lastIndex]= (hltObjMap[lastIndex] | 0x80000000); // set top bit for final
        }

	if (!foundHltObject && isMuonTrigger)
	  cout<<"Warning: canAna::triggerSelection: matching trigger module not found! "
	      <<a<<endl;
	} // end if fGoodHLT
    } // if passed
  } // end for loop hlt

  // This is just to mark the trigger type and number of triggers for this event
  // the type only valid for the last type if there are more than 1
  if (fhltType>999) {cout<<fhltType<<endl; fhltType -= 1000;}
  fhltType = fhltType + (1000 * (foundNumHlts-1));

  // Diagnostics printout
  if ( (fVerbose>9) || (fVerbose==-32))
    cout<<" number of found matching hlt objects: "<<foundNumHltObjects<<endl;

  // Diagnostics printout
  if (pdTrigger && isMuonTrigger && (foundNumHltObjects==0)) {
    cout<<" No matching hlt objects found:  "<<endl;
    //cout<< " Passed triggers: ";
    for (int i = 0; i < NHLT; ++i) {
      result = wasRun = false;
      a = fpEvt->fHLTNames[i];
      //ps = fpEvt->fHLTPrescale[i];
      wasRun = fpEvt->fHLTWasRun[i];
      result = fpEvt->fHLTResult[i];
      //error  = fpEvt->fHLTError[i];
      if (wasRun && result) cout << a << " ";
    }
    cout<<endl;
  }

  // TESTS
  if (fVerbose>999) {  // Just testing
    cout << " ----------------------------------------------------------------------" << endl;

    TTrgObjv2 *pTO;
    cout<<" Dump TTrgObjv2 "<<fpEvt->nTrgObjv2()<<endl;
    for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {
      pTO = fpEvt->getTrgObjv2(i);
      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      cout<<i<<" hlt "<<pTO->fHltPath<<" hlt-index "<<pTO->fHltIndex<<" module label "
	  <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<" "<<num<<endl;

      for(int n=0;n<num;++n) {
	int index = muonIndex[n];
	int id = muonID[n];
	TLorentzVector p = muonP[n];
	cout<<n<<" index "<<index<<" id "<<id<<" pt/eta/phi "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;
      }

    }

    // print the hlt object map
    for(map<unsigned int, unsigned int, less<unsigned int> >::iterator iter=hltObjMap.begin();
	iter!=hltObjMap.end(); ++iter) {
      cout<<hex<<iter->first<<" "<<iter->second<<" "<<dec<<endl;
    }

  } // end testing

  // if (0) {  THIS IS SOMETHING OLD FROM URS?
  //   TTrgObj *p;
  //   for (int i = 0; i < fpEvt->nTrgObj(); ++i) {
  //     p = fpEvt->getTrgObj(i);
  //     //      cout << p->fLabel << endl;
  //     //      if (!p->fLabel.CompareTo("hltL1sL1DoubleMu33HighQ:HLT::")) cout << "= " << p->fLabel << endl;
  //     if (!p->fLabel.CompareTo("hltL1sL1DoubleMuOpen:HLT::")) {
  //       if (0) cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta()
  // 		    << " phi = " << p->fP.Phi() <<  endl;
  //     }

  //     if (1 && !p->fLabel.CompareTo("hltL1sL1DoubleMu33HighQ:HLT::")) {
  //       //cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
  //       ((TH1D*)fHistDir->Get("L1_0"))->Fill(p->fP.Eta());
  //     }

  //     if (1 && !p->fLabel.CompareTo("hltL1sL1DoubleMu0or33HighQ:HLT::")) {
  //       //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
  //       ((TH1D*)fHistDir->Get("L1_1"))->Fill(p->fP.Eta());
  //     }

  //     if (1 && !p->fLabel.CompareTo("hltDimuon33L1Filtered0:HLT::")) {
  //       //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
  //       ((TH1D*)fHistDir->Get("L1_2"))->Fill(p->fP.Eta());
  //     }

  //     if (1 && !p->fLabel.CompareTo("hltDimuonL2PreFiltered0:HLT::")) {
  //       //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
  //       ((TH1D*)fHistDir->Get("L1_3"))->Fill(p->fP.Eta());
  //     }

  //     if (1 && !p->fLabel.CompareTo("hltDimuonL2PreFiltered0:HLT::")) {
  //       //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
  //       ((TH1D*)fHistDir->Get("L1_4"))->Fill(p->fP.Eta());
  //     }

  //     if (1 && !p->fLabel.CompareTo("hltDoubleDisplacedMu4L3PreFiltered:HLT::")) {
  //       //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
  //       ((TH1D*)fHistDir->Get("L1_5"))->Fill(p->fP.Eta());
  //     }

  //     if (1 && !p->fLabel.CompareTo("hltDoubleMu3p5LowMassDisplacedL3Filtered:HLT::")) {
  //       //      cout << p->fLabel << " pT = " << p->fP.Perp() << " eta = " << p->fP.Eta() << " phi = " << p->fP.Phi() <<  endl;
  //       ((TH1D*)fHistDir->Get("L1_6"))->Fill(p->fP.Eta());
  //     }
  //   }

  // }

  if (false == fGoodHLT) {
    if (fVerbose > 1) cout << "------->  event NOT triggered!" << endl;
  }

}

// ----------------------------------------------------------------------
void candAna::bookHist() {

  fHistDir->cd();

  TH1D *h11(0);
  (void)h11;
  h11 = new TH1D(Form("mon%s", fName.c_str()), Form("mon%s", fName.c_str()), 50, 0., 50.);

  h11 = new TH1D(Form("dr_%s", fName.c_str()), Form("dr(track, trigger bject) %s", fName.c_str()), 501, -0.001, 0.5);

  h11 = new TH1D("l1pt",  "(pT(L3) - pT(L1))/(pT(L3)", 40, -2., 2.);
  h11 = new TH1D("l1pt_10_11",  "(pT(L3) - pT(L1))/(pT(L3) for 10 < pT < 11", 40, -2., 2.);
  h11 = new TH1D("l1pt_15_17",  "(pT(L3) - pT(L1))/(pT(L3) for 15 < pT < 17", 40, -2., 2.);
  h11 = new TH1D("l1ptec_5_6",  "(pT(L3) - pT(L1))/(pT(L3) for 5 < pT < 6 (endcap)", 40, -2., 2.);
  h11 = new TH1D("l1ptec_10_11",  "(pT(L3) - pT(L1))/(pT(L3) for 10 < pT < 11 (endcap)", 40, -2., 2.);

  h11 = new TH1D("l1eta", "(eta(L3) - eta(L1))", 40, -2., 2.);
  h11 = new TH1D("l1phi", "(phi(L3) - phi(L1))", 40, -2., 2.);
  new TProfile("pl1pt","Profile of (pT(L3) - pT(L1))/pT(L3)", 30, 0., 30., -2., 2., "S");
  new TProfile("pl1eta","Profile of L1 eta resolution", 25, -2.5, 2.5, -2., 2., "S");
  new TProfile("pl1phi","Profile of L1 phi resolution", 32, -3.15, 3.15, -2., 2., "S");
  new TProfile("pl1ptec","Profile of (pT(L3) - pT(L1))/pT(L3) (endcap)", 30, 0., 30., -2., 2., "S");


  h11 = new TH1D("L1_0", "hltL1sL1DoubleMu33HighQ", 50, -2.5, 2.5);
  h11 = new TH1D("L1_1", "hltL1sL1DoubleMu0or33HighQ", 50, -2.5, 2.5);
  h11 = new TH1D("L1_2", "hltDimuon33L1Filtered0", 50, -2.5, 2.5);
  h11 = new TH1D("L1_3", "hltDimuon33L2PreFiltered0", 50, -2.5, 2.5);
  h11 = new TH1D("L1_4", "hltDimuonL2PreFiltered0", 50, -2.5, 2.5);
  h11 = new TH1D("L1_5", "hltDoubleDisplacedMu4L3PreFiltered", 50, -2.5, 2.5);
  h11 = new TH1D("L1_6", "hltDoubleMu3p5LowMassDisplacedL3Filtered", 50, -2.5, 2.5);

  const double firstRun = 251000.5, lastRun = 271000.5; const int numRuns = 20000;
  h11 = new TH1D("run1", "runs json&hltmatched", numRuns, firstRun, lastRun);

  //h11 = new TH1D("test1", "test1",50, 0., 50.);
  h11 = new TH1D("test2", "test2",1000, 0., 0.1);
  h11 = new TH1D("test3", "test3",50, 0., 50.);
  //h11 = new TH1D("test4", "test4",100, 0., 100.);
  //h11 = new TH1D("test5", "test5",100, 0., 100.);
  h11 = new TH1D("test6", "test6",200, 0., 2.);
  h11 = new TH1D("test7", "test7", 400, 0., 4.);
  h11 = new TH1D("test8", "test8",200, -1., 1.);

  h11 = new TH1D("test11", "dr", 1000, 0., 1.);
  h11 = new TH1D("test12", "dr", 1000, 0., 1.);
  h11 = new TH1D("test13", "dpt", 200, -1., 1.);
  h11 = new TH1D("test14", "dpt", 200, -1., 1.);

  h11 = new TH1D("gp1cms", "p1cms", 50, 0, 10.);
  h11 = new TH1D("gp2cms", "p2cms", 50, 0, 10.);
  h11 = new TH1D("gt1cms", "t1cms", 50, -1, 1.);
  h11 = new TH1D("gt2cms", "t2cms", 50, -1, 1.);

  h11 = new TH1D("gp1cmsg", "p1cms (with photons)", 50, 0, 10.);
  h11 = new TH1D("gp2cmsg", "p2cms (with photons)", 50, 0, 10.);
  h11 = new TH1D("gt1cmsg", "t1cms (with photons)", 50, 0, 1.);
  h11 = new TH1D("gt2cmsg", "t2cms (with photons)", 50, 0, 1.);

  h11 = new TH1D("rp1cms", "p1cms", 50, 0, 10.);
  h11 = new TH1D("rp2cms", "p2cms", 50, 0, 10.);
  h11 = new TH1D("rt1cms", "t1cms", 50, -1, 1.);
  h11 = new TH1D("rt2cms", "t2cms", 50, -1, 1.);
  h11 = new TH1D("rt3cms", "t3cms", 50, -1, 1.);
  TH2D *h22(0);
  (void)h22;
  h22 = new TH2D("tvsm",   "tvsm", 50, 4.9, 5.9, 50, -1., 1.);
  h22 = new TH2D("h2dtest",   "h2dtest", 50, 0.0, 0.05, 50, 0.,0.5);

  h11 = new TH1D("rp1cmsg", "p1cms (with photons)", 50, 0, 10.);
  h11 = new TH1D("rp2cmsg", "p2cms (with photons)", 50, 0, 10.);
  h11 = new TH1D("rt1cmsg", "t1cms (with photons)", 50, 0, 1.);
  h11 = new TH1D("rt2cmsg", "t2cms (with photons)", 50, 0, 1.);
  h11 = new TH1D("rt3cmsg", "t2cms (with photons)", 50, 0, 1.);

  h11 = new TH1D("gt1", "gt1", 50, -1, 1.);
  h11 = new TH1D("gt2", "gt2", 50, -1, 1.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  setupReducedTree(fTree);

  // -- tree for muon id MVA
  if (0) {
    fMuonIdTree = new TTree("muonidtree", "muonidtree");
    setupMuonIdTree(fMuonIdTree);
  }

  // -- Efficiency/Acceptance Tree
  fEffTree = new TTree("effTree", "effTree");
  fEffTree->Branch("run",    &fRun,               "run/L");
  fEffTree->Branch("evt",    &fEvt,               "evt/L");
  fEffTree->Branch("hlt",    &fGoodHLT,           "hlt/O");
  fEffTree->Branch("hlt1",   &fGoodHLT1,          "hlt1/O");
  fEffTree->Branch("procid", &fProcessType,       "procid/I");
  fEffTree->Branch("bidx",   &fGenBTmi,           "bidx/I");

  fEffTree->Branch("gm",     &fETgm,              "gm/F");
  fEffTree->Branch("gpt",    &fETgpt,             "gpt/F");
  fEffTree->Branch("geta",   &fETgeta,            "geta/F");
  fEffTree->Branch("gtau",   &fETgtau,            "gtau/F");

  fEffTree->Branch("m1pt",   &fETm1pt,            "m1pt/F");
  fEffTree->Branch("g1pt",   &fETg1pt,            "g1pt/F");
  fEffTree->Branch("m1eta",  &fETm1eta,           "m1eta/F");
  fEffTree->Branch("g1eta",  &fETg1eta,           "g1eta/F");
  fEffTree->Branch("m1q",    &fETm1q,             "m1q/I");
  fEffTree->Branch("m1gt",   &fETm1gt,            "m1gt/O");
  fEffTree->Branch("m1id",   &fETm1id,            "m1id/O");
  fEffTree->Branch("m1tmid", &fETm1tmid,          "m1tmid/O");
  fEffTree->Branch("m1mvaid",&fETm1mvaid,         "m1mvaid/O");

  fEffTree->Branch("m2pt",   &fETm2pt,            "m2pt/F");
  fEffTree->Branch("g2pt",   &fETg2pt,            "g2pt/F");
  fEffTree->Branch("m2eta",  &fETm2eta,           "m2eta/F");
  fEffTree->Branch("g2eta",  &fETg2eta,           "g2eta/F");
  fEffTree->Branch("m2q",    &fETm2q,             "m2q/I");
  fEffTree->Branch("m2gt",   &fETm2gt,            "m2gt/O");
  fEffTree->Branch("m2id",   &fETm2id,            "m2id/O");
  fEffTree->Branch("m2tmid", &fETm2tmid,          "m2tmid/O");
  fEffTree->Branch("m2mvaid",&fETm2mvaid,         "m2mvaid/O");

  fEffTree->Branch("m",      &fETcandMass,        "m/F");
  fEffTree->Branch("tau",    &fETtau,             "tau/F");


  // -- Analysis distributions
  TH1D *h = new TH1D("analysisDistributions", "analysisDistributions", 10000, 0., 10000.);
  (void)h;

}


// ----------------------------------------------------------------------
void candAna::setupMuonIdTree(TTree *t) {

  t->Branch("pt",                   &fMuonData.pt,  "pt/F");
  t->Branch("eta",                  &fMuonData.eta, "eta/F");

  t->Branch("intvalidmuonhits",     &fMuonData.validMuonHits, "validmuonhits/I");
  t->Branch("intnmatchedstations",  &fMuonData.nMatchedStations, "nmatchedstations/I");
  t->Branch("intvalidpixelhits",    &fMuonData.validPixelHits, "validpixelhits/I");
  t->Branch("inttrklayerswithhits", &fMuonData.trkLayerWithHits, "trklayerswithhits/I");

  t->Branch("gchi2",                &fMuonData.glbNChi2, "gchi2/F");

  t->Branch("itrkvalidfraction",    &fMuonData.trkValidFract, "itrkvalidfraction/F");
  t->Branch("segcomp",              &fMuonData.segComp, "segcomp/F");
  t->Branch("chi2lmom",             &fMuonData.chi2LocMom, "chi2lmom/F");
  t->Branch("chi2lpos",             &fMuonData.chi2LocPos, "chi2lpos/F");
  t->Branch("gtrkprob",             &fMuonData.glbTrackProb, "gtrkprob/F");
  t->Branch("ntrkvhits",            &fMuonData.NTrkVHits, "ntrkvhits/F");
  t->Branch("ntrkehitsout",         &fMuonData.NTrkEHitsOut, "ntrkehitsout/F");

  t->Branch("kink",                 &fMuonData.kink, "kink/F");
  t->Branch("dpt",                  &fMuonData.dpt, "dpt/F");
  t->Branch("dptrel",               &fMuonData.dptrel, "dptrel/F");
  t->Branch("deta",                 &fMuonData.deta, "deta/F");
  t->Branch("dphi",                 &fMuonData.dphi, "dphi/F");
  t->Branch("dr",                   &fMuonData.dr, "dr/F");

}

// ----------------------------------------------------------------------
void candAna::setupReducedTree(TTree *t) {

  t->Branch("run",     &fRun,               "run/L");
  t->Branch("evt",     &fEvt,               "evt/L");
  t->Branch("ls",      &fLS,                "ls/I");
  t->Branch("rlumi",   &fLumi,              "rlumi/D");
  t->Branch("json",    &fJSON,              "json/O");
  t->Branch("tm",      &fCandTM,            "tm/I");
  t->Branch("pr",      &fGenBpartial,       "pr/I");
  t->Branch("procid",  &fProcessType,       "procid/I");
  t->Branch("hlt",     &fGoodHLT,           "hlt/O");
  t->Branch("hlt1",    &fGoodHLT1,          "hlt1/O");
  t->Branch("l1s",     &fL1Seeds,           "l1s/I");
  t->Branch("ps",      &fHltPrescale,       "ps/I");
  t->Branch("tos",     &fTOS,               "tos/O");
  t->Branch("tis",     &fTIS,               "tis/O");
  t->Branch("reftrg",  &fRefTrigger,        "reftrg/O");
  t->Branch("pvidx",   &fPvIdx,             "pvidx/I");
  t->Branch("pvz",     &fPvZ,               "pvz/D");
  t->Branch("pvn",     &fPvN,               "pvn/I");
  t->Branch("cb",      &fCowboy,            "cb/O");
  t->Branch("rr",      &fRunRange,          "rr/I");
  t->Branch("bdt",     &fBDT,               "bdt/D");
  t->Branch("npv",     &fPvN,               "npv/I");
  t->Branch("pvw8",    &fPvAveW8,           "pvw8/D");
  t->Branch("presel",  &fPreselection,      "presel/O");

  // -- global cuts and weights
  t->Branch("gmuid",   &fGoodMuonsID,       "gmuid/O");
  t->Branch("gmutmid", &fGoodMuonsTmID,     "gmutmid/O");
  t->Branch("gmumvaid",&fGoodMuonsMvaID,    "gmumvaid/O");
  t->Branch("gmupt",   &fGoodMuonsPt,       "gmupt/O");
  t->Branch("gmueta",  &fGoodMuonsEta,      "gmueta/O");
  t->Branch("gtqual",  &fGoodTracks,        "gtqual/O");
  t->Branch("gtpt",    &fGoodTracksPt,      "gtpt/O");
  t->Branch("gteta",   &fGoodTracksEta,     "gteta/O");

  // -- PV
  t->Branch("pvlip",    &fCandPvLip,        "pvlip/D");
  t->Branch("pvlips",   &fCandPvLipS,       "pvlips/D");
  t->Branch("pv2lip",   &fCandPv2Lip,       "pv2lip/D");
  t->Branch("pv2lips",  &fCandPv2LipS,      "pv2lips/D");
  t->Branch("pvip",     &fCandPvIp,         "pvip/D");
  t->Branch("pvips",    &fCandPvIpS,        "pvips/D");
  t->Branch("pvip3d",   &fCandPvIp3D,       "pvip3d/D");
  t->Branch("pvips3d",  &fCandPvIpS3D,      "pvips3d/D");

  // -- cand
  t->Branch("chan",    &fChan,              "chan/I");
  t->Branch("q",       &fCandQ,             "q/I");
  t->Branch("type",    &fCandType,          "type/I");
  t->Branch("pt",      &fCandPt,            "pt/D");
  t->Branch("eta",     &fCandEta,           "eta/D");
  t->Branch("phi",     &fCandPhi,           "phi/D");
  t->Branch("tau",     &fCandTau,           "tau/D");
  t->Branch("taue",    &fCandTauE,          "taue/D");
  t->Branch("tauxy",   &fCandTauxy,         "tauxy/D");
  t->Branch("tauxye",  &fCandTauxyE,        "tauxye/D");
  t->Branch("m",       &fCandM,             "m/D");
  t->Branch("me",      &fCandME,            "me/D");
  t->Branch("cm",      &fCandM2,            "cm/D");
  t->Branch("m3",      &fCandM3,            "m3/D");
  t->Branch("m4",      &fCandM4,            "m4/D");
  t->Branch("cosa",    &fCandCosA,          "cosa/D");
  t->Branch("alpha",   &fCandA,             "alpha/D");
  t->Branch("iso",     &fCandIso,           "iso/D");
  t->Branch("isotrk",  &fCandIsoTrk,        "isotrk/I");
  t->Branch("closetrk",&fCandCloseTrk,      "closetrk/I");
  t->Branch("chi2",    &fCandChi2,          "chi2/D");
  t->Branch("dof",     &fCandDof,           "dof/D");
  t->Branch("chi2dof", &fCandChi2Dof,       "chi2dof/D");
  t->Branch("prob",    &fCandProb,          "prob/D");
  t->Branch("fls3d",   &fCandFLS3d,         "fls3d/D");
  t->Branch("fl3d",    &fCandFL3d,          "fl3d/D");
  t->Branch("flxy",    &fCandFLxy,          "flxy/D");
  t->Branch("fl3dE",   &fCandFL3dE,         "fl3dE/D");
  t->Branch("flsxy",   &fCandFLSxy,         "flsxy/D");
  t->Branch("docatrk", &fCandDocaTrk,       "docatrk/D");
  t->Branch("docatrkbdt", &fCandDocaTrkBdt, "docatrkbdt/D");
  t->Branch("maxdoca", &fCandDoca,          "maxdoca/D");
  t->Branch("lip",     &fCandPvLip,         "lip/D");
  t->Branch("lipE",    &fCandPvLipE,        "lipE/D");
  t->Branch("tip",     &fCandPvTip,         "tip/D");
  t->Branch("tipE",    &fCandPvTipE,        "tipE/D");
  t->Branch("dr",      &fDeltaR,            "dr/D");

  t->Branch("pvdchi2",   &fCandPvDeltaChi2, "pvdchi2/D");
  t->Branch("closetrks1", &fCandCloseTrkS1,   "closetrks1/I");
  t->Branch("closetrks2", &fCandCloseTrkS2,   "closetrks2/I");
  t->Branch("closetrks3", &fCandCloseTrkS3,   "closetrks3/I");
  t->Branch("othervtx", &fCandOtherVtx, "othervtx/D");

  t->Branch("osiso",   &fOsIso,             "osiso/D");
  t->Branch("osreliso",&fOsRelIso,          "osreliso/D");
  t->Branch("osmpt",   &fOsMuonPt,          "osmpt/D");
  t->Branch("osmptrel",&fOsMuonPtRel,       "osmptrel/D");
  t->Branch("osmdr",   &fOsMuonDeltaR,      "osmdr/D");

  // -- muons
  t->Branch("m1q",     &fMu1Q,              "m1q/I");
  t->Branch("m1id",    &fMu1Id,             "m1id/O");
  t->Branch("m1gmid",  &fMu1GmId,           "m1gmid/O");
  t->Branch("m1tmid",  &fMu1TmId,           "m1tmid/O");
  t->Branch("m1mvaid", &fMu1MvaId,          "m1mvaid/O");
  t->Branch("m1rtmid", &fMu1rTmId,          "m1rtmid/O");
  t->Branch("m1rmvaid",&fMu1rMvaId,          "m1rmvaid/O");
  t->Branch("m1mvabdt",&fMu1BDT,            "m1mvabdt/D");
  t->Branch("m1rmvabdt",&fMu1rBDT,          "m1rmvabdt/D");

  t->Branch("m1mvaidlm",   &fMu1MvaIdLM,    "m1mvaidlm/O");
  t->Branch("m1rmvaidlm",  &fMu1rMvaIdLM,   "m1rmvaidlm/O");
  t->Branch("m1mvabdtlm",  &fMu1BDTLM,      "m1mvabdtlm/D");
  t->Branch("m1rmvabdtlm", &fMu1rBDTLM,     "m1rmvabdtlm/D");

  t->Branch("m1trigm", &fMu1TrigM,          "m1trigm/D");

  t->Branch("m1pt",    &fMu1Pt,             "m1pt/D");
  t->Branch("m1eta",   &fMu1Eta,            "m1eta/D");
  t->Branch("m1phi",   &fMu1Phi,            "m1phi/D");
  t->Branch("m1ip",    &fMu1IP,             "m1ip/D");
  t->Branch("m1ips",   &fMu1IPS,            "m1ips/D");
  t->Branch("m1gt",    &fMu1TkQuality,      "m1gt/I");
  t->Branch("m1pix",   &fMu1Pix,            "m1pix/I");
  t->Branch("m1bpix",  &fMu1BPix,           "m1bpix/I");
  t->Branch("m1bpixl1",&fMu1BPixL1,         "m1bpixl1/I");
  t->Branch("m1chi2",  &fMu1Chi2,           "m1chi2/D");
  t->Branch("m1pv",    &fMu1PV,             "m1pv/I");
  t->Branch("m1vtxprob",&fMu1VtxProb,       "m1vtxprob/D");
  t->Branch("m1xpdist",&fMu1XpDist,         "m1xpdist/D");
  t->Branch("m1iso",   &fMu1Iso,            "m1iso/D");

  t->Branch("m2q",     &fMu2Q,              "m2q/I");
  t->Branch("m2id",    &fMu2Id,             "m2id/O");
  t->Branch("m2gmid",  &fMu2GmId,           "m2gmid/O");
  t->Branch("m2tmid",  &fMu2TmId,           "m2tmid/O");
  t->Branch("m2mvaid", &fMu2MvaId,          "m2mvaid/O");
  t->Branch("m2rtmid", &fMu2rTmId,          "m2rtmid/O");
  t->Branch("m2rmvaid",&fMu2rMvaId,         "m2rmvaid/O");
  t->Branch("m2mvabdt",&fMu2BDT,            "m2mvabdt/D");
  t->Branch("m2rmvabdt",&fMu2rBDT,          "m2rmvabdt/D");

  t->Branch("m2mvaidlm",   &fMu2MvaIdLM,    "m2mvaidlm/O");
  t->Branch("m2rmvaidlm",  &fMu2rMvaIdLM,   "m2rmvaidlm/O");
  t->Branch("m2mvabdtlm",  &fMu2BDTLM,      "m2mvabdtlm/D");
  t->Branch("m2rmvabdtlm", &fMu2rBDTLM,     "m2rmvabdtlm/D");

  t->Branch("m2trigm", &fMu2TrigM,          "m2trigm/D");

  t->Branch("m2pt",    &fMu2Pt,             "m2pt/D");
  t->Branch("m2eta",   &fMu2Eta,            "m2eta/D");
  t->Branch("m2phi",   &fMu2Phi,            "m2phi/D");
  t->Branch("m2ip",    &fMu2IP,             "m2ip/D");
  t->Branch("m2ips",   &fMu2IPS,            "m2ips/D");
  t->Branch("m2gt",    &fMu2TkQuality,      "m2gt/I");
  t->Branch("m2pix",   &fMu2Pix,            "m2pix/I");
  t->Branch("m2bpix",  &fMu2BPix,           "m2bpix/I");
  t->Branch("m2bpixl1",&fMu2BPixL1,         "m2bpixl1/I");
  t->Branch("m2chi2",  &fMu2Chi2,           "m2chi2/D");
  t->Branch("m2pv",    &fMu2PV,             "m2pv/I");
  t->Branch("m2vtxprob",&fMu2VtxProb,       "m2vtxprob/D");
  t->Branch("m2xpdist",&fMu2XpDist,         "m2xpdist/D");
  t->Branch("m2iso",   &fMu2Iso,            "m2iso/D");

  t->Branch("mudist",  &fMuDist,            "mudist/D");
  t->Branch("mudeltar",&fMuDeltaR,          "mudeltar/D");
  t->Branch("hltm",    &fHLTmatch,          "hltm/O");
  t->Branch("hltt",    &fhltType,           "hltt/I");

  t->Branch("g1pt",    &fMu1PtGen,          "g1pt/D");
  t->Branch("g2pt",    &fMu2PtGen,          "g2pt/D");
  t->Branch("g1eta",   &fMu1EtaGen,         "g1eta/D");
  t->Branch("g2eta",   &fMu2EtaGen,         "g2eta/D");
  t->Branch("g1phi",   &fMu1PhiGen,         "g1phi/D");
  t->Branch("g2phi",   &fMu2PhiGen,         "g2phi/D");
  t->Branch("gmass",   &fGenMass,           "gmass/D");
  t->Branch("gtau",    &fGenLifeTime,       "gtau/D");
  t->Branch("g1id",    &fMu1GenID,          "g1id/I");
  t->Branch("g2id",    &fMu2GenID,          "g2id/I");

  t->Branch("t1pt",    &fMu1PtNrf,          "t1pt/D");
  t->Branch("t1eta",   &fMu1EtaNrf,         "t1eta/D");
  t->Branch("t2pt",    &fMu2PtNrf,          "t2pt/D");
  t->Branch("t2eta",   &fMu2EtaNrf,         "t2eta/D");

  t->Branch("hm1pt",  &fHltMu1Pt,  "hm1pt/D");
  t->Branch("hm1eta", &fHltMu1Eta, "hm1eta/D");
  t->Branch("hm1phi", &fHltMu1Phi, "hm1phi/D");
  t->Branch("hm2pt",  &fHltMu2Pt,  "hm2pt/D");
  t->Branch("hm2eta", &fHltMu2Eta, "hm2eta/D");
  t->Branch("hm2phi", &fHltMu2Phi, "hm2phi/D");


  // -- all muon variables for Marco
  t->Branch("m1bdt", &fMu1Data.mbdt, "m1bdt/F");
  t->Branch("m1validMuonHits", &fMu1Data.validMuonHits, "m1validMuonHits/F");
  t->Branch("m1glbNChi2", &fMu1Data.glbNChi2, "m1glbNChi2/F");

  t->Branch("m1nMatchedStations", &fMu1Data.nMatchedStations, "m1nMatchedStations/I");
  t->Branch("m1validPixelHits", &fMu1Data.validPixelHits, "m1validPixelHits/I");
  t->Branch("m1validPixelHits2", &fMu1Data.validPixelHits2, "m1validPixelHits2/I");
  t->Branch("m1trkLayerWithHits", &fMu1Data.trkLayerWithHits, "m1trkLayerWithHits/I");

  t->Branch("m1trkValidFract", &fMu1Data.trkValidFract, "m1trkValidFract/F");
  t->Branch("m1pt", &fMu1Data.pt, "m1pt/F");
  t->Branch("m1eta", &fMu1Data.eta, "m1eta/F");
  t->Branch("m1segComp", &fMu1Data.segComp, "m1segComp/F");
  t->Branch("m1chi2LocMom", &fMu1Data.chi2LocMom, "m1chi2LocMom/F");
  t->Branch("m1chi2LocPos", &fMu1Data.chi2LocPos, "m1chi2LocPos/F");
  t->Branch("m1glbTrackProb", &fMu1Data.glbTrackProb, "m1glbTrackProb/F");
  t->Branch("m1NTrkVHits", &fMu1Data.NTrkVHits, "m1NTrkVHits/F");
  t->Branch("m1NTrkEHitsOut", &fMu1Data.NTrkEHitsOut, "m1NTrkEHitsOut/F");
  t->Branch("m1kink", &fMu1Data.kink, "m1kink/F");
  t->Branch("m1dpt", &fMu1Data.dpt, "m1dpt/F");
  t->Branch("m1dptrel", &fMu1Data.dptrel, "m1dptrel/F");
  t->Branch("m1deta", &fMu1Data.deta, "m1deta/F");
  t->Branch("m1dphi", &fMu1Data.dphi, "m1dphi/F");
  t->Branch("m1dr", &fMu1Data.dr, "m1dr/F");


  t->Branch("m2bdt", &fMu2Data.mbdt, "m2bdt/F");
  t->Branch("m2validMuonHits", &fMu2Data.validMuonHits, "m2validMuonHits/F");
  t->Branch("m2glbNChi2", &fMu2Data.glbNChi2, "m2glbNChi2/F");

  t->Branch("m2nMatchedStations", &fMu2Data.nMatchedStations, "m2nMatchedStations/I");
  t->Branch("m2validPixelHits", &fMu2Data.validPixelHits, "m2validPixelHits/I");
  t->Branch("m2validPixelHits2", &fMu2Data.validPixelHits2, "m2validPixelHits2/I");
  t->Branch("m2trkLayerWithHits", &fMu2Data.trkLayerWithHits, "m2trkLayerWithHits/I");

  t->Branch("m2trkValidFract", &fMu2Data.trkValidFract, "m2trkValidFract/F");
  t->Branch("m2pt", &fMu2Data.pt, "m2pt/F");
  t->Branch("m2eta", &fMu2Data.eta, "m2eta/F");
  t->Branch("m2segComp", &fMu2Data.segComp, "m2segComp/F");
  t->Branch("m2chi2LocMom", &fMu2Data.chi2LocMom, "m2chi2LocMom/F");
  t->Branch("m2chi2LocPos", &fMu2Data.chi2LocPos, "m2chi2LocPos/F");
  t->Branch("m2glbTrackProb", &fMu2Data.glbTrackProb, "m2glbTrackProb/F");
  t->Branch("m2NTrkVHits", &fMu2Data.NTrkVHits, "m2NTrkVHits/F");
  t->Branch("m2NTrkEHitsOut", &fMu2Data.NTrkEHitsOut, "m2NTrkEHitsOut/F");
  t->Branch("m2kink", &fMu2Data.kink, "m2kink/F");
  t->Branch("m2dpt", &fMu2Data.dpt, "m2dpt/F");
  t->Branch("m2dptrel", &fMu2Data.dptrel, "m2dptrel/F");
  t->Branch("m2deta", &fMu2Data.deta, "m2deta/F");
  t->Branch("m2dphi", &fMu2Data.dphi, "m2dphi/F");
  t->Branch("m2dr", &fMu2Data.dr, "m2dr/F");

}

// ----------------------------------------------------------------------
void candAna::readCuts(string fileName, int dump) {

  // -- define default values for some cuts
  NOPRESELECTION = 0;
  IGNORETRIGGER  = 0;
  DSNAME = "Charmonium";

  // -- set up cut sequence for analysis
  basicCuts();
  moreBasicCuts();
  candidateCuts();
  moreCandidateCuts();

  fCutFile = fileName;

  if (dump) cout << "==> candAna: Reading " << fCutFile << " for cut settings" << endl;
  vector<string> cutLines;
  readFile(fCutFile, cutLines);

  char CutName[100];
  float CutValue;

  char  buffer[1000]; //  XmlName[1000];
  fHistDir->cd();
  if (dump) cout << "gDirectory: "; fHistDir->pwd();
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
    if ((fNchan + 1) == lineItems.size()) {
      for (unsigned int j = 1; j < lineItems.size(); ++j) {
	if (cutname == "metaMin") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->metaMin = cutvalue; ok = 1;
	}

	if (cutname == "metaMax") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->metaMax = cutvalue; ok = 1;
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

	if (cutname == "pv2lip") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pv2lip = cutvalue; ok = 1;
	}

	if (cutname == "pv2lips") {
	  cutvalue = atof(lineItems[j].c_str());
	  fCuts[j-1]->pv2lips = cutvalue; ok = 1;
	}

	if (cutname == "l1seeds") {
	  vector<string> vl1seeds = split(lineItems[j], ',');
	  for (unsigned int is = 0; is < vl1seeds.size(); ++is) {
	    fCuts[j-1]->l1seeds.push_back(atoi(vl1seeds[is].c_str()));
	  }
	}

      }

    }

    // -- now back to the original cut reading
    sprintf(buffer, "%s", cutLines[i].c_str());

    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue);
      if (dump) cout << "TYPE:           " << TYPE << endl;
      if (1313 == TYPE) cstring = "#mu^{+}#mu^{-}";
      if (301313 == TYPE) cstring = "#mu^{+}#mu^{-}";
      if (200521 == TYPE) cstring = "#mu^{+}#mu^{-}K^{+}";
      if (300521 == TYPE) cstring = "#mu^{+}#mu^{-}K^{+}";
      ibin = 1;
      hcuts->SetBinContent(ibin, TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
    }

    if (!strcmp(CutName, "SELMODE")) {
      SELMODE = int(CutValue);
      if (dump) cout << "SELMODE:           " << SELMODE << endl;
      ibin = 2;
      hcuts->SetBinContent(ibin, SELMODE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Selection Mode :: %i", CutName, SELMODE));
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

    if (!strcmp(CutName, "TRUTHCAND")) {
      TRUTHCAND = int(CutValue);
      if (dump) cout << "TRUTHCAND:           " << TRUTHCAND << endl;
      ibin = 4;
      hcuts->SetBinContent(ibin, TYPE);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: Candidate type", CutName));
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

    if (!strcmp(CutName, "PVAVEW8")) {
      PVAVEW8 = CutValue;
      if (dump) cout << "PVAVEW8:           " << PVAVEW8 << endl;
      ibin = 28;
      hcuts->SetBinContent(ibin, PVAVEW8);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: <w^{PV}_{trk}> :: %4.3f", CutName, PVAVEW8));
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
    if (!strcmp(CutName, "TRACKQUALITY")) {
      TRACKQUALITY = static_cast<int>(CutValue);
      if (dump) cout << "TRACKQUALITY:           " << TRACKQUALITY << " " << endl;
      ibin = 100;
      hcuts->SetBinContent(ibin, TRACKQUALITY);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: track quality :: %d", CutName, TRACKQUALITY));
    }

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

    // -- Muons
    if (!strcmp(CutName, "MUIDMASK")) {
      MUIDMASK = int(CutValue);
      if (dump) cout << "MUIDMASK:           " << MUIDMASK << endl;
      ibin = 200;
      hcuts->SetBinContent(ibin, MUIDMASK);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MuIDMask :: %d", CutName, MUIDMASK));
    }

    if (!strcmp(CutName, "MUIDRESULT")) {
      // MUIDRESULT == 0: compare result of & with ">=0"
      // MUIDRESULT != 0: compare result of & with "==MUIDRESULT"
      MUIDRESULT = int(CutValue);
      if (dump) cout << "MUIDRESULT:           " << MUIDRESULT << endl;
      ibin = 201;
      hcuts->SetBinContent(ibin, MUIDRESULT);
      hcuts->GetXaxis()->SetBinLabel(ibin, Form("%s :: MuIDResult :: %d", CutName, MUIDRESULT));
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

  if (dump) {
    for (unsigned int i = 0; i < fNchan; ++i) {
      fCuts[i]->dump();
    }
    cout << "------------------------------------" << endl;
  }
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
bool candAna::tightMuon(TAnaTrack *pT, bool hadronsPass) {

  const int verbose(0);

  if (verbose) cout << fYear << " --------- pT = " << pT->fPlab.Perp() << endl;

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    //cout << "NOTRIGGER requested... " << endl;
    return true;
  }

  //             654 3210
  // 80 = 0x50 = 0101 0000
  // global muon
  //  bool muflag = ((pT->fMuID & 2) == 2);
  // GMPT&&TMA:
  //  bool muflag = ((pT->fMuID & 80) == 80);
  // GMPT: 0100 0000 = 0x40
  bool muflag = ((pT->fMuID & 0x40) == 0x40);
  if (verbose) cout << "muflag: " << hex << pT->fMuID << dec << " -> " << muflag << endl;

  bool mucuts(true);
  if (verbose) cout << "mu index: " << pT->fMuIndex << " track index: " << pT->fIndex << endl;
  if (pT->fMuIndex > -1) {
    TAnaMuon *pM = fpEvt->getMuon(pT->fMuIndex);
    if (pM->fGtrkNormChi2 > 10.) mucuts = false;
    if (pM->fNvalidMuonHits < 1) mucuts = false;
    if (pM->fNmatchedStations < 2) mucuts = false;
    if (verbose) cout << "matched muon stations: " << pM->fNmatchedStations << " -> " << mucuts << endl;
  } else {
    mucuts = false;
  }

  bool trackcuts(true);

  if (fpReader->numberOfPixLayers(pT) < 1) trackcuts = false;
  if (verbose)  cout << "pixel layers: " << fpReader->numberOfPixLayers(pT) << " -> " << trackcuts << endl;

  //int trkHitsOld = fpReader->numberOfTrackerLayers(pT); UNUSED
  int trkHits    = pT->fLayersWithHits;

  if (fYear == 2011) {
    if (trkHits < 9) trackcuts = false;
    if (verbose)  cout << "valid hits: " << pT->fValidHits << " trackHist: " << trkHits << " -> " << trackcuts << endl;
  } else if (fYear == 2012) {
    if (trkHits < 6) trackcuts = false;
    if (verbose)  cout << "number of tracker layers: " << trkHits << " -> " << trackcuts << endl;
  } else {
    if (pT->fValidHits < 11) trackcuts = false;
    if (verbose)  cout << "valid hits: " << pT->fValidHits << " -> " << trackcuts << endl;
  }


  if (muflag && mucuts && trackcuts) {
    if (verbose) cout << " +++ passed "<<endl;
    return true;
  } else {
    //cout<<" failed "<<endl;
    return false;
  }



}


// ----------------------------------------------------------------------
bool candAna::tightMuon(TSimpleTrack *pT, bool hadronsPass) {
  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    return true;
  }
  if (0 == pT->getMuonID()) {
    return false;
  }

  TAnaMuon *pM(0);
  int idx = pT->getIndex();
  for (int i = 0; i < fpEvt->nMuons(); ++i) {
    pM = fpEvt->getMuon(i);
    if (idx == pM->fIndex) {
      return tightMuon(pM, hadronsPass);
    }
  }
  return false;
}




// ----------------------------------------------------------------------
bool candAna::mvaMuon(TAnaMuon *pt, double &result, bool hadronsPass) {

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    return true;
  }

  result = pt->fBarrelBDTresponse;
  return (result > 0);

}




// ----------------------------------------------------------------------
bool candAna::mvaMuon(TSimpleTrack *pt, double &result, bool hadronsPass) {

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    result = 99.;
    return true;
  }

  if (0 == pt->getMuonID()) {
    result = -3.;
    return false;
  }

  TAnaMuon *pM(0);
  int idx = pt->getIndex();
  for (int i = 0; i < fpEvt->nMuons(); ++i) {
    pM = fpEvt->getMuon(i);
    if (idx == pM->fIndex) {
      return mvaMuon(pM, result, hadronsPass);
    }
  }
  return false;

}


// ----------------------------------------------------------------------
bool candAna::mvaMuon(TAnaTrack *pt, double &result, bool hadronsPass) {

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    result = 99.;
    return true;
  }

  if (0 == pt->fMuID) {
    result = -3.;
    return false;
  }


  int idx = pt->fMuIndex;
  if (idx > -1 && idx < fpEvt->nMuons()) {
    TAnaMuon *pM = fpEvt->getMuon(idx);
    return mvaMuon(pM, result, hadronsPass);
  } else {
    cout << "muon index out of range!!!!!" << endl;
  }
  return false;

}


// ----------------------------------------------------------------------
bool candAna::mvaMuonLM(TAnaMuon *pt, double &result, bool hadronsPass) {

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    return true;
  }

  if (!tightMuon(pt)) {
    result = -2.;
    return false;
  }

  mrd.trkValidFract    = pt->fItrkValidFraction;
  mrd.glbNChi2         = pt->fGtrkNormChi2;
  mrd.pt               = pt->fPlab.Perp();
  mrd.eta              = pt->fPlab.Eta();
  mrd.segComp          = pt->fSegmentComp;
  mrd.chi2LocMom       = pt->fChi2LocalMomentum;
  mrd.chi2LocPos       = pt->fChi2LocalPosition;
  mrd.glbTrackProb     = pt->fGtrkProb;
  mrd.NTrkVHits        = static_cast<float>(pt->fNumberOfValidTrkHits);
  mrd.NTrkEHitsOut     = static_cast<float>(pt->fNumberOfLostTrkHits);

  mrd.dpt                  = pt->fNmatchedStations;
  mrd.intvalidpixelhits    = fpReader->numberOfPixLayers(pt); // FIXME, kind of correct
  mrd.inttrklayerswithhits = fpReader->numberOfTrackerLayers(pt);
  mrd.intnmatchedstations  = pt->fNmatchedStations;

  mrd.kink             = pt->fMuonChi2;

  mrd.dpt              = pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag();
  mrd.dptrel           = TMath::Abs(pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag())/pt->fInnerPlab.Mag();
  if (pt->fOuterPlab.Mag() > 3.) {
    mrd.deta             = pt->fInnerPlab.Eta() - pt->fOuterPlab.Eta();
    mrd.dphi             = pt->fInnerPlab.DeltaPhi(pt->fOuterPlab);
    mrd.dr               = pt->fInnerPlab.DeltaR(pt->fOuterPlab);
  } else {
    mrd.deta             = -99.;
    mrd.dphi             = -99.;
    mrd.dr               = -99.;
  }


  result = fMvaMuonID->EvaluateMVA("BDT");
  if (result > MUBDT) return true;
  return false;
}


// ----------------------------------------------------------------------
bool candAna::mvaMuonLM(TSimpleTrack *pt, double &result, bool hadronsPass) {

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    result = 99.;
    return true;
  }

  if (0 == pt->getMuonID()) {
    result = -3.;
    return false;
  }

  TAnaMuon *pM(0);
  int idx = pt->getIndex();
  for (int i = 0; i < fpEvt->nMuons(); ++i) {
    pM = fpEvt->getMuon(i);
    if (idx == pM->fIndex) {
      return mvaMuonLM(pM, result, hadronsPass);
    }
  }
  return false;

}


// ----------------------------------------------------------------------
bool candAna::mvaMuonLM(TAnaTrack *pt, double &result, bool hadronsPass) {

  if (hadronsPass && HLTRANGE.begin()->first == "NOTRIGGER") {
    result = 99.;
    return true;
  }

  if (0 == pt->fMuID) {
    result = -3.;
    return false;
  }


  int idx = pt->fMuIndex;
  if (idx > -1 && idx < fpEvt->nMuons()) {
    TAnaMuon *pM = fpEvt->getMuon(idx);
    return mvaMuonLM(pM, result, hadronsPass);
  } else {
    cout << "muon index out of range!!!!!" << endl;
  }
  return false;

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
pair<int, int> candAna::nCloseTracks(TAnaCand *pC, double dcaCut, double dcaCutS, double ptCut) {
  int cnt(0), cnts(0);
  int nsize = pC->fNstTracks.size();
  int pvIdx = pC->fPvIdx;
  //int pvIdx2= nearestPV(pvIdx, 0.1);  UNUSED d.k.
  //if (TMath::Abs(fCandPv2LipS) > 2) pv2Idx = -1;

  TSimpleTrack *pT;
  double pt(0.);
  if (nsize > 0) {
    for (int i = 0; i<nsize; ++i) {
      int trkId = pC->fNstTracks[i].first;
      double doca = pC->fNstTracks[i].second.first;

      if (doca > dcaCut) continue; // check the doca cut

      pT = fpEvt->getSimpleTrack(trkId);
      // -- check that any track associated with a definitive vertex is from the same or the closest (compatible) other PV
      if ((pT->getPvIndex() > -1) && (pT->getPvIndex() != pvIdx)) continue;

      pt = pT->getP().Perp();
      if (pt < ptCut) continue;

      ++cnt;
    }
  }

  if (nsize > 0) {
    for (int i = 0; i<nsize; ++i) {
      int trkId = pC->fNstTracks[i].first;
      double docas = pC->fNstTracks[i].second.first/pC->fNstTracks[i].second.second;

      if (docas > dcaCutS) continue; // check the doca cut

      pT = fpEvt->getSimpleTrack(trkId);
      // -- check that any track associated with a definitive vertex is from the same or the closest (compatible) other PV
      if ((pT->getPvIndex() > -1) && (pT->getPvIndex() != pvIdx)) continue;

      pt = pT->getP().Perp();
      if (pt < ptCut) continue;

      ++cnts;
    }
  }

  return make_pair(cnt, cnts);
}


// ----------------------------------------------------------------------
TAnaCand* candAna::osCand(TAnaCand *pC) {
  TAnaCand *a = new TAnaCand();
  a->fPvIdx = pC->fPvIdx;
  a->fType = pC->fType;
  a->fSig1 = a->fSig2 = -1;
  a->fPlab = TVector3(-pC->fPlab.X(), -pC->fPlab.Y(), -pC->fPlab.Z()); //???
  return a;
}


// ----------------------------------------------------------------------
double candAna::osIsolation(TAnaCand *pC, double r, double ptmin) {
  double iso(0.);
  TSimpleTrack *pT(0);
  int overlap(0), verbose(0);

  // -- look at all tracks that are associated to the same vertex
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i);
    if (verbose) {
      cout << "   track " << i
     	   << " with pT = " << pT->getP().Perp()
     	   << " eta = " << pT->getP().Eta()
     	   << " pointing at PV " << pT->getPvIndex();
    }

    // -- check against overlap to primary candidate
    overlap = 0;
    for (int j = pC->fSig1; j <= pC->fSig2; ++j) {
      if (i == fpEvt->getSigTrack(j)->fIndex) {
	overlap = 1;
	break;
      }
    }
    if (1 == overlap) continue;


// ???
//     if ((pT->fPvIdx > -1) && (pT->fPvIdx != pvIdx)) {
//       if (verbose) cout << " doca track " << trkId << " skipped because it is from a different PV " << pT->fPvIdx <<endl;
//       continue;
//     }


    if (pT->getPvIndex() != pC->fPvIdx) {
      if (verbose) cout << " skipped because of PV index mismatch" << endl; 	     //FIXME
      continue;
    }

    double pt = pT->getP().Perp();
    if (pt < ptmin) {
      if (verbose) cout << " skipped because of pt = " << pt << endl;
      continue;
    }

    if (pT->getP().DeltaR(pC->fPlab) > r) {
      iso += pt;
      if (verbose) cout << endl;
    }
    else {
      if (verbose) cout << " skipped because of deltaR = " << pT->getP().DeltaR(pC->fPlab) << endl;
    }
  }

  return iso;

}

// ----------------------------------------------------------------------
int candAna::osMuon(TAnaCand *pC, double r) {

  double mpt(-1.);
  int idx (-1);
  TSimpleTrack *pT(0);
  int overlap(0), verbose(0);

  // -- look at all tracks that are associated to the same vertex
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    pT = fpEvt->getSimpleTrack(i);
    if (!tightMuon(pT)) continue;

    // -- check against overlap to primary candidate
    overlap = 0;
    for (int j = pC->fSig1; j <= pC->fSig2; ++j) {
      if (i == fpEvt->getSigTrack(j)->fIndex) {
	overlap = 1;
	break;
      }
    }
    if (1 == overlap) continue;

    if (verbose) {
      cout << "   track " << i
     	   << " with pT = " << pT->getP().Perp()
     	   << " eta = " << pT->getP().Eta()
     	   << " pointing at PV " << pT->getPvIndex();
    }

    // ???
    //     if (pT->fPvIdx != pvIdx) {
    //       if (verbose) cout << " skipped because of PV index mismatch" << endl; 	     //FIXME
    //       continue;
    //     }

    if ((pT->getPvIndex() > -1) && (pT->getPvIndex() != pC->fPvIdx)) {
      if (verbose) cout << " track " << i << " skipped because it is from a different PV " << pT->getPvIndex() <<endl;
      continue;
    }

    if (pT->getP().DeltaR(pC->fPlab) > r) {
      if (pT->getP().Perp() > mpt) {
	mpt = pT->getP().Perp();
	idx = i;
      }
      if (verbose) cout << endl;
    }
    else {
      if (verbose) cout << " skipped because of deltaR = " << pT->getP().DeltaR(pC->fPlab) << endl;
    }
  }

  return idx;


}




// ----------------------------------------------------------------------
double candAna::isoClassicWithDOCA(TAnaCand *pC, double docaCut, double r, double ptmin) {
  const double ptCut(ptmin), coneSize(r);
  const bool verbose(false);

  double iso(-1.), pt(0.), sumPt(0.), candPt(0.), candPtScalar(0.);
  TSimpleTrack *ps;
  vector<int> cIdx, pIdx;
  int pvIdx = pC->fPvIdx;

  fCandI0trk = 0;
  fCandI1trk = 0;
  fCandI2trk = 0;

  getSigTracks(cIdx, pC);
  for (unsigned int i = 0; i < cIdx.size(); ++i) {
    //ps = fpEvt->getSimpleTrack(i);
    ps = fpEvt->getSimpleTrack(cIdx[i]);

    if (verbose) cout << " track idx = " << ps->getIndex() << " with ID = " << fpEvt->getSimpleTrackMCID(ps->getIndex()) << endl;
    candPtScalar += ps->getP().Perp();
    if (verbose) {
      int tIdx = fpEvt->getSimpleTrack(ps->getIndex())->getPvIndex();
      if (pvIdx != tIdx) {
    	cout << "Signal track pointing to PV " << tIdx << " instead of " << pvIdx << endl;
      }
    }
  }

  candPt = pC->fPlab.Perp();

  // -- look at all tracks that are associated to the same vertex
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    ps = fpEvt->getSimpleTrack(i);
    if (verbose) {
      cout << "   track " << i
     	   << " with pT = " << ps->getP().Perp()
     	   << " eta = " << ps->getP().Eta()
     	   << " pointing at PV " << ps->getPvIndex();
    }


    if (ps->getPvIndex() != pvIdx) {
      if (verbose) cout << " skipped because of PV index mismatch" << endl; 	     //FIXME
      continue;
    }


    pt = ps->getP().Perp();
    if (pt < ptCut) {
      if (verbose) cout << " skipped because of pt = " << pt << endl;
      continue;
    }
    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  {
      if (verbose) cout << " skipped because it is a sig track " << endl;
      continue;
    }
    if (ps->getP().DeltaR(pC->fPlab) < coneSize) {
      pIdx.push_back(i);
      ++fCandI0trk;
      sumPt += pt;
      if (verbose) cout << endl;
    }
    else {
      if (verbose) cout << " skipped because of deltaR = " << ps->getP().DeltaR(pC->fPlab) << endl;
    }
  }

  // -- Now consider the DOCA tracks
  int nsize = pC->fNstTracks.size();
  if (nsize>0) {
    for(int i = 0; i<nsize; ++i) {
      int trkId = pC->fNstTracks[i].first;
      double doca = pC->fNstTracks[i].second.first;
      // double docaE = pC->fNstTracks[i].second.second;

      if (doca > docaCut) continue; // check the doca cut

      ps = fpEvt->getSimpleTrack(trkId);


      if ((ps->getPvIndex() > -1) && (ps->getPvIndex() != pvIdx)) {
	if (verbose) cout << " doca track " << trkId << " skipped because it is from a different PV " << ps->getPvIndex() <<endl;
	continue;
      }

      pt = ps->getP().Perp();
      if (pt < ptCut) {
	if (verbose) cout << " doca track " << trkId << " skipped because of pt = " << pt << endl;
	continue;
      }

      if (ps->getP().DeltaR(pC->fPlab) > coneSize) {
	if (verbose) cout << " doca track " << trkId << " skipped because of deltaR = " << ps->getP().DeltaR(pC->fPlab) << endl;
	continue;
      }

      // -- Skip tracks already included above
      if (pIdx.end() != find(pIdx.begin(), pIdx.end(), trkId))  continue;
      if (cIdx.end() != find(cIdx.begin(), cIdx.end(), trkId))  continue;

      ++fCandI1trk;
      sumPt += pt;
      if (verbose) cout << " doca track " << trkId << " included "<<doca<<" "<<pt<<endl;

    } // for loop over tracks
  } // end if

  fCandI2trk = fCandI0trk + fCandI1trk;

  iso = candPt/(candPt + sumPt);

  //   if (verbose) cout << "--> iso = " << candPt << " .. " << sumPt << " = " << iso << endl;
  //   if (verbose) cout << "--> iso = " << pC->fPlab.Perp() << " .. " << sumPt << " = " << pC->fPlab.Perp()/(pC->fPlab.Perp() + sumPt) << endl;

  return iso;
}


// ----------------------------------------------------------------------
void candAna::xpDistMuons() {

  fMu1XpDist = -99.;
  fMu2XpDist = -99.;

  if (fpMuon1->fMuIndex < 0) return;
  if (fpMuon2->fMuIndex < 0) return;

  TAnaMuon *pm[2];
  pm[0] = fpEvt->getMuon(fpMuon1->fMuIndex);
  pm[1] = fpEvt->getMuon(fpMuon2->fMuIndex);

  //int set0(-99), set1(-99);
  for (int m = 0; m < 2; ++m) {
    for (int i = 0; i < 10; ++i) {

      int trackIdx = pm[m]->fXpTracks[i].idx;
      if (trackIdx < 0 || trackIdx > fpEvt->nSimpleTracks()) continue;

      // -- ignore this track if it is the other muon of this decay
      if (trackIdx == pm[1-m]->fIndex) continue;

      // -- check that the track's PV is the current one of it is set
      int pvIdx = fpEvt->getSimpleTrack(trackIdx)->getPvIndex();
      if (pvIdx > -1) {
	if (pvIdx != fpCand->fPvIdx) continue;
      }

      // -- now set the mu i xp dist and break for this one
      if (0 == m) {
	fMu1XpDist = pm[m]->fXpTracks[i].dist;
	//set0 = i;
	break;
      }
      if (1 == m) {
	fMu2XpDist = pm[m]->fXpTracks[i].dist;
	//set1 = i;
	break;
      }
    }
  }

  //  cout << "used tracks " << set0 << " and " << set1 << endl;

}



// ----------------------------------------------------------------------
double candAna::isoMuon(TAnaCand *pCand, TAnaMuon *pMuon) {

  TSimpleTrack *sTrack;
  std::vector<near_track_t> nearTracks;
  std::map<int,float>::const_iterator it;
  std::map<int,int> cand_tracks;
  near_track_t nt;
  size_t k;

  // get the candidate structure
  findAllTrackIndices(pCand, &cand_tracks);

  TVector3 plabMu = pMuon->fPlab;
  for (it = pMuon->fNstTracks.begin(); it != pMuon->fNstTracks.end(); ++it) {

    // no tracks from candidate...
    if (cand_tracks.count(it->first) > 0)
      continue;

    sTrack = fpEvt->getSimpleTrack(it->first);

    // no tracks from foreign primary vertex
    if (sTrack->getPvIndex() >= 0 && sTrack->getPvIndex() != pCand->fPvIdx)
      continue;

    nt.ix = it->first;
    nt.doca = it->second;
    nt.p = sTrack->getP().Mag();
    nt.pt = sTrack->getP().Perp();
    TVector3 recTrack = sTrack->getP();
    nt.pt_rel = (recTrack - (recTrack * plabMu) * plabMu).Mag() / plabMu.Mag2();
    nt.deltaR = plabMu.DeltaR(recTrack);
    nearTracks.push_back(nt);
  }

  // compute isolation variable
  double result = 0;
  for (k = 0; k < nearTracks.size(); k++) {
    if (nearTracks[k].deltaR < 0.5 && nearTracks[k].pt > 0.5 && nearTracks[k].doca < 0.1) // 1 mm
      result += nearTracks[k].p;
  }
  result = plabMu.Mag()/(plabMu.Mag() + result);

  return result;
}


// ----------------------------------------------------------------------
void candAna::findAllTrackIndices(TAnaCand* pCand, map<int,int> *indices) {
  int j;

  // iterate through all own tracks. has to be done first, so the duplicate signal tracks
  // won't be added in the daughter anymore
  for (j = pCand->fSig1; j <= pCand->fSig2 && j>=0; j++)
    indices->insert(make_pair(fpEvt->getSigTrack(j)->fIndex,j));

  for (j = pCand->fDau1; j <= pCand->fDau2 && j>=0; j++)
    findAllTrackIndices(fpEvt->getCand(j),indices);
}


// ----------------------------------------------------------------------
void candAna::muScaleCorrectedMasses() {
  fCandM3 = fCandM4 = -99.;
  TLorentzVector myNegMuon, myPosMuon;
  if (fpMuon1->fQ < 0) {
    myNegMuon.SetXYZM(fpMuon1->fPlab.X(), fpMuon1->fPlab.Y(), fpMuon1->fPlab.Z(), MMUON);
    myPosMuon.SetXYZM(fpMuon2->fPlab.X(), fpMuon2->fPlab.Y(), fpMuon2->fPlab.Z(), MMUON);
  } else {
    myPosMuon.SetXYZM(fpMuon1->fPlab.X(), fpMuon1->fPlab.Y(), fpMuon1->fPlab.Z(), MMUON);
    myNegMuon.SetXYZM(fpMuon2->fPlab.X(), fpMuon2->fPlab.Y(), fpMuon2->fPlab.Z(), MMUON);
  }

  double mass0 = (myNegMuon + myPosMuon).Mag();
  double mass1 = (myNegMuon + myPosMuon).Mag();
  //  cout << "candMass: " << fCandM << " TLV mass = " << mass0 << " msc mass = " << mass1 << endl;
  fCandM3 = mass1;
  fCandM4 = (mass1/mass0)*fCandM;

}


// ----------------------------------------------------------------------
double candAna::constrainedMass() {
  vector<int> rectracks;
  int bla;
  for (int it = fpCand->fSig1; it <= fpCand->fSig2; ++it) {
    bla = fpEvt->getSigTrack(it)->fIndex;
    rectracks.push_back(bla);
  }

  TAnaCand *pC(0), *pDau(0), *pMcCand(0);
  int mctype = TYPE + 100000;
  unsigned int nmatch(0);
  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pC = fpEvt->getCand(iC);
    if (pC->fType != mctype) continue;
    nmatch = 0;
    for (int it = pC->fSig1; it <= pC->fSig2; ++it) {
      bla = fpEvt->getSigTrack(it)->fIndex;
      if (rectracks.end() != find(rectracks.begin(), rectracks.end(), bla)) {
	++nmatch;
      }
    }
    if (pC->fDau1 > 0) {
      pDau = fpEvt->getCand(pC->fDau1);
      for (int it = pDau->fSig1; it <= pDau->fSig2; ++it) {
	bla = fpEvt->getSigTrack(it)->fIndex;
	if (rectracks.end() != find(rectracks.begin(), rectracks.end(), bla)) {
	  ++nmatch;
	}
      }
    }
    if (nmatch == rectracks.size()) {
      pMcCand = pC;
      break;
    }
  }
  if (pMcCand) {
    return pMcCand->fMass;
  } else {
    return -2.;
  }
}


// ----------------------------------------------------------------------
void candAna::runRange() {
  fRunRange = 6;
  if ((fRun >= 160329) && (fRun <= 163261)) {
    fRunRange = 0; // the beginning
  }
  if ((fRun >= 163269) && (fRun <= 163869)) {
    fRunRange = 1; // displaced J/psi
  }
  if ((fRun >= 165088) && (fRun <= 167913)) {
    fRunRange = 2; // HLTDimuon7
  }
  if ((fRun >= 170249) && (fRun <= 173198)) {
    fRunRange = 3; // 2e33
  }
  if ((fRun >= 173236) && (fRun <= 178380)) {
    fRunRange = 4; // 3e33 WITH the ETA cut!
  }
  if ((fRun >= 178420) && (fRun <= 999999)) {
    fRunRange = 5; // 5e33
  }

  if (fIsMC) {
    if (string::npos != fHLTPath.find("HLT_Dimuon7_Jpsi_Displaced_v1")) {
      fRunRange = 2;
    }
    if (string::npos != fHLTPath.find("HLT_DoubleMu3p5_Jpsi_Displaced_v2_Bs")) {
      fRunRange = 3;
    }
    if (string::npos != fHLTPath.find("HLT_DoubleMu4_Jpsi_Displaced_v1")) {
      fRunRange = 4;
    }
  }

}




// ----------------------------------------------------------------------
int candAna::nearestPV(int pvIdx, double maxDist) {

  if (pvIdx==-1) return -1;  // add protection d.k. 9/6/12

  TAnaVertex *v0 = fpEvt->getPV(pvIdx);

  if (0 == v0) return -1;  // add protection

  double zV0 = v0->fPoint.Z();

  int idx(-1);
  double z(0.), zabs(0.), zmin(99.);
  for (int i = 0; i < fpEvt->nPV(); ++i) {
    if (i == pvIdx) continue;
    z = fpEvt->getPV(i)->fPoint.Z();
    zabs = TMath::Abs(zV0 - z);
    if (zabs < zmin) {
      idx = i;
      zmin = zabs;
    }
  }

  if (zmin < maxDist) {
    //     cout << "pvIdx = " << pvIdx << " at = " << zV0
    // 	 << ", nearest other PV with idx = " << idx << " at z = " << fpEvt->getPV(idx)->fPoint.Z()
    // 	 << " and delta(z) = " << zmin << endl;
    return idx;
  } else {
    return -1;
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
void candAna::calcBDT() {
  fBDT = -99.;
  //??  if (5 == mode && 5.2 < mass && mass < 5.45 && fb.iso < 0.7) continue;
  if (fChan < 0) return;
  if (0 == fReaderEvents0.size()) {
    cout << "no BDT defined" << endl;
    return;
  }

  if (!preselection(fRTD, fChan)) return;

  //   if (fCandPt > 100) return;
  //   if (fCandPt < 6) return;
  //   if (fMu1Pt < 4) return;
  //   if (fMu2Pt < 4) return;
  //   if (fCandFL3d > 1.5) return;
  //   if (fCandFL3d < 0.) return;
  //   if (fCandM > 5.9) return;
  //   if (fCandM < 4.9) return;

  //   if (!fb.hlt) return;
  //   if (!fb.gmuid) return;

  frd.pt = fCandPt;
  frd.eta = fCandEta;
  frd.m1eta = fMu1Eta;
  frd.m2eta = fMu2Eta;
  frd.m1pt = fMu1Pt;
  frd.m2pt = fMu2Pt;
  frd.fls3d = fCandFLS3d;
  frd.alpha = fCandA;
  frd.maxdoca = fCandDoca;
  frd.pvip = fCandPvIp;
  frd.pvips = fCandPvIpS;
  frd.iso = fCandIso;
  frd.docatrk = fCandDocaTrk;
  frd.chi2dof = fCandChi2/fCandDof;
  frd.closetrk = fCandCloseTrk;

  frd.m  = fCandM;
  //  cout << "Evt = " << fEvt << " %3 = " << fEvt%3 << " chan = " << fChan << " " << " etas = " << fMu1Eta << " " << fMu2Eta;

  fBDT = 0.;
  return;
  if (0 == fEvt%3) {
    fBDT   = fReaderEvents0[fChan]->EvaluateMVA("BDT");
  } else if (1 == fEvt%3) {
    fBDT   = fReaderEvents1[fChan]->EvaluateMVA("BDT");
  } else if (2 == fEvt%3) {
    fBDT   = fReaderEvents2[fChan]->EvaluateMVA("BDT");
  } else {
    cout << "all hell break loose" << endl;
  }
  //  cout << " bdt = " << fBDT << endl;
}


// // ----------------------------------------------------------------------
// int candAna::detChan(double m1eta, double m2eta) {
//   // -- simple two channel analysis: channel 0 if both muons in barrel, channel 1 else
//   if (TMath::Abs(m1eta) < 1.4 && TMath::Abs(m2eta) < 1.4) return 0;
//   if (TMath::Abs(m1eta) < 2.4 && TMath::Abs(m2eta) < 2.4) return 1;
//   return -1;
// }


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
	if (stype == "trkValidFract") {
	  cout << "  adding trkValidFract" << endl;
	  reader->AddVariable( "trkValidFract", &d.trkValidFract);
	  continue;
	}
	if (stype == "glbNChi2") {
	  reader->AddVariable("glbNChi2", &d.glbNChi2);
	  cout << "  adding glbNChi2" << endl;
	  continue;
	}
	if (stype == "eta") {
	  cout << "  adding eta" << endl;
	  reader->AddVariable("eta", &d.eta);
	  continue;
	}
	if (stype == "pt") {
	  cout << "  adding pt" << endl;
	  reader->AddVariable("pt", &d.pt);
	  continue;
	}
	if (stype == "segComp") {
	  cout << "  adding segComp" << endl;
	  reader->AddVariable("segComp", &d.segComp);
	  continue;
	}
	if (stype == "chi2LocMom") {
	  cout << "  adding chi2LocMom" << endl;
	  reader->AddVariable("chi2LocMom", &d.chi2LocMom);
	  continue;
	}
	if (stype == "chi2LocPos") {
	  cout << "  adding chi2LocPos" << endl;
	  reader->AddVariable("chi2LocPos", &d.chi2LocPos);
	  continue;
	}
	if (stype == "glbTrackProb") {
	  cout << "  adding glbTrackProb" << endl;
	  reader->AddVariable("glbTrackProb", &d.glbTrackProb);
	  continue;
	}
	if (stype == "NTrkVHits") {
	  cout << "  adding NTrkVHits" << endl;
	  reader->AddVariable("NTrkVHits", &d.NTrkVHits);
	  continue;
	}
	if (stype == "NTrkEHitsOut") {
	  cout << "  adding NTrkEHitsOut" << endl;
	  reader->AddVariable("NTrkEHitsOut", &d.NTrkEHitsOut);
	  continue;
	}
	if (stype == "NTrkEHitsOut") {
	  cout << "  adding NTrkEHitsOut" << endl;
	  reader->AddVariable("NTrkEHitsOut", &d.NTrkEHitsOut);
	  continue;
	}

	// -- new convention for UL's BDT
	if (stype == "intnmatchedstations") {
	  cout << "  adding intnmatchedstations" << endl;
	  reader->AddVariable("intnmatchedstations", &d.intnmatchedstations);
	  continue;
	}
	if (stype == "intvalidpixelhits") {
	  cout << "  adding intvalidpixelhits" << endl;
	  reader->AddVariable("intvalidpixelhits", &d.intvalidpixelhits);
	  continue;
	}
	if (stype == "inttrklayerswithhits") {
	  cout << "  adding inttrklayerswithhits" << endl;
	  reader->AddVariable("inttrklayerswithhits", &d.inttrklayerswithhits);
	  continue;
	}
	if (stype == "gchi2") {
	  cout << "  adding gchi2" << endl;
	  reader->AddVariable("gchi2", &d.glbNChi2);
	  continue;
	}
	if (stype == "itrkvalidfraction") {
	  cout << "  adding itrkvalidfraction" << endl;
	  reader->AddVariable("itrkvalidfraction", &d.trkValidFract);
	  continue;
	}
	if (stype == "segcomp") {
	  cout << "  adding segcomp" << endl;
	  reader->AddVariable("segcomp", &d.segComp);
	  continue;
	}
	if (stype == "chi2lmom") {
	  cout << "  adding chi2lmom" << endl;
	  reader->AddVariable("chi2lmom", &d.chi2LocMom);
	  continue;
	}
	if (stype == "chi2lpos") {
	  cout << "  adding chi2lpos" << endl;
	  reader->AddVariable("chi2lpos", &d.chi2LocPos);
	  continue;
	}
	if (stype == "gtrkprob") {
	  cout << "  adding gtrkprob" << endl;
	  reader->AddVariable("gtrkprob", &d.glbTrackProb);
	  continue;
	}
	if (stype == "ntrkvhits") {
	  cout << "  adding ntrkvhits" << endl;
	  reader->AddVariable("ntrkvhits", &d.NTrkVHits);
	  continue;
	}
	if (stype == "inttrklayerswithhits") {
	  cout << "  adding inttrklayerswithhits" << endl;
	  reader->AddVariable("inttrklayerswithhits", &d.inttrklayerswithhits);
	  continue;
	}
	if (stype == "ntrkehitsout") {
	  cout << "  adding ntrkehitsout" << endl;
	  reader->AddVariable("ntrkehitsout", &d.NTrkEHitsOut);
	  continue;
	}
	if (stype == "kink") {
	  cout << "  adding kink" << endl;
	  reader->AddVariable("kink", &d.kink);
	  continue;
	}
	if (stype == "dpt") {
	  cout << "  adding dpt" << endl;
	  reader->AddVariable("dpt", &d.dpt);
	  continue;
	}
	if (stype == "deta") {
	  cout << "  adding deta" << endl;
	  reader->AddVariable("deta", &d.deta);
	  continue;
	}
	if (stype == "dphi") {
	  cout << "  adding dphi" << endl;
	  reader->AddVariable("dphi", &d.dphi);
	  continue;
	}
	if (stype == "dr") {
	  cout << "  adding dr" << endl;
	  reader->AddVariable("dr", &d.dr);
	  continue;
	}
	if (stype == "dptrel") {
	  cout << "  adding dptrel" << endl;
	  reader->AddVariable("dptrel", &d.dptrel);
	  continue;
	}


      }
      break;
    }
  }

  reader->BookMVA("BDT", TString(weightFile.c_str()));
  return reader;
}


// ----------------------------------------------------------------------
TMVA::Reader* candAna::setupReader(string xmlFile, readerData &rd) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

  TString dir    = "weights/";
  TString methodNameprefix = "BDT";
  //  TString methodName = TString(fBdt) + TString(" method");
  //  TString weightfile = dir + fBdt + "_" + methodNameprefix + TString(".weights.xml");
  TString weightfile = xmlFile;

  // -- read in variables from weight file
  vector<string> allLines;
  char  buffer[2000];
  cout << "setupReader, open file " << weightfile << endl;
  ifstream is(weightfile);
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
	//	cout << "ivar " << j-i << " variable string: ->" << stype << "<-" << endl;
	if (stype == "m1pt") {
	  cout << "  adding m1pt" << endl;
	  reader->AddVariable( "m1pt", &rd.m1pt);
	}
	if (stype == "m2pt") {
	  cout << "  adding m2pt" << endl;
	  reader->AddVariable( "m2pt", &rd.m2pt);
	}
	if (stype == "m1eta") {
	  cout << "  adding m1eta" << endl;
	  reader->AddVariable( "m1eta", &rd.m1eta);
	}
	if (stype == "m2eta") {
	  reader->AddVariable( "m2eta", &rd.m2eta);
	  cout << "  adding m2eta" << endl;
	}
	if (stype == "pt") {
	  cout << "  adding pt" << endl;
	  reader->AddVariable( "pt", &rd.pt);
	}
	if (stype == "eta") {
	  cout << "  adding eta" << endl;
	  reader->AddVariable( "eta", &rd.eta);
	}
	if (stype == "fls3d") {
	  cout << "  adding fls3d" << endl;
	  reader->AddVariable( "fls3d", &rd.fls3d);
	}
	if (stype == "alpha") {
	  cout << "  adding alpha" << endl;
	  reader->AddVariable( "alpha", &rd.alpha);
	}
	if (stype == "maxdoca") {
	  cout << "  adding maxdoca" << endl;
	  reader->AddVariable( "maxdoca", &rd.maxdoca);
	}
	if (stype == "pvip") {
	  cout << "  adding pvip" << endl;
	  reader->AddVariable( "pvip", &rd.pvip);
	}
	if (stype == "pvips") {
	  cout << "  adding pvips" << endl;
	  reader->AddVariable( "pvips", &rd.pvips);
	}
	if (stype == "iso") {
	  cout << "  adding iso" << endl;
	  reader->AddVariable( "iso", &rd.iso);
	}
	if (stype == "docatrk") {
	  cout << "  adding docatrk" << endl;
	  reader->AddVariable( "docatrk", &rd.docatrk);
	}
	if (stype == "closetrk") {
	  cout << "  adding closetrk" << endl;
	  reader->AddVariable( "closetrk", &rd.closetrk);
	}
	if (stype == "closetrks1") {
	  cout << "  adding closetrks1" << endl;
	  reader->AddVariable( "closetrks1", &rd.closetrks1);
	}
	if (stype == "closetrks2") {
	  cout << "  adding closetrks2" << endl;
	  reader->AddVariable( "closetrks2", &rd.closetrks2);
	}
	if (stype == "closetrks3") {
	  cout << "  adding closetrks3" << endl;
	  reader->AddVariable( "closetrks3", &rd.closetrks3);
	}
	if (stype == "chi2dof") {
	  cout << "  adding chi2dof" << endl;
	  reader->AddVariable( "chi2dof", &rd.chi2dof);
	}
	if (stype == "m1iso") {
	  cout << "  adding m1iso" << endl;
	  reader->AddVariable( "m1iso", &rd.m1iso);
	}
	if (stype == "m2iso") {
	  cout << "  adding m2iso" << endl;
	  reader->AddVariable( "m2iso", &rd.m2iso);
	}
	if (stype == "pvdchi2") {
	  cout << "  adding pvdchi2" << endl;
	  reader->AddVariable( "pvdchi2", &rd.pvdchi2);
	}
	if (stype == "othervtx") {
	  cout << "  adding othervtx" << endl;
	  reader->AddVariable( "othervtx", &rd.othervtx);
	}
	if (stype == "pv2lip") {
	  cout << "  adding pv2lip" << endl;
	  reader->AddVariable( "pv2lip", &rd.pv2lip);
	}
	if (stype == "pv2lips") {
	  cout << "  adding pv2lips" << endl;
	  reader->AddVariable( "pv2lips", &rd.pv2lips);
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
	cout << "ivar " << j-i << " spectator string: ->" << stype << "<-" << endl;
	if (stype == "m") {
	  cout << "  adding m as spectator" << endl;
	  reader->AddSpectator( "m", &rd.m);
	}
      }
      break;
    }
  }

  // --- Book the MVA methods
  reader->BookMVA("BDT", weightfile);
  return reader;
}


// ----------------------------------------------------------------------
void candAna::replaceAll(std::string &s, std::string a, std::string b) {

  TString ts(s.c_str());
  ts.ReplaceAll(a.c_str(), b.c_str());
  s = ts.Data();

}



// ----------------------------------------------------------------------
void candAna::fillRedTreeData() {
  // -- this only fills the variables that are needed for the preselection() function
  fRTD.hlt       = fGoodHLT;
  fRTD.gmuid     = fGoodMuonsID;

  fRTD.pt        = fCandPt;
  fRTD.eta       = fCandEta;
  fRTD.m         = fCandM;

  fRTD.m1pt      = fMu1Pt;
  fRTD.m2pt      = fMu2Pt;

  fRTD.m1eta     = fMu1Eta;
  fRTD.m2eta     = fMu2Eta;

  fRTD.m1q       = fMu1Q;
  fRTD.m2q       = fMu2Q;

  fRTD.pvip      = fCandPvIp;
  fRTD.pvips     = fCandPvIpS;

  fRTD.pvlip     = fCandPvLip;
  fRTD.pvlips    = fCandPvLipS;

  fRTD.closetrk  = fCandCloseTrk;
  fRTD.iso       = fCandIso;

  fRTD.flsxy     = fCandFLSxy;
  fRTD.fl3d      = fCandFL3d;
  fRTD.fls3d     = fCandFLS3d;

  fRTD.chi2dof   = fCandChi2Dof;
  fRTD.chi2      = fCandChi2;
  fRTD.dof       = fCandDof;

  fRTD.alpha     = fCandA;

  fRTD.docatrk   = fCandDocaTrk;
  fRTD.maxdoca   = fCandDoca;

}
//-------------------------------------------------------------------------------

// A trigger matcher based on deltaR (from Frank) + pt matching.
// check 2 muons, use only the selected hlt objects which correspond to triggers
// which passed and were on out trigger list.
// Only consider trig objects which match our trigger list.
// The main cuts are: deltaRthr for DR and deltaPtMatch for pt
// uses TTrgObjv2
bool candAna::doTriggerMatching(TAnaTrack *fp1, TAnaTrack *fp2) { // call the normal version with (true)
  int indx1=-1, indx2=-1;
  const double deltaRthr(0.02); // final cut, Frank had 0.5, change 0.020
  const double deltaPtMatch(0.15); // the pt matching cut
  const int verboseThr = 30;
  //const bool localPrint = false;
  bool localPrint = (fVerbose==-32) || (fVerbose > verboseThr);
  int mu1match(-1), mu2match(-1);
  string hlt1, hlt2;
  double deltaRmin1(100),deltaRmin2(100);
  double trigMatchDeltaPt1 = 99., trigMatchDeltaPt2 = 99.;
  bool match=false;
  TTrgObjv2 *pTO;
  TLorentzVector tlvMu1, tlvMu2;

  if (localPrint) {
    cout << "mu1: pt,eta,phi: " << fp1->fPlab.Perp() << " " << fp1->fPlab.Eta() << " " << fp1->fPlab.Phi()<< endl;
    cout << "mu2: pt,eta,phi: " << fp2->fPlab.Perp() << " " << fp2->fPlab.Eta() << " " << fp2->fPlab.Phi()<< endl;
  }

  tlvMu1.SetPtEtaPhiM(fp1->fPlab.Perp(),fp1->fPlab.Eta(),fp1->fPlab.Phi(),MMUON); // assume a muon
  tlvMu2.SetPtEtaPhiM(fp2->fPlab.Perp(),fp2->fPlab.Eta(),fp2->fPlab.Phi(),MMUON); // assume a muon

  map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    int hltIndex = pTO->fHltIndex;

    bool lastModule= false;
    bool activeModule = false;
    ix=hltObjMap.find(i);
    //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
    //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
    if (ix!=hltObjMap.end()) {
      activeModule=true; // signal the teh module belongs a to selected/active HLT
      lastModule = ( (ix->second & 0x80000000) != 0);
    }

    if (lastModule) { // this object was selected, use last module
      if (localPrint) cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
			 <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<" "<<activeModule<<endl;

      bool match1=false, match2=false;
      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
	int index = muonIndex[n];
	int id = muonID[n];
	TLorentzVector p = muonP[n];

	if ( abs(id) != 13 ) { // if not muon trigger skip
	  if (localPrint) cout<<" a none hlt-muon found in a trigger object, skip it, id= "
			     <<id<<" "<<pTO->fHltPath<<" "<<pTO->fLabel<<" "<<pTO->fType<<endl;
	  continue;  // skip checking non-muon objects
	}

	// check direction matching
	double deltaR1 = p.DeltaR(tlvMu1);
	double deltaR2 = p.DeltaR(tlvMu2);

	if (localPrint) {
	  cout<<" particle"<<n<<" index "<<index<<" id "<<id
	      <<" pt/eta/phi "<<p.Pt()<<"/"<<p.Eta()<<"/"<<p.Phi()<<" i/n "<<i<<"/"<<n
	      <<" dr "<<deltaR1 <<" "<<deltaR2<<endl;
	}

	// muon 1
	if (deltaR1<deltaRmin1) {
	  deltaRmin1=deltaR1;  // best match until now
	  if (fVerbose > verboseThr || localPrint) {cout << " mu1 selected "<< deltaR1 <<endl;}
	    // check now the pt matching
	  double trigMatchDeltaPt=999.;
	  if (fp1->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp1->fPlab.Mag())/fp1->fPlab.Mag();
	  if ( trigMatchDeltaPt < deltaPtMatch ) {  // check if it is good enough
	    if (deltaR1<deltaRthr) {
	      trigMatchDeltaPt1=trigMatchDeltaPt;
	      mu1match = n;
	      hlt1 = pTO->fLabel;
	      indx1=i;
	      match1=true;
	    } // if delta
	  } // if pt match
	} // if direction match

	// muon 2
	if (deltaR2<deltaRmin2) {
	  deltaRmin2=deltaR2;
	  if (localPrint) {cout << " mu2 selected "<< deltaR2 <<endl;}
	    // check now the pt matching
	  double trigMatchDeltaPt=999.;
	  if (fp2->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp2->fPlab.Mag())/fp2->fPlab.Mag();
	  if ( trigMatchDeltaPt < deltaPtMatch ) {
	    if (deltaR2<deltaRthr) {
	      trigMatchDeltaPt2=trigMatchDeltaPt;
	      mu2match = n;
	      hlt2 = pTO->fLabel;
	      indx2=i;
	      match2=true;
	    } // if delta
	  } // if pt match
	} // if direction match
      } // end for loop n

      // check that at least one module matched both
      match = match || (match1&&match2);

    } // end if valid module

  } // loop over all modules

  if (localPrint)
    cout << " best match "
	 <<indx1<<" "<< deltaRmin1 << " "<<mu1match<<" "<<hlt1<<" "<<trigMatchDeltaPt1<<" "
	 <<indx2<<" "<< deltaRmin2 << " "<<mu2match<<" "<<hlt2<<" "<<trigMatchDeltaPt2<<endl;

  ((TH1D*)fHistDir->Get("test8"))->Fill(trigMatchDeltaPt1);
  ((TH1D*)fHistDir->Get("test8"))->Fill(trigMatchDeltaPt2);
  ((TH1D*)fHistDir->Get("test2"))->Fill(deltaRmin1);
  ((TH1D*)fHistDir->Get("test2"))->Fill(deltaRmin2);

  if (mu1match>-1) {
    double tmp=fMu1TrigM;
    fMu1TrigM = deltaRmin1;
    if (tmp!=fMu1TrigM) cout<<"Warning:  two methods inconsistent-mu1 "<<tmp<<" "<<fMu1TrigM<<endl;
  }
  if (mu2match>-1) {
    double tmp=fMu2TrigM;
    fMu2TrigM = deltaRmin2;
    if (tmp!=fMu2TrigM) cout<<"Warning:  two methods inconsistent-mu2 "<<tmp<<" "<<fMu2TrigM<<endl;
  }

  bool HLTmatch = false;
  if (match && mu1match>-1 && mu2match>-1) {
    if (indx1!=indx2) { // should never happen since usually we only have one selected trigger
      cout<<"Warning:  best match for the two muons is to two different modules "<<indx1<<" "<<indx2<<endl;
      HLTmatch=true;
    } else if (mu1match==mu2match) { // matched to 2 same tracks, what to do? skip it?
      cout<<"Error:  two muons matched to same particle "<<indx1<<" "<<indx2<<" "<<mu1match<<endl;
    } else { // ok
      if (localPrint) cout<<" matching OK"<<indx1<<endl;
      HLTmatch=true;
    }
  }

  return HLTmatch;
}
//-------------------------------------------------------------------------------
// match track to trigger, return DR
// return the best, smalles DR
// anyTrig = true - match to any triggered which fired & created a muon object (might be from another DS)
//         = false - match only to the trigger selected from our list
// muonsOnly = true - uae only muon trigger particles to match
//           = false - use all trigger particles
double candAna::doTriggerMatchingR(TAnaTrack *fp1, bool anyTrig, bool muonsOnly, bool anyModule) {
  int indx1=-1;
  //const double deltaRthr(0.02); // final cut, Frank had 0.5, change 0.020
  const double deltaPtMatch(0.15); // the pt matching cut
  const int verboseThr = 30;
  //const bool localPrint = false;
  bool localPrint = (fVerbose==-32) || (fVerbose > verboseThr);
  int mu1match(-1);
  string hlt1;
  double deltaRmin1(100);
  double trigMatchDeltaPt1 = 99.;

  TTrgObjv2 *pTO;
  TLorentzVector tlvMu1;

  if (localPrint) {
    cout << " doTriggerMatchingR pt,eta,phi: " << fp1->fPlab.Perp() << " " << fp1->fPlab.Eta() << " " << fp1->fPlab.Phi()<< endl;
  }

  tlvMu1.SetPtEtaPhiM(fp1->fPlab.Perp(),fp1->fPlab.Eta(),fp1->fPlab.Phi(),MMUON); // assume a muon

  //cout<<" size "<<hltObjMap.size()<<endl;
  map<unsigned int, unsigned int, less<unsigned int> >::iterator  ix;
  for(int i=0; i!=fpEvt->nTrgObjv2(); i++) { // loop over all objects
    pTO = fpEvt->getTrgObjv2(i);
    int hltIndex = pTO->fHltIndex;

    bool lastModule= false;
    bool activeModule = false;
    ix=hltObjMap.find(i);
    //cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
    //	<<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber;
    if (ix!=hltObjMap.end()) {
      activeModule=true; // signal the teh module belongs a to selected/active HLT
      int num = ix->first;  // index of the module
      int hltN = (ix->second & 0x7FFFFFFF); // index of the HLT path the module belongs to
      lastModule = ( (ix->second & 0x80000000) != 0); // last module flag
      //cout<<" "<<num<<" "<<hltN<<" "<<lastN;
      if (num != i) cout<<" very very wrong1 "<<num<<" "<<i<<endl;
      if (hltN != hltIndex ) cout<<" very very wrong2 "<<hltN<<" "<<hltIndex<<endl;

    }

    if (anyTrig || lastModule || (anyModule&&activeModule) ) { // this object was selected, or use all
      if (localPrint) cout<<i<<" selected hlt "<<pTO->fHltPath<<" hlt-index "<<hltIndex<<" module label "
			 <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<endl;

      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      for(int n=0;n<num;++n) {  // loop over particles in this module, usually 2
	int index = muonIndex[n];
	int id = muonID[n];
	TLorentzVector p = muonP[n];

	if ( abs(id) != 13 ) { // if not muon trigger skip
	  //cout<<" match track "<<id<<" "<<pTO->fHltPath<<" "<<pTO->fLabel<<" "<<pTO->fType<<endl;
	  if (muonsOnly) continue;  // skip checking non-muon objects
	}

	// check direction matching
	double deltaR1 = p.DeltaR(tlvMu1);

	if (localPrint) {
	  cout<<" particle"<<n<<" index "<<index<<" id "<<id
	      <<" pt/eta/phi "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;
	  cout <<" mu1 "<< i<<" "<<pTO->fLabel << " "<<n <<" "<<deltaR1 <<endl;
	}

	// muon 1
	if (deltaR1<deltaRmin1) {
	  deltaRmin1=deltaR1;
	  mu1match = n;
	  hlt1 = pTO->fLabel;
	  indx1=i;

	  if (fVerbose > verboseThr || localPrint) {cout << " mu1 selected "<< deltaR1 <<endl;}

	  // check now the pt matching
	  double trigMatchDeltaPt=999.;
	  if (fp1->fPlab.Mag() > 0.) trigMatchDeltaPt = TMath::Abs(p.Rho()  - fp1->fPlab.Mag())/fp1->fPlab.Mag();
	  if ( trigMatchDeltaPt < deltaPtMatch ) trigMatchDeltaPt1=trigMatchDeltaPt;

	} // if direction match

      } // end for loop n
    } // end if valid module

  } // loop over all modules


  if (localPrint)
    cout << " best match "
	 <<indx1<<" "<< deltaRmin1 << " "<<mu1match<<" "<<hlt1<<" "<<trigMatchDeltaPt1<<endl;

  fTrigMatchDeltaPt=trigMatchDeltaPt1;

  if (mu1match>-1) {
    if (localPrint) cout<<" matching OK "<<deltaRmin1<<" "<<hlt1<<endl;
    //double tmp=fMu1TrigM;
    //fMu1TrigM = deltaRmin1;
    //if (tmp!=fMu1TrigM) cout<<" two methods inconsistent "<<tmp<<" "<<fMu1TrigM<<endl;
  }

  return deltaRmin1;
}

// ---------------------------------------------------------------------------------
// To match a single track to a trigger object (selected or all)
// pt - track
// anyTrig - if true use all trigger objects
// calls doTriggerMatchingR()
  bool candAna::doTriggerMatching(TAnaTrack *pt, bool anyTrig, bool muonsOnly, bool anyModule) {

  const double deltaRthrsh(0.02); // final cut, Frank had 0.5, change 0.020

  double dR = doTriggerMatchingR(pt,anyTrig,muonsOnly,anyModule);
  ((TH1D*)fHistDir->Get("test6"))->Fill(dR);

  bool HLTmatch = (dR<deltaRthrsh );
  return HLTmatch;

}

//--------------------------------------------------------------
bool candAna::doTriggerVeto(TAnaTrack *fp, bool muonsOnly, bool matchPt,
			    bool anyModule, float deltaRthr, int histoOffset) {

  cout<<" OBSOLETE, DO NOT USE "<<endl;

  // The valid code is in CandAnaDstar.cc
  return false;

}


// ----------------------------------------------------------------------
// -- check whether the reftrigger's objects are matched to the candidate's tracks
bool candAna::refTrigger(TAnaCand *pC, string refTriggerPath) {
  bool result(false);
  int verbose(0);

  // -- get list of indices of tracks making up candidate
  vector<int> sigIdx;
  getSigTracks(sigIdx, pC);

  if (verbose) {
    cout << "==> candAna::refTrigger> in DS = " << DSNAME
	 << " HLT = " << fGoodHLT << ", JSON = " << fJSON
	 << ", candidate " <<  pC->fType << " with tracks " << endl;
    for (unsigned int i = 0; i < sigIdx.size(); ++i) {
      cout << "muon = " << fpEvt->getSimpleTrack(sigIdx[i])->getMuonID()
	   << " " << Form(" %4d ", sigIdx[i])
	   << " pT/eta/phi = "
	   << fpEvt->getSimpleTrack(sigIdx[i])->getP().Perp() << " "
	   << fpEvt->getSimpleTrack(sigIdx[i])->getP().Eta() << " "
	   << fpEvt->getSimpleTrack(sigIdx[i])->getP().Phi() << " "
	   << endl;
    }
  }


  // -- determine trigger objects for reference path
  TTrgObjv2 *pTO(0);
  set<int> trgTrkIdx;
  for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {  // loop over all saved hlt objects
    pTO = fpEvt->getTrgObjv2(i);
    if (refTriggerPath == pTO->fHltPath) {
      if (!triggerFired(refTriggerPath)) {
	cout << "%^&^%&^%&%&^%&%&^%&^%&^%&^%&^%&  refTrigger in trigger objects, but not fired!" << endl;
      }
      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      // -- skip L1 and L2 objects (bad resolution for matching)
      if (pTO->fType.Contains("L1T")) continue;
      if (pTO->fType.Contains("L1Filter")) continue;
      if (pTO->fType.Contains("L2")) continue;
      if (verbose) cout << "  " << pTO->fHltPath << ": " << pTO->fType << " .. " << pTO->fLabel << "  " << " with n(particles) = " << num << endl;
      for (int j = 0; j < num; ++j) {
	double dr(0.);
	int trkIdx = matchTrgObj2Trk(muonP[j].Vect(), dr);
	if (trkIdx < 0) {
	  if (verbose) cout << "XXXXXXXXX NO MATCHING TRACK FOUND" << endl;
	  continue;
	}
	if (verbose) cout << "        " << muonP[j].Perp() << "/" << muonP[j].Eta() << "/" << muonP[j].Phi() << " muon? " << muonID[j]
			  << " matched to track idx " << trkIdx << " pt/eta/phi = "
			  << fpEvt->getSimpleTrack(trkIdx)->getP().Perp() << "/"
			  << fpEvt->getSimpleTrack(trkIdx)->getP().Eta() << "/"
			  << fpEvt->getSimpleTrack(trkIdx)->getP().Phi()
			  << " with dr = " << dr
			  << endl;
	trgTrkIdx.insert(trkIdx);
      }

    }
  }

  // -- now check that all objects of ref trigger are matched to the cand's muons
  set<int>::iterator it;
  int nmatch(0);
  for (it = trgTrkIdx.begin(); it != trgTrkIdx.end(); ++it) {
    for (unsigned int i = 0; i < sigIdx.size(); ++i) {
      if (*it == sigIdx[i]) {
	++nmatch;
	break;
      }
    }
  }

  if (verbose) cout << "REF TRIGGER nmatch = " << nmatch << endl;
  result = (nmatch == 2);
  return result;
}

// ----------------------------------------------------------------------
// -- calculate closest distance (deltaR) for a simpleTrack to a PD trigger object
// -- return by reference: deltaR, distance at M1, distance at M2
// -- FIXME: Implement calculating M1 and M2 distances (is this beneficial?)
void candAna::dist2PdTrigger(TSimpleTrack *pS, double &dr, double &dm1, double &dm2) {

  TTrgObjv2 *pTO(0);
  double result(99.);
  dr = dm1 = dm2 = -1.;
  int verbose(0);
  TVector3 pT = pS->getP();
  // TAnaMuon *pM = fpEvt->getSimpleTrackMuon(pS->getIndex());
  // TVector3 rm1, rm2;
  // if (pM) {
  //   rm1 = pM->fPositionAtM2;
  // } else {
  //   rm1.SetXYZ(0., 0., 0.);
  // }

  // cout << "--- trg objects -------------------------------------------------------------------" << endl;
  // for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {  // loop over all saved hlt objects
  //   pTO = fpEvt->getTrgObjv2(i);
  //   cout << "trgobjv2 type = " << pTO->fType << ", label = " << pTO->fLabel << " num = " << pTO->fP.size()
  // 	 << " pt/eta/phi = " << pTO->fP[0].Perp() << "/" << pTO->fP[0].Eta() << "/" << pTO->fP[0].Phi()
  // 	 << endl;
  //   continue;
  //   if (pTO->fType.Contains("l3muon")) {
  //     vector<int> muonIndex = pTO->fIndex;
  //     vector<int> muonID = pTO->fID;
  //     vector<TLorentzVector> muonP = pTO->fP;
  //     int num = muonIndex.size();
  //     cout << "  " << pTO->fHltPath << ": " << pTO->fType << " .. " << pTO->fLabel << "  " << " with n(particles) = " << num << endl;
  //     for (int j = 0; j < num; ++j) {
  // 	dr = muonP[j].Vect().DeltaR(pT);
  // 	cout << "        " << muonP[j].Perp() << "/" << muonP[j].Eta() << "/" << muonP[j].Phi() << " muon? " << muonID[j] << endl;
  //     }
  //   }
  // }


  // -- get list of PD triggers from histogram  e.g. triggers_Charmonium_run273730
  TH1D* ht = (TH1D*)fpReader->getFile()->Get(Form("triggers_%s_run%d", DSNAME.c_str(), static_cast<int>(fRun)));
  if (!ht) return;
  string hltPath("nada");
  if (verbose) cout << "==> candAna::tis> trigger objects for these paths" << endl;
  for (int j = 1; j <= ht->GetNbinsX(); ++j) {
    hltPath =  ht->GetXaxis()->GetBinLabel(j);
    // -- determine trigger objects for this path
    for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {  // loop over all saved hlt objects
      pTO = fpEvt->getTrgObjv2(i);
      if (hltPath == pTO->fHltPath) {
	vector<int> muonIndex = pTO->fIndex;
	vector<int> muonID = pTO->fID;
	vector<TLorentzVector> muonP = pTO->fP;
	int num = muonIndex.size();
	// -- skip L1 and L2 objects (bad resolution for matching)
	if (pTO->fType.Contains("L1Filter")) continue;
	if (pTO->fType.Contains("L1T")) continue;
	if (pTO->fType.Contains("L2")) continue;
	if (verbose) cout << "  " << pTO->fHltPath << ": " << pTO->fType << " .. " << pTO->fLabel << "  " << " with n(particles) = " << num << endl;
	double dr(99.);
	for (int j = 0; j < num; ++j) {
	  dr = muonP[j].Vect().DeltaR(pT);
	  if (dr < result) result = dr;
	  if (verbose > 0) {
	    cout << "hlt path =  " << hltPath << endl;
	    cout << "  " << pTO->fHltPath << ": " << pTO->fType << " .. " << pTO->fLabel << "  " << " with n(particles) = " << num << "/" << j << endl;
	    cout << "        " << muonP[j].Perp() << "/" << muonP[j].Eta() << "/" << muonP[j].Phi() << " muon? " << muonID[j]
		 << " distance to track " << " pt/eta/phi = "
		 << pT.Perp() << "/"
		 << pT.Eta() << "/"
		 << pT.Phi()
		 << " with dr = " << dr
		 << endl;
	  }
	}

      }
    }

  }
  dr = result;
}

// ----------------------------------------------------------------------
// -- search for a PD trigger that has no overlap with the tracks of the candidate
bool candAna::tis(TAnaCand *pC) {
  bool result(false);
  int verbose(0);

  // -- get list of indices of tracks making up candidate
  vector<int> sigIdx;
  getSigTracks(sigIdx, pC);

  // -- get list of PD triggers from histogram  e.g. triggers_Charmonium_run273730
  TH1D* ht = (TH1D*)fpReader->getFile()->Get(Form("triggers_%s_run%d", DSNAME.c_str(), static_cast<int>(fRun)));
  if (!ht) return false;
  string hltPath("nada");
  TTrgObjv2 *pTO(0);
  if (verbose) cout << "==> candAna::tis> trigger objects for these paths" << endl;
  map<string, set<int> > trgTrkIdx;
  TH1D *h1 = (TH1D*)(fHistDir->Get(Form("dr_%s", fName.c_str())));
  for (int j = 1; j <= ht->GetNbinsX(); ++j) {
    hltPath =  ht->GetXaxis()->GetBinLabel(j);
    // -- determine trigger objects for this path
    for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {  // loop over all saved hlt objects
      pTO = fpEvt->getTrgObjv2(i);
      if (hltPath == pTO->fHltPath) {
	vector<int> muonIndex = pTO->fIndex;
	vector<int> muonID = pTO->fID;
	vector<TLorentzVector> muonP = pTO->fP;
	int num = muonIndex.size();
	// -- skip L1 and L2 objects (bad resolution for matching)
	if (pTO->fType.Contains("L1Filter")) continue;
	if (pTO->fType.Contains("L1T")) continue;
	if (pTO->fType.Contains("L2")) continue;
	if (verbose) cout << "  " << pTO->fHltPath << ": " << pTO->fType << " .. " << pTO->fLabel << "  " << " with n(particles) = " << num << endl;
	for (int j = 0; j < num; ++j) {
	  double dr(0.);
	  int trkIdx = matchTrgObj2Trk(muonP[j].Vect(), dr);
	  h1->Fill(dr);
	  if (trkIdx < 0) {
	    if (verbose) cout << "XXXXXXXXX NO MATCHING TRACK FOUND" << endl;
	    continue;
	  }
	  if (verbose) cout << "        " << muonP[j].Perp() << "/" << muonP[j].Eta() << "/" << muonP[j].Phi() << " muon? " << muonID[j]
			    << " matched to track idx " << trkIdx << " pt/eta/phi = "
			    << fpEvt->getSimpleTrack(trkIdx)->getP().Perp() << "/"
			    << fpEvt->getSimpleTrack(trkIdx)->getP().Eta() << "/"
			    << fpEvt->getSimpleTrack(trkIdx)->getP().Phi()
			    << " with dr = " << dr
			    << endl;
	  trgTrkIdx[hltPath].insert(trkIdx);
	}

      }
    }

  }

  if (verbose) cout << "==> searching for non-overlapping trigger" << endl;
  map<string, set<int> >::iterator it;
  for (it = trgTrkIdx.begin(); it != trgTrkIdx.end(); ++it) {
    set<int>::iterator is;
    if (verbose) cout << it->first << " size = " << it->second.size() << ": ";
    bool overlap(false);
    for (is = it->second.begin(); is != it->second.end(); ++is) {
      for (unsigned int i = 0; i < sigIdx.size(); ++i) {
	if (*is == sigIdx[i]) {
	  if (verbose) cout << "(" << *is << ") ";
	  overlap = true;
	} else {
	  if (verbose) cout << " " << *is << " ";
	}
      }
    }
    if (!overlap) {
      result = true;
      if (verbose) {
	cout << " TIS trigger: NOT overlapping!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
      }
    }
    if (verbose) cout << endl;

  }



  return result;
}


// ----------------------------------------------------------------------
// -- check whether all trigger primitives are matched to the candidate's tracks
bool candAna::tos(TAnaCand *pC) {
  bool result(false);
  int verbose(0);

  if (fHLT1Path == "nada" || !fGoodHLT1) {
    if (verbose) cout << "event not triggered: fGoodHLT1 = " << fGoodHLT1 << " fHLT1Path = " << fHLT1Path << endl;
    return false;
  }
  // -- get list of indices of tracks making up candidate
  vector<int> sigIdx;
  getSigTracks(sigIdx, pC);

  string hltPath(fHLT1Path);
  TTrgObjv2 *pTO(0);
  if (verbose) cout << "==> candAna::tos> trigger objects for this path ->" << hltPath << "<-  for cand type = " << pC->fType << endl;
  // cout << "cand tracks = ";
  // for (unsigned int i = 0; i < sigIdx.size(); ++i) {
  //   cout << sigIdx[i] << " ";
  // }
  // cout << endl;

  map<string, set<int> > trgTrkIdx;
  // -- determine trigger objects for this path
  for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {  // loop over all saved hlt objects
    pTO = fpEvt->getTrgObjv2(i);
    if (hltPath == pTO->fHltPath) {
      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      // -- skip L1 and L2 objects (bad resolution for matching)
      if (pTO->fType.Contains("L1Filter")) continue;
      if (pTO->fType.Contains("L1T")) continue;
      if (pTO->fType.Contains("L2")) continue;
      if (verbose) cout << "  " << pTO->fHltPath << ": " << pTO->fType << " .. " << pTO->fLabel << "  " << " with n(particles) = " << num << endl;
      for (int j = 0; j < num; ++j) {
	double dr(0.);
	int trkIdx = matchTrgObj2Trk(muonP[j].Vect(), dr);
	if (trkIdx < 0) {
	  if (verbose) cout << "XXXXXXXXX NO MATCHING TRACK FOUND" << endl;
	  continue;
	}
	if (verbose) cout << "        " << muonP[j].Perp() << "/" << muonP[j].Eta() << "/" << muonP[j].Phi() << " muon? " << muonID[j]
			  << " matched to track idx " << trkIdx << " pt/eta/phi = "
			  << fpEvt->getSimpleTrack(trkIdx)->getP().Perp() << "/"
			  << fpEvt->getSimpleTrack(trkIdx)->getP().Eta() << "/"
			  << fpEvt->getSimpleTrack(trkIdx)->getP().Phi()
			  << " with dr = " << dr
			  << endl;
	trgTrkIdx[hltPath].insert(trkIdx);
      }

    }
  }

  if (verbose) cout << "==> searching for completely overlapping trigger, trgTrkIdx.size() = " << trgTrkIdx.size() << endl;
  map<string, set<int> >::iterator it;
  for (it = trgTrkIdx.begin(); it != trgTrkIdx.end(); ++it) {
    set<int>::iterator is;
    if (verbose) cout << it->first << " size = " << it->second.size() << ": ";
    int overlaps(0);
    for (is = it->second.begin(); is != it->second.end(); ++is) {
      for (unsigned int i = 0; i < sigIdx.size(); ++i) {
	if (*is == sigIdx[i]) {
	  ++overlaps;
	  if (verbose) cout << "(" << *is << ") ";
	} else {
	  if (verbose) cout << " " << *is << " ";
	}
      }
    }
    if (overlaps == it->second.size()) {
      result = true;
      if (verbose) {
	cout << " TOS trigger: COMPLETELY  overlapping!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
      }
    }
    if (verbose) cout << endl;

  }
  return result;
}


// ----------------------------------------------------------------------
int candAna::matchTrgObj2Trk(TVector3 t, double &dr) {
  double dRthrsh(0.3), dRmin(99.);
  int dRminIdx(-1);
  TVector3 p3;
  dr = -0.001;
  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
    p3 = fpEvt->getSimpleTrack(i)->getP();
    double dR = p3.DeltaR(t);
    double tPt = t.Perp();
    double pPt = p3.Perp();
    if (pPt > 0) pPt = tPt/pPt;
    if ((dR < dRthrsh) && (dR < dRmin) && (pPt > 0.5) && (pPt < 1.5)) {
      dRmin = dR;
      dRminIdx = i;
    }
  }
  if (dRmin < 98.) {
    dr = dRmin;
  } else {
    dr = -0.001;
  }
  return dRminIdx;
}


// ----------------------------------------------------------------------
bool candAna::triggerFired(std::string triggerPath) {
  for (int i = 0; i < NHLT; ++i) {
    if (fpEvt->fHLTResult[i]) {
      if (fpEvt->fHLTNames[i] == triggerPath) {
	return true;
      }
    }
  }
  return false;
}


// ----------------------------------------------------------------------
void candAna::boostGames() {

  if (fpMuon1 == NULL || fpMuon2 == NULL) return; // protection for DSTAR d.k. 15/1/2013

  double gcosTheta, gcosTheta2, rcosTheta, rcosTheta2;
  TVector3 pvec = TVector3(0., 0., 1.);
  if ((fGenBTmi > -1) && (fGenM1Tmi > -1) && (fGenM2Tmi > -1 )) {
    TGenCand *pB(0), *pM1(0), *pM2(0);

    pB  = fpEvt->getGenTWithIndex(fGenBTmi);
    pM1 = fpEvt->getGenTWithIndex(fGenM1Tmi);
    pM2 = fpEvt->getGenTWithIndex(fGenM2Tmi);

    TVector3 boost = pB->fP.Vect();
    boost.SetMag(boost.Mag()/pB->fP.E());
    TLorentzVector pM1Cms = pM1->fP;
    pM1Cms.Boost(-boost);
    TLorentzVector pM2Cms = pM2->fP;
    pM2Cms.Boost(-boost);

    //    N = P_{beam} x P_b / | P_{beam} x P_b |
    TVector3 bvec = pB->fP.Vect();
    TVector3 nvec = pvec.Cross(bvec);
    TLorentzVector nvec4; nvec4.SetXYZM(nvec.X(), nvec.Y(), nvec.Z(), 0);
    nvec4.Boost(-boost);

    if (pM1->fQ > 0) {
      gcosTheta = pM1Cms.CosTheta();
      gcosTheta2 = TMath::Cos(nvec4.Vect().Angle(pM1Cms.Vect()));

    } else {
      gcosTheta = pM2Cms.CosTheta();
      gcosTheta2 = TMath::Cos(nvec4.Vect().Angle(pM2Cms.Vect()));
    }

    if (0 == fNGenPhotons) {
      ((TH1D*)fHistDir->Get("gp1cms"))->Fill(pM1Cms.Rho());
      ((TH1D*)fHistDir->Get("gp2cms"))->Fill(pM2Cms.Rho());
      ((TH1D*)fHistDir->Get("gt1cms"))->Fill(gcosTheta);
      ((TH1D*)fHistDir->Get("gt2cms"))->Fill(gcosTheta2);
    } else {
      ((TH1D*)fHistDir->Get("gp1cmsg"))->Fill(pM1Cms.Rho());
      ((TH1D*)fHistDir->Get("gp2cmsg"))->Fill(pM2Cms.Rho());

      ((TH1D*)fHistDir->Get("gt1cmsg"))->Fill(gcosTheta);
      ((TH1D*)fHistDir->Get("gt2cmsg"))->Fill(gcosTheta2);
    }
  }




  // -- reco version
  TVector3 boost = fpCand->fPlab;
  double eboost  = TMath::Sqrt(fpCand->fPlab*fpCand->fPlab + fpCand->fMass*fpCand->fMass);
  boost.SetMag(boost.Mag()/eboost);
  TLorentzVector pM1Cms; pM1Cms.SetXYZM(fpMuon1->fPlab.X(), fpMuon1->fPlab.Y(), fpMuon1->fPlab.Z(), MMUON);
  pM1Cms.Boost(-boost);
  TLorentzVector pM2Cms; pM2Cms.SetXYZM(fpMuon2->fPlab.X(), fpMuon2->fPlab.Y(), fpMuon2->fPlab.Z(), MMUON);
  pM2Cms.Boost(-boost);

  //    N = P_{beam} x P_b / | P_{beam} x P_b |
  TVector3 bvec = fpCand->fPlab;
  TVector3 nvec = pvec.Cross(bvec);
  TLorentzVector nvec4;
  //  nvec4.SetXYZM(nvec.X(), nvec.Y(), nvec.Z(), 0);
  nvec4.SetXYZT(nvec.X(), nvec.Y(), nvec.Z(), 0);
  nvec4.Boost(-boost);


  if (fMu1Q > 0) {
    rcosTheta  = pM1Cms.CosTheta();
    rcosTheta2 = TMath::Cos(nvec4.Vect().Angle(pM1Cms.Vect()));
  } else {
    rcosTheta = pM2Cms.CosTheta();
    rcosTheta2 = TMath::Cos(nvec4.Vect().Angle(pM2Cms.Vect()));
  }



  ((TH2D*)fHistDir->Get("tvsm"))->Fill(fpCand->fMass, rcosTheta2);

  ((TH1D*)fHistDir->Get("rp1cms"))->Fill(pM1Cms.Rho());
  ((TH1D*)fHistDir->Get("rp2cms"))->Fill(pM2Cms.Rho());

  ((TH1D*)fHistDir->Get("rt1cms"))->Fill(rcosTheta);
  ((TH1D*)fHistDir->Get("rt2cms"))->Fill(rcosTheta2);
  if (fGoodHLT) {
    ((TH1D*)fHistDir->Get("rt3cms"))->Fill(rcosTheta2);
  }
  ((TH1D*)fHistDir->Get("gt1"))->Fill(rcosTheta-gcosTheta);
  ((TH1D*)fHistDir->Get("gt2"))->Fill(rcosTheta2-gcosTheta2);

  if (0)  cout << "muon 1:  p = " << fpMuon1->fPlab.Mag() << " " << pM1Cms.Rho()
	       << " muon 2: p = " << fpMuon2->fPlab.Mag() << " " << pM2Cms.Rho()
	       << " theta g = " << gcosTheta  << " r = " << rcosTheta
	       << " theta g2 = " << gcosTheta2  << " r = " << rcosTheta2
	       << endl;


}


// ----------------------------------------------------------------------
double candAna::distToMuon(TSimpleTrack *ps) {

  int numMuons(fpEvt->nMuons());

  TVector3 trackMom = ps->getP();
  int psIdx = ps->getIndex();

  TVector3 muonMom;
  TAnaMuon *pM(0);

  if (0) cout << "simple track index = " << ps->getIndex() << " with pt/eta/phi = "
	      << ps->getP().Perp() << "/" << ps->getP().Eta() << "/" << ps->getP().Phi()
	      << endl;

  int pmIdx(-1), bestIdx(-1);
  double ptMuon(0.), drMin(9999.), dr(0.);
  for (int im = 0; im < numMuons; ++im) {
    pM = fpEvt->getMuon(im);
    pmIdx = pM->fIndex;
    // skip if not global muon
    if ((pM->fMuID & 2) != 2) {
      if (0) cout << "muon with track index " << pmIdx << " is not a GM" << endl;
      continue;
    }
    // skip if same track index
    if (pmIdx == psIdx) {
      if (0) cout << "skipping muon with same track index " << pmIdx << endl;
      continue;
    }
    dr = pM->fPlab.DeltaR(trackMom);
    if (dr < drMin) {
      if (0) cout << "muon with track index " << pmIdx << " has smaller dr = " << dr << endl;
      drMin = dr;
      bestIdx = pmIdx;
    }
  }

  TSimpleTrack *s = fpEvt->getSimpleTrack(bestIdx);

  if (0) cout << " muon " << bestIdx
	      << " with pt/eta/phi = "
	      << s->getP().Perp() << "/" << s->getP().Eta() << "/" << s->getP().Phi()
	      << " has dr = " << drMin << endl;
  return drMin;
}



//-----------------------------------------------------------------------------------
// Loops over all muons, returns dR of the closests muon, excluding the
// same track muon (if exists)
double candAna::matchToMuon(TAnaTrack *pt, bool skipSame) {
  //bool print = true;
  bool print = false;

  int numMuons = fpEvt->nMuons();
  if (print) cout << "Found " << numMuons << " rec muons in event" << endl;

  TVector3 trackMom = pt->fPlab;  // test track momentum
  int it0 = pt->fIndex;
  if (print) cout<<" check track "<<it0<<endl;

  TVector3 muonMom;
  TAnaMuon * muon = 0;
  //TSimpleTrack * pTrack = 0;
  double ptMuon=0., dRMin=9999.;
  int select = -1;
  for (int it = 0; it<numMuons; ++it) { // loop over muons
    muon = fpEvt->getMuon(it);
    //if (print) muon->dump();

    // check if this is a nice  muon, accept only global and tracker muons
    int muonId = muon->fMuID;
    //     if ( (muonId & 0x6) == 0 ) continue;  // skip muons which are not global/tracker

    int itrk = muon->fIndex;
    // Eliminate pure standalone muons and calo muons. Skip same track comparion (only of skipSame=true)
    if (itrk<0 || (skipSame && itrk == it0)) {
      if (print) {
	if (itrk<0) cout<<"standalone only or calo muon? "<<hex<<muonId<<dec<<" "<<(muon->fPlab).Perp()<<endl;
	else       cout<<"skip, same track "<<endl;
      }
      continue;  // skip same track comparion and standalone/calo muons
    }

    //if (itrk>0) cout<<" tracker muon? "<<hex<<(muonId&0x6)<<dec<<" "<<(muon->fPlab).Perp()<<endl;

    // Use direct access, withour going through SimpleTracks
    muonMom = muon->fPlab;
    ptMuon  = muonMom.Perp();
    double dR = muonMom.DeltaR(trackMom);
    //double etaMuon = muonMom.Eta();
    //double phiMuon = muonMom.Phi();
    if (print) cout<<it<<" "<<ptMuon<<" "<<dR<<endl;

    // Go through reco track, Find the reco track  NOT NEEDED
//     if (itrk>=0 && itrk< (fpEvt->nSimpleTracks()) ) {  // if the simple track exists
//       //pTrack = fpEvt->getRecTrack(itrk);
//       pTrack = fpEvt->getSimpleTrack(itrk);
//       //cout<<it<<" "<<itrk<<" "<<pTrack<<endl;
//       if (pTrack != 0) {
// 	muonMom = pTrack->getP();
//  	ptMuon  = muonMom.Perp();
// 	dR = muonMom.DeltaR(trackMom);
// 	if (print) cout<<it<<" "<<ptMuon<<" "<<dR<<endl;
//       } //if ptrack
//     } // if track

    if (dR<dRMin) {dRMin=dR; select=it;} // select the best fit

  } // loop over muons

  if (select>-1) {
    int idx =  (fpEvt->getMuon(select))->fIndex;  // find the track index
    if (print) cout<<"final muon match "<<dRMin<<" "<<select<<" "<<idx<<endl;

    ((TH1D*)fHistDir->Get("test7"))->Fill(dRMin);
  }

  return dRMin;
}


// ----------------------------------------------------------------------
void candAna::play3() {

  if (!fJSON) return;
  cout << "Evt: " << fEvt << " ----------------------------------------------------------------------" << endl;
  static int first(1);
  if (first) {
    first = 0;
    fHistDir->cd();

    vector<string> particles;
    particles.push_back("bu");
    particles.push_back("bs");
    particles.push_back("ks");
    particles.push_back("mm");
    particles.push_back("dstar");

    for (unsigned int i = 0; i < particles.size(); ++i) {
      new TH1D(Form("%sflsxy", particles[i].c_str()), Form("%s flsxy", particles[i].c_str()), 100, 0., 50.);
      if ("ks" == particles[i]) {
	new TH1D(Form("%smass", particles[i].c_str()),  Form("%s mass", particles[i].c_str()), 100, 0.4, 0.6);
      } else if ("dstar" == particles[i]) {
	new TH1D(Form("%smass", particles[i].c_str()),  Form("%s mass", particles[i].c_str()), 100, 1.8, 2.2);
	new TH1D(Form("%sdm", particles[i].c_str()),  Form("%s dm", particles[i].c_str()), 100, 0.13, 0.16);
	new TH1D(Form("%smd0", particles[i].c_str()),  Form("%s md0", particles[i].c_str()), 100, 1.6, 2.0);
      } else if ("bs" == particles[i]) {
	new TH1D(Form("%smass", particles[i].c_str()),  Form("%s mass", particles[i].c_str()), 100, 4.8, 6.0);
      } else if ("mm" == particles[i]) {
	new TH1D(Form("%smass", particles[i].c_str()),  Form("%s mass", particles[i].c_str()), 100, 4.2, 6.7);
      } else if ("bu" == particles[i]) {
	new TH1D(Form("%smass", particles[i].c_str()),  Form("%s mass", particles[i].c_str()), 100, 4.8, 6.0);
      }
      new TH1D(Form("%schi2", particles[i].c_str()),  Form("%schi2", particles[i].c_str()), 100, 0., 20.);
      new TH1D(Form("%sschi2", particles[i].c_str()),  Form("%sschi2", particles[i].c_str()), 100, 0., 20.);
      new TH1D(Form("%spvips", particles[i].c_str()), Form("%spvips", particles[i].c_str()), 100, 0., 5.);
      new TH1D(Form("%sfls3d", particles[i].c_str()), Form("%sfls3d", particles[i].c_str()), 100, 0., 50.);
      new TH1D(Form("%salpha", particles[i].c_str()), Form("%salpha", particles[i].c_str()), 100, 0., .4);
      new TH1D(Form("%sdocatrk", particles[i].c_str()), Form("%sdocatrk", particles[i].c_str()), 100, 0., .1);
      new TH1D(Form("%stip", particles[i].c_str()),   Form("%stip", particles[i].c_str()), 100, 0., .02);
      new TH1D(Form("%slip", particles[i].c_str()),   Form("%slip", particles[i].c_str()), 100, -0.02, .02);
      new TH1D(Form("%spvw8", particles[i].c_str()),   Form("%spvw8", particles[i].c_str()), 100, 0., 1.0);
    }
  }


  TAnaCand *pCand(0), *sCand(0);
  // cout << "----------------------------------------------------------------------" << endl;
  // cout << "-- Event " << fEvt << endl;
  // cout << "----------------------------------------------------------------------" << endl;

  for (int iC = 0; iC < fpEvt->nCands(); ++iC) {
    pCand = fpEvt->getCand(iC);

    double PvAveW8(0.);
    if (pCand->fPvIdx > -1 && pCand->fPvIdx < fpEvt->nPV()) {
      TAnaVertex *pv = fpEvt->getPV(pCand->fPvIdx);
      PvAveW8 = ((pv->fNdof+2.)/2.)/pv->getNtracks();
    }

    if (300521 == pCand->fType) {
      if (pCand->fMaxDoca > 0.06) continue;
      if (pCand->fPvIP3d/pCand->fPvIP3dE > 5) continue;
      if (pCand->fVtx.fDxy/pCand->fVtx.fDxyE < 3) continue;
      if ((pCand->fMass < 4.8) || (pCand->fMass > 6.0)) continue;
      sCand = fpEvt->getCand(pCand->fDau1);
      if (sCand->fMaxDoca > 0.06) continue;
      if ((sCand->fMass < 2.9) || (sCand->fMass > 3.2)) continue;
      if (sCand->fPlab.Perp() < 6.9) continue;

      ((TH1D*)fHistDir->Get("buflsxy"))->Fill(pCand->fVtx.fDxy/pCand->fVtx.fDxyE);
      ((TH1D*)fHistDir->Get("bumass"))->Fill(pCand->fMass);
      ((TH1D*)fHistDir->Get("bupvips"))->Fill(pCand->fPvIP3d/pCand->fPvIP3dE);
      ((TH1D*)fHistDir->Get("buchi2"))->Fill(pCand->fVtx.fChi2/pCand->fVtx.fNdof);
      ((TH1D*)fHistDir->Get("buschi2"))->Fill(sCand->fVtx.fChi2/sCand->fVtx.fNdof);
      ((TH1D*)fHistDir->Get("bufls3d"))->Fill(pCand->fVtx.fD3d/pCand->fVtx.fD3dE);
      ((TH1D*)fHistDir->Get("bualpha"))->Fill(pCand->fAlpha);
      ((TH1D*)fHistDir->Get("budocatrk"))->Fill(pCand->fNstTracks[0].second.first);
      ((TH1D*)fHistDir->Get("butip"))->Fill(pCand->fPvTip);
      ((TH1D*)fHistDir->Get("bulip"))->Fill(pCand->fPvLip);
      ((TH1D*)fHistDir->Get("bupvw8"))->Fill(PvAveW8);

    } else if (300531 == pCand->fType) {
      if (pCand->fMaxDoca > 0.06) continue;
      if (pCand->fPvIP3d/pCand->fPvIP3dE > 5) continue;
      if (pCand->fVtx.fDxy/pCand->fVtx.fDxyE < 2.5) continue;
      if (pCand->fVtx.fChi2 > 10.0) continue;
      if ((pCand->fMass < 4.8) || (pCand->fMass > 6.0)) continue;
      sCand = fpEvt->getCand(pCand->fDau1);
      if (sCand->fMaxDoca > 0.06) continue;
      if ((sCand->fMass < 2.9) || (sCand->fMass > 3.2)) continue;
      if (sCand->fPlab.Perp() < 6.9) continue;
      double schi2 = pCand->fVtx.fChi2/pCand->fVtx.fNdof;
      sCand = fpEvt->getCand(pCand->fDau2);
      if ((sCand->fMass < 0.98) || (sCand->fMass > 1.06)) continue;

      ((TH1D*)fHistDir->Get("bsflsxy"))->Fill(pCand->fVtx.fDxy/pCand->fVtx.fDxyE);
      ((TH1D*)fHistDir->Get("bsmass"))->Fill(pCand->fMass);
      ((TH1D*)fHistDir->Get("bspvips"))->Fill(pCand->fPvIP3d/pCand->fPvIP3dE);
      ((TH1D*)fHistDir->Get("bschi2"))->Fill(pCand->fVtx.fChi2/pCand->fVtx.fNdof);
      ((TH1D*)fHistDir->Get("bsschi2"))->Fill(schi2);
      ((TH1D*)fHistDir->Get("bsfls3d"))->Fill(pCand->fVtx.fD3d/pCand->fVtx.fD3dE);
      ((TH1D*)fHistDir->Get("bsalpha"))->Fill(pCand->fAlpha);
      ((TH1D*)fHistDir->Get("bsdocatrk"))->Fill(pCand->fNstTracks[0].second.first);
      ((TH1D*)fHistDir->Get("bstip"))->Fill(pCand->fPvTip);
      ((TH1D*)fHistDir->Get("bslip"))->Fill(pCand->fPvLip);
      ((TH1D*)fHistDir->Get("bspvw8"))->Fill(PvAveW8);
    } else if (11310 == pCand->fType) {
      if (pCand->fMaxDoca > 0.1) continue;
      if (pCand->fPvIP3d/pCand->fPvIP3dE > 5) continue;
      if (pCand->fVtx.fDxy/pCand->fVtx.fDxyE < 5.) continue;
      if (pCand->fVtx.fDxy > 4.) continue;
      if (pCand->fVtx.fChi2 > 20.) continue;
      if ((pCand->fMass < 0.44) || (pCand->fMass > 0.56)) continue;
      double pt1 = fpEvt->getSigTrack(pCand->fSig1)->fPlab.Perp();
      double pt2 = fpEvt->getSigTrack(pCand->fSig2)->fPlab.Perp();
      if (pt2 > pt1) {
	double pt0 = pt1;
	pt1 = pt2;
	pt2 = pt0;
      }
      if (pt1 < 3.5) continue;
      if (pt2 < 2.5) continue;

      ((TH1D*)fHistDir->Get("ksflsxy"))->Fill(pCand->fVtx.fDxy/pCand->fVtx.fDxyE);
      ((TH1D*)fHistDir->Get("ksmass"))->Fill(pCand->fMass);
      ((TH1D*)fHistDir->Get("kspvips"))->Fill(pCand->fPvIP3d/pCand->fPvIP3dE);
      ((TH1D*)fHistDir->Get("kschi2"))->Fill(pCand->fVtx.fChi2/pCand->fVtx.fNdof);
      ((TH1D*)fHistDir->Get("ksfls3d"))->Fill(pCand->fVtx.fD3d/pCand->fVtx.fD3dE);
      ((TH1D*)fHistDir->Get("ksalpha"))->Fill(pCand->fAlpha);
      ((TH1D*)fHistDir->Get("ksdocatrk"))->Fill(pCand->fNstTracks[0].second.first);
      ((TH1D*)fHistDir->Get("kstip"))->Fill(pCand->fPvTip);
      ((TH1D*)fHistDir->Get("kslip"))->Fill(pCand->fPvLip);
      ((TH1D*)fHistDir->Get("kspvw8"))->Fill(PvAveW8);

      cout << Form("==> Cand: idx=%3d type=%d m=%5.3f+/-%5.3f, pT=%6.2f eta=%+4.3f f=%+4.3f maxdoca=%5.4f",
		   pCand->fIndex, pCand->fType, pCand->fMass, pCand->fMassE,
		   pCand->fPlab.Perp(), pCand->fPlab.Eta(), pCand->fPlab.Phi(), pCand->fMaxDoca) << endl;
      cout << "daughter cands: " << pCand->fDau1 << " .. " << pCand->fDau2
	   << " sig tracks: ";
      if (pCand->fSig1 > -1) cout << fpEvt->getSigTrack(pCand->fSig1)->fIndex;
      if (pCand->fSig2 > -1) cout << " .. " << fpEvt->getSigTrack(pCand->fSig2)->fIndex;
      cout << " pt= " << pt1 << "/" << pt2;
      cout << " flxy=" << pCand->fVtx.fDxy << " flsxy=" << pCand->fVtx.fDxy/pCand->fVtx.fDxyE;
      cout << " pvips= " << pCand->fPvIP3d/pCand->fPvIP3dE;
      cout << endl;
      pCand->fVtx.dump();
      cout << " referring to PV " << pCand->fPvIdx << " at "
	   << fpEvt->getPV(pCand->fPvIdx)->fPoint.X() << "/"
	   << fpEvt->getPV(pCand->fPvIdx)->fPoint.Y() << "/"
	   << fpEvt->getPV(pCand->fPvIdx)->fPoint.Z()
	   << " with pvw8 = " << ((fpEvt->getPV(pCand->fPvIdx)->fNdof+2.)/2.)/fpEvt->getPV(pCand->fPvIdx)->getNtracks()
	   << endl;
      // fPvNtrk = pv->getNtracks();
      // fPvNdof = pv->fNdof;


      // sCand = fpEvt->getCand(pCand->fDau1);
      // cout << Form("=> J/psi cand: idx=%3d type=%d m=%5.3f+/-%5.3f, pT=%6.2f f=%+4.3f eta=%+4.3f maxdoca=%5.4f",
      // 		   sCand->fIndex, sCand->fType, sCand->fMass, sCand->fMassE,
      // 		   sCand->fPlab.Perp(), sCand->fPlab.Phi(), sCand->fPlab.Eta(), sCand->fMaxDoca) << endl;
      // cout << " sig tracks: ";
      // if (sCand->fSig1 > -1) cout << fpEvt->getSigTrack(sCand->fSig1)->fIndex;
      // if (sCand->fSig2 > -1) cout << " .. " << fpEvt->getSigTrack(sCand->fSig2)->fIndex;
      // cout << " flxy=" << sCand->fVtx.fDxy << " flsxy=" << sCand->fVtx.fDxy/sCand->fVtx.fDxyE;
      // cout << " pvips= " << sCand->fPvIP3d/sCand->fPvIP3dE;
      // cout << endl;
      // sCand->fVtx.dump();

      // sCand = fpEvt->getCand(pCand->fDau2);
      // cout << Form("=> phi cand: idx=%3d type=%d m=%5.3f+/-%5.3f, pT=%6.2f f=%+4.3f eta=%+4.3f maxdoca=%5.4f",
      // 		   sCand->fIndex, sCand->fType, sCand->fMass, sCand->fMassE,
      // 		   sCand->fPlab.Perp(), sCand->fPlab.Phi(), sCand->fPlab.Eta(), sCand->fMaxDoca) << endl;
      // cout << " sig tracks: ";
      // if (sCand->fSig1 > -1) cout << fpEvt->getSigTrack(sCand->fSig1)->fIndex;
      // if (sCand->fSig2 > -1) cout << " .. " << fpEvt->getSigTrack(sCand->fSig2)->fIndex;
      // cout << " flxy=" << sCand->fVtx.fDxy << " flsxy=" << sCand->fVtx.fDxy/sCand->fVtx.fDxyE;
      // cout << " pvips= " << sCand->fPvIP3d/sCand->fPvIP3dE;
      // cout << endl;
      // sCand->fVtx.dump();

    } else if (1313 == pCand->fType) {
      if (pCand->fPlab.Perp() < 5.0) continue;
      if (pCand->fMaxDoca > 0.06) continue;
      if (pCand->fPvIP3d/pCand->fPvIP3dE > 5) continue;
      if (pCand->fVtx.fDxy/pCand->fVtx.fDxyE < 2.5) continue;
      if ((pCand->fMass < 4.2) || (pCand->fMass > 6.7)) continue;
      double pt1 = fpEvt->getSigTrack(pCand->fSig1)->fPlab.Perp();
      double pt2 = fpEvt->getSigTrack(pCand->fSig2)->fPlab.Perp();
      if (pt2 > pt1) {
	double pt0 = pt1;
	pt1 = pt2;
	pt2 = pt0;
      }
      if (pt1 < 4.0) continue;
      if (pt2 < 4.0) continue;
      ((TH1D*)fHistDir->Get("mmflsxy"))->Fill(pCand->fVtx.fDxy/pCand->fVtx.fDxyE);
      ((TH1D*)fHistDir->Get("mmmass"))->Fill(pCand->fMass);
      ((TH1D*)fHistDir->Get("mmpvips"))->Fill(pCand->fPvIP3d/pCand->fPvIP3dE);
      ((TH1D*)fHistDir->Get("mmchi2"))->Fill(pCand->fVtx.fChi2/pCand->fVtx.fNdof);
      ((TH1D*)fHistDir->Get("mmfls3d"))->Fill(pCand->fVtx.fD3d/pCand->fVtx.fD3dE);
      ((TH1D*)fHistDir->Get("mmalpha"))->Fill(pCand->fAlpha);
      ((TH1D*)fHistDir->Get("mmdocatrk"))->Fill(pCand->fNstTracks[0].second.first);
      ((TH1D*)fHistDir->Get("mmtip"))->Fill(pCand->fPvTip);
      ((TH1D*)fHistDir->Get("mmlip"))->Fill(pCand->fPvLip);
      ((TH1D*)fHistDir->Get("mmpvw8"))->Fill(PvAveW8);
    } else if (300054 == pCand->fType) {
      if (fEvt == 32579099) {
	sCand = fpEvt->getCand(pCand->fDau1);  // D0 candidate
	double dm = pCand->fMass - sCand->fMass;
	cout << "dm = " << dm << endl;
	cout << Form("Cand: idx=%3d type=%d m=%5.3f+/-%5.3f, pT=%6.2f eta=%+4.3f f=%+4.3f maxdoca=%5.4f",
		     pCand->fIndex, pCand->fType, pCand->fMass, pCand->fMassE,
		     pCand->fPlab.Perp(), pCand->fPlab.Eta(), pCand->fPlab.Phi(), pCand->fMaxDoca) << endl;
	cout << "daughter cands: " << pCand->fDau1 << " .. " << pCand->fDau2
	     << " sig tracks: ";
	if (pCand->fSig1 > -1) cout << fpEvt->getSigTrack(pCand->fSig1)->fIndex;
	if (pCand->fSig2 > -1) cout << " .. " << fpEvt->getSigTrack(pCand->fSig2)->fIndex;
	cout << " flxy=" << pCand->fVtx.fDxy << " flsxy=" << pCand->fVtx.fDxy/pCand->fVtx.fDxyE;
	cout << " pvips= " << pCand->fPvIP3d/pCand->fPvIP3dE;
	cout << endl;
	pCand->fVtx.dump();
	cout << " referring to PV " << pCand->fPvIdx << " at "
	     << fpEvt->getPV(pCand->fPvIdx)->fPoint.X() << "/"
	     << fpEvt->getPV(pCand->fPvIdx)->fPoint.Y() << "/"
	     << fpEvt->getPV(pCand->fPvIdx)->fPoint.Z()
	     << endl;
	cout << Form(" D0 cand: idx=%3d type=%d m=%5.3f+/-%5.3f, pT=%6.2f f=%+4.3f eta=%+4.3f maxdoca=%5.4f",
		     sCand->fIndex, sCand->fType, sCand->fMass, sCand->fMassE,
		     sCand->fPlab.Perp(), sCand->fPlab.Phi(), sCand->fPlab.Eta(), sCand->fMaxDoca) << endl;
	cout << " sig tracks: ";
	if (sCand->fSig1 > -1) cout << fpEvt->getSigTrack(sCand->fSig1)->fIndex;
	if (sCand->fSig2 > -1) cout << " .. " << fpEvt->getSigTrack(sCand->fSig2)->fIndex;
	cout << " flxy=" << sCand->fVtx.fDxy << " flsxy=" << sCand->fVtx.fDxy/sCand->fVtx.fDxyE;
	cout << " pvips= " << sCand->fPvIP3d/sCand->fPvIP3dE;
	cout << endl;
	sCand->fVtx.dump();
      }
      if (pCand->fMaxDoca > 0.1) continue;
      if ((pCand->fMass < 1.9) || (pCand->fMass > 2.2)) continue;
      sCand = fpEvt->getCand(pCand->fDau1);  // D0 candidate
      if ((sCand->fMass < 1.75) || (sCand->fMass > 2.0)) continue;
      double dm = pCand->fMass - sCand->fMass;
      if (dm > 0.17) continue;
      ((TH1D*)fHistDir->Get("dstarflsxy"))->Fill(pCand->fVtx.fDxy/pCand->fVtx.fDxyE);
      ((TH1D*)fHistDir->Get("dstarmass"))->Fill(pCand->fMass);
      ((TH1D*)fHistDir->Get("dstardm"))->Fill(dm);
      ((TH1D*)fHistDir->Get("dstarmd0"))->Fill(sCand->fMass);
      ((TH1D*)fHistDir->Get("dstarpvips"))->Fill(pCand->fPvIP3d/pCand->fPvIP3dE);
      ((TH1D*)fHistDir->Get("dstarchi2"))->Fill(pCand->fVtx.fChi2/pCand->fVtx.fNdof);
      ((TH1D*)fHistDir->Get("dstarfls3d"))->Fill(pCand->fVtx.fD3d/pCand->fVtx.fD3dE);
      ((TH1D*)fHistDir->Get("dstaralpha"))->Fill(pCand->fAlpha);
      ((TH1D*)fHistDir->Get("dstardocatrk"))->Fill(pCand->fNstTracks[0].second.first);
      ((TH1D*)fHistDir->Get("dstartip"))->Fill(pCand->fPvTip);
      ((TH1D*)fHistDir->Get("dstarlip"))->Fill(pCand->fPvLip);
      ((TH1D*)fHistDir->Get("dstarpvw8"))->Fill(PvAveW8);

      cout << "===> event " << fEvt << "  ";
      cout << Form("Cand: idx=%3d type=%d m=%5.3f+/-%5.3f, pT=%6.2f f=%+4.3f eta=%+4.3f maxdoca=%5.4f",
		   pCand->fIndex, pCand->fType, pCand->fMass, pCand->fMassE,
		   pCand->fPlab.Perp(), pCand->fPlab.Phi(), pCand->fPlab.Eta(), pCand->fMaxDoca) << endl;
      cout << "daughter cands: " << pCand->fDau1 << " .. " << pCand->fDau2
	   << " sig tracks: ";
      if (pCand->fSig1 > -1) cout << fpEvt->getSigTrack(pCand->fSig1)->fIndex;
      if (pCand->fSig2 > -1) cout << " .. " << fpEvt->getSigTrack(pCand->fSig2)->fIndex;
      cout << " flxy=" << pCand->fVtx.fDxy << " flsxy=" << pCand->fVtx.fDxy/pCand->fVtx.fDxyE;
      cout << " pvips= " << pCand->fPvIP3d/pCand->fPvIP3dE;
      cout << endl;
      pCand->fVtx.dump();
      cout << Form(" D0 cand: idx=%3d type=%d m=%5.3f+/-%5.3f, pT=%6.2f f=%+4.3f eta=%+4.3f maxdoca=%5.4f",
		   sCand->fIndex, sCand->fType, sCand->fMass, sCand->fMassE,
		   sCand->fPlab.Perp(), sCand->fPlab.Phi(), sCand->fPlab.Eta(), sCand->fMaxDoca) << endl;
      cout << " sig tracks: ";
      if (sCand->fSig1 > -1) cout << fpEvt->getSigTrack(sCand->fSig1)->fIndex;
      if (sCand->fSig2 > -1) cout << " .. " << fpEvt->getSigTrack(sCand->fSig2)->fIndex;
      cout << endl;
    } else if (300443 == pCand->fType) {
    } else if (300333 == pCand->fType) {
    } else if (400521 == pCand->fType) {
    } else if (400531 == pCand->fType) {
    } else if (400443 == pCand->fType) {
    } else if (400333 == pCand->fType) {
    } else {
    }
  }
}


// ----------------------------------------------------------------------
void candAna::triggerEff(std::string ref, std::string os, int mode) {

  TH1D *h1(0);
  (void*)h1;
  string tname = Form("os_%d", mode);
  if (0 == ((TH1D*)gFile->Get(Form("%s_ptp", tname.c_str())))) {
    TDirectory *pDir = gDirectory;
    gFile->cd();
    cout << "triggerEff booking hists for mode = " << mode << endl;
    h1 = new TH1D(Form("%s_ptp", tname.c_str()), "pt (pass)", 50, 0., 50.);
    h1 = new TH1D(Form("%s_pta", tname.c_str()), "pt (all)", 50, 0., 50.);
    h1 = new TH1D(Form("%s_ptp1", tname.c_str()), "pt muon1 (pass)", 50, 0., 50.);
    h1 = new TH1D(Form("%s_pta1", tname.c_str()), "pt muon1 (all)",  50, 0., 50.);
    h1 = new TH1D(Form("%s_ptp2", tname.c_str()), "pt muon2 (pass)", 50, 0., 50.);
    h1 = new TH1D(Form("%s_pta2", tname.c_str()), "pt muon2 (all)",  50, 0., 50.);

    h1 = new TH1D(Form("%s_etap", tname.c_str()), "eta (pass)", 50, -2.5, 2.5);
    h1 = new TH1D(Form("%s_etaa", tname.c_str()), "eta (all)",  50, -2.5, 2.5);
    h1 = new TH1D(Form("%s_etap1", tname.c_str()), "eta muon1 (pass)", 50, -2.5, 2.5);
    h1 = new TH1D(Form("%s_etaa1", tname.c_str()), "eta muon1 (all)",  50, -2.5, 2.5);
    h1 = new TH1D(Form("%s_etap2", tname.c_str()), "eta muon2 (pass)", 50, -2.5, 2.5);
    h1 = new TH1D(Form("%s_etaa2", tname.c_str()), "eta muon2 (all)",  50, -2.5, 2.5);

    h1 = new TH1D(Form("%s_phip", tname.c_str()), "phi (pass)", 50, -3.15, 3.15);
    h1 = new TH1D(Form("%s_phia", tname.c_str()), "phi (all)",  50, -3.15, 3.15);
    h1 = new TH1D(Form("%s_phip1", tname.c_str()), "phi muon1 (pass)", 50, -3.15, 3.15);
    h1 = new TH1D(Form("%s_phia1", tname.c_str()), "phi muon1 (all)",  50, -3.15, 3.15);
    h1 = new TH1D(Form("%s_phip2", tname.c_str()), "phi muon2 (pass)", 50, -3.15, 3.15);
    h1 = new TH1D(Form("%s_phia2", tname.c_str()), "phi muon2 (all)",  50, -3.15, 3.15);

    pDir->cd();
  }

  bool refTrigger(false), osTrigger(false);
  for (int i = 0; i < NHLT; ++i) {
    if ((fpEvt->fHLTNames[i] == ref) && (fpEvt->fHLTResult[i])) {
      refTrigger = true;
    }
    if ((fpEvt->fHLTNames[i] == os) && (fpEvt->fHLTResult[i])) {
      osTrigger = true;
    }
  }

  if (!refTrigger) return;

  TAnaCand *pC(0);
  int m1(-1), m2(-1);
  double m1Pt(0.), m2Pt(0.), m1Eta(-99.), m2Eta(-66.), dEta(99.), m1Phi(-66.), m2Phi(-66.);
  for (int i = 0; i < fpEvt->nCands(); ++i) {
    pC = fpEvt->getCand(i);
    if ((pC->fType == 300511) || (pC->fType == 300531) ||(pC->fType == 300521)) {
      vector<int> sigIdx;
      getSigTracks(sigIdx, pC);
      m1 = m2 = -1;
      for (unsigned int i = 0; i < sigIdx.size(); ++i) {
	if (fpEvt->getSimpleTrack(sigIdx[i])->getMuonID()) {
	  if (m1 == -1) {
	    m1 = sigIdx[i];
	    m1Pt  = fpEvt->getSimpleTrack(sigIdx[i])->getP().Perp();
	    m1Eta = fpEvt->getSimpleTrack(sigIdx[i])->getP().Eta();
	    m1Phi = fpEvt->getSimpleTrack(sigIdx[i])->getP().Phi();
	  } else {
	    m2 = sigIdx[i];
	    m2Pt  = fpEvt->getSimpleTrack(sigIdx[i])->getP().Perp();
	    m2Eta = fpEvt->getSimpleTrack(sigIdx[i])->getP().Eta();
	    m2Phi = fpEvt->getSimpleTrack(sigIdx[i])->getP().Phi();
	  }
	}
      }

      if (m1Pt < m2Pt) {
	double bla = m1Pt;
	m1Pt = m2Pt;
	m2Pt = bla;

	bla = m1Eta;
	m1Eta = m2Eta;
	m2Eta = bla;

	bla = m1Phi;
	m1Phi = m2Phi;
	m2Phi = bla;
      }

      dEta = TMath::Abs(m1Eta - m2Eta);

      double flsxy = pC->fVtx.fDxy/pC->fVtx.fDxyE;
      if (1 == mode) {
	if (!(flsxy > 3.0 && (m1Pt > 4.0) && (m2Pt > 3.0) && (dEta < 1.8) && (TMath::Abs(m1Eta) < 1.6) && (TMath::Abs(m2Eta) < 1.6))) continue;
      }

      if (2 == mode) {
	if (!(flsxy > 3.0 && (m1Pt > 4.0) && (m2Pt > 3.0) && (dEta < 1.8) && (TMath::Abs(m1Eta) < 1.6) && (TMath::Abs(m2Eta) < 1.6))) continue;
      }

      tname = Form("os_%d", mode);
      ((TH1D*)gFile->Get(Form("%s_pta", tname.c_str())))->Fill(pC->fPlab.Perp());
      ((TH1D*)gFile->Get(Form("%s_pta1", tname.c_str())))->Fill(m1Pt);
      ((TH1D*)gFile->Get(Form("%s_pta2", tname.c_str())))->Fill(m2Pt);
      ((TH1D*)gFile->Get(Form("%s_etaa", tname.c_str())))->Fill(pC->fPlab.Eta());
      ((TH1D*)gFile->Get(Form("%s_etaa1", tname.c_str())))->Fill(m1Eta);
      ((TH1D*)gFile->Get(Form("%s_etaa2", tname.c_str())))->Fill(m2Eta);
      ((TH1D*)gFile->Get(Form("%s_phia", tname.c_str())))->Fill(pC->fPlab.Phi());
      ((TH1D*)gFile->Get(Form("%s_phia1", tname.c_str())))->Fill(m1Phi);
      ((TH1D*)gFile->Get(Form("%s_phia2", tname.c_str())))->Fill(m2Phi);

      if (osTrigger) {
	((TH1D*)gFile->Get(Form("%s_ptp", tname.c_str())))->Fill(pC->fPlab.Perp());
	((TH1D*)gFile->Get(Form("%s_ptp1", tname.c_str())))->Fill(m1Pt);
	((TH1D*)gFile->Get(Form("%s_ptp2", tname.c_str())))->Fill(m2Pt);
	((TH1D*)gFile->Get(Form("%s_etap", tname.c_str())))->Fill(pC->fPlab.Eta());
	((TH1D*)gFile->Get(Form("%s_etap1", tname.c_str())))->Fill(m1Eta);
	((TH1D*)gFile->Get(Form("%s_etap2", tname.c_str())))->Fill(m2Eta);
	((TH1D*)gFile->Get(Form("%s_phip", tname.c_str())))->Fill(pC->fPlab.Phi());
	((TH1D*)gFile->Get(Form("%s_phip1", tname.c_str())))->Fill(m1Phi);
	((TH1D*)gFile->Get(Form("%s_phip2", tname.c_str())))->Fill(m2Phi);
      }

      break; // fill only for one candidate with fulfilling J/psi muons and displacement
    }
  }

}


// ----------------------------------------------------------------------
void candAna::play() {
  static int first(1);
  TH1D *h1(0);
  (void*)h1;
  vector<string> tnames;
  tnames.push_back("rf");
  tnames.push_back("os");
  tnames.push_back("c0");
  tnames.push_back("c1");
  string tname("nada");
  if (first) {
    first = 0;
    TDirectory *pDir = gDirectory;
    gFile->cd();
    for (unsigned int i = 0; i < tnames.size(); ++i) {
      h1 = new TH1D(Form("%s_m1", tnames[i].c_str()), "pt (muon1)", 50, 0., 50.);
      h1 = new TH1D(Form("%s_m2", tnames[i].c_str()), "pt (muon2)", 50, 0., 50.);
      h1 = new TH1D(Form("%s_e1", tnames[i].c_str()), "eta (muon1)", 50, -2.5, 2.5);
      h1 = new TH1D(Form("%s_e2", tnames[i].c_str()), "eta (muon2)", 50, -2.5, 2.5);
      h1 = new TH1D(Form("%s_f1", tnames[i].c_str()), "phi (muon1)", 50, -3.15, 3.15);
      h1 = new TH1D(Form("%s_f2", tnames[i].c_str()), "phi (muon2)", 50, -3.15, 3.15);
      h1 = new TH1D(Form("%s_dEta", tnames[i].c_str()), "dEta", 50, -2.5, 2.5);

      h1 = new TH1D(Form("%s_ptp", tnames[i].c_str()), "pt (pass)", 50, 0., 50.);
      h1 = new TH1D(Form("%s_pta", tnames[i].c_str()), "pt (all)", 50, 0., 50.);
      h1 = new TH1D(Form("%s_ptp1", tnames[i].c_str()), "pt muon1 (pass)", 50, 0., 50.);
      h1 = new TH1D(Form("%s_pta1", tnames[i].c_str()), "pt muon1 (all)",  50, 0., 50.);
      h1 = new TH1D(Form("%s_ptp2", tnames[i].c_str()), "pt muon2 (pass)", 50, 0., 50.);
      h1 = new TH1D(Form("%s_pta2", tnames[i].c_str()), "pt muon2 (all)",  50, 0., 50.);

      h1 = new TH1D(Form("%s_etap", tnames[i].c_str()), "eta (pass)", 50, -2.5, 2.5);
      h1 = new TH1D(Form("%s_etaa", tnames[i].c_str()), "eta (all)",  50, -2.5, 2.5);
      h1 = new TH1D(Form("%s_etap1", tnames[i].c_str()), "eta muon1 (pass)", 50, -2.5, 2.5);
      h1 = new TH1D(Form("%s_etaa1", tnames[i].c_str()), "eta muon1 (all)",  50, -2.5, 2.5);
      h1 = new TH1D(Form("%s_etap2", tnames[i].c_str()), "eta muon2 (pass)", 50, -2.5, 2.5);
      h1 = new TH1D(Form("%s_etaa2", tnames[i].c_str()), "eta muon2 (all)",  50, -2.5, 2.5);

      h1 = new TH1D(Form("%s_phip", tnames[i].c_str()), "phi (pass)", 50, -3.15, 3.15);
      h1 = new TH1D(Form("%s_phia", tnames[i].c_str()), "phi (all)",  50, -3.15, 3.15);
      h1 = new TH1D(Form("%s_phip1", tnames[i].c_str()), "phi muon1 (pass)", 50, -3.15, 3.15);
      h1 = new TH1D(Form("%s_phia1", tnames[i].c_str()), "phi muon1 (all)",  50, -3.15, 3.15);
      h1 = new TH1D(Form("%s_phip2", tnames[i].c_str()), "phi muon2 (pass)", 50, -3.15, 3.15);
      h1 = new TH1D(Form("%s_phia2", tnames[i].c_str()), "phi muon2 (all)",  50, -3.15, 3.15);
    }
    pDir->cd();
  }
  vector<string> refTriggerList;
  //  refTriggerList.push_back("HLT_Mu7p5_Track2_Jpsi_v2");
  refTriggerList.push_back("HLT_Dimuon16_Jpsi_v2");
  //NO has OS seed!  refTriggerList.push_back("HLT_Dimuon10_Jpsi_Barrel_v2");
  bool refTrigger(false), osTrigger(false), ctrl0Trigger(false), ctrl1Trigger(false);
  for (int i = 0; i < NHLT; ++i) {
    for (unsigned int j = 0; j < refTriggerList.size(); ++j) {
      if ((fpEvt->fHLTNames[i] == refTriggerList[j]) && (fpEvt->fHLTResult[i])) {
	//	cout << "event triggered by " << fpEvt->fHLTNames[i] << endl;
	refTrigger = true;
      }
    }
    if ((fpEvt->fHLTNames[i] == "HLT_DoubleMu4_3_Jpsi_Displaced_v2") && (fpEvt->fHLTResult[i])) {
      //      cout << "event triggered by " << fpEvt->fHLTNames[i] << endl;
      osTrigger = true;
    }
    if ((fpEvt->fHLTNames[i] == "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_v2") && (fpEvt->fHLTResult[i])) {
      //      cout << "event triggered by " << fpEvt->fHLTNames[i] << endl;
      ctrl0Trigger = true;
    }
    if ((fpEvt->fHLTNames[i] == "HLT_Dimuon0er16_Jpsi_NoVertexing_v2") && (fpEvt->fHLTResult[i])) {
      //      cout << "event triggered by " << fpEvt->fHLTNames[i] << endl;
      ctrl1Trigger = true;
    }

  }

  if (!refTrigger) return;

  TAnaCand *pC(0);
  int m1(-1), m2(-1);
  double m1Pt(0.), m2Pt(0.), m1Eta(-99.), m2Eta(-66.), dEta(99.), m1Phi(-66.), m2Phi(-66.);
  for (int i = 0; i < fpEvt->nCands(); ++i) {
    pC = fpEvt->getCand(i);
    if ((pC->fType == 300511) || (pC->fType == 300531) ||(pC->fType == 300521)) {
      vector<int> sigIdx;
      getSigTracks(sigIdx, pC);
      m1 = m2 = -1;
      for (unsigned int i = 0; i < sigIdx.size(); ++i) {
	// cout << "muon = " << fpEvt->getSimpleTrack(sigIdx[i])->getMuonID()
	//      << " " << Form(" %4d ", sigIdx[i])
	//      << " pT/eta/phi = "
	//      << fpEvt->getSimpleTrack(sigIdx[i])->getP().Perp() << " "
	//      << fpEvt->getSimpleTrack(sigIdx[i])->getP().Eta() << " "
	//      << fpEvt->getSimpleTrack(sigIdx[i])->getP().Phi() << " "
	//      << endl;
	if (fpEvt->getSimpleTrack(sigIdx[i])->getMuonID()) {
	  if (m1 == -1) {
	    m1 = sigIdx[i];
	    m1Pt  = fpEvt->getSimpleTrack(sigIdx[i])->getP().Perp();
	    m1Eta = fpEvt->getSimpleTrack(sigIdx[i])->getP().Eta();
	    m1Phi = fpEvt->getSimpleTrack(sigIdx[i])->getP().Phi();
	  } else {
	    m2 = sigIdx[i];
	    m2Pt  = fpEvt->getSimpleTrack(sigIdx[i])->getP().Perp();
	    m2Eta = fpEvt->getSimpleTrack(sigIdx[i])->getP().Eta();
	    m2Phi = fpEvt->getSimpleTrack(sigIdx[i])->getP().Phi();
	  }
	}
      }

      if (m1Pt < m2Pt) {
	double bla = m1Pt;
	m1Pt = m2Pt;
	m2Pt = bla;

	bla = m1Eta;
	m1Eta = m2Eta;
	m2Eta = bla;

	bla = m1Phi;
	m1Phi = m2Phi;
	m2Phi = bla;
      }

      dEta = TMath::Abs(m1Eta - m2Eta);
      ((TH1D*)gFile->Get("rf_dEta"))->Fill(dEta);
      ((TH1D*)gFile->Get("rf_m1"))->Fill(m1Pt);
      ((TH1D*)gFile->Get("rf_m2"))->Fill(m2Pt);
      ((TH1D*)gFile->Get("rf_e1"))->Fill(m1Eta);
      ((TH1D*)gFile->Get("rf_e2"))->Fill(m2Eta);
      ((TH1D*)gFile->Get("rf_f1"))->Fill(m1Phi);
      ((TH1D*)gFile->Get("rf_f2"))->Fill(m2Phi);

      // cout << "m1Pt = " << m1Pt << " m1Eta = " << m1Eta << " m2Pt = " << m2Pt << " m2Eta = " << m2Eta << " dEta = " << dEta << endl;

      double flsxy = pC->fVtx.fDxy/pC->fVtx.fDxyE;
      if (flsxy > 3.0 && (m1Pt > 4.0) && (m2Pt > 3.0) && (dEta < 1.8) && (TMath::Abs(m1Eta) < 1.6) && (TMath::Abs(m2Eta) < 1.6)) {
	tname = "os";
	((TH1D*)gFile->Get(Form("%s_pta", tname.c_str())))->Fill(pC->fPlab.Perp());
	((TH1D*)gFile->Get(Form("%s_pta1", tname.c_str())))->Fill(m1Pt);
	((TH1D*)gFile->Get(Form("%s_pta2", tname.c_str())))->Fill(m2Pt);
	((TH1D*)gFile->Get(Form("%s_etaa", tname.c_str())))->Fill(pC->fPlab.Eta());
	((TH1D*)gFile->Get(Form("%s_etaa1", tname.c_str())))->Fill(m1Eta);
	((TH1D*)gFile->Get(Form("%s_etaa2", tname.c_str())))->Fill(m2Eta);
	((TH1D*)gFile->Get(Form("%s_phia", tname.c_str())))->Fill(pC->fPlab.Phi());
	((TH1D*)gFile->Get(Form("%s_phia1", tname.c_str())))->Fill(m1Phi);
	((TH1D*)gFile->Get(Form("%s_phia2", tname.c_str())))->Fill(m2Phi);

	if (osTrigger) {
	  tname = "os";
	  ((TH1D*)gFile->Get(Form("%s_ptp", tname.c_str())))->Fill(pC->fPlab.Perp());
	  ((TH1D*)gFile->Get(Form("%s_ptp1", tname.c_str())))->Fill(m1Pt);
	  ((TH1D*)gFile->Get(Form("%s_ptp2", tname.c_str())))->Fill(m2Pt);
	  ((TH1D*)gFile->Get(Form("%s_etap", tname.c_str())))->Fill(pC->fPlab.Eta());
	  ((TH1D*)gFile->Get(Form("%s_etap1", tname.c_str())))->Fill(m1Eta);
	  ((TH1D*)gFile->Get(Form("%s_etap2", tname.c_str())))->Fill(m2Eta);
	  ((TH1D*)gFile->Get(Form("%s_phip", tname.c_str())))->Fill(pC->fPlab.Phi());
	  ((TH1D*)gFile->Get(Form("%s_phip1", tname.c_str())))->Fill(m1Phi);
	  ((TH1D*)gFile->Get(Form("%s_phip2", tname.c_str())))->Fill(m2Phi);
	}

	if (ctrl0Trigger) {
	  tname = "c0";
	  ((TH1D*)gFile->Get(Form("%s_ptp", tname.c_str())))->Fill(pC->fPlab.Perp());
	  ((TH1D*)gFile->Get(Form("%s_ptp1", tname.c_str())))->Fill(m1Pt);
	  ((TH1D*)gFile->Get(Form("%s_ptp2", tname.c_str())))->Fill(m2Pt);
	  ((TH1D*)gFile->Get(Form("%s_etap", tname.c_str())))->Fill(pC->fPlab.Eta());
	  ((TH1D*)gFile->Get(Form("%s_etap1", tname.c_str())))->Fill(m1Eta);
	  ((TH1D*)gFile->Get(Form("%s_etap2", tname.c_str())))->Fill(m2Eta);
	  ((TH1D*)gFile->Get(Form("%s_phip", tname.c_str())))->Fill(pC->fPlab.Phi());
	  ((TH1D*)gFile->Get(Form("%s_phip1", tname.c_str())))->Fill(m1Phi);
	  ((TH1D*)gFile->Get(Form("%s_phip2", tname.c_str())))->Fill(m2Phi);
	}

	if (ctrl1Trigger) {
	  tname = "c1";
	  ((TH1D*)gFile->Get(Form("%s_ptp", tname.c_str())))->Fill(pC->fPlab.Perp());
	  ((TH1D*)gFile->Get(Form("%s_ptp1", tname.c_str())))->Fill(m1Pt);
	  ((TH1D*)gFile->Get(Form("%s_ptp2", tname.c_str())))->Fill(m2Pt);
	  ((TH1D*)gFile->Get(Form("%s_etap", tname.c_str())))->Fill(pC->fPlab.Eta());
	  ((TH1D*)gFile->Get(Form("%s_etap1", tname.c_str())))->Fill(m1Eta);
	  ((TH1D*)gFile->Get(Form("%s_etap2", tname.c_str())))->Fill(m2Eta);
	  ((TH1D*)gFile->Get(Form("%s_phip", tname.c_str())))->Fill(pC->fPlab.Phi());
	  ((TH1D*)gFile->Get(Form("%s_phip1", tname.c_str())))->Fill(m1Phi);
	  ((TH1D*)gFile->Get(Form("%s_phip2", tname.c_str())))->Fill(m2Phi);
	}


	break; // fill only for one candidate with fulfilling J/psi muons and displacement
      }
    }
  }


  // // -- get list of indices of tracks making up candidate

}

void candAna::play2() {

  const int PRINT = 1;
  //const float dRMin = 0.02;
  //int muTightHlt = 0, muHlt = 0, muPairHlt=0, muPairTightHlt=0;
  static int cand=0;
  cand++;

  if (PRINT>0) cout<<" PLAY2: Candidate "<<cand<<endl<<" Muons "<<fpEvt->nMuons()<<endl;
  for (int it = 0; it< fpEvt->nMuons(); ++it) { // loop over muons
    TAnaMuon *muon = fpEvt->getMuon(it); // get muon






    // some direct muon methods
    //if (PRINT>0) muon->dump(); // big printout

    // some direct muon methods
    //int muonId = muon->fMuID; // muon ID bits
    //int muonIdx = muon->fMuIndex;  // index in the muon list = it
    //int genidx = muon->fGenIndex; // gen index

    // get track index
    int itrk = muon->fIndex; // index of the SimpleTrack
    TVector3 muonMom = muon->fPlab;
    //double ptMuon  = muonMom.Perp();

    if (itrk<0) continue; // skip muons without inner tracker info
    bool muid  = tightMuon(muon);  // tight muons

    if (PRINT) cout<<it<<" "<<itrk<<" "<<muid<<" "<<muonMom.Perp()<<" "
		  <<muonMom.Eta()<<" "<<muonMom.Phi()<<endl;

    //if (pM->fMuID & 1) {

    // // now check all pairs
    // for (int it2 = (it+1); it2< fpEvt->nMuons(); ++it2) { // loop over muons
    //   TAnaMuon *muon2 = fpEvt->getMuon(it2); // get muon
    //   if (muon2->fIndex<0) continue; // skip muons without inner tracker info
    //   bool matched = doTriggerMatching(muon,muon2); // see if a pair matches HLT dimuon object
    //   if (matched) {
    // 	muPairHlt++;
    // 	bool muid2  = tightMuon(muon2);  // tight muons
    // 	if (muid&&muid2) muPairTightHlt++;
    //   } // if matched
    // } // 2ns muons


  } // END MUON LOOP

  //if (PRINT) cout<<" HLT matched muons in this event "<<muHlt<<" "<<muTightHlt<<" "<<muPairHlt<<" "<<muPairTightHlt<<endl;


  // ((TH1D*)fHistDir->Get("test5"))->Fill(muHlt);
  // ((TH1D*)fHistDir->Get("test4"))->Fill(muTightHlt);

  int numSTracks = fpEvt->nSimpleTracks();
  if (PRINT>0) cout<<" Number of simple tracks "<<numSTracks<<endl;
  if (PRINT>9) {
    TSimpleTrack *pT(0);
    for (int i = 0; i < numSTracks; ++i) {
      pT = fpEvt->getSimpleTrack(i);
      double pt  = pT->getP().Perp();
      double eta = pT->getP().Eta();
      double phi = pT->getP().Phi();
      cout<<i<<" "<<pT->getIndex()<<" "<<pt<<" "<<eta<<" "<<phi<<endl;
    }
  }

  cout << "Found " << fpEvt->nCands() << " cands in event" << endl;
  cout << "Found " << fpEvt->nGenCands() << " gen cands in event" << endl;
  cout << "Found " << fpEvt->nSigTracks() << " sig tracks in event" << endl;
  cout << "Found " << fpEvt->nRecTracks() << " rec tracks in event" << endl;

  // cout << " Tracks ------------------------------" << endl;
  // TAnaTrack *pTrack;
  // for (int it = 0; it < fpEvt->nRecTracks(); ++it) {
  //   pTrack = fpEvt->getRecTrack(it);
  //   //double pt = pTrack->fPlab.Perp();
  //   cout << "R: "; pTrack->dump();
  // }


  if (PRINT>8) {
    cout << " Candidates ------------------------------" <<endl;
    TAnaCand *pCand;
    for (int it = 0; it < fpEvt->nCands(); ++it) {
      pCand = fpEvt->getCand(it);
      cout << "C: " << pCand->fType << " "; pCand->dump();
      //int s1=pCand->fSig1;
      //int s2=pCand->fSig2;
      //int i1 =fpEvt->getSigTrack(pCand->fSig1)->fIndex;
      //int i2 =fpEvt->getSigTrack(pCand->fSig2)->fIndex;
      //cout<<s1<<"-"<<s2<<" "<<i1<<"-"<<i2<<endl;
      TAnaTrack *p0;
      // look at signal tracks
      for (int it = pCand->fSig1; it <= pCand->fSig2; ++it) {
	p0 = fpEvt->getSigTrack(it);
	//int i = p0->fIndex;
	p0->dump();
	//cout<<it<<" "<<i<<" "<<p0->fGenIndex<<" "<<p0->fPlab.Pt()<<endl;
      }

      // another way
      if (0) {
	vector<int> cIdx;
	TSimpleTrack *ps;
	getSigTracks(cIdx, pCand);  // get vector if SigTrack indices
	for (unsigned int i = 0; i < cIdx.size(); ++i) {
	  int idx = cIdx[i];
	  ps = fpEvt->getSimpleTrack(idx); // get simple track with index i
	  cout << " track idx = " << i<<" "<<ps->getIndex()
	    //<< " with ID = " << fpEvt->getSimpleTrackMCID(ps->getIndex())
	       <<" "<<cIdx[i]<<" "<<(ps->getP()).Perp()
	       <<  endl;
	  //ps->dump();
	}
      } // if
    }
  }

  cout<<" Num PVs "<<fpEvt->nPV()<<" : ";
  for (int i = 0; i < fpEvt->nPV(); ++i) {
    double z = fpEvt->getPV(i)->fPoint.Z();
    cout<<z<<" ";
  }
  cout<<endl;

  // Look at Trig Object v2
  TTrgObjv2 *pTO;
   cout<<" Dump TTrgObjv2 "<<fpEvt->nTrgObjv2()<<endl;
    for (int i = 0; i < fpEvt->nTrgObjv2(); ++i) {
      pTO = fpEvt->getTrgObjv2(i);
      //pTO->dump();

      cout<<i<<" hlt "<<pTO->fHltPath<<" hlt-index "<<pTO->fHltIndex<<" module label "
	  <<pTO->fLabel<<" type "<<pTO->fType<<" number "<<pTO->fNumber<<endl;

      vector<int> muonIndex = pTO->fIndex;
      vector<int> muonID = pTO->fID;
      vector<TLorentzVector> muonP = pTO->fP;
      int num = muonIndex.size();
      for(int n=0;n<num;++n) {
	int index = muonIndex[n];
	int id = muonID[n];
	TLorentzVector p = muonP[n];
	cout<<n<<" index "<<index<<" id "<<id<<" pt/eta/phi "<<p.Pt()<<" "<<p.Eta()<<" "<<p.Phi()<<endl;
      }

    }

}


// ----------------------------------------------------------------------
void candAna::fillMuonData(muonData &a, TAnaMuon *pt) {

  if (0 == pt) {

    a.mbdt          = -99.;
    a.pt            = -99.;
    a.eta           = -99.;
    a.validMuonHits    = -99.;
    a.glbNChi2         = -99.;
    a.nMatchedStations = -99;
    a.validPixelHits   = -99;
    a.validPixelHits2  = -99;
    a.trkLayerWithHits = -99.;

    a.trkValidFract = -99.;
    a.segComp       = -99.;
    a.chi2LocMom    = -99.;
    a.chi2LocPos    = -99.;
    a.glbTrackProb  = -99.;
    a.NTrkVHits     = -99.;
    a.NTrkEHitsOut  = -99.;

    a.kink          = -99.;

    a.dpt           = -99.;
    a.dptrel        = -99.;
    a.deta          = -99.;
    a.dphi          = -99.;
    a.dr            = -99.;

  } else {

    a.mbdt          = -99.;
    a.pt            = pt->fPlab.Perp();
    a.eta           = pt->fPlab.Eta();
    a.validMuonHits    = pt->fNvalidMuonHits;
    a.glbNChi2         = pt->fGtrkNormChi2;
    a.nMatchedStations = pt->fNmatchedStations;
    a.validPixelHits   = fpReader->numberOfPixLayers(pt);
    a.validPixelHits2  = fpReader->numberOfPixelHits(pt);
    a.trkLayerWithHits = fpReader->numberOfTrackerLayers(pt);

    a.trkValidFract = pt->fItrkValidFraction;
    a.segComp       = pt->fSegmentComp;
    a.chi2LocMom    = pt->fChi2LocalMomentum;
    a.chi2LocPos    = pt->fChi2LocalPosition;
    a.glbTrackProb  = pt->fGtrkProb;
    a.NTrkVHits     = static_cast<float>(pt->fNumberOfValidTrkHits);
    a.NTrkEHitsOut  = static_cast<float>(pt->fNumberOfLostTrkHits);

    a.kink          = pt->fMuonChi2;

    a.dpt           = pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag();
    a.dptrel        = TMath::Abs(pt->fInnerPlab.Mag() - pt->fOuterPlab.Mag())/pt->fInnerPlab.Mag();
    if (pt->fOuterPlab.Mag() > 3.) {
      a.deta          = pt->fInnerPlab.Eta() - pt->fOuterPlab.Eta();
      a.dphi          = pt->fInnerPlab.DeltaPhi(pt->fOuterPlab);
      a.dr            = pt->fInnerPlab.DeltaR(pt->fOuterPlab);
    } else {
      a.deta          = -99.;
      a.dphi          = -99.;
      a.dr            = -99.;
    }

  }
}


// ----------------------------------------------------------------------
void candAna::print1() {

  cout << "--------------------------------------------------" << endl;
  cout << "event : " << fEvt << " " << fEvent << endl;
  cout << "muon 1: " << fpMuon1->fRefPlab.X() << " " << fpMuon1->fRefPlab.Y() << " " << fpMuon1->fRefPlab.Z() << endl;
  cout << "muon 2: " << fpMuon2->fRefPlab.X() << " " << fpMuon2->fRefPlab.Y() << " " << fpMuon2->fRefPlab.Z() << endl;
  int pvidx = (fpCand->fPvIdx > -1? fpCand->fPvIdx : 0);
  cout << "PV: x = " << fpEvt->getPV(pvidx)->fPoint.X()
       << " y = " << fpEvt->getPV(pvidx)->fPoint.Y()
       << " z = " << fpEvt->getPV(pvidx)->fPoint.Z()
       << " c = " << fpEvt->getPV(pvidx)->fChi2 << "/" << fpEvt->getPV(pvidx)->fNdof
       << endl;
  TAnaVertex sv = fpCand->fVtx;
  cout << "SV: x = " << sv.fPoint.X() << " y = " << sv.fPoint.Y() << " z = " << sv.fPoint.Z()
       << " c = " << sv.fChi2 << "/" << sv.fNdof
       << endl;
  cout << "fl3d = " << fCandFL3d << " fl3de = " << fCandFL3dE << " fls3d = " << fCandFLS3d << endl;
  cout << "pvip = " << fCandPvIp3D << " pvipe = " << fCandPvIpE3D << " pvips = " << fCandPvIpS3D << endl;
  cout << "alpha = " << fCandA << " cosalpha = " <<  fCandCosA << endl;
}


// ----------------------------------------------------------------------
void candAna::pvStudy(bool bookHist) {
  static int   chan;
  static int   idx1, idx2, idx3;
  static float gx, gy, gz;
  static float sx, sy, sz;
  static float p1x, p1y, p1z, p1d;
  static float p2x, p2y, p2z, p2d;
  static float p3x, p3y, p3z, p3d;
  static float pt, eta, phi, m;
  static float lz1, lz2;
  static float d1, d2, d3; // distance between reco PV and gen PV
  static float dsv; // distance between reco SV and gen SV
  static float a1, a2, a3; // pointing angle between
  static float gfl, fl1, fl2, fl3; // flight length
  static float fl3d, fls3d; // flight length significance from the real candidate
  static float gt, t1, t2, t3; // (3D) lifetime
  static float gs, s1, s2, s3; // (2D) lifetime
  static float mult1, mult2; // PV track multiplicity
  static float stat1, stat2; // PV status
  static float chi1, chi2; // PV chi2
  static float prob1, prob2; // PV prob
  static float dz12, dzmin;
  static float npv;
  if (bookHist) {
    TDirectory *pDir = gDirectory;
    fHistDir->cd();

    fPvStudyTree = new TTree("pvstudy", "pvstudy");
    fPvStudyTree->Branch("pt",   &pt,  "pt/F");
    fPvStudyTree->Branch("eta",  &eta, "eta/F");
    fPvStudyTree->Branch("phi",  &phi, "phi/F");
    fPvStudyTree->Branch("m",    &m, "m/F");
    fPvStudyTree->Branch("chan", &chan, "chan/I");
    fPvStudyTree->Branch("npv",  &npv, "npv/F");
    fPvStudyTree->Branch("idx1", &idx1, "idx1/I");
    fPvStudyTree->Branch("idx2", &idx2, "idx2/I");
    fPvStudyTree->Branch("idx3", &idx3, "idx3/I");

    fPvStudyTree->Branch("lz1",  &lz1, "lz1/F");
    fPvStudyTree->Branch("lz2",  &lz2, "lz2/F");
    fPvStudyTree->Branch("mult1",&mult1, "mult1/F");
    fPvStudyTree->Branch("mult2",&mult2, "mult2/F");
    fPvStudyTree->Branch("prob1",&prob1, "prob1/F");
    fPvStudyTree->Branch("prob2",&prob2, "prob2/F");
    fPvStudyTree->Branch("chi1", &chi1, "chi1/F");
    fPvStudyTree->Branch("chi2", &chi2, "chi2/F");
    fPvStudyTree->Branch("dz12", &dz12,  "dz12/F");
    fPvStudyTree->Branch("dzmin",&dzmin,  "dzmin/F");

    fPvStudyTree->Branch("gfl",   &gfl,   "gfl/F");
    fPvStudyTree->Branch("gt",    &gt,    "gt/F");
    fPvStudyTree->Branch("fls3d", &fls3d, "fls3d/F");
    fPvStudyTree->Branch("fl3d",  &fl3d,  "fl3d/F");
    fPvStudyTree->Branch("fl1",   &fl1,   "fl1/F");
    fPvStudyTree->Branch("fl2",   &fl2,   "fl2/F");
    fPvStudyTree->Branch("fl3",   &fl3,   "fl3/F");
    fPvStudyTree->Branch("t1",    &t1,    "t1/F");
    fPvStudyTree->Branch("t2",    &t2,    "t2/F");
    fPvStudyTree->Branch("t3",    &t3,    "t3/F");

    // -- 2D version
    fPvStudyTree->Branch("gs",    &gs,    "gs/F");
    fPvStudyTree->Branch("s1",    &s1,    "s1/F");
    fPvStudyTree->Branch("s2",    &s2,    "s2/F");
    fPvStudyTree->Branch("s3",    &s3,    "s3/F");

    fPvStudyTree->Branch("gx",  &gx, "gx/F");
    fPvStudyTree->Branch("gy",  &gy, "gy/F");
    fPvStudyTree->Branch("gz",  &gz, "gz/F");

    fPvStudyTree->Branch("p1x",  &p1x, "p1x/F");
    fPvStudyTree->Branch("p1y",  &p1y, "p1y/F");
    fPvStudyTree->Branch("p1z",  &p1z, "p1z/F");
    fPvStudyTree->Branch("p1d",  &p1d, "p1d/F");

    fPvStudyTree->Branch("p2x",  &p2x, "p2x/F");
    fPvStudyTree->Branch("p2y",  &p2y, "p2y/F");
    fPvStudyTree->Branch("p2z",  &p2z, "p2z/F");
    fPvStudyTree->Branch("p2d",  &p2d, "p2d/F");

    fPvStudyTree->Branch("p3x",  &p3x, "p3x/F");
    fPvStudyTree->Branch("p3y",  &p3y, "p3y/F");
    fPvStudyTree->Branch("p3z",  &p3z, "p3z/F");

    fPvStudyTree->Branch("sx",  &sx,  "sx/F");
    fPvStudyTree->Branch("sy",  &sy,  "sy/F");
    fPvStudyTree->Branch("sz",  &sz,  "sz/F");
    fPvStudyTree->Branch("dsv", &dsv, "dsv/F");

    fPvStudyTree->Branch("d1",  &d1, "d1/F");
    fPvStudyTree->Branch("a1",  &a1, "a1/F");
    fPvStudyTree->Branch("d2",  &d2, "d2/F");
    fPvStudyTree->Branch("a2",  &a2, "a2/F");
    fPvStudyTree->Branch("d3",  &d3, "d3/F");
    fPvStudyTree->Branch("a3",  &a3, "a3/F");
    pDir->cd();
    return;
  }

  // -- this only works for truth-matched cands (no other cand available at this point)
  if (fGenBTmi < 0) return;
  if (fCandTmi < 0) return;
  TAnaCand *pCand = fpEvt->getCand(fCandTmi);
  if (0 == pCand) return;
  if (pCand->fPv2Idx < 0) return;

  fls3d = pCand->fVtx.fD3d/pCand->fVtx.fD3dE;
  fl3d  = pCand->fVtx.fD3d;

  TGenCand *pB = fpEvt->getGenCand(fGenBTmi);
  TGenCand *pM1 = fpEvt->getGenCand(fGenM1Tmi);
  TGenCand *pM2 = fpEvt->getGenCand(fGenM2Tmi);
  chan = detChan(pM1->fP.Eta(), pM2->fP.Eta());

  TVector3 sv = pCand->fVtx.fPoint;
  TVector3 fl, fl2d;
  TVector3 genPV, genSV, distPV(99., 99., 99.), distPV2(99., 99., 99.),  distPVL(99., 99., 99.);

  genSV = TVector3(fpEvt->getGenCand(fGenM1Tmi)->fV.X(), fpEvt->getGenCand(fGenM1Tmi)->fV.y(), fpEvt->getGenCand(fGenM1Tmi)->fV.Z());
  sx = genSV.X();
  sy = genSV.Y();
  sz = genSV.Z();
  dsv = (genSV - sv).Mag();

  double gm = pB->fP.Mag();
  double gp = pB->fP.P();
  // mother pointer
  TGenCand *pM = fpEvt->getGenTWithIndex(pB->fMom1);
  // use the mother if it has the same PDGID (it oscillated)
  if (TMath::Abs(pB->fID) != TMath::Abs(pM->fID)) pM = pB;
  double x = (pM1->fV - pM->fV).Mag();
  TVector3 v2d = pM1->fV - pM->fV;
  v2d.SetZ(0.);
  double x2d = v2d.Mag();
  double gp2d = pB->fP.Perp();

  genPV = TVector3(pM->fV.X(), pM->fV.y(), pM->fV.Z());
  gx = genPV.X();
  gy = genPV.Y();
  gz = genPV.Z();

  if ((pM->fV - genPV).Mag() > 0.0001) {
    cout << "pM->fV = (" << pM->fV.X() << "," << pM->fV.Y() << "," << pM->fV.Z() << " gen PV = (" <<  genPV.X() << "," << genPV.Y() << "," << genPV.Z() << ")" << endl;
    fpEvt->dumpGenBlock();
  }

  gt = x*gm/gp/TMath::Ccgs();
  gs = x2d*gm/gp2d/TMath::Ccgs();
  gfl = x;

  p1x = fpEvt->getPV(pCand->fPvIdx)->fPoint.X();
  p1y = fpEvt->getPV(pCand->fPvIdx)->fPoint.Y();
  p1z = fpEvt->getPV(pCand->fPvIdx)->fPoint.Z();
  lz1   = pCand->fPvLip;
  mult1 = fpEvt->getPV(pCand->fPvIdx)->fNtracks;
  chi1  = fpEvt->getPV(pCand->fPvIdx)->fChi2;
  prob1 = fpEvt->getPV(pCand->fPvIdx)->fProb;

  p2x = fpEvt->getPV(pCand->fPv2Idx)->fPoint.X();
  p2y = fpEvt->getPV(pCand->fPv2Idx)->fPoint.Y();
  p2z = fpEvt->getPV(pCand->fPv2Idx)->fPoint.Z();
  lz2 = pCand->fPv2Lip;
  mult2 = fpEvt->getPV(pCand->fPv2Idx)->fNtracks;
  chi2  = fpEvt->getPV(pCand->fPv2Idx)->fChi2;
  prob2 = fpEvt->getPV(pCand->fPv2Idx)->fProb;

  dz12  = p2z-p1z;
  dzmin = 99.;
  for (int ipv = 0; ipv < fpEvt->nPV(); ++ipv) {
    if (ipv == pCand->fPvIdx) continue;
    double delta = fpEvt->getPV(ipv)->fPoint.Z() - p1z;
    if (TMath::Abs(delta) < TMath::Abs(dzmin)) {
      dzmin = delta;
    }
  }
  distPV = fpEvt->getPV(pCand->fPvIdx)->fPoint - genPV;
  distPV2 = fpEvt->getPV(pCand->fPv2Idx)->fPoint - genPV;

  d1 = distPV.Mag();
  d2 = distPV2.Mag();

  TVector3 plab   = pCand->fPlab;
  m =  pCand->fMass;
  pt = plab.Perp();
  eta = plab.Eta();
  phi = plab.Phi();

  const double massOverC = m/TMath::Ccgs();

  fl = sv - fpEvt->getPV(pCand->fPvIdx)->fPoint;
  fl1 = fl.Mag();
  a1 = TMath::ACos(plab.Dot(fl) / (plab.Mag() * fl.Mag()));
  t1 = fl1 * TMath::Cos(a1) / plab.Mag() * massOverC;

  TVector3 plab2d = plab;
  plab2d.SetZ(0.);
  fl2d = fl;
  fl2d.SetZ(0.);
  double fl12d = fl2d.Mag();
  double a12d = TMath::ACos(plab2d.Dot(fl2d) / (plab2d.Mag() * fl2d.Mag()));
  s1 = fl12d * TMath::Cos(a12d) / plab2d.Mag() * massOverC;


  fl = sv - fpEvt->getPV(pCand->fPv2Idx)->fPoint;
  fl2 = fl.Mag();
  a2 = TMath::ACos(plab.Dot(fl) / (plab.Mag() * fl.Mag()));
  t2 = fl2 * TMath::Cos(a2) / plab.Mag() * massOverC;

  fl2d = fl;
  fl2d.SetZ(0.);
  double fl22d = fl2d.Mag();
  double a22d = TMath::ACos(plab2d.Dot(fl2d) / (plab2d.Mag() * fl2d.Mag()));
  s2 = fl22d * TMath::Cos(a22d) / plab2d.Mag() * massOverC;


  // -- loop over all PV in event and try out other approaches
  int minAlphaIdx(-1);
  double minAlpha(99.), alphaL;
  int minPvIpIdx(-1);
  double minPvIp(99.), pvipL;
  npv = fpEvt->nPV();

  for (int ipv = 0; ipv < fpEvt->nPV(); ++ipv) {
    fl = sv - fpEvt->getPV(ipv)->fPoint;
    distPVL = fpEvt->getPV(ipv)->fPoint - genPV;
    // -- best pointing angle
    alphaL = TMath::ACos(plab.Dot(fl) / (plab.Mag() * fl.Mag()));
    if (alphaL < minAlpha) {
      minAlphaIdx = ipv;
      minAlpha = alphaL;
    }
  }

  distPVL = fpEvt->getPV(minAlphaIdx)->fPoint - genPV;
  d3 = distPVL.Mag();
  fl = sv - fpEvt->getPV(minAlphaIdx)->fPoint;
  fl3 = fl.Mag();
  a3 = TMath::ACos(plab.Dot(fl) / (plab.Mag() * fl.Mag()));
  t3 = fl3 * TMath::Cos(a3) / plab.Mag() * massOverC;

  fl2d = fl;
  fl2d.SetZ(0.);
  double fl32d = fl2d.Mag();
  double a32d = TMath::ACos(plab2d.Dot(fl2d) / (plab2d.Mag() * fl2d.Mag()));
  s3 = fl32d * TMath::Cos(a32d) / plab2d.Mag() * massOverC;


  p3x = fpEvt->getPV(minAlphaIdx)->fPoint.X();
  p3y = fpEvt->getPV(minAlphaIdx)->fPoint.Y();
  p3z = fpEvt->getPV(minAlphaIdx)->fPoint.Z();

  idx1 = pCand->fPvIdx;
  idx2 = pCand->fPv2Idx;
  idx3 = minAlphaIdx;

  // -- fill delta(z) to closest other PV
  p1d = p2d = p3d = 99.;
  for (int ipv = 0; ipv < fpEvt->nPV(); ++ipv) {
    if (ipv != pCand->fPvIdx) {
      double dz = TMath::Abs(fpEvt->getPV(ipv)->fPoint.Z() - fpEvt->getPV(pCand->fPvIdx)->fPoint.Z());
      if (dz < p1d) p1d = dz;
    }

    if (ipv != pCand->fPv2Idx) {
      double dz = TMath::Abs(fpEvt->getPV(ipv)->fPoint.Z() - fpEvt->getPV(pCand->fPv2Idx)->fPoint.Z());
      if (dz < p2d) p2d = dz;
    }

    if (ipv != pCand->fPv2Idx) {
      double dz = TMath::Abs(fpEvt->getPV(ipv)->fPoint.Z() - fpEvt->getPV(minAlphaIdx)->fPoint.Z());
      if (dz < p3d) p3d = dz;
    }
  }

  fPvStudyTree->Fill();

}


// ----------------------------------------------------------------------
double candAna::getDetVarComb(TAnaMuon *mu) {
  double combination(0);
  if ((mu->fvDThits.size() < 4)
      || (mu->fvRPChits.size() < 4)
      || (mu->fvCSChits.size() < 4)) {
    // cout << "TAnaMuon without enough entries in hit vectors: "
    // 	 << mu->fvDThits.size() << " / "
    //   	 << mu->fvRPChits.size() << " / "
    //   	 << mu->fvCSChits.size() << " / "
    // 	 << endl;
    return 0.;
  }
  for (unsigned int station = 0; station < 4; ++station) {
    combination += ((mu->fvDThits)[station])/2.;
    combination += (mu->fvRPChits)[station];
    if ((mu->fvCSChits)[station] > 6) {
      combination += 6;
    } else {
      combination += (mu->fvCSChits)[station];
    }
  }
  return combination;
}
