#include "plotResults.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>


#include "TROOT.h"
#include "TStyle.h"
#include "TKey.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TF1.h"
#include "TFitResult.h"

#include "common/dataset.hh"
#include "common/util.hh"
#include "common/Lumi.hh"

ClassImp(plotResults)

using namespace std;

// ----------------------------------------------------------------------
plotResults::plotResults(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotResults::loadFiles(files);

  changeSetup(dir, "plotResults", setup);
  init();

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fMassLo = 4.5;
  fMassHi = 6.5;

  fNoLo = 5.10;
  fNoHi = 5.40;

  fCsLo = 5.27;
  fCsHi = 5.47;

  fBgLo = 4.9;
  fBgHi = 5.9;

  fSgLo = 5.20;
  fSgHi = 5.45;


  fChan = 0;

  TH2D *h2(0);
  TH1D *h(0);
  for (int i = 0; i < fNchan; ++i) {
    h2 = new TH2D(Form("hAccAll%d", i), Form("hAccAll%d", i), 25, 0., 2.5, 25, 0., 50.);
    fhAccAll.push_back(h2);
    h2 = new TH2D(Form("hAccPass%d", i), Form("hAccPass%d", i), 25, 0., 2.5, 25, 0., 50.);
    fhAccPass.push_back(h2);

    h = new TH1D(Form("hAccPtAll%d", i), Form("hAccPtAll%d", i), 25, 0., 50.);
    fhAccPtAll.push_back(h);
    h = new TH1D(Form("hAccPtPass%d", i), Form("hAccPtPass%d", i), 25, 0., 50.);
    fhAccPtPass.push_back(h);

    h = new TH1D(Form("hAccEtaAll%d", i), Form("hAccEtaAll%d", i), 25, 0., 2.5);
    fhAccEtaAll.push_back(h);
    h = new TH1D(Form("hAccEtaPass%d", i), Form("hAccEtaPass%d", i), 25, 0., 2.5);
    fhAccEtaPass.push_back(h);

    h = new TH1D(Form("hGenAndAccNumbers%d", i), Form("hGenAndAccNumbers%d", i), 100, 0., 100.);
    fhGenAndAccNumbers.push_back(h);

    h = new TH1D(Form("hMassAbsNoCuts%d", i), Form("hMassAbsNoCuts%d", i), 100, 0, 10);
    fhMassAbsNoCuts.push_back(h);

    h = new TH1D(Form("hMassNoCuts%d", i), Form("hMassNoCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassNoCuts.push_back(h);

    h = new TH1D(Form("fhMassWithAnaCuts%d", i), Form("hMassChan%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAnaCuts.push_back(h);

    h = new TH1D(Form("hMassWithMuonCuts%d", i), Form("hMassWithMuonCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMuonCuts.push_back(h);

    h = new TH1D(Form("hMassWithTriggerCuts%d", i), Form("hMassWithTriggerCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithTriggerCuts.push_back(h);

    h = new TH1D(Form("hMassWithAllCuts%d", i), Form("hMassWithAllCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCuts.push_back(h);

    h = new TH1D(Form("hMassWithAllCutsBlind%d", i), Form("hMassWithAllCutsBlind%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithAllCutsBlind.push_back(h);

    h = new TH1D(Form("hW8MassWithAllCuts%d", i), Form("hW8MassWithAllCuts%d", i), NBINS, fMassLo, fMassHi);
    fhW8MassWithAllCuts.push_back(h);

    h = new TH1D(Form("hW8MassWithAllCutsBlind%d", i), Form("hW8MassWithAllCutsBlind%d", i), NBINS, fMassLo, fMassHi);
    fhW8MassWithAllCutsBlind.push_back(h);

    h = new TH1D(Form("hMassWithMassCuts%d", i), Form("hW8MassWithMassCuts%d", i), NBINS, fMassLo, fMassHi);
    fhMassWithMassCuts.push_back(h);

    h = new TH1D(Form("hNo%d", i), Form("hNo%d", i), 100, 4.9, 5.9);
    fhNo.push_back(h);

    h = new TH1D(Form("hNoC%d", i), Form("hNoC%d", i), 200, 4.9, 5.9);
    fhNoC.push_back(h);

    h = new TH1D(Form("hCs%d", i), Form("hCs%d", i), 100, 4.9, 5.9);
    fhCs.push_back(h);

    h = new TH1D(Form("hCsC%d", i), Form("hCsC%d", i), 200, 4.9, 5.9);
    fhCsC.push_back(h);

    h = new TH1D(Form("hB0%d", i), Form("hB0%d", i), 100, 4.9, 5.9);
    fhB0.push_back(h);

    h = new TH1D(Form("hB0C%d", i), Form("hB0C%d", i), 200, 4.9, 5.9);
    fhB0C.push_back(h);

  }

  if (1) {
    int year(fYear);
    year = 2012;
    fptFakePosKaons     = new PidTable(Form("../common/pidtables/%d-kaonPosFakeRate-mvaMuon.dat", year));
    fptFakeNegKaons     = new PidTable(Form("../common/pidtables/%d-kaonNegFakeRate-mvaMuon.dat", year));

    fptFakePosPions     = new PidTable(Form("../common/pidtables/%d-pionPosFakeRate-mvaMuon.dat", year));
    fptFakeNegPions     = new PidTable(Form("../common/pidtables/%d-pionNegFakeRate-mvaMuon.dat", year));

    fptFakePosProtons   = new PidTable(Form("../common/pidtables/%d-protonPosFakeRate-mvaMuon.dat", year));
    fptFakeNegProtons   = new PidTable(Form("../common/pidtables/%d-protonNegFakeRate-mvaMuon.dat", year));
  }

}


// ----------------------------------------------------------------------
plotResults::~plotResults() {

}


// ----------------------------------------------------------------------
void plotResults::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotResults::makeAll(string what) {
  if (what == "all" || what == "dumpdatasets") {
    dumpDatasets();
  }

  if (what == "all" || what == "genvalidation") {
    genSummary("bdmmMcOffAcc", "candAnaMuMu");
    genSummary("bdmmMcOff", "candAnaMuMu");
    genSummary("bdmmMc", "candAnaMuMu");
    genSummary("bsmmMcOffAcc", "candAnaMuMu");
    genSummary("bsmmMcOff", "candAnaMuMu");
    genSummary("bsmmMc", "candAnaMuMu");
    genSummary("bsmm0Mc", "candAnaMuMu");
    genSummary("bsmm2Mc", "candAnaMuMu");
    genSummary("bsmm3Mc", "candAnaMuMu");
    genSummary("bsmm4Mc", "candAnaMuMu");
    genSummary("bsmm5Mc", "candAnaMuMu");

    genSummary("bupsikMcOffAcc", "candAnaBu2JpsiK");
    genSummary("bupsikMcOff", "candAnaBu2JpsiK");
    genSummary("bupsikMc", "candAnaBu2JpsiK");
    genSummary("bspsiphiMcOffAcc", "candAnaBs2JpsiPhi");
    genSummary("bspsiphiMcOff", "candAnaBs2JpsiPhi");
    genSummary("bspsiphiMc", "candAnaBs2JpsiPhi");
    genSummary("bdpsikstarMcAcc", "candAnaBd2JpsiKstar");
    genSummary("bdpsikstarMc", "candAnaBd2JpsiKstar");


    // -- loop over all (effective) two-body backgrounds
    for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
      if (string::npos == it->first.find("Bg")) continue;
      cout << "===> Create genSummary for " << it->first << endl;
      genSummary(it->first, "candAnaMuMu");
    }

  }

  // -- this will recreate fHistFile!
  if (what == "dbx") {
    fillAndSaveHistograms(100000);
  }



}



// ----------------------------------------------------------------------
void plotResults::bookHist(string dsname) {


}


// ----------------------------------------------------------------------
void plotResults::dumpDatasets() {
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << formatTex(fBfPsiMuMu, Form("%s:BfPsiMuMu:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(fBfPsiMuMuE, Form("%s:BfPsiMuMu:err", fSuffix.c_str()), 5) << endl;

  fTEX << formatTex(fBfPhiKpKm, Form("%s:BfPhiKpKm:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(fBfPhiKpKmE, Form("%s:BfPhiKpKm:err", fSuffix.c_str()), 5) << endl;

  fTEX << formatTex(fBfKstarKpPim, Form("%s:BfKstarKpPim:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(fBfKstarKpPim, Form("%s:BfKstarKpPim:err", fSuffix.c_str()), 5) << endl;

  fTEX << formatTexErrSci(fCrossSection, 0., Form("%s:PythiaCrossSection:val", fSuffix.c_str()), -1) << endl;

  fTEX << "% ----------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    TH1D *h1 = it->second->getHist("monEvents", false);
    int nEvtFile = static_cast<int>(h1->GetBinContent(1));
    int nCands   = static_cast<int>(h1->GetBinContent(2));
    double epscand  = static_cast<double>(nCands)/nEvtFile;
    double epscandE = dEff(nCands, nEvtFile);
    fTEX << Form("\\vdef{%s:%s:name} {%s}", fSuffix.c_str(), it->first.c_str(), it->first.c_str()) << endl;
    fTEX << Form("\\vdef{%s:%s:decay} {%s}", fSuffix.c_str(), it->first.c_str(), it->second->fLatexName.c_str()) << endl;
    fTEX << formatTex(nEvtFile, Form("%s:%s:nEvtFile", fSuffix.c_str(), it->first.c_str()), 0) << endl;
    fTEX << formatTex(nCands, Form("%s:%s:nCands", fSuffix.c_str(), it->first.c_str()), 0) << endl;
    fTEX << formatTex(epscand, Form("%s:%s:epsCand:val", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    fTEX << formatTex(epscandE, Form("%s:%s:epsCand:err", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    fTEX << formatTex(it->second->fFilterEff, Form("%s:%s:filterEff:val", fSuffix.c_str(), it->first.c_str()), 6) << endl;
    fTEX << formatTex(0., Form("%s:%s:filterEff:err", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    if (it->second->fBf > 0.) {
      fTEX << formatTexErrSci(it->second->fBf, it->second->fBfE, Form("%s:%s:bf", fSuffix.c_str(), it->first.c_str()), 2) << endl;
      double eqLumi(0.);
      if (it->second->fFilterEff > 0) {
	eqLumi = nEvtFile/fCrossSection/it->second->fBf/it->second->fFilterEff;
	fTEX << formatTex(eqLumi, Form("%s:%s:eqLumi:val", fSuffix.c_str(), it->first.c_str()), 1) << endl;
	cout << it->first << ": eqLumi = " << eqLumi << " filterEff: " << it->second->fFilterEff << endl;
      }
    }
    if (it->second->fLumi > 0.) {
      fTEX << formatTex(it->second->fLumi, Form("%s:%s:lumi", fSuffix.c_str(), it->first.c_str()), 1) << endl;
    }
  }

}


// ----------------------------------------------------------------------
void plotResults::genSummary(std::string dsname, std::string dir) {
  TH1D *hpt    = new TH1D("pt", "pt", 50, 0., 50.0);
  TH1D *heta   = new TH1D("eta", "eta", 40, -4., 4.0);
  TH1D *tpt    = new TH1D("tpt", "pt (HLT)", 50, 0., 50.0); setFilledHist(tpt, kBlue, kYellow, 1000);
  TH1D *teta   = new TH1D("teta", "eta (HLT)", 40, -4., 4.0); setFilledHist(teta, kBlue, kYellow, 1000);
  TH1D *hm1eta = new TH1D("m1eta", "m1 eta", 40, -4., 4.0);
  TH1D *hm2eta = new TH1D("m2eta", "m2 eta", 40, -4., 4.0);
  TH1D *tm1eta = new TH1D("tm1eta", "m1 eta (HLT)", 40, -4., 4.0); setFilledHist(tm1eta, kBlue, kYellow, 1000);
  TH1D *tm2eta = new TH1D("tm2eta", "m2 eta (HLT)", 40, -4., 4.0); setFilledHist(tm2eta, kBlue, kYellow, 1000);
  TH1D *hketa  = new TH1D("keta", "kaon eta", 40, -4., 4.0);
  TH1D *hm1pt  = new TH1D("m1pt", "m1 pt", 100, 0., 10.0);
  TH1D *hm2pt  = new TH1D("m2pt", "m2 pt", 100, 0., 10.0);
  TH1D *tm1pt  = new TH1D("tm1pt", "m1 pt (HLT)", 100, 0., 10.0); setFilledHist(tm1pt, kBlue, kYellow, 1000);
  TH1D *tm2pt  = new TH1D("tm2pt", "m2 pt (HLT)", 100, 0., 10.0); setFilledHist(tm2pt, kBlue, kYellow, 1000);
  TH1D *hkpt   = new TH1D("kpt", "kaon pt", 50, 0., 10.0);
  TH1D *htau   = new TH1D("tau", "tau", 100, 0., 15.e-12);

  TTree *T = getTree(dsname, dir, "effTree");
  T->Draw("gtau>>tau");
  T->Draw("gpt>>pt");
  T->Draw("geta>>eta");

  T->Draw("g1pt>>m1pt");
  T->Draw("g2pt>>m2pt");
  T->Draw("g1eta>>m1eta");
  T->Draw("g2eta>>m2eta");

  T->Draw("gpt>>tpt", "hlt");
  T->Draw("g1pt>>tm1pt", "hlt");
  T->Draw("g2pt>>tm2pt", "hlt");

  T->Draw("geta>>teta", "hlt");
  T->Draw("g1eta>>tm1eta", "hlt");
  T->Draw("g2eta>>tm2eta", "hlt");


  bool addKaon(false);
  bool addHLT(true);
  if (string::npos != dsname.find("Bg")) {
    addHLT = false;
  }

  if (string::npos != dsname.find("bupsik")) {
    T->Draw("g3eta>>keta");
    T->Draw("g3pt>>kpt");
    addKaon = true;
  }
  if (string::npos != dsname.find("bspsiphi")) {
    T->Draw("g3eta>>keta");
    T->Draw("g4eta>>keta");
    T->Draw("g3pt>>kpt");
    T->Draw("g4pt>>kpt");
    addKaon = true;
  }

  tl->SetTextSize(0.05);
  makeCanvas(1);
  int ncol(4);
  if (addKaon) {
    c1->Divide(5,2);
    ncol = 5;
  } else {
    c1->Divide(4,2);
  }

  c1->cd(1);
  setTitles(hpt, "p_{T} [GeV]", "entries/bin");
  hpt->Draw();
  if (addHLT) {
    tpt->Draw("same");
    hpt->Draw("sameaxis");
  }

  c1->cd(2);
  hm1pt->Draw();
  setTitles(hm1pt, "p_{T} [GeV]", "entries/bin");
  if (addHLT) {
    tm1pt->Draw("same");
    hm1pt->Draw("sameaxis");
  }

  c1->cd(3);
  hm2pt->Draw();
  setTitles(hm2pt, "p_{T} [GeV]", "entries/bin");
  if (addHLT) {
    tm2pt->Draw("same");
    hm2pt->Draw("sameaxis");
  }

  if (addKaon) {
    c1->cd(ncol-1);
    setTitles(hkpt, "p_{T} [GeV]", "entries/bin");
    hkpt->Draw();
  }

  c1->cd(ncol);
  gPad->SetLogy(1);
  htau->Fit("expo", "l");
  TF1 *f = (TF1*)htau->GetFunction("expo");
  double chi2 = f->GetChisquare();
  int    ndf  = f->GetNDF();
  double t    = -1./f->GetParameter(1);
  double tE   = -t*f->GetParError(1)/f->GetParameter(1);
  t  *= 1.e12;
  tE *= 1.e12;

  c1->cd(ncol+1);
  setTitles(heta, "#eta", "entries/bin");
  heta->Draw();
  if (addHLT) {
    teta->Draw("same");
    heta->Draw("axissame");
  }
  tl->DrawLatexNDC(0.55, 0.3, "B");

  c1->cd(ncol+2);
  setTitles(hm1eta, "#eta", "entries/bin");
  hm1eta->Draw();
  if (addHLT) {
    tm1eta->Draw("same");
    hm1eta->Draw("sameaxis");
  }
  tl->DrawLatexNDC(0.40, 0.3, "leading muon");

  c1->cd(ncol+3);
  setTitles(hm2eta, "#eta", "entries/bin");
  hm2eta->Draw();
  if (addHLT) {
    tm2eta->Draw("same");
    hm2eta->Draw("sameaxis");
  }
  tl->DrawLatexNDC(0.35, 0.3, "subleading muon");

  if (addKaon) {
    c1->cd(2*ncol-1);
    setTitles(hketa, "#eta", "entries/bin");
    hketa->Draw();
    tl->DrawLatexNDC(0.35, 0.3, "kaon(s)");
  }

  c1->cd(2*ncol);
  tl->SetTextSize(0.07);
  if (addHLT) {
    tl->DrawLatexNDC(0.2, 0.9, Form("%s (HLT)", dsname.c_str()));
  } else {
    tl->DrawLatexNDC(0.2, 0.9, Form("%s ", dsname.c_str()));
  }
  tl->DrawLatexNDC(0.2, 0.8, Form("Events: %d", T->GetEntries()));

  tl->DrawLatexNDC(0.2, 0.35, Form("#tau"));
  tl->DrawLatexNDC(0.25, 0.35, Form("= (%3.2f #pm %5.2f)ps", t, tE));

  tl->DrawLatexNDC(0.2, 0.25, Form("#varepsilon"));
  tl->DrawLatexNDC(0.25, 0.25, Form("= %4.3f ", teta->GetSumOfWeights()/heta->GetSumOfWeights()));

  c1->SaveAs(Form("%s/genSummary-%s.pdf", fDirectory.c_str(), dsname.c_str()));

}



// ----------------------------------------------------------------------
void plotResults::fillAndSaveHistograms(int nevents) {
  // -- dump histograms
  fHistFile = TFile::Open(fHistFileName.c_str(), "RECREATE");
  cout << " opened, running on " << nevents << " entries" << endl;

  TTree *t(0);


  // -- rare backgrounds
  if (1) {
    resetHistograms();
    //    rareBgHists("nada", nevents);
  }

  // -- normalization
  if (1) {
    resetHistograms();
    setup("bupsikData");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents);
    saveHistograms(fSetup);

    resetHistograms();
    fSetup = "bupsikMcOff";
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents);
    //    otherNumbers(fSetup);
    saveHistograms(fSetup);

    resetHistograms();
    setup("bspsiphiData");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents);
    saveHistograms(fSetup);

    resetHistograms();
    fSetup = "bspsiphiMcOff";
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents);
    //    otherNumbers(fSetup);
    saveHistograms(fSetup);

    resetHistograms();
    setup("bmmData");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents);
    saveHistograms(fSetup);

    resetHistograms();
    setup("bdmmMcOff");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents);
    saveHistograms(fSetup);

    resetHistograms();
    setup("bsmmMcOff");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents);
    saveHistograms(fSetup);

  }

  fHistFile->Close();
}






// ----------------------------------------------------------------------
void plotResults::loopFunction1() {

  if (fChan < 0) return;

  double mass = fb.m;

  fhMassAbsNoCuts[fChan]->Fill(mass);
  fhAccAll[fChan]->Fill(TMath::Abs(fb.eta), fb.pt);
  fhAccPtAll[fChan]->Fill(fb.pt);
  fhAccEtaAll[fChan]->Fill(TMath::Abs(fb.eta));

  if (!fGoodAcceptance) return;
  fhAccPass[fChan]->Fill(TMath::Abs(fb.eta), fb.pt);
  fhAccPtPass[fChan]->Fill(fb.pt);
  fhAccEtaPass[fChan]->Fill(TMath::Abs(fb.eta));

  // -- this is the base, after the raw acceptance cuts
  fhMassNoCuts[fChan]->Fill(mass);

  if (fDoUseBDT) {
    if (!fGoodQ) return;
    if (!fGoodPvAveW8) return;
    if (!fGoodTracks) return;
    if (!fGoodTracksPt) return;
    if (!fGoodTracksEta) return;
    if (!fGoodBdtPt) return;
    if (!fGoodMuonsEta) return;
    if (!fGoodJpsiCuts) return;
    cout << "confusione assoluta" << endl;
  } else {
    if (!fGoodQ) return;
    if (!fGoodMuonsPt) return;
    if (!fGoodMuonsEta) return;
    if (!fGoodJpsiCuts) return;
    if (!fGoodPvAveW8) return;
    if (!fGoodMaxDoca) return;
    if (!fGoodLip) return;
    if (!fGoodLipS) return;
    if (!fGoodIp) return;
    if (!fGoodIpS) return;
    if (!fGoodPt) return;
    if (!fGoodEta) return;
    if (!fGoodAlpha) return;
    if (!fGoodChi2) return;
    if (!fGoodFLS) return;
    if (!fGoodCloseTrack) return;
    if (!fGoodIso) return;
    if (!fGoodDocaTrk) return;
  }

  fhMassWithAnaCuts[fChan]->Fill(mass);

  // -- weighted with fake rate
  fhW8MassWithAllCuts[fChan]->Fill(mass, fW8MisId);

  // -- blind version
  if (BMM == fMode && !(5.2 < mass && mass < 5.45)) {
    fhW8MassWithAllCutsBlind[fChan]->Fill(mass, fW8MisId);
  }

  // -- MUON ID
  if (false == fGoodMuonsID) return;
  fhMassWithMuonCuts[fChan]->Fill(mass);

  // -- Trigger
  if (false == fGoodHLT) return;
  fhMassWithTriggerCuts[fChan]->Fill(mass);

  fhMassWithAllCuts[fChan]->Fill(mass);
  if (BMM == fMode && !(5.2 < mass && mass < 5.45)) {
    fhMassWithAllCutsBlind[fChan]->Fill(mass);
  }

  if (fMode == BU2JPSIKP) {
    fhNo[fChan]->Fill(mass);
    fhNoC[fChan]->Fill(fb.cm);
  }

  if (fMode == BS2JPSIPHI) {
    fhCs[fChan]->Fill(mass);
    fhCsC[fChan]->Fill(fb.cm);
  }

  if (fMode == BD2JPSIKSTAR) {
    fhB0[fChan]->Fill(mass);
    fhB0C[fChan]->Fill(fb.cm);
  }


  if (fMode == BSMM && mass < fCuts[fChan]->mBsLo) return;
  if (fMode == BSMM && mass > fCuts[fChan]->mBsHi) return;
  if (fMode == BDMM && mass < fCuts[fChan]->mBdLo) return;
  if (fMode == BDMM && mass > fCuts[fChan]->mBdHi) return;
  if (fMode == BU2JPSIKP && mass < fNoLo) return;
  if (fMode == BU2JPSIKP && mass > fNoHi) return;
  if (fMode == BS2JPSIPHI && mass < fCsLo) return;
  if (fMode == BS2JPSIPHI && mass > fCsHi) return;
  if (fMode == BD2JPSIKSTAR && mass < fNoLo) return;
  if (fMode == BD2JPSIKSTAR && mass > fNoHi) return;

  fhMassWithMassCuts[fChan]->Fill(mass);






}

// ----------------------------------------------------------------------
void plotResults::loopFunction2() {

}


// ----------------------------------------------------------------------
void plotResults::loopFunction3() {


}



// ----------------------------------------------------------------------
void plotResults::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
  int nentries = Int_t(t->GetEntries());
  int nbegin(0), nend(nentries);
  if (nevts > 0 && nentries > nevts) {
    nentries = nevts;
    nbegin = 0;
    nend = nevts;
  }
  if (nevts > 0 && nstart > 0) {
    nentries = nstart + nevts;
    nbegin = nstart;
    if (nstart + nevts < t->GetEntries()) {
      nend = nstart + nevts;
    } else {
      nend = t->GetEntries();
    }
  }

  nentries = nend - nstart;

  int step(1000000);
  if (nentries < 5000000)  step = 500000;
  if (nentries < 1000000)  step = 100000;
  if (nentries < 100000)   step = 10000;
  if (nentries < 10000)    step = 1000;
  if (nentries < 1000)     step = 100;
  step = 500000;
  cout << "==> plotResults::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotResults::*pF)(void);
  if (ifunc == 1) pF = &plotResults::loopFunction1;
  if (ifunc == 2) pF = &plotResults::loopFunction2;
  if (ifunc == 3) pF = &plotResults::loopFunction3;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotResults::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotResults::loadFile loading files listed in " << files << endl;

  char buffer[1000];
  ifstream is(files.c_str());
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}

    string sbuffer = string(buffer);
    replaceAll(sbuffer, " ", "");
    replaceAll(sbuffer, "\t", "");
    if (sbuffer.size() < 1) continue;

    string::size_type m1 = sbuffer.find("lumi=");
    string stype = sbuffer.substr(5, m1-5);

    string::size_type m2 = sbuffer.find("file=");
    string slumi = sbuffer.substr(m1+5, m2-m1-6);
    string sfile = sbuffer.substr(m2+5);
    string sname("nada"), sdecay("nada");

    TFile *pF(0);
    dataset *ds(0);

    if (string::npos != stype.find("Data")) {
      // -- Charmonium
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("blablabla")) {
        sname = "bmmCharmonium";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

    } else if (string::npos != stype.find("mc")) {
      // -- MC
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2.;

      if (string::npos != stype.find("wrongreco,")) {
        sname = "wrongReco";
        sdecay = "wrongReco";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

      if (string::npos != stype.find("bctopsimunu,")) {
        sname = "bcpsimunuMc";
        sdecay = "bcpsimunu";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

    }

    if (sname != "nada") {
      ds->fLcolor = ds->fColor;
      ds->fFcolor = ds->fColor;
      ds->fName   = sdecay;
      ds->fFullName = sname;
      insertDataset(sname, ds);
    } else {
      delete ds;
    }
  }

  is.close();

  int cnt(0);
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << Form("   %20s: %70s %15s %9s %8s %8s", "Dataset name", "Filename", "LaTeX", "Lumi", "Eff", "BF") << endl;
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    // cout << it->first << endl;
    // cout << it->second->fName << endl;
    // cout << it->second->fF->GetName() << endl;
    cout << Form("%2d %20s: %70s %15s %8.1f %.2e %.2e",
		 cnt,
		 it->first.c_str(),
		 it->second->fF->GetName(),
		 it->second->fLatexName.c_str(),
		 it->second->fLumi,
		 it->second->fFilterEff,
		 it->second->fBf
		 )
	 << endl;
    ++cnt;
  }
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void plotResults::saveHistograms(string smode) {

  fHistFile->cd();
  fHistFile->mkdir(fSetup.c_str());
  fHistFile->cd(fSetup.c_str());
  TDirectory *dir = gDirectory;

  // int mode(0);
  // if (string::npos != smode.find("bsmmMc"))          mode = 0;
  // if (string::npos != smode.find("bdmmMc"))          mode = 1;
  // if (string::npos != smode.find("bsmmMcOffAcc"))    mode = 2;

  // if (string::npos != smode.find("bmmData"))         mode = 5;
  // if (string::npos != smode.find("SgDataAMS"))       mode = 6;

  // if (string::npos != smode.find("bupsikMc"))        mode = 10;
  // if (string::npos != smode.find("bupsikData"))      mode = 15;

  // if (string::npos != smode.find("bspsiphiMc"))      mode = 20;
  // if (string::npos != smode.find("bspsiphiData"))    mode = 25;

  TH1D *h1(0);
  TH2D *h2(0);

  for (unsigned int i = 0; i < fNchan; ++i) {
    string modifier = "cnc";
    h1 = (TH1D*)(fhGenAndAccNumbers[i]->Clone(Form("hGenAndAccNumbers_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetDirectory(dir);
    h1->Write();

    h2 = (TH2D*)(fhAccAll[i]->Clone(Form("hAccAll_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h2->SetDirectory(dir);
    h2->Write();
    h2 = (TH2D*)(fhAccPass[i]->Clone(Form("hAccPass_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h2->SetDirectory(dir);
    h2->Write();

    h1 = (TH1D*)(fhAccEtaAll[i]->Clone(Form("hAccEtaAll_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetDirectory(dir);
    h1->Write();
    h1 = (TH1D*)(fhAccEtaPass[i]->Clone(Form("hAccEtaPass_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhAccPtAll[i]->Clone(Form("hAccPtAll_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetDirectory(dir);
    h1->Write();
    h1 = (TH1D*)(fhAccPtPass[i]->Clone(Form("hAccPtPass_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhMassAbsNoCuts[i]->Clone(Form("hMassAbsNoCuts_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetTitle(Form("hMassAbsNoCuts_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhMassNoCuts[i]->Clone(Form("hMassNoCuts_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetTitle(Form("hMassNoCuts_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhMassWithAnaCuts[i]->Clone(Form("hMassWithAnaCuts_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetTitle(Form("hMassWithAnaCuts_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhMassWithMuonCuts[i]->Clone(Form("hMassWithMuonCuts_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetTitle(Form("hMassWithMuonCuts_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhMassWithTriggerCuts[i]->Clone(Form("hMassWithTriggerCuts_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetTitle(Form("hMassWithTriggerCuts_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhMassWithAllCuts[i]->Clone(Form("hMassWithAllCuts_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetTitle(Form("hMassWithAllCuts_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhMassWithMassCuts[i]->Clone(Form("hMassWithMassCuts_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetTitle(Form("hMassWithMassCuts_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();

    h1 = (TH1D*)(fhMassWithMassCuts[i]->Clone(Form("hMassWithAllCutsBlind_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
    h1->SetTitle(Form("hMassWithAllCutsBlind_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();


    if (string::npos != fSetup.find("bupsik")) {
      h1 = (TH1D*)(fhNo[i]->Clone(Form("hNo_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hNo_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhNoC[i]->Clone(Form("hNoC_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hNoC_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();
    }

    if (string::npos != fSetup.find("bspsiphi")) {
      h1 = (TH1D*)(fhCs[i]->Clone(Form("hCs_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hCs_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhCsC[i]->Clone(Form("hCsC_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hCsC_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

    }

    if (string::npos != fSetup.find("bdpsikstar")) {
      h1 = (TH1D*)(fhB0[i]->Clone(Form("hB0_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hB0_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhB0C[i]->Clone(Form("hB0C_%s_%s_chan%d", modifier.c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hB0C_%s_%s_chan%d %s", modifier.c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

    }



  }



}


// ----------------------------------------------------------------------
void plotResults::resetHistograms() {

  for (unsigned int i = 0; i < fNchan; ++i) {
    fhGenAndAccNumbers[i]->Reset();

    fhMassAbsNoCuts[i]->Reset();
    fhMassNoCuts[i]->Reset();

    fhMassWithAnaCuts[i]->Reset();

    fhMassWithMuonCuts[i]->Reset();

    fhMassWithTriggerCuts[i]->Reset();

    fhMassWithAllCuts[i]->Reset();
    fhMassWithAllCutsBlind[i]->Reset();

    fhW8MassWithAllCuts[i]->Reset();
    fhW8MassWithAllCutsBlind[i]->Reset();

    fhMassWithMassCuts[i]->Reset();

    fhNo[i]->Reset();
    fhNoC[i]->Reset();
    fhCs[i]->Reset();
    fhCsC[i]->Reset();
    fhB0[i]->Reset();
    fhB0C[i]->Reset();

    fhAccAll[i]->Reset();
    fhAccPass[i]->Reset();
    fhAccPtAll[i]->Reset();
    fhAccPtPass[i]->Reset();
    fhAccEtaAll[i]->Reset();
    fhAccEtaPass[i]->Reset();



  }

}
