#include "plotWork.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
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

ClassImp(plotWork)

using namespace std;

// ----------------------------------------------------------------------
plotWork::plotWork(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  loadFiles(files);

  if (setup == "") {
    fHistFileName = Form("%s/plotWork.root", dir.c_str());
  } else {
    fHistFileName = Form("%s/plotWork-%s.root", dir.c_str(), setup.c_str());
  }

  fTexFileName = fHistFileName;
  replaceAll(fTexFileName, ".root", ".tex");
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  fTEX.open(fTexFileName.c_str(), ios::app);


  if (cuts == "") {
    cout << "==> no cuts provided!?" << endl;
  } else {
    cout << Form("==> reading cuts/%s", cuts.c_str()) << endl;
    readCuts(Form("cuts/%s", cuts.c_str()));
    printCuts(cout);
  }

  fAnaCuts.addCut("fGoodHLT", "HLT", fGoodHLT);
  fAnaCuts.addCut("fGoodMuonsID", "lepton ID", fGoodMuonsID);
  fAnaCuts.addCut("fGoodMuonsPt", "p_{T,#mu} [GeV]", fGoodMuonsPt);
  fAnaCuts.addCut("fGoodMuonsEta", "#eta_{#mu}", fGoodMuonsEta);
  fAnaCuts.addCut("fGoodTracks", "good tracks", fGoodTracks);
  fAnaCuts.addCut("fGoodTracksPt", "p_{T,trk} [GeV]", fGoodTracksPt);
  fAnaCuts.addCut("fGoodTracksEta", "#eta_{trk} ", fGoodTracksEta);

  fAnaCuts.addCut("fGoodQ", "q_{1} 1_{2}", fGoodQ);
  fAnaCuts.addCut("fGoodPvAveW8", "<w8>", fGoodPvAveW8);
  fAnaCuts.addCut("fGoodIp", "IP", fGoodIp);
  fAnaCuts.addCut("fGoodIpS", "IPS", fGoodIpS);
  fAnaCuts.addCut("fGoodLip", "LIP", fGoodLip);
  fAnaCuts.addCut("fGoodLipS", "LIPS", fGoodLipS);
  fAnaCuts.addCut("fGoodMaxDoca", "MAXDOCA", fGoodMaxDoca);
  fAnaCuts.addCut("fGoodPt", "p_{T,B}", fGoodPt);
  fAnaCuts.addCut("fGoodEta", "#eta_{B}", fGoodEta);
  fAnaCuts.addCut("fGoodAlpha", "#alpha", fGoodAlpha);
  fAnaCuts.addCut("fGoodFLS", "l/#sigma(l)", fGoodFLS);
  fAnaCuts.addCut("fGoodChi2", "#chi^{2}", fGoodChi2);
  fAnaCuts.addCut("fGoodIso", "I_{trk}", fGoodIso);
  fAnaCuts.addCut("fGoodCloseTrack", "close track veto", fGoodCloseTrack);
  fAnaCuts.addCut("fGoodDocaTrk", "d_{ca}(trk)", fGoodDocaTrk);
  fAnaCuts.addCut("fGoodBDT", "bdt", fGoodBDT);
  fAnaCuts.addCut("fGoodLastCut", "lastCut", fGoodLastCut);


  fChan = 0;

}


// ----------------------------------------------------------------------
plotWork::~plotWork() {

}


// ----------------------------------------------------------------------
void plotWork::makeAll(int bitmask) {

  if (bitmask & 0x1) {
  }

}


// ----------------------------------------------------------------------
void plotWork::privateVsOfficial(string mode, int nevt) {
  MASSMIN = 4.5;
  MASSMAX = 6.5;
  string dir("candAnaMuMu");
  if (string::npos != mode.find("bdmm")) {
    fSetup = "SgMc";
    fIsSignal = true;
    BGLBOXMIN = 4.80;
    BGLBOXMAX = 5.20;
    SIGBOXMIN = 5.20;
    SIGBOXMAX = 5.45;
    BGHBOXMIN = 5.45;
    BGHBOXMAX = 6.00;
  }

  if (string::npos != mode.find("bujpsikp")) {
    BGLBOXMIN = 5.00;
    BGLBOXMAX = 5.18;
    SIGBOXMIN = 5.23;
    SIGBOXMAX = 5.33;
    BGHBOXMIN = 5.40;
    BGHBOXMAX = 5.50;
  }

  if (string::npos != mode.find("bsjpsiphi")) {
    BGLBOXMIN = 5.10;
    BGLBOXMAX = 5.29;
    SIGBOXMIN = 5.34;
    SIGBOXMAX = 5.40;
    BGHBOXMIN = 5.45;
    BGHBOXMAX = 5.70;
  }

  TTree *t(0);
  string lmode = mode;
  fDS[lmode]->fName = "private";
  t = getTree(lmode, dir);
  setupTree(t, fSetup);
  fOffset = 0;
  bookDistributions(lmode);
  int nevents(-1);
  if (nevt > -1) nevents = nevt;
  loopOverTree(t, 1, nevents);

  lmode = mode + "_official";
  fDS[lmode]->fName = "official";
  fDS[lmode]->fColor = kRed+2;
  fDS[lmode]->fLcolor = kRed+2;
  fDS[lmode]->fFcolor = kRed+2;
  fDS[lmode]->fFillStyle = 3356;
  t = getTree(lmode, dir);
  setupTree(t, fSetup);
  fOffset = 1;
  bookDistributions(lmode);
  loopOverTree(t, 1, nevents);

  // a0->hPresel[2]->Draw("e");
  // a1->hPresel[2]->Scale(a0->hPresel[2]->GetSumOfWeights()/a1->hPresel[2]->GetSumOfWeights());
  // a1->hPresel[2]->Draw("samehist");

  //  gStyle->SetOptTitle(0);


  TH1* h0 = fpMuon1Pt[0]->hPresel[2];
  TH1* h1 = fpMuon1Pt[1]->hPresel[2];  h0->SetTitle("");
  fDS[mode]->setHistStyle(h0);  fDS[lmode]->setHistStyle(h1);   overlay(h0, mode, h1, lmode, 0, "", UNITY, false, true, 0.4, 0.7);
  savePad(Form("mcval-m1pt-%s-%s.pdf", mode.c_str(), lmode.c_str()));

  h0 = fpMuon2Pt[0]->hPresel[2];  h1 = fpMuon2Pt[1]->hPresel[2];  h0->SetTitle("");
  fDS[mode]->setHistStyle(h0);  fDS[lmode]->setHistStyle(h1);   overlay(h0, mode, h1, lmode, 0, "", UNITY, false, true, 0.4, 0.7);
  savePad(Form("mcval-m2pt-%s-%s.pdf", mode.c_str(), lmode.c_str()));

  h0 = fpPt[0]->hPresel[2];  h1 = fpPt[1]->hPresel[2]; h0->SetTitle("");
  fDS[mode]->setHistStyle(h0);  fDS[lmode]->setHistStyle(h1);   overlay(h0, mode, h1, lmode, 0, "", UNITY, false, true, 0.4, 0.7);
  savePad(Form("mcval-pt-%s-%s.pdf", mode.c_str(), lmode.c_str()));

  h0 = fpIso[0]->hPresel[2];  h1 = fpIso[1]->hPresel[2]; h0->SetTitle("");
  fDS[mode]->setHistStyle(h0);  fDS[lmode]->setHistStyle(h1);   overlay(h0, mode, h1, lmode, 0, "", UNITY, false, true, 0.4, 0.7);
  savePad(Form("mcval-iso-%s-%s.pdf", mode.c_str(), lmode.c_str()));

  h0 = fpFLS3d[0]->hPresel[2];  h1 = fpFLS3d[1]->hPresel[2]; h0->SetTitle("");
  fDS[mode]->setHistStyle(h0);  fDS[lmode]->setHistStyle(h1);   overlay(h0, mode, h1, lmode, 0, "", UNITY, false, true, 0.4, 0.7);
  savePad(Form("mcval-fls3d-%s-%s.pdf", mode.c_str(), lmode.c_str()));

  h0 = fpIpS[0]->hPresel[2];  h1 = fpIpS[1]->hPresel[2]; h0->SetTitle("");
  fDS[mode]->setHistStyle(h0);  fDS[lmode]->setHistStyle(h1);   overlay(h0, mode, h1, lmode, 0, "", UNITY, false, true, 0.4, 0.7);
  savePad(Form("mcval-ips-%s-%s.pdf", mode.c_str(), lmode.c_str()));

  h0 = fpIp[0]->hPresel[2];  h1 = fpIp[1]->hPresel[2];  h0->SetTitle("");
  fDS[mode]->setHistStyle(h0);  fDS[lmode]->setHistStyle(h1);   overlay(h0, mode, h1, lmode, 0, "", UNITY, false, true, 0.4, 0.7);
  savePad(Form("mcval-ip-%s-%s.pdf", mode.c_str(), lmode.c_str()));

  h0 = fpPvN[0]->hPresel[2];  h1 = fpPvN[1]->hPresel[2]; h0->SetTitle("");
  fDS[mode]->setHistStyle(h0);  fDS[lmode]->setHistStyle(h1);   overlay(h0, mode, h1, lmode, 0, "", UNITY, false, true, 0.4, 0.7);
  savePad(Form("mcval-pvn-%s-%s.pdf", mode.c_str(), lmode.c_str()));
}


// ----------------------------------------------------------------------
void plotWork::loopFunction1() {

  double mass = fb.cm;
  if (fIsMC) mass = fb.m;
  if (fIsSignal) mass = fb.m;

  TLorentzVector a;
  a.SetPtEtaPhiM(fb.pt,fb.eta,fb.phi,mass);

  fpMuon1Pt[fOffset]->fill(fb.m1pt, mass);
  fpMuon2Pt[fOffset]->fill(fb.m2pt, mass);

  fpMuonsEta[fOffset]->fill(fb.m1eta, mass);
  fpMuonsEta[fOffset]->fill(fb.m2eta, mass);
  fpPt[fOffset]->fill(fb.pt, mass);
  fpP[fOffset]->fill(a.P(), mass);
  fpPz[fOffset]->fill(a.Pz(), mass);
  fpEta[fOffset]->fill(fb.eta, mass);
  fpAlpha[fOffset]->fill(fb.alpha, mass);

  fpIso[fOffset]->fill(fb.iso, mass);
  fpCloseTrk[fOffset]->fill(fb.closetrk, mass);
  fpDocaTrk[fOffset]->fill(fb.docatrk, mass);

  fpChi2Dof[fOffset]->fill(fb.chi2/fb.dof, mass);
  fpPChi2Dof[fOffset]->fill(fb.pchi2dof, mass);

  fpFLS3d[fOffset]->fill(fb.fls3d, mass);
  fpFL3d[fOffset]->fill(fb.fl3d, mass);
  fpFL3dE[fOffset]->fill(fb.fl3dE, mass);

  fpMaxDoca[fOffset]->fill(fb.maxdoca, mass);
  fpIp[fOffset]->fill(fb.pvip, mass);
  fpIpS[fOffset]->fill(fb.pvips, mass);
  fpPvZ[fOffset]->fill(fb.pvz, mass);
  fpPvN[fOffset]->fill(fb.pvn, mass);
  fpPvAveW8[fOffset]->fill(fb.pvw8, mass);


}



// ----------------------------------------------------------------------
void plotWork::bookDistributions(string mode) {
  string name = Form("%s_%s_", fSetup.c_str(), mode.c_str());
  fpMuon1Pt[fOffset]   = bookDistribution(Form("%smuon1pt", name.c_str()), "p_{T, #mu1} [GeV]", "fGoodMuonsPt", 60, 0., 30.);
  fpMuon2Pt[fOffset]   = bookDistribution(Form("%smuon2pt", name.c_str()), "p_{T, #mu2} [GeV]", "fGoodMuonsPt", 40, 0., 20.);
  fpMuonsEta[fOffset]  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsEta", 40, -2.5, 2.5);
  fpPt[fOffset]        = bookDistribution(Form("%spt", name.c_str()), "p_{T}(B) [GeV]", "fGoodPt", 60, 0., 60.);
  fpP[fOffset]         = bookDistribution(Form("%sp", name.c_str()), "p(B) [GeV]", "fGoodPt", 50, 0., 100.);
  fpPz[fOffset]        = bookDistribution(Form("%spz", name.c_str()), "p_{z}(B) [GeV]", "fGoodPt", 50, 0., 100.);
  fpEta[fOffset]       = bookDistribution(Form("%seta", name.c_str()), "#eta(B)", "fGoodEta", 40, -2.5, 2.5);
  fpAlpha[fOffset]     = bookDistribution(Form("%salpha", name.c_str()), "#alpha_{3D}", "fGoodAlpha", 50, 0., 0.1);
  fpIso[fOffset]       = bookDistribution(Form("%siso", name.c_str()),  "isolation", "fGoodIso", 52, 0., 1.04);
  fpCloseTrk[fOffset]  = bookDistribution(Form("%sclosetrk", name.c_str()),  "N_{trk}^{close}", "fGoodCloseTrack", 10, 0., 10.);
  fpDocaTrk[fOffset]   = bookDistribution(Form("%sdocatrk", name.c_str()), "d_{ca}^{0} [cm]", "fGoodDocaTrk", 50, 0., 0.20);

  fpChi2Dof[fOffset]   = bookDistribution(Form("%schi2dof", name.c_str()),  "#chi^{2}/dof", "fGoodChi2", 40, 0., 4.);
  fpPChi2Dof[fOffset]  = bookDistribution(Form("%spchi2dof", name.c_str()),  "P(#chi^{2},dof)", "fGoodChi2", 50, 0., 1.0);

  fpFLS3d[fOffset]     = bookDistribution(Form("%sfls3d", name.c_str()), "l_{3D}/#sigma(l_{3D})", "fGoodFLS", 60, 0., 120.);
  fpFL3d[fOffset]      = bookDistribution(Form("%sfl3d", name.c_str()),  "l_{3D} [cm]", "fGoodFLS", 60, 0., 1.5);
  fpFL3dE[fOffset]     = bookDistribution(Form("%sfl3de", name.c_str()), "#sigma(l_{3D}) [cm]", "fGoodFLS", 50, 0., 0.05);

  fpMaxDoca[fOffset]   = bookDistribution(Form("%smaxdoca", name.c_str()), "d^{max} [cm]", "fGoodMaxDoca", 60, 0., 0.03);
  fpIp[fOffset]        = bookDistribution(Form("%sip", name.c_str()), "#delta_{3D} [cm]", "fGoodIp", 50, 0., 0.020);
  fpIpS[fOffset]       = bookDistribution(Form("%sips", name.c_str()), "#delta_{3D}/#sigma(#delta_{3D})", "fGoodIpS", 50, 0., 4);

  fpPvZ[fOffset]       = bookDistribution(Form("%spvz", name.c_str()), "z_{PV} [cm]", "fGoodHLT", 40, -20., 20.);
  fpPvN[fOffset]       = bookDistribution(Form("%spvn", name.c_str()), "N(PV) ", "fGoodHLT", 40, 0., 40.);
  fpPvAveW8[fOffset]   = bookDistribution(Form("%spvavew8", name.c_str()), "<w^{PV}>", "fGoodPvAveW8", 50, 0.5, 1.);


}



// ----------------------------------------------------------------------
AnalysisDistribution* plotWork::bookDistribution(string hn, string ht, string hc, int nbins, double lo, double hi) {
  AnalysisDistribution *p = new AnalysisDistribution(hn.c_str(), ht.c_str(), nbins, lo, hi);
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX);
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX);
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX);
  p->setAnalysisCuts(&fAnaCuts, hc.c_str());
  p->setPreselCut(&fPreselection);

  return p;
}


// ----------------------------------------------------------------------
void plotWork::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotWork::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotWork::*pF)(void);
  if (ifunc == 1) pF = &plotWork::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
// NOTE: This works with the output of genAnalysis!!!!
void plotWork::prodSummary(string ds1, int year) {
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

  static const double aparticles[] = {511, 521, 531, 5122};
  vector<int> particles(aparticles, aparticles + sizeof(aparticles)/sizeof(aparticles[0]));

  if (fDS[ds1]) {
    fDS[ds1]->cd("");
  } else {
    cout << "fDS[" << ds1 << "] not found" << endl;
    return;
  }
  // -- Masses
  cout << "Masses" << endl;
  cout << Form("B0: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
	       ((TH1D*)gFile->Get("m511"))->GetMean(),
	       M511,
	       ((TH1D*)gFile->Get("m511"))->GetMean() - M511
	       )
       << endl;

  cout << Form("B+: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
	       ((TH1D*)gFile->Get("m521"))->GetMean(),
	       M521,
	       ((TH1D*)gFile->Get("m521"))->GetMean() - M521
	       )
       << endl;

  cout << Form("Bs: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
	       ((TH1D*)gFile->Get("m531"))->GetMean(),
	       M531,
	       ((TH1D*)gFile->Get("m531"))->GetMean() - M531
	       )
       << endl;

  cout << Form("Lb: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
	       ((TH1D*)gFile->Get("m5122"))->GetMean(),
	       M5122,
	       ((TH1D*)gFile->Get("m5122"))->GetMean() - M5122
	       )
       << endl;

  if (((TH1D*)gFile->Get("m541"))->GetEntries() > 1)
    cout << Form("Bc: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
		 ((TH1D*)gFile->Get("m541"))->GetMean(),
		 M541,
		 ((TH1D*)gFile->Get("m541"))->GetMean() - M541
		 )
	 << endl;



  // -- Lifetimes
  cout << "Lifetime" << endl;
  double t(1.), tE(1.), chi2(0);
  int ndf(0);
  TH1D *h = (TH1D*)gFile->Get("t521");
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
  h = (TH1D*)gFile->Get("t511");
  if (h->GetEntries() > 100) {
    h->Fit("expo", "ql");
    //      gPad->SaveAs("t511.pdf");
    f = (TF1*)h->GetFunction("expo");
    chi2 = f->GetChisquare()/f->GetNDF();
    t    = -1./f->GetParameter(1);
    tE   = -t*f->GetParError(1)/f->GetParameter(1);
    cout << Form("B0: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L511, (t-L511)/tE, chi2) << endl;
  }

  h = (TH1D*)gFile->Get("t531");
  if (h->GetEntries() > 100) {
    h->Fit("expo", "ql");
    //      gPad->SaveAs("t531.pdf");
    f = (TF1*)h->GetFunction("expo");
    chi2 = f->GetChisquare()/f->GetNDF();
    t    = -1./f->GetParameter(1);
    tE   = -t*f->GetParError(1)/f->GetParameter(1);
    cout << Form("Bs: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L531, (t-L531)/tE, chi2) << endl;
  }

  h = (TH1D*)gFile->Get("t5122");
  if (h->GetEntries() > 100) {
    h->Fit("expo", "ql");
    //      gPad->SaveAs("t5122.pdf");
    f = (TF1*)h->GetFunction("expo");
    chi2 = f->GetChisquare()/f->GetNDF();
    t    = -1./f->GetParameter(1);
    tE   = -t*f->GetParError(1)/f->GetParameter(1);
    cout << Form("Lb: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L5122, (t-L5122)/tE, chi2) << endl;
  }

  h = (TH1D*)gFile->Get("t541");
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
  tl->SetTextSize(0.03);
  gPad->SetLogy(1);
  h = ((TH1D*)gFile->Get("t521"));
  f = (TF1*)h->GetFunction("expo");
  chi2 = f->GetChisquare();
  ndf  = f->GetNDF();
  t    = -1./f->GetParameter(1);
  tE   = -t*f->GetParError(1)/f->GetParameter(1);
  h->SetTitle("");
  setTitles(h, "c#tau [#mum]", "Entries/bin");
  h->Draw();
  setItalic(); tl->DrawLatexNDC(0.2, 0.45, "B^{+}"); setRoman();
  tl->DrawLatexNDC(0.2, 0.40, Form("#chi^{2}/ndf")); tl->DrawLatexNDC(0.3, 0.40, Form("= %3.1f/%d", chi2, ndf));
  tl->DrawLatexNDC(0.2, 0.35, Form("#tau(fit)"));    tl->DrawLatexNDC(0.3, 0.35, Form("= %5.2f #pm %5.2f", t, tE));
  tl->DrawLatexNDC(0.2, 0.30, Form("#tau(pdg)"));    tl->DrawLatexNDC(0.3, 0.30, Form("= %5.2f ", L521));
  tl->DrawLatexNDC(0.2, 0.25, Form("pull"));         tl->DrawLatexNDC(0.3, 0.25, Form("= %+5.2f ", (t-L521)/tE));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  savePad("t521.pdf");

  h = ((TH1D*)gFile->Get("t511"));
  f = (TF1*)h->GetFunction("expo");
  chi2 = f->GetChisquare();
  ndf  = f->GetNDF();
  t    = -1./f->GetParameter(1);
  tE   = -t*f->GetParError(1)/f->GetParameter(1);
  h->SetTitle("");
  setTitles(h, "c#tau [#mum]", "Entries/bin");
  h->Draw();
  setItalic(); tl->DrawLatexNDC(0.2, 0.45, "B^{0}"); setRoman();
  tl->DrawLatexNDC(0.2, 0.40, Form("#chi^{2}/ndf")); tl->DrawLatexNDC(0.3, 0.40, Form("= %3.1f/%d", chi2, ndf));
  tl->DrawLatexNDC(0.2, 0.35, Form("#tau(fit)"));    tl->DrawLatexNDC(0.3, 0.35, Form("= %5.2f #pm %5.2f", t, tE));
  tl->DrawLatexNDC(0.2, 0.30, Form("#tau(pdg)"));    tl->DrawLatexNDC(0.3, 0.30, Form("= %5.2f ", L511));
  tl->DrawLatexNDC(0.2, 0.25, Form("pull"));         tl->DrawLatexNDC(0.3, 0.25, Form("= %+5.2f ", (t-L511)/tE));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  savePad("t511.pdf");

  h = ((TH1D*)gFile->Get("t531"));
  f = (TF1*)h->GetFunction("expo");
  chi2 = f->GetChisquare();
  ndf  = f->GetNDF();
  t    = -1./f->GetParameter(1);
  tE   = -t*f->GetParError(1)/f->GetParameter(1);
  h->SetTitle("");
  setTitles(h, "c#tau [#mum]", "Entries/bin");
  h->Draw();
  setItalic(); tl->DrawLatexNDC(0.2, 0.45, "B_{s}"); setRoman();
  tl->DrawLatexNDC(0.2, 0.40, Form("#chi^{2}/ndf")); tl->DrawLatexNDC(0.3, 0.40, Form("= %3.1f/%d", chi2, ndf));
  tl->DrawLatexNDC(0.2, 0.35, Form("#tau(fit)"));    tl->DrawLatexNDC(0.3, 0.35, Form("= %5.2f #pm %5.2f", t, tE));
  tl->DrawLatexNDC(0.2, 0.30, Form("#tau(pdg)"));    tl->DrawLatexNDC(0.3, 0.30, Form("= %5.2f ", L531));
  tl->DrawLatexNDC(0.2, 0.25, Form("pull"));         tl->DrawLatexNDC(0.3, 0.25, Form("= %+5.2f ", (t-L531)/tE));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  savePad("t531.pdf");

  h = ((TH1D*)gFile->Get("t5122"));
  f = (TF1*)h->GetFunction("expo");
  chi2 = f->GetChisquare();
  ndf  = f->GetNDF();
  t    = -1./f->GetParameter(1);
  tE   = -t*f->GetParError(1)/f->GetParameter(1);
  h->SetTitle("");
  setTitles(h, "c#tau [#mum]", "Entries/bin");
  h->Draw();
  setItalic(); tl->DrawLatexNDC(0.2, 0.45, "#Lambda_{b}"); setRoman();
  tl->DrawLatexNDC(0.2, 0.40, Form("#chi^{2}/ndf")); tl->DrawLatexNDC(0.3, 0.40, Form("= %3.1f/%d", chi2, ndf));
  tl->DrawLatexNDC(0.2, 0.35, Form("#tau(fit)"));    tl->DrawLatexNDC(0.3, 0.35, Form("= %5.2f #pm %5.2f", t, tE));
  tl->DrawLatexNDC(0.2, 0.30, Form("#tau(pdg)"));    tl->DrawLatexNDC(0.3, 0.30, Form("= %5.2f ", L5122));
  tl->DrawLatexNDC(0.2, 0.25, Form("pull"));         tl->DrawLatexNDC(0.3, 0.25, Form("= %+5.2f ", (t-L5122)/tE));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  savePad("t5122.pdf");
}


// ----------------------------------------------------------------------
void plotWork::setCuts(string cuts) {
  cout << "==> plotWork::setCuts: " << cuts << endl;

  istringstream ss(cuts);
  string token, name, sval;

  while (getline(ss, token, ',')) {

    string::size_type m1 = token.find("=");
    name = token.substr(0, m1);
    sval = token.substr(m1+1);

    if (string::npos != name.find("PTLO")) {
      float val;
      val = atof(sval.c_str());
      PTLO = val;
    }

  }
}


// ----------------------------------------------------------------------
void plotWork::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> Loading files listed in " << files << endl;

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
    string sname, sdecay;

    //    cout << "stype: ->" << stype << "<-" << endl;

    TFile *pF(0);
    if (string::npos != stype.find("data")) {
      // -- DATA
      pF = loadFile(sfile);

      dataset *ds = new dataset();
      ds->fSize = 1;
      ds->fWidth = 2;

      if (string::npos != stype.find("bmm")) {
        sname = "data_bmm";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bu2jpsik")) {
        sname = "data_bu2jpsik";
        sdecay = "bu2jpsik";
	ds->fColor = kBlack;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      ds->fLcolor = ds->fColor;
      ds->fFcolor = ds->fColor;
      ds->fName   = sdecay;
      ds->fFullName = sname;
      fDS.insert(make_pair(sname, ds));


    } else {
      // -- MC
      pF = loadFile(sfile);
      //      cout << "  " << sfile << ": " << pF << endl;

      dataset *ds = new dataset();
      ds->fSize = 1;
      ds->fWidth = 2;

      string filter = "";
      if (string::npos != stype.find("acc")) filter = "_acc";

      string prod = "";
      if (string::npos != stype.find("official")) filter = "_official";

      // -------------------------
      // -- genAnalysis files below
      // -------------------------
      if (string::npos != stype.find("bdmm,")) {
        sname = "bdmm" + filter + prod;
        sdecay = "bdmm";
	ds->fColor = kBlue-7;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      //      cout << "  inserting as name ->" << sname << "<- and decay = " << sdecay << endl;
      ds->fLcolor = ds->fColor;
      ds->fFcolor = ds->fColor;
      ds->fName   = sdecay;
      ds->fFullName = sname;
      fDS.insert(make_pair(sname, ds));



    }


  }

  is.close();
  cout << "Summary of files loaded: " << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    cout << "    " << it->first << ": " << it->second->fName << ", " << it->second->fF->GetName() << endl;
  }
}
