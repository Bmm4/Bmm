#include "plotClass.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"

#include "common/util.hh"

#include "preselection.hh"

ClassImp(plotClass)

using namespace std;

// ----------------------------------------------------------------------
plotClass::plotClass(string dir, string files, string cuts, string setup) {

  setTdrStyle();

  gStyle->SetHatchesSpacing(2);

  fDBX = true;
  fDoUseBDT = false;
  fVerbose = true;

  fDirectory = dir;
  fSetup = setup;
  fSuffix = setup;
  fMode = UNSET;

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

  fEpsilon = 0.00001;
  fLumi = 20.;

  legg = 0;
  c0 = c1 = c2 = c3 = c4 = c5 =0;
  tl = new TLatex();
  box = new TBox();
  pa = new TArrow();
  pl = new TLine();
  legge = 0;

  fAccPt = 3.5;
  fAccEtaGen = 2.5;
  fAccEtaRec = 2.4;

  c0 = (TCanvas*)gROOT->FindObject("c0");
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);

  fHistFile = 0; // this must be opened in a derived class!

  fStampString = "preliminary";
  if (fDoUseBDT) {
    fStampString = "BDT";
  } else {
    fStampString = "CNC";
  }
  fStampCms = "BMM4";
  fStampLumi = "2.2 fb^{-1}";

  string sfiles(files);
  if (string::npos != sfiles.find("2011")) {
    fYear = 2011;
    fStampLumi = "L = 5 fb^{-1} (#sqrt{s} = 7 TeV)";
  }
  if (string::npos != sfiles.find("2012")) {
    fYear = 2012;
    fStampLumi = "L = 20 fb^{-1} (#sqrt{s} = 8 TeV)";
  }
  if (string::npos != sfiles.find("2015")) {
    fYear = 2015;
    fStampLumi = "L = 2.5 fb^{-1} (#sqrt{s} = 13 TeV)";
  }
  if (string::npos != sfiles.find("2016")) {
    fYear = 2016;
    fStampLumi = "L = 13.1 fb^{-1} (#sqrt{s} = 13 TeV)";
  }
  if (setup == "") fSuffix = Form("%d", fYear);

  fIF = new initFunc();

  if (1) {
    int year(2012);
    string directory("../common/pidtables/");
    string name("");
    name = directory + Form("%d-kaonPosFakeRate-mvaMuon.dat", year); fptFakePosKaons     = new PidTable(Form(name.c_str()));
    name = directory + Form("%d-kaonNegFakeRate-mvaMuon.dat", year); fptFakeNegKaons     = new PidTable(Form(name.c_str()));
    name = directory + Form("%d-pionPosFakeRate-mvaMuon.dat", year); fptFakePosPions     = new PidTable(Form(name.c_str()));
    name = directory + Form("%d-pionNegFakeRate-mvaMuon.dat", year); fptFakeNegPions     = new PidTable(Form(name.c_str()));
    name = directory + Form("%d-protonPosFakeRate-mvaMuon.dat", year); fptFakePosProtons = new PidTable(Form(name.c_str()));
    name = directory + Form("%d-protonNegFakeRate-mvaMuon.dat", year); fptFakeNegProtons = new PidTable(Form(name.c_str()));

    name = directory + Form("%d-L1L2_data_all.dat", year);        fptT1     = new PidTable(name.c_str());
    name = directory + Form("%d-L3_data_all.dat", year);          fptT2     = new PidTable(name.c_str());
    name = directory + Form("%d-MuonID_data_all.dat", year);      fptM      = new PidTable(name.c_str());
  }

  fAnaCuts.addCut("fGoodHLT", "HLT", fGoodHLT);
  fAnaCuts.addCut("fGoodMuonsID", "lepton ID", fGoodMuonsID);
  fAnaCuts.addCut("fGoodMuonsPt", "p_{T,#mu} [GeV]", fGoodMuonsPt);
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

  // -- NOTE: This should be synchronized to AN-16-178/trunk/symbols.tex
  fVarToTex.insert(make_pair("mpt", "p_{T_{#mu}} #it{[GeV]}"));
  fVarToTex.insert(make_pair("m1pt", "p_{T_{#mu,1}} #it{[GeV]}"));
  fVarToTex.insert(make_pair("m2pt", "p_{T_{#mu,2}} #it{[GeV]}"));
  fVarToTex.insert(make_pair("meta", "#eta_{#mu}"));
  fVarToTex.insert(make_pair("m1eta", "#eta_{#mu,1}"));
  fVarToTex.insert(make_pair("m2eta", "#eta_{#mu,2}"));
  fVarToTex.insert(make_pair("fls3d", "l_{3D}/#sigma(l_{3D})"));
  fVarToTex.insert(make_pair("alpha", "#alpha_{3D}"));
  fVarToTex.insert(make_pair("chi2dof", "#chi^{2}/dof"));

  fVarToTex.insert(make_pair("iso", "isolation"));
  fVarToTex.insert(make_pair("m1iso", "#mu_{1} isolation"));
  fVarToTex.insert(make_pair("m2iso", "#mu_{2} isolation"));
  fVarToTex.insert(make_pair("docatrk", "d_{ca}^{0} #it{[cm]}"));
  fVarToTex.insert(make_pair("closetrk", "N_{trk}^{close}"));
  fVarToTex.insert(make_pair("closetrks1", "N_{trk}^{close, 1#sigma}"));
  fVarToTex.insert(make_pair("closetrks2", "N_{trk}^{close, 2#sigma}"));
  fVarToTex.insert(make_pair("closetrks3", "N_{trk}^{close, 3#sigma}"));

  fVarToTex.insert(make_pair("maxdoca", "d^{max} #it{[cm]}"));
  fVarToTex.insert(make_pair("pvip", "#delta_{3D} #it{[cm]}"));
  fVarToTex.insert(make_pair("pvips", "#delta_{3D}/#sigma(#delta_{3D})"));

  fVarToTex.insert(make_pair("pt", "p_{T_{B}} #it{[GeV]}"));
  fVarToTex.insert(make_pair("eta", "#eta_{B}"));

  // -- initialize cuts
  cout << "==> Reading cuts from " << Form("%s", cuts.c_str()) << endl;
  readCuts(Form("%s/%s", fDirectory.c_str(), cuts.c_str()));
  fNchan = fCuts.size();
}

// ----------------------------------------------------------------------
plotClass::~plotClass() {
}


// ----------------------------------------------------------------------
void plotClass::changeSetup(string dir, string name, string setup) {
  if (setup == "") {
    fHistFileName = Form("%s/%s.%d.root", dir.c_str(), name.c_str(), fYear);
    fNumbersFileName = fDirectory + Form("/%s.%d.txt", name.c_str(), fYear);
  } else {
    fHistFileName = Form("%s/%s-%s.%d.root", dir.c_str(), name.c_str(), setup.c_str(), fYear);
    fNumbersFileName = fDirectory + Form("/%s.%s.%d.txt", name.c_str(), setup.c_str(), fYear);
  }

  fTexFileName = fNumbersFileName;
  replaceAll(fTexFileName, ".txt", ".tex");
  fTEX.open(fTexFileName.c_str(), ios::app);
}


// ----------------------------------------------------------------------
void plotClass::init() {
}

// ----------------------------------------------------------------------
// see http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=15054
void plotClass::closeHistFile() {
  fHistFile->Write();
}

// ----------------------------------------------------------------------
void plotClass::cd(std::string dataset, std::string dir) {
  if (0 == fDS.count(dataset)) {
    cout << "unknown dataset: " << dataset << endl;
  } else {
    fDS[dataset]->cd(dir.c_str());
  }
}

// ----------------------------------------------------------------------
void plotClass::bookHist(string name) {
  cout << "==> plotClass: bookHist " << name << endl;
}

// ----------------------------------------------------------------------
void plotClass::makeAll(int bitmask) {
  cout << "==> plotClass: makeAll " << bitmask << endl;
}

// ----------------------------------------------------------------------
void plotClass::treeAnalysis() {
  cout << "==> plotClass: treeAnalysis " << endl;
}

// ----------------------------------------------------------------------
void plotClass::normHist(TH1 *h, string ds, int method) {
  double scale(1.);
  string smethod("");
  // -- normalize to 1
  if (method == UNITY) {
    smethod = "unity";
    scale = (h->Integral() > 0 ? 1./h->Integral() : 1.);
    h->GetYaxis()->SetTitle("normalized to 1");
  } else if (method == SOMETHING) {
    smethod = "something";
    scale = fNorm * (h->Integral() > 0 ? fNorm/h->Integral() : 1.);
    h->GetYaxis()->SetTitle("weighted events");
  } else if (method == XSECTION) {
    smethod = "xsection";
    // -- normalize to EFFECTIVE xsec*bf (EFFECTIVE to account for cuts)
    //    the cross section is known for ds
    //    ds corresponds to know lumi
    //
    //    n = xsec * L
    //    "integral" over histogram should be EFFECTIVE xsec
    scale = (h->Integral() > 0 ? fDS[ds]->fXsec*fDS[ds]->fBf/h->Integral() : 1.);
    h->GetYaxis()->SetTitle("pb");
  } else if (method == LUMI) {
    smethod = "lumi";
    // -- normalize to xsec*bf
    //    n = xsec * L
    //    "integral" over histogram should be events expected in fLumi
    scale = (h->Integral() > 0 ? fLumi/fDS[ds]->fLumi : 1.);
    h->GetYaxis()->SetTitle(Form("events in %4.0f/fb", fLumi));
  } else if (method == NONORM) {
    smethod = "nonorm";
    scale = 1.;
  } else {
    scale = 1.;
  }

  cout << "==>plotClass:  normHist scaling by " << scale << ", based on method " << smethod << endl;

  double c(0.), e(0.);
  for (int i = 0; i <= h->GetNbinsX(); ++i) {
    c = h->GetBinContent(i);
    e = h->GetBinError(i);
    h->SetBinContent(i, c*scale);
    h->SetBinError(i, e*scale);
  }
}


// ----------------------------------------------------------------------
void plotClass::overlay(TH1* h1, string f1, TH1* h2, string f2, TH1* h3, string f3, int method, bool loga, bool legend, double xleg, double yleg) {
  const bool verbose(false);

  showOverflow(h1);
  showOverflow(h2);
  if (h3) showOverflow(h3);

  normHist(h1, f1, method);
  normHist(h2, f2, method);
  if (h3) normHist(h3, f3, method);
  double ymin(0.0001);
  double h1max(h1->GetBinContent(h1->GetMaximumBin()));
  double h2max(h2->GetBinContent(h2->GetMaximumBin()));
  double h3max(h3?h3->GetBinContent(h3->GetMaximumBin()):-1);
  double hmax(h1max);
  int imax(1);
  if (h2max > h1max) {
    hmax = h2max;
    imax = 2;
  }
  if (h3max > h2max) {
    hmax = h3max;
    imax = 3;
  }
  hmax *= 1.2;
  if (verbose)  {
    cout << "hmax = " << hmax << " from imax = " << imax;
    if (1 == imax) cout << " bin " << h1->GetMaximumBin() << " with maximum " << h1->GetBinContent(h1->GetMaximumBin()) << endl;
    if (2 == imax) cout << " bin " << h2->GetMaximumBin() << " with maximum " << h2->GetBinContent(h2->GetMaximumBin()) << endl;
    if (3 == imax) cout << " bin " << h3->GetMaximumBin() << " with maximum " << h3->GetBinContent(h3->GetMaximumBin()) << endl;
  }
  if (loga) {
    gPad->SetLogy(1);
    hmax *= 2.;
    double hmin(h1->GetMinimum(ymin));
    if (verbose) cout << "hmin1 = " << hmin << endl;
    if (h2->GetMinimum(ymin) < hmin) {
      hmin = h2->GetMinimum(ymin);
      if (verbose) cout << "hmin2 = " << hmin << endl;
    }
    if (h3 && h3->GetMinimum(ymin) < hmin) {
      hmin = h3->GetMinimum(ymin);
      if (verbose) cout << "hmin3 = " << hmin << endl;
    }
    h1->SetMinimum(0.1*hmin);
    if (verbose) cout << "hmin = " << hmin << endl;
  } else {
    gPad->SetLogy(0);
    h1->SetMinimum(0.);
  }

  h1->SetMaximum(hmax);

  TH1* h0 = h1->DrawCopy("hist");
  setTitles(h1, h0);
  h2->DrawCopy("histsame");
  if (h3) h3->DrawCopy("histsame");

  if (legend) {
    newLegend(xleg, yleg, xleg+0.25, yleg+0.15);
    legg->SetTextSize(0.03);
    string text;
    text = fDS[f1]->fName.c_str();
    if (fDBX) {
      text = Form("%s: %4.3f#pm%4.3f, %4.3f", fDS[f1]->fName.c_str(), h1->GetMean(), h1->GetMeanError(), h1->GetRMS());
    }
    legg->AddEntry(h1, text.c_str(), "f");

    text = fDS[f2]->fName.c_str();
    if (fDBX) {
      text = Form("%s: %4.3f#pm%4.3f, %4.3f", fDS[f2]->fName.c_str(), h2->GetMean(), h2->GetMeanError(), h2->GetRMS());
    }
    legg->AddEntry(h2, text.c_str(), "f");

    if (h3) {
      text = fDS[f3]->fName.c_str();
      if (fDBX) {
	text = Form("%s: %4.3f#pm%4.3f, %4.3f", fDS[f3]->fName.c_str(), h3->GetMean(), h3->GetMeanError(), h3->GetRMS());
      }
      legg->AddEntry(h3, text.c_str(), "f");
    }

    legg->Draw();
  }


  if (verbose) cout << "==>plotClass: overlay(" << f1 << ", " << h1->GetName() << " integral= " << h1->Integral()
		    << ", " << f2 << ", " << h2->GetName() << " integral= " << h2->Integral()
		    << (h3? Form(", %s, %s, %f", f3.c_str(), h3->GetName(), h3->Integral()) : "")
		    << ")  log: " << loga << " legend = " << legend
		    << endl;
}

// ----------------------------------------------------------------------
void plotClass::overlay(string h1name, string f1, string h2name, string f2, string h3name, string f3, int method, bool loga,
			bool legend, double xleg, double yleg) {

  cout << h1name << " from " << f1 << " vs. " << h2name << " from " << f2 << " vs. " << h3name << " from " << f3 << endl;

  TH1D *h1 = fDS[f1]->getHist(Form("%s", h1name.c_str()), true);
  TH1D *h2 = fDS[f2]->getHist(Form("%s", h2name.c_str()), true);
  TH1D *h3(0);
  if (h3name != "") {
    h3 = fDS[f3]->getHist(Form("%s", h3name.c_str()), true);
  } else {
    h3 = 0;
  }

  overlay(h1, f1, h2, f2, h3, f3, method, loga, legend, xleg, yleg);
}


// ----------------------------------------------------------------------
void plotClass::loopFunction1() {
}


// ----------------------------------------------------------------------
void plotClass::loopFunction2() {
}

// ----------------------------------------------------------------------
void plotClass::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  if (2 == ifunc)          step = 10000;
  cout << "==> plotClass::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"  << " looping from  " << nbegin << " .. " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  void (plotClass::*pF)(void);
  if (ifunc == 1) pF = &plotClass::loopFunction1;
  if (ifunc == 2) pF = &plotClass::loopFunction2;

  cout << "pF: " << pF << endl;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotClass::setupTree(TTree *t, string mode) {

  if (string::npos != mode.find("Mc")) {
    fIsMC = true;
  } else {
    fIsMC = false;
  }

  t->SetBranchAddress("pt", &fb.pt);
  t->SetBranchAddress("q", &fb.q);

  t->SetBranchAddress("tis", &fb.tis);
  t->SetBranchAddress("reftrg", &fb.reftrg);

  t->SetBranchAddress("tau", &fb.tau);
  t->SetBranchAddress("gtau", &fb.gtau);

  t->SetBranchAddress("bdt",&fb.bdt);
  t->SetBranchAddress("lip",&fb.lip);
  t->SetBranchAddress("lipE",&fb.lipE);
  t->SetBranchAddress("tip",&fb.tip);
  t->SetBranchAddress("tipE",&fb.tipE);

  t->SetBranchAddress("closetrk",&fb.closetrk);
  t->SetBranchAddress("pvlip",   &fb.pvlip);
  t->SetBranchAddress("pvlips",  &fb.pvlips);
  t->SetBranchAddress("pv2lip",  &fb.pv2lip);
  t->SetBranchAddress("pv2lips", &fb.pv2lips);
  t->SetBranchAddress("maxdoca", &fb.maxdoca);
  t->SetBranchAddress("pvip",    &fb.pvip);
  t->SetBranchAddress("pvips",   &fb.pvips);
  t->SetBranchAddress("pvip3d",  &fb.pvip3d);
  t->SetBranchAddress("pvips3d", &fb.pvips3d);
  t->SetBranchAddress("pvw8",    &fb.pvw8);
  t->SetBranchAddress("pvz",     &fb.pvz);

  t->SetBranchAddress("m1pix",    &fb.m1pix);
  t->SetBranchAddress("m2pix",    &fb.m2pix);
  t->SetBranchAddress("m1bpix",   &fb.m1bpix);
  t->SetBranchAddress("m2bpix",   &fb.m2bpix);
  t->SetBranchAddress("m1bpixl1", &fb.m1bpixl1);
  t->SetBranchAddress("m2bpixl1", &fb.m2bpixl1);

  t->SetBranchAddress("rr",     &fb.rr);
  t->SetBranchAddress("pvn",    &fb.pvn);
  t->SetBranchAddress("run",    &fb.run);
  t->SetBranchAddress("l1s",    &fb.l1s);
  t->SetBranchAddress("evt",    &fb.evt);
  t->SetBranchAddress("hlt",    &fb.hlt);
  t->SetBranchAddress("hltm",   &fb.hltm);
  t->SetBranchAddress("ls",     &fb.ls);
  t->SetBranchAddress("ps",     &fb.ps);
  t->SetBranchAddress("chan",   &fb.chan);
  t->SetBranchAddress("cb",     &fb.cb);
  t->SetBranchAddress("json",   &fb.json);
  t->SetBranchAddress("gmuid",  &fb.gmuid);
  t->SetBranchAddress("gmutmid", &fb.gmutmid);
  t->SetBranchAddress("gmumvaid", &fb.gmumvaid);
  t->SetBranchAddress("gtqual", &fb.gtqual);
  t->SetBranchAddress("tm",     &fb.tm);
  t->SetBranchAddress("procid", &fb.procid);
  t->SetBranchAddress("m",      &fb.m);
  t->SetBranchAddress("m3",     &fb.m3);
  t->SetBranchAddress("m4",     &fb.m4);
  t->SetBranchAddress("me",     &fb.me);
  t->SetBranchAddress("cm",     &fb.cm);
  t->SetBranchAddress("pt",     &fb.pt);
  t->SetBranchAddress("phi",    &fb.phi);
  t->SetBranchAddress("eta",    &fb.eta);
  t->SetBranchAddress("cosa",   &fb.cosa);
  t->SetBranchAddress("alpha",  &fb.alpha);
  t->SetBranchAddress("iso",    &fb.iso);
  t->SetBranchAddress("chi2",   &fb.chi2);
  t->SetBranchAddress("dof",    &fb.dof);
  t->SetBranchAddress("prob",   &fb.pchi2dof);
  t->SetBranchAddress("chi2dof",&fb.chi2dof);
  t->SetBranchAddress("flsxy",  &fb.flsxy);
  t->SetBranchAddress("fls3d",  &fb.fls3d);
  t->SetBranchAddress("fl3d",   &fb.fl3d);
  t->SetBranchAddress("fl3dE",  &fb.fl3dE);
  t->SetBranchAddress("m1pt",   &fb.m1pt);
  t->SetBranchAddress("m1gt",   &fb.m1gt);
  t->SetBranchAddress("m1eta",  &fb.m1eta);
  t->SetBranchAddress("m1phi",  &fb.m1phi);
  t->SetBranchAddress("m1q",    &fb.m1q);
  t->SetBranchAddress("m2pt",   &fb.m2pt);
  t->SetBranchAddress("m2gt",   &fb.m2gt);
  t->SetBranchAddress("m2eta",  &fb.m2eta);
  t->SetBranchAddress("m2phi",  &fb.m2phi);
  t->SetBranchAddress("m2q",    &fb.m2q);
  t->SetBranchAddress("docatrk",&fb.docatrk);

  t->SetBranchAddress("m1id",     &fb.m1id);
  t->SetBranchAddress("m1rmvaid", &fb.m1rmvaid);
  t->SetBranchAddress("m1trigm",  &fb.m1trigm);
  t->SetBranchAddress("m1rmvabdt",&fb.m1rmvabdt);
  t->SetBranchAddress("m1tmid",   &fb.m1tmid);
  t->SetBranchAddress("m1gmid",   &fb.m1gmid);

  t->SetBranchAddress("m2id",     &fb.m2id);
  t->SetBranchAddress("m2rmvaid", &fb.m2rmvaid);
  t->SetBranchAddress("m2trigm",  &fb.m2trigm);
  t->SetBranchAddress("m2rmvabdt",&fb.m2rmvabdt);
  t->SetBranchAddress("m2tmid",   &fb.m2tmid);
  t->SetBranchAddress("m2gmid",   &fb.m2gmid);

  t->SetBranchAddress("m1iso",     &fb.m1iso);
  t->SetBranchAddress("m2iso",     &fb.m2iso);
  t->SetBranchAddress("closetrks1",&fb.closetrks1);
  t->SetBranchAddress("closetrks2",&fb.closetrks2);
  t->SetBranchAddress("closetrks3",&fb.closetrks3);
  t->SetBranchAddress("othervtx",  &fb.othervtx);
  t->SetBranchAddress("pvdchi2",   &fb.pvdchi2);

  t->SetBranchAddress("g1pt",   &fb.g1pt);
  t->SetBranchAddress("g2pt",   &fb.g2pt);
  t->SetBranchAddress("g1eta",  &fb.g1eta);
  t->SetBranchAddress("g2eta",  &fb.g2eta);
  t->SetBranchAddress("g1id",   &fb.g1id);
  t->SetBranchAddress("g2id",   &fb.g2id);
  if (string::npos != mode.find("psi")) {
    t->SetBranchAddress("mpsi", &fb.mpsi);
    t->SetBranchAddress("psipt", &fb.psipt);
    t->SetBranchAddress("psicosa", &fb.psicosa);
    t->SetBranchAddress("psiprob", &fb.psiprob);
    t->SetBranchAddress("psiflsxy", &fb.psiflsxy);
    t->SetBranchAddress("psimaxdoca", &fb.psimaxdoca);

  }
  if (string::npos != mode.find("bupsik")) {
    if (string::npos != mode.find("Mc")) {
      t->SetBranchAddress("g3pt", &fb.g3pt);
      t->SetBranchAddress("g3eta",&fb.g3eta);
    }
    t->SetBranchAddress("kpt",  &fb.k1pt);
    t->SetBranchAddress("kgt",  &fb.k1gt);
    t->SetBranchAddress("keta", &fb.k1eta);
  }

  if (string::npos != mode.find("bspsiphi")) {
    if (string::npos != mode.find("Mc")) {
      t->SetBranchAddress("g3pt", &fb.g3pt);
      t->SetBranchAddress("g3eta",&fb.g3eta);
      t->SetBranchAddress("g4pt", &fb.g4pt);
      t->SetBranchAddress("g4eta",&fb.g4eta);
    }
    t->SetBranchAddress("mpsi", &fb.mpsi);
    t->SetBranchAddress("psipt", &fb.psipt);
    t->SetBranchAddress("mkk",  &fb.mkk);
    t->SetBranchAddress("dr",   &fb.dr);
    t->SetBranchAddress("k1pt", &fb.k1pt);
    t->SetBranchAddress("k1gt", &fb.k1gt);
    t->SetBranchAddress("k1eta",&fb.k1eta);
    t->SetBranchAddress("k2pt", &fb.k2pt);
    t->SetBranchAddress("k2gt", &fb.k2gt);
    t->SetBranchAddress("k2eta",&fb.k2eta);
  } else {
    fb.mkk = 999.;
    fb.dr = 999.;
  }

  if (string::npos != mode.find("dstarpi")) {
    t->SetBranchAddress("md0",&fb.md0);
    t->SetBranchAddress("dm",&fb.dm);
    t->SetBranchAddress("ptd0",&fb.ptd0);
  }
}





// ----------------------------------------------------------------------
void plotClass::candAnalysis(/*int mode*/) {

  cuts *pCuts(0);
  fChan = detChan(fb.m1eta, fb.m2eta);
  if (fChan < 0) {
    if (0) cout << "plotClass::candAnalysis: " << fb.run << " " << fb.evt
		<< " could not determine channel: " << fb.m1eta << " " << fb.m2eta << endl;
    fBDT = -99.;
    fGoodHLT = fGoodMuonsID = false;
    fGoodQ = fGoodPvAveW8 = fGoodMaxDoca = fGoodIp = fGoodIpS = fGoodPt = fGoodEta = fGoodAlpha =  fGoodChi2 = fGoodFLS = false;
    fGoodCloseTrack = fGoodIso = fGoodDocaTrk = fGoodLastCut = fPreselection = false;
    fGoodAcceptance = fGoodBdtPt = fGoodMuonsPt = fGoodMuonsEta = fGoodTracks =  fGoodTracksPt = fGoodTracksEta = false;
    return;
  }
  pCuts = fCuts[fChan];

  bool bp2jpsikp(false), bs2jpsiphi(false);
  if (BU2JPSIKP == fMode)  bp2jpsikp = true;
  if (BS2JPSIPHI == fMode)  bs2jpsiphi = true;

  // -- reset all
  fBDT = -99.;
  fGoodHLT = fGoodMuonsID = false;
  fGoodQ = fGoodPvAveW8 = fGoodMaxDoca = fGoodIp = fGoodIpS = fGoodPt = fGoodEta = fGoodAlpha =  fGoodChi2 = fGoodFLS = false;
  fGoodCloseTrack = fGoodIso = fGoodDocaTrk = fGoodLastCut = fPreselection = false;

  fGoodJpsiCuts = true;

  fGoodAcceptance = true;
  fGoodBdtPt      = true;
  fGoodMuonsPt    = true;
  fGoodMuonsEta   = true;
  fGoodTracks     = fb.gtqual;
  fGoodTracksPt   = true;
  fGoodTracksEta  = true;

  fIsCowboy = fb.cb;

  if (fIsMC) {
    if (fb.g1pt < fAccPt) fGoodAcceptance = false;
    if (fb.g2pt < fAccPt) fGoodAcceptance = false;
    if (TMath::Abs(fb.g1eta) > 2.5) fGoodAcceptance = false;
    if (TMath::Abs(fb.g2eta) > 2.5) fGoodAcceptance = false;
  } else {
    static int runComplained(-1);
    if (!fb.json) {
      if (fb.run != runComplained) {
	cout << "json failure for run " << fb.run << " LS " << fb.ls << endl;
	runComplained = fb.run;
      }
      return;
    }
  }

  if (fb.m1pt < fAccPt) fGoodAcceptance = false;
  if (fb.m2pt < fAccPt) fGoodAcceptance = false;
  if (0 == fb.m1gt)  fGoodAcceptance = false;
  if (0 == fb.m2gt)  fGoodAcceptance = false;

  // if (fb.m1pt < pCuts->bdtPt) {
  //   fGoodBdtPt = false;
  // }
  // if (fb.m2pt < pCuts->bdtPt) {
  //   fGoodBdtPt = false;
  // }

  if (fb.m1pt < pCuts->m1pt) {
    fGoodMuonsPt = false;
  }
  if (fb.m2pt < pCuts->m2pt) {
    fGoodMuonsPt = false;
  }
  if (TMath::Abs(fb.m1eta) > 2.4) {
    fGoodAcceptance = false;
  }
  if (TMath::Abs(fb.m2eta) > 2.4) {
    fGoodAcceptance = false;
  }

  // -- by virtue of having a fChan this is true
  fGoodMuonsEta = true;

  if (bp2jpsikp) {
    if (fIsMC) {
      // gen-level cuts for Bu2JpsiKp
      if (fb.g1pt < fAccPt) fGoodAcceptance = false; // FIXME?
      if (fb.g2pt < fAccPt) fGoodAcceptance = false; // FIXME?
      if (TMath::Abs(fb.g3eta) > 2.5) fGoodAcceptance = false;
      if (fb.g3pt < 0.4) fGoodAcceptance = false;
    }
    if (TMath::Abs(fb.k1eta) > 2.4) {
      fGoodAcceptance = false;
      fGoodTracksEta = false;
    }
    if (fb.k1pt < 0.5) {
      fGoodAcceptance = false;
      fGoodTracksPt = false;
    }
    if (0 == fb.k1gt)  fGoodAcceptance = false;
  }

  if (bs2jpsiphi) {
    if (fIsMC) {
      if (TMath::Abs(fb.g3eta) > 2.5) fGoodAcceptance = false;
      if (TMath::Abs(fb.g4eta) > 2.5) fGoodAcceptance = false;
      // gen-level cuts for Bs2JpsiPhi
      if (fb.g1pt < fAccPt) fGoodAcceptance = false; // FIXME?
      if (fb.g2pt < fAccPt) fGoodAcceptance = false; // FIXME?
      if (fb.g3pt < 0.4) fGoodAcceptance = false;
      if (fb.g4pt < 0.4) fGoodAcceptance = false;
    }
    if (TMath::Abs(fb.k1eta) > 2.4) {
      fGoodAcceptance = false;
      fGoodTracksEta = false;
    }
    if (TMath::Abs(fb.k2eta) > 2.4) {
      fGoodAcceptance = false;
      fGoodTracksEta = false;
    }
    if (fb.k1pt < 0.5) {
      fGoodAcceptance = false;
      fGoodTracksPt = false;
    }
    if (fb.k2pt < 0.5) {
      fGoodAcceptance = false;
      fGoodTracksPt = false;
    }
    if (0 == fb.k1gt)  fGoodAcceptance = false;
    if (0 == fb.k2gt)  fGoodAcceptance = false;

    if (fb.dr   > 0.3) fGoodJpsiCuts = false;
    if (fb.mkk  < 1.01) fGoodJpsiCuts = false;
    if (fb.mkk  > 1.03) fGoodJpsiCuts = false;
  }

  if (bs2jpsiphi || bp2jpsikp) {
    if (fb.mpsi > 3.2) fGoodJpsiCuts = false;
    if (fb.mpsi < 3.0) fGoodJpsiCuts = false;
    if (fb.psipt < 7.0) fGoodJpsiCuts = false;
  } else {
    fGoodJpsiCuts = true;
  }

  if (fDoUseBDT) {
    if (fGoodAcceptance
        && fGoodTracks
        && fGoodTracksPt
        && fGoodTracksEta
        && fGoodBdtPt
        && fGoodMuonsEta
        && fGoodJpsiCuts
        ) {
      calcBDT();
      fb.bdt = fBDT;
    }
//     else {
//       cout << "acceptance:    " << fGoodAcceptance  << endl;
//       cout << "goodtracks:    " << fGoodTracks  << endl;
//       cout << "goodtrackspt:  " << fGoodTracksPt  << endl;
//       cout << "goodtrackseta: " << fGoodTracksEta  << endl;
//       cout << "goodBdtPt:     " << fGoodBdtPt  << endl;
//       cout << "goodMuonsEta:  " << fGoodMuonsEta  << endl;
//       cout << "goodJpsiCuts:  " << fGoodJpsiCuts  << endl;
//     }
  }

  fGoodMuonsID  = fb.m1gmid && fb.m2gmid;

  fW8 = 1.;
  fW8MmuID = fW8Mtrig = fW8DmuID = fW8Dtrig = -1.;
  double w1(-1.), w2(-1.);

  if (fIsMC) {
    PidTable *pT, *pT1, *pT2;

    // -- Weights with data PidTables
    // if (fIsCowboy) {
    //   pT  = fptCbM;
    //   pT1 = fptCbT1;
    //   pT2 = fptCbT2;
    // } else {
    //   pT  = fptSgM;
    //   pT1 = fptSgT1;
    //   pT2 = fptSgT2;
    // }

    pT  = fptM;
    pT1 = fptT1;
    pT2 = fptT2;


    double am1eta = TMath::Abs(fb.m1eta);
    double am2eta = TMath::Abs(fb.m2eta);

    w1       = pT->effD(fb.m1pt, am1eta, fb.m1phi);
    w2       = pT->effD(fb.m2pt, am2eta, fb.m2phi);
    fW8DmuID = w1*w2;

    w1       = pT1->effD(fb.m1pt, am1eta, fb.m1phi) * pT2->effD(fb.m1pt, am1eta, fb.m1phi);
    w2       = pT1->effD(fb.m2pt, am2eta, fb.m2phi) * pT2->effD(fb.m2pt, am2eta, fb.m2phi);
    fW8Dtrig = w1*w2;

    // -- Weights with MC PidTables
    // if (fIsCowboy) {
    //   pT  = fptCbMMC;
    //   pT1 = fptCbT1MC;
    //   pT2 = fptCbT2MC;
    // } else {
    //   pT  = fptSgMMC;
    //   pT1 = fptSgT1MC;
    //   pT2 = fptSgT2MC;
    // }
    pT  = fptM;
    pT1 = fptT1;
    pT2 = fptT2;

    w1       = pT->effD(fb.m1pt, am1eta, fb.m1phi);
    w2       = pT->effD(fb.m2pt, am2eta, fb.m2phi);
    fW8MmuID = w1*w2;

    w1       = pT1->effD(fb.m1pt, am1eta, fb.m1phi) * pT2->effD(fb.m1pt, am1eta, fb.m1phi);
    w2       = pT1->effD(fb.m2pt, am2eta, fb.m2phi) * pT2->effD(fb.m2pt, am2eta, fb.m2phi);
    fW8Mtrig = w1*w2;

    if (RARE == fMode) {
      w1 = w2 = -1.;
      // -- track 1
      if (  13 == fb.g1id) w1 = pT->effD(fb.m1pt, am1eta, 1.);
      if ( -13 == fb.g1id) w1 = pT->effD(fb.m1pt, am1eta, 1.);

      if (  321 == fb.g1id) w1 = fptFakePosKaons->effD(fb.m1pt, am1eta, 1.);
      if ( -321 == fb.g1id) w1 = fptFakeNegKaons->effD(fb.m1pt, am1eta, 1.);

      if (  211 == fb.g1id) w1 = fptFakePosPions->effD(fb.m1pt, am1eta, 1.);
      if ( -211 == fb.g1id) w1 = fptFakeNegPions->effD(fb.m1pt, am1eta, 1.);

      if ( 2212 == fb.g1id) w1 = fptFakePosProtons->effD(fb.m1pt, am1eta, 1.);
      if (-2212 == fb.g1id) w1 = fptFakeNegProtons->effD(fb.m1pt, am1eta, 1.);

      // -- track 2
      if (  13 == fb.g2id) w2 = pT->effD(fb.m2pt, am2eta, 1.);
      if ( -13 == fb.g2id) w2 = pT->effD(fb.m2pt, am2eta, 1.);

      if (  321 == fb.g2id) w2 = fptFakePosKaons->effD(fb.m2pt, am2eta, 1.);
      if ( -321 == fb.g2id) w2 = fptFakeNegKaons->effD(fb.m2pt, am2eta, 1.);

      if (  211 == fb.g2id) w2 = fptFakePosPions->effD(fb.m2pt, am2eta, 1.);
      if ( -211 == fb.g2id) w2 = fptFakeNegPions->effD(fb.m2pt, am2eta, 1.);

      if ( 2212 == fb.g2id) w2 = fptFakePosProtons->effD(fb.m2pt, am2eta, 1.);
      if (-2212 == fb.g2id) w2 = fptFakeNegProtons->effD(fb.m2pt, am2eta, 1.);

      fW8MisId = w1*w2;
    }

  }

  fGoodQ          = (fb.m1q*fb.m2q < 0);
  fGoodPvAveW8    = (fb.pvw8 > 0.7);
  fGoodMaxDoca    = (TMath::Abs(fb.maxdoca) < pCuts->maxdoca);
  fGoodIp         = (TMath::Abs(fb.pvip) < pCuts->pvip);
  fGoodIpS        = (TMath::Abs(fb.pvips) < pCuts->pvips);

  fGoodLip        = (TMath::Abs(fb.pvlip) < pCuts->pvlip);
  fGoodLipS       = (TMath::Abs(fb.pvlips) < pCuts->pvlips);

  fGoodPt         = (fb.pt > pCuts->pt);
  fGoodEta        = ((fb.eta > -2.40) && (fb.eta < 2.40));
  fGoodAlpha      = (fb.alpha < pCuts->alpha);
  fGoodChi2       = (fb.chi2/fb.dof < pCuts->chi2dof);
  fGoodFLS        = (fb.fls3d > pCuts->fls3d);
  if (TMath::IsNaN(fb.fls3d)) fGoodFLS = false;

  fGoodCloseTrack = (fb.closetrk < pCuts->closetrk);
  fGoodIso        = (fb.iso > pCuts->iso);
  fGoodDocaTrk    = (fb.docatrk > pCuts->docatrk);
  fGoodLastCut    = true;

  //FIXME  fGoodBDT        = (fBDT > pCuts->bdtMin);
  fGoodBDT        = true;
  //FIXME  fGoodHLT        = fb.hlt &&fb.hltm;
  fGoodHLT        = fb.hlt;

  // -- no trigger matching for rare decays!
  if (RARE == fMode) fGoodHLT = fb.hlt;

  fPreselection   = (fGoodHLT && fGoodMuonsID && fGoodMuonsPt && fGoodMuonsEta && (fb.alpha < 0.2) && (fb.fls3d > 5));
  if (bs2jpsiphi || bp2jpsikp) {
    fPreselection = fPreselection && fGoodJpsiCuts;

    if (!fGoodJpsiCuts) {
      //      cout << "dr: " << fb.dr  << " mkk: " << fb.mkk << " mpsi: " << fb.mpsi << " psipt: " << fb.psipt << endl;
    }
  }

  if (!fPreselection) {
    if (0)     cout << "HLT: " << fGoodHLT << " muons: " << fGoodMuonsID << " muonspt: " << fGoodMuonsPt << " muonseta: " << fGoodMuonsEta
		    << " a: " << fb.alpha << " fls3d: " << fb.fls3d
		    << " mpsi: " << fb.mpsi << " psipt: " << fb.psipt
		    << " dr: " << fb.dr  << " mkk: " << fb.mkk
		    << endl;
  }

  fAnaCuts.update();

}



// ----------------------------------------------------------------------
TTree* plotClass::getTree(string ds, string dir, string tree) {
  if (!fDS[ds]) {
    cout << "xx> plotClass::getTree: dataset " << ds << " not found" << endl;
    return 0;
  }
  TTree *t(0);
  if (!dir.compare("")) {
    t = (TTree*)fDS[ds]->fF->Get(tree.c_str());
  } else {
    t = (TTree*)fDS[ds]->fF->Get(Form("%s/%s", dir.c_str(), tree.c_str()));
  }
  cout << "plotClass::getTree(" << ds << ", " << dir << ", " << tree << "): " << t << endl;
  return t;
}

// ----------------------------------------------------------------------
TFile* plotClass::loadFile(string file) {
  TFile *f = TFile::Open(file.c_str());
  return f;
}



// ----------------------------------------------------------------------
void plotClass::replaceAll(string &sInput, const string &oldString, const string &newString) {
  string::size_type foundpos = sInput.find(oldString);
  while (foundpos != string::npos)  {
    sInput.replace(sInput.begin() + foundpos, sInput.begin() + foundpos + oldString.length(), newString);
    foundpos = sInput.find(oldString);
  }
}

// ----------------------------------------------------------------------
void plotClass::newLegend(double x1, double y1, double x2, double y2, string title) {
  //  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2, title.c_str());
  legg->SetFillStyle(0);
  legg->SetBorderSize(0);
  legg->SetTextSize(0.04);
  legg->SetFillColor(0);
  legg->SetTextFont(52);
}

// ----------------------------------------------------------------------
void plotClass::makeCanvas(int i) {
  if (i & 16) {
    c5 = new TCanvas("c5", "c5", 210,   0, 800, 900);
    c5->ToggleEventStatus();
  }
  if (i & 8) {
    c4 = new TCanvas("c4", "c4", 210,   0, 800, 600);
    c4->ToggleEventStatus();
  }
  if (i & 4) {
    c3 = new TCanvas("c3", "c3", 200,  20, 800, 800);
    c3->ToggleEventStatus();
  }
  if (i & 1) {
    //    c1 = new TCanvas("c1", "c1", 20,  60, 1200, 400);
    c1 = new TCanvas("c1", "c1", 20,  60, 1000, 400);
    c1->ToggleEventStatus();
  }
  if (i & 2) {
    c2 = new TCanvas("c2", "c2", 300, 200, 400, 800);
    c2->ToggleEventStatus();
  }
}


// ----------------------------------------------------------------------
void plotClass::calcBDT() {
  fBDT = -99.;

  if (!preselection(fb, fChan)) return;

  //??  if (5 == mode && 5.2 < mass && mass < 5.45 && fb.iso < 0.7) continue;
  //  if (rejectInvIso && 5.2 < fb.m && fb.m < 5.45 && fb.iso < 0.7) return;
  //   if (fb.pt > 100) return;
  //   if (fb.pt < 6) return;
  //   if (fb.m1pt < 4) return;
  //   if (fb.m2pt < 4) return;
  //   if (fb.fl3d > 1.5) return;
  //   if (fb.m > 5.9) return;
  //   if (fb.m < 4.9) return;

  //   if (!fb.hlt) return;
  //   if (!fb.gmuid) return;

  frd.pt = fb.pt;
  frd.eta = fb.eta;
  frd.m1eta = fb.m1eta;
  frd.m2eta = fb.m2eta;
  frd.m1pt = fb.m1pt;
  frd.m2pt = fb.m2pt;
  frd.fls3d = fb.fls3d;
  frd.alpha = fb.alpha;
  frd.maxdoca = fb.maxdoca;
  frd.pvip = fb.pvip;
  frd.pvips = fb.pvips;
  frd.iso = fb.iso;
  frd.docatrk = fb.docatrk;
  frd.chi2dof = fb.chi2dof;
  frd.closetrk = fb.closetrk;

  frd.m1iso = fb.m1iso;
  frd.m2iso = fb.m2iso;

  frd.closetrks1 = fb.closetrks1;
  frd.closetrks2 = fb.closetrks2;
  frd.closetrks3 = fb.closetrks3;

  frd.pv2lip  = fb.pv2lip;
  frd.pv2lips = fb.pv2lips;

  frd.m  = fb.m;
  int remainder = TMath::Abs(fb.evt%3);
  if (0 == remainder) {
    fBDT   = fReaderEvents0[fChan]->EvaluateMVA("BDT");
  } else if (1 == remainder) {
    fBDT   = fReaderEvents1[fChan]->EvaluateMVA("BDT");
  } else if (2 == remainder) {
    fBDT   = fReaderEvents2[fChan]->EvaluateMVA("BDT");
  } else {
    cout << "all hell break loose" << endl;
  }
}


// ----------------------------------------------------------------------
TMVA::Reader* plotClass::setupReader(string xmlFile, readerData &rd) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );

  TString dir    = "weights/";
  TString methodNameprefix = "BDT";
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
  cout << "adding variables: ";
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\"");
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2);
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
        m1 = allLines[j].find("Expression=\"")+10;
        m2 = allLines[j].find("\" Label=\"");
        stype = allLines[j].substr(m1+2, m2-m1-2);
        //      cout << "ivar " << j-i << " variable string: ->" << stype << "<-" << endl;
        if (stype == "m1pt") {
          cout << " m1pt";
          reader->AddVariable( "m1pt", &rd.m1pt);
        }
        if (stype == "m2pt") {
          cout << " m2pt";
          reader->AddVariable( "m2pt", &rd.m2pt);
        }
        if (stype == "m1eta") {
          cout << " m1eta";
          reader->AddVariable( "m1eta", &rd.m1eta);
        }
        if (stype == "m2eta") {
          reader->AddVariable( "m2eta", &rd.m2eta);
          cout << " m2eta";
        }
        if (stype == "pt") {
          cout << " pt";
          reader->AddVariable( "pt", &rd.pt);
        }
        if (stype == "eta") {
          cout << " eta";
          reader->AddVariable( "eta", &rd.eta);
        }
        if (stype == "fls3d") {
          cout << " fls3d";
          reader->AddVariable( "fls3d", &rd.fls3d);
        }
        if (stype == "alpha") {
          cout << " alpha";
          reader->AddVariable( "alpha", &rd.alpha);
        }
        if (stype == "maxdoca") {
          cout << " maxdoca";
          reader->AddVariable( "maxdoca", &rd.maxdoca);
        }
        if (stype == "pvip") {
          cout << " pvip";
          reader->AddVariable( "pvip", &rd.pvip);
        }
        if (stype == "pvips") {
          cout << " pvips";
          reader->AddVariable( "pvips", &rd.pvips);
        }
        if (stype == "iso") {
          cout << " iso";
          reader->AddVariable( "iso", &rd.iso);
        }
        if (stype == "docatrk") {
          cout << " docatrk";
          reader->AddVariable( "docatrk", &rd.docatrk);
        }
        if (stype == "closetrk") {
          cout << " closetrk";
          reader->AddVariable( "closetrk", &rd.closetrk);
        }
        if (stype == "chi2dof") {
          cout << " chi2dof";
          reader->AddVariable( "chi2dof", &rd.chi2dof);
        }
        if (stype == "closetrks1") {
          cout << " closetrks1";
          reader->AddVariable( "closetrks1", &rd.closetrks1);
        }
        if (stype == "closetrks2") {
          cout << " closetrks2";
          reader->AddVariable( "closetrks2", &rd.closetrks2);
        }
        if (stype == "closetrks3") {
          cout << " closetrks3";
          reader->AddVariable( "closetrks3", &rd.closetrks3);
        }
        if (stype == "m1iso") {
          cout << " m1iso";
          reader->AddVariable( "m1iso", &rd.m1iso);
        }
        if (stype == "m2iso") {
          cout << " m2iso";
          reader->AddVariable( "m2iso", &rd.m2iso);
        }
        if (stype == "othervtx") {
          cout << " othervtx";
          reader->AddVariable( "othervtx", &rd.othervtx);
        }
        if (stype == "pvdchi2") {
          cout << " pvdchi2";
          reader->AddVariable( "pvdchi2", &rd.pvdchi2);
        }
        if (stype == "pv2lip") {
          cout << " pv2lip";
          reader->AddVariable( "pv2lip", &rd.pv2lip);
        }
        if (stype == "pv2lips") {
          cout << " pv2lips";
          reader->AddVariable( "pv2lips", &rd.pv2lips);
        }
      }
      break;
    }
  }
  cout << endl;

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
	//        cout << "ivar " << j-i << " spectator string: ->" << stype << "<-" << endl;
        if (stype == "m") {
	  // cout << "  adding m as spectator" << endl;
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
void plotClass::readFile(string filename, vector<string> &lines) {
  cout << "    readFile " << filename << endl;
  char  buffer[200];
  ifstream is(filename.c_str());
  if (!is) {
    cout << "file ->" << filename << "<- not found, exit(1)" << endl;
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
void plotClass::readCuts(string filename) {
  cout << "==> plotClass: Reading " << filename << " for cut settings" << endl;
  vector<string> cutLines;
  readFile(filename, cutLines);

  float cutvalue;
  int dump(0), ok(0);
  string cutname("nada");

  cuts *a = 0;

  fCuts.clear();

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
    cleanupString(cutLines[i]);
    vector<string> lineItems = split(cutLines[i], ' ');
    if (lineItems.size() < 2) {
      //      cout << "no complete cut line provided ->" << cutLines[i] << "<-" << endl;
      continue;
    }
    cutname  = lineItems[0];

    for (unsigned int j = 1; j < lineItems.size(); ++j) {
      a = fCuts[j-1];

      cutvalue = atof(lineItems[j].c_str());
      if (cutname == "nchan") {
	ok = 1;
      }

      if (cutname == "index") {
	a->index = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "index:                " << cutvalue << endl;
      }

      if (cutname == "metaMin") {
	a->metaMin = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "metaMin:              " << cutvalue << endl;
      }

      if (cutname == "metaMax") {
	a->metaMax = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "metaMax:              " << cutvalue << endl;
      }

      if (cutname == "mBdLo") {
	a->mBdLo = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "mBdLo:                " << cutvalue << endl;
      }

      if (cutname == "mBdLo") {
	a->mBdLo = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "mBdLo:                " << cutvalue << endl;
      }

      if (cutname == "mBdHi") {
	a->mBdHi = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "mBdHi:                " << cutvalue << endl;
      }

      if (cutname == "mBsLo") {
	a->mBsLo = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "mBsLo:                " << cutvalue << endl;
      }

      if (cutname == "mBsHi") {
	a->mBsHi = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "mBsHi:                " << cutvalue << endl;
      }

      if (cutname == "etaMin") {
	a->etaMin = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "etaMin:               " << cutvalue << endl;
      }

      if (cutname == "etaMax") {
	a->etaMax = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "etaMax:               " << cutvalue << endl;
      }

      if (cutname == "pt") {
	a->pt = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "pt:                   " << cutvalue << endl;
      }

      if (cutname == "m1pt") {
	a->m1pt = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "m1pt:                 " << cutvalue << endl;
      }

      if (cutname == "m2pt") {
	a->m2pt = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "m2pt:                 " << cutvalue << endl;
      }

      if (cutname == "iso") {
	a->iso = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "iso:                   " << cutvalue << endl;
      }

      if (cutname == "chi2dof") {
	a->chi2dof = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "chi2dof:               " << cutvalue << endl;
      }

      if (cutname == "alpha") {
	a->alpha = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "alpha:                 " << cutvalue << endl;
      }

      if (cutname == "fls3d") {
	a->fls3d = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "fls3d:                 " << cutvalue << endl;
      }

      if (cutname == "docatrk") {
	a->docatrk = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "docatrk:               " << cutvalue << endl;
      }

      if (cutname == "maxdoca") {
	a->maxdoca = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "maxdoca:               " << cutvalue << endl;
      }

      if (cutname == "closetrk") {
	a->closetrk = static_cast<int>(cutvalue); ok = 1;
	if (dump) cout << j-1 << " " << "closetrk:              " << cutvalue << endl;
      }

      if (cutname == "pvip") {
	a->pvip = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "pvip:                  " << cutvalue << endl;
      }

      if (cutname == "pvips") {
	a->pvips = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "pvips:                 " << cutvalue << endl;
      }


      if (cutname == "pvlip") {
	a->pvlip = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "pvlip:                 " << cutvalue << endl;
      }

      if (cutname == "pvlips") {
	a->pvlips = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "pvlips:                " << cutvalue << endl;
      }

      if (cutname == "pv2lip") {
	a->pv2lip = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "pv2lip:                " << cutvalue << endl;
      }

      if (cutname == "pv2lips") {
	a->pv2lips = cutvalue; ok = 1;
	if (dump) cout << j-1 << " " << "pv2lips:               " << cutvalue << endl;
      }

      if (cutname == "l1seeds") {
	vector<string> vl1seeds = split(lineItems[j], ',');
	for (unsigned int is = 0; is < vl1seeds.size(); ++is) {
	  fCuts[j-1]->l1seeds.push_back(atoi(vl1seeds[is].c_str()));
	}
      }


    }

    if (!ok) cout << "==> what about " << cutname << endl;
  }

  cout << "==> finished reading cut setting, fCuts.size() =  " << fCuts.size() << endl;

}

// ----------------------------------------------------------------------
int plotClass::detChan(double m1eta, double m2eta) {

  double m1 = TMath::Abs(m1eta);
  double m2 = TMath::Abs(m2eta);

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
}


// ----------------------------------------------------------------------
void plotClass::printCuts(ostream &OUT) {

  OUT << "----------------------------------------------------------------------" << endl;
  OUT << "channel    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  OUT << Form("%10d", fCuts[i]->index);
  OUT << endl;

  OUT << "etaMin     ";
  fTEX << Form("\\vdef{%s:etaB:var}  {\\ensuremath{{|\\eta| } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->etaMin);
    fTEX <<  Form("\\vdef{%s:etaMin:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->etaMin) << endl;
  }
  OUT << endl;

  OUT << "etaMax     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->etaMax);
    fTEX <<  Form("\\vdef{%s:etaMax:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->etaMax) << endl;
  }
  OUT << endl;

  OUT << "l1seeds     ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    for (unsigned is = 0; is < fCuts[i]->l1seeds.size(); ++is) OUT << Form("%d ", fCuts[i]->l1seeds[is]);
    for (int is = fCuts[i]->l1seeds.size()*2; is < 10; ++is) OUT << " ";
    for (unsigned is = 0; is < fCuts[i]->l1seeds.size(); ++is)
      fTEX <<  Form("\\vdef{%s:l1seeds:%d}   {\\ensuremath{{%d } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->l1seeds[is]) << endl;
  }
  OUT << endl;

  OUT << "mBdLo      ";
  fTEX << Form("\\vdef{%s:mBd:var}  {\\ensuremath{{m(\\Bz) } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    OUT << Form("%10.3f", fCuts[i]->mBdLo);
    fTEX << Form("\\vdef{%s:mBdLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->mBdLo) << endl;
  }
  OUT << endl;

  OUT << "mBdHi      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    OUT << Form("%10.3f", fCuts[i]->mBdHi);
    fTEX << Form("\\vdef{%s:mBdHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->mBdHi) << endl;
  }
  OUT << endl;

  OUT << "mBsLo      ";
  fTEX << Form("\\vdef{%s:mBs:var}  {\\ensuremath{{m(\\Bs) } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->mBsLo);
    fTEX <<  Form("\\vdef{%s:mBsLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->mBsLo) << endl;
  }
  OUT << endl;

  OUT << "mBsHi      ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->mBsHi);
    fTEX <<  Form("\\vdef{%s:mBsHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->mBsHi) << endl;
  }
  OUT << endl;

  OUT << "pt         ";
  fTEX << Form("\\vdef{%s:ptb:var}  {\\ensuremath{{\\ptb } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pt);
    fTEX <<  Form("\\vdef{%s:pt:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->pt) << endl;
  }
  OUT << endl;

  OUT << "m1pt       ";
  fTEX << Form("\\vdef{%s:ptmuone:var}  {\\ensuremath{{\\ptmuone } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->m1pt);
    fTEX <<  Form("\\vdef{%s:m1pt:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->m1pt) << endl;
  }
  OUT << endl;

  OUT << "m2pt       ";
  fTEX << Form("\\vdef{%s:ptmutwo:var}  {\\ensuremath{{\\ptmutwo } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->m2pt);
    fTEX <<  Form("\\vdef{%s:m2pt:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->m2pt) << endl;
  }
  OUT << endl;

  OUT << "metaMin    ";
  fTEX << Form("\\vdef{%s:etamuf:var}  {\\ensuremath{{|\\eta(\\mu_f)| } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->metaMin);
    fTEX <<  Form("\\vdef{%s:metaMin:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->metaMin) << endl;
  }
  OUT << endl;

  OUT << "metaMax    ";
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->metaMax);
    fTEX <<  Form("\\vdef{%s:metaMax:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->metaMax) << endl;
  }
  OUT << endl;

  OUT << "iso        ";
  fTEX << Form("\\vdef{%s:iso:var}  {\\ensuremath{{\\iso } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->iso);
    fTEX <<  Form("\\vdef{%s:iso:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->iso) << endl;
  }
  OUT << endl;

  OUT << "chi2dof    ";
  fTEX << Form("\\vdef{%s:chidof:var}  {\\ensuremath{{\\chidof } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->chi2dof);
    fTEX <<  Form("\\vdef{%s:chi2dof:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->chi2dof) << endl;
  }
  OUT << endl;

  OUT << "alpha      ";
  fTEX << Form("\\vdef{%s:alpha:var}  {\\ensuremath{{\\alpha } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->alpha);
    fTEX <<  Form("\\vdef{%s:alpha:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->alpha) << endl;
  }
  OUT << endl;

  OUT << "fls3d      ";
  fTEX << Form("\\vdef{%s:fls:var}  {\\ensuremath{{\\fls } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->fls3d);
    fTEX <<  Form("\\vdef{%s:fls3d:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->fls3d) << endl;
  }
  OUT << endl;

  OUT << "docatrk    ";
  fTEX << Form("\\vdef{%s:docatrk:var}  {\\ensuremath{{\\docatrk } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->docatrk);
    fTEX <<  Form("\\vdef{%s:docatrk:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->docatrk) << endl;
  }
  OUT << endl;

  OUT << "closetrk   ";
  fTEX << Form("\\vdef{%s:closetrk:var}  {\\ensuremath{{\\closetrk } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->closetrk);
    fTEX <<  Form("\\vdef{%s:closetrk:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->closetrk) << endl;
  }
  OUT << endl;

  OUT << "maxdoca    ";
  fTEX << Form("\\vdef{%s:maxdoca:var}  {\\ensuremath{{\\maxdoca } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->maxdoca);
    fTEX <<  Form("\\vdef{%s:maxdoca:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->maxdoca) << endl;
  }
  OUT << endl;

  OUT << "pvip       ";
  fTEX << Form("\\vdef{%s:ip:var}  {\\ensuremath{{\\pvip } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pvip);
    fTEX <<  Form("\\vdef{%s:pvip:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->pvip) << endl;
  }
  OUT << endl;

  OUT << "pvips      ";
  fTEX << Form("\\vdef{%s:ips:var}  {\\ensuremath{{\\pvips } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pvips);
    fTEX <<  Form("\\vdef{%s:pvips:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->pvips) << endl;
  }
  OUT << endl;

  OUT << "pvlip      ";
  fTEX << Form("\\vdef{%s:lip:var}  {\\ensuremath{{\\pvlip } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pvlip);
    fTEX <<  Form("\\vdef{%s:pvlip:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->pvlip) << endl;
  }
  OUT << endl;

  OUT << "pvlips     ";
  fTEX << Form("\\vdef{%s:lips:var}  {\\ensuremath{{\\pvlips } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pvlips);
    fTEX <<  Form("\\vdef{%s:pvlips:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->pvlips) << endl;
  }
  OUT << endl;

  OUT << "pv2lip     ";
  fTEX << Form("\\vdef{%s:liptwo:var}  {\\ensuremath{{\\pvliptwo } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pv2lip);
    fTEX <<  Form("\\vdef{%s:pv2lip:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->pv2lip) << endl;
  }
  OUT << endl;

  OUT << "pv2lips    ";
  fTEX << Form("\\vdef{%s:lipstwo:var}  {\\ensuremath{{\\pvlipstwo } } }", fSuffix.c_str()) << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i)  {
    OUT << Form("%10.3f", fCuts[i]->pv2lips);
    fTEX <<  Form("\\vdef{%s:pv2lips:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), fCuts[i]->index, fCuts[i]->pv2lips) << endl;
  }
  OUT << endl;

  OUT.flush();

  return;
}



// ----------------------------------------------------------------------
void plotClass::setItalic() {
  tl->SetTextFont(52);
}


// ----------------------------------------------------------------------
void plotClass::setRoman() {
  tl->SetTextFont(42);
}


// ----------------------------------------------------------------------
void plotClass::savePad(string name, TCanvas *c) {
  if (0 == c) {
    gPad->SaveAs(Form("%s/%s", fDirectory.c_str(), name.c_str()));
  } else {
    c->SaveAs(Form("%s/%s", fDirectory.c_str(), name.c_str()));
  }
}



// ----------------------------------------------------------------------
void plotClass::insertDataset(std::string dsname, dataset *ds) {
  if (fDS.find(dsname) != fDS.end()) {
    cout << "######## Error: " << dsname  << " already present in fDS, NOT inserting again" << endl;
  } else {
    cout << "     inserting: " << dsname  << " " << endl;
    fDS.insert(make_pair(dsname, ds));
  }
}


// ----------------------------------------------------------------------
void plotClass::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotClass::loadFile loading files listed in " << files << endl;

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
    string::size_type m2 = sbuffer.find("file=");
    string stype("nada");
    bool useBf(false);
    cout << "sbuffer: " << sbuffer  << " m1: " << m1 << endl;
    if (m1 > sbuffer.size()) {
      m1 = sbuffer.find("bf=");
      useBf = true;
    }
    stype = sbuffer.substr(5, m1-5);
    string slumi("nada"), sbf("nada");
    if (useBf) {
      sbf = sbuffer.substr(m1+3, m2-m1-3);
    } else {
      slumi = sbuffer.substr(m1+5, m2-m1-5);
    }
    string sfile = sbuffer.substr(m2+5);

    string sname("nada"), sdecay("nada"), ldecay("");

    double bf(0.), bfE(0.);
    if (useBf) {
      cout << "sbf = " << sbf  << endl;
      float val, err;
      int expo;
      sscanf(sbf.c_str(), "(%f,%f)e%d", &val, &err, &expo);
      bf = val*TMath::Power(10., expo);
      bfE = err*TMath::Power(10., expo);
    } else {
      bf = 0.;
      bfE = 0.;
    }

    //    if (useBf) cout << " -> BF = " << bf << " +/- " << bfE << endl;

    TFile *pF(0);
    dataset *ds(0);

    if (string::npos != stype.find("data")) {
      // -- DATA
      pF = loadFile(sfile);

      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("bmm,")) {
        sname = "bmmData";
        sdecay = "dimuon";
	ldecay = "dimuon";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	ds->fLumi   = atof(slumi.c_str());
      }

      if (string::npos != stype.find("bupsik,")) {
        sname = "bupsikData";
        sdecay = "B^{+} #rightarrow J/#kern[-0.2]{#it{#psi}}K^{+}";
        ldecay = "\\bupsik";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	ds->fLumi   = atof(slumi.c_str());
      }

      if (string::npos != stype.find("bspsiphi,")) {
        sname = "bspsiphiData";
        sdecay = "B^{0}_{s} #rightarrow J/#kern[-0.2]{#it{#psi}}#it{#phi}";
	ldecay = "\\bspsiphi";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	ds->fLumi   = atof(slumi.c_str());
      }

      if (string::npos != stype.find("bdpsikstar,")) {
        sname = "bdpsikstarData";
        sdecay = "B^{0} #rightarrow J/#kern[-0.2]{#it{#psi}}K^{*0}";
	ldecay = "\\bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	ds->fLumi   = atof(slumi.c_str());
      }

    } else if (string::npos != stype.find("mc")) {
      // -- MC
      pF = loadFile(sfile);

      ds = new dataset();
      ds->fSize = 0.1;
      ds->fWidth = 2.;

      if (string::npos != stype.find("bupsik,")) {
        sname = "bupsikMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
	sdecay = "B^{+} #rightarrow J/#kern[-0.2]{#it{#psi}}K^{+}";
        ldecay = "\\bupsik";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

      if (string::npos != stype.find("bspsiphi,")) {
        sname = "bspsiphiMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
        sdecay = "B^{0}_{s} #rightarrow J/#kern[-0.2]{#it{#psi}}#it{#phi}";
        ldecay = "\\bspsiphi";
	ds->fColor = kRed;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bsmm,")) {
        sname = "bsmmMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
        sdecay = "B^{0}_{s} #rightarrow #it{#mu#mu}";
        ldecay = "\\bsmm";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }


      if (string::npos != stype.find("bdpsikstar,")) {
        sname = "bdpsikstarMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
        sdecay = "B^{0} #rightarrow J/#kern[-0.2]{#it{#psi}}K^{*0}";
        ldecay = "\\bdpsikstar";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bdmm,")) {
        sname = "bdmmMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
        sdecay = "B^{0} #rightarrow #it{#mu#mu}";
        ldecay = "\\bdmm";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bspipi,")) {
        sname = "bspipiMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
	if (string::npos != stype.find("bg")) sname += "Bg";
        sdecay = "B^{s} #rightarrow #it{#pi#pi}";
        ldecay = "\\bspipi";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bskk,")) {
        sname = "bskkMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
	if (string::npos != stype.find("bg")) sname += "Bg";
        sdecay = "B^{s} #rightarrow KK";
        ldecay = "\\bskk";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bskmunu,")) {
        sname = "bskmunuMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
	if (string::npos != stype.find("bg")) sname += "Bg";
        sdecay = "B^{s} #rightarrow K#it{#mu}#it{#nu}";
        ldecay = "\\bskmunu";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("lbppi,")) {
        sname = "lbppiMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
	if (string::npos != stype.find("bg")) sname += "Bg";
        sdecay = "#it{#Lambda}_{b} #rightarrow p #it{#pi}";
        ldecay = "\\lbppi";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("lbpmunu,")) {
	sname = "lbpmunuMc";
	if (string::npos != stype.find("mcOff")) sname += "Off";
	if (string::npos != stype.find("acc")) sname += "Acc";
	if (string::npos != stype.find("bg")) sname += "Bg";
        sdecay = "#it{#Lambda}_{b} #rightarrow p #it{#mu}  #it{#nu}";
        ldecay = "\\lbpmunu";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = bf;
	ds->fBfE    = bfE;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }


    }
    // -- insert if with a "valid" sname
    if (sname != "nada") {
      ds->fLcolor    = ds->fColor;
      ds->fFcolor    = ds->fColor;
      ds->fName      = sdecay;
      ds->fLatexName = ldecay;
      ds->fFullName  = sname;
      insertDataset(sname, ds);
    } else {
      delete ds;
    }


  }

  is.close();
}


// ----------------------------------------------------------------------
// downloaded on 2016/06/6 from https://ghm.web.cern.ch/ghm/plots/
TStyle * plotClass::setTdrStyle() {

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  // For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

  // For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

  // For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  // For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  // tdrStyle->SetErrorMarker(20);
  //tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

  //For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

  //For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

  // For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(52);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

  // Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

  // For the Global title:
  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(52);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

  // For the axis titles:
  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(52, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

  // For the axis labels:
  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(52, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

  // For the axis:
  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  // Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  // Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->SetHatchesLineWidth(1);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->cd();

  return tdrStyle;
}
