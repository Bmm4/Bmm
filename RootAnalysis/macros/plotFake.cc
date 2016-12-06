#include "plotFake.hh"

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
#include "TPad.h"
#include "TF1.h"
#include "TFitResult.h"

#include "common/dataset.hh"
#include "common/util.hh"

ClassImp(plotFake)

using namespace std;

// ----------------------------------------------------------------------
plotFake::plotFake(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotFake::loadFiles(files);

  changeSetup(dir, "plotFake", setup);

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

  fChannelList.clear();
  for (unsigned int i = 0; i < 3; ++i) {
    fChannelList.push_back(Form("%d", i));
  }
  fNchan = fChannelList.size();


  fDoList.clear();
  fDoList.push_back("FakePt");
  fDoList.push_back("FakeEta");

  fDoList.push_back("AllPt");
  fDoList.push_back("AllEta");

  fDoList.push_back("FakeTisFakePt");
  fDoList.push_back("FakeTisFakeEta");
  fDoList.push_back("FakeTisAllPt");
  fDoList.push_back("FakeTisAllEta");

  fDoList.push_back("FakeTisDtFakePt");
  fDoList.push_back("FakeTisDtFakeEta");
  fDoList.push_back("FakeTisDtAllPt");
  fDoList.push_back("FakeTisDtAllEta");

  fDoList.push_back("FakeTisDtDmFakePt");
  fDoList.push_back("FakeTisDtDmFakeEta");
  fDoList.push_back("FakeTisDtDmAllPt");
  fDoList.push_back("FakeTisDtDmAllEta");

  fDoList.push_back("FakeBdt");
  fDoList.push_back("FakeTip");
  fDoList.push_back("FakeLip");
  fDoList.push_back("FakeInnerChi2");
  fDoList.push_back("FakeOuterChi2");
  fDoList.push_back("FakeChi2LocalPosition");
  fDoList.push_back("FakeChi2LocalMomentum");
  fDoList.push_back("FakeStaTrkMult");
  fDoList.push_back("FakeTmTrkMult");
  fDoList.push_back("FakeDeltaR");
  fDoList.push_back("FakeItrkValidFraction");
  fDoList.push_back("FakeSegmentComp");
  fDoList.push_back("FakeGtrkNormChi2");
  fDoList.push_back("FakeDzRef");
  fDoList.push_back("FakeDxyRef");
  fDoList.push_back("FakeGtrkProb");
  fDoList.push_back("FakeMuonChi2");
  fDoList.push_back("FakeGlbKinkFinder");
  fDoList.push_back("FakeStaRelChi2");
  fDoList.push_back("FakeTrkRelChi2");
  fDoList.push_back("FakeGlbDeltaEtaPhi");
  fDoList.push_back("FakeTimeInOut");
  fDoList.push_back("FakeTimeInOutE");
  fDoList.push_back("FakeNvalidMuonHits");
  fDoList.push_back("FakeNmatchedStations");
  fDoList.push_back("FakeLayersWithHits");
  fDoList.push_back("FakeNumberOfValidTrkHits");
  fDoList.push_back("FakeNumberOfLostTrkHits");
  fDoList.push_back("FakeNumberOfValidPixHits");
  fDoList.push_back("FakeRPChits1");
  fDoList.push_back("FakeRPChits2");
  fDoList.push_back("FakeRPChits3");
  fDoList.push_back("FakeRPChits4");
  fDoList.push_back("FakeCombHits");

  fAnaCuts.clear();
  fAnaCuts.addCut("GoodCand", "good cand", fGoodCand);
  fAnaCuts.addCut("GoodPt", "good pt", fGoodPt);
  fAnaCuts.addCut("GlobalMuon", "global muon ID", fGlobalMuon);

  fAnaCuts.addCut("Good", "all, triggered independently of signal", fGood);
  fAnaCuts.addCut("GoodFake", "fake, triggered independently of signal", fGoodFake);

  fAnaCuts.addCut("TIS", "all, triggered independently of signal", fGoodTIS);
  fAnaCuts.addCut("TISFAKE", "fake, triggered independently of signal", fGoodTISFake);

  fAnaCuts.addCut("TISDT", "all, triggered independently of signal, dist(trigger)", fGoodTISDT);
  fAnaCuts.addCut("TISDTFAKE", "fake, triggered independently of signal, dist(trigger)", fGoodTISDTFake);

  fAnaCuts.addCut("TISDTDM", "all, triggered independently of signal, dist(trigger), dist(muon)", fGoodTISDTDM);
  fAnaCuts.addCut("TISDTDMFAKE", "fake, triggered independently of signal, dist(trigger), dist(muon)", fGoodTISDTDMFake);

  fAnaCuts.dumpAll();
}


// ----------------------------------------------------------------------
plotFake::~plotFake() {

}



// ----------------------------------------------------------------------
void plotFake::init() {
  cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotFake::makeAll(string what) {

  if (what == "dbx") {
    //      fakeRate("fakeData_lambda", "fakeMc_lambda", "FakeTisDtDmFakePt", "FakeTisDtDmAllPt");
    makeSample("fakeData", "ks");
    makeSample("fakeMc", "ks");
    makeOverlay("fakeData_ks", "fakeMc_ks", "Cu");
  }

  if (what == "all" || string::npos != what.find("sample")) {
    init();
    if ((what == "all") || (what == "sample") || ((string::npos != what.find("sample") && string::npos != what.find("ks")))) {
      makeSample("fakeData", "ks");
      makeSample("fakeMc", "ks");
    }
    if ((what == "all") || (what == "sample") || ((string::npos != what.find("sample") && string::npos != what.find("psi")))) {
       makeSample("fakeData", "psi");
       makeSample("fakeMc", "psi");
    }
    if ((what == "all") || (what == "sample") || ((string::npos != what.find("sample") && string::npos != what.find("phi")))) {
      makeSample("fakeData", "phi");
      makeSample("fakeMc", "phi");
    }

    if ((what == "all") || (what == "sample") || ((string::npos != what.find("sample") && string::npos != what.find("lambda")))) {
      makeSample("fakeData", "lambda");
      makeSample("fakeMc", "lambda");
    }

  }

  if (what == "all" || string::npos != what.find("plot")) {
    fTEX.close();
    system(Form("/bin/rm -f %s", fTexFileName.c_str()));
    fTEX.open(fTexFileName.c_str(), ios::app);
    system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
    system(Form("/bin/rm -f %s/sbsctrl_ad*_fake*.pdf", fDirectory.c_str()));
    system(Form("/bin/rm -f %s/ad*_fake*.pdf", fDirectory.c_str()));

    if ((what == "all") || (what == "plot") || (string::npos != what.find("plot") && string::npos != what.find("ks"))) {
      makeOverlay("fakeData_ks", "fakeMc_ks", "Cu");
    }
    if ((what == "all") || (what == "plot") || (string::npos != what.find("plot") && string::npos != what.find("psi"))) {
      makeOverlay("fakeData_psi", "fakeMc_psi", "Cu");
    }
    if ((what == "all") || (what == "plot") || (string::npos != what.find("plot") && string::npos != what.find("phi"))) {
      makeOverlay("fakeData_phi", "fakeMc_phi", "Cu");
    }

    if ((what == "all") || (what == "plot") || (string::npos != what.find("plot") && string::npos != what.find("lambda"))) {
      makeOverlay("fakeData_lambda", "fakeMc_lambda", "Cu");
    }
  }

  if (what == "all" || string::npos != what.find("fakerate")) {
    if ((what == "all") || (what == "fakerate") || (string::npos != what.find("ks"))) {
      fakeRate("fakeData_ks", "fakeMc_ks", "FakePt", "AllPt", 0.1);
      fakeRate("fakeData_ks", "fakeMc_ks", "FakeEta", "AllEta", 0.1);

      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisFakePt", "FakeTisAllPt");
      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisFakeEta", "FakeTisAllEta");

      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisDtFakePt", "FakeTisDtAllPt");
      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisDtFakeEta", "FakeTisDtAllEta");

      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisDtDmFakePt", "FakeTisDtDmAllPt");
      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisDtDmFakeEta", "FakeTisDtDmAllEta");
    }

    if ((what == "all") || (what == "fakerate") || (string::npos != what.find("phi"))) {
      fakeRate("fakeData_phi", "fakeMc_phi", "FakeTisDtDmFakePt", "FakeTisDtDmAllPt");
      fakeRate("fakeData_phi", "fakeMc_phi", "FakeTisDtDmFakeEta", "FakeTisDtDmAllEta");
    }

    if ((what == "all") || (what == "fakerate") || (string::npos != what.find("lambda"))) {
      fakeRate("fakeData_lambda", "fakeMc_lambda", "FakeTisDtDmFakePt", "FakeTisDtDmAllPt");
      fakeRate("fakeData_lambda", "fakeMc_lambda", "FakeTisDtDmFakeEta", "FakeTisDtDmAllEta");
    }

  }

  // plotMass("fakeData_ks", "Cu");

}


// ----------------------------------------------------------------------
void  plotFake::massPlots(std::string varname) {

  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  vector<string> modes;
  modes.push_back("ks");
  modes.push_back("phi");
  modes.push_back("lambda");
  modes.push_back("psi");
  string hname(""), fname("");
  TH1D *h1(0);
  gStyle->SetOptStat(0);
  for (unsigned int im = 0; im < modes.size(); ++im) {
    for (unsigned int ic = 0; ic < 3; ++ic) {
      hname = Form("ad%d_fakeData_%s_%s", ic, modes[im].c_str(), varname.c_str());
      h1 = (TH1D*)fHistFile->Get(hname.c_str());
      cout << "hname = " << hname << endl;
      if (string::npos != modes[im].find("ks")) {
	setTitles(h1, "m_{#pi#pi} [GeV]", "Candidates", 0.05, 1.1, 1.8);
      } else if (string::npos != modes[im].find("phi")) {
	setTitles(h1, "m_{KK} [GeV]", "Candidates", 0.05, 1.1, 1.8);
      } else if (string::npos != modes[im].find("lambda")) {
	setTitles(h1, "m_{p#pi} [GeV]", "Candidates", 0.05, 1.1, 1.8);
      } else if (string::npos != modes[im].find("psi")) {
	setTitles(h1, "m_{#mu#mu} [GeV]", "Candidates", 0.05, 1.1, 1.8);
      }
      shrinkPad(0.15, 0.2);
      if (h1) h1->Draw();
      savePad(Form("fakemass_ad%d_%s_%s.pdf", ic, modes[im].c_str(), varname.c_str()));
    }
  }
}

// // ----------------------------------------------------------------------
// void plotFake::plotMass(string sample, string selection) {

//   cout << "fHistFileName: " << fHistFileName;
//   fHistFile = TFile::Open(fHistFileName.c_str());
//   cout << " opened " << endl;

//   gStyle->SetOptFit(0);
//   gStyle->SetOptStat(0);
//   gStyle->SetOptTitle(0);

//   TH1D *h(0);
//   fIF->fLo = 5.0;
//   fIF->fHi = 5.5;
//   TF1 *lBg = fIF->expoErr(fIF->fLo, fIF->fHi);
//   string header;
//   if (string::npos != sample.find("psi")) header = "J/#psi #rightarrow #mu #mu";
//   else if (string::npos != sample.find("ks")) header = "K_{S} #rightarrow #pi^{+} #pi^{-}";
//   else if (string::npos != sample.find("phi")) header = "#phi #rightarrow K{+} K^{-}";
//   else if (string::npos != sample.find("lambda")) header = "#Lambda #rightarrow p #pi^{+}";
//   tl->SetTextFont(42);

//   double xPos(0.58);

//   for (unsigned int i = 0; i < fNchan; ++i) {
//     h = (TH1D*)gDirectory->Get(Form("ad%d_%s_FakePtMass%s", i, sample.c_str(), selection.c_str()));
//     if (!h) break;
//     TF1 *f1 = fIF->expoErrGauss(h, 5.28, 0.04);
//     setTitles(h, "m [GeV]", Form("Entries / %4.3f GeV", h->GetBinWidth(1)), 0.05, 1.1, 1.3);
//     if (header != "Dimuon") {
//       h->Fit(f1, "", "e");

//       for (int i = 0; i < lBg->GetNpar(); ++i) {
// 	lBg->SetParameter(i, f1->GetParameter(3+i));
//       }

//       double c  = f1->GetParameter(0);
//       cout << "OVERALL INTEGRAL: " << f1->Integral(5.15, 5.45) << " BACKGROUND INTEGRAL: " << lBg->Integral(5.15, 5.45) << endl;
//       c = f1->Integral(5.15, 5.45) - lBg->Integral(5.15, 5.45);
//       double cE = f1->GetParError(0);
//       double ierr = f1->IntegralError(5.15, 5.45)/h->GetBinWidth(1);

//       double signal = c/h->GetBinWidth(1);
//       double signalE(0.);
//       if (ierr > TMath::Sqrt(signal)) {
// 	signalE = ierr;
//       } else {
// 	signalE = cE/c*signal;
//       }
//       tl->SetTextSize(0.025);
//       tl->DrawLatexNDC(xPos, 0.70, Form("Signal: %5.1f  #pm %5.1f", signal, signalE));
//       tl->DrawLatexNDC(xPos, 0.66, Form("Mass:   %5.4f  #pm %5.4f GeV", f1->GetParameter(1), f1->GetParError(1)));
//       tl->DrawLatexNDC(xPos, 0.62, Form("Width:  %5.4f  #pm %5.4f GeV", f1->GetParameter(2), f1->GetParError(2)));
//     } else {
//       h->Draw("hist");
//     }

//     tl->SetTextSize(0.05);
//     tl->DrawLatexNDC(xPos, 0.80, header.c_str());
//     tl->SetTextSize(0.025);
//     tl->DrawLatexNDC(xPos, 0.75, Form("%2.1f < |#eta(#mu_{f})| < %2.1f", fCuts[i]->metaMin, fCuts[i]->metaMax));

//     stamp(0., fStampCms, fStampString, 0., fStampLumi);

//     c0->SaveAs(Form("%s/mass_ad%d_%d_%s_%s.pdf", fDirectory.c_str(), i, fYear, sample.c_str(), selection.c_str()));
//   }
//   fHistFile->Close();
// }


// ----------------------------------------------------------------------
void plotFake::makeSample(std::string dataset, std::string sample, int nevents, int nstart) {

  string dir("");

  // -- dump histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;

  fSample = dataset + "_" + sample;
  cout << "fSample  = " << fSample << endl;

  fChannelSample.clear();
  for (int ic = 0; ic < fNchan; ++ic) {
    fChannelSample.push_back(Form("ad%s_%s", fChannelList[ic].c_str(), fSample.c_str()));
  }

  if (string::npos != fSample.find("Mc")) {
    fIsMC = true;
  } else {
    fIsMC = false;
  }

  MASSMIN = 0.45;
  MASSMAX = 0.55;
  if (string::npos != sample.find("ks")) {
    fMode = FAKEKS;
    dir = "candAnaFake310";
    // BGLBOXMIN = 0.450;
    // BGLBOXMAX = 0.465;
    // SIGBOXMIN = 0.480;
    // SIGBOXMAX = 0.515;
    // BGHBOXMIN = 0.530;
    // BGHBOXMAX = 0.550;
  }

  if (string::npos != sample.find("psi")) {
    fMode = FAKEPSI;
    dir = "candAnaFake443";
    // BGLBOXMIN = 2.800;
    // BGLBOXMAX = 2.900;
    // SIGBOXMIN = 3.050;
    // SIGBOXMAX = 3.150;
    // BGHBOXMIN = 3.200;
    // BGHBOXMAX = 3.300;
  }

  if (string::npos != sample.find("phi")) {
    fMode = FAKEPHI;
    dir = "candAnaFake333";
    // BGLBOXMIN = 0.990;
    // BGLBOXMAX = 1.005;
    // SIGBOXMIN = 1.010;
    // SIGBOXMAX = 1.030;
    // BGHBOXMIN = 1.035;
    // BGHBOXMAX = 1.045;
  }

  if (string::npos != sample.find("lambda")) {
    fMode = FAKELAMBDA;
    dir = "candAnaFake3122";
    // BGLBOXMIN = 1.095;
    // BGLBOXMAX = 1.105;
    // SIGBOXMIN = 1.110;
    // SIGBOXMAX = 1.122;
    // BGHBOXMIN = 1.130;
    // BGHBOXMAX = 1.140;
  }


  if (fIsMC) dir = "candAnaFakeMC";


  // -- mass box definitions are in there!
  bookDistributions();

  TTree *t = getTree(dataset, dir, "fakeTree");
  if (0 == t) {
    cout << "fakeTree for sample = " << sample << " not found in dataset = " << dataset << endl;
    return;
  }
  setupTree(t);
  fCds = fSample;
  loopOverTree(t, 1, nevents, nstart);

  fHistFile->Write();
  fHistFile->Close();

}


// ----------------------------------------------------------------------
void plotFake::makeOverlay(string what1, string what2, string selection) {

  cout << "fHistFileName: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  for (unsigned int i = 0; i < fChannelList.size(); ++i) {
    cout << "===> sbsDistributions " << Form("ad%s_%s", fChannelList[i].c_str(), what1.c_str()) << endl;
    sbsDistributions(Form("ad%s_%s", fChannelList[i].c_str(), what1.c_str()), selection);
    sbsDistributions(Form("ad%s_%s", fChannelList[i].c_str(), what2.c_str()), selection);

    for (unsigned int id = 0; id < fDoList.size(); ++id) {
      c0->cd();
      overlay(Form("sbs_ad%s_%s_%s%s", fChannelList[i].c_str(), what1.c_str(), fDoList[id].c_str(), selection.c_str()),
	      Form("sbs_ad%s_%s_%s%s", fChannelList[i].c_str(), what2.c_str(), fDoList[id].c_str(), selection.c_str())
	      );
      savePad(Form("fakeoverlay_%s_ad%s_%s_ad%s_%s.pdf", fDoList[id].c_str(), fChannelList[i].c_str(), what1.c_str(), fChannelList[i].c_str(), what2.c_str()));
    }
  }
  fHistFile->Close();

}

// ----------------------------------------------------------------------
void plotFake::fillDistributions() {



}


// ----------------------------------------------------------------------
void plotFake::bookDistributions() {

  adsetFake *a(0);
  for (unsigned int i = 0; i < fChannelList.size(); ++i) {
    string mapname = Form("ad%s_%s", fChannelList[i].c_str(), fSample.c_str());
    string name = Form("%s_", mapname.c_str());

    if (string::npos != fSample.find("ks")) {
      if (0 == i) {
	SIGBOXMIN = 0.482;
	SIGBOXMAX = 0.510;
	BGLBOXMIN = 0.450;
	BGLBOXMAX = 0.465;
	BGHBOXMIN = 0.530;
	BGHBOXMAX = 0.550;
      } else if (1 == i) {
	SIGBOXMIN = 0.480;
	SIGBOXMAX = 0.515;
	BGLBOXMIN = 0.450;
	BGLBOXMAX = 0.465;
	BGHBOXMIN = 0.530;
	BGHBOXMAX = 0.550;
      } else {
	SIGBOXMIN = 0.475;
	SIGBOXMAX = 0.520;
	BGLBOXMIN = 0.450;
	BGLBOXMAX = 0.465;
	BGHBOXMIN = 0.530;
	BGHBOXMAX = 0.550;
      }
    }

    cout << "fSample: " << fSample << endl;
    cout << "SIG: " << SIGBOXMIN << " .. " << SIGBOXMAX << endl;
    cout << "BGL: " << BGLBOXMIN << " .. " << BGLBOXMAX << endl;
    cout << "BGH: " << BGHBOXMIN << " .. " << BGHBOXMAX << endl;

    a = new adsetFake();
    a->fpFakeEta  = bookDistribution(Form("%sFakeEta", name.c_str()), "#eta", "GoodFake", 40, -2.4, 2.4);
    a->fpFakePt   = bookDistribution(Form("%sFakePt", name.c_str()), "p_{T} [GeV]", "GoodFake", 40, 0., 20.);
    a->fpAllEta  = bookDistribution(Form("%sAllEta", name.c_str()), "#eta", "Good", 40, -2.4, 2.4);
    a->fpAllPt   = bookDistribution(Form("%sAllPt", name.c_str()), "p_{T} [GeV]", "Good", 40, 0., 20.);

    a->fpFakeBdt       = bookDistribution(Form("%sFakeBdt", name.c_str()), "BDT", "GlobalMuon", 50, -1., 1.);
    a->fpFakeTip       = bookDistribution(Form("%sFakeTip", name.c_str()), "TIP [cm]", "GlobalMuon", 50, 0., 2.);
    a->fpFakeLip       = bookDistribution(Form("%sFakeLip", name.c_str()), "LIP [cm]", "GlobalMuon", 50, 0., 2.);
    a->fpFakeInnerChi2 = bookDistribution(Form("%sFakeInnerChi2", name.c_str()), "inner track #chi^{2}", "GlobalMuon", 51, 0., 20.);
    a->fpFakeOuterChi2 = bookDistribution(Form("%sFakeOuterChi2", name.c_str()), "outer track #chi^{2}", "GlobalMuon", 51, 0., 02.);

    a->fpFakeChi2LocalPosition = bookDistribution(Form("%sFakeChi2LocalPosition", name.c_str()), "local position #chi^{2}", "GlobalMuon", 51, 0., 102.);
    a->fpFakeChi2LocalMomentum = bookDistribution(Form("%sFakeChi2LocalMomentum", name.c_str()), "local momentum #chi^{2}", "GlobalMuon", 51, 0., 102.);
    a->fpFakeStaTrkMult = bookDistribution(Form("%sFakeStaTrkMult", name.c_str()), "STA trk multipicity", "GlobalMuon", 12, -2., 10.);
    a->fpFakeTmTrkMult = bookDistribution(Form("%sFakeTmTrkMult", name.c_str()), "TM trk multiplicity", "GlobalMuon", 20, 0., 20.);

    a->fpFakeDeltaR = bookDistribution(Form("%sFakeDeltaR", name.c_str()), "deltaR", "GlobalMuon", 40, 0., 1.0);
    a->fpFakeItrkValidFraction = bookDistribution(Form("%sFakeItrkValidFraction", name.c_str()), "inner track valid fraction", "GlobalMuon", 50, 0., 1.02);
    a->fpFakeSegmentComp = bookDistribution(Form("%sFakeSegmentComp", name.c_str()), "segment compatibility", "GlobalMuon", 50, 0., 1.02);
    a->fpFakeGtrkNormChi2 = bookDistribution(Form("%sFakeGtrkNormChi2", name.c_str()), "global track norm. #chi^{2}", "GlobalMuon", 40, 0., 12.);
    a->fpFakeDzRef = bookDistribution(Form("%sFakeDzRef", name.c_str()), "dzrf", "GlobalMuon", 40, -20., 20.);
    a->fpFakeDxyRef = bookDistribution(Form("%sFakeDxyRef", name.c_str()), "dxyref", "GlobalMuon", 100, -1., 1.);
    a->fpFakeGtrkProb = bookDistribution(Form("%sFakeGtrkProb", name.c_str()), "global track prob", "GlobalMuon", 51, 0., 102.);
    a->fpFakeMuonChi2 = bookDistribution(Form("%sFakeMuonChi2", name.c_str()), "muon #chi^{2}", "GlobalMuon", 51, 0., 1020.);
    a->fpFakeGlbKinkFinder = bookDistribution(Form("%sFakeGlbKinkFinder", name.c_str()), "log10(GlbKinkFinder)", "GlobalMuon", 50, -5., 15.);
    a->fpFakeStaRelChi2 = bookDistribution(Form("%sFakeStaRelChi2", name.c_str()), "StaRelChi2", "GlobalMuon", 51, 0., 102.);
    a->fpFakeTrkRelChi2 = bookDistribution(Form("%sFakeTrkRelChi2", name.c_str()), "TrkRelChi2", "GlobalMuon", 50, 0., 20.);
    a->fpFakeGlbDeltaEtaPhi = bookDistribution(Form("%sFakeGlbDeltaEtaPhi", name.c_str()), "GlbDeltaEtaPhi", "GlobalMuon", 50, -2., 4.);
    a->fpFakeTimeInOut = bookDistribution(Form("%sFakeTimeInOut", name.c_str()), "TimeInOut", "GlobalMuon", 50, -200., 200.);
    a->fpFakeTimeInOutE = bookDistribution(Form("%sFakeTimeInOutE", name.c_str()), "TimeInOutE", "GlobalMuon", 50, 0., 10.);
    a->fpFakeTimeInOutS = bookDistribution(Form("%sFakeTimeInOutS", name.c_str()), "TimeInOut/TimeInOutE", "GlobalMuon", 50, -100., 100.);
    a->fpFakeNvalidMuonHits = bookDistribution(Form("%sFakeNvalidMuonHits", name.c_str()), "NvalidMuonHits", "GlobalMuon", 51, 0., 51.);
    a->fpFakeNmatchedStations = bookDistribution(Form("%sFakeNmatchedStations", name.c_str()), "NmatchedStations", "GlobalMuon", 10, 0., 10.);
    a->fpFakeLayersWithHits = bookDistribution(Form("%sFakeLayersWithHits", name.c_str()), "LayersWithHits", "GlobalMuon", 20, 0., 20.);
    a->fpFakeNumberOfValidTrkHits = bookDistribution(Form("%sFakeNumberOfValidTrkHits", name.c_str()), "NumberOfValidTrkHits", "GlobalMuon", 30, 0., 30.);
    a->fpFakeNumberOfLostTrkHits = bookDistribution(Form("%sFakeNumberOfLostTrkHits", name.c_str()), "NumberOfLostTrkHits", "GlobalMuon", 30, 0., 30.);
    a->fpFakeNumberOfValidPixHits = bookDistribution(Form("%sFakeNumberOfValidPixHits", name.c_str()), "NumberOfValidPixHits", "GlobalMuon", 12, 0., 12.);
    a->fpFakeRPChits1 = bookDistribution(Form("%sFakeRPChits1", name.c_str()), "RPChits1", "GlobalMuon", 10, 0., 10.);
    a->fpFakeRPChits2 = bookDistribution(Form("%sFakeRPChits2", name.c_str()), "RPChits2", "GlobalMuon", 10, 0., 10.);
    a->fpFakeRPChits3 = bookDistribution(Form("%sFakeRPChits3", name.c_str()), "RPChits3", "GlobalMuon", 10, 0., 10.);
    a->fpFakeRPChits4 = bookDistribution(Form("%sFakeRPChits4", name.c_str()), "RPChits4", "GlobalMuon", 10, 0., 10.);
    a->fpFakeCombHits = bookDistribution(Form("%sFakeCombHits", name.c_str()), "CombHits", "GlobalMuon", 35, 0., 35.);

    a->fpFakeTisAllEta  = bookDistribution(Form("%sFakeTisAllEta", name.c_str()), "#eta", "TIS", 40, -2.4, 2.4);
    a->fpFakeTisAllPt   = bookDistribution(Form("%sFakeTisAllPt", name.c_str()), "p_{T} [GeV]", "TIS", 40, 0., 20.);
    a->fpFakeTisFakeEta = bookDistribution(Form("%sFakeTisFakeEta", name.c_str()), "#eta", "TISFAKE", 40, -2.4, 2.4);
    a->fpFakeTisFakePt  = bookDistribution(Form("%sFakeTisFakePt", name.c_str()), "p_{T} [GeV]", "TISFAKE", 40, 0., 20.);

    a->fpFakeTisDtAllEta  = bookDistribution(Form("%sFakeTisDtAllEta", name.c_str()), "#eta", "TISDT", 40, -2.4, 2.4);
    a->fpFakeTisDtAllPt   = bookDistribution(Form("%sFakeTisDtAllPt", name.c_str()), "p_{T} [GeV]", "TISDT", 40, 0., 20.);
    a->fpFakeTisDtFakeEta = bookDistribution(Form("%sFakeTisDtFakeEta", name.c_str()), "#eta", "TISDTFAKE", 40, -2.4, 2.4);
    a->fpFakeTisDtFakePt  = bookDistribution(Form("%sFakeTisDtFakePt", name.c_str()), "p_{T} [GeV]", "TISDTFAKE", 40, 0., 20.);

    a->fpFakeTisDtDmAllEta  = bookDistribution(Form("%sFakeTisDtDmAllEta", name.c_str()), "#eta", "TISDTDM", 40, -2.4, 2.4);
    a->fpFakeTisDtDmAllPt   = bookDistribution(Form("%sFakeTisDtDmAllPt", name.c_str()), "p_{T} [GeV]", "TISDTDM", 40, 0., 20.);
    a->fpFakeTisDtDmFakeEta = bookDistribution(Form("%sFakeTisDtDmFakeEta", name.c_str()), "#eta", "TISDTDMFAKE", 40, -2.4, 2.4);
    a->fpFakeTisDtDmFakePt  = bookDistribution(Form("%sFakeTisDtDmFakePt", name.c_str()), "p_{T} [GeV]", "TISDTDMFAKE", 40, 0., 20.);


    fAdMap.insert(make_pair(mapname, a));
    cout << "bookDistributions: mapname = " << mapname << endl;
  }

}

// ----------------------------------------------------------------------
AnalysisDistribution* plotFake::bookDistribution(string hn, string ht, std::string hc, int nbins, double lo, double hi) {

  double masslo(0.0), masshi(0.0);
  if (string::npos != fSample.find("ks")) {
    masslo = 0.45;
    masshi = 0.55;
  }
  if (string::npos != fSample.find("phi")) {
    masslo = 0.99;
    masshi = 1.05;
  }
  if (string::npos != fSample.find("lambda")) {
    masslo = 1.08;
    masshi = 1.15;
  }
  if (string::npos != fSample.find("psi")) {
    masslo = 2.75;
    masshi = 3.35;
  }

  AnalysisDistribution *p = new AnalysisDistribution(hn.c_str(), ht.c_str(), nbins, lo, hi, masslo, masshi);
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX);
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX);
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX);
  p->setAnalysisCuts(&fAnaCuts, hc.c_str());
  p->setPreselCut(&fGoodCand);
  return p;
}


// ----------------------------------------------------------------------
void plotFake::sbsDistributions(string sample, string selection, std::string what) {
  cout << "plotFake::sbsDistributions(" << sample << ", " << selection << ", " << what << ")" << endl;
  string sbsControlPlotsFileName = Form("sbsctrl");
  AnalysisDistribution a(Form("%s_FakePt", sample.c_str()));
  a.fVerbose = 1;
  a.fControlPlotsFileName = sbsControlPlotsFileName;
  a.fDirectory = fDirectory;

  a.fpIF->resetLimits();

  int type(0);
  if (string::npos != sample.find("ks")) {
    type = 1;
    a.fMassPeak  = 0.498;
    a.fMassSigma = 0.005;
    a.fMassLo    = 0.450;
    a.fMassHi    = 0.550;
  } else if (string::npos != sample.find("psi")) {
    type = 1;
    a.fMassPeak  = 3.097;
    a.fMassSigma = 0.030;
    a.fMassLo    = 2.930;
    a.fMassHi    = 3.280;
  } else if (string::npos != sample.find("lambda")) {
    type = 1;
    a.fMassPeak  = 1.116;
    a.fMassSigma = 0.002;
    a.fMassLo    = 1.095;
    a.fMassHi    = 1.140;
    if (0 && string::npos == sample.find("Mc")) {
      string bla =  Form("%s_FakePtMassNm", sample.c_str());
      TH1D *h = (TH1D*)gDirectory->Get(Form("%s", bla.c_str()));
      cout << "=> Looking for prefit histogram " << bla.c_str() << ", at h = " << h << " with nentries = " << h->GetSumOfWeights() << endl;
      fIF->fLo = a.fMassLo;
      fIF->fHi = a.fMassHi;
      fIF->limitPar(0, 0., 1.e7);
      TF1 *f1 = fIF->pol1gauss(h, a.fMassPeak, a.fMassSigma);
      cout << "now prefitting: h = " << h << " f1 = " << f1 << endl;
      h->Fit(f1, "lsr", "", a.fMassLo, a.fMassHi);
      savePad(Form("prefit-%s-%s-%s.pdf", sample.c_str(), selection.c_str(), what.c_str()));
      a.fpIF->limitPar(1, f1->GetParameter(1) - 3.*f1->GetParError(1), f1->GetParameter(1) + 3.*f1->GetParError(1));
      a.fpIF->limitPar(2, f1->GetParameter(2) - 3.*f1->GetParError(2), f1->GetParameter(2) + 3.*f1->GetParError(2));
      cout << "done prefitting: h = " << h << " f1 = " << f1 << endl;
      cout << "done prefitting: h = " << h << " f1 = " << f1 << endl;
    }
  } else if (string::npos != sample.find("phi")) {
    type = 2;
    a.fMassPeak  = 1.019;
    a.fMassSigma = 0.008;
    a.fMassLo    = 0.990;
    a.fMassHi    = 1.045;
  } else {
    a.fMassPeak = 5.27;
  }

  // -- override the above choice in case of MC
  if (string::npos != sample.find("Mc")) {
    type = 0; // signal window
  }

  cout << "----------------------------------------------------------------------" << endl;
  cout << "type = " << type << endl;
  cout << "----------------------------------------------------------------------" << endl;
  TH1D *h(0);
  bool restricted = (what != "");
  string bla;
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    bla =  Form("%s_%s", sample.c_str(), fDoList[i].c_str());
    if (0 == type) {
      cout << "=> Looking for signal histogram " << Form("%s%s0", bla.c_str(), selection.c_str()) << endl;
      h = (TH1D*)gDirectory->Get(Form("%s%s0", bla.c_str(), selection.c_str()));
      cout << "=> cloning into  "
	   << Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), selection.c_str()) << endl;
      h = (TH1D*)h->Clone(Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), selection.c_str()));
    } else if (1 == type) {
      cout << "=> sbsDistribution histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      h = a.sbsDistribution(bla.c_str(), selection.c_str());
    } else if (2 == type) {
      cout << "=> sbsDistributionPhiKK histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      h = a.sbsDistributionPhiKK(bla.c_str(), selection.c_str());
    }

    cout << "  Title: " << h->GetTitle() << " with integral: " << h->GetSumOfWeights() << endl;

  }

  cout << "sbsDistributions: created histogram with name: " << h->GetName() << " in directory " << h->GetDirectory()->GetName() << endl;

}


// ----------------------------------------------------------------------
void plotFake::overlay(string sample1, string sample2, string what) {

  cout << "overlay gDirectory = " << gDirectory->GetName() << endl;

  TH1D *h1 = (TH1D*)gDirectory->Get(sample1.c_str());
  TH1D *h2 = (TH1D*)gDirectory->Get(sample2.c_str());

  bool doLegend(true);
  bool leftLegend(false);

  cout << "plotFake::overlay samples " << sample1 << " and " << sample2 << endl;
  cout << "plotFake::overlay h1 " << h1 << " and " << h2 << endl;

  bool drawGrid(true);

  c0->SetTopMargin(0.1);
  c0->SetBottomMargin(0.15);
  c0->SetRightMargin(0.05);
  if (drawGrid) c0->SetGridx();

  if (h2->GetSumOfWeights() > 0.) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

  h1->SetTitle("");
  h1->SetMinimum(0.01);
  double ymax = h1->GetMaximum();
  if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();

  c0->SetLogy(0);
  if ((string::npos != sample1.find("Chi2"))
      || (string::npos != sample1.find("ItrkValidFraction"))
      || (string::npos != sample1.find("GtrkProb"))
      || (string::npos != sample1.find("GlbDeltaEtaPhi"))
      || (string::npos != sample1.find("TimeInOut"))
      || (string::npos != sample1.find("Tip"))
      || (string::npos != sample1.find("Lip"))
      ) {
    h1->SetMinimum(0.5);
    ymax *= 5.;
    c0->SetLogy(1);
  }

  if (string::npos != sample1.find("NvalidMuonHits")) leftLegend = true;
  if (string::npos != sample1.find("LayersWithHits")) leftLegend = true;

  h1->SetStats(0);
  string hname = h1->GetName();
  if (string::npos != hname.find("Data")) {
    h1->Draw("e");
  } else {
    h1->Draw();
  }

  string label1("bla"), label2("bla");
  if (string::npos != sample1.find("psi")) {
    label1 = "muons";
  } else if (string::npos != sample2.find("ks")) {
    label1 = "pions";
  } else if (string::npos != sample2.find("phi")) {
    label1 = "kaons";
  } else if (string::npos != sample2.find("lambda")) {
    label1 = "protons";
  }

  if (string::npos != sample2.find("psi")) {
    label2 = "muons";
    setFilledHist(h2, kBlue, kBlue, 3365);
  } else if (string::npos != sample2.find("ks")) {
    label2 = "pions";
    setFilledHist(h2, kRed, kRed, 3365);
  } else if (string::npos != sample2.find("phi")) {
    label2 = "kaons";
    setFilledHist(h2, kRed, kRed, 3365);
  } else if (string::npos != sample2.find("lambda")) {
    label2 = "protons";
    setFilledHist(h2, kRed, kRed, 3365);
  }

  string header("");
  if (label1 == label2) {
    header = label1;
    label1 = "";
    label2 = "";
  }

  char loption1[100], loption2[100];
  if (string::npos != sample1.find("Data")) {
    label1 += " data";
    sprintf(loption1, "ep");
  }

  if (string::npos != sample1.find("Mc")) {
    label1 += " MC";
    sprintf(loption1, "f");
  }

  if (string::npos != sample2.find("Data")) {
    label2 += " data";
    sprintf(loption2, "ep");
  }

  if (string::npos != sample2.find("Mc")) {
    label2 += " MC";
    sprintf(loption2, "f");
  }


  h2->Draw("samehist");
  h1->SetMaximum(1.2*ymax);

  if (doLegend) {
    if (leftLegend) {
      newLegend(0.26, 0.7, 0.46, 0.87);
    } else {
      newLegend(0.70, 0.7, 0.90, 0.87);
    }

    legg->SetHeader(header.c_str());
    legg->SetTextSize(0.05);
    legg->AddEntry(h1, label1.c_str(), loption1);
    legg->AddEntry(h2, label2.c_str(), loption2);
    legg->Draw();
  }

  //  stamp(0.18, fStampCms, fStampString, 0.4, fStampLumi);

}



// ----------------------------------------------------------------------
void plotFake::fakeRate(string dataset1, string dataset2, string varF, string varA, double ymax) {

  cout << "fHistFileName: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  string label1("bla"), label2("bla");
  if (string::npos != dataset1.find("psi")) {
    label1 = "muons";
  } else if (string::npos != dataset1.find("ks")) {
    label1 = "pions";
  } else if (string::npos != dataset1.find("phi")) {
    label1 = "kaons";
  } else if (string::npos != dataset1.find("lambda")) {
    label1 = "protons";
  }

  if (string::npos != dataset2.find("psi")) {
    label2 = "muons";
  } else if (string::npos != dataset2.find("ks")) {
    label2 = "pions";
  } else if (string::npos != dataset2.find("phi")) {
    label2 = "kaons";
  } else if (string::npos != dataset2.find("lambda")) {
    label2 = "protons";
  }

  string header("");
  if (label1 == label2) {
    header = label1;
    label1 = "";
    label2 = "";
  }

  char loption1[100], loption2[100];
  if (string::npos != dataset1.find("Data")) {
    label1 += " data";
    sprintf(loption1, "ep");
  }

  if (string::npos != dataset1.find("Mc")) {
    label1 += " MC";
    sprintf(loption1, "l");
  }

  if (string::npos != dataset2.find("Data")) {
    label2 += " data";
    sprintf(loption2, "ep");
  }

  if (string::npos != dataset2.find("Mc")) {
    label2 += " MC";
    sprintf(loption2, "l");
  }


  for (unsigned int i = 0; i < fChannelList.size(); ++i) {
    cout << "===> sbsDistributions " << Form("ad%s_%s", fChannelList[i].c_str(), dataset1.c_str()) << "Si" << varF << endl;
    sbsDistributions(Form("ad%s_%s", fChannelList[i].c_str(), dataset1.c_str()), "Si", varF);
    cout << "===> sbsDistributions " << Form("ad%s_%s", fChannelList[i].c_str(), dataset2.c_str()) << "Si" << varF << endl;
    sbsDistributions(Form("ad%s_%s", fChannelList[i].c_str(), dataset2.c_str()), "Si", varF);

    cout << "===> sbsDistributions " << Form("ad%s_%s", fChannelList[i].c_str(), dataset1.c_str()) << "Si" << varA << endl;
    sbsDistributions(Form("ad%s_%s", fChannelList[i].c_str(), dataset1.c_str()), "Si", varA);
    cout << "===> sbsDistributions " << Form("ad%s_%s", fChannelList[i].c_str(), dataset2.c_str()) << "Si" << varA << endl;
    sbsDistributions(Form("ad%s_%s", fChannelList[i].c_str(), dataset2.c_str()), "Si", varA);

    c0->Clear();

    TH1D *h1p = (TH1D*)gDirectory->Get(Form("sbs_ad%s_%s_%sSi", fChannelList[i].c_str(), dataset1.c_str(), varF.c_str()));
    TH1D *h1a = (TH1D*)gDirectory->Get(Form("sbs_ad%s_%s_%sSi", fChannelList[i].c_str(), dataset1.c_str(), varA.c_str()));

    TH1D *h2p = (TH1D*)gDirectory->Get(Form("sbs_ad%s_%s_%sSi", fChannelList[i].c_str(), dataset2.c_str(), varF.c_str()));
    TH1D *h2a = (TH1D*)gDirectory->Get(Form("sbs_ad%s_%s_%sSi", fChannelList[i].c_str(), dataset2.c_str(), varA.c_str()));


    h1p->Divide(h1a);
    h2p->Divide(h2a);

    h1p->SetMinimum(0.);
    h1p->SetMaximum(ymax);
    h1p->SetTitle("");
    h1p->Draw();
    setHist(h2p, kBlue);
    h2p->Draw("histsame");


    newLegend(0.21, 0.7, 0.41, 0.87);
    legg->SetHeader(header.c_str());
    legg->SetTextSize(0.05);
    legg->AddEntry(h1p, label1.c_str(), loption1);
    legg->AddEntry(h2p, label2.c_str(), loption2);
    legg->Draw();

    savePad(Form("fakerate_%s_ad%s_%s_ad%s_%s.pdf", varA.c_str(), fChannelList[i].c_str(), dataset1.c_str(), fChannelList[i].c_str(), dataset2.c_str()));


  }





}

// ----------------------------------------------------------------------
void plotFake::playKs(string cuts, string name) {
  string dsname = "fakeData";
  string dir = "candAnaFake310";
  TTree *T = getTree(dsname, dir, "events");

  TH1D *h = (TH1D*)gDirectory->Get("h1");
  if (h) delete h;
  TH1D *h1 = new TH1D("h1", "", 100, 0.45, 0.55); h1->Sumw2();
  h1->SetMinimum(0.);
  h1->SetNdivisions(505, "X");

  T->Draw("m>>h1", cuts.c_str());
  cout << "==> h1 entries: " << h1->GetSumOfWeights() << endl;
  fitKs(h1);
  cout << cuts << endl;
  savePad(Form("playKs-mass-%s.pdf", name.c_str()));
}

// ----------------------------------------------------------------------
void plotFake::fitKs(TH1D *h) {
  fIF->fVerbose = true;
  fIF->fVerbose = false;
  fIF->resetLimits();
  fIF->limitPar(1, 0.495, 0.503);
  fIF->limitPar(2, 0.004, 0.008);
  TF1* f1 = fIF->pol1gauss2c(h, 0.498, 0.006);
  //  f1->FixParameter(4, 0.);
  TFitResultPtr r = h->Fit(f1, "lsq", "e");
  double bwidth = h->GetBinWidth(h->FindBin(1));
  double ypeak = f1->GetParameter(0);
  double xpeak = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);
  double nSig(5.0);
  double xmin = xpeak - nSig*sigma;
  double xmax = xpeak + nSig*sigma;
  xmin = 0.485;
  xmax = 0.515;
  double aintegral  = f1->Integral(xmin, xmax)/bwidth;
  // -- (slight?) overestimate of signal integral error by including the background
  double fintegralE = f1->IntegralError(xmin, xmax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/bwidth;
  // -- now set constant of pol0 to zero and integrate over signal only
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  double fintegral  = f1->Integral(xmin, xmax)/bwidth;
  if (fintegral > 0) {
    double err1 = f1->GetParError(0)/f1->GetParameter(0);
    double err2 = f1->GetParError(2)/f1->GetParameter(2);
    double errT = TMath::Sqrt(err1*err1 + err2*err2);
    fYieldE = errT*fintegral;
    fYieldE = fintegralE;
    fYield  = fintegral;
  } else {
    double fallbackXmin(0.485);
    double fallbackXmax(0.511);
    int nsgbins = h->FindBin(fallbackXmax) - h->FindBin(fallbackXmin) + 1;
    double sg = h->Integral(h->FindBin(fallbackXmin), h->FindBin(fallbackXmax));
    int nbgbins = h->FindBin(fallbackXmax) - 1 - 1 + 1;
    nbgbins    += h->GetNbinsX() - (h->FindBin(fallbackXmax) + 1) + 1;
    double bg = h->Integral(1, h->FindBin(fallbackXmin)-1) +  h->Integral(h->FindBin(fallbackXmax)+1, h->GetNbinsX());

    fYield  = sg - nsgbins*(bg/nbgbins);
    fYieldE = (sg > 1? TMath::Sqrt(sg): 1.);
  }

  pa->DrawArrow(xmin, 0.4*ypeak, xmin, 0.);
  pa->DrawArrow(xmax, 0.4*ypeak, xmax, 0.);
  cout << "S = " << fYield << " +/- " << fYieldE << " background: " << (aintegral-fintegral)
       << " S/B = " << fYield/(aintegral-fintegral)
       << endl;
  tl->DrawLatexNDC(0.2, 0.92, Form("N_{sig} = %5.0f #pm %5.0f", fYield, fYieldE));

}

// ----------------------------------------------------------------------
void plotFake::playPhi(string cuts, string name) {
  string dsname = "fakeData";
  string dir = "candAnaFake333";
  TTree *T = getTree(dsname, dir, "events");
  TH1D *h = (TH1D*)gDirectory->Get("h1");
  if (h) delete h;

  TH1D *h1 = new TH1D("h1", "", 100, 1.00, 1.05); h1->Sumw2();
  h1->SetNdivisions(505, "X");
  T->Draw("m>>h1", cuts.c_str());
  fitPhi(h1);
  cout << cuts << endl;
  savePad(Form("playPhi-mass-%s.pdf", name.c_str()));
}


// ----------------------------------------------------------------------
void plotFake::fitPhi(TH1D *h) {
  fIF->fVerbose = false;
  fIF->resetLimits();
  fIF->limitPar(1, 1.01, 1.03);
  fIF->limitPar(2, 0.002, 0.005);
  fIF->limitPar(3, 0.001, 0.150);
  fIF->limitPar(4, 0.0051, 0.030);

  TF1 *f1 = fIF->phiKK(h);
  h->SetMinimum(0.);
  TFitResultPtr r =  h->Fit(f1, "lsq", "e");

  TF1 *f2 = fIF->argus(h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1));
  f2->SetLineColor(kBlue);
  f2->SetLineStyle(kDashed);
  f2->SetParameter(0, f1->GetParameter(5));
  f2->SetParameter(1, f1->GetParameter(6));
  f2->SetParameter(2, f1->GetParameter(7));
  f2->Draw("same");

  double bwidth = h->GetBinWidth(h->FindBin(1));
  double ypeak = f1->GetParameter(0);
  double xpeak = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);
  double nSig(3.0);
  double xmin = xpeak - nSig*sigma;
  double xmax = xpeak + nSig*sigma;
  double aintegral  = f1->Integral(xmin, xmax)/bwidth;
  // -- (slight?) overestimate of signal integral error by including the background
  double fintegralE = f1->IntegralError(xmin, xmax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/bwidth;
  // -- now set constant of bg to zero and integrate over signal only
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  double fintegral  = f1->Integral(xmin, xmax)/bwidth;
  if (fintegral > 0) {
    double err1 = f1->GetParError(0)/f1->GetParameter(0);
    double err2 = f1->GetParError(2)/f1->GetParameter(2);
    double errT = TMath::Sqrt(err1*err1 + err2*err2);
    fYieldE = errT*fintegral;
    fYieldE = fintegralE;
    fYield  = fintegral;
  } else {
    double fallbackXmin(1.000);
    double fallbackXmax(1.050);
    int nsgbins = h->FindBin(fallbackXmax) - h->FindBin(fallbackXmin) + 1;
    double sg = h->Integral(h->FindBin(fallbackXmin), h->FindBin(fallbackXmax));
    int nbgbins = h->FindBin(fallbackXmax) - 1 - 1 + 1;
    nbgbins    += h->GetNbinsX() - (h->FindBin(fallbackXmax) + 1) + 1;
    double bg = h->Integral(1, h->FindBin(fallbackXmin)-1) +  h->Integral(h->FindBin(fallbackXmax)+1, h->GetNbinsX());

    fYield  = sg - nsgbins*(bg/nbgbins);
    fYieldE = (sg > 1? TMath::Sqrt(sg): 1.);
  }

  pa->DrawArrow(xmin, ypeak, xmin, 0.);
  pa->DrawArrow(xmax, ypeak, xmax, 0.);
  cout << "S = " << fYield << " +/- " << fYieldE << " background: " << (aintegral-fintegral)
       << " S/B = " << fYield/(aintegral-fintegral)
       << endl;
  tl->DrawLatexNDC(0.2, 0.92, Form("N_{sig} = %5.0f #pm %5.0f", fYield, fYieldE));

}


// ----------------------------------------------------------------------
void plotFake::playLambda(string cuts, string name) {
  string dsname = "fakeData";
  string dir = "candAnaFake3122";
  TTree *T = getTree(dsname, dir, "events");
  TH1D *h = (TH1D*)gDirectory->Get("h1");
  if (h) delete h;

  TH1D *h1 = new TH1D("h1", "", 100, 1.105, 1.125); h1->Sumw2();
  h1->SetNdivisions(505, "X");
  h1->SetMinimum(0.);
  T->Draw("m>>h1", cuts.c_str());
  fitLambda(h1);
  cout << cuts << endl;

  savePad(Form("playLambda-mass-%s.pdf", name.c_str()));

}



// ----------------------------------------------------------------------
void plotFake::fitLambda(TH1D *h) {

  fIF->fVerbose = true;
  fIF->fVerbose = false;
  fIF->resetLimits();
  fIF->limitPar(0, 0, 1e8);
  fIF->limitPar(1, 1.110, 1.120);
  fIF->limitPar(2, 0.001, 0.003);
  fIF->limitPar(3, -0.1,   0.4);
  fIF->limitPar(4, 0.0031, 0.01);
  TF1* f1 = fIF->pol1gauss2c(h, 1.1158, 0.001);
  //  fIF->dumpParameters(f1);
  TFitResultPtr r = h->Fit(f1, "lsq", "e");
  double bwidth = h->GetBinWidth(h->FindBin(1));
  double ypeak = f1->GetParameter(0);
  double xpeak = f1->GetParameter(1);
  double sigma = f1->GetParameter(2);
  double nSig(5.0);
  double xmin = xpeak - nSig*sigma;
  double xmax = xpeak + nSig*sigma;
  xmin = 1.112;
  xmax = 1.120;
  double aintegral  = f1->Integral(xmin, xmax)/bwidth;
  // -- (slight?) overestimate of signal integral error by including the background
  double fintegralE = f1->IntegralError(xmin, xmax, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/bwidth;

  // -- draw background
  TF1* fb = fIF->pol1(h->GetBinLowEdge(1), h->GetBinLowEdge(h->GetNbinsX()+1));
  fb->SetParameter(0, f1->GetParameter(5));
  fb->SetParameter(1, f1->GetParameter(6));
  fb->SetLineStyle(kDashed);
  fb->Draw("same");

  // -- now set constant of pol0 to zero and integrate over signal only
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);

  double fintegral  = f1->Integral(xmin, xmax)/bwidth;
  if (fintegral > 0) {
    double err1 = f1->GetParError(0)/f1->GetParameter(0);
    double err2 = f1->GetParError(2)/f1->GetParameter(2);
    double errT = TMath::Sqrt(err1*err1 + err2*err2);
    fYieldE = errT*fintegral;
    fYieldE = fintegralE;
    fYield  = fintegral;
  } else {
    double fallbackXmin(1.110);
    double fallbackXmax(1.120);
    int nsgbins = h->FindBin(fallbackXmax) - h->FindBin(fallbackXmin) + 1;
    double sg = h->Integral(h->FindBin(fallbackXmin), h->FindBin(fallbackXmax));
    int nbgbins = h->FindBin(fallbackXmax) - 1 - 1 + 1;
    nbgbins    += h->GetNbinsX() - (h->FindBin(fallbackXmax) + 1) + 1;
    double bg = h->Integral(1, h->FindBin(fallbackXmin)-1) +  h->Integral(h->FindBin(fallbackXmax)+1, h->GetNbinsX());

    fYield  = sg - nsgbins*(bg/nbgbins);
    fYieldE = (sg > 1? TMath::Sqrt(sg): 1.);
  }

  pa->DrawArrow(xmin, 0.7*ypeak, xmin, 0.);
  pa->DrawArrow(xmax, 0.7*ypeak, xmax, 0.);
  cout << "S = " << fYield << " +/- " << fYieldE << " background: " << (aintegral-fintegral)
       << " S/B = " << fYield/(aintegral-fintegral)
       << endl;
  tl->DrawLatexNDC(0.2, 0.92, Form("N_{sig} = %5.0f #pm %5.0f", fYield, fYieldE));

}


// ----------------------------------------------------------------------
void plotFake::bookHist(int mode) {

}



// ----------------------------------------------------------------------
void plotFake::loopFunction1() {

  fGoodCand = true;
  if (fMode == FAKELAMBDA) {
    if (fCandPvIp > 0.01) fGoodCand = false;
    if (fCandPvIpS > 2)   fGoodCand = false;
    if (fCandFLSxy < 15)  fGoodCand = false;
    if (fCandFLS3d < 15)  fGoodCand = false;
  }


  double mass = fCandM;

  if (fIsMC) {
    mass = SIGBOXMIN + 0.5 * (SIGBOXMAX - SIGBOXMIN);
  }

  string mapname("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");

  //  cout << "cand m = " << fCandM << endl;

  for (int i = 0; i < fFakeNtrk; ++i) {
    if (TMath::Abs(fFakeEta[i]) < 0.8) {
      fChan = 0;
    } else if (0.8 < TMath::Abs(fFakeEta[i])  && TMath::Abs(fFakeEta[i]) < 1.3) {
      fChan = 1;
    } else {
      fChan = 2;
    }

    if (fIsMC) {
      if (fMode == FAKEKS) {
	fGoodCand = (fFakeId[i] == 211);
      } else if (fMode == FAKEPHI) {
	fGoodCand = (fFakeId[i] == 321);
      } else if (fMode == FAKELAMBDA) {
	fGoodCand = (fFakeId[i] == 2212);
      } else if (fMode == FAKEPSI) {
	fGoodCand = (fFakeId[i] == 13);
      }
    }

    mapname = fChannelSample[fChan];

    fGlobalMuon  = (fFakeGm[i]>0);
    fGoodPt      = (fFakePt[i] > 4.);
    fGoodDtrig   = (fFakeDtrig[i] > 0.01);
    fGoodDmuon   = (fFakeDmuon[i] > 0.5);
    if (fIsMC) {
      fTIS       = true;
      fGoodDtrig = true; // FIXME?
    }

    fGood     = fGoodCand && fGoodPt;
    fGoodFake = fGood && fGlobalMuon;

    fGoodTIS         = fTIS       && fGoodCand && fGoodPt;
    fGoodTISFake     = fGoodTIS   && fGlobalMuon;

    fGoodTISDT       = fGoodTIS   && fGoodDtrig;
    fGoodTISDTFake   = fGoodTISDT && fGlobalMuon;

    fGoodTISDTDM     = fGoodTISDT   && fGoodDmuon;
    fGoodTISDTDMFake = fGoodTISDTDM && fGlobalMuon;

    fAnaCuts.update();


    fAdMap[mapname]->fpAllPt->fill(fFakePt[i], mass);
    fAdMap[mapname]->fpAllEta->fill(fFakeEta[i], mass);

    fAdMap[mapname]->fpFakePt->fill(fFakePt[i], mass);
    fAdMap[mapname]->fpFakeEta->fill(fFakeEta[i], mass);

    fAdMap[mapname]->fpFakeTisAllPt->fill(fFakePt[i], mass);
    fAdMap[mapname]->fpFakeTisAllEta->fill(fFakeEta[i], mass);
    fAdMap[mapname]->fpFakeTisFakePt->fill(fFakePt[i], mass);
    fAdMap[mapname]->fpFakeTisFakeEta->fill(fFakeEta[i], mass);

    fAdMap[mapname]->fpFakeTisDtAllEta->fill(fFakeEta[i], mass);
    fAdMap[mapname]->fpFakeTisDtAllPt->fill(fFakePt[i], mass);
    fAdMap[mapname]->fpFakeTisDtFakeEta->fill(fFakeEta[i], mass);
    fAdMap[mapname]->fpFakeTisDtFakePt->fill(fFakePt[i], mass);

    fAdMap[mapname]->fpFakeTisDtDmAllEta->fill(fFakeEta[i], mass);
    fAdMap[mapname]->fpFakeTisDtDmAllPt->fill(fFakePt[i], mass);
    fAdMap[mapname]->fpFakeTisDtDmFakeEta->fill(fFakeEta[i], mass);
    fAdMap[mapname]->fpFakeTisDtDmFakePt->fill(fFakePt[i], mass);

    fAdMap[mapname]->fpFakeBdt->fill(fFakeBdt[i], mass);
    fAdMap[mapname]->fpFakeTip->fill(fFakeTip[i], mass);
    fAdMap[mapname]->fpFakeLip->fill(fFakeLip[i], mass);

    fAdMap[mapname]->fpFakeInnerChi2->fill(fFakeInnerChi2[i], mass);
    fAdMap[mapname]->fpFakeOuterChi2->fill(fFakeOuterChi2[i], mass);

    fAdMap[mapname]->fpFakeChi2LocalPosition->fill(fFakeChi2LocalPosition[i], mass);
    fAdMap[mapname]->fpFakeChi2LocalMomentum->fill(fFakeChi2LocalMomentum[i], mass);
    fAdMap[mapname]->fpFakeStaTrkMult->fill(fFakeStaTrkMult[i], mass);
    fAdMap[mapname]->fpFakeTmTrkMult->fill(fFakeTmTrkMult[i], mass);

    fAdMap[mapname]->fpFakeDeltaR->fill(fFakeDeltaR[i], mass);
    fAdMap[mapname]->fpFakeItrkValidFraction->fill(fFakeItrkValidFraction[i], mass);
    fAdMap[mapname]->fpFakeSegmentComp->fill(fFakeSegmentComp[i], mass);
    fAdMap[mapname]->fpFakeGtrkNormChi2->fill(fFakeGtrkNormChi2[i], mass);
    fAdMap[mapname]->fpFakeDzRef->fill(fFakeDzRef[i], mass);
    fAdMap[mapname]->fpFakeDxyRef->fill(fFakeDxyRef[i], mass);
    fAdMap[mapname]->fpFakeGtrkProb->fill(fFakeGtrkProb[i], mass);
    fAdMap[mapname]->fpFakeMuonChi2->fill(fFakeMuonChi2[i], mass);
    fAdMap[mapname]->fpFakeGlbKinkFinder->fill(TMath::Log10(fFakeGlbKinkFinder[i]), mass);
    fAdMap[mapname]->fpFakeStaRelChi2->fill(fFakeStaRelChi2[i], mass);
    fAdMap[mapname]->fpFakeTrkRelChi2->fill(fFakeTrkRelChi2[i], mass);
    fAdMap[mapname]->fpFakeGlbDeltaEtaPhi->fill(fFakeGlbDeltaEtaPhi[i], mass);
    fAdMap[mapname]->fpFakeTimeInOut->fill(fFakeTimeInOut[i], mass);
    fAdMap[mapname]->fpFakeTimeInOutE->fill(fFakeTimeInOutE[i], mass);
    if (fFakeTimeInOutE[i] > 0.) fAdMap[mapname]->fpFakeTimeInOutS->fill(fFakeTimeInOut[i]/fFakeTimeInOutE[i], mass);
    fAdMap[mapname]->fpFakeNvalidMuonHits->fill(fFakeNvalidMuonHits[i], mass);
    fAdMap[mapname]->fpFakeNmatchedStations->fill(fFakeNmatchedStations[i], mass);
    fAdMap[mapname]->fpFakeLayersWithHits->fill(fFakeLayersWithHits[i], mass);
    fAdMap[mapname]->fpFakeNumberOfValidTrkHits->fill(fFakeNumberOfValidTrkHits[i], mass);
    fAdMap[mapname]->fpFakeNumberOfLostTrkHits->fill(fFakeNumberOfLostTrkHits[i], mass);
    fAdMap[mapname]->fpFakeNumberOfValidPixHits->fill(fFakeNumberOfValidPixHits[i], mass);
    fAdMap[mapname]->fpFakeRPChits1->fill(fFakeRPChits1[i], mass);
    fAdMap[mapname]->fpFakeRPChits2->fill(fFakeRPChits2[i], mass);
    fAdMap[mapname]->fpFakeRPChits3->fill(fFakeRPChits3[i], mass);
    fAdMap[mapname]->fpFakeRPChits4->fill(fFakeRPChits4[i], mass);
    fAdMap[mapname]->fpFakeCombHits->fill(fFakeCombHits[i], mass);

  }


}


// ----------------------------------------------------------------------
void plotFake::analysis() {


}


// ----------------------------------------------------------------------
void plotFake::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotFake::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotFake::*pF)(void);
  if (ifunc == 1) pF = &plotFake::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotFake::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotFake::loadFiles> Loading files listed in " << files << endl;

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

    //    cout << "stype: ->" << stype << "<-" << endl;

    TFile *pF(0);
    dataset *ds(0);
    if (string::npos != stype.find("data")) {
      // -- DATA
      pF = loadFile(sfile);

      ds = new dataset();
      if (string::npos != stype.find("XXXX")) {
	ds->fSize = 1;
	ds->fWidth = 2;

	ds->fLcolor = ds->fColor;
	ds->fFcolor = ds->fColor;
	ds->fName   = sdecay;
	ds->fFullName = sname;
      }
    } else if (string::npos != stype.find("mc")) {
      // -- MC
      pF = loadFile(sfile);
      //      cout << "  " << sfile << ": " << pF << endl;

      dataset *ds = new dataset();
      ds->fSize = 1;
      ds->fWidth = 2;

      if (string::npos != stype.find("XXX")) {
        sname = "XXX";
        sdecay = "bu2jpsik";
	ds->fColor = kBlue-7;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }
    }

    if (sname != "nada") {
      //      cout << "  inserting as " << sname << " and " << sdecay << endl;
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
  cout << "------------------------------------------------------------------------------------------" << endl;
  cout << Form("   %30s: %20s: ", "Dataset name", "Decay mode name") << "Filename:" << endl;
  cout << "------------------------------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    // cout << it->first << endl;
    // cout << it->second->fName << endl;
    // cout << it->second->fF->GetName() << endl;
    cout << Form("%2d %30s: %20s: ", cnt, it->first.c_str(), it->second->fName.c_str()) << it->second->fF->GetName() << endl;
    ++cnt;
  }
  cout << "------------------------------------------------------------------------------------------" << endl;

}


// ----------------------------------------------------------------------
void plotFake::setupTree(TTree *t) {

  t->SetBranchAddress("m",       &fCandM);
  t->SetBranchAddress("pvip",    &fCandPvIp);
  t->SetBranchAddress("pvips",   &fCandPvIpS);
  t->SetBranchAddress("chi2dof", &fCandChi2Dof);
  t->SetBranchAddress("fls3d",   &fCandFLS3d);
  t->SetBranchAddress("fl3d",    &fCandFL3d);
  t->SetBranchAddress("flxy",    &fCandFLxy);
  t->SetBranchAddress("fl3dE",   &fCandFL3dE);
  t->SetBranchAddress("flsxy",   &fCandFLSxy);
  t->SetBranchAddress("maxdoca", &fCandDoca);
  t->SetBranchAddress("tis",     &fTIS);
  t->SetBranchAddress("cowboy",  &fCowboy);

  t->SetBranchAddress("ntrk",    &fFakeNtrk);
  t->SetBranchAddress("id",      fFakeId);
  t->SetBranchAddress("q",       fFakeQ);
  t->SetBranchAddress("gm",      fFakeGm);
  t->SetBranchAddress("pt",      fFakePt);
  t->SetBranchAddress("eta",     fFakeEta);
  t->SetBranchAddress("phi",     fFakePhi);
  t->SetBranchAddress("dtrig",   fFakeDtrig);
  t->SetBranchAddress("dmuon",   fFakeDmuon);
  t->SetBranchAddress("bdt",     fFakeBdt);

  t->SetBranchAddress("tip", fFakeTip);
  t->SetBranchAddress("lip", fFakeLip);

  t->SetBranchAddress("innerchi2", fFakeInnerChi2);
  t->SetBranchAddress("outerchi2", fFakeOuterChi2);
  t->SetBranchAddress("chi2localposition", fFakeChi2LocalPosition);
  t->SetBranchAddress("chi2localmomentum", fFakeChi2LocalMomentum);
  t->SetBranchAddress("statrkmult", fFakeStaTrkMult);
  t->SetBranchAddress("tmtrkmult", fFakeTmTrkMult);
  t->SetBranchAddress("deltar", fFakeDeltaR);
  t->SetBranchAddress("itrkvalidfraction", fFakeItrkValidFraction);
  t->SetBranchAddress("segmentcomp", fFakeSegmentComp);
  t->SetBranchAddress("gtrknormchi2", fFakeGtrkNormChi2);
  t->SetBranchAddress("dzref", fFakeDzRef);
  t->SetBranchAddress("dxyref", fFakeDxyRef);
  t->SetBranchAddress("gtrktailprob", fFakeGtrkProb);
  t->SetBranchAddress("numberofvalidtrkhits", fFakeNumberOfValidTrkHits);
  t->SetBranchAddress("numberoflosttrkhits", fFakeNumberOfLostTrkHits);
  t->SetBranchAddress("muonchi2", fFakeMuonChi2);
  t->SetBranchAddress("glbkinkfinder", fFakeGlbKinkFinder);
  t->SetBranchAddress("starelchi2", fFakeStaRelChi2);
  t->SetBranchAddress("trkrelchi2", fFakeTrkRelChi2);
  t->SetBranchAddress("glbdeltaetaphi", fFakeGlbDeltaEtaPhi);
  t->SetBranchAddress("timeinout", fFakeTimeInOut);
  t->SetBranchAddress("timeinoute", fFakeTimeInOutE);
  t->SetBranchAddress("nvalidmuonhits", fFakeNvalidMuonHits);
  t->SetBranchAddress("nvalidmuonhits", fFakeNvalidMuonHits);
  t->SetBranchAddress("nmatchedstations", fFakeNmatchedStations);
  t->SetBranchAddress("layerswithhits", fFakeLayersWithHits);
  t->SetBranchAddress("numberofvalidpixhits", fFakeNumberOfValidPixHits);
  t->SetBranchAddress("rpchits1", fFakeRPChits1);
  t->SetBranchAddress("rpchits2", fFakeRPChits2);
  t->SetBranchAddress("rpchits3", fFakeRPChits3);
  t->SetBranchAddress("rpchits4", fFakeRPChits4);
  t->SetBranchAddress("mudethitscomb", fFakeCombHits);



}
