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
plotFake::plotFake(string dir, string files, string cuts, string setup, int year): plotClass(dir, files, cuts, setup, year) {
  plotClass::loadFiles(files);
  plotFake::loadFiles(files);

  changeSetup(dir, "plotFake", setup);

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  fChan = 0;

  fChannelList.clear();
  // -- four eta regions!
  for (unsigned int i = 0; i < 4; ++i) {
    fChannelList.push_back(Form("%d", i));
  }

  fillDoList("all");

  fCncCuts.clear();
  fCncCuts.addCut("GoodCand", "good cand", fGoodCand);
  fCncCuts.addCut("GoodPt", "good pt", fGoodPt);
  fCncCuts.addCut("GlobalMuon", "global muon ID", fGlobalMuon);

  fCncCuts.addCut("Good", "all, good pt", fGood);
  fCncCuts.addCut("GoodFake", "fake, good pt", fGoodFake);

  fCncCuts.addCut("TIS", "all, triggered independently of signal", fGoodTIS);
  fCncCuts.addCut("TISFAKE", "fake, triggered independently of signal", fGoodTISFake);

  fCncCuts.addCut("TISDT", "all, triggered independently of signal, dist(trigger)", fGoodTISDT);
  fCncCuts.addCut("TISDTFAKE", "fake, triggered independently of signal, dist(trigger)", fGoodTISDTFake);

  fCncCuts.addCut("TISDTDM", "all, triggered independently of signal, dist(trigger), dist(muon)", fGoodTISDTDM);
  fCncCuts.addCut("TISDTDMFAKE", "fake, triggered independently of signal, dist(trigger), dist(muon)", fGoodTISDTDMFake);

  fCncCuts.dumpAll();
}


// ----------------------------------------------------------------------
plotFake::~plotFake() {

}

// ----------------------------------------------------------------------
void plotFake::fillDoList(string what) {

  fDoList.clear();

  if (what == "fakerate") {
    fDoList.push_back("FakeTisDtDmFakeEta");
    fDoList.push_back("FakeTisDtDmFakePt");
    fDoList.push_back("FakeTisDtDmAllEta");
    fDoList.push_back("FakeTisDtDmAllPt");
    fDoList.push_back("FakeTisDtFakeEta");
    fDoList.push_back("FakeTisDtFakePt");
    fDoList.push_back("FakeTisDtAllEta");
    fDoList.push_back("FakeTisDtAllPt");
    fDoList.push_back("FakeTisFakeEta");
    fDoList.push_back("FakeTisFakePt");
    fDoList.push_back("FakeTisAllEta");
    fDoList.push_back("FakeTisAllPt");
    fDoList.push_back("FakePt");
    fDoList.push_back("FakeEta");
    fDoList.push_back("AllEta");
    fDoList.push_back("AllPt");
  }

  if (what == "all") {
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
    fDoList.push_back("FakeQprod");
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
  }


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

  if (what == "dbx1") {
    makeOverlay("fakeData_ks", "fakeMc_ks", "Cu");
    makeOverlay("fakeData_phi", "fakeMc_phi", "Cu");
    makeOverlay("fakeData_lambda", "fakeMc_lambda", "Cu");
    makeOverlay("fakeData_psi", "fakeMc_psi", "Cu");
    for (unsigned int ic = 0; ic < fNchan; ++ic) {
      fTEX << formatTex(fCuts[ic]->muonbdt, Form("%s:muonidBdtCut_responseCut_chan%i:val", fSuffix.c_str(), ic), 3) << endl;
    }
  }

  if (what == "all" || string::npos != what.find("sample")) {
    init();
    if ((what == "all") || (what == "sample") || ((string::npos != what.find("sample") && string::npos != what.find("ks")))) {
      makeSample("fakeData", "ks");
      makeSample("fakeMc", "ks");
    }
    if ((what == "all") || (what == "sample") || ((string::npos != what.find("sample") && string::npos != what.find("psi")))) {
      makeSample("fakeData", "psi", 5e6);
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
    system(Form("/bin/rm -f %s/%d%s*sbsctrl_ad*_fake*.pdf", fDirectory.c_str(), fYear, fSetup.c_str()));
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

    for (unsigned int ic = 0; ic < fNchan; ++ic) {
      fTEX << formatTex(fCuts[ic]->muonbdt, Form("%s:muonidBdtCut_responseCut_chan%i:val", fSuffix.c_str(), ic), 3) << endl;
    }
  }

  if (what == "all" || string::npos != what.find("fakerate") || string::npos != what.find("plot")) {
    fillDoList("fakerate");
    if ((what == "all") || (what == "fakerate") || (what == "plot") || (string::npos != what.find("ks"))) {
      fakeRate("fakeData_ks", "fakeMc_ks", "FakePt", "AllPt", 0.1);
      fakeRate("fakeData_ks", "fakeMc_ks", "FakeEta", "AllEta", 0.1);

      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisFakePt", "FakeTisAllPt");
      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisFakeEta", "FakeTisAllEta");

      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisDtFakePt", "FakeTisDtAllPt");
      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisDtFakeEta", "FakeTisDtAllEta");

      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisDtDmFakePt", "FakeTisDtDmAllPt");
      fakeRate("fakeData_ks", "fakeMc_ks", "FakeTisDtDmFakeEta", "FakeTisDtDmAllEta");
    }

    if ((what == "all") || (what == "fakerate") || (what == "plot") || (string::npos != what.find("phi"))) {
      fakeRate("fakeData_phi", "fakeMc_phi", "FakeTisDtDmFakePt", "FakeTisDtDmAllPt");
      fakeRate("fakeData_phi", "fakeMc_phi", "FakeTisDtDmFakeEta", "FakeTisDtDmAllEta");
    }

    if ((what == "all") || (what == "fakerate") || (what == "plot") || (string::npos != what.find("lambda"))) {
      fakeRate("fakeData_lambda", "fakeMc_lambda", "FakeTisDtDmFakePt", "FakeTisDtDmAllPt");
      fakeRate("fakeData_lambda", "fakeMc_lambda", "FakeTisDtDmFakeEta", "FakeTisDtDmAllEta");
    }

    if ((what == "all") || (what == "fakerate") || (what == "plot") || (string::npos != what.find("psi"))) {
      fakeRate("fakeData_psi", "fakeMc_psi", "FakeTisDtDmFakePt", "FakeTisDtDmAllPt");
      fakeRate("fakeData_psi", "fakeMc_psi", "FakeTisDtDmFakeEta", "FakeTisDtDmAllEta");
    }

  }

  if ((what == "all") || string::npos != what.find("pidtables") || string::npos != what.find("plot")) {
    if (2016 == fYear) {
      mkPidTables("bmm4-25");
      plotPidTables("");
    }

    if (2012 == fYear) {
      mkPidTables("bmm4-42");
    }

    if (2011 == fYear) {
      mkPidTables("bmm4-42");
    }
    plotPidTables("");
  }
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
	setTitles(h1, "m_{#pi#pi} #it{[GeV]}", "Candidates", 0.05, 1.1, 1.8);
      } else if (string::npos != modes[im].find("phi")) {
	setTitles(h1, "m_{KK} #it{[GeV]}", "Candidates", 0.05, 1.1, 1.8);
      } else if (string::npos != modes[im].find("lambda")) {
	setTitles(h1, "m_{p#pi} #it{[GeV]}", "Candidates", 0.05, 1.1, 1.8);
      } else if (string::npos != modes[im].find("psi")) {
	setTitles(h1, "m_{#mu#mu} #it{[GeV]}", "Candidates", 0.05, 1.1, 1.8);
      }
      shrinkPad(0.15, 0.2);
      if (h1) h1->Draw();
      savePad(Form("%s-fakemass_ad%d_%s_%s.pdf", fSetup.c_str(), ic, modes[im].c_str(), varname.c_str()));
    }
  }
}


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
  for (int ic = 0; ic < fChannelList.size(); ++ic) {
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
  }

  if (string::npos != sample.find("psi")) {
    fMode = FAKEPSI;
    dir = "candAnaFake443";
  }

  if (string::npos != sample.find("phi")) {
    fMode = FAKEPHI;
    dir = "candAnaFake333";
  }

  if (string::npos != sample.find("lambda")) {
    fMode = FAKELAMBDA;
    dir = "candAnaFake3122";
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
  fCds = fDS[dataset];
  loopOverTree(t, 1, nevents, nstart);

  fHistFile->Write();
  fHistFile->Close();

}


// ----------------------------------------------------------------------
void plotFake::makeOverlay(string what1, string what2, string selection, string what) {

  cout << "fHistFileName: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  string label1("bla");
  if (string::npos != what1.find("psi")) {
    label1 = "muons";
  } else if (string::npos != what1.find("ks")) {
    label1 = "pions";
  } else if (string::npos != what1.find("phi")) {
    label1 = "kaons";
  } else if (string::npos != what1.find("lambda")) {
    label1 = "protons";
  }


  bool restricted = (what != "");
  c0->cd();
  c0->Clear();
  c0->SetCanvasSize(700, 700);
  shrinkPad(0.15, 0.18, 0.1);
  for (unsigned int ic = 0; ic < fChannelList.size(); ++ic) {
    //  for (unsigned int i = 0; i < 1; ++i) {
    cout << "===> sbsDistributions " << Form("ad%s_%s", fChannelList[ic].c_str(), what1.c_str()) << endl;
    sbsDistributions(Form("ad%s_%s", fChannelList[ic].c_str(), what1.c_str()), selection, what);
    sbsDistributions(Form("ad%s_%s", fChannelList[ic].c_str(), what2.c_str()), selection, what);

    for (unsigned int id = 0; id < fDoList.size(); ++id) {
      if (restricted) {
	if (string::npos == fDoList[id].find(what)) continue;
      }
      c0->cd();
      overlay(Form("sbs_ad%s_%s_%s%s", fChannelList[ic].c_str(), what1.c_str(), fDoList[id].c_str(), selection.c_str()),
	      Form("sbs_ad%s_%s_%s%s", fChannelList[ic].c_str(), what2.c_str(), fDoList[id].c_str(), selection.c_str())
	      );
      savePad(Form("%s-fakeoverlay_ad%s_%s_ad%s_%s_%s-%s.pdf",
		   fSetup.c_str(),
		   fChannelList[ic].c_str(), what1.c_str(),
		   fChannelList[ic].c_str(), what2.c_str(),
		   fDoList[id].c_str(), selection.c_str()));

      // -- determine muon id/misid systematics // FAIL: I don't think this works in any form for the misid rate!!!!
      if (string::npos != fDoList[id].find("FakeBdt")) {
	TH1D *h1 = (TH1D*)gDirectory->Get(Form("sbs_ad%s_%s_%s%s", fChannelList[ic].c_str(), what1.c_str(), fDoList[id].c_str(), selection.c_str()));
	TH1D *h2 = (TH1D*)gDirectory->Get(Form("sbs_ad%s_%s_%s%s", fChannelList[ic].c_str(), what2.c_str(), fDoList[id].c_str(), selection.c_str()));
	TH1D *hd = (TH1D*)h1->Clone("hd"); hd->Reset();
	int nmax = h1->GetNbinsX();
	double itot = h1->Integral();
	double imax(0.), icut(0.);
	for (int ib = 1; ib <= h1->GetNbinsX(); ++ib) {
	  double int1 = h1->Integral(ib, nmax);
	  if (int1 < 0.) int1 *= -1.;
	  double int2 = h2->Integral(ib, nmax);
	  if (int2 < 0.) int2 *= -1.;
	  double iave = 0.5*(int1+int2);
	  double err1 = dEff(static_cast<int>(int1), static_cast<int>(itot));
	  double err2 = dEff(static_cast<int>(int2), static_cast<int>(itot));
	  double sgn1 = TMath::Abs(int1/itot - int2/itot)/TMath::Sqrt(err1*err1 + err2*err2);
	  // -- this would be the absolute efficiency error
	  hd->SetBinContent(ib, (int1-int2)/itot);
	  if (TMath::Abs((int1-int2)/itot) > imax) imax = TMath::Abs((int1-int2)/itot);
	  if (ib == h1->FindBin(fCuts[0]->muonbdt)) icut = TMath::Abs((int1-int2)/itot);
	  // -- this is the relative efficiency error (turns out to be too large to be practical, for misid rates!)
	  // if (iave > 0.) {
	  //   hd->SetBinContent(ib, (int1-int2)/iave);
	  // } else {
	  //   hd->SetBinContent(ib, 0.);
	  // }
	  // if (TMath::Abs((int1-int2)/itot) > imax) imax = TMath::Abs((int1-int2)/iave);
	  // if (ib == h1->FindBin(fCuts[0]->muonbdt)) icut = TMath::Abs((int1-int2)/iave);
	  // -- difference-significance scaled control sample error
	  // double error(0.);
	  // if (int1 > 0) {
	  //   error = err1/(int1/itot); // relative error on data efficiency
	  // } else {
	  //   error = 0.;
	  // }
	  // if (sgn1 > 1.) {
	  //   error = error * sgn1;
	  // }
	  // hd->SetBinContent(ib, error);
	  // if (error > imax) imax = error;
	  // if (ib == h1->FindBin(fCuts[0]->muonbdt)) icut = error;
	  hd->SetBinError(ib, TMath::Sqrt(err1*err1 + err2*err2));
	  cout << "bin " << ib << " center = " << h1->GetBinCenter(ib)
	       << " itot = " << itot
	       << " eps1 = " << int1/itot << "+/-" << err1
	       << " eps2 = " << int2/itot << "+/-" << err2
	       << " diff = " << (int1-int2)/itot
	       << endl;
	}
	hd->SetMinimum(-0.50);
	hd->SetMaximum(0.50);
	hd->GetXaxis()->SetTitle("BDT > ");
	hd->GetYaxis()->SetTitle("(#varepsilon(cut; data) - #varepsilon(cut; MC)");
	// hd->GetYaxis()->SetTitle("(#varepsilon(cut; data) - #varepsilon(cut; MC)/mean");
	// hd->GetYaxis()->SetTitle("#Delta#varepsilon(cut; data) (#times f)");
	hd->SetTitleOffset(1.5, "y");
	hd->SetNdivisions(505, "X");
	hd->SetNdivisions(505, "Y");
	gPad->SetGridy(true);
	hd->Draw();
	tl->SetTextSize(0.04);
	tl->DrawLatexNDC(0.20, 0.92, Form("#Delta max: %4.3f cut: %4.3f", imax, icut));
	tl->DrawLatexNDC(0.73, 0.92, Form("%s chan %d", label1.c_str(), ic));
	tl->DrawLatexNDC(0.20, 0.84, Form("(cut BDT > %4.3f", fCuts[0]->muonbdt));

	double err(0.);
	if (imax < 0.02) err = 0.02;
	else if (imax < 0.05) err = 0.05;
	else if (imax < 0.1) err = 0.1;
	else if (imax < 0.15) err = 0.15;
	else if (imax < 0.20) err = 0.20;
	else if (imax < 0.25) err = 0.25;
	else if (imax < 0.30) err = 0.30;
	else  err = 0.50;
	fTEX << formatTex(imax, Form("%s:muonidBdtMax_%s_chan%i:val", fSuffix.c_str(), label1.c_str(), ic), 3) << endl;
	fTEX << formatTex(err, Form("%s:muonidBdtMax_%s_chan%i:err", fSuffix.c_str(), label1.c_str(), ic), 3) << endl;

	if (icut < 0.02) err = 0.02;
	else if (icut < 0.05) err = 0.05;
	else if (icut < 0.1) err = 0.1;
	else if (icut < 0.15) err = 0.15;
	else if (icut < 0.20) err = 0.20;
	else if (icut < 0.25) err = 0.25;
	else if (icut < 0.30) err = 0.30;
	else  err = 0.50;
	fTEX << formatTex(icut, Form("%s:muonidBdtCut_%s_chan%i:val", fSuffix.c_str(), label1.c_str(), ic), 3) << endl;
	fTEX << formatTex(err, Form("%s:muonidBdtCut_%s_chan%i:err", fSuffix.c_str(), label1.c_str(), ic), 3) << endl;

	savePad(Form("%s-systematics_ad%s_%s_ad%s_%s_%s-%s.pdf",
		     fSetup.c_str(),
		     fChannelList[ic].c_str(),
		     what1.c_str(), fChannelList[ic].c_str(),
		     what2.c_str(), fDoList[id].c_str(),
		     selection.c_str()));
	gPad->SetGridy(false);
      }

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
	BGLBOXMAX = 0.475;
	BGHBOXMIN = 0.515;
	BGHBOXMAX = 0.550;
      } else if (1 == i) {
	SIGBOXMIN = 0.480;
	SIGBOXMAX = 0.515;
	BGLBOXMIN = 0.450;
	BGLBOXMAX = 0.475;
	BGHBOXMIN = 0.520;
	BGHBOXMAX = 0.550;
      } else {
	SIGBOXMIN = 0.475;
	SIGBOXMAX = 0.520;
	BGLBOXMIN = 0.450;
	BGLBOXMAX = 0.470;
	BGHBOXMIN = 0.520;
	BGHBOXMAX = 0.550;
      }
    }

    if (string::npos != fSample.find("lambda")) {
      // BGLBOXMIN = 1.095;
      // BGLBOXMAX = 1.105;
      // SIGBOXMIN = 1.110;
      // SIGBOXMAX = 1.122;
      // BGHBOXMIN = 1.130;
      // BGHBOXMAX = 1.140;
      if (0 == i) {
	BGLBOXMIN = 1.095;
	BGLBOXMAX = 1.105;
	SIGBOXMIN = 1.110;
	SIGBOXMAX = 1.122;
	BGHBOXMIN = 1.130;
	BGHBOXMAX = 1.140;
      } else if (1 == i) {
	BGLBOXMIN = 1.095;
	BGLBOXMAX = 1.105;
	SIGBOXMIN = 1.110;
	SIGBOXMAX = 1.122;
	BGHBOXMIN = 1.130;
	BGHBOXMAX = 1.140;
      } else {
	BGLBOXMIN = 1.095;
	BGLBOXMAX = 1.105;
	SIGBOXMIN = 1.110;
	SIGBOXMAX = 1.122;
	BGHBOXMIN = 1.130;
	BGHBOXMAX = 1.140;
      }
    }

    if (string::npos != fSample.find("psi")) {
      // BGLBOXMIN = 2.800;
      // BGLBOXMAX = 2.900;
      // SIGBOXMIN = 3.050;
      // SIGBOXMAX = 3.150;
      // BGHBOXMIN = 3.200;
      // BGHBOXMAX = 3.300;
      if (0 == i) {
	BGLBOXMIN = 2.920;
	BGLBOXMAX = 2.960;
	SIGBOXMIN = 3.000;
	SIGBOXMAX = 3.200;
	BGHBOXMIN = 3.240;
	BGHBOXMAX = 3.300;
      } else if (1 == i) {
	BGLBOXMIN = 2.920;
	BGLBOXMAX = 2.960;
	SIGBOXMIN = 3.000;
	SIGBOXMAX = 3.200;
	BGHBOXMIN = 3.240;
	BGHBOXMAX = 3.300;
      } else {
	BGLBOXMIN = 2.920;
	BGLBOXMAX = 2.960;
	SIGBOXMIN = 3.000;
	SIGBOXMAX = 3.200;
	BGHBOXMIN = 3.240;
	BGHBOXMAX = 3.300;
      }
    }

    if (string::npos != fSample.find("phi")) {
      // BGLBOXMIN = 0.990;
      // BGLBOXMAX = 1.005;
      // SIGBOXMIN = 1.010;
      // SIGBOXMAX = 1.030;
      // BGHBOXMIN = 1.035;
      // BGHBOXMAX = 1.045;
      if (0 == i) {
	BGLBOXMIN = 0.990;
	BGLBOXMAX = 1.005;
	SIGBOXMIN = 1.010;
	SIGBOXMAX = 1.030;
	BGHBOXMIN = 1.035;
	BGHBOXMAX = 1.045;
      } else if (1 == i) {
	BGLBOXMIN = 0.990;
	BGLBOXMAX = 1.005;
	SIGBOXMIN = 1.010;
	SIGBOXMAX = 1.030;
	BGHBOXMIN = 1.035;
	BGHBOXMAX = 1.045;
      } else {
	BGLBOXMIN = 0.990;
	BGLBOXMAX = 1.005;
	SIGBOXMIN = 1.010;
	SIGBOXMAX = 1.030;
	BGHBOXMIN = 1.035;
	BGHBOXMAX = 1.045;
      }
    }

    cout << "fSample: " << fSample << endl;
    cout << "SIG: " << SIGBOXMIN << " .. " << SIGBOXMAX << endl;
    cout << "BGL: " << BGLBOXMIN << " .. " << BGLBOXMAX << endl;
    cout << "BGH: " << BGHBOXMIN << " .. " << BGHBOXMAX << endl;

    a = new adsetFake();
    a->fpFakeEta  = bookDistribution(Form("%sFakeEta", name.c_str()), "#eta", "GoodFake", 48, -2.4, 2.4);
    a->fpFakePt   = bookDistribution(Form("%sFakePt", name.c_str()), "#it{p_{T}} [GeV]", "GoodFake", 10, 0., 20.);
    a->fpAllEta  = bookDistribution(Form("%sAllEta", name.c_str()), "#eta", "Good", 48, -2.4, 2.4);
    a->fpAllPt   = bookDistribution(Form("%sAllPt", name.c_str()), "#it{p_{T}} [GeV]", "Good", 10, 0., 20.);

    a->fpFakeBdt       = bookDistribution(Form("%sFakeBdt", name.c_str()), "BDT", "GlobalMuon", 40, -1., 1.0);
    a->fpFakeTip       = bookDistribution(Form("%sFakeTip", name.c_str()), "TIP [cm]", "GlobalMuon", 20, 0., 2.);
    a->fpFakeLip       = bookDistribution(Form("%sFakeLip", name.c_str()), "LIP [cm]", "GlobalMuon", 20, 0., 2.);
    a->fpFakeQprod     = bookDistribution(Form("%sFakeQprod", name.c_str()), "qprod", "GlobalMuon", 3, -1., 2.);
    a->fpFakeInnerChi2 = bookDistribution(Form("%sFakeInnerChi2", name.c_str()), "inner track #chi^{2}", "GlobalMuon", 20, 0., 20.);
    a->fpFakeOuterChi2 = bookDistribution(Form("%sFakeOuterChi2", name.c_str()), "outer track #chi^{2}", "GlobalMuon", 40, 0., 20.);

    a->fpFakeChi2LocalPosition = bookDistribution(Form("%sFakeChi2LocalPosition", name.c_str()), "local position #chi^{2}", "GlobalMuon", 21, 0., 102.);
    a->fpFakeChi2LocalMomentum = bookDistribution(Form("%sFakeChi2LocalMomentum", name.c_str()), "local momentum #chi^{2}", "GlobalMuon", 21, 0., 102.);
    a->fpFakeStaTrkMult = bookDistribution(Form("%sFakeStaTrkMult", name.c_str()), "STA trk multipicity", "GlobalMuon", 12, -2., 10.);
    a->fpFakeTmTrkMult = bookDistribution(Form("%sFakeTmTrkMult", name.c_str()), "TM trk multiplicity", "GlobalMuon", 20, 0., 20.);

    a->fpFakeDeltaR = bookDistribution(Form("%sFakeDeltaR", name.c_str()), "deltaR", "GlobalMuon", 20, 0., 1.0);
    a->fpFakeItrkValidFraction = bookDistribution(Form("%sFakeItrkValidFraction", name.c_str()), "inner track valid fraction", "GlobalMuon", 21, 0., 1.02);
    a->fpFakeSegmentComp = bookDistribution(Form("%sFakeSegmentComp", name.c_str()), "segment compatibility", "GlobalMuon", 21, 0., 1.02);
    a->fpFakeGtrkNormChi2 = bookDistribution(Form("%sFakeGtrkNormChi2", name.c_str()), "global track norm. #chi^{2}", "GlobalMuon", 24, 0., 12.);
    a->fpFakeDzRef = bookDistribution(Form("%sFakeDzRef", name.c_str()), "dzrf", "GlobalMuon", 20, -20., 20.);
    a->fpFakeDxyRef = bookDistribution(Form("%sFakeDxyRef", name.c_str()), "dxyref", "GlobalMuon", 50, -1., 1.);
    a->fpFakeGtrkProb = bookDistribution(Form("%sFakeGtrkProb", name.c_str()), "global track prob", "GlobalMuon", 21, 0., 102.);
    a->fpFakeMuonChi2 = bookDistribution(Form("%sFakeMuonChi2", name.c_str()), "muon #chi^{2}", "GlobalMuon", 21, 0., 1020.);
    a->fpFakeGlbKinkFinder = bookDistribution(Form("%sFakeGlbKinkFinder", name.c_str()), "log10(GlbKinkFinder)", "GlobalMuon", 20, -5., 15.);
    a->fpFakeStaRelChi2 = bookDistribution(Form("%sFakeStaRelChi2", name.c_str()), "StaRelChi2", "GlobalMuon", 21, 0., 102.);
    a->fpFakeTrkRelChi2 = bookDistribution(Form("%sFakeTrkRelChi2", name.c_str()), "TrkRelChi2", "GlobalMuon", 20, 0., 20.);
    a->fpFakeGlbDeltaEtaPhi = bookDistribution(Form("%sFakeGlbDeltaEtaPhi", name.c_str()), "GlbDeltaEtaPhi", "GlobalMuon", 20, -2., 4.);
    a->fpFakeTimeInOut = bookDistribution(Form("%sFakeTimeInOut", name.c_str()), "TimeInOut", "GlobalMuon", 20, -200., 200.);
    a->fpFakeTimeInOutE = bookDistribution(Form("%sFakeTimeInOutE", name.c_str()), "TimeInOutE", "GlobalMuon", 20, 0., 10.);
    a->fpFakeTimeInOutS = bookDistribution(Form("%sFakeTimeInOutS", name.c_str()), "TimeInOut/TimeInOutE", "GlobalMuon", 20, -100., 100.);
    a->fpFakeNvalidMuonHits = bookDistribution(Form("%sFakeNvalidMuonHits", name.c_str()), "NvalidMuonHits", "GlobalMuon", 21, 0., 51.);
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

    a->fpFakeTisAllEta  = bookDistribution(Form("%sFakeTisAllEta", name.c_str()), "#eta", "TIS", 12, -2.4, 2.4);
    a->fpFakeTisAllPt   = bookDistribution(Form("%sFakeTisAllPt", name.c_str()), "#it{p_{T}} [GeV]", "TIS", 5, 0., 20.);
    a->fpFakeTisFakeEta = bookDistribution(Form("%sFakeTisFakeEta", name.c_str()), "#eta", "TISFAKE", 12, -2.4, 2.4);
    a->fpFakeTisFakePt  = bookDistribution(Form("%sFakeTisFakePt", name.c_str()), "#it{p_{T}} [GeV]", "TISFAKE", 5, 0., 20.);

    a->fpFakeTisDtAllEta  = bookDistribution(Form("%sFakeTisDtAllEta", name.c_str()), "#eta", "TISDT", 12, -2.4, 2.4);
    a->fpFakeTisDtAllPt   = bookDistribution(Form("%sFakeTisDtAllPt", name.c_str()), "#it{p_{T}} [GeV]", "TISDT", 5, 0., 20.);
    a->fpFakeTisDtFakeEta = bookDistribution(Form("%sFakeTisDtFakeEta", name.c_str()), "#eta", "TISDTFAKE", 12, -2.4, 2.4);
    a->fpFakeTisDtFakePt  = bookDistribution(Form("%sFakeTisDtFakePt", name.c_str()), "#it{p_{T}} [GeV]", "TISDTFAKE", 5, 0., 20.);

    a->fpFakeTisDtDmAllEta  = bookDistribution(Form("%sFakeTisDtDmAllEta", name.c_str()), "#eta", "TISDTDM", 12, -2.4, 2.4);
    a->fpFakeTisDtDmAllPt   = bookDistribution(Form("%sFakeTisDtDmAllPt", name.c_str()), "#it{p_{T}} [GeV]", "TISDTDM", 5, 0., 20.);
    a->fpFakeTisDtDmFakeEta = bookDistribution(Form("%sFakeTisDtDmFakeEta", name.c_str()), "#eta", "TISDTDMFAKE", 12, -2.4, 2.4);
    a->fpFakeTisDtDmFakePt  = bookDistribution(Form("%sFakeTisDtDmFakePt", name.c_str()), "#it{p_{T}} [GeV]", "TISDTDMFAKE", 5, 0., 20.);


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
    masslo = 1.095;
    masshi = 1.145;
  }
  if (string::npos != fSample.find("psi")) {
    masslo = 2.90;
    masshi = 3.30;
  }

  AnalysisDistribution *p = new AnalysisDistribution(hn.c_str(), ht.c_str(), nbins, lo, hi, masslo, masshi);
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX);
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX);
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX);
  p->setAnalysisCuts(&fCncCuts, hc.c_str());
  p->setPreselCut(&fGoodCand);
  return p;
}


// ----------------------------------------------------------------------
void plotFake::sbsDistributions(string sample, string selection, std::string what) {
  cout << "plotFake::sbsDistributions(" << sample << ", " << selection << ", " << what << ")" << endl;
  string sbsControlPlotsFileName = Form("%s-sbsctrl", fSetup.c_str());
  AnalysisDistribution a(Form("%s_FakePt", sample.c_str()));
  a.fVerbose = 0;
  a.fControlPlotsFileName = sbsControlPlotsFileName;
  a.fDirectory = fDirectory;

  a.fpIF->resetLimits();

  int type(0), subtype(0);
  if (string::npos != sample.find("ks")) {
    type = 1;
    subtype = 211;
    a.fMassPeak  = 0.498;
    a.fMassSigma = 0.005;
    a.fMassLo    = 0.450;
    a.fMassHi    = 0.550;
  } else if (string::npos != sample.find("psi")) {
    type = 1;
    subtype = 13;
    a.fMassPeak  = 3.097;
    a.fMassSigma = 0.030;
    a.fMassLo    = 2.920;
    a.fMassHi    = 3.300;
  } else if (string::npos != sample.find("lambda")) {
    type = 1;
    subtype = 2212;
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
      savePad(Form("%s-prefit-%s-%s-%s.pdf", fSetup.c_str(), sample.c_str(), selection.c_str(), what.c_str()));
      a.fpIF->limitPar(1, f1->GetParameter(1) - 3.*f1->GetParError(1), f1->GetParameter(1) + 3.*f1->GetParError(1));
      a.fpIF->limitPar(2, f1->GetParameter(2) - 3.*f1->GetParError(2), f1->GetParameter(2) + 3.*f1->GetParError(2));
      cout << "done prefitting: h = " << h << " f1 = " << f1 << endl;
      cout << "done prefitting: h = " << h << " f1 = " << f1 << endl;
    }
  } else if (string::npos != sample.find("phi")) {
    type = 2;
    subtype = 321;
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

  bool restricted = (what != "");
  cout << "----------------------------------------------------------------------" << endl;
  cout << "type = " << type << " fDoList.size() = " << fDoList.size() << " restricted = " << restricted << " what = " << what << endl;
  cout << "----------------------------------------------------------------------" << endl;
  TH1D *h(0);
  string bla;
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    bla =  Form("%s_%s", sample.c_str(), fDoList[i].c_str());
    cout << "bla = " << bla << " selection = " << selection << endl;
    if (0 == type) {
      cout << "=> Looking for signal histogram " << Form("%s%s0", bla.c_str(), selection.c_str()) << endl;
      h = (TH1D*)gDirectory->Get(Form("%s%s0", bla.c_str(), selection.c_str()));
      cout << "=> cloning into  "
	   << Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), selection.c_str()) << endl;
      h = (TH1D*)h->Clone(Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), selection.c_str()));
    } else if (1 == type) {
      cout << "=> sbsDistribution histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      TH1D *hm = (TH1D*)gDirectory->Get(Form("%sMass%s", bla.c_str(), selection.c_str()));
      if (13 == subtype)   hm->SetXTitle("#it{m_{#mu#mu}} [GeV]");
      if (211 == subtype)  hm->SetXTitle("#it{m_{#pi#pi}} [GeV]");
      if (2212 == subtype) hm->SetXTitle("#it{m_{p#pi}} [GeV]");
      h = a.sbsDistribution(bla.c_str(), selection.c_str());
    } else if (2 == type) {
      cout << "=> sbsDistributionPhiKK histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      TH1D *hm = (TH1D*)gDirectory->Get(Form("%sMass%s", bla.c_str(), selection.c_str()));
      if (321 == subtype)   hm->SetXTitle("#it{m_{KK}} [GeV]");
      h = a.sbsDistributionPhiKK(bla.c_str(), selection.c_str());
    }
    cout << "  Title: " << h->GetTitle() << " with integral: " << h->GetSumOfWeights() << endl;
  }
  if (h) {
    cout << "sbsDistributions: created histogram with name: " << h->GetName() << " in directory " << h->GetDirectory()->GetName() << endl;
  } else {
    cout << "what = " << what << " not in doList: " << endl;
    for (unsigned int i = 0; i < fDoList.size(); ++i) {
      cout << "  " << fDoList[i] << endl;
    }
  }

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

  if (string::npos != sample1.find("NvalidMuonHits"))        leftLegend = true;
  if (string::npos != sample1.find("LayersWithHits"))        leftLegend = true;
  if (string::npos != sample1.find("FakeItrkValidFraction")) leftLegend = true;

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

  // label1 += Form(" (%d%s)", fYear, fSetup.c_str());
  // label2 += Form(" (%d%s)", fYear, fSetup.c_str());
  label1 += Form(" (%d)", fYear);
  label2 += Form(" (%d)", fYear);

  if (string::npos != sample1.find("OuterChi2") && (string::npos != label1.find("muons")))    doLegend = false;
  if (string::npos != sample1.find("FakeCombHits") && (string::npos != label1.find("muons"))) leftLegend = true;
  if (string::npos != sample1.find("SegmentComp") && (string::npos != label1.find("muons")))  leftLegend = true;
  if (string::npos != sample1.find("Bdt") && (string::npos != label1.find("muons")))          leftLegend = true;
  if (string::npos != sample1.find("Qprod")) {
    leftLegend = true;
    h1->SetMinimum(0.5);
    // ymax *= 5.;
    //    c0->SetLogy(1);
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
      newLegend(0.60, 0.7, 0.90, 0.87);
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

  c0->cd();
  c0->Clear();
  c0->SetCanvasSize(700, 700);

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

  int nBins;
  double int1V, int1E, int2V, int2E;
  vector<double> vval, verr;
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
    int1V = int2V = int1E = int2E = 0.;
    nBins = 1;
    double xlo(3.9), xhi(23.9);
    if (string::npos != varF.find("Eta")) {
      xlo = -2.5;
      xhi = 2.5;
    }
    for (int ibin = 1; ibin <= h1p->GetNbinsX(); ++ibin) {
      if (h1p->GetBinLowEdge(ibin) < xlo) continue;
      if (h1p->GetBinLowEdge(ibin+1) > xhi) break;
      //      ++nBins;
      cout << "looking at bin " << i << " with bin center: " << h1p->GetBinCenter(ibin) << ": ";
      if (h1p->GetBinContent(ibin) > 0) {
	int1V += h1p->GetBinContent(ibin);
        int1E += h1p->GetBinError(ibin)*h1p->GetBinError(ibin);
      }
      int2V += h2p->GetBinContent(ibin);
      int2E += h2p->GetBinError(ibin)*h2p->GetBinError(ibin);
      cout << h1p->GetBinContent(ibin)    << ", " << h2p->GetBinContent(ibin) << endl;
    }
    int1V = int1V/nBins;
    int1E = TMath::Sqrt(int1E)/nBins;
    int2V = int2V/nBins;
    int2E = TMath::Sqrt(int2E)/nBins;
    cout << "SYSTEMATIC " << dataset1 << " integral 1: " << int1V << " +/- " << int1E << " Nbins = " << nBins << endl;
    setHist(h2p, kBlue);
    h2p->Draw("histsame");
    cout << "SYSTEMATIC " << dataset2 << " integral 2: " << int2V << " +/- " << int2E << endl;
    if (string::npos != varA.find("FakeTisDtDmAllPt") && (i < 2)) {
      vval.push_back(int1V);
      verr.push_back(int1E);
      vval.push_back(int2V);
      verr.push_back(int2E);
    }
    newLegend(0.21, 0.7, 0.41, 0.87);

    legg->SetHeader(header.c_str());
    legg->SetTextSize(0.05);
    // legg->AddEntry(h1p, Form("%s (%5.4f#pm%5.4f)", label1.c_str(), int1V, int1E), loption1);
    // legg->AddEntry(h2p, Form("%s (%5.4f#pm%5.4f)", label2.c_str(), int2V, int2E), loption2);
    legg->AddEntry(h1p, Form("%s", label1.c_str()), loption1);
    legg->AddEntry(h2p, Form("%s", label2.c_str()), loption2);
    legg->Draw();

    tl->DrawLatexNDC(0.60, 0.78, Form("%6.5f#pm%6.5f", int1V, int1E));
    tl->DrawLatexNDC(0.60, 0.72, Form("%6.5f#pm%6.5f", int2V, int2E));

    double dV = int1V - int2V;
    double dE = TMath::Sqrt(int1E*int1E + int2E*int2E);
    tl->DrawLatexNDC(0.49, 0.66, Form("#Delta = %+6.5f#pm%6.5f", dV, dE));


    savePad(Form("%s-fakerate_%s_ad%s_%s_ad%s_%s.pdf",
		 fSetup.c_str(),
		 varA.c_str(), fChannelList[i].c_str(), dataset1.c_str(), fChannelList[i].c_str(), dataset2.c_str()));
  }

  if (varA == "FakeTisDtDmAllPt") {
    double sysVal(0.), sysErr(0.);
    average(sysVal, sysErr, vval, verr);
    fTEX << formatTex(sysVal, Form("%s:muonidSystematics_%s:val", fSetup.c_str(), dataset1.c_str()), 5) << endl;
    fTEX << formatTex(sysErr, Form("%s:muonidSystematics_%s:err", fSetup.c_str(), dataset1.c_str()), 5) << endl;
    fTEX << formatTex(sysErr/sysVal, Form("%s:muonidSystematics_%s:sys", fSetup.c_str(), dataset1.c_str()), 4) << endl;
  }

  fHistFile->Close();
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
  h1->SetNdivisions(504, "X");
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
    if (fCandPvIp  > 0.006) fGoodCand = false;
    if (fCandFLS3d < 25)  fGoodCand = false;
    if (fCandFLxy  > 3.9)  fGoodCand = false;
  }
  if (fMode == FAKEKS) {
    if (fCandPvIp > 0.006) fGoodCand = false;
    if (fCandFLS3d < 20)   fGoodCand = false;
  }

  double mass = fCandM;

  if (fIsMC) {
    mass = SIGBOXMIN + 0.5 * (SIGBOXMAX - SIGBOXMIN);
  }

  string mapname("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");


  int nfakes(0);
  for (int i = 0; i < fFakeNtrk; ++i) {
    if (fFakeGm[i]>0) ++nfakes;
  }

  bool singleFake = (nfakes == 1);
  if (fMode == FAKEPSI) singleFake = true;

  for (int i = 0; i < fFakeNtrk; ++i) {
    if (TMath::Abs(fFakeEta[i]) < 0.7) {
      fChan = 0;
    } else if (0.7 < TMath::Abs(fFakeEta[i])  && TMath::Abs(fFakeEta[i]) < 1.4) {
      fChan = 1;
    } else if (1.4 < TMath::Abs(fFakeEta[i])  && TMath::Abs(fFakeEta[i]) < 2.0) {
      fChan = 2;
    } else if (2.0 < TMath::Abs(fFakeEta[i])  && TMath::Abs(fFakeEta[i]) < 2.4) {
      fChan = 3;
    } else {
      continue;
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

    if (fYear > 2012 || fIsMC) {
      fFakeSsVeto = false;
    }

    fGoodPt      = (fFakePt[i] > 4.);
    fGoodDtrig   = (fFakeDtrig[i] > 0.01) && fTIS && !fFakeSsVeto;
    fGoodDmuon   = (fFakeDmuon[i] > 0.5);
    fGoodDNmuon  = (fFakeDNmuon[i] > 0.5);

    if (fIsMC) {
      fTIS        = true;
      fGoodDtrig  = true; // does not work on MC (PD not well defined there)
    } else {
      fGoodDmuon =  fGoodDNmuon;
    }

    if (fMode == FAKEPSI) {
      fGlobalMuon  = (fFakeGm[i] > 0) && singleFake;
    } else {
      // -- there are no light resonances -> muon triggers in the Charmonium PD!
      fGlobalMuon  = (fFakeGm[i] > 0) && singleFake && fGoodDtrig;
    }

    fGood     = fGoodCand && fGoodPt;
    fGoodFake = fGood && fGlobalMuon;
    fGoodFake = fGood && fGlobalMuon && (fFakeBdt[i] > fCuts[0]->muonbdt);

    fGoodTIS         = fTIS       && fGood;
    fGoodTISFake     = fGoodTIS   && fGoodFake;

    fGoodTISDT       = fGoodTIS   && fGoodDtrig;
    fGoodTISDTFake   = fGoodTISDT && fGoodFake;

    fGoodTISDTDM     = fGoodTISDT   && fGoodDmuon;
    fGoodTISDTDMFake = fGoodTISDTDM && fGoodFake;

    if ((fMode == FAKELAMBDA) && fFakeId[i] != 2212) {
      //      cout << "discarding i = " << i << " fFakeId[i] = " << fFakeId[i] << endl;
      continue;
    }

    fCncCuts.update();


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

    fAdMap[mapname]->fpFakeQprod->fill(fFakeQprod[i], mass);
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
void plotFake::mkPidTables(string prefix) {

  if ("" == prefix) prefix = "bmm4-25";

  fHists.clear();
  vector<int> vIds;
  vIds.push_back(13);
  vIds.push_back(211);
  vIds.push_back(321);
  vIds.push_back(2212);
  vector<string> q;
  q.push_back("Pos");
  q.push_back("Neg");

  cout << "fHistFileName: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;

  TH2D *h2(0);
  string name, title;
  name = Form("allNegId13");
  h2 = (TH2D*)fHistFile->Get(name.c_str());

  fMuBdtCutB = fCuts[0]->muonbdt;
  fMuBdtCutE = fCuts[1]->muonbdt;

  // -- create histograms first if not already present
  if (0 == h2) {
    for (unsigned int i = 0; i < vIds.size(); ++i) {
      for (unsigned int iq = 0; iq < q.size(); ++iq) {
	name = Form("all%sId%d", q[iq].c_str(), vIds[i]);
	title = Form("bdt > %4.2f,%4.2f", fMuBdtCutB, fMuBdtCutE);
	h2 = new TH2D(name.c_str(), title.c_str(), 21, 0., 2.1, 100, 0., 50.);
	setTitles(h2, "|#eta|", "#it{p_{T}} [GeV]");
	fHists[name] = h2;
	name = Form("pass%sId%d", q[iq].c_str(), vIds[i]);
	h2 = new TH2D(name.c_str(), title.c_str(), 21, 0., 2.1, 100, 0., 50.);
	setTitles(h2, "|#eta|", "#it{p_{T}} [GeV]");
	fHists[name] = h2;
      }
    }
    TTree *t = getTree("fakeMc", "candAnaFakeMC", "fakeTree");
    fCds = fDS["fakeMc"];
    if (0 == t) {
      return;
    }
    setupTree(t);
    fIsMC = true;

    int nevents(0), nstart(0);
    cout << "Running loopOverTree with fMuBdtCutB = " << fMuBdtCutB << " and fMuBdtCutE = " << fMuBdtCutE << endl;
    loopOverTree(t, 2, nevents, nstart);

    for (map<string, TH1*>::iterator it = fHists.begin(); it != fHists.end(); ++it) {
      it->second->Write();
    }
  }

  // -- analyze histograms and write PidTables
  PidTable a("fakeTemplate.dat");
  PidTable a5("fakeTemplatePt5.dat");
  PidTable a10("fakeTemplatePt10.dat");
  PidTable A("effTemplate.dat");
  PidTable *pa = &a;

  PidTable b;
  PidTable c;
  string aname, pname;

  for (unsigned int i = 0; i < vIds.size(); ++i) {
    if (13 == vIds[i]) continue; // use special template for muons
    for (unsigned int iq = 0; iq < q.size(); ++iq) {
      if (2011 == fYear) {
	if (2212 == vIds[i]) {
	  pa = &a5;
	} else {
	  pa = &a10;
	}
      } else if (2012 == fYear) {
	if (2212 == vIds[i]) {
	  pa = &a5;
	} else {
	  pa = &a10;
	}
      } else if (2016 == fYear) {
	if (2212 == vIds[i]) {
	  pa = &a10;
	} else {
	  pa = &a10;
	}
      } else {
	pa = &a;
      }

      aname = Form("all%sId%d", q[iq].c_str(), vIds[i]);
      pname = Form("pass%sId%d", q[iq].c_str(), vIds[i]);
      pa->flush();
      b.flush();
      b.readFromHist(fHistFile, pname.c_str(), aname.c_str());
      pa->fillEff(b);
      h2 = (TH2D*)fHistFile->Get(pname.c_str());
      pa->setComment(h2->GetTitle());
      name = Form("weights/pidtables/%d-%d%s-%s.dat", fYear, vIds[i], q[iq].c_str(), prefix.c_str());
      cout << name << endl;
      pa->dumpToFile(name.c_str());
    }
  }


  for (unsigned int iq = 0; iq < q.size(); ++iq) {
    aname = Form("all%sId%d", q[iq].c_str(), 13);
    pname = Form("pass%sId%d", q[iq].c_str(), 13);
    A.flush();
    b.flush();
    b.readFromHist(fHistFile, pname.c_str(), aname.c_str());
    A.fillEff(b);
    h2 = (TH2D*)fHistFile->Get(pname.c_str());
    A.setComment(h2->GetTitle());
    name = Form("weights/pidtables/%d-%d%s-%s.dat", fYear, 13, q[iq].c_str(), prefix.c_str());
    cout << name << endl;
    A.dumpToFile(name.c_str());
  }


  fHistFile->Close();
}

// ----------------------------------------------------------------------
void plotFake::plotPidTables(string prefix) {

  // -- Dump exactly what is used in the analysis
  PidTable *a;
  gStyle->SetOptTitle(0);
  tl->SetTextSize(0.07);

  c0->cd();
  c0->Clear();
  c0->SetCanvasSize(700, 700);

  // -- hadrons
  if ("" == prefix) {
    gPad->SetLogy(0);
    gStyle->SetPaintTextFormat("5.4f");
    double xbins[] = {0., 0.7, 1.4, 2.1};
    double ybins[] = {0., 4., 5., 7., 10., 20., 30.};
    TH2D *h2 = new TH2D("h2", "", 3, xbins, 6, ybins);
    setTitles(h2, "|#eta|", "#it{p_{T}} [GeV]", 0.05, 1.1, 1.3);
    h2->SetMinimum(0.0);
    h2->SetMaximum(0.002);
    h2->SetMarkerSize(1.3);
    h2->SetMarkerColor(kBlack);

    gStyle->SetOptStat(0);
    shrinkPad(0.15, 0.15, 0.1);
    tl->SetTextSize(0.04);

    a = fptFakePosKaons;  h2->Reset(); h2->SetTitle(Form("pos. kaons (%d, %s)", fYear, a->getComment().Data()));
    cout << "Hallo" << endl;
    a->print(cout);
    a->eff2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-effFakePosKaons.pdf", fDirectory.c_str(), fSetup.c_str()));
    a->err2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-errFakePosKaons.pdf", fDirectory.c_str(), fSetup.c_str()));

    a = fptFakeNegKaons;  h2->Reset(); h2->SetTitle(Form("neg. kaons (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-effFakeNegKaons.pdf", fDirectory.c_str(), fSetup.c_str()));
    a->err2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-errFakeNegKaons.pdf", fDirectory.c_str(), fSetup.c_str()));

    a = fptFakePosPions;  h2->Reset(); h2->SetTitle(Form("pos. pions (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-effFakePosPions.pdf", fDirectory.c_str(), fSetup.c_str()));
    a->err2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-errFakePosPions.pdf", fDirectory.c_str(), fSetup.c_str()));

    a = fptFakeNegPions;  h2->Reset(); h2->SetTitle(Form("neg. pions (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-effFakeNegPions.pdf", fDirectory.c_str(), fSetup.c_str()));
    a->err2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-errFakeNegPions.pdf", fDirectory.c_str(), fSetup.c_str()));

    a = fptFakePosProtons;  h2->Reset(); h2->SetTitle(Form("pos. protons (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-effFakePosProtons.pdf", fDirectory.c_str(), fSetup.c_str()));
    a->err2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-errFakePosProtons.pdf", fDirectory.c_str(), fSetup.c_str()));

    a = fptFakeNegProtons;  h2->Reset(); h2->SetTitle(Form("neg. protons (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-effFakeNegProtons.pdf", fDirectory.c_str(), fSetup.c_str()));
    a->err2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-errFakeNegProtons.pdf", fDirectory.c_str(), fSetup.c_str()));
  }


  // -- muons
  if ("" == prefix) {
    gStyle->SetPaintTextFormat("3.2f");
    double xbins[] = {0.0, 0.3, 0.60, 0.90, 1.2, 1.5, 1.8, 2.1};
    double ybins[] = {3.9, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 15., 20., 30.};
    TH2D *h2 = new TH2D("h2", "", 7, xbins, 10, ybins);
    setTitles(h2, "|#eta|", "#it{p_{T}} [GeV]", 0.05, 1.1, 1.3);
    h2->SetMinimum(0.0);
    h2->SetMaximum(1.0);
    h2->SetMarkerSize(1.3);
    h2->GetYaxis()->SetMoreLogLabels();
    //    h2->SetMarkerColor(kWhite);

    gStyle->SetOptStat(0);

    shrinkPad(0.15, 0.15, 0.1);
    gPad->SetLogy(1);

    // -- muon id
    a = fptPosMuons;  h2->Reset(); h2->SetTitle(Form("pos. muons (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2);
    h2->Draw("coltext");
    tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-effPosMuons.pdf", fDirectory.c_str(), fSetup.c_str()));

    a = fptNegMuons;  h2->Reset(); h2->SetTitle(Form("neg. muons (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2);
    h2->Draw("coltext");
    tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-effNegMuons.pdf", fDirectory.c_str(), fSetup.c_str()));
  }




  // -- and now the custom display of non-standard PidTables
  vector<int> vIds;
  vIds.push_back(13);
  vIds.push_back(211);
  vIds.push_back(321);
  vIds.push_back(2212);
  vector<string> sIds;
  sIds.push_back("muons");
  sIds.push_back("pions");
  sIds.push_back("kaons");
  sIds.push_back("protons");
  vector<string> q;
  q.push_back("Pos");
  q.push_back("Neg");
  // -- hadrons
  if ("" != prefix) {
    gStyle->SetPaintTextFormat("5.4f");
    double xbins[] = {0., 0.7, 1.4, 2.1};
    double ybins[] = {0., 4., 5., 7., 10., 20., 30.};
    TH2D *h2 = new TH2D("h2", "", 3, xbins, 6, ybins);
    setTitles(h2, "|#eta|", "#it{p_{T}} [GeV]");
    h2->SetMinimum(0.0);
    h2->SetMaximum(0.002);
    h2->SetMarkerSize(1.3);
    h2->SetMarkerColor(kBlack);

    gStyle->SetOptStat(0);
    shrinkPad(0.15, 0.15, 0.25);
    tl->SetTextSize(0.04);
    string name("");
    for (unsigned int i = 0; i < vIds.size(); ++i) {
      if (13 == vIds[i]) continue; // use special template for muons
      for (unsigned int iq = 0; iq < q.size(); ++iq) {
	name = Form("weights/pidtables/%d-%d%s-%s.dat", fYear, vIds[i], q[iq].c_str(), prefix.c_str());
	cout << "pidtable with name = " << name << endl;
	PidTable *a = new PidTable(Form(name.c_str()));

	h2->Reset(); h2->SetTitle(Form("%s %s (%d, %s)", q[iq].c_str(), sIds[i].c_str(), fYear, a->getComment().Data()));
	a->eff2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
	c0->SaveAs(Form("%s/%s-eff-%d%s-%s.pdf", fDirectory.c_str(), fSetup.c_str(), vIds[i], q[iq].c_str(), prefix.c_str()));
	a->err2d(h2); h2->Draw("coltext"); tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
	c0->SaveAs(Form("%s/%s-err-%d%s-%s.pdf", fDirectory.c_str(), fSetup.c_str(), vIds[i], q[iq].c_str(), prefix.c_str()));
      }
    }
  }


  // -- muons
  if ("" != prefix) {
    gStyle->SetPaintTextFormat("3.2f");
    double xbins[] = {0.0, 0.3, 0.60, 0.90, 1.2, 1.5, 1.8, 2.1};
    double ybins[] = {4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10., 15., 20., 30.};
    TH2D *h2 = new TH2D("h2", "", 7, xbins, 10, ybins);
    setTitles(h2, "|#eta|", "#it{p_{T}} [GeV]");
    h2->SetMinimum(0.0);
    h2->SetMaximum(1.0);
    h2->SetMarkerSize(1.3);
    //    h2->SetMarkerColor(kWhite);

    gStyle->SetOptStat(0);

    shrinkPad(0.15, 0.15, 0.25);
    gPad->SetLogy(1);

    // -- muon id
    string name = Form("weights/pidtables/%d-%d%s-%s.dat", fYear, 13, "Pos", prefix.c_str());
    PidTable *a = new PidTable(Form(name.c_str()));
    h2->Reset(); h2->SetTitle(Form("Pos muons (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2);  h2->Draw("coltext");  tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-eff-%d%s-%s.pdf", fDirectory.c_str(), fSetup.c_str(), 13, "Pos", prefix.c_str()));

    name = Form("weights/pidtables/%d-%d%s-%s.dat", fYear, 13, "Neg", prefix.c_str());
    a = new PidTable(Form(name.c_str()));
    h2->Reset(); h2->SetTitle(Form("Neg muons (%d, %s)", fYear, a->getComment().Data()));
    a->eff2d(h2);  h2->Draw("coltext");  tl->DrawLatexNDC(0.20, 0.92, h2->GetTitle());
    c0->SaveAs(Form("%s/%s-eff-%d%s-%s.pdf", fDirectory.c_str(), fSetup.c_str(), 13, "Neg", prefix.c_str()));
  }





}


// ----------------------------------------------------------------------
void plotFake::loopFunction2() {
  string aname("all"), pname("pass");
  string name(""), h0name(""), h1name("");

  for (int i = 0; i < fFakeNtrk; ++i) {
    if (fFakePt[i] < 4.0) continue;
    if (TMath::Abs(fFakeEta[i]) > 2.1) continue;
    if (0 == fFakeHP[i]) continue;
    if (13 == fFakeId[i]) {
      if (fFakeQ[i] < 0) {
	name = "NegId13";
      } else {
	name = "PosId13";
      }
    } else if (211   == fFakeId[i]) {
      if (fFakeQ[i] < 0) {
	name = "NegId211";
      } else {
	name = "PosId211";
      }
    } else if (321  == fFakeId[i]) {
      if (fFakeQ[i] < 0) {
	name = "NegId321";
      } else {
	name = "PosId321";
      }
    } else if (2212  == fFakeId[i]) {
      if (fFakeQ[i] < 0) {
	name = "NegId2212";
      } else {
	name = "PosId2212";
      }
    } else {
      //      cout << "unknown particle: " << fFakeId[i] << endl;
      continue;
    }
    h0name = aname + name;
    //    cout << h0name << endl;
    fHists[h0name]->Fill(TMath::Abs(fFakeEta[i]), fFakePt[i]);
    h1name = pname + name;
    //    cout << h1name << endl;
    if (TMath::Abs(fFakeEta[i]) < 1.4) {
      if (fFakeBdt[i] > fMuBdtCutB) fHists[h1name]->Fill(TMath::Abs(fFakeEta[i]), fFakePt[i]);
    } else {
      if (fFakeBdt[i] > fMuBdtCutE) fHists[h1name]->Fill(TMath::Abs(fFakeEta[i]), fFakePt[i]);
    }
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
  cout << "==> plotFake::loopOverTree> loop over dataset " << fCds->fName << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotFake::*pF)(void);
  if (ifunc == 1) pF = &plotFake::loopFunction1;
  if (ifunc == 2) pF = &plotFake::loopFunction2;

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
  t->SetBranchAddress("vetoSameSign", &fFakeSsVeto);
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
  t->SetBranchAddress("dNmuon",  fFakeDNmuon);

  t->SetBranchAddress("hp",      fFakeHP);
  t->SetBranchAddress("bdt",     fFakeBdt);

  t->SetBranchAddress("tip", fFakeTip);
  t->SetBranchAddress("lip", fFakeLip);

  t->SetBranchAddress("qprod", fFakeQprod);
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
