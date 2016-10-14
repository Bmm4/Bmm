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
  init();

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

  fChannelList.clear();
  for (unsigned int i = 0; i < 1; ++i) {
    fChannelList.push_back(Form("%d", i));
  }

  fDoList.clear();
  fDoList.push_back("FakePt");
  fDoList.push_back("FakeEta");

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
  fDoList.push_back("FakeDz");
  fDoList.push_back("FakeLip");
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

  fAnaCuts.clear();
  fAnaCuts.addCut("GoodCand", "good cand", fGoodCand);
  fAnaCuts.addCut("GoodPt", "good pt", fGoodPt);
  fAnaCuts.addCut("GlobalMuon", "global muon ID", fGlobalMuon);

  fAnaCuts.dumpAll();
}


// ----------------------------------------------------------------------
plotFake::~plotFake() {

}



// ----------------------------------------------------------------------
void plotFake::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotFake::makeAll(string what) {


  if (what == "all") {
    init();
    makeSample("fakeData", "ks");
    makeSample("fakeMc", "ks");

    makeSample("fakeData", "psi");
    makeSample("fakeMc", "psi");

    makeOverlay("fakeData_ks", "fakeMc_ks");
  }

  if (string::npos != what.find("plot")) {
    fTEX.close();
    system(Form("/bin/rm -f %s", fTexFileName.c_str()));
    fTEX.open(fTexFileName.c_str(), ios::app);
    system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
    system(Form("/bin/rm -f %s/sbsctrl_ad*_fake*.pdf", fDirectory.c_str()));
    system(Form("/bin/rm -f %s/ad*_fake*.pdf", fDirectory.c_str()));

    makeOverlay("fakeData_ks", "fakeMc_ks", "Cu");
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
    BGLBOXMIN = 0.450;
    BGLBOXMAX = 0.465;
    SIGBOXMIN = 0.480;
    SIGBOXMAX = 0.515;
    BGHBOXMIN = 0.530;
    BGHBOXMAX = 0.550;
  }

  if (string::npos != sample.find("psi")) {
    fMode = FAKEPSI;
    dir = "candAnaFake443";
    BGLBOXMIN = 2.800;
    BGLBOXMAX = 2.900;
    SIGBOXMIN = 3.050;
    SIGBOXMAX = 3.150;
    BGHBOXMIN = 3.200;
    BGHBOXMAX = 3.300;
  }

  if (fIsMC) dir = "candAnaFakeMC";


  // -- must be after the mass box definitions!
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
      overlay(Form("sbs_ad%s_%s_%s%s", fChannelList[i].c_str(), what1.c_str(), fDoList[id].c_str(), selection.c_str()),
	      Form("sbs_ad%s_%s_%s%s", fChannelList[i].c_str(), what2.c_str(), fDoList[id].c_str(), selection.c_str())
	      );
      savePad(Form("fakeoverlay_%s_ad%s_%s_ad%s_%s.pdf", fDoList[id].c_str(), fChannelList[i].c_str(), what1.c_str(), fChannelList[i].c_str(), what1.c_str()));
    }
  }
  fHistFile->Close();

}

// ----------------------------------------------------------------------
void plotFake::fillDistributions() {



}


// ----------------------------------------------------------------------
void plotFake::bookDistributions() {

  cout << "SIG: " << SIGBOXMIN << " .. " << SIGBOXMAX << endl;
  cout << "BGL: " << BGLBOXMIN << " .. " << BGLBOXMAX << endl;
  cout << "BGH: " << BGHBOXMIN << " .. " << BGHBOXMAX << endl;

  adset *a(0);
  for (unsigned int i = 0; i < fChannelList.size(); ++i) {
    string mapname = Form("ad%s_%s", fChannelList[i].c_str(), fSample.c_str());
    string name = Form("%s_", mapname.c_str());

    a = new adset();
    a->fpFakeEta  = bookDistribution(Form("%sFakeEta", name.c_str()), "#eta", "GlobalMuon", 40, -2.4, 2.4);
    a->fpFakePt   = bookDistribution(Form("%sFakePt", name.c_str()), "p_{T} [GeV]", "GlobalMuon", 40, 0., 20.);
    a->fpFakeInnerChi2 = bookDistribution(Form("%sFakeInnerChi2", name.c_str()), "inner track #chi^{2}", "GlobalMuon", 51, 0., 102.);
    a->fpFakeOuterChi2 = bookDistribution(Form("%sFakeOuterChi2", name.c_str()), "outer track #chi^{2}", "GlobalMuon", 51, 0., 102.);

    a->fpFakeChi2LocalPosition = bookDistribution(Form("%sFakeChi2LocalPosition", name.c_str()), "local position #chi^{2}", "GlobalMuon", 51, 0., 102.);
    a->fpFakeChi2LocalMomentum = bookDistribution(Form("%sFakeChi2LocalMomentum", name.c_str()), "local momentum #chi^{2}", "GlobalMuon", 51, 0., 102.);
    a->fpFakeStaTrkMult = bookDistribution(Form("%sFakeStaTrkMult", name.c_str()), "STA trk multipicity", "GlobalMuon", 12, -2., 10.);
    a->fpFakeTmTrkMult = bookDistribution(Form("%sFakeTmTrkMult", name.c_str()), "TM trk multiplicity", "GlobalMuon", 20, 0., 20.);

    a->fpFakeDeltaR = bookDistribution(Form("%sFakeDeltaR", name.c_str()), "deltaR", "GlobalMuon", 40, 0., 1.0);
    a->fpFakeItrkValidFraction = bookDistribution(Form("%sFakeItrkValidFraction", name.c_str()), "inner track valid fraction", "GlobalMuon", 50, 0., 1.02);
    a->fpFakeSegmentComp = bookDistribution(Form("%sFakeSegmentComp", name.c_str()), "segment compatibility", "GlobalMuon", 50, 0., 1.02);
    a->fpFakeGtrkNormChi2 = bookDistribution(Form("%sFakeGtrkNormChi2", name.c_str()), "global track norm. #chi^{2}", "GlobalMuon", 40, 0., 12.);
    a->fpFakeDz = bookDistribution(Form("%sFakeDz", name.c_str()), "dz", "GlobalMuon", 40, -20., 20.);
    a->fpFakeLip = bookDistribution(Form("%sFakeLip", name.c_str()), "lip", "GlobalMuon", 100, -1., 1.);
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
    masslo = 1.10;
    masshi = 1.50;
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
  // cout << "fHistFileName: " << fHistFileName;
  // fHistFile = TFile::Open(fHistFileName.c_str());
  // cout << " opened " << endl;

  string sbsControlPlotsFileName = Form("sbsctrl");

  AnalysisDistribution a(Form("%s_FakePt", sample.c_str()));
  a.fVerbose = 1;
  a.fControlPlotsFileName = sbsControlPlotsFileName;
  a.fDirectory = fDirectory;

  int type(0);
  if (string::npos != sample.find("ks"))             type = 1; // pol1
  if (string::npos != sample.find("psi"))            type = 1; // pol1

  // -- override the above choice in case of MC
  if (string::npos != sample.find("Mc"))             type = 0; // signal window

  if (string::npos != sample.find("ks")) {
    a.fMassPeak  = 0.498;
    a.fMassSigma = 0.005;
    a.fMassLo    = 0.450;
    a.fMassLo    = 0.550;
  } else if (string::npos != sample.find("psi")) {
    a.fMassPeak  = 3.097;
    a.fMassSigma = 0.030;
    a.fMassLo    = 2.800;
    a.fMassLo    = 3.300;
  } else {
    a.fMassPeak = 5.27;
    if (string::npos != fChannel.find("0")) {
      a.fMassSigma = 0.02;
    } else {
      a.fMassSigma = 0.03;
    }
  }

  cout << "gDIRECTORY: "; gDirectory->pwd();
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
      cout << "=> sbsDistributionPol1ErrGauss histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      h = a.sbsDistribution(bla.c_str(), selection.c_str());
    }

    cout << "  Title: " << h->GetTitle() << " with integral: " << h->GetSumOfWeights() << endl;

    c0->cd();
    h->Draw();
    savePad(Form("%s.pdf", bla.c_str()));
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
  h1->SetStats(0);
  string hname = h1->GetName();
  if (string::npos != hname.find("Data")) {
    h1->Draw("e");
  } else {
    h1->Draw();
  }

  string label1("bla"), label2("bla");
  if (string::npos != sample2.find("psi")) {
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
      newLegend(0.21, 0.7, 0.41, 0.87);
    } else {
      newLegend(0.50, 0.7, 0.75, 0.87);
    }

    legg->SetHeader(header.c_str());
    legg->SetTextSize(0.05);
    legg->AddEntry(h1, label1.c_str(), loption1);
    legg->AddEntry(h2, label2.c_str(), loption2);
    legg->Draw();
  }

  stamp(0.18, fStampCms, fStampString, 0.4, fStampLumi);

}



// ----------------------------------------------------------------------
void plotFake::fakeRate(string var, string dataset, string particle) {

  tl->SetNDC(kTRUE);

  c0->Clear();
  //  c0->Divide(4, 4);
  c0->Divide(2, 3);
  int ipad(0);

  vector<string> mode;
  mode.push_back("nmu");
  //  mode.push_back("muo");

  map<string, TH1D*> hmode;

  vector<string> histos;
  if (particle == "pion") {
    histos.push_back(Form("candAnaFake310/%s1", var.c_str()));
    histos.push_back(Form("candAnaFake310/%s2", var.c_str()));
  } else if (particle == "kaon") {
    histos.push_back(Form("candAnaFake333/%s1", var.c_str()));
    //    histos.push_back(Form("candAnaFake333/%s2", var.c_str()));
  }  else if (particle == "proton") {
    histos.push_back(Form("candAnaFake3122/%s1", var.c_str()));
  } else {
    cout << "particle " << particle << " not known, returning" << endl;
  }

  // -- get "default" to properly initialize results histograms
  string hname = histos[0] + mode[0];
  TH2D *h2 = fDS[dataset]->getHist2(hname, false);

  for (unsigned int imode = 0; imode < mode.size(); ++imode) {
    hname = particle + "_" + mode[imode];
    TH1D *hr = new TH1D(hname.c_str(), hname.c_str(), h2->GetNbinsX(), h2->GetXaxis()->GetXbins()->GetArray());
    hr->Sumw2();
    hmode[hname] = hr;
    for (unsigned int ihist = 0; ihist < histos.size(); ++ihist) {
      hname = histos[ihist] + mode[imode];
      cout << "====> getting " << hname << " from dataset " << dataset << endl;
      TH2D *h2 = fDS[dataset]->getHist2(hname, false);
      int nbins(h2->GetNbinsX());
      cout << "x bins: " << nbins  << endl;
      TH1D *h1(0);
      //      for (int i = 1; i <= nbins; ++i) {
      for (int i = 2; i <= 3; ++i) {
	c0->cd(++ipad);
	h1 = h2->ProjectionY(Form("hist%dpt%d", ihist, i), i, i);
	h1->SetTitle(Form("hist%dpt%d %s", ihist, i, h1->GetTitle()));
	if (string::npos != hname.find("Fake310")) {
	  fitKs(h1);
	} else if (string::npos != hname.find("Fake333")) {
	  fitPhi(h1);
	} else if (string::npos != hname.find("Fake3122")) {
	  fitLambda(h1);
	}
	tl->DrawLatex(0.2, 0.7, Form("%4.2f#pm%4.2f", fYield, fYieldE));
	hr->SetBinContent(i, hr->GetBinContent(i) + fYield);
	hr->SetBinError(i, TMath::Sqrt(hr->GetBinError(i)*hr->GetBinError(i) + fYieldE*fYieldE));
      }
      c0->cd(++ipad);
      h2->Draw("colz");
    }
  }
  // c0->cd(++ipad);
  // hmode[Form("%s_nmu", particle.c_str())]->Draw("e1");
  // for (int i = 1; i <= hmode["pion_nmu"]->GetNbinsX(); ++i) {
  //   cout << hmode["pion_nmu"]->GetBinLowEdge(i) << ": " << hmode["pion_nmu"]->GetBinContent(i)
  // 	 << " +/- " << hmode["pion_nmu"]->GetBinError(i)
  // 	 << endl;
  // }

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
  double NSIG(5.0);
  double xmin = xpeak - NSIG*sigma;
  double xmax = xpeak + NSIG*sigma;
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
  double NSIG(3.0);
  double xmin = xpeak - NSIG*sigma;
  double xmax = xpeak + NSIG*sigma;
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
  double NSIG(5.0);
  double xmin = xpeak - NSIG*sigma;
  double xmax = xpeak + NSIG*sigma;
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
  fChannel = "0";

  string mapname = Form("ad%s_%s", fChannel.c_str(), fSample.c_str());
  double mass = fCandM;

  if (fIsMC) {
    mass = SIGBOXMIN + 0.5 * (SIGBOXMAX - SIGBOXMIN);
  }

  //  cout << "cand m = " << fCandM << endl;

  for (int i = 0; i < fFakeNtrk; ++i) {
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

    fGlobalMuon = (fFakeGm[i]>0);
    fGoodPt = fFakePt[i] > 4.;
    fAnaCuts.update();
    fAdMap[mapname]->fpFakePt->fill(fFakePt[i], mass);
    fAdMap[mapname]->fpFakeEta->fill(fFakeEta[i], mass);


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
    fAdMap[mapname]->fpFakeDz->fill(fFakeDz[i], mass);
    fAdMap[mapname]->fpFakeLip->fill(fFakeLip[i], mass);
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
  //  t->SetBranchAddress("maxdoca", &fCandDoca);

  t->SetBranchAddress("ntrk",    &fFakeNtrk);
  t->SetBranchAddress("id",      fFakeId);
  t->SetBranchAddress("q",       fFakeQ);
  t->SetBranchAddress("gm",      fFakeGm);
  t->SetBranchAddress("pt",      fFakePt);
  t->SetBranchAddress("eta",     fFakeEta);
  t->SetBranchAddress("phi",     fFakePhi);

  t->SetBranchAddress("bdt",     fFakeBdt);

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
  t->SetBranchAddress("dz", fFakeDz);
  t->SetBranchAddress("lip", fFakeLip);
  t->SetBranchAddress("gtrkprob", fFakeGtrkProb);
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



}
