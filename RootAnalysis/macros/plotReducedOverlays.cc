#include "plotReducedOverlays.hh"

#include "common/AnalysisDistribution.hh"
#include "common/HistCutEfficiency.hh"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVirtualFitter.h"

#include "common/util.hh"

using namespace std;

ClassImp(plotReducedOverlays)

// ----------------------------------------------------------------------
plotReducedOverlays::plotReducedOverlays(string dir, string files, string cuts, string setup) : plotClass(dir, files, cuts, setup) {

  changeSetup(dir, "plotReducedOverlays", setup);

  TVirtualFitter::SetMaxIterations(50000);

  cout << "==> plotReducedOverlays files: " << files << " dir: " << dir << " cuts: " << cuts << " setup: " << setup << endl;

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();
  fNchan = 2;
  fDoCNC = true;
  //  fDoCNC = false;

  printCuts(cout);

  plotClass::loadFiles(files);
  plotReducedOverlays::loadFiles(files);

  //  fIncludeOverflowInLastBin = true;
  fIncludeOverflowInLastBin = false;

  fIsMC = false;
  fIsSignal = false;

  fSel0 = false;
  fSel1 = false;
  fSel2 = false;

  fDoList.clear();
  fDoList.push_back("bdtsel2");
  if (1) {
    fDoList.push_back("fls3d");
    fDoList.push_back("muon1pt");
    fDoList.push_back("muon2pt");

    fDoList.push_back("muonseta");
    fDoList.push_back("pt");
    fDoList.push_back("p");
    fDoList.push_back("pz");
    fDoList.push_back("eta");
    fDoList.push_back("alpha");

    fDoList.push_back("iso");
    fDoList.push_back("closetrk");
    fDoList.push_back("docatrk");

    fDoList.push_back("chi2dof");
    fDoList.push_back("pchi2dof");
    fDoList.push_back("fl3d");
    fDoList.push_back("fl3de");
    fDoList.push_back("flsxy");
    fDoList.push_back("flxy");
    fDoList.push_back("flxye");

    fDoList.push_back("maxdoca");
    fDoList.push_back("ip");
    fDoList.push_back("ips");
    fDoList.push_back("pvn");
    fDoList.push_back("pvntrk");
    fDoList.push_back("pv2ntrk");
    fDoList.push_back("pvz");
    fDoList.push_back("dzmin");
    fDoList.push_back("dz12");
    fDoList.push_back("pvavew8");

    fDoList.push_back("lip");
    fDoList.push_back("lips");

    fDoList.push_back("lip2");
    fDoList.push_back("lips2");

    fDoList.push_back("m1iso");
    fDoList.push_back("m2iso");
    fDoList.push_back("othervtx");
    fDoList.push_back("pvdchi2");
    // fDoList.push_back("closetrks1");
    // fDoList.push_back("closetrks2");
    // fDoList.push_back("closetrks3");

    fDoList.push_back("tau");
    fDoList.push_back("bdt");
    fDoList.push_back("bdtsel0");
    fDoList.push_back("bdtsel1");
  }

  fChannelList.clear();
  // -- no restrictions, just normal analysis channels
  for (unsigned int i = 0; i < fNchan; ++i) {
    fChannelList.push_back(Form("%d", i));
  }

  // -- small N(PV)
  fChannelList.push_back("0lopu");
  //  fChannelList.push_back("1lopu");

  // -- high N(PV)
  fChannelList.push_back("0hipu");
  //  fChannelList.push_back("1hipu");

  // -- close to other PV
  fChannelList.push_back("0cpv");
  //  fChannelList.push_back("1cpv");

  // -- far to other PV
  fChannelList.push_back("0fpv");
  //  fChannelList.push_back("1fpv");

  // -- small fls3d
  fChannelList.push_back("0sfl");
  //  fChannelList.push_back("1sfl");

  // -- big fls3d
  fChannelList.push_back("0bfl");
  //  fChannelList.push_back("1bfl");

  // -- large dzmin
  fChannelList.push_back("0ldz");
  //  fChannelList.push_back("1bfl");

  // -- small dzmin
  fChannelList.push_back("0sdz");
  //  fChannelList.push_back("1bfl");


}


// ----------------------------------------------------------------------
plotReducedOverlays::~plotReducedOverlays() {

}


// ----------------------------------------------------------------------
void plotReducedOverlays::init() {

  system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  system(Form("/bin/rm -f %s/plotSbsHistograms-%d%s.root", fDirectory.c_str(), fYear, fSetup.c_str()));
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  // system(Form("/bin/rm -f %s/sbsctrl*.pdf", fDirectory.c_str()));
  // system(Form("/bin/rm -f %s/sbsctrl/sbsctrl%d%s_ad*_*.pdf", fDirectory.c_str(), fYear, fSetup.c_str()));
  // system(Form("/bin/rm -f %s/adfpy/adfpy*.pdf", fDirectory.c_str()));
  // system(Form("/bin/rm -f %s/adfpy*.pdf", fDirectory.c_str()));
  system(Form("/bin/rm -f %s/overlay%d%s_ad*_*.pdf", fDirectory.c_str(), fYear, fSetup.c_str()));
  system(Form("/bin/rm -f %s/mass_ad*_*.pdf", fDirectory.c_str()));

}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeAll(string what) {

  if (what == "dbx") {
    //    init();
    //    makeSampleOverlay("bdpsikstarData", "bdpsikstarMcComb");
    //    makeSampleOverlay("bupsikData", "bupsikMcComb", "bdt");
    //    makeSampleOverlay("bmmData", "bdmmMcComb", "bdt");
    makeSampleOverlay("bupsikData", "bupsikMcComb", "bdt");
    return;
  }

  if (what == "official") {
    comparePrivateAndOfficial();
  }
  if (string::npos != what.find("sys")) {
    allSystematics();
    return;
  }

  if (string::npos != what.find("years")) {
    overlay3Samples("bdmmMcComb", "plotSbsHistograms-2016BF.root",
		    "bdmmMcComb", "plotSbsHistograms-2012.root",
		    "bdmmMcComb", "plotSbsHistograms-2011.root",
		    "bdt"
		    );
    return;
  }


  if (what == "all") {
    init();

    printCuts(cout);
    // -- data vs combined MC
    makeSampleOverlay("bmmData", "bdmmMcComb");
    makeSampleOverlay("bupsikData", "bupsikMcComb");
    makeSampleOverlay("bspsiphiData", "bspsiphiMcComb");
    //    makeSampleOverlay("bdpsikstarData", "bdpsikstarMcComb");

    allSystematics();

    // -- validation of private MC vs official MC
    if (2016 == fYear) makeSampleOverlay("bupsikMc", "bupsikMcOff");
  }

  if (string::npos != what.find("plot")) {
    string hfname = Form("%s/plotSbsHistograms-%d%s.root", fDirectory.c_str(), fYear, fSetup.c_str());
    system(Form("/bin/rm -f %s", hfname.c_str()));

    fTEX.close();
    system(Form("/bin/rm -f %s", fTexFileName.c_str()));
    fTEX.open(fTexFileName.c_str(), ios::app);
    system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
    system(Form("/bin/rm -f %s/sbsctrl%d_ad*_*.pdf", fDirectory.c_str(), fYear));
    system(Form("/bin/rm -f %s/overlay%d_ad*_*.pdf", fDirectory.c_str(), fYear));
    system(Form("/bin/rm -f %s/mass%d_ad*_*.pdf", fDirectory.c_str(), fYear));

    printCuts(cout);
    if (fDoCNC) {
      makeOverlay("bspsiphiData", "bspsiphiMcComb", "cnc");
      makeOverlay("bmmData", "bdmmMcComb", "cnc");
      makeOverlay("bupsikData", "bupsikMcComb", "cnc");
      makeOverlay("bdpsikstarData", "bdpsikstarMcComb", "cnc");
    }

    makeOverlay("bspsiphiData", "bspsiphiMcComb", "bdt");
    makeOverlay("bmmData", "bdmmMcComb", "bdt");
    makeOverlay("bupsikData", "bupsikMcComb", "bdt");
    makeOverlay("bdpsikstarData", "bdpsikstarMcComb", "bdt");
  }



  if (0) {
    plotMass("bspsiphiData", "Cu");
    plotMass("bupsikData", "Cu");
    plotMass("bdpsikstarData", "Cu");
    plotMass("bmmData", "Presel");
  }


}


// ----------------------------------------------------------------------
void plotReducedOverlays::comparePrivateAndOfficial() {

  changeSetup(fDirectory, "plotReducedOverlays", "comparePrivateAndOfficial");
  makeSampleOverlay("bupsikMc", "bupsikMcOff", "Ao");

}

// ----------------------------------------------------------------------
void plotReducedOverlays::makeSampleOverlay(string sample1, string sample2, string selection) {
  makeSample(sample1);
  makeSample(sample2);

  // makeSample(sample1, 1.e6);
  // makeSample(sample2, 1.e6);

  fStampString = "nada";
  if (fDoCNC) makeOverlay(sample1, sample2, "cnc");
  makeOverlay(sample1, sample2, "bdt");
}


// ----------------------------------------------------------------------
void plotReducedOverlays::plotMass(string sample, string selection) {

  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << "fHistFileName: " << fHistFileName << " opened " << endl;

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TH1D *h(0);
  fIF->fLo = 5.0;
  fIF->fHi = 5.5;
  TF1 *lBg = fIF->expoErr(fIF->fLo, fIF->fHi);
  string header;
  if (string::npos != sample.find("bspsiphi")) header = "B_{s} #rightarrow J/#psi #phi";
  else if (string::npos != sample.find("bupsik")) header = "B^{+} #rightarrow J/#psi K^{+}";
  else if (string::npos != sample.find("bdpsikstar")) header = "B^{0} #rightarrow J/#psi K^{*}";
  else if (string::npos != sample.find("mm")) header = "Dimuon";
  tl->SetTextFont(42);

  double xPos(0.58);

  vector<string> cuts;
  cuts.push_back("bdt");
  if (fDoCNC) cuts.push_back("cnc");

  double preco(5.16), sigma(0.03);

  string hname("");
  for (unsigned int iv = 0; iv < cuts.size(); ++iv) {
    if (string::npos != cuts[iv].find("cnc")) fStampString = "CNC";
    if (string::npos != cuts[iv].find("bdt")) fStampString = "BDT";

    for (unsigned int i = 0; i < fNchan; ++i) {
      hname = Form("ad%d%s_%s_tauMass%s", i, cuts[iv].c_str(), sample.c_str(), selection.c_str());
      h = (TH1D*)gDirectory->Get(hname.c_str());
      if (!h) break;
      TF1 *f1 = fIF->expoErrGauss(h, 5.28, sigma, preco);
      setTitles(h, "#it{m} [GeV]", Form("Entries / %4.3f GeV", h->GetBinWidth(1)), 0.05, 1.1, 1.5);
      if (header != "Dimuon") {
	h->Fit(f1, "", "e");

	for (int i = 0; i < lBg->GetNpar(); ++i) {
	  lBg->SetParameter(i, f1->GetParameter(3+i));
	}

	double c  = f1->GetParameter(0);
	cout << "OVERALL INTEGRAL: " << f1->Integral(5.15, 5.45) << " BACKGROUND INTEGRAL: " << lBg->Integral(5.15, 5.45) << endl;
	c = f1->Integral(5.15, 5.45) - lBg->Integral(5.15, 5.45);
	double cE = f1->GetParError(0);
	double ierr = f1->IntegralError(5.15, 5.45)/h->GetBinWidth(1);

	double signal = c/h->GetBinWidth(1);
	double signalE(0.);
	if (ierr > TMath::Sqrt(signal)) {
	  signalE = ierr;
	} else {
	  signalE = cE/c*signal;
	}
	tl->SetTextSize(0.025);
	tl->DrawLatexNDC(xPos, 0.70, Form("Signal: %5.1f  #pm %5.1f", signal, signalE));
	tl->DrawLatexNDC(xPos, 0.66, Form("Mass:   %5.4f  #pm %5.4f GeV", f1->GetParameter(1), f1->GetParError(1)));
	tl->DrawLatexNDC(xPos, 0.62, Form("Width:  %5.4f  #pm %5.4f GeV", f1->GetParameter(2), f1->GetParError(2)));
      } else {
	h->Draw("hist");
      }

      tl->SetTextSize(0.05);
      tl->DrawLatexNDC(xPos, 0.80, header.c_str());
      tl->SetTextSize(0.025);
      tl->DrawLatexNDC(xPos, 0.75, Form("%2.1f < |#eta(#mu_{f})| < %2.1f", fCuts[i]->metaMin, fCuts[i]->metaMax));

      if (1) {
	TH1D *hMassBGL    = (TH1D*)gDirectory->Get(Form("ad%d%s_%s_tauMassBGL", i, cuts[iv].c_str(), sample.c_str()));
	TH1D *hMassBGH    = (TH1D*)gDirectory->Get(Form("ad%d%s_%s_tauMassBGH", i, cuts[iv].c_str(), sample.c_str()));
	TH1D *hMassSG     = (TH1D*)gDirectory->Get(Form("ad%d%s_%s_tauMassSG", i, cuts[iv].c_str(), sample.c_str()));
	TArrow aa;
	double ymax = h->GetMaximum();
	double x0 = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
	double y0 = 0.2*ymax;
	double x1 = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
	double y1 = 0.2*ymax;

	double x2 = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
	double y2 = 0.2*ymax;
	double x3 = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
	double y3 = 0.2*ymax;

	double x4 = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
	double y4 = 0.2*ymax;

	double x5 = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
	double y5 = 0.2*ymax;

	aa.SetLineWidth(2);
	aa.SetLineColor(kRed);
	aa.DrawArrow(x0, y0, x0, 0.);
	aa.DrawArrow(x1, y1, x1, 0.);
	aa.SetLineColor(kBlack);
	aa.DrawArrow(x2, y2, x2, 0.);
	aa.DrawArrow(x3, y3, x3, 0.);

	aa.SetLineColor(kBlue);
	aa.DrawArrow(x4, y4, x4, 0.);
	aa.DrawArrow(x5, y5, x5, 0.);
      }


      stamp(0., fStampCms, fStampString, 0., fStampLumi);

      c0->SaveAs(Form("%s/mass%d%s_ad%d%s_%s_%s.pdf", fDirectory.c_str(), fYear, fSetup.c_str(),
		      i, cuts[iv].c_str(), sample.c_str(), selection.c_str()));
    }
  }
  fHistFile->Close();
  cout << "fHistFileName: " << fHistFileName << " closed " << endl;
}



// ----------------------------------------------------------------------
void plotReducedOverlays::makeOverlay(string sample1, string sample2, string selection) {
  vector<string> cuts;
  if (string::npos != selection.find("bdt")) {
    fStampString = "BDT";
    cuts.push_back("HLT");
    cuts.push_back("Presel");
  }
  if (string::npos != selection.find("cnc")) {
    fStampString = "CNC";
    cuts.push_back("Ao");
    cuts.push_back("Presel");
  }

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "plotReducedOverlays::makeOverlay(" << sample1 << ", " << sample2 << ", " << selection << ") fStampString = " << fStampString << endl;
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "fHistFileName: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  for (unsigned int iv = 0; iv < cuts.size(); ++iv) {
    for (unsigned int i = 0; i < fChannelList.size(); ++i) {
      cout << "===> sbsDistributions(" << Form("ad%s%s_%s", fChannelList[i].c_str(), selection.c_str(), sample1.c_str()) << ", " << cuts[iv] << ")" << endl;
      sbsDistributions(Form("ad%s%s_%s", fChannelList[i].c_str(), selection.c_str(), sample1.c_str()), cuts[iv]);
      sbsDistributions(Form("ad%s%s_%s", fChannelList[i].c_str(), selection.c_str(), sample2.c_str()), cuts[iv]);

      overlay(Form("ad%s%s_%s", fChannelList[i].c_str(), selection.c_str(), sample1.c_str()),
	      Form("ad%s%s_%s", fChannelList[i].c_str(), selection.c_str(), sample2.c_str()),
	      cuts[iv]);
    }
  }


  // -- dump the histograms
  string hfname = Form("%s/plotSbsHistograms-%d%s.root", fDirectory.c_str(), fYear, fSetup.c_str());
  TFile *fl = TFile::Open(hfname.c_str(), "UPDATE");
  cout << "+++ saving sbs histograms into " << hfname << endl;
  TH1D *h1(0), *h2(0);
  string n1, n2;
  for (unsigned int iv = 0; iv < cuts.size(); ++iv) {
    for (unsigned int ic = 0; ic < fChannelList.size(); ++ic) {
      for (unsigned int i = 0; i < fDoList.size(); ++i) {
	n1 = Form("sbs_%s_%s%s",
		  Form("ad%s%s_%s", fChannelList[ic].c_str(), selection.c_str(), sample1.c_str()),
		  fDoList[i].c_str(),
		  cuts[iv].c_str());
	n2 = Form("sbs_%s_%s%s",
		  Form("ad%s%s_%s", fChannelList[ic].c_str(), selection.c_str(), sample2.c_str()),
		  fDoList[i].c_str(),
		  cuts[iv].c_str());
	h1 = (TH1D*)fHistFile->Get(n1.c_str());
	cout << "n1: " << n1 << " -> " << h1 << endl;

	h2 = (TH1D*)fHistFile->Get(n2.c_str());
	cout << "n2: " << n2 << " -> " << h2 << endl;
	if (0 == h1 || 0 == h2) {
	  cout << "  histograms not found" << endl;
	  continue;
	}
	h1->SetDirectory(fl);
	h1->Write();
	h2->SetDirectory(fl);
	h2->Write();
      }
    }
  }
  fl->Close();

  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeOverlay2Channels(string sample, string channel1, string channel2, string selection) {
}


// ----------------------------------------------------------------------
void plotReducedOverlays::compareBsAndBp(string file) {
  overlay2Files(file.c_str(), "NoMc", file.c_str(), "CsMc", "B", "B");
  overlay2Files(file.c_str(), "NoData", file.c_str(), "CsData", "B", "B");
}



// ----------------------------------------------------------------------
void plotReducedOverlays::makeSample(string sample, int nevents, int nstart) {

  string dir("");

  // -- dump histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;

  fSample = sample;

  if (string::npos != fSample.find("Mc")) {
    fIsMC = true;
  } else {
    fIsMC = false;
  }

  MASSMIN = 4.5;
  MASSMAX = 6.5;
  if ((string::npos != fSample.find("bdmm")) || (string::npos != fSample.find("bmm")) || (string::npos != fSample.find("bsmm"))) {
    dir = "candAnaMuMu";
    fMode = BMM;
    fIsSignal = true;
    BGLBOXMIN = 4.90;
    BGLBOXMAX = 5.20;
    SIGBOXMIN = 5.20;
    SIGBOXMAX = 5.45;
    BGHBOXMIN = 5.45;
    BGHBOXMAX = 6.00;
  }

  if (string::npos != fSample.find("bupsik")) {
    dir = "candAnaBu2JpsiK";
    fMode = BU2JPSIKP;

    BGLBOXMIN = 5.00;
    BGLBOXMAX = 5.10;

    SIGBOXMIN = 5.18;
    SIGBOXMAX = 5.38;

    BGHBOXMIN = 5.60;
    BGHBOXMAX = 6.00;
  }

    if (string::npos != fSample.find("bdpsikstar")) {
    dir = "candAnaBd2JpsiKstar";
    fMode = BD2JPSIKSTAR;
    BGLBOXMIN = 5.00;
    BGLBOXMAX = 5.10;

    SIGBOXMIN = 5.18;
    SIGBOXMAX = 5.38;

    BGHBOXMIN = 5.60;
    BGHBOXMAX = 6.00;
  }

  if (string::npos != fSample.find("bspsiphi")) {
    dir = "candAnaBs2JpsiPhi";
    fMode = BS2JPSIPHI;
    BGLBOXMIN = 5.00;
    BGLBOXMAX = 5.20;

    SIGBOXMIN = 5.27;
    SIGBOXMAX = 5.47;

    BGHBOXMIN = 5.60;
    BGHBOXMAX = 6.00;
  }

  // -- must be after the mass box definitions!
  if (fDoCNC) bookDistributions("cnc");
  bookDistributions("bdt");

  TTree *t = getTree(fSample, dir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  }

  // nevents = 2e6;
  // nstart = 0;
  setupTree(t, fSample);
  fCds = fDS[fSample];
  loopOverTree(t, 1, nevents, nstart);

  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotReducedOverlays::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotReducedOverlays::loopOverTree> loop over dataset " << fCds->fName << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotReducedOverlays::*pF)(void);
  if (ifunc == 1) pF = &plotReducedOverlays::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}



// ----------------------------------------------------------------------
void plotReducedOverlays::loopFunction1() {

  // -- modify here the fGoodHLT to accomodate the preselections in HFBmm_cff.py, but not in the truth-based MC candidates!
  bool cmsswPresel =  (fb.m1pt>4.) && (fb.m2pt>4.)  && (fb.pvips < 5.) && (fb.flsxy > 4) && (fb.maxdoca < 0.08) && (fb.chi2 < 10.);

  if (fMode == BU2JPSIKP) {
    cmsswPresel = cmsswPresel && (fb.kpt > 0.6);
  }

  if (fMode == BS2JPSIPHI) {
    cmsswPresel = cmsswPresel && (fb.k1pt > 0.6) && (fb.k2pt > 0.6);
  }

  if (fMode == BD2JPSIKSTAR) {
    cmsswPresel = cmsswPresel && (fb.kpt > 0.6) && (fb.pipt > 0.6);
  }

  fGoodHLT = fb.hlt1 && fb.tos && fb.l1t
    && fGoodAcceptance
    && fGoodMuonsID && fGoodMuonsPt && fGoodMuonsEta && fGoodTracksPt && fGoodTracksEta
    && cmsswPresel
    && fGoodJpsiCuts
    && fGoodDcand;

  // -- update ana cuts!
  double bdtCut(0.1);
  if (2011 == fYear) bdtCut = 0.15;
  if (2012 == fYear) bdtCut = 0.20;
  //  fPreselection    = (fGoodHLT && (fb.alpha < 0.1) && (fb.fls3d > 6) && (fb.pvips < 4) && (fb.docatrk < 0.2));
  //  fPreselectionBDT = (fGoodHLT && (fBDT > 0.1));
  fPreselection    = (fGoodHLT && fBDT > fCuts[fChan]->bdtCut);
  fPreselectionBDT = fPreselection;

  fCncCuts.update();
  fBdtCuts.update();

  bool loPU = (fb.pvn <  9);
  bool hiPU = (fb.pvn > 24);

  bool closePV = (TMath::Abs(fb.pv2lip) < 0.1);
  bool farPV   = (TMath::Abs(fb.pv2lip) > 1.2);

  bool sfl = (fb.fls3d < 12.);
  bool bfl = (fb.fls3d > 50.);

  bool ldz = (fb.dzmin > 1.2);
  bool sdz = (fb.dzmin < 0.08);

  if (0) {
    cout
      << " presel = " << fPreselection
      << " preselBDT = " << fPreselectionBDT
      << " cmsswPresel = " << cmsswPresel
      << " goodHLT = " << fGoodHLT
      << " hlt1 = " << fb.hlt1
      << " tos " << fb.tos
      << " l1t = " << fb.l1t
      << " acc = " << fGoodAcceptance
      << " gt = " << fGoodTracks
      << " trkPt = " << fGoodTracksPt
      << " trkEta = " << fGoodTracksEta
      << " gmID = " << fGoodMuonsID
      << " gmPt = " << fGoodMuonsPt
      << " gmEta = " << fGoodMuonsEta
      << " gbdtPt = " << fGoodBdtPt
      << " Q = " << fGoodQ
      << " j/psi = " << fGoodJpsiCuts
      << " dcand = " << fGoodDcand
      << " pv = " << fGoodPvAveW8
      << " BDT = " << fBDT
      << endl;
    if (!cmsswPresel) {
      cout << " m1pt = " << fb.m1pt
	   << " m2pt = " << fb.m2pt
	   << " pvips = " << fb.pvips
	   << " flsxy = " << fb.flsxy
	   << " maxdoca = " << fb.maxdoca
	   << " chi2 = " << fb.chi2
	   << endl;
    }
  }


  if (fGoodHLT
      && (fBDT > -1.) && (fb.alpha < 0.1) && (fb.fls3d > 7)
      ) {
    fSel0 = true;
  } else {
    fSel0 = false;
  }

  if (fGoodHLT
      && (fBDT > -1.) && (fb.alpha < 0.07) && (fb.fls3d > 10) && (fb.pvips < 2.)
      ) {
    fSel1 = true;
  } else {
    fSel1 = false;
  }

  if (fGoodHLT
      && (fBDT > -1.) && (fb.alpha < 0.05) && (fb.fls3d > 12) && (fb.pvips < 2.)
      ) {
    fSel2 = true;
  } else {
    fSel2 = false;
  }

  if (-1 < fChan && fChan < fNchan) {
    fChannel = Form("%d", fChan);
    if (fDoCNC) fillDistributions("cnc");
    fillDistributions("bdt");

    if (0 == fChan) {
      if (loPU) {
      	fChannel = Form("%dlopu", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }
      if (hiPU) {
      	fChannel = Form("%dhipu", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }

      if (closePV) {
      	fChannel = Form("%dcpv", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }
      if (farPV) {
      	fChannel = Form("%dfpv", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }

      if (sfl) {
      	fChannel = Form("%dsfl", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }
      if (bfl) {
      	fChannel = Form("%dbfl", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }

      if (sdz) {
      	fChannel = Form("%dsdz", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }
      if (ldz) {
      	fChannel = Form("%dldz", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }

    }
  }
}


// ----------------------------------------------------------------------
void plotReducedOverlays::bookDistributions(std::string selmode) {

  cout << "SIG: " << SIGBOXMIN << " .. " << SIGBOXMAX << endl;
  cout << "BGL: " << BGLBOXMIN << " .. " << BGLBOXMAX << endl;
  cout << "BGH: " << BGHBOXMIN << " .. " << BGHBOXMAX << endl;

  adset *a(0);
  AnalysisCuts *pCuts = &fCncCuts;
  bool bdt(false);
  bool *p = &fPreselection;
  if (string::npos != selmode.find("bdt")) {
    pCuts = &fBdtCuts;
    p     = &fPreselectionBDT;
    bdt = true;
  }
  for (unsigned int i = 0; i < fChannelList.size(); ++i) {
    string mapname = Form("ad%s%s_%s", fChannelList[i].c_str(), selmode.c_str(), fSample.c_str());
    string name = Form("%s_", mapname.c_str());

    a = new adset();
    a->fpPvN      = bookDistribution(Form("%spvn", name.c_str()), "N(PV) ", "fGoodHLT", pCuts, 40, 0., 40., p);
    a->fpPvNtrk   = bookDistribution(Form("%spvntrk", name.c_str()), "N_{trk}(PV) ", "fGoodHLT", pCuts, 40, 0., 120., p);
    a->fpPv2Ntrk  = bookDistribution(Form("%spv2ntrk", name.c_str()), "N_{trk}(PV2) ", "fGoodHLT", pCuts, 40, 0., 120., p);
    a->fpDzmin    = bookDistribution(Form("%sdzmin", name.c_str()), "min(#Delta z) [cm] ", "fGoodHLT", pCuts, 50, -2., 2., p);
    a->fpDz12     = bookDistribution(Form("%sdz12", name.c_str()), "#Delta z [cm] ", "fGoodHLT", pCuts, 50, -2., 2., p);

    a->fpPvZ      = bookDistribution(Form("%spvz", name.c_str()), "z_{PV} [cm]", "fGoodHLT", pCuts, 40, -20., 20., p);
    a->fpPvAveW8  = bookDistribution(Form("%spvavew8", name.c_str()), "<w^{PV}>", "fGoodPvAveW8", pCuts, 50, 0.5, 1., p);

    if (bdt) {
      a->fpBDT      = bookDistribution(Form("%sbdt", name.c_str()), "BDT", "fGoodEta", pCuts, 100, -1.0, 1.0, p);
      a->fpBDTSel0  = bookDistribution(Form("%sbdtsel0", name.c_str()), "BDTsel0", "fGoodEta", pCuts, 100, -1.0, 1.0, &fSel0);
      a->fpBDTSel1  = bookDistribution(Form("%sbdtsel1", name.c_str()), "BDTsel1", "fGoodEta", pCuts, 100, -1.0, 1.0, &fSel1);
      a->fpBDTSel2  = bookDistribution(Form("%sbdtsel2", name.c_str()), "BDTsel2", "fGoodEta", pCuts, 100, -1.0, 1.0, &fSel2);
    } else {
      a->fpBDT      = bookDistribution(Form("%sbdt", name.c_str()), "BDT", "fGoodCloseTrack", pCuts, 100, -1.0, 1.0, p);
      a->fpBDTSel0  = bookDistribution(Form("%sbdtsel0", name.c_str()), "BDTsel0", "fGoodMaxDoca", pCuts, 100, -1.0, 1.0, p);
      a->fpBDTSel1  = bookDistribution(Form("%sbdtsel1", name.c_str()), "BDTsel1", "fGoodFLS", pCuts, 100, -1.0, 1.0, p);
      a->fpBDTSel2  = bookDistribution(Form("%sbdtsel2", name.c_str()), "BDTsel2", "fGoodIpS", pCuts, 100, -1.0, 1.0, p);
    }

    a->fpMuon1Pt  = bookDistribution(Form("%smuon1pt", name.c_str()), "#it{p}_{T, #mu1} [GeV]", "fGoodMuonsPt", pCuts, 60, 0., 30., p);
    a->fpMuon2Pt   = bookDistribution(Form("%smuon2pt", name.c_str()), "#it{p}_{T, #mu2} [GeV]", "fGoodMuonsPt", pCuts, 40, 0., 20., p);
    a->fpMuonsEta  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsEta", pCuts, 40, -2.5, 2.5, p);

    a->fpPt        = bookDistribution(Form("%spt", name.c_str()), "#it{p}_{T}#it{(B)} [GeV]", "fGoodPt", pCuts, 60, 0., 60., p);
    a->fpP         = bookDistribution(Form("%sp", name.c_str()), "#it{p(B)} [GeV]", "fGoodPt", pCuts, 50, 0., 100., p);
    a->fpPz        = bookDistribution(Form("%spz", name.c_str()), "#it{p_{z}(B)} [GeV]", "fGoodPt", pCuts, 50, 0., 100., p);
    a->fpEta       = bookDistribution(Form("%seta", name.c_str()), "#eta#it{(B)}", "fGoodEta", pCuts, 40, -2.5, 2.5, p);

    a->fpChi2Dof   = bookDistribution(Form("%schi2dof", name.c_str()),  "#chi^{2}/dof", (bdt? "fGoodBDT": "fGoodChi2"), pCuts, 30, 0., 3., p);
    a->fpPChi2Dof  = bookDistribution(Form("%spchi2dof", name.c_str()),  "#it{P(#chi^{2},dof)}", (bdt? "fGoodBDT": "fGoodChi2"), pCuts, 50, 0., 1.0, p);
    a->fpPvDchi2   = bookDistribution(Form("%spvdchi2", name.c_str()),  "#Delta(#chi^{2})", (bdt? "fGoodBDT": "fGoodChi2"), pCuts, 100, 0., 2000., p);
    a->fpOtherVtx  = bookDistribution(Form("%sothervtx", name.c_str()),  "othervtx", (bdt? "fGoodBDT": "fGoodChi2"), pCuts, 40, 0., 1., p);
    a->fpMaxDoca   = bookDistribution(Form("%smaxdoca", name.c_str()), "#it{d}^{max}_{ca} [cm]", (bdt? "fGoodBDT": "fGoodMaxDoca"), pCuts, 60, 0., 0.03, p);

    a->fpAlpha     = bookDistribution(Form("%salpha", name.c_str()), "#alpha_{3D}", (bdt? "fGoodBDT": "fGoodAlpha"), pCuts, 50, 0., 0.1, p);

    a->fpFLS3d     = bookDistribution(Form("%sfls3d", name.c_str()), "#it{l_{3D}}/#sigma(#it{l_{3D}})", (bdt? "fGoodBDT": "fGoodFLS"), pCuts, 60, 0., 120., p);
    a->fpFL3d      = bookDistribution(Form("%sfl3d", name.c_str()),  "#it{l_{3D}} [cm]", (bdt? "fGoodBDT": "fGoodFLS"), pCuts, 60, 0., 1.5, p);
    a->fpFL3dE     = bookDistribution(Form("%sfl3de", name.c_str()), "#sigma(#it{l_{3D}}) [cm]", (bdt? "fGoodBDT": "fGoodFLS"), pCuts, 50, 0., 0.05, p);

    a->fpFLSxy     = bookDistribution(Form("%sflsxy", name.c_str()), "#it{l_{xy}}/#sigma(#it{l_{xy}})", (bdt? "fGoodBDT": "fGoodFLS"),
				      pCuts, 60, 0., 120., p);
    a->fpFLxy      = bookDistribution(Form("%sflxy", name.c_str()),  "#it{l_{xy}} [cm]", (bdt? "fGoodBDT": "fGoodFLS"),
				      pCuts, 60, 0., 1.5, p);
    a->fpFLxyE     = bookDistribution(Form("%sflxye", name.c_str()), "#sigma(#it{l_{xy}}) [cm]", (bdt? "fGoodBDT": "fGoodFLS"),
				      pCuts, 50, 0., 0.05, p);

    a->fpIp        = bookDistribution(Form("%sip", name.c_str()), "#delta_{#it{3D}} [cm]", (bdt? "fGoodBDT": "fGoodIp"),
				      pCuts, 60, 0., 0.012, p);
    a->fpIpS       = bookDistribution(Form("%sips", name.c_str()), "#delta_{#it{3D}}/#sigma(#delta_{3D})", (bdt? "fGoodBDT": "fGoodIpS"),
				      pCuts, 40, 0., 4, p);

    a->fpLip       = bookDistribution(Form("%slip", name.c_str()), "|#it{l_{z}}| [cm]", (bdt? "fGoodBDT": "fGoodLip"),
				      pCuts, 50, 0., 0.015, p);
    a->fpLip2      = bookDistribution(Form("%slip2", name.c_str()), "|#it{l_{z}^{(2)}}| [cm]", (bdt? "fGoodBDT": "fGoodLip"),
				      pCuts, 100, 0., 4.0, p);

    a->fpLipS      = bookDistribution(Form("%slips", name.c_str()), "|#it{l_{z}}|/#sigma(#it{l_{z}})", (bdt? "fGoodBDT": "fGoodLipS"),
				      pCuts, 50, 0., 4, p);
    a->fpLipS2     = bookDistribution(Form("%slips2", name.c_str()), "|#it{l_{z}}^{(2)}|/#sigma(#it{l_{z}}^{(2)})", (bdt? "fGoodBDT": "fGoodLipS"),
				      pCuts, 50, 0., 20., p);


    a->fpDocaTrk   = bookDistribution(Form("%sdocatrk", name.c_str()), "#it{d}_{ca}^{0} [cm]", (bdt? "fGoodBDT": "fGoodDocaTrk"), pCuts, 50, 0., 0.20, p);
    a->fpIso       = bookDistribution(Form("%siso", name.c_str()),  "isolation", (bdt? "fGoodBDT": "fGoodIso"), pCuts, 26, 0., 1.04, p);
    a->fpM1Iso     = bookDistribution(Form("%sm1iso", name.c_str()),  "m1 isolation", (bdt? "fGoodBDT": "fGoodM1Iso"), pCuts, 26, 0., 1.04, p);
    a->fpM2Iso     = bookDistribution(Form("%sm2iso", name.c_str()),  "m2 isolation", (bdt? "fGoodBDT": "fGoodM2Iso"), pCuts, 26, 0., 1.04, p);

    a->fpCloseTrk  = bookDistribution(Form("%sclosetrk", name.c_str()),  "#it{N}_{trk}^{close}", (bdt? "fGoodBDT": "fGoodCloseTrack"),
				      pCuts, 10, 0., 10., p);
    // a->fpCloseTrkS1= bookDistribution(Form("%sclosetrks1", name.c_str()),  "N_{trk}^{close, 1#sigma}",
    // (bdt? "fGoodBDT": "fGoodCloseTrackS1"), pCuts, 10, 0., 10., p);
    // a->fpCloseTrkS2= bookDistribution(Form("%sclosetrks2", name.c_str()),  "N_{trk}^{close, 2#sigma}",
    // (bdt? "fGoodBDT": "fGoodCloseTrackS2"), pCuts, 10, 0., 10., p);
    // a->fpCloseTrkS3= bookDistribution(Form("%sclosetrks3", name.c_str()),  "N_{trk}^{close, 3#sigma}",
    // (bdt? "fGoodBDT": "fGoodCloseTrackS3"), pCuts, 10, 0., 10., p);

    a->fpTau       = bookDistribution(Form("%stau", name.c_str()), "#tau [ps]", (bdt? "fGoodBDT": "fGoodCloseTrack"), pCuts, 50, 0., 10., p);

    fAdMap.insert(make_pair(mapname, a));
    cout << "bookDistributions: mapname = " << mapname << endl;

  }

}


// ----------------------------------------------------------------------
void plotReducedOverlays::sbsDistributions(string sample, string selection, string what) {
  string sbsControlPlotsFileName = Form("sbsctrl%d%s", fYear, fSetup.c_str());
  string sid = sample.substr(0, sample.find("_"));
  replaceAll(sid, "bdt", "");
  replaceAll(sid, "cnc", "");
  sid = sid.substr(sid.find("ad"));

  AnalysisDistribution a(Form("%s_muon1pt", sample.c_str()));
  a.fVerbose = 1;
  a.fControlPlotsFileName = sbsControlPlotsFileName;
  a.fDirectory = fDirectory + string("/sbsctrl");

  int type(0);
  if (string::npos != sample.find("bmmData"))        type = 0; // sidebands
  if (string::npos != sample.find("Mc"))             type = 1; // signal window
  if (string::npos != sample.find("bupsikData"))     type = 2; // expo+err2
  if (string::npos != sample.find("bspsiphiData"))   type = 3; // dedicated
  if (string::npos != sample.find("bdpsikstarData")) type = 4; // dedicated

  cout << "gDIRECTORY: "; gDirectory->pwd();
  TH1D *h(0);
  bool restricted = (what != "");
  string bla;
  string lselection(selection);
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    bla =  Form("%s_%s", sample.c_str(), fDoList[i].c_str());
    if (0 == type) {
      cout << "=> Looking for sideband histogram " << Form("%s%s1", bla.c_str(), lselection.c_str()) << endl;
      //      h = (TH1D*)gDirectory->Get(Form("%s%s1", bla.c_str(), lselection.c_str()));
      h = a.sbDistribution(bla.c_str(), lselection.c_str(), 2);
      cout << "=> cloning into  "
	   << Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), lselection.c_str()) << endl;
      h = (TH1D*)h->Clone(Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), lselection.c_str()));
    } else if (1 == type) {
      cout << "=> Looking for signal histogram " << Form("%s%s0", bla.c_str(), lselection.c_str()) << endl;
      h = (TH1D*)gDirectory->Get(Form("%s%s0", bla.c_str(), lselection.c_str()));
      cout << "=> cloning into  "
	   << Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), lselection.c_str()) << endl;
      h = (TH1D*)h->Clone(Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), lselection.c_str()));
    } else if (2 == type) {
      cout << "=> sbsDistributionExpoErrGauss histogram " << Form("%s for lselection %s", bla.c_str(), lselection.c_str())
	   << endl;
      h = a.sbsDistributionExpoErrGauss(bla.c_str(), lselection.c_str());
    } else if (3 == type) {
      cout << "=> sbsDistributionBs2JpsiPhi histogram " << Form("%s for lselection %s", bla.c_str(), lselection.c_str()) << endl;
      h = a.sbsDistributionBs2JpsiPhi(bla.c_str(), lselection.c_str());
    } else if (4 == type) {
      cout << "=> sbsDistributionBd2JpsiKstar histogram " << Form("%s for lselection %s", bla.c_str(), lselection.c_str()) << endl;
      h = a.sbsDistributionBd2JpsiKstar(bla.c_str(), lselection.c_str());
    }
    cout << "  Title: " << h->GetTitle() << " with integral: " << h->GetSumOfWeights() << endl;
  }
}


// ----------------------------------------------------------------------
void plotReducedOverlays::allSystematics() {

  // -- BDT cut systematics
  sysBdtCut("bupsikData", "bupsikMcComb", "bdt");
  sysBdtCut("bspsiphiData", "bspsiphiMcComb", "bdt");
  sysBdtCut("bdpsikstarData", "bdpsikstarMcComb", "bdt");
}


// ----------------------------------------------------------------------
void plotReducedOverlays::sysBdtCut(string sample1, string sample2, string selection, string file2) {

  c0->SetCanvasSize(700, 700);
  c0->cd();
  shrinkPad(0.15, 0.18);

  string hfname = Form("%s/plotSbsHistograms-%d%s.root", fDirectory.c_str(), fYear, fSetup.c_str());
  TFile *f1 = TFile::Open(hfname.c_str());
  TFile *f2 = f1;
  if (file2 != "nada") {
    Form("%s/plotSbsHistograms-%s.root", fDirectory.c_str(), file2.c_str());
    f2 = TFile::Open(hfname.c_str());
  }

  vector<string> dolist;
  dolist.push_back("bdt");
  dolist.push_back("bdtsel0");
  dolist.push_back("bdtsel1");
  dolist.push_back("bdtsel2");
  vector<string> cutlevel;
  cutlevel.push_back("Presel");
  //  cutlevel.push_back("HLT");

  TH1D *h1(0), *h2(0);
  double n1(0.), n2(0.), N1(0.), N2(0.);
  double eps1(0.), eps1E(0.), eps2(0.), eps2E(0.), deps(0.), depsE(0.);
  int lo(-1), hi(-1);
  string hname("");

  TArrow aa;
  for (int id = 0; id < dolist.size(); ++id) {
    for (int ic = 0; ic < cutlevel.size(); ++ic) {
      for (int i = 0; i < fChannelList.size(); ++i) {
	hname = Form("sbs_ad%s%s_%s_%s%s", fChannelList[i].c_str(), selection.c_str(), sample1.c_str(), dolist[id].c_str(), cutlevel[ic].c_str());
	h1 = (TH1D*)f1->Get(hname.c_str());
	hname = Form("sbs_ad%s%s_%s_%s%s", fChannelList[i].c_str(), selection.c_str(), sample2.c_str(), dolist[id].c_str(), cutlevel[ic].c_str());
	h2 = (TH1D*)f2->Get(hname.c_str());
	if (!h1 || !h2) {
	  cout << "histogram(s) " << hname << " not found, h1/h2 = " << h1 << "/" << h2 << endl;
	  continue;
	}
	int ichan(0);
	if (string::npos != fChannelList[i].find("1")) ichan = 1;
	lo    = h1->FindBin(fCuts[ichan]->bdtCut);
	hi    = h1->GetNbinsX();
	n1    = h1->Integral(lo, hi);
	N1    = h1->Integral(0, hi);
	eps1  = n1/N1;
	n2    = h2->Integral(lo, hi);
	N2    = h2->Integral(0, hi);
	eps2  = n2/N2;
	deps  = 1. - (eps1/eps2);
	eps1E = dRatio(n1, N1);
	eps2E = dRatio(n2, N2);
	depsE = dRatio(eps1, eps1E, eps2, eps2E);

	h1->Scale(1./h1->GetSumOfWeights());
	h2->Scale(1./h2->GetSumOfWeights());
	h1->SetMinimum(0.);
	setHist(h1, fDS[sample1]);
	h1->Draw();
	setHist(h2, fDS[sample2]);
	h2->Draw("samehist");
	cout << " chan " << ichan << " bdt cut: " << fCuts[ichan]->bdtCut
	     << " eps1 = " << eps1 << " eps2 = " << eps2
	     << " -> del diff  = " << deps
	     << " +/- " << depsE
	     << endl;

	aa.DrawArrow(fCuts[ichan]->bdtCut, 0.3*h1->GetMaximum(), fCuts[ichan]->bdtCut, 0.);
	tl->SetTextSize(0.025);
	tl->SetTextAngle(0.);
	tl->DrawLatexNDC(0.15, 0.92,
			 Form("#varepsilon(%s) = %4.3f, #varepsilon(%s) = %4.3f, #Delta = %4.3f #pm %4.3f",
			      sample1.c_str(), eps1, sample2.c_str(), eps2, deps, depsE));
	tl->SetTextSize(0.03);
	tl->DrawLatexNDC(0.25, 0.80, Form("BDT > %4.3f", fCuts[ichan]->bdtCut));

	hname = Form("%d%s:ad%s%s_%s_%s_%s%s", fYear, fSetup.c_str(), fChannelList[i].c_str(), selection.c_str(),
		     sample1.c_str(), sample2.c_str(), dolist[id].c_str(), cutlevel[ic].c_str());
	fTEX << formatTex(eps1,  Form("%s:sysBdt:eps1V", hname.c_str()), 3) << endl;
	fTEX << formatTex(eps1E, Form("%s:sysBdt:eps1E", hname.c_str()), 3) << endl;
	fTEX << formatTex(eps2,  Form("%s:sysBdt:eps2V", hname.c_str()), 3) << endl;
	fTEX << formatTex(eps2E, Form("%s:sysBdt:eps2E", hname.c_str()), 3) << endl;
	fTEX << formatTex(deps,  Form("%s:sysBdt:depsV", hname.c_str()), 3) << endl;
	fTEX << formatTex(depsE, Form("%s:sysBdt:depsE", hname.c_str()), 3) << endl;

	hname = Form("sbso%d%s_ad%s%s_%s_%s_%s%s", fYear, fSetup.c_str(), fChannelList[i].c_str(), selection.c_str(),
		     sample1.c_str(), sample2.c_str(), dolist[id].c_str(), cutlevel[ic].c_str());

	tl->SetTextAngle(90.);
	tl->SetTextSize(0.025);
	tl->DrawLatexNDC(0.92, 0.20, hname.c_str());

	savePad(Form("sbso/%s.pdf", hname.c_str()));
      }
    }
  }
}


// ----------------------------------------------------------------------
void plotReducedOverlays::overlay3Samples(string sample1, string file1,
					  string sample2, string file2,
					  string sample3, string file3,
					  string selection) {

  c0->SetCanvasSize(700, 700);
  c0->cd();
  shrinkPad(0.15, 0.18);

  string label1 = file1;
  replaceAll(label1, ".root", "");
  replaceAll(label1, "plotSbsHistograms-", "");
  string label2 = file2;
  replaceAll(label2, ".root", "");
  replaceAll(label2, "plotSbsHistograms-", "");
  string label3 = file3;
  replaceAll(label3, ".root", "");
  replaceAll(label3, "plotSbsHistograms-", "");

  TFile *f1 = TFile::Open((fDirectory + "/" + file1).c_str());
  TFile *f2 = TFile::Open((fDirectory + "/" + file2).c_str());
  TFile *f3 = TFile::Open((fDirectory + "/" + file3).c_str());

  vector<string> cutlevel;
  cutlevel.push_back("Presel");

  TH1D *h1(0), *h2(0), *h3(0);
  int lo(-1), hi(-1);
  string hname("");

  for (int id = 0; id < fDoList.size(); ++id) {
    for (int ic = 0; ic < cutlevel.size(); ++ic) {
      for (int i = 0; i < fChannelList.size(); ++i) {
	hname = Form("sbs_ad%s%s_%s_%s%s", fChannelList[i].c_str(), selection.c_str(), sample1.c_str(), fDoList[id].c_str(), cutlevel[ic].c_str());
	h1 = (TH1D*)f1->Get(hname.c_str());
	hname = Form("sbs_ad%s%s_%s_%s%s", fChannelList[i].c_str(), selection.c_str(), sample2.c_str(), fDoList[id].c_str(), cutlevel[ic].c_str());
	h2 = (TH1D*)f2->Get(hname.c_str());
	hname = Form("sbs_ad%s%s_%s_%s%s", fChannelList[i].c_str(), selection.c_str(), sample3.c_str(), fDoList[id].c_str(), cutlevel[ic].c_str());
	h3 = (TH1D*)f3->Get(hname.c_str());
	if (!h1 || !h2 || !h3) {
	  cout << "histogram(s) " << hname << " not found, h1/h2 = " << h1 << "/" << h2 << endl;
	  continue;
	}

	h1->Scale(1./h1->GetSumOfWeights());
	h2->Scale(1./h2->GetSumOfWeights());
	h3->Scale(1./h3->GetSumOfWeights());
	h1->SetMinimum(0.);
	double ymax(h1->GetMaximum());
	if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();
	if (h3->GetMaximum() > ymax) ymax = h3->GetMaximum();
	h1->SetMaximum(1.2*ymax);
	setHist(h1, kBlack);
	h1->Draw("hist");
	setHist(h2, kBlue);
	h2->Draw("samehist");
	setHist(h3, kRed);
	h3->Draw("samehist");

	newLegend(0.50, 0.7, 0.75, 0.87);
	legg->SetHeader((selection + " (" + cutlevel[ic] + ")").c_str());
	legg->SetTextSize(0.05);
	legg->AddEntry(h1, label1.c_str(), "l");
	legg->AddEntry(h2, label2.c_str(), "l");
	legg->AddEntry(h3, label3.c_str(), "l");
	legg->Draw();
	hname = Form("mco%d%s_ad%s%s_%s_%s_%s%s", fYear, fSetup.c_str(), fChannelList[i].c_str(), selection.c_str(),
		     sample1.c_str(), sample2.c_str(), fDoList[id].c_str(), cutlevel[ic].c_str());

	savePad(Form("sbso/%s.pdf", hname.c_str()));
      }
    }
  }
}



// ----------------------------------------------------------------------
void plotReducedOverlays::overlay(string sample1, string sample2, string selection, string what) {

  gStyle->SetOptTitle(0);
  c0->SetCanvasSize(700, 700);
  c0->cd();
  shrinkPad(0.15, 0.18);

  string ds1 = sample1.substr(sample1.find("_")+1);
  string ds2 = sample2.substr(sample2.find("_")+1);

  TH1D *h1(0), *h2(0);
  string n1, n2;
  bool restricted = (what != "");
  bool doLegend(true);
  bool leftLegend(false);
  string lselection(selection);
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    // if ((string::npos != fDoList[i].find("bdt")) && (string::npos != sample1.find("cnc"))) {
    //   lselection = "Cu";
    // } else {
    //   lselection = selection;
    // }
    n1 =  Form("sbs_%s_%s%s", sample1.c_str(), fDoList[i].c_str(), lselection.c_str());
    n2 =  Form("sbs_%s_%s%s", sample2.c_str(), fDoList[i].c_str(), lselection.c_str());
    doLegend = true;
    if (string::npos != fDoList[i].find("eta")) doLegend = false;
    if (string::npos != fDoList[i].find("bdt")) doLegend = false;
    h1 = (TH1D*)gDirectory->Get(n1.c_str());
    if (fIncludeOverflowInLastBin) addOverflow(h1);
    cout << "n1: " << n1 << " -> " << h1 << endl;

    h2 = (TH1D*)gDirectory->Get(n2.c_str());
    if (fIncludeOverflowInLastBin) addOverflow(h2);
    cout << "n2: " << n2 << " -> " << h2 << endl;
    if (0 == h1 || 0 == h2) {
      cout << "  histograms not found" << endl;
      continue;
    }
    if (h2->GetSumOfWeights() > 0) {
      h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
    }

    cout << "setHist for " << ds1 << " and " << ds2 << endl;
    cout << "fStampString = " << fStampString << endl;

    h1->SetNdivisions(505, "X");


    overlayAndRatio(c0, h1, h2);
    setHist(h1, fDS[ds1]);
    setHist(h2, fDS[ds2]);

    //    setHist(h2, fDS[ds2]);

    if (doLegend) {
      if (leftLegend) {
	newLegend(0.21, 0.7, 0.41, 0.87);
      } else {
	newLegend(0.50, 0.7, 0.75, 0.87);
      }

      char loption1[100], loption2[100];
      string header, h1string, h2string;
      if (string::npos != sample1.find("bspsiphi") && string::npos != sample2.find("bspsiphi")) header = "B_{s} #rightarrow J/#psi #phi";
      else if (string::npos != sample1.find("bupsik") && string::npos != sample2.find("bupsik")) header = "B^{+} #rightarrow J/#psi K^{+}";
      else if (string::npos != sample1.find("bdpsikstar") && string::npos != sample2.find("bdpsikstar")) header = "B^{0} #rightarrow J/#psi K^{*}";
      else if (string::npos != sample1.find("mm") && string::npos != sample2.find("mm")) header = "Dimuon";
      else header = "Zoge am Boge";

      if (string::npos != sample1.find("Mc")) {
	sprintf(loption1, "f");
	if (string::npos != sample1.find("Sg")) {
	  h1string = "B_{s} #rightarrow #mu^{+} #mu^{-} (MC)";
	} else {
	  h1string = "MC simulation";
	  // -- change things for the private/official validation
	  if (string::npos == sample1.find("Off") && string::npos != sample2.find("Off")) {
	    h1string += " (private)";
	    sprintf(loption1, "p");
	    h1->SetMarkerStyle(24);
	    h1->SetMarkerSize(1.5);
	    h1->Draw("esame");
	  }
	}
      } else if (string::npos != sample1.find("Data")) {
	sprintf(loption1, "p");
	if (string::npos != sample1.find("mm")) {
	  h1string = "data sidebands";
	} else {
	  h1string = "data";
	}
      } else {
	h1string = "??";
      }

      if (string::npos != sample2.find("Mc")) {
	sprintf(loption2, "f");
	if (string::npos != sample2.find("bdmm")) {
	  h2string = "B^{0} #rightarrow #mu^{+} #mu^{-}";
	} else if (string::npos != sample2.find("bsmm")) {
	  h2string = "B^{0}_{s} #rightarrow #mu^{+} #mu^{-}";
	} else {
	  h2string = "MC simulation";
	  // -- change things for the private/official validation
	  if (string::npos == sample1.find("Off") && string::npos != sample2.find("Off")) {
	    h2string += " (official)";
	  }
	}
      } else if (string::npos != sample2.find("Data")) {
	sprintf(loption2, "p");
	if (string::npos != sample2.find("bmm")) {
	  h2string = "data sidebands";
	} else {
	  h2string = "data";
	}
      } else {
	h2string = "??";
      }

      legg->SetHeader(header.c_str());
      legg->SetTextSize(0.05);
      legg->AddEntry(h1, h1string.c_str(), loption1);
      legg->AddEntry(h2, h2string.c_str(), loption2);

      legg->Draw();
    } else {
      legg = 0;
    }

    stamp(0.18, fStampCms, fStampString, 0.4, fStampLumi);

    if (1) {
      TLatex ll;
      ll.SetTextAngle(90.);
      ll.SetTextSize(0.03);
      ll.DrawLatexNDC(0.97, 0., Form("%s/%s/%s/%s", sample1.c_str(), sample2.c_str(), lselection.c_str(), fDoList[i].c_str()));
    }

    c0->Modified();
    c0->Update();
    c0->SaveAs(Form("%s/overlays/%d/overlay%d%s_%s_%s_%s_%s.pdf",
		    fDirectory.c_str(), fYear, fYear, fSetup.c_str(), sample1.c_str(), sample2.c_str(), fDoList[i].c_str(), lselection.c_str()));



    // -- clean up
    TPad *p = (TPad*)c0->GetPrimitive("pad1");
    delete p;
    p = (TPad*)c0->GetPrimitive("pad2");
    TH1D *h = (TH1D*)(p->GetPrimitive("hr"));
    delete h;
    delete p;
    if (legg) delete legg;


  }
}

// ----------------------------------------------------------------------
void plotReducedOverlays::overlay2Files(std::string file1, std::string sample1,
					std::string file2, std::string sample2,
					std::string chan1, std::string chan2,
					std::string selection, std::string what) {

  string hfname1 = Form("%s/%s.root", fDirectory.c_str(), file1.c_str());
  TFile *f1 = TFile::Open(hfname1.c_str());
  string hfname2 = Form("%s/%s.root", fDirectory.c_str(), file2.c_str());
  TFile *f2 = TFile::Open(hfname2.c_str());

  string ds1 = sample1.substr(sample1.find("_")+1);
  string ds2 = sample2.substr(sample2.find("_")+1);

  string n1(""), n2("");
  TH1D *h1(0), *h2(0);
  bool doLegend(false), leftLegend(false);
  for (int id = 0; id < fDoList.size(); ++id) {
    for (int i = 0; i < fChannelList.size(); ++i) {
      n1 =  Form("sbs_%s_%s%s", sample1.c_str(), fDoList[i].c_str(), selection.c_str());
      n2 =  Form("sbs_%s_%s%s", sample2.c_str(), fDoList[i].c_str(), selection.c_str());
      doLegend = true;
      if (string::npos != fDoList[i].find("eta")) doLegend = false;
      if (string::npos != fDoList[i].find("bdt")) doLegend = false;
      h1 = (TH1D*)f1->Get(n1.c_str());
      if (fIncludeOverflowInLastBin) addOverflow(h1);
      cout << "n1: " << n1 << " -> " << h1 << endl;

      h2 = (TH1D*)f2->Get(n2.c_str());
      if (fIncludeOverflowInLastBin) addOverflow(h2);
      cout << "n2: " << n2 << " -> " << h2 << endl;
      if (0 == h1 || 0 == h2) {
	cout << "  histograms not found" << endl;
	continue;
      }
      if (h2->GetSumOfWeights() > 0) {
	h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
      }

      cout << "setHist for " << ds1 << " and " << ds2 << endl;
      cout << "fStampString = " << fStampString << endl;

      h1->SetNdivisions(505, "X");

      overlayAndRatio(c0, h1, h2);
      setHist(h1, fDS[ds1]);
      setHist(h2, fDS[ds2]);

      if (doLegend) {
	if (leftLegend) {
	  newLegend(0.21, 0.7, 0.41, 0.87);
	} else {
	  newLegend(0.50, 0.7, 0.75, 0.87);
	}

	char loption1[100], loption2[100];
	string header, h1string, h2string;
	if (string::npos != sample1.find("bspsiphi") && string::npos != sample2.find("bspsiphi")) header = "B_{s} #rightarrow J/#psi #phi";
	else if (string::npos != sample1.find("bupsik") && string::npos != sample2.find("bupsik")) header = "B^{+} #rightarrow J/#psi K^{+}";
	else if (string::npos != sample1.find("bdpsikstar") && string::npos != sample2.find("bdpsikstar")) header = "B^{0} #rightarrow J/#psi K^{*}";
	else if (string::npos != sample1.find("mm") && string::npos != sample2.find("mm")) header = "Dimuon";
	else header = "Zoge am Boge";

	if (string::npos != sample1.find("Mc")) {
	  sprintf(loption1, "f");
	  if (string::npos != sample1.find("Sg")) {
	    h1string = "B_{s} #rightarrow #mu^{+} #mu^{-} (MC)";
	  } else {
	    h1string = "MC simulation";
	    // -- change things for the private/official validation
	    if (string::npos == sample1.find("Off") && string::npos != sample2.find("Off")) {
	      h1string += " (private)";
	      sprintf(loption1, "p");
	      h1->SetMarkerStyle(24);
	      h1->SetMarkerSize(1.5);
	      h1->Draw("esame");
	    }
	  }
	} else if (string::npos != sample1.find("Data")) {
	  sprintf(loption1, "p");
	  if (string::npos != sample1.find("mm")) {
	    h1string = "data sidebands";
	  } else {
	    h1string = "data";
	  }
	} else {
	  h1string = "??";
	}

	if (string::npos != sample2.find("Mc")) {
	  sprintf(loption2, "f");
	  if (string::npos != sample2.find("bdmm")) {
	    h2string = "B^{0} #rightarrow #mu^{+} #mu^{-}";
	  } else if (string::npos != sample2.find("bsmm")) {
	    h2string = "B^{0}_{s} #rightarrow #mu^{+} #mu^{-}";
	  } else {
	    h2string = "MC simulation";
	    // -- change things for the private/official validation
	    if (string::npos == sample1.find("Off") && string::npos != sample2.find("Off")) {
	      h2string += " (official)";
	    }
	  }
	} else if (string::npos != sample2.find("Data")) {
	  sprintf(loption2, "p");
	  if (string::npos != sample2.find("bmm")) {
	    h2string = "data sidebands";
	  } else {
	    h2string = "data";
	  }
	} else {
	  h2string = "??";
	}

	legg->SetHeader(header.c_str());
	legg->SetTextSize(0.05);
	legg->AddEntry(h1, h1string.c_str(), loption1);
	legg->AddEntry(h2, h2string.c_str(), loption2);

	legg->Draw();
      } else {
	legg = 0;
      }

      stamp(0.18, fStampCms, fStampString, 0.4, fStampLumi);

      if (1) {
	TLatex ll;
	ll.SetTextAngle(90.);
	ll.SetTextSize(0.03);
	ll.DrawLatexNDC(0.97, 0., Form("%s/%s/%s/%s", sample1.c_str(), sample2.c_str(), selection.c_str(), fDoList[i].c_str()));
      }

      c0->Modified();
      c0->Update();
      c0->SaveAs(Form("%s/overlay%d%s_%s_%s_%s_%s.pdf",
		      fDirectory.c_str(), fYear, fSetup.c_str(), sample1.c_str(), sample2.c_str(), fDoList[i].c_str(), selection.c_str()));



      // -- clean up
      TPad *p = (TPad*)c0->GetPrimitive("pad1");
      delete p;
      p = (TPad*)c0->GetPrimitive("pad2");
      TH1D *h = (TH1D*)(p->GetPrimitive("hr"));
      delete h;
      delete p;
      if (legg) delete legg;
    }
  }

}


// ----------------------------------------------------------------------
void plotReducedOverlays::loopFunction2() {

}


// ----------------------------------------------------------------------
void plotReducedOverlays::fillDistributions(string selmode) {

  string mapname = Form("ad%s%s_%s", fChannel.c_str(), selmode.c_str(), fSample.c_str());
  //  cout << "fillDistributions: mapname = " << mapname << endl;
  double mass = fb.cm;
  if (fIsMC) mass = fb.m;
  if (fIsSignal) mass = fb.m;
  //  mass = fb.m;

  TLorentzVector a;
  a.SetPtEtaPhiM(fb.pt,fb.eta,fb.phi,mass);

  fAdMap[mapname]->fpMuon1Pt->fill(fb.m1pt, mass);
  fAdMap[mapname]->fpMuon2Pt->fill(fb.m2pt, mass);

  fAdMap[mapname]->fpMuonsEta->fill(fb.m1eta, mass);
  fAdMap[mapname]->fpMuonsEta->fill(fb.m2eta, mass);
  fAdMap[mapname]->fpPt->fill(fb.pt, mass);
  fAdMap[mapname]->fpP->fill(a.P(), mass);
  fAdMap[mapname]->fpPz->fill(a.Pz(), mass);
  fAdMap[mapname]->fpEta->fill(fb.eta, mass);
  fAdMap[mapname]->fpAlpha->fill(fb.alpha, mass);

  fAdMap[mapname]->fpIso->fill(fb.iso, mass);
  fAdMap[mapname]->fpCloseTrk->fill(fb.closetrk, mass);
  fAdMap[mapname]->fpDocaTrk->fill(fb.docatrk, mass);

  fAdMap[mapname]->fpChi2Dof->fill(fb.chi2/fb.dof, mass);
  fAdMap[mapname]->fpPChi2Dof->fill(fb.pchi2dof, mass);

  fAdMap[mapname]->fpFLS3d->fill(fb.fls3d, mass);
  fAdMap[mapname]->fpFL3d->fill(fb.fl3d, mass);
  fAdMap[mapname]->fpFL3dE->fill(fb.fl3dE, mass);

  fAdMap[mapname]->fpFLSxy->fill(fb.flsxy, mass);
  fAdMap[mapname]->fpFLxy->fill(fb.flxy, mass);
  fAdMap[mapname]->fpFLxyE->fill(fb.flxy/fb.flsxy, mass);

  fAdMap[mapname]->fpTau->fill(1.e12*fb.tau, mass);

  fAdMap[mapname]->fpMaxDoca->fill(fb.maxdoca, mass);
  fAdMap[mapname]->fpIp->fill(fb.pvip, mass);
  fAdMap[mapname]->fpIpS->fill(fb.pvips, mass);
  fAdMap[mapname]->fpPvZ->fill(fb.pvz, mass);
  fAdMap[mapname]->fpDzmin->fill(fb.dzmin, mass);
  fAdMap[mapname]->fpDz12->fill(fb.dz12, mass);
  fAdMap[mapname]->fpPvN->fill(fb.pvn, mass);
  fAdMap[mapname]->fpPvNtrk->fill(fb.pvntrk, mass);
  fAdMap[mapname]->fpPv2Ntrk->fill(fb.pv2ntrk, mass);
  fAdMap[mapname]->fpPvAveW8->fill(fb.pvw8, mass);

  fAdMap[mapname]->fpM1Iso->fill(fb.m1iso, mass);
  fAdMap[mapname]->fpM2Iso->fill(fb.m2iso, mass);

  fAdMap[mapname]->fpLip->fill(TMath::Abs(fb.pvlip), mass);
  fAdMap[mapname]->fpLipS->fill(TMath::Abs(fb.pvlips), mass);

  fAdMap[mapname]->fpLip2->fill(TMath::Abs(fb.pv2lip), mass);
  fAdMap[mapname]->fpLipS2->fill(TMath::Abs(fb.pv2lips), mass);
  // fAdMap[mapname]->fpLastCut->fill(mass, mass);

  fAdMap[mapname]->fpOtherVtx->fill(fb.othervtx, mass);
  fAdMap[mapname]->fpPvDchi2->fill(fb.pvdchi2, mass);

  fAdMap[mapname]->fpBDT->fill(fBDT, mass);
  fAdMap[mapname]->fpBDTSel0->fill(fBDT, mass);
  fAdMap[mapname]->fpBDTSel1->fill(fBDT, mass);
  fAdMap[mapname]->fpBDTSel2->fill(fBDT, mass);
}

// ----------------------------------------------------------------------
AnalysisDistribution* plotReducedOverlays::bookDistribution(string hn, string ht, string hc, AnalysisCuts *pCuts,
							    int nbins, double lo, double hi, bool *presel) {
  AnalysisDistribution *p = new AnalysisDistribution(hn.c_str(), ht.c_str(), nbins, lo, hi);
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX);
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX);
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX);
  p->setAnalysisCuts(pCuts, hc.c_str());
  if (0 == presel) {
    p->setPreselCut(&fPreselection);
  } else {
    p->setPreselCut(presel);
  }

  return p;
}

// ----------------------------------------------------------------------
void plotReducedOverlays::overlayAndRatio(TCanvas *c, TH1D *h1, TH1D *h2) {
  bool drawGrid(true);
  bool fitRatio(false);
  string sname = string(h1->GetName());
  bool isBmmData = (string::npos != sname.find("bmmData"));

  // -- Upper plot
  double splity(0.3);
  c->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0.0, splity, 1.0, 1.0);
  pad1->SetTopMargin(0.1);
  pad1->SetBottomMargin(0.);
  pad1->SetRightMargin(0.05);
  if (drawGrid) pad1->SetGridx();
  pad1->Draw();
  pad1->cd();

  h1->SetTitle("");
  h1->SetMinimum(0.001);
  double ymax = h1->GetMaximum();
  if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();
  h1->SetStats(0);
  string hname = h1->GetName();
  if (string::npos != hname.find("Data")) {
    h1->Draw("e");
  } else {
    h1->Draw();
  }
  h2->Draw("samehist");
  h1->SetMaximum(1.2*ymax);

  // -- Lower plot
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1.0, splity);
  pad2->SetTopMargin(0.);
  pad2->SetBottomMargin(0.4);
  pad2->SetRightMargin(0.05);
  pad2->Draw();
  if (drawGrid) pad2->SetGridy();
  if (drawGrid) pad2->SetGridx();
  pad2->cd();

  TH1D *hr = (TH1D*)h1->Clone("hr");
  hr->SetLineColor(kBlack);
  hr->SetMinimum(0.4);
  hr->SetMaximum(1.6);
  hr->SetStats(0);
  hr->Divide(h2);
  hr->SetMarkerStyle(24);
  hr->Draw("e0");
  pl->DrawLine(hr->GetBinLowEdge(1), 1., hr->GetBinLowEdge(hr->GetNbinsX()), 1.0);
  if (isBmmData) {
    hr->Reset();
    double ts(tl->GetTextSize());
    tl->SetTextSize(0.1);
    tl->DrawLatexNDC(0.3, 0.8, "ratio not meaningful");
  }


  if (fitRatio) {
    hr->Fit("pol1", "q");
    tl->SetTextAngle(90.);
    double ts(tl->GetTextSize());
    tl->SetTextSize(0.09);
    tl->SetTextColor(kRed);
    tl->DrawLatexNDC(0.91, 0.3, Form("p1 = %+4.3f #pm %4.3f",
				     hr->GetFunction("pol1")->GetParameter(1),
				     hr->GetFunction("pol1")->GetParError(1)));
    tl->SetTextSize(ts);
    tl->SetTextColor(kBlack);
    tl->SetTextAngle(0.);
  }

  double psize = 0.07;
  double pratio = (1-splity)/splity;
  h1->GetYaxis()->SetTitleSize(psize);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitleOffset(1.55);

  h1->GetYaxis()->SetLabelSize(psize);
  h1->GetXaxis()->SetLabelSize(psize);


  // hr settings
  hr->SetTitle("");
  hr->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  hr->GetXaxis()->SetTitleOffset(1.0);
  hr->GetYaxis()->SetTitle("ratio");
  hr->GetYaxis()->SetTitleOffset(0.4);
  hr->GetYaxis()->CenterTitle();

  hr->GetYaxis()->SetNdivisions(204);
  hr->GetYaxis()->SetTitleFont(42);
  hr->GetXaxis()->SetTitleFont(42);
  hr->GetYaxis()->SetTitleSize(pratio*psize);
  hr->GetXaxis()->SetTitleSize(pratio*psize);
  hr->GetYaxis()->SetLabelSize(pratio*psize);
  hr->GetXaxis()->SetLabelSize(pratio*psize);

  pad1->cd();
}


// ----------------------------------------------------------------------
void plotReducedOverlays::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotReducedOverlays::loadFile loading files listed in " << files << endl;

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

    if (string::npos != stype.find("data")) {
      // -- DATA
      pF = loadFile(sfile);

      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("XXX")) {
        sname = "bdpsikstarData";
        sdecay = "bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
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

      if ((string::npos != stype.find("bdmm")) && (string::npos != stype.find("official"))) {
        sname = "bdmmOfficial";
        sdecay = "bdmm";
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
      fDS.insert(make_pair(sname, ds));
    } else {
      delete ds;
    }


  }

  is.close();

  cout << "------------------------------------------------------------------------------------------" << endl;
  cout << Form("%30s: %20s: ", "Dataset name", "decay mode name") << "Filename" << endl;
  cout << "------------------------------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    cout << Form("%30s: %20s: ", it->first.c_str(), it->second->fName.c_str()) << it->second->fF->GetName() << endl;
  }
  cout << "------------------------------------------------------------------------------------------" << endl;
}
