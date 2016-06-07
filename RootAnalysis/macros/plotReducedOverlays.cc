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

  TVirtualFitter::SetMaxIterations(50000);

  cout << "==> plotReducedOverlays files: " << files << " dir: " << dir << " cuts: " << cuts << " setup: " << setup << endl;

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  loadFiles(files);

  if (setup == "") {
    fHistFileName = Form("%s/plotReducedOverlays.root", dir.c_str());
    fNumbersFileName = fDirectory + "/plotReducedOverlays.txt";
  } else {
    fHistFileName = Form("%s/plotReducedOverlays-%s.root", dir.c_str(), setup.c_str());
    fNumbersFileName = fDirectory + "/plotReducedOverlays." + setup + ".txt";
  }

  fTexFileName = fNumbersFileName;
  replaceAll(fTexFileName, ".txt", ".tex");
  fTEX.open(fNumbersFileName.c_str(), ios::app);

  fOffset = 0;

  fIsMC = false;
  fIsSignal = false;

  fSel0 = false;
  fSel1 = false;
  fSel2 = false;

  fDoList.clear();
  fDoList.push_back("muon1pt");
  if (1) {
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
    fDoList.push_back("fls3d");
    fDoList.push_back("fl3d");
    fDoList.push_back("fl3de");

    fDoList.push_back("maxdoca");
    fDoList.push_back("ip");
    fDoList.push_back("ips");
    fDoList.push_back("pvn");
    fDoList.push_back("pvavew8");

    fDoList.push_back("lip");
    fDoList.push_back("lips");

    fDoList.push_back("lip2");
    fDoList.push_back("lips2");

    fDoList.push_back("m1iso");
    fDoList.push_back("m2iso");
    fDoList.push_back("othervtx");
    fDoList.push_back("pvdchi2");
    fDoList.push_back("closetrks1");
    fDoList.push_back("closetrks2");
    fDoList.push_back("closetrks3");
  }

  fChannelList.clear();
  fChannelList.push_back("0");
  // fChannelList.push_back("1");
  // fChannelList.push_back("2");

  // fChannelList.push_back("0lopu");
  // fChannelList.push_back("1lopu");
  // fChannelList.push_back("2lopu");

  // fChannelList.push_back("0hipu");
  // fChannelList.push_back("1hipu");
  // fChannelList.push_back("2hipu");

}


// ----------------------------------------------------------------------
plotReducedOverlays::~plotReducedOverlays() {

}


// ----------------------------------------------------------------------
void plotReducedOverlays::init() {

  system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  system(Form("/bin/rm -f %s/sbsctrl_ad*_*.pdf", fDirectory.c_str()));
  system(Form("/bin/rm -f %s/overlay_ad*_*.pdf", fDirectory.c_str()));
  system(Form("/bin/rm -f %s/mass_ad*_*.pdf", fDirectory.c_str()));

}

// ----------------------------------------------------------------------
void plotReducedOverlays::makeAll(string selection) {

  init();

  makeSampleOverlay("bspsiphiData", "bspsiphiMc", "Ao");
  makeSampleOverlay("bmmData", "bdmmMc", "Presel");
  makeSampleOverlay("bupsikData", "bupsikMc", "Ao");
  makeSampleOverlay("bdpsikstarData", "bdpsikstarMc", "Ao");

  plotMass("bspsiphiData", "Cu");
  plotMass("bupsikData", "Cu");
  plotMass("bdpsikstarData", "Cu");
  plotMass("bmmData", "Presel");

}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeSampleOverlay(string sample1, string sample2, string selection) {

  makeSample(sample1, selection);
  makeSample(sample2, selection);

  makeOverlay(sample1,
	      sample2,
	      selection);


  //  overlay(sample1, sample2, selection);
}


// ----------------------------------------------------------------------
void plotReducedOverlays::plotMass(string sample, string selection) {

  cout << "fHistFileName: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

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

  for (unsigned int i = 0; i < 5; ++i) {
    h = (TH1D*)gDirectory->Get(Form("ad%d_%s_lastcutMass%s", i, sample.c_str(), selection.c_str()));
    if (!h) break;
    TF1 *f1 = fIF->expoErrGauss(h, 5.28, 0.04);
    setTitles(h, "m [GeV]", Form("Entries / %4.3f GeV", h->GetBinWidth(1)), 0.05, 1.1, 1.3);
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
      tl->DrawLatexNDC(0.55, 0.70, Form("Signal: %5.1f  #pm %5.1f", signal, signalE));
      tl->DrawLatexNDC(0.55, 0.66, Form("Mass:   %5.4f  #pm %5.4f GeV", f1->GetParameter(1), f1->GetParError(1)));
      tl->DrawLatexNDC(0.55, 0.62, Form("Width:  %5.4f  #pm %5.4f GeV", f1->GetParameter(2), f1->GetParError(2)));
    } else {
      h->Draw("hist");
    }

    tl->SetTextSize(0.05);
    tl->DrawLatexNDC(0.55, 0.80, header.c_str());
    tl->SetTextSize(0.025);
    tl->DrawLatexNDC(0.55, 0.75, Form("%2.1f < #eta(#mu_{f}) < %2.1f", fCuts[i]->metaMin, fCuts[i]->metaMax));

    stamp(0., fStampCms, fStampString, 0., fStampLumi);

    c0->SaveAs(Form("%s/mass_ad%d_%s_%s.pdf", fDirectory.c_str(), i, sample.c_str(), selection.c_str()));
  }
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotReducedOverlays::sbsSingleFile(string file1, string sample1, string channel, string selection) {

  string hfname  = file1;
  cout << "fHistFile: " << hfname;
  fHistFile = TFile::Open(hfname.c_str());
  cout << " opened " << endl;

  fChannel = channel;

  sbsDistributions(sample1, selection);
  sbsDistributions(sample1, "HLT", "bdtMin");

  TDirectory *pD = gDirectory;

  hfname  = fDirectory + "/anaBmm.plotReducedOverlaysSbs." + fSuffix + ".root";
  TFile *fl = TFile::Open(hfname.c_str(), "UPDATE");
  TH1D *h1(0);

  // -- and now all histograms as well
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    h1 = (TH1D*)pD->Get(Form("sbs_%s_%s_%s%s", fChannel.c_str(), sample1.c_str(), fDoList[i].c_str(), selection.c_str()));
    h1->SetDirectory(fl);
    h1->Write();
  }

  fl->Close();

  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotReducedOverlays::makeOverlay(string sample1, string sample2, string selection) {

  cout << "fHistFileName: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  for (unsigned int i = 0; i < fChannelList.size(); ++i) {
    sbsDistributions(Form("ad%s_%s", fChannelList[i].c_str(), sample1.c_str()), selection);
    sbsDistributions(Form("ad%s_%s", fChannelList[i].c_str(), sample2.c_str()), selection);

    overlay(Form("ad%s_%s", fChannelList[i].c_str(), sample1.c_str()),
	    Form("ad%s_%s", fChannelList[i].c_str(), sample2.c_str()),
	    selection);
  }
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
void plotReducedOverlays::makeSample(string sample, string selection, int nevents, int nstart) {

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
    BGLBOXMIN = 4.80;
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
    BGLBOXMAX = 5.18;
    if (string::npos != fChannel.find("0")) {
      SIGBOXMIN = 5.23;
      SIGBOXMAX = 5.33;
    } else {
      SIGBOXMIN = 5.21;
      SIGBOXMAX = 5.34;
    }
    BGHBOXMIN = 5.40;
    BGHBOXMAX = 5.50;
  }

    if (string::npos != fSample.find("bdpsikstar")) {
    dir = "candAnaBd2JpsiKstar";
    fMode = BD2JPSIKSTAR;
    BGLBOXMIN = 5.00;
    BGLBOXMAX = 5.18;
    if (string::npos != fChannel.find("0")) {
      SIGBOXMIN = 5.23;
      SIGBOXMAX = 5.33;
    } else {
      SIGBOXMIN = 5.21;
      SIGBOXMAX = 5.34;
    }
    BGHBOXMIN = 5.40;
    BGHBOXMAX = 5.50;
  }

  if (string::npos != fSample.find("bspsiphi")) {
    dir = "candAnaBs2JpsiPhi";
    fMode = BS2JPSIPHI;
    BGLBOXMIN = 5.10;
    BGLBOXMAX = 5.29;
    if (string::npos != fChannel.find("0")) {
      SIGBOXMIN = 5.34;
      SIGBOXMAX = 5.40;
    } else {
      SIGBOXMIN = 5.32;
      SIGBOXMAX = 5.41;
    }
    BGHBOXMIN = 5.45;
    BGHBOXMAX = 5.70;
  }

  // -- must be after the mass box definitions!
  bookDistributions();

  TTree *t = getTree(fSample, dir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  }
  setupTree(t, fSample);
  fCds = fSample;
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
  cout << "==> plotReducedOverlays::loopOverTree> loop over dataset " << fCds << " in file "
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

  // -- modify here the fGoodHLT to increase S/B for the BDT distribution
  fGoodHLT        = fb.hlt;


  // -- update ana cuts!
  fAnaCuts.update();

  // bool loPU = (fb.pvn <  6);
  // bool hiPU = (fb.pvn > 15);

  // bool closePV(false);
  // if (1 == fChan) {
  //   closePV = (fb.pvlip2 < 0.15 && fb.pvlip2 > 0.04);
  // } else {
  //   closePV = (fb.pvlip2 < 0.15 && fb.pvlip2 > 0.02);
  // }

  // bool farPV   = (fb.pvlip2 > 0.5);


  if (fb.hlt && fGoodMuonsID && (fBDT > -1.) && fb.fls3d > 10) {
    fSel0 = true;
  } else {
    fSel0 = false;
  }

  if (fb.hlt && fGoodMuonsID && (fBDT > -1.) && fb.fls3d > 15) {
    fSel1 = true;
  } else {
    fSel1 = false;
  }

  if (fb.hlt && fGoodMuonsID && (fBDT > -1.) && fb.fls3d > 20) {
    fSel2 = true;
  } else {
    fSel2 = false;
  }



  if (0 == fChan) {
    fChannel = "0";
    fillDistributions();
    if (0) {
      // if (loPU) {
      // 	fChannel = "0lopu";
      // 	fillDistributions();
      // } else {
      // 	fChannel = "0hipu";
      // 	fillDistributions();
      // }
    }
  }

  if (1 == fChan) {
    fChannel = "1";
    fillDistributions();
    if (0) {
      // if (loPU) {
      // 	fChannel = "1lopu";
      // 	fillDistributions();
      // } else {
      // 	fChannel = "1hipu";
      // 	fillDistributions();
      // }
    }
  }

  if (2 == fChan) {
    fChannel = "2";
    fillDistributions();
    if (0) {
      //   if (loPU) {
      // 	fChannel = "2lopu";
      // 	fillDistributions();
      //   } else {
      // 	fChannel = "2hipu";
      // 	fillDistributions();
      //   }
    }
  }

  //  fillDistributions();

}


// ----------------------------------------------------------------------
void plotReducedOverlays::bookDistributions() {

  cout << "SIG: " << SIGBOXMIN << " .. " << SIGBOXMAX << endl;
  cout << "BGL: " << BGLBOXMIN << " .. " << BGLBOXMAX << endl;
  cout << "BGH: " << BGHBOXMIN << " .. " << BGHBOXMAX << endl;

  adset *a(0);
  for (unsigned int i = 0; i < fChannelList.size(); ++i) {
    string mapname = Form("ad%s_%s", fChannelList[i].c_str(), fSample.c_str());
    string name = Form("%s_", mapname.c_str());

    a = new adset();
    a->fpMuon1Pt  = bookDistribution(Form("%smuon1pt", name.c_str()), "p_{T, #mu1} [GeV]", "fGoodMuonsPt", 60, 0., 30.);

    a->fpMuon2Pt   = bookDistribution(Form("%smuon2pt", name.c_str()), "p_{T, #mu2} [GeV]", "fGoodMuonsPt", 40, 0., 20.);
    a->fpMuonsEta  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsPt", 40, -2.5, 2.5);
    a->fpPt        = bookDistribution(Form("%spt", name.c_str()), "p_{T}(B) [GeV]", "fGoodPt", 60, 0., 60.);
    a->fpP         = bookDistribution(Form("%sp", name.c_str()), "p(B) [GeV]", "fGoodPt", 50, 0., 100.);
    a->fpPz        = bookDistribution(Form("%spz", name.c_str()), "p_{z}(B) [GeV]", "fGoodPt", 50, 0., 100.);
    a->fpEta       = bookDistribution(Form("%seta", name.c_str()), "#eta(B)", "fGoodEta", 40, -2.5, 2.5);
    a->fpAlpha     = bookDistribution(Form("%salpha", name.c_str()), "#alpha_{3D}", "fGoodAlpha", 50, 0., 0.1);
    a->fpIso       = bookDistribution(Form("%siso", name.c_str()),  "isolation", "fGoodIso", 52, 0., 1.04);
    a->fpCloseTrk  = bookDistribution(Form("%sclosetrk", name.c_str()),  "N_{trk}^{close}", "fGoodCloseTrack", 10, 0., 10.);
    a->fpDocaTrk   = bookDistribution(Form("%sdocatrk", name.c_str()), "d_{ca}^{0} [cm]", "fGoodDocaTrk", 50, 0., 0.20);

    a->fpChi2Dof   = bookDistribution(Form("%schi2dof", name.c_str()),  "#chi^{2}/dof", "fGoodChi2", 40, 0., 4.);
    a->fpPChi2Dof  = bookDistribution(Form("%spchi2dof", name.c_str()),  "P(#chi^{2},dof)", "fGoodChi2", 50, 0., 1.0);

    a->fpFLS3d     = bookDistribution(Form("%sfls3d", name.c_str()), "l_{3D}/#sigma(l_{3D})", "fGoodFLS", 60, 0., 120.);
    a->fpFL3d      = bookDistribution(Form("%sfl3d", name.c_str()),  "l_{3D} [cm]", "fGoodFLS", 60, 0., 1.5);
    a->fpFL3dE     = bookDistribution(Form("%sfl3de", name.c_str()), "#sigma(l_{3D}) [cm]", "fGoodFLS", 50, 0., 0.05);

    a->fpMaxDoca   = bookDistribution(Form("%smaxdoca", name.c_str()), "d^{max} [cm]", "fGoodMaxDoca", 60, 0., 0.03);
    a->fpIp        = bookDistribution(Form("%sip", name.c_str()), "#delta_{3D} [cm]", "fGoodIp", 50, 0., 0.015);
    a->fpIpS       = bookDistribution(Form("%sips", name.c_str()), "#delta_{3D}/#sigma(#delta_{3D})", "fGoodIpS", 50, 0., 4);
    // a->fpPvZ      = bookDistribution(Form("%spvz", name.c_str()), "z_{PV} [cm]", "fGoodHLT", 40, -20., 20.);
    a->fpPvN       = bookDistribution(Form("%spvn", name.c_str()), "N(PV) ", "fGoodHLT", 40, 0., 40.);
    a->fpPvAveW8   = bookDistribution(Form("%spvavew8", name.c_str()), "<w^{PV}>", "fGoodPvAveW8", 50, 0.5, 1.);

    a->fpBDT       = bookDistribution(Form("%sbdt", name.c_str()), "BDT", "fGoodHLT", 200, -1.0, 1.0);

    a->fpCloseTrkS1= bookDistribution(Form("%sclosetrks1", name.c_str()),  "N_{trk}^{close, 1#sigma}", "fGoodCloseTrack", 10, 0., 10.);
    a->fpCloseTrkS2= bookDistribution(Form("%sclosetrks2", name.c_str()),  "N_{trk}^{close, 2#sigma}", "fGoodCloseTrack", 10, 0., 10.);
    a->fpCloseTrkS3= bookDistribution(Form("%sclosetrks3", name.c_str()),  "N_{trk}^{close, 3#sigma}", "fGoodCloseTrack", 10, 0., 10.);
    a->fpM1Iso     = bookDistribution(Form("%sm1iso", name.c_str()),  "m1 isolation", "fGoodIso", 52, 0., 1.04);
    a->fpM2Iso     = bookDistribution(Form("%sm2iso", name.c_str()),  "m2 isolation", "fGoodIso", 52, 0., 1.04);

    a->fpPvDchi2   = bookDistribution(Form("%spvdchi2", name.c_str()),  "#Delta(#chi^{2})", "fGoodChi2", 100, 0., 2000.);
    a->fpOtherVtx  = bookDistribution(Form("%sothervtx", name.c_str()),  "othervtx", "fGoodChi2", 40, 0., 1.);

    a->fpLip       = bookDistribution(Form("%slip", name.c_str()), "l_{z} [cm]", "fGoodLip", 50, 0., 0.015);
    a->fpLipS      = bookDistribution(Form("%slips", name.c_str()), "l_{z}/#sigma(l_{z})", "fGoodLipS", 50, 0., 4);

    a->fpLip2      = bookDistribution(Form("%slip2", name.c_str()), "l_{z}^{(2)} [cm]", "fGoodLip", 100, 0., 4.0);
    a->fpLipS2     = bookDistribution(Form("%slips2", name.c_str()), "l_{z}^{(2)}/#sigma(l_{z}^{(2)})", "fGoodLipS", 50, 0., 20.);

    a->fpLastCut   = bookDistribution(Form("%slastcut", name.c_str()), "lastcut", "fGoodLastCut", 50, 4., 6.);

    a->fpBDTSel0   = bookSpecialDistribution(Form("%sbdtsel0", name.c_str()), "BDTsel0", "fGoodHLT", 200, -1.0, 1.0, &fSel0);
    a->fpBDTSel1   = bookSpecialDistribution(Form("%sbdtsel1", name.c_str()), "BDTsel1", "fGoodHLT", 200, -1.0, 1.0, &fSel1);
    a->fpBDTSel2   = bookSpecialDistribution(Form("%sbdtsel2", name.c_str()), "BDTsel2", "fGoodHLT", 200, -1.0, 1.0, &fSel2);

    fAdMap.insert(make_pair(mapname, a));
    cout << "bookDistributions: mapname = " << mapname << endl;

  }

}


// ----------------------------------------------------------------------
void plotReducedOverlays::bookDistributions0() {

  string name = Form("ad%s_%s_", fChannel.c_str(), fSample.c_str());

  cout << "fOffset: " << fOffset << " name = " << name << endl;
  fpMuon1Pt[fOffset]   = bookDistribution(Form("%smuon1pt", name.c_str()), "p_{T, #mu1} [GeV]", "fGoodMuonsPt", 60, 0., 30.);
  fpMuon2Pt[fOffset]   = bookDistribution(Form("%smuon2pt", name.c_str()), "p_{T, #mu2} [GeV]", "fGoodMuonsPt", 40, 0., 20.);
  fpMuonsEta[fOffset]  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsPt", 40, -2.5, 2.5);
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
  fpIp[fOffset]        = bookDistribution(Form("%sip", name.c_str()), "#delta_{3D} [cm]", "fGoodIp", 50, 0., 0.015);
  fpIpS[fOffset]       = bookDistribution(Form("%sips", name.c_str()), "#delta_{3D}/#sigma(#delta_{3D})", "fGoodIpS", 50, 0., 4);
  //  fpPvZ[fOffset]       = bookDistribution(Form("%spvz", name.c_str()), "z_{PV} [cm]", "fGoodHLT", 40, -20., 20.);
  fpPvN[fOffset]       = bookDistribution(Form("%spvn", name.c_str()), "N(PV) ", "fGoodHLT", 40, 0., 40.);
  fpPvAveW8[fOffset]   = bookDistribution(Form("%spvavew8", name.c_str()), "<w^{PV}>", "fGoodPvAveW8", 50, 0.5, 1.);

  fpBDT[fOffset]       = bookDistribution(Form("%sbdt", name.c_str()), "BDT", "fGoodHLT", 200, -1.0, 1.0);

  fpCloseTrkS1[fOffset]= bookDistribution(Form("%sclosetrks1", name.c_str()),  "N_{trk}^{close, 1#sigma}", "fGoodCloseTrack", 10, 0., 10.);
  fpCloseTrkS2[fOffset]= bookDistribution(Form("%sclosetrks2", name.c_str()),  "N_{trk}^{close, 2#sigma}", "fGoodCloseTrack", 10, 0., 10.);
  fpCloseTrkS3[fOffset]= bookDistribution(Form("%sclosetrks3", name.c_str()),  "N_{trk}^{close, 3#sigma}", "fGoodCloseTrack", 10, 0., 10.);
  fpM1Iso[fOffset]     = bookDistribution(Form("%sm1iso", name.c_str()),  "m1 isolation", "fGoodIso", 52, 0., 1.04);
  fpM2Iso[fOffset]     = bookDistribution(Form("%sm2iso", name.c_str()),  "m2 isolation", "fGoodIso", 52, 0., 1.04);

  fpPvDchi2[fOffset]   = bookDistribution(Form("%spvdchi2", name.c_str()),  "#Delta(#chi^{2})", "fGoodChi2", 100, 0., 2000.);
  fpOtherVtx[fOffset]  = bookDistribution(Form("%sothervtx", name.c_str()),  "othervtx", "fGoodChi2", 40, 0., 1.);

  fpLip[fOffset]       = bookDistribution(Form("%slip", name.c_str()), "l_{z} [cm]", "fGoodLip", 50, 0., 0.015);
  fpLipS[fOffset]      = bookDistribution(Form("%slips", name.c_str()), "l_{z}/#sigma(l_{z})", "fGoodLipS", 50, 0., 4);

  fpLip2[fOffset]      = bookDistribution(Form("%slip2", name.c_str()), "l_{z}^{(2)} [cm]", "fGoodLip", 100, 0., 4.0);
  fpLipS2[fOffset]     = bookDistribution(Form("%slips2", name.c_str()), "l_{z}^{(2)}/#sigma(l_{z}^{(2)})", "fGoodLipS", 50, 0., 20.);


  fpBDTSel0[fOffset]   = bookSpecialDistribution(Form("%sbdtsel0", name.c_str()), "BDTsel0", "fGoodHLT", 200, -1.0, 1.0, &fSel0);
  fpBDTSel1[fOffset]   = bookSpecialDistribution(Form("%sbdtsel1", name.c_str()), "BDTsel1", "fGoodHLT", 200, -1.0, 1.0, &fSel1);
  fpBDTSel2[fOffset]   = bookSpecialDistribution(Form("%sbdtsel2", name.c_str()), "BDTsel2", "fGoodHLT", 200, -1.0, 1.0, &fSel2);

}


// ----------------------------------------------------------------------
void plotReducedOverlays::sbsDistributions(string sample, string selection, string what) {

  string sbsControlPlotsFileName = Form("sbsctrl");

  AnalysisDistribution a(Form("%s_muon1pt", sample.c_str()));
  a.fVerbose = 1;
  a.fControlPlotsFileName = sbsControlPlotsFileName;
  a.fDirectory = fDirectory;

  int type(0);
  if (string::npos != sample.find("bmmData"))      type = 0; // sidebands
  if (string::npos != sample.find("Mc"))           type = 1; // signal window
  if (string::npos != sample.find("bupsikData"))   type = 2; // pol1+err
  if (string::npos != sample.find("bspsiphiData")) type = 3; // expo
  double preco(5.15);
  if (string::npos != sample.find("bspsiphiData")) {
    preco = 5.2;
    if (string::npos != fChannel.find("0")) {
      a.fMassPeak = 5.37;
      a.fMassSigma = 0.06;
    } else {
      a.fMassPeak = 5.37;
      a.fMassSigma = 0.06;
    }
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
      cout << "=> Looking for sideband histogram " << Form("%s%s1", bla.c_str(), selection.c_str()) << endl;
      h = (TH1D*)gDirectory->Get(Form("%s%s1", bla.c_str(), selection.c_str()));
      cout << "=> cloning into  "
	   << Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), selection.c_str()) << endl;
      h = (TH1D*)h->Clone(Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), selection.c_str()));
    } else if (1 == type) {
      cout << "=> Looking for signal histogram " << Form("%s%s0", bla.c_str(), selection.c_str()) << endl;
      h = (TH1D*)gDirectory->Get(Form("%s%s0", bla.c_str(), selection.c_str()));
      cout << "=> cloning into  "
	   << Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), selection.c_str()) << endl;
      h = (TH1D*)h->Clone(Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), selection.c_str()));
    } else if (2 == type) {
      cout << "=> sbsDistributionPol1ErrGauss histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      h = a.sbsDistributionPol1ErrGauss(bla.c_str(), selection.c_str(), preco);
    } else if (3 == type) {
      cout << "=> sbsDistributionExpoGauss histogram " << Form("%s for selection %s", bla.c_str(), selection.c_str()) << endl;
      h = a.sbsDistributionExpoGauss(bla.c_str(), selection.c_str());
    }
    cout << "  Title: " << h->GetTitle() << " with integral: " << h->GetSumOfWeights() << endl;
  }

}


// ----------------------------------------------------------------------
void plotReducedOverlays::allSystematics() {
  for (int i = 0; i < 2; ++i) {
    systematics("CsData", "CsMc", i);
    systematics("NoData", "NoMc", i);
    systematics("CsData", "SgMc", i);
    systematics("SgMc",   "CsMc", i);
  }

}


// ----------------------------------------------------------------------
void plotReducedOverlays::systematics(string sample1, string sample2, int chan) {
  gStyle->SetOptTitle(0);
  zone(1);
  c0->cd();

  string sChan = (chan == 0? "B": "E");

  double bdtCut = fCuts[chan]->bdtMin;
  cout << "bdtCut = " << bdtCut << endl;

  string hfname  = fDirectory + "/anaBmm.plotReducedOverlaysSystematics." + fSuffix + ".root";
  cout << "fHistFile: " << hfname << endl;
  fHistFile = TFile::Open(hfname.c_str(), "");

  // -- extract here the means of the NPV distributions
  TH1D *hpv1 = (TH1D*)fHistFile->Get(Form("sbs_%s_%s_pvnPresel", sChan.c_str(), sample1.c_str()));
  TH1D *hpv2 = (TH1D*)fHistFile->Get(Form("sbs_%s_%s_pvnPresel", sChan.c_str(), sample2.c_str()));

  if (string::npos != sample2.find("CsMc")) {
    setHist(hpv2, kRed, 20, 1.5);
    setFilledHist(hpv2, kRed, kRed, 3365);
  } else if (string::npos != sample2.find("SgMc")) {
    setHist(hpv2, kRed, 20, 1.5);
    setFilledHist(hpv2, kRed, kRed, 3344);
  } else {
    setHist(hpv2, kBlue, 20, 1.5);
    setFilledHist(hpv2, kBlue, kBlue, 3365);
  }

  shrinkPad(0.15, 0.18);
  hpv1->SetMinimum(0.);
  hpv1->Draw();
  hpv2->Scale(hpv1->GetSumOfWeights()/hpv2->GetSumOfWeights());
  hpv2->Draw("samehist");
  tl->DrawLatex(0.15, 0.92, Form("Data: %3.2f#pm%3.2f", hpv1->GetMean(), hpv1->GetMeanError()));
  tl->DrawLatex(0.6, 0.92, Form("MC: %3.2f#pm%3.2f", hpv2->GetMean(), hpv2->GetMeanError()));

  c0->SaveAs(Form("%s/%s-systematics-npv_%s_%s_chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan));


  TH1D *h1 = (TH1D*)fHistFile->Get(Form("sbs_%s_%s_bdtHLT", sChan.c_str(), sample1.c_str()));
  for (int i = 1; i <= h1->GetNbinsX(); ++i) if (h1->GetBinContent(i) < 0) h1->SetBinContent(i, -0.0001);
  HistCutEfficiency a1(h1, bdtCut, 0);
  double eps1 = a1.hiEff;
  double eps1E = a1.hiErr;
  cout << "eps1 = " << eps1 << " lo eff = " << a1.loEff << endl;

  TH1D *h2 = (TH1D*)fHistFile->Get(Form("sbs_%s_%s_bdtHLT", sChan.c_str(), sample2.c_str()));
  for (int i = 1; i <= h2->GetNbinsX(); ++i) if (h2->GetBinContent(i) < 0) h2->SetBinContent(i, -0.0001);
  HistCutEfficiency a2(h2, bdtCut, 0);
  double eps2 = a2.hiEff;
  double eps2E = a2.hiErr;

  cout << "eps2 = " << eps2 << " lo eff = " << a2.loEff << endl;
  double deltaEps = eps1-eps2;
  double adeltaEps = TMath::Abs(eps1-eps2);
  double rdeltaEps = 2.*TMath::Abs(eps1-eps2)/(eps1+eps2);

  fTEX << formatTex(fCuts[chan]->bdtMin, Form("%s:sysCutOnBdtchan%i:val", fSuffix.c_str(), chan), 3) << endl;
  fTEX << formatTex(deltaEps, Form("%s:deltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;
  fTEX << formatTex(adeltaEps, Form("%s:absDeltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;
  fTEX << formatTex(rdeltaEps, Form("%s:relDeltaEps%s%schan%i:val", fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan), 3) << endl;

  if (string::npos != sample2.find("CsMc")) {
    setHist(h2, kRed, 20, 1.5);
    setFilledHist(h2, kRed, kRed, 3365);
  } else if (string::npos != sample2.find("SgMc")) {
    setHist(h2, kRed, 20, 1.5);
    setFilledHist(h2, kRed, kRed, 3344);
  } else {
    setHist(h2, kBlue, 20, 1.5);
    setFilledHist(h2, kBlue, kBlue, 3365);
  }

  h1->SetMinimum(0.);
  h1->SetMaximum(1.3*h1->GetMaximum());
  h1->Draw("e");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  h2->Draw("histsame");

  newLegend(0.2, 0.77, 0.45, 0.87);
  legg->SetTextSize(0.035);
  string text1, text2;
  if (string::npos != sample1.find("CsData")) text1 = "B_{s}^{0} #rightarrow J/#psi #phi (data)";
  if (string::npos != sample1.find("NoData")) text1 = "B^{+} #rightarrow J/#psi K (data)";
  if (string::npos != sample1.find("SgMc"))   text1 = "B_{s}^{0} #rightarrow #mu #mu (MC)";

  if (string::npos != sample2.find("SgMc")) text2 = "B_{s}^{0} #rightarrow #mu #mu (MC)";
  if (string::npos != sample2.find("CsMc")) text2 = "B_{s}^{0} #rightarrow J/#psi #phi (MC)";
  if (string::npos != sample2.find("NoMc")) text2 = "B^{+} #rightarrow J/#psi K (MC)";


  legg->AddEntry(h1, Form("#varepsilon = %4.3f#pm%4.3f, %s", eps1, eps1E, text1.c_str()), "p");
  legg->AddEntry(h2, Form("#varepsilon = %4.3f#pm%4.3f, %s", eps2, eps2E, text2.c_str()), "f");
  legg->Draw();

  tl->SetTextSize(0.035);
  tl->DrawLatex(0.21, 0.70, Form("rel difference: %4.3f", rdeltaEps));
  tl->DrawLatex(0.21, 0.65, Form("b> %4.3f", bdtCut));

  double yhi = 0.3*h1->GetMaximum();
  double ylo = 0.;
  pa->DrawArrow(bdtCut, yhi, bdtCut, ylo);

  c0->SaveAs(Form("%s/%s-systematics_%s_%s_chan%d.pdf", fDirectory.c_str(), fSuffix.c_str(), sample1.c_str(), sample2.c_str(), chan));

}




// ----------------------------------------------------------------------
void plotReducedOverlays::overlay(string sample1, string sample2, string selection, string what) {

  gStyle->SetOptTitle(0);
  c0->cd();
  shrinkPad(0.15, 0.18);

  string ds1 = sample1.substr(sample1.find("_")+1);
  string ds2 = sample2.substr(sample2.find("_")+1);

  TH1D *h1(0), *h2(0);
  string n1, n2;
  bool restricted = (what != "");
  bool doLegend(true);
  bool leftLegend(false);
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    n1 =  Form("sbs_%s_%s%s", sample1.c_str(), fDoList[i].c_str(), selection.c_str());
    n2 =  Form("sbs_%s_%s%s", sample2.c_str(), fDoList[i].c_str(), selection.c_str());
    if (string::npos != fDoList[i].find("eta")) doLegend = false; else doLegend = true;
    if (string::npos != fDoList[i].find("bdt")) leftLegend = true; else leftLegend = false;
    h1 = (TH1D*)gDirectory->Get(n1.c_str());
    cout << "n1: " << n1 << " -> " << h1 << endl;
    h2 = (TH1D*)gDirectory->Get(n2.c_str());
    cout << "n2: " << n2 << " -> " << h2 << endl;
    if (0 == h1 || 0 == h2) {
      cout << "  histograms not found" << endl;
      continue;
    }
    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    cout << "setHist for " << ds1 << " and " << ds2 << endl;

    overlayAndRatio(c0, h1, h2);
    setHist(h1, fDS[ds1]);
    setHist(h2, fDS[ds2]);
    setHist(h2, fDS[ds2]);

    if (doLegend) {
      if (leftLegend) {
	newLegend(0.21, 0.7, 0.41, 0.85);
      } else {
	newLegend(0.50, 0.7, 0.75, 0.85);
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
      legg->AddEntry(h1, h1string.c_str(), loption1);
      legg->AddEntry(h2, h2string.c_str(), loption2);

      legg->Draw();
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
    c0->SaveAs(Form("%s/overlay_%s_%s_%s_%s.pdf",
		    fDirectory.c_str(), sample1.c_str(), sample2.c_str(), fDoList[i].c_str(), selection.c_str()));
  }


}

// ----------------------------------------------------------------------
void plotReducedOverlays::overlayOld(string sample1, string sample2, string selection, string what) {

  gStyle->SetOptTitle(0);
  c0->cd();
  shrinkPad(0.15, 0.18);

  TH1D *h1(0), *h2(0);
  string n1, n2;
  bool restricted = (what != "");
  bool doLegend(true);
  bool leftLegend(false);
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    n1 =  Form("sbs_%s_%s%s", sample1.c_str(), fDoList[i].c_str(), selection.c_str());
    n2 =  Form("sbs_%s_%s%s", sample2.c_str(), fDoList[i].c_str(), selection.c_str());
    if (string::npos != fDoList[i].find("eta")) doLegend = false; else doLegend = true;
    if (string::npos != fDoList[i].find("bdt")) leftLegend = true; else leftLegend = false;
    h1 = (TH1D*)gDirectory->Get(n1.c_str());
    cout << "n1: " << n1 << " -> " << h1 << endl;
    h2 = (TH1D*)gDirectory->Get(n2.c_str());
    cout << "n2: " << n2 << " -> " << h2 << endl;
    if (0 == h1 || 0 == h2) {
      cout << "  histograms not found" << endl;
      continue;
    }
    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    h1->Draw();
    h2->Draw("samehist");
    setHist(h1, kBlack, 20, 1.5);
    if (string::npos != sample2.find("CsMc")) {
      setHist(h2, kRed, 20, 1.5);
      setFilledHist(h2, kRed, kRed, 3365);
    } else {
      setHist(h2, kBlue, 20, 1.5);
      setFilledHist(h2, kBlue, kBlue, 3365);
    }

    double ymax = (h1->GetMaximum() > h2->GetMaximum()? 1.2*h1->GetMaximum() : 1.2*h2->GetMaximum());
    h1->SetMinimum(0);

    //    h1->SetNdivisions(504, "Y");
    h1->SetTitleOffset(1.0, "Y");
    h1->SetTitleSize(0.06, "Y");
    h1->SetLabelSize(0.055, "Y");

    h1->SetNdivisions(504, "X");
    h1->SetTitleOffset(1.0, "X");
    h1->SetTitleSize(0.06, "X");
    h1->SetLabelSize(0.055, "X");
    h1->SetMinimum(0.01);
    h1->SetMaximum(ymax);
    h1->Draw("e");
    h2->Draw("samehist");

    if (doLegend) {
      if (leftLegend) {
	newLegend(0.21, 0.7, 0.41, 0.85);
      } else {
	newLegend(0.50, 0.7, 0.75, 0.85);
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
      legg->AddEntry(h1, h1string.c_str(), loption1);
      legg->AddEntry(h2, h2string.c_str(), loption2);

      legg->Draw();
    }

    stamp(0.18, fStampCms, fStampString, 0.4, fStampLumi);

    if (1) {
      TLatex ll;
      ll.SetTextAngle(90.);
      ll.SetTextSize(0.03);
      ll.DrawLatexNDC(0.93, 0.17, Form("%s/%s/%s/%s", sample1.c_str(), sample2.c_str(), selection.c_str(), fDoList[i].c_str()));
    }

    if (string::npos != fDoList[i].find("npv")) {
      tl->DrawLatex(0.2, 0.92, Form("means: MC(%4.3f) Data(%4.3f)", h1->GetMean(), h2->GetMean()));
    }

    c0->Modified();
    c0->Update();
    c0->SaveAs(Form("%s/overlay_%s_%s_%s_%s.pdf",
		    fDirectory.c_str(), sample1.c_str(), sample2.c_str(), fDoList[i].c_str(), selection.c_str()));
  }


}


// ----------------------------------------------------------------------
void plotReducedOverlays::overlay2Files(std::string file1, std::string sample1,
					std::string file2, std::string sample2,
					std::string chan1, std::string chan2,
					std::string selection, std::string what) {

  TFile *f1 =  TFile::Open(file1.c_str(), "");
  string fn1 = file1;
  rmPath(fn1);
  rmSubString(fn1, ".root");
  cout << "fn1: " << fn1 << endl;
  TFile *f2 =  TFile::Open(file2.c_str(), "");
  string fn2 = file2;
  rmPath(fn2);
  rmSubString(fn2, ".root");
  cout << "fn2: " << fn2 << endl;

  gStyle->SetOptTitle(0);
  c0->cd();
  TH1D *h1(0), *h2(0);
  string n1, n2;
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    n1 =  Form("sbs_%s_%s_%s%s", chan1.c_str(), sample1.c_str(), fDoList[i].c_str(), selection.c_str());
    n2 =  Form("sbs_%s_%s_%s%s", chan2.c_str(), sample2.c_str(), fDoList[i].c_str(), selection.c_str());
    h1 = (TH1D*)f1->Get(n1.c_str());
    cout << "n1: " << n1 << " -> " << h1 << endl;
    h2 = (TH1D*)f2->Get(n2.c_str());
    cout << "n2: " << n2 << " -> " << h2 << endl;
    if (0 == h1 || 0 == h2) {
      cout << "  histograms not found" << endl;
      continue;
    }
    if (h2->GetSumOfWeights() > 0) h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());

    h1->Draw();
    h2->Draw("samehist");
    setHist(h1, kBlack, 20, 1.5);
    if (string::npos != sample1.find("Mc")) {
      setFilledHist(h1, kBlack, kBlack, 3356);
    }

    if (string::npos != sample2.find("CsMc")) {
      setHist(h2, kRed, 20, 1.5);
      setFilledHist(h2, kRed, kRed, 3365);
    } else {
      setHist(h2, kBlue, 20, 1.5);
      setFilledHist(h2, kBlue, kBlue, 3365);
    }

    double ymax = (h1->GetMaximum() > h2->GetMaximum()? 1.2*h1->GetMaximum() : 1.2*h2->GetMaximum());
    h1->SetMinimum(0.01);
    h1->SetMaximum(ymax);
    if (string::npos != sample1.find("Mc")) {
      h1->Draw("hist");
    } else {
      h1->Draw("e");
    }
    h2->Draw("samehist");
    //    double ks = h1->KolmogorovTest(h2);
    //    tl->DrawLatex(0.2, 0.85, Form("P(KS)= %4.3f", ks));

    newLegend(0.0, 0.91, 0.75, 0.98);
    legg->SetTextSize(0.022);
    string text1, text2;

    if (!what.compare("multichan")) {
      if (string::npos != sample1.find("Mc")) {
	//	legg->AddEntry(h1, Form("%s:%s/%s", file1.c_str(), sample1.c_str(), chan1.c_str()), "f");
	legg->AddEntry(h1, Form("%s:%s", sample1.c_str(), chan1.c_str()), "f");
      } else {
	//	legg->AddEntry(h1, Form("%s:%s/%s", file1.c_str(), sample1.c_str(), chan1.c_str()), "p");
	legg->AddEntry(h1, Form("%s:%s", sample1.c_str(), chan1.c_str()), "p");
      }
      if (string::npos != sample2.find("Mc")) {
	//	legg->AddEntry(h2, Form("%s:%s/%s", file2.c_str(), sample2.c_str(), chan2.c_str()), "f");
	legg->AddEntry(h2, Form("%s:%s", sample2.c_str(), chan2.c_str()), "f");
      } else {
	//	legg->AddEntry(h2, Form("%s:%s/%s", file2.c_str(), sample2.c_str(), chan2.c_str()), "p");
	legg->AddEntry(h2, Form("%s:%s", sample2.c_str(), chan2.c_str()), "p");
      }
    } else {
      if (string::npos != sample1.find("Mc")) {
	//	legg->AddEntry(h1, Form("%s:%s", file1.c_str(), fn1.c_str()), "f");
	legg->AddEntry(h1, Form("%s:%s", sample1.c_str(), chan1.c_str()), "f");
      } else {
	//	legg->AddEntry(h1, Form("%s:%s", file1.c_str(), fn1.c_str()), "p");
	legg->AddEntry(h1, Form("%s:%s", sample1.c_str(), chan1.c_str()), "p");
      }
      //      legg->AddEntry(h2, Form("%s:%s", file2.c_str(), fn2.c_str()), "f");
      legg->AddEntry(h2, Form("%s:%s", sample2.c_str(), chan2.c_str()), "f");
    }
    legg->Draw();


    c0->Modified();
    c0->Update();
    c0->SaveAs(Form("%s/%s-overlay2files-%s-%s-%s_%s-%s_%s-%s.pdf",
		    fDirectory.c_str(), fSuffix.c_str(), what.c_str(), sample1.c_str(), chan1.c_str(), sample2.c_str(), chan2.c_str(),
		    fDoList[i].c_str(), selection.c_str()));
  }



}


// ----------------------------------------------------------------------
void plotReducedOverlays::loopFunction2() {

}


// ----------------------------------------------------------------------
void plotReducedOverlays::fillDistributions() {

  string mapname = Form("ad%s_%s", fChannel.c_str(), fSample.c_str());
  //  cout << "fillDistributions: mapname = " << mapname << endl;
  double mass = fb.cm;
  if (fIsMC) mass = fb.m;
  if (fIsSignal) mass = fb.m;
  mass = fb.m;

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

  fAdMap[mapname]->fpMaxDoca->fill(fb.maxdoca, mass);
  fAdMap[mapname]->fpIp->fill(fb.pvip, mass);
  fAdMap[mapname]->fpIpS->fill(fb.pvips, mass);
  //  fpPvZ->fill(fb.pvz, mass);
  fAdMap[mapname]->fpPvN->fill(fb.pvn, mass);
  fAdMap[mapname]->fpPvAveW8->fill(fb.pvw8, mass);

  fAdMap[mapname]->fpCloseTrkS1->fill(fb.closetrks1, mass);
  fAdMap[mapname]->fpCloseTrkS2->fill(fb.closetrks2, mass);
  fAdMap[mapname]->fpCloseTrkS3->fill(fb.closetrks3, mass);

  fAdMap[mapname]->fpM1Iso->fill(fb.m1iso, mass);
  fAdMap[mapname]->fpM2Iso->fill(fb.m2iso, mass);

  fAdMap[mapname]->fpLip->fill(fb.pvlip, mass);
  fAdMap[mapname]->fpLipS->fill(fb.pvlips, mass);

  fAdMap[mapname]->fpLip2->fill(fb.pvlip2, mass);
  fAdMap[mapname]->fpLipS2->fill(fb.pvlips2, mass);
  fAdMap[mapname]->fpLastCut->fill(mass, mass);

  fAdMap[mapname]->fpOtherVtx->fill(fb.othervtx, mass);
  fAdMap[mapname]->fpPvDchi2->fill(fb.pvdchi2, mass);

  fAdMap[mapname]->fpBDT->fill(fBDT, mass);

  fAdMap[mapname]->fpBDTSel0->fill(fBDT, mass);
  fAdMap[mapname]->fpBDTSel1->fill(fBDT, mass);
  fAdMap[mapname]->fpBDTSel2->fill(fBDT, mass);
}


// ----------------------------------------------------------------------
void plotReducedOverlays::fillDistributions0() {

  //  cout << "BDT: " << fBDT << " mass =  " << fb.m << " pt = " << fb.pt << endl;

  double mass = fb.cm;
  if (fIsMC) mass = fb.m;
  if (fIsSignal) mass = fb.m;
  // FIXME remove this once the constrained mass has propagated back into the reduced tree!
  mass = fb.m;

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
  //  fpPvZ[fOffset]->fill(fb.pvz, mass);
  fpPvN[fOffset]->fill(fb.pvn, mass);
  fpPvAveW8[fOffset]->fill(fb.pvw8, mass);

  fpCloseTrkS1[fOffset]->fill(fb.closetrks1, mass);
  fpCloseTrkS2[fOffset]->fill(fb.closetrks2, mass);
  fpCloseTrkS3[fOffset]->fill(fb.closetrks3, mass);

  fpM1Iso[fOffset]->fill(fb.m1iso, mass);
  fpM2Iso[fOffset]->fill(fb.m2iso, mass);

  fpLip[fOffset]->fill(fb.pvlip, mass);
  fpLipS[fOffset]->fill(fb.pvlips, mass);

  fpLip2[fOffset]->fill(fb.pvlip2, mass);
  fpLipS2[fOffset]->fill(fb.pvlips2, mass);

  fpOtherVtx[fOffset]->fill(fb.othervtx, mass);
  fpPvDchi2[fOffset]->fill(fb.pvdchi2, mass);

  fpBDT[fOffset]->fill(fBDT, mass);

  fpBDTSel0[fOffset]->fill(fBDT, mass);
  fpBDTSel1[fOffset]->fill(fBDT, mass);
  fpBDTSel2[fOffset]->fill(fBDT, mass);
}




// ----------------------------------------------------------------------
AnalysisDistribution* plotReducedOverlays::bookDistribution(string hn, string ht, string hc, int nbins, double lo, double hi) {
  AnalysisDistribution *p = new AnalysisDistribution(hn.c_str(), ht.c_str(), nbins, lo, hi);
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX);
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX);
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX);
  p->setAnalysisCuts(&fAnaCuts, hc.c_str());
  p->setPreselCut(&fPreselection);

  return p;
}


// ----------------------------------------------------------------------
AnalysisDistribution* plotReducedOverlays::bookSpecialDistribution(string hn, string ht, string hc, int nbins, double lo, double hi, bool *presel) {
  AnalysisDistribution *p = new AnalysisDistribution(hn.c_str(), ht.c_str(), nbins, lo, hi);
  p->setSigWindow(SIGBOXMIN, SIGBOXMAX);
  p->setBg1Window(BGLBOXMIN, BGLBOXMAX);
  p->setBg2Window(BGHBOXMIN, BGHBOXMAX);
  p->setAnalysisCuts(&fAnaCuts, hc.c_str());
  p->setPreselCut(presel);

  return p;
}

// ----------------------------------------------------------------------
void plotReducedOverlays::overlayAndRatio(TCanvas *c, TH1D *h1, TH1D *h2) {
  bool drawGrid(true);
  bool fitRatio(false);

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
  h1->SetMinimum(0.01);
  double ymax = h1->GetMaximum();
  if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();
  h1->SetStats(0);
  h1->Draw();
  h2->Draw("same");
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
  hr->Sumw2();
  hr->SetStats(0);
  hr->Divide(h2);
  hr->SetMarkerStyle(24);
  //  hr->Draw("e0"); // this will draw error bars if the marker is out of range
  hr->Draw("e0");
  pl->DrawLine(hr->GetBinLowEdge(1), 1., hr->GetBinLowEdge(hr->GetNbinsX()), 1.0);


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

    TFile *pF(0);
    dataset *ds(0);

    if (string::npos != stype.find("data")) {
      // -- DATA
      pF = loadFile(sfile);

      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("bmm")) {
        sname = "bmmData";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }

      if (string::npos != stype.find("bupsik")) {
        sname = "bupsikData";
        sdecay = "bupsik";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }

      if (string::npos != stype.find("bspsiphi")) {
        sname = "bspsiphiData";
        sdecay = "bspsiphi";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }

      if (string::npos != stype.find("bdpsikstar")) {
        sname = "bdpsikstarData";
        sdecay = "bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }

    } else {
      // -- MC
      pF = loadFile(sfile);

      ds = new dataset();
      ds->fSize = 0.1;
      ds->fWidth = 2.;

      if (string::npos != stype.find("bupsik")) {
        sname = "bupsikMc";
        sdecay = "bupsik";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
	fDS.insert(make_pair(sname, ds));
      }

      if (string::npos != stype.find("bspsiphi")) {
        sname = "bspsiphiMc";
        sdecay = "bspsiphi";
	ds->fColor = kRed;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }

      if (string::npos != stype.find("bsmm")) {
        sname = "bsmmMc";
        sdecay = "bsmm";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }


      if (string::npos != stype.find("bdpsikstar")) {
        sname = "bdpsikstarMc";
        sdecay = "bdpsikstar";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }

      if (string::npos != stype.find("bdmm")) {
        sname = "bdmmMc";
        sdecay = "bdmm";
	ds->fColor = kBlue;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }

    }
    ds->fLcolor = ds->fColor;
    ds->fFcolor = ds->fColor;
    ds->fName   = sdecay;
    ds->fFullName = sname;
    fDS.insert(make_pair(sname, ds));



  }

  is.close();

  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    cout << it->first << ": " << it->second->fName << ", " << it->second->fF->GetName() << endl;
  }
}
