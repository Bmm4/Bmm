#include "plotReducedOverlays.hh"

#include "common/AnalysisDistribution.hh"
#include "common/HistCutEfficiency.hh"
#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVirtualFitter.h"
#include <TTimeStamp.h>

#include "common/util.hh"
#include "common/fitPsYield.hh"

using namespace std;

ClassImp(plotReducedOverlays)

// ----------------------------------------------------------------------
plotReducedOverlays::plotReducedOverlays(string dir, string files, string cuts, string setup, int year) : plotClass(dir, files, cuts, setup, year) {

  changeSetup(dir, "plotReducedOverlays", setup);

  fTEX.open(fTexFileName.c_str(), ios::app);

  TVirtualFitter::SetMaxIterations(50000);

  cout << "==> plotReducedOverlays files: " << files << " dir: " << dir << " cuts: " << cuts << " setup: " << setup << endl;

  // -- initialize cuts
  // string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  // cout << "===> Reading cuts from " << cutfile << endl;
  // readCuts(cutfile);
  // fNchan = fCuts.size();
  fNchan = 2;
  fDoCNC = true;
  //  fDoCNC = false;

  plotClass::loadFiles(files);
  plotReducedOverlays::loadFiles(files);

  //  fIncludeOverflowInLastBin = true;
  fIncludeOverflowInLastBin = false;

  fStampCms = fSuffix;

  fIsMC = false;
  fIsSignal = false;

  fSel0 = false;
  fSel1 = false;
  fSel2 = false;
  fSel3 = false;

  fDoList.clear();
  fDoList.push_back("bdtsel2");
  if (1) {
    fDoList.push_back("fls3d");
    fDoList.push_back("muon1pt");
    fDoList.push_back("muon2pt");

    fDoList.push_back("muonseta");
    fDoList.push_back("muonsphi");

    fDoList.push_back("muon1bdt");
    fDoList.push_back("muon2bdt");
    fDoList.push_back("muonsbdt");
    fDoList.push_back("muonmbdt");

    fDoList.push_back("pt");
    fDoList.push_back("p");
    fDoList.push_back("pz");
    fDoList.push_back("eta");
    fDoList.push_back("phi");
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

    fDoList.push_back("osmdr");

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

    fDoList.push_back("kaonspt");
    fDoList.push_back("kaonseta");
    fDoList.push_back("kaonsphi");

    fDoList.push_back("phidr");
    fDoList.push_back("mkk");
    // -- START for tracking studies
    if (0) {
      fDoList.push_back("muontkqual");
      fDoList.push_back("muonalg");
      fDoList.push_back("muonvalhits");
      fDoList.push_back("muonpixhits");
      fDoList.push_back("muontrkhits");
      fDoList.push_back("muonvalhitfraction");
      fDoList.push_back("muonlayerswithhits");
      fDoList.push_back("muonchi2");
      fDoList.push_back("muondz");
      fDoList.push_back("muondzE");
      fDoList.push_back("muond0");
      fDoList.push_back("muond0E");
      fDoList.push_back("muondsz");
      fDoList.push_back("muondszE");
      fDoList.push_back("muondxy");
      fDoList.push_back("muondxyE");
      fDoList.push_back("muonptE");
      //      fDoList.push_back("muonptEpt");
      fDoList.push_back("muonetaE");
      fDoList.push_back("muonphiE");
      fDoList.push_back("kaontkqual");
      fDoList.push_back("kaonalg");
      fDoList.push_back("kaonvalhits");
      fDoList.push_back("kaonpixhits");
      fDoList.push_back("kaontrkhits");
      fDoList.push_back("kaonvalhitfraction");
      fDoList.push_back("kaonlayerswithhits");
      fDoList.push_back("kaonchi2");
      fDoList.push_back("kaondz");
      fDoList.push_back("kaondzE");
      fDoList.push_back("kaond0");
      fDoList.push_back("kaond0E");
      fDoList.push_back("kaondsz");
      fDoList.push_back("kaondszE");
      fDoList.push_back("kaondxy");
      fDoList.push_back("kaondxyE");
      fDoList.push_back("kaonptE");
      fDoList.push_back("kaonptEpt");
      fDoList.push_back("kaonetaE");
      fDoList.push_back("kaonphiE");
      // -- END for tracking studies
    }
  }

  fChannelList.clear();
  // -- no restrictions, just normal analysis channels
  for (unsigned int i = 0; i < fNchan; ++i) {
    fChannelList.push_back(Form("%d", i));
  }


  if (0) {
    // -- small N(PV)
    fChannelList.push_back("0lopu");
    //  fChannelList.push_back("1lopu");

    // -- high N(PV)
    fChannelList.push_back("0hipu");
    //  fChannelList.push_back("1hipu");

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

    // -- top left
    fChannelList.push_back("0tl");
    // -- top right
    fChannelList.push_back("0tr");
    // -- bottom left
    fChannelList.push_back("0bl");
    // -- bottom right
    fChannelList.push_back("0br");
  }

}


// ----------------------------------------------------------------------
plotReducedOverlays::~plotReducedOverlays() {

}


// ----------------------------------------------------------------------
void plotReducedOverlays::init() {

  system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  system(Form("/bin/rm -f %s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str()));
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

  if (what == "bdtopt") {
    fChannelList.clear();
    for (unsigned int i = 0; i < fNchan; ++i) {
      fChannelList.push_back(Form("%d", i));
    }
    init();
    // makeSample("bupsikData", 100000);
    // makeSample("bupsikMcComb", 10000);
    makeSample("bupsikData");
    makeSample("bupsikMcComb");
    fStampString = "nada";
    makeOverlay("bupsikData", "bupsikMcComb", "bdt");
    makeOverlay("bupsikData", "bupsikMcComb", "cnc");
    return;
  }

  if (what == "bdtoptplot") {
    fChannelList.clear();
    for (unsigned int i = 0; i < fNchan; ++i) {
      fChannelList.push_back(Form("%d", i));
    }
    system(Form("/bin/rm -f %s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str()));
    system(Form("/bin/rm -f %s/overlay%s_ad*_*.pdf", fDirectory.c_str(), fSuffix.c_str()));
    system(Form("/bin/rm -f %s/mass_ad*_*.pdf", fDirectory.c_str()));

    fStampString = "nada";
    makeOverlay("bupsikData", "bupsikMcComb", "bdt");
    sysBdtCut("bupsikData", "bupsikMcComb", "bdt");
    sysComparison("bupsikData", "bupsikMcComb", "bdt");
    return;
  }

  if (what == "dbx") {
    // sysDoubleRatioFromMassFits("bspsiphi", "bupsik", "ad0bdt");
    // sysDoubleRatioFromMassFits("bspsiphi", "bupsik", "ad1bdt");
    // sysDoubleRatio("bdmm", "bupsik", "muonmbdt");
  }

  if (what == "dbx1") {
    fChannelList.clear();
    for (unsigned int i = 0; i < fNchan; ++i) {
      fChannelList.push_back(Form("%d", i));
    }

    system(Form("/bin/rm -f %s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str()));
    makeSampleOverlay("bspsiphiData", "bspsiphiMcComb", 100000);
    return;
  }

  if (what == "dbx2") {
    fDoCNC = false;
    fChannelList.clear();
    for (unsigned int i = 0; i < fNchan; ++i) {
      fChannelList.push_back(Form("%d", i));
    }

    system(Form("/bin/rm -f %s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str()));
    makeSampleOverlay("bspsiphiData", "bspsiphiMcComb");
    makeSampleOverlay("bupsikData", "bupsikMcComb");
    //    sysDoubleRatioFromMassFits("bspsiphi", "bupsik", "ad0bdt");
    return;
  }

  if (what == "nobs") {
    if (0) {
      fChannelList.clear();
      for (unsigned int i = 0; i < fNchan; ++i) {
	fChannelList.push_back(Form("%d", i));
      }
    }
    system(Form("/bin/rm -f %s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str()));
    makeSampleOverlay("bupsikData", "bupsikMcV09");
    return;
  }

  if (what == "wbs") {
    if (0) {
      fChannelList.clear();
      for (unsigned int i = 0; i < fNchan; ++i) {
	fChannelList.push_back(Form("%d", i));
      }
    }
    system(Form("/bin/rm -f %s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str()));
    makeSampleOverlay("bupsikData", "bupsikMcV08");
    return;
  }

  if (what == "mcbs") {
    if (1) {
      fChannelList.clear();
      for (unsigned int i = 0; i < fNchan; ++i) {
	fChannelList.push_back(Form("%d", i));
      }
    }
    system(Form("/bin/rm -f %s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str()));
    makeSampleOverlay("bupsikMcV08", "bupsikMcV09");
    return;
  }

  if (what == "dbxplot") {
    fChannelList.clear();
    for (unsigned int i = 0; i < fNchan; ++i) {
      fChannelList.push_back(Form("%d", i));
    }

    system(Form("/bin/rm -f %s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str()));
    makeOverlay("bupsikData", "bupsikMcComb", "bdt");
    if (fDoCNC) makeOverlay("bupsikData", "bupsikMcComb", "cnc");

    return;
  }

  if (what == "dbxo2f") {
    overlay2Files("2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad0bdt", "Presel",
 		  "2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad0cnc", "Ao");
    overlay2Files("2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad1bdt", "Presel",
 		  "2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad1cnc", "Ao");
    overlay2Files("2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad0bdt", "Presel",
 		  "2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad0cnc", "Ao");
    overlay2Files("2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad1bdt", "Presel",
 		  "2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad1cnc", "Ao");

    overlay2Files("2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad0bdt", "Presel",
      		  "2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad0bdt", "Presel");
    overlay2Files("2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad1bdt", "Presel",
		  "2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad1bdt", "Presel");
    overlay2Files("2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad0cnc", "Presel",
		  "2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad0cnc", "Presel");
    overlay2Files("2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad1cnc", "Presel",
		  "2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad1cnc", "Presel");
    overlay2Files("2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad0cnc", "Ao",
		  "2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad0cnc", "Ao");
    overlay2Files("2016BF-00/plotSbsHistograms-2016BF-00.root", "bupsikData", "ad1cnc", "Ao",
		  "2016GH-00/plotSbsHistograms-2016GH-00.root", "bupsikData", "ad1cnc", "Ao");

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

  if (what == "2017") {
    init();
    // -- data vs combined MC
    makeSampleOverlay("bmmData", "bdmmMcComb");
    makeSampleOverlay("bupsikData", "bupsikMcComb");
  }


  if (what == "all") {
    if (1) {
      fChannelList.clear();
      for (unsigned int i = 0; i < fNchan; ++i) {
	fChannelList.push_back(Form("%d", i));
      }
    }

    init();

    // -- data vs combined MC
    makeSampleOverlay("bspsiphiData", "bspsiphiMcComb");
    makeSampleOverlay("bupsikData", "bupsikMcComb");
    makeSampleOverlay("bmmData", "bdmmMcComb");
    makeSampleOverlay("bdpsikstarData", "bdpsikstarMcComb");

    allSystematics();

    // -- validation of private MC vs official MC
  }

  if (string::npos != what.find("validatemc")) {
    if (2016 == fYear) makeSampleOverlay("bupsikMcComb", "bupsikMcOff");
  }

  if (string::npos != what.find("plot")) {
    string hfname = Form("%s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str());
    system(Form("/bin/rm -f %s", hfname.c_str()));

    system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
    system(Form("/bin/rm -f %s/sbsctrl%d_ad*_*.pdf", fDirectory.c_str(), fYear));
    system(Form("/bin/rm -f %s/overlay%d_ad*_*.pdf", fDirectory.c_str(), fYear));
    system(Form("/bin/rm -f %s/mass%d_ad*_*.pdf", fDirectory.c_str(), fYear));

    if (fDoCNC) {
      makeOverlay("bspsiphiData", "bspsiphiMcComb", "cnc");
      makeOverlay("bmmData", "bdmmMcComb", "cnc");
      makeOverlay("bupsikData", "bupsikMcComb", "cnc");
      //      makeOverlay("bdpsikstarData", "bdpsikstarMcComb", "cnc");
    }

    makeOverlay("bspsiphiData", "bspsiphiMcComb", "bdt");
    makeOverlay("bmmData", "bdmmMcComb", "bdt");
    makeOverlay("bupsikData", "bupsikMcComb", "bdt");
    //    makeOverlay("bdpsikstarData", "bdpsikstarMcComb", "bdt");
    allSystematics();
  }

}




// ----------------------------------------------------------------------
void plotReducedOverlays::scaleFactors2016(string filename) {
  TFile *f1 = TFile::Open(filename.c_str());

  vector<string> vars;
  vars.push_back("maxdoca");
  vars.push_back("fl3de");
  vars.push_back("fl3d");
  string sel = "Ao";
  string sample = "bupsik";
  TH1D *hm(0), *hd(0);
  string hname("");
  vector<string> channels;
  channels.push_back("ad0cnc");
  channels.push_back("ad1cnc");
  channels.push_back("ad0bdt");
  channels.push_back("ad1bdt");
  for (unsigned int ichan = 0; ichan < channels.size(); ++ichan) {
    if (string::npos != channels[ichan].find("bdt")) {
      sel = "Presel";
    } else {
      sel = "Ao";
    }
    for (unsigned int iv = 0; iv < vars.size(); ++iv) {
      hname = Form("sbs_%s_%sMcComb_%s%s", channels[ichan].c_str(), sample.c_str(), vars[iv].c_str(), sel.c_str());
      hm = (TH1D*)f1->Get(hname.c_str());
      if (!hm) {
	cout << "histogram " << hname << " not found!?" << endl;
	continue;
      }

      hname = Form("sbs_%s_%sData_%s%s", channels[ichan].c_str(), sample.c_str(), vars[iv].c_str(), sel.c_str());
      hd = (TH1D*)f1->Get(hname.c_str());
      if (!hd) {
	cout << "histogram " << hname << " not found!?" << endl;
	continue;
      }
      cout << "chan " << channels[ichan] << " " << vars[iv] << " data/mc: " << hd->GetMean()/hm->GetMean() << endl;
    }
  }
}


// ----------------------------------------------------------------------
void plotReducedOverlays::comparePrivateAndOfficial() {

  changeSetup(fDirectory, "plotReducedOverlays", "comparePrivateAndOfficial");
  makeSampleOverlay("bupsikMc", "bupsikMcOff");

}

// ----------------------------------------------------------------------
void plotReducedOverlays::makeSampleOverlay(string sample1, string sample2, int nevents) {
  makeSample(sample1, nevents);
  makeSample(sample2, nevents);

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
    cuts.push_back("Cu");
  }
  if (string::npos != selection.find("cnc")) {
    fStampString = "CNC";
    cuts.push_back("Ao");
    cuts.push_back("Cu");
    cuts.push_back("Nm");
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
  string hfname = Form("%s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str());
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
  }

  if (string::npos != fSample.find("bupsik")) {
    dir = "candAnaBu2JpsiK";
    fMode = BU2JPSIKP;
    fIsSignal = false;
  }

  if (string::npos != fSample.find("bdpsikstar")) {
    dir = "candAnaBd2JpsiKstar";
    fMode = BD2JPSIKSTAR;
    fIsSignal = false;
  }

  if (string::npos != fSample.find("bspsiphi")) {
    dir = "candAnaBs2JpsiPhi";
    fMode = BS2JPSIPHI;
    fIsSignal = false;
  }

  // -- define boxes in bookDistributions!
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
  cout << "loopOverTree: nevts = " << nevts << endl;
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
    if (jentry%step == 0) {
      TTimeStamp ts;
      cout << Form(" .. evt = %d", jentry)
	   << ", time now: " << ts.AsString("lc")
	   << endl;
    }

    candAnalysis();
    (this->*pF)();
  }

}



// ----------------------------------------------------------------------
void plotReducedOverlays::loopFunction1() {

  bool cmsswPresel =  (fb.m1pt>4.) && (fb.m2pt>4.)  && (fb.pvips < 5.) && (fb.flsxy > 4.) && (fb.maxdoca < 0.08) && (fb.chi2/fb.dof < 5.);
  bool candanaPresel =  (fb.alpha<0.2) && (fb.fls3d>4.) && (TMath::Abs(fb.pvips) < 4.) && (TMath::Abs(fb.pvip) < 0.02);

  if (fMode == BU2JPSIKP) {
    cmsswPresel = cmsswPresel && (fb.kpt > 0.6);
  }

  if (fMode == BS2JPSIPHI) {
    cmsswPresel = cmsswPresel && (fb.k1pt > 0.6) && (fb.k2pt > 0.6);
  }

  if (fMode == BD2JPSIKSTAR) {
    cmsswPresel = cmsswPresel && (fb.kpt > 0.6) && (fb.pipt > 0.6);
  }

  // -- tests
  //  cmsswPresel = cmsswPresel && (fb.m1pix  >=1) && (fb.m2pix >= 1);

  // // -- copy from plotClass
  // fPreselection = fb.hlt1 && fb.tos && fb.l1t
  //   && fGoodAcceptance
  //   && fGoodGlobalMuons && fGoodMuonsPt && fGoodMuonsEta && fGoodTracksPt && fGoodTracksEta
  //   && cmsswPresel
  //   && fGoodJpsiCuts
  //   && fGoodDcand
  //   && (fb.alpha < 0.1) && (fb.fls3d > 5) && (fb.pvips < 4.) ;


  // -- local redefinition
  fGoodHLT = fb.hlt1 && fb.tos && fb.l1t
    && fGoodAcceptance
    && fGoodGlobalMuonsKin
    && fGoodTracksPt && fGoodTracksEta
    && cmsswPresel
    && candanaPresel
    && fGoodJpsiCuts
    && fGoodDcand
    && fGoodBdtPresel
    ;



  double bdtCut(0.1);
  if (-1 < fChan && fChan < fNchan) {
    bdtCut = fCuts[fChan]->bdtCut;
  }
  // -- remember: cmsswPresel (etc) is already part of fGoodHLT!
  fPreselection    = (fGoodHLT && (fb.alpha < 0.1) && (fb.fls3d > 5) && (fb.pvips < 4));
  fPreselectionBDT = (fGoodHLT && (fBDT > bdtCut));

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

  if (fGoodHLT
      && fGoodBDT
      ) {
    fSel3 = true;
  } else {
    fSel3 = false;
  }


  // -- update ana cuts!
  fCncCuts.update();
  fBdtCuts.update();

  bool loPU = (fb.pvn <  9);
  bool hiPU = (fb.pvn > 24);

  bool sfl = (fb.fls3d < 12.);
  bool bfl = (fb.fls3d > 50.);

  bool ldz = (TMath::Abs(fb.dzmin) > 1.2);
  bool sdz = (TMath::Abs(fb.dzmin) < 0.08);

  bool topr = fb.phi > 0.5*TMath::Pi();
  bool topl = fb.phi < 0.5*TMath::Pi() && fb.phi > 0.;
  bool botl = fb.phi < 0.              && fb.phi > -0.5*TMath::Pi();
  bool botr = fb.phi < -0.5*TMath::Pi();

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

  if (-1 < fChan && fChan < fNchan) {
    fChannel = Form("%d", fChan);
    if (fDoCNC) fillDistributions("cnc");
    fillDistributions("bdt");

    // -- skip the rest if we are on a 'normal'/'fast' running schedule
    if (fChannelList.size() == fNchan) return;

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

      if (botl) {
      	fChannel = Form("%dbl", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }
      if (botr) {
      	fChannel = Form("%dbr", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }
      if (topr) {
      	fChannel = Form("%dtr", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }
      if (topl) {

      	fChannel = Form("%dtl", fChan);
	if (fDoCNC) fillDistributions("cnc");
	fillDistributions("bdt");
      }

    }
  }
}


// ----------------------------------------------------------------------
void plotReducedOverlays::bookDistributions(std::string selmode) {
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

    if (fMode == BMM) {
      BGLBOXMIN = 4.90;
      BGLBOXMAX = 5.20;
      SIGBOXMIN = 5.20;
      SIGBOXMAX = 5.45;
      BGHBOXMIN = 5.45;
      BGHBOXMAX = 6.00;
    }

    if (fMode == BU2JPSIKP) {
      BGLBOXMIN = 5.02;
      BGLBOXMAX = 5.10;
      if (fCuts[i]->metaMax < 1.0) {
	SIGBOXMIN = 5.23;
	SIGBOXMAX = 5.33;
	BGHBOXMIN = 5.35;
	BGHBOXMAX = 5.45;
      } else if (fCuts[i]->metaMax < 1.5) {
	SIGBOXMIN = 5.23;
	SIGBOXMAX = 5.34;
	BGHBOXMIN = 5.37;
	BGHBOXMAX = 5.47;
      } else {
	SIGBOXMIN = 5.24; // must be tighter to reduce more strongly contribution from partially reco'ed background!
	SIGBOXMAX = 5.35;
	BGHBOXMIN = 5.42;
	BGHBOXMAX = 5.52;
      }
    }

    if (fMode == BD2JPSIKSTAR) {
      BGLBOXMIN = 5.00;
      BGLBOXMAX = 5.10;
      if (fCuts[i]->metaMax < 1.0) {
	SIGBOXMIN = 5.23;
	SIGBOXMAX = 5.33;
	BGHBOXMIN = 5.36;
      } else if (fCuts[i]->metaMax < 1.5) {
	SIGBOXMIN = 5.22;
	SIGBOXMAX = 5.34;
	BGHBOXMIN = 5.38;
      } else {
	SIGBOXMIN = 5.20;
	SIGBOXMAX = 5.38;
	BGHBOXMIN = 5.45;
      }
      BGHBOXMAX = 5.90;
    }

    if (fMode == BS2JPSIPHI) {
      if (fCuts[i]->metaMax < 1.0) {
	BGLBOXMIN = 5.20;
	BGLBOXMAX = 5.25;
	SIGBOXMIN = 5.35;
	SIGBOXMAX = 5.40;
	BGHBOXMIN = 5.45;
	BGHBOXMAX = 5.50;
      } else if (fCuts[i]->metaMax < 1.5) {
	BGLBOXMIN = 5.20;
	BGLBOXMAX = 5.25;
	SIGBOXMIN = 5.35;
	SIGBOXMAX = 5.40;
	BGHBOXMIN = 5.45;
	BGHBOXMAX = 5.50;
      } else {
	BGLBOXMIN = 5.20;
	BGLBOXMAX = 5.25;
	SIGBOXMIN = 5.32;
	SIGBOXMAX = 5.40;
	BGHBOXMIN = 5.48;
	BGHBOXMAX = 5.53;
      }
    }

    cout << "SIG: " << SIGBOXMIN << " .. " << SIGBOXMAX << endl;
    cout << "BGL: " << BGLBOXMIN << " .. " << BGLBOXMAX << endl;
    cout << "BGH: " << BGHBOXMIN << " .. " << BGHBOXMAX << endl;



    a = new adset();
    a->fpPvN      = bookDistribution(Form("%spvn", name.c_str()), "N(PV) ", "fGoodHLT", pCuts, 40, 0., 40., p);
    a->fpPvNtrk   = bookDistribution(Form("%spvntrk", name.c_str()), "N_{trk}(PV) ", "fGoodHLT", pCuts, 40, 0., 120., p);
    a->fpPv2Ntrk  = bookDistribution(Form("%spv2ntrk", name.c_str()), "N_{trk}(PV2) ", "fGoodHLT", pCuts, 40, 0., 120., p);
    a->fpDzmin    = bookDistribution(Form("%sdzmin", name.c_str()), "min(#Delta z) [cm] ", "fGoodHLT", pCuts, 50, -2., 2., p);
    a->fpDz12     = bookDistribution(Form("%sdz12", name.c_str()), "#Delta z [cm] ", "fGoodHLT", pCuts, 50, -2., 2., p);

    a->fpPvZ      = bookDistribution(Form("%spvz", name.c_str()), "z_{PV} [cm]", "fGoodHLT", pCuts, 40, -20., 20., p);
    a->fpPvAveW8  = bookDistribution(Form("%spvavew8", name.c_str()), "<w^{PV}>", "fGoodPvAveW8", pCuts, 50, 0.5, 1., p);

    if (bdt) {
      a->fpBDT      = bookDistribution(Form("%sbdt", name.c_str()), "BDT", "fGoodGlobalMuonsKin", pCuts, 100, -1.0, 1.0, &fPreselection); // else you have bdt>cut from fPreselectionBDT
      a->fpBDTSel0  = bookDistribution(Form("%sbdtsel0", name.c_str()), "BDTsel0", "fGoodEta", pCuts, 100, -1.0, 1.0, &fSel0);
      a->fpBDTSel1  = bookDistribution(Form("%sbdtsel1", name.c_str()), "BDTsel1", "fGoodEta", pCuts, 100, -1.0, 1.0, &fSel1);
      a->fpBDTSel2  = bookDistribution(Form("%sbdtsel2", name.c_str()), "BDTsel2", "fGoodEta", pCuts, 100, -1.0, 1.0, &fSel2);
      a->fpBDTSel3  = bookDistribution(Form("%sbdtsel3", name.c_str()), "BDTsel3", "fGoodGlobalMuonsKin", pCuts, 100, -1.0, 1.0, &fSel3);
    } else {
      a->fpBDT      = bookDistribution(Form("%sbdt", name.c_str()), "BDT", "fGoodCloseTrack", pCuts, 100, -1.0, 1.0, p);
      a->fpBDTSel0  = bookDistribution(Form("%sbdtsel0", name.c_str()), "BDTsel0", "fGoodMaxDoca", pCuts, 100, -1.0, 1.0, p);
      a->fpBDTSel1  = bookDistribution(Form("%sbdtsel1", name.c_str()), "BDTsel1", "fGoodFLS", pCuts, 100, -1.0, 1.0, p);
      a->fpBDTSel2  = bookDistribution(Form("%sbdtsel2", name.c_str()), "BDTsel2", "fGoodIpS", pCuts, 100, -1.0, 1.0, p);
      a->fpBDTSel3  = bookDistribution(Form("%sbdtsel3", name.c_str()), "BDTsel3", "fGoodIpS", pCuts, 100, -1.0, 1.0, p);
    }

    a->fpMuon1Pt   = bookDistribution(Form("%smuon1pt", name.c_str()), "#it{p}_{T, #mu1} [GeV]", "fGoodMuonsPt", pCuts, 60, 0., 30., p);
    a->fpMuon2Pt   = bookDistribution(Form("%smuon2pt", name.c_str()), "#it{p}_{T, #mu2} [GeV]", "fGoodMuonsPt", pCuts, 40, 0., 20., p);
    a->fpMuonsEta  = bookDistribution(Form("%smuonseta", name.c_str()), "#eta_{#mu}", "fGoodMuonsEta", pCuts, 40, -2.5, 2.5, p);
    a->fpMuonsPhi  = bookDistribution(Form("%smuonsphi", name.c_str()), "#phi_{#mu}", "fGoodMuonsEta", pCuts, 40, -3.15, 3.15, p);

    a->fpMuon1Bdt   = bookDistribution(Form("%smuon1bdt", name.c_str()), "BDT_{T, #mu1} ", "fGoodGlobalMuonsKin", pCuts, 100, -1., 1., p);
    a->fpMuon2Bdt   = bookDistribution(Form("%smuon2bdt", name.c_str()), "BDT_{T, #mu2} ", "fGoodGlobalMuonsKin", pCuts, 100, -1., 1., p);
    a->fpMuonsBdt   = bookDistribution(Form("%smuonsbdt", name.c_str()), "BDT_{T, #mu} ",  "fGoodGlobalMuonsKin", pCuts, 100, -1., 1., p);
    a->fpMuonmBdt   = bookDistribution(Form("%smuonmbdt", name.c_str()), "BDT_{T, #mu}^{min} ",  "fGoodGlobalMuonsKin", pCuts, 100, -1., 1., p);

    a->fpPt        = bookDistribution(Form("%spt", name.c_str()), "#it{p}_{T}#it{(B)} [GeV]", "fGoodPt", pCuts, 60, 0., 60., p);
    a->fpP         = bookDistribution(Form("%sp", name.c_str()), "#it{p(B)} [GeV]", "fGoodPt", pCuts, 50, 0., 100., p);
    a->fpPz        = bookDistribution(Form("%spz", name.c_str()), "#it{p_{z}(B)} [GeV]", "fGoodPt", pCuts, 50, 0., 100., p);
    a->fpEta       = bookDistribution(Form("%seta", name.c_str()), "#eta#it{(B)}", "fGoodEta", pCuts, 40, -2.5, 2.5, p);
    a->fpPhi       = bookDistribution(Form("%sphi", name.c_str()), "#phi#it{(B)}", "fGoodEta", pCuts, 40, -3.15, 3.15, p);

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
				      pCuts, 50, 0., 100., p);
    a->fpOsmDr     = bookDistribution(Form("%sosmdr", name.c_str()), "osmdr", (bdt? "fGoodBDT": "fGoodFLS"),
				      pCuts, 50, 0., 5., p);


    a->fpDocaTrk   = bookDistribution(Form("%sdocatrk", name.c_str()), "#it{d}_{ca}^{0} [cm]", (bdt? "fGoodBDT": "fGoodDocaTrk"), pCuts, 50, 0., 0.20, p);
    a->fpIso       = bookDistribution(Form("%siso", name.c_str()),  "isolation", (bdt? "fGoodBDT": "fGoodIso"), pCuts, 26, 0., 1.04, p);
    a->fpM1Iso     = bookDistribution(Form("%sm1iso", name.c_str()),  "m1 isolation", (bdt? "fGoodBDT": "fGoodM1Iso"), pCuts, 26, 0., 1.04, p);
    a->fpM2Iso     = bookDistribution(Form("%sm2iso", name.c_str()),  "m2 isolation", (bdt? "fGoodBDT": "fGoodM2Iso"), pCuts, 26, 0., 1.04, p);

    a->fpCloseTrk  = bookDistribution(Form("%sclosetrk", name.c_str()),  "#it{N}_{trk}^{close}", (bdt? "fGoodBDT": "fGoodCloseTrack"),
				      pCuts, 10, 0., 10., p);
    a->fpTau       = bookDistribution(Form("%stau", name.c_str()), "#tau [ps]", (bdt? "fGoodBDT": "fGoodCloseTrack"), pCuts, 50, 0., 10., p);

    if (fMode == BU2JPSIKP || fMode == BD2JPSIKSTAR || fMode == BS2JPSIPHI) {
      a->fpKaonsPt   = bookDistribution(Form("%skaonspt", name.c_str()), "#it{p}_{T, K} [GeV]", "fGoodMuonsPt", pCuts, 40, 0., 20., p);
      a->fpKaonsEta  = bookDistribution(Form("%skaonseta", name.c_str()), "#eta_{K}", "fGoodMuonsEta", pCuts, 40, -2.5, 2.5, p);
      a->fpKaonsPhi  = bookDistribution(Form("%skaonsphi", name.c_str()), "#phi_{K}", "fGoodMuonsEta", pCuts, 40, -3.15, 3.15, p);
      if (fMode == BS2JPSIPHI) {
	a->fpPhiDeltaR  = bookDistribution(Form("%sphidr", name.c_str()), "deltaR(#phi)", "fGoodGlobalMuonsKin", pCuts, 40, 0., 1.0, p);
	a->fpMkk  = bookDistribution(Form("%smkk", name.c_str()), "m(K,K)", "fGoodGlobalMuonsKin", pCuts, 40, 0.96, 1.06, p);
      }
      if (0) {
	// -- START for tracking studies
	a->fpMutkqual           = bookDistribution(Form("%smuontkqual", name.c_str()), "tkqual", "fGoodMuonsPt", pCuts, 32, 0., 32., p);
	a->fpMualg              = bookDistribution(Form("%smuonalg", name.c_str()), "alg", "fGoodMuonsPt", pCuts, 15, 0., 15., p);
	a->fpMuvalhits          = bookDistribution(Form("%smuonvalhits", name.c_str()), "valhits", "fGoodMuonsPt", pCuts, 40, 0., 40., p);
	a->fpMutrkhits          = bookDistribution(Form("%smuontrkhits", name.c_str()), "trkhits", "fGoodMuonsPt", pCuts, 20, 0., 20., p);
	a->fpMupixhits          = bookDistribution(Form("%smuonpixhits", name.c_str()), "pixhits", "fGoodMuonsPt", pCuts, 10, 0., 10., p);
	a->fpMuvalhitfraction   = bookDistribution(Form("%smuonvalhitfraction", name.c_str()), "valhitfraction", "fGoodMuonsPt", pCuts, 51, 0., 1.02, p);
	a->fpMulayerswithhits   = bookDistribution(Form("%smuonlayerswithhits", name.c_str()), "layerswithhits", "fGoodMuonsPt", pCuts, 25, 0., 25., p);
	a->fpMuchi2             = bookDistribution(Form("%smuonchi2", name.c_str()), "chi2", "fGoodMuonsPt", pCuts, 40, 0., 40., p);
	a->fpMudz               = bookDistribution(Form("%smuondz", name.c_str()), "dz", "fGoodMuonsPt", pCuts, 60, -15., 15., p);
	a->fpMudzE              = bookDistribution(Form("%smuondzE", name.c_str()), "log10(dzE)", "fGoodMuonsPt", pCuts, 40, -3., -1., p);
	a->fpMud0               = bookDistribution(Form("%smuond0", name.c_str()), "d0", "fGoodMuonsPt", pCuts, 100, -0.5, 0.5, p);
	a->fpMud0E              = bookDistribution(Form("%smuond0E", name.c_str()), "log10(d0E)", "fGoodMuonsPt", pCuts, 50, -3., -1., p);
	a->fpMudsz              = bookDistribution(Form("%smuondsz", name.c_str()), "dsz", "fGoodMuonsPt", pCuts, 60, -15., 15., p);
	a->fpMudszE             = bookDistribution(Form("%smuondszE", name.c_str()), "log10(dszE)", "fGoodMuonsPt", pCuts, 50, -3., -1., p);
	a->fpMudxy              = bookDistribution(Form("%smuondxy", name.c_str()), "dxy", "fGoodMuonsPt", pCuts, 50, -1., 1., p);
	a->fpMudxyE             = bookDistribution(Form("%smuondxyE", name.c_str()), "log10(dxyE)", "fGoodMuonsPt", pCuts, 50, -3., -1., p);
	a->fpMuptE              = bookDistribution(Form("%smuonptE", name.c_str()), "ptE", "fGoodMuonsPt", pCuts, 40, 0., 0.4, p);
	a->fpMuptEpt            = bookDistribution(Form("%smuonptEpt", name.c_str()), "ptE/pt", "fGoodMuonsPt", pCuts, 40, 0., 0.02, p);
	a->fpMuetaE             = bookDistribution(Form("%smuonetaE", name.c_str()), "etaE", "fGoodMuonsPt", pCuts, 50, 0., 0.01, p);
	a->fpMuphiE             = bookDistribution(Form("%smuonphiE", name.c_str()), "phiE", "fGoodMuonsPt", pCuts, 50, 0., 0.01, p);
	a->fpKatkqual           = bookDistribution(Form("%skaontkqual", name.c_str()), "tkqual", "fGoodMuonsPt", pCuts, 32, 0., 32., p);
	a->fpKaalg              = bookDistribution(Form("%skaonalg", name.c_str()), "alg", "fGoodMuonsPt", pCuts, 15, 0., 15., p);
	a->fpKavalhits          = bookDistribution(Form("%skaonvalhits", name.c_str()), "valhits", "fGoodMuonsPt", pCuts, 40, 0., 40., p);
	a->fpKatrkhits          = bookDistribution(Form("%skaontrkhits", name.c_str()), "trkhits", "fGoodMuonsPt", pCuts, 20, 0., 20., p);
	a->fpKapixhits          = bookDistribution(Form("%skaonpixhits", name.c_str()), "pixhits", "fGoodMuonsPt", pCuts, 10, 0., 10., p);
	a->fpKavalhitfraction   = bookDistribution(Form("%skaonvalhitfraction", name.c_str()), "valhitfraction", "fGoodMuonsPt", pCuts, 51, 0., 1.02, p);
	a->fpKalayerswithhits   = bookDistribution(Form("%skaonlayerswithhits", name.c_str()), "layerswithhits", "fGoodMuonsPt", pCuts, 25, 0., 25., p);
	a->fpKachi2             = bookDistribution(Form("%skaonchi2", name.c_str()), "chi2", "fGoodMuonsPt", pCuts, 40, 0., 40., p);
	a->fpKadz               = bookDistribution(Form("%skaondz", name.c_str()), "dz", "fGoodMuonsPt", pCuts, 60, -15., 15., p);
	a->fpKadzE              = bookDistribution(Form("%skaondzE", name.c_str()), "log10(dzE)", "fGoodMuonsPt", pCuts, 40, -3., -1., p);
	a->fpKad0               = bookDistribution(Form("%skaond0", name.c_str()), "d0", "fGoodMuonsPt", pCuts, 100, -0.5, 0.5, p);
	a->fpKad0E              = bookDistribution(Form("%skaond0E", name.c_str()), "log10(d0E)", "fGoodMuonsPt", pCuts, 50, -3., -1., p);
	a->fpKadsz              = bookDistribution(Form("%skaondsz", name.c_str()), "dsz", "fGoodMuonsPt", pCuts, 60, -15., 15., p);
	a->fpKadszE             = bookDistribution(Form("%skaondszE", name.c_str()), "log10(dszE)", "fGoodMuonsPt", pCuts, 50, -3., -1., p);
	a->fpKadxy              = bookDistribution(Form("%skaondxy", name.c_str()), "dxy", "fGoodMuonsPt", pCuts, 50, -1., 1., p);
	a->fpKadxyE             = bookDistribution(Form("%skaondxyE", name.c_str()), "log10(dxyE)", "fGoodMuonsPt", pCuts, 50, -3., -1., p);
	a->fpKaptE              = bookDistribution(Form("%skaonptE", name.c_str()), "ptE", "fGoodMuonsPt", pCuts, 40, 0., 0.3, p);
	a->fpKaptEpt            = bookDistribution(Form("%skaonptEpt", name.c_str()), "ptE/pt", "fGoodMuonsPt", pCuts, 40, 0., 0.02, p);
	a->fpKaetaE             = bookDistribution(Form("%skaonetaE", name.c_str()), "etaE", "fGoodMuonsPt", pCuts, 50, 0., 0.005, p);
	a->fpKaphiE             = bookDistribution(Form("%skaonphiE", name.c_str()), "phiE", "fGoodMuonsPt", pCuts, 50, 0., 0.005, p);
	// -- END for tracking studies
      }
    }


    fAdMap.insert(make_pair(mapname, a));
    cout << "bookDistributions: mapname = " << mapname << endl;

  }

}


// ----------------------------------------------------------------------
void plotReducedOverlays::sbsDistributions(string sample, string selection, string what) {
  string sbsControlPlotsFileName = Form("sbsctrl%s", fSuffix.c_str());
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

  cout << "gDIRECTORY: "; gDirectory->pwd(); cout << endl;
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
      if (!h) {
	cout << "=> not found: " << bla << endl;
      } else {
	cout << "=> cloning into  "
	     << Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), lselection.c_str()) << endl;
	h = (TH1D*)h->Clone(Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), lselection.c_str()));
      }
    } else if (1 == type) {
      cout << "=> Looking for signal histogram " << Form("%s%s0", bla.c_str(), lselection.c_str()) << endl;
      h = (TH1D*)gDirectory->Get(Form("%s%s0", bla.c_str(), lselection.c_str()));
      if (!h) {
	cout << "=> not found: " << bla << endl;
      } else {
	cout << "=> cloning into  "
	     << Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), lselection.c_str()) << endl;
	h = (TH1D*)h->Clone(Form("sbs_%s_%s%s", sample.c_str(), fDoList[i].c_str(), lselection.c_str()));
      }
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
    if (h) {
      cout << "  Title: " << h->GetTitle() << " with integral: " << h->GetSumOfWeights() << endl;
      //      zeroNegativeEntries(h);
      //      cout << " after zeroNegativeEntries(h) integral: "  << h->GetSumOfWeights() << endl;
    }
  }
}


// ----------------------------------------------------------------------
void plotReducedOverlays::allSystematics() {

  // -- BDT cut systematics
  // sysBdtCut("bupsikData", "bupsikMcComb", "bdt");
  // sysBdtCut("bspsiphiData", "bspsiphiMcComb", "bdt");
  // sysBdtCut("bdpsikstarData", "bdpsikstarMcComb", "bdt");

  for (int i = 0; i < fNchan; ++i) {
    sysDoubleRatio("bspsiphi", "bupsik", Form("ad%dbdt", i), "muonmbdt", "HLT", fCuts[i]->muonbdt);
    sysDoubleRatio("bdmm", "bupsik", Form("ad%dbdt", i), "muonmbdt", "HLT", fCuts[i]->muonbdt);
    sysDoubleRatio("bspsiphi", "bupsik", Form("ad%dbdt", i), "bdt", "Presel", fCuts[i]->bdtCut);
    sysDoubleRatio("bspsiphi", "bupsik", Form("ad%dbdt", i), "bdt", "Cu", fCuts[i]->bdtCut);
    sysDoubleRatio("bspsiphi", "bupsik", Form("ad%dbdt", i), "bdt", "HLT", fCuts[i]->bdtCut);
  }

}

// ----------------------------------------------------------------------
void plotReducedOverlays::sysComparison(string sample1, string sample2, string selection, string file2) {
  c0->SetCanvasSize(700, 700);
  c0->cd();
  shrinkPad(0.15, 0.18);

  string hfname = Form("%s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str());
  TFile *f1 = TFile::Open(hfname.c_str());
  TFile *f2 = f1;
  if (file2 != "nada") {
    Form("%s/plotSbsHistograms-%s.root", fDirectory.c_str(), file2.c_str());
    f2 = TFile::Open(hfname.c_str());
  }

  vector<string> dolist;
  dolist.push_back("tau");
  dolist.push_back("fls3d");
  dolist.push_back("bdtsel0");
  dolist.push_back("iso");
  vector<string> cutlevel;
  cutlevel.push_back("Presel");
  //  cutlevel.push_back("HLT");

  TH1D *h1(0), *h2(0);
  string hname("");

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

	h1->Scale(1./h1->GetSumOfWeights());
	h2->Scale(1./h2->GetSumOfWeights());
	//FIXME h1->SetMinimum(0.);
	setHist(h1, fDS[sample1]);
	h1->Draw();
	setHist(h2, fDS[sample2]);
	h2->Draw("samehist");
	double ks = h1->KolmogorovTest(h2);
	double ch = h1->Chi2Test(h2, "NORM");
	double chi2, ndof;
	double ch2 = chi2TestErr(h1, h2, chi2, ndof);

	cout << "tests:  KS = " << ks << " p(chi2) = " << ch << " and the second version = " << ch2 << endl;
	tl->SetTextAngle(0.);
	tl->SetTextSize(0.02);
	tl->DrawLatexNDC(0.2, 0.92, Form("KS: %3.2f (%3.2e)", ks, ks));
	tl->DrawLatexNDC(0.5, 0.92, Form("P_{#chi}: %3.2f (%3.2e)", ch, ch));
	tl->DrawLatexNDC(0.7, 0.92, Form("P_{#chi}: %3.2f (%3.2e)", ch2, ch2));

	hname = Form("%d%s:ad%s%s_%s_%s_%s%s", fYear, fSetup.c_str(), fChannelList[i].c_str(), selection.c_str(),
		     sample1.c_str(), sample2.c_str(), dolist[id].c_str(), cutlevel[ic].c_str());
	//\vdef{2016GH-1709:ad0bdt_bupsikData_bupsikMcComb_tauPresel:sysKS:val}   {\ensuremath{{0.000 } } }
	fTEX << formatTex(ks, Form("%s:sysKS:val", hname.c_str()), 3) << endl;
	fTEX << "\\vdef{" << Form("%s:sysKS:val2", hname.c_str()) << "} {" <<  Form("%3.2e", ks) << "}" << endl;
	fTEX << formatTex(ch, Form("%s:sysCH:val", hname.c_str()), 3) << endl;
	fTEX << "\\vdef{" << Form("%s:sysCH:val2", hname.c_str()) << "} {" <<  Form("%3.2e", ch) << "}" << endl;

	hname = Form("sysComp%d%s_ad%s%s_%s_%s_%s%s", fYear, fSetup.c_str(), fChannelList[i].c_str(), selection.c_str(),
		     sample1.c_str(), sample2.c_str(), dolist[id].c_str(), cutlevel[ic].c_str());
	savePad(Form("sbso/%s.pdf", hname.c_str()));
      }
    }
  }

}


// ----------------------------------------------------------------------
// -- calculate double ratio [eps(mu|bspsiphiMc)/eps(mu|bupsiMc)] / [eps(mu|bspsiphiData)/eps(mu|bupsiData)] from sideband-subtracted distribution
void plotReducedOverlays::sysDoubleRatio(string sample1, string sample2, string chansel, string var, string cutlevel, double cut) {

  // -- add an epsilon
  cut += 1.e-6;

  c0->SetCanvasSize(700, 700);
  c0->cd();
  shrinkPad(0.15, 0.18);

  string hfname = Form("%s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str());
  cout << "open " << hfname << endl;
  TFile *f1 = TFile::Open(hfname.c_str());

  string h1mName = Form("%sMcComb", sample1.c_str());
  string h1dName = Form("%sData", sample1.c_str());
  string h2mName = Form("%sMcComb", sample2.c_str());
  string h2dName = Form("%sData", sample2.c_str());

  // -- data dimuons have no 's' or 'd' in name
  if (string::npos != h1dName.find("bsmm")) replaceAll(h1dName, "bsmm", "bmm");
  if (string::npos != h1dName.find("bdmm")) replaceAll(h1dName, "bdmm", "bmm");

  double nd1(0.), nd2(0.), Nd1(0.), Nd2(0.);
  double nm1(0.), nm2(0.), Nm1(0.), Nm2(0.);
  double epsd1(0.), epsd1E(0.), epsd2(0.), epsd2E(0.), depsd(0.), depsdE(0.);
  double epsm1(0.), epsm1E(0.), epsm2(0.), epsm2E(0.), depsm(0.), depsmE(0.);
  int lo(-1), hi(-1);

  string sysname("");
  string d1name, d2name, m1name, m2name;
  TH1D *d1(0), *d2(0), *m1(0), *m2(0);
  d1name = Form("sbs_%s_%s_%s%s", chansel.c_str(), h1dName.c_str(), var.c_str(), cutlevel.c_str());
  d1 = (TH1D*)f1->Get(d1name.c_str());
  d2name = Form("sbs_%s_%s_%s%s", chansel.c_str(), h2dName.c_str(), var.c_str(), cutlevel.c_str());
  d2 = (TH1D*)f1->Get(d2name.c_str());

  m1name = Form("sbs_%s_%s_%s%s", chansel.c_str(), h1mName.c_str(), var.c_str(), cutlevel.c_str());
  m1 = (TH1D*)f1->Get(m1name.c_str());
  m2name = Form("sbs_%s_%s_%s%s", chansel.c_str(), h2mName.c_str(), var.c_str(), cutlevel.c_str());
  m2 = (TH1D*)f1->Get(m2name.c_str());

  if (!d1 || !d2 || !m1 || !m2) {
    cout << "histogram(s) not found: " << endl;
    cout << d1name << ": " << d1 << endl;
    cout << d2name << ": " << d2 << endl;
    cout << m1name << ": " << m1 << endl;
    cout << m2name << ": " << m2 << endl;
    return;
  }

  int INTLO = d1->FindBin(0.);
  INTLO = 1;
  int ichan(0);
  if (string::npos != chansel.find("1")) ichan = 1;
  lo    = d1->FindBin(cut);
  hi    = d1->GetNbinsX();
  nd1   = d1->Integral(lo, hi);
  Nd1   = d1->Integral(INTLO, hi);
  epsd1 = nd1/Nd1;
  nd2   = d2->Integral(lo, hi);
  Nd2   = d2->Integral(INTLO, hi);
  epsd2 = nd2/Nd2;
  depsd = (epsd1/epsd2);
  epsd1E = dRatio(nd1, Nd1);
  epsd2E = dRatio(nd2, Nd2);
  depsdE = dRatio(epsd1, epsd1E, epsd2, epsd2E);

  INTLO = m1->FindBin(0.);
  INTLO = 1;
  lo    = m1->FindBin(cut);
  hi    = m1->GetNbinsX();
  nm1   = m1->Integral(lo, hi);
  Nm1   = m1->Integral(INTLO, hi);
  epsm1 = nm1/Nm1;
  nm2   = m2->Integral(lo, hi);
  Nm2   = m2->Integral(INTLO, hi);
  epsm2 = nm2/Nm2;
  depsm = (epsm1/epsm2);
  epsm1E = dRatio(nm1, Nm1);
  epsm2E = dRatio(nm2, Nm2);
  depsmE = dRatio(epsm1, epsm1E, epsm2, epsm2E);
  cout << "lo: " << lo << " hi: " << hi << " INTLO = " << INTLO << " epsd1: " << epsd1 << " epsd2 = " << epsd2 << endl;
  cout << "data: depsd = " << depsd << " +/- " << depsdE << " nd1: " << nd1 << "/" << Nd1 << " nd2: " << nd2 << "/" << Nd2 << " " << d1name << "/" << d2name << endl;
  cout << "mc:   depsm = " << depsm << " +/- " << depsmE << " nm1: " << nm1 << "/" << Nm1 << " nm2: " << nm2 << "/" << Nm2 << " " << m1name << "/" << m2name << endl;
  cout << d1->GetName() << "/" << d2->GetName() << endl;
  cout << m1->GetName() << "/" << m2->GetName() << endl;
  d1->Scale(1./d1->GetSumOfWeights());
  d2->Scale(1./d2->GetSumOfWeights());
  m1->Scale(1./m1->GetSumOfWeights());
  m2->Scale(1./m2->GetSumOfWeights());
  //FIXME	h1->SetMinimum(0.);
  setHist(d1, fDS[h1dName]);
  d1->SetMarkerColor(fDS[h1mName]->fColor);
  setHist(d2, fDS[h2dName]);
  d2->SetMarkerColor(fDS[h2mName]->fColor);
  d2->SetMarkerStyle(20);
  setHist(m1, fDS[h1mName]);
  setHist(m2, fDS[h2mName]);

  double ymax(d1->GetMaximum());
  if (d2->GetMaximum() > ymax) ymax = d2->GetMaximum();
  if (m1->GetMaximum() > ymax) ymax = m1->GetMaximum();
  if (m2->GetMaximum() > ymax) ymax = m2->GetMaximum();
  d1->SetMaximum(1.2*ymax);
  d1->SetMinimum(-0.02*ymax);
  d1->Draw();
  d2->Draw("samee");
  m1->Draw("histsame");
  m2->Draw("histsame");


  newLegend(0.2, 0.65, 0.5, 0.87);
  legg->SetHeader((chansel + "/" + cutlevel + "/" + fSetup).c_str());
  legg->SetTextSize(0.03);
  legg->AddEntry(d1, Form("%s #varepsilon = %4.3f#pm%4.3f", sample1.c_str(), epsd1, epsd1E), "p");
  legg->AddEntry(d2, Form("%s #varepsilon = %4.3f#pm%4.3f", sample2.c_str(), epsd2, epsd2E), "p");
  legg->AddEntry(m1, Form("%s #varepsilon = %4.3f#pm%4.3f", sample1.c_str(), epsm1, epsm1E), "l");
  legg->AddEntry(m2, Form("%s #varepsilon = %4.3f#pm%4.3f", sample2.c_str(), epsm2, epsm2E), "l");
  legg->Draw();
  tl->SetTextSize(0.03);
  tl->DrawLatexNDC(0.2, 0.60, Form("ratio(data): %5.4f #pm %5.4f", depsd, depsdE));
  tl->DrawLatexNDC(0.2, 0.55, Form("ratio(MC):  %5.4f #pm %5.4f", depsm, depsmE));
  tl->DrawLatexNDC(0.2, 0.50, Form("Systematics: %5.4f ", 1.-(depsd/depsm)));
  pl->DrawLine(cut, 0., cut, 1.2*ymax);
  tl->DrawLatex(cut, 1.3*ymax, Form("cut = %5.4f", cut));

  sysname = Form("%s:%s_%s_%s_%s_%s", fSetup.c_str(), chansel.c_str(), sample1.c_str(),  sample2.c_str(), var.c_str(), cutlevel.c_str());
  fTEX << formatTex(epsm1,  Form("%s_Eps1Mc:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(epsm2,  Form("%s_Eps2Mc:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(epsd1,  Form("%s_Eps1Dt:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(epsd2,  Form("%s_Eps2Dt:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(depsm,  Form("%s_ratioEpsMc:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(depsmE, Form("%s_ratioEpsMc:err", sysname.c_str()), 3) << endl;
  fTEX << formatTex(depsd,  Form("%s_ratioEpsDt:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(depsdE, Form("%s_ratioEpsDt:err", sysname.c_str()), 3) << endl;
  fTEX << formatTex(1.-(depsd/depsm),  Form("%s_ratioEps:sys", sysname.c_str()), 3) << endl;

  savePad(Form("sbso/sysDoubleRatio%s_%s_%s_%s_%s_chan%d.pdf", fSetup.c_str(), var.c_str(), sample1.c_str(), sample2.c_str(), cutlevel.c_str(), ichan));


}


// ----------------------------------------------------------------------
void plotReducedOverlays::sysBdtCut(string sample1, string sample2, string selection, string file2, bool bdtscan) {

  c0->SetCanvasSize(700, 700);
  c0->cd();
  shrinkPad(0.15, 0.18);

  string hfname = Form("%s/plotSbsHistograms-%s.root", fDirectory.c_str(), fSuffix.c_str());
  TFile *f1 = TFile::Open(hfname.c_str());
  TFile *f2 = f1;
  if (file2 != "nada") {
    Form("%s/plotSbsHistograms-%s.root", fDirectory.c_str(), file2.c_str());
    f2 = TFile::Open(hfname.c_str());
  }

  // -- bdtscan: update the value of fCuts[ichan]->bdtCut with the one found in the scanBDT.root file
  if (bdtscan) {
    string rname = Form("%s/scanBDT-%d%s.root", fDirectory.c_str(), fYear, fSetup.c_str());
    TFile *fr = TFile::Open(rname.c_str());
    if (fr && fr->IsOpen()) {
      TH1D *hcuts = (TH1D*)fr->Get("bdtCuts");
      if (hcuts) {
	for (int ic = 0; ic < fNchan; ++ic) {
	  cout << "Overriding BDT cut for chan " << ic << " from " << fCuts[ic]->bdtCut << " to " << hcuts->GetBinContent(ic+1) << endl;
	  fCuts[ic]->bdtCut = hcuts->GetBinContent(ic+1);
	}
      }
      fr->Close();
    }
  }
  vector<string> dolist;
  dolist.push_back("bdt");
  dolist.push_back("bdtsel0");
  dolist.push_back("bdtsel1");
  dolist.push_back("bdtsel2");
  vector<string> cutlevel;
  cutlevel.push_back("Presel");
  cutlevel.push_back("HLT");

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
	//FIXME	h1->SetMinimum(0.);
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
	tl->DrawLatexNDC(0.25, 0.85, Form("%s", fSetup.c_str()));
	tl->DrawLatexNDC(0.25, 0.80, Form("BDT > %4.3f", fCuts[ichan]->bdtCut));

	hname = Form("%s:ad%s%s_%s_%s_%s%s", fSetup.c_str(), fChannelList[i].c_str(), selection.c_str(),
		     sample1.c_str(), sample2.c_str(), dolist[id].c_str(), cutlevel[ic].c_str());
	fTEX << formatTex(eps1,  Form("%s:sysBdt:eps1V", hname.c_str()), 3) << endl;
	fTEX << formatTex(eps1E, Form("%s:sysBdt:eps1E", hname.c_str()), 3) << endl;
	fTEX << formatTex(eps2,  Form("%s:sysBdt:eps2V", hname.c_str()), 3) << endl;
	fTEX << formatTex(eps2E, Form("%s:sysBdt:eps2E", hname.c_str()), 3) << endl;
	fTEX << formatTex(deps,  Form("%s:sysBdt:depsV", hname.c_str()), 3) << endl;
	fTEX << formatTex(depsE, Form("%s:sysBdt:depsE", hname.c_str()), 3) << endl;

	hname = Form("sbso%s_ad%s%s_%s_%s_%s%s", fSetup.c_str(), fChannelList[i].c_str(), selection.c_str(),
		     sample1.c_str(), sample2.c_str(), dolist[id].c_str(), cutlevel[ic].c_str());

	tl->SetTextAngle(90.);
	tl->SetTextSize(0.025);
	tl->DrawLatexNDC(0.92, 0.150, hname.c_str());
	tl->SetTextAngle(0.);

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
	//FIXME h1->SetMinimum(0.);
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
  TH1D *h1o(0), *h2o(0);
  string n1, n2;
  bool restricted = (what != "");
  bool doLegend(true);
  bool leftLegend(false);
  string lselection(selection);
  for (unsigned int i = 0; i < fDoList.size(); ++i) {
    if (restricted) {
      if (string::npos == fDoList[i].find(what)) continue;
    }
    n1 =  Form("sbs_%s_%s%s", sample1.c_str(), fDoList[i].c_str(), lselection.c_str());
    n2 =  Form("sbs_%s_%s%s", sample2.c_str(), fDoList[i].c_str(), lselection.c_str());
    doLegend = true;
    if (string::npos != fDoList[i].find("eta")) doLegend = false;
    if (string::npos != fDoList[i].find("bdt")) doLegend = false;
    h1o = (TH1D*)gDirectory->Get(n1.c_str());
    h2o = (TH1D*)gDirectory->Get(n2.c_str());
    if (0 == h1o || 0 == h2o) {
      cout << "  histograms not found" << endl;
      continue;
    }

    h1 = (TH1D*)h1o->Clone(Form("clone_%s", n1.c_str()));
    if (fIncludeOverflowInLastBin) addOverflow(h1);
    cout << "n1: " << n1 << " -> h1: " << h1 << " -> h1o: " << h1o << " sow = " << h1o->GetSumOfWeights() << endl;
    h2 = (TH1D*)h2o->Clone(Form("clone_%s", n2.c_str()));
    if (fIncludeOverflowInLastBin) addOverflow(h2);
    cout << "n2: " << n2 << " -> h2: " << h2 << " -> h2o: " << h2o << " sow = " << h1o->GetSumOfWeights() << endl;
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
	  // -- change things for the V08/09 validation
	  if (string::npos != sample1.find("V08") && string::npos != sample2.find("V09")) {
	    h1string += " (w/ BS)";
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
	  // -- change things for the V08/V09 comparison
	  if (string::npos != sample1.find("V08") && string::npos != sample2.find("V09")) {
	    h2string += " (no BS)";
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
    c0->SaveAs(Form("%s/overlays/%d/overlay%s_%s_%s_%s_%s.pdf",
		    fDirectory.c_str(), fYear, fSuffix.c_str(), sample1.c_str(), sample2.c_str(), fDoList[i].c_str(), lselection.c_str()));



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
void plotReducedOverlays::overlay2Files(std::string file1, std::string sample1, std::string chan1, std::string selection1,
					std::string file2, std::string sample2, std::string chan2, std::string selection2,
					std::string what) {

  bool restricted = (what != "");
  string hfname1 = Form("%s", file1.c_str());
  TFile *f1 = TFile::Open(hfname1.c_str());
  string hfname2 = Form("%s", file2.c_str());
  TFile *f2 = TFile::Open(hfname2.c_str());

  rmPath(hfname1);
  replaceAll(hfname1, ".root", "");
  replaceAll(hfname1, "plotSbsHistograms-", "");
  rmPath(hfname2);
  replaceAll(hfname2, ".root", "");
  replaceAll(hfname2, "plotSbsHistograms-", "");

  string ds1 = sample1.substr(sample1.find("_")+1);
  string ds2 = sample2.substr(sample2.find("_")+1);

  string n1(""), n2("");
  TH1D *h1o(0), *h2o(0);
  TH1D *h1(0), *h2(0);
  bool doLegend(false), leftLegend(false);
  for (int id = 0; id < fDoList.size(); ++id) {
    if (restricted) {
      if (string::npos == fDoList[id].find(what)) continue;
    }
    n1 =  Form("sbs_%s_%s_%s%s", chan1.c_str(), sample1.c_str(), fDoList[id].c_str(), selection1.c_str());
    n2 =  Form("sbs_%s_%s_%s%s", chan2.c_str(), sample2.c_str(), fDoList[id].c_str(), selection2.c_str());
    doLegend = true;
    if (string::npos != fDoList[id].find("eta")) doLegend = false;
    if (string::npos != fDoList[id].find("bdt")) doLegend = false;
    h1o = (TH1D*)f1->Get(n1.c_str());
    h2o = (TH1D*)f2->Get(n2.c_str());
    if (0 == h1o || 0 == h2o) {
      cout << "  histograms not found" << endl;
      continue;
    }
    h1 = (TH1D*)h1o->Clone(Form("clone_%s", n1.c_str()));
    if (fIncludeOverflowInLastBin) addOverflow(h1);
    // cout << "n1: " << n1 << " -> " << h1 << endl;
    h2 = (TH1D*)h2o->Clone(Form("clone_%s", n2.c_str()));
    if (fIncludeOverflowInLastBin) addOverflow(h2);
    // cout << "n2: " << n2 << " -> " << h2 << endl;
    if (h2->GetSumOfWeights() > 0) {
      h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
    }

    // cout << "setHist for " << ds1 << " and " << ds2 << endl;
    // cout << "fStampString = " << fStampString << endl;

    h1->SetNdivisions(505, "X");

    setHist(h1, fDS[ds1]);
    setHist(h2, fDS[ds2]);

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

    // -- change things for the legacy/rereco validation
    if (string::npos != sample1.find("bupsikData") && string::npos != sample2.find("bupsikDataLegacy")) {
      h1string += " (rereco)";
      sprintf(loption1, "p");
      h1->SetMarkerStyle(24);
      h1->SetMarkerSize(1.5);
      h1->Draw("samehist");
      string bla = h2->GetName();
      replaceAll(bla, "Data", "");
      h2->SetMarkerStyle(25);
      h2->SetMarkerSize(1.5);
      h2->SetName(bla.c_str());
      h2string += " (legacy)";
      sprintf(loption2, "p");
    }
    // -- change things for the 2016BF/GH comparison
    if (string::npos != file1.find("2016BF") && string::npos != file2.find("2016GH")) {
      h1string += Form(" (2016BF/%s)", chan1.c_str());
      sprintf(loption1, "p");
      h1->SetMarkerStyle(24);
      h1->SetMarkerSize(1.5);
      //dbx h1->Draw("samee");
      string bla = h2->GetName();
      replaceAll(bla, "Data", "");
      h2->SetMarkerStyle(25);
      h2->SetMarkerSize(1.5);
      h2->SetName(bla.c_str());
      h2string += Form(" (2016GH/%s)", chan2.c_str());
      //dbx h2->Draw("samehist");
      sprintf(loption2, "f");
    }
    // -- change things for the cnc/bdt comparison
    if (((string::npos != file1.find("2016BF") && string::npos != file2.find("2016BF"))
	 || ((string::npos != file1.find("2016GH") && string::npos != file2.find("2016GH"))))
	&& ((string::npos != chan1.find("bdt")) && (string::npos != chan2.find("cnc")))) {
      h1string = Form(" (%s/%s)", hfname1.c_str(), chan1.c_str());
      h1->SetMarkerStyle(24);
      h1->SetMarkerSize(1.5);
      // dbx h1->Draw("samee");
      sprintf(loption1, "p");
      string bla = h2->GetName();
      replaceAll(bla, "Data", "");
      h2->SetMarkerStyle(25);
      h2->SetMarkerSize(1.5);
      h2->SetName(bla.c_str());
      h2string = Form(" (%s/%s)", hfname2.c_str(), chan2.c_str());
      //dbx h2->Draw("samehist");
      sprintf(loption2, "f");
    }

    overlayAndRatio(c0, h1, h2);

    if (doLegend) {
      if (leftLegend) {
	newLegend(0.21, 0.7, 0.41, 0.87);
      } else {
	newLegend(0.50, 0.7, 0.75, 0.87);
      }

      legg->SetHeader(header.c_str());
      legg->SetTextSize(0.05);
      legg->AddEntry(h1, h1string.c_str(), loption1);
      legg->AddEntry(h2, h2string.c_str(), loption2);

      legg->Draw();
    } else {
      legg = 0;
    }


    if (1) {
      TLatex ll;
      ll.SetTextAngle(90.);
      ll.SetTextSize(0.03);
      ll.DrawLatexNDC(0.97, 0., Form("%s/%s/%s/%s/%s/%s/%s", sample1.c_str(), hfname1.c_str(), selection1.c_str(),
				     sample2.c_str(), hfname2.c_str(), selection2.c_str(), fDoList[id].c_str()));
    }

    c0->Modified();
    c0->Update();
    c0->SaveAs(Form("%s/o2f_%s_%s_%s_%s_%s_%s_%s_%s_%s.pdf",
		    fDirectory.c_str(),
		    hfname1.c_str(), sample1.c_str(), chan1.c_str(), selection1.c_str(),
		    hfname2.c_str(), sample2.c_str(), chan2.c_str(), selection2.c_str(),
		    fDoList[id].c_str()));



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
void plotReducedOverlays::loopFunction2() {

}


// ----------------------------------------------------------------------
void plotReducedOverlays::fillDistributions(string selmode) {

  string mapname = Form("ad%s%s_%s", fChannel.c_str(), selmode.c_str(), fSample.c_str());
  //  cout << "fillDistributions: mapname = " << mapname << endl;
  double mass = fb.cm;
  //if (fIsMC) mass = fb.m;
  if (fIsSignal) mass = fb.m;
  //  mass = fb.m;

  double ps   = static_cast<double>(fb.ps);
  if (fYear < 2016) {
    ps = 1.;
  }
  double w8 = fb.corrW8*ps;

  TLorentzVector a;
  a.SetPtEtaPhiM(fb.pt, fb.eta, fb.phi, mass);

  fAdMap[mapname]->fpMuon1Pt->fill(fb.m1pt, mass, w8);
  fAdMap[mapname]->fpMuon2Pt->fill(fb.m2pt, mass, w8);
  fAdMap[mapname]->fpMuonsEta->fill(fb.m1eta, mass, w8);
  fAdMap[mapname]->fpMuonsEta->fill(fb.m2eta, mass, w8);
  fAdMap[mapname]->fpMuonsPhi->fill(fb.m1phi, mass, w8);
  fAdMap[mapname]->fpMuonsPhi->fill(fb.m2phi, mass, w8);

  fAdMap[mapname]->fpMuon1Bdt->fill(fb.m1mvabdt, mass, w8);
  fAdMap[mapname]->fpMuon2Bdt->fill(fb.m2mvabdt, mass, w8);
  fAdMap[mapname]->fpMuonsBdt->fill(fb.m1mvabdt, mass, w8);
  fAdMap[mapname]->fpMuonsBdt->fill(fb.m2mvabdt, mass, w8);
  if (fb.m1mvabdt < fb.m2mvabdt) {
    fAdMap[mapname]->fpMuonmBdt->fill(fb.m1mvabdt, mass, w8);
  } else {
    fAdMap[mapname]->fpMuonmBdt->fill(fb.m2mvabdt, mass, w8);
  }

  fAdMap[mapname]->fpPt->fill(fb.pt, mass, w8);
  fAdMap[mapname]->fpP->fill(a.P(), mass, w8);
  fAdMap[mapname]->fpPz->fill(a.Pz(), mass, w8);
  fAdMap[mapname]->fpEta->fill(fb.eta, mass, w8);
  fAdMap[mapname]->fpPhi->fill(fb.phi, mass, w8);
  fAdMap[mapname]->fpAlpha->fill(fb.alpha, mass, w8);

  fAdMap[mapname]->fpIso->fill(fb.iso, mass, w8);
  fAdMap[mapname]->fpCloseTrk->fill(fb.closetrk, mass, w8);
  fAdMap[mapname]->fpDocaTrk->fill(fb.docatrk, mass, w8);

  fAdMap[mapname]->fpChi2Dof->fill(fb.chi2/fb.dof, mass, w8);
  fAdMap[mapname]->fpPChi2Dof->fill(fb.pchi2dof, mass, w8);

  fAdMap[mapname]->fpFLS3d->fill(fb.fls3d, mass, w8);
  fAdMap[mapname]->fpFL3d->fill(fb.fl3d, mass, w8);
  fAdMap[mapname]->fpFL3dE->fill(fb.fl3dE, mass, w8);

  fAdMap[mapname]->fpFLSxy->fill(fb.flsxy, mass, w8);
  fAdMap[mapname]->fpFLxy->fill(fb.flxy, mass, w8);
  fAdMap[mapname]->fpFLxyE->fill(fb.flxy/fb.flsxy, mass, w8);

  fAdMap[mapname]->fpTau->fill(1.e12*fb.tau, mass, w8);

  fAdMap[mapname]->fpMaxDoca->fill(fb.maxdoca, mass, w8);
  fAdMap[mapname]->fpIp->fill(fb.pvip, mass, w8);
  fAdMap[mapname]->fpIpS->fill(fb.pvips, mass, w8);
  fAdMap[mapname]->fpPvZ->fill(fb.pvz, mass, w8);
  fAdMap[mapname]->fpDzmin->fill(fb.dzmin, mass, w8);
  fAdMap[mapname]->fpDz12->fill(fb.dz12, mass, w8);
  fAdMap[mapname]->fpPvN->fill(fb.pvn, mass, w8);
  fAdMap[mapname]->fpPvNtrk->fill(fb.pvntrk, mass, w8);
  fAdMap[mapname]->fpPv2Ntrk->fill(fb.pv2ntrk, mass, w8);
  fAdMap[mapname]->fpPvAveW8->fill(fb.pvw8, mass, w8);

  fAdMap[mapname]->fpM1Iso->fill(fb.m1iso, mass, w8);
  fAdMap[mapname]->fpM2Iso->fill(fb.m2iso, mass, w8);

  fAdMap[mapname]->fpLip->fill(TMath::Abs(fb.pvlip), mass, w8);
  fAdMap[mapname]->fpLipS->fill(TMath::Abs(fb.pvlips), mass, w8);

  fAdMap[mapname]->fpLip2->fill(TMath::Abs(fb.pv2lip), mass, w8);
  fAdMap[mapname]->fpOsmDr->fill(fb.osmdr, mass, w8);

  // fAdMap[mapname]->fpLastCut->fill(mass, mass);

  fAdMap[mapname]->fpOtherVtx->fill(fb.othervtx, mass, w8);
  fAdMap[mapname]->fpPvDchi2->fill(fb.pvdchi2, mass, w8);

  fAdMap[mapname]->fpBDT->fill(fBDT, mass, w8);
  fAdMap[mapname]->fpBDTSel0->fill(fBDT, mass, w8);
  fAdMap[mapname]->fpBDTSel1->fill(fBDT, mass, w8);
  fAdMap[mapname]->fpBDTSel2->fill(fBDT, mass, w8);
  fAdMap[mapname]->fpBDTSel3->fill(fBDT, mass, w8);


  if (fMode == BU2JPSIKP || fMode == BD2JPSIKSTAR || fMode == BS2JPSIPHI) {
    if (fAdMap[mapname]->fpMutkqual) {
      fAdMap[mapname]->fpMutkqual->fill(fb.m1tkqual, mass, w8);
      fAdMap[mapname]->fpMutkqual->fill(fb.m2tkqual, mass, w8);

      fAdMap[mapname]->fpMualg->fill(fb.m1alg, mass, w8);
      fAdMap[mapname]->fpMualg->fill(fb.m2alg, mass, w8);

      fAdMap[mapname]->fpMuvalhits->fill(fb.m1valhits, mass, w8);
      fAdMap[mapname]->fpMuvalhits->fill(fb.m2valhits, mass, w8);

      fAdMap[mapname]->fpMupixhits->fill(fb.m1pix, mass, w8);
      fAdMap[mapname]->fpMupixhits->fill(fb.m2pix, mass, w8);

      fAdMap[mapname]->fpMutrkhits->fill(fb.m1trk, mass, w8);
      fAdMap[mapname]->fpMutrkhits->fill(fb.m2trk, mass, w8);

      fAdMap[mapname]->fpMuvalhitfraction->fill(fb.m1valhitfraction, mass, w8);
      fAdMap[mapname]->fpMuvalhitfraction->fill(fb.m2valhitfraction, mass, w8);

      fAdMap[mapname]->fpMulayerswithhits->fill(fb.m1layerswithhits, mass, w8);
      fAdMap[mapname]->fpMulayerswithhits->fill(fb.m2layerswithhits, mass, w8);

      fAdMap[mapname]->fpMuchi2->fill(fb.m1chi2, mass, w8);
      fAdMap[mapname]->fpMuchi2->fill(fb.m2chi2, mass, w8);

      fAdMap[mapname]->fpMudz->fill(fb.m1dz, mass, w8);
      fAdMap[mapname]->fpMudz->fill(fb.m2dz, mass, w8);

      fAdMap[mapname]->fpMudzE->fill(TMath::Log10(fb.m1dzE), mass, w8);
      fAdMap[mapname]->fpMudzE->fill(TMath::Log10(fb.m2dzE), mass, w8);

      fAdMap[mapname]->fpMud0->fill(fb.m1d0, mass, w8);
      fAdMap[mapname]->fpMud0->fill(fb.m2d0, mass, w8);

      fAdMap[mapname]->fpMud0E->fill(TMath::Log10(fb.m1d0E), mass, w8);
      fAdMap[mapname]->fpMud0E->fill(TMath::Log10(fb.m2d0E), mass, w8);

      fAdMap[mapname]->fpMudsz->fill(fb.m1dsz, mass, w8);
      fAdMap[mapname]->fpMudsz->fill(fb.m2dsz, mass, w8);

      fAdMap[mapname]->fpMudszE->fill(TMath::Log10(fb.m1dszE), mass, w8);
      fAdMap[mapname]->fpMudszE->fill(TMath::Log10(fb.m2dszE), mass, w8);

      fAdMap[mapname]->fpMudxy->fill(fb.m1dxy, mass, w8);
      fAdMap[mapname]->fpMudxy->fill(fb.m2dxy, mass, w8);

      fAdMap[mapname]->fpMudxyE->fill(TMath::Log10(fb.m1dxyE), mass, w8);
      fAdMap[mapname]->fpMudxyE->fill(TMath::Log10(fb.m2dxyE), mass, w8);

      fAdMap[mapname]->fpMuptE->fill(fb.m1ptE, mass, w8);
      fAdMap[mapname]->fpMuptE->fill(fb.m2ptE, mass, w8);
      fAdMap[mapname]->fpMuptEpt->fill(fb.m1ptE/fb.m1pt, mass, w8);
      fAdMap[mapname]->fpMuptEpt->fill(fb.m2ptE/fb.m2pt, mass, w8);
      fAdMap[mapname]->fpMuetaE->fill(fb.m1etaE, mass, w8);
      fAdMap[mapname]->fpMuetaE->fill(fb.m2etaE, mass, w8);
      fAdMap[mapname]->fpMuphiE->fill(fb.m1phiE, mass, w8);
      fAdMap[mapname]->fpMuphiE->fill(fb.m2phiE, mass, w8);

      if (fMode == BU2JPSIKP) {
	fAdMap[mapname]->fpKatkqual->fill(fb.ktkqual, mass, w8);
	fAdMap[mapname]->fpKaalg->fill(fb.kalg, mass, w8);
	fAdMap[mapname]->fpKavalhits->fill(fb.kvalhits, mass, w8);
	fAdMap[mapname]->fpKapixhits->fill(fb.kpix, mass, w8);
	fAdMap[mapname]->fpKatrkhits->fill(fb.ktrk, mass, w8);
	fAdMap[mapname]->fpKavalhitfraction->fill(fb.kvalhitfraction, mass, w8);
	fAdMap[mapname]->fpKalayerswithhits->fill(fb.klayerswithhits, mass, w8);
	fAdMap[mapname]->fpKachi2->fill(fb.kchi2, mass, w8);
	fAdMap[mapname]->fpKadz->fill(fb.kdz, mass, w8);
	fAdMap[mapname]->fpKadzE->fill(TMath::Log10(fb.kdzE), mass, w8);
	fAdMap[mapname]->fpKad0->fill(fb.kd0, mass, w8);
	fAdMap[mapname]->fpKad0E->fill(TMath::Log10(fb.kd0E), mass, w8);
	fAdMap[mapname]->fpKadsz->fill(fb.kdsz, mass, w8);
	fAdMap[mapname]->fpKadszE->fill(TMath::Log10(fb.kdszE), mass, w8);
	fAdMap[mapname]->fpKadxy->fill(fb.kdxy, mass, w8);
	fAdMap[mapname]->fpKadxyE->fill(TMath::Log10(fb.kdxyE), mass, w8);
	fAdMap[mapname]->fpKaptE->fill(fb.kptE, mass, w8);
	fAdMap[mapname]->fpKaptEpt->fill(fb.kptE/fb.kpt, mass, w8);
	fAdMap[mapname]->fpKaetaE->fill(fb.ketaE, mass, w8);
	fAdMap[mapname]->fpKaphiE->fill(fb.kphiE, mass, w8);
      }
    }

    if (fMode == BU2JPSIKP) {
      fAdMap[mapname]->fpKaonsPt->fill(fb.kpt, mass, w8);
      fAdMap[mapname]->fpKaonsEta->fill(fb.keta, mass, w8);
      fAdMap[mapname]->fpKaonsPhi->fill(fb.kphi, mass, w8);
    } else if (fMode == BD2JPSIKSTAR) {
      fAdMap[mapname]->fpKaonsPt->fill(fb.kpt, mass, w8);
      fAdMap[mapname]->fpKaonsEta->fill(fb.keta, mass, w8);
      fAdMap[mapname]->fpKaonsPhi->fill(fb.kphi, mass, w8);
    } else if (fMode == BS2JPSIPHI) {
      fAdMap[mapname]->fpKaonsPt->fill(fb.k1pt, mass, w8);
      fAdMap[mapname]->fpKaonsEta->fill(fb.k1eta, mass, w8);
      fAdMap[mapname]->fpKaonsPhi->fill(fb.k1phi, mass, w8);
      fAdMap[mapname]->fpKaonsPt->fill(fb.k2pt, mass, w8);
      fAdMap[mapname]->fpKaonsEta->fill(fb.k2eta, mass, w8);
      fAdMap[mapname]->fpKaonsPhi->fill(fb.k2phi, mass, w8);

      fAdMap[mapname]->fpPhiDeltaR->fill(fb.phidr, mass, w8);
      fAdMap[mapname]->fpMkk->fill(fb.mkk, mass, w8);
    }

  }
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
  double ymax = h1->GetMaximum();
  if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();
  h1->SetStats(0);
  string hname = h1->GetName();
  if (string::npos != hname.find("Data")) {
    h1->Draw("e");
  } else {
    h1->Draw();
  }
  hname = h2->GetName();
  if (string::npos != hname.find("Data")) {
    h2->Draw("samee");
  } else {
    h2->Draw("samehist");
  }
  h1->SetMaximum(1.2*ymax);
  h1->SetMinimum(-0.02*ymax);

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

      if (string::npos != stype.find("bupsikV08")) {
        sname = "bupsikMcV08";
        sdecay = "bupsik";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

      if (string::npos != stype.find("bupsikV09")) {
        sname = "bupsikMcV09";
        sdecay = "bupsik";
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
  cout << Form("%-25s: ", "Dataset name") << "Filename" << endl;
  cout << "------------------------------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    cout << Form("%-25s: ", it->first.c_str());
    if (it->second->fF) {
      cout << it->second->fF->GetName() << endl;
    } else {
      cout << "##### file not opened?!?!?! #####" << endl;
    }
  }
  cout << "------------------------------------------------------------------------------------------" << endl;
}
