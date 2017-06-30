#include "plotTrigger.hh"

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

ClassImp(plotTrigger)

using namespace std;

// ----------------------------------------------------------------------
plotTrigger::plotTrigger(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotTrigger::loadFiles(files);

  changeSetup(dir, "plotTrigger", setup);
  init();

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

}


// ----------------------------------------------------------------------
plotTrigger::~plotTrigger() {

}


// ----------------------------------------------------------------------
void plotTrigger::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotTrigger::makeAll(string what) {
}


// ----------------------------------------------------------------------
void plotTrigger::plotL1Seeds(std::string dsname) {
  if (string::npos != dsname.find("bupsik")) fMode = BU2JPSIKP;
  if (string::npos != dsname.find("bspsiphi")) fMode = BS2JPSIPHI;
  if (string::npos != dsname.find("bdpsikstar")) fMode = BD2JPSIKSTAR;

  fSample = dsname;
  string dir = "candAnaMuMu";
  if (string::npos != fSample.find("bspsiphi")) {
    fMode = BS2JPSIPHI;
    dir  = "candAnaBs2JpsiPhi";
  } else if (string::npos != fSample.find("bupsik")) {
    fMode = BU2JPSIKP;
    dir  = "candAnaBu2JpsiK";
  } else if (string::npos != fSample.find("bdpsikstar")) {
    fMode = BD2JPSIKSTAR;
    dir  = "candAnaBd2JpsiKstar";
  }

  TH1D *h1(0);
  fHistFile = TFile::Open(fHistFileName.c_str(), "");
  if (fHistFile) {
    h1 = (TH1D*)fHistFile->Get(Form("h_%s_%s", "s0", dsname.c_str()));
    fHistFile->Close();
  }
  if (!h1) {
    // -- create histograms
    cout << "fHistFile: " << fHistFileName;
    fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
    cout << " opened " << endl;


    TTree *t = getTree(fSample, dir);
    if (0 == t) {
      cout << "tree for sample = " << fSample << " not found" << endl;
      return;
    }
    //  bookHist(fSample);
    fpHL1s0 = new TH1D(Form("h_%s_%s", "s0", dsname.c_str()), Form("h_%s_%s", "s0", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s1 = new TH1D(Form("h_%s_%s", "s1", dsname.c_str()), Form("h_%s_%s", "s1", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s2 = new TH1D(Form("h_%s_%s", "s2", dsname.c_str()), Form("h_%s_%s", "s2", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s3 = new TH1D(Form("h_%s_%s", "s3", dsname.c_str()), Form("h_%s_%s", "s3", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s4 = new TH1D(Form("h_%s_%s", "s4", dsname.c_str()), Form("h_%s_%s", "s4", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s5 = new TH1D(Form("h_%s_%s", "s5", dsname.c_str()), Form("h_%s_%s", "s5", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s6 = new TH1D(Form("h_%s_%s", "s6", dsname.c_str()), Form("h_%s_%s", "s6", dsname.c_str()), 7236, 273150, 280385);
    fpHL1All = new TH1D(Form("h_%s_%s", "All", dsname.c_str()), Form("h_%s_%s", "All", dsname.c_str()), 7236, 273150, 280385);

    setupTree(t, fSample);
    fCds = fDS[fSample];
    loopOverTree(t, 4);

    fpHL1s0->Write();
    fpHL1s1->Write();
    fpHL1s2->Write();
    fpHL1s3->Write();
    fpHL1s4->Write();
    fpHL1s5->Write();
    fpHL1s6->Write();
    fpHL1All->Write();

    fHistFile->Close();
  }

  fHistFile = TFile::Open(fHistFileName.c_str(), "");
  TH1D *hAll = (TH1D*)fHistFile->Get(Form("h_%s_%s", "All", dsname.c_str()));
  TH1D *hSeed(0), *hEff(0);
  gPad->SetGridx(1);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  vector<string> sseed;
  sseed.push_back("L1_DoubleMu0er1p6_dEtaMax1p8");
  sseed.push_back("L1_DoubleMu0er1p6_dEta_Max1p8_OS");
  sseed.push_back("L1_DoubleMu0er1p4_dEta_Max1p8_OS");
  sseed.push_back("L1_DoubleMu_10_0_dEta_Max1p8");
  sseed.push_back("L1_DoubleMu_11_4");
  sseed.push_back("L1_DoubleMu_12_5");

  for (int iseed = 0; iseed < 6; ++iseed) {
    hSeed = (TH1D*)fHistFile->Get(Form("h_%s_%s", Form("s%d", iseed), dsname.c_str()));
    hEff = (TH1D*)hSeed->Clone(Form("eff_%s", hSeed->GetName())); hEff->Sumw2(); hEff->Reset();
    hEff->SetTitle((sseed[iseed] + Form("(%s)", fDS[dsname]->fName.c_str())).c_str());
    cout << "hEff title = " << hEff->GetTitle() << endl;
    hEff->Divide(hSeed, hAll, 1., 1., "b");
    setTitles(hEff, "run", "fraction");
    hEff->Draw();
    savePad(Form("effL1Seed_%d_%s.pdf", iseed, dsname.c_str()));

  }
}





// ----------------------------------------------------------------------
void plotTrigger::plotTisEfficiency(string dsname) {

  // -- read histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  vector<string> vds;

  TIter next(fHistFile->GetListOfKeys());
  TKey *key(0);
  TH1D *hpass(0), *hnorm(0);
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;

    if (TString(key->GetName()).Contains("_norm_")) {
      string hname = key->GetName();
      replaceAll(hname, "h_norm_", "");
      vds.push_back(hname);
    }
  }

  string name;
  for (unsigned int i = 0; i < vds.size(); ++i) {
    if (dsname != "all" && vds[i] != dsname) continue;
    hpass = (TH1D*)(fHistFile->Get(Form("h_pass_%s", vds[i].c_str())));
    hnorm = (TH1D*)(fHistFile->Get(Form("h_norm_%s", vds[i].c_str())));

    hnorm->Draw();
    hpass->Draw("samee");

    name = vds[i];
    replaceAll(name, "Data", " Charmonium");
    replaceAll(name, "SingleMuon", " BLA");
    replaceAll(name, "BLA", " SingleMuon");

    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.16, 0.92, name.c_str());
    tl->DrawLatexNDC(0.60, 0.92, Form("#varepsilon = %d/%d = %3.2f #pm %3.2f",
				     static_cast<int>(hpass->GetSumOfWeights()),
				     static_cast<int>(hnorm->GetSumOfWeights()),
				     hpass->GetSumOfWeights()/hnorm->GetSumOfWeights(),
				     dEff(static_cast<int>(hpass->GetSumOfWeights()),
					  static_cast<int>(hnorm->GetSumOfWeights())
					  ))
		     );

    savePad(Form("eff-%s.pdf", vds[i].c_str()));
  }



}


// ----------------------------------------------------------------------
void plotTrigger::runTisEfficiency(string dsname) {

  if (string::npos != dsname.find("bupsik")) fMode = BU2JPSIKP;
  if (string::npos != dsname.find("bspsiphi")) fMode = BS2JPSIPHI;
  if (string::npos != dsname.find("bdpsikstar")) fMode = BD2JPSIKSTAR;

  // -- dump histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;

  fSample = dsname;
  string dir = "candAnaBu2JpsiK";

  TTree *t = getTree(fSample, dir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  }
  //  bookHist(fSample);
  fpHnorm = new TH1D(Form("h_%s_%s", "norm", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);
  fpHpass = new TH1D(Form("h_%s_%s", "pass", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);


  setupTree(t, fSample);
  fCds = fDS[fSample];
  loopOverTree(t, 1);


  fpHnorm->Write();
  fpHpass->Write();
  fHistFile->Close();

}


// ----------------------------------------------------------------------
void plotTrigger::refTrgEfficiency(string selection, string dsname) {


  setup(dsname);
  fSample = dsname;

  zone(2,2);

  TTree *t = getTree(fSample, fTreeDir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  }

  TH1D *h1 = new TH1D(Form("h_0"), Form("h_0"), 40, 4.8, 6.0);
  TH1D *h2 = new TH1D(Form("h_0_hlt"), Form("h_0_hlt"), 40, 4.8, 6.0);

  double nNorm(0.), nPass(0.), effMc(0.), effRt(0.), effMcE(0.), effRtE
    (0.);
  double mBp(5.28), sBp(0.05);

  string fitopt("lm");

  // -- basic HLT efficiency
  string tselection = selection;
  t->Draw("m >> h_0", tselection.c_str());
  tselection = selection + " && hlt";
  t->Draw("m >> h_0_hlt", tselection.c_str());

  TF1 *f1 = fIF->bupsik(h1);
  c0->cd(1);
  h1->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  f1->SetParameter(7, 0.);
  f1->SetParameter(8, 0.);
  f1->SetParameter(9, 0.);
  nNorm = f1->Integral(5.1, 5.4)/h1->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "sel. 0");
  savePad("h_0.pdf");
  delete f1;

  c0->cd(2);
  fIF->limitPar(1, 5.1, 5.5);
  f1 = fIF->bupsik(h2);
  h2->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  f1->SetParameter(7, 0.);
  f1->SetParameter(8, 0.);
  f1->SetParameter(9, 0.);
  nPass = f1->Integral(5.1, 5.4)/h2->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "sel. 0 && HLT");
  savePad("h_0_hlt.pdf");
  delete f1;
  effMc = nPass/nNorm;
  effMcE = dEff(static_cast<int>(nPass), static_cast<int>(nNorm));
  cout << "==> efficiency = " << nPass << "/" << nNorm << " = " << nPass/nNorm << endl;


  // -- HLT efficiency derived from reference trigger
  TH1D *h3 = new TH1D(Form("h_rt"), Form("h_rt"), 40, 4.8, 6.0);
  TH1D *h4 = new TH1D(Form("h_rt_hlt"), Form("h_rt_hlt"), 40, 4.8, 6.0);

  c0->cd(3);
  tselection = selection + " && reftrg";
  t->Draw("m >> h_rt", tselection.c_str());
  tselection = selection + " && reftrg && hlt";
  t->Draw("m >> h_rt_hlt", tselection.c_str());

  f1 = fIF->bupsik(h3);
  h3->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  f1->SetParameter(7, 0.);
  f1->SetParameter(8, 0.);
  f1->SetParameter(9, 0.);
  nNorm = f1->Integral(5.1, 5.4)/h3->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "sel. 0 && reftrg");
  savePad("h_rt.pdf");
  delete f1;

  c0->cd(4);
  f1 = fIF->bupsik(h4);
  h4->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  f1->SetParameter(7, 0.);
  f1->SetParameter(8, 0.);
  f1->SetParameter(9, 0.);
  nPass = f1->Integral(5.1, 5.4)/h4->GetBinWidth(1);
  effRt = nPass/nNorm;
  effRtE = dEff(static_cast<int>(nPass), static_cast<int>(nNorm));
  tl->DrawLatexNDC(0.2, 0.96, "sel. 0 && reftrg && HLT");
  savePad("h_rt_hlt.pdf");
  delete f1;
  cout << "==> efficiency = " << nPass << "/" << nNorm << " = " << nPass/nNorm << endl;
  cout << "==> cuts: " << selection << endl;
  cout << "==> efficiencies: MC = " << Form("%4.2f +/- %4.2f", effMc, effMcE)
       << "; ref trigger = " << Form("%4.2f +/- %4.2f", effRt, effRtE)
       << endl;

}



// ----------------------------------------------------------------------
void plotTrigger::loopFunction1() {

}


// ----------------------------------------------------------------------
void plotTrigger::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotTrigger::loopOverTree> loop over dataset " << (fCds?fCds->fName:"undefined") << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotTrigger::*pF)(void);
  if (ifunc == 1) pF = &plotTrigger::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotTrigger::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotTrigger::loadFile loading files listed in " << files << endl;

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

    if (string::npos != stype.find("SingleMuon")) {
      // -- SingleMuon
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("bmm")) {
        sname = "bmmSingleMuon";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bupsik")) {
        sname = "bupsikSingleMuon";
        sdecay = "bupsik";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bspsiphi")) {
        sname = "bspsiphiSingleMuon";
        sdecay = "bspsiphi";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bdpsikstar")) {
        sname = "bdpsikstarSingleMuon";
        sdecay = "bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

    } else if (string::npos != stype.find("Charmonium")) {
      // -- Charmonium
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("bmm")) {
        sname = "bmmCharmonium";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bupsik")) {
        sname = "bupsikCharmonium";
        sdecay = "bupsik";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bspsiphi")) {
        sname = "bspsiphiCharmonium";
        sdecay = "bspsiphi";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bdpsikstar")) {
        sname = "bdpsikstarCharmonium";
        sdecay = "bdpsikstar";
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

    } else if (string::npos != stype.find("relval")) {
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2.;
      if (string::npos != stype.find("bsmm,relval")) {
	ds->fF      = pF;
	sname   = "bsmmrelval";
      }

      if (string::npos != stype.find("bdmm,relval")) {
	ds->fF      = pF;
	sname   = "bdmmrelval";
      }

      if (string::npos != stype.find("bupsik,relval")) {
	ds->fF      = pF;
	sname   = "bupsikrelval";
      }

      if (string::npos != stype.find("bspsiphi,relval")) {
	ds->fF      = pF;
	sname   = "bspsiphirelval";
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
    if (!it->second->fF) {
      cout << "missing " << it->first << endl;
    }
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
