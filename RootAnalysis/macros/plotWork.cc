#include "plotWork.hh"

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

ClassImp(plotWork)

using namespace std;

// ----------------------------------------------------------------------
plotWork::plotWork(string dir, string files, string cuts, string setup, int year): plotClass(dir, files, cuts, setup, year) {
  plotClass::loadFiles(files);
  plotWork::loadFiles(files);

  changeSetup(dir, "plotWork", setup);
  init();

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

  MKKLO = 0.9;
  MKKHI = 1.2;
  DR    = 99.;
  PTK1  = 0.0;
  PTK2  = 0.0;
  PTPSI = 0.0;
}


// ----------------------------------------------------------------------
plotWork::~plotWork() {

}

// ----------------------------------------------------------------------
string plotWork::removeVarFromSelection(string var, string selection) {
  // this will split on &, not on &&. So every second element in cuts is empty
  cout << "remove var ->" << var << "<-" << endl;
  vector<string> cuts = split(selection, '&');
  string redcuts("");
  for (unsigned int i = 0; i < cuts.size(); ++i) {
    if (cuts[i] == "") continue;
    cout << "cuts[i]= ->" << cuts[i] << "<-" << endl;
    if ((string::npos == cuts[i].find(" "+var+"<")) && (string::npos == cuts[i].find(" "+var+">"))) {
      redcuts += cuts[i];
      if (i < cuts.size()-1) redcuts += " && ";
    }
  }
  return redcuts;

}

// ----------------------------------------------------------------------
// for removeVarFromSelection(...) to work, there MUST be
//   - a space in front of the variable
//   - no space between the variable and the <> operator
// !!
// Variables not adhering to this will not be removed.
string plotWork::selectionString(int imode, int itrig) {
  string selection("");
  if (10 == imode) {
    selection = "abs(m1eta)<1.6 && abs(m2eta)<1.6";
    selection += " && m1pt>4 && m2pt>4 && m1q*m2q<0";
    selection += " && pt>6.5";
    selection += " && fls3d>10 && alpha<0.05 && pvips<2 && pvip<0.008 && chi2dof<2.2";
    selection += " && iso>0.8 && docatrk>0.015 && closetrk<2";
  } else if (11 == imode) {
    selection = "abs(m1eta)<1.6 && abs(m2eta)<1.6";
    selection += " && m1pt>8 && m2pt>8 && m1q*m2q<0";
    selection += " && pt>6.5";
    selection += " && fls3d>10 && alpha<0.05 && pvips<2 && pvip<0.008 && chi2dof<2.2";
    selection += " && iso>0.8 && docatrk>0.015 && closetrk<2";
  } else if (12 == imode) {
    selection = "abs(m1eta)<1.6 && abs(m2eta)<1.6";
    selection += " && m1pt>8 && m2pt>8 && m1q*m2q<0";
    selection += " && pt>6.5";
    selection += " && reftrg";
    selection += " && fls3d>10 && alpha<0.05 && pvips<2 && pvip<0.008 && chi2dof<2.2";
    selection += " && iso>0.8 && docatrk>0.015 && closetrk<2";
  }

  if (0 == itrig) {
    selection += " ";
  } else if (1 == itrig) {
    selection += " && reftrg";
  } else if (2 == itrig) {
    selection += " && tis";
  } else if (3 == itrig) {
    selection += " && hlt1 && tos";
  }

  if (string::npos != fSample.find("psi")) {
    selection += " && psimaxdoca<0.5 && mpsi>2.9 && mpsi<3.3 && psipt>6.9";
  }

  return selection;
}


// ----------------------------------------------------------------------
void plotWork::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotWork::makeAll(string what) {

  if (what == "work" || string::npos != what.find("runtis")) {
    runTisEfficiency("bupsikData");
    runTisEfficiency("bspsiphiData");
    runTisEfficiency("bdpsikstarData");
    runTisEfficiency("bmmData");

    runTisEfficiency("bupsikSingleMuon");
    runTisEfficiency("bspsiphiSingleMuon");
    runTisEfficiency("bdpsikstarSingleMuon");
    runTisEfficiency("bmmSingleMuon");
  }

  if (what == "work" || string::npos != what.find("ups2")) {
    ups2("", "");
  }

  if (what == "work" || string::npos != what.find("ups1")) {
    ups1("", "");
  }

  if (what == "work" || string::npos != what.find("plottis")) {
    plotTisEfficiency("all");
  }

  if (what == "all" || string::npos != what.find("effvar")) {
    efficiencyVariable("all", "hlt", 10, 0, 0, 0, "bupsikMc");
    efficiencyVariable("all", "reftrg", 10, 0, 0, 0, "bupsikMc");

    efficiencyVariable("all", "hlt", 11, 0, 0, 0, "bupsikMc");
    efficiencyVariable("all", "reftrg", 11, 0, 0, 0, "bupsikMc");

    efficiencyVariable("all", "hlt", 12, 0, 0, 0, "bupsikMc");

  }

  if (what == "work" || string::npos != what.find("fitstudies")) {
    fS = -1.;
    int ndata(-1), nmc(-1);
    MKKLO = 0.9;
    MKKHI = 1.3;
    DR    = 99.;
    PTK1  = 0.0;
    PTK2  = 0.0;

    fitStudies("bspsiphiData", Form("norm"), ndata);
    if (fS > -1.) return;
    fitStudies("bspsiphiMc", Form("norm"), nmc);

    for (int i = 0; i < 6; ++i) {
      MKKLO = 1.015 - i*0.005;
      MKKHI = 1.025 + i*0.005;
      fitStudies("bspsiphiData", Form("mkk%d", i), ndata);
      fitStudies("bspsiphiMc", Form("mkk%d", i), nmc);
    }

    MKKLO = 1.01;
    MKKHI = 1.03;
    DR    = 2.0;
    PTK1  = 0.0;
    PTK2  = 0.0;
    for (int i = 0; i < 8; ++i) {
      DR = 0.5 - i*0.05;
      fitStudies("bspsiphiData", Form("dr%d", i), ndata);
      fitStudies("bspsiphiMc", Form("dr%d", i), nmc);
    }


    MKKLO = 1.01;
    MKKHI = 1.03;
    DR    = 99.0;
    PTK1  = 0.0;
    PTK2  = 0.0;
    for (int i = 0; i < 8; ++i) {
      PTK1 = 0.6 + i*0.1;
      fitStudies("bspsiphiData", Form("ptk1%d", i), ndata);
      fitStudies("bspsiphiMc", Form("ptk1%d", i), nmc);
    }

    MKKLO = 1.01;
    MKKHI = 1.03;
    DR    = 99.0;
    PTK1  = 0.0;
    PTK2  = 0.0;
    for (int i = 0; i < 8; ++i) {
      PTK2 = 0.6 + i*0.1;
      fitStudies("bspsiphiData", Form("ptk2%d", i), ndata);
      fitStudies("bspsiphiMc", Form("ptk2%d", i), nmc);
    }
  }

  if (what == "relval") {
    genSummary("bsmmrelval", "candAnaMuMu");
    genSummary("bdmmrelval", "candAnaMuMu");
    genSummary("bupsikrelval", "candAnaBu2JpsiK");
    genSummary("bspsiphirelval", "candAnaBs2JpsiPhi");
  }

  if (what == "all" || what == "l1seeds") {
    plotL1Seeds("bupsikData");
    plotL1Seeds("bmmData");
  }

}



// ----------------------------------------------------------------------
void plotWork::bookHist(string dsname) {

  fpHnorm = new TH1D(Form("h_%s_%s", "norm", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);
  fpHpass = new TH1D(Form("h_%s_%s", "pass", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);

}



// ----------------------------------------------------------------------
void plotWork::genSummary(std::string dsname, std::string dir) {
  TH1D *pt    = new TH1D("pt", "pt", 50, 0., 50.0); setTitles(pt, fVarToTex["pt"].c_str(), "", 0.07, 0.9, 1.1, 0.06);
  TH1D *eta   = new TH1D("eta", "eta", 30, 0., 9.0); setTitles(eta, fVarToTex["eta"].c_str(), "", 0.07, 0.9, 1.1, 0.06);
  TH1D *tpt    = new TH1D("tpt", "pt (HLT)", 50, 0., 50.0); setFilledHist(tpt, kBlue, kYellow, 1000);
  TH1D *teta   = new TH1D("teta", "eta (HLT)", 30, 0., 9.0); setFilledHist(teta, kBlue, kYellow, 1000);
  TH1D *meta  = new TH1D("meta", "muon eta", 30, 0., 9.0); setHist(meta, kBlue); setTitles(meta, fVarToTex["meta"].c_str(), "", 0.07, 0.9, 1.1, 0.06);
  TH1D *tmeta  = new TH1D("tmeta", "muon eta (HLT)", 30, 0., 9.0); setFilledHist(tmeta, kBlue, kYellow, 1000);
  TH1D *keta  = new TH1D("keta", "kaon eta", 30, 0., 9.0); setHist(keta, kRed);
  TH1D *mpt   = new TH1D("mpt", "muon pt", 50, 0., 10.0); setHist(mpt, kBlue); setTitles(mpt, fVarToTex["mpt"].c_str(), "", 0.07, 0.9, 1.1, 0.06);
  TH1D *tmpt   = new TH1D("tmpt", "muon pt (HLT)", 50, 0., 10.0); setFilledHist(tmpt, kBlue, kYellow, 1000);
  TH1D *kpt   = new TH1D("kpt", "kaon pt", 50, 0., 10.0); setHist(kpt, kRed);
  TH1D *tau   = new TH1D("tau", "tau", 100, 0., 20.e-12); setTitles(tau, "#tau [ps]", "", 0.07, 0.9, 1.1, 0.06);

  TTree *T = getTree(dsname, dir, "effTree");
  T->Draw("gtau>>tau");
  T->Draw("gpt>>pt");
  T->Draw("TMath::Abs(geta)>>eta");

  T->Draw("g1pt>>mpt");
  T->Draw("g2pt>>mpt");
  T->Draw("TMath::Abs(g1eta)>>meta");
  T->Draw("TMath::Abs(g2eta)>>meta");

  if (string::npos == dsname.find("Bg")) {
    T->Draw("gpt>>tpt", "hlt");
    T->Draw("g1pt>>tmpt", "hlt");
    T->Draw("g2pt>>tmpt", "hlt");

    T->Draw("TMath::Abs(geta)>>teta", "hlt");
    T->Draw("TMath::Abs(g1eta)>>tmeta", "hlt");
    T->Draw("TMath::Abs(g2eta)>>tmeta", "hlt");
  }

  bool addKaon(false);
  if (string::npos != dsname.find("bupsik")) {
    T->Draw("TMath::Abs(g3eta)>>keta");
    T->Draw("g3pt>>kpt");
    addKaon = true;
  }
  if (string::npos != dsname.find("bspsiphi")) {
    T->Draw("TMath::Abs(g3eta)>>keta");
    T->Draw("TMath::Abs(g4eta)>>keta");
    T->Draw("g3pt>>kpt");
    T->Draw("g4pt>>kpt");
    addKaon = true;
  }
  if (string::npos != dsname.find("bdpsikstar")) {
    T->Draw("TMath::Abs(g3eta)>>keta");
    T->Draw("TMath::Abs(g4eta)>>keta");
    T->Draw("g3pt>>kpt");
    T->Draw("g4pt>>kpt");
    addKaon = true;
  }

  tl->SetTextSize(0.05);
  makeCanvas(1);
  int ncol(4);
  c1->Divide(5,1);
  ncol = 5;

  tl->SetTextSize(0.07);

  c1->cd(1);
  gPad->SetTopMargin(0.12);
  pt->Draw();
  tpt->Draw("same");
  tl->DrawLatexNDC(0.15, 0.93, Form("%s ", dsname.c_str()));
  pt->Draw("sameaxis");

  newLegend(0.5, 0.7, 0.8, 0.85);
  legg->SetTextSize(0.07);
  legg->AddEntry(pt,  "B", "l");
  legg->AddEntry(tpt, "HLT passed", "f");
  legg->Draw();


  c1->cd(2);
  gPad->SetTopMargin(0.12);
  eta->SetMinimum(0.);
  eta->Draw();
  teta->Draw("same");
  eta->Draw("sameaxis");
  tl->SetTextSize(0.07);
  tl->DrawLatexNDC(0.75, 0.75, "B");
  tl->DrawLatexNDC(0.15, 0.93, Form("Events: %d", T->GetEntries()));

  c1->cd(3);
  gPad->SetTopMargin(0.12);
  mpt->Draw();
  tmpt->Draw("same");
  mpt->Draw("axissame");
  if (addKaon) kpt->Draw("same");
  tl->SetTextColor(kBlue);
  tl->DrawLatexNDC(0.65, 0.75, "muon(s)");
  tl->SetTextColor(kRed);
  if (addKaon)  tl->DrawLatexNDC(0.65, 0.65, "hadron(s)");
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.15, 0.93, Form("#varepsilon"));
  tl->DrawLatexNDC(0.20, 0.93, Form("= %4.3f ", teta->GetSumOfWeights()/eta->GetSumOfWeights()));

  c1->cd(4);
  gPad->SetTopMargin(0.12);
  meta->Draw();
  meta->SetMinimum(0.);
  tmeta->Draw("same");
  if (addKaon) keta->Draw("same");
  meta->Draw("axissame");

  c1->cd(5);
  gPad->SetTopMargin(0.12);
  gPad->SetLogy(1);
  gStyle->SetOptFit(0);
  tau->Fit("expo", "lr", "", 0., 15.e-12);
  TF1 *f = (TF1*)tau->GetFunction("expo");
  double chi2 = f->GetChisquare();
  int    ndf  = f->GetNDF();
  double t    = -1./f->GetParameter(1);
  double tE   = -t*f->GetParError(1)/f->GetParameter(1);
  t  *= 1.e12;
  tE *= 1.e12;

  tl->DrawLatexNDC(0.15, 0.93, Form("#tau"));
  tl->DrawLatexNDC(0.20, 0.93, Form("= (%5.4f #pm %5.4f)ps", t, tE));

  c1->SaveAs(Form("%s/genSummary-%s.pdf", fDirectory.c_str(), dsname.c_str()));

}

// ----------------------------------------------------------------------
void plotWork::plotL1Seeds(std::string dsname) {
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
void plotWork::plotTisEfficiency(string dsname) {

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
void plotWork::runTisEfficiency(string dsname) {

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
void plotWork::refTrgEfficiency(string selection, string dsname) {

  if (string::npos != dsname.find("bupsik")) fMode = BU2JPSIKP;
  if (string::npos != dsname.find("bspsiphi")) fMode = BS2JPSIPHI;
  if (string::npos != dsname.find("bdpsikstar")) fMode = BD2JPSIKSTAR;

  zone(2,2);

  string dir = "candAnaBu2JpsiK";
  fSample = dsname;

  TTree *t = getTree(fSample, dir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  }

  TH1D *h1 = new TH1D(Form("h_mc"), Form("h_mc"), 40, 4.8, 6.0);
  TH1D *h2 = new TH1D(Form("h_mc_hlt"), Form("h_mc_hlt"), 40, 4.8, 6.0);

  double nNorm(0.), nPass(0.), effMc(0.), effRt(0.), effMcE(0.), effRtE
    (0.);
  double mBp(5.28), sBp(0.05);

  string fitopt("lm");

  // -- basic HLT efficiency derived from MC
  string tselection = selection;
  t->Draw("m >> h_mc", tselection.c_str());
  tselection = selection + " && hlt";
  t->Draw("m >> h_mc_hlt", tselection.c_str());

  fIF->limitPar(1, 5.1, 5.5);
  TF1 *f1 = fIF->pol1gauss2c(h1, mBp, sBp);
  c0->cd(1);
  h1->Fit(f1, fitopt.c_str());

  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  nNorm = f1->Integral(5.1, 5.4)/h1->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "MC");
  savePad("h_mc.pdf");
  delete f1;

  c0->cd(2);
  fIF->limitPar(1, 5.1, 5.5);
  f1 = fIF->pol1gauss2c(h2, mBp, sBp);
  h2->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  nPass = f1->Integral(5.1, 5.4)/h2->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "MC && HLT");
  savePad("h_mc_hlt.pdf");
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

  fIF->limitPar(1, 5.1, 5.5);
  f1 = fIF->pol1gauss2c(h3, mBp, sBp);
  h3->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  nNorm = f1->Integral(5.1, 5.4)/h3->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "ref trigger");
  savePad("h_rt.pdf");
  delete f1;

  c0->cd(4);
  fIF->limitPar(1, 5.1, 5.5);
  f1 = fIF->pol1gauss2c(h4, mBp, sBp);
  h4->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  nPass = f1->Integral(5.1, 5.4)/h4->GetBinWidth(1);
  effRt = nPass/nNorm;
  effRtE = dEff(static_cast<int>(nPass), static_cast<int>(nNorm));
  tl->DrawLatexNDC(0.2, 0.96, "ref trigger && HLT");
  savePad("h_rt_hlt.pdf");
  delete f1;
  cout << "==> efficiency = " << nPass << "/" << nNorm << " = " << nPass/nNorm << endl;
  cout << "==> cuts: " << selection << endl;
  cout << "==> efficiencies: MC = " << Form("%4.2f +/- %4.2f", effMc, effMcE)
       << "; ref trigger = " << Form("%4.2f +/- %4.2f", effRt, effRtE)
       << endl;

}


// ----------------------------------------------------------------------
void plotWork::efficiencyVariable(string var, string effvar, int iselection, int nbin, double xmin, double xmax, string dsname) {

  fSample = dsname;
  cout << "==> plotWork::efficiencyVariable> sample = " << fSample << endl;;
  cout << "==> plotWork::efficiencyVariable> orig selection =  " << selectionString(iselection, 0) << endl;
  string selection = removeVarFromSelection(var, selectionString(iselection, 0));
  cout << "==> plotWork::efficiencyVariable> adap selection = " << selection << endl;

  if (var == "all") {
    var = "m1pt";     efficiencyVariable(var, effvar, iselection, 40,   0.,  40., dsname);
    var = "m2pt";     efficiencyVariable(var, effvar, iselection, 20,   0.,  20., dsname);
    var = "m1eta";    efficiencyVariable(var, effvar, iselection, 20, -2.0,  2.0, dsname);
    var = "m2eta";    efficiencyVariable(var, effvar, iselection, 20, -2.0,  2.0, dsname);
    var = "pt";       efficiencyVariable(var, effvar, iselection, 20,   0.,  40., dsname);
    var = "eta";      efficiencyVariable(var, effvar, iselection, 20, -2.0,  2.0, dsname);
    var = "fls3d";    efficiencyVariable(var, effvar, iselection, 20,   0., 120., dsname);
    var = "chi2dof";  efficiencyVariable(var, effvar, iselection, 25,   0.,  5.0, dsname);
    var = "iso";      efficiencyVariable(var, effvar, iselection, 51,   0., 1.02, dsname);
    var = "m1iso";    efficiencyVariable(var, effvar, iselection, 51,   0., 1.02, dsname);
    var = "m2iso";    efficiencyVariable(var, effvar, iselection, 51,   0., 1.02, dsname);
    var = "closetrk"; efficiencyVariable(var, effvar, iselection,  6,   0.,   6., dsname);
    var = "docatrk";  efficiencyVariable(var, effvar, iselection, 20,   0.,  0.2, dsname);
    var = "pvip";     efficiencyVariable(var, effvar, iselection, 20,   0., 0.02, dsname);
    var = "pvips";    efficiencyVariable(var, effvar, iselection, 20,   0.,  4.0, dsname);
    var = "maxdoca";  efficiencyVariable(var, effvar, iselection, 20,   0., 0.06, dsname);
    var = "alpha";    efficiencyVariable(var, effvar, iselection, 20,   0.,  0.1, dsname);
    return;
  }

  // -- dump histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;

  string dir("candAnaMuMu");
  if (string::npos != dsname.find("bupsik")) {
    dir = "candAnaBu2JpsiK";
  }

  TTree *t = getTree(fSample, dir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  } else {
    cout << "tree for sample = " << fSample << " found" << endl;
  }

  gStyle->SetHatchesLineWidth(2);

  string normName = Form("effVar_%s_%s_%s_norm", fSample.c_str(), var.c_str(), effvar.c_str());
  string passName = Form("effVar_%s_%s_%s_pass", fSample.c_str(), var.c_str(), effvar.c_str());
  string effName  = Form("effVar_%s_%s_%s_eff", fSample.c_str(), var.c_str(), effvar.c_str());
  TH1D *h1 = new TH1D(normName.c_str(), normName.c_str(), nbin, xmin, xmax);
  h1->Sumw2();
  setTitles(h1, fVarToTex[var].c_str(), "");
  setFilledHist(h1, kBlue, kYellow, 1000, 2);
  TH1D *h2 = new TH1D(passName.c_str(), passName.c_str(), nbin, xmin, xmax);
  h2->Sumw2();
  setFilledHist(h2, kBlue, kBlue, 3354, 2);

  // -- basic HLT efficiency derived from MC
  string tselection = selection;
  t->Draw(Form("%s >> %s", var.c_str(), normName.c_str()), tselection.c_str());
  cout << "==> " << var << " SEL histogram contents =        " << h1->Integral(1, h1->GetNbinsX()+1) << endl;
  tselection = selection + " && " + effvar;
  t->Draw(Form("%s >> %s", var.c_str(), passName.c_str()), tselection.c_str());
  cout << "==> " << var << " SEL && " << effvar << " histogram contents = " << h2->Integral(1, h1->GetNbinsX()+1) << endl;

  TH1D *h3 = (TH1D*)(h1->Clone(effName.c_str())); h3->Reset();
  setHist(h3);
  h3->Divide(h2, h1, 1., 1., "b");
  setTitles(h3, fVarToTex[var].c_str(), "Efficiency");

  zone(1,2);
  h1->Draw("hist");
  h2->Draw("histsame");
  h1->Draw("axissame");

  c0->cd(2);
  h3->SetMinimum(0.);
  h3->SetMaximum(1.);
  h3->Draw("e");

  c0->cd();
  savePad(Form("trgEfficiency_%s_%s_%s_sel%d.pdf", fSample.c_str(), var.c_str(), effvar.c_str(), iselection));

  h1->Write();
  h2->Write();
  h3->Write();
  fHistFile->Close();

}



// ----------------------------------------------------------------------
void plotWork::bmm5Trigger(std::string cuts, std::string texname) {


  if (cuts == "all") {
    bmm5Trigger("gmugmid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>4", "m2pt4eta14");
    bmm5Trigger("gmugmid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>5", "m2pt5eta14");
    bmm5Trigger("gmugmid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>6", "m2pt6eta14");
    bmm5Trigger("gmugmid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>7", "m2pt7eta14");
    bmm5Trigger("gmugmid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>4&&fls3d>10&&alpha<0.05&&iso>0.8", "m2pt4eta14alphafls3diso");
    bmm5Trigger("gmugmid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>6&&fls3d>10&&alpha<0.05&&iso>0.8", "m2pt6eta14alphafls3diso");

    bmm5Trigger("gmuid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>4", "gmuidm2pt4eta14");
    bmm5Trigger("gmuid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>5", "gmuidm2pt5eta14");
    bmm5Trigger("gmuid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>6", "gmuidm2pt6eta14");
    bmm5Trigger("gmuid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>7", "gmuidm2pt7eta14");
    bmm5Trigger("gmuid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>4&&fls3d>10&&alpha<0.05&&iso>0.8", "gmuidm2pt4eta14alphafls3diso");
    bmm5Trigger("gmuid&&TMath::Abs(m1eta)<1.4&&TMath::Abs(m2eta)<1.4&&m2pt>6&&fls3d>10&&alpha<0.05&&iso>0.8", "gmuidm2pt6eta14alphafls3diso");


    return;
  }

  std::ofstream TEX;
  TEX.open("bmm5Trigger.tex", ios::app);

  vector<string> seeds;
  if (0) {
    seeds.push_back("L1_DoubleMu7_OS");
    seeds.push_back("L1_DoubleMu_13_6");
    seeds.push_back("L1_DoubleMu0er1p6_dEta_Max1p8");
    seeds.push_back("L1_DoubleMu0er1p6_dEta_Max1p8_OS");
    seeds.push_back("L1_DoubleMu0er1p2_dEta_Max1p8_OS");
    seeds.push_back("L1_DoubleMu0er1p2_OS_MASS0to10");
    seeds.push_back("L1_DoubleMu0er1p4_dR0to1p8_OS");
    seeds.push_back("L1_DoubleMu6_OS_MASS0to10");

    seeds.push_back("L1_DoubleMu6_OS_dR0to1");
    seeds.push_back("L1_DoubleMu0er1p4_OS_MASS0to10");
    seeds.push_back("L1_DoubleMu0er1p4_dR_0to1p4_OS");
    seeds.push_back("L1_DoubleMu0er1p4_dR_0to1p2_OS");
    seeds.push_back("L1_DoubleMu0er1p4_dR_0to1p0_OS");

    seeds.push_back("L1_DoubleMu0");
    seeds.push_back("L1_DoubleMu0_bph1");
    seeds.push_back("L1_DoubleMu0_bph2");
    seeds.push_back("L1_DoubleMu0_bph3");
    seeds.push_back("L1_DoubleMu0_bph4");
    seeds.push_back("L1_DoubleMu0_bph5");
    seeds.push_back("L1_DoubleMu0_bph6");
    seeds.push_back("L1_DoubleMu0_bph7");
    seeds.push_back("L1_DoubleMu0_bph8");
  }

  seeds.push_back("L1_DoubleMu0");
  seeds.push_back("L1_DoubleMu0_bph1");

  TFile *f(0);
  TTree *t(0);
  for (unsigned int i = 0; i < seeds.size(); ++i) {
    string fname = Form("T1-rerunHLT-BsToMuMu-%s.bmmReader.mix-Bs2MuMu.root", seeds[i].c_str());
    f = TFile::Open(fname.c_str());
    t = (TTree*)(f->Get("candAnaMuMu/events"));
    if (0 == t) {
      cout << "tree not found in " << fname;
      continue;
    }
    double norm = t->Draw("m", cuts.c_str(), "goff");
    string pcuts = cuts + "&&hlt1&&tos";
    double pass = t->Draw("m", pcuts.c_str(), "goff");
    double effV = pass/norm;
    double effE = dEff(static_cast<int>(pass), static_cast<int>(norm));
    cout << "seed: " << seeds[i] << " eff: " << effV << " +/- " <<  effE << " " << texname << endl;
    string useed = seeds[i];
    replaceAll(useed, "_", "\\%");
    replaceAll(useed, "%", "_");
    TEX << Form("\\vdef{%s:%s:seed}   {\\tt %s } ",  seeds[i].c_str(), texname.c_str(), useed.c_str()) << endl;
    TEX << formatTex(effV, Form("%s:%s:val", seeds[i].c_str(), texname.c_str()), 3) << endl;
    TEX << formatTex(effE, Form("%s:%s:err", seeds[i].c_str(), texname.c_str()), 3) << endl;
    f->Close();
  }

  TEX.close();

}



// ----------------------------------------------------------------------
void plotWork::fitStudies(string dsname, string tag, int nevt, int nstart) {
  fSample = dsname;
  string dir = "candAnaBs2JpsiPhi";
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

  TH1D *h1(0), *h2(0);
  // -- check for data histogram
  fHistFile = TFile::Open(fHistFileName.c_str(), "");
  if (fHistFile) {
    h1 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", dsname.c_str(), 0, tag.c_str()));
    fHistFile->Close();
  }
  if (!h1) {
    cout << "fHistFile: " << fHistFileName;
    fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
    cout << " opened " << endl;

    fHma.clear();
    fHmc.clear();
    string cuts = Form("%4.3f<mkk<%4.3f, dr<%3.2f, ptk1>%3.2f, ptk2>%3.2f", MKKLO, MKKHI, DR, PTK1, PTK2);
    for (int i = 0; i < fNchan; ++i) {
      h1 = new TH1D(Form("hma_%s_chan%d_%s", dsname.c_str(), i, tag.c_str()),
		    Form("m %s_chan%d %s cuts:%s", dsname.c_str(), i, tag.c_str(), cuts.c_str()),
		    90, 5.0, 5.9);
      fHma.push_back(h1);
      h1 = new TH1D(Form("hmc_%s_chan%d_%s", dsname.c_str(), i, tag.c_str()),
		    Form("m (constrained) %s_chan%d cuts:%s", dsname.c_str(), i, cuts.c_str()),
		    90, 5.0, 5.9);
      fHmc.push_back(h1);
    }

    TTree *t = getTree(fSample, dir);
    if (0 == t) {
      cout << "tree for sample = " << fSample << " not found" << endl;
      return;
    }

    setupTree(t, fSample);
    fCds = fDS[fSample];
    loopOverTree(t, 3, nevt, nstart);

    fHistFile->Write();
    fHistFile->Close();
    fS = -1.;
    return;
  }

  // --------------------------------
  // -- Fit and analyze the histograms
  // --------------------------------
  fHistFile = TFile::Open(fHistFileName.c_str(), "");

  zone(2, 4);
  gStyle->SetOptFit(0);
  string mcname = dsname;
  replaceAll(mcname, "Charmonium", "Mc");
  replaceAll(mcname, "Data", "Mc");
  tag = "norm";
  string cuts;
  vector<double> normMc, normDa;
  for (int i = 0; i < fNchan; ++i) {
    h1 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", dsname.c_str(), i, tag.c_str()));
    cuts = h1->GetTitle(); cuts = cuts.substr(cuts.find("cuts:") + 5);
    h1->SetMinimum(0.);
    h2 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", mcname.c_str(), i, tag.c_str()));
    h2->SetMinimum(0.);
    fIF->fLo = 5.0;
    fIF->fHi = 5.9;
    c0->cd(2*i+1);
    fitStudiesFit0(h1, i);
    normDa.push_back(fS);
    tl->DrawLatexNDC(0.60, 0.80, Form("N   = %5.1f", fN));
    tl->DrawLatexNDC(0.60, 0.70, Form("W/S = %5.1f / %5.1f", fW, fS));
    tl->DrawLatexNDC(0.67, 0.60, Form("= %4.3f", fW/fS));
    tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
    fTEX << Form("\\vdef{%s:fitstudies_chan%d_norm:var} {%s}", fSuffix.c_str(), i, tag.c_str()) << endl;
    fTEX << Form("\\vdef{%s:fitstudies_chan%d_norm:cuts} {%s}", fSuffix.c_str(), i, cuts.c_str()) << endl;
    fTEX << formatTex(fS, Form("%s:fitstudies_chan%d_norm:daS", fSuffix.c_str(), i), 1)	 << endl;
    fTEX << formatTex(fSsigma, Form("%s:fitstudies_chan%d_norm:daSsigma", fSuffix.c_str(), i), 3)	 << endl;
    fTEX << formatTex(fW/fS, Form("%s:fitstudies_chan%d_norm:W/S", fSuffix.c_str(), i), 2) << endl;
    fTEX << formatTex(fS/fB,Form("%s:fitstudies_chan%d_norm:S/B", fSuffix.c_str(), i), 2)	 << endl;
    c0->cd(2*i+2);
    fitStudiesFit0(h2, i);
    normMc.push_back(fEntries);
    tl->DrawLatexNDC(0.60, 0.80, Form("N   = %5.1f", fN));
    tl->DrawLatexNDC(0.60, 0.70, Form("W/S = %5.1f / %5.1f", fW, fS));
    tl->DrawLatexNDC(0.67, 0.60, Form("= %4.3f", fW/fS));
    tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
    tl->DrawLatexNDC(0.62, 0.40, Form("RMS = %4.3f", fSRMS));
    fTEX << formatTex(fS,Form("%s:fitstudies_chan%d_norm:mcS", fSuffix.c_str(), i), 2)	 << endl;
    fTEX << formatTex(fSsigma,Form("%s:fitstudies_chan%d_norm:mcSsigma", fSuffix.c_str(), i), 3)	 << endl;
    fTEX << formatTex(fSRMS,Form("%s:fitstudies_chan%d_norm:mcSrms", fSuffix.c_str(), i), 3)	 << endl;
 }
  c0->cd(2*(fNchan-1)+1);
  tl->DrawLatexNDC(0.80, 0.40, tag.c_str());
  savePad(Form("fitStudies_normalization.pdf"), c0);

  tag = "mkk";
  for (int imod = 0; imod < 6; ++imod) {
    string tg = Form("%s%d", tag.c_str(), imod);
    for (int i = 0; i < fNchan; ++i) {
      h1 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", dsname.c_str(), i, tg.c_str()));
      cuts = h1->GetTitle(); cuts = cuts.substr(cuts.find("cuts:") + 5);
      h1->SetMinimum(0.);
      h2 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", mcname.c_str(), i, tg.c_str()));
      h2->SetMinimum(0.);
      fIF->fLo = 5.0;
      fIF->fHi = 5.9;
      c0->cd(2*i+1);
      fitStudiesFit0(h1, i);
      tl->DrawLatexNDC(0.60, 0.80, Form("N   = %5.1f", fN));
      tl->DrawLatexNDC(0.60, 0.70, Form("W/S = %5.1f / %5.1f", fW, fS));
      tl->DrawLatexNDC(0.67, 0.60, Form("= %4.3f", fW/fS));
      tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
      fTEX << Form("\\vdef{%s:fitstudies_chan%d_%s:var} {%s}", fSuffix.c_str(), i, tg.c_str(), tg.c_str()) << endl;
      fTEX << Form("\\vdef{%s:fitstudies_chan%d_%s:cuts} {%s}", fSuffix.c_str(), i, tg.c_str(), cuts.c_str()) << endl;
      fTEX << formatTex(fS, Form("%s:fitstudies_chan%d_%s:daS", fSuffix.c_str(), i, tg.c_str()), 1) << endl;
      fTEX << formatTex(fSsigma, Form("%s:fitstudies_chan%d_%s:daSsigma", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fW/fS, Form("%s:fitstudies_chan%d_%s:W/S", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      fTEX << formatTex(fS/fB, Form("%s:fitstudies_chan%d_%s:S/B", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      c0->cd(2*i+2);
      fitStudiesFit0(h2, i);
      tl->DrawLatexNDC(0.60, 0.70, Form("#varepsilon = %4.3f", fEntries/normMc[i]));
      tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
      tl->DrawLatexNDC(0.62, 0.40, Form("RMS = %4.3f", fSRMS));
      fTEX << formatTex(fS,Form("%s:fitstudies_chan%d_%s:mcS", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      fTEX << formatTex(fSsigma,Form("%s:fitstudies_chan%d_%s:mcSsigma", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fSRMS,Form("%s:fitstudies_chan%d_%s:mcSrms", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fEntries/normMc[i],Form("%s:fitstudies_chan%d_%s:eff", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
    }
    c0->cd(2*(fNchan-1)+1);
    tl->DrawLatexNDC(0.80, 0.40, Form("%s", tg.c_str()));
    savePad(Form("fitStudies_%s.pdf", tg.c_str()), c0);
  }

  tag = "dr";
  for (int imod = 0; imod < 8; ++imod) {
    string tg = Form("%s%d", tag.c_str(), imod);
    for (int i = 0; i < fNchan; ++i) {
      h1 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", dsname.c_str(), i, tg.c_str()));
      cuts = h1->GetTitle(); cuts = cuts.substr(cuts.find("cuts:") + 5);
      h1->SetMinimum(0.);
      h2 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", mcname.c_str(), i, tg.c_str()));
      h2->SetMinimum(0.);
      fIF->fLo = 5.0;
      fIF->fHi = 5.9;
      c0->cd(2*i+1);
      fitStudiesFit0(h1, i);
      tl->DrawLatexNDC(0.60, 0.80, Form("N   = %5.1f", fN));
      tl->DrawLatexNDC(0.60, 0.70, Form("W/S = %5.1f / %5.1f", fW, fS));
      tl->DrawLatexNDC(0.67, 0.60, Form("= %4.3f", fW/fS));
      tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
      fTEX << Form("\\vdef{%s:fitstudies_chan%d_%s:var} {%s}", fSuffix.c_str(), i, tg.c_str(), tg.c_str()) << endl;
      fTEX << Form("\\vdef{%s:fitstudies_chan%d_%s:cuts} {%s}", fSuffix.c_str(), i, tg.c_str(), cuts.c_str()) << endl;
      fTEX << formatTex(fS, Form("%s:fitstudies_chan%d_%s:daS", fSuffix.c_str(), i, tg.c_str()), 1) << endl;
      fTEX << formatTex(fSsigma, Form("%s:fitstudies_chan%d_%s:daSsigma", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fW/fS, Form("%s:fitstudies_chan%d_%s:W/S", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      fTEX << formatTex(fS/fB, Form("%s:fitstudies_chan%d_%s:S/B", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      c0->cd(2*i+2);
      fitStudiesFit0(h2, i);
      tl->DrawLatexNDC(0.60, 0.70, Form("#varepsilon = %4.3f", fEntries/normMc[i]));
      tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
      tl->DrawLatexNDC(0.62, 0.40, Form("RMS = %4.3f", fSRMS));
      fTEX << formatTex(fS,Form("%s:fitstudies_chan%d_%s:mcS", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      fTEX << formatTex(fSsigma, Form("%s:fitstudies_chan%d_%s:mcSsigma", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fSRMS,Form("%s:fitstudies_chan%d_%s:mcSrms", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fEntries/normMc[i],Form("%s:fitstudies_chan%d_%s:eff", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
    }
    c0->cd(2*(fNchan-1)+1);
    tl->DrawLatexNDC(0.80, 0.40, Form("%s", tg.c_str()));
    savePad(Form("fitStudies_%s.pdf", tg.c_str()), c0);
  }


  tag = "ptk1";
  for (int imod = 0; imod < 8; ++imod) {
    string tg = Form("%s%d", tag.c_str(), imod);
    for (int i = 0; i < fNchan; ++i) {
      h1 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", dsname.c_str(), i, tg.c_str()));
      cuts = h1->GetTitle(); cuts = cuts.substr(cuts.find("cuts:") + 5);
      h1->SetMinimum(0.);
      h2 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", mcname.c_str(), i, tg.c_str()));
      h2->SetMinimum(0.);
      fIF->fLo = 5.0;
      fIF->fHi = 5.9;
      c0->cd(2*i+1);
      fitStudiesFit0(h1, i);
      tl->DrawLatexNDC(0.60, 0.80, Form("N   = %5.1f", fN));
      tl->DrawLatexNDC(0.60, 0.70, Form("W/S = %5.1f / %5.1f", fW, fS));
      tl->DrawLatexNDC(0.67, 0.60, Form("= %4.3f", fW/fS));
      tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
      fTEX << Form("\\vdef{%s:fitstudies_chan%d_%s:var} {%s}", fSuffix.c_str(), i, tg.c_str(), tg.c_str()) << endl;
      fTEX << Form("\\vdef{%s:fitstudies_chan%d_%s:cuts} {%s}", fSuffix.c_str(), i, tg.c_str(), cuts.c_str()) << endl;
      fTEX << formatTex(fS, Form("%s:fitstudies_chan%d_%s:daS", fSuffix.c_str(), i, tg.c_str()), 1) << endl;
      fTEX << formatTex(fSsigma, Form("%s:fitstudies_chan%d_%s:daSsigma", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fW/fS, Form("%s:fitstudies_chan%d_%s:W/S", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      fTEX << formatTex(fS/fB, Form("%s:fitstudies_chan%d_%s:S/B", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      c0->cd(2*i+2);
      fitStudiesFit0(h2, i);
      tl->DrawLatexNDC(0.60, 0.70, Form("#varepsilon = %4.3f", fEntries/normMc[i]));
      tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
      tl->DrawLatexNDC(0.62, 0.40, Form("RMS = %4.3f", fSRMS));
      fTEX << formatTex(fS,Form("%s:fitstudies_chan%d_%s:mcS", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      fTEX << formatTex(fSsigma, Form("%s:fitstudies_chan%d_%s:mcSsigma", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fSRMS,Form("%s:fitstudies_chan%d_%s:mcSrms", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fEntries/normMc[i],Form("%s:fitstudies_chan%d_%s:eff", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
    }
    c0->cd(2*(fNchan-1)+1);
    tl->DrawLatexNDC(0.80, 0.40, Form("%s", tg.c_str()));
    savePad(Form("fitStudies_%s.pdf", tg.c_str()), c0);
  }

  tag = "ptk2";
  for (int imod = 0; imod < 8; ++imod) {
    string tg = Form("%s%d", tag.c_str(), imod);
    for (int i = 0; i < fNchan; ++i) {
      h1 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", dsname.c_str(), i, tg.c_str()));
      cuts = h1->GetTitle(); cuts = cuts.substr(cuts.find("cuts:") + 5);
      h1->SetMinimum(0.);
      h2 = (TH1D*)fHistFile->Get(Form("hma_%s_chan%d_%s", mcname.c_str(), i, tg.c_str()));
      h2->SetMinimum(0.);
      fIF->fLo = 5.0;
      fIF->fHi = 5.9;
      c0->cd(2*i+1);
      fitStudiesFit0(h1, i);
      tl->DrawLatexNDC(0.60, 0.80, Form("N   = %5.1f", fN));
      tl->DrawLatexNDC(0.60, 0.70, Form("W/S = %5.1f / %5.1f", fW, fS));
      tl->DrawLatexNDC(0.67, 0.60, Form("= %4.3f", fW/fS));
      tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
      fTEX << Form("\\vdef{%s:fitstudies_chan%d_%s:var} {%s}", fSuffix.c_str(), i, tg.c_str(), tg.c_str()) << endl;
      fTEX << Form("\\vdef{%s:fitstudies_chan%d_%s:cuts} {%s}", fSuffix.c_str(), i, tg.c_str(), cuts.c_str()) << endl;
      fTEX << formatTex(fS, Form("%s:fitstudies_chan%d_%s:daS", fSuffix.c_str(), i, tg.c_str()), 1) << endl;
      fTEX << formatTex(fSsigma, Form("%s:fitstudies_chan%d_%s:daSsigma", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fW/fS, Form("%s:fitstudies_chan%d_%s:W/S", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      fTEX << formatTex(fS/fB, Form("%s:fitstudies_chan%d_%s:S/B", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      c0->cd(2*i+2);
      fitStudiesFit0(h2, i);
      tl->DrawLatexNDC(0.60, 0.70, Form("#varepsilon = %4.3f", fEntries/normMc[i]));
      tl->DrawLatexNDC(0.67, 0.50, Form("#sigma = %4.3f", fSsigma));
      tl->DrawLatexNDC(0.62, 0.40, Form("RMS = %4.3f", fSRMS));
      fTEX << formatTex(fS,Form("%s:fitstudies_chan%d_%s:mcS", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
      fTEX << formatTex(fSsigma, Form("%s:fitstudies_chan%d_%s:mcSsigma", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fSRMS,Form("%s:fitstudies_chan%d_%s:mcSrms", fSuffix.c_str(), i, tg.c_str()), 3) << endl;
      fTEX << formatTex(fEntries/normMc[i],Form("%s:fitstudies_chan%d_%s:eff", fSuffix.c_str(), i, tg.c_str()), 2) << endl;
    }
    c0->cd(2*(fNchan-1)+1);
    tl->DrawLatexNDC(0.80, 0.40, Form("%s", tg.c_str()));
    savePad(Form("fitStudies_%s.pdf", tg.c_str()), c0);
  }


}


// ----------------------------------------------------------------------
void  plotWork::fitStudiesFit0(TH1D *h1, int i) {

  TF1 *f1(0);
  double xmin(5.0), xmax(6.0);
  double mBs(5.37), sBs(0.04), stepBs(5.21);

  TF1 *fb = new TF1(Form("b%s_%d", h1->GetName(), i), "[0] * TMath::Exp([1]*x)", 5.0, 6.0);
  TF1 *fs = new TF1(Form("s%s_%d", h1->GetName(), i), "[0] * TMath::Gaus(x, [1], [2], false)", 5.0, 6.0);
  TF1 *fw = new TF1(Form("w%s_%d", h1->GetName(), i), "[0] * TMath::Gaus(x, [1], [2], false)", 5.0, 6.0);
  TH1D *hr = (TH1D*)h1->Clone("hr"); hr->Reset();
  fSRMS = h1->GetRMS();
  fIF->fLo = xmin;
  fIF->fHi = xmax;
  //  f1 = fIF->pol1gauss2(h1, mBs, sBs, 0.1, 0.1);
  f1 = fIF->expogauss2(h1, mBs, sBs, 0.1, 0.1);
  h1->Fit(f1, "lr", "", xmin, 5.9);
  f1 = (TF1*)f1->Clone(Form("f%d", i));
  fb->SetParameters(f1->GetParameter(6),
		    f1->GetParameter(7));
  fb->SetLineStyle(kDashed);
  fb->Draw("same");
  fs->SetParameters(f1->GetParameter(0),
		    f1->GetParameter(1),
		    f1->GetParameter(2));
  fs->SetLineColor(kBlue);
  fs->SetLineStyle(kDashed);
  fs->Draw("same");
  fw->SetParameters(f1->GetParameter(3)*f1->GetParameter(0),
		    f1->GetParameter(4),
		    f1->GetParameter(5));
  fw->SetLineColor(kMagenta);
  fw->SetLineStyle(kDashed);
  fw->Draw("same");

  xmin = f1->GetParameter(1) - 2.*f1->GetParameter(2);
  xmax = f1->GetParameter(1) + 2.*f1->GetParameter(2);
  double S = fs->Integral(xmin, xmax)/h1->GetBinWidth(1);
  double W = fw->Integral(xmin, xmax)/h1->GetBinWidth(1);
  double B = fb->Integral(xmin, xmax)/h1->GetBinWidth(1);
  if (B < 1.) B = 1.;
  // -- estimate RMS of combined Gaussian
  for (int i = 0; i < S; ++i) {
    hr->Fill(fs->GetRandom());
  }
  // -- if the two Gaussians are close together, count both of them as signal
  fN  = S;
  fNE = S * f1->GetParError(0) / f1->GetParameter(0);;
  if (TMath::Abs(f1->GetParameter(1) - f1->GetParameter(4)) < 0.7*f1->GetParameter(2)) {
    S = S + W;
    for (int i = 0; i < W; ++i) {
      hr->Fill(fw->GetRandom());
    }
  }



  tl->SetTextSize(0.06);
  tl->DrawLatexNDC(0.2, 0.95, Form("S/B = %5.1f / %5.1f", S, B));
  tl->DrawLatexNDC(0.7, 0.95, Form("fit: %3.1f ", f1->GetChisquare()/f1->GetNDF()));
  tl->DrawLatexNDC(0.2, 0.88, Form("S/B = %3.1f", S/B));
  tl->DrawLatexNDC(0.6, 0.88, Form("%4.3f < x < %4.3f ", xmin, xmax));
  fS  = S;
  fSE = S * f1->GetParError(0) / f1->GetParameter(0);
  fW  = W;
  fWE = W * f1->GetParError(3) / f1->GetParameter(3);
  fB  = B;
  fBE = B * f1->GetParError(7) / f1->GetParameter(7);
  fChi2Dof = f1->GetChisquare()/f1->GetNDF();
  fSsigma  = f1->GetParameter(2);
  fSsigma  = hr->GetRMS();
  fEntries = h1->GetSumOfWeights();
}




// ----------------------------------------------------------------------
void plotWork::loopFunction1() {

  if (!fb.tis) return;

  if (!fGoodQ) return;
  if (!fGoodPvAveW8) return;
  if (!fGoodMaxDoca) return;
  if (!fGoodIp) return;
  if (!fGoodIpS) return;

  if (!fGoodLip) return;
  if (!fGoodLipS) return;

  if (!fGoodAlpha) return;
  if (!fGoodChi2) return;

  if (!fGoodCloseTrack) return;
  if (!fGoodIso) return;
  if (!fGoodDocaTrk) return;

  if (TMath::Abs(fb.flsxy) < 3.0) return;
  if (TMath::Abs(fb.fls3d) < 10.0) return;

  if (TMath::Abs(fb.m1eta) > 1.6) return;
  if (TMath::Abs(fb.m2eta) > 1.6) return;

  if (fb.m1pt < 4.0) return;
  if (fb.m2pt < 3.0) return;


  if ((fMode == BU2JPSIKP) || (fMode == BD2JPSIKSTAR) || (fMode == BS2JPSIPHI)) {
    if (TMath::Abs(fb.mpsi) < 2.9) return;
    if (TMath::Abs(fb.mpsi) > 3.3) return;

    if (TMath::Abs(fb.psipt) < 7.0) return;
    if (TMath::Abs(fb.psiprob) < 0.1) return;
  }

  fpHnorm->Fill(fb.m);

  if (!fb.hlt1) return;
  fpHpass->Fill(fb.m);

}


// ----------------------------------------------------------------------
void plotWork::ups1(std::string file1, std::string file2) {

  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;

  TH1D *h = (TH1D*)fHistFile->Get("Hmass0_1");
  if (0 == h) {
    for (int i = 0; i < 13; ++i) {
	fHmass0.push_back(new TH1D(Form("Hmass0_%d", i), Form(" %3.1f < |#it{#eta}_{f}| < %3.1f", i*0.2, (i+1)*0.2), 80, 5.0, 5.8));
	fHmass1.push_back(new TH1D(Form("Hmass1_%d", i), Form(" %3.1f < |#it{#eta}_{f}| < %3.1f", i*0.2, (i+1)*0.2), 80, 5.0, 5.8));
      }

    TFile *f0 = TFile::Open("/scratch/ursl/bmm4/v06/bmm-mc-RunIISpring16DR80-BsToMuMu_BMuonFilter-v06.root");
    TTree *t0 = (TTree*)f0->Get("candAnaMuMu/events");

    setupTree(t0);
    fSetupInt = 0;
    loopOverTree(t0, 2);

    //  TFile *f1 = TFile::Open("/scratch/ursl/bmm4/ups/bmm-mc-ups-bsmm-v04-0000.bmmReader.mix-Bs2MuMu.root");
    TFile *f1 = TFile::Open("/scratch/ursl/bmm4/ups/bmm-mc-ups-bsmmExtended.root");
    cout << "f1 = " << f1 << endl;
    TTree *t1 = (TTree*)f1->Get("candAnaMuMu/events");
    cout << "t1 = " << t1 << endl;

    fSetupInt = 1;
    setupTree(t1);
    loopOverTree(t1, 2);

    fHistFile->cd();
    for (int i = 0; i < 13; ++i) {
      fHmass0[i]->Write();
      fHmass1[i]->Write();
    }
    fHistFile->Close();

    fHmass0.clear();
    fHmass1.clear();
  }

  fHistFile = TFile::Open(fHistFileName.c_str());
  gStyle->SetOptFit(0);

  shrinkPad(0.15, 0.19);
  double eps(0.01);
  TH1D *h1 = new TH1D("mass0", "Run 2", 13, 0., 2.6); h1->Sumw2();
  TH1D *h2 = new TH1D("mass1", "Phase 2", 13, 0.+eps, 2.6+eps); h2->Sumw2();
  TH1D *s1 = new TH1D("rms0", "Run 2", 13, 0., 2.6); h1->Sumw2();
  TH1D *s2 = new TH1D("rms1", "Phase 2", 13, 0.+eps, 2.6+eps); h2->Sumw2();
  TH1D *p1 = new TH1D("peak0", "Run 2", 13, 0., 2.6); p1->Sumw2();
  TH1D *p2 = new TH1D("peak1", "Phase 2", 13, 0.+eps, 2.6+eps); p2->Sumw2();
  TH1D *w1 = new TH1D("sigma0", "Run 2", 13, 0., 2.6); w1->Sumw2();
  TH1D *w2 = new TH1D("sigma1", "Phase 2", 13, 0.+eps, 2.6+eps); w2->Sumw2();
  TH1D *hmass0(0), *hmass1(0);
  double FITRMS(2.0);
  for (int i = 0; i < 13; ++i) {
    hmass0 = (TH1D*)fHistFile->Get(Form("Hmass0_%d", i));
    hmass1 = (TH1D*)fHistFile->Get(Form("Hmass1_%d", i));
    cout << "mass = " << hmass0->GetMean() << " RMS = " << hmass0->GetRMS() << endl;
    h1->SetBinContent(i+1, hmass0->GetMean());
    h1->SetBinError(i+1, hmass0->GetMeanError());
    s1->SetBinContent(i+1, hmass0->GetRMS());
    s1->SetBinError(i+1, hmass0->GetRMSError());

    cout << "mass = " << hmass1->GetMean() << " RMS = " << hmass1->GetRMS() << endl;
    h2->SetBinContent(i+1, hmass1->GetMean());
    h2->SetBinError(i+1, hmass1->GetMeanError());
    s2->SetBinContent(i+1, hmass1->GetRMS());
    s2->SetBinError(i+1, hmass1->GetRMSError());

    setFilledHist(hmass0, kBlue, kBlue, 3365);
    setFilledHist(hmass1, kRed, kRed, 3354);

    // -- do the fitting before the scaling
    double peak0V(0.),  peak0E(0.),  peak1V(0.),  peak1E(0.);
    double sigma0V(0.), sigma0E(0.), sigma1V(0.), sigma1E(0.);
    if (hmass1->GetSumOfWeights() > 0.) {
      hmass1->Fit("gaus", "r0", "", 5.37 - FITRMS*hmass1->GetRMS(), 5.37 + FITRMS*hmass1->GetRMS());
      peak1V  = hmass1->GetFunction("gaus")->GetParameter(1);
      peak1E  = hmass1->GetFunction("gaus")->GetParError(1);
      sigma1V = hmass1->GetFunction("gaus")->GetParameter(2);
      sigma1E = hmass1->GetFunction("gaus")->GetParError(2);
      hmass1->GetFunction("gaus")->SetLineColor(kRed);
      p2->SetBinContent(i+1, peak1V);
      p2->SetBinError(i+1, peak1E);
      w2->SetBinContent(i+1, sigma1V);
      w2->SetBinError(i+1, sigma1E);
    }
    if (hmass0->GetSumOfWeights() > 0.) {
      hmass0->Fit("gaus", "r0", "", 5.37 - FITRMS*hmass0->GetRMS(), 5.37 + FITRMS*hmass0->GetRMS());
      peak0V  = hmass0->GetFunction("gaus")->GetParameter(1);
      peak0E  = hmass0->GetFunction("gaus")->GetParError(1);
      sigma0V = hmass0->GetFunction("gaus")->GetParameter(2);
      sigma0E = hmass0->GetFunction("gaus")->GetParError(2);
      hmass0->GetFunction("gaus")->SetLineColor(kBlue);
      p1->SetBinContent(i+1, peak0V);
      p1->SetBinError(i+1, peak0E);
      w1->SetBinContent(i+1, sigma0V);
      w1->SetBinError(i+1, sigma0E);
    }

    if (hmass0->GetSumOfWeights() > 0.) hmass0->Scale(1./hmass0->GetSumOfWeights());
    if (hmass1->GetSumOfWeights() > 0.) hmass1->Scale(1./hmass1->GetSumOfWeights());

    hmass1->SetMaximum(1.2*hmass1->GetMaximum()/hmass1->GetSumOfWeights());
    hmass1->Draw();
    hmass0->Draw("same");

    tl->SetTextSize(0.04); tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.2, 0.92, Form("%s", hmass0->GetTitle()));
    tl->SetTextSize(0.03); tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.65, 0.85, "Run-2");
    tl->SetTextSize(0.03); tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.65, 0.80, Form("RMS: %4.3f MeV", hmass0->GetRMS()));
    tl->SetTextSize(0.03); tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.65, 0.76, Form("peak: %5.4f MeV", peak0V));

    tl->SetTextSize(0.03); tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.65, 0.65, "Phase-2");
    tl->SetTextSize(0.03); tl->SetTextColor(kRed);   tl->DrawLatexNDC(0.65, 0.60, Form("RMS: %4.3f MeV", hmass1->GetRMS()));
    tl->SetTextSize(0.03); tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.65, 0.56, Form("peak: %5.4f MeV", peak1V));
    savePad(Form("ups1-mass-bin%d.pdf", i));
  }

  c0->Clear();
  setHist(p1, kBlue, 24, 1.2);
  setHist(p2, kRed, 25, 1.2);
  p1->SetMinimum(5.3);
  p1->SetMaximum(5.4);
  setTitles(p1, "|#it{#eta}_{f}|", "MPV(#it{m_{#mu#mu}}) [GeV]", 0.05, 1.1, 1.6);
  p1->Draw("e");
  p2->Draw("esame");

  tl->SetTextColor(kBlack);
  //  stamp(0.25, "CMS Phase-2 ", "Simulation ", 0., "");
  newLegend(0.25, 0.28, 0.45, 0.42);
  legg->SetTextSize(0.05);
  legg->AddEntry(p1,  "Run 2", "p");
  legg->AddEntry(p2, "Phase-2", "p");
  legg->Draw();

  savePad(Form("ups1-massMPV-vsEta.pdf"));


  c0->Clear();
  setHist(w1, kBlue, 24, 1.7);
  setHist(w2, kRed, 25, 1.7);
  w1->SetMinimum(0.);
  w1->SetMaximum(0.22);
  setTitles(w1, "|#it{#eta}_{f}|", "#sigma(#it{m_{#mu#mu}}) [GeV]", 0.05, 1.1, 1.6);
  w1->Draw("e");
  w2->Draw("esame");

  stamp(0.25, "CMS Phase-2", "Simulation ", 0., "");
  newLegend(0.25, 0.48, 0.45, 0.62);
  legg->SetTextSize(0.05);
  legg->AddEntry(w1,  "Run 2", "p");
  legg->AddEntry(w2, "Phase-2", "p");
  legg->Draw();

  savePad("ups1-massSigma-vsEta.pdf");


  c0->Clear();
  setHist(h1, kBlue, 24, 1.2);
  setHist(h2, kRed, 25, 1.2);
  h1->SetMinimum(5.3);
  h1->SetMaximum(5.4);
  setTitles(h1, "|#it{#eta}_{f}|", "mean(#it{m_{#mu#mu}}) [GeV]", 0.05, 1.1, 1.6);
  h1->Draw("e");
  h2->Draw("esame");

  //  stamp(0.25, "CMS Phase-2", "Simulation ", 0., "");
  newLegend(0.25, 0.28, 0.45, 0.42);
  legg->SetTextSize(0.05);
  legg->AddEntry(h1,  "Run 2", "p");
  legg->AddEntry(h2, "Phase-2", "p");
  legg->Draw();

  savePad(Form("ups1-massMean-vsEta.pdf"));


  c0->Clear();
  setHist(s1, kBlue, 24, 1.7);
  setHist(s2, kRed, 25, 1.7);
  setTitles(s1, "|#it{#eta}_{f}|", "RMS(#it{m_{#mu#mu}}) [GeV]", 0.05, 1.1, 1.7);
  s1->SetMinimum(0.);
  s1->SetMaximum(0.22);
  s1->Draw("e");
  s2->Draw("esame");

  stamp(0.25, "CMS Phase-2", "Simulation ", 0., "");
  newLegend(0.25, 0.48, 0.45, 0.62);
  legg->SetTextSize(0.05);
  legg->AddEntry(s1,  "Run 2", "p");
  legg->AddEntry(s2, "Phase-2", "p");
  legg->Draw();

  savePad("ups1-massRms-vsEta.pdf");

}



// ----------------------------------------------------------------------
void plotWork::ups2(std::string file1, std::string file2) {

  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;


  TH1D *h = (TH1D*)fHistFile->Get("HBd0_1");
  if (0 == h) {
    cout << "hist not found, unning over files" << endl;

    for (int i = 0; i < 13; ++i) {
      fHBd0.push_back(new TH1D(Form("HBd0_%d", i), Form(" %3.1f < |#eta_{f}| < %3.1f", i*0.2, (i+1)*0.2), 80, 5.0, 5.8));
      fHBd1.push_back(new TH1D(Form("HBd1_%d", i), Form(" %3.1f < |#eta_{f}| < %3.1f", i*0.2, (i+1)*0.2), 80, 5.0, 5.8));

      fHBs0.push_back(new TH1D(Form("HBs0_%d", i), Form(" %3.1f < |#eta_{f}| < %3.1f", i*0.2, (i+1)*0.2), 80, 5.0, 5.8));
      fHBs1.push_back(new TH1D(Form("HBs1_%d", i), Form(" %3.1f < |#eta_{f}| < %3.1f", i*0.2, (i+1)*0.2), 80, 5.0, 5.8));

      fHBg0.push_back(new TH1D(Form("HBg0_%d", i), Form(" %3.1f < |#eta_{f}| < %3.1f", i*0.2, (i+1)*0.2), 80, 5.0, 5.8));
      fHBg1.push_back(new TH1D(Form("HBg1_%d", i), Form(" %3.1f < |#eta_{f}| < %3.1f", i*0.2, (i+1)*0.2), 80, 5.0, 5.8));

    }

    int NEVT(-1);

    //    TFile *f0 = TFile::Open("/scratch/ursl/bmm4/v06/bmm-mc-RunIISpring16DR80-BsToMuMu_BMuonFilter-v03.root");
    TFile *f0 = TFile::Open("/scratch/ursl/bmm4/v06/bmm-mc-RunIISpring16DR80-BsToMuMu_BMuonFilter-v06.root");
    TTree *t0 = (TTree*)f0->Get("candAnaMuMu/events");
    setupTree(t0);
    fSetupInt = 2;
    loopOverTree(t0, 2, NEVT);
    f0->Close();

    TFile *f1 = TFile::Open("/scratch/ursl/bmm4/ups/bmm-mc-ups-bsmmExtended.root");
    TTree *t1 = (TTree*)f1->Get("candAnaMuMu/events");
    fSetupInt = 3;
    setupTree(t1);
    loopOverTree(t1, 2, NEVT);
    f1->Close();

    TFile *f2 = TFile::Open("/scratch/ursl/bmm4/v03/bmm-mc-RunIISpring16DR80-BdToMuMu_BMuonFilter-v03.root");
    TTree *t2 = (TTree*)f2->Get("candAnaMuMu/events");
    setupTree(t2);
    fSetupInt = 4;
    loopOverTree(t2, 2, NEVT);
    f2->Close();

    TFile *f3 = TFile::Open("/scratch/ursl/bmm4/ups/bmm-mc-ups-bdmmExtended.root");
    TTree *t3 = (TTree*)f3->Get("candAnaMuMu/events");
    fSetupInt = 5;
    setupTree(t3);
    loopOverTree(t3, 2, NEVT);
    f3->Close();

    TFile *f4 = TFile::Open("/scratch/ursl/bmm4/v03/bmm-mc-RunIISpring16DR80-BdToPiMuNu_BMuonFilter-v03.root");
    TTree *t4 = (TTree*)f4->Get("candAnaMuMu/events");
    setupTree(t4);
    fSetupInt = 6;
    loopOverTree(t4, 2, NEVT);
    f4->Close();

    TFile *f5 = TFile::Open("/scratch/ursl/bmm4/ups/bmm-mc-ups-bdpimunuExtended.root");
    TTree *t5 = (TTree*)f5->Get("candAnaMuMu/events");
    fSetupInt = 7;
    setupTree(t5);
    loopOverTree(t5, 2, NEVT);
    f5->Close();

    fHistFile->cd();
    for (int i = 0; i < 13; ++i) {
      fHBd0[i]->Write();
      fHBd1[i]->Write();

      fHBs0[i]->Write();
      fHBs1[i]->Write();

      fHBg0[i]->Write();
      fHBg1[i]->Write();
    }
    fHistFile->Close();
  }

  fHistFile = TFile::Open(fHistFileName.c_str());
  gStyle->SetOptFit(0);

  TH1D *hbs0Comb = new TH1D("hbs0Comb", "hbs0Comb", 80, 5.0, 5.8);
  TH1D *hbs1Comb = new TH1D("hbs1Comb", "hbs1Comb", 80, 5.0, 5.8);
  TH1D *hbd0Comb = new TH1D("hbd0Comb", "hbd0Comb", 80, 5.0, 5.8);
  TH1D *hbd1Comb = new TH1D("hbd1Comb", "hbd1Comb", 80, 5.0, 5.8);

  shrinkPad(0.15, 0.19);
  double eps(0.01);
  // -- N(B0) / N(Bs) in B0 signal window
  TH1D *h0 = new TH1D("dmass0", "Run 2", 13, 0., 2.6); h0->Sumw2();
  TH1D *h1 = new TH1D("dmass1", "Phase 2", 13, 0.+eps, 2.6+eps); h1->Sumw2();

  // -- Separation of Bs and B0 in terms of peak / TMath::Sqrt(sig1*sig1 + sig2*sig2)
  TH1D *hSep0 = new TH1D("hsep0", "Run 2", 13, 0., 2.6); hSep0->Sumw2();
  TH1D *hSep1 = new TH1D("hsep1", "Phase 2", 13, 0.+eps, 2.6+eps); hSep1->Sumw2();

  // -- EXPO slope of background
  TH1D *hBg0 = new TH1D("bgslope0", "Run 2", 13, 0., 2.6); hBg0->Sumw2();
  TH1D *hBg1 = new TH1D("bgslope1", "Phase 2", 13, 0.+eps, 2.6+eps); hBg1->Sumw2();
  double nbs, nbd;
  double nbs0E, nbs1E, nbd0E, nbd1E;
  double FITRMS(2.0);
  gStyle->SetOptFit(1);
  for (int i = 0; i < 13; ++i) {
    TH1D *hbs0 = (TH1D*)fHistFile->Get(Form("HBs0_%d", i));
    TH1D *hbs1 = (TH1D*)fHistFile->Get(Form("HBs1_%d", i));
    TH1D *hbd0 = (TH1D*)fHistFile->Get(Form("HBd0_%d", i));
    TH1D *hbd1 = (TH1D*)fHistFile->Get(Form("HBd1_%d", i));
    TH1D *hbg0 = (TH1D*)fHistFile->Get(Form("HBg0_%d", i));
    TH1D *hbg1 = (TH1D*)fHistFile->Get(Form("HBg1_%d", i));
    int lobin(hbs0->FindBin(5.2)), hibin(hbs0->FindBin(5.3));

    if (hbs0->Integral() > 0 && hbd0->Integral() > 0) {
      hbs0->Fit("gaus", "r", "", 5.37 - FITRMS*hbs0->GetRMS(), 5.37 + FITRMS*hbs0->GetRMS());
      savePad(Form("ups2-temp-hbs0-%d.pdf", i));
      hbd0->Fit("gaus", "r", "", 5.28 - FITRMS*hbd0->GetRMS(), 5.28 + FITRMS*hbd0->GetRMS());
      savePad(Form("ups2-temp-hbd0-%d.pdf", i));

      double peakA = hbs0->GetFunction("gaus")->GetParameter(1);
      double sigA  = hbs0->GetFunction("gaus")->GetParameter(2);
      double sigAe = hbs0->GetFunction("gaus")->GetParError(2);
      double peakB = hbd0->GetFunction("gaus")->GetParameter(1);
      double sigB  = hbd0->GetFunction("gaus")->GetParameter(2);
      double sigBe = hbd0->GetFunction("gaus")->GetParError(2);
      double sigC  = TMath::Sqrt(sigA*sigA + sigB*sigB);
      double sigCe = TMath::Sqrt(4*sigA*sigA*sigAe*sigAe  + 4*sigB*sigB*sigBe*sigBe)/(2.*sigC);
      hSep0->SetBinContent(i+1, (peakA-peakB)/sigC);
      double err = (peakA-peakB)/sigC - (peakA-peakB)/(sigC+sigCe);
      hSep0->SetBinError(i+1, err);
      cout << hbs0->GetTitle() << ": peakA-peakB = " << peakA << " - " << peakB << " = " << peakA-peakB
	   << " sigA = " << sigA << " sigB = " << sigB << ", sigC = " << sigC << "+/-" << sigCe << " -> "
	   << (peakA-peakB)/sigC << " +/- " << err
	   << endl;
    }
    if (hbs1->Integral() > 0 && hbd1->Integral() > 0) {
      hbs1->Fit("gaus", "r", "", 5.37 - FITRMS*hbs1->GetRMS(), 5.37 + FITRMS*hbs1->GetRMS());
      savePad(Form("ups2-temp-hbs1-%d.pdf", i));
      hbd1->Fit("gaus", "r", "", 5.28 - FITRMS*hbd1->GetRMS(), 5.28 + FITRMS*hbd1->GetRMS());
      savePad(Form("ups2-temp-hbd1-%d.pdf", i));

      double peakA = hbs1->GetFunction("gaus")->GetParameter(1);
      double sigA  = hbs1->GetFunction("gaus")->GetParameter(2);
      double sigAe = hbs1->GetFunction("gaus")->GetParError(2);
      double peakB = hbd1->GetFunction("gaus")->GetParameter(1);
      double sigB  = hbd1->GetFunction("gaus")->GetParameter(2);
      double sigBe = hbd1->GetFunction("gaus")->GetParError(2);
      double sigC  = TMath::Sqrt(sigA*sigA + sigB*sigB);
      double sigCe = TMath::Sqrt(4*sigA*sigA*sigAe*sigAe  + 4*sigB*sigB*sigBe*sigBe)/(2.*sigC);
      hSep1->SetBinContent(i+1, (peakA-peakB)/sigC);
      double err = (peakA-peakB)/sigC - (peakA-peakB)/(sigC+sigCe);
      hSep1->SetBinError(i+1, err);
      cout << hbs1->GetTitle() << ": peakA-peakB = " << peakA << " - " << peakB << " = " << peakA-peakB
	   << " sigA = " << sigA << " sigB = " << sigB << ", sigC = " << sigC << "+/-" << sigCe << " -> "
	   << (peakA-peakB)/sigC << " +/- " << err
	   << endl;
    }

    if (i < 8) {
      hbs0Comb->Add(hbs0);
      hbs1Comb->Add(hbs1);

      hbd0Comb->Add(hbd0);
      hbd1Comb->Add(hbd1);
    }

    // -- scale complete histogram to run-2 expectations
    //    this normalization is necessary because else the event counts for Bs and B0 are 'given' by the MC production numbers
    if (hbs0->GetEntries() > 0) {
      nbs0E = TMath::Sqrt(hbs0->Integral(lobin, hibin))/hbs0->Integral(lobin, hibin);
      hbs0->Scale((6.878 + 12.304)/hbs0->Integral());
    }
    if (hbs1->GetEntries() > 0) {
      nbs1E = TMath::Sqrt(hbs1->Integral(lobin, hibin))/hbs1->Integral(lobin, hibin);
      hbs1->Scale((6.878 + 12.304)/hbs1->Integral());
    }
    if (hbd0->GetEntries() > 0) {
      nbd0E = TMath::Sqrt(hbd0->Integral(lobin, hibin))/hbd0->Integral(lobin, hibin);
      hbd0->Scale((0.838 + 1.460)/hbd0->Integral());
    }
    if (hbd1->GetEntries() > 0) {
      nbd1E = TMath::Sqrt(hbd1->Integral(lobin, hibin))/hbd1->Integral(lobin, hibin);
      hbd1->Scale((0.838 + 1.460)/hbd1->Integral());
    }

    if (hbs0->Integral() > 0 && hbd0->Integral() > 0) {
      nbs = hbs0->Integral(lobin, hibin);
      nbd = hbd0->Integral(lobin, hibin);
      cout << i << " nbs0 = " << nbs << "/" << nbd << " = " << nbd/nbs << ", err2 = " <<  TMath::Abs(nbd0E*nbd0E + nbs0E*nbs0E) << endl;
      h0->SetBinContent(i+1, nbd/nbs);
      h0->SetBinError(i+1, (nbd/nbs)*TMath::Abs(nbd0E*nbd0E + nbs0E*nbs0E));
    }
    if (hbs1->Integral() > 0 && hbd1->Integral() > 0) {
      nbs = hbs1->Integral(lobin, hibin);
      nbd = hbd1->Integral(lobin, hibin);
      cout << i << " nbs1 = " << nbs << "/" << nbd << " = " << nbd/nbs << ", err2 = " <<  TMath::Abs(nbd0E*nbd0E + nbs0E*nbs0E) << endl;
      h1->SetBinContent(i+1, nbd/nbs);
      h1->SetBinError(i+1, (nbd/nbs)*TMath::Abs(nbd1E*nbd1E + nbs1E*nbs1E));
    }

  }

  c0->Clear();
  setHist(h0, kBlue, 24, 1.7);
  setHist(h1, kRed, 25, 1.7);
  h0->SetMinimum(0.0);
  h0->SetMaximum(10.);
  setTitles(h0, "|#it{#eta}_{f}|", "#it{N(B^{0})}/#it{N(B^{0}_{s})} in [5.2,5.3] GeV", 0.05, 1.1, 1.6);
  h0->Draw("e");
  h1->Draw("esame");

  stamp(0.25, "CMS Phase-2", "Simulation ", 0., "");
  newLegend(0.6, 0.58, 0.8, 0.72);
  legg->SetTextSize(0.05);
  legg->AddEntry(h0,  "Run 2", "p");
  legg->AddEntry(h1, "Phase-2", "p");
  legg->Draw();

  savePad("ups2-b0-bs-xfeed-vsEta.pdf");

  c0->Clear();
  setHist(hSep0, kBlue, 24, 1.7);
  setHist(hSep1, kRed, 25, 1.7);
  hSep0->SetMinimum(0.0);
  hSep0->SetMaximum(4.);
  setTitles(hSep0, "|#it{#eta}_{f}|", "#it{B^{0}} and #it{B^{0}_{s}} separation significance", 0.05, 1.1, 1.6);
  hSep0->Draw("e");
  hSep1->Draw("esame");

  stamp(0.25, "CMS Phase-2", "Simulation ", 0., "");
  newLegend(0.6, 0.58, 0.8, 0.72);
  legg->SetTextSize(0.05);
  legg->AddEntry(hSep0,  "Run 2", "p");
  legg->AddEntry(hSep1, "Phase-2", "p");
  legg->Draw();

  savePad("ups2-b0-bs-sep-vsEta.pdf");

  c0->Clear();
  //Bs: 6.878 + 12.304
  //Bz: 0.838 + 1.460
  double bdratio = (0.838 + 1.460)/(6.878 + 12.304);
  hbd0Comb->Scale(bdratio/hbd0Comb->Integral(1, hbd0Comb->GetNbinsX()+1));
  hbd1Comb->Scale(bdratio/hbd1Comb->Integral(1, hbd1Comb->GetNbinsX()+1));

  hbs0Comb->Scale(1./hbs0Comb->Integral(1, hbs0Comb->GetNbinsX()+1));
  hbs1Comb->Scale(1./hbs1Comb->Integral(1, hbs1Comb->GetNbinsX()+1));

  setFilledHist(hbs0Comb, kBlue+2, kBlue+2, 3365);
  setFilledHist(hbs1Comb, kBlue+2, kBlue+2, 3365);

  setFilledHist(hbd0Comb, kRed+2, kRed+2, 3354);
  setFilledHist(hbd1Comb, kRed+2, kRed+2, 3354);

  setTitles(hbs0Comb, "#it{m_{#mu#mu}} [GeV]", "a.u.", 0.05, 1.1, 1.6);
  hbs0Comb->SetMaximum(1.3*hbs0Comb->GetMaximum());
  hbs0Comb->Draw();
  hbd0Comb->Draw("same");

  stamp(0.25, "CMS Phase-2", "Simulation ", 0., "");
  newLegend(0.65, 0.48, 0.85, 0.72);
  legg->SetTextSize(0.05);
  //  legg->SetHeader("|#it{#eta}_{f}| < 1.4");
  legg->SetHeader("Run 2");
  legg->AddEntry(hbs0Comb, "#it{B^{0}_{s}}", "f");
  legg->AddEntry(hbd0Comb, "#it{B^{0}}", "f");
  legg->Draw();

  savePad(Form("ups2-b0-bs-run2-overlay.pdf"));

  setTitles(hbs1Comb, "#it{m_{#mu#mu}} [GeV]", "a.u.", 0.05, 1.1, 1.6);
  hbs1Comb->SetMaximum(1.3*hbs1Comb->GetMaximum());
  hbs1Comb->Draw();
  hbd1Comb->Draw("same");
  stamp(0.25, "CMS Phase-2", "Simulation ", 0., "");
  newLegend(0.65, 0.48, 0.85, 0.72);
  legg->SetTextSize(0.05);
  //  legg->SetHeader("|#it{#eta}_{f}| < 1.4");
  legg->SetHeader("Phase-2");
  legg->AddEntry(hbs1Comb, "#it{B^{0}_{s}}", "f");
  legg->AddEntry(hbd1Comb, "#it{B^{0}}", "f");
  legg->Draw();
  savePad("ups2-b0-bs-phase2-overlay.pdf");



}

// ----------------------------------------------------------------------
void plotWork::loopFunction2() {
  if ((TMath::Abs(fb.m2eta) < 1.4) && fb.m2pt < 4.) return;
  if ((TMath::Abs(fb.m2eta) > 1.4) && fb.m2pt < 2.) return;
  if (!fb.gmugmid) return;

  //  if (TMath::Abs(fb.fls3d) < 10.) return;

  // if (TMath::Abs(fb.maxdoca) > 0.08) return;
  // if (TMath::Abs(fb.pvips) > 5.) return;



  double meta = fb.m1eta;
  if (TMath::Abs(meta) < TMath::Abs(fb.m2eta)) meta = fb.m2eta;

  int ieta = TMath::Abs(meta)/0.2;
  //  cout << "eta = " << TMath::Abs(meta) << " -> " << ieta << endl;
  if (ieta > 12) return;

  if (0 == fSetupInt) {
    fHmass0[ieta]->Fill(fb.m);
  }
  if (1 == fSetupInt) {
    fHmass1[ieta]->Fill(fb.m);
  }

  if (2 == fSetupInt) {
    fHBs0[ieta]->Fill(fb.m);
  }
  if (3 == fSetupInt) {
    fHBs1[ieta]->Fill(fb.m);
  }

  if (4 == fSetupInt) {
    fHBd0[ieta]->Fill(fb.m);
  }
  if (5 == fSetupInt) {
    fHBd1[ieta]->Fill(fb.m);
  }

  if (6 == fSetupInt) {
    fHBg0[ieta]->Fill(fb.m);
  }
  if (7 == fSetupInt) {
    fHBg1[ieta]->Fill(fb.m);
  }


}



// ----------------------------------------------------------------------
void plotWork::loopFunction3() {

  if (fChan < 0) return;

  if (TMath::Abs(fb.mkk) < MKKLO) return;
  if (TMath::Abs(fb.mkk) > MKKHI) return;
  if (TMath::Abs(fb.phidr) > DR) return;
  if (fb.k1pt > fb.k2pt) {
    if (TMath::Abs(fb.k1pt) < PTK1) return;
    if (TMath::Abs(fb.k2pt) < PTK2) return;
  } else {
    if (TMath::Abs(fb.k1pt) < PTK2) return;
    if (TMath::Abs(fb.k2pt) < PTK1) return;
  }
  if (TMath::Abs(fb.psipt) < PTPSI) return;


  if (TMath::Abs(fb.flsxy) < fCuts[fChan]->flsxy) return;
  if (TMath::Abs(fb.fls3d) < fCuts[fChan]->fls3d) return;

  if (TMath::Abs(fb.chi2dof) > fCuts[fChan]->chi2dof) return;
  if (TMath::Abs(fb.alpha) > fCuts[fChan]->alpha) return;
  if (TMath::Abs(fb.pvip) > fCuts[fChan]->pvip) return;
  if (TMath::Abs(fb.pvips) > fCuts[fChan]->pvips) return;

  if (TMath::Abs(fb.iso) < fCuts[fChan]->iso) return;
  if (TMath::Abs(fb.docatrk) < fCuts[fChan]->docatrk) return;
  if (TMath::Abs(fb.closetrk) > fCuts[fChan]->closetrk) return;


  //  if (fChan == 3) cout << "eta: " << TMath::Abs(fb.m1eta) << "/" << TMath::Abs(fb.m2eta) << " -> chan = " << fChan << endl;

  fHma[fChan]->Fill(fb.m);
  fHmc[fChan]->Fill(fb.cm);

}


// ----------------------------------------------------------------------
void plotWork::loopFunction4() {

  if (!fb.hlt1) return;
  static int runPrinted(-1);
  if (fb.run != runPrinted) {
    cout << "fb.run = " << fb.run << endl;
    runPrinted = fb.run;
  }

  fpHL1All->Fill(fb.run);
  if (fb.l1s & (0x1<<0)) fpHL1s0->Fill(fb.run);
  if (fb.l1s & (0x1<<1)) fpHL1s1->Fill(fb.run);
  if (fb.l1s & (0x1<<2)) fpHL1s2->Fill(fb.run);
  if (fb.l1s & (0x1<<3)) fpHL1s3->Fill(fb.run);
  if (fb.l1s & (0x1<<4)) fpHL1s4->Fill(fb.run);
  if (fb.l1s & (0x1<<5)) fpHL1s5->Fill(fb.run);
  if (fb.l1s & (0x1<<6)) fpHL1s6->Fill(fb.run);

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
  cout << "==> plotWork::loopOverTree> loop over dataset " << (fCds?fCds->fName:"undefined") << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotWork::*pF)(void);
  if (ifunc == 1) pF = &plotWork::loopFunction1;
  if (ifunc == 2) pF = &plotWork::loopFunction2;
  if (ifunc == 3) pF = &plotWork::loopFunction3;
  if (ifunc == 4) pF = &plotWork::loopFunction4;

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
void plotWork::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotWork::loadFile loading files listed in " << files << endl;

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
