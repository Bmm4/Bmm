#include "plotStuff.hh"

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

ClassImp(plotStuff)

using namespace std;

// ----------------------------------------------------------------------
plotStuff::plotStuff(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotStuff::loadFiles(files);

  changeSetup(dir, "plotStuff", setup);
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
plotStuff::~plotStuff() {

}



// ----------------------------------------------------------------------
void plotStuff::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotStuff::makeAll(string what) {

  if (what == "dbx") {
    //    changeSetup("results", "yieldstability", "");
    yieldStability("bupsikData", "HLT");
    yieldStability("bmmData", "HLT");
    yieldStability("bspsiphiData", "HLT");
    yieldStability("bdpsikstarData", "HLT");
  }

  if (what == "all" || what == "yieldstability") {
    yieldStability("bupsikData", "HLT");
    yieldStability("bmmData", "HLT");
    yieldStability("bspsiphiData", "HLT");
    yieldStability("bdpsikstarData", "HLT");
    yieldStabilityRatios("HLT");
  }

  if (what == "all" || what == "pvstudy") {
    pvStudy("bdmmMcOff", "&& (fl1>0.01)", "fl1");
    pvStudy("bupsikMcOff", "&&(fl1>0.01)", "fl1");

    pvStudy("bdmmMcOff", "&& (fl1>0.01) && (idx1!=idx3)", "fl1_idx");
    pvStudy("bupsikMcOff", "&& (fl1>0.01) && (idx1!=idx3)", "fl1_idx");
  }

  if (what == "all" || what == "pustudy") {
    puStudy("bsmmMcOff");
    puStudy("bupsikMcOff");
  }
}



// ----------------------------------------------------------------------
void plotStuff::bookHist(string dsname) {


}

// ----------------------------------------------------------------------
void plotStuff::puStudy(string dsname) {

  char mode[200]; sprintf(mode, "%s", dsname.c_str());
  // -- book histograms and profiles
  for (int i = 0; i < NCHAN; ++i) {
    fpHmultFar[i]   = new TH1D(Form("pvmultfar_%s_chan%d", mode, i), "PV multiplicity dzmin > 2 cm", 40, 0., 120.);
    setFilledHist(fpHmultFar[i], kBlue, kBlue, 3354);
    setTitles(fpHmultFar[i], "PV track multiplicity", "a.u.", 0.05, 1.1, 1.8);

    fpHmultClose05[i] = new TH1D(Form("pvmultclose05_%s_chan%d", mode, i), "PV multiplicity dzmin < 0.05 cm", 40, 0., 120.);
    setFilledHist(fpHmultClose05[i], kGreen+2, kGreen+2, 3365);

    fpHmultFar2[i]  = new TH1D(Form("pvmultfar2_%s_chan%d", mode, i), "PV multiplicity | l_{z}^{2} | > 1 cm", 40, 0., 120.);
    setFilledHist(fpHmultFar2[i], kBlue, kBlue, 3354);
    setTitles(fpHmultFar2[i], "PV track multiplicity", "a.u.", 0.05, 1.1, 1.8);
    fpHmultClose2[i] = new TH1D(Form("pvmultclose2_%s_chan%d", mode, i), "PV multiplicity  | l_{z}^{2} | < 0.05 cm", 40, 0., 120.);
    setFilledHist(fpHmultClose2[i], kGreen+2, kGreen+2, 3365);


    fpHdzmin[i] = new TH1D(Form("dzmin_%s_chan%d", mode, i), "dzmin", 100, 0., 1.0);
    setFilledHist(fpHdzmin[i], kBlue, kYellow, 1000);
    setTitles(fpHdzmin[i], "minimum |#Delta z|", "a.u.", 0.05, 1.1, 1.8);

    fpHlz1[i] = new TH1D(Form("lz1_%s_chan%d", mode, i), "lz1", 100, 0., 0.2);
    setTitles(fpHlz1[i], "| l_{z}^{(1)} | [cm]", "a.u.", 0.05, 1.1, 1.8);

    fpHlz2[i] = new TH1D(Form("lz2_%s_chan%d", mode, i), "lz2", 100, 0., 0.2);
    setTitles(fpHlz2[i], "| l_{z}^{(2)} | [cm]", "a.u.", 0.05, 1.1, 1.8);

    fpHtmlz1[i] = new TH1D(Form("tmlz1_%s_chan%d", mode, i), "lz1 (correct)", 100, 0., 0.2);
    setTitles(fpHtmlz1[i], "| l_{z}^{(1)} | [cm]", "a.u.", 0.05, 1.1, 1.8);

    // -- vs dzmin
    fpP1Mult[i] = new TProfile(Form("p1mult_%s_chan%d", mode, i), "mult", 40, 0., 2.0, 0., 80., "");
    fpP1Mult[i]->SetMinimum(0.);     fpP1Mult[i]->SetMaximum(80.);
    setTitles(fpP1Mult[i], "minimum |#Delta z| [cm]", "PV trk mult", 0.05, 1.1, 1.8, 0.04);
    fpP1flsxy[i] = new TProfile(Form("p1flsxy_%s_chan%d", mode, i), "flsxy", 40, 0., 2.0, 0., 120., "");
    fpP1flsxy[i]->SetMinimum(0.);     fpP1flsxy[i]->SetMaximum(50.);
    setTitles(fpP1flsxy[i], "minimum |#Delta z| [cm]", "mean flsxy", 0.05, 1.1, 1.8, 0.04);
    fpP1fls3d[i] = new TProfile(Form("p1fls3d_%s_chan%d", mode, i), "fls3d", 40, 0., 2.0, 0., 120., "");
    setTitles(fpP1fls3d[i], "minimum |#Delta z| [cm]", "mean fls3d", 0.05, 1.1, 1.8, 0.04);
    fpP1fls3d[i]->SetMinimum(0.);     fpP1fls3d[i]->SetMaximum(50.);
    fpP1fl3d[i] = new TProfile(Form("p1fl3d_%s_chan%d", mode, i), "fl3d", 40, 0., 2.0, 0., 120., "");
    setTitles(fpP1fl3d[i], "minimum |#Delta z| [cm]", "mean fl3d", 0.05, 1.1, 1.8, 0.04);
    fpP1fl3d[i]->SetMinimum(0.);     fpP1fl3d[i]->SetMaximum(0.5);

    fpP1dfl3d[i] = new TProfile(Form("p1dfl3d_%s_chan%d", mode, i), "dfl3d", 40, 0., 2.0, 0., 50., "");
    fpP1dfl3d[i]->SetMinimum(-0.03); fpP1dfl3d[i]->SetMaximum(0.03);
    setTitles(fpP1dfl3d[i], "minimum |#Delta z| [cm]", "mean #Delta fl3d", 0.05, 1.1, 1.8, 0.04);

    fpP1tau[i] = new TProfile(Form("p1tau_%s_chan%d", mode, i), "tau", 40, 0., 2.0, -1.e-10, 1.e-10, "");
    setTitles(fpP1tau[i], "minimum |#Delta z| [cm]", "mean #tau_{eff} ", 0.05, 1.1, 1.8, 0.04);
    fpP1dtau[i] = new TProfile(Form("p1dtau_%s_chan%d", mode, i), "dtau", 40, 0., 2.0, -1.e-10, 1.e-10, "");
    fpP1dtau[i]->SetMinimum(-0.2e-12); fpP1dtau[i]->SetMaximum(+0.2e-12);
    setTitles(fpP1dtau[i], "minimum |#Delta z| [cm]", "mean #Delta#tau_{eff}", 0.05, 1.1, 1.8, 0.04);

    // -- vs lz2
    fpP2Mult[i] = new TProfile(Form("p2mult_%s_chan%d", mode, i), "mult", 40, 0., 2.0, 0., 80., "");
    fpP2Mult[i]->SetMinimum(0.);     fpP2Mult[i]->SetMaximum(80.);
    setTitles(fpP2Mult[i], "| l_{z}^{(2)} | [cm]", "PV trk mult", 0.05, 1.1, 1.8, 0.04);
    fpP2flsxy[i] = new TProfile(Form("p2flsxy_%s_chan%d", mode, i), "flsxy", 40, 0., 1.0, 0., 120., "");
    fpP2flsxy[i]->SetMinimum(0.);     fpP2flsxy[i]->SetMaximum(50.);
    setTitles(fpP2flsxy[i], "| l_{z}^{(2)} | [cm]", "mean flsxy", 0.05, 1.1, 1.8, 0.04);
    fpP2fls3d[i] = new TProfile(Form("p2fls3d_%s_chan%d", mode, i), "fls3d", 40, 0., 1.0, 0., 120., "");
    setTitles(fpP2fls3d[i], "| l_{z}^{(2)} | [cm]", "mean fls3d", 0.05, 1.1, 1.8, 0.04);
    fpP2fls3d[i]->SetMinimum(0.);     fpP2fls3d[i]->SetMaximum(50.);
    fpP2fl3d[i] = new TProfile(Form("p2fl3d_%s_chan%d", mode, i), "fl3d", 40, 0., 2.0, 0., 120., "");
    setTitles(fpP2fl3d[i], "minimum |#Delta z| [cm]", "mean fl3d", 0.05, 1.1, 1.8, 0.04);
    fpP2fl3d[i]->SetMinimum(0.);     fpP2fl3d[i]->SetMaximum(0.5);

    fpP2dfl3d[i] = new TProfile(Form("p2dfl3d_%s_chan%d", mode, i), "dfl3d", 40, 0., 1.0, 0., 50., "");
    fpP2dfl3d[i]->SetMinimum(-0.03); fpP2dfl3d[i]->SetMaximum(0.03);
    setTitles(fpP2dfl3d[i], "| l_{z}^{(2)} | [cm]", "mean #Delta fl3d", 0.05, 1.1, 1.8, 0.04);

    fpP2tau[i] = new TProfile(Form("p2tau_%s_chan%d", mode, i), "tau", 40, 0., 1.0, -1.e-10, 1.e-10, "");
    setTitles(fpP2tau[i], "| l_{z}^{(2)} | [cm]", "mean #tau_{eff} ", 0.05, 1.1, 1.8, 0.04);
    fpP2dtau[i] = new TProfile(Form("p2dtau_%s_chan%d", mode, i), "dtau", 40, 0., 1.0, -1.e-10, 1.e-10, "");
    fpP2dtau[i]->SetMinimum(-0.2e-12); fpP2dtau[i]->SetMaximum(+0.2e-12);
    setTitles(fpP2dtau[i], "| l_{z}^{(2)} | [cm]", "mean #Delta#tau_{eff}", 0.05, 1.1, 1.8, 0.04);

    // -- with a cut of fls3d, vs dzmin
    fpP3tau[i] = new TProfile(Form("p3tau_%s_chan%d", mode, i), "tau (fls3d > 5)", 40, 0., 1.0, -1.e-10, 1.e-10, "");
    setTitles(fpP3tau[i], "minimum |#Delta z| [cm]", "mean #tau_{eff} ", 0.05, 1.1, 1.8, 0.04);
    fpP3dtau[i] = new TProfile(Form("p3dtau_%s_chan%d", mode, i), "dtau (fls3d > 5)", 40, 0., 1.0, -1.e-10, 1.e-10, "");
    fpP3dtau[i]->SetMinimum(-0.2e-12); fpP3dtau[i]->SetMaximum(+0.2e-12);
    setTitles(fpP3dtau[i], "minimum |#Delta z| [cm]", "mean #Delta#tau_{eff}", 0.05, 1.1, 1.8, 0.04);

  }

  fSample = dsname;
  setup(fSample);

  TTree *t = getTree(dsname, fTreeDir, "pvstudy");
  setupPvTree(t);

  int nevts(-1);
  int nstart(-1);

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
  cout << "==> plotStuff::loopOverPvTree> loop over dataset " << fCds->fName << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  void (plotStuff::*pF)(void);
  pF = &plotStuff::loopFunction2;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
    (this->*pF)();
  }


  // ----------------------------------------------------------------------
  // -- Plot things
  // ----------------------------------------------------------------------

  for (int i = 0; i < NCHAN; ++i) {
    // -- minimum z separation
    shrinkPad(0.13, 0.20);
    fpHdzmin[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-dzmin.pdf", i, mode));

    // -- profiles vs minimum z separation
    fpP1Mult[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof1-mult.pdf", i, mode));
    fpP1flsxy[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof1-flsxy.pdf", i, mode));
    fpP1fls3d[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof1-fls3d.pdf", i, mode));
    fpP1fl3d[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof1-fl3d.pdf", i, mode));

    fpP1dfl3d[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof1-dfl3d.pdf", i, mode));

    fpP1tau[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof1-tau.pdf", i, mode));
    fpP1dtau[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof1-dtau.pdf", i, mode));

    // -- minimum z separation
    shrinkPad(0.13, 0.20);
    gPad->SetLogy(1);
    setFilledHist(fpHlz1[i], kBlue, kYellow, 1000);
    fpHlz1[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-lz1.pdf", i, mode));

    gPad->SetLogy(1);
    fpHlz2[i]->Draw();
    setFilledHist(fpHlz2[i], kBlue, kYellow, 1000);
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-lz2.pdf", i, mode));

    setFilledHist(fpHlz1[i], kBlue, kBlue, 3354);
    setFilledHist(fpHlz2[i], kRed, kRed, 3365);
    setTitles(fpHlz1[i], "| l_{z} | [cm]", "a.u.", 0.05, 1.1, 1.8);
    fpHlz1[i]->Draw();
    fpHlz2[i]->Draw("same");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));

    newLegend(0.41, 0.7, 0.71, 0.90);
    legg->SetHeader("longitudinal IP to PV");
    legg->SetTextSize(0.05);
    legg->AddEntry(fpHlz1[i], "best PV", "f");
    legg->AddEntry(fpHlz2[i], "second-best PV", "f");
    legg->Draw();

    savePad(Form("pustudyChan%d-%s-overlay-lz1-lz2.pdf", i, mode));

    // -- profiles vs lz2
    gPad->SetLogy(0);
    fpP2Mult[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof2-mult.pdf", i, mode));
    fpP2flsxy[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof2-flsxy.pdf", i, mode));
    fpP2fls3d[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof2-fls3d.pdf", i, mode));
    fpP2fl3d[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof2-fl3d.pdf", i, mode));

    fpP2dfl3d[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof2-dfl3d.pdf", i, mode));

    fpP2tau[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof2-tau.pdf", i, mode));
    fpP2dtau[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof2-dtau.pdf", i, mode));

    fpP3tau[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof3-tau.pdf", i, mode));
    fpP3dtau[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof3-dtau.pdf", i, mode));

    // -- PV multiplicity
    shrinkPad(0.13, 0.20);
    overlay(fpHmultFar[i], "close", fpHmultClose05[i], "far", 0, "", UNITY, false, false);
    newLegend(0.41, 0.7, 0.71, 0.90);
    legg->SetHeader("other PV");
    legg->SetTextSize(0.05);
    legg->AddEntry(fpHmultFar[i], "#it{|d_{z}^{min}| > 1cm}", "f");
    legg->AddEntry(fpHmultClose05[i], "#it{|d_{z}^{min}| < 0.05cm}", "f");
    legg->Draw();

    tl->SetTextSize(0.04);
    tl->SetTextColor(kGreen+2);
    tl->DrawLatexNDC(0.40, 0.62, Form("#mu/RMS = %3.1f/%3.1f", fpHmultClose05[i]->GetMean(), fpHmultClose05[i]->GetRMS()));
    tl->SetTextColor(kBlue);
    tl->DrawLatexNDC(0.5, 0.5, Form("#mu/RMS = %3.1f/%3.1f", fpHmultFar[i]->GetMean(), fpHmultFar[i]->GetRMS()));
    tl->SetTextColor(kBlack);

    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-farCloseMultiplicity.pdf", i, mode));


    // -- PV multiplicity lz2
    shrinkPad(0.13, 0.20);
    overlay(fpHmultFar2[i], "close", fpHmultClose2[i], "far", 0, "", UNITY, false, false);
    newLegend(0.41, 0.7, 0.71, 0.90);
    legg->SetHeader("other PV");
    legg->SetTextSize(0.05);
    legg->AddEntry(fpHmultFar2[i], "#it{| l_{z}^{2} | > 1cm}", "f");
    legg->AddEntry(fpHmultClose2[i], "#it{| l_{z}^{2} | < 0.05cm}", "f");
    legg->Draw();

    tl->SetTextSize(0.03);
    tl->SetTextColor(kGreen+2);
    tl->DrawLatexNDC(0.60, 0.62, Form("#mu/RMS = %3.1f/%3.1f", fpHmultClose2[i]->GetMean(), fpHmultClose2[i]->GetRMS()));
    tl->SetTextColor(kBlue);
    tl->DrawLatexNDC(0.65, 0.5, Form("#mu/RMS = %3.1f/%3.1f", fpHmultFar2[i]->GetMean(), fpHmultFar2[i]->GetRMS()));
    tl->SetTextColor(kBlack);

    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.75, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-farCloseMultiplicity2.pdf", i, mode));
  }
}

// ----------------------------------------------------------------------
void plotStuff::pvStudy(string dsname, string selection, string fmod) {
  fSample = dsname;
  string dir = "candAnaBu2JpsiK";
  if (string::npos != fSample.find("mm")) {
    dir = "candAnaMuMu";
  } else if (string::npos != fSample.find("bspsiphi")) {
    dir = "candAnaBs2JpsiPhi";
  } else if (string::npos != fSample.find("bdpsikstar")) {
    dir = "candAnaBd2JpsiKstar";
  }
  TTree *T = getTree(dsname, dir, "pvstudy");

  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");

  TH1D *h1(0);
  for (int i = 0; i < fNchan; ++i) {
    h1  = new TH1D(Form("dd1_ch%d", i), Form("dist(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "dist(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dd2_ch%d", i), Form("dist(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "dist(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dd3_ch%d", i), Form("dist(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "dist(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("dx1_ch%d", i), Form("dx(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta x(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dx2_ch%d", i), Form("dx(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta x(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dx3_ch%d", i), Form("dx(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta x(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("dy1_ch%d", i), Form("dy(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta y(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dy2_ch%d", i), Form("dy(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta y(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dy3_ch%d", i), Form("dy(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta y(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("dz1_ch%d", i), Form("dz(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "#Delta z(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dz2_ch%d", i), Form("dz(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "#Delta z(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dz3_ch%d", i), Form("dz(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "#Delta z(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("dt1_ch%d", i), Form("#delta t(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dt2_ch%d", i), Form("#delta t(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dt3_ch%d", i), Form("#delta t(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("ds1_ch%d", i), Form("#delta t^{2D}(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t^{2D}(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("ds2_ch%d", i), Form("#delta t^{2D}(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t^{2D}(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("ds3_ch%d", i), Form("#delta t^{2D}(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t^{2D}(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
  }

  for (int i = 0; i < 4; ++i) {
    T->Draw(Form("d1 >> dd1_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("d2 >> dd2_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("d3 >> dd3_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));

    T->Draw(Form("p1z-gz >> dz1_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("p2z-gz >> dz2_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("p3z-gz >> dz3_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));

    T->Draw(Form("p1x-gx >> dx1_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("p2x-gx >> dx2_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("p3x-gx >> dx3_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));

    T->Draw(Form("p1y-gy >> dy1_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("p2y-gy >> dy2_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("p3y-gy >> dy3_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));

    T->Draw(Form("t1-gt >> dt1_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("t2-gt >> dt2_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("t3-gt >> dt3_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));

    T->Draw(Form("s1-gs >> ds1_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("s2-gs >> ds2_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
    T->Draw(Form("s3-gs >> ds3_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
  }

  TH1D *h2(0);
  tl->SetTextSize(0.03);
  makeCanvas(1);
  zone(4, 1, c1);
  tl->SetTextSize(0.06);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dd1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dd3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.25, 0.80, Form("r1: %5.3f/%5.3f [cm]", h1->GetMean(), h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.25, 0.73, Form("r3: %5.3f/%5.3f [cm]", h2->GetMean(), h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  savePad(Form("pvstudy-dd-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dx1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dx3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.25, 0.80, Form("r1: %5.3f/%5.3f [cm]", h1->GetMean(), h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.25, 0.73, Form("r3: %5.3f/%5.3f [cm]", h2->GetMean(), h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dx.pdf"), c1);
  savePad(Form("pvstudy-dx-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dy1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dy3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.25, 0.80, Form("r1: %5.3f/%5.3f [cm]", h1->GetMean(), h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.25, 0.73, Form("r3: %5.3f/%5.3f [cm]", h2->GetMean(), h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dy.pdf"), c1);
  savePad(Form("pvstudy-dy-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dz1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dz3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.25, 0.80, Form("r1: %5.3f/%5.3f [cm]", h1->GetMean(), h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.25, 0.73, Form("r3: %5.3f/%5.3f [cm]", h2->GetMean(), h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dz.pdf"), c1);
  savePad(Form("pvstudy-dz-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);


  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dt1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dt3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->SetMaximum(13.*h1->GetMaximum());
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.2, 0.84, Form("r1: %5.3f/%5.3f [ps]", 1.e12*h1->GetMean(), 1.e12*h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.2, 0.79, Form("r3: %5.3f/%5.3f [ps]", 1.e12*h2->GetMean(), 1.e12*h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dt.pdf"), c1);
  savePad(Form("pvstudy-dt-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("ds1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("ds3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->SetMaximum(13.*h1->GetMaximum());
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.2, 0.84, Form("r1: %5.3f/%5.3f [ps]", 1.e12*h1->GetMean(), 1.e12*h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.2, 0.79, Form("r3: %5.3f/%5.3f [ps]", 1.e12*h2->GetMean(), 1.e12*h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dt2d.pdf"), c1);
  savePad(Form("pvstudy-dt2d-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  fHistFile->Close();

}


// ----------------------------------------------------------------------
void plotStuff::yieldStability(string dsname, string trg) {
  double MINLUMI(1000.);
  double mBp(5.28), sBp(0.015), stepBp(5.15);
  double xmin(5.0), xmax(5.9), ymax(0.), expoLo(5.16), expoHi(5.85);

  fSample = dsname;
  fMode = BMM;
  setup(dsname);
  if (string::npos != fSample.find("bspsiphi")) {
    mBp    = 5.369;
    sBp    = 0.015;
    stepBp = 5.15;
  }


  // -- check whether there are any histograms pre-produced already
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  bool ok = fHistFile->cd(fTreeDir.c_str());
  cout << "OK = " << ok << endl;
  TH2D *h2(0), *hBlock(0);
  if (!ok) {
    TDirectory *hDir = fHistFile->mkdir(fTreeDir.c_str());
    fHistFile->cd(fTreeDir.c_str());
    hDir = gDirectory;
    cout << "created " << hDir->GetName() << endl;

    TTree *t = getTree(fSample, fTreeDir);
    if (0 == t) {
      cout << "tree for sample = " << fSample << " not found" << endl;
      return;
    }
    setupTree(t, fSample);
    fCds = fDS[fSample];
    loopOverTree(t, 1);

    cout << "writing output histograms: " << fYieldHLT.size() << endl;
    for (map<string, TH2D*>::iterator it = fYieldHLT.begin(); it != fYieldHLT.end(); ++it) {
      cout << "run " << it->first << endl;
      it->second->Draw("colz");
      it->second->SetDirectory(hDir);
      it->second->Write();
    }
    for (map<string, TH2D*>::iterator it = fYieldRTR.begin(); it != fYieldRTR.end(); ++it) {
      it->second->Draw("colz");
      it->second->SetDirectory(hDir);
      it->second->Write();
    }

    fYieldHLT.clear();
    fYieldRTR.clear();
  } else {
    cout << "histograms exist already, looping over them" << endl;
    TIter next(gDirectory->GetListOfKeys());
    TKey *key(0);
    int run(-1), runMin(9999999), runMax(0), firstLumiRun(99), lastLumiRun(99);
    vector<int> vruns;
    while ((key = (TKey*)next())) {
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
      if (TString(key->GetName()).Contains("_HLT_")) {
	string hname = key->GetName();
	if (0 == hBlock) {
	  hBlock = (TH2D*)((TH2D*)gDirectory->Get(hname.c_str()))->Clone("hBlock");
	  hBlock->Reset();
	}
	replaceAll(hname, "h_HLT_", "");
	run = atoi(hname.c_str());
	if (run > runMax) runMax = run;
	if (run < runMin) runMin = run;
	if (find(vruns.begin(), vruns.end(), run) == vruns.end()) {
	  vruns.push_back(run);
	  //if (run == 283946) vruns.push_back(run);
	}
	//	if (run < 278800) continue;
	//	if (run > 274000) break;
      }
    }
    fHistFile->Close();

    if (vruns.size() > 0) {
      cout << "analyzing runs " << runMin << " .. " <<  runMax << endl;

      // -- create run blocks based on integrated lumi
      Lumi lumi("../common/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
      firstLumiRun = lumi.firstRun();
      lastLumiRun  = lumi.lastRun();
      cout << "lumiRuns = " << firstLumiRun << " .. " << lastLumiRun << endl;
      double intLumi(0.);
      map<pair<int, double>, vector<int> > runBlocks;
      vector<int> segment;
      for (unsigned int irun = 0; irun < vruns.size(); ++irun) {
	intLumi += lumi.lumi(vruns[irun]);
	segment.push_back(vruns[irun]);
	if (intLumi > MINLUMI) {
	  runBlocks.insert(make_pair(make_pair(segment[0], intLumi), segment));
	  intLumi = 0.;
	  segment.clear();
	}
      }

      // -- get the histograms
      fHistFile = TFile::Open(fHistFileName.c_str());
      string hname("");
      for (int ichan = 0; ichan < fNchan; ++ichan) {
	cout << "--> chan " << ichan << endl;
	for (map<pair<int, double>, vector<int> >::iterator it = runBlocks.begin(); it != runBlocks.end(); ++it) {
	  cout << Form("new block: %d %4.1f: ", it->first.first, it->first.second) << endl;
	  for (unsigned int i = 0; i < it->second.size(); ++i) {
	    hname = Form("%s/h_%s_%d_chan%d", fTreeDir.c_str(), trg.c_str(), it->second[i], ichan);
	    h2 = (TH2D*)(fHistFile->Get(hname.c_str()));
	    cout << it->second[i] << " (" << hname << ": " << h2 << ") ";
	    if (0 == h2) continue;
	    h2->Draw("colz");
	    c0->Modified();
	    c0->Update();
	  }
	  cout << endl;
	}
      }
    }
  }


}

// ----------------------------------------------------------------------
void plotStuff::yieldStabilityOld(string dsname, string trg) {
  int MAXPS(20);
  double MINLUMI(1000.);
  double mBp(5.28), sBp(0.015), stepBp(5.15);
  double xmin(5.0), xmax(5.9), ymax(0.), expoLo(5.16), expoHi(5.85);
  fIF->fName = "fit";
  TF1 *f1 = fIF->pol1Err2gauss2c(xmin, xmax);
  fIF->fName = "comp";
  TF1 *fg  = fIF->gauss2c(xmin, xmax);
  fg->SetLineColor(kBlue+1);
  TF1 *fe = fIF->err2(xmin, xmax);
  fe->SetLineColor(kRed+2);
  fe->SetLineStyle(kSolid);
  TF1 *fp = fIF->pol1(xmin, xmax);
  fp->SetLineColor(kRed+2);
  fp->SetLineStyle(kSolid);

  gStyle->SetOptStat(11111);

  fIF->fVerbose = false;

  // -- dump histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;


  fSample = dsname;
  fMode = BMM;
  string dir = "candAnaMuMu";
  if (string::npos != fSample.find("bupsik")) {
    fMode = BU2JPSIKP;
    dir = "candAnaBu2JpsiK";
  }
  if (string::npos != fSample.find("bdpsikstar")) {
    fMode = BD2JPSIKSTAR;
    dir = "candAnaBd2JpsiKstar";
  }
  if (string::npos != fSample.find("bspsiphi")) {
    fMode = BS2JPSIPHI;
    mBp    = 5.369;
    sBp    = 0.015;
    stepBp = 5.15;
    dir    = "candAnaBs2JpsiPhi";
  }


  // -- check whether there are any histograms pre-produced already
  bool ok = fHistFile->cd(dir.c_str());
  cout << "OK = " << ok << endl;
  if (ok) {
    cout << "histograms exist already, looping over them" << endl;
    TIter next(gDirectory->GetListOfKeys());
    TKey *key(0);
    int run(-1), runMin(9999999), runMax(0), firstLumiRun(99), lastLumiRun(99);
    vector<int> vruns;
    while ((key = (TKey*)next())) {
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
      if (TString(key->GetName()).Contains("_HLT_")) {
	string hname = key->GetName();
	replaceAll(hname, "h_HLT_", "");
	run = atoi(hname.c_str());
	if (run > runMax) runMax = run;
	if (run < runMin) runMin = run;
	if (find(vruns.begin(), vruns.end(), run) == vruns.end()) {
	  vruns.push_back(run);
	  //if (run == 283946) vruns.push_back(run);
	}
	//	if (run < 278800) continue;
	//	if (run > 274000) break;
      }
    }

    if (vruns.size() > 0) {
      cout << "analyzing runs " << runMin << " .. " <<  runMax << endl;

      // -- create run blocks based on integrated lumi
      Lumi lumi("../common/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
      firstLumiRun = lumi.firstRun();
      lastLumiRun  = lumi.lastRun();
      cout << "lumiRuns = " << firstLumiRun << " .. " << lastLumiRun << endl;
      double intLumi(0.);
      map<pair<int, double>, vector<int> > runBlocks;
      vector<int> segment;
      for (unsigned int irun = 0; irun < vruns.size(); ++irun) {
	intLumi += lumi.lumi(vruns[irun]);
	segment.push_back(vruns[irun]);
	if (intLumi > MINLUMI) {
	  runBlocks.insert(make_pair(make_pair(segment[0], intLumi), segment));
	  intLumi = 0.;
	  segment.clear();
	}
      }

      // -- print the blocks
      for (map<pair<int, double>, vector<int> >::iterator it = runBlocks.begin(); it != runBlocks.end(); ++it) {
	cout << Form("%d %4.1f: ", it->first.first, it->first.second);
	for (unsigned int i = 0; i < it->second.size(); ++i) {
	  cout << it->second[i] << " ";
	}
	cout << endl;
      }

      //FIXME      return;

      // -- the result histograms
      vector<TH1D *> vRunHLT;
      for (unsigned int ichan = 0; ichan < fNchan; ++ichan) {
	vRunHLT.push_back(new TH1D(Form("hRun%s_%s_chan%d", trg.c_str(), dsname.c_str(), ichan),
				   Form("hRun%s_%s_chan%d", trg.c_str(), dsname.c_str(), ichan),
				   lastLumiRun-firstLumiRun+1, firstLumiRun, lastLumiRun));
	vRunHLT[ichan]->Sumw2();
      }

      // -- book map of histograms for 'all' and individual prescales
      fIF->fLo = xmin;
      fIF->fHi = xmax;
      TH2D *h2(0);
      for (int i = 0; i < MAXPS; ++i) {
	h2 = (TH2D*)(gDirectory->Get(Form("h_%s_%d_%d", trg.c_str(), vruns[0], i)));
	if (h2) {
	  cout << "histogram " << Form("h_%s_%d_%d", trg.c_str(), vruns[0], i) << " found -- breaking" << endl;
	  break;
	}
      }
      if (0 == h2) {
	cout << "did not find any histogram???" << endl;
	return;
      }
      TH2D *h2Sum = (TH2D*)h2->Clone("h2sum"); h2Sum->Reset();
      TH1D *h1 = h2->ProjectionX("chan_0", 1, 1);
      h1->SetName("chan_0"); h1->Reset();
      map<string, TH1D*> vBlockHist;
      for (unsigned int ichan = 0; ichan < fNchan; ++ichan) {
	for (unsigned int ips = 1; ips < MAXPS; ++ips) {
	  vBlockHist.insert(make_pair(Form("chan%d_ps%d", ichan, ips), (TH1D*)h1->Clone(Form("chan%d_ps%d", ichan, ips))));
	  vBlockHist[Form("chan%d_ps%d", ichan, ips)]->Reset();
	  vBlockHist.insert(make_pair(Form("chan%d", ichan), (TH1D*)h1->Clone(Form("chan%d", ichan))));
	  vBlockHist[Form("chan%d", ichan)]->Reset();
	}
      }


      // -- now loop over blocks
      for (map<pair<int, double>, vector<int> >::iterator it = runBlocks.begin(); it != runBlocks.end(); ++it) {
	int iblock = it->first.first;
	double blockLumi = it->first.second;
	cout << Form("block %d %4.1f: ", iblock, blockLumi);
	// -- clear the block histograms
	for (unsigned int ichan = 0; ichan < fNchan; ++ichan) {
	  for (unsigned int ips = 1; ips < MAXPS; ++ips) {
	    vBlockHist[Form("chan%d_ps%d", ichan, ips)]->Reset();
	    vBlockHist[Form("chan%d", ichan)]->Reset();
	  }
	}
	// -- add up all runs in block into per-PS and combined histograms
	for (unsigned int i = 0; i < it->second.size(); ++i) {
	  int irun = it->second[i];
	  // -- get all prescales
	  for (int ips = 1; ips < MAXPS; ++ips) {
	    h2 = (TH2D*)(gDirectory->Get(Form("h_%s_%d_%d", trg.c_str(), irun, ips)));
	    if (h2) {
	      for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
		h1 = h2->ProjectionX("bla", ichan+1, ichan+1);
		vBlockHist[Form("chan%d_ps%d", ichan, ips)]->Add(h1);
		vBlockHist[Form("chan%d", ichan)]->Add(h1);
		delete h1;
	      }
	    }
	  }
	}

	// -- fit the block histograms
	cout << "all: " << vBlockHist[Form("chan%d", 0)]->GetEntries() << endl;
	for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
	  if (fMode != BMM) {
	    // -- fit combined ps hist to determine signal and error function parameters
	    h1 = vBlockHist[Form("chan%d", ichan)];
	    for (int ipar = 0; ipar < f1->GetNpar(); ++ipar) {
	      f1->ReleaseParameter(ipar);
	    }

	    cout << "========> Fitting combined ps for channel " << ichan << " h1->GetSumOfWeights() = " << h1->GetSumOfWeights() << endl;
	    double p0, p1;
	    fIF->fLo = expoLo;
	    fIF->fHi = expoHi;
	    fIF->initPol1(p0, p1, h1);
	    fIF->fLo = xmin;
	    fIF->fHi = xmax;
	    double A   = 0.5*p1*(expoHi*expoHi - expoLo*expoLo) + p0*(expoHi - expoLo);
	    double g0 = (h1->Integral(h1->FindBin(expoLo), h1->FindBin(expoHi))*h1->GetBinWidth(1) - A);
	    double errN = 0.6*(h1->GetBinContent(h1->FindBin(expoLo - 0.1)) - h1->GetBinContent(h1->FindBin(expoLo)));
	    f1->SetParameter(0, 0.8*g0); f1->SetParLimits(0, 0., h1->GetMaximum());
	    f1->SetParameter(1, mBp);    f1->SetParLimits(1, mBp - 1.*sBp, mBp + 1.*sBp);
	    f1->SetParameter(2, sBp);    f1->SetParLimits(2, 0.010, 0.040);
	    f1->SetParameter(3, 0.2);    f1->SetParLimits(3, 0.05, 0.60);
	    f1->SetParameter(4, 5*sBp);  f1->SetParLimits(4, 0.050, 0.150);
	    f1->SetParameter(5, p0);     //f1->SetParLimits(3, 0., 1.e10);
	    f1->SetParameter(6, p1);     //f1->SetParLimits(4, -1.e10, 2.);
	    f1->SetParameter(7, stepBp); f1->SetParLimits(7, stepBp - 2.*sBp, stepBp + 2.*sBp);
	    f1->SetParameter(8, 2.*sBp); f1->SetParLimits(8, 2.*sBp - sBp, 2.*sBp + 1.5*sBp);
	    f1->SetParameter(9, errN);   f1->SetParLimits(9, 0., h1->GetMaximum());

	    h1->Fit(f1, "lr", "", xmin, xmax);
	    fg->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4));
	    fg->Draw("same");
	    for (int ipar = 0; ipar < fe->GetNpar(); ++ipar) fe->SetParameter(ipar, f1->GetParameter(ipar+7));
	    fe->Draw("same");
	    for (int ipar = 0; ipar < fp->GetNpar(); ++ipar) fp->SetParameter(ipar, f1->GetParameter(ipar+5));
	    fp->Draw("same");
	    tl->SetTextSize(0.04);
	    tl->DrawLatexNDC(0.5, 0.5, Form("Signal: %5.1f", fg->Integral(5.15, 5.4)/h1->GetBinWidth(1)));
	    tl->SetTextSize(0.03);
	    tl->DrawLatexNDC(0.2, 0.92, Form("yield-%s-%d-chan%d-allps.pdf", trg.c_str(), iblock, ichan));
	    savePad(Form("yield-%s-%s-%d-chan%d-allps.pdf", trg.c_str(), dsname.c_str(), iblock, ichan));

	    double A0(f1->GetParameter(0));
	    double peak(f1->GetParameter(1));
	    double peakE(f1->GetParError(1));
	    double sigma(f1->GetParameter(2));
	    double sigmaE(f1->GetParError(2));
	    double frac2(f1->GetParameter(3));
	    double frac2E(f1->GetParError(3));
	    double sigma2(f1->GetParameter(4));
	    double sigma2E(f1->GetParError(4));
	    double pol0(f1->GetParameter(5));
	    double pol0E(f1->GetParError(5));
	    double pol1(f1->GetParameter(6));
	    double pol1E(f1->GetParError(6));
	    double step(f1->GetParameter(7));
	    double stepE(f1->GetParError(7));
	    double res(f1->GetParameter(8));
	    double resE(f1->GetParError(8));
	    double level(f1->GetParameter(9));
	    double levelE(f1->GetParError(9));

	    double NSG   = fg->Integral(5.1, 5.5)/h1->GetBinWidth(1);
	    double NTOT  = h1->GetSumOfWeights();
	    double SALL  = NSG/NTOT;
	    double SALLE = TMath::Sqrt(1./NSG + 1./NTOT)*SALL;
	    // -- now fit all histograms for the different prescales
	    double nAll = h1->GetMaximum();
	    double norm(0.), normE(0.);
	    for (int ips = 1; ips < MAXPS; ++ips) {
	      h1 = vBlockHist[Form("chan%d_ps%d", ichan, ips)];
	      //	    h1->SetMinimum(-0.2*h1->GetMaximum());
	      if (0 == h1) {
		cout << "xxx no histogram " << Form("chan%d_ps%d", ichan, ips) << endl;
		continue;
	      }
	      double psNorm(0.), psNormE(0.), par0E(0.), intError(0.);
	      if (h1 && h1->GetSumOfWeights() < 200) {
		// -- FIXME replace with S/B scaled nentries!
		if (h1->GetSumOfWeights() > 0) {
		  h1->Draw();
		  psNorm  = SALL*h1->GetSumOfWeights();
		  psNormE = TMath::Sqrt(1./h1->GetSumOfWeights() + SALLE*SALLE/SALL/SALL)*psNorm;
		  cout << "using S/All scaled histogram entries,  h1->GetSumOfWeights() = " << h1->GetSumOfWeights()
		       << " SALL = " << SALL << " +/- " << SALLE << " -> psNorm = " << psNorm << " +/- " << psNormE
		       << endl;
		} else {
		  continue;
		}
	      } else{
		for (int ipar = 0; ipar < f1->GetNpar(); ++ipar) {
		  f1->ReleaseParameter(ipar);
		}
		double scale = h1->GetMaximum()/nAll;
		cout << "========> Fitting ps  = " << ips << " for channel " << ichan << " h1->GetSumOfWeights() = " << h1->GetSumOfWeights() << endl;
		cout << "SCALE = " << scale << endl;
		f1->SetParameter(0, scale*A0);      f1->SetParLimits(0, 0., 1.e10);
		f1->SetParameter(1, peak);	  f1->SetParLimits(1, peak - peakE,     peak + peakE);
		f1->SetParameter(2, sigma);	  f1->SetParLimits(2, sigma - sigmaE,   sigma + sigmaE);
		f1->SetParameter(3, frac2);         f1->SetParLimits(3, 0., 1.);
		f1->SetParameter(4, sigma2);	  f1->SetParLimits(4, sigma2 - sigma2E, sigma2 + sigma2E);
		f1->SetParameter(5, scale*pol0);  //  f1->SetParLimits(5, pol0 - pol0E, pol0 + pol0E);
		f1->SetParameter(6, scale*pol1);  //  f1->SetParLimits(6, pol1 - pol1E, pol1 + pol1E);
		f1->FixParameter(7, step);
		f1->FixParameter(8, res);
		f1->SetParameter(9, scale*level);	  f1->SetParLimits(9, scale*level*0.8, scale*level*1.2);
		fIF->dumpParameters(f1);
		h1->Fit(f1, "lr", "", xmin, xmax);
		par0E = f1->GetParError(0);
		intError = f1->IntegralError(peak-3*sigma, peak+3*sigma);
		psNormE = intError;
		for (int ipar = 0; ipar < fg->GetNpar(); ++ipar) {
		  fg->SetParameter(ipar, f1->GetParameter(ipar));
		  fg->SetParError(ipar, f1->GetParameter(ipar));
		}
		for (int ipar = 0; ipar < fe->GetNpar(); ++ipar) fe->SetParameter(ipar, f1->GetParameter(ipar+7));
		for (int ipar = 0; ipar < fp->GetNpar(); ++ipar) fp->SetParameter(ipar, f1->GetParameter(ipar+5));
		fg->Draw("same");
		f1->Draw("same");
		fe->Draw("same");
		fp->Draw("same");
		psNorm = fg->Integral(peak-3*sigma, peak+3*sigma);
		psNorm  /= h1->GetBinWidth(1);
	      }
	      psNormE /= h1->GetBinWidth(1);
	      if (psNormE < 0.01*psNorm) psNormE = par0E/h1->GetBinWidth(1);
	      if (psNormE > psNorm) psNormE = TMath::Sqrt(psNorm);
	      norm += ips*psNorm;
	      normE += (ips*psNormE)*(ips*psNormE);
	      cout << "prescale: " << ips << " Nsig = " << psNorm << " (area = " << fg->GetParameter(0)
		   << ") -> running sum: " << norm << " running error: " << TMath::Sqrt(normE) << endl;
	      if (TMath::Sqrt(normE) < 0.001) {
		cout << "XXXXXXXXX psNormE            = " << psNormE << endl;
		cout << "XXXXXXXXX fg->GetParError(0) = " << par0E << endl;
		cout << "XXXXXXXXX intError           = " << intError << endl;
		cout << "XXXXXXXXX sqrt(psNorm)       = " << TMath::Sqrt(psNorm)*h1->GetBinWidth(1) << endl;
	      }

	      tl->SetTextSize(0.04);
	      tl->DrawLatexNDC(0.5, 0.5, Form("Signal: %5.1f", psNorm));
	      tl->SetTextSize(0.03);
	      tl->DrawLatexNDC(0.2, 0.92, Form("yield-%s-%d-chan%d-ps%d.pdf", trg.c_str(), iblock, ichan, ips));
	      savePad(Form("yield-%s-%s-%d-chan%d-ps%d.pdf", trg.c_str(), dsname.c_str(), iblock, ichan, ips));
	    }

	    double result  = norm/blockLumi;
	    double resultE = TMath::Sqrt(normE)/blockLumi;
	    // -- normalize yield to 1/pb
	    cout << "==> Filling for chan = " << ichan << " into bin " << static_cast<double>(iblock)
		 << " result = " << result << " +/- " << resultE
		 << " for blockLumi = " << blockLumi
		 << endl;
	    vRunHLT[ichan]->SetBinContent(vRunHLT[ichan]->FindBin(static_cast<double>(iblock)), result);
	    vRunHLT[ichan]->SetBinError(vRunHLT[ichan]->FindBin(static_cast<double>(iblock)), resultE);
	  } else {
	    double allCnt(0.);
	    bool SKIP(true);
	    for (int ips = 1; ips < MAXPS; ++ips) {
	      h1 = vBlockHist[Form("chan%d_ps%d", ichan, ips)];
	      if (0 == h1) {
		cout << "xxx no histogram " << Form("chan%d_ps%d", ichan, ips) << endl;
		continue;
	      }
	      if (ips > 1) {
		if (h1->GetSumOfWeights() > 0) {
		  cout << "XXXXXXXXXXX BMM with prescale > 1 XXXXXXXXXXXXXXXXX" << endl;
		  SKIP = false;
		} else {
		  SKIP = true;
		}
	      } else {
		SKIP = false;
	      }
	      if (!SKIP) {
		double loCnt = h1->Integral(h1->FindBin(5.0), h1->FindBin(5.2));
		double hiCnt = h1->Integral(h1->FindBin(5.5), h1->FindBin(5.9));
		allCnt += ips*(loCnt + hiCnt);
		h1->Draw();
		tl->DrawLatexNDC(0.5, 0.5, Form("Signal: %5.1f+%5.1f = %5.1f", loCnt, hiCnt, allCnt));
		savePad(Form("yield-%s-%s-%d-chan%d-ps%d.pdf", trg.c_str(), dsname.c_str(), iblock, ichan, ips));
	      }
	    }
	    double result  = allCnt/blockLumi;
	    double resultE = TMath::Sqrt(allCnt)/blockLumi;
	    vRunHLT[ichan]->SetBinContent(vRunHLT[ichan]->FindBin(static_cast<double>(iblock)), result);
	    vRunHLT[ichan]->SetBinError(vRunHLT[ichan]->FindBin(static_cast<double>(iblock)), resultE);

	  }
	}
      }
      gStyle->SetOptStat(0);
      gStyle->SetOptFit(0);
      for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
	setTitles(vRunHLT[ichan], "run", Form("N(%s)", fDS[dsname]->fName.c_str()), 0.05, 1.1, 1.9);
	vRunHLT[ichan]->Draw();
	savePad(Form("yieldVsBlock-%s-%s-chan%d.pdf", trg.c_str(), dsname.c_str(), ichan));
	if (1) {
	  vRunHLT[ichan]->SetDirectory(gDirectory);
	  vRunHLT[ichan]->Write();
	}
      }
    }

  } else {
    TDirectory *hDir = fHistFile->mkdir(dir.c_str());
    fHistFile->cd(dir.c_str());
    hDir = gDirectory;
    cout << "created " << hDir->GetName() << endl;

    TTree *t = getTree(fSample, dir);
    if (0 == t) {
      cout << "tree for sample = " << fSample << " not found" << endl;
      return;
    }
    setupTree(t, fSample);
    fCds = fDS[fSample];
    loopOverTree(t, 1);

    cout << "writing output histograms: " << fYieldHLT.size() << endl;
    for (map<string, TH2D*>::iterator it = fYieldHLT.begin(); it != fYieldHLT.end(); ++it) {
      cout << "run " << it->first << endl;
      it->second->Draw("colz");
      it->second->SetDirectory(hDir);
      it->second->Write();
    }
    for (map<string, TH2D*>::iterator it = fYieldRTR.begin(); it != fYieldRTR.end(); ++it) {
      it->second->Draw("colz");
      it->second->SetDirectory(hDir);
      it->second->Write();
    }

    fYieldHLT.clear();
    fYieldRTR.clear();
  }

  //  fHistFile->Write();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotStuff::yieldStabilityRatios(string trgname) {

  vector<pair<string, string> > overlays;
  overlays.push_back(make_pair("bupsikData", "bspsiphiData"));
  overlays.push_back(make_pair("bupsikData", "bdpsikstarData"));
  overlays.push_back(make_pair("bmmData",    "bupsikData"));
  overlays.push_back(make_pair("bmmData",    "bspsiphiData"));

  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;

  string dir1(""), dir2("");
  TH1D *h1(0), *h2(0), *hR(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  tl->SetTextSize(0.03);
  c0->Clear();
  shrinkPad(0.15, 0.2);
  for (vector<pair<string, string> >::iterator it = overlays.begin(); it != overlays.end(); ++it) {
    if (string::npos != it->first.find("bupsik"))      dir1 = "candAnaBu2JpsiK";
    if (string::npos != it->first.find("bspsiphi"))    dir1 = "candAnaBs2JpsiPhi";
    if (string::npos != it->first.find("bdpsikstar"))  dir1 = "candAnaBd2JpsiKstar";
    if (string::npos != it->first.find("bmm"))         dir1 = "candAnaMuMu";

    if (string::npos != it->second.find("bupsik"))     dir2 = "candAnaBu2JpsiK";
    if (string::npos != it->second.find("bspsiphi"))   dir2 = "candAnaBs2JpsiPhi";
    if (string::npos != it->second.find("bdpsikstar")) dir2 = "candAnaBd2JpsiKstar";
    if (string::npos != it->second.find("bmm"))        dir2 = "candAnaMuMu";

    for (int ichan = 0; ichan < fNchan; ++ichan) {
      cout << "overlay " << it->first << " and " << it->second << " chan " << ichan << endl;
      h1 = (TH1D*)fHistFile->Get(Form("%s/hRun%s_%s_chan%d", dir1.c_str(), trgname.c_str(), it->first.c_str(), ichan));
      h2 = (TH1D*)fHistFile->Get(Form("%s/hRun%s_%s_chan%d", dir2.c_str(), trgname.c_str(), it->second.c_str(), ichan));
      cout << "h1 = " << h1 << " h2 = " << h2 << endl;
      if (h1 && h2) {
	hR = (TH1D*)h1->Clone(Form("hr_%s_%s", h1->GetName(), h2->GetName()));
	hR->Clear();
	setTitles(hR, "run", Form("N(%s) / N(%s)", fDS[it->first]->fName.c_str(), fDS[it->second]->fName.c_str()), 0.05, 1.1, 1.9);
	hR->Divide(h1, h2);
	hR->Fit("pol1");
	tl->DrawLatexNDC(0.2, 0.92, Form("p0 = %4.3f#pm%4.3f ",
					 hR->GetFunction("pol1")->GetParameter(0),
					 hR->GetFunction("pol1")->GetParError(0)
					 ));
	tl->DrawLatexNDC(0.45, 0.92, Form("p1 = %4.3f#pm%4.3f (x 1e6)",
					  1.e6*hR->GetFunction("pol1")->GetParameter(1),
					  1.e6*hR->GetFunction("pol1")->GetParError(1)
					  ));
	tl->DrawLatexNDC(0.80, 0.92, Form("chan %d", ichan));

	savePad(Form("yieldStabilityRatio-%s-%s-chan%d.pdf", it->first.c_str(), it->second.c_str(), ichan));
      }
    }

  }



}



// ----------------------------------------------------------------------
void plotStuff::loopFunction1() {

  if (!fGoodMuonsID) return;

  if (!fGoodQ) return;
  if (!fGoodPvAveW8) return;
  if (!fGoodMaxDoca) return;

  if (!fb.json) return;

  if (fb.flsxy    < fCuts[fChan]->flsxy) return;
  if (fb.fls3d    < fCuts[fChan]->fls3d) return;

  if (fb.chi2dof  > fCuts[fChan]->chi2dof) return;
  if (fb.alpha    > fCuts[fChan]->alpha) return;
  if (fb.pvip     > fCuts[fChan]->pvip) return;
  if (fb.pvips    > fCuts[fChan]->pvips) return;

  if (fb.iso      < fCuts[fChan]->iso) return;
  if (fb.docatrk  < fCuts[fChan]->docatrk) return;
  if (fb.closetrk > fCuts[fChan]->closetrk) return;

  if (fb.m1pt < 4.0) return;
  if (fb.m2pt < 4.0) return;

  double m = fb.m;
  if ((fMode == BU2JPSIKP) || (fMode == BD2JPSIKSTAR) || (fMode == BS2JPSIPHI)) {
    if (TMath::Abs(fb.mpsi) < 2.9) return;
    if (TMath::Abs(fb.mpsi) > 3.3) return;

    if (TMath::Abs(fb.psipt) < 6.9) return;
    if (TMath::Abs(fb.psicosa) < 0.9) return;
    if (TMath::Abs(fb.psiprob) < 0.1) return;
    if (TMath::Abs(fb.psiflsxy) < 4) return;
    if (fMode == BS2JPSIPHI) {
      if (fb.mkk   < 1.01) return;
      if (fb.mkk   > 1.03) return;
      if (fb.phidr > 0.30) return;
      if (fb.k1pt  < 0.80) return;
      if (fb.k2pt  < 0.80) return;
    }

    if (fMode == BU2JPSIKP) {
      if (fb.kpt < 0.70) return;
    }

    if (fMode == BD2JPSIKSTAR) {
      if (fb.kpt < 0.70) return;
      if (fb.pipt < 0.70) return;
      if (fb.mkpi < 0.86) return;
      if (fb.mkpi > 0.94) return;
    }

    m = fb.cm;
  }


  if (0 == fYieldHLT.count(Form("%d_chan%d", static_cast<int>(fb.run), fChan))) {
    TH2D *h = new TH2D(Form("h_HLT_%d_chan%d", static_cast<int>(fb.run), fChan), Form("run%d chan%d", fb.run, fChan), 90, 5.0, 5.9, MAXPS+1, -1., MAXPS);
    fYieldHLT.insert(make_pair(Form("%d_chan%d", static_cast<int>(fb.run), fChan), h));

    h = new TH2D(Form("h_RTR_%d_chan%d", static_cast<int>(fb.run), fChan), Form("run%d chan%d", fb.run, fChan), 90, 5.0, 5.9, MAXPS+1, -1., MAXPS);
    fYieldRTR.insert(make_pair(Form("%d_chan%d", static_cast<int>(fb.run), fChan), h));
  }


  if (fb.hlt1 && fb.tos) {
    fYieldHLT[Form("%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, -0.1, static_cast<double>(fb.ps));
    fYieldHLT[Form("%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, 0.1);
    fYieldHLT[Form("%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, fb.ps+0.1);
  }

  if (fb.reftrg) {
    fYieldRTR[Form("%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, -0.1, static_cast<double>(fb.ps));
    fYieldRTR[Form("%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, 0.1);
    fYieldRTR[Form("%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, fb.ps+0.1);
  }

}


// ----------------------------------------------------------------------
// -- this is for the puStudy!!!
void plotStuff::loopFunction2() {

  int ic = fpv.chan;
  int idx[] = {4, ic};
  int imax = (ic > -1?2:1);
  // -- we are only interested in 'real' candidates (else you get strange effects, e.g. back-to-back low-pT muons with large eta-B)
  if (!fpv.hlt1) return;

  if (fpv.m1pt < 4.0) return;
  if (fpv.m2pt < 4.0) return;

  // -- fill per-channel and combined histograms
  for (int i = 0; i < imax; ++i) {
    if (TMath::Abs(fpv.dzmin) > 1.0) fpHmultFar[idx[i]]->Fill(fpv.mult1);
    if (TMath::Abs(fpv.dzmin) < 0.05) fpHmultClose05[idx[i]]->Fill(fpv.mult1);

    // -- vs dzmin
    fpHdzmin[idx[i]]->Fill(TMath::Abs(fpv.dzmin));
    fpP1Mult[idx[i]]->Fill(fpv.dzmin, fpv.mult1);
    fpP1flsxy[idx[i]]->Fill(fpv.dzmin, fpv.flsxy);
    fpP1fls3d[idx[i]]->Fill(fpv.dzmin, fpv.fls3d);
    fpP1fl3d[idx[i]]->Fill(fpv.dzmin, fpv.fl3d);
    fpP1dfl3d[idx[i]]->Fill(fpv.dzmin, fpv.fl3d - fpv.gfl);

    fpP1tau[idx[i]]->Fill(fpv.dzmin,  fpv.t1);
    fpP1dtau[idx[i]]->Fill(fpv.dzmin, fpv.t1 - fpv.gt);

    // -- vs lz2
    if (TMath::Abs(fpv.lz2) > 1.0) fpHmultFar2[idx[i]]->Fill(fpv.mult1);
    if (TMath::Abs(fpv.lz2) < 0.05) fpHmultClose2[idx[i]]->Fill(fpv.mult1);

    fpHlz1[idx[i]]->Fill(TMath::Abs(fpv.lz1));
    if (1) {
      fpHtmlz1[idx[i]]->Fill(TMath::Abs(fpv.lz1));
    }
    fpHlz2[idx[i]]->Fill(TMath::Abs(fpv.lz2));
    fpP2Mult[idx[i]]->Fill(fpv.lz2, fpv.mult1);
    fpP2flsxy[idx[i]]->Fill(fpv.lz2, fpv.flsxy);
    fpP2fls3d[idx[i]]->Fill(fpv.lz2, fpv.fls3d);
    fpP2fl3d[idx[i]]->Fill(fpv.lz2, fpv.fl3d);
    fpP2dfl3d[idx[i]]->Fill(fpv.lz2, fpv.fl3d - fpv.gfl);

    fpP2tau[idx[i]]->Fill(fpv.lz2,  fpv.t1);
    fpP2dtau[idx[i]]->Fill(fpv.lz2, fpv.t1 - fpv.gt);

    // -- cutting on fls3d > 10
    if (fpv.fls3d > 5) {
      fpP3tau[idx[i]]->Fill(fpv.dzmin,  fpv.t1);
      fpP3dtau[idx[i]]->Fill(fpv.dzmin,  fpv.t1 - fpv.gt);
    }
  }
}



void plotStuff::loopFunction3() { }
void plotStuff::loopFunction4() { }

// ----------------------------------------------------------------------
void plotStuff::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotStuff::loopOverTree> loop over dataset " << fCds->fName << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotStuff::*pF)(void);
  if (ifunc == 1) pF = &plotStuff::loopFunction1;
  if (ifunc == 2) pF = &plotStuff::loopFunction2;
  if (ifunc == 3) pF = &plotStuff::loopFunction3;
  if (ifunc == 4) pF = &plotStuff::loopFunction4;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotStuff::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotStuff::loadFile loading files listed in " << files << endl;

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

      if (string::npos != stype.find("wrongreco,")) {
        sname = "wrongReco";
        sdecay = "wrongReco";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

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


// ----------------------------------------------------------------------
void plotStuff::setupPvTree(TTree *t) {

  t->SetBranchAddress("m1pt", &fpv.m1pt);
  t->SetBranchAddress("m2pt", &fpv.m2pt);
  t->SetBranchAddress("pt",   &fpv.pt);
  t->SetBranchAddress("eta",  &fpv.eta);
  t->SetBranchAddress("phi",  &fpv.phi);
  t->SetBranchAddress("m",    &fpv.m);
  t->SetBranchAddress("chan", &fpv.chan);
  t->SetBranchAddress("hlt1", &fpv.hlt1);
  t->SetBranchAddress("npv",  &fpv.npv);
  t->SetBranchAddress("idx1", &fpv.idx1);
  t->SetBranchAddress("idx2", &fpv.idx2);
  t->SetBranchAddress("idx3", &fpv.idx3);
  t->SetBranchAddress("lz1",  &fpv.lz1);
  t->SetBranchAddress("lz2",  &fpv.lz2);
  t->SetBranchAddress("mult1",&fpv.mult1);
  t->SetBranchAddress("mult2",&fpv.mult2);
  t->SetBranchAddress("prob1",&fpv.prob1);
  t->SetBranchAddress("prob2",&fpv.prob2);
  t->SetBranchAddress("chi1", &fpv.chi1);
  t->SetBranchAddress("chi2", &fpv.chi2);
  t->SetBranchAddress("dz12", &fpv.dz12);
  t->SetBranchAddress("dzmin",&fpv.dzmin);
  t->SetBranchAddress("gfl",   &fpv.gfl);
  t->SetBranchAddress("gt",    &fpv.gt);
  t->SetBranchAddress("flsxy", &fpv.flsxy);
  t->SetBranchAddress("flxy",  &fpv.flxy);
  t->SetBranchAddress("fls3d", &fpv.fls3d);
  t->SetBranchAddress("fl3d",  &fpv.fl3d);
  t->SetBranchAddress("fl1",   &fpv.fl1);
  t->SetBranchAddress("fl2",   &fpv.fl2);
  t->SetBranchAddress("fl3",   &fpv.fl3);
  t->SetBranchAddress("t1",    &fpv.t1);
  t->SetBranchAddress("t2",    &fpv.t2);
  t->SetBranchAddress("t3",    &fpv.t3);
  t->SetBranchAddress("gs",    &fpv.gs);
  t->SetBranchAddress("s1",    &fpv.s1);
  t->SetBranchAddress("s2",    &fpv.s2);
  t->SetBranchAddress("s3",    &fpv.s3);
  t->SetBranchAddress("gx",  &fpv.gx);
  t->SetBranchAddress("gy",  &fpv.gy);
  t->SetBranchAddress("gz",  &fpv.gz);
  t->SetBranchAddress("p1x",  &fpv.p1x);
  t->SetBranchAddress("p1y",  &fpv.p1y);
  t->SetBranchAddress("p1z",  &fpv.p1z);
  t->SetBranchAddress("p1d",  &fpv.p1d);
  t->SetBranchAddress("p2x",  &fpv.p2x);
  t->SetBranchAddress("p2y",  &fpv.p2y);
  t->SetBranchAddress("p2z",  &fpv.p2z);
  t->SetBranchAddress("p2d",  &fpv.p2d);

  t->SetBranchAddress("p3x",  &fpv.p3x);
  t->SetBranchAddress("p3y",  &fpv.p3y);
  t->SetBranchAddress("p3z",  &fpv.p3z);

  t->SetBranchAddress("sx",  &fpv.sx);
  t->SetBranchAddress("sy",  &fpv.sy);
  t->SetBranchAddress("sz",  &fpv.sz);
  t->SetBranchAddress("dsv", &fpv.dsv);

  t->SetBranchAddress("d1",  &fpv.d1);
  t->SetBranchAddress("a1",  &fpv.a1);
  t->SetBranchAddress("d2",  &fpv.d2);
  t->SetBranchAddress("a2",  &fpv.a2);
  t->SetBranchAddress("d3",  &fpv.d3);
  t->SetBranchAddress("a3",  &fpv.a3);


}
