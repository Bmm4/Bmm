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
#include "TGaxis.h"
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
#include "common/fitPsYield.hh"

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

  if (2011 == fYear) {
    fLargeRuns.push_back(165993);
    fLargeRuns.push_back(166950);
    fLargeRuns.push_back(167898);
    fLargeRuns.push_back(172822);
    fLargeRuns.push_back(172868);
    fLargeRuns.push_back(177730);
    fLargeRuns.push_back(178100);
  } else if (2012 == fYear) {
    fLargeRuns.push_back(194050);
    fLargeRuns.push_back(194912);
    fLargeRuns.push_back(195552);
    fLargeRuns.push_back(199435);
    fLargeRuns.push_back(201191);
    fLargeRuns.push_back(201278);
    fLargeRuns.push_back(202178);
    fLargeRuns.push_back(203002);
    fLargeRuns.push_back(205344);
    fLargeRuns.push_back(206745);
    fLargeRuns.push_back(207099);
    fLargeRuns.push_back(207231);
    fLargeRuns.push_back(207454);
    fLargeRuns.push_back(207905);
  } else if (2016 == fYear) {
    fLargeRuns.push_back(273725); // B
    fLargeRuns.push_back(274241);
    fLargeRuns.push_back(274316);
    fLargeRuns.push_back(274335);
    fLargeRuns.push_back(274388);
    fLargeRuns.push_back(274422);
    fLargeRuns.push_back(274968);
    fLargeRuns.push_back(275001);
    fLargeRuns.push_back(275125);
    fLargeRuns.push_back(275310);
    fLargeRuns.push_back(275376);
    fLargeRuns.push_back(275782); // C
    fLargeRuns.push_back(275836);
    fLargeRuns.push_back(275847);         // very high
    fLargeRuns.push_back(275890);
    fLargeRuns.push_back(276242);         // very low
    fLargeRuns.push_back(276282);
    fLargeRuns.push_back(276363);
    fLargeRuns.push_back(276437); // D
    fLargeRuns.push_back(276501);
    fLargeRuns.push_back(276525);
    fLargeRuns.push_back(276542);
    fLargeRuns.push_back(276581);
    fLargeRuns.push_back(276655);
    fLargeRuns.push_back(276776);
    fLargeRuns.push_back(276811);
    fLargeRuns.push_back(276831);
    fLargeRuns.push_back(276870); // E
    fLargeRuns.push_back(276950);
    fLargeRuns.push_back(277096);
    fLargeRuns.push_back(277194);
    fLargeRuns.push_back(277305);        // catastrophic run
    fLargeRuns.push_back(278167); // F
    fLargeRuns.push_back(278308);
    fLargeRuns.push_back(278406);
    fLargeRuns.push_back(278509);
    fLargeRuns.push_back(278808);
    fLargeRuns.push_back(278820); // G
    fLargeRuns.push_back(278822);
    fLargeRuns.push_back(278969);
    fLargeRuns.push_back(279694);
    fLargeRuns.push_back(279931);
    fLargeRuns.push_back(279975);
    fLargeRuns.push_back(280249);
    fLargeRuns.push_back(280385);
    fLargeRuns.push_back(281693); // H
    fLargeRuns.push_back(281797);
    fLargeRuns.push_back(281976);
    fLargeRuns.push_back(282037);
    fLargeRuns.push_back(282092);
    fLargeRuns.push_back(282735);
    fLargeRuns.push_back(282814);
    fLargeRuns.push_back(283270);
    fLargeRuns.push_back(283408);
    fLargeRuns.push_back(283478);
    fLargeRuns.push_back(283865);
    fLargeRuns.push_back(283946);
  }

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

  if (what == "runs") {
    runStudy("bupsikData");
    runStudy("bmmData");
  }

  if (what == "ysfill") {
    yieldStability("bupsikData", "ysfill");
    yieldStability("bmmData", "ysfill");
    yieldStability("bspsiphiData", "ysfill");
  }
  if (what == "ysplot") {
    yieldStability("bupsikData", "NOC");
    yieldStability("bupsikData", "TOS");
    yieldStability("bmmData", "NOC");
    yieldStability("bmmData", "TOS");
    yieldStabilityRatios("NOC");
    yieldStabilityRatios("TOS");
  }

  if (what == "all" || what == "massresolution") {
    massResolution("bsmmMcComb", "bsmmMcRun1");
  }
  if (what == "all" || what == "tauefficiency") {
    tauEfficiency("all", "", "", "");
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
    puStudy("bmmData");
    puStudy("bsmmMcOff", "bmmData");
    puStudy("bupsikMcOff");
    puStudy("bupsikData");
  }


  if (string::npos != what.find("wrongreco")) {
    wrongReco("wrongReco", "candAnaBd2JpsiKstarAsBu", "hlt");
    wrongReco("wrongReco", "candAnaBd2JpsiKstarAsBs", "1.01 < mkk && mkk < 1.03 && k1pt > 0.7 && k2pt > 0.7");
    wrongReco("bcpsimunuMcCombBg", "candAnaMuMu", "");
    wrongReco("bupsipiMcRun1", "candAnaBu2JpsiK", "gmuid&&abs(m1eta)<1.4&&fabs(m2eta)<1.4&&m2pt>4&&fls3d>10&&iso>0.8&&alpha<0.05");

    plotWrongReco("alpha", 50, 0., 0.1, "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
    plotWrongReco("closetrk", 20, 0., 20., "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
    plotWrongReco("iso", 51, 0., 1.01, "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
    plotWrongReco("pvip", 50, 0., 0.02, "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
    plotWrongReco("pvips", 50, 0., 4., "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
    plotWrongReco("closetrks1", 20, 0., 20., "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
    plotWrongReco("closetrks2", 20, 0., 20., "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
    plotWrongReco("closetrks3", 20, 0., 20., "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
    plotWrongReco("docatrk", 100, 0., 0.02, "", "wrongReco", "candAnaBd2JpsiKstarAsBu", "bupsikMcComb", "candAnaBu2JpsiK");
  }

}



// ----------------------------------------------------------------------
void plotStuff::bookHist(string dsname) {


}


// ----------------------------------------------------------------------
void plotStuff::puStudy(string dsname, std::string dsname2) {

  char mode[200]; sprintf(mode, "%s", dsname.c_str());
  bool fillHist(false);
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");

  // ----------------------------------------------------------------------
  // -- Get or book histograms and profiles
  // ----------------------------------------------------------------------
  for (int i = 0; i < NCHAN; ++i) {
    if (!(fpHmultFar[i] = (TH1D*)fHistFile->Get(Form("pvmultfar_%s_chan%d", mode, i)))) {
      fillHist = true;
      fpHmultFar[i]   = new TH1D(Form("pvmultfar_%s_chan%d", mode, i), "PV multiplicity dzmin > 2 cm", 40, 0., 120.);
    }
    setFilledHist(fpHmultFar[i], kBlue, kBlue, 3354);
    setTitles(fpHmultFar[i], "PV track multiplicity", "a.u.", 0.05, 1.1, 1.8);

    if (!(fpHmultClose05[i] = (TH1D*)fHistFile->Get(Form("pvmultclose05_%s_chan%d", mode, i))))
      fpHmultClose05[i] = new TH1D(Form("pvmultclose05_%s_chan%d", mode, i), "PV multiplicity dzmin < 0.05 cm", 40, 0., 120.);
    setFilledHist(fpHmultClose05[i], kGreen+2, kGreen+2, 3365);

    if (!(fpHmultFar2[i] = (TH1D*)fHistFile->Get(Form("pvmultfar2_%s_chan%d", mode, i))))
      fpHmultFar2[i]  = new TH1D(Form("pvmultfar2_%s_chan%d", mode, i), "PV multiplicity | l_{z}^{2} | > 1 cm", 40, 0., 120.);
    setFilledHist(fpHmultFar2[i], kBlue, kBlue, 3354);
    setTitles(fpHmultFar2[i], "PV track multiplicity", "a.u.", 0.05, 1.1, 1.8);
    if (!(fpHmultClose2[i] = (TH1D*)fHistFile->Get(Form("pvmultclose2_%s_chan%d", mode, i))))
      fpHmultClose2[i] = new TH1D(Form("pvmultclose2_%s_chan%d", mode, i), "PV multiplicity  | l_{z}^{2} | < 0.05 cm", 40, 0., 120.);
    setFilledHist(fpHmultClose2[i], kGreen+2, kGreen+2, 3365);

    if (!(fpHpvn[i] = (TH1D*)fHistFile->Get(Form("pvn_%s_chan%d", mode, i))))
      fpHpvn[i]   = new TH1D(Form("pvn_%s_chan%d ", mode, i), "PV multiplicity", 40, 0., 40.);
    setFilledHist(fpHpvn[i], kBlue, kYellow, 1000);
    setTitles(fpHpvn[i], "PV multiplicity", "a.u.", 0.05, 1.1, 1.8);

    if (!(fpHdzmin[i] = (TH1D*)fHistFile->Get(Form("dzmin_%s_chan%d", mode, i))))
      fpHdzmin[i] = new TH1D(Form("dzmin_%s_chan%d", mode, i), "dzmin", 200, -2., 2.0);
    setFilledHist(fpHdzmin[i], kBlue, kYellow, 1000);
    setTitles(fpHdzmin[i], "minimum #Delta z #it{[cm]}", "a.u.", 0.05, 1.1, 1.8);

    if (!(fpHlz1[i] = (TH1D*)fHistFile->Get(Form("lz1_%s_chan%d", mode, i))))
      fpHlz1[i] = new TH1D(Form("lz1_%s_chan%d", mode, i), "lz1", 400, -2.0, 2.0);
    setTitles(fpHlz1[i], "#it{|} l_{z}^{(1)} #it{| [cm]}", "a.u.", 0.05, 1.1, 1.8);
    if (!(fpHlzs1[i] = (TH1D*)fHistFile->Get(Form("lzs1_%s_chan%d", mode, i))))
      fpHlzs1[i] = new TH1D(Form("lzs1_%s_chan%d", mode, i), "lzs1", 200, -100., 100);
    setTitles(fpHlzs1[i], "l_{z}^{(1)} #sigma(l_{z}^{(1)}) ", "a.u.", 0.05, 1.1, 1.8);

    if (!(fpHlz2[i] = (TH1D*)fHistFile->Get(Form("lz2_%s_chan%d", mode, i))))
      fpHlz2[i] = new TH1D(Form("lz2_%s_chan%d", mode, i), "lz2", 400, -2.0, 2.0);
    setTitles(fpHlz2[i], "#it{|} l_{z}^{(2)} #it{| [cm]}", "a.u.", 0.05, 1.1, 1.8);
    if (!(fpHlzs2[i] = (TH1D*)fHistFile->Get(Form("lzs2_%s_chan%d", mode, i))))
      fpHlzs2[i] = new TH1D(Form("lzs2_%s_chan%d", mode, i), "lzs2", 200, -100., 100);
    setTitles(fpHlzs2[i], "l_{z}^{(2)} / #sigma(l_{z}^{(2)})", "a.u.", 0.05, 1.1, 1.8);


    // -- vs dzmin
    if (!(fpP0Mult[i] = (TProfile*)fHistFile->Get(Form("p0mult_%s_chan%d", mode, i))))
      fpP0Mult[i] = new TProfile(Form("p0mult_%s_chan%d", mode, i), "mult", 40, -1.0, 1.0, "");
    fpP0Mult[i]->SetMinimum(0.);     fpP0Mult[i]->SetMaximum(80.);
    setTitles(fpP0Mult[i], "minimum #it{|}#Delta z#it{| [cm]}", "PV trk mult", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP0fls3d[i] = (TProfile*)fHistFile->Get(Form("p0fls3d_%s_chan%d", mode, i))))
      fpP0fls3d[i] = new TProfile(Form("p0fls3d_%s_chan%d", mode, i), "fls3d", 40, -1.0, 1.0, "");
    setTitles(fpP0fls3d[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean fls3d", 0.05, 1.1, 1.8, 0.04);
    fpP0fls3d[i]->SetMinimum(0.);     fpP0fls3d[i]->SetMaximum(50.);
    if (!(fpP0fl3d[i] = (TProfile*)fHistFile->Get(Form("p0fl3d_%s_chan%d", mode, i))))
      fpP0fl3d[i] = new TProfile(Form("p0fl3d_%s_chan%d", mode, i), "fl3d", 40, -1.0, 1.0, "");
    setTitles(fpP0fl3d[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean fl3d", 0.05, 1.1, 1.8, 0.04);
    fpP0fl3d[i]->SetMinimum(0.);     fpP0fl3d[i]->SetMaximum(0.5);
    if (!(fpP0dfl3d[i] = (TProfile*)fHistFile->Get(Form("p0dfl3d_%s_chan%d", mode, i))))
      fpP0dfl3d[i] = new TProfile(Form("p0dfl3d_%s_chan%d", mode, i), "dfl3d", 40, -1.0, 1.0, "");
    fpP0dfl3d[i]->SetMinimum(-0.03); fpP0dfl3d[i]->SetMaximum(0.03);
    setTitles(fpP0dfl3d[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean #Delta fl3d", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP0Iso[i] = (TProfile*)fHistFile->Get(Form("p0Iso_%s_chan%d", mode, i))))
      fpP0Iso[i] = new TProfile(Form("p0Iso_%s_chan%d", mode, i), "Iso", 40, -1.0, 1.0, "");
    fpP0Iso[i]->SetMinimum(0.3); fpP0Iso[i]->SetMaximum(1.01);
    setTitles(fpP0Iso[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean Isolation", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP0Ntrk[i] = (TProfile*)fHistFile->Get(Form("p0Ntrk_%s_chan%d", mode, i))))
      fpP0Ntrk[i] = new TProfile(Form("p0Ntrk_%s_chan%d", mode, i), "Ntrk", 40, -1.0, 1.0, "");
    fpP0Ntrk[i]->SetMinimum(0.); fpP0Ntrk[i]->SetMaximum(20.);
    setTitles(fpP0Ntrk[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean closetrk", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP0PvIps[i] = (TProfile*)fHistFile->Get(Form("p0PvIps_%s_chan%d", mode, i))))
      fpP0PvIps[i] = new TProfile(Form("p0PvIps_%s_chan%d", mode, i), "PvIps", 40, -1.0, 1.0, "");
    fpP0PvIps[i]->SetMinimum(0.); fpP0PvIps[i]->SetMaximum(5.);
    setTitles(fpP0PvIps[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean PvIps", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP0Bdt[i] = (TProfile*)fHistFile->Get(Form("p0Bdt_%s_chan%d", mode, i))))
      fpP0Bdt[i] = new TProfile(Form("p0Bdt_%s_chan%d", mode, i), "Bdt", 40, -1.0, 1.0, "");
    fpP0Bdt[i]->SetMinimum(-1.); fpP0Bdt[i]->SetMaximum(1.);
    setTitles(fpP0Bdt[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean Bdt", 0.05, 1.1, 1.8, 0.04);

    if (!(fpP0l1[i] = (TProfile*)fHistFile->Get(Form("p0l1_%s_chan%d", mode, i))))
      fpP0l1[i] = new TProfile(Form("p0l1_%s_chan%d", mode, i), "lz1", 40, -1.0, 1.0, "");
    fpP0l1[i]->SetMinimum(-0.1); fpP0l1[i]->SetMaximum(0.1);
    setTitles(fpP0l1[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean l_{z}^{(1)} [cm]", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP0l2[i] = (TProfile*)fHistFile->Get(Form("p0l2_%s_chan%d", mode, i))))
      fpP0l2[i] = new TProfile(Form("p0l2_%s_chan%d", mode, i), "lz2", 40, -1.0, 1.0, "");
    fpP0l2[i]->SetMinimum(-2.); fpP0l2[i]->SetMaximum(2.);
    setTitles(fpP0l2[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean l_{z}^{(2)} [cm]", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP0npv[i] = (TProfile*)fHistFile->Get(Form("p0npv_%s_chan%d", mode, i))))
      fpP0npv[i] = new TProfile(Form("p0npv_%s_chan%d", mode, i), "npv", 40, -1.0, 1.0, "");
    fpP0npv[i]->SetMinimum(0.); fpP0npv[i]->SetMaximum(40.);
    setTitles(fpP0npv[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean N_{PV}", 0.05, 1.1, 1.8, 0.04);

    if (!(fpP0tau[i] = (TProfile*)fHistFile->Get(Form("p0tau_%s_chan%d", mode, i))))
      fpP0tau[i] = new TProfile(Form("p0tau_%s_chan%d", mode, i), "tau", 40, -1.0, 1.0, "");
    fpP0tau[i]->SetMinimum(0.e-12); fpP0tau[i]->SetMaximum(+2.e-12);
    setTitles(fpP0tau[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean #tau_{eff} ", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP0dtau[i] = (TProfile*)fHistFile->Get(Form("p0dtau_%s_chan%d", mode, i))))
      fpP0dtau[i] = new TProfile(Form("p0dtau_%s_chan%d", mode, i), "dtau", 40, -1.0, 1.0, "");
    fpP0dtau[i]->SetMinimum(-0.2e-12); fpP0dtau[i]->SetMaximum(+0.2e-12);
    setTitles(fpP0dtau[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean #Delta#tau_{eff}", 0.05, 1.1, 1.8, 0.04);

    // -- vs lz
    if (!(fpP1Mult[i] = (TProfile*)fHistFile->Get(Form("p1Mult_%s_chan%d", mode, i))))
      fpP1Mult[i] = new TProfile(Form("p1mult_%s_chan%d", mode, i), "mult", 40, -0.03, 0.030, "");
    fpP1Mult[i]->SetMinimum(0.);     fpP1Mult[i]->SetMaximum(80.);
    setTitles(fpP1Mult[i], "l_{z}^{(1)} #it{ [cm]}", "PV trk mult", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP1fls3d[i] = (TProfile*)fHistFile->Get(Form("p1fls3d_%s_chan%d", mode, i))))
      fpP1fls3d[i] = new TProfile(Form("p1fls3d_%s_chan%d", mode, i), "fls3d", 40, -0.03, 0.03, "");
    setTitles(fpP1fls3d[i], "l_{z}^{(1)} #it{ [cm]}", "mean fls3d", 0.05, 1.1, 1.8, 0.04);
    fpP1fls3d[i]->SetMinimum(0.);     fpP1fls3d[i]->SetMaximum(50.);
    if (!(fpP1fl3d[i] = (TProfile*)fHistFile->Get(Form("p1fl3d_%s_chan%d", mode, i))))
      fpP1fl3d[i] = new TProfile(Form("p1fl3d_%s_chan%d", mode, i), "fl3d", 40, -0.03, 0.03, "");
    setTitles(fpP1fl3d[i], "l_{z}^{(1)} #it{ [cm]}", "mean fl3d", 0.05, 1.1, 1.8, 0.04);
    fpP1fl3d[i]->SetMinimum(0.);     fpP1fl3d[i]->SetMaximum(0.5);
    if (!(fpP1dfl3d[i] = (TProfile*)fHistFile->Get(Form("p1dfl3d_%s_chan%d", mode, i))))
      fpP1dfl3d[i] = new TProfile(Form("p1dfl3d_%s_chan%d", mode, i), "dfl3d", 40, -0.03, 0.03, "");
    fpP1dfl3d[i]->SetMinimum(-0.03); fpP1dfl3d[i]->SetMaximum(0.03);
    setTitles(fpP1dfl3d[i], "l_{z}^{(1)} #it{ [cm]}", "mean #Delta fl3d", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP1Iso[i] = (TProfile*)fHistFile->Get(Form("p1Iso_%s_chan%d", mode, i))))
      fpP1Iso[i] = new TProfile(Form("p1Iso_%s_chan%d", mode, i), "Iso", 40, -0.03, 0.03, "");
    fpP1Iso[i]->SetMinimum(0.3); fpP1Iso[i]->SetMaximum(1.01);
    setTitles(fpP1Iso[i], "l_{z}^{(1)} #it{ [cm]}", "mean Isolation", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP1Ntrk[i] = (TProfile*)fHistFile->Get(Form("p1Ntrk_%s_chan%d", mode, i))))
      fpP1Ntrk[i] = new TProfile(Form("p1Ntrk_%s_chan%d", mode, i), "Ntrk", 40, -0.03, 0.03, "");
    fpP1Ntrk[i]->SetMinimum(0.); fpP1Ntrk[i]->SetMaximum(20.);
    setTitles(fpP1Ntrk[i], "l_{z}^{(1)} #it{ [cm]}", "mean closetrk", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP1PvIps[i] = (TProfile*)fHistFile->Get(Form("p1PvIps_%s_chan%d", mode, i))))
      fpP1PvIps[i] = new TProfile(Form("p1PvIps_%s_chan%d", mode, i), "PvIps", 40, -0.03, 0.03, "");
    fpP1PvIps[i]->SetMinimum(0.); fpP1PvIps[i]->SetMaximum(5.);
    setTitles(fpP1PvIps[i], "l_{z}^{(1)} #it{ [cm]}", "mean PvIps", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP1Bdt[i] = (TProfile*)fHistFile->Get(Form("p1Bdt_%s_chan%d", mode, i))))
      fpP1Bdt[i] = new TProfile(Form("p1Bdt_%s_chan%d", mode, i), "Bdt", 40, -0.03, 0.03, "");
    fpP1Bdt[i]->SetMinimum(-1.); fpP1Bdt[i]->SetMaximum(1.);
    setTitles(fpP1Bdt[i], "l_{z}^{(1)} #it{ [cm]}", "mean Bdt", 0.05, 1.1, 1.8, 0.04);

    if (!(fpP1dzmin[i] = (TProfile*)fHistFile->Get(Form("p1dzmin_%s_chan%d", mode, i))))
      fpP1dzmin[i] = new TProfile(Form("p1dzmin_%s_chan%d", mode, i), "lz1", 40, -0.03, 0.03, "");
    fpP1dzmin[i]->SetMinimum(-1.); fpP1dzmin[i]->SetMaximum(1.);
    setTitles(fpP1dzmin[i], "l_{z}^{(1)} #it{[cm]}", "mean minimum #Delta z [cm]", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP1l2[i] = (TProfile*)fHistFile->Get(Form("p1l2_%s_chan%d", mode, i))))
      fpP1l2[i] = new TProfile(Form("p1l2_%s_chan%d", mode, i), "lz2", 40, -0.03, 0.03, "");
    fpP1l2[i]->SetMinimum(-1.); fpP1l2[i]->SetMaximum(1.);
    setTitles(fpP1l2[i], "l_{z}^{(1)} #it{[cm]}", "mean l_{z}^{(2)} [cm]", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP1npv[i] = (TProfile*)fHistFile->Get(Form("p1npv_%s_chan%d", mode, i))))
      fpP1npv[i] = new TProfile(Form("p1npv_%s_chan%d", mode, i), "npv", 40, -0.03, 0.03, "");
    fpP1npv[i]->SetMinimum(0.); fpP1npv[i]->SetMaximum(40.);
    setTitles(fpP1npv[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean N_{PV}", 0.05, 1.1, 1.8, 0.04);

    if (!(fpP1tau[i] = (TProfile*)fHistFile->Get(Form("p1tau_%s_chan%d", mode, i))))
      fpP1tau[i] = new TProfile(Form("p1tau_%s_chan%d", mode, i), "tau", 40, -0.03, 0.03, "");
    fpP1tau[i]->SetMinimum(0.e-12); fpP1tau[i]->SetMaximum(+2.e-12);
    setTitles(fpP1tau[i], "l_{z}^{(1)} #it{ [cm]}", "mean #tau_{eff} ", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP1dtau[i] = (TProfile*)fHistFile->Get(Form("p1dtau_%s_chan%d", mode, i))))
      fpP1dtau[i] = new TProfile(Form("p1dtau_%s_chan%d", mode, i), "dtau", 40, -0.03, 0.03, "");
    fpP1dtau[i]->SetMinimum(-0.2e-12); fpP1dtau[i]->SetMaximum(+0.2e-12);
    setTitles(fpP1dtau[i], "l_{z}^{(1)} #it{ [cm]}", "mean #Delta#tau_{eff}", 0.05, 1.1, 1.8, 0.04);

    // -- vs lz2
    if (!(fpP2Mult[i] = (TProfile*)fHistFile->Get(Form("p2mult_%s_chan%d", mode, i))))
      fpP2Mult[i] = new TProfile(Form("p2mult_%s_chan%d", mode, i), "mult", 40, -2., 2., "");
    fpP2Mult[i]->SetMinimum(0.);     fpP2Mult[i]->SetMaximum(80.);
    setTitles(fpP2Mult[i], "l_{z}^{(2)} #it{ [cm]}", "PV trk mult", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP2fls3d[i] = (TProfile*)fHistFile->Get(Form("p2fls3d_%s_chan%d", mode, i))))
      fpP2fls3d[i] = new TProfile(Form("p2fls3d_%s_chan%d", mode, i), "fls3d", 40, -2., 2.0, "");
    setTitles(fpP2fls3d[i], "l_{z}^{(2)} #it{ [cm]}", "mean fls3d", 0.05, 1.1, 1.8, 0.04);
    fpP2fls3d[i]->SetMinimum(0.);     fpP2fls3d[i]->SetMaximum(50.);
    if (!(fpP2fl3d[i] = (TProfile*)fHistFile->Get(Form("p2fl3d_%s_chan%d", mode, i))))
      fpP2fl3d[i] = new TProfile(Form("p2fl3d_%s_chan%d", mode, i), "fl3d", 40, -2., 2.0, "");
    setTitles(fpP2fl3d[i], "l_{z}^{(2)} #it{ [cm]}", "mean fl3d", 0.05, 1.1, 1.8, 0.04);
    fpP2fl3d[i]->SetMinimum(0.);     fpP2fl3d[i]->SetMaximum(0.5);
    if (!(fpP2dfl3d[i] = (TProfile*)fHistFile->Get(Form("p2dfl3d_%s_chan%d", mode, i))))
      fpP2dfl3d[i] = new TProfile(Form("p2dfl3d_%s_chan%d", mode, i), "dfl3d", 40, -2., 2.0, "");
    fpP2dfl3d[i]->SetMinimum(-0.03); fpP2dfl3d[i]->SetMaximum(0.03);
    setTitles(fpP2dfl3d[i], "l_{z}^{(2)} #it{ [cm]}", "mean #Delta fl3d", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP2Iso[i] = (TProfile*)fHistFile->Get(Form("p2Iso_%s_chan%d", mode, i))))
      fpP2Iso[i] = new TProfile(Form("p2Iso_%s_chan%d", mode, i), "Iso", 40, -2., 2.0, "");
    fpP2Iso[i]->SetMinimum(0.3); fpP2Iso[i]->SetMaximum(1.01);
    setTitles(fpP2Iso[i], "l_{z}^{(2)} #it{ [cm]}", "mean Isolation", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP2Ntrk[i] = (TProfile*)fHistFile->Get(Form("p2Ntrk_%s_chan%d", mode, i))))
      fpP2Ntrk[i] = new TProfile(Form("p2Ntrk_%s_chan%d", mode, i), "Ntrk", 40, -2., 2.0, "");
    fpP2Ntrk[i]->SetMinimum(0.); fpP2Ntrk[i]->SetMaximum(20.);
    setTitles(fpP2Ntrk[i], "l_{z}^{(2)} #it{ [cm]}", "mean closetrk", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP2PvIps[i] = (TProfile*)fHistFile->Get(Form("p2PvIps_%s_chan%d", mode, i))))
      fpP2PvIps[i] = new TProfile(Form("p2PvIps_%s_chan%d", mode, i), "PvIps", 40, -2., 2.0, "");
    fpP2PvIps[i]->SetMinimum(0.); fpP2PvIps[i]->SetMaximum(5.);
    setTitles(fpP2PvIps[i], "l_{z}^{(2)} #it{ [cm]}", "mean PvIps", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP2Bdt[i] = (TProfile*)fHistFile->Get(Form("p2Bdt_%s_chan%d", mode, i))))
      fpP2Bdt[i] = new TProfile(Form("p2Bdt_%s_chan%d", mode, i), "Bdt", 40, -2., 2.0, "");
    fpP2Bdt[i]->SetMinimum(-1.); fpP2Bdt[i]->SetMaximum(1.);
    setTitles(fpP2Bdt[i], "l_{z}^{(2)} #it{ [cm]}", "mean Bdt", 0.05, 1.1, 1.8, 0.04);

    if (!(fpP2l1[i] = (TProfile*)fHistFile->Get(Form("p2l1_%s_chan%d", mode, i))))
      fpP2l1[i] = new TProfile(Form("p2l1_%s_chan%d", mode, i), "lz1", 40, -2., 2.0, "");
    fpP2l1[i]->SetMinimum(-0.1); fpP2l1[i]->SetMaximum(0.1);
    setTitles(fpP2l1[i], "l_{z}^{(2)} #it{ [cm]}", "mean l_{z}^{(1)} [cm]", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP2dzmin[i] = (TProfile*)fHistFile->Get(Form("p2dzmin_%s_chan%d", mode, i))))
      fpP2dzmin[i] = new TProfile(Form("pP2dzmin_%s_chan%d", mode, i), "lz2", 40, -2., 2.0, "");
    fpP2dzmin[i]->SetMinimum(-1.); fpP2dzmin[i]->SetMaximum(1.);
    setTitles(fpP2dzmin[i], "l_{z}^{(2)} #it{ [cm]}", "mean minimum #Delta z [cm]", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP2npv[i] = (TProfile*)fHistFile->Get(Form("p2npv_%s_chan%d", mode, i))))
      fpP2npv[i] = new TProfile(Form("p2npv_%s_chan%d", mode, i), "npv", 40, -2.0, 2.0, "");
    fpP2npv[i]->SetMinimum(0.); fpP2npv[i]->SetMaximum(40.);
    setTitles(fpP2npv[i], "minimum #it{|}#Delta z#it{| [cm]}", "mean N_{PV}", 0.05, 1.1, 1.8, 0.04);

    if (!(fpP2tau[i] = (TProfile*)fHistFile->Get(Form("p2tau_%s_chan%d", mode, i))))
      fpP2tau[i] = new TProfile(Form("p2tau_%s_chan%d", mode, i), "tau", 40, -2., 2.0, "");
    fpP2tau[i]->SetMinimum(0.e-12); fpP2tau[i]->SetMaximum(+2.e-12);
    setTitles(fpP2tau[i], "l_{z}^{(2)} #it{ [cm]}", "mean #tau_{eff} ", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP2dtau[i] = (TProfile*)fHistFile->Get(Form("p2dtau_%s_chan%d", mode, i))))
      fpP2dtau[i] = new TProfile(Form("p2dtau_%s_chan%d", mode, i), "dtau", 40, -2., 2.0, "");
    fpP2dtau[i]->SetMinimum(-0.2e-12); fpP2dtau[i]->SetMaximum(+0.2e-12);
    setTitles(fpP2dtau[i], "l_{z}^{(2)} #it{ [cm]}", "mean #Delta#tau_{eff}", 0.05, 1.1, 1.8, 0.04);


    // -- vs NPV
    if (!(fpP3Mult[i] = (TProfile*)fHistFile->Get(Form("p3mult_%s_chan%d", mode, i))))
      fpP3Mult[i] = new TProfile(Form("p3mult_%s_chan%d", mode, i), "mult", 40, 0., 40.0, "");
    fpP3Mult[i]->SetMinimum(0.);     fpP3Mult[i]->SetMaximum(80.);
    setTitles(fpP3Mult[i], "N_{PV}", "PV trk mult", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP3fls3d[i] = (TProfile*)fHistFile->Get(Form("p3fls3d_%s_chan%d", mode, i))))
      fpP3fls3d[i] = new TProfile(Form("p3fls3d_%s_chan%d", mode, i), "fls3d", 40, 0., 40.0, "");
    setTitles(fpP3fls3d[i], "N_{PV}", "mean fls3d", 0.05, 1.1, 1.8, 0.04);
    fpP3fls3d[i]->SetMinimum(0.);     fpP3fls3d[i]->SetMaximum(50.);
    if (!(fpP3fl3d[i] = (TProfile*)fHistFile->Get(Form("p3fl3d_%s_chan%d", mode, i))))
      fpP3fl3d[i] = new TProfile(Form("p3fl3d_%s_chan%d", mode, i), "fl3d", 40, 0., 40.0, "");
    setTitles(fpP3fl3d[i], "N_{PV}", "mean fl3d", 0.05, 1.1, 1.8, 0.04);
    fpP3fl3d[i]->SetMinimum(0.);     fpP3fl3d[i]->SetMaximum(0.5);
    if (!(fpP3dfl3d[i] = (TProfile*)fHistFile->Get(Form("p3dfl3d_%s_chan%d", mode, i))))
      fpP3dfl3d[i] = new TProfile(Form("p3dfl3d_%s_chan%d", mode, i), "dfl3d", 40, 0., 40.0, "");
    fpP3dfl3d[i]->SetMinimum(-0.03); fpP3dfl3d[i]->SetMaximum(0.03);
    setTitles(fpP3dfl3d[i], "N_{PV}", "mean #Delta fl3d", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP3Iso[i] = (TProfile*)fHistFile->Get(Form("p3Iso_%s_chan%d", mode, i))))
      fpP3Iso[i] = new TProfile(Form("p3Iso_%s_chan%d", mode, i), "Iso", 40, 0., 40.0, "");
    fpP3Iso[i]->SetMinimum(0.3); fpP3Iso[i]->SetMaximum(1.01);
    setTitles(fpP3Iso[i], "N_{PV}", "mean Isolation", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP3Ntrk[i] = (TProfile*)fHistFile->Get(Form("p3Ntrk_%s_chan%d", mode, i))))
      fpP3Ntrk[i] = new TProfile(Form("p3Ntrk_%s_chan%d", mode, i), "Ntrk", 40, 0., 40.0, "");
    fpP3Ntrk[i]->SetMinimum(0.); fpP3Ntrk[i]->SetMaximum(20.);
    setTitles(fpP3Ntrk[i], "N_{PV}", "mean closetrk", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP3PvIps[i] = (TProfile*)fHistFile->Get(Form("p3PvIps_%s_chan%d", mode, i))))
      fpP3PvIps[i] = new TProfile(Form("p3PvIps_%s_chan%d", mode, i), "PvIps", 40, 0., 40.0, "");
    fpP3PvIps[i]->SetMinimum(0.); fpP3PvIps[i]->SetMaximum(5.);
    setTitles(fpP3PvIps[i], "N_{PV}", "mean PvIps", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP3Bdt[i] = (TProfile*)fHistFile->Get(Form("p3Bdt_%s_chan%d", mode, i))))
      fpP3Bdt[i] = new TProfile(Form("p3Bdt_%s_chan%d", mode, i), "Bdt", 40, 0., 40.0, "");
    fpP3Bdt[i]->SetMinimum(-1.); fpP3Bdt[i]->SetMaximum(1.);
    setTitles(fpP3Bdt[i], "N_{PV}", "mean Bdt", 0.05, 1.1, 1.8, 0.04);

    if (!(fpP3l1[i] = (TProfile*)fHistFile->Get(Form("p3l1_%s_chan%d", mode, i))))
      fpP3l1[i] = new TProfile(Form("p3l1_%s_chan%d", mode, i), "lz1", 40, 0., 40.0, "");
    fpP3l1[i]->SetMinimum(-0.1); fpP3l1[i]->SetMaximum(0.1);
    setTitles(fpP3l1[i], "N_{PV}", "mean l_{z}^{(1)} [cm]", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP3l2[i] = (TProfile*)fHistFile->Get(Form("p3l2_%s_chan%d", mode, i))))
      fpP3l2[i] = new TProfile(Form("p3l2_%s_chan%d", mode, i), "npv", 40, 0.0, 40.0, "S");
    fpP3l2[i]->SetMinimum(0.); fpP3l2[i]->SetMaximum(1.);
    setTitles(fpP3l2[i], "N_{PV}", "mean l_{z}^{(2)} [cm]", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP3dzmin[i] = (TProfile*)fHistFile->Get(Form("p3dzmin_%s_chan%d", mode, i))))
      fpP3dzmin[i] = new TProfile(Form("p3dzmin_%s_chan%d", mode, i), "lz2", 40, 0., 40.0, "S");
    fpP3dzmin[i]->SetMinimum(0.); fpP3dzmin[i]->SetMaximum(1.);
    setTitles(fpP3dzmin[i], "N_{PV}", "mean minimum #Delta z [cm]", 0.05, 1.1, 1.8, 0.04);

    if (!(fpP3tau[i] = (TProfile*)fHistFile->Get(Form("p3tau_%s_chan%d", mode, i))))
      fpP3tau[i] = new TProfile(Form("p3tau_%s_chan%d", mode, i), "tau", 40, 0., 40.0, "");
    fpP3tau[i]->SetMinimum(0.e-12); fpP3tau[i]->SetMaximum(+2.e-12);
    setTitles(fpP3tau[i], "N_{PV}", "mean #tau_{eff} ", 0.05, 1.1, 1.8, 0.04);
    if (!(fpP3dtau[i] = (TProfile*)fHistFile->Get(Form("p3dtau_%s_chan%d", mode, i))))
      fpP3dtau[i] = new TProfile(Form("p3dtau_%s_chan%d", mode, i), "dtau", 40, 0., 40.0, "");
    fpP3dtau[i]->SetMinimum(-0.2e-12); fpP3dtau[i]->SetMaximum(+0.2e-12);
    setTitles(fpP3dtau[i], "N_{PV}", "mean #Delta#tau_{eff}", 0.05, 1.1, 1.8, 0.04);

  }

  // ----------------------------------------------------------------------
  // -- Loop over tree if necessary
  // ----------------------------------------------------------------------
  if (fillHist) {
    fSample = dsname;
    setup(fSample);

    TTree *t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    fCds = fDS[fSample];
    loopOverTree(t, 3);
    fHistFile->Write();
  }

  // ----------------------------------------------------------------------
  // -- Plot things
  // ----------------------------------------------------------------------
  gStyle->SetOptFit(0);

  TF1 *f1 = fIF->pol1(-999., 999.);
  f1->SetLineColor(kBlack);
  for (int i = 0; i < 2; ++i) {
    // -- PVN
    shrinkPad(0.13, 0.20);
    gPad->SetGrid(1,1);
    fpHpvn[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-pvn.pdf", i, mode));

    // -- minimum z separation
    shrinkPad(0.13, 0.20);
    gPad->SetGrid(1,1);
    fpHdzmin[i]->Draw();
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-dzmin.pdf", i, mode));

    // -- overlay pvlip and pv2lip (full range)
    gPad->SetLogy(1);
    setFilledHist(fpHlz1[i], kBlue, kBlue, 3354);
    setFilledHist(fpHlz2[i], kRed, kRed, 3365);
    setTitles(fpHlz1[i], "l_{z} #it{[cm]}", "a.u.", 0.05, 1.1, 1.8);
    fpHlz1[i]->SetMaximum(35.*fpHlz1[i]->GetMaximum());
    fpHlz1[i]->Draw();
    fpHlz2[i]->Draw("same");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));

    newLegend(0.41, 0.7, 0.71, 0.90);
    legg->SetHeader("longitudinal IP to PV");
    legg->SetTextSize(0.05);
    legg->AddEntry(fpHlz1[i], "best PV", "f");
    legg->AddEntry(fpHlz2[i], "second-best PV", "f");
    legg->Draw();

    savePad(Form("pustudyChan%d-%s-overlay-lz1-lz2-full.pdf", i, mode));

    // -- overlay pvlip and pv2lip (reduced range)
    gPad->SetLogy(1);
    fpHlz1[i]->SetAxisRange(-0.15, 0.15, "X");
    fpHlz1[i]->Draw();
    fpHlz2[i]->Draw("same");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));

    newLegend(0.41, 0.7, 0.71, 0.90);
    legg->SetHeader("longitudinal IP to PV");
    legg->SetTextSize(0.05);
    legg->AddEntry(fpHlz1[i], "best PV", "f");
    legg->AddEntry(fpHlz2[i], "second-best PV", "f");
    legg->Draw();

    savePad(Form("pustudyChan%d-%s-overlay-lz1-lz2-reduced.pdf", i, mode));


    // -- overlay pvlips and pv2lips
    gPad->SetLogy(1);
    setFilledHist(fpHlzs1[i], kBlue, kBlue, 3354);
    setFilledHist(fpHlzs2[i], kRed, kRed, 3365);
    setTitles(fpHlzs1[i], "#it{|} l_{z}/#sigma(l_{z}) #it{|}", "a.u.", 0.05, 1.1, 1.8);
    fpHlzs1[i]->Draw();
    fpHlzs2[i]->Draw("same");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));

    newLegend(0.41, 0.7, 0.71, 0.90);
    legg->SetHeader("long. IP sig. to PV");
    legg->SetTextSize(0.05);
    legg->AddEntry(fpHlzs1[i], "best PV", "f");
    legg->AddEntry(fpHlzs2[i], "second-best PV", "f");
    legg->Draw();

    savePad(Form("pustudyChan%d-%s-overlay-lzs1-lzs2.pdf", i, mode));



    // -- PV multiplicity
    shrinkPad(0.13, 0.20);
    gPad->SetLogy(0);
    overlay(fpHmultFar[i], "close", fpHmultClose05[i], "far", 0, "", UNITY, false, false);
    newLegend(0.41, 0.7, 0.71, 0.90);
    legg->SetHeader("other PV");
    legg->SetTextSize(0.05);
    legg->AddEntry(fpHmultFar[i], "#it{|}d_{z}^{min}#it{|} > 1cm", "f");
    legg->AddEntry(fpHmultClose05[i], "#it{|}d_{z}^{min}#it{|} < 0.05cm", "f");
    legg->Draw();

    tl->SetTextSize(0.03);
    tl->SetTextColor(kGreen+2);
    tl->DrawLatexNDC(0.60, 0.62, Form("#mu/RMS = %3.1f/%3.1f", fpHmultClose05[i]->GetMean(), fpHmultClose05[i]->GetRMS()));
    tl->SetTextColor(kBlue);
    tl->DrawLatexNDC(0.65, 0.5, Form("#mu/RMS = %3.1f/%3.1f", fpHmultFar[i]->GetMean(), fpHmultFar[i]->GetRMS()));
    tl->SetTextColor(kBlack);

    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-farCloseMultiplicity.pdf", i, mode));


    // -- PV multiplicity lz2
    shrinkPad(0.13, 0.20);
    overlay(fpHmultFar2[i], "close", fpHmultClose2[i], "far", 0, "", UNITY, false, false);
    newLegend(0.41, 0.7, 0.71, 0.90);
    legg->SetHeader("other PV");
    legg->SetTextSize(0.05);
    legg->AddEntry(fpHmultFar2[i], "#it{|} l_{z}^{2} #it{|} > 1cm", "f");
    legg->AddEntry(fpHmultClose2[i], "#it{|} l_{z}^{2} #it{|} < 0.05cm", "f");
    legg->Draw();

    tl->SetTextSize(0.03);
    tl->SetTextColor(kGreen+2);
    tl->DrawLatexNDC(0.60, 0.62, Form("#mu/RMS = %3.1f/%3.1f", fpHmultClose2[i]->GetMean(), fpHmultClose2[i]->GetRMS()));
    tl->SetTextColor(kBlue);
    tl->DrawLatexNDC(0.65, 0.5, Form("#mu/RMS = %3.1f/%3.1f", fpHmultFar2[i]->GetMean(), fpHmultFar2[i]->GetRMS()));
    tl->SetTextColor(kBlack);

    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-farCloseMultiplicity2.pdf", i, mode));


    // -- profiles vs dzmin
    fpP0Mult[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP0Mult[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof0-mult.pdf", i, mode));

    fpP0fls3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP0fls3d[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof0-fls3d.pdf", i, mode));
    fpP0fl3d[i]->Fit(f1, "q");
    if (dsname2 != "") o2Profile(fpP0fl3d[i], dsname, dsname2, i);
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-fl3d.pdf", i, mode));
    fpP0dfl3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-dfl3d.pdf", i, mode));

    fpP0l1[i]->Fit(f1, "q");
    if (dsname2 != "") o2Profile(fpP0l1[i], dsname, dsname2, i);
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-l1z.pdf", i, mode));
    fpP0l2[i]->Fit(f1, "q");
    if (dsname2 != "") o2Profile(fpP0l2[i], dsname, dsname2, i);
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-l2z.pdf", i, mode));
    fpP0npv[i]->Fit(f1, "q");
    if (dsname2 != "") o2Profile(fpP0npv[i], dsname, dsname2, i);
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-npv.pdf", i, mode));

    fpP0tau[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    if (dsname2 != "") o2Profile(fpP0tau[i], dsname, dsname2, i);
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %9.8f #pm %9.8f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-tau.pdf", i, mode));
    fpP0dtau[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-dtau.pdf", i, mode));

    fpP0Iso[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP0Iso[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof0-iso.pdf", i, mode));
    fpP0Ntrk[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP0Ntrk[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof0-ntrk.pdf", i, mode));
    fpP0PvIps[i]->Fit(f1, "q");
    if (dsname2 != "") o2Profile(fpP0PvIps[i], dsname, dsname2, i);
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-pvips.pdf", i, mode));
    fpP0Bdt[i]->Fit(f1, "q");
    if (dsname2 != "") o2Profile(fpP0Bdt[i], dsname, dsname2, i);
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof0-bdt.pdf", i, mode));


    // -- profiles vs lz1
    fpP1Mult[i]->Fit(f1, "q");
    savePad(Form("pustudyChan%d-%s-prof1-mult.pdf", i, mode));
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP1Mult[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof1-mult.pdf", i, mode));
    fpP1fls3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP1fls3d[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof1-fls3d.pdf", i, mode));
    fpP1fl3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP1fl3d[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof1-fl3d.pdf", i, mode));
    fpP1dfl3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof1-dfl3d.pdf", i, mode));

    fpP1dzmin[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    if (dsname2 != "") o2Profile(fpP1dzmin[i], dsname, dsname2, i);
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof1-dzmin.pdf", i, mode));
    fpP1l2[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP1l2[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof1-l2z.pdf", i, mode));
    fpP1npv[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP1npv[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof1-npv.pdf", i, mode));

    fpP1tau[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %9.8f #pm %9.8f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP1tau[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof1-tau.pdf", i, mode));
    fpP1dtau[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof1-dtau.pdf", i, mode));

    fpP1Iso[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    if (dsname2 != "") o2Profile(fpP1Iso[i], dsname, dsname2, i);
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof1-iso.pdf", i, mode));
    fpP1Ntrk[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    if (dsname2 != "") o2Profile(fpP1Ntrk[i], dsname, dsname2, i);
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof1-ntrk.pdf", i, mode));
    fpP1PvIps[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP1PvIps[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof1-pvips.pdf", i, mode));
    fpP1Bdt[i]->Fit(f1, "q");
    if (dsname2 != "") o2Profile(fpP1Bdt[i], dsname, dsname2, i);
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    savePad(Form("pustudyChan%d-%s-prof1-bdt.pdf", i, mode));

    // -- profiles vs lz2
    fpP2Mult[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    if (dsname2 != "") o2Profile(fpP2Mult[i], dsname, dsname2, i);
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof2-mult.pdf", i, mode));
    fpP2fls3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2fls3d[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-fls3d.pdf", i, mode));
    fpP2fl3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2fl3d[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-fl3d.pdf", i, mode));
    fpP2dfl3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof2-dfl3d.pdf", i, mode));

    fpP2dzmin[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2dzmin[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-dzmin.pdf", i, mode));
    fpP2l1[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2l1[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-l1z.pdf", i, mode));
    fpP2npv[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2npv[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-npv.pdf", i, mode));

    fpP2tau[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %9.8f #pm %9.8f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2tau[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-tau.pdf", i, mode));
    fpP2dtau[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof2-dtau.pdf", i, mode));

    fpP2Iso[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2Iso[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-iso.pdf", i, mode));
    fpP2Ntrk[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2Ntrk[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-ntrk.pdf", i, mode));
    fpP2PvIps[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2PvIps[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-pvips.pdf", i, mode));
    fpP2Bdt[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP2Bdt[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof2-bdt.pdf", i, mode));


    // -- profiles vs pvn
    fpP3Mult[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3Mult[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-mult.pdf", i, mode));
    fpP3fls3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3fls3d[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-fls3d.pdf", i, mode));
    fpP3fl3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3fl3d[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-fl3d.pdf", i, mode));
    fpP3dfl3d[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof3-dfl3d.pdf", i, mode));

    fpP3dzmin[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3dzmin[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-dzmin.pdf", i, mode));
    fpP3l1[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3l1[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-l1z.pdf", i, mode));
    fpP3l2[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3l2[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-l2z.pdf", i, mode));

    fpP3tau[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %9.8f #pm %9.8f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3tau[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-tau.pdf", i, mode));
    fpP3dtau[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    savePad(Form("pustudyChan%d-%s-prof3-dtau.pdf", i, mode));

    fpP3Iso[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3Iso[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-iso.pdf", i, mode));
    fpP3Ntrk[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3Ntrk[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-ntrk.pdf", i, mode));
    fpP3PvIps[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3PvIps[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-pvips.pdf", i, mode));
    fpP3Bdt[i]->Fit(f1, "q");
    tl->SetTextSize(0.04); tl->DrawLatexNDC(0.25, 0.92, Form("Chan %d", i));
    tl->SetTextSize(0.03); tl->DrawLatexNDC(0.5, 0.92, Form("p1 = %5.4f #pm %5.4f", f1->GetParameter(1), f1->GetParError(1)));
    if (dsname2 != "") o2Profile(fpP3Bdt[i], dsname, dsname2, i);
    savePad(Form("pustudyChan%d-%s-prof3-bdt.pdf", i, mode));
  }

  fHistFile->Close();
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
    T->Draw(Form("p3z-gz >> dz2_ch%d", i), Form("chan == %d && m1pt>4 && m2pt>4 %s", i, selection.c_str()));
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
  cout << "yieldStability> dsname: " << dsname << " for setup: " << trg << endl;
  double MINLUMI(2.);
  double mBp(5.28), sBp(0.015), stepBp(5.15);
  double xmin(5.0), xmax(5.9), ymax(0.), expoLo(5.16), expoHi(5.85);

  int nchan = fNchan;
  nchan = 2;

  fSample = trg;
  fMode = BMM;
  setup(dsname);
  if (string::npos != fSetup.find("bspsiphi")) {
    mBp    = 5.369;
    sBp    = 0.015;
    stepBp = 5.15;
  }


  // -- check whether there are any histograms pre-produced already
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  bool ok = fHistFile->cd(fTreeDir.c_str());
  cout << "OK = " << ok << endl;
  TH2D *h2(0), *hBlock(0);

  TDirectory *hDir(0);
  if (!ok) {
    hDir = fHistFile->mkdir(fTreeDir.c_str());
    fHistFile->cd(fTreeDir.c_str());
    hDir = gDirectory;
  }

  if (ok) {
    ok = false;
    TIter next(gDirectory->GetListOfKeys());
    TKey *key(0);
    while ((key = (TKey*)next())) {
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
      if (TString(key->GetName()).Contains(Form("%s_", trg.c_str()))) {
	ok = true;
      }
    }
  }
  if (!ok) {
    TTree *t = getTree(fSetup, fTreeDir);
    if (0 == t) {
      cout << "tree for sample = " << fSetup << " not found" << endl;
      return;
    }
    setupTree(t, fSetup);
    fCds = fDS[fSetup];
    //    loopOverTree(t, 1, 8e6, 17500000); // FIXME
    loopOverTree(t, 1);

    cout << "writing output histograms: " << fYieldHists.size() << endl;
    for (map<string, TH2D*>::iterator it = fYieldHists.begin(); it != fYieldHists.end(); ++it) {
      cout << "run " << it->first << endl;
      it->second->Draw("colz");
      it->second->SetDirectory(hDir);
      it->second->Write();
    }
    fYieldHists.clear();
  } else {
  }

  cout << "histograms exist already, looping over them" << endl;
  TIter next(gDirectory->GetListOfKeys());
  TKey *key(0);
  int run(-1), runMin(9999999), runMax(0), firstLumiRun(99), lastLumiRun(99);
  vector<int> vruns;
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
    if (TString(key->GetName()).Contains(Form("%s_", trg.c_str()))) {
      string hname = key->GetName();
      if (0 == hBlock) {
	hBlock = (TH2D*)((TH2D*)gDirectory->Get(hname.c_str()))->Clone("hBlock");
	hBlock->Reset();
      }
      replaceAll(hname, Form("%s_", trg.c_str()), "");
      run = atoi(hname.c_str());
      if (run > runMax) runMax = run;
      if (run < runMin) runMin = run;
      if (find(vruns.begin(), vruns.end(), run) == vruns.end()) {
	vruns.push_back(run);
      }
    }
  }

  if (vruns.size() > 0) {
    cout << "analyzing runs " << runMin << " .. " <<  runMax << endl;

    // -- create run blocks based on integrated lumi
    Lumi *lumi(0);
    if (2016 == fYear) {
      if (string::npos != dsname.find("psi")) {
	lumi = new Lumi("../common/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
	//this is useless, includes DiMuonX_Y seeds: lumi = new Lumi("../common/json/2016-HLT_DoubleMu4_3_Jpsi_Displaced.lumi");
      } else {
	lumi = new Lumi("../common/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
	//this is useless, includes DiMuonX_Y seeds: lumi = new Lumi("../common/json/2016-HLT_DoubleMu4_3_Bs.lumi");
      }
    } else if (2012 == fYear) {
      lumi = new Lumi("../common/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.lumi");
    } else if (2011 == fYear) {
      lumi = new Lumi("../common/json/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_MuonPhys_v2.lumi");
    }
    firstLumiRun = lumi->firstRun();
    lastLumiRun  = lumi->lastRun();
    cout << "lumiRuns = " << firstLumiRun << " .. " << lastLumiRun << endl;
    double intLumi(0.);
    map<pair<int, double>, vector<int> > runBlocks;
    vector<int> segment;
    for (unsigned int irun = 0; irun < vruns.size(); ++irun) {
      intLumi += lumi->lumi(vruns[irun]);
      segment.push_back(vruns[irun]);
      if (intLumi > MINLUMI) {
	cout << "Adding " << segment[0] << ": " << intLumi << endl;
	runBlocks.insert(make_pair(make_pair(segment[0], intLumi), segment));
	intLumi = 0.;
	segment.clear();
      }
    }
    // -- add the last block as well
    if (segment[0] > 0) {
      cout << "Final adding " << segment[0] << ": " << intLumi << endl;
      runBlocks.insert(make_pair(make_pair(segment[0], intLumi), segment));
    }

    // -- the result histograms
    vector<TH1D *> vRunHLT, vRunHLTLumi;
    for (unsigned int ichan = 0; ichan < nchan; ++ichan) {
      vRunHLT.push_back(new TH1D(Form("hRun%s_%s_chan%d", trg.c_str(), dsname.c_str(), ichan),
				 Form("hRun%s_%s_chan%d", trg.c_str(), dsname.c_str(), ichan),
				 lastLumiRun-firstLumiRun+1, firstLumiRun, lastLumiRun));
      vRunHLT[ichan]->Sumw2();

      vRunHLTLumi.push_back(new TH1D(Form("hRunLumi%s_%s_chan%d", trg.c_str(), dsname.c_str(), ichan),
				 Form("hRunLumi%s_%s_chan%d", trg.c_str(), dsname.c_str(), ichan),
				 lastLumiRun-firstLumiRun+1, firstLumiRun, lastLumiRun));
      vRunHLTLumi[ichan]->Sumw2();
    }

    // -- get the histograms
    //    fHistFile = TFile::Open(fHistFileName.c_str());
    string hname("");
    TDirectory *pDir = fHistFile->GetDirectory(fTreeDir.c_str());
    for (int ichan = 0; ichan < nchan; ++ichan) {
      double lumi(0.), totalLumi(0.);
      cout << "--> chan " << ichan << endl;
      for (map<pair<int, double>, vector<int> >::iterator it = runBlocks.begin(); it != runBlocks.end(); ++it) {
	int iblock = it->first.first;
	lumi = it->first.second;
	totalLumi += lumi;
	cout << endl;
	cout << Form("new block: %d (%d) Lumi = %4.1f/%4.1f ", it->first.first, iblock, lumi, totalLumi) << endl;
	hBlock->Reset();
	hBlock->SetName(Form("hBlock_%s_%d_chan%d", dsname.c_str(), iblock, ichan));
	for (unsigned int i = 0; i < it->second.size(); ++i) {
	  hname = Form("%s/%s_%d_chan%d", fTreeDir.c_str(), trg.c_str(), it->second[i], ichan);
	  h2 = (TH2D*)(fHistFile->Get(hname.c_str()));
	  cout << it->second[i] << " (" << hname << ": " << h2 << ") ";
	  if (0 == h2) continue;
	  cout << "adding " << hname << " h2 = " << h2 <<  " with lumi: " << lumi << endl;
	  hBlock->Add(h2);
	}
	double result(0.), resultE(0.);
	if (string::npos != dsname.find("bupsik") || string::npos != dsname.find("bspsiphi")) {
	  fitPsYield fpy(hBlock, 0);
	  if (string::npos != dsname.find("bupsik")) {
	    fpy.fitBu2JpsiKp(0, fDirectory + "/ys-" + trg + "-");
	  } else if (string::npos != dsname.find("bspsiphi")) {
	    fpy.fitBs2JpsiPhi(0, fDirectory + "/");
	  }
	  if (1) {
	    result  = fpy.getSignalYield();
	    resultE = fpy.getSignalError();
	  }
	  if (0) {
	    result = hBlock->Integral(1, hBlock->GetNbinsX(), 2, 2);
	    resultE = TMath::Sqrt(result);
	  }
	  if (0) {
	    result  = fpy.getSignalW8Yield();
	    resultE = fpy.getSignalW8Error();
	  }
	  if (0) {
	    result  = fpy.getSignalUnW8Yield();
	    resultE = fpy.getSignalUnW8Error();
	  }
	} else {
	  result = hBlock->Integral(1, hBlock->GetNbinsX(), 2, 2);
	  resultE = TMath::Sqrt(result);
	  c0->Clear();
	  hBlock->Draw("colz");
	  savePad(Form("ys-hBlock_%s_%d-chan%d.pdf", dsname.c_str(), iblock, ichan));
	}
	cout << "result = " << result << " +/- " << resultE
	     << " lumi-normalized = " << result/lumi << " +/- " << resultE/lumi
	     << " (lumi = " << lumi << ")"
	     << " filling into bin " << vRunHLT[ichan]->FindBin(static_cast<double>(iblock)) << endl;
	vRunHLTLumi[ichan]->SetBinContent(vRunHLTLumi[ichan]->FindBin(static_cast<double>(iblock)), result/lumi);
	vRunHLTLumi[ichan]->SetBinError(vRunHLTLumi[ichan]->FindBin(static_cast<double>(iblock)), resultE/lumi);

	vRunHLT[ichan]->SetBinContent(vRunHLT[ichan]->FindBin(static_cast<double>(iblock)), result);
	vRunHLT[ichan]->SetBinError(vRunHLT[ichan]->FindBin(static_cast<double>(iblock)), resultE);
      }
    }

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gPad->SetGridx();
    gPad->SetGridy();
    for (unsigned ichan = 0; ichan < nchan; ++ichan) {
      setTitles(vRunHLT[ichan], "run", Form("N(%s)", fDS[dsname]->fName.c_str()), 0.05, 1.1, 1.9);
      setTitles(vRunHLTLumi[ichan], "run", Form("N(%s)", fDS[dsname]->fName.c_str()), 0.05, 1.1, 1.9);
      vRunHLTLumi[ichan]->Draw();
      if (2016 == fYear) {
	double ymax(vRunHLTLumi[ichan]->GetMaximum());
	plot2016Eras(ymax);
      }
      savePad(Form("ys-yieldPerLumi-%s%d-%s-chan%d.pdf", trg.c_str(), fYear, dsname.c_str(), ichan));

      vRunHLT[ichan]->Draw();
      if (2016 == fYear) {
	double ymax(vRunHLTLumi[ichan]->GetMaximum());
	plot2016Eras(ymax);
      }
      savePad(Form("ys-yield-%s%d-%s-chan%d.pdf", trg.c_str(), fYear, dsname.c_str(), ichan));

      if (1) {
	vRunHLTLumi[ichan]->SetDirectory(gDirectory);
	vRunHLTLumi[ichan]->Write();

	vRunHLT[ichan]->SetDirectory(gDirectory);
	vRunHLT[ichan]->Write();
      }
    }

    delete lumi;
  }

  fHistFile->Close();

}



// ----------------------------------------------------------------------
void plotStuff::massResolution(std::string file1, std::string file2) {
  fHmass0.clear();
  fHmass1.clear();
  for (int i = 0; i < 25; ++i) {
    fHmass0.push_back(new TH1D(Form("Hmass0_%d", i), Form(" %3.1f < |#eta_{f}| < %3.1f", i*0.1, (i+1)*0.1), 80, 5.0, 5.8));
    fHmass1.push_back(new TH1D(Form("Hmass1_%d", i), Form(" %3.1f < |#eta_{f}| < %3.1f", i*0.1, (i+1)*0.1), 80, 5.0, 5.8));
  }

  // -- sample 1
  fSample = file1;
  setup(fSample);

  TTree *t0 = getTree(fSample, fTreeDir);
  setupTree(t0, fSample);
  fCds = fDS[fSample];
  fDBX = 0;
  loopOverTree(t0, 5);


  // -- sample 2
  fSample = file2;
  setup(fSample);

  TTree *t1 = getTree(fSample, fTreeDir);
  setupTree(t1, fSample);
  fCds = fDS[fSample];
  fDBX = 1;
  loopOverTree(t1, 5);

  gStyle->SetOptFit(0);

  shrinkPad(0.15, 0.19);
  double eps(0.01);
  TH1D *h1 = new TH1D("mass0", "Run 1", 25, 0., 2.5); h1->Sumw2();
  TH1D *h2 = new TH1D("mass1", "Run 2", 25, 0.+eps, 2.5+eps); h2->Sumw2();
  TH1D *s1 = new TH1D("rms0", "Run 1", 25, 0., 2.5); h1->Sumw2();
  TH1D *s2 = new TH1D("rms1", "Run 2", 25, 0.+eps, 2.5+eps); h2->Sumw2();
  TH1D *p1 = new TH1D("peak0", "Run 1", 25, 0., 2.5); p1->Sumw2();
  TH1D *p2 = new TH1D("peak1", "Run 2", 25, 0.+eps, 2.5+eps); p2->Sumw2();
  TH1D *w1 = new TH1D("sigma0", "Run 1", 25, 0., 2.5); w1->Sumw2();
  TH1D *w2 = new TH1D("sigma1", "Run 2", 25, 0.+eps, 2.5+eps); w2->Sumw2();
  for (int i = 0; i < 25; ++i) {
    cout << "mass = " << fHmass0[i]->GetMean() << " RMS = " << fHmass0[i]->GetRMS() << endl;
    h1->SetBinContent(i+1, fHmass0[i]->GetMean());
    h1->SetBinError(i+1, fHmass0[i]->GetMeanError());
    s1->SetBinContent(i+1, fHmass0[i]->GetRMS());
    s1->SetBinError(i+1, fHmass0[i]->GetRMSError());

    cout << "mass = " << fHmass1[i]->GetMean() << " RMS = " << fHmass1[i]->GetRMS() << endl;
    h2->SetBinContent(i+1, fHmass1[i]->GetMean());
    h2->SetBinError(i+1, fHmass1[i]->GetMeanError());
    s2->SetBinContent(i+1, fHmass1[i]->GetRMS());
    s2->SetBinError(i+1, fHmass1[i]->GetRMSError());

    setFilledHist(fHmass0[i], kBlue, kBlue, 3365);
    setFilledHist(fHmass1[i], kRed, kRed, 3354);

    // -- do the fitting before the scaling
    double peak0V(0.),  peak0E(0.),  peak1V(0.),  peak1E(0.);
    double sigma0V(0.), sigma0E(0.), sigma1V(0.), sigma1E(0.);
    if (fHmass1[i]->GetSumOfWeights() > 0.) {
      fHmass1[i]->Fit("gaus", "r0", "", 5.37 - fHmass1[i]->GetRMS(), 5.37 + fHmass1[i]->GetRMS());
      peak1V  = fHmass1[i]->GetFunction("gaus")->GetParameter(1);
      peak1E  = fHmass1[i]->GetFunction("gaus")->GetParError(1);
      sigma1V = fHmass1[i]->GetFunction("gaus")->GetParameter(2);
      sigma1E = fHmass1[i]->GetFunction("gaus")->GetParError(2);
      fHmass1[i]->GetFunction("gaus")->SetLineColor(kRed);
      p2->SetBinContent(i+1, peak1V);
      p2->SetBinError(i+1, peak1E);
      w2->SetBinContent(i+1, sigma1V);
      w2->SetBinError(i+1, sigma1E);
    }
    if (fHmass0[i]->GetSumOfWeights() > 0.) {
      fHmass0[i]->Fit("gaus", "r0", "", 5.37 - fHmass0[i]->GetRMS(), 5.37 + fHmass0[i]->GetRMS());
      peak0V  = fHmass0[i]->GetFunction("gaus")->GetParameter(1);
      peak0E  = fHmass0[i]->GetFunction("gaus")->GetParError(1);
      sigma0V = fHmass0[i]->GetFunction("gaus")->GetParameter(2);
      sigma0E = fHmass0[i]->GetFunction("gaus")->GetParError(2);
      fHmass0[i]->GetFunction("gaus")->SetLineColor(kBlue);
      p1->SetBinContent(i+1, peak0V);
      p1->SetBinError(i+1, peak0E);
      w1->SetBinContent(i+1, sigma0V);
      w1->SetBinError(i+1, sigma0E);
    }

    if (fHmass0[i]->GetSumOfWeights() > 0.) fHmass0[i]->Scale(1./fHmass0[i]->GetSumOfWeights());
    if (fHmass1[i]->GetSumOfWeights() > 0.) fHmass1[i]->Scale(1./fHmass1[i]->GetSumOfWeights());

    fHmass1[i]->SetMaximum(1.2*fHmass1[i]->GetMaximum()/fHmass1[i]->GetSumOfWeights());
    fHmass1[i]->Draw();
    fHmass0[i]->Draw("same");

    tl->SetTextSize(0.04); tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.2, 0.92, Form("%s", fHmass0[i]->GetTitle()));
    tl->SetTextSize(0.03); tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.65, 0.80, Form("RMS: %4.3f MeV", fHmass0[i]->GetRMS()));
    tl->SetTextSize(0.03); tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.65, 0.76, Form("peak: %5.4f MeV", peak0V));
    tl->SetTextSize(0.03); tl->SetTextColor(kRed);   tl->DrawLatexNDC(0.65, 0.70, Form("RMS: %4.3f MeV", fHmass1[i]->GetRMS()));
    tl->SetTextSize(0.03); tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.65, 0.66, Form("peak: %5.4f MeV", peak1V));
    savePad(Form("massRes-mass-bin%d.pdf", i));
  }

  c0->Clear();
  setHist(p1, kBlue, 24, 1.2);
  setHist(p2, kRed, 25, 1.2);
  p1->SetMinimum(5.3);
  p1->SetMaximum(5.4);
  setTitles(p1, "#it{|}#eta_{#it{f}}#it{|}", "#it{MPV}(m_{#it{#mu #mu}}) #it{[GeV]}", 0.05, 1.1, 1.6);
  p1->Draw("e");
  p2->Draw("esame");

  newLegend(0.25, 0.7, 0.45, 0.85);
  legg->SetTextSize(0.05);
  legg->AddEntry(p1, "Run 1", "p");
  legg->AddEntry(p2, "Run 2", "p");
  legg->Draw();

  savePad(Form("massRes-massMPV-vsEta.pdf"));


  c0->Clear();
  setHist(w1, kBlue, 24, 1.2);
  setHist(w2, kRed, 25, 1.2);
  w1->SetMinimum(0.);
  w1->SetMaximum(0.22);
  setTitles(w1, "#it{|}#eta_{#it{f}}#it{|}", "#sigma(m_{#it{#mu #mu}}) #it{[GeV]}", 0.05, 1.1, 1.6);
  w1->Draw("e");
  w2->Draw("esame");

  newLegend(0.25, 0.7, 0.45, 0.85);
  legg->SetTextSize(0.05);
  legg->AddEntry(w1, "Run 1", "p");
  legg->AddEntry(w2, "Run 2", "p");
  legg->Draw();

  savePad(Form("massRes-massSigma-vsEta.pdf"));


  c0->Clear();
  setHist(h1, kBlue, 24, 1.2);
  setHist(h2, kRed, 25, 1.2);
  h1->SetMinimum(5.3);
  h1->SetMaximum(5.4);
  setTitles(h1, "#it{|}#eta_{#it{f}}#it{|}", "#it{mean}(m_{#it{#mu #mu}}) #it{[GeV]}", 0.05, 1.1, 1.6);
  h1->Draw("e");
  h2->Draw("esame");

  newLegend(0.25, 0.7, 0.45, 0.85);
  legg->SetTextSize(0.05);
  legg->AddEntry(h1, "Run 1", "p");
  legg->AddEntry(h2, "Run 2", "p");
  legg->Draw();

  savePad(Form("massRes-massMean-vsEta.pdf"));


  c0->Clear();
  setHist(s1, kBlue, 24, 1.2);
  setHist(s2, kRed, 25, 1.2);
  setTitles(s1, "#it{|}#eta_{#it{f}}#it{|}", "#it{RMS}(m_{#it{#mu #mu}}) #it{[GeV]}", 0.05, 1.1, 1.7);
  s1->SetMinimum(0.);
  s1->SetMaximum(0.22);
  s1->Draw("e");
  s2->Draw("esame");

  newLegend(0.25, 0.7, 0.45, 0.85);
  legg->SetTextSize(0.05);
  legg->AddEntry(s1, "Run 1", "p");
  legg->AddEntry(s2, "Run 2", "p");
  legg->Draw();

  savePad(Form("massRes-massRms-vsEta.pdf"));
}


// ----------------------------------------------------------------------
// massResolution
void plotStuff::loopFunction5() {
  if (fb.m2pt < 4.) return;
  if (fb.pvn  > 10) return;

  double meta = fb.m1eta;
  if (TMath::Abs(meta) < TMath::Abs(fb.m2eta)) meta = fb.m2eta;

  int ieta = TMath::Abs(meta)/0.1;
  if (ieta > 24) {
    //    cout << "eta = " << TMath::Abs(meta) << " -> " << ieta << endl;
    return;
  }
  if (0 == fDBX) {
    fHmass0[ieta]->Fill(fb.m);
  } else {
    fHmass1[ieta]->Fill(fb.m);
  }

}


// ----------------------------------------------------------------------
void plotStuff::yieldStabilityRatios(string trgname) {

  vector<pair<string, string> > overlays;
  overlays.push_back(make_pair("bmmData",    "bupsikData"));

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

	savePad(Form("ys-ratio-%s-%s-%s-chan%d.pdf", trgname.c_str(), it->first.c_str(), it->second.c_str(), ichan));
      }
    }

  }



}



// ----------------------------------------------------------------------
void plotStuff::loopFunction1() {

  bool goodRun(false);

  if (fLargeRuns.end() != find(fLargeRuns.begin(), fLargeRuns.end(), fb.run)) {
    goodRun = true;
    //    cout << "found " << fb.run << " in fLargeRuns" << endl;
  } else {
    //    cout << "NOT found " << fb.run << " in fLargeRuns" << endl;
  }

  if (!goodRun) return;

  if (fb.m1pt < 4.0) return;
  if (fb.m2pt < 4.0) return;

  double m = fb.m;
  if ((fMode == BU2JPSIKP) || (fMode == BD2JPSIKSTAR) || (fMode == BS2JPSIPHI)) {
    if (fb.mpsi < 3.04) return;
    if (fb.mpsi > 3.15) return;

    if (fb.psipt   < 7.0) return;
    if (fb.psiprob < 0.1) return;
    if (fMode == BU2JPSIKP) {
      if (fb.kpt < 0.70) return;
    }

    if (fMode == BS2JPSIPHI) {
      if (fb.mkk   < 1.01) return;
      if (fb.mkk   > 1.03) return;
      if (fb.phidr > 0.30) return;
      if (fb.k1pt  < 0.80) return;
      if (fb.k2pt  < 0.80) return;
    }

    if (fMode == BD2JPSIKSTAR) {
      if (!fb.kstarfail) return;
      if (fb.kpt < 0.70) return;
      if (fb.pipt < 0.70) return;
      if (fb.mkpi < 0.86) return;
      if (fb.mkpi > 0.94) return;
    }

  }

  char hname[200];
  if (fYear < 2013.) fb.ps = 1;

  if (fb.json && fb.m1q*fb.m2q<0 && fb.pvw8 > 0.7 && fb.gmuid) {
    // do nothing
  } else {
    return;
  }

  if ((fSample == "ysfill") || (fSample == "NOC")) {
    // no cut
    sprintf(hname, "NOC_%d_chan%d", static_cast<int>(fb.run), fChan);
    if (0 == fYieldHists.count(hname)) fYieldHists.insert(make_pair(hname, new TH2D(hname, hname, 90, 5.0, 5.9, MAXPS+1, -1., MAXPS)));
    fYieldHists[Form("NOC_%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, -0.1, static_cast<double>(fb.ps));
    fYieldHists[Form("NOC_%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, 0.1);
    fYieldHists[Form("NOC_%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, fb.ps+0.1);
  }

  if (fb.fls3d < 4.) return;
  if (fb.chi2dof  > fCuts[fChan]->chi2dof) return;
  if (fb.alpha    > 0.2) return;
  if (fb.iso      < fCuts[fChan]->iso) return;

  if ((fSample == "ysfill") || (fSample == "TOS")) {
    if (fb.hlt1 && fb.tos) {
      sprintf(hname, "TOS_%d_chan%d", static_cast<int>(fb.run), fChan);
      if (0 == fYieldHists.count(hname)) fYieldHists.insert(make_pair(hname, new TH2D(hname, hname, 90, 5.0, 5.9, MAXPS+1, -1., MAXPS)));
      fYieldHists[Form("TOS_%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, -0.1, static_cast<double>(fb.ps));
      fYieldHists[Form("TOS_%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, 0.1);
      fYieldHists[Form("TOS_%d_chan%d", static_cast<int>(fb.run), fChan)]->Fill(m, fb.ps+0.1);
    }
  }

}


// ----------------------------------------------------------------------
void plotStuff::loopFunction2() {

  if (fLargeRuns.end() != find(fLargeRuns.begin(), fLargeRuns.end(), fb.run)) {
  } else {
    // -- comment the following line if ALL runs should be looked at
    //    return;
  }

  char hname0[200];
  char hname1[200];
  sprintf(hname0, "ls0_%s_%d", fSample.c_str(), static_cast<int>(fb.run));
  sprintf(hname1, "ls1_%s_%d", fSample.c_str(), static_cast<int>(fb.run));

  if (0 == fHLs0.count(hname0)) {
    fHLs0.insert(make_pair(hname0, new TH2D(hname0, hname0, 3500, 0., 3500., 6, -1., 5.)));
    fHLs1.insert(make_pair(hname1, new TH2D(hname1, hname1, 3500, 0., 3500., 6, -1., 5.)));
  }


  if (!fb.json) return;
  if (!fb.hlt1) return;
  if (!fb.tos) return;
  if (fb.m1pt < 4.0) return;
  if (fb.m2pt < 4.0) return;
  if (TMath::Abs(fb.m1eta) > 1.4) return;
  if (TMath::Abs(fb.m2eta) > 1.4) return;
  if (fb.m1q * fb.m2q > 0) return;

  bool doubleMu0(false), highPtMu(false);
  if (((fb.l1s & 0x1) == 1) || ((fb.l1s & 0x2) == 2) || ((fb.l1s & 0x4) == 4)) {
    fHLs0[hname0]->Fill(fb.ls, -0.5);
    if (fb.chan > -1) fHLs0[hname0]->Fill(fb.ls, fb.chan);
  }

  if (((fb.l1s & 8) == 8) || ((fb.l1s & 16) == 16) || ((fb.l1s & 32) == 32) || ((fb.l1s & 64) == 64)) {
    fHLs1[hname1]->Fill(fb.ls, -0.5);
    if (fb.chan > -1) fHLs1[hname1]->Fill(fb.ls, fb.chan);
  }

  fProf["npv"]->Fill(fb.run, fb.pvn);


}

// ----------------------------------------------------------------------
// -- this is the puStudy on the events tree!!!
void plotStuff::loopFunction3() {

  int ic = fb.chan;
  int idx[] = {4, ic};
  int imax = (ic > -1?2:1);
  // -- we are only interested in 'real' candidates (else you get strange effects, e.g. back-to-back low-pT muons with large eta-B)
  if (!fb.hlt1) return;
  if (fb.m1pt < 4.0) return;
  if (fb.m2pt < 4.0) return;

  // -- check for data BMM sidebands
  if (fMode == BSMM && fb.procid < 0) {
    if (fb.m < 5.45) return;
  }

  // -- check for data BU2JPSIKP sidebands
  if (fMode == BU2JPSIKP && fb.procid < 0) {
    if (fb.m > 5.35) return;
    if (fb.m < 5.20) return;
  }

  // -- fill per-channel and combined histograms
  for (int i = 0; i < imax; ++i) {
    if (TMath::Abs(fb.dzmin) > 1.0) fpHmultFar[idx[i]]->Fill(fb.pvntrk);
    if (TMath::Abs(fb.dzmin) < 0.05) fpHmultClose05[idx[i]]->Fill(fb.pvntrk);

    fpHpvn[idx[i]]->Fill(fb.pvn);
    fpHdzmin[idx[i]]->Fill(fb.dzmin);

    if (TMath::Abs(fb.pv2lip) > 1.0) fpHmultFar2[idx[i]]->Fill(fb.pvntrk);
    if (TMath::Abs(fb.pv2lip) < 0.05) fpHmultClose2[idx[i]]->Fill(fb.pvntrk);

    fpHlz1[idx[i]]->Fill(fb.pvlip);
    fpHlzs1[idx[i]]->Fill(fb.pvlips);
    fpHlz2[idx[i]]->Fill(fb.pv2lip);
    fpHlzs2[idx[i]]->Fill(fb.pv2lips);

    // -- vs dzmin
    fpP0Mult[idx[i]]->Fill(fb.dzmin, fb.pvntrk);
    fpP0Iso[idx[i]]->Fill(fb.dzmin, fb.iso);
    fpP0Ntrk[idx[i]]->Fill(fb.dzmin, fb.closetrk);
    fpP0PvIps[idx[i]]->Fill(fb.dzmin, fb.pvips);
    if (fb.bdt > -10.) fpP0Bdt[idx[i]]->Fill(fb.dzmin, fb.bdt);
    fpP0fls3d[idx[i]]->Fill(fb.dzmin, fb.fls3d);
    fpP0fl3d[idx[i]]->Fill(fb.dzmin, fb.fl3d);
    fpP0dfl3d[idx[i]]->Fill(fb.dzmin, fb.fl3d - fb.gfl3d);
    fpP0l1[idx[i]]->Fill(fb.dzmin, fb.pvlip);
    fpP0l2[idx[i]]->Fill(fb.dzmin, fb.pv2lip);
    fpP0npv[idx[i]]->Fill(fb.dzmin, fb.pvn);

    fpP0tau[idx[i]]->Fill(fb.dzmin,  fb.tau);
    fpP0dtau[idx[i]]->Fill(fb.dzmin, fb.tau - fb.gtau);

    // -- vs lz1
    fpP1Mult[idx[i]]->Fill(fb.pvlip, fb.pvntrk);
    fpP1Iso[idx[i]]->Fill(fb.pvlip, fb.iso);
    fpP1Ntrk[idx[i]]->Fill(fb.pvlip, fb.closetrk);
    fpP1PvIps[idx[i]]->Fill(fb.pvlip, fb.pvips);
    if (fb.bdt > -10.) fpP1Bdt[idx[i]]->Fill(fb.pvlip, fb.bdt);
    fpP1fls3d[idx[i]]->Fill(fb.pvlip, fb.fls3d);
    fpP1fl3d[idx[i]]->Fill(fb.pvlip, fb.fl3d);
    fpP1dfl3d[idx[i]]->Fill(fb.pvlip, fb.fl3d - fb.gfl3d);
    fpP1dzmin[idx[i]]->Fill(fb.pvlip, fb.dzmin);
    fpP1l2[idx[i]]->Fill(fb.pvlip, fb.pv2lip);
    fpP1npv[idx[i]]->Fill(fb.pvlip, fb.pvn);

    fpP1tau[idx[i]]->Fill(fb.pvlip,  fb.tau);
    fpP1dtau[idx[i]]->Fill(fb.pvlip, fb.tau - fb.gtau);

    // -- vs lz2
    fpP2Mult[idx[i]]->Fill(fb.pv2lip, fb.pvntrk);
    fpP2Iso[idx[i]]->Fill(fb.pv2lip, fb.iso);
    fpP2Ntrk[idx[i]]->Fill(fb.pv2lip, fb.closetrk);
    if (fb.bdt > -10.) fpP2Bdt[idx[i]]->Fill(fb.pv2lip, fb.bdt);
    fpP2PvIps[idx[i]]->Fill(fb.pv2lip, fb.pvips);
    fpP2fls3d[idx[i]]->Fill(fb.pv2lip, fb.fls3d);
    fpP2fl3d[idx[i]]->Fill(fb.pv2lip, fb.fl3d);
    fpP2dfl3d[idx[i]]->Fill(fb.pv2lip, fb.fl3d - fb.gfl3d);
    fpP2dzmin[idx[i]]->Fill(fb.pv2lip, fb.dzmin);
    fpP2l1[idx[i]]->Fill(fb.pv2lip, fb.pvlip);
    fpP2npv[idx[i]]->Fill(fb.pv2lip, fb.pvn);

    fpP2tau[idx[i]]->Fill(fb.pv2lip,  fb.tau);
    fpP2dtau[idx[i]]->Fill(fb.pv2lip, fb.tau - fb.gtau);


    // -- vs NPV
    fpP3Mult[idx[i]]->Fill(fb.pvn, fb.pvntrk);
    fpP3Iso[idx[i]]->Fill(fb.pvn, fb.iso);
    fpP3Ntrk[idx[i]]->Fill(fb.pvn, fb.closetrk);
    if (fb.bdt > -10.) fpP3Bdt[idx[i]]->Fill(fb.pvn, fb.bdt);
    fpP3PvIps[idx[i]]->Fill(fb.pvn, fb.pvips);
    fpP3fls3d[idx[i]]->Fill(fb.pvn, fb.fls3d);
    fpP3fl3d[idx[i]]->Fill(fb.pvn, fb.fl3d);
    fpP3dfl3d[idx[i]]->Fill(fb.pvn, fb.fl3d - fb.gfl3d);
    fpP3dzmin[idx[i]]->Fill(fb.pvn, TMath::Abs(fb.dzmin));
    fpP3l1[idx[i]]->Fill(fb.pvn, fb.pvlip);
    fpP3l2[idx[i]]->Fill(fb.pvn, TMath::Abs(fb.pv2lip));

    fpP3tau[idx[i]]->Fill(fb.pvn,  fb.tau);
    fpP3dtau[idx[i]]->Fill(fb.pvn, fb.tau - fb.gtau);


  }
}


// ----------------------------------------------------------------------
void plotStuff::loopFunction4() {
  if (!fb.hlt1) return;
  if (!fb.tos) return;
  if (0 != fb.chan) return;

  if (fb.flsxy < 4.) return;
  if (fb.m < 4.9) return;
  if (fb.m > 5.9) return;

  string prf("RRA_");

  if (fb.run > fSplitRun) {
    prf = "RRB_";
  }
  fPlots[prf + "run"]->Fill(0.001*fb.run);

  fPlots[prf + "m1pt"]->Fill(fb.m1pt);
  fPlots[prf + "m1eta"]->Fill(fb.m1eta);
  fPlots[prf + "m1phi"]->Fill(fb.m1phi);
  fPlots[prf + "m2pt"]->Fill(fb.m2pt);
  fPlots[prf + "m2eta"]->Fill(fb.m2eta);
  fPlots[prf + "m2phi"]->Fill(fb.m2phi);
  fPlots[prf + "m1pix"]->Fill(fb.m1pix);
  fPlots[prf + "m1bpix"]->Fill(fb.m1bpix);
  fPlots[prf + "m1bpixl1"]->Fill(fb.m1bpixl1);
  fPlots[prf + "m2pix"]->Fill(fb.m2pix);
  fPlots[prf + "m2bpix"]->Fill(fb.m2bpix);
  fPlots[prf + "m2bpixl1"]->Fill(fb.m2bpixl1);

  fPlots[prf + "m"]->Fill(fb.m);
  fPlots[prf + "pt"]->Fill(fb.pt);
  fPlots[prf + "eta"]->Fill(fb.eta);
  fPlots[prf + "phi"]->Fill(fb.phi);
  fPlots[prf + "tau"]->Fill(fb.tau);
  fPlots[prf + "taue"]->Fill(fb.taue);
  fPlots[prf + "bdt"]->Fill(fb.bdt);

  fPlots[prf + "lip"]->Fill(fb.lip);
  fPlots[prf + "lipe"]->Fill(fb.lipE);
  fPlots[prf + "tip"]->Fill(fb.tip);
  fPlots[prf + "tipe"]->Fill(fb.tipE);

  fPlots[prf + "closetrk"]->Fill(fb.closetrk);
  fPlots[prf + "pvlip"]->Fill(fb.pvlip);
  fPlots[prf + "pvlips"]->Fill(fb.pvlips);
  fPlots[prf + "maxdoca"]->Fill(fb.maxdoca);
  fPlots[prf + "pvip"]->Fill(fb.pvip);
  fPlots[prf + "pvips"]->Fill(fb.pvips);
  fPlots[prf + "alpha"]->Fill(fb.alpha);
  fPlots[prf + "fls3d"]->Fill(fb.fls3d);
  fPlots[prf + "fl3d"]->Fill(fb.fl3d);
  fPlots[prf + "fl3de"]->Fill(fb.fl3dE);
  fPlots[prf + "docatrk"]->Fill(fb.docatrk);
  fPlots[prf + "iso"]->Fill(fb.iso);
  fPlots[prf + "m1iso"]->Fill(fb.m1iso);
  fPlots[prf + "m2iso"]->Fill(fb.m2iso);
  fPlots[prf + "chi2dof"]->Fill(fb.chi2dof);

}

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
  if (ifunc == 5) pF = &plotStuff::loopFunction5;

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


// ----------------------------------------------------------------------
void plotStuff::runStudy(string ds) {

  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");

  TDirectory *hDir(0);
  hDir = gDirectory;

  bool ok = false;

  int runStart(273150), runEnd(284044);
  int runBins = runEnd - runStart + 1;

  TIter next(gDirectory->GetListOfKeys());
  TKey *key(0);
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH2D")) continue;
    if (TString(key->GetName()).Contains(Form("ls0_%s", ds.c_str()))) {
      ok = true;
      break;
    }
  }

  if (!ok) {
    cout << "Histograms not found" << endl;
    fProf.insert(make_pair("npv", new TProfile("npv", "npv", runBins, runStart, runEnd)));

    fSample = ds;
    setup(fSample);
    fCds = fDS[fSample];
    TTree *t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    loopOverTree(t, 2);

    if (1) {
      cout << "writing output histograms: " << fYieldHists.size() << endl;

      fProf["npv"]->SetDirectory(hDir);
      fProf["npv"]->Write();

      for (map<string, TH2D*>::iterator it = fHLs0.begin(); it != fHLs0.end(); ++it) {
	cout << "name = " << it->first << " hist " << it->second->GetName() << endl;
	it->second->SetDirectory(hDir);
	it->second->Write();
      }

      for (map<string, TH2D*>::iterator it = fHLs1.begin(); it != fHLs1.end(); ++it) {
	cout << "name = " << it->first << " hist " << it->second->GetName() << endl;
	it->second->SetDirectory(hDir);
	it->second->Write();
      }
    }
  }

  Lumi *lumi(0);
  if (2016 == fYear) {
    if (string::npos != ds.find("psi")) {
      //	lumi = new Lumi("../common/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
      lumi = new Lumi("../common/json/2016-HLT_DoubleMu4_3_Jpsi_Displaced.lumi");
    } else {
      lumi = new Lumi("../common/json/2016-HLT_DoubleMu4_3_Bs.lumi");
    }
  } else if (2012 == fYear) {
    lumi = new Lumi("../common/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.lumi");
  } else if (2011 == fYear) {
    lumi = new Lumi("../common/json/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_MuonPhys_v2.lumi");
  }


  shrinkPad(0.1, 0.1, 0.07, 0.12);
  TH1D *hl = new TH1D("lostLs", "lost lumi sections", runBins, runStart, runEnd);
  TH1D *h0y = new TH1D("yields", "yields low pT", runBins, runStart, runEnd);
  TH1D *h1y = new TH1D("yields", "yields high pT", runBins, runStart, runEnd);
  TH2D *hc0b = new TH2D("hc0b", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc0c = new TH2D("hc0c", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc0d = new TH2D("hc0d", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc0e = new TH2D("hc0e", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc0f = new TH2D("hc0f", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc0g = new TH2D("hc0g", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc0h = new TH2D("hc0h", "", 50, 0., 50., 50, 0., 1200.);

  TH2D *hc1b = new TH2D("hc1b", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc1c = new TH2D("hc1c", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc1d = new TH2D("hc1d", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc1e = new TH2D("hc1e", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc1f = new TH2D("hc1f", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc1g = new TH2D("hc1g", "", 50, 0., 50., 50, 0., 1200.);
  TH2D *hc1h = new TH2D("hc1h", "", 50, 0., 50., 50, 0., 1200.);

  TProfile *hp = (TProfile*)gDirectory->Get("npv");
  next = gDirectory->GetListOfKeys();
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH2D")) continue;
    if (TString(key->GetName()).Contains(Form("ls0_%s_", ds.c_str()))) {
      string hname = key->GetName();
      TH2D *h0 = (TH2D*)((TH2D*)gDirectory->Get(hname.c_str()));
      replaceAll(hname, "ls0", "ls1");
      TH2D *h1 = (TH2D*)((TH2D*)gDirectory->Get(hname.c_str()));
      TH1D *h0l = h0->ProjectionX(Form("%s_x0", hname.c_str()), 1,1); h0l->SetLineColor(kBlue);
      TH1D *h1l = h1->ProjectionX(Form("%s_x1", hname.c_str()), 1,1); h1l->SetLineColor(kRed);
      h0l->Rebin(10);
      h1l->Rebin(10);
      replaceAll(hname, Form("ls1_%s_", ds.c_str()), "");
      int run = atoi(hname.c_str());
      string rr = runRange(run);
      double hltLumi = lumi->lumi(run);
      if (hltLumi > 10.) {
	if ("B" == rr) {
	  hc0b->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1b->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("C" == rr) {
	  hc0c->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1c->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("D" == rr) {
	  hc0d->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1d->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("E" == rr) {
	  hc0e->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1e->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("F" == rr) {
	  hc0f->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1f->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("G" == rr) {
	  hc0g->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1g->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("H" == rr) {
	  hc0h->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1h->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	}
	h0y->SetBinContent(h1y->FindBin(run), h0l->Integral()/hltLumi);
	h1y->SetBinContent(h1y->FindBin(run), h1l->Integral()/hltLumi);

	cout << "run = " << run
	     << " lumi = " << hltLumi
	     << " nPV: " << hp->GetBinContent(hp->FindBin(run))
	     << " yield(0) = " <<  h0l->Integral()
	     << " yield(0)*pb = " <<  h0l->Integral()/hltLumi
	     << " yield(1) = " <<  h1l->Integral()
	     << " yield(1)*pb = " <<  h1l->Integral()/hltLumi
	     << endl;
      }
      // -- skip initial runs for the missed ls count
      if (run > 274000) {
	int lsCnt(0), lsTot(0);
	double ymax = h1l->GetMaximum();
	for (int ib = 0; ib < h0l->GetNbinsX(); ++ib) {
	  if (h1l->GetBinContent(ib) > 0.05*ymax) {
	    ++lsTot;
	    if (h0l->GetBinContent(ib) < h1l->GetBinContent(ib)) ++lsCnt;
	  }
	}
	if (lsCnt>1) {
	  cout << "XXXX Run " << run << " missed ls fraction: " << static_cast<double>(lsCnt)/lsTot << endl;
	  hl->SetBinContent(hl->FindBin(run), static_cast<double>(lsCnt)/lsTot);
	}
      }
      TH1D *hr = (TH1D*)h1l->Clone("h1r");
      hr->Divide(h0l, h1l, 1., 1.);
      hr->Scale(0.6*h0l->GetMaximum()/hr->GetMaximum());
      hr->SetLineColor(kBlack);
      h0l->Draw("hist");
      h1l->Draw("samehist");
      hr->Draw("samehist");
      tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.2, 0.93, h0l->GetName());
      tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.2, 0.90, h1l->GetName());
      savePad(Form("runStudy_%s_l1seeds_%s.pdf", ds.c_str(), hname.c_str()));
    }
  }

  hl->SetMinimum(0.);
  hl->Draw();
  plot2016Eras(hl->GetMaximum());
  savePad(Form("runStudy_%s_missedLS.pdf", ds.c_str()));
  h0y->SetMinimum(0.);
  h0y->Draw("e");
  plot2016Eras(h0y->GetMaximum());
  savePad(Form("runStudy_%s_yields_DoubleMu0.pdf", ds.c_str()));
  h1y->SetMinimum(0.);
  h1y->Draw("e");
  plot2016Eras(h1y->GetMaximum());
  savePad(Form("runStudy_%s_yields_highPt.pdf", ds.c_str()));

  cout << "MAXIMUM: " << hc0b->GetMaximum() << endl;
  //  hc0b->SetMaximum(1.2*hc0b->GetMaximum());
  hc0b->SetMarkerColor(kBlack);     hc0b->Draw();
  hc0c->SetMarkerColor(kGreen+1);   hc0c->Draw("same");
  hc0d->SetMarkerColor(kRed+1);     hc0d->Draw("same");
  hc0e->SetMarkerColor(kRed+2);     hc0e->Draw("same");
  hc0f->SetMarkerColor(kRed+3);     hc0f->Draw("same");
  hc0g->SetMarkerColor(kBlue+1);    hc0g->Draw("same");
  hc0h->SetMarkerColor(kBlue+3);    hc0h->Draw("same");

  tl->SetTextColor(kBlack);   tl->DrawLatexNDC(0.3, 0.92, "DoubleMu0 seeds");
  tl->SetTextColor(kBlack);   tl->DrawLatexNDC(0.7, 0.80, "2016B");
  tl->SetTextColor(kGreen+1); tl->DrawLatexNDC(0.7, 0.72, "2016C");
  tl->SetTextColor(kRed+1);   tl->DrawLatexNDC(0.7, 0.64, "2016D");
  tl->SetTextColor(kRed+2);   tl->DrawLatexNDC(0.7, 0.56, "2016E");
  tl->SetTextColor(kRed+3);   tl->DrawLatexNDC(0.7, 0.48, "2016F");
  tl->SetTextColor(kBlue+1);  tl->DrawLatexNDC(0.7, 0.40, "2016G");
  tl->SetTextColor(kBlue+3);  tl->DrawLatexNDC(0.7, 0.32, "2016H");

  savePad(Form("runStudy_%s_yieldVsNpv_Seed0_allEras.pdf", ds.c_str()));

  hc1b->SetMarkerColor(kBlack);     hc1b->Draw();
  hc1c->SetMarkerColor(kGreen+1);   hc1c->Draw("same");
  hc1d->SetMarkerColor(kRed+1);     hc1d->Draw("same");
  hc1e->SetMarkerColor(kRed+2);     hc1e->Draw("same");
  hc1f->SetMarkerColor(kRed+3);     hc1f->Draw("same");
  hc1g->SetMarkerColor(kBlue+1);    hc1g->Draw("same");
  hc1h->SetMarkerColor(kBlue+3);    hc1h->Draw("same");

  tl->SetTextColor(kBlack);   tl->DrawLatexNDC(0.3, 0.92, "DoubleMu_1X_Y seeds");
  tl->SetTextColor(kBlack);   tl->DrawLatexNDC(0.7, 0.80, "2016B");
  tl->SetTextColor(kGreen+1); tl->DrawLatexNDC(0.7, 0.72, "2016C");
  tl->SetTextColor(kRed+1);   tl->DrawLatexNDC(0.7, 0.64, "2016D");
  tl->SetTextColor(kRed+2);   tl->DrawLatexNDC(0.7, 0.56, "2016E");
  tl->SetTextColor(kRed+3);   tl->DrawLatexNDC(0.7, 0.48, "2016F");
  tl->SetTextColor(kBlue+1);  tl->DrawLatexNDC(0.7, 0.40, "2016G");
  tl->SetTextColor(kBlue+3);  tl->DrawLatexNDC(0.7, 0.32, "2016H");

  savePad(Form("runStudy_%s_yieldVsNpv_Seed1_allEras.pdf", ds.c_str()));

  delete hl;
  delete h0y;
  delete h1y;
  delete hc0b;
  delete hc0c;
  delete hc0d;
  delete hc0e;
  delete hc0f;
  delete hc0g;
  delete hc0h;

  delete hc1b;
  delete hc1c;
  delete hc1d;
  delete hc1e;
  delete hc1f;
  delete hc1g;
  delete hc1h;

  fHistFile->Close();

}


// ----------------------------------------------------------------------
void plotStuff::tauEfficiency(string varname, string cut, string otherSelection, string dsname) {

  if (varname == "all") {
    for (int i = 0; i < 2; ++i) {
      tauEfficiency(Form("flsxy_chan%d", i), "flsxy>10", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("fls3d_chan%d", i), "fls3d>15", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("alpha_chan%d", i), "alpha<0.05", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");

      tauEfficiency(Form("chi2dof_chan%d", i), "chi2dof<2", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("maxdoca_chan%d", i), "maxdoca<0.02", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");

      tauEfficiency(Form("docatrk_chan%d", i), "docatrk>0.015", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("closetrk_chan%d", i), "closetrk<3", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");

      tauEfficiency(Form("iso_chan%d", i), "iso>0.8", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("m1iso_chan%d", i), "m1iso>0.8", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("m2iso_chan%d", i), "m2iso>0.8", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");

      tauEfficiency(Form("pvip_chan%d", i), "pvip<0.01", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("pvips_chan%d", i), "pvips<5", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");


      tauEfficiency("hlt1_chan0", "hlt1", "m2pt>4.&&(chan==0)", "bsmmMcComb");
      tauEfficiency("hlt1_chan1", "hlt1", "m2pt>4.&&(chan==1)", "bsmmMcComb");
      tauEfficiency("hlt1_chan2", "hlt1", "m2pt>4.&&(chan==2)", "bsmmMcComb");
      tauEfficiency("hlt1_chan3", "hlt1", "m2pt>4.&&(chan==3)", "bsmmMcComb");


      tauEfficiency(Form("m1pt_chan%d", i), "m1pt>8", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("m2pt_chan%d", i), "m2pt>6", Form("m1pt>6.&&(chan==%d)", i), "bsmmMcComb");

      tauEfficiency(Form("cnc_chan%d", i), "cnc", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("bdt_chan%d", i), "bdt>0.34", Form("m2pt>4.&&(chan==%d)", i), "bsmmMcComb");

      tauEfficiency(Form("gmuid_chan%d", i), "gmuid", Form("m1pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("gtqual_chan%d", i), "gtqual", Form("m1pt>4.&&(chan==%d)", i), "bsmmMcComb");

      tauEfficiency(Form("m1bpix_chan%d", i), "m1bpix>0", Form("m1pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("m2bpix_chan%d", i), "m2bpix>0", Form("m1pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("m1bpixl1_chan%d", i), "m1bpixl1>0", Form("m1pt>4.&&(chan==%d)", i), "bsmmMcComb");
      tauEfficiency(Form("m2bpixl1_chan%d", i), "m2bpixl1>0", Form("m1pt>4.&&(chan==%d)", i), "bsmmMcComb");
    }
    return;
  }


  fSample = dsname;
  cout << "==> plotWork::efficiencyVariable> sample = " << fSample << endl;;

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

  string normName = Form("tau_%s_norm", fSample.c_str());
  string passName = Form("tau_%s_pass", fSample.c_str());
  string effName  = Form("tau_%s_eff", fSample.c_str());
  TH1D *h1 = new TH1D(normName.c_str(), normName.c_str(), 48, 0., 12.e-12);
  h1->Sumw2();
  setTitles(h1, "#tau #it{[ps]}", "");
  setFilledHist(h1, kBlue, kYellow, 1000, 2);
  TH1D *h2 = new TH1D(passName.c_str(), passName.c_str(), 48, 0., 12.e-12);
  h2->Sumw2();
  setFilledHist(h2, kBlue, kBlue, 3354, 2);

  // -- basic HLT efficiency derived from MC
  string tselection = otherSelection;
  t->Draw(Form("tau >> %s", normName.c_str()), tselection.c_str());
  cout << "==> " << tselection << " histogram contents =        " << h1->Integral(1, h1->GetNbinsX()+1) << endl;
  tselection = cut + " && " + otherSelection;
  t->Draw(Form("tau >> %s", passName.c_str()), tselection.c_str());
  cout << "==> " << tselection << " histogram contents = " << h2->Integral(1, h1->GetNbinsX()+1) << endl;

  TH1D *h3 = (TH1D*)(h1->Clone(effName.c_str())); h3->Reset();
  setHist(h3);
  h3->Divide(h2, h1, 1., 1., "b");
  setTitles(h3, "#tau #it{[ps]}", "Efficiency");

  gPad->SetTicks(1,0);
  gPad->SetTopMargin(0.11);
  gPad->SetRightMargin(0.15);
  h1->Draw("hist");
  h2->Draw("histsame");
  h1->Draw("axissame");
  tl->SetTextSize(0.027);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.5, 0.84, Form("cuts: %s", cut.c_str()));
  tl->DrawLatexNDC(0.5, 0.80, Form("base: %s", otherSelection.c_str()));
  tl->SetTextSize(0.05);
  tl->DrawLatexNDC(0.2, 0.92, varname.c_str());

  double xmin(0.e-12);
  if (h3->GetBinContent(2) < 0.4) xmin = 4.e-12;
  double sf = 0.7*h1->GetMaximum();
  h3->Scale(sf);
  h3->SetMarkerColor(kBlue);
  h3->SetLineColor(kBlue);
  h3->Fit("pol1", "r", "same", xmin, 12.e-12);
  h3->GetFunction("pol1")->SetLineWidth(2);
  h3->GetFunction("pol1")->SetLineColor(kBlue);

  tl->SetTextSize(0.027);
  tl->SetTextColor(kBlue);
  tl->DrawLatexNDC(0.5, 0.76,
		   Form("slope = (%3.1f #pm %3.1f)e9", 1.e-9/sf*h3->GetFunction("pol1")->GetParameter(1), 1.e-9/sf*h3->GetFunction("pol1")->GetParError(1)));

  // draw axis on the right side of the pad
  TGaxis *axis = new TGaxis(1.0*h3->GetXaxis()->GetXmax(), 0., 1.0*h3->GetXaxis()->GetXmax(), h1->GetMaximum(), 0.01, 1./0.7, 510, "+L");
  axis->SetLabelColor(kBlue);
  axis->SetLineColor(kBlue);
  axis->Draw();


  c0->cd();
  savePad(Form("tauEff%d_%s_%s.pdf", fYear, varname.c_str(), fSample.c_str()));


}


// ----------------------------------------------------------------------
void plotStuff::yieldStudy(int run, string ds) {
  fSplitRun = run;

  TH1::SetDefaultSumw2(kTRUE);
  vector<string> range;
  range.push_back("RRA_");
  range.push_back("RRB_");
  string var("");
  for (unsigned int i = 0; i < range.size(); ++i) {
    var = range[i] + "run"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 13, 271., 284.)));
    var = range[i] + "m1pt"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 40, 0., 40.)));
    var = range[i] + "m1eta"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, -2.5, 2.5)));
    var = range[i] + "m1phi"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, -3.15, 3.15)));
    var = range[i] + "m2pt"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 40, 0., 20.)));
    var = range[i] + "m2eta"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, -2.5, 2.5)));
    var = range[i] + "m2phi"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, -3.15, 3.15)));
    var = range[i] + "m1pix"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 5, 0., 5)));
    var = range[i] + "m1bpix"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 5, 0., 5)));
    var = range[i] + "m1bpixl1"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 5, 0., 5)));
    var = range[i] + "m2pix"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 5, 0., 5)));
    var = range[i] + "m2bpix"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 5, 0., 5)));
    var = range[i] + "m2bpixl1"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 5, 0., 5)));


    var = range[i] + "m"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 4.9, 5.9)));
    var = range[i] + "pt"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 50.)));
    var = range[i] + "eta"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, -2.5, 2.5)));
    var = range[i] + "phi"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, -3.15, 3.15)));
    var = range[i] + "tau"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 100, 0., 1.5e-12)));
    var = range[i] + "taue"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 100, 0., 0.2e-12)));
    var = range[i] + "bdt"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 40, -1., 1.)));

    var = range[i] + "lip"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 40, 0., 0.1)));
    var = range[i] + "lipe"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 40, 0., 0.01)));
    var = range[i] + "tip"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 40, 0., 0.04)));
    var = range[i] + "tipe"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 40, 0., 0.005)));

    var = range[i] + "closetrk"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 20, 0., 20.)));
    var = range[i] + "pvlip"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 0.02)));
    var = range[i] + "pvlips"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 5.)));
    var = range[i] + "maxdoca"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 0.05)));
    var = range[i] + "pvip"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 0.04)));
    var = range[i] + "pvips"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 5)));

    var = range[i] + "alpha"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 0.2)));
    var = range[i] + "fls3d"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 60, 0., 40.)));
    var = range[i] + "fl3d"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 0.8)));
    var = range[i] + "fl3de"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 0.04)));
    var = range[i] + "docatrk"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 0.08)));
    var = range[i] + "iso"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 1.02)));
    var = range[i] + "m1iso"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 1.02)));
    var = range[i] + "m2iso"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 1.02)));
    var = range[i] + "chi2dof"; fPlots.insert(make_pair(var, new TH1D(Form("%s", var.c_str()), Form("%s", var.c_str()), 50, 0., 5.)));

  }

  fSample = ds;
  setup(fSample);
  fCds = fDS[fSample];
  TTree *t = getTree(fSample, fTreeDir);
  setupTree(t, fSample);
  loopOverTree(t, 4);

  c0->Clear();
  c0->Divide(1,2);

  for (map<string, TH1D*>::iterator it = fPlots.begin(); it != fPlots.end(); ++it) {
    string A = it->first;
    if (string::npos != A.find("RRB_")) continue;
    cout << "A = " << A << endl;
    string B = A;
    replaceAll(B, "RRA_", "RRB_");
    cout << "B = " << B << endl;
    TH1D *ha = it->second;
    TH1D *hb = fPlots[B];
    cout << "ha = " << ha << " hb = " << hb << endl;
    setFilledHist(hb, kBlue, kBlue, 3354);
    setFilledHist(ha, kRed, kRed, 3365);

    ha->Scale(1./ha->GetSumOfWeights());
    hb->Scale(1./hb->GetSumOfWeights());
    if (hb->GetMaximum() > ha->GetMaximum()) ha->SetMaximum(1.2*hb->GetMaximum());

    c0->cd(1);
    ha->SetMinimum(0.);
    ha->Draw("hist");
    hb->Draw("histsame");
    tl->DrawLatexNDC(0.24, 0.02, fSample.c_str());
    tl->DrawLatexNDC(0.60, 0.02, B.c_str());

    newLegend(0.71, 0.7, 0.81, 0.90);
    legg->SetTextSize(0.05);
    legg->AddEntry(ha, Form("< %d", fSplitRun), "f");
    legg->AddEntry(hb, Form("> %d", fSplitRun), "f");
    legg->Draw();


    TH1D *hratio = (TH1D*)ha->Clone("ratio");
    hratio->Reset();
    hratio->Divide(hb, ha, 1., 1.);
    setHist(hratio, kBlack, kBlack);
    c0->cd(2);
    hratio->SetMinimum(0.);
    hratio->SetMaximum(1.5);
    hratio->Draw("e");
    replaceAll(B, "RRB_", "");
    savePad(Form("yieldStudy-%s-%d-%s.pdf", B.c_str(), fSplitRun, fSample.c_str()), c0);
  }
}



// ----------------------------------------------------------------------
void plotStuff::plotWrongReco(string var, int nbin, double min, double max, string selection,
			     string wds, string wdir, string cds, string cdir) {

  string name = var + wds + wdir;
  TH1D *h1 = new TH1D(name.c_str(), name.c_str(), nbin, min, max);
  setFilledHist(h1, kRed, kRed, 3365);
  setTitles(h1, fVarToTex[var].c_str(), "Entries / Bin", 0.05, 1.2, 1.5, 0.05, 52);
  TTree *t = getTree(wds, wdir);
  t->Draw(Form("%s>>%s", var.c_str(), name.c_str()), selection.c_str());

  name = var + cds + cdir;
  TH1D *h2 = new TH1D(name.c_str(), name.c_str(), nbin, min, max);
  setFilledHist(h2, kBlue, kBlue, 3354);
  setTitles(h2, fVarToTex[var].c_str(), "Entries / Bin", 0.05, 1.2, 1.5, 0.05, 52);
  t = getTree(cds, cdir);
  t->Draw(Form("%s>>%s", var.c_str(), name.c_str()), selection.c_str());

  double int1 = h1->Integral();
  double int2 = h2->Integral();
  h1->Scale(1./int1);
  h2->Scale(1./int2);

  double ymax(h1->GetMaximum());
  if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();
  ymax *= 1.2;

  shrinkPad(0.15, 0.18);
  h1->SetMaximum(ymax);
  h1->Draw();
  h2->Draw("samehist");

  // setRoman();
  setItalic();
  tl->DrawLatexNDC(0.2, 0.92, Form("%s candidates", fDS[cds]->fName.c_str()));

  string wrg("");
  if (wdir == "candAnaBd2JpsiKstarAsBu") {
    wrg = "bdpsikstarMcComb";
  }
  newLegend(0.55, 0.7, 0.75, 0.87);
  legg->SetHeader("true decay");
  legg->SetTextSize(0.04);
  legg->AddEntry(h1, fDS[wrg]->fName.c_str(), "f");
  legg->AddEntry(h2, fDS[cds]->fName.c_str(), "f");
  legg->Draw();

  savePad(Form("plotWrongReco_%s_%s_%s.pdf", var.c_str(), wds.c_str(), cds.c_str()));
}



// ----------------------------------------------------------------------
void plotStuff::wrongReco(string ds1, string mode, string selection) {

  setItalic();

  string mapname = ds1 + "-" + mode.substr(string("candAna").size());
  cout << "==> mapname: " << mapname << endl;

  shrinkPad(0.15, 0.2);
  string name(Form("H1_%s", mapname.c_str()));
  string xtitle("m_{#it{#mu#mu} K^{+}} #it{[GeV]}");
  if (string::npos != mapname.find("AsBs")) xtitle = "m_{#it{#mu#mu} K^{+}K^{-}} #it{[GeV]}";
  if (string::npos != mapname.find("bc")) xtitle = "m_{#it{#mu#mu}} #it{[GeV]}";
  TH1D *h1(0);
  if (string::npos != ds1.find("bupsipi")) {
    h1 = new TH1D(name.c_str(), name.c_str(), 100, 4.8, 6.8);
  } else {
    h1 = new TH1D(name.c_str(), name.c_str(), 100, 5.0, 6.0);
  }
  setTitles(h1, xtitle.c_str(), "Entries / Bin", 0.05, 1.2, 2.0, 0.05, 52);
  TTree *t = getTree(ds1, mode);
  t->Draw(Form("m>>%s", name.c_str()), selection.c_str());

  if (mode == "candAnaBd2JpsiKstarAsBu") {
    tl->SetTextSize(0.05);  tl->DrawLatexNDC(0.2, 0.92, fDS["bdpsikstarMcComb"]->fName.c_str());
    TF1 *f1 = fIF->expoErr(5.0, 6.0);
    double preco(5.145);
    double e0(preco),  e0Min(preco-0.01), e0Max(preco+0.01);
    double e1(0.075),  e1Min(0.050), e1Max(0.100);
    double e2(1.15), e2Min(1.05),  e2Max(1.25);
    double e3(h1->GetMaximum());
    double p0, p1;
    fIF->fLo = 5.25;
    fIF->fHi = 5.6;
    fIF->initExpo(p0, p1, h1);

    f1->SetParameters(p0, p1, e0, e1, e2, e3);
    f1->FixParameter(0, p0);
    f1->SetLineWidth(2);
    fIF->fLo = 5.0;
    fIF->fHi = 6.0;
    h1->Fit(f1, "lr", "", 5.02, 6.0);
  }

  if (mode == "candAnaBd2JpsiKstarAsBs") {
    tl->SetTextSize(0.05);  tl->DrawLatexNDC(0.2, 0.92, fDS["bdpsikstarMcComb"]->fName.c_str());
    fIF->fLo = 5.0;
    fIF->fHi = 6.0;
    //      TF1 *f1 = fIF->pol1Landau(h1, 5.4, 0.1);
    //      TF1 *f1 = fIF->pol1gauss(h1, 5.4, 0.1);
    TF1 *f1 = fIF->pol1gauss2(h1, 5.4, 0.1, 0.05, 0.05);

    f1->SetLineWidth(2);
    h1->Fit(f1, "lr", "", 5.0, 6.0);
  }

  if (mode == "candAnaMuMu") {
    tl->SetTextSize(0.05);  tl->DrawLatexNDC(0.2, 0.92, fDS["bcpsimunuMcCombBg"]->fName.c_str());
    fIF->fLo = 5.0;
    fIF->fHi = 6.0;
    //      TF1 *f1 = fIF->pol1Landau(h1, 5.4, 0.1);
    //      TF1 *f1 = fIF->pol1gauss(h1, 5.4, 0.1);
    TF1 *f1 = fIF->expo(h1);

    f1->SetLineWidth(2);
    h1->Fit(f1, "lr", "", 5.0, 6.0);
  }

  if (mode == "candAnaBu2JpsiK") {
    tl->SetTextSize(0.05);  tl->DrawLatexNDC(0.2, 0.92, fDS[ds1]->fName.c_str());
    fIF->fLo = 4.9;
    fIF->fHi = 6.4;
    TF1 *f1 = fIF->gauss3(h1);

    f1->SetLineWidth(2);
    h1->Fit(f1, "lr", "", 4.9, 6.4);
  }

  tl->SetTextAngle(90.);
  tl->SetTextSize(0.02);
  setRoman();
  tl->DrawLatexNDC(0.95, 0.17, selection.c_str());
  tl->SetTextAngle(0.);

  savePad(Form("wrongReco-%s.pdf", mapname.c_str()));
}



// ----------------------------------------------------------------------
void plotStuff::o2Profile(TProfile *p, string dsname, string dsname2, int i) {
  string sname = p->GetName();
  sname = sname.substr(0, sname.find("_"));
  TProfile *s = (TProfile*)fHistFile->Get(Form("%s_%s_chan%d", sname.c_str(), dsname2.c_str(), i));
  s->SetMarkerStyle(24);
  s->SetMarkerColor(kRed);
  s->SetLineColor(kRed);
  s->Draw("same");
  newLegend(0.6, 0.8, 0.9, 0.90);
  legg->SetTextSize(0.03);
  legg->AddEntry(p, Form("%s", dsname.c_str()), "p");
  legg->AddEntry(s, Form("%s", dsname2.c_str()), "p");
  legg->Draw();
}


// ----------------------------------------------------------------------
void plotStuff::plot2016Eras(double ymax) {
  tl->SetNDC(false);
  pl->SetLineColor(kRed);
  tl->DrawLatex(273150., 1.0, "B"); pl->DrawLine(273150., 0., 273150., ymax);
  tl->DrawLatex(275657., 1.0, "C"); pl->DrawLine(275657., 0., 275657., ymax);
  tl->DrawLatex(276315., 1.0, "D"); pl->DrawLine(276315., 0., 276315., ymax);
  tl->DrawLatex(276831., 1.0, "E"); pl->DrawLine(276831., 0., 276831., ymax);
  tl->DrawLatex(277772., 1.0, "F"); pl->DrawLine(277772., 0., 277772., ymax);
  tl->DrawLatex(278820., 1.0, "G"); pl->DrawLine(278820., 0., 278820., ymax);
  tl->DrawLatex(280919., 1.0, "H"); pl->DrawLine(280919., 0., 280919., ymax);
}
