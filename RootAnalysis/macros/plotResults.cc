#include "plotResults.hh"

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
#include "THStack.h"
#include "TFitResult.h"

#include "common/dataset.hh"
#include "common/util.hh"
#include "common/Lumi.hh"
#include "common/fitPsYield.hh"

ClassImp(plotResults)

using namespace std;

// ----------------------------------------------------------------------
plotResults::plotResults(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup),
									       fNoNumbers(5),
									       fCsNumbers(5),
									       fB0Numbers(5),
									       fBsmmNumbers(5),
									       fBdmmNumbers(5),
									       fHhNumbers(5),
									       fSlNumbers(5),
									       fCombNumbers(5),
									       fNpNumbers(5),
									       fBgNumbers(5),
									       fSgAndBgNumbers(5)
{
  plotClass::loadFiles(files);
  plotResults::loadFiles(files);

  changeSetup(dir, "plotResults", setup);
  init();

  fSaveSmallTree = false;

  // -- initialize cuts
  // string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  // cout << "===> Reading cuts from " << cutfile << endl;
  // readCuts(cutfile);
  // fNchan = fCuts.size();

  printCuts(cout);

  fMassLo = 4.9;
  fMassHi = 5.9;

  fNoLo = 5.10;
  fNoHi = 5.40;

  fCsLo = 5.27;
  fCsHi = 5.47;

  fBgLo = 4.9;
  fBgHi = 5.9;

  fSgLo = 5.20;
  fSgHi = 5.45;


  fChan = 0;

  TH2D *h2(0);
  TH1D *h(0);
  string mode("cnc");
  for (int imode = 0; imode < 2; ++imode) {
    if (0 == imode) {
      mode =  "cnc" + fSuffix;
    } else {
      mode = "bdt" + fSuffix;
    }
    for (int i = 0; i < fNchan; ++i) {
      h2 = new TH2D(Form("h%sAccAll%d", mode.c_str(), i), Form("h%sAccAll%d", mode.c_str(), i), 25, 0., 2.5, 25, 0., 50.);
      fhAccAll[mode].push_back(h2);
      h2 = new TH2D(Form("h%sAccPass%d", mode.c_str(), i), Form("h%sAccPass%d", mode.c_str(), i), 25, 0., 2.5, 25, 0., 50.);
      fhAccPass[mode].push_back(h2);

      h = new TH1D(Form("h%sAccPtAll%d", mode.c_str(), i), Form("h%sAccPtAll%d", mode.c_str(), i), 25, 0., 50.);
      fhAccPtAll[mode].push_back(h);
      h = new TH1D(Form("h%sAccPtPass%d", mode.c_str(), i), Form("h%sAccPtPass%d", mode.c_str(), i), 25, 0., 50.);
      fhAccPtPass[mode].push_back(h);

      h = new TH1D(Form("h%sAccEtaAll%d", mode.c_str(), i), Form("h%sAccEtaAll%d", mode.c_str(), i), 25, 0., 2.5);
      fhAccEtaAll[mode].push_back(h);
      h = new TH1D(Form("h%sAccEtaPass%d", mode.c_str(), i), Form("h%sAccEtaPass%d", mode.c_str(), i), 25, 0., 2.5);
      fhAccEtaPass[mode].push_back(h);

      h = new TH1D(Form("h%sGenAndAccNumbers%d", mode.c_str(), i), Form("h%sGenAndAccNumbers%d", mode.c_str(), i), 100, 0., 100.);
      fhGenAndAccNumbers[mode].push_back(h);

      h = new TH1D(Form("h%sMassAbsNoCuts%d", mode.c_str(), i), Form("h%sMassAbsNoCuts%d", mode.c_str(), i), 400, 2., 6.);
      fhMassAbsNoCuts[mode].push_back(h);

      h = new TH1D(Form("h%sMassNoCuts%d", mode.c_str(), i), Form("h%sMassNoCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassNoCuts[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithAnaCuts%d", mode.c_str(), i), Form("h%sMassChan%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithAnaCuts[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithMuonCuts%d", mode.c_str(), i), Form("h%sMassWithMuonCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithMuonCuts[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithTriggerCuts%d", mode.c_str(), i), Form("h%sMassWithTriggerCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithTriggerCuts[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithAllCuts%d", mode.c_str(), i), Form("h%sMassWithAllCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithAllCuts[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithAllCutsBlind%d", mode.c_str(), i), Form("h%sMassWithAllCutsBlind%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithAllCutsBlind[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithAllCutsSeagull%d", mode.c_str(), i), Form("h%sMassWithAllCutsSeagull%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithAllCutsSeagull[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithAllCutsSeagullBlind%d", mode.c_str(), i), Form("h%sMassWithAllCutsSeagullBlind%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithAllCutsSeagullBlind[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithAllCutsCowboy%d", mode.c_str(), i), Form("h%sMassWithAllCutsCowboy%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithAllCutsCowboy[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithAllCutsCowboyBlind%d", mode.c_str(), i), Form("h%sMassWithAllCutsCowboyBlind%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithAllCutsCowboyBlind[mode].push_back(h);

      // -- weighted with misid
      h = new TH1D(Form("h%sW8MassWithAllCuts%d", mode.c_str(), i), Form("h%sW8MassWithAllCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhW8MassWithAllCuts[mode].push_back(h);

      h = new TH1D(Form("h%sW8MassWithAllCutsSeagull%d", mode.c_str(), i), Form("h%sW8MassWithAllCutsSeagull%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhW8MassWithAllCutsSeagull[mode].push_back(h);

      h = new TH1D(Form("h%sW8MassWithAllCutsCowboy%d", mode.c_str(), i), Form("h%sW8MassWithAllCutsCowboy%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhW8MassWithAllCutsCowboy[mode].push_back(h);

      h = new TH1D(Form("h%sMassWithMassCuts%d", mode.c_str(), i), Form("h%sW8MassWithMassCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhMassWithMassCuts[mode].push_back(h);

      h2 = new TH2D(Form("h%sNorm%d", mode.c_str(), i), Form("h%sNorm%d", mode.c_str(), i), 100, 4.9, 5.9, MAXPS+1, -1., MAXPS);
      fhNorm[mode].push_back(h2);
      h2 = new TH2D(Form("h%sNormC%d", mode.c_str(), i), Form("h%sNormC%d", mode.c_str(), i), 200, 4.9, 5.9, MAXPS+1, -1., MAXPS);
      fhNormC[mode].push_back(h2);
    }
  }

  // -- define names and channels of all anaNumbers
  for (int i = 0; i < fNchan; ++i) {
    fNoNumbers[i].fChan   = i;
    fNoNumbers[i].fName   = "bupsik";
    fNoNumbers[i].fNameMc = "bupsikMcComb";
    fNoNumbers[i].fNameDa = "bupsikData";
    fNoNumbers[i].fMcYield.clear();
    fNoNumbers[i].fFrac.clear();

    fCsNumbers[i].fChan = i;
    fCsNumbers[i].fName = "bspsiphi";
    fCsNumbers[i].fNameMc = "bspsiphiMcComb";
    fCsNumbers[i].fNameDa = "bspsiphiData";
    fCsNumbers[i].fMcYield.clear();
    fCsNumbers[i].fFrac.clear();

    fB0Numbers[i].fChan = i;
    fB0Numbers[i].fName = "bdpsikstar";
    fB0Numbers[i].fNameMc = "bdpsikstarMc";
    fB0Numbers[i].fNameDa = "bdpsikstarData";
    fB0Numbers[i].fMcYield.clear();
    fB0Numbers[i].fFrac.clear();

    fBsmmNumbers[i].fChan = i;
    fBsmmNumbers[i].fName = "bsmm";
    fBsmmNumbers[i].fNameMc = "bsmmMcComb";
    fBsmmNumbers[i].fNameDa = "bmmData";
    fBsmmNumbers[i].fMcYield.clear();
    fBsmmNumbers[i].fFrac.clear();

    fBdmmNumbers[i].fChan = i;
    fBdmmNumbers[i].fName = "bdmm";
    fBdmmNumbers[i].fNameMc = "bdmmMcComb";
    fBdmmNumbers[i].fNameDa = "bmmData";
    fBdmmNumbers[i].fMcYield.clear();
    fBdmmNumbers[i].fFrac.clear();

    fSlNumbers[i].fChan = i;
    fSlNumbers[i].fName = "sl";
    fSlNumbers[i].fNameMc = "nada";
    fSlNumbers[i].fNameDa = "nada";

    fHhNumbers[i].fChan = i;
    fHhNumbers[i].fName = "hh";
    fHhNumbers[i].fNameMc = "nada";
    fHhNumbers[i].fNameDa = "nada";
    fHhNumbers[i].fMcYield.clear();

    fNpNumbers[i].fChan = i;
    fNpNumbers[i].fName = "np";
    fNpNumbers[i].fNameMc = "nada";
    fNpNumbers[i].fNameDa = "nada";

    fBgNumbers[i].fChan = i;
    fBgNumbers[i].fName = "bg";
    fBgNumbers[i].fNameMc = "nada";
    fBgNumbers[i].fNameDa = "nada";

    fSgAndBgNumbers[i].fChan = i;
    fSgAndBgNumbers[i].fName = "sgandbg";
    fSgAndBgNumbers[i].fNameMc = "nada";
    fSgAndBgNumbers[i].fNameDa = "nada";


    fHhNumbers[i].fMcYield.clear();
    fSlNumbers[i].fMcYield.clear();
    fNpNumbers[i].fMcYield.clear();
    fBgNumbers[i].fMcYield.clear();
    fBgNumbers[i].fObsYield.clear();
    fSgAndBgNumbers[i].fMcYield.clear();
    fSgAndBgNumbers[i].fObsYield.clear();
    for (int j = 0; j < NWIN; ++j) {
      number aaa;
      fHhNumbers[i].fMcYield.push_back(aaa);
      fSlNumbers[i].fMcYield.push_back(aaa);

      fNpNumbers[i].fMcYield.push_back(aaa);

      fBgNumbers[i].fMcYield.push_back(aaa);
      fBgNumbers[i].fObsYield.push_back(aaa);

      fSgAndBgNumbers[i].fMcYield.push_back(aaa);
      fSgAndBgNumbers[i].fObsYield.push_back(aaa);

      fNoNumbers[i].fFrac.push_back(aaa);
      fNoNumbers[i].fMcYield.push_back(aaa);

      fCsNumbers[i].fFrac.push_back(aaa);
      fCsNumbers[i].fMcYield.push_back(aaa);

      fB0Numbers[i].fFrac.push_back(aaa);
      fB0Numbers[i].fMcYield.push_back(aaa);

      fBdmmNumbers[i].fFrac.push_back(aaa);
      fBdmmNumbers[i].fMcYield.push_back(aaa);

      fBsmmNumbers[i].fFrac.push_back(aaa);
      fBsmmNumbers[i].fMcYield.push_back(aaa);
    }

    fCombNumbers[i].fChan = i;
    fCombNumbers[i].fName = "comb";
    fCombNumbers[i].fNameMc = "nada";
    fCombNumbers[i].fNameDa = "bmmData";
    fCombNumbers[i].fObsYield.clear();
    fCombNumbers[i].fFitYield.clear();
    for (int j = 0; j < NWIN; ++j) {
      number aaa;
      fCombNumbers[i].fObsYield.push_back(aaa);
      fCombNumbers[i].fFitYield.push_back(aaa);
    }

  }

  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    if (skipThisBg(it->first)) continue;
    vector<anaNumbers*> va;
    for (int i = 0; i < fNchan; ++i) {
      anaNumbers *a = new anaNumbers(it->first, i);
      a->fNameMc = it->first;
      number aaa;
      a->fName = it->first;
      a->fFrac.clear();
      a->fMcYield.clear();
      a->fObsYield.clear();
      a->fFitYield.clear();
      for (int j = 0; j < 5; ++j) {
	a->fFrac.push_back(aaa);
	a->fMcYield.push_back(aaa);
	a->fObsYield.push_back(aaa);
	a->fFitYield.push_back(aaa);
      }
      replaceAll(a->fName, "McOff", "");
      replaceAll(a->fName, "McComb", "");
      replaceAll(a->fName, "Mc", "");
      va.push_back(a);
    }
    fRareNumbers.insert(make_pair(it->first, va));
  }


  // -- define systematics, per channel!
  double fsfu[] = {0.078, 0.078, 0.078, 0.078};
  fSystematics["fsfu"] = vector<double>(fsfu, fsfu + sizeof(fsfu)/sizeof(fsfu[0]));
  double bfbupsik[] = {0.031, 0.031, 0.031, 0.031};
  fSystematics["bfbupsik"] = vector<double>(bfbupsik, bfbupsik + sizeof(bfbupsik)/sizeof(bfbupsik[0]));
  double bfbspsiphi[] = {0.076, 0.076, 0.076, 0.076};
  fSystematics["bfbspsiphi"] = vector<double>(bfbspsiphi, bfbspsiphi + sizeof(bfbspsiphi)/sizeof(bfbspsiphi[0]));

  double pxy[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["pxy"] = vector<double>(pxy, pxy + sizeof(pxy)/sizeof(pxy[0]));

  double acceptance[] = {0.035, 0.035, 0.050, 0.050};
  fSystematics["acceptance"] = vector<double>(acceptance, acceptance + sizeof(acceptance)/sizeof(acceptance[0]));
  double effanabupsik[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effanabupsik"] = vector<double>(effanabupsik, effanabupsik + sizeof(effanabupsik)/sizeof(effanabupsik[0]));
  double effanabspsiphi[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effanabspsiphi"] = vector<double>(effanabspsiphi, effanabspsiphi + sizeof(effanabspsiphi)/sizeof(effanabspsiphi[0]));
  double effanabsmm[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effanabsmm"] = vector<double>(effanabsmm, effanabsmm + sizeof(effanabsmm)/sizeof(effanabsmm[0]));
  double effanabdmm[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effanabdmm"] = vector<double>(effanabdmm, effanabdmm + sizeof(effanabdmm)/sizeof(effanabdmm[0]));


  double effcandbupsik[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effcandbupsik"] = vector<double>(effcandbupsik, effcandbupsik + sizeof(effcandbupsik)/sizeof(effcandbupsik[0]));
  double effcandbspsiphi[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effcandbspsiphi"] = vector<double>(effcandbspsiphi, effcandbspsiphi + sizeof(effcandbspsiphi)/sizeof(effcandbspsiphi[0]));
  double effcandbsmm[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effcandbsmm"] = vector<double>(effcandbsmm, effcandbsmm + sizeof(effcandbsmm)/sizeof(effcandbsmm[0]));
  double effcandbdmm[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effcandbdmm"] = vector<double>(effcandbdmm, effcandbdmm + sizeof(effcandbdmm)/sizeof(effcandbdmm[0]));

  double effmuid[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["effmuid"] = vector<double>(effmuid, effmuid + sizeof(effmuid)/sizeof(effmuid[0]));
  double efftrig[] = {0.05, 0.05, 0.050, 0.050};
  fSystematics["efftrig"] = vector<double>(efftrig, efftrig + sizeof(efftrig)/sizeof(efftrig[0]));
  double fakemuid[] = {0.2, 0.2, 0.2, 0.2};
  fSystematics["fakemuid"] = vector<double>(effmuid, effmuid + sizeof(effmuid)/sizeof(effmuid[0]));

  double bunorm[] = {0.03, 0.04, 0.05, 0.10};
  fSystematics["normbupsik"] = vector<double>(bunorm, bunorm + sizeof(bunorm)/sizeof(bunorm[0]));
  double bsnorm[] = {0.04, 0.04, 0.05, 0.06};
  fSystematics["normbspsiphi"] = vector<double>(bsnorm, bsnorm + sizeof(bsnorm)/sizeof(bsnorm[0]));


  cout << "-----------------------------------------------------" << endl;
  cout << "Systematic contribution: ";
  for (int i = 0; i < fNchan; ++i) {
    cout << " chan " << i;
  }
  cout << endl;
  cout << "-----------------------------------------------------" << endl;
  for (std::map<std::string, vector<double> >::iterator it = fSystematics.begin(); it != fSystematics.end(); ++it) {
    cout << Form("%23s: ", it->first.c_str());
    for (unsigned int i = 0; i < it->second.size(); ++i) {
      cout << Form("  %4.3f", it->second.at(i));
    }
    cout << endl;
  }
  cout << "-----------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
plotResults::~plotResults() {
  cout << "plotResults destructor" << endl;
  resetHistograms(true);
  cout << "done with histogram deletion" << endl;

}


// ----------------------------------------------------------------------
void plotResults::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotResults::makeAll(string what) {
  if (what == "all" || what == "dumpdatasets" || what == "ana" || what == "genvalidation") {
    dumpDatasets();
  }

  if (what == "all" || what == "genvalidation") {
    if (2016 == fYear) {
      genSummary("bdmmMcCombAcc", "candAnaMuMu");
      genSummary("bdmmMcOff", "candAnaMuMu");
      genSummary("bdmmMc", "candAnaMuMu");
      genSummary("bsmmMcCombAcc", "candAnaMuMu");
      genSummary("bsmmMcOff", "candAnaMuMu");
      genSummary("bsmmMc", "candAnaMuMu");
      genSummary("bsmm80Mc", "candAnaMuMu");
      genSummary("bsmm75Mc", "candAnaMuMu");
      genSummary("bsmm70Mc", "candAnaMuMu");
      genSummary("bsmm69Mc", "candAnaMuMu");
      genSummary("bsmm68Mc", "candAnaMuMu");
      genSummary("bsmm67Mc", "candAnaMuMu");
      genSummary("bsmm66Mc", "candAnaMuMu");
      genSummary("bsmm65Mc", "candAnaMuMu");
      genSummary("bsmm60Mc", "candAnaMuMu");
      genSummary("bsmm55Mc", "candAnaMuMu");
      genSummary("bsmm50Mc", "candAnaMuMu");
      genSummary("bsmm45Mc", "candAnaMuMu");
      genSummary("bsmm40Mc", "candAnaMuMu");
      genSummary("bsmm35Mc", "candAnaMuMu");

      genSummary("bupsikMcCombAcc", "candAnaBu2JpsiK");
      genSummary("bupsikMcOff", "candAnaBu2JpsiK");
      genSummary("bupsikMc", "candAnaBu2JpsiK");
      genSummary("bupsikMcComb", "candAnaBu2JpsiK");
      genSummary("bspsiphiMcCombAcc", "candAnaBs2JpsiPhi");
      genSummary("bspsiphiMcOff", "candAnaBs2JpsiPhi");
      genSummary("bspsiphiMc", "candAnaBs2JpsiPhi");
      genSummary("bspsiphiMcComb", "candAnaBs2JpsiPhi");
      genSummary("bspsifMcComb", "candAnaBs2Jpsif0");
      genSummary("bspsifMcCombAcc", "candAnaBs2Jpsif0");
      genSummary("bdpsikstarMcCombAcc", "candAnaBd2JpsiKstar");
      genSummary("bdpsikstarMcComb", "candAnaBd2JpsiKstar");


      // -- loop over all (effective) two-body backgrounds
      for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
	if (string::npos == it->first.find("Bg")) continue;
	cout << "===> Create genSummary for " << it->first << endl;
	genSummary(it->first, "candAnaMuMu");
      }
    } else if (2012 == fYear) {
      genSummary("bupsikMcComb", "candAnaBu2JpsiK");
      genSummary("bupsikMcCombAcc", "candAnaBu2JpsiK");

      genSummary("bspsiphiMcComb", "candAnaBs2JpsiPhi");
      genSummary("bspsiphiMcCombAcc", "candAnaBs2JpsiPhi");

      genSummary("bdpsikstarMcComb", "candAnaBd2JpsiKstar");
      genSummary("bdpsikstarMcCombAcc", "candAnaBd2JpsiKstar");

      genSummary("bspsifMcComb", "candAnaBs2Jpsif0");
      genSummary("bspsifMcCombAcc", "candAnaBs2Jpsif0");

      genSummary("bdmmMcComb", "candAnaMuMu");
      genSummary("bdmmMcCombAcc", "candAnaMuMu");

      genSummary("bsmmMcComb", "candAnaMuMu");
      genSummary("bsmmMcCombAcc", "candAnaMuMu");

    }

  }

  // -- this will recreate fHistFile!
  if (what == "all" || string::npos != what.find("fill")) {
    fillAndSaveHistograms();
  }

  if (what == "dbx") {
    fillAndSaveHistograms(0, 1e5);
  }

  if (what == "all" || string::npos != what.find("ana")) {
    dumpDatasets();
    fHistWithAllCuts = "hMassWithAllCuts";
    calculateNumbers("cnc" + fSuffix);
    calculateNumbers("bdt" + fSuffix);
  }

  if (what == "all" || string::npos != what.find("cnc")) {
    fHistWithAllCuts = "hMassWithAllCuts";
    calculateNumbers("cnc" + fSuffix);
  }

  if (what == "all" || string::npos != what.find("bdt")) {
    fHistWithAllCuts = "hMassWithAllCuts";
    calculateNumbers("bdt" + fSuffix);
  }



}



// ----------------------------------------------------------------------
void plotResults::bookHist(string dsname) {


}


// ----------------------------------------------------------------------
void plotResults::dumpDatasets() {

  std::ofstream TEX;
  string dsname = Form("%s/%d-datasets.tex", fDirectory.c_str(), fYear);
  system(Form("/bin/rm -f %s", dsname.c_str()));
  TEX.open(dsname.c_str(), ios::app);

  TEX << "% ----------------------------------------------------------------------" << endl;
  TEX << formatTex(fBfPsiMuMu, Form("%s:BfPsiMuMu:val", fSuffix.c_str()), 5) << endl;
  TEX << formatTex(fBfPsiMuMuE, Form("%s:BfPsiMuMu:err", fSuffix.c_str()), 5) << endl;

  TEX << formatTex(fBfPhiKpKm, Form("%s:BfPhiKpKm:val", fSuffix.c_str()), 5) << endl;
  TEX << formatTex(fBfPhiKpKmE, Form("%s:BfPhiKpKm:err", fSuffix.c_str()), 5) << endl;
  TEX << formatTex(fBfKstarKpPim, Form("%s:BfKstarKpPim:val", fSuffix.c_str()), 5) << endl;
  TEX << formatTex(fBfKstarKpPimE, Form("%s:BfKstarKpPim:err", fSuffix.c_str()), 5) << endl;

  TEX << formatTexErrSci(fCrossSection, 0., Form("%s:PythiaCrossSection:val", fSuffix.c_str()), -1) << endl;

  TEX << "% ----------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    TH1D *h1 = it->second->getHist("monEvents", false);
    int nEvtFile = static_cast<int>(h1->GetBinContent(1));
    int nCands   = static_cast<int>(h1->GetBinContent(2));
    double epscand  = static_cast<double>(nCands)/nEvtFile;
    double epscandE = dEff(nCands, nEvtFile);
    TEX << Form("\\vdef{%s:%s:name} {%s}", fSuffix.c_str(), it->first.c_str(), it->first.c_str()) << endl;
    TEX << Form("\\vdef{%s:%s:decay} {%s}", fSuffix.c_str(), it->first.c_str(), it->second->fLatexName.c_str()) << endl;
    TEX << formatTex(nEvtFile, Form("%s:%s:nEvtFile", fSuffix.c_str(), it->first.c_str()), 0) << endl;
    TEX << formatTex(nCands, Form("%s:%s:nCands", fSuffix.c_str(), it->first.c_str()), 0) << endl;
    TEX << formatTex(epscand, Form("%s:%s:epsCand:val", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    TEX << formatTex(epscandE, Form("%s:%s:epsCand:err", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    TEX << formatTex(it->second->fFilterEff, Form("%s:%s:filterEff:val", fSuffix.c_str(), it->first.c_str()), 6) << endl;
    TEX << formatTex(0., Form("%s:%s:filterEff:err", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    if (it->second->fBf > 0.) {
      TEX << formatTexErrSci(it->second->fBf, it->second->fBfE, Form("%s:%s:bf", fSuffix.c_str(), it->first.c_str()), 2) << endl;
      double eqLumi(0.);
      if (it->second->fFilterEff > 0) {
	eqLumi = nEvtFile/fCrossSection/it->second->fBf/it->second->fFilterEff;
	TEX << formatTex(eqLumi, Form("%s:%s:eqLumi:val", fSuffix.c_str(), it->first.c_str()), 1) << endl;
	cout << it->first << ": eqLumi = " << eqLumi << " filterEff: " << it->second->fFilterEff << endl;
      }
    }
    if (it->second->fLumi > 0.) {
      TEX << formatTex(it->second->fLumi, Form("%s:%s:lumi", fSuffix.c_str(), it->first.c_str()), 1) << endl;
    }
  }
  TEX.close();
}


// ----------------------------------------------------------------------
void plotResults::genSummary(std::string dsname, std::string dir) {
  TH1D *hpt(0), *heta(0), *tpt(0), *teta(0),
    *hm1eta(0), *hm2eta(0), *tm1eta(0), *tm2eta (0), *hketa(0),
    *hm1pt(0), *hm2pt(0),
    *tm1pt(0), *tm2pt(0), *hkpt(0),
    *htau(0);

  TH1D *h1 = (TH1D*)gDirectory->Get("pt");
  if (0 == h1) {
    hpt    = new TH1D("pt", "pt", 50, 0., 50.0);
    heta   = new TH1D("eta", "eta", 40, -4., 4.0);
    tpt    = new TH1D("tpt", "pt (HLT)", 50, 0., 50.0); setFilledHist(tpt, kBlue, kYellow, 1000);
    teta   = new TH1D("teta", "eta (HLT)", 40, -4., 4.0); setFilledHist(teta, kBlue, kYellow, 1000);
    hm1eta = new TH1D("m1eta", "m1 eta", 40, -4., 4.0);
    hm2eta = new TH1D("m2eta", "m2 eta", 40, -4., 4.0);
    tm1eta = new TH1D("tm1eta", "m1 eta (HLT)", 40, -4., 4.0); setFilledHist(tm1eta, kBlue, kYellow, 1000);
    tm2eta = new TH1D("tm2eta", "m2 eta (HLT)", 40, -4., 4.0); setFilledHist(tm2eta, kBlue, kYellow, 1000);
    hketa  = new TH1D("keta", "kaon eta", 40, -4., 4.0);
    hm1pt  = new TH1D("m1pt", "m1 pt", 100, 0., 10.0);
    hm2pt  = new TH1D("m2pt", "m2 pt", 100, 0., 10.0);
    tm1pt  = new TH1D("tm1pt", "m1 pt (HLT)", 100, 0., 10.0); setFilledHist(tm1pt, kBlue, kYellow, 1000);
    tm2pt  = new TH1D("tm2pt", "m2 pt (HLT)", 100, 0., 10.0); setFilledHist(tm2pt, kBlue, kYellow, 1000);
    hkpt   = new TH1D("kpt", "kaon pt", 50, 0., 10.0);
    htau   = new TH1D("tau", "tau", 100, 0., 15.e-12);
    TTree *T = getTree(dsname, dir, "effTree");
  } else {
    hpt    = (TH1D*)gDirectory->Get("pt");
    heta   = (TH1D*)gDirectory->Get("eta");
    tpt    = (TH1D*)gDirectory->Get("tpt");
    teta   = (TH1D*)gDirectory->Get("teta");
    hm1eta = (TH1D*)gDirectory->Get("m1eta");
    hm2eta = (TH1D*)gDirectory->Get("m2eta");
    tm1eta = (TH1D*)gDirectory->Get("tm1eta");
    tm2eta = (TH1D*)gDirectory->Get("tm2eta");
    hketa  = (TH1D*)gDirectory->Get("keta");
    hm1pt  = (TH1D*)gDirectory->Get("m1pt");
    hm2pt  = (TH1D*)gDirectory->Get("m2pt");
    tm1pt  = (TH1D*)gDirectory->Get("tm1pt");
    tm2pt  = (TH1D*)gDirectory->Get("tm2pt");
    hkpt   = (TH1D*)gDirectory->Get("kpt");
    htau   = (TH1D*)gDirectory->Get("tau");
  }
  TTree *T = getTree(dsname, dir, "effTree");

  T->Draw("gtau>>tau");
  T->Draw("gpt>>pt");
  T->Draw("geta>>eta");

  T->Draw("g1pt>>m1pt");
  T->Draw("g2pt>>m2pt");
  T->Draw("g1eta>>m1eta");
  T->Draw("g2eta>>m2eta");

  T->Draw("gpt>>tpt", "hlt");
  T->Draw("g1pt>>tm1pt", "hlt");
  T->Draw("g2pt>>tm2pt", "hlt");

  T->Draw("geta>>teta", "hlt");
  T->Draw("g1eta>>tm1eta", "hlt");
  T->Draw("g2eta>>tm2eta", "hlt");


  bool addKaon(false);
  bool addHLT(true);
  if (string::npos != dsname.find("Bg")) {
    addHLT = false;
  }

  if (string::npos != dsname.find("bupsik")) {
    T->Draw("g3eta>>keta");
    T->Draw("g3pt>>kpt");
    addKaon = true;
  }
  if (string::npos != dsname.find("bspsiphi")) {
    T->Draw("g3eta>>keta");
    T->Draw("g4eta>>keta");
    T->Draw("g3pt>>kpt");
    T->Draw("g4pt>>kpt");
    addKaon = true;
  }

  tl->SetTextSize(0.05);
  makeCanvas(1);
  int ncol(4);
  if (addKaon) {
    c1->Divide(5,2);
    ncol = 5;
  } else {
    c1->Divide(4,2);
  }

  c1->cd(1);
  setTitles(hpt, "p_{T} [GeV]", "entries/bin");
  hpt->Draw();
  if (addHLT) {
    tpt->Draw("same");
    hpt->Draw("sameaxis");
  }

  c1->cd(2);
  hm1pt->Draw();
  setTitles(hm1pt, "p_{T} [GeV]", "entries/bin");
  if (addHLT) {
    tm1pt->Draw("same");
    hm1pt->Draw("sameaxis");
  }

  c1->cd(3);
  hm2pt->Draw();
  setTitles(hm2pt, "p_{T} [GeV]", "entries/bin");
  if (addHLT) {
    tm2pt->Draw("same");
    hm2pt->Draw("sameaxis");
  }

  if (addKaon) {
    c1->cd(ncol-1);
    setTitles(hkpt, "p_{T} [GeV]", "entries/bin");
    hkpt->Draw();
  }

  c1->cd(ncol);
  gPad->SetLogy(1);
  htau->Fit("expo", "l");
  TF1 *f = (TF1*)htau->GetFunction("expo");
  double t    = -1./f->GetParameter(1);
  double tE   = -t*f->GetParError(1)/f->GetParameter(1);
  t  *= 1.e12;
  tE *= 1.e12;

  c1->cd(ncol+1);
  setTitles(heta, "#eta", "entries/bin");
  heta->Draw();
  if (addHLT) {
    teta->Draw("same");
    heta->Draw("axissame");
  }
  tl->DrawLatexNDC(0.55, 0.3, "B");

  c1->cd(ncol+2);
  setTitles(hm1eta, "#eta", "entries/bin");
  hm1eta->Draw();
  if (addHLT) {
    tm1eta->Draw("same");
    hm1eta->Draw("sameaxis");
  }
  tl->DrawLatexNDC(0.40, 0.3, "leading muon");

  c1->cd(ncol+3);
  setTitles(hm2eta, "#eta", "entries/bin");
  hm2eta->Draw();
  if (addHLT) {
    tm2eta->Draw("same");
    hm2eta->Draw("sameaxis");
  }
  tl->DrawLatexNDC(0.35, 0.3, "subleading muon");

  if (addKaon) {
    c1->cd(2*ncol-1);
    setTitles(hketa, "#eta", "entries/bin");
    hketa->Draw();
    tl->DrawLatexNDC(0.35, 0.3, "kaon(s)");
  }

  c1->cd(2*ncol);
  tl->SetTextSize(0.07);
  if (addHLT) {
    tl->DrawLatexNDC(0.2, 0.9, Form("%s (HLT)", dsname.c_str()));
  } else {
    tl->DrawLatexNDC(0.2, 0.9, Form("%s ", dsname.c_str()));
  }
  tl->DrawLatexNDC(0.2, 0.8, Form("Events: %d", T->GetEntries()));

  tl->DrawLatexNDC(0.2, 0.35, Form("#tau"));
  tl->DrawLatexNDC(0.25, 0.35, Form("= (%4.3f #pm %5.3f)ps", t, tE));

  tl->DrawLatexNDC(0.2, 0.25, Form("#varepsilon"));
  tl->DrawLatexNDC(0.25, 0.25, Form("= %4.3f ", teta->GetSumOfWeights()/heta->GetSumOfWeights()));

  c1->SaveAs(Form("%s/genSummary-%d-%s.pdf", fDirectory.c_str(), fYear, dsname.c_str()));

}



// ----------------------------------------------------------------------
void plotResults::fillAndSaveHistograms(int start, int nevents) {
  // -- dump histograms
  fHistFile = TFile::Open(fHistFileName.c_str(), "RECREATE");
  cout << " opened, running on " << nevents << " entries" << endl;

  TTree *t(0);

  fSaveSmallTree = true;

  if (0) {
    // -- for debugging
    resetHistograms();
    setup("bdmmMc");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, -1, start);
    saveHistograms(fSetup);
    fHistFile->Close();
    return;
  }

  // -- rare backgrounds
  if (1) {
    resetHistograms();
    rareBgHists("nada", nevents);
  }

  // -- normalization modes
  if (1) {
    resetHistograms();
    setup("bupsikData");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents, start);
    fSetup = "bupsikData";
    saveHistograms(fSetup);

    resetHistograms();
    setup("bupsikMcComb");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents, start);
    cout << "done with loopovertree" << endl;
    otherNumbers(fSetup);
    saveHistograms(fSetup);

    resetHistograms();
    setup("bmmData");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents, start);
    fSetup = "bmmData";
    saveHistograms(fSetup);

    resetHistograms();
    setup("bdmmMcComb");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents, start);
    otherNumbers(fSetup);
    saveHistograms(fSetup);

    resetHistograms();
    setup("bsmmMcComb");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents, start);
    otherNumbers(fSetup);
    saveHistograms(fSetup);

    if (1) {
    resetHistograms();
    setup("bspsiphiData");
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents, start);
    fSetup = "bspsiphiData";
    saveHistograms(fSetup);

    resetHistograms();
    fSetup = "bspsiphiMcComb";
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    loopOverTree(t, 1, nevents, start);
    otherNumbers(fSetup);
    saveHistograms(fSetup);
    }
  }

  fHistFile->Close();

  fSaveSmallTree = false;
}


// ----------------------------------------------------------------------
void plotResults::initNumbers(anaNumbers &a) {
  a.clear();
}


// ----------------------------------------------------------------------
void plotResults::rareBgHists(string smode, int nevents) {
  int nloop(0);
  // -- loop over all (effective) two-body backgrounds
  int start(0);
  TTree *t(0);
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    if (skipThisBg(it->first)) continue;

    resetHistograms();
    setup(it->first);
    t = getTree(fSetup, fTreeDir);
    setupTree(t, fSetup);
    cout << "==============================================================" << endl;
    cout << "==> rareBgHists for " << it->first << " and setup = " << fSetup << endl;
    cout << "==============================================================" << endl;
    //    loopOverTree(t, 1, 100000, start);
    loopOverTree(t, 1, nevents, start);
    otherNumbers(fSetup);
    saveHistograms(fSetup);
    ++nloop;
    //    if (nloop > 3) break;
  }
}


// ----------------------------------------------------------------------
void plotResults::otherNumbers(string smode) {
  int ibin(0);
  string accname = smode + "Acc";
  //  replaceAll(accname, "Comb", "Off");
  // -- just in case you are running a cross-check on the acceptance sample:
  replaceAll(accname, "AccAcc", "Acc");
  // -- For Bg samples, the situation is more complex:
  if (string::npos != smode.find("Bg")) {
    accname = rareAccName(smode);
  }

  anaNumbers* aa[fNchan];
  if (string::npos != fSetup.find("bupsik"))  {
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = &fNoNumbers[ichan];
  }
  if (string::npos != fSetup.find("bspsiphi"))  {
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = &fCsNumbers[ichan];
  }

  if (string::npos != fSetup.find("bdmm"))  {
    cout << "setting aa to fBdmmNumbers" << endl;
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = &fBdmmNumbers[ichan];
  }

  if (string::npos != fSetup.find("bsmm"))  {
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = &fBsmmNumbers[ichan];
  }

  if (string::npos != fSetup.find("Bg"))  {
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = fRareNumbers[fSetup][ichan];
  }

  double effGenSel(0.), effGenSelE(0.);

  if (string::npos != fSetup.find("Mc"))  {
    effGenSel = fDS[smode]->fFilterEff/fDS[accname]->fFilterEff;
    effGenSelE = dRatio(fDS[smode]->fFilterEff, fDS[smode]->fFilterEffE, fDS[accname]->fFilterEff, fDS[accname]->fFilterEffE);
  }

  cout << "smode = " << smode << " fSetup = " << fSetup
       << " accname: " << accname << " directory: " << fTreeDir
       << " numbers: " << aa[0]->fName << " chan = " << aa[0]->fChan
       << " effGenSel = " << effGenSel << " +/- " << effGenSelE
       << endl;

  vector<string> modifier;
  modifier.push_back("cnc" + fSuffix);
  modifier.push_back("bdt" + fSuffix);

  for (unsigned int i = 0; i < fNchan; ++i) {
    fChan = i;
    // -- fill the numbers into anaNumbers
    getAccAndEffFromEffTree(accname, *aa[i], *fCuts[i], -1);
    // -- for tests on the acc samples, overwrite fEffGenSel
    if (string::npos != smode.find("Acc")) {
      aa[i]->fEffGenSel.val = 1.;
      aa[i]->fEffGenSel.estat = 0.;
    }
    fDS[fSetup]->cd(fTreeDir.c_str());
    double effFilter  = fDS[fSetup]->fFilterEff;
    double effFilterE = fDS[fSetup]->fFilterEffE;
    if (effFilter < 1e-6) effFilter = 1.0;
    double genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();

    for (unsigned int im = 0; im < modifier.size(); ++im) {
      // -- fill numbers from anaNumbers into histogram
      ibin = 1;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, effFilter);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "effFilter");
      ibin = 2;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, effFilterE);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "effFilterE");
      ibin = 3;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, effGenSel);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "effGenSel");
      ibin = 4;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, effGenSelE);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "effGenSelE");

      ibin = 11;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, genFileYield);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "genFileYield");
      ibin = 12;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, genFileYield/effGenSel);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "genYield");
      ibin = 13;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, dRatio(genFileYield, TMath::Sqrt(genFileYield), effGenSel, effGenSelE));
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "genYieldE");

      ibin = 20;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, aa[i]->fAcc.val);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "acc");
      ibin = 21;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, aa[i]->fAcc.estat);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "accE");

      ibin = 22;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, aa[i]->fAccGenFileYield.val);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "accGenFileYield");
      ibin = 23;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, aa[i]->fAccGenYield.val);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "accGenYield");

      ibin = 24;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, aa[i]->fAccRecoYield.val);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "accRecoYield");
      ibin = 25;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, aa[i]->fAccCandYield.val);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "accCandYield");


      ibin = 30;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, aa[i]->fEffCand.val);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "effCand");
      ibin = 31;
      fhGenAndAccNumbers[modifier[im]][i]->SetBinContent(ibin, aa[i]->fEffCand.estat);
      fhGenAndAccNumbers[modifier[im]][i]->GetXaxis()->SetBinLabel(ibin, "effCandE");
    }
  }
}

// ----------------------------------------------------------------------
void plotResults::getAccAndEffFromEffTree(string ds, anaNumbers &a, cuts &b, int proc) {

  TTree *t  = getTree(ds, fTreeDir, "effTree");
  if (!t) {
    cout << "plotResults::getAccAndEffFromEffTree(" << a.fName << "): no tree `effTree' found ????????????"
         << endl;
    return;
  } else {
    a.fEffFilter.val = fDS[ds]->fFilterEff;
    cout << "plotResults::getAccAndEffFromEffTree(" << a.fName << ")" << endl
         << " get acceptance from file " << fDS[ds]->fFullName << " and dir = " << Form("%s/effTree", ds.c_str())
	 << endl
         << " with a.fEffFilter = " << a.fEffFilter.val
         << " with a.fEffGenSel = " << a.fEffGenSel.val
         << endl;
  }

  bool sg(false), no(false), cs(false);

  int   bprocid, bidx, bchan;
  bool  bhlt1;
  float bg1pt, bg2pt, bg1eta, bg2eta;
  float bm1pt, bm1eta, bm2pt, bm2eta;
  bool  bm1gt, bm2gt;
  float bm;

  float bg3pt, bg4pt, bg3eta, bg4eta;
  float bk1pt, bk2pt, bk1eta, bk2eta;
  bool  bk1gt, bk2gt;

  t->SetBranchAddress("hlt1", &bhlt1);
  t->SetBranchAddress("procid",&bprocid);
  t->SetBranchAddress("bidx",&bidx);
  t->SetBranchAddress("chan",&bchan);

  t->SetBranchAddress("g1pt",&bg1pt);
  t->SetBranchAddress("g2pt",&bg2pt);
  t->SetBranchAddress("g1eta",&bg1eta);
  t->SetBranchAddress("g2eta",&bg2eta);

  t->SetBranchAddress("m1pt",&bm1pt);
  t->SetBranchAddress("m2pt",&bm2pt);
  t->SetBranchAddress("m1eta",&bm1eta);
  t->SetBranchAddress("m2eta",&bm2eta);

  t->SetBranchAddress("m1gt",&bm1gt);
  t->SetBranchAddress("m2gt",&bm2gt);

  if (string::npos != a.fName.find("bdmm")) {
    cout << "anaBmm::getAccAndEffFromEffTree(" << a.fName << "): SIGNAL " << endl;
    sg = true;
  }
  if (string::npos != a.fName.find("bsmm")) {
    cout << "anaBmm::getAccAndEffFromEffTree(" << a.fName << "): SIGNAL " << endl;
    sg = true;
  }
  if (string::npos != a.fName.find("Bg")) {
    cout << "anaBmm::getAccAndEffFromEffTree(" << a.fName << "): RARE BACKGROUND " << endl;
    sg = true;
  }
  if (string::npos != a.fName.find("bupsik")) {
    cout << "anaBmm::getAccAndEffFromEffTree(" << a.fName << "): NORMALIZATION " << endl;
    no = true;
    t->SetBranchAddress("g3pt", &bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
    t->SetBranchAddress("kpt", &bk1pt);
    t->SetBranchAddress("keta",&bk1eta);
    t->SetBranchAddress("kgt", &bk1gt);
  }
  if (string::npos != a.fName.find("bspsiphi")) {
    cout << "anaBmm::getAccAndEffFromEffTree(" << a.fName << "): CONTROL SAMPLE " << endl;
    cs = true;
    t->SetBranchAddress("g3pt", &bg3pt);
    t->SetBranchAddress("g3eta",&bg3eta);
    t->SetBranchAddress("k1pt", &bk1pt);
    t->SetBranchAddress("k1eta",&bk1eta);
    t->SetBranchAddress("k1gt", &bk1gt);

    t->SetBranchAddress("g4pt", &bg4pt);
    t->SetBranchAddress("g4eta",&bg4eta);
    t->SetBranchAddress("k2pt", &bk2pt);
    t->SetBranchAddress("k2eta",&bk2eta);
    t->SetBranchAddress("k2gt", &bk2gt);
  }
  if (string::npos != a.fName.find("bdpsikstar")) {
    cout << "FFFFFFIIIIIIIIIIIXXXXXXXXXMMMMMMMMMEEEEEEEE!!!!!!!" << endl;
    return;
  }
  t->SetBranchAddress("m",&bm);


  int nentries = Int_t(t->GetEntries());
  int ngen(0), nreco(0), ncand(0);
  cout << "channel = " << a.fChan << endl;
  for (int jentry = 0; jentry < nentries; jentry++) {
    t->GetEntry(jentry);
    if (bidx < 0) continue;
    ++ngen;
    if (proc > 0 && bprocid != proc) continue;
    if (sg) {
      // -- Signal
      if (TMath::Abs(bg1eta) < fAccEtaGen && TMath::Abs(bg2eta) < fAccEtaGen
          && bg1pt > fAccPt && bg2pt > fAccPt
          && bm1pt > fAccPt && bm2pt > fAccPt
          && TMath::Abs(bm1eta) < fAccEtaRec && TMath::Abs(bm2eta) < fAccEtaRec
          && bm1gt && bm2gt
          ) {
        if (bchan == a.fChan) {
          ++nreco; // for acceptance
          if (bm > 0) {
            ++ncand; // for cand efficiency
          }
        }
      }
    } else if (no) {
      // -- Normalization
      if (TMath::Abs(bg1eta) < fAccEtaGen && TMath::Abs(bg2eta) < fAccEtaGen  && TMath::Abs(bg3eta) < fAccEtaGen
          && bg1pt > fAccPt && bg2pt > fAccPt && bg3pt > 0.4
          && bm1pt > fAccPt && bm2pt > fAccPt && bk1pt > 0.5
          && TMath::Abs(bm1eta) < fAccEtaRec && TMath::Abs(bm2eta) < fAccEtaRec && TMath::Abs(bk1eta) < fAccEtaRec
          && bm1gt && bm2gt && bk1gt
          ) {
        if (bchan == a.fChan) {
          ++nreco; // for acceptance
          if (bm > 0) {
            ++ncand; // for cand efficiency
          }
        }
      }
    } else if (cs) {
      // -- control sample
      if (TMath::Abs(bg1eta) < fAccEtaGen && TMath::Abs(bg2eta) < fAccEtaGen && TMath::Abs(bg3eta) < fAccEtaGen && TMath::Abs(bg4eta) < fAccEtaGen
          && bg1pt > fAccPt && bg2pt > fAccPt && bg3pt > 0.4 && bg4pt > 0.4
          && bm1pt > fAccPt && bm2pt > fAccPt && bk1pt > 0.5 && bk2pt > 0.5
          && TMath::Abs(bm1eta) < fAccEtaRec && TMath::Abs(bm2eta) < fAccEtaRec && TMath::Abs(bk1eta) < fAccEtaRec && TMath::Abs(bk2eta) < fAccEtaRec
          && bm1gt && bm2gt && bk1gt && bk2gt
          ) {
        if (bchan == a.fChan) {
          ++nreco; // for acceptance
          if (bm > 0) {
            ++ncand; // for cand efficiency
          }
        }
      }
    }
  }

  a.fAccGenFileYield.val = ngen;
  a.fAccGenYield.val     = ngen;
  a.fAccRecoYield.val    = nreco;
  a.fAccCandYield.val    = ncand;

  if (a.fAccGenYield.val > 0) {
    a.fAcc.val   = a.fAccRecoYield.val/a.fAccGenYield.val;
    a.fAcc.estat = dEff(static_cast<int>(a.fAccRecoYield.val), static_cast<int>(a.fAccGenYield.val));
  }

  if (a.fAccRecoYield.val > 0) {
    a.fEffCand.val   = a.fAccCandYield.val/a.fAccRecoYield.val;
    a.fEffCand.estat = dEff(static_cast<int>(a.fAccCandYield.val), static_cast<int>(a.fAccRecoYield.val));
  }

  // -- that's it. The remaining numbers will be filled from histograms, based on the events tree

  cout << "AccGenFileYield: " << a.fAccGenFileYield.val << endl;
  cout << "AccGenYield:     " << a.fAccGenYield.val << endl;
  cout << "AccRecoYield:    " << a.fAccRecoYield.val << endl;
  cout << "candYield:       " << a.fAccCandYield.val << endl;
  cout << "acc:             " << a.fAcc.val << endl;
  cout << "effCand:         " << a.fEffCand.val << endl;
}


// ----------------------------------------------------------------------
// fills the number of events expected (total integral over entire histogram!)
// (corresponding to the efftot given, careful about mass cuts!)
// For example:
//          BF(Bs -> mu mu)  epstot(Bs)              ( fs)
//   n_s = ----------------- ----------  N(B+) pRatio(=--)
//          BF(B+ -> mu muK) epstot(B+)		     ( fu)
// ----------------------------------------------------------------------
void plotResults::scaleYield(anaNumbers &aSig, anaNumbers &aNorm, double pRatio) {
  cout << "+++ scaleYield: " << aSig.fNameMc << " wrt " << aNorm.fNameMc << endl;
  double sgBf  = fDS[aSig.fNameMc]->fBf;
  double sgBfE = fDS[aSig.fNameMc]->fBfE;

  double noBf  = fDS[aNorm.fNameMc]->fBf;
  double noBfE = fDS[aNorm.fNameMc]->fBfE;

  double yield  = (sgBf/noBf) * pRatio * (aSig.fEffTot.val/aNorm.fEffTot.val) * aNorm.fSignalFit.val;
  //  -- FIXME: error propagation!
  aSig.fScaledYield.val = yield;
  aSig.fScaledYield.setErrors(0.05*yield, 0.05*yield);

  // -- scale fMcYield
  if (aSig.fMcYield[4].val > 0.) {
    double scaleFactor = yield/aSig.fMcYield[4].val;
    cout << "+++ fMcYield[4] = " << aSig.fMcYield[4].val << " new scaled signal yield = " << yield << " -> scale factor = " << scaleFactor << endl;
    for (int i = 0; i < NWIN; ++i) {
      cout << "[" << i << "] " << aSig.fMcYield[i].val << " -> ";
      aSig.fMcYield[i].val = aSig.fMcYield[i].val*scaleFactor;
      cout << aSig.fMcYield[i].val;
    }
    cout << endl;
  } else {
    cout << "+++ fMcYield[4] = " << aSig.fMcYield[4].val << " (not larger than 0), not scaling fMcYield" << endl;
  }

  cout << "+++ sgBf = " << sgBf << " noBf = " << noBf << " pRatio = " << pRatio
       << " sgEffTot = " << aSig.fEffTot.val
       << " noEffTot = " << aNorm.fEffTot.val
       << " noSignalFit = " << aNorm.fSignalFit.val
       << " ->  scaled Yield = " << aSig.fScaledYield.val
       << endl;
}




// ----------------------------------------------------------------------
void plotResults::calculateNumbers(string mode) {
  cout << "==> calculateNumbers for mode: " << mode << endl;

  if (string::npos != mode.find("cnc")) {
    fDoUseBDT = false;
  } else {
    fDoUseBDT = true;
  }
  // -- open histogram file
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  for (unsigned int chan = 0; chan < fNchan; ++chan) {
    cout << "calculateNumbers for channel " << chan << endl;
    fChan = chan;
    initNumbers(fNoNumbers[chan]);
    initNumbers(fCsNumbers[chan]);
    initNumbers(fB0Numbers[chan]);
    initNumbers(fBsmmNumbers[chan]);
    initNumbers(fBdmmNumbers[chan]);
    initNumbers(fHhNumbers[chan]);
    initNumbers(fSlNumbers[chan]);
    initNumbers(fCombNumbers[chan]);
    initNumbers(fNpNumbers[chan]);
    initNumbers(fBgNumbers[chan]);
    initNumbers(fSgAndBgNumbers[chan]);

    // -- first to provide scaleYield base
    calculateB2JpsiNumbers(fNoNumbers[chan]);
    calculateB2JpsiNumbers(fCsNumbers[chan]);

    // -- before rare backgrounds to provide trigger efficiency per channel!
    calculateSgNumbers(fBsmmNumbers[chan]);
    calculateSgNumbers(fBdmmNumbers[chan]);
    calculateCombBgNumbers(fCombNumbers[chan]);

    // -- and finally the rare backgrounds
    calculateRareBgNumbers(chan);

    // -- do two channels only
    if (chan == 1) break;
  }

  fHistFile->Close();


}


// ----------------------------------------------------------------------
void plotResults::numbersFromHist(anaNumbers &aa, string syst) {
  // previously mode had been used to differentiate between bdmm, bsmm, no, cs, ...
  int chan = aa.fChan;

  // -- efficiency and acceptance
  string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  modifier += "_" + aa.fNameMc;

  cout << " numbersFromHist for name: " << aa.fName << ", chan: " << aa.fChan
       << " all cuts hist: " << Form("%s_%s_chan%d", fHistWithAllCuts.c_str(), modifier.c_str(), chan)
       << endl;

  fHistFile->cd(aa.fNameMc.c_str());
  TH1D *hAcceptance              = (TH1D*)gDirectory->Get(Form("hGenAndAccNumbers_%s_chan%d", modifier.c_str(), chan));

  TH1D *hMassAbsNoCuts           = (TH1D*)gDirectory->Get(Form("hMassAbsNoCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassNoCuts              = (TH1D*)gDirectory->Get(Form("hMassNoCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassWithAnaCuts         = (TH1D*)gDirectory->Get(Form("hMassWithAnaCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassWithMuonCuts        = (TH1D*)gDirectory->Get(Form("hMassWithMuonCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassWithTriggerCuts     = (TH1D*)gDirectory->Get(Form("hMassWithTriggerCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassWithAllCuts         = (TH1D*)gDirectory->Get(Form("%s_%s_chan%d", fHistWithAllCuts.c_str(), modifier.c_str(), chan));
  TH1D *hMassWithMassCuts        = (TH1D*)gDirectory->Get(Form("hMassWithMassCuts_%s_chan%d", modifier.c_str(), chan));

  TH1D *hMuId                    = (TH1D*)gDirectory->Get(Form("hMuId_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMuIdMC                  = (TH1D*)gDirectory->Get(Form("hMuIdMC_%s_chan%d", modifier.c_str(), chan));

  TH1D *hMuTr                    = (TH1D*)gDirectory->Get(Form("hMuTr_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMuTrMC                  = (TH1D*)gDirectory->Get(Form("hMuTrMC_%s_chan%d", modifier.c_str(), chan));
  fHistFile->cd();

  double effFilter    = getValueByLabel(hAcceptance, "effFilter");
  double effFilterE   = getValueByLabel(hAcceptance, "effFilterE");
  if (effFilter < 1e-6) {
    cout << "resetting effFilter to 1" << endl;
    effFilter = 1.0;
  }
  aa.fEffFilter.val = effFilter;
  aa.fEffFilter.setErrors(effFilterE, effFilter*0.001);


  double effGenSel    = getValueByLabel(hAcceptance, "effGenSel");
  double effGenSelE   = getValueByLabel(hAcceptance, "effGenSelE");
  if (effGenSel < 1e-6) {
    cout << "resetting effGenSel to 1" << endl;
    effGenSel = 1.0;
  }
  aa.fEffGenSel.val = effGenSel;
  aa.fEffGenSel.setErrors(effGenSelE, effGenSelE);


  double genFileYield = getValueByLabel(hAcceptance, "genFileYield");

  aa.fAccGenFileYield.val = getValueByLabel(hAcceptance, "accGenFileYield");
  aa.fAccGenYield.val     = getValueByLabel(hAcceptance, "accGenYield");
  aa.fGenFileYield.val    = getValueByLabel(hAcceptance, "genFileYield");
  aa.fGenYield.val        = getValueByLabel(hAcceptance, "genYield");
  aa.fGenYield.estat      = getValueByLabel(hAcceptance, "genYieldE");
  aa.fAccRecoYield.val    = getValueByLabel(hAcceptance, "accRecoYield");
  aa.fAccCandYield.val    = getValueByLabel(hAcceptance, "accCandYield");

  aa.fAcc.val             = getValueByLabel(hAcceptance, "acc");
  aa.fAcc.setErrors(getValueByLabel(hAcceptance, "accE"), aa.fAcc.val * fSystematics["acceptance"][chan]);

  aa.fEffCand.val         = getValueByLabel(hAcceptance, "effCand");
  aa.fEffCand.setErrors(getValueByLabel(hAcceptance, "effCandE"), aa.fEffCand.val * fSystematics["effcand" + syst][chan]);

  double a = massIntegral(hMassNoCuts, ALL, chan);
  double b = massIntegral(hMassWithAnaCuts, ALL, chan);
  double c = massIntegral(hMassWithMuonCuts, ALL, chan);
  double d = massIntegral(hMassWithTriggerCuts, ALL, chan);
  double e = massIntegral(hMassWithAllCuts, ALL, chan);
  double f = massIntegral(hMassWithMassCuts, ALL, chan);

  aa.fAccCandYield.val     = getValueByLabel(hAcceptance, "accCandYield");

  aa.fCandYield.val        = a;
  aa.fAnaYield.val         = b;
  aa.fEffAna.val           = b/a;
  aa.fEffAna.setErrors(dEff(static_cast<int>(b), static_cast<int>(a)), aa.fEffAna.val * fSystematics["effana" + syst][chan]);

  aa.fMuidYield.val        = c;
  aa.fEffMuidMC.val        = c/b;
  aa.fEffMuidMC.setErrors(dEff(static_cast<int>(c), static_cast<int>(b)), aa.fEffMuidMC.val * fSystematics["effmuid"][chan]);

  aa.fTrigYield.val        = d;
  aa.fEffTrigMC.val        = d/c;
  aa.fEffTrigMC.setErrors(dEff(static_cast<int>(c), static_cast<int>(b)),  aa.fEffTrigMC.val * fSystematics["effmuid"][chan]);

  aa.fEffTot.val           = e/aa.fGenYield.val;
  aa.fEffTot.setErrors(dEff(e, TMath::Sqrt(e), aa.fGenYield.val, aa.fGenYield.estat), quadraticSum(4, fSystematics["acceptance"][chan]
											      , fSystematics["effcand" + syst][chan]
											      , fSystematics["effana" + syst][chan]
											      , fSystematics["effmuid"][chan]
											      , fSystematics["efftrig"][chan])  * aa.fEffTot.val);

  // -- prod efficiencies
  aa.fEffProdMC.val        = aa.fAcc.val * aa.fEffCand.val * aa.fEffAna.val * aa.fEffMuidMC.val * aa.fEffTrigMC.val;
  aa.fEffProdMC.setErrors(aa.fEffProdMC.val * quadraticSum(5,
							   aa.fAcc.estat/aa.fAcc.val,
							   aa.fEffAna.estat/aa.fEffAna.val,
							   aa.fEffCand.estat/aa.fEffCand.val,
							   aa.fEffMuidMC.estat/aa.fEffMuidMC.val,
							   aa.fEffTrigMC.estat/aa.fEffTrigMC.val),
			  aa.fEffProdMC.val * quadraticSum(5,
							   aa.fAcc.estat/aa.fAcc.val,
							   aa.fEffAna.estat/aa.fEffAna.val,
							   aa.fEffCand.estat/aa.fEffCand.val,
							   aa.fEffMuidMC.estat/aa.fEffMuidMC.val,
							   aa.fEffTrigMC.estat/aa.fEffTrigMC.val)
			  );

  aa.fTotGenYield.val     = e/(aa.fEffProdMC.val);
  aa.fTotGenYield.setErrors(0.001*aa.fTotGenYield.val, 0.001*aa.fTotGenYield.val); //FIXME
  aa.fProdGenYield.val     = e/(aa.fEffTot.val);
  aa.fProdGenYield.setErrors(0.001*aa.fProdGenYield.val, 0.001*aa.fProdGenYield.val); //FIXME

  double tot   = massIntegral(hMassWithAllCuts, ALL, chan);
  double lo    = massIntegral(hMassWithAllCuts, LO, chan);
  double hi    = massIntegral(hMassWithAllCuts, HI, chan);
  double bd    = massIntegral(hMassWithAllCuts, BD, chan);
  double bs    = massIntegral(hMassWithAllCuts, BS, chan);

  aa.fFrac[0].val   = lo/tot;
  aa.fFrac[0].estat = dEff(static_cast<int>(lo), static_cast<int>(tot));
  aa.fFrac[1].val   = bd/tot;
  aa.fFrac[1].estat = dEff(static_cast<int>(bd), static_cast<int>(tot));
  aa.fFrac[2].val   = bs/tot;
  aa.fFrac[2].estat = dEff(static_cast<int>(bs), static_cast<int>(tot));
  aa.fFrac[3].val   = hi/tot;
  aa.fFrac[3].estat = dEff(static_cast<int>(hi), static_cast<int>(tot));

  aa.fMcYield[0].val   = lo;
  aa.fMcYield[1].val   = bd;
  aa.fMcYield[2].val   = bs;
  aa.fMcYield[3].val   = hi;
  aa.fMcYield[4].val   = tot;
}


// ----------------------------------------------------------------------
void plotResults::calculateB2JpsiNumbers(anaNumbers &a) {
  c0->Clear();
  string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  fSuffixSel = modifier;
  cout << "==> calculateB2JpsiNumbers for name: " << a.fName << ", chan: " << a.fChan << " fSuffixSel: " << fSuffixSel << endl;

  // -- MC: efficiency and acceptance
  char mode[200];
  sprintf(mode, "%s", a.fName.c_str());
  fSetup = a.fNameMc;
  int chan = a.fChan;
  if (string::npos != a.fName.find("bupsik")) {
    numbersFromHist(a, "bupsik");
  } else if (string::npos != a.fName.find("bspsiphi")) {
    numbersFromHist(a, "bspsiphi");
  }
  // -- cross check the scaled yield normalized to B+ -> J/psi K+
  if (string::npos != a.fName.find("bspsiphi")) {
    double pRatio(fFsfu.val);
    scaleYield(a, fNoNumbers[chan], pRatio);
  }
  // -- data: fit yields
  fSetup = a.fNameDa;
  string  name = Form("hNorm_%s_%s_chan%d", modifier.c_str(), fSetup.c_str(), chan);
  bool ok = fHistFile->cd(fSetup.c_str());
  cout << "cd to " << fSetup << ": " << ok << endl;
  fitPsYield fpy(name, 0);
  if (string::npos != a.fName.find("bupsik")) {
    fpy.fitBu2JpsiKp(5, fDirectory + "/");
  } else if (string::npos != a.fName.find("bspsiphi")) {
    fpy.fitBs2JpsiPhi(5, fDirectory + "/");
  }
  a.fSignalFit.val   = fpy.getSignalYield();
  a.fSignalFit.estat = fpy.getSignalError();
  a.fSignalFit.esyst = fSystematics["norm" + a.fName][chan] * fpy.getSignalYield();
  a.fSignalFit.etot  = TMath::Sqrt(a.fSignalFit.estat*a.fSignalFit.estat + a.fSignalFit.esyst*a.fSignalFit.esyst);

  printNumbers(a, cout);

  cout << "chan " << a.fChan << " total yield: " << fpy.getSignalYield() << " +/- "  << fpy.getSignalError() << endl;
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- NORMALIZATION: " << mode << " chan: " << chan << endl;
  dumpTex(a.fEffGenSel, Form("%s:GENSEL-%s-chan%d", fSuffixSel.c_str(), mode, chan), 6);
  dumpTex(a.fAcc, Form("%s:ACC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 6);
  dumpTex(a.fEffCand, Form("%s:EFF-CAND-%s-chan%d", fSuffixSel.c_str(), mode, chan), 6);
  dumpTex(a.fEffMuidMC, Form("%s:EFF-MU-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 6);
  dumpTex(a.fEffTrigMC, Form("%s:EFF-TRIG-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 6);
  dumpTex(a.fEffAna, Form("%s:EFF-ANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), 6);
  dumpTex(a.fEffProdMC, Form("%s:EFF-PRODMC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 6);
  dumpTex(a.fEffTot, Form("%s:EFF-TOT-%s-chan%d", fSuffixSel.c_str(), mode, chan), 6);

  dumpTex(a.fGenFileYield, Form("%s:N-GENFILEYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fGenYield, Form("%s:N-GENYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fTotGenYield, Form("%s:N-TOTYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fProdGenYield, Form("%s:N-PRODYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);

  if (string::npos != a.fName.find("bspsiphi")) {
    dumpTex(a.fScaledYield, Form("%s:N-SCALEDYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  }
  dumpTex(a.fSignalFit, Form("%s:N-OBS-%s-chan%d", fSuffixSel.c_str(), mode, chan), 1);

  c0->Modified();
  c0->Update();
  fHistFile->cd();
}


// ----------------------------------------------------------------------
void plotResults::calculateCombBgNumbers(anaNumbers &a, int mode, double lo, double hi) {
  string hname = fHistWithAllCuts;

  c0->Clear();
  string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  fSuffixSel = modifier;
  cout << "==> calculateCombBgNumbers for name: " << a.fName << ", chan: " << a.fChan << " fSuffixSel = " << fSuffixSel << endl;

  // -- get the histogram
  fSetup = a.fNameDa;
  string  name = Form("%s_%s_%s_chan%d", hname.c_str(), modifier.c_str(), fSetup.c_str(), a.fChan);
  TH1D *h1 = (TH1D*)gDirectory->Get(Form("%s/%s", fSetup.c_str(), name.c_str()));
  cout << "getting  histogram ->" << Form("%s/%s", fSetup.c_str(), name.c_str()) << "<-" << endl;
  TF1 *lF1(0), *lF2(0);

  if (0 == mode) {
    lF1 = fIF->pol0BsBlind(h1);
    lF2 = fIF->pol0(h1);
  }

  lF2->SetLineStyle(kDashed);
  double binw = h1->GetBinWidth(1);
  TFitResultPtr r;
  bool notFit(false);
  if (massIntegral(h1, HI, a.fChan) > 3) {
    h1->Fit(lF1, "rl", "", lo, hi);
  } else {
    notFit = true;
    lF1->SetParameter(0, 3.*binw/(hi-lo));              // A = p0 * (hi-lo), N = A/binw -> p0 = A/(hi-lo) = N*binw/(hi-lo)
    lF1->SetParError(0, TMath::Sqrt(3.)*binw/(hi-lo));  // A = p0 * (hi-lo), N = A/binw -> p0 = A/(hi-lo) = N*binw/(hi-lo)
  }
  setTitles(h1, "#it{m}_{#it{#mu #mu}} [GeV]", Form("Candidates / %4.3f GeV", h1->GetBinWidth(1)));

  h1->DrawCopy();
  lF2->SetParameters(lF1->GetParameters());
  lF2->SetParErrors(lF1->GetParErrors());
  lF2->SetLineColor(kBlue);
  if (notFit) lF1->Draw("same");
  lF2->Draw("same");

  a.fObsYield[0].val = massIntegral(h1, LO, a.fChan);
  a.fFitYield[0].val = lF2->Integral(fBgLo, fCuts[a.fChan]->mBdLo)/binw;

  a.fObsYield[1].val = massIntegral(h1, BD, a.fChan);
  a.fFitYield[1].val = lF2->Integral(fCuts[a.fChan]->mBdLo, fCuts[a.fChan]->mBdHi)/binw;

  a.fObsYield[2].val = massIntegral(h1, BS, a.fChan);
  a.fFitYield[2].val = lF2->Integral(fCuts[a.fChan]->mBsLo, fCuts[a.fChan]->mBsHi)/binw;

  a.fObsYield[3].val = massIntegral(h1, HI, a.fChan);
  a.fFitYield[3].val = lF2->Integral(fCuts[a.fChan]->mBsHi, fBgHi)/binw;

  a.fObsYield[4].val = massIntegral(h1, ALL, a.fChan);
  a.fFitYield[4].val = lF2->Integral(fBgLo, fBgHi)/binw;

  // -- dump numbers
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- DATA/COMBINATORIAL " << mode << " chan " << a.fChan << endl;
  for (unsigned i = 0; i < NWIN; ++i) {
    dumpTex(a.fFitYield[i], Form("%s:N-FIT-MBIN%d-CB-chan%d", fSuffixSel.c_str(), i, a.fChan), 3);
    dumpTex(a.fObsYield[i], Form("%s:N-OBS-MBIN%d-DATA-chan%d", fSuffixSel.c_str(), i, a.fChan), 0);
  }

  double tsize = tl->GetTextSize();
  tl->SetTextSize(0.03);
  tl->DrawLatexNDC(0.6, 0.80, Form("Blind %s", fSuffixSel.c_str()));
  tl->DrawLatexNDC(0.6, 0.76, Form("chan: %d", a.fChan));
  tl->DrawLatexNDC(0.6, 0.72, Form("mode: %d", mode));
  tl->SetTextSize(tsize);
  c0->Modified();
  c0->Update();
  savePad(Form("%s-combBg-mode%d-chan%d.pdf", fSuffixSel.c_str(), mode, a.fChan));
}

// ----------------------------------------------------------------------
void plotResults::calculateSgNumbers(anaNumbers &a) {
  cout << "==> calculateSgNumbers for name: " << a.fName << ", chan: " << a.fChan << endl;
  c0->Clear();
  string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  fSuffixSel = modifier;

  // -- MC: efficiency and acceptance
  char mode[200];
  sprintf(mode, "%s", a.fName.c_str());
  fSetup = a.fNameMc;
  int chan = a.fChan;
  numbersFromHist(a, "bsmm");
  double pRatio(fFsfu.val);
  if (string::npos != a.fName.find("bdmm")) pRatio = 1.;
  scaleYield(a, fNoNumbers[chan], pRatio);

  // -- data: fitted/interpolated yields
  fSetup = a.fNameDa;

  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- SIGNAL " << mode << " chan " << chan << endl;
  dumpTex(a.fEffGenSel, Form("%s:GENSEL-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
  dumpTex(a.fAcc, Form("%s:ACC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
  dumpTex(a.fEffCand, Form("%s:EFF-CAND-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
  dumpTex(a.fEffMuidMC, Form("%s:EFF-MU-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
  dumpTex(a.fEffTrigMC, Form("%s:EFF-TRIG-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
  dumpTex(a.fEffAna, Form("%s:EFF-ANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
  dumpTex(a.fEffProdMC, Form("%s:EFF-PRODMC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
  dumpTex(a.fEffTot, Form("%s:EFF-TOT-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);

  dumpTex(a.fGenFileYield, Form("%s:N-GENFILEYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fGenYield, Form("%s:N-GENYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fTotGenYield, Form("%s:N-TOTYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fProdGenYield, Form("%s:N-PRODYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fScaledYield, Form("%s:N-SCALEDYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);

  for (unsigned i = 0; i < a.fMcYield.size(); ++i) {
    dumpTex(a.fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-%s-chan%d", fSuffixSel.c_str(), i, mode, chan), 3);
  }

  //  dumpTex(a.fSignalFit, Form("%s:N-OBS-%s-chan%d", fSuffixSel.c_str(), mode, chan), 1);

  c0->Modified();
  c0->Update();
  fHistFile->cd();

}

// ----------------------------------------------------------------------
void plotResults::calculateRareBgNumbers(int chan) {
  string hname = fHistWithAllCuts; // "hMassWithAllCutsBlind";
  string uname = fHistWithAllCuts; // "hMassWithAllCuts";
  string wname = fHistWithAllCuts; replaceAll(wname, "hMass", "hW8Mass");

  cout << "==> calculateRareBgNumbers for chan: " << chan << endl;
  int nloop(0);
  string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  fSuffixSel = modifier;
  char mode[200];
  // -- book summed histograms
  THStack hSl("hSl", "rare semileptonic decays");
  TH1D *h1Sl = new TH1D("h1Sl", "h1Sl", NBINS/5, fMassLo, fMassHi); h1Sl->Sumw2();
  vector<TH1*> vSl;
  vector<string> vSlnames, vSloptions;
  THStack hHh("hHh", "rare peaking decays");
  TH1D *h1Hh = new TH1D("h1Hh", "h1Hh", NBINS/5, fMassLo, fMassHi); h1Hh->Sumw2();
  vector<TH1*> vHh;
  vector<string> vHhnames, vHhoptions;
  THStack hBg("hBg", "rare decays");
  vector<TH1*> vBg;
  vector<string> vBgnames, vBgoptions;

  // -- loop over all (effective) two-body backgrounds
  TTree *t(0);
  string u8name(""), w8name(""), accname("bs");
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    if (skipThisBg(it->first)) continue;
    cout << "calculateRareBgNumbers for " << it->first << " in chan = " << chan << endl;
    anaNumbers *a = fRareNumbers[it->first][chan];
    string accname = rareAccName(it->first);
    cout << "rec: " << Form("%s/hGenAndAccNumbers_%s_%s_chan%d", it->first.c_str(), modifier.c_str(), it->first.c_str(), chan) << endl;
    TH1D *h1 = (TH1D*)fHistFile->Get(Form("%s/hGenAndAccNumbers_%s_%s_chan%d", it->first.c_str(), modifier.c_str(), it->first.c_str(), chan));
    double effGenSel = getValueByLabel(h1,  "effGenSel");
    double genYield  = getValueByLabel(h1,  "genYield");
    double acc       = getValueByLabel(h1,  "acc");

    u8name = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), uname.c_str(), modifier.c_str(), it->first.c_str(), chan);
    w8name = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), wname.c_str(), modifier.c_str(), it->first.c_str(), chan);
    TH1D *hu = (TH1D*)fHistFile->Get(u8name.c_str());
    TH1D *hw = (TH1D*)fHistFile->Get(w8name.c_str());

    double utot     = massIntegral(hu, ALL, chan);
    double ulo      = massIntegral(hu, LO, chan);
    double ubd      = massIntegral(hu, BD, chan);
    double ubs      = massIntegral(hu, BS, chan);
    double uhi      = massIntegral(hu, HI, chan);
    double eAccAna  = utot/genYield; // this does NOT yet include misid and trigger efficiency, will be applied later on
    double eTot     = eAccAna * (massIntegral(hw, ALL, chan) / utot) * fBsmmNumbers[chan].fEffTrigMC.val;

    fSetup = a->fNameMc;
    a->fEffAccAna.val = eAccAna;
    a->fEffTot.val = eTot;

    if (utot > 0) {
      a->fFrac[0].val   = ulo/utot;
      a->fFrac[0].estat = dEff(static_cast<int>(ulo), static_cast<int>(utot));
      a->fFrac[1].val   = ubd/utot;
      a->fFrac[1].estat = dEff(static_cast<int>(ubd), static_cast<int>(utot));
      a->fFrac[2].val   = ubs/utot;
      a->fFrac[2].estat = dEff(static_cast<int>(ubs), static_cast<int>(utot));
      a->fFrac[3].val   = uhi/utot;
      a->fFrac[3].estat = dEff(static_cast<int>(uhi), static_cast<int>(utot));
    } else {
      a->fFrac[0].val = 0.;
      a->fFrac[1].val = 0.;
      a->fFrac[2].val = 0.;
      a->fFrac[3].val = 0.;
      a->fFrac[4].val = 0.;
    }

    double pRatio(0.);
    if (string::npos != it->first.find("bs")) pRatio = fFsfu.val;
    if (string::npos != it->first.find("lb")) pRatio = fFsfu.val; // FIXME!
    if (string::npos != it->first.find("bd")) pRatio = 1.;
    if (string::npos != it->first.find("bu")) pRatio = 1.;

    scaleYield(*a, fNoNumbers[chan], pRatio);
    hw->Scale(a->fScaledYield.val/massIntegral(hw, ALL, chan));

    int nmuons = ((string::npos != it->first.find("mu")) ? 1 : 0);
    if (string::npos != it->first.find("mumu")) nmuons = 2;
    double muonSys(0.), fsfuSys(0.);
    if (0 == nmuons) {
      muonSys = 2.*fSystematics["fakemuid"][chan];
    } else if (1 == nmuons) {
      muonSys = quadraticSum(2, fSystematics["fakemuid"][chan], fSystematics["effmuid"][chan]);
    } else {
      muonSys = 2.*fSystematics["effmuid"][chan];
    }
    if (pRatio < 0.9) {
      fsfuSys = fSystematics["fsfu"][chan];
    }
    double esystRel = quadraticSum(7,
				   fsfuSys,
				   fSystematics["acceptance"][chan],
				   fSystematics["effcandbsmm"][chan],
				   fSystematics["effanabsmm"][chan],
				   fSystematics["efftrig"][chan],
				   muonSys,
				   fSystematics["normbupsik"][chan]
				   );

    a->fMcYield[0].val = massIntegral(hw, LO, chan);
    double integral = massIntegral(hu, LO, chan);
    double estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fMcYield[0].setErrors(a->fMcYield[0].val*estatRel, a->fMcYield[0].val*esystRel);

    a->fMcYield[1].val = massIntegral(hw, BD, chan);
    integral = massIntegral(hu, BD, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fMcYield[1].setErrors(a->fMcYield[1].val*estatRel, a->fMcYield[1].val*esystRel);

    a->fMcYield[2].val = massIntegral(hw, BS, chan);
    integral = massIntegral(hu, BS, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    // FIXME??
    if (a->fName == "bdpimunuBg") {
      cout << "XXX a->fMcYield[2].val = " << a->fMcYield[2].val << " integral = " << integral << " estatRel = " << estatRel << endl;
    }
    a->fMcYield[2].setErrors(a->fMcYield[2].val*estatRel, a->fMcYield[2].val*esystRel);

    a->fMcYield[3].val = massIntegral(hw, HI, chan);
    integral = massIntegral(hu, HI, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fMcYield[3].setErrors(a->fMcYield[3].val*estatRel, a->fMcYield[3].val*esystRel);

    a->fMcYield[4].val = massIntegral(hw, ALL, chan);
    integral = massIntegral(hu, ALL, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fMcYield[4].setErrors(a->fMcYield[4].val*estatRel, a->fMcYield[4].val*esystRel);

    setFilledHist(hw, it->second->fLcolor, it->second->fFcolor, it->second->fFillStyle, 1);
    TH1D *hwrb = (TH1D*)hw->Clone(Form("hwrb_%s", hw->GetName()));
    hwrb->Rebin(5);
    hwrb->SetNdivisions(510, "XYZ");
    setTitles(hwrb, "m_{#it{#mu #mu}} [GeV]", Form("Candidates / %4.3f GeV", hwrb->GetBinWidth(1)));
    if (0 == nmuons) {
      for (unsigned im = 0; im < a->fMcYield.size(); ++im) {
	fHhNumbers[chan].fMcYield[im].val += a->fMcYield[im].val;
	fHhNumbers[chan].fMcYield[im].add2Errors(a->fMcYield[im]);
      }
      h1Hh->Add(hwrb);
      hHh.Add(hwrb);
      vHh.insert(vHh.begin(), hw);
      vHhnames.insert(vHhnames.begin(), it->second->fName);
      vHhoptions.insert(vHhoptions.begin(), "f");
      hBg.Add(hwrb);
      vBg.insert(vBg.begin(), hw);
      vBgnames.insert(vBgnames.begin(), it->second->fName);
      vBgoptions.insert(vBgoptions.begin(), "f");
    } else if (nmuons >= 1) {
      for (unsigned im = 0; im < a->fMcYield.size(); ++im) {
	fSlNumbers[chan].fMcYield[im].val += a->fMcYield[im].val;
	fSlNumbers[chan].fMcYield[im].add2Errors(a->fMcYield[im]);
      }
      h1Sl->Add(hwrb);
      hSl.Add(hwrb);
      vSl.insert(vSl.begin(), hw);
      vSlnames.insert(vSlnames.begin(), it->second->fName);
      vSloptions.insert(vSloptions.begin(), "f");
      hBg.Add(hwrb);
      vBg.insert(vBg.begin(), hw);
      vBgnames.insert(vBgnames.begin(), it->second->fName);
      vBgoptions.insert(vBgoptions.begin(), "f");
    }
    sprintf(mode, "%s", a->fName.c_str());
    fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    fTEX << "% -- RARE " << mode << " chan " << chan << endl;
    dumpTex(a->fEffGenSel, Form("%s:GENSEL-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
    dumpTex(a->fAcc, Form("%s:ACC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
    dumpTex(a->fEffCand, Form("%s:EFF-CAND-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
    dumpTex(a->fEffMuidMC, Form("%s:EFF-MU-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
    dumpTex(a->fEffTrigMC, Form("%s:EFF-TRIG-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
    dumpTex(a->fEffAna, Form("%s:EFF-ANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
    dumpTex(a->fEffProdMC, Form("%s:EFF-PRODMC-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
    dumpTex(a->fEffAccAna, Form("%s:EFF-ACCANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);
    dumpTex(a->fEffTot, Form("%s:EFF-TOT-%s-chan%d", fSuffixSel.c_str(), mode, chan), 5);

    dumpTex(a->fGenFileYield, Form("%s:N-GENFILEYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
    dumpTex(a->fGenYield, Form("%s:N-GENYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
    dumpTex(a->fTotGenYield, Form("%s:N-TOTYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
    dumpTex(a->fProdGenYield, Form("%s:N-PRODYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);

    dumpTex(a->fScaledYield, Form("%s:N-SCALEDYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);

    for (unsigned i = 0; i < a->fMcYield.size(); ++i) {
      dumpTex(a->fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-%s-chan%d", fSuffixSel.c_str(), i, mode, chan), 3);
    }
    hw->SetMinimum(0.);
    hw->Draw("hist");
    c0->Modified();
    c0->Update();
    if (1) {
      savePad(Form("%s.pdf", hw->GetName()));
    }
    ++nloop;
    //    if (nloop > 3) break;
  }

  // -- add combinatorial background for more combined sums
  for (unsigned im = 0; im < NWIN; ++im) {
    fNpNumbers[chan].fMcYield[im].val +=  fSlNumbers[chan].fMcYield[im].val;
    fNpNumbers[chan].fMcYield[im].val +=  fCombNumbers[chan].fFitYield[im].val;

    fBgNumbers[chan].fMcYield[im].val +=  fHhNumbers[chan].fMcYield[im].val;
    fBgNumbers[chan].fMcYield[im].val +=  fNpNumbers[chan].fMcYield[im].val;

    fSgAndBgNumbers[chan].fMcYield[im].val +=  fNpNumbers[chan].fMcYield[im].val;
  }

  c0->Clear();
  hSl.Draw("hist");
  hSl.GetXaxis()->SetTitle("#it{m}_{#it{#mu #mu}} [GeV]");
  hSl.GetYaxis()->SetTitle("Candidates / Bin");
  TLegend *lSl = ::newLegend("semileptonic decays", 0.50, 0.6, 0.85, 0.85, vSl, vSlnames, vSloptions);
  lSl->Draw();
  c0->Modified();
  c0->Update();
  if (1) {
    savePad(Form("%s-rare-sl-chan%d.pdf", fSuffixSel.c_str(), chan));
  }


  c0->Clear();
  hHh.Draw("hist");
  hHh.GetXaxis()->SetTitle("#it{m}_{#it{#mu #mu}} [GeV]");
  hHh.GetYaxis()->SetTitle("Candidates / Bin");
  TLegend *lHh = ::newLegend("hadronic decays", 0.56, 0.4, 0.85, 0.85, vHh, vHhnames, vHhoptions);
  lHh->Draw();
  c0->Modified();
  c0->Update();
  if (1) {
    savePad(Form("%s-rare-hh-chan%d.pdf", fSuffixSel.c_str(), chan));
  }


  c0->Clear();
  hBg.Draw("hist");
  hBg.GetXaxis()->SetTitle("#it{m}_{#it{#mu #mu}} [GeV]");
  hBg.GetYaxis()->SetTitle("Candidates / Bin");
  TLegend *lBg = ::newLegend("rare decays", 0.50, 0.3, 0.85, 0.85, vBg, vBgnames, vBgoptions);
  lBg->Draw();
  c0->Modified();
  c0->Update();
  if (1) {
    savePad(Form("%s-rare-bg-chan%d.pdf", fSuffixSel.c_str(), chan));
  }

  // -- create overlay of data and the stacked (and scaled) backgrounds
  fSetup = fCombNumbers[chan].fNameDa;
  string  name = Form("%s_%s_%s_chan%d", hname.c_str(), modifier.c_str(), fSetup.c_str(), chan);
  cout << "getting  histogram ->" << Form("%s/%s", fSetup.c_str(), name.c_str()) << "<-" << endl;
  TH1D *h1e = (TH1D*)((TH1D*)gDirectory->Get(Form("%s/%s", fSetup.c_str(), name.c_str())))->Clone("h1erebin");
  h1e->Rebin(5);
  TH1D *h1c = (TH1D*)((TH1D*)gDirectory->Get(Form("%s/%s", fSetup.c_str(), name.c_str())))->Clone("h1crebin");
  h1c->Clear();
  h1c->Rebin(5);
  double bgtot(0.);
  for (int i = 0; i < 4; ++i) {
    bgtot += fCombNumbers[chan].fFitYield[i].val;
  }
  bgtot = bgtot/h1c->GetNbinsX();
  for (int i = 1; i <= h1c->GetNbinsX(); ++i) {
    h1c->SetBinContent(i, bgtot);
  }

  setHist(h1e, kBlack, 0, 0, 2.5);
  h1e->Draw("hist");
  THStack hCombBg("hCombBg", "comb + rare decays");
  setFilledHist(h1c, kGray, kGray, 1000);
  hCombBg.Add(h1c);
  setFilledHist(h1Sl, kYellow-3, kYellow-3, 1000);
  hCombBg.Add(h1Sl);
  setFilledHist(h1Hh, kMagenta-2, kMagenta-2, 1000);
  hCombBg.Add(h1Hh);

  newLegend(0.4, 0.6, 0.86, 0.86, "Background");
  legg->AddEntry(h1e, "Data", "l");
  legg->AddEntry(h1Hh, "rare hadronic", "f");
  legg->AddEntry(h1Sl, "rare semileptonic", "f");
  legg->AddEntry(h1c, "combinatorial", "f");
  legg->Draw();

  hCombBg.Draw("histsame");
  h1e->Draw("histsame");
  h1e->Draw("axissame");
  c0->Modified();
  c0->Update();
  if (1) {
    savePad(Form("%s-data-background-chan%d.pdf", fSuffixSel.c_str(), chan));
  }

  // -- and now the scaled version FIXME definition of scale
  double bgLo   = massIntegral(h1c, LO, chan) + massIntegral(h1Sl, LO, chan) + massIntegral(h1Hh, LO, chan);
  double diffLo = fCombNumbers[chan].fObsYield[0].val - bgLo;
  double scale  = diffLo/massIntegral(h1Sl, LO, chan);

  setHist(h1e, kBlack);
  h1e->Draw("hist");
  THStack hCombScBg("hCombScBg", "comb + scaled rare decays");
  h1Sl->Add(h1Sl, scale);
  setFilledHist(h1c, kGray, kGray, 1000);
  hCombScBg.Add(h1c);
  setFilledHist(h1Sl, kYellow-3, kYellow-3, 1000);
  hCombScBg.Add(h1Sl);
  setFilledHist(h1Hh, kMagenta-2, kMagenta-2, 1000);
  hCombScBg.Add(h1Hh);

  newLegend(0.4, 0.6, 0.86, 0.86, "Background");
  legg->AddEntry(h1e, "Data", "l");
  legg->AddEntry(h1Hh, "rare hadronic", "f");
  legg->AddEntry(h1Sl, "scal. rare semileptonic", "f");
  legg->AddEntry(h1c, "combinatorial", "f");
  legg->Draw();

  tl->SetTextSize(0.035);
  tl->DrawLatexNDC(0.65, 0.92, Form("scale = %3.1f", scale));

  hCombScBg.Draw("histsame");
  h1e->Draw("histsame");
  h1e->Draw("axissame");
  c0->Modified();
  c0->Update();
  if (1) {
    savePad(Form("%s-data-scaled-background-chan%d.pdf", fSuffixSel.c_str(), chan));
  }


  // -- calculate scaled numbers
  fSlNumbers[chan].fScaleFactor = 1. + scale;
  for (int iw = 0; iw < NWIN; ++iw) {
    fSlNumbers[chan].fFitYield[iw].val = fSlNumbers[chan].fMcYield[iw].val * fSlNumbers[chan].fScaleFactor;
    fNpNumbers[chan].fFitYield[iw].val = fNpNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
    fBgNumbers[chan].fFitYield[iw].val = fBgNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
    fSgAndBgNumbers[chan].fFitYield[iw].val = fSgAndBgNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
  }
  // -- add signal
  for (int iw = 0; iw < NWIN; ++iw) {
    fSgAndBgNumbers[chan].fFitYield[iw].val = fSgAndBgNumbers[chan].fMcYield[iw].val
      + fBsmmNumbers[chan].fMcYield[iw].val
      + fBdmmNumbers[chan].fMcYield[iw].val;
  }

  // -- dump combined/summed numbers: rare sl decays and rare hadronic (peaking) decays
    fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    fTEX << "% -- SUMMARY BACKGROUND " << mode << " chan " << chan << endl;
  for (unsigned i = 0; i < fHhNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fHhNumbers[chan].fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-HH-chan%d", fSuffixSel.c_str(), i, chan), 3);
  }

  fTEX << formatTex(fSlNumbers[chan].fScaleFactor, Form("%s:SCALEFACTOR-SL-chan%d:val", fSuffixSel.c_str(), chan), 3) << endl;
  for (unsigned i = 0; i < fSlNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fSlNumbers[chan].fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-SL-chan%d", fSuffixSel.c_str(), i, chan), 3);
    dumpTex(fSlNumbers[chan].fFitYield[i], Form("%s:N-FITSCLYIELD-MBIN%d-SL-chan%d", fSuffixSel.c_str(), i, chan), 3);
  }
  for (unsigned i = 0; i < fNpNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fNpNumbers[chan].fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-NP-chan%d", fSuffixSel.c_str(), i, chan), 3);
    dumpTex(fNpNumbers[chan].fFitYield[i], Form("%s:N-FITSCLYIELD-MBIN%d-NP-chan%d", fSuffixSel.c_str(), i, chan), 3);
  }
  for (unsigned i = 0; i < fBgNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fBgNumbers[chan].fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-BG-chan%d", fSuffixSel.c_str(), i, chan), 3);
    dumpTex(fBgNumbers[chan].fFitYield[i], Form("%s:N-FITSCLYIELD-MBIN%d-BG-chan%d", fSuffixSel.c_str(), i, chan), 3);
  }

  for (unsigned i = 0; i < fBgNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fSgAndBgNumbers[chan].fFitYield[i], Form("%s:N-FITSCLYIELD-MBIN%d-SGANDBG-chan%d", fSuffixSel.c_str(), i, chan), 3);
  }



  delete h1e;
  delete h1c;

}



// ----------------------------------------------------------------------
void plotResults::loopFunction1() {

  if (fChan < 0) return;
  double mass = fb.m;

  vector<string> modifier;
  modifier.push_back("cnc" + fSuffix);
  modifier.push_back("bdt" + fSuffix);

  for (unsigned int im = 0; im < modifier.size(); ++im) {
    fhMassAbsNoCuts[modifier[im]][fChan]->Fill(mass);
    fhAccAll[modifier[im]][fChan]->Fill(TMath::Abs(fb.eta), fb.pt);
    fhAccPtAll[modifier[im]][fChan]->Fill(fb.pt);
    fhAccEtaAll[modifier[im]][fChan]->Fill(TMath::Abs(fb.eta));
  }

  if (!fGoodAcceptance) return;

  for (unsigned int im = 0; im < modifier.size(); ++im) {
    fhAccPass[modifier[im]][fChan]->Fill(TMath::Abs(fb.eta), fb.pt);
    fhAccPtPass[modifier[im]][fChan]->Fill(fb.pt);
    fhAccEtaPass[modifier[im]][fChan]->Fill(TMath::Abs(fb.eta));
    // -- this is the base, after the raw acceptance cuts
    fhMassNoCuts[modifier[im]][fChan]->Fill(mass);
  }


  // -----------------
  // -- CNC histograms
  // -----------------
  if (fGoodQ
      && fGoodPvAveW8
      && fGoodTracks
      && fGoodTracksPt
      && fGoodTracksEta
      && fGoodMuonsPt
      && fGoodMuonsEta
      && fGoodJpsiCuts
      && fGoodMaxDoca
      && fGoodLip
      && fGoodLipS
      && fGoodIp
      && fGoodIpS
      && fGoodPt
      && fGoodEta
      && fGoodAlpha
      && fGoodChi2
      && fGoodFLS
      && fGoodCloseTrack
      && fGoodIso
      && fGoodDocaTrk
      ) {
    fhMassWithAnaCuts[modifier[0]][fChan]->Fill(mass);

    if (fGoodMuonsID) {
      fhMassWithMuonCuts[modifier[0]][fChan]->Fill(mass);
      if (fGoodHLT) {
	fhMassWithTriggerCuts[modifier[0]][fChan]->Fill(mass);
	fhMassWithAllCuts[modifier[0]][fChan]->Fill(mass);
	if (fIsCowboy) {
	  fhMassWithAllCutsCowboy[modifier[0]][fChan]->Fill(mass);
	} else {
	  fhMassWithAllCutsSeagull[modifier[0]][fChan]->Fill(mass);
	}

	// -- blind version (to have the possibility to make a blind plot also after unblinding)
	if ((5.2 < mass) && (mass < 5.45)) {
	} else {
	  fhMassWithAllCutsBlind[modifier[0]][fChan]->Fill(mass);
	  if (fIsCowboy) {
	    fhMassWithAllCutsCowboyBlind[modifier[0]][fChan]->Fill(mass);
	  } else {
	    fhMassWithAllCutsSeagullBlind[modifier[0]][fChan]->Fill(mass);
	  }
	}

	// -- weighted with fake rate
	fhW8MassWithAllCuts[modifier[0]][fChan]->Fill(mass, fW8MisId);
	if (fIsCowboy) {
	  fhW8MassWithAllCutsCowboy[modifier[0]][fChan]->Fill(mass, fW8MisId);
	} else {
	  fhW8MassWithAllCutsSeagull[modifier[0]][fChan]->Fill(mass, fW8MisId);
	}

	if (fYear < 2016) {
	  if (fb.ps == 0) fb.ps = 1;
	}

	// - include prescale values on y axis
	if ((fMode == BS2JPSIPHI)
	    || (fMode == BD2JPSIKSTAR)
	    || (fMode == BU2JPSIKP)) {
	  fhNorm[modifier[0]][fChan]->Fill(mass, -0.1, static_cast<double>(fb.ps));
	  fhNorm[modifier[0]][fChan]->Fill(mass, 0.1);
	  fhNorm[modifier[0]][fChan]->Fill(mass, fb.ps+0.1);
	  fhNormC[modifier[0]][fChan]->Fill(fb.cm, -0.1, static_cast<double>(fb.ps));
	  fhNormC[modifier[0]][fChan]->Fill(fb.cm, 0.1);
	  fhNormC[modifier[0]][fChan]->Fill(fb.cm, fb.ps+0.1);
	}


	if (fMode == BSMM && fCuts[fChan]->mBsLo < mass && mass < fCuts[fChan]->mBsHi) {
	  fhMassWithMassCuts[modifier[0]][fChan]->Fill(mass);
	}

	if (fMode == BDMM && fCuts[fChan]->mBdLo < mass && mass < fCuts[fChan]->mBdHi) {
	  fhMassWithMassCuts[modifier[0]][fChan]->Fill(mass);
	}

	if (fMode == BU2JPSIKP && fNoLo < mass && mass < fNoHi) {
	  fhMassWithMassCuts[modifier[0]][fChan]->Fill(mass);
	}

	if (fMode == BS2JPSIPHI && fCsLo < mass && mass < fCsHi) {
	  fhMassWithMassCuts[modifier[0]][fChan]->Fill(mass);
	}

	if (fMode == BD2JPSIKSTAR && fNoLo < mass && mass < fNoHi) {
	  fhMassWithMassCuts[modifier[0]][fChan]->Fill(mass);
	}

      }
    }
  }

  // -----------------
  // -- BDT histograms
  // -----------------
  if (fGoodQ
      && fGoodPvAveW8
      && fGoodTracks
      && fGoodTracksPt
      && fGoodTracksEta
      && fGoodMuonsPt  // PidTables do not really work below 4 GeV!!
      && fGoodMuonsEta
      && fGoodJpsiCuts
      && fGoodBDT
      ) {
    fhMassWithAnaCuts[modifier[1]][fChan]->Fill(mass);

    if (fGoodMuonsID) {
      fhMassWithMuonCuts[modifier[1]][fChan]->Fill(mass);
      if (fGoodHLT) {
	fhMassWithTriggerCuts[modifier[1]][fChan]->Fill(mass);
	fhMassWithAllCuts[modifier[1]][fChan]->Fill(mass);
	if (fIsCowboy) {
	  fhMassWithAllCutsCowboy[modifier[1]][fChan]->Fill(mass);
	} else {
	  fhMassWithAllCutsSeagull[modifier[1]][fChan]->Fill(mass);
	}
	// -- blind version
	if ((5.2 < mass) && (mass < 5.45)) {
	} else {
	  fhMassWithAllCutsBlind[modifier[1]][fChan]->Fill(mass);
	  if (fIsCowboy) {
	    fhMassWithAllCutsCowboyBlind[modifier[0]][fChan]->Fill(mass);
	  } else {
	    fhMassWithAllCutsSeagullBlind[modifier[0]][fChan]->Fill(mass);
	  }
	}

	// -- weighted with fake rate
	if (0) cout << "fW8MisId = " << fW8MisId
		    << " m1eta =  " << fb.m1eta << " m2eta =  " << fb.m2eta
		    << " m1pt =  " << fb.m1pt << " m2pt =  " << fb.m2pt
		    << " g1id =  " << fb.g1id << " g2id =  " << fb.g2id
		    << endl;

	fhW8MassWithAllCuts[modifier[1]][fChan]->Fill(mass, fW8MisId);
	if (fIsCowboy) {
	  fhW8MassWithAllCutsCowboy[modifier[1]][fChan]->Fill(mass, fW8MisId);
	} else {
	  fhW8MassWithAllCutsSeagull[modifier[1]][fChan]->Fill(mass, fW8MisId);
	}

	// - include prescale values on y axis
	if ((fMode == BS2JPSIPHI)
	    || (fMode == BD2JPSIKSTAR)
	    || (fMode == BU2JPSIKP)) {
	  fhNorm[modifier[1]][fChan]->Fill(mass, -0.1, static_cast<double>(fb.ps));
	  fhNorm[modifier[1]][fChan]->Fill(mass, 0.1);
	  fhNorm[modifier[1]][fChan]->Fill(mass, fb.ps+0.1);
	  fhNormC[modifier[1]][fChan]->Fill(fb.cm, -0.1, static_cast<double>(fb.ps));
	  fhNormC[modifier[1]][fChan]->Fill(fb.cm, 0.1);
	  fhNormC[modifier[1]][fChan]->Fill(fb.cm, fb.ps+0.1);
	}

	if (fMode == BSMM && fCuts[fChan]->mBsLo < mass && mass < fCuts[fChan]->mBsHi) {
	  fhMassWithMassCuts[modifier[1]][fChan]->Fill(mass);
	}

	if (fMode == BDMM && fCuts[fChan]->mBdLo < mass && mass < fCuts[fChan]->mBdHi) {
	  fhMassWithMassCuts[modifier[1]][fChan]->Fill(mass);
	}

	if (fMode == BU2JPSIKP && fNoLo < mass && mass < fNoHi) {
	  fhMassWithMassCuts[modifier[1]][fChan]->Fill(mass);
	}

	if (fMode == BS2JPSIPHI && fCsLo < mass && mass < fCsHi) {
	  fhMassWithMassCuts[modifier[1]][fChan]->Fill(mass);
	}

	if (fMode == BD2JPSIKSTAR && fNoLo < mass && mass < fNoHi) {
	  fhMassWithMassCuts[modifier[1]][fChan]->Fill(mass);
	}
      }
    }
  }

}

// ----------------------------------------------------------------------
void plotResults::loopFunction2() {

}


// ----------------------------------------------------------------------
void plotResults::loopFunction3() {


}



// ----------------------------------------------------------------------
void plotResults::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotResults::loopOverTree> loop over dataset " << fCds->fName << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotResults::*pF)(void);
  if (ifunc == 1) pF = &plotResults::loopFunction1;
  if (ifunc == 2) pF = &plotResults::loopFunction2;
  if (ifunc == 3) pF = &plotResults::loopFunction3;

  // -- setup small tree
  TDirectory *dir(0);
  TTree *small(0);
  TFile *fLocal(0);
  if (fSaveSmallTree) {
    string tname(fSetup);
    replaceAll(tname, "Off", "");
    replaceAll(tname, "Comb", "");

    dir = gDirectory;
    fLocal = TFile::Open(Form("%s/small%d-%s.root", fDirectory.c_str(), fYear, tname.c_str()), "RECREATE");
    small = new TTree(Form("%s", tname.c_str()), Form("%s", tname.c_str()));
    small->SetDirectory(fLocal);
    small->Branch("run",    &fb.run,       "run/I");
    small->Branch("evt",    &fb.evt,       "evt/I");
    small->Branch("ls",     &fb.ls,        "ls/I");
    small->Branch("ps",     &fb.ps,        "ps/I");

    small->Branch("chan",   &fb.chan,      "chan/I");
    small->Branch("muid",   &fb.gmuid,     "muid/O");
    small->Branch("bdt",    &fBDT,         "bdt/D");
    small->Branch("cnc",    &fGoodCNC,     "cnc/O");
    small->Branch("pt",     &fb.pt,        "pt/D");
    small->Branch("eta",    &fb.eta,       "eta/D");
    small->Branch("phi",    &fb.phi,       "phi/D");
    small->Branch("m",      &fb.m,         "m/D");
    small->Branch("me",     &fb.me,        "me/D");

    small->Branch("tau",    &fb.tau,       "tau/D");
    small->Branch("taue",   &fb.taue,      "taue/D");
    small->Branch("gtau",   &fb.gtau,      "gtau/D");

    small->Branch("m1pt",   &fb.m1pt,      "m1pt/D");
    small->Branch("m1phi",  &fb.m1phi,     "m1phi/D");
    small->Branch("m1eta",  &fb.m1eta,     "m1eta/D");
    small->Branch("m1bdt",  &fb.m1mvabdt,  "m1bdt/D");
    small->Branch("m1rbdt", &fb.m1rmvabdt, "m1rbdt/D");
    small->Branch("m1q",    &fb.m1q,       "m1q/I");

    small->Branch("m2pt",   &fb.m2pt,      "m2pt/D");
    small->Branch("m2phi",  &fb.m2phi,     "m2phi/D");
    small->Branch("m2eta",  &fb.m2eta,     "m2eta/D");
    small->Branch("m2bdt",  &fb.m2mvabdt,  "m2bdt/D");
    small->Branch("m2rbdt", &fb.m2rmvabdt, "m2rbdt/D");
    small->Branch("m2q",    &fb.m2q,       "m2q/I");
  }

  // ----------------------------
  // -- the real loop starts here
  // ----------------------------
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();

    // -- filling of small tree
    if (fSaveSmallTree
	&& fGoodHLT && fGoodGlobalMuons
	&& (fGoodCNC || fBDT > -2.)
	) {
      small->Fill();
    }
  }

  if (fSaveSmallTree) {
    small->Write();
    fLocal->Close();
    dir->cd();
  }

}


// ----------------------------------------------------------------------
void plotResults::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotResults::loadFile loading files listed in " << files << endl;

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

    if (string::npos != stype.find("Data")) {
      // -- Charmonium
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("blablabla")) {
        sname = "bmmCharmonium";
        sdecay = "bmm";
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
	ds->fBf     = -1.;
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


// ----------------------------------------------------------------------
void plotResults::saveHistograms(string smode) {

  fHistFile->cd();
  fHistFile->mkdir(fSetup.c_str());
  fHistFile->cd(fSetup.c_str());
  TDirectory *dir = gDirectory;

  TH1D *h1(0);
  TH2D *h2(0);

  vector<string> modifier;
  modifier.push_back("cnc" + fSuffix);
  modifier.push_back("bdt" + fSuffix);

  for (unsigned int im = 0; im < modifier.size(); ++im) {
    for (unsigned int i = 0; i < fNchan; ++i) {
      h1 = (TH1D*)(fhGenAndAccNumbers[modifier[im]][i]->Clone(Form("hGenAndAccNumbers_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();

      h2 = (TH2D*)(fhAccAll[modifier[im]][i]->Clone(Form("hAccAll_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h2->SetDirectory(dir);
      h2->Write();
      h2 = (TH2D*)(fhAccPass[modifier[im]][i]->Clone(Form("hAccPass_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h2->SetDirectory(dir);
      h2->Write();

      h1 = (TH1D*)(fhAccEtaAll[modifier[im]][i]->Clone(Form("hAccEtaAll_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();
      h1 = (TH1D*)(fhAccEtaPass[modifier[im]][i]->Clone(Form("hAccEtaPass_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhAccPtAll[modifier[im]][i]->Clone(Form("hAccPtAll_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();
      h1 = (TH1D*)(fhAccPtPass[modifier[im]][i]->Clone(Form("hAccPtPass_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassAbsNoCuts[modifier[im]][i]->Clone(Form("hMassAbsNoCuts_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassAbsNoCuts_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassNoCuts[modifier[im]][i]->Clone(Form("hMassNoCuts_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassNoCuts_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAnaCuts[modifier[im]][i]->Clone(Form("hMassWithAnaCuts_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAnaCuts_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithMuonCuts[modifier[im]][i]->Clone(Form("hMassWithMuonCuts_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithMuonCuts_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithTriggerCuts[modifier[im]][i]->Clone(Form("hMassWithTriggerCuts_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithTriggerCuts_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCuts[modifier[im]][i]->Clone(Form("hMassWithAllCuts_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCuts_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsBlind[modifier[im]][i]->Clone(Form("hMassWithAllCutsBlind_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsBlind_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsSeagull[modifier[im]][i]->Clone(Form("hMassWithAllCutsSeagull_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsSeagull_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsSeagullBlind[modifier[im]][i]->Clone(Form("hMassWithAllCutsSeagullBlind_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsSeagullBlind_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsCowboy[modifier[im]][i]->Clone(Form("hMassWithAllCutsCowboy_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsCowboy_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsCowboyBlind[modifier[im]][i]->Clone(Form("hMassWithAllCutsCowboyBlind_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsCowboyBlind_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhW8MassWithAllCuts[modifier[im]][i]->Clone(Form("hW8MassWithAllCuts_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hW8MassWithAllCuts_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();


      h1 = (TH1D*)(fhW8MassWithAllCutsSeagull[modifier[im]][i]->Clone(Form("hW8MassWithAllCutsSeagull_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hW8MassWithAllCutsSeagull_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhW8MassWithAllCutsCowboy[modifier[im]][i]->Clone(Form("hW8MassWithAllCutsCowboy_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hW8MassWithAllCutsCowboy_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();


      h1 = (TH1D*)(fhMassWithMassCuts[modifier[im]][i]->Clone(Form("hMassWithMassCuts_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithMassCuts_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();



      if ((string::npos != fSetup.find("bupsik"))
	  || (string::npos != fSetup.find("bspsiphi"))
	  || (string::npos != fSetup.find("bdpsikstar"))
	  ) {
	h2 = (TH2D*)(fhNorm[modifier[im]][i]->Clone(Form("hNorm_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
	h2->SetTitle(Form("hNorm_%s_%s_chan%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
	h2->SetDirectory(dir);
	h2->Write();

	h2 = (TH2D*)(fhNormC[modifier[im]][i]->Clone(Form("hNormC_%s_%s_chan%d", modifier[im].c_str(), smode.c_str(), i)));
	h2->SetTitle(Form("hNormC_%s_%s_%d %s", modifier[im].c_str(), smode.c_str(), i, smode.c_str()));
	h2->SetDirectory(dir);
	h2->Write();
      }
    }
  }

  TH1D *hcuts = fDS["fakeMc"]->getHist("candAnaFakeMC/hcuts");
  hcuts->SetName("hcuts");
  hcuts->SetDirectory(fHistFile);
  cout << "writing hcuts " << hcuts << " name: " << hcuts->GetName() << " to " << fHistFile->GetName() << endl;
  hcuts->Write();

}

// ----------------------------------------------------------------------
void plotResults::resetHistograms(bool deleteThem) {

  vector<string> modifier;
  modifier.push_back("cnc" + fSuffix);
  modifier.push_back("bdt" + fSuffix);

  for (unsigned int im = 0; im < modifier.size(); ++im) {
    for (int i = 0; i < fNchan; ++i) {
      fhAccAll[modifier[im]][i]->Reset();
      if (deleteThem) delete fhAccAll[modifier[im]][i];
      fhAccPass[modifier[im]][i]->Reset();
      if (deleteThem) delete fhAccPass[modifier[im]][i];

      fhAccPtAll[modifier[im]][i]->Reset();
      if (deleteThem) delete fhAccPtAll[modifier[im]][i];
      fhAccPtPass[modifier[im]][i]->Reset();
      if (deleteThem) delete fhAccPtPass[modifier[im]][i];

      fhAccEtaAll[modifier[im]][i]->Reset();
      if (deleteThem) delete fhAccEtaAll[modifier[im]][i];
      fhAccEtaPass[modifier[im]][i]->Reset();
      if (deleteThem) delete fhAccEtaPass[modifier[im]][i];

      fhGenAndAccNumbers[modifier[im]][i]->Reset();
      if (deleteThem) delete fhGenAndAccNumbers[modifier[im]][i];

      fhMassAbsNoCuts[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassAbsNoCuts[modifier[im]][i];
      fhMassNoCuts[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassNoCuts[modifier[im]][i];

      fhMassWithAnaCuts[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAnaCuts[modifier[im]][i];

      fhMassWithMuonCuts[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithMuonCuts[modifier[im]][i];

      fhMassWithTriggerCuts[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithTriggerCuts[modifier[im]][i];

      fhMassWithAllCuts[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCuts[modifier[im]][i];
      fhMassWithAllCutsBlind[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsBlind[modifier[im]][i];

      fhMassWithAllCutsSeagull[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsSeagull[modifier[im]][i];
      fhMassWithAllCutsSeagullBlind[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsSeagullBlind[modifier[im]][i];

      fhMassWithAllCutsCowboy[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsCowboy[modifier[im]][i];
      fhMassWithAllCutsCowboyBlind[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsCowboyBlind[modifier[im]][i];

      fhW8MassWithAllCuts[modifier[im]][i]->Reset();
      if (deleteThem) delete fhW8MassWithAllCuts[modifier[im]][i];

      fhW8MassWithAllCutsSeagull[modifier[im]][i]->Reset();
      if (deleteThem) delete fhW8MassWithAllCutsSeagull[modifier[im]][i];

      fhW8MassWithAllCutsCowboy[modifier[im]][i]->Reset();
      if (deleteThem) delete fhW8MassWithAllCutsCowboy[modifier[im]][i];

      fhMassWithMassCuts[modifier[im]][i]->Reset();
      if (deleteThem) delete fhMassWithMassCuts[modifier[im]][i];

      fhNorm[modifier[im]][i]->Reset();
      if (deleteThem) delete fhNorm[modifier[im]][i];
      fhNormC[modifier[im]][i] ->Reset();
      if (deleteThem) delete fhNormC[modifier[im]][i];
    }
  }
}


// ----------------------------------------------------------------------
void plotResults::dumpTex(number &a, const string &s, int ndigits, int sgn) {
  fTEX << formatTex(a.val,   Form("%s:val", s.c_str()), ndigits) << endl;
  fTEX << formatTex(a.estat, Form("%s:estat", s.c_str()), ndigits) << endl;
  fTEX << formatTex(a.esyst, Form("%s:esyst", s.c_str()), ndigits) << endl;
  fTEX << formatTex(a.etot,  Form("%s:etot", s.c_str()), ndigits) << endl;
  fTEX << formatTexErrSci(a.val, a.etot, Form("%s:all", s.c_str()), ndigits, sgn) << endl;
}


// ----------------------------------------------------------------------
void plotResults::printNumbers(anaNumbers &a, ostream &OUT) {
  number n;
  OUT << "==> anaNumbers for " << a.fName << " chan " << a.fChan << endl;
  n = a.fAccGenFileYield; OUT << "    fAccGenFileYield = " << Form("%6.0f +/- %6.1f +/- %6.1f", n.val,  n.estat, n.esyst) << endl;
  n = a.fAccGenYield;     OUT << "    fAccGenYield =     " << Form("%6.0f +/- %6.1f +/- %6.1f", n.val,  n.estat, n.esyst) << endl;
  n = a.fGenFileYield;    OUT << "    fGenFileYield =    " << Form("%6.0f +/- %6.1f +/- %6.1f", n.val,  n.estat, n.esyst) << endl;
  n = a.fGenYield;        OUT << "    fGenYield =        " << Form("%6.0f +/- %6.1f +/- %6.1f", n.val,  n.estat, n.esyst) << endl;
  n = a.fEffGenSel;       OUT << "    fEffGenSel =       " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fEffFilter;       OUT << "    fEffFilter =       " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fAccRecoYield;    OUT << "    fAccRecoYield =    " << Form("%7.0f +/- %7.0f +/- %7.0f", n.val,  n.estat, n.esyst) << endl;
  n = a.fAccCandYield;    OUT << "    fAccCandYield =    " << Form("%7.0f +/- %7.0f +/- %7.0f", n.val,  n.estat, n.esyst) << endl;
  n = a.fCandYield;       OUT << "    fCandYield =       " << Form("%7.0f +/- %7.0f +/- %7.0f", n.val,  n.estat, n.esyst) << endl;
  n = a.fAnaYield;        OUT << "    fAnaYield =        " << Form("%7.0f +/- %7.0f +/- %7.0f", n.val,  n.estat, n.esyst) << endl;
  n = a.fMuidYield;       OUT << "    fMuidYield =       " << Form("%7.0f +/- %7.0f +/- %7.0f", n.val,  n.estat, n.esyst) << endl;
  n = a.fTrigYield;       OUT << "    fTrigYield =       " << Form("%7.0f +/- %7.0f +/- %7.0f", n.val,  n.estat, n.esyst) << endl;
  n = a.fAcc;             OUT << "    fAcc =             " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fEffCand;         OUT << "    fEffCand =         " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fEffAna;          OUT << "    fEffAna =          " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fEffMuidMC;       OUT << "    fEffMuidMC =       " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fEffTrigMC;       OUT << "    fEffTrigMC =       " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fEffTot;          OUT << "    fEffTot =          " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fEffProdMC;       OUT << "    fEffProdMC =       " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
  n = a.fSignalFit;       OUT << "    fSignalFit =       " << Form("%7.0f +/- %7.0f +/- %7.0f", n.val,  n.estat, n.esyst) << endl;
}

// ----------------------------------------------------------------------
bool plotResults::skipThisBg(string name) {
  if (string::npos == name.find("Bg")) return true;
  if (string::npos != name.find("Acc")) return true;
  if (string::npos != name.find("bdpimunuMcOffBg")) return true; // use the private sample instead
  if (string::npos != name.find("bcpsimunu")) return true;
  if (string::npos != name.find("McOffBg") || string::npos != name.find("McBg")) {
    string tname = name;
    replaceAll(tname, "McBg", "McCombBg");
    if ((tname != name) && fDS.count(tname) > 0) {
      cout << "instead of " << name << " there is also " << tname << ", which will be used" << endl;
      return true;
    }
    tname = name;
    replaceAll(tname, "McOffBg", "McCombBg");
    if ((tname != name) && fDS.count(tname) > 0) {
      cout << "instead of " << name << " there is also " << tname << ", which will be used" << endl;
      return true;
    }
  }
  return false;
}


// ----------------------------------------------------------------------
string plotResults::rareAccName(string sname) {
  string accname("");
  // -- the acceptance sample cannot be determined trivially
  if (string::npos != sname.find("bdpimumuMcCombBg")) {
    accname = "bdpimumuMcCombAccBg";
    return accname;
  }
  if (string::npos != sname.find("bupimumuMcCombBg")) {
    accname = "bupimumuMcCombAccBg";
    return accname;
  }

  if (string::npos != sname.find("bs")) {
    if (string::npos == sname.find("mu")) {
      accname = "bskkMcCombAccBg";
    } else {
      accname = "bskmunuMcCombAccBg";
    }
  }
  if (string::npos != sname.find("bd")) {
    if (string::npos == sname.find("mu")) {
      accname = "bdkpiMcCombAccBg";
    } else {
      accname = "bdpimunuMcCombAccBg";
    }
  }
  if (string::npos != sname.find("lb")) {
    if (string::npos == sname.find("mu")) {
      accname = "lbpkMcCombAccBg";
    } else {
      accname = "lbpmunuMcCombAccBg";
    }
  }
  return accname;
}


// ----------------------------------------------------------------------
double plotResults::massIntegral(TH1* h, INTMODE imode, int ichan) {
  int lo(0), hi(0);
  double eps(1.e-6);
  if (imode == LO) {
    lo = h->FindBin(fMassLo + eps);
    hi = h->FindBin(fCuts[ichan]->mBdLo - eps);
  }
  if (imode == BD) {
    lo = h->FindBin(fCuts[ichan]->mBdLo + eps);
    hi = h->FindBin(fCuts[ichan]->mBdHi - eps);
  }
  if (imode == BS) {
    lo = h->FindBin(fCuts[ichan]->mBsLo + eps);
    hi = h->FindBin(fCuts[ichan]->mBsHi - eps);
  }
  if (imode == HI) {
    lo = h->FindBin(fCuts[ichan]->mBsHi + eps);
    hi = h->FindBin(fMassHi - eps);
  }
  if (imode == ALL) {
    lo = 0;
    hi = h->GetNbinsX()+1;
  }

  return h->Integral(lo, hi);
}
