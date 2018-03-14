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
#include "TMarker.h"
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
plotResults::plotResults(string dir, string files, string cuts, string setup, int year): plotClass(dir, files, cuts, setup, year),
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

  fHistStrings.clear();
  fHistStrings.push_back("cnc" + fSuffix);
  fHistStrings.push_back("bdt" + fSuffix);
  for (int i = 0; i <= 80; ++i) {
    fHistStrings.push_back(Form("bdt_%d_", i) + fSuffix);
  }

  TH2D *h2(0);
  TH1D *h(0);
  string mode("cnc");
  for (int i = 0; i < fHistStrings.size(); ++i) {
    mode = fHistStrings[i];
    for (int i = 0; i < fNchan; ++i) {
      h2 = new TH2D(Form("h%sAccAll%d", mode.c_str(), i), Form("h%sAccAll%d", mode.c_str(), i), 25, 0., 2.5, 25, 0., 50.);
      fhAccAll[mode].push_back(h2);
      h2 = new TH2D(Form("h%sAccPass%d", mode.c_str(), i), Form("h%sAccPass%d", mode.c_str(), i), 25, 0., 2.5, 25, 0., 50.);
      fhAccPass[mode].push_back(h2);

      h = new TH1D(Form("h%sBdtCrossCheck%d", mode.c_str(), i), Form("h%sBdtCrossCheck%d", mode.c_str(), i), 200, -1., 1.);
      fhBdtCrossCheck[mode].push_back(h);

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

      h2 = new TH2D(Form("h%sW8Norm%d", mode.c_str(), i), Form("h%sW8Norm%d", mode.c_str(), i), 100, 4.9, 5.9, MAXPS+1, -1., MAXPS);
      fhW8Norm[mode].push_back(h2);
      h2 = new TH2D(Form("h%sW8NormC%d", mode.c_str(), i), Form("h%sW8NormC%d", mode.c_str(), i), 200, 4.9, 5.9, MAXPS+1, -1., MAXPS);
      fhW8NormC[mode].push_back(h2);
    }
  }

  for (int i = 0; i < fNchan; ++i) {
    h = new TH1D(Form("h%sBdt%d", mode.c_str(), i), Form("h%sBdt%d", mode.c_str(), i), 200, -1., 1.);
    fhBdt.push_back(h);
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
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotResults::makeAll(string what) {

  if (what == "bdtopt") {
    fHistWithAllCuts = "hMassWithAllCuts";
    //    fillAndSaveHistograms(0, 10000);
    fillAndSaveHistograms();
    fSuffixSel = "bdt" + fSuffix;
    for (int i = 0; i < fNchan; ++i) {
      calculateNumbers("bdt" + fSuffix, i);
    }
    scanBDT(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), true);
    return;
  }

  if (what == "bdtoptplot") {
    scanBDT(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), false);
    return;
  }

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
    } else if ((2012 == fYear) || (2011 == fYear)) {
      genSummary("bupsikMcComb", "candAnaBu2JpsiK");
      genSummary("bupsikMcCombAcc", "candAnaBu2JpsiK");

      genSummary("bspsiphiMcComb", "candAnaBs2JpsiPhi");
      genSummary("bspsiphiMcCombAcc", "candAnaBs2JpsiPhi");

      genSummary("bdpsikstarMcComb", "candAnaBd2JpsiKstar");
      genSummary("bdpsikstarMcCombAcc", "candAnaBd2JpsiKstar");

      genSummary("bdmmMcComb", "candAnaMuMu");
      genSummary("bdmmMcCombAcc", "candAnaMuMu");

      genSummary("bsmmMcComb", "candAnaMuMu");
      genSummary("bsmmMcCombAcc", "candAnaMuMu");

      // -- loop over all (effective) two-body backgrounds
      for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
	if (string::npos == it->first.find("Bg")) continue;
	cout << "===> Create genSummary for " << it->first << endl;
	genSummary(it->first, "candAnaMuMu");
      }
    }

  }

  // -- this will recreate fHistFile!
  if (what == "small") {
    string mode = "bupsikMcComb";
    fHistFile = TFile::Open("bla.root", "RECREATE");
    fSaveSmallTree = true;
    resetHistograms();
    setup(mode);
    TTree *t = getTree(mode, fTreeDir);
    setupTree(t, mode);
    loopOverTree(t, 1);
    saveHistograms(mode);
    fHistFile->Close();
    return;
  }

  // -- this will recreate fHistFile!
  if ((what == "all") || (string::npos != what.find("fill"))  || (string::npos != what.find("ana"))) {
    //    fillAndSaveHistograms(0, 50000);
    fillAndSaveHistograms();
  }

  if ((what == "all") || (string::npos != what.find("ana"))) {
    dumpDatasets();
    fHistWithAllCuts = "hMassWithAllCuts";
    fSuffixSel = "cnc" + fSuffix;
    for (int i = 0; i < fNchan; ++i) {
      calculateNumbers("cnc" + fSuffix, i);
      fSuffixSel = "bdt" + fSuffix;
      calculateNumbers("bdt" + fSuffix, i);
    }
    scanBDT(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), true);
  }

  if (what == "all" || string::npos != what.find("cnc")) {
    fHistWithAllCuts = "hMassWithAllCuts";
    fSuffixSel = "cnc" + fSuffix;
    for (int i = 0; i < fNchan; ++i) {
      calculateNumbers("cnc" + fSuffix, i);
    }
  }

  if (what == "all" || string::npos != what.find("bdt")) {
    fHistWithAllCuts = "hMassWithAllCuts";
    fSuffixSel = "bdt" + fSuffix;
    for (int i = 0; i < fNchan; ++i) {
      calculateNumbers("bdt" + fSuffix, i);
    }
    scanBDT(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), true);
  }

  if (what == "dbx2") {

    frd.pt = 16.016298;
    frd.eta = -0.337338;
    frd.fls3d = 5.8128776;
    frd.alpha = 0.1133278;
    frd.maxdoca = 0.001;
    frd.pvip = 0.0034266;
    frd.pvips = 1.0469902;
    frd.iso = 0.9435766;
    frd.docatrk = 0.0077934;
    frd.chi2dof = 1.335e-06;
    frd.closetrk = 10;

    frd.fl3d = 0.0321190;
    frd.flsxy = 6.6656122;
    frd.m1iso = 0.8524752;
    frd.m2iso = 0.7810096;

    frd.closetrks1 = 0.;
    frd.closetrks2 = 0.;
    frd.closetrks3 = 0.;

    frd.pv2lip  = 0.;
    frd.pv2lips = 0.;

    frd.m  = fb.m;
    cout << fPresel.preselection() << endl;
    cout << fReaderEvents0[0]->EvaluateMVA("BDT") << endl;
    cout << fReaderEvents1[0]->EvaluateMVA("BDT") << endl;
    cout << fReaderEvents2[0]->EvaluateMVA("BDT") << endl;
  }



}



// ----------------------------------------------------------------------
void plotResults::bookHist(string dsname) {
  fHistWithAllCuts = "hMassWithAllCuts";
  fSuffixSel = "bdt" + fSuffix;
  for (int i = 0; i < fNchan; ++i) {
    calculateNumbers("bdt" + fSuffix, i);
  }
}


// ----------------------------------------------------------------------
void plotResults::scanBDT(string fname, bool createTexFile) {
  fTEX.close();
  cout << "starting scanBDT, createTexFile = " << createTexFile << endl;
  int BDTMIN(0);
  if (createTexFile) {
    vector<int> BDTMAX;
    for (int ichan = 0; ichan < fNchan; ++ichan) {BDTMAX.push_back(-1);}
    cout << "creating tex file" << endl;
    string dsname = Form("%s", fname.c_str());
    fHistWithAllCuts = "hMassWithAllCuts";
    // -- check over which range you have to run. Find maximum PER CHANNEL
    fHistFile = TFile::Open(fHistFileName.c_str());
    for (int ichan = 0; ichan < fNchan; ++ichan) {
      TH1D *h0 = (TH1D*)fHistFile->Get(Form("bdmmMcComb/hMassWithAllCuts_bdt_%d_%s_bdmmMcComb_chan%d", 0, fSuffix.c_str(), ichan));
      if (!h0) {
	cout << "histogram " << Form("bdmmMcComb/hMassWithAllCuts_bdt_%d_%s_bdmmMcComb_chan%d", 0, fSuffix.c_str(), ichan)
	     << " not found" << endl;
	return;
      }
      double total = h0->GetSumOfWeights();
      cout << "total events " << total << " for " << h0->GetName() << endl;
      for (int ib = BDTMIN; ib <= 80; ++ib) {
	h0 = (TH1D*)fHistFile->Get(Form("bdmmMcComb/hMassWithAllCuts_bdt_%d_%s_bdmmMcComb_chan%d", ib, fSuffix.c_str(), ichan));
	double integral = h0->GetSumOfWeights();
	cout << "   integral(" << ib << ") = " << integral << endl;
	if (integral < 0.05*total) {
	  BDTMAX[ichan] = ib;
	  cout << " chan " << ichan << " going up to bdt < " << BDTMAX[ichan]
	       << ", because there events = " << h0->GetSumOfWeights()
	       << endl;
	  break;
	}
      }
    }
    fHistFile->Close();

    // -- and now loop over all requested
    fTEX.open(dsname.c_str());
    for (int ichan = 0; ichan < 2; ++ichan) {
      cout << "producing numbers for chan = " << ichan << " and BDTMAX = " << BDTMAX[ichan] << endl;
      for (int ib = BDTMIN; ib <= BDTMAX[ichan]; ++ib) {
	string idx = Form("bdt_%d_", ib);
	fSuffixSel = idx + fSuffix;
	cout << "  bdt > " << ib*0.01 << " fSuffixSel = " << fSuffixSel << endl;
	calculateNumbers(idx + fSuffix, ichan);
      }
    }
    fTEX.close();
  }

  vector<string> plots;
  plots.push_back("CSBF");
  plots.push_back("BSMMBFS");
  plots.push_back("BSMMBFU");
  plots.push_back("ZAD");
  plots.push_back("ZAS");
  plots.push_back("SSB");
  plots.push_back("SOB");
  plots.push_back("SOBS");
  plots.push_back("SOBB");
  plots.push_back("alphaS");
  plots.push_back("alphaU");

  fHistFile = TFile::Open(fHistFileName.c_str());
  vector<TH1D*> hbdt;
  TIter next(fHistFile->GetListOfKeys());
  TKey *key(0);
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;
    fHistFile->cd(key->GetName());
    cout << gDirectory->GetName() << endl;
    for (int ic = 0; ic < fNchan; ++ic) {
      hbdt.push_back((TH1D*)((TH1D*)gDirectory->Get(Form("hBdt_%s_chan%d", gDirectory->GetName(), ic)))->Clone());
    }
  }

  string histfilename = Form("%s", fname.c_str());
  replaceAll(histfilename, ".tex", ".root");
  cout << "fHistFile: " << histfilename;
  TFile *HistFile = TFile::Open(histfilename.c_str(), "RECREATE");
  cout << " opened " << endl;

  int iMax(-1), bMax(-1);
  for (unsigned int i = 0; i < hbdt.size(); ++i) {
    string sname = hbdt[i]->GetName();
    if (string::npos != sname.find("bdmm") && (string::npos != sname.find("chan0") || string::npos != sname.find("chan1"))) {
      double tot = hbdt[i]->Integral();
      double rs = 0.;
      cout << "hist " << hbdt[i]->GetName() << " tot = " << tot << endl;
      for (int ibin = 1; ibin <= hbdt[i]->GetNbinsX(); ++ibin) {
	rs += hbdt[i]->GetBinContent(ibin);
	cout << "ibin = " << ibin
	     << " rs = " << rs << " tot = " << tot << " rs/tot = " << rs/tot
	     << " iMax = " << iMax
	     << endl;
	if (rs/tot > 0.95) {
	  if (ibin > iMax) {
	    iMax = ibin;
	    bMax = 100 * hbdt[i]->GetBinLowEdge(ibin);
	    cout << "STOP tot =  " << tot << " rs = " << rs << " iMax = " << iMax  << " bMax = " << bMax
		 << " absolute BDT < " << hbdt[i]->GetBinLowEdge(ibin)
		 << endl;
	  }
	  break;
	}
      }
    }
    // cout << "write out " << hbdt[i]->GetName() << endl;
    hbdt[i]->SetDirectory(HistFile);
    hbdt[i]->Write();
    delete hbdt[i];
  }

  cout << " go up to bMax = " << bMax << endl;

  TH1D *h1 = new TH1D("bdtCuts", "bdtCuts", 100, 0., 100.);
  for (int i = 0; i < fNchan; ++i) {
    h1->SetBinContent(i+1, fCuts[i]->bdtCut);
    h1->GetXaxis()->SetBinLabel(i+1, Form("chan%d", i));
  }
  h1->SetDirectory(HistFile);
  h1->Write();

  // -- now read in tex file
  vector<string> allLines;
  ifstream is(fname);
  char buffer[2000];
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  is.close();
  string name;
  for (unsigned int i = 0; i < plots.size(); ++i) {
    for (unsigned int ic = 0; ic < fNchan; ++ic) {
      int imax(-1);
      double maximum(-1.);
      h1 = new TH1D(Form("bdtScan_%s_chan%d", plots[i].c_str(), ic),
		    Form("bdtScan_%s_chan%d", plots[i].c_str(), ic),
		    80, 0., 80.);
      for (unsigned int ib = BDTMIN; ib <= bMax; ++ib) {
	string idx = Form("bdt_%d_", ib);
	fSuffixSel = idx + fSuffix;
	name = Form("%s:%s:chan%d:val", fSuffixSel.c_str(), plots[i].c_str(), ic);
	double val = findVarValue(name, allLines);
	if (val > maximum) {
	  maximum = val;
	  imax = i;
	}
	if (0) cout << "  " << name << ": " << val << " or: " << Form("%8.6f", val) << " idx ->" << idx << "<-"
		    << " filling into bin " << h1->FindBin(ib)
		    << endl;
	h1->SetBinContent(h1->FindBin(ib), val);
      }
      if (ic < 2) {
	if (string::npos != plots[i].find("ZAS")) cout << Form("maximum chan%d%s: ", ic, plots[i].c_str()) << maximum << endl;
	if (string::npos != plots[i].find("ZAD")) cout << Form("maximum chan%d%s: ", ic, plots[i].c_str()) << maximum << endl;
	if (string::npos != plots[i].find("SSB")) cout << Form("maximum chan%d%s: ", ic, plots[i].c_str()) << maximum << endl;
	if (string::npos != plots[i].find("SOB")) cout << Form("maximum chan%d%s: ", ic, plots[i].c_str()) << maximum << endl;
      }
      h1->SetDirectory(HistFile);
      h1->Write();
      delete h1;
    }
  }
  HistFile->Close();
  fHistFile->Close();
}


// ----------------------------------------------------------------------
void plotResults::displayScanBDT(string what, int mode, int chan) {

  c0->Clear();

  vector<string> inputFiles;
  vector<Color_t> colors;

  if ("all" == what) {
    displayScanBDT("SSB", mode, chan);
    displayScanBDT("ZAD", mode, chan);
    displayScanBDT("ZAS", mode, chan);
    displayScanBDT("CSBF", mode, chan);
    displayScanBDT("BSMMBFS", mode, chan);
    return;
  }

  if (0 == mode) {
    inputFiles.push_back("2011/scanBDT-2011.root");   colors.push_back(kRed);
    inputFiles.push_back("2012/scanBDT-2012.root");   colors.push_back(kBlack);
    // inputFiles.push_back("results/scanBDT-2016BF.root"); colors.push_back(kGreen+2);
    // inputFiles.push_back("results/scanBDT-2016GH.root"); colors.push_back(kBlue);
    inputFiles.push_back("2016BF-00/scanBDT-2016BF-00.root"); colors.push_back(kGreen+2);
    inputFiles.push_back("2016GH-00/scanBDT-2016GH-00.root"); colors.push_back(kBlue);
    inputFiles.push_back("2016BF-01/scanBDT-2016BF-01.root"); colors.push_back(kYellow+2);
    inputFiles.push_back("2016GH-01/scanBDT-2016GH-01.root"); colors.push_back(kCyan+2);
  } else if (1 == mode) {
    inputFiles.push_back("results/scanBDT-2016BF.root");              colors.push_back(kBlack);
    inputFiles.push_back("results/scanBDT-2016BF-389-23.root");   colors.push_back(kRed);
    inputFiles.push_back("results/scanBDT-2016BF-409-23.root");   colors.push_back(kGreen+1);
    inputFiles.push_back("results/scanBDT-2016BF-429-23.root");   colors.push_back(kGreen+2);
    inputFiles.push_back("results/scanBDT-2016BF-419-23.root");   colors.push_back(kBlue);
  } else if (2 == mode) {
    inputFiles.push_back("results/scanBDT-2016GH.root");    colors.push_back(kBlack);
    inputFiles.push_back("se/scanBDT-2016GH-5559.root");    colors.push_back(kRed);
    inputFiles.push_back("se/scanBDT-2016GH-2329.root");    colors.push_back(kGreen+1);
    inputFiles.push_back("se/scanBDT-2016GH-2489.root");    colors.push_back(kBlue);
  } else if (3 == mode) {
    int i = kRed;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-9619.root");    colors.push_back(i);  i = i -1;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-7639.root");    colors.push_back(i); i = i -1;
    i = kBlue;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-16309.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-2029.root");    colors.push_back(i); i = i -1;
    i = kMagenta;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-9999.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-4459.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-11549.root");    colors.push_back(i);  i = i -1;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-6769.root");    colors.push_back(i); i = i -1;
    i = kGreen;
    inputFiles.push_back("/scratch/ursl/bmm4/se/abdt-4/scanBDT-2016GH-10019.root");    colors.push_back(i); i = i -1;
  } else if (4 == mode) {
    int i = kRed;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-2329.root");    colors.push_back(i);  i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-2489.root");    colors.push_back(i); i = i -1;
    i = kBlue;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-2569.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-2809.root");    colors.push_back(i); i = i -1;
    i = kMagenta;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-2969.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-3289.root");    colors.push_back(i); i = i -1;
    i = kYellow-2;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-3369.root");    colors.push_back(i);  i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-3529.root");    colors.push_back(i); i = i -1;
    i = kGreen;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-4909.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-5209.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-5849.root");    colors.push_back(i); i = i -1;
  } else if (5 == mode) {
    int i = kRed;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-2399.root");    colors.push_back(i);  i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-2799.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-2959.root");    colors.push_back(i); i = i -1;
    i = kBlue;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-3279.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-3439.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-4849.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-5809.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-5839.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-7409.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-8049.root");    colors.push_back(i); i = i -1;
    i = kYellow-2;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-36559.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-37029.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-38159.root");    colors.push_back(i); i = i -1;
    i = kGreen;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-38319.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-38799.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-39119.root");    colors.push_back(i); i = i -1;
    inputFiles.push_back("se/abdt-5/scanBDT-2016BF-41329.root");    colors.push_back(i); i = i -1;
    //    i = kMagenta;
    // inputFiles.push_back("se/abdt-5/scanBDT-2016BF-3369.root");    colors.push_back(i);  i = i -1;
    // inputFiles.push_back("se/abdt-5/scanBDT-2016BF-3529.root");    colors.push_back(i); i = i -1;
  } else if (6 == mode) {
    int i = kCyan;
    inputFiles.push_back("2016BF-00/scanBDT-2016BF-00.root");    colors.push_back(i);
    i = kBlue;
    inputFiles.push_back("2016BF-10/scanBDT-2016BF-10.root");    colors.push_back(i);
    i = kRed;
    inputFiles.push_back("2016BF-11/scanBDT-2016BF-11.root");    colors.push_back(i);
    i = kOrange;
    inputFiles.push_back("2016BF-12/scanBDT-2016BF-12.root");    colors.push_back(i);
    i = kSpring;
    inputFiles.push_back("2016BF-13/scanBDT-2016BF-13.root");    colors.push_back(i);
    i = kGreen+1;
    inputFiles.push_back("2016BF-14/scanBDT-2016BF-14.root");    colors.push_back(i);
    i = kBlack;
    inputFiles.push_back("2016BF-20/scanBDT-2016BF-20.root");    colors.push_back(i);
    i = kCyan+1;
  } else if (7 == mode) {
    int i = kCyan;
    inputFiles.push_back("2016GH-00/scanBDT-2016GH-00.root");    colors.push_back(i);
    i = kBlue;
    inputFiles.push_back("2016GH-10/scanBDT-2016GH-10.root");    colors.push_back(i);
    i = kRed;
    inputFiles.push_back("2016GH-11/scanBDT-2016GH-11.root");    colors.push_back(i);
    i = kOrange;
    inputFiles.push_back("2016GH-12/scanBDT-2016GH-12.root");    colors.push_back(i);
    i = kSpring;
    inputFiles.push_back("2016GH-13/scanBDT-2016GH-13.root");    colors.push_back(i);
    i = kGreen+1;
    inputFiles.push_back("2016GH-14/scanBDT-2016GH-14.root");    colors.push_back(i);
    i = kBlack;
    inputFiles.push_back("2016GH-20/scanBDT-2016GH-20.root");    colors.push_back(i);
  } else if (8 == mode) {
    int i = kBlue;
    inputFiles.push_back("s00/scanBDT-2016BFs00.root");    colors.push_back(i);
    i = kGreen+2;
    inputFiles.push_back("s01/scanBDT-2016BFs01.root");    colors.push_back(i);
    i = kRed;
    inputFiles.push_back("s02/scanBDT-2016BFs02.root");    colors.push_back(i);
    i = kMagenta;
    inputFiles.push_back("s03/scanBDT-2016BFs03.root");    colors.push_back(i);
  } else if (9 == mode) {
    int i = kBlue;
    inputFiles.push_back("s00/scanBDT-2016GHs00.root");    colors.push_back(i);
    i = kGreen+2;
    inputFiles.push_back("s01/scanBDT-2016GHs01.root");    colors.push_back(i);
    i = kRed;
    inputFiles.push_back("s02/scanBDT-2016GHs02.root");    colors.push_back(i);
    i = kMagenta;
    inputFiles.push_back("s03/scanBDT-2016GHs03.root");    colors.push_back(i);
  } else if (11 == mode) {
    int i = kBlue;
    inputFiles.push_back("BFs00_1479/scanBDT-2016BFs00_1479.root");    colors.push_back(i);
    i = kGreen+2;
    inputFiles.push_back("BFs01_3539/scanBDT-2016BFs01_3539.root");    colors.push_back(i);
    i = kRed;
    inputFiles.push_back("BFs02_1779/scanBDT-2016BFs02_1779.root");    colors.push_back(i);
    i = kMagenta;
    inputFiles.push_back("s03/scanBDT-2016BFs03.root");    colors.push_back(i);
  } else if (12 == mode) {
    int i = kBlue;
    inputFiles.push_back("GHs00_1159/scanBDT-2016GHs00_1159.root");    colors.push_back(i);
    i = kGreen+2;
    //    inputFiles.push_back("GHs01_7319/scanBDT-2016GHs01_7319.root");    colors.push_back(i);
    i = kRed;
    inputFiles.push_back("GHs02_2719/scanBDT-2016GHs02_2719.root");    colors.push_back(i);
    i = kMagenta;
    inputFiles.push_back("s03/scanBDT-2016GHs03.root");    colors.push_back(i);
  }

  string bname("hBdt_bsmmMcComb");
  if (string::npos != what.find("CSBF")) {
    bname = "hBdt_bspsiphiMcComb";
  }

  TFile *f(0);

  for (unsigned int ifile = 0; ifile < inputFiles.size(); ++ifile) {
    cout << "open " << inputFiles[ifile] << endl;
    f = TFile::Open(inputFiles[ifile].c_str());
    string hname(Form("bdtScan_%s", what.c_str()));
    string pdfname(inputFiles[ifile].c_str());
    replaceAll(pdfname, "results/", "");
    replaceAll(pdfname, "scanBDT-", "");
    replaceAll(pdfname, ".root", "");

    TH1D *h0(0), *h1(0);
    TH1D *b0(0), *b1(0);
    string s("");

    TH1D *hbdt = (TH1D*)f->Get("bdtCuts");
    cout << "hbdt = " << hbdt << endl;
    double bdtCut = 100.*hbdt->GetBinContent(chan+1);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    s = hname + Form("_chan%d", chan);
    h0 = (TH1D*)f->Get(s.c_str());
    if (!h0) {
      cout << "Did not find histogram ->" << s << "<-" << endl;
      return;
    }
    double bdty   = h0->GetBinContent(h0->FindBin(bdtCut));
    double bdtx(0.);
    for (int ic = 0; ic < h0->GetNbinsX(); ++ic) {
      if (h0->GetBinContent(ic) > 0.) bdtx = h0->GetBinCenter(ic);
    }
    bdtx /= 100.;

    h0->SetMinimum(0.0);
    if (what == "CSBF") {
      h0->SetMaximum(5.e-5);
    } else if (what == "BSMMBFS") {
      h0->SetMaximum(7.e-9);
    } else if (what == "SSB") {
      h0->SetMaximum(4.);
    } else if (what == "ZAD") {
      h0->SetMaximum(1.);
    } else if (what == "ZAS") {
      h0->SetMaximum(4.5);
    } else {
      h0->SetMaximum(1.4*h0->GetMaximum());
    }
    h0->SetLineColor(colors[ifile]);
    h0->SetMarkerColor(colors[ifile]);
    h0->SetMarkerStyle(2);
    h0->SetMarkerSize(0.6);
    setTitles(h0, "100 #times BDT >", what.c_str(), 0.05, 1.2, 1.6);
    if (0 == ifile) {
      h0->DrawCopy("p"); // "p"
    } else {
      h0->DrawCopy("psame");
    }
    TMarker *pm = new TMarker(bdtCut, bdty, 28);
    pm->SetMarkerColor(colors[ifile]);
    pm->SetMarkerSize(2.);
    pm->Draw();
    h0->GetListOfFunctions()->Add(pm);

    tl->SetTextSize(0.015);
    string bla = inputFiles[ifile];
    rmPath(bla);
    replaceAll(bla, "scanBDT-", "");
    replaceAll(bla, ".root", "");
    tl->SetTextColor(colors[ifile]); tl->DrawLatexNDC(0.75, 0.15 + ifile*0.016, bla.c_str());
    //    f->Close();
  }

  tl->SetTextSize(0.03);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.7, 0.92, Form("Chan %d", chan));
  c0->SaveAs(Form("%s/%s-%d-chan%d.pdf", fDirectory.c_str(), what.c_str(), mode, chan));

}



// ----------------------------------------------------------------------
void plotResults::dumpDatasets() {

  ofstream TEX;
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

  T->Draw("gpt>>tpt", "hlt1");
  T->Draw("g1pt>>tm1pt", "hlt1");
  T->Draw("g2pt>>tm2pt", "hlt1");

  T->Draw("geta>>teta", "hlt1");
  T->Draw("g1eta>>tm1eta", "hlt1");
  T->Draw("g2eta>>tm2eta", "hlt1");


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
  string mode("");

  if (0) {
    // -- ONLY for debugging
    resetHistograms();
    mode = "bdmmMcComb";
    setup(mode);
    t = getTree(mode, fTreeDir);
    setupTree(t, mode);
    //    loopOverTree(t, 1, 1000000, start);
    loopOverTree(t, 1, -1, start);
    saveHistograms(mode);
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
    mode = "bupsikData";
    setup(mode);
    t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    loopOverTree(t, 1, nevents, start);
    saveHistograms(fSample);

    resetHistograms();
    setup("bupsikMcComb");
    t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    loopOverTree(t, 1, nevents, start);
    cout << "done with loopovertree" << endl;
    otherNumbers(fSample);
    saveHistograms(fSample);

    resetHistograms();
    setup("bmmData");
    t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    loopOverTree(t, 1, nevents, start);
    fSample = "bmmData";
    saveHistograms(fSample);

    resetHistograms();
    setup("bdmmMcComb");
    t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    loopOverTree(t, 1, nevents, start);
    otherNumbers(fSample);
    saveHistograms(fSample);

    resetHistograms();
    setup("bsmmMcComb");
    t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    loopOverTree(t, 1, nevents, start);
    otherNumbers(fSample);
    saveHistograms(fSample);

    if (1) {
      resetHistograms();
      setup("bspsiphiData");
      t = getTree(fSample, fTreeDir);
      setupTree(t, fSample);
      loopOverTree(t, 1, nevents, start);
      fSample = "bspsiphiData";
      saveHistograms(fSample);

      resetHistograms();
      fSample = "bspsiphiMcComb";
      t = getTree(fSample, fTreeDir);
      setupTree(t, fSample);
      loopOverTree(t, 1, nevents, start);
      otherNumbers(fSample);
      saveHistograms(fSample);
    }
  }

  fHistFile->Close();

  fSaveSmallTree = false;
}


// ----------------------------------------------------------------------
void plotResults::initNumbers(anaNumbers &a) {
  a.clear();
  a.fSel = (fDoUseBDT?"bdt":"cnc");
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
    t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    cout << "==============================================================" << endl;
    cout << "==> rareBgHists for " << it->first << " and fSample = " << fSample << endl;
    cout << "==============================================================" << endl;
    //    loopOverTree(t, 1, 100000, start);
    loopOverTree(t, 1, nevents, start);
    otherNumbers(fSample);
    saveHistograms(fSample);
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
  if (string::npos != fSample.find("bupsik"))  {
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = &fNoNumbers[ichan];
  }
  if (string::npos != fSample.find("bspsiphi"))  {
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = &fCsNumbers[ichan];
  }

  if (string::npos != fSample.find("bdmm"))  {
    cout << "setting aa to fBdmmNumbers" << endl;
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = &fBdmmNumbers[ichan];
  }

  if (string::npos != fSample.find("bsmm"))  {
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = &fBsmmNumbers[ichan];
  }

  if (string::npos != fSample.find("Bg"))  {
    for (int ichan = 0; ichan < fNchan; ++ichan) aa[ichan] = fRareNumbers[fSample][ichan];
  }

  double effGenSel(0.), effGenSelE(0.);

  if (string::npos != fSample.find("Mc"))  {
    effGenSel = fDS[smode]->fFilterEff/fDS[accname]->fFilterEff;
    effGenSelE = dRatio(fDS[smode]->fFilterEff, fDS[smode]->fFilterEffE, fDS[accname]->fFilterEff, fDS[accname]->fFilterEffE);
  }

  cout << "smode = " << smode << " fSample = " << fSample
       << " accname: " << accname << " directory: " << fTreeDir
       << " numbers: " << aa[0]->fName << " chan = " << aa[0]->fChan
       << " effGenSel = " << effGenSel << " +/- " << effGenSelE
       << endl;

  // vector<string> modifier;
  // modifier.push_back("cnc" + fSuffix);
  // modifier.push_back("bdt" + fSuffix);

  for (unsigned int i = 0; i < fNchan; ++i) {
    fChan = i;
    // -- fill the numbers into anaNumbers
    getAccAndEffFromEffTree(accname, *aa[i], *fCuts[i], -1);
    // -- for tests on the acc samples, overwrite fEffGenSel
    if (string::npos != smode.find("Acc")) {
      aa[i]->fEffGenSel.val = 1.;
      aa[i]->fEffGenSel.estat = 0.;
    }
    fDS[fSample]->cd(fTreeDir.c_str());
    double effFilter  = fDS[fSample]->fFilterEff;
    double effFilterE = fDS[fSample]->fFilterEffE;
    if (effFilter < 1e-6) effFilter = 1.0;
    double genFileYield = ((TTree*)gDirectory->Get("effTree"))->GetEntries();

    for (unsigned int im = 0; im < fHistStrings.size(); ++im) {
      // -- fill numbers from anaNumbers into histogram
      ibin = 1;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, effFilter);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "effFilter");
      ibin = 2;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, effFilterE);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "effFilterE");
      ibin = 3;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, effGenSel);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "effGenSel");
      ibin = 4;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, effGenSelE);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "effGenSelE");

      ibin = 11;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, genFileYield);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "genFileYield");
      ibin = 12;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, genFileYield/effGenSel);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "genYield");
      ibin = 13;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, dRatio(genFileYield, TMath::Sqrt(genFileYield), effGenSel, effGenSelE));
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "genYieldE");

      ibin = 20;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, aa[i]->fAcc.val);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "acc");
      ibin = 21;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, aa[i]->fAcc.estat);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "accE");

      ibin = 22;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, aa[i]->fAccGenFileYield.val);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "accGenFileYield");
      ibin = 23;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, aa[i]->fAccGenYield.val);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "accGenYield");

      ibin = 24;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, aa[i]->fAccRecoYield.val);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "accRecoYield");
      ibin = 25;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, aa[i]->fAccCandYield.val);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "accCandYield");


      ibin = 30;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, aa[i]->fEffCand.val);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "effCand");
      ibin = 31;
      fhGenAndAccNumbers[fHistStrings[im]][i]->SetBinContent(ibin, aa[i]->fEffCand.estat);
      fhGenAndAccNumbers[fHistStrings[im]][i]->GetXaxis()->SetBinLabel(ibin, "effCandE");
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
void plotResults::scaleYield(anaNumbers &aSig, anaNumbers &aNorm, double pRatio, bool useW8) {
  cout << "+++ scaleYield: " << aSig.fNameMc << " wrt " << aNorm.fNameMc << endl;
  double sgBf  = fDS[aSig.fNameMc]->fBf;
  double sgBfE = fDS[aSig.fNameMc]->fBfE;

  double noBf  = fDS[aNorm.fNameMc]->fBf;
  double noBfE = fDS[aNorm.fNameMc]->fBfE;

  double yield(0.);
  if (useW8) {
    yield = (sgBf/noBf) * pRatio * (aSig.fEffTot.val/aNorm.fEffTot.val) * aNorm.fW8SignalFit.val;
  } else {
    yield = (sgBf/noBf) * pRatio * (aSig.fEffTot.val/aNorm.fEffTot.val) * aNorm.fSignalFit.val;
  }
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
      cout << aSig.fMcYield[i].val << " ";
      // FIXME error scaling
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
void plotResults::calculateNumbers(string mode, int chan) {
  cout << "==> calculateNumbers for mode: " << mode << " chan = " << chan << endl;

  // -- do two channels only
  if (chan > 1) {
    cout << "chan >1 requested, refusing to do so" << endl;
    return;
  }

  if (string::npos != mode.find("cnc")) {
    fDoUseBDT = false;
  } else {
    fDoUseBDT = true;
  }
  // -- open histogram file
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

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

  // -- calculate Bs-> J/psi phi BF
  calculateBs2Bu(chan);
  calculatePerformance(chan);

  fHistFile->Close();
}

// ----------------------------------------------------------------------
void plotResults::calculateBs2Bu(int ichan) {

  double bsBf  = fDS[fCsNumbers[ichan].fNameMc]->fBf;
  double bsBfE = fDS[fCsNumbers[ichan].fNameMc]->fBfE;

  double buBf  = fDS[fNoNumbers[ichan].fNameMc]->fBf;
  double buBfE = fDS[fNoNumbers[ichan].fNameMc]->fBfE;


  number bf;
  double alpha = (1./fFsfu.val)
    * (fNoNumbers[ichan].fEffTot.val / fCsNumbers[ichan].fEffTot.val)
    * buBf;

  bf.val   = (fCsNumbers[ichan].fW8SignalFit.val / fNoNumbers[ichan].fW8SignalFit.val) * alpha;


  cout << "*******************************" << endl;
  cout << "psiphi w8 yield: " << fCsNumbers[ichan].fW8SignalFit.val << " +/- " << fCsNumbers[ichan].fW8SignalFit.estat << endl;
  cout << "psik w8 yield:   " << fNoNumbers[ichan].fW8SignalFit.val << " +/- " << fNoNumbers[ichan].fW8SignalFit.estat << endl;
  cout << "BF:              " << bf.val << " +/- " << bf.estat << " for selection " << fNoNumbers[ichan].fSel << endl;
  cout << "*******************************" << endl;

  int NDIG(6);
  // string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  // fSuffixSel = modifier;
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- BF(Bs -> J/psi phi):  chan: " << ichan << " BFCS/CSBF" << endl;
  dumpTex(bf, Form("%s:CSBF:chan%d", fSuffixSel.c_str(), ichan), 7);
  dumpTex(fNoNumbers[ichan].fW8SignalFit, Form("%s:BUNW8:chan%d", fSuffixSel.c_str(), ichan), 7);
  dumpTex(fNoNumbers[ichan].fSignalFit, Form("%s:BUNW0:chan%d", fSuffixSel.c_str(), ichan), 7);
  dumpTex(fCsNumbers[ichan].fW8SignalFit, Form("%s:BSNW8:chan%d", fSuffixSel.c_str(), ichan), 7);
  dumpTex(fCsNumbers[ichan].fSignalFit, Form("%s:BSNW0:chan%d", fSuffixSel.c_str(), ichan), 7);
  }


// ----------------------------------------------------------------------
void plotResults::calculatePerformance(int ichan) {

  double bsBf  = fDS[fCsNumbers[ichan].fNameMc]->fBf;
  double bsBfE = fDS[fCsNumbers[ichan].fNameMc]->fBfE;

  double buBf  = fDS[fNoNumbers[ichan].fNameMc]->fBf;
  double buBfE = fDS[fNoNumbers[ichan].fNameMc]->fBfE;

  number bfU, bfS;
  double alphaU = (1./fFsfu.val)
    * (fNoNumbers[ichan].fEffTot.val / fBsmmNumbers[ichan].fEffTot.val)
    * buBf
    / fNoNumbers[ichan].fW8SignalFit.val;
  double alphaUE = TMath::Sqrt(
			       (fFsfu.etot/fFsfu.val)*(fFsfu.etot/fFsfu.val)
			       + (fNoNumbers[ichan].fEffTot.etot/fNoNumbers[ichan].fEffTot.val)
			       * (fNoNumbers[ichan].fEffTot.etot/fNoNumbers[ichan].fEffTot.val)
			       + (fBsmmNumbers[ichan].fEffTot.etot/fBsmmNumbers[ichan].fEffTot.val)
			       * (fBsmmNumbers[ichan].fEffTot.etot/fBsmmNumbers[ichan].fEffTot.val)
			       + (buBfE/buBf)*(buBfE/buBf)
			       + (fNoNumbers[ichan].fW8SignalFit.etot/fNoNumbers[ichan].fW8SignalFit.val)
			       * (fNoNumbers[ichan].fW8SignalFit.etot/fNoNumbers[ichan].fW8SignalFit.val)
			       );
  alphaUE = alphaUE * alphaU;

  double alphaS =
    (fCsNumbers[ichan].fEffTot.val / fBsmmNumbers[ichan].fEffTot.val)
    * bsBf
    / fCsNumbers[ichan].fW8SignalFit.val;
  double alphaSE = TMath::Sqrt(
			       + (fCsNumbers[ichan].fEffTot.etot/fCsNumbers[ichan].fEffTot.val)
			       * (fCsNumbers[ichan].fEffTot.etot/fCsNumbers[ichan].fEffTot.val)
			       + (fBsmmNumbers[ichan].fEffTot.etot/fBsmmNumbers[ichan].fEffTot.val)
			       * (fBsmmNumbers[ichan].fEffTot.etot/fBsmmNumbers[ichan].fEffTot.val)
			       + (bsBfE/bsBf)*(bsBfE/bsBf)
			       + (fCsNumbers[ichan].fW8SignalFit.etot/fCsNumbers[ichan].fW8SignalFit.val)
			       * (fCsNumbers[ichan].fW8SignalFit.etot/fCsNumbers[ichan].fW8SignalFit.val)
			       );
  alphaSE   = alphaSE * alphaS;

  bfU.val   = fBsmmNumbers[ichan].fScaledYield.val  * alphaU;
  bfU.estat = 0.;
  bfU.esyst = 0.;
  bfU.etot  = TMath::Sqrt(alphaUE*alphaUE*fBsmmNumbers[ichan].fScaledYield.val*fBsmmNumbers[ichan].fScaledYield.val
			  + alphaU*alphaU*fBsmmNumbers[ichan].fScaledYield.etot*fBsmmNumbers[ichan].fScaledYield.etot);

  bfS.val   = fBsmmNumbers[ichan].fScaledYield.val  * alphaS;
  bfS.estat = 0.;
  bfS.esyst = 0.;
  bfS.etot  = TMath::Sqrt(alphaSE*alphaSE*fBsmmNumbers[ichan].fScaledYield.val*fBsmmNumbers[ichan].fScaledYield.val
			  + alphaS*alphaS*fBsmmNumbers[ichan].fScaledYield.etot*fBsmmNumbers[ichan].fScaledYield.etot);

  cout << "*******************************" << endl;
  cout << "expected Bsmm yield:           " << fBsmmNumbers[ichan].fScaledYield.val << " +/- " << fBsmmNumbers[ichan].fScaledYield.etot << endl;
  cout << "psik w8 yield:                 " << fNoNumbers[ichan].fW8SignalFit.val << " +/- " << fNoNumbers[ichan].fW8SignalFit.etot << endl;
  cout << "alphaU:                        " << alphaU << " +/- " << alphaUE  << endl;
  cout << "BF (B+ norm):                  " << bfU.val << " +/- " << bfU.etot << " for selection " << fNoNumbers[ichan].fSel << endl;
  cout << "*******************************" << endl;

  cout << "*******************************" << endl;
  cout << "expected Bsmm yield:           " << fBsmmNumbers[ichan].fScaledYield.val << " +/- " << fBsmmNumbers[ichan].fScaledYield.etot << endl;
  cout << "psiphi w8 yield:               " << fCsNumbers[ichan].fW8SignalFit.val << " +/- " << fCsNumbers[ichan].fW8SignalFit.etot << endl;
  cout << "alphaS:                        " << alphaS << " +/- " << alphaSE << endl;
  cout << "BF (Bs norm):                  " << bfS.val << " +/- " << bfS.etot << " for selection " << fCsNumbers[ichan].fSel << endl;
  cout << "*******************************" << endl;

  int NDIG(6);
  // string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  // fSuffixSel = modifier;
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- BF(Bs -> mumu; Bu -> psi K):  chan: " << ichan << endl;
  fTEX << formatTexErrSci(alphaU, alphaUE, Form("%s:alphaU:chan%d", fSuffixSel.c_str(), ichan), 4) << endl;
  fTEX << formatTex(alphaU, Form("%s:alphaU:chan%d:val", fSuffixSel.c_str(), ichan), 11) << endl;
  dumpTex(bfU, Form("%s:BSMMBFU:chan%d", fSuffixSel.c_str(), ichan), 10);
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- BF(Bs -> mumu; Bs -> psi phi):  chan: " << ichan << endl;
  fTEX << formatTexErrSci(alphaS, alphaSE, Form("%s:alphaS:chan%d", fSuffixSel.c_str(), ichan), 4) << endl;
  dumpTex(bfS, Form("%s:BSMMBFS:chan%d", fSuffixSel.c_str(), ichan), 10);
  fTEX << formatTex(alphaS, Form("%s:alphaS:chan%d:val", fSuffixSel.c_str(), ichan), 11) << endl;

  // -- calculate S/B here (for Bs -> mu mu)
  double s   = fBsmmNumbers[ichan].fScaledYield.val;
  double b   = fBgNumbers[ichan].fFitYield[2].val;
  double sob = s/b;
  double ssb = s/TMath::Sqrt(s+b);
  if (b < 0.01) b = 0.01;
  double be  = fBgNumbers[ichan].fFitYield[2].etot;
  if (be/b < 0.01) be = 0.1*b;
  // -- this is from Glen Cowan "Discovery sensitivity for a counting experiment with background uncertainty", eq (20)
  //    https://www.pp.rhul.ac.uk/~cowan/stat/medsig/medsigNote.pdf
  double zas = TMath::Sqrt(2.*(
			       (s+b)*TMath::Log(((s+b)*(b+be*be))/(b*b + (s+b)*be*be))
			       - ((b*b)/(be*be))*TMath::Log(1. + ((be*be*s)/(b*(b+be*be))))
			       )
			   );
  s   = fBdmmNumbers[ichan].fScaledYield.val;
  b   = fBgNumbers[ichan].fFitYield[1].val;
  if (b < 0.01) b = 0.01;
  be  = fBgNumbers[ichan].fFitYield[1].etot;
  if (be/b < 0.01) be = 0.1*b;
  double zad = TMath::Sqrt(2.*(
			       (s+b)*TMath::Log(((s+b)*(b+be*be))/(b*b + (s+b)*be*be))
			       - ((b*b)/(be*be))*TMath::Log(1. + ((be*be*s)/(b*(b+be*be))))
			       )
			   );
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- SOB chan: " << ichan << endl;
  fTEX << formatTex(sob, Form("%s:SOB:chan%d:val", fSuffixSel.c_str(), ichan), 4)
       << endl;
  fTEX << formatTex(ssb, Form("%s:SSB:chan%d:val", fSuffixSel.c_str(), ichan), 4)
       << endl;
  fTEX << formatTex(zad, Form("%s:ZAD:chan%d:val", fSuffixSel.c_str(), ichan), 4)
       << endl;
  fTEX << formatTex(zas, Form("%s:ZAS:chan%d:val", fSuffixSel.c_str(), ichan), 4)
       << endl;
  fTEX << formatTex(s,
		    Form("%s:SOBS:chan%d:val", fSuffixSel.c_str(), ichan), 4)
       << endl;
  fTEX << formatTex(b,
		    Form("%s:SOBB:chan%d:val", fSuffixSel.c_str(), ichan), 4)
       << endl;

}



// ----------------------------------------------------------------------
void plotResults::numbersFromHist(anaNumbers &aa, string syst) {
  // previously mode had been used to differentiate between bdmm, bsmm, no, cs, ...
  int chan = aa.fChan;

  // -- efficiency and acceptance
  // string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  string modifier = fSuffixSel + string("_") + aa.fNameMc;

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
  aa.fMcYield[0].estat = TMath::Sqrt(lo);
  aa.fMcYield[1].val   = bd;
  aa.fMcYield[1].estat = TMath::Sqrt(bd);
  aa.fMcYield[2].val   = bs;
  aa.fMcYield[2].estat = TMath::Sqrt(bs);
  aa.fMcYield[3].val   = hi;
  aa.fMcYield[3].estat = TMath::Sqrt(hi);
  aa.fMcYield[4].val   = tot;
  aa.fMcYield[4].estat = TMath::Sqrt(tot);
}


// ----------------------------------------------------------------------
void plotResults::calculateB2JpsiNumbers(anaNumbers &a) {
  c0->Clear();
  c0->SetCanvasSize(700, 700);
  // string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  // fSuffixSel = modifier;
  cout << "==> calculateB2JpsiNumbers for name: " << a.fName << ", chan: " << a.fChan << " fSuffixSel: " << fSuffixSel << endl;

  // -- MC: efficiency and acceptance
  char mode[200];
  sprintf(mode, "%s", a.fName.c_str());
  fSample = a.fNameMc;
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
  fSample = a.fNameDa;
  //  string  name = Form("hNorm_%s_%s_chan%d", modifier.c_str(), fSample.c_str(), chan);
  string  name = Form("hNorm_%s_%s_chan%d", fSuffixSel.c_str(), fSample.c_str(), chan);
  bool ok = fHistFile->cd(fSample.c_str());
  cout << "cd to " << fSample << ": " << ok << ", normalization fitting: " << name << endl;
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

  // -- fit also weighted yields (with correction weights)
  //  string name2 = Form("hW8Norm_%s_%s_chan%d", modifier.c_str(), fSample.c_str(), chan);
  string name2 = Form("hW8Norm_%s_%s_chan%d", fSuffixSel.c_str(), fSample.c_str(), chan);
  ok = fHistFile->cd(fSample.c_str());
  cout << "cd to " << fSample << ": " << ok << ", W8 normalization fitting: " << name << endl;
  fitPsYield fpy2(name2, 0);
  if (string::npos != a.fName.find("bupsik")) {
    fpy2.fitBu2JpsiKp(5, fDirectory + "/");
  } else if (string::npos != a.fName.find("bspsiphi")) {
    fpy2.fitBs2JpsiPhi(5, fDirectory + "/");
  }
  a.fW8SignalFit.val   = fpy2.getSignalYield();
  a.fW8SignalFit.estat = fpy2.getSignalError();
  a.fW8SignalFit.esyst = fSystematics["norm" + a.fName][chan] * fpy2.getSignalYield();
  a.fW8SignalFit.etot  = TMath::Sqrt(a.fW8SignalFit.estat*a.fW8SignalFit.estat + a.fW8SignalFit.esyst*a.fW8SignalFit.esyst);

  printNumbers(a, cout);

  cout << "chan " << a.fChan << " total yield: " << fpy.getSignalYield() << " +/- "  << fpy.getSignalError() << endl;
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- NORMALIZATION: " << mode << " chan: " << chan << endl;
  int NDIG(5);
  dumpTex(a.fEffGenSel, Form("%s:GENSEL-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fAcc, Form("%s:ACC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffCand, Form("%s:EFF-CAND-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffMuidMC, Form("%s:EFF-MU-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTrigMC, Form("%s:EFF-TRIG-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffAna, Form("%s:EFF-ANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffProdMC, Form("%s:EFF-PRODMC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTot, Form("%s:EFF-TOT-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);

  NDIG = 3;
  dumpTex(a.fGenFileYield, Form("%s:N-GENFILEYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fGenYield, Form("%s:N-GENYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fTotGenYield, Form("%s:N-TOTYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fProdGenYield, Form("%s:N-PRODYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);

  if (string::npos != a.fName.find("bspsiphi")) {
    dumpTex(a.fScaledYield, Form("%s:N-SCALEDYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  }
  dumpTex(a.fSignalFit,   Form("%s:N-OBS-%s-chan%d", fSuffixSel.c_str(), mode, chan), 1);
  dumpTex(a.fW8SignalFit, Form("%s:N-W8OBS-%s-chan%d", fSuffixSel.c_str(), mode, chan), 1);

  c0->Modified();
  c0->Update();
  fHistFile->cd();
}


// ----------------------------------------------------------------------
void plotResults::calculateCombBgNumbers(anaNumbers &a, int mode, double lo, double hi) {
  string hname = fHistWithAllCuts;
  c0->SetCanvasSize(700, 700);

  c0->Clear();
  // string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  // fSuffixSel = modifier;
  cout << "==> calculateCombBgNumbers for name: " << a.fName << ", chan: " << a.fChan << " fSuffixSel = " << fSuffixSel << endl;

  // -- get the histogram
  fSample = a.fNameDa;
  //  string  name = Form("%s_%s_%s_chan%d", hname.c_str(), modifier.c_str(), fSample.c_str(), a.fChan);
  string  name = Form("%s_%s_%s_chan%d", hname.c_str(), fSuffixSel.c_str(), fSample.c_str(), a.fChan);
  TH1D *h1 = (TH1D*)gDirectory->Get(Form("%s/%s", fSample.c_str(), name.c_str()));
  cout << "getting  histogram ->" << Form("%s/%s", fSample.c_str(), name.c_str()) << "<-" << endl;
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
  c0->SetCanvasSize(700, 700);
  // string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  // fSuffixSel = modifier;

  // -- MC: efficiency and acceptance
  char mode[200];
  sprintf(mode, "%s", a.fName.c_str());
  fSample = a.fNameMc;
  int chan = a.fChan;
  numbersFromHist(a, "bsmm");
  double pRatio(fFsfu.val);
  if (string::npos != a.fName.find("bdmm")) pRatio = 1.;
  scaleYield(a, fNoNumbers[chan], pRatio);

  // -- data: fitted/interpolated yields
  fSample = a.fNameDa;

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
  c0->Clear();
  c0->SetCanvasSize(700, 700);
  c0->cd();
  shrinkPad(0.12, 0.15, 0.10);

  cout << "==> calculateRareBgNumbers for chan: " << chan << endl;
  int nloop(0);
  // string modifier = (fDoUseBDT?"bdt":"cnc") + fSuffix;
  // fSuffixSel = modifier;
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
    cout << "rec: " << Form("%s/hGenAndAccNumbers_%s_%s_chan%d", it->first.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan) << endl;
    // TH1D *h1 = (TH1D*)fHistFile->Get(Form("%s/hGenAndAccNumbers_%s_%s_chan%d", it->first.c_str(), modifier.c_str(), it->first.c_str(), chan));
    TH1D *h1 = (TH1D*)fHistFile->Get(Form("%s/hGenAndAccNumbers_%s_%s_chan%d", it->first.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan));
    double effGenSel = getValueByLabel(h1,  "effGenSel");
    double genYield  = getValueByLabel(h1,  "genYield");
    double acc       = getValueByLabel(h1,  "acc");

    // u8name = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), uname.c_str(), modifier.c_str(), it->first.c_str(), chan);
    // w8name = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), wname.c_str(), modifier.c_str(), it->first.c_str(), chan);
    u8name = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), uname.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan);
    w8name = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), wname.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan);
    TH1D *hu = (TH1D*)fHistFile->Get(u8name.c_str());
    TH1D *hw = (TH1D*)fHistFile->Get(w8name.c_str());

    double utot     = massIntegral(hu, ALL, chan);
    double ulo      = massIntegral(hu, LO, chan);
    double ubd      = massIntegral(hu, BD, chan);
    double ubs      = massIntegral(hu, BS, chan);
    double uhi      = massIntegral(hu, HI, chan);
    double eAccAna  = utot/genYield; // this does NOT yet include misid and trigger efficiency, will be applied later on
    double eTot     = eAccAna * (massIntegral(hw, ALL, chan) / utot) * fBsmmNumbers[chan].fEffTrigMC.val;

    fSample = a->fNameMc;
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

  // -- add combinatorial background for more combined sums FIXME add errors
  for (unsigned im = 0; im < NWIN; ++im) {
    fNpNumbers[chan].fMcYield[im].val +=  fSlNumbers[chan].fMcYield[im].val;
    fNpNumbers[chan].fMcYield[im].add2Errors(fSlNumbers[chan].fMcYield[im]);
    fNpNumbers[chan].fMcYield[im].val +=  fCombNumbers[chan].fFitYield[im].val;
    fNpNumbers[chan].fMcYield[im].add2Errors(fCombNumbers[chan].fFitYield[im]);

    fBgNumbers[chan].fMcYield[im].val +=  fHhNumbers[chan].fMcYield[im].val;
    fBgNumbers[chan].fMcYield[im].add2Errors(fHhNumbers[chan].fMcYield[im]);
    fBgNumbers[chan].fMcYield[im].val +=  fNpNumbers[chan].fMcYield[im].val;
    fBgNumbers[chan].fMcYield[im].add2Errors(fNpNumbers[chan].fMcYield[im]);

    // -- the sl numbers are added later on after "fitting" TRUE?? FIXME
    fSgAndBgNumbers[chan].fMcYield[im].val +=  fNpNumbers[chan].fMcYield[im].val;
    fSgAndBgNumbers[chan].fMcYield[im].add2Errors(fNpNumbers[chan].fMcYield[im]);
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
  TLegend *lBg = ::newLegend("rare decays", 0.50, 0.2, 0.85, 0.85, vBg, vBgnames, vBgoptions);
  lBg->Draw();
  c0->Modified();
  c0->Update();
  if (1) {
    savePad(Form("%s-rare-bg-chan%d.pdf", fSuffixSel.c_str(), chan));
  }

  // -- create overlay of data and the stacked (and scaled) backgrounds
  fSample = fCombNumbers[chan].fNameDa;
  //  string  name = Form("%s_%s_%s_chan%d", hname.c_str(), modifier.c_str(), fSample.c_str(), chan);
  string  name = Form("%s_%s_%s_chan%d", hname.c_str(), fSuffixSel.c_str(), fSample.c_str(), chan);
  cout << "getting  histogram ->" << Form("%s/%s", fSample.c_str(), name.c_str()) << "<-" << endl;
  TH1D *h1e = (TH1D*)((TH1D*)gDirectory->Get(Form("%s/%s", fSample.c_str(), name.c_str())))->Clone("h1erebin");
  h1e->Rebin(5);
  TH1D *h1c = (TH1D*)((TH1D*)gDirectory->Get(Form("%s/%s", fSample.c_str(), name.c_str())))->Clone("h1crebin");
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
  if (scale < -1.) scale = -1.;
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
  tl->DrawLatexNDC(0.55, 0.92, Form("sl correction: = %+3.1f %%", 100*scale));

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
    fSlNumbers[chan].fFitYield[iw].val      = fSlNumbers[chan].fMcYield[iw].val * fSlNumbers[chan].fScaleFactor;
    fNpNumbers[chan].fFitYield[iw].val      = fNpNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
    fBgNumbers[chan].fFitYield[iw].val      = fBgNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
    fSgAndBgNumbers[chan].fFitYield[iw].val = fSgAndBgNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
  }
  // -- add signal
  for (int iw = 0; iw < NWIN; ++iw) {
    fSgAndBgNumbers[chan].fFitYield[iw].val = fSgAndBgNumbers[chan].fMcYield[iw].val + fBsmmNumbers[chan].fMcYield[iw].val + fBdmmNumbers[chan].fMcYield[iw].val;
  }

  // -- dump combined/summed numbers: rare sl decays and rare hadronic (peaking) decays
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- SUMMARY BACKGROUND " << mode << " chan " << chan << endl;
  fTEX << formatTex(100*scale, Form("%s:SCALEPERCENT-SL-chan%d", fSuffixSel.c_str(), chan), 1) << endl;
  fTEX << formatTex(scale + 1.0, Form("%s:SCALEPLUSONE-SL-chan%d", fSuffixSel.c_str(), chan), 3) << endl;
  fTEX << formatTex(scale, Form("%s:SCALE-SL-chan%d", fSuffixSel.c_str(), chan), 3) << endl;

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
  if (fChan >= fNchan) return;
  double mass = fb.m;
  double ps   = static_cast<double>(fb.ps);
  if (fYear < 2016) {
    ps = 1.;
  }

  for (unsigned int im = 0; im < fHistStrings.size(); ++im) {
    fhMassAbsNoCuts[fHistStrings[im]][fChan]->Fill(mass);
    fhAccAll[fHistStrings[im]][fChan]->Fill(TMath::Abs(fb.eta), fb.pt);
    fhAccPtAll[fHistStrings[im]][fChan]->Fill(fb.pt);
    fhAccEtaAll[fHistStrings[im]][fChan]->Fill(TMath::Abs(fb.eta));
  }

  if (!fGoodAcceptance) return;

  for (unsigned int im = 0; im < fHistStrings.size(); ++im) {
    fhAccPass[fHistStrings[im]][fChan]->Fill(TMath::Abs(fb.eta), fb.pt);
    fhAccPtPass[fHistStrings[im]][fChan]->Fill(fb.pt);
    fhAccEtaPass[fHistStrings[im]][fChan]->Fill(TMath::Abs(fb.eta));
    // -- this is the base, after the raw acceptance cuts
    fhMassNoCuts[fHistStrings[im]][fChan]->Fill(mass);
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
      && fGoodPhi
      && fGoodEta
      && fGoodAlpha
      && fGoodChi2
      && fGoodFLS
      && fGoodCloseTrack
      && fGoodIso
      && fGoodDocaTrk
      ) {
    fhMassWithAnaCuts[fHistStrings[0]][fChan]->Fill(mass);

    if (fGoodMuonsID) {
      fhMassWithMuonCuts[fHistStrings[0]][fChan]->Fill(mass);
      if (fGoodHLT) {
	fhMassWithTriggerCuts[fHistStrings[0]][fChan]->Fill(mass);
	fhMassWithAllCuts[fHistStrings[0]][fChan]->Fill(mass);
	if (fIsCowboy) {
	  fhMassWithAllCutsCowboy[fHistStrings[0]][fChan]->Fill(mass);
	} else {
	  fhMassWithAllCutsSeagull[fHistStrings[0]][fChan]->Fill(mass);
	}

	// -- blind version (to have the possibility to make a blind plot also after unblinding)
	if ((5.2 < mass) && (mass < 5.45)) {
	} else {
	  fhMassWithAllCutsBlind[fHistStrings[0]][fChan]->Fill(mass);
	  if (fIsCowboy) {
	    fhMassWithAllCutsCowboyBlind[fHistStrings[0]][fChan]->Fill(mass);
	  } else {
	    fhMassWithAllCutsSeagullBlind[fHistStrings[0]][fChan]->Fill(mass);
	  }
	}

	// -- weighted with fake rate
	fhW8MassWithAllCuts[fHistStrings[0]][fChan]->Fill(mass, fW8MisId);
	if (fIsCowboy) {
	  fhW8MassWithAllCutsCowboy[fHistStrings[0]][fChan]->Fill(mass, fW8MisId);
	} else {
	  fhW8MassWithAllCutsSeagull[fHistStrings[0]][fChan]->Fill(mass, fW8MisId);
	}


	// - include prescale values on y axis
	if ((fMode == BS2JPSIPHI)
	    || (fMode == BD2JPSIKSTAR)
	    || (fMode == BU2JPSIKP)) {
	  fhNorm[fHistStrings[0]][fChan]->Fill(mass, -0.1, ps);
	  fhNorm[fHistStrings[0]][fChan]->Fill(mass, 0.1);
	  fhNorm[fHistStrings[0]][fChan]->Fill(mass, ps+0.1);
	  fhNormC[fHistStrings[0]][fChan]->Fill(fb.cm, -0.1, ps);
	  fhNormC[fHistStrings[0]][fChan]->Fill(fb.cm, 0.1);
	  fhNormC[fHistStrings[0]][fChan]->Fill(fb.cm, ps+0.1);

	  fhW8Norm[fHistStrings[0]][fChan]->Fill(mass, -0.1, fb.corrW8*ps);
	  fhW8Norm[fHistStrings[0]][fChan]->Fill(mass, 0.1, fb.corrW8);
	  fhW8Norm[fHistStrings[0]][fChan]->Fill(mass, ps+0.1, fb.corrW8);
	  fhW8NormC[fHistStrings[0]][fChan]->Fill(fb.cm, -0.1, fb.corrW8*ps);
	  fhW8NormC[fHistStrings[0]][fChan]->Fill(fb.cm, 0.1, fb.corrW8);
	  fhW8NormC[fHistStrings[0]][fChan]->Fill(fb.cm, ps+0.1, fb.corrW8);

	}


	if (fMode == BSMM && fCuts[fChan]->mBsLo < mass && mass < fCuts[fChan]->mBsHi) {
	  fhMassWithMassCuts[fHistStrings[0]][fChan]->Fill(mass);
	}

	if (fMode == BDMM && fCuts[fChan]->mBdLo < mass && mass < fCuts[fChan]->mBdHi) {
	  fhMassWithMassCuts[fHistStrings[0]][fChan]->Fill(mass);
	}

	if (fMode == BU2JPSIKP && fNoLo < mass && mass < fNoHi) {
	  fhMassWithMassCuts[fHistStrings[0]][fChan]->Fill(mass);
	}

	if (fMode == BS2JPSIPHI && fCsLo < mass && mass < fCsHi) {
	  fhMassWithMassCuts[fHistStrings[0]][fChan]->Fill(mass);
	}

	if (fMode == BD2JPSIKSTAR && fNoLo < mass && mass < fNoHi) {
	  fhMassWithMassCuts[fHistStrings[0]][fChan]->Fill(mass);
	}

      }
    }
  }

  // -----------------
  // -- BDT histograms
  // -----------------
  bool goodBDT(false);
  goodBDT = fGoodBDT;

  for (unsigned int i = 1; i < fHistStrings.size(); ++i) {
    if (1 == i) {
      goodBDT = fGoodBDT;
    } else {
      goodBDT = (fBDT > (i-2)*0.01);
    }

    if (fGoodQ
	&& fGoodPvAveW8
	&& fGoodTracks
	&& fGoodTracksPt
	&& fGoodTracksEta
	&& fGoodMuonsPt  // PidTables do not really work below 4 GeV!!
	&& fGoodMuonsEta
	&& fGoodJpsiCuts
	&& fGoodMuonsID
	&& fGoodDcand
	&& fGoodHLT
	) {
      fhBdt[fChan]->Fill(fBDT);
    }

    if (fGoodQ
	&& fGoodPvAveW8
	&& fGoodTracks
	&& fGoodTracksPt
	&& fGoodTracksEta
	&& fGoodMuonsPt  // PidTables do not really work below 4 GeV!!
	&& fGoodMuonsEta
	&& fGoodJpsiCuts
	&& goodBDT
	) {

      fhBdtCrossCheck[fHistStrings[i]][fChan]->Fill(fBDT);

      fhMassWithAnaCuts[fHistStrings[i]][fChan]->Fill(mass);

      if (fGoodMuonsID && fGoodDcand) {
	fhMassWithMuonCuts[fHistStrings[i]][fChan]->Fill(mass);
	if (fGoodHLT) {
	  fhMassWithTriggerCuts[fHistStrings[i]][fChan]->Fill(mass);
	  fhMassWithAllCuts[fHistStrings[i]][fChan]->Fill(mass);
	  if (fIsCowboy) {
	    fhMassWithAllCutsCowboy[fHistStrings[i]][fChan]->Fill(mass);
	  } else {
	    fhMassWithAllCutsSeagull[fHistStrings[i]][fChan]->Fill(mass);
	  }
	  // -- blind version
	  if ((5.2 < mass) && (mass < 5.45)) {
	  } else {
	    fhMassWithAllCutsBlind[fHistStrings[i]][fChan]->Fill(mass);
	    if (fIsCowboy) {
	      fhMassWithAllCutsCowboyBlind[fHistStrings[i]][fChan]->Fill(mass);
	    } else {
	      fhMassWithAllCutsSeagullBlind[fHistStrings[i]][fChan]->Fill(mass);
	    }
	  }

	  // -- weighted with fake rate
	  if (0) cout << "fW8MisId = " << fW8MisId
		      << " m1eta =  " << fb.m1eta << " m2eta =  " << fb.m2eta
		      << " m1pt =  " << fb.m1pt << " m2pt =  " << fb.m2pt
		      << " g1id =  " << fb.g1id << " g2id =  " << fb.g2id
		      << endl;

	  fhW8MassWithAllCuts[fHistStrings[i]][fChan]->Fill(mass, fW8MisId);
	  if (fIsCowboy) {
	    fhW8MassWithAllCutsCowboy[fHistStrings[i]][fChan]->Fill(mass, fW8MisId);
	  } else {
	    fhW8MassWithAllCutsSeagull[fHistStrings[i]][fChan]->Fill(mass, fW8MisId);
	  }

	  // - include prescale values on y axis
	  if ((fMode == BS2JPSIPHI)
	      || (fMode == BD2JPSIKSTAR)
	      || (fMode == BU2JPSIKP)) {
	    fhNorm[fHistStrings[i]][fChan]->Fill(mass, -0.1, ps);
	    fhNorm[fHistStrings[i]][fChan]->Fill(mass, 0.1);
	    fhNorm[fHistStrings[i]][fChan]->Fill(mass, ps+0.1);
	    fhNormC[fHistStrings[i]][fChan]->Fill(fb.cm, -0.1, ps);
	    fhNormC[fHistStrings[i]][fChan]->Fill(fb.cm, 0.1);
	    fhNormC[fHistStrings[i]][fChan]->Fill(fb.cm, ps+0.1);

	    fhW8Norm[fHistStrings[i]][fChan]->Fill(mass, -0.1, fb.corrW8*ps);
	    fhW8Norm[fHistStrings[i]][fChan]->Fill(mass, 0.1, fb.corrW8);
	    fhW8Norm[fHistStrings[i]][fChan]->Fill(mass, ps+0.1, fb.corrW8);
	    fhW8NormC[fHistStrings[i]][fChan]->Fill(fb.cm, -0.1, fb.corrW8*ps);
	    fhW8NormC[fHistStrings[i]][fChan]->Fill(fb.cm, 0.1, fb.corrW8);
	    fhW8NormC[fHistStrings[i]][fChan]->Fill(fb.cm, ps+0.1, fb.corrW8);
	  }

	  if (fMode == BSMM && fCuts[fChan]->mBsLo < mass && mass < fCuts[fChan]->mBsHi) {
	    fhMassWithMassCuts[fHistStrings[i]][fChan]->Fill(mass);
	  } else if (fMode == BDMM && fCuts[fChan]->mBdLo < mass && mass < fCuts[fChan]->mBdHi) {
	    fhMassWithMassCuts[fHistStrings[i]][fChan]->Fill(mass);
	  } else if (fMode == BU2JPSIKP && fNoLo < mass && mass < fNoHi) {
	    fhMassWithMassCuts[fHistStrings[i]][fChan]->Fill(mass);
	  } else if (fMode == BS2JPSIPHI && fCsLo < mass && mass < fCsHi) {
	    fhMassWithMassCuts[fHistStrings[i]][fChan]->Fill(mass);
	  } else if (fMode == BD2JPSIKSTAR && fNoLo < mass && mass < fNoHi) {
	    fhMassWithMassCuts[fHistStrings[i]][fChan]->Fill(mass);
	  }
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
    string tname(fSample);
    replaceAll(tname, "Off", "");
    replaceAll(tname, "Comb", "");

    dir = gDirectory;
    fLocal = TFile::Open(Form("%s/small%s-%s.root", fDirectory.c_str(), fSetup.c_str(), tname.c_str()), "RECREATE");
    small = new TTree(Form("%s", tname.c_str()), Form("%s", tname.c_str()));
    small->SetDirectory(fLocal);
    small->Branch("run",    &fb.run,       "run/I");
    small->Branch("evt",    &fb.evt,       "evt/I");
    small->Branch("ls",     &fb.ls,        "ls/I");
    small->Branch("ps",     &fb.ps,        "ps/I");
    small->Branch("cw8",    &fb.corrW8,    "cw8/D");

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
    small->Branch("tauxy",  &fb.tauxy,     "tauxy/D");
    small->Branch("tauxye", &fb.tauxye,    "tauxye/D");
    small->Branch("gtau",   &fb.gtau,      "gtau/D");

    small->Branch("m1pt",   &fb.m1pt,      "m1pt/D");
    small->Branch("m1eta",  &fb.m1eta,     "m1eta/D");
    small->Branch("m1phi",  &fb.m1phi,     "m1phi/D");
    small->Branch("m1bdt",  &fb.m1mvabdt,  "m1bdt/D");
    small->Branch("m1rbdt", &fb.m1rmvabdt, "m1rbdt/D");
    small->Branch("m1q",    &fb.m1q,       "m1q/I");

    small->Branch("m2pt",   &fb.m2pt,      "m2pt/D");
    small->Branch("m2eta",  &fb.m2eta,     "m2eta/D");
    small->Branch("m2phi",  &fb.m2phi,     "m2phi/D");
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
  fHistFile->mkdir(fSample.c_str());
  fHistFile->cd(fSample.c_str());
  TDirectory *dir = gDirectory;

  TH1D *h1(0);
  TH2D *h2(0);

  // vector<string> modifier;
  // modifier.push_back("cnc" + fSuffix);
  // modifier.push_back("bdt" + fSuffix);

  for (unsigned int i = 0; i < fNchan; ++i) {
    h1 = (TH1D*)(fhBdt[i]->Clone(Form("hBdt_%s_chan%d", smode.c_str(), i)));
    h1->SetTitle(Form("hBdt_%s_chan%d %s", smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();
  }

  for (unsigned int im = 0; im < fHistStrings.size(); ++im) {
    for (unsigned int i = 0; i < fNchan; ++i) {
      h1 = (TH1D*)(fhGenAndAccNumbers[fHistStrings[im]][i]->Clone(Form("hGenAndAccNumbers_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();

      h2 = (TH2D*)(fhAccAll[fHistStrings[im]][i]->Clone(Form("hAccAll_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h2->SetDirectory(dir);
      h2->Write();
      h2 = (TH2D*)(fhAccPass[fHistStrings[im]][i]->Clone(Form("hAccPass_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h2->SetDirectory(dir);
      h2->Write();

      h1 = (TH1D*)(fhAccEtaAll[fHistStrings[im]][i]->Clone(Form("hAccEtaAll_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();
      h1 = (TH1D*)(fhAccEtaPass[fHistStrings[im]][i]->Clone(Form("hAccEtaPass_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();


      h1 = (TH1D*)(fhBdtCrossCheck[fHistStrings[im]][i]->Clone(Form("hBdtCrossCheck_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhAccPtAll[fHistStrings[im]][i]->Clone(Form("hAccPtAll_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();
      h1 = (TH1D*)(fhAccPtPass[fHistStrings[im]][i]->Clone(Form("hAccPtPass_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassAbsNoCuts[fHistStrings[im]][i]->Clone(Form("hMassAbsNoCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassAbsNoCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassNoCuts[fHistStrings[im]][i]->Clone(Form("hMassNoCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassNoCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAnaCuts[fHistStrings[im]][i]->Clone(Form("hMassWithAnaCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAnaCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithMuonCuts[fHistStrings[im]][i]->Clone(Form("hMassWithMuonCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithMuonCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithTriggerCuts[fHistStrings[im]][i]->Clone(Form("hMassWithTriggerCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithTriggerCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCuts[fHistStrings[im]][i]->Clone(Form("hMassWithAllCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsBlind[fHistStrings[im]][i]->Clone(Form("hMassWithAllCutsBlind_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsBlind_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsSeagull[fHistStrings[im]][i]->Clone(Form("hMassWithAllCutsSeagull_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsSeagull_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsSeagullBlind[fHistStrings[im]][i]->Clone(Form("hMassWithAllCutsSeagullBlind_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsSeagullBlind_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsCowboy[fHistStrings[im]][i]->Clone(Form("hMassWithAllCutsCowboy_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsCowboy_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhMassWithAllCutsCowboyBlind[fHistStrings[im]][i]->Clone(Form("hMassWithAllCutsCowboyBlind_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithAllCutsCowboyBlind_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhW8MassWithAllCuts[fHistStrings[im]][i]->Clone(Form("hW8MassWithAllCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hW8MassWithAllCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();


      h1 = (TH1D*)(fhW8MassWithAllCutsSeagull[fHistStrings[im]][i]->Clone(Form("hW8MassWithAllCutsSeagull_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hW8MassWithAllCutsSeagull_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhW8MassWithAllCutsCowboy[fHistStrings[im]][i]->Clone(Form("hW8MassWithAllCutsCowboy_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hW8MassWithAllCutsCowboy_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();


      h1 = (TH1D*)(fhMassWithMassCuts[fHistStrings[im]][i]->Clone(Form("hMassWithMassCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hMassWithMassCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      if ((string::npos != fSample.find("bupsik"))
	  || (string::npos != fSample.find("bspsiphi"))
	  || (string::npos != fSample.find("bdpsikstar"))
	  ) {
	h2 = (TH2D*)(fhNorm[fHistStrings[im]][i]->Clone(Form("hNorm_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
	h2->SetTitle(Form("hNorm_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
	h2->SetDirectory(dir);
	h2->Write();

	h2 = (TH2D*)(fhNormC[fHistStrings[im]][i]->Clone(Form("hNormC_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
	h2->SetTitle(Form("hNormC_%s_%s_%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
	h2->SetDirectory(dir);
	h2->Write();

	h2 = (TH2D*)(fhW8Norm[fHistStrings[im]][i]->Clone(Form("hW8Norm_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
	h2->SetTitle(Form("hW8Norm_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
	h2->SetDirectory(dir);
	h2->Write();

	h2 = (TH2D*)(fhW8NormC[fHistStrings[im]][i]->Clone(Form("hW8NormC_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
	h2->SetTitle(Form("hW8NormC_%s_%s_%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
	h2->SetDirectory(dir);
	h2->Write();

      }
    }
  }

  TH1D *hcuts = fDS["bmmData"]->getHist("candAnaMuMu/hcuts");
  hcuts->SetName("hcuts");
  hcuts->SetDirectory(fHistFile);
  cout << "writing hcuts " << hcuts << " name: " << hcuts->GetName() << " to " << fHistFile->GetName() << endl;
  hcuts->Write();

}

// ----------------------------------------------------------------------
void plotResults::resetHistograms(bool deleteThem) {

  // vector<string> modifier;
  // modifier.push_back("cnc" + fSuffix);
  // modifier.push_back("bdt" + fSuffix);
  // for (unsigned int im = 0; im < modifier.size(); ++im) {

  for (unsigned int i = 0; i < fhBdt.size(); ++i) {
    fhBdt[i]->Reset();
    if (deleteThem) delete fhBdt[i];
  }


  for (unsigned int im = 0; im < fHistStrings.size(); ++im) {
    for (int i = 0; i < fNchan; ++i) {
      fhAccAll[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAccAll[fHistStrings[im]][i];
      fhAccPass[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAccPass[fHistStrings[im]][i];

      fhBdtCrossCheck[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhBdtCrossCheck[fHistStrings[im]][i];

      fhAccPtAll[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAccPtAll[fHistStrings[im]][i];
      fhAccPtPass[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAccPtPass[fHistStrings[im]][i];

      fhAccEtaAll[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAccEtaAll[fHistStrings[im]][i];
      fhAccEtaPass[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAccEtaPass[fHistStrings[im]][i];

      fhGenAndAccNumbers[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhGenAndAccNumbers[fHistStrings[im]][i];

      fhMassAbsNoCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassAbsNoCuts[fHistStrings[im]][i];
      fhMassNoCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassNoCuts[fHistStrings[im]][i];

      fhMassWithAnaCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAnaCuts[fHistStrings[im]][i];

      fhMassWithMuonCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithMuonCuts[fHistStrings[im]][i];

      fhMassWithTriggerCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithTriggerCuts[fHistStrings[im]][i];

      fhMassWithAllCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCuts[fHistStrings[im]][i];
      fhMassWithAllCutsBlind[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsBlind[fHistStrings[im]][i];

      fhMassWithAllCutsSeagull[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsSeagull[fHistStrings[im]][i];
      fhMassWithAllCutsSeagullBlind[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsSeagullBlind[fHistStrings[im]][i];

      fhMassWithAllCutsCowboy[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsCowboy[fHistStrings[im]][i];
      fhMassWithAllCutsCowboyBlind[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithAllCutsCowboyBlind[fHistStrings[im]][i];

      fhW8MassWithAllCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhW8MassWithAllCuts[fHistStrings[im]][i];

      fhW8MassWithAllCutsSeagull[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhW8MassWithAllCutsSeagull[fHistStrings[im]][i];

      fhW8MassWithAllCutsCowboy[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhW8MassWithAllCutsCowboy[fHistStrings[im]][i];

      fhMassWithMassCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhMassWithMassCuts[fHistStrings[im]][i];

      fhNorm[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhNorm[fHistStrings[im]][i];
      fhNormC[fHistStrings[im]][i] ->Reset();
      if (deleteThem) delete fhNormC[fHistStrings[im]][i];
      fhW8Norm[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhW8Norm[fHistStrings[im]][i];
      fhW8NormC[fHistStrings[im]][i] ->Reset();
      if (deleteThem) delete fhW8NormC[fHistStrings[im]][i];
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
  n = a.fW8SignalFit;     OUT << "    fW8SignalFit =     " << Form("%7.0f +/- %7.0f +/- %7.0f", n.val,  n.estat, n.esyst) << endl;
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


// ----------------------------------------------------------------------
double plotResults::findVarValue(string varName, vector<string> &lines) {
  for (unsigned int i = 0; i < lines.size(); ++i) {
    if (string::npos != lines[i].find(varName) && string::npos != lines[i].find("val")) {
      string::size_type m1 = lines[i].find("ensuremath{");
      string::size_type m2 = lines[i].find("}", m1);
      string snum = lines[i].substr(m1+12, m2-m1-12-1);
      return atof(snum.c_str());
    }
  }
  return 0.;
}
