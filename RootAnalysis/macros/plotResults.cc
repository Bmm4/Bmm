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
#include "TTimeStamp.h"
#include "ROOT/TProcessExecutor.hxx"
#include "Math/QuantFuncMathMore.h"

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

  fSaveSmallTree = false;

  // -- initialize cuts
  // string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  // cout << "===> Reading cuts from " << cutfile << endl;
  // readCuts(cutfile);
  // fNchan = fCuts.size();

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
  for (unsigned int ih = 0; ih < fHistStrings.size(); ++ih) {
    mode = fHistStrings[ih];
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

      // -- ams
      h = new TH1D(Form("h%sAmMassWithAllCuts%d", mode.c_str(), i), Form("h%sAmMassWithAllCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhAmMassWithAllCuts[mode].push_back(h);
      h = new TH1D(Form("h%sAmW8MassWithAllCuts%d", mode.c_str(), i), Form("h%sAmW8MassWithAllCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhAmW8MassWithAllCuts[mode].push_back(h);
      h = new TH1D(Form("h%sAmMassWithMuonCuts%d", mode.c_str(), i), Form("h%sAmMassWithMuonCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhAmMassWithMuonCuts[mode].push_back(h);

      // -- real (anti-)muon sample
      h = new TH1D(Form("h%sRMMassWithMuonCuts%d", mode.c_str(), i), Form("h%sRMMassWithMuonCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhRMMassWithMuonCuts[mode].push_back(h);
      h = new TH1D(Form("h%sRMMassWithAllCuts%d", mode.c_str(), i), Form("h%sRMMassWithAllCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhRMMassWithAllCuts[mode].push_back(h);

      h = new TH1D(Form("h%sRAmMassWithMuonCuts%d", mode.c_str(), i), Form("h%sRAmMassWithMuonCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhRAmMassWithMuonCuts[mode].push_back(h);
      h = new TH1D(Form("h%sRAmMassWithAllCuts%d", mode.c_str(), i), Form("h%sRAmMassWithAllCuts%d", mode.c_str(), i), NBINS, fMassLo, fMassHi);
      fhRAmMassWithAllCuts[mode].push_back(h);


      if (ih < 2) {
	h2 = new TH2D(Form("h%sMassSysStep0%d", mode.c_str(), i), Form("h%sMassChan%d", mode.c_str(), i), 100, 4.9, 5.9, MAXPS+1, -1., MAXPS);
	fhMassSysStep0[mode].push_back(h2);

	h2 = new TH2D(Form("h%sMassSysStep0C%d", mode.c_str(), i), Form("h%sMassChan%d", mode.c_str(), i), 200, 4.9, 5.9, MAXPS+1, -1., MAXPS);
	fhMassSysStep0C[mode].push_back(h2);

	h2 = new TH2D(Form("h%sMassSysStep1%d", mode.c_str(), i), Form("h%sMassChan%d", mode.c_str(), i), 100, 4.9, 5.9, MAXPS+1, -1., MAXPS);
	fhMassSysStep1[mode].push_back(h2);

	h2 = new TH2D(Form("h%sMassSysStep1C%d", mode.c_str(), i), Form("h%sMassChan%d", mode.c_str(), i), 200, 4.9, 5.9, MAXPS+1, -1., MAXPS);
	fhMassSysStep1C[mode].push_back(h2);
      }
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
    fNoNumbers[i].fAmYield.clear();
    fNoNumbers[i].fFrac.clear();

    fCsNumbers[i].fChan = i;
    fCsNumbers[i].fName = "bspsiphi";
    fCsNumbers[i].fNameMc = "bspsiphiMcComb";
    fCsNumbers[i].fNameDa = "bspsiphiData";
    fCsNumbers[i].fMcYield.clear();
    fCsNumbers[i].fAmYield.clear();
    fCsNumbers[i].fFrac.clear();

    fB0Numbers[i].fChan = i;
    fB0Numbers[i].fName = "bdpsikstar";
    fB0Numbers[i].fNameMc = "bdpsikstarMc";
    fB0Numbers[i].fNameDa = "bdpsikstarData";
    fB0Numbers[i].fMcYield.clear();
    fB0Numbers[i].fAmYield.clear();
    fB0Numbers[i].fFrac.clear();

    fBsmmNumbers[i].fChan = i;
    fBsmmNumbers[i].fName = "bsmm";
    fBsmmNumbers[i].fNameMc = "bsmmMcComb";
    fBsmmNumbers[i].fNameDa = "bmmData";
    fBsmmNumbers[i].fMcYield.clear();
    fBsmmNumbers[i].fAmYield.clear();
    fBsmmNumbers[i].fFrac.clear();

    fBdmmNumbers[i].fChan = i;
    fBdmmNumbers[i].fName = "bdmm";
    fBdmmNumbers[i].fNameMc = "bdmmMcComb";
    fBdmmNumbers[i].fNameDa = "bmmData";
    fBdmmNumbers[i].fMcYield.clear();
    fBdmmNumbers[i].fAmYield.clear();
    fBdmmNumbers[i].fFrac.clear();

    fSlNumbers[i].fChan = i;
    fSlNumbers[i].fName = "sl";
    fSlNumbers[i].fNameMc = "nada";
    fSlNumbers[i].fNameDa = "nada";

    fHhNumbers[i].fChan = i;
    fHhNumbers[i].fName = "hh";
    fHhNumbers[i].fNameMc = "nada";
    fHhNumbers[i].fNameDa = "nada";

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
    fBgNumbers[i].fObsYieldAm.clear();
    fSgAndBgNumbers[i].fMcYield.clear();
    fSgAndBgNumbers[i].fObsYield.clear();
    fSgAndBgNumbers[i].fObsYieldAm.clear();

    fHhNumbers[i].fAmYield.clear();
    fSlNumbers[i].fAmYield.clear();
    fNpNumbers[i].fAmYield.clear();
    fBgNumbers[i].fAmYield.clear();
    fSgAndBgNumbers[i].fAmYield.clear();

    for (int j = 0; j < NWIN; ++j) {
      number aaa;
      fHhNumbers[i].fMcYield.push_back(aaa);
      fHhNumbers[i].fAmYield.push_back(aaa);

      fSlNumbers[i].fMcYield.push_back(aaa);
      fSlNumbers[i].fAmYield.push_back(aaa);

      fNpNumbers[i].fMcYield.push_back(aaa);
      fNpNumbers[i].fAmYield.push_back(aaa);

      fBgNumbers[i].fMcYield.push_back(aaa);
      fBgNumbers[i].fAmYield.push_back(aaa);
      fBgNumbers[i].fObsYield.push_back(aaa);
      fBgNumbers[i].fObsYieldAm.push_back(aaa);

      fSgAndBgNumbers[i].fMcYield.push_back(aaa);
      fSgAndBgNumbers[i].fAmYield.push_back(aaa);
      fSgAndBgNumbers[i].fObsYield.push_back(aaa);
      fSgAndBgNumbers[i].fObsYieldAm.push_back(aaa);

      fNoNumbers[i].fFrac.push_back(aaa);
      fNoNumbers[i].fMcYield.push_back(aaa);
      fNoNumbers[i].fAmYield.push_back(aaa);

      fCsNumbers[i].fFrac.push_back(aaa);
      fCsNumbers[i].fMcYield.push_back(aaa);
      fCsNumbers[i].fAmYield.push_back(aaa);

      fB0Numbers[i].fFrac.push_back(aaa);
      fB0Numbers[i].fMcYield.push_back(aaa);
      fB0Numbers[i].fAmYield.push_back(aaa);

      fBdmmNumbers[i].fFrac.push_back(aaa);
      fBdmmNumbers[i].fMcYield.push_back(aaa);
      fBdmmNumbers[i].fAmYield.push_back(aaa);

      fBsmmNumbers[i].fFrac.push_back(aaa);
      fBsmmNumbers[i].fMcYield.push_back(aaa);
      fBsmmNumbers[i].fAmYield.push_back(aaa);
    }

    fCombNumbers[i].fChan = i;
    fCombNumbers[i].fName = "comb";
    fCombNumbers[i].fNameMc = "nada";
    fCombNumbers[i].fNameDa = "bmmData";
    fCombNumbers[i].fObsYield.clear();
    fCombNumbers[i].fObsYieldAm.clear();
    fCombNumbers[i].fFitYield.clear();
    fCombNumbers[i].fFitYieldAm.clear();
    for (int j = 0; j < NWIN; ++j) {
      number aaa;
      fCombNumbers[i].fObsYield.push_back(aaa);
      fCombNumbers[i].fObsYieldAm.push_back(aaa);
      fCombNumbers[i].fFitYield.push_back(aaa);
      fCombNumbers[i].fFitYieldAm.push_back(aaa);
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
      a->fAmYield.clear();
      a->fObsYield.clear();
      a->fObsYieldAm.clear();
      a->fFitYield.clear();
      a->fFitYieldAm.clear();
      for (int j = 0; j < 5; ++j) {
	a->fFrac.push_back(aaa);
	a->fMcYield.push_back(aaa);
	a->fAmYield.push_back(aaa);
	a->fObsYield.push_back(aaa);
	a->fObsYieldAm.push_back(aaa);
	a->fFitYield.push_back(aaa);
	a->fFitYieldAm.push_back(aaa);
      }
      replaceAll(a->fName, "McOff", "");
      replaceAll(a->fName, "McComb", "");
      replaceAll(a->fName, "Mc", "");
      va.push_back(a);
    }
    fRareNumbers.insert(make_pair(it->first, va));
  }


  // -- define systematics, per channel!
  double fsfu[] = {0.048, 0.048};
  fSystematics["fsfu"] = vector<double>(fsfu, fsfu + sizeof(fsfu)/sizeof(fsfu[0]));
  double bfbupsik[] = {0.030, 0.030};
  fSystematics["bfbupsik"] = vector<double>(bfbupsik, bfbupsik + sizeof(bfbupsik)/sizeof(bfbupsik[0]));
  double bfbspsiphi[] = {0.074, 0.074};
  fSystematics["bfbspsiphi"] = vector<double>(bfbspsiphi, bfbspsiphi + sizeof(bfbspsiphi)/sizeof(bfbspsiphi[0]));

  // -- not needed anymore with UML
  // double pxy[] = {0.05, 0.05};
  // fSystematics["pxy"] = vector<double>(pxy, pxy + sizeof(pxy)/sizeof(pxy[0]));

  double acceptance[] = {0.030, 0.030};
  fSystematics["acceptance"] = vector<double>(acceptance, acceptance + sizeof(acceptance)/sizeof(acceptance[0]));
  double effanabupsik[] = {0.06, 0.07};
  fSystematics["effanabupsik"] = vector<double>(effanabupsik, effanabupsik + sizeof(effanabupsik)/sizeof(effanabupsik[0]));
  double effanabspsiphi[] = {0.05, 0.06};
  fSystematics["effanabspsiphi"] = vector<double>(effanabspsiphi, effanabspsiphi + sizeof(effanabspsiphi)/sizeof(effanabspsiphi[0]));
  double effanabsmm[] = {0.05, 0.06};
  fSystematics["effanabsmm"] = vector<double>(effanabsmm, effanabsmm + sizeof(effanabsmm)/sizeof(effanabsmm[0]));
  double effanabdmm[] = {0.05, 0.06};
  fSystematics["effanabdmm"] = vector<double>(effanabdmm, effanabdmm + sizeof(effanabdmm)/sizeof(effanabdmm[0]));

  double effana4r[fNchan];
  if (fYear == 2016) {
    effana4r[0] = 0.05;
    effana4r[1] = 0.05;
  } else if (fYear == 2012) {
    effana4r[0] = 0.10;
    effana4r[1] = 0.10;
  } else if (fYear == 2011) {
    effana4r[0] = 0.10;
    effana4r[1] = 0.07;
  }
  fSystematics["effana4r"] = vector<double>(effana4r, effana4r + sizeof(effana4r)/sizeof(effana4r[0]));

  double effcandbupsik[] = {0.01, 0.01};
  fSystematics["effcandbupsik"] = vector<double>(effcandbupsik, effcandbupsik + sizeof(effcandbupsik)/sizeof(effcandbupsik[0]));
  double effcandbspsiphi[] = {0.01, 0.01};
  fSystematics["effcandbspsiphi"] = vector<double>(effcandbspsiphi, effcandbspsiphi + sizeof(effcandbspsiphi)/sizeof(effcandbspsiphi[0]));
  double effcandbsmm[] = {0.01, 0.01};
  fSystematics["effcandbsmm"] = vector<double>(effcandbsmm, effcandbsmm + sizeof(effcandbsmm)/sizeof(effcandbsmm[0]));
  double effcandbdmm[] = {0.01, 0.01};
  fSystematics["effcandbdmm"] = vector<double>(effcandbdmm, effcandbdmm + sizeof(effcandbdmm)/sizeof(effcandbdmm[0]));

  double effmuid[fNchan];
  if (fYear == 2016) {
    for (int i = 0; i < fNchan; ++i) effmuid[i] = 0.03;
  } else if (fYear <= 2012) {
    for (int i = 0; i < fNchan; ++i) effmuid[i] = 0.03;
  }
  fSystematics["effmuid"] = vector<double>(effmuid, effmuid + sizeof(effmuid)/sizeof(effmuid[0]));
  double efftrig[] = {0.03, 0.03};
  fSystematics["efftrig"] = vector<double>(efftrig, efftrig + sizeof(efftrig)/sizeof(efftrig[0]));
  double efflifetime[] = {0.03, 0.04};
  fSystematics["efflifetime"] = vector<double>(efflifetime, efflifetime + sizeof(efflifetime)/sizeof(efflifetime[0]));

  // -- pions and kaons
  for (int i = 0; i < fNchan; ++i) effmuid[i] = 0.10;
  fSystematics["fakemuid"] = vector<double>(effmuid, effmuid + sizeof(effmuid)/sizeof(effmuid[0]));
  //  -- protons
  for (int i = 0; i < fNchan; ++i) effmuid[i] = 0.60;
  fSystematics["fakemuidProtons"] = vector<double>(effmuid, effmuid + sizeof(effmuid)/sizeof(effmuid[0]));

  double bunorm[] = {0.04, 0.04};
  fSystematics["normbupsik"] = vector<double>(bunorm, bunorm + sizeof(bunorm)/sizeof(bunorm[0]));
  double bsnorm[] = {0.03, 0.08};
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
  cout << "***************************************************************" << endl;
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);
}


// ----------------------------------------------------------------------
void plotResults::makeAll(string what) {
  if (what == "display") {
    int istart = 20;
    for (int ichan = 0; ichan < 2; ++ichan) {
      for (int i = istart; i < istart+2; ++i) {
	displayScanBDT("all", i, ichan);
      }
    }
  }

  if (what == "show") {
    showScanBDT("all");
  }

  if (what == "ds" ) {
    dumpDatasets();
    printCuts(cout);
  }

  if (what == "bdtopt") {
    init();
    fHistWithAllCuts = "hMassWithAllCuts";
    fillAndSaveHistograms();
    //    fillAndSaveHistograms(0, 10000);
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

  if (what == "dumpdatasets" || what == "genvalidation") {
    init();
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
    init();
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
  if ((what == "all") || (string::npos != what.find("fill"))) {
    //    fillAndSaveHistograms(0, 50000);
    fillAndSaveHistograms();
  }

  if ((what == "all") || (string::npos != what.find("ana"))) {
    init();

    dumpDatasets();
    printCuts(cout);

    fHistWithAllCuts = "hMassWithAllCuts";
    for (int i = 0; i < fNchan; ++i) {
      if (1) {
	fSuffixSel = "cnc" + fSuffix;
	calculateNumbers("cnc" + fSuffix, i);
      }
      fSuffixSel = "bdt" + fSuffix;
      calculateNumbers("bdt" + fSuffix, i);
    }
    if (1) scanBDT(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), true);
  }

  if (string::npos != what.find("cnc")) {
    init();
    fHistWithAllCuts = "hMassWithAllCuts";
    fSuffixSel = "cnc" + fSuffix;
    for (int i = 0; i < fNchan; ++i) {
      calculateNumbers("cnc" + fSuffix, i);
    }
  }

  if (string::npos != what.find("bdt")) {
    fHistWithAllCuts = "hMassWithAllCuts";
    fSuffixSel = "bdt" + fSuffix;
    for (int i = 0; i < fNchan; ++i) {
      calculateNumbers("bdt" + fSuffix, i);
    }
    scanBDT(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), true);
  }


  if ((what == "all") || (string::npos != what.find("ana")) || (string::npos != what.find("sys"))) {

    if (2016 == fYear) {
      sysEffLifetime("WithAnaCuts");
      sysEffLifetime("WithAllCuts");
    }

    // -- now the systematics part
    replaceAll(fTexFileName, "Results", "Systematics");
    cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
    system(Form("/bin/rm -f %s", fTexFileName.c_str()));
    fTEX.open(fTexFileName.c_str(), ios::app);

    for (std::map<std::string, vector<double> >::iterator it = fSystematics.begin(); it != fSystematics.end(); ++it) {
      for (unsigned int i = 0; i < it->second.size(); ++i) {
	fTEX << formatTex(it->second.at(i),  Form("%s:%s_systematics_chan%d:sys", fSetup.c_str(), it->first.c_str(), i), 3) << endl;
      }
    }


    for (int i = 0; i < fNchan; ++i) {
      sysDoubleRatio("bspsiphi", "bupsik", Form("ad%dbdt", i), "muonmbdt", "HLT", fCuts[i]->muonbdt);
      sysDoubleRatio("bdmm", "bupsik", Form("ad%dbdt", i), "muonmbdt", "HLT", fCuts[i]->muonbdt);

      sysNorm("bupsikData", i, 0);
      sysNorm("bspsiphiData", i, 1);

      sysAna("bupsikData", "bspsiphiData", "C", i, 0, 1);
      sysAna("bupsikData", "bspsiphiData", "", i, 0, 1);
    }

    sysBsVsBu();

    scanBDTEffRatio(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), "EFF-TRIG-MC-bupsik", "EFF-TRIG-MC-bsmm");

    fTEX.close();
  }

  if (what == "argh") {
    fSaveSmallTree = false;

    string bac = fHistFileName;
    replaceAll(fHistFileName, "Results", "Lifetime");
    fHistFile = TFile::Open(fHistFileName.c_str(), "RECREATE");
    vector<int> v;
    v.push_back(35);
    v.push_back(40);
    v.push_back(45);
    v.push_back(50);
    v.push_back(55);
    v.push_back(60);
    v.push_back(65);
    v.push_back(66);
    v.push_back(67);
    v.push_back(68);
    v.push_back(69);
    v.push_back(70);
    v.push_back(75);
    v.push_back(80);
    TTree *t(0);

    for (unsigned int i = 0; i < v.size(); ++i) {
      resetHistograms();
      setup(Form("bsmm%dMc", v[i]));
      t = getTree(fSample, fTreeDir);
      setupTree(t, fSample);
      loopOverTree(t, 1);
      saveHistograms(fSample);
    }
    fHistFile->Close();
    fHistFileName = bac;
    return;
  }

  if (what == "dbx1") {
    fHistFile = TFile::Open(fHistFileName.c_str());
    fSuffixSel = "bdt" + fSuffix;
    //    sysBsVsBu();
    //    scanBDTEffRatio(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), "EFF-TRIG-MC-bupsik", "EFF-TRIG-MC-bsmm");


    // scanBDTEffRatio(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), "EFF-CAND-bupsik", "EFF-ANA-bsmm");
    // scanBDTEffRatio(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), "EFF-PRODMC-bupsik", "EFF-ANA-bsmm");
    // scanBDTEffRatio(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), "EFF-TOT-bupsik", "EFF-ANA-bsmm");
    // scanBDTEffRatio(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), "EFF-ANA-bupsik", "EFF-ANA-bsmm");
    // scanBDTEffRatio(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), "EFF-MU-MC-bupsik", "EFF-MU-MC-bsmm");
    // scanBDTEffRatio(Form("%s/scanBDT-%s.tex", fDirectory.c_str(), fSuffix.c_str()), "EFF-TRIG-MC-bupsik", "EFF-TRIG-MC-bsmm");
  }

  if (what == "dbx2") {
    fHistWithAllCuts = "hMassWithAllCuts";
    fSuffixSel = "bdt" + fSuffix;
    for (int i = 0; i < fNchan; ++i) {
      calculateNumbers("bdt" + fSuffix, i);
    }
  }

  if (what == "dbx3") {
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

    // -- and now loop over all requested: invoke calculateNumbers!!
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
    cout << "write out " << hbdt[i]->GetName() << endl;
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
    for (int ic = 0; ic < fNchan; ++ic) {
      double maximum(-1.);
      h1 = new TH1D(Form("bdtScan_%s_chan%d", plots[i].c_str(), ic),
		    Form("bdtScan_%s_chan%d", plots[i].c_str(), ic),
		    80, 0., 80.);
      for (int ib = BDTMIN; ib <= bMax; ++ib) {
	string idx = Form("bdt_%d_", ib);
	fSuffixSel = idx + fSuffix;
	name = Form("%s:%s:chan%d:val", fSuffixSel.c_str(), plots[i].c_str(), ic);
	double val = findVarValue(name, allLines);
	if (val > maximum) {
	  maximum = val;
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

  showScanBDT("all");
}


// ----------------------------------------------------------------------
void plotResults::sysEffLifetime(string cutname) {
  fTEX.open(fTexFileName.c_str(), ios::app);
  string bac = fHistFileName;
  vector<int> v;
  v.push_back(35);
  v.push_back(40);
  v.push_back(45);
  v.push_back(50);
  v.push_back(55);
  v.push_back(60);
  v.push_back(65);
  v.push_back(66);
  v.push_back(67);
  v.push_back(68);
  v.push_back(69);
  v.push_back(70);
  v.push_back(75);
  v.push_back(80);

  double tsm(1.615);
  TH1D *hAll(0), *hCut(0);
  double nAll(0.), nCut(0.);
  string hname("");
  for (int ic = 0; ic < fNchan; ++ic) {
    TH1D *hEff = new TH1D("heff", "", 60, 1.3, 1.9); hEff->Sumw2();
    setTitles(hEff, "#tau_{gen} [ps]", "#varepsilon", 0.05, 1.2, 1.1);
    replaceAll(fHistFileName, "Results", "Lifetime");
    fHistFile = TFile::Open(fHistFileName.c_str());

    for (unsigned int i = 0; i < v.size(); ++i) {
      hname = Form("bsmm%dMc/hMassNoCuts_bdt%s_bsmm%dMc_chan%d", v[i], fSetup.c_str(), v[i], ic);
      //	    cout << "hname = " << hname << endl;
      hAll = (TH1D*)fHistFile->Get(hname.c_str());
      nAll = hAll->Integral();
      hname = Form("bsmm%dMc/hMass%s_bdt%s_bsmm%dMc_chan%d", v[i], cutname.c_str(), fSetup.c_str(), v[i], ic);
      //	    cout << "hname = " << hname << endl;
      hCut = (TH1D*)fHistFile->Get(hname.c_str());
      nCut = hCut->Integral();
      cout << "Efficiency(" << v[i] << ", ic=" << ic << ") = " << nCut/nAll << endl;
      double effV = nCut/nAll;
      double effE = dEff(static_cast<int>(nCut), static_cast<int>(nAll));
      hEff->SetBinContent(hEff->FindBin(1+static_cast<double>(v[i])/100), effV);
      hEff->SetBinError(hEff->FindBin(1+static_cast<double>(v[i])/100), effE);
      fTEX << formatTex(effV,  Form("%s:%s_efftau%d_chan%d:val", cutname.c_str(), fSetup.c_str(), v[i], ic), 3) << endl;
      fTEX << formatTex(effE,  Form("%s:%s_efftau%d_chan%d:err", cutname.c_str(), fSetup.c_str(), v[i], ic), 3) << endl;
    }

    fHistFile->Close();
    fHistFileName = bac;
    fHistFile = TFile::Open(fHistFileName.c_str());
    hAll = (TH1D*)fHistFile->Get(Form("bsmmMcComb/hMassNoCuts_bdt%s_bsmmMcComb_chan%d", fSetup.c_str(), ic));
    nAll = hAll->Integral();
    hCut = (TH1D*)fHistFile->Get(Form("bsmmMcComb/hMass%s_bdt%s_bsmmMcComb_chan%d", cutname.c_str(), fSetup.c_str(), ic));
    nCut = hCut->Integral();
    cout << "Efficiency(default) = " << nCut/nAll << endl;
    double effV = nCut/nAll;
    double effE = dEff(static_cast<int>(nCut), static_cast<int>(nAll));
    hEff->SetBinContent(hEff->FindBin(1.471), effV);
    hEff->SetBinError(hEff->FindBin(1.471), effE);
    fTEX << formatTex(effV,  Form("%s:%s_efftauDef_chan%d:val", cutname.c_str(), fSetup.c_str(), ic), 3) << endl;
    fTEX << formatTex(effE,  Form("%s:%s_efftauDef_chan%d:err", cutname.c_str(), fSetup.c_str(), ic), 3) << endl;
    fHistFile->Close();

    //      hEff->Draw();
    double ymax = 1.1*hEff->GetBinContent(hEff->FindBin(1.8));
    hEff->SetMaximum(ymax);
    hEff->SetMinimum(0.001);
    hEff->Fit("pol1");
    double p0 = hEff->GetFunction("pol1")->GetParameter(0);
    double p1 = hEff->GetFunction("pol1")->GetParameter(1);
    tl->SetTextSize(0.05);
    tl->SetTextColor(kBlack);
    tl->SetTextAngle(90.);
    pl->SetLineColor(kBlack);
    pl->SetLineStyle(kDotted);
    pl->DrawLine(1.415, 0., 1.415, ymax);
    tl->DrawLatex(1.405, 0.2*ymax, "(#Gamma_{L}, A = -1)");
    pl->SetLineStyle(kSolid);
    pl->DrawLine(tsm, 0., tsm, ymax);
    tl->DrawLatex(1.605, 0.2*ymax, "SM (#Gamma_{H}, A = +1)");
    pl->SetLineStyle(kDotted);
    pl->DrawLine(0.5*(tsm+1.415), 0., 0.5*(tsm+1.415), ymax);
    tl->DrawLatex(1.505, 0.2*ymax, "(A = 0)");
    pl->SetLineColor(kRed);
    tl->SetTextColor(kRed);
    pl->SetLineStyle(kSolid);
    pl->DrawLine(1.475, 0., 1.475, ymax);
    tl->DrawLatex(1.47, 0.05*ymax, "CMS MC");
    tl->SetTextAngle(0.);
    tl->SetTextColor(kBlack);
    double ecms = hEff->GetBinContent(hEff->FindBin(1.475));
    double esmh = p0 + p1*tsm;
    double esml = p0 + p1*1.415;
    double esm0 = p0 + p1*1.515;
    tl->SetTextSize(0.03);
    tl->DrawLatex(1.34, 0.15*ymax, Form("#varepsilon(fit) = %4.3f", esml));
    tl->SetTextColor(kRed);
    tl->DrawLatex(1.48, 0.03*ymax, Form("#varepsilon(CMS) = %4.3f", ecms));
    tl->SetTextColor(kBlack);
    tl->DrawLatex(1.52, 0.15*ymax, Form("#varepsilon(fit)  = %4.3f", esm0));
    tl->DrawLatex(1.62, 0.15*ymax, Form("#varepsilon(fit)  = %4.3f", esmh));
    tl->DrawLatex(1.62, 0.11*ymax, Form("#Delta(rel) = %4.3f (compared to CMS MC)", (esmh-ecms)/esmh));
    tl->DrawLatex(1.62, 0.07*ymax, Form("#Delta(rel) = %4.3f (compared to #Gamma_{L}, A=-1)", (esmh-esml)/esmh));
    tl->DrawLatex(1.62, 0.03*ymax, Form("#Delta(rel) = %4.3f (compared to A=0)", (esmh-esm0)/esmh));

    fTEX << formatTex(esmh/ecms, Form("%s:%s_efftauScale_chan%d:val", cutname.c_str(), fSetup.c_str(), ic), 3) << endl;
    fTEX << formatTex(ecms,  Form("%s:%s_efftauCMS_chan%d:val", cutname.c_str(), fSetup.c_str(), ic), 3) << endl;
    fTEX << formatTex(esmh,  Form("%s:%s_efftauHeavy_chan%d:val", cutname.c_str(), fSetup.c_str(), ic), 3) << endl;
    fTEX << formatTex(esml,  Form("%s:%s_efftauLight_chan%d:val", cutname.c_str(), fSetup.c_str(), ic), 3) << endl;
    fTEX << formatTex(esm0,  Form("%s:%s_efftauA0_chan%d:val", cutname.c_str(), fSetup.c_str(), ic), 3) << endl;

    fTEX << formatTex((esmh-esml)/esmh/TMath::Sqrt(12),  Form("%s:%s_efftau_chan%d:sys", cutname.c_str(), fSetup.c_str(), ic), 3) << endl;

    savePad(Form("sys/sysEffLifetime_%s_%s_chan%d.pdf", fSetup.c_str(), cutname.c_str(), ic));
    delete hEff;
  }
  fTEX.close();

}




// ----------------------------------------------------------------------
void plotResults::sysBsVsBu() {
  // -- read in all tex files
  map<string, string> files;
  files.insert(make_pair("2011s01", Form("%s/plotResults.2011s01.tex", fDirectory.c_str())));
  files.insert(make_pair("2012s01", Form("%s/plotResults.2012s01.tex", fDirectory.c_str())));
  files.insert(make_pair("2016BFs01", Form("%s/plotResults.2016BFs01.tex", fDirectory.c_str())));
  files.insert(make_pair("2016GHs01", Form("%s/plotResults.2016GHs01.tex", fDirectory.c_str())));
  vector<string> allLines;
  char buffer[2000];

  map<string, pair<double, double> > bfbs;
  string name;
  for (map<string, string>::iterator it = files.begin(); it != files.end(); ++it) {
    allLines.clear();
    ifstream is(it->second);
    while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
    is.close();
    for (int ichan = 0; ichan < fNchan; ++ichan) {
      name = Form("%s:CSBF:chan%d", fSuffixSel.c_str(), ichan);
      double val = findVarValue(name, allLines);
      double err = findVarEstat(name, allLines);
      name = Form("%s_chan%d", it->first.c_str(), ichan);
      cout << name << ": " << val << " +/- " << err << endl;
      bfbs.insert(make_pair(name, make_pair(val, err)));
    }
  }

  // -- make chanX vectors and overall vector
  vector<double> chan0V, chan0E, chan1V, chan1E, allV, allE;
  for (map<string, pair<double, double> >::iterator it = bfbs.begin(); it != bfbs.end(); ++it) {
    if (string::npos != it->first.find("chan0")) {
      chan0V.push_back(it->second.first);
      chan0E.push_back(it->second.second);
      allV.push_back(it->second.first);
      allE.push_back(it->second.second);
    }
    if (string::npos != it->first.find("chan1")) {
      chan1V.push_back(it->second.first);
      chan1E.push_back(it->second.second);
      allV.push_back(it->second.first);
      allE.push_back(it->second.second);
    }
  }

  double val0(0.), err0(0.), val1(0.), err1(0.), valA(0.), errA(0.), chi0(0.), chi1(0.), chiA(0.);
  average(val0, err0, chan0V, chan0E, chi0);
  average(val1, err1, chan1V, chan1E, chi1);
  average(valA, errA, allV, allE, chiA);


  double temp(0.);
  for (unsigned int i = 0; i < allV.size(); ++i) {
    temp += (allV[i] - valA)*(allV[i] - valA);
  }
  double rms = TMath::Sqrt(temp/(allV.size() - 1));


  string setup("meta");
  cout << formatTex(val0, Form("%s:CSBF_COMB_chan0:val", setup.c_str()), 7) << endl;
  cout << formatTex(err0, Form("%s:CSBF_COMB_chan0:err", setup.c_str()), 7) << endl;
  cout << formatTex(chi0, Form("%s:CSBF_COMB_chan0:chi", setup.c_str()), 7) << endl;
  cout << formatTex(val1, Form("%s:CSBF_COMB_chan1:val", setup.c_str()), 7) << endl;
  cout << formatTex(err1, Form("%s:CSBF_COMB_chan1:err", setup.c_str()), 7) << endl;
  cout << formatTex(chi1, Form("%s:CSBF_COMB_chan1:chi", setup.c_str()), 7) << endl;
  cout << formatTex(valA, Form("%s:CSBF_COMB_combi:val", setup.c_str()), 7) << endl;
  cout << formatTex(errA, Form("%s:CSBF_COMB_combi:err", setup.c_str()), 7) << endl;
  cout << formatTex(chiA, Form("%s:CSBF_COMB_combi:chi", setup.c_str()), 7) << endl;
  cout << formatTex(rms,  Form("%s:CSBF_COMB_combi:rms", setup.c_str()), 7) << endl;
  cout << formatTex(rms/valA,  Form("%s:CSBF_COMB_combi:sys", setup.c_str()), 4) << endl;

}


// ----------------------------------------------------------------------
// -- calculate double ratio [eps(mu|bspsiphiMc)/eps(mu|bupsiMc)] / [eps(mu|bspsiphiData)/eps(mu|bupsiData)] from integrals of signal mass fits
void plotResults::sysNorm(std::string sample1, int ichan, int version) {

  string hfname = fHistFileName;
  cout << "open " << hfname << endl;
  TFile *f1 = TFile::Open(hfname.c_str());

  string nc(Form("hMassSysStep1_bdt%s", fSetup.c_str()));
  string wc(Form("hMassSysStep1C_bdt%s", fSetup.c_str()));

  string d1name, d2name, m1name, m2name;
  TH2D *h2NC(0), *h2WC(0);
  TH1D *d1NC(0), *d1WC(0);
  h2NC = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample1.c_str(), nc.c_str(), sample1.c_str(), ichan));
  cout << Form("%s/%s_%s_chan%d", sample1.c_str(), nc.c_str(), sample1.c_str(), ichan) << ": " << h2NC << endl;
  d1NC = h2NC->ProjectionX(Form("proj_%s_%s_chan%d", nc.c_str(), sample1.c_str(), ichan), 2, 2);
  h2WC = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample1.c_str(), wc.c_str(), sample1.c_str(), ichan));
  cout << Form("%s/%s_%s_chan%d", sample1.c_str(), wc.c_str(), sample1.c_str(), ichan) << ": " << h2WC << endl;
  d1WC = h2WC->ProjectionX(Form("proj_%s_%s_chan%d", wc.c_str(), sample1.c_str(), ichan), 2, 2);


  double xmin(5.2), sigma(0.01);

  double ncSgVal(0.), wcSgVal(0.),  ncSgErr(0.), wcSgErr(0.);
  fitPsYield a(0, fVerbose);
  string sname = fDirectory + string("/sys/sysNorm") + fSetup;
  if (string::npos != sample1.find("bspsiphi")) {
    xmin = 5.1;
    sigma = 0.030;
    if (0 == version) {
      a.fit0_Bs2JpsiPhi(d1NC, -1, sname, xmin, 5.80, sigma + ichan*0.015);
    } else if (1 == version) {
      a.fit1_Bs2JpsiPhi(d1NC, -1, sname, xmin, 5.80, sigma + ichan*0.015);
    }
    ncSgVal = a.getSignalUnW8Yield();
    ncSgErr = a.getSignalUnW8Error();
    xmin = 5.2;
    sigma = 0.012;
    if (1 == ichan) {
      if (fYear < 2015) {
	sigma = 0.020;
      } else {
	sigma = 0.015;
      }
    }
    a.fit1_Bs2JpsiPhi(d1WC, -1, sname, xmin, 5.80, sigma);
    wcSgVal = a.getSignalUnW8Yield();
    wcSgErr = a.getSignalUnW8Error();
  }

  if (string::npos != sample1.find("bupsik")) {
    xmin = 4.95; // 4% OK
    xmin = 5.00; // 4% OK
    sigma = 0.030;
    a.fit0_Bu2JpsiKp(d1NC, -1, sname, xmin, 5.80, sigma + ichan*0.015);
    ncSgVal = a.getSignalUnW8Yield();
    ncSgErr = a.getSignalUnW8Error();
    xmin = 5.1; // 5.5%?
    xmin = 5.15; // 5%
    xmin = 5.05; // not a good choice
    xmin = 5.00; // 4% OK but bad chi2
    sigma = 0.010;
    if (1 == ichan) {
      if (fYear < 2015) {
	sigma = 0.020;
      } else {
	sigma = 0.015;
      }
    }
    a.fit0_Bu2JpsiKp(d1WC, -1, sname, xmin, 5.80, sigma);
    wcSgVal = a.getSignalUnW8Yield();
    wcSgErr = a.getSignalUnW8Error();
  }

  double sysV = (ncSgVal-wcSgVal)/ncSgVal;
  double sysE = dRatio(wcSgVal, wcSgErr, ncSgVal, ncSgErr);

  replaceAll(sample1, "Data", "");
  string sysname = Form("%s:sysnorm_chan%d_%s", fSetup.c_str(), ichan, sample1.c_str());
  fTEX << formatTex(wcSgVal,  Form("%s_wconst:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(wcSgErr,  Form("%s_wconst:err", sysname.c_str()), 3) << endl;
  fTEX << formatTex(ncSgVal,  Form("%s_nconst:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(ncSgErr,  Form("%s_nconst:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(sysV,  Form("%s_nconst:sys", sysname.c_str()), 3) << endl;
  fTEX << formatTex(sysE,  Form("%s_nconst:sysErr", sysname.c_str()), 3) << endl;


}


// ----------------------------------------------------------------------
// -- calculate double ratio [eps(mu|bspsiphiMc)/eps(mu|bupsiMc)] / [eps(mu|bspsiphiData)/eps(mu|bupsiData)] from integrals of signal mass fits
void plotResults::sysAna(std::string sample1, std::string sample2, string massc, int ichan, int buversion, int bsversion) {
  c0->SetCanvasSize(700, 700);
  c0->cd();
  shrinkPad(0.15, 0.18);

  string hfname = fHistFileName;
  cout << "open " << hfname << endl;
  TFile *f1 = TFile::Open(hfname.c_str());

  string strAll(Form("hMassSysStep0%s_bdt%s", massc.c_str(), fSetup.c_str()));
  string strCut(Form("hMassSysStep1%s_bdt%s", massc.c_str(), fSetup.c_str()));

  string d1name, d2name, m1name, m2name;
  TH2D *h2d1A(0), *h2d1C(0), *h2d2A(0), *h2d2C(0);
  TH1D *d1A(0), *d1C(0), *d2A(0), *d2C(0);
  TH2D *h2m1A(0), *h2m1C(0), *h2m2A(0), *h2m2C(0);
  TH1D *m1A(0), *m1C(0), *m2A(0), *m2C(0);
  // hMassWithAnaCuts_bdt2011s01_bupsikData_chan1;1
  h2d1A = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample1.c_str(), strAll.c_str(), sample1.c_str(), ichan));
  cout << Form("%s/%s_%s_chan%d", sample1.c_str(), strAll.c_str(), sample1.c_str(), ichan) << ": " << h2d1A << endl;
  d1A = h2d1A->ProjectionX(Form("proj_%s_%s_chan%d", strAll.c_str(), sample1.c_str(), ichan), 2, 2);
  h2d1C = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample1.c_str(), strCut.c_str(), sample1.c_str(), ichan));
  cout << Form("%s/%s_%s_chan%d", sample1.c_str(), strCut.c_str(), sample1.c_str(), ichan) << ": " << h2d1C << endl;
  d1C = h2d1C->ProjectionX(Form("proj_%s_%s_chan%d", strCut.c_str(), sample1.c_str(), ichan), 2, 2);

  h2d2A = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample2.c_str(), strAll.c_str(), sample2.c_str(), ichan));
  cout << Form("%s/%s_%s_chan%d", sample2.c_str(), strAll.c_str(), sample2.c_str(), ichan) << ": " << h2d2A << endl;
  d2A = h2d2A->ProjectionX(Form("proj_%s_%s_chan%d", strAll.c_str(), sample2.c_str(), ichan), 2, 2);
  h2d2C = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample2.c_str(), strCut.c_str(), sample2.c_str(), ichan));
  cout << Form("%s/%s_%s_chan%d", sample2.c_str(), strCut.c_str(), sample2.c_str(), ichan) << ": " << h2d2C << endl;
  d2C = h2d2C->ProjectionX(Form("proj_%s_%s_chan%d", strCut.c_str(), sample2.c_str(), ichan), 2, 2);

  replaceAll(sample1, "Data", "McComb");
  replaceAll(sample2, "Data", "McComb");
  h2m1A = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample1.c_str(), strAll.c_str(), sample1.c_str(), ichan));
  m1A = h2m1A->ProjectionX(Form("proj_%s_%s_chan%d", strAll.c_str(), sample1.c_str(), ichan), 2, 2);
  double m1ASgVal = m1A->GetSumOfWeights();
  h2m1C = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample1.c_str(), strCut.c_str(), sample1.c_str(), ichan));
  m1C = h2m1C->ProjectionX(Form("proj_%s_%s_chan%d", strAll.c_str(), sample1.c_str(), ichan), 2, 2);
  double m1CSgVal = m1C->GetSumOfWeights();
  h2m2A = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample2.c_str(), strAll.c_str(), sample2.c_str(), ichan));
  m2A = h2m2A->ProjectionX(Form("proj_%s_%s_chan%d", strAll.c_str(), sample2.c_str(), ichan), 2, 2);
  double m2ASgVal = m2A->GetSumOfWeights();
  h2m2C = (TH2D*)f1->Get(Form("%s/%s_%s_chan%d", sample2.c_str(), strCut.c_str(), sample2.c_str(), ichan));
  m2C = h2m2C->ProjectionX(Form("proj_%s_%s_chan%d", strCut.c_str(), sample2.c_str(), ichan), 2, 2);
  double m2CSgVal = m2C->GetSumOfWeights();

  if (1) {
    c0->Divide(2,2);
    c0->cd(1);  d1A->Draw(); d1C->Draw("samehist");
    c0->cd(2);  d2A->Draw(); d2C->Draw("samehist");
    c0->cd(3);  m1A->Draw(); m1C->Draw("samehist");
    c0->cd(4);  m2A->Draw(); m2C->Draw("samehist");
    savePad("bla.pdf", c0);
  }

  double xmin(5.2), sigma(0.01);
  if (string::npos == massc.find("C")) {
    // NOT the mass-constrained version
    sigma = 0.025;
    xmin = 5.1;
  }

  fitPsYield a(0, fVerbose);
  string version(Form("s%du%d", bsversion, buversion));
  string sname = fDirectory + string("/sys/sysAna") + fSetup;
  if (0 == bsversion) {
    a.fit0_Bs2JpsiPhi(d2A, -1, sname, xmin, 5.80, sigma + ichan*0.015);
  } else if (1 == bsversion) {
    a.fit1_Bs2JpsiPhi(d2A, -1, sname, xmin, 5.80, sigma + ichan*0.015);
  }
  double d2ASgVal = a.getSignalUnW8Yield();
  double chi2dof1 = a.getSignalUnW8Chi2();
  if (0 == bsversion) {
    a.fit0_Bs2JpsiPhi(d2C, -1, sname, xmin, 5.80, sigma + ichan*0.015);
  } else if (1 == bsversion) {
    a.fit1_Bs2JpsiPhi(d2C, -1, sname, xmin, 5.80, sigma + ichan*0.015);
  }
  double d2CSgVal = a.getSignalUnW8Yield();
  double chi2dof2 = a.getSignalUnW8Chi2();
  double ndof     = a.getSignalUnW8Ndof();

  // a.fit1_Bs2JpsiPhi(m1C, -1, sname, xmin, 5.80, sigma + ichan*0.005);
  // double m1CSgFit = a.getSignalUnW8Yield();
  // a.fit1_Bs2JpsiPhi(m1A, -1, sname, xmin, 5.80, sigma + ichan*0.005);
  // double m1ASgFit = a.getSignalUnW8Yield();

  xmin = 5.05;
  if (string::npos == massc.find("C")) {
    // NOT the mass-constrained version
    sigma = 0.025;
    xmin = 5.0;
  }

  fitPsYield b(0, fVerbose);
  if (0 == buversion) {
    b.fit0_Bu2JpsiKp(d1A, -1, sname, xmin, 5.80, sigma + ichan*0.005);
  } else if (1 == buversion) {
    b.fit1_Bu2JpsiKp(d1A, -1, sname, xmin, 5.80, sigma + ichan*0.005);
  }
  double d1ASgVal = b.getSignalUnW8Yield();
  if (0 == buversion) {
    b.fit0_Bu2JpsiKp(d1C, -1, sname, xmin, 5.80, sigma + ichan*0.005);
  } else if (1 == buversion) {
    b.fit1_Bu2JpsiKp(d1C, -1, sname, xmin, 5.80, sigma + ichan*0.005);
  }
  double d1CSgVal = b.getSignalUnW8Yield();

  // b.fit0_Bu2JpsiKp(m2C, -1, sname, xmin, 5.80, sigma + ichan*0.005);
  // b.fit0_Bu2JpsiKp(m2A, -1, sname, xmin, 5.80, sigma + ichan*0.005);


  double eps1Data  = d1CSgVal/d1ASgVal;
  double eps1DataE = dEff(static_cast<int>(d1CSgVal), static_cast<int>(d1ASgVal));
  double eps2Data  = d2CSgVal/d2ASgVal;
  double eps2DataE = dEff(static_cast<int>(d2CSgVal), static_cast<int>(d2ASgVal));
  double ratioData = eps1Data/eps2Data;
  double ratioDataE= dRatio(eps1Data, eps1DataE, eps2Data, eps2DataE);
  double eps1Mc    = m1CSgVal/m1ASgVal;
  double eps1McE   = dEff(static_cast<int>(m1CSgVal), static_cast<int>(m1ASgVal));
  double eps2Mc    = m2CSgVal/m2ASgVal;
  double eps2McE   = dEff(static_cast<int>(m2CSgVal), static_cast<int>(m2ASgVal));
  double ratioMc   = eps1Mc/eps2Mc;
  double ratioMcE  = dRatio(eps1Mc, eps1McE, eps2Mc, eps2McE);

  cout << "Data: bs: " << d1CSgVal << "/" << d1ASgVal << " = " << d1CSgVal/d1ASgVal  << " bu: " << d2CSgVal << "/" << d2ASgVal << " = " << d2CSgVal/d2ASgVal << endl;
  cout << "MC:   bs: " << m1CSgVal << "/" << m1ASgVal << " = " << m1CSgVal/m1ASgVal  << " bu: " << m2CSgVal << "/" << m2ASgVal << " = " << m2CSgVal/m2ASgVal << endl;
  cout << "Data eps1/eps2 = " << eps1Data << "/" << eps2Data << " = " << ratioData << endl;
  cout << "MC   eps1/eps2 = " << eps1Mc   << "/" << eps2Mc   << " = " << ratioMc << endl;
  cout << "systematics: " << 1. - (ratioData/ratioMc) << " chi2/ndf = " << chi2dof2/ndof << " _ " << chi2dof1/ndof << endl;

  replaceAll(sample1, "McComb", "");
  replaceAll(sample2, "McComb", "");
  string vmassc = massc + version;
  string sysname = Form("%s:sysana%s_chan%d_%s_%s", fSetup.c_str(), vmassc.c_str(), ichan, sample1.c_str(),  sample2.c_str());
  fTEX << formatTex(eps1Mc,  Form("%s_Eps1Mc:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(eps2Mc,  Form("%s_Eps2Mc:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(eps1Data,  Form("%s_Eps1Dt:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(eps2Data,  Form("%s_Eps2Dt:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(ratioMc,  Form("%s_ratioEpsMc:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(ratioMcE, Form("%s_ratioEpsMc:err", sysname.c_str()), 3) << endl;
  fTEX << formatTex(ratioData,  Form("%s_ratioEpsDt:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(ratioDataE, Form("%s_ratioEpsDt:err", sysname.c_str()), 3) << endl;

  // -- write out the systematic on the double ratio
  fTEX << formatTex((ratioData/ratioMc),  Form("%s_doubleRatioEps:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(dRatio(ratioData, ratioDataE, ratioMc, ratioMcE),  Form("%s_doubleRatioEps:err", sysname.c_str()), 3) << endl;
  fTEX << formatTex(1.-(ratioData/ratioMc),  Form("%s_doubleRatioEps:sys", sysname.c_str()), 3) << endl;
  fTEX << formatTex(dRatio(ratioData, ratioDataE, ratioMc, ratioMcE),  Form("%s_doubleRatioEps:sysErr", sysname.c_str()), 3) << endl;

  // -- write out the systematic on the single ratios data/MC
  sysname = Form("%s:sysana%s_chan%d_%s", fSetup.c_str(), vmassc.c_str(), ichan, sample1.c_str());
  fTEX << formatTex(1.-(eps1Data/eps1Mc),  Form("%s_singleRatioEps:sys", sysname.c_str()), 3) << endl;
  fTEX << formatTex(dRatio(eps1Data, eps1DataE, eps1Mc, eps1McE),  Form("%s_singleRatioEps:sysErr", sysname.c_str()), 3) << endl;

  sysname = Form("%s:sysana%s_chan%d_%s", fSetup.c_str(), vmassc.c_str(), ichan, sample2.c_str());
  fTEX << formatTex(1.-(eps2Data/eps2Mc),  Form("%s_singleRatioEps:sys", sysname.c_str()), 3) << endl;
  fTEX << formatTex(dRatio(eps2Data, eps2DataE, eps2Mc, eps2McE),  Form("%s_singleRatioEps:sysErr", sysname.c_str()), 3) << endl;

}


// ----------------------------------------------------------------------
// -- calculate double ratio [eps(mu|bupsikMc)/eps(mu|bspsiphiMc)] / [eps(mu|bupsikData)/eps(mu|bspsiphiData)] from sideband-subtracted distribution
void plotResults::sysDoubleRatio(string sample1, string sample2, string chansel, string var, string cutlevel, double cut) {

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
  fTEX << formatTex((depsd/depsm),  Form("%s_ratioEps:val", sysname.c_str()), 3) << endl;
  fTEX << formatTex(1.-(depsd/depsm),  Form("%s_ratioEps:sys", sysname.c_str()), 3) << endl;

  savePad(Form("sys/sysDoubleRatio%s_%s_%s_%s_%s_chan%d.pdf", fSetup.c_str(), var.c_str(), sample1.c_str(), sample2.c_str(), cutlevel.c_str(), ichan));


}




//----------------------------------------------------------------------
void plotResults::scanBDTEffRatio(string fname, string string1, string string2) {
  cout << "starting scanBDTEffRatio, string1 = " << string1 << " string2 = " << string2 << endl;
  int BDTMIN(0);
  double bdtCutOffHi(-1.), bdtCutOffLo(-1.);
  fHistFile = TFile::Open(fHistFileName.c_str());
  TH1D *hbdt(0);

  // -- now read in tex file
  vector<string> allLines;
  ifstream is(fname);
  char buffer[2000];
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  is.close();
  string name1, name2;
  TH1D *h1(0);
  for (int ic = 0; ic < fNchan; ++ic) {
    hbdt = (TH1D*)fHistFile->Get(Form("bsmmMcComb/hBdt_bsmmMcComb_chan%d", ic));
    double btot = hbdt->Integral(hbdt->FindBin(fCuts[ic]->bdtCut), hbdt->GetNbinsX());
    for (int ibin = hbdt->FindBin(fCuts[ic]->bdtCut); ibin < hbdt->GetNbinsX(); ++ibin) {
      if (hbdt->Integral(ibin, hbdt->GetNbinsX()) < 0.5*btot) {
	bdtCutOffHi = hbdt->GetBinCenter(ibin);
	bdtCutOffLo = hbdt->GetBinCenter(hbdt->FindBin(fCuts[ic]->bdtCut) - (ibin - hbdt->FindBin(fCuts[ic]->bdtCut)));
	break;
      }
    }
    cout << "btot = " << btot << ", bdtCutOffHi = " << bdtCutOffHi << ", bdtCutOffLo = " << bdtCutOffLo << endl;
    double ratio(-1.), val0(-1.), valMin(999.), valMax(-999.), valCutHi(-999.), valCutLo(-999.);
    h1 = new TH1D(Form("bdtScanEffRatio_%s_%s_chan%d", string1.c_str(), string2.c_str(), ic),
		  Form("bdtScanEffRatio_%s_%s_chan%d", string1.c_str(), string2.c_str(), ic),
		  80, 0., 80.);
    for (unsigned int ib = BDTMIN; ib <= 80; ++ib) {

      string idx = Form("bdt_%d_", ib);
      fSuffixSel = idx + fSuffix;
      name1 = Form("%s:%s-chan%d", fSuffixSel.c_str(), string1.c_str(), ic);
      double err1 = findVarEstat(name1, allLines);
      name1 = Form("%s:%s-chan%d", fSuffixSel.c_str(), string1.c_str(), ic);
      double val1 = findVarValue(name1, allLines);
      name2 = Form("%s:%s-chan%d", fSuffixSel.c_str(), string2.c_str(), ic);
      double err2 = findVarEstat(name2, allLines);
      name2 = Form("%s:%s-chan%d", fSuffixSel.c_str(), string2.c_str(), ic);
      double val2 = findVarValue(name2, allLines);
      if (0) cout << "  " << name1 << ": " << val1 << " +/- " << err1
		  << "  " << name2 << ": " << val2 << " +/- " << err2
		  << " filling into bin " << h1->FindBin(ib)
		  << endl;
      if (val2 > 0) {
	ratio = val1/val2;
      } else {
	ratio = 0.;
      }
      h1->SetBinContent(h1->FindBin(ib), ratio);
      if (val1 < 1.e-6 && val2 < 1.e-6) break;
      if (h1->GetBinLowEdge(ib) < 100*fCuts[ic]->bdtCut &&  100*fCuts[ic]->bdtCut < h1->GetBinLowEdge(ib+1)) {
	val0 = ratio;
      }
      if (h1->GetBinLowEdge(ib) < 100*bdtCutOffHi &&  100*bdtCutOffHi < h1->GetBinLowEdge(ib+1)) {
	valCutHi = ratio;
      }
      if (h1->GetBinLowEdge(ib) < 100*bdtCutOffLo &&  100*bdtCutOffLo < h1->GetBinLowEdge(ib+1)) {
	valCutLo = ratio;
      }
      if (ratio < valMin) valMin = ratio;
      if (ratio > valMax) valMax = ratio;
    }
    setTitles(h1, "100 #times BDT >", Form("%s/%s", string1.c_str(), string2.c_str()), 0.05, 1.2, 1.1);
    h1->Draw();
    pl->SetLineWidth(2); pl->SetLineColor(kRed);
    pl->DrawLine(100*fCuts[ic]->bdtCut, 1.05*h1->GetMaximum(), 100*fCuts[ic]->bdtCut, 0.0);
    pl->SetLineWidth(2); pl->SetLineColor(kMagenta);
    pl->DrawLine(100*bdtCutOffHi, 1.05*h1->GetMaximum(), 100*bdtCutOffHi, 0.0);
    pl->SetLineWidth(2); pl->SetLineColor(kCyan);
    pl->DrawLine(100*bdtCutOffLo, 1.05*h1->GetMaximum(), 100*bdtCutOffLo, 0.0);
    tl->SetTextSize(0.041);
    tl->SetTextColor(kBlack);
    tl->DrawLatexNDC(0.25, 0.85, Form("r_{max}")); tl->DrawLatexNDC(0.30, 0.85, Form("= %4.3f", valMax));
    tl->DrawLatexNDC(0.25, 0.80, Form("r_{def}")); tl->DrawLatexNDC(0.30, 0.80, Form("= %4.3f", val0));
    tl->DrawLatexNDC(0.25, 0.75, Form("r_{min}")); tl->DrawLatexNDC(0.30, 0.75, Form("= %4.3f", valMin));
    tl->DrawLatexNDC(0.25, 0.70, Form("r_{+50%%}")); tl->DrawLatexNDC(0.30, 0.70, Form("= %4.3f", valCutHi));
    tl->DrawLatexNDC(0.25, 0.65, Form("r_{-50%%}")); tl->DrawLatexNDC(0.30, 0.65, Form("= %4.3f", valCutLo));

    tl->DrawLatexNDC(0.25, 0.55, Form("#Delta(+)")); tl->DrawLatexNDC(0.45, 0.55, Form("= %4.3f", (valMax-val0)/val0));
    tl->DrawLatexNDC(0.25, 0.50, Form("#Delta(-)")); tl->DrawLatexNDC(0.45, 0.50, Form("= %4.3f", (val0-valMin)/val0));

    tl->DrawLatexNDC(0.25, 0.45, Form("#Delta(+50%%)")); tl->DrawLatexNDC(0.45, 0.45, Form("= %4.3f", TMath::Abs((valCutHi-val0))/val0));
    tl->DrawLatexNDC(0.25, 0.40, Form("#Delta(-50%%)")); tl->DrawLatexNDC(0.45, 0.40, Form("= %4.3f", TMath::Abs((valCutLo-val0))/val0));

    tl->DrawLatexNDC(0.55, 0.92, Form("%s/chan%d", fSetup.c_str(), ic));

    savePad(Form("sys/sysScanBDTEffRatio_%s-%s-%s-chan%d.pdf", fSetup.c_str(), string1.c_str(), string2.c_str(), ic));

    // -- write out the (maximum) systematics
    string sysname = Form("%s:sysTriggerScanHi_%s_%s", fSetup.c_str(), string1.c_str(), string2.c_str());
    fTEX << formatTex(TMath::Abs((valCutHi-val0))/val0,  Form("%s_chan%d:sys", sysname.c_str(), ic), 3) << endl;
    sysname = Form("%s:sysTriggerScanLo_%s_%s", fSetup.c_str(), string1.c_str(), string2.c_str());
    fTEX << formatTex(TMath::Abs((valCutLo-val0))/val0,  Form("%s_chan%d:sys", sysname.c_str(), ic), 3) << endl;
    sysname = Form("%s:sysTriggerScanMax_%s_%s", fSetup.c_str(), string1.c_str(), string2.c_str());
    fTEX << formatTex(TMath::Abs((valMax-val0))/val0,  Form("%s_chan%d:sys", sysname.c_str(), ic), 3) << endl;
    sysname = Form("%s:sysTriggerScanMin_%s_%s", fSetup.c_str(), string1.c_str(), string2.c_str());
    fTEX << formatTex(TMath::Abs((val0-valMin))/val0,  Form("%s_chan%d:sys", sysname.c_str(), ic), 3) << endl;
    sysname = Form("%s:sysTriggerScanTot_%s_%s", fSetup.c_str(), string1.c_str(), string2.c_str());
    if (TMath::Abs((valMax-val0))/val0 > TMath::Abs((val0-valMin))/val0) {
      fTEX << formatTex(TMath::Abs((valMax-val0))/val0,  Form("%s_chan%d:sys", sysname.c_str(), ic), 3) << endl;
    } else {
      fTEX << formatTex(TMath::Abs((val0-valMin))/val0,  Form("%s_chan%d:sys", sysname.c_str(), ic), 3) << endl;
    }
  }

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

  if (20 == mode) {
    int i = kCyan+1;
    inputFiles.push_back("BFs00_2/scanBDT-2016BFs00_2.root");    colors.push_back(i);
    i = kCyan+2;
    inputFiles.push_back("BFs00_219/scanBDT-2016BFs00_219.root");    colors.push_back(i);
    i = kGreen+1;
    inputFiles.push_back("BFs00_8669/scanBDT-2016BFs00_8669.root");    colors.push_back(i);
    i = kBlue+1;
    inputFiles.push_back("BFs01_100/scanBDT-2016BFs01_100.root");    colors.push_back(i);
    i = kBlue+2;
    inputFiles.push_back("BFs01_7849/scanBDT-2016BFs01_7849.root");    colors.push_back(i);
    i = kRed+1;
    inputFiles.push_back("BFs02_1919/scanBDT-2016BFs02_1919.root");    colors.push_back(i);
  } else if (21 == mode) {
    int i = kCyan+1;
    inputFiles.push_back("GHs00_2/scanBDT-2016GHs00_2.root");    colors.push_back(i);
    i = kCyan+2;
    inputFiles.push_back("GHs00_219/scanBDT-2016GHs00_219.root");    colors.push_back(i);
    i = kBlue+1;
    inputFiles.push_back("GHs01_101/scanBDT-2016GHs01_101.root");    colors.push_back(i);
    i = kRed+1;
    inputFiles.push_back("GHs02_32959/scanBDT-2016GHs02_32959.root");    colors.push_back(i);
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

    TH1D *h0(0);
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
void plotResults::showScanBDT(string what) {
  c0->Clear();
  c0->SetCanvasSize(700, 500);

  if ("all" == what) {
    showScanBDT("ZAS");
    showScanBDT("CSBF");
    showScanBDT("BSMMBFS");
    return;
  }

  string bname("hBdt_bsmmMcComb");
  if (string::npos != what.find("CSBF")) {
    bname = "hBdt_bspsiphiMcComb";
  }

  vector<Color_t> colors;
  colors.push_back(kBlue+2);
  colors.push_back(kRed+2);

  string inputFile = Form("%s/scanBDT-%s.root", fDirectory.c_str(), fSuffix.c_str());
  cout << "open " << inputFile << endl;
  TFile *f = TFile::Open(inputFile.c_str());
  string hname(Form("bdtScan_%s", what.c_str()));
  string pdfname(inputFile.c_str());
  replaceAll(pdfname, Form("%s/", fDirectory.c_str()), "");
  replaceAll(pdfname, "scanBDT-", "");
  replaceAll(pdfname, ".root", "");

  TH1D *h0(0);
  string s("");

  TH1D *hbdt = (TH1D*)f->Get("bdtCuts");
  cout << "hbdt = " << hbdt << endl;
  for (int chan = 0; chan < 2; ++chan) {
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
    h0->SetLineColor(colors[chan]);
    h0->SetMarkerColor(colors[chan]);
    h0->SetMarkerStyle(2);
    h0->SetMarkerSize(0.6);
    setTitles(h0, "100 #times BDT >", what.c_str(), 0.05, 1.2, 1.6);
    if (0 == chan) {
      h0->DrawCopy("p"); // "p"
    } else {
      h0->DrawCopy("psame");
    }
    TMarker *pm = new TMarker(bdtCut, bdty, 28);
    pm->SetMarkerColor(colors[chan]);
    pm->SetMarkerSize(2.);
    pm->Draw();
    h0->GetListOfFunctions()->Add(pm);

    tl->SetTextSize(0.05);
    tl->SetTextColor(colors[chan]); tl->DrawLatexNDC(0.25, 0.72 - chan*0.055, Form("chan %d", chan));
  }

  tl->SetTextSize(0.05);
  tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.25, 0.80, Form("%s", fSetup.c_str()));
  c0->SaveAs(Form("%s/showScanBDT-%s-%s.pdf", fDirectory.c_str(), fSetup.c_str(), what.c_str()));

}



// ----------------------------------------------------------------------
void plotResults::dumpDatasets() {

  printCuts(cout);

  ofstream TEX;
  string dsname = Form("%s/%s-datasets.tex", fDirectory.c_str(), fSuffix.c_str());
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
      TEX << formatTexErrSci(it->second->fBf, it->second->fBfE, Form("%s:%s:bf", fSuffix.c_str(), it->first.c_str()), 3) << endl;
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
  tl->DrawLatexNDC(0.2, 0.8, Form("Events: %d", static_cast<int>(T->GetEntries())));

  tl->DrawLatexNDC(0.2, 0.35, Form("#tau"));
  tl->DrawLatexNDC(0.25, 0.35, Form("= (%4.3f #pm %5.3f)ps", t, tE));

  tl->DrawLatexNDC(0.2, 0.25, Form("#varepsilon"));
  tl->DrawLatexNDC(0.25, 0.25, Form("= %4.3f ", teta->GetSumOfWeights()/heta->GetSumOfWeights()));

  c1->SaveAs(Form("%s/genSummary-%d-%s.pdf", fDirectory.c_str(), fYear, dsname.c_str()));

}

// ----------------------------------------------------------------------
void plotResults::fillAndSaveHistograms(int start, int nevents) {


  //  ROOT::EnableImplicitMT(4); // does not work. Maybe not enabled in ROOT compilation?
  // /cvmfs/cms.cern.ch/slc6_amd64_gcc530/lcg/root/6.08.07/

  // -- dump histograms
  fHistFile = TFile::Open(fHistFileName.c_str(), "RECREATE");
  cout << " opened, running on " << nevents << " entries" << endl;

  TTree *t(0);

  fSaveSmallTree = true;
  string mode("");

  if (0) {
    // resetHistograms();
    rareBgHists("nada", 10000);

    resetHistograms();
    setup("bsmmMcComb");
    t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    loopOverTree(t, 1, 10000, start);
    otherNumbers(fSample);
    saveHistograms(fSample);

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
  if (2016 == fYear) {
    string bac = fHistFileName;
    replaceAll(fHistFileName, "Results", "Lifetime");
    fHistFile = TFile::Open(fHistFileName.c_str(), "RECREATE");
    vector<int> v;
    v.push_back(35);
    v.push_back(40);
    v.push_back(45);
    v.push_back(50);
    v.push_back(55);
    v.push_back(60);
    v.push_back(65);
    v.push_back(66);
    v.push_back(67);
    v.push_back(68);
    v.push_back(69);
    v.push_back(70);
    v.push_back(75);
    v.push_back(80);

    for (unsigned int i = 0; i < v.size(); ++i) {
      resetHistograms();
      setup(Form("bsmm%dMc", v[i]));
      t = getTree(fSample, fTreeDir);
      setupTree(t, fSample);
      loopOverTree(t, 1);
      saveHistograms(fSample);
    }
    fHistFileName = bac;
  }

}

// ----------------------------------------------------------------------
void plotResults::initNumbers(anaNumbers &a) {
  a.clear();
  a.fSel = (fDoUseBDT?"bdt":"cnc");
}


// ----------------------------------------------------------------------
void plotResults::rareBgHists(string smode, int nevents, int nbg) {
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
    if (nbg > -1 && (nloop > nbg)) break;
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

  for (int i = 0; i < fNchan; ++i) {
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

  double noBf  = fDS[aNorm.fNameMc]->fBf;

  double yield(0.);
  if (useW8) {
    yield = (sgBf/noBf) * pRatio * (aSig.fEffTot.val/aNorm.fEffTot.val) * aNorm.fW8SignalFit.val;
  } else {
    yield = (sgBf/noBf) * pRatio * (aSig.fEffTot.val/aNorm.fEffTot.val) * aNorm.fSignalFit.val;
  }
  aSig.fScaledYield.val = yield;

  // cout << "===scaleYield: sgBf = " << sgBf << " noBf = " << noBf << " eff(sg) = " << aSig.fEffTot.val
  //      << " eff(no) = " << aNorm.fEffTot.val
  //      << " N(no) = " << aNorm.fSignalFit.val
  //      << " yield = " << yield
  //      << endl;

  // -- scale fMcYield
  if (aSig.fMcYield[4].val > 0.) {
    double scaleFactor = yield/aSig.fMcYield[4].val;
    for (int i = 0; i < NWIN; ++i) {
      cout << "[" << i << "] " << aSig.fMcYield[i].val << " -> ";
      aSig.fMcYield[i].val = aSig.fMcYield[i].val*scaleFactor;
      cout << aSig.fMcYield[i].val << " ";
      // FIXME error scaling
    }
    cout << endl;
  } else {
    cout << "+++ fMcYield[4] = " << aSig.fMcYield[4].val << " (not larger than 0), not scaling fMcYield (that's OK in background)" << endl;
  }

  // -- anti-muon selection (but normalized to the normal B+)
  if (useW8) {
    yield = (sgBf/noBf) * pRatio * (aSig.fEffTotAm.val/aNorm.fEffTot.val) * aNorm.fW8SignalFit.val;
  } else {
    yield = (sgBf/noBf) * pRatio * (aSig.fEffTotAm.val/aNorm.fEffTot.val) * aNorm.fSignalFit.val;
  }

  // cout << "===scaleYieldAMS: sgBf = " << sgBf << " noBf = " << noBf << " eff(sg) = " << aSig.fEffTotAm.val
  //      << " eff(no) = " << aNorm.fEffTot.val
  //      << " N(no) = " << aNorm.fSignalFit.val
  //      << " yield = " << yield
  //      << endl;

  aSig.fScaledYieldAm.val = yield;
  if (aSig.fAmYield[4].val > 0.) {
    double scaleFactor = yield/aSig.fAmYield[4].val;
    cout << "+++ fAmYield[4] = " << aSig.fAmYield[4].val << " new scaled signal yield = " << yield << " -> scale factor = " << scaleFactor << endl;
    for (int i = 0; i < NWIN; ++i) {
      cout << "[" << i << "] " << aSig.fAmYield[i].val << " -> ";
      aSig.fAmYield[i].val = aSig.fAmYield[i].val*scaleFactor;
      cout << aSig.fAmYield[i].val << " ";
      // FIXME error scaling
    }
    cout << endl;
  } else {
    cout << "+++ fAmYield[4] = " << aSig.fAmYield[4].val << " (not larger than 0), not scaling fAmYield (that's OK in background)" << endl;
  }

  cout << "===scaleYield= " << aSig.fNameMc << " sgBf = " << sgBf << " noBf = " << noBf << " pRatio = " << pRatio
       << " sgEffTot = " << aSig.fEffTot.val
       << " sgEffTotAm = " << aSig.fEffTotAm.val
       << " noEffTot = " << aNorm.fEffTot.val
       << endl
       << "   noSignalFit = " << aNorm.fSignalFit.val
       << " noW8SignalFit = " << aNorm.fW8SignalFit.val
       << " ->  scaled Yield = " << aSig.fScaledYield.val
       << " ->  scaled AM Yield = " << aSig.fScaledYieldAm.val
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

  double buBf  = fDS[fNoNumbers[ichan].fNameMc]->fBf;
  double buBfE = fDS[fNoNumbers[ichan].fNameMc]->fBfE;


  number bf;
  double alpha = (1./fFsfu.val)
    * (fNoNumbers[ichan].fEffTot4R.val / fCsNumbers[ichan].fEffTot4R.val)
    * buBf;

  double alphaEstat = alpha * TMath::Sqrt((fFsfu.estat/fFsfu.val)*(fFsfu.estat/fFsfu.val)
					  + (fNoNumbers[ichan].fEffTot4R.estat/fNoNumbers[ichan].fEffTot4R.val)*(fNoNumbers[ichan].fEffTot4R.estat/fNoNumbers[ichan].fEffTot4R.val)
					  + (fCsNumbers[ichan].fEffTot4R.estat/fCsNumbers[ichan].fEffTot4R.val)*(fCsNumbers[ichan].fEffTot4R.estat/fCsNumbers[ichan].fEffTot4R.val));
  double alphaEsyst = alpha * TMath::Sqrt((fFsfu.esyst/fFsfu.val)*(fFsfu.esyst/fFsfu.val)
					  + (fNoNumbers[ichan].fEffTot4R.esyst/fNoNumbers[ichan].fEffTot4R.val)*(fNoNumbers[ichan].fEffTot4R.esyst/fNoNumbers[ichan].fEffTot4R.val)
					  + (fCsNumbers[ichan].fEffTot4R.esyst/fCsNumbers[ichan].fEffTot4R.val)*(fCsNumbers[ichan].fEffTot4R.esyst/fCsNumbers[ichan].fEffTot4R.val)
					  + (buBfE/buBf)*(buBfE/buBf));

  bf.val   = (fCsNumbers[ichan].fW8SignalFit.val / fNoNumbers[ichan].fW8SignalFit.val) * alpha;
  bf.estat = bf.val * TMath::Sqrt((fCsNumbers[ichan].fW8SignalFit.estat/fCsNumbers[ichan].fW8SignalFit.val)*(fCsNumbers[ichan].fW8SignalFit.estat/fCsNumbers[ichan].fW8SignalFit.val)
				  + (fNoNumbers[ichan].fW8SignalFit.estat/fNoNumbers[ichan].fW8SignalFit.val)*(fNoNumbers[ichan].fW8SignalFit.estat/fNoNumbers[ichan].fW8SignalFit.val)
				  + (alphaEstat/alpha)*(alphaEstat/alpha));
  bf.esyst = bf.val * TMath::Sqrt((fCsNumbers[ichan].fW8SignalFit.esyst/fCsNumbers[ichan].fW8SignalFit.val)*(fCsNumbers[ichan].fW8SignalFit.esyst/fCsNumbers[ichan].fW8SignalFit.val)
				  + (fNoNumbers[ichan].fW8SignalFit.esyst/fNoNumbers[ichan].fW8SignalFit.val)*(fNoNumbers[ichan].fW8SignalFit.esyst/fNoNumbers[ichan].fW8SignalFit.val)
				  + (alphaEsyst/alpha)*(alphaEsyst/alpha));
  bf.calcEtot();


  cout << "*******************************" << endl;
  cout << "psiphi w8 yield: " << fCsNumbers[ichan].fW8SignalFit.val
       << " +/- " << fCsNumbers[ichan].fW8SignalFit.estat
       << " +/- " << fCsNumbers[ichan].fW8SignalFit.esyst
       << " -> " << fCsNumbers[ichan].fW8SignalFit.esyst/fCsNumbers[ichan].fW8SignalFit.val
       << endl;
  cout << "psik w8 yield:   " << fNoNumbers[ichan].fW8SignalFit.val
       << " +/- " << fNoNumbers[ichan].fW8SignalFit.estat
       << " +/- " << fNoNumbers[ichan].fW8SignalFit.esyst
       << " -> " << fNoNumbers[ichan].fW8SignalFit.esyst/fNoNumbers[ichan].fW8SignalFit.val
       << endl;
  cout << "No efftot:       " << fNoNumbers[ichan].fEffTot4R.val << " +/- " << fNoNumbers[ichan].fEffTot4R.esyst
       << " -> " << fNoNumbers[ichan].fEffTot4R.esyst/fNoNumbers[ichan].fEffTot4R.val << endl;
  cout << "Cs efftot:       " << fCsNumbers[ichan].fEffTot4R.val << " +/- " << fCsNumbers[ichan].fEffTot4R.esyst
       << " -> " << fCsNumbers[ichan].fEffTot4R.esyst/fCsNumbers[ichan].fEffTot4R.val       << endl;
  cout << "buBf:            " << buBf << " +/- " << buBfE
       << " -> " << buBfE/buBf
       << endl;
  cout << "fsfu:            " << fFsfu.val << " +/- " << fFsfu.esyst
       << " -> " << fFsfu.esyst/fFsfu.val
       << endl;
  cout << "alpha:           " << alpha << " +/- " << alphaEstat << " +/- " << alphaEsyst
       << " -> "  << alphaEsyst/alpha
       << endl;
  cout << "BF:              " << bf.val << " +/- " << bf.estat << " +/- " << bf.esyst
       << " -> " << bf.esyst/bf.val
       << endl;
  cout << "*******************************" << endl;

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
    * (fNoNumbers[ichan].fEffTot4R.val / fBsmmNumbers[ichan].fEffTot4R.val)
    * buBf
    / fNoNumbers[ichan].fW8SignalFit.val;
  double alphaUE = TMath::Sqrt(
			       (fFsfu.etot/fFsfu.val)*(fFsfu.etot/fFsfu.val)
			       + (fNoNumbers[ichan].fEffTot4R.etot/fNoNumbers[ichan].fEffTot4R.val)
			       * (fNoNumbers[ichan].fEffTot4R.etot/fNoNumbers[ichan].fEffTot4R.val)
			       + (fBsmmNumbers[ichan].fEffTot4R.etot/fBsmmNumbers[ichan].fEffTot4R.val)
			       * (fBsmmNumbers[ichan].fEffTot4R.etot/fBsmmNumbers[ichan].fEffTot4R.val)
			       + (buBfE/buBf)*(buBfE/buBf)
			       + (fNoNumbers[ichan].fW8SignalFit.etot/fNoNumbers[ichan].fW8SignalFit.val)
			       * (fNoNumbers[ichan].fW8SignalFit.etot/fNoNumbers[ichan].fW8SignalFit.val)
			       );
  alphaUE = alphaUE * alphaU;

  double alphaS =
    (fCsNumbers[ichan].fEffTot4R.val / fBsmmNumbers[ichan].fEffTot4R.val)
    * bsBf
    / fCsNumbers[ichan].fW8SignalFit.val;
  double alphaSE = TMath::Sqrt(
			       + (fCsNumbers[ichan].fEffTot4R.etot/fCsNumbers[ichan].fEffTot4R.val)
			       * (fCsNumbers[ichan].fEffTot4R.etot/fCsNumbers[ichan].fEffTot4R.val)
			       + (fBsmmNumbers[ichan].fEffTot4R.etot/fBsmmNumbers[ichan].fEffTot4R.val)
			       * (fBsmmNumbers[ichan].fEffTot4R.etot/fBsmmNumbers[ichan].fEffTot4R.val)
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

  TH1D *hMassNoCuts              = (TH1D*)gDirectory->Get(Form("hMassNoCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassWithAnaCuts         = (TH1D*)gDirectory->Get(Form("hMassWithAnaCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassWithMuonCuts        = (TH1D*)gDirectory->Get(Form("hMassWithMuonCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassWithTriggerCuts     = (TH1D*)gDirectory->Get(Form("hMassWithTriggerCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hMassWithAllCuts         = (TH1D*)gDirectory->Get(Form("%s_%s_chan%d", fHistWithAllCuts.c_str(), modifier.c_str(), chan));
  TH1D *hAmMassWithMuonCuts      = (TH1D*)gDirectory->Get(Form("hAmMassWithMuonCuts_%s_chan%d", modifier.c_str(), chan));
  TH1D *hAmMassWithAllCuts       = (TH1D*)gDirectory->Get(Form("hAmMassWithAllCuts_%s_chan%d", modifier.c_str(), chan));

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
  double cam = massIntegral(hAmMassWithMuonCuts, ALL, chan);
  double d = massIntegral(hMassWithTriggerCuts, ALL, chan);
  double e = massIntegral(hMassWithAllCuts, ALL, chan);
  double g = massIntegral(hAmMassWithAllCuts, ALL, chan);

  aa.fAccCandYield.val     = getValueByLabel(hAcceptance, "accCandYield");

  aa.fCandYield.val        = a;
  aa.fAnaYield.val         = b;
  aa.fEffAna.val           = b/a;
  //DECREPIT  aa.fEffAna.setErrors(dEff(static_cast<int>(b), static_cast<int>(a)), aa.fEffAna.val * fSystematics["effana" + syst][chan]);
  aa.fEffAna.setErrors(dEff(static_cast<int>(b), static_cast<int>(a)), aa.fEffAna.val * fSystematics["effana4r"][chan]);

  aa.fMuidYield.val        = c;
  aa.fEffMuidMC.val        = c/b;
  aa.fEffMuidMC.setErrors(dEff(static_cast<int>(c), static_cast<int>(b)), aa.fEffMuidMC.val * fSystematics["effmuid"][chan]);

  aa.fTrigYield.val        = d;
  aa.fEffTrigMC.val        = d/c;
  aa.fEffTrigMC.setErrors(dEff(static_cast<int>(d), static_cast<int>(c)),  aa.fEffTrigMC.val * fSystematics["effmuid"][chan]);

  aa.fEffTot.val           = e/aa.fGenYield.val;
  aa.fEffTot.setErrors(dEff(e, TMath::Sqrt(e), aa.fGenYield.val, aa.fGenYield.estat), quadraticSum(5, fSystematics["acceptance"][chan]
											      , fSystematics["effcand" + syst][chan]
											   // , fSystematics["effana" + syst][chan]
											      , fSystematics["effana4r"][chan]
											      , fSystematics["effmuid"][chan]
											      , fSystematics["efftrig"][chan])  * aa.fEffTot.val);

  aa.fEffTrigMCAm.val        = g/cam;
  aa.fEffTrigMCAm.setErrors(dEff(static_cast<int>(g), static_cast<int>(cam)),  aa.fEffTrigMCAm.val * fSystematics["effmuid"][chan]);

  aa.fEffTotAm.val         = g/aa.fGenYield.val;
  aa.fEffTotAm.setErrors(dEff(g, TMath::Sqrt(g), aa.fGenYield.val, aa.fGenYield.estat), quadraticSum(5, fSystematics["acceptance"][chan]
											      , fSystematics["effcand" + syst][chan]
											   // , fSystematics["effana" + syst][chan]
											      , fSystematics["effana4r"][chan]
											      , fSystematics["effmuid"][chan]
											      , fSystematics["efftrig"][chan])  * aa.fEffTotAm.val);

  // -- the following is "for ratios", excluding the correlated (canceling, to some order) systematics
  aa.fEffTot4R.val         = aa.fEffTot.val;
  aa.fEffTot4R.setErrors(dEff(e, TMath::Sqrt(e), aa.fGenYield.val, aa.fGenYield.estat), quadraticSum(2
												     , fSystematics["effcand" + syst][chan]
												     , fSystematics["effana4r"][chan])  * aa.fEffTot4R.val);

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


  tot   = massIntegral(hAmMassWithAllCuts, ALL, chan);
  lo    = massIntegral(hAmMassWithAllCuts, LO, chan);
  hi    = massIntegral(hAmMassWithAllCuts, HI, chan);
  bd    = massIntegral(hAmMassWithAllCuts, BD, chan);
  bs    = massIntegral(hAmMassWithAllCuts, BS, chan);

  aa.fAmYield[0].val   = lo;
  aa.fAmYield[1].val   = bd;
  aa.fAmYield[2].val   = bs;
  aa.fAmYield[3].val   = hi;
  aa.fAmYield[4].val   = tot;



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
  int chan(a.fChan);
  int imode(0);
  if (string::npos != a.fName.find("bupsik")) {
    numbersFromHist(a, "bupsik");
  } else if (string::npos != a.fName.find("bspsiphi")) {
    numbersFromHist(a, "bspsiphi");
    imode = 1;
  }
  // -- cross check the scaled yield normalized to B+ -> J/psi K+
  if (string::npos != a.fName.find("bspsiphi")) {
    double pRatio(fFsfu.val);
    scaleYield(a, fNoNumbers[chan], pRatio);
    // -- FIXME errors
    a.fScaledYield.setErrors(0.05*a.fScaledYield.val, 0.05*a.fScaledYield.val);
  }

  // -- need some finetuning for Bs2JpsiPhi
  double xmin = 5.2;
  double sigma = 0.030;
  if (1 == chan) {
    if (fYear < 2015) {
      sigma = 0.020;
    } else {
      sigma = 0.015;
    }
  }

  // -- data: fit yields
  fSample = a.fNameDa;
  //  string  name = Form("hNorm_%s_%s_chan%d", modifier.c_str(), fSample.c_str(), chan);
  string  name = Form("hNorm_%s_%s_chan%d", fSuffixSel.c_str(), fSample.c_str(), chan);
  bool ok = fHistFile->cd(fSample.c_str());
  cout << "cd to " << fSample << ": " << ok << ", normalization fitting: " << name << endl;
  fitPsYield fpy(name);
  if (string::npos != a.fName.find("bupsik")) {
    fpy.fitBu2JpsiKp(5, fDirectory + "/norm", imode, 5.0, 5.8, 0.025);
  } else if (string::npos != a.fName.find("bspsiphi")) {
    fpy.fitBs2JpsiPhi(5, fDirectory + "/norm", imode, xmin, 5.8, sigma + chan*0.015);
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
    fpy2.fitBu2JpsiKp(5, fDirectory + "/norm", imode, 5.0, 5.8, 0.025);
  } else if (string::npos != a.fName.find("bspsiphi")) {
    fpy2.fitBs2JpsiPhi(5, fDirectory + "/norm", imode, xmin, 5.8, sigma + chan*0.015);
  }
  a.fW8SignalFit.val   = fpy2.getSignalYield();
  a.fW8SignalFit.estat = fpy2.getSignalError();
  a.fW8SignalFit.esyst = fSystematics["norm" + a.fName][chan] * fpy2.getSignalYield();
  a.fW8SignalFit.etot  = TMath::Sqrt(a.fW8SignalFit.estat*a.fW8SignalFit.estat + a.fW8SignalFit.esyst*a.fW8SignalFit.esyst);

  printNumbers(a, cout);

  cout << "chan " << a.fChan << " total yield: " << fpy.getSignalYield() << " +/- "  << fpy.getSignalError() << endl;
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- NORMALIZATION: " << mode << " chan: " << chan << endl;
  int NDIG(6);
  dumpTex(a.fEffGenSel, Form("%s:GENSEL-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fAcc, Form("%s:ACC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffCand, Form("%s:EFF-CAND-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffMuidMC, Form("%s:EFF-MU-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTrigMC, Form("%s:EFF-TRIG-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTrigMCAm, Form("%s:EFF-TRIG-MC-AM-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffAna, Form("%s:EFF-ANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffProdMC, Form("%s:EFF-PRODMC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTot, Form("%s:EFF-TOT-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTotAm, Form("%s:EFF-TOTAM-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);

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

  string aname = name; replaceAll(aname, "hMass", "hAmMass");
  TH1D *h1am = (TH1D*)gDirectory->Get(Form("%s/%s", fSample.c_str(), aname.c_str()));
  cout << "getting  histogram ->" << Form("%s/%s", fSample.c_str(), aname.c_str()) << "<-" << endl;

  TF1 *lF1(0), *lF2(0);
  TF1 *lamF1(0), *lamF2(0);

  if (0 == mode) {
    lF1   = fIF->pol0BsBlind(h1);
    lF2   = fIF->pol0(h1);
    lamF1 = fIF->pol0BsBlind(h1am);
    lamF2 = fIF->pol0(h1am);
  }
  double massInt = massIntegral(h1, HI, a.fChan);
  double massIntam = massIntegral(h1am, HI, a.fChan);

  cout << "lF1   par0 = " << lF1->GetParameter(0) << " from hist " << h1->GetName() << " with " << h1->GetEntries() << " entries "
       << " integral = " << massInt << endl;
  cout << "lamF1 par0 = " << lamF1->GetParameter(0) << " from hist " << h1am->GetName() << " with " << h1am->GetEntries() << " entries "
       << " integral = " << massIntam << endl;


  lF2->SetLineStyle(kDashed);
  lamF2->SetLineStyle(kDashed);
  double binw = h1->GetBinWidth(1);
  TFitResultPtr r;
  bool notFit(false),  notFitAm(false);

  double elo(0.), ehi(0.);
  double erel = poissonError(static_cast<int>(massInt), ehi, elo);
  // // -- choose the smaller lower error
  // erel = elo;
  cout << "combErrors: errors " << static_cast<int>(massInt) << " +" << ehi << " -" << elo << " -> +/- " << erel
       << " -> relative error: " << erel/massInt << endl;
  erel = erel/massInt;
  double erelLo = elo/massInt;
  double erelHi = ehi/massInt;

  double erelAm = poissonError(static_cast<int>(massIntam), ehi, elo);
  erelAm = erelAm/massIntam;
  double erelAmLo = elo/massIntam;
  double erelAmHi = ehi/massIntam;

  //  erelAm = elo;

  if (massInt > 3) {
    h1->Fit(lF1, "rl", "", lo, hi);
  } else {
    notFit = true;
    lF1->SetParameter(0, massInt*binw/(hi-lo));              // A = p0 * (hi-lo), N = A/binw -> p0 = A/(hi-lo) = N*binw/(hi-lo)
    //    lF1->SetParError(0, TMath::Sqrt(massInt)*binw/(hi-lo));  // A = p0 * (hi-lo), N = A/binw -> p0 = A/(hi-lo) = N*binw/(hi-lo)
    lF1->SetParError(0, erel*lF1->GetParameter(0));
  }

  if (massIntam > 3) {
    h1am->Fit(lamF1, "rl", "", lo, hi);
  } else {
    notFitAm = true;
    lamF1->SetParameter(0, massIntam*binw/(hi-lo));              // A = p0 * (hi-lo), N = A/binw -> p0 = A/(hi-lo) = N*binw/(hi-lo)
    lamF1->SetParError(0, TMath::Sqrt(massIntam)*binw/(hi-lo));  // A = p0 * (hi-lo), N = A/binw -> p0 = A/(hi-lo) = N*binw/(hi-lo)
    lamF1->SetParError(0, erelAm*lamF1->GetParameter(0));
  }

  setTitles(h1, "#it{m}_{#it{#mu #mu}} [GeV]", Form("Candidates / %4.3f GeV", h1->GetBinWidth(1)));
  setTitles(h1am, "#it{m}_{#it{#mu #mu}}^{AMS} [GeV]", Form("Candidates / %4.3f GeV", h1am->GetBinWidth(1)));

  h1->DrawCopy();
  lF2->SetParameters(lF1->GetParameters());
  lF2->SetParErrors(lF1->GetParErrors());
  lF2->SetLineColor(kBlue);
  //  if (notFit) lF1->Draw("same");
  lF1->Draw("same");
  lF2->Draw("same");

  double tsize = tl->GetTextSize();
  tl->SetTextSize(0.03);
  tl->DrawLatexNDC(0.6, 0.80, Form("Blind %s", fSuffixSel.c_str()));
  tl->DrawLatexNDC(0.6, 0.76, Form("chan: %d", a.fChan));
  tl->DrawLatexNDC(0.6, 0.72, Form("mode: %d", mode));
  tl->SetTextSize(tsize);
  c0->Modified();
  c0->Update();
  savePad(Form("%s-combBg-mode%d-chan%d.pdf", fSuffixSel.c_str(), mode, a.fChan));

  h1am->DrawCopy();
  lamF2->SetParameters(lamF1->GetParameters());
  lamF2->SetParErrors(lamF1->GetParErrors());
  lamF2->SetLineColor(kBlue);
  //  if (notFit) lamF1->Draw("same");
  lamF1->Draw("same");
  lamF2->Draw("same");
  tl->DrawLatexNDC(0.6, 0.80, Form("Blind %s", fSuffixSel.c_str()));
  tl->DrawLatexNDC(0.6, 0.76, Form("chan: %d", a.fChan));
  tl->DrawLatexNDC(0.6, 0.72, Form("mode: %d", mode));
  tl->SetTextSize(tsize);
  c0->Modified();
  c0->Update();
  savePad(Form("%s-combBg-AM-mode%d-chan%d.pdf", fSuffixSel.c_str(), mode, a.fChan));

  // -- low sideband
  a.fObsYield[0].val = massIntegral(h1, LO, a.fChan);
  a.fObsYield[0].estat = TMath::Sqrt(a.fObsYield[0].val);
  a.fObsYield[0].calcEtot();
  a.fFitYield[0].val = lF2->Integral(fBgLo, fCuts[a.fChan]->mBdLo)/binw;
  a.fFitYield[0].estat = erel * a.fFitYield[0].val;
  a.fFitYield[0].estatLo = erelLo * a.fFitYield[0].val;
  a.fFitYield[0].estatHi = erelHi * a.fFitYield[0].val;
  a.fFitYield[0].calcEtot();

  a.fObsYieldAm[0].val = massIntegral(h1am, LO, a.fChan);
  //  cout << "lo sideband a.fObsYieldAm[1].val = " << a.fObsYieldAm[0].val << endl;
  a.fObsYieldAm[0].estat = TMath::Sqrt(a.fObsYieldAm[0].val);
  a.fObsYieldAm[0].calcEtot();
  a.fFitYieldAm[0].val = lamF2->Integral(fBgLo, fCuts[a.fChan]->mBdLo)/binw;
  a.fFitYieldAm[0].estat = erelAm * a.fFitYieldAm[0].val;
  a.fFitYieldAm[0].calcEtot();


  // -- B0 window
  a.fObsYield[1].val = massIntegral(h1, BD, a.fChan);
  a.fObsYield[1].estat = TMath::Sqrt(a.fObsYield[1].val);
  a.fObsYield[1].calcEtot();
  a.fFitYield[1].val = lF2->Integral(fCuts[a.fChan]->mBdLo, fCuts[a.fChan]->mBdHi)/binw;
  a.fFitYield[1].estat = erel * a.fFitYield[1].val;
  a.fFitYield[1].estatLo = erelLo * a.fFitYield[1].val;
  a.fFitYield[1].estatHi = erelHi * a.fFitYield[1].val;
  a.fFitYield[1].calcEtot();

  a.fObsYieldAm[1].val = massIntegral(h1am, BD, a.fChan);
  //  cout << "B0 window a.fObsYieldAm[1].val = " << a.fObsYieldAm[1].val << endl;
  a.fObsYieldAm[1].estat = TMath::Sqrt(a.fObsYieldAm[1].val);
  a.fObsYieldAm[1].calcEtot();
  a.fFitYieldAm[1].val = lamF2->Integral(fCuts[a.fChan]->mBdLo, fCuts[a.fChan]->mBdHi)/binw;
  a.fFitYieldAm[1].estat = erelAm * a.fFitYieldAm[1].val;
  a.fFitYieldAm[1].calcEtot();


  // -- Bs window
  a.fObsYield[2].val = massIntegral(h1, BS, a.fChan);
  a.fObsYield[2].estat = TMath::Sqrt(a.fObsYield[2].val);
  a.fObsYield[2].calcEtot();
  a.fFitYield[2].val = lF2->Integral(fCuts[a.fChan]->mBsLo, fCuts[a.fChan]->mBsHi)/binw;
  a.fFitYield[2].estat = erel * a.fFitYield[2].val;
  a.fFitYield[2].estatLo = erelLo * a.fFitYield[2].val;
  a.fFitYield[2].estatHi = erelHi * a.fFitYield[2].val;
  a.fFitYield[2].calcEtot();

  a.fObsYieldAm[2].val = massIntegral(h1am, BS, a.fChan);
  //  cout << "Bs window a.fObsYieldAm[2].val = " << a.fObsYieldAm[2].val << endl;
  a.fObsYieldAm[2].estat = TMath::Sqrt(a.fObsYieldAm[2].val);
  a.fObsYieldAm[2].calcEtot();
  a.fFitYieldAm[2].val = lamF2->Integral(fCuts[a.fChan]->mBsLo, fCuts[a.fChan]->mBsHi)/binw;
  a.fFitYieldAm[2].estat = erelAm * a.fFitYieldAm[2].val;
  a.fFitYieldAm[2].calcEtot();

  // -- high sideband
  a.fObsYield[3].val = massIntegral(h1, HI, a.fChan);
  a.fObsYield[3].estat = TMath::Sqrt(a.fObsYield[3].val);
  a.fObsYield[3].calcEtot();
  a.fFitYield[3].val = lF2->Integral(fCuts[a.fChan]->mBsHi, fBgHi)/binw;
  a.fFitYield[3].estat = erel * a.fFitYield[3].val;
  a.fFitYield[3].estatLo = erelLo * a.fFitYield[3].val;
  a.fFitYield[3].estatHi = erelHi * a.fFitYield[3].val;
  a.fFitYield[3].calcEtot();

  a.fObsYieldAm[3].val = massIntegral(h1am, HI, a.fChan);
  //  cout << "hi sideband a.fObsYieldAm[1].val = " << a.fObsYieldAm[3].val << endl;
  a.fObsYieldAm[3].estat = TMath::Sqrt(a.fObsYieldAm[3].val);
  a.fObsYieldAm[3].calcEtot();
  a.fFitYieldAm[3].val = lamF2->Integral(fCuts[a.fChan]->mBsHi, fBgHi)/binw;
  a.fFitYieldAm[3].estat = erelAm * a.fFitYieldAm[3].val;
  a.fFitYieldAm[3].calcEtot();


  // -- everything
  a.fObsYield[4].val = massIntegral(h1, ALL, a.fChan);
  a.fObsYield[4].estat = TMath::Sqrt(a.fObsYield[4].val);
  a.fObsYield[4].calcEtot();
  a.fFitYield[4].val = lF2->Integral(fBgLo, fBgHi)/binw;
  a.fFitYield[4].estat = erel * a.fFitYield[4].val;
  a.fFitYield[4].estatLo = erelLo * a.fFitYield[4].val;
  a.fFitYield[4].estatHi = erelHi * a.fFitYield[4].val;
  a.fFitYield[4].calcEtot();

  a.fObsYieldAm[4].val = massIntegral(h1am, ALL, a.fChan);
  a.fObsYieldAm[4].estat = TMath::Sqrt(a.fObsYieldAm[4].val);
  a.fObsYieldAm[4].calcEtot();
  a.fFitYieldAm[4].val = lamF2->Integral(fBgLo, fBgHi)/binw;
  a.fFitYieldAm[4].estat = erelAm * a.fFitYieldAm[4].val;
  a.fFitYieldAm[4].calcEtot();

  // -- dump numbers
  int NDIG(6);
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- DATA/COMBINATORIAL " << mode << " chan " << a.fChan << endl;
  for (unsigned i = 0; i < NWIN; ++i) {
    dumpTex(a.fFitYield[i], Form("%s:N-FIT-MBIN%d-CB-chan%d", fSuffixSel.c_str(), i, a.fChan), NDIG);
    dumpTex(a.fObsYield[i], Form("%s:N-OBS-MBIN%d-DATA-chan%d", fSuffixSel.c_str(), i, a.fChan), 0);
    dumpTex(a.fFitYieldAm[i], Form("%s:N-FIT-AM-MBIN%d-DATA-chan%d", fSuffixSel.c_str(), i, a.fChan), NDIG);
    dumpTex(a.fObsYieldAm[i], Form("%s:N-OBS-AM-MBIN%d-DATA-chan%d", fSuffixSel.c_str(), i, a.fChan), 0);
  }

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
  double relBfE = fDS[a.fNameMc]->fBfE/fDS[a.fNameMc]->fBf;
  if (string::npos != a.fName.find("bdmm")) {
    pRatio = 1.;
  }
  scaleYield(a, fNoNumbers[chan], pRatio);

  double    fsfuSys(0.);
  if (pRatio < 0.9) {
    fsfuSys = fSystematics["fsfu"][chan];
  }
  double esystRel = quadraticSum(8,
				 fsfuSys,
				 fSystematics["acceptance"][chan],
				 fSystematics["effcandbsmm"][chan],
				 //fSystematics["effanabsmm"][chan],
				 fSystematics["effana4r"][chan],
				 fSystematics["efftrig"][chan],
				 2.*fSystematics["effmuid"][chan],
				 fSystematics["normbupsik"][chan],
				 relBfE
				 );
  cout << "DBXDBXSG:  esystRel = " << esystRel
       << " from: "
       << " fsfsu: " << fsfuSys << " "
       << " acc: "   << fSystematics["acceptance"][chan] << " "
       << " cand: "  << fSystematics["effcandbsmm"][chan] << " "
    //<< " ana: "   << fSystematics["effanabsmm"][chan] << " "
       << " ana: "   << fSystematics["effana4r"][chan] << " "
       << " trig: "  << fSystematics["efftrig"][chan] << " "
       << " muid: "  << 2.*fSystematics["effmuid"][chan] << " "
       << " norm: "  << fSystematics["normbupsik"][chan] << " "
       << " bf: "  << relBfE << " "
       << endl;

  double estatRel = TMath::Sqrt(a.fProdGenYield.val)/a.fProdGenYield.val;
  a.fScaledYield.setErrors(estatRel*a.fScaledYield.val, a.fScaledYield.val*esystRel);
  a.fScaledYieldAm.setErrors(estatRel*a.fScaledYieldAm.val, a.fScaledYieldAm.val*esystRel);

  for (int i = 0; i < NWIN; ++i)  {
    a.fMcYield[i].setErrors(estatRel*a.fMcYield[i].val, a.fMcYield[i].val*esystRel);
    a.fAmYield[i].setErrors(estatRel*a.fAmYield[i].val, a.fAmYield[i].val*esystRel);
  }


  // -- data: fitted/interpolated yields
  fSample = a.fNameDa;

  int NDIG(6);
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- SIGNAL " << mode << " chan " << chan << endl;
  dumpTex(a.fEffGenSel, Form("%s:GENSEL-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fAcc, Form("%s:ACC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffCand, Form("%s:EFF-CAND-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffMuidMC, Form("%s:EFF-MU-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTrigMC, Form("%s:EFF-TRIG-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTrigMCAm, Form("%s:EFF-TRIG-MC-AM-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffAna, Form("%s:EFF-ANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffProdMC, Form("%s:EFF-PRODMC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTot, Form("%s:EFF-TOT-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
  dumpTex(a.fEffTotAm, Form("%s:EFF-TOTAM-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);

  dumpTex(a.fGenFileYield, Form("%s:N-GENFILEYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fGenYield, Form("%s:N-GENYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fTotGenYield, Form("%s:N-TOTYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fProdGenYield, Form("%s:N-PRODYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
  dumpTex(a.fScaledYield, Form("%s:N-SCALEDYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);

  for (unsigned i = 0; i < a.fMcYield.size(); ++i) {
    dumpTex(a.fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-%s-chan%d", fSuffixSel.c_str(), i, mode, chan), NDIG);
    dumpTex(a.fAmYield[i], Form("%s:N-AM-SCALEDYIELD-MBIN%d-%s-chan%d", fSuffixSel.c_str(), i, mode, chan), NDIG);
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
  string aname = fHistWithAllCuts; replaceAll(aname, "hMass", "hAmW8Mass");
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

  // -- determine trigger efficiency for normal muon id and AMS:
  double eff(0.), effE(0.);
  double triggerEffAms(0.), triggerEffAmsE(0.), triggerEff(0.), triggerEffE(0.);
  double triggerEffAmsAve(0.), triggerEffAmsAveE(0.), triggerEffAve(0.), triggerEffAveE(0.);
  double triggerEffSignal    = fBsmmNumbers[chan].fEffTrigMC.val;
  double triggerEffAmsSignal = fBsmmNumbers[chan].fEffTrigMCAm.val;

  triggerEfficiency(chan, eff, effE, "ams", "bg");
  if (eff > 0) {
    triggerEffAms = eff;
    triggerEffAmsE = effE;
    triggerEffAmsAve  = 0.5*(triggerEffAms + triggerEffAmsSignal);
    triggerEffAmsAveE = TMath::Abs(triggerEffAmsAve - triggerEffAmsSignal);
  }

  triggerEfficiency(chan, eff, effE, "muid", "bg");
  if (eff > 0) {
    triggerEff = eff;
    triggerEffE = effE;
    triggerEffAve  = 0.5*(triggerEff + triggerEffSignal);
    triggerEffAveE = TMath::Abs(triggerEffAve - triggerEffSignal);
  } else {
    triggerEff = triggerEffAms;
    triggerEffE = triggerEffAmsE;
    triggerEffAve  = triggerEffAmsAve;
    triggerEffAveE = triggerEffAmsAveE;
  }

  // -- temporary fix (background trigger efficiency is NOT the HLT trigger efficiency, simply L1T)
  // -- ASSUME background efficiency as 50% of signal trigger efficiency,
  // -- ASSUME background efficiency error  as 100% of background trigger efficiency

  triggerEffAmsAve  = 0.5*triggerEffAmsSignal;
  triggerEffAmsAveE = 0.5*triggerEffAmsSignal;

  triggerEffAve     = 0.5*triggerEffSignal;
  triggerEffAveE    = 0.5*triggerEffSignal;

  cout << "XXX Trigger efficiencies: "
       << " sg = "  << triggerEffSignal
       << " sgAMS = "  << triggerEffAmsSignal
       << " bg = "  << triggerEffAve
       << " bgAMS = "  << triggerEffAmsAve
       << endl;

  // -- loop over all (effective) two-body backgrounds
  string u8name(""), w8name(""), amname(""), accname("bs");
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    if (skipThisBg(it->first)) continue;
    cout << "calculateRareBgNumbers for " << it->first << " in chan = " << chan << endl;
    anaNumbers *a = fRareNumbers[it->first][chan];
    string accname = rareAccName(it->first);
    cout << "rec: " << Form("%s/hGenAndAccNumbers_%s_%s_chan%d", it->first.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan) << endl;
    TH1D *h1 = (TH1D*)fHistFile->Get(Form("%s/hGenAndAccNumbers_%s_%s_chan%d", it->first.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan));
    double genYield  = getValueByLabel(h1,  "genYield");

    u8name = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), uname.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan);
    w8name = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), wname.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan);
    amname = Form("%s/%s_%s_%s_chan%d", it->first.c_str(), aname.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan);
    TH1D *hu = (TH1D*)fHistFile->Get(u8name.c_str());
    TH1D *hw = (TH1D*)fHistFile->Get(w8name.c_str());
    TH1D *ha = (TH1D*)fHistFile->Get(amname.c_str());

    double utot     = massIntegral(hu, ALL, chan);
    double ulo      = massIntegral(hu, LO, chan);
    double ubd      = massIntegral(hu, BD, chan);
    double ubs      = massIntegral(hu, BS, chan);
    double uhi      = massIntegral(hu, HI, chan);
    double eAccAna  = utot/genYield;
    // -- now include trigger efficiency and misid weight:
    double eTot     = eAccAna * (massIntegral(hw, ALL, chan) / utot) * fBsmmNumbers[chan].fEffTrigMC.val;
    double eTotAm   = eAccAna * (massIntegral(ha, ALL, chan) / utot) * fBsmmNumbers[chan].fEffTrigMCAm.val;
    // -- switch to background estimates:
    eTot     = eAccAna * (massIntegral(hw, ALL, chan) / utot) * triggerEffAve;
    eTotAm   = eAccAna * (massIntegral(ha, ALL, chan) / utot) * triggerEffAmsAve;
    if (1) cout << "acc = " << eAccAna << " integral = " << massIntegral(ha, ALL, chan) << " utot = " << utot << " eff = " << triggerEffAmsAve << endl;

    fSample = a->fNameMc;
    a->fEffAccAna.val = eAccAna;
    a->fEffAccAna.setErrors(0.01*eAccAna, 0.05*eAccAna);
    a->fEffTot.val = eTot;
    a->fEffTot.setErrors(0.01*eTot, 0.5*eTot);
    a->fEffTotAm.val = eTotAm;
    a->fEffTotAm.setErrors(0.01*eTotAm, 0.5*eTotAm);
    a->fEffTot4R.val = eTot;
    a->fEffTot4R.setErrors(0.01*eTot, 0.02*eTot);

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

    // -- scale total expected yield to measured B+ with epsTot
    //    NB: At this point, fMcYield[] and fAmYield[] have not been filled yet!
    scaleYield(*a, fNoNumbers[chan], pRatio);
    cout << "^^^^^^^^^^^^^^^ a->fScaledYield.val/massIntegral(hw, ALL, chan) = "
	 << a->fScaledYield.val << " / " << massIntegral(hw, ALL, chan) << " = "
	 << a->fScaledYield.val/massIntegral(hw, ALL, chan)
	 << endl;
    // -- scale the existing histogram to the expected yield (this accounts for the equivalent MC lumi)
    hw->Scale(a->fScaledYield.val/massIntegral(hw, ALL, chan));
    ha->Scale(a->fScaledYieldAm.val/massIntegral(ha, ALL, chan));

    int nprotons = ((string::npos != it->first.find("lb")) ? 1 : 0);
    int nmuons = ((string::npos != it->first.find("mu")) ? 1 : 0);
    if (string::npos != it->first.find("mumu")) nmuons = 2;
    double muonSys(0.), fsfuSys(0.);
    if (0 == nmuons) {
      if (0 == nprotons) {
	muonSys = 2.*fSystematics["fakemuid"][chan];
      } else {
	muonSys = quadraticSum(2, fSystematics["fakemuid"][chan], fSystematics["fakemuidProtons"][chan]);
      }
    } else if (1 == nmuons) {
      if (0 == nprotons) {
	muonSys = quadraticSum(2, fSystematics["fakemuid"][chan], fSystematics["effmuid"][chan]);
      } else {
	muonSys = quadraticSum(2, fSystematics["fakemuidProtons"][chan], fSystematics["effmuid"][chan]);
      }
    } else {
      muonSys = 2.*fSystematics["effmuid"][chan];
    }
    if (pRatio < 0.9) {
      fsfuSys = fSystematics["fsfu"][chan];
    }
    double esystRel = quadraticSum(8,
				   fsfuSys,
				   fSystematics["acceptance"][chan],
				   fSystematics["effcandbsmm"][chan],
				   //  fSystematics["effanabsmm"][chan],
				   fSystematics["effana4r"][chan],
				   //fSystematics["efftrig"][chan],
				   triggerEffAveE/triggerEffAve, //FIXME: from trigger uncertainty above
				   muonSys,
				   fSystematics["normbupsik"][chan],
				   it->second->fBfE/it->second->fBf
				   );

    // -- MC yield with default muon ID
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

    // -- MC yield with anti-muon ID
    a->fAmYield[0].val = massIntegral(ha, LO, chan);
    integral = massIntegral(hu, LO, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fAmYield[0].setErrors(a->fAmYield[0].val*estatRel, a->fAmYield[0].val*esystRel);

    a->fAmYield[1].val = massIntegral(ha, BD, chan);
    integral = massIntegral(hu, BD, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fAmYield[1].setErrors(a->fAmYield[1].val*estatRel, a->fAmYield[1].val*esystRel);

    a->fAmYield[2].val = massIntegral(ha, BS, chan);
    integral = massIntegral(hu, BS, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fAmYield[2].setErrors(a->fAmYield[2].val*estatRel, a->fAmYield[2].val*esystRel);

    a->fAmYield[3].val = massIntegral(ha, HI, chan);
    integral = massIntegral(hu, HI, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fAmYield[3].setErrors(a->fAmYield[3].val*estatRel, a->fAmYield[3].val*esystRel);

    a->fAmYield[4].val = massIntegral(ha, ALL, chan);
    integral = massIntegral(hu, ALL, chan);
    estatRel = (integral > 0.?TMath::Sqrt(integral)/integral : 0.);
    a->fAmYield[4].setErrors(a->fAmYield[4].val*estatRel, a->fAmYield[4].val*esystRel);

    setFilledHist(hw, it->second->fLcolor, it->second->fFcolor, it->second->fFillStyle, 1);
    TH1D *hwrb = (TH1D*)hw->Clone(Form("hwrb_%s", hw->GetName()));
    int rebin(5);
    hwrb->Rebin(rebin);
    hwrb->SetNdivisions(510, "XYZ");
    setTitles(hwrb, "m_{#it{#mu #mu}} [GeV]", Form("Candidates / %4.3f GeV", rebin*hwrb->GetBinWidth(1)));
    if (0 == nmuons) {
      for (unsigned im = 0; im < a->fMcYield.size(); ++im) {
	fHhNumbers[chan].fMcYield[im].val += a->fMcYield[im].val;
	fHhNumbers[chan].fMcYield[im].add2Errors(a->fMcYield[im]);

	fHhNumbers[chan].fAmYield[im].val += a->fAmYield[im].val;
	fHhNumbers[chan].fAmYield[im].add2Errors(a->fAmYield[im]);
	if (1 == im) cout << "hh: "
			  << " summed Mc = " << fHhNumbers[chan].fMcYield[im].val
			  << " with " << a->fMcYield[im].val
			  << " summed Am = " << fHhNumbers[chan].fAmYield[im].val
			  << " with " << a->fAmYield[im].val
			  << " from " << a->fName
			  << endl;
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

	fSlNumbers[chan].fAmYield[im].val += a->fAmYield[im].val;
	fSlNumbers[chan].fAmYield[im].add2Errors(a->fAmYield[im]);
      }
      cout << "sl: " << fSlNumbers[chan].fMcYield[3].esyst << endl;
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
    int NDIG(6);
    fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
    fTEX << "% -- RARE " << mode << " chan " << chan << endl;
    dumpTex(a->fEffGenSel, Form("%s:GENSEL-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fAcc, Form("%s:ACC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffCand, Form("%s:EFF-CAND-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffMuidMC, Form("%s:EFF-MU-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffTrigMC, Form("%s:EFF-TRIG-MC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffTrigMCAm, Form("%s:EFF-TRIG-MC-AM-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffAna, Form("%s:EFF-ANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffProdMC, Form("%s:EFF-PRODMC-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffAccAna, Form("%s:EFF-ACCANA-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffTot, Form("%s:EFF-TOT-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);
    dumpTex(a->fEffTotAm, Form("%s:EFF-TOTAM-%s-chan%d", fSuffixSel.c_str(), mode, chan), NDIG);

    dumpTex(a->fGenFileYield, Form("%s:N-GENFILEYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
    dumpTex(a->fGenYield, Form("%s:N-GENYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
    dumpTex(a->fTotGenYield, Form("%s:N-TOTYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);
    dumpTex(a->fProdGenYield, Form("%s:N-PRODYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);

    dumpTex(a->fScaledYield, Form("%s:N-SCALEDYIELD-%s-chan%d", fSuffixSel.c_str(), mode, chan), 3);

    for (unsigned i = 0; i < a->fMcYield.size(); ++i) {
      dumpTex(a->fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-%s-chan%d", fSuffixSel.c_str(), i, mode, chan), NDIG);
      dumpTex(a->fAmYield[i], Form("%s:N-AM-MBIN%d-%s-chan%d", fSuffixSel.c_str(), i, mode, chan), NDIG);
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
    fNpNumbers[chan].fMcYield[im].add2Errors(fSlNumbers[chan].fMcYield[im]);
    fNpNumbers[chan].fMcYield[im].val +=  fCombNumbers[chan].fFitYield[im].val;
    fNpNumbers[chan].fMcYield[im].add2Errors(fCombNumbers[chan].fFitYield[im]);

    fBgNumbers[chan].fMcYield[im].val +=  fHhNumbers[chan].fMcYield[im].val;
    fBgNumbers[chan].fMcYield[im].add2Errors(fHhNumbers[chan].fMcYield[im]);
    fBgNumbers[chan].fMcYield[im].val +=  fNpNumbers[chan].fMcYield[im].val;
    fBgNumbers[chan].fMcYield[im].add2Errors(fNpNumbers[chan].fMcYield[im]);

    // -- the sl numbers are corrected later on after "fitting" for the scale factor
    fSgAndBgNumbers[chan].fMcYield[im].val +=  fNpNumbers[chan].fMcYield[im].val;
    fSgAndBgNumbers[chan].fMcYield[im].add2Errors(fNpNumbers[chan].fMcYield[im]);
    fSgAndBgNumbers[chan].fMcYield[im].val +=  fHhNumbers[chan].fMcYield[im].val;
    fSgAndBgNumbers[chan].fMcYield[im].add2Errors(fHhNumbers[chan].fMcYield[im]);
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
  h1e->GetYaxis()->SetTitle(Form("Candidates / %4.3f GeV", h1e->GetBinWidth(1)));
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
  h1e->SetMaximum(9.);
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
    fSlNumbers[chan].fFitYield[iw].setErrors(fSlNumbers[chan].fScaleFactor*fSlNumbers[chan].fMcYield[iw].estat,
					     fSlNumbers[chan].fScaleFactor*fSlNumbers[chan].fMcYield[iw].esyst);
    fNpNumbers[chan].fFitYield[iw].val      = fNpNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
    fNpNumbers[chan].fFitYield[iw].setErrors(fNpNumbers[chan].fMcYield[iw].estat+fSlNumbers[chan].fMcYield[iw].estat,
					     fNpNumbers[chan].fMcYield[iw].esyst+fSlNumbers[chan].fMcYield[iw].esyst);
    fBgNumbers[chan].fFitYield[iw].val      = fBgNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
    fBgNumbers[chan].fFitYield[iw].setErrors(fBgNumbers[chan].fMcYield[iw].estat+fSlNumbers[chan].fMcYield[iw].estat,
					     fBgNumbers[chan].fMcYield[iw].esyst+fSlNumbers[chan].fMcYield[iw].esyst);
    fSgAndBgNumbers[chan].fFitYield[iw].val = fSgAndBgNumbers[chan].fMcYield[iw].val + (fSlNumbers[chan].fMcYield[iw].val * (fSlNumbers[chan].fScaleFactor - 1.));
    fSgAndBgNumbers[chan].fFitYield[iw].setErrors(fSgAndBgNumbers[chan].fMcYield[iw].estat+fSlNumbers[chan].fMcYield[iw].estat,
						  fSgAndBgNumbers[chan].fMcYield[iw].esyst+fSlNumbers[chan].fMcYield[iw].esyst);
  }
  // -- add signal
  for (int iw = 0; iw < NWIN; ++iw) {
    fSgAndBgNumbers[chan].fFitYield[iw].val = fSgAndBgNumbers[chan].fFitYield[iw].val + fBsmmNumbers[chan].fMcYield[iw].val + fBdmmNumbers[chan].fMcYield[iw].val;
    fSgAndBgNumbers[chan].fFitYield[iw].add2Errors(fBsmmNumbers[chan].fMcYield[iw]);
    fSgAndBgNumbers[chan].fFitYield[iw].add2Errors(fBdmmNumbers[chan].fMcYield[iw]);
  }

  // -- dump combined/summed numbers: rare sl decays and rare hadronic (peaking) decays
  int NDIG(6);
  fTEX << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  fTEX << "% -- SUMMARY BACKGROUND " << mode << " chan " << chan << endl;
  fTEX << formatTex(100*scale, Form("%s:SCALEPERCENT-SL-chan%d", fSuffixSel.c_str(), chan), 1) << endl;
  fTEX << formatTex(scale + 1.0, Form("%s:SCALEPLUSONE-SL-chan%d", fSuffixSel.c_str(), chan), 3) << endl;
  fTEX << formatTex(scale, Form("%s:SCALE-SL-chan%d", fSuffixSel.c_str(), chan), 3) << endl;

  for (unsigned i = 0; i < fHhNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fHhNumbers[chan].fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-HH-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
    dumpTex(fHhNumbers[chan].fAmYield[i], Form("%s:N-AM-SCALEDYIELD-MBIN%d-HH-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
  }

  fTEX << formatTex(fSlNumbers[chan].fScaleFactor, Form("%s:SCALEFACTOR-SL-chan%d:val", fSuffixSel.c_str(), chan), NDIG) << endl;
  for (unsigned i = 0; i < fSlNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fSlNumbers[chan].fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-SL-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
    dumpTex(fSlNumbers[chan].fFitYield[i], Form("%s:N-FITSCLYIELD-MBIN%d-SL-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
  }
  for (unsigned i = 0; i < fNpNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fNpNumbers[chan].fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-NP-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
    dumpTex(fNpNumbers[chan].fFitYield[i], Form("%s:N-FITSCLYIELD-MBIN%d-NP-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
  }
  for (unsigned i = 0; i < fBgNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fBgNumbers[chan].fMcYield[i], Form("%s:N-SCALEDYIELD-MBIN%d-BG-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
    dumpTex(fBgNumbers[chan].fFitYield[i], Form("%s:N-FITSCLYIELD-MBIN%d-BG-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
  }

  for (unsigned i = 0; i < fBgNumbers[chan].fMcYield.size(); ++i) {
    dumpTex(fSgAndBgNumbers[chan].fFitYield[i], Form("%s:N-FITSCLYIELD-MBIN%d-SGANDBG-chan%d", fSuffixSel.c_str(), i, chan), NDIG);
  }



  delete h1e;
  delete h1c;

}


// ----------------------------------------------------------------------
void plotResults::triggerEfficiency(int chan, double &eff, double &effE, string mode, string sample) {
  eff = effE = -1.;
  string mname(""), pname("MassWithAllCuts"), aname("MassWithMuonCuts");
  if (string::npos != mode.find("ams")) {
    mname = "hRAm";
  }
  if (string::npos != mode.find("muid")) {
    mname = "hRM";
  }

  vector<double> vEff, vEffE;
  double passTot(0.), allTot(0.);
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    TH1D *hp(0), *ha(0);

    // double pRatio(0.);
    // if (string::npos != it->first.find("bs")) pRatio = fFsfu.val;
    // if (string::npos != it->first.find("lb")) pRatio = fFsfu.val; // FIXME!
    // if (string::npos != it->first.find("bd")) pRatio = 1.;
    // if (string::npos != it->first.find("bu")) pRatio = 1.;

    // -- muons and dimuons allowed only with "sg"
    int nmuons = ((string::npos != it->first.find("mu")) ? 1 : 0);
    if (string::npos != it->first.find("mumu")) nmuons = 2;

    if (string::npos != sample.find("bg")) {
      if (nmuons > 0) continue;
      if (skipThisBg(it->first)) continue;
    }

    if (string::npos != sample.find("sg")) {
      if (string::npos == it->first.find("mmMc")) continue;
    }

    string hpname = Form("%s/%s%s_%s_%s_chan%d", it->first.c_str(), mname.c_str(), pname.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan);
    string haname = Form("%s/%s%s_%s_%s_chan%d", it->first.c_str(), mname.c_str(), aname.c_str(), fSuffixSel.c_str(), it->first.c_str(), chan);
    cout << " triggerEfficiency(" << mode << "," << sample << "): " << haname << ": ";
    hp = (TH1D*)fHistFile->Get(hpname.c_str());
    ha = (TH1D*)fHistFile->Get(haname.c_str());

    double pass(0.), all(0.), lEff(0.), lEffE(0.);
    if (hp) pass = hp->GetEntries();
    if (ha) all  = ha->GetEntries();
    passTot += pass;
    allTot  += all;
    if (all < 1) {
      cout << "skipping " << it->first << " because of no statistics: " << pass << "/" << all << endl;
      continue;
    }
    lEff = pass/all;
    lEffE  = dEff(static_cast<int>(pass), static_cast<int>(all));
    cout << "Efficiency for " << it->first << ": " << pass << "/" << all << " = " << lEff << " +/- " << lEffE << endl;
    vEff.push_back(lEff);
    vEffE.push_back(lEffE);
  }

  double val0(0.), err0(0.), chi0(0.);
  if (passTot > 0.) {
    average(val0, err0, vEff, vEffE, chi0);
  } else {
    val0 = (passTot+1)/(allTot+2);
    err0 = TMath::Sqrt((passTot+1)*(allTot-passTot+1)) / TMath::Sqrt((allTot+3)*(allTot+2)*(allTot+2));
    cout << "No statistics in all, using d'Agostini p86: " << val0 << " +/- " << err0 << endl;
  }
  eff  = val0;
  effE = err0;
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

    // -- ams anti-muon selection
    if (fGoodAmsID) {
      fhAmMassWithMuonCuts[fHistStrings[0]][fChan]->Fill(mass);
      if (fGoodHLT) {
	fhAmMassWithAllCuts[fHistStrings[0]][fChan]->Fill(mass);
	fhAmW8MassWithAllCuts[fHistStrings[0]][fChan]->Fill(mass, fW8AmMisId);
      }
    }
    // -- now the "real" muon ID (to extract the "real" trigger efficiency also for backgrounds)
    if (fGoodRealMuonsID && fGoodDcand) {
      fhRMMassWithMuonCuts[fHistStrings[0]][fChan]->Fill(mass);
      if (fb.hlt1 && fb.tos && fb.l1t) {
	fhRMMassWithAllCuts[fHistStrings[0]][fChan]->Fill(mass);
      }
    }
    if (fGoodRealAmsID && fGoodDcand) {
      fhRAmMassWithMuonCuts[fHistStrings[0]][fChan]->Fill(mass);
      if (fb.hlt1 && fb.tos && fb.l1t) {
	fhRAmMassWithAllCuts[fHistStrings[0]][fChan]->Fill(mass);
      }
    }


  }

  // -----------------
  // -- BDT histograms
  // -----------------
  bool goodBDT(false);

  if ((fMode == BS2JPSIPHI)
      || (fMode == BD2JPSIKSTAR)
      || (fMode == BU2JPSIKP)) {
    if (fGoodHLT
	&& fGoodAcceptance
	&& fGoodGlobalMuonsKin
	&& fGoodTracksPt && fGoodTracksEta
	&& fGoodCmssw
	&& fGoodCandAna
	&& fGoodJpsiCuts
	&& fGoodDcand
	&& fGoodBdtPresel
	&& fGoodGlobalMuonsKin
	) {
      fhMassSysStep0[fHistStrings[1]][fChan]->Fill(mass, -0.1, ps);
      fhMassSysStep0[fHistStrings[1]][fChan]->Fill(mass, 0.1);
      fhMassSysStep0[fHistStrings[1]][fChan]->Fill(mass, ps+0.1);
      fhMassSysStep0C[fHistStrings[1]][fChan]->Fill(fb.cm, -0.1, ps);
      fhMassSysStep0C[fHistStrings[1]][fChan]->Fill(fb.cm, 0.1);
      fhMassSysStep0C[fHistStrings[1]][fChan]->Fill(fb.cm, ps+0.1);

      if (fGoodBDT) {
	fhMassSysStep1[fHistStrings[1]][fChan]->Fill(mass, -0.1, ps);
	fhMassSysStep1[fHistStrings[1]][fChan]->Fill(mass, 0.1);
	fhMassSysStep1[fHistStrings[1]][fChan]->Fill(mass, ps+0.1);
	fhMassSysStep1C[fHistStrings[1]][fChan]->Fill(fb.cm, -0.1, ps);
	fhMassSysStep1C[fHistStrings[1]][fChan]->Fill(fb.cm, 0.1);
	fhMassSysStep1C[fHistStrings[1]][fChan]->Fill(fb.cm, ps+0.1);
      }
    }
  }


  if (fGoodQ
      && fGoodPvAveW8
      && fGoodTracks
      && fGoodTracksPt
      && fGoodTracksEta
      && fGoodMuonsPt  // PidTables do not really work below 4 GeV!!
      && fGoodMuonsEta
      && fGoodJpsiCuts
      && fGoodCmssw
      && fGoodCandAna
      && fGoodBdtPresel
      ) {

    if (fGoodMuonsID && fGoodDcand && fGoodHLT) {
      fhBdt[fChan]->Fill(fBDT);
    }

    for (unsigned int i = 1; i < fHistStrings.size(); ++i) {
      if (1 == i) {
	goodBDT = fGoodBDT;
      } else {
	goodBDT = (fBDT > (i-2)*0.01);
      }

      if (goodBDT) {
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
			<< " fW8AmMisId = " << fW8AmMisId
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
	// -- ams anti-muon selection
	if (fGoodAmsID && fGoodDcand) {
	  fhAmMassWithMuonCuts[fHistStrings[i]][fChan]->Fill(mass);
	  if (fGoodHLT) {
	    fhAmMassWithAllCuts[fHistStrings[i]][fChan]->Fill(mass);
	    fhAmW8MassWithAllCuts[fHistStrings[i]][fChan]->Fill(mass, fW8AmMisId);
	  }
	}
	// -- now the "real" muon ID (to extract the "real" trigger efficiency also for backgrounds)
	if (fGoodRealMuonsID && fGoodDcand) {
	  fhRMMassWithMuonCuts[fHistStrings[i]][fChan]->Fill(mass);
	  if (fb.hlt1 && fb.tos && fb.l1t) {
	    fhRMMassWithAllCuts[fHistStrings[i]][fChan]->Fill(mass);
	  }
	}

	if (fGoodRealAmsID && fGoodDcand) {
	  fhRAmMassWithMuonCuts[fHistStrings[i]][fChan]->Fill(mass);
	  if (fb.hlt1 && fb.tos && fb.l1t) {
	    fhRAmMassWithAllCuts[fHistStrings[i]][fChan]->Fill(mass);
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
    if (jentry%step == 0) {
      TTimeStamp ts;
      cout << Form(" .. evt = %d", jentry)
	   << ", time now: " << ts.AsString("lc")
	   << endl;
    }

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
    cout << Form("%2d %20s: %70s %15s %8.1f %.2e %.2e+/-%.2e",
		 cnt,
		 it->first.c_str(),
		 it->second->fF->GetName(),
		 it->second->fLatexName.c_str(),
		 it->second->fLumi,
		 it->second->fFilterEff,
		 it->second->fBf,
		 it->second->fBfE
		 )
	 << endl;
    ++cnt;
  }
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void plotResults::saveHistograms(string smode, TFile *f) {
  cout << "plotResults::saveHistograms> saving " << smode << " to " << fHistFile->GetName() << endl;
  TFile *fFile(0);
  if (0 == f) {
    fFile = fHistFile;
  } else {
    fFile = f;
  }
  fFile->cd();
  fFile->mkdir(fSample.c_str());
  fFile->cd(fSample.c_str());
  TDirectory *dir = gDirectory;

  TH1D *h1(0);
  TH2D *h2(0);

  // vector<string> modifier;
  // modifier.push_back("cnc" + fSuffix);
  // modifier.push_back("bdt" + fSuffix);

  for (int i = 0; i < fNchan; ++i) {
    h1 = (TH1D*)(fhBdt[i]->Clone(Form("hBdt_%s_chan%d", smode.c_str(), i)));
    h1->SetTitle(Form("hBdt_%s_chan%d %s", smode.c_str(), i, smode.c_str()));
    h1->SetDirectory(dir);
    h1->Write();
  }

  for (unsigned int im = 0; im < fHistStrings.size(); ++im) {
    for (int i = 0; i < fNchan; ++i) {
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

      if (im < 2) {
	h1 = (TH1D*)(fhMassSysStep0[fHistStrings[im]][i]->Clone(Form("hMassSysStep0_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
	h1->SetTitle(Form("hMassSysStep0_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
	h1->SetDirectory(dir);
	h1->Write();

	h1 = (TH1D*)(fhMassSysStep0C[fHistStrings[im]][i]->Clone(Form("hMassSysStep0C_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
	h1->SetTitle(Form("hMassSysStep0C_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
	h1->SetDirectory(dir);
	h1->Write();

	h1 = (TH1D*)(fhMassSysStep1[fHistStrings[im]][i]->Clone(Form("hMassSysStep1_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
	h1->SetTitle(Form("hMassSysStep1_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
	h1->SetDirectory(dir);
	h1->Write();

	h1 = (TH1D*)(fhMassSysStep1C[fHistStrings[im]][i]->Clone(Form("hMassSysStep1C_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
	h1->SetTitle(Form("hMassSysStep1C_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
	h1->SetDirectory(dir);
	h1->Write();
      }

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

      h1 = (TH1D*)(fhAmMassWithAllCuts[fHistStrings[im]][i]->Clone(Form("hAmMassWithAllCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hAmMassWithAllCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhAmW8MassWithAllCuts[fHistStrings[im]][i]->Clone(Form("hAmW8MassWithAllCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hAmW8MassWithAllCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhAmMassWithMuonCuts[fHistStrings[im]][i]->Clone(Form("hAmMassWithMuonCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hAmMassWithMuonCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();


      h1 = (TH1D*)(fhRMMassWithMuonCuts[fHistStrings[im]][i]->Clone(Form("hRMMassWithMuonCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hRMMassWithMuonCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhRMMassWithAllCuts[fHistStrings[im]][i]->Clone(Form("hRMMassWithAllCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hRMMassWithAllCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhRAmMassWithMuonCuts[fHistStrings[im]][i]->Clone(Form("hRAmMassWithMuonCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hRAmMassWithMuonCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
      h1->SetDirectory(dir);
      h1->Write();

      h1 = (TH1D*)(fhRAmMassWithAllCuts[fHistStrings[im]][i]->Clone(Form("hRAmMassWithAllCuts_%s_%s_chan%d", fHistStrings[im].c_str(), smode.c_str(), i)));
      h1->SetTitle(Form("hRAmMassWithAllCuts_%s_%s_chan%d %s", fHistStrings[im].c_str(), smode.c_str(), i, smode.c_str()));
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
  hcuts->SetDirectory(fFile);
  cout << "writing hcuts " << hcuts << " name: " << hcuts->GetName() << " to " << fFile->GetName() << endl;
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

      if (im < 2) {
	fhMassSysStep0[fHistStrings[im]][i]->Reset();
	if (deleteThem) delete fhMassSysStep0[fHistStrings[im]][i];

	fhMassSysStep0C[fHistStrings[im]][i]->Reset();
	if (deleteThem) delete fhMassSysStep0C[fHistStrings[im]][i];

	fhMassSysStep1[fHistStrings[im]][i]->Reset();
	if (deleteThem) delete fhMassSysStep1[fHistStrings[im]][i];

	fhMassSysStep1C[fHistStrings[im]][i]->Reset();
	if (deleteThem) delete fhMassSysStep1C[fHistStrings[im]][i];
      }

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

      fhAmMassWithAllCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAmMassWithAllCuts[fHistStrings[im]][i];
      fhAmW8MassWithAllCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAmW8MassWithAllCuts[fHistStrings[im]][i];
      fhAmMassWithMuonCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhAmMassWithMuonCuts[fHistStrings[im]][i];


      fhRAmMassWithAllCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhRAmMassWithAllCuts[fHistStrings[im]][i];
      fhRAmMassWithMuonCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhRAmMassWithMuonCuts[fHistStrings[im]][i];
      fhRMMassWithAllCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhRMMassWithAllCuts[fHistStrings[im]][i];
      fhRMMassWithMuonCuts[fHistStrings[im]][i]->Reset();
      if (deleteThem) delete fhRMMassWithMuonCuts[fHistStrings[im]][i];


    }
  }
}


// ----------------------------------------------------------------------
void plotResults::dumpTex(number &a, const string &s, int ndigits, int sgn) {
  fTEX << formatTex(a.val,   Form("%s:val", s.c_str()), ndigits) << endl;
  fTEX << formatTex(a.estat, Form("%s:estat", s.c_str()), ndigits) << endl;
  if (a.estatLo > 0.) {
    fTEX << formatTex(a.estatLo, Form("%s:estlo", s.c_str()), ndigits) << endl;
    fTEX << formatTex(a.estatHi, Form("%s:esthi", s.c_str()), ndigits) << endl;
  }
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
  n = a.fEffTrigMCAm;     OUT << "    fEffTrigMCAm =     " << Form("%8.6f +/- %8.6f +/- %8.6f", n.val,  n.estat, n.esyst) << endl;
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


// ----------------------------------------------------------------------
double plotResults::findVarEtot(string varName, vector<string> &lines) {
  for (unsigned int i = 0; i < lines.size(); ++i) {
    if (string::npos != lines[i].find(varName) && string::npos != lines[i].find("etot")) {
      string::size_type m1 = lines[i].find("ensuremath{");
      string::size_type m2 = lines[i].find("}", m1);
      string snum = lines[i].substr(m1+12, m2-m1-12-1);
      return atof(snum.c_str());
    }
  }
  return 0.;
}

// ----------------------------------------------------------------------
double plotResults::findVarEstat(string varName, vector<string> &lines) {
  for (unsigned int i = 0; i < lines.size(); ++i) {
    if (string::npos != lines[i].find(varName) && string::npos != lines[i].find("estat")) {
      string::size_type m1 = lines[i].find("ensuremath{");
      string::size_type m2 = lines[i].find("}", m1);
      string snum = lines[i].substr(m1+12, m2-m1-12-1);
      return atof(snum.c_str());
    }
  }
  return 0.;
}
