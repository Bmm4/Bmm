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
plotWork::plotWork(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotWork::loadFiles(files);

  changeSetup(dir, "plotWork", setup);

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

}


// ----------------------------------------------------------------------
plotWork::~plotWork() {

}

// ----------------------------------------------------------------------
string plotWork::removeVarFromSelection(string var, string selection) {
  // this will split on &, not on &&. So every second element in cuts is empty
  vector<string> cuts = split(selection, '&');
  string redcuts("");
  for (unsigned int i = 0; i < cuts.size(); ++i) {
    if (cuts[i] == "") continue;
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
    selection += " && hlt";
  }

  if (fSample.find("psi")) {
    selection += " && psimaxdoca<0.5 && mpsi>2.9 && mpsi<3.3 && psipt>6.9";
  }

  return selection;
}


// ----------------------------------------------------------------------
void plotWork::init() {
  system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
}


// ----------------------------------------------------------------------
void plotWork::makeAll(int bitmask) {



  if (bitmask & 0x1) {
    runTisEfficiency("bupsikData");
    runTisEfficiency("bspsiphiData");
    runTisEfficiency("bdpsikstarData");
    runTisEfficiency("bmmData");

    runTisEfficiency("bupsikSingleMuon");
    runTisEfficiency("bspsiphiSingleMuon");
    runTisEfficiency("bdpsikstarSingleMuon");
    runTisEfficiency("bmmSingleMuon");
  }

  if (bitmask & 0x2) {
    plotTisEfficiency("all");
  }

  if (bitmask & 0x4) {
    efficiencyVariable("all", "hlt", 10, 0, 0, 0, "bupsikMc");
    efficiencyVariable("all", "reftrg", 10, 0, 0, 0, "bupsikMc");

    efficiencyVariable("all", "hlt", 11, 0, 0, 0, "bupsikMc");
    efficiencyVariable("all", "reftrg", 11, 0, 0, 0, "bupsikMc");

    efficiencyVariable("all", "hlt", 12, 0, 0, 0, "bupsikMc");

  }
}


// ----------------------------------------------------------------------
void plotWork::bookHist(string dsname) {

  fpHnorm = new TH1D(Form("h_%s_%s", "norm", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);
  fpHpass = new TH1D(Form("h_%s_%s", "pass", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);

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
  fCds = fSample;
  loopOverTree(t, 1);


  fpHnorm->Write();
  fpHpass->Write();
  fHistFile->Close();

}


// ----------------------------------------------------------------------
void plotWork::refTrgEfficiency(string selection, string dsname) {

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
  string selection = removeVarFromSelection(var, selectionString(iselection, 0));
  cout << "==> plotWork::efficiencyVariable> selection = " << selection << endl;

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

  string dir("");
  if (string::npos != dsname.find("bupsik")) {
    dir = "candAnaBu2JpsiK";
  }

  TTree *t = getTree(fSample, dir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
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
void plotWork::yieldStability(string dsname, string trg) {

  // -- dump histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;


  fSample = dsname;
  string dir = "candAnaBu2JpsiK";
  if (string::npos != fSample.find("bmm")) {
    dir = "candAnaMuMu";
  } else if (string::npos != fSample.find("bdpsikstar")) {
    dir = "candAnaBd2JpsiKstar";
  } else if (string::npos != fSample.find("bdpsikstar")) {
    dir = "candAnaBd2JpsiKstar";
  }

  // -- check whether there are any histograms pre-produced already
  bool ok = fHistFile->cd(dir.c_str());
  cout << "OK? " << ok << endl;
  if (ok) {
    cout << "histograms exist already, looping over them" << endl;
    TIter next(gDirectory->GetListOfKeys());
    TKey *key(0);
    TH1D *hHLT(0), *hRTR(0);
    vector<int> vds;
    int run(-1), runMin(9999999), runMax(0);
    while ((key = (TKey*)next())) {
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
      if (TString(key->GetName()).Contains("_HLT_")) {
	string hname = key->GetName();
	replaceAll(hname, "h_HLT_", "");
	run = atoi(hname.c_str());
	if (run > runMax) runMax = run;
	if (run < runMin) runMin = run;
	vds.push_back(run);
      }
    }

    if (vds.size() > 0) {
      Lumi lumi("../common/json/Cert_271036-275125_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
      cout << "runs " << runMin << " .. " <<  runMax << endl;
      TH1D *hRunHLT = new TH1D(Form("hRun%s", trg.c_str()), "", runMax-runMin+1, runMin, runMax); hRunHLT->Sumw2();

      double mBp(5.28), sBp(0.04), stepBp(5.145);
      double xmin(5.0), xmax(5.8), ymax(0.);
      fIF->fLo = xmin;
      fIF->fHi = xmax;
      for (unsigned int i = 0; i < vds.size(); ++i) {
	TH1D *h1 = (TH1D*)(gDirectory->Get(Form("h_%s_%d", trg.c_str(), vds[i])));
	if (h1->Integral(1, h1->GetNbinsX()+1) < 100) continue;
	TF1 *f1 = fIF->expoErrGauss(h1, mBp, sBp, stepBp);
	h1->Fit(f1, "lr", "", xmin, xmax);
	double nNormE = f1->IntegralError(5.1, 5.4)/f1->Integral(5.1, 5.4);
	if (nNormE < 0.1) nNormE = f1->GetParError(0)/f1->GetParameter(0);
	f1->SetParameter(3, 0.);
	f1->SetParameter(4, 0.);
	f1->SetParameter(5, 0.);
	f1->SetParameter(6, 0.);
	double nNorm = f1->Integral(5.1, 5.4)/h1->GetBinWidth(1);
	double result = nNorm/lumi.lumi(vds[i]);
	cout << "Run " << vds[i] << " bin: " << hRunHLT->FindBin(vds[i])
	     << " nNorm: " << nNorm << " +/- " << nNormE*nNorm
	     << " lumi: " << lumi.lumi(vds[i])
	     <<  " result: " << result
	     << endl;
	// normalize yield to 1/pb
	hRunHLT->SetBinContent(hRunHLT->FindBin(static_cast<double>(vds[i])), result);
	hRunHLT->SetBinError(hRunHLT->FindBin(static_cast<double>(vds[i])), nNormE*nNorm/lumi.lumi(vds[i]));
	if (result > ymax) ymax = result;
	tl->DrawLatexNDC(0.2, 0.92, Form("Nsig = %4.1f +/- %4.1f", nNorm, nNormE*nNorm));
	tl->DrawLatexNDC(0.2, 0.85, Form("%d", vds[i]));
	c0->Modified();
	c0->Update();
	savePad(Form("yield-%s-%d.pdf", trg.c_str(), vds[i]));
      }
      hRunHLT->SetMinimum(0.);
      hRunHLT->SetMaximum(1.3*ymax);
      hRunHLT->Draw("e");
      savePad(Form("yield-vs-runs-%d-%s-%s.pdf", fYear, fSample.c_str(), trg.c_str()));
      return;
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
    fCds = fSample;
    loopOverTree(t, 2);

    cout << "writing output histograms: " << fYieldHLT.size() << endl;
    for (map<int, TH1D*>::iterator it = fYieldHLT.begin(); it != fYieldHLT.end(); ++it) {
      cout << "run " << it->first << endl;
      it->second->Draw();
      it->second->SetDirectory(hDir);
      it->second->Write();
    }
    for (map<int, TH1D*>::iterator it = fYieldRTR.begin(); it != fYieldRTR.end(); ++it) {
      it->second->Draw();
      it->second->SetDirectory(hDir);
      it->second->Write();
    }
  }

  fHistFile->Write();
  fHistFile->Close();
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


  if ((fMode == BU2JPSIKP) || (fMode = BD2JPSIKSTAR) || (fMode = BS2JPSIPHI)) {
    if (TMath::Abs(fb.mpsi) < 2.9) return;
    if (TMath::Abs(fb.mpsi) > 3.3) return;

    if (TMath::Abs(fb.psipt) < 6.9) return;
    if (TMath::Abs(fb.psicosa) < 0.9) return;
    if (TMath::Abs(fb.psiprob) < 0.1) return;
    if (TMath::Abs(fb.psiflsxy) < 3) return;
  }

  fpHnorm->Fill(fb.m);

  if (!fb.hlt) return;
  fpHpass->Fill(fb.m);

}

// ----------------------------------------------------------------------
void plotWork::loopFunction2() {


  //  if (!fGoodMuonsID) return;

  if (!fGoodQ) return;
  if (!fGoodPvAveW8) return;
  if (!fGoodMaxDoca) return;

  if (!fb.json) return;

  if (!fGoodLip) return;
  if (!fGoodLipS) return;

  if (fb.docatrk > 0.15) return;
  if (fb.closetrk > 3) return;
  if (fb.iso < 0.7) return;
  if (fb.fls3d < 5) return;
  if (fb.chi2dof > 5.0) return;
  if (fb.alpha > 0.2) return;

  if (fb.pvip > 0.01) return;
  if (fb.pvips > 3) return;


  if (TMath::Abs(fb.flsxy) < 3.0) return;

  if (TMath::Abs(fb.m1eta) > 1.6) return;
  if (TMath::Abs(fb.m2eta) > 1.6) return;

  if (fb.m1pt < 4.0) return;
  if (fb.m2pt < 4.0) return;

  double m = fb.m;
  if ((fMode == BU2JPSIKP) || (fMode = BD2JPSIKSTAR) || (fMode = BS2JPSIPHI)) {
    if (TMath::Abs(fb.mpsi) < 2.9) return;
    if (TMath::Abs(fb.mpsi) > 3.3) return;

    if (TMath::Abs(fb.psipt) < 6.9) return;
    if (TMath::Abs(fb.psicosa) < 0.9) return;
    if (TMath::Abs(fb.psiprob) < 0.1) return;
    if (TMath::Abs(fb.psiflsxy) < 3) return;
    m = fb.cm;
  }


  if (0 == fYieldHLT.count(fb.run)) {
    TH1D *h = new TH1D(Form("h_HLT_%d", static_cast<int>(fb.run)), Form("h_HLT_%d", static_cast<int>(fb.run)), 60, 4.8, 6.0);
    fYieldHLT.insert(make_pair(fb.run, h));

    h = new TH1D(Form("h_RTR_%d", static_cast<int>(fb.run)), Form("h_RTR_%d", static_cast<int>(fb.run)), 60, 4.8, 6.0);
    fYieldRTR.insert(make_pair(fb.run, h));
  }


  if (fb.hlt) {
    fYieldHLT[fb.run]->Fill(m);
  }

  if (fb.reftrg) {
    fYieldRTR[fb.run]->Fill(m);
  }

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
  cout << "==> plotWork::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotWork::*pF)(void);
  if (ifunc == 1) pF = &plotWork::loopFunction1;
  if (ifunc == 2) pF = &plotWork::loopFunction2;

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

      if (string::npos != stype.find("bupsikBla,")) {
        sname = "bupsikMcBla";
        sdecay = "bupsik";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

      if (string::npos != stype.find("YYY")) {
        sname = "bupsikMc";
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
