#include "fitPsYield.hh"
#include "util.hh"
#include "initFunc.hh"

#include <sstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include "TClass.h"
#include "TMath.h"
#include "TStyle.h"
#include "TKey.h"

// ----------------------------------------------------------------------
// -- debugging example for TH2D
// -----------------------------
// TFile *_file0 = TFile::Open("results/plotResults-BF.2016.root")
// gFile->cd("bupsikData")
// string name("hNorm_bdt2016BF_bupsikData_chan0")
// fitPsYield a(name, 0)
// a.fitBu2JpsiKp(5, "results")
//
// -- debugging example for TH1D
// -----------------------------
// TFile *_file0 = TFile::Open("results/plotReducedOverlays-BF.2016.root");
// TH1D *h = (TH1D*)_file0->Get("ad1cnc_bspsiphiData_tauMassAo");
// fitPsYield a(0, 1);
// a.fit0_Bu2JpsiKp(h, -1, ".");
// ----------------------------------------------------------------------



using namespace std;

// ----------------------------------------------------------------------
fitPsYield::fitPsYield(string hname, TDirectory *pD, int verbose): fVerbose(verbose), fBaseName(hname), fpIF(new initFunc),
								   fCombined(0), fCombinedW8(0) {
  fData.clear();
  TDirectory *pDir = pD;
  if (0 == pDir) pDir = gDirectory;
  TIter next(pDir->GetListOfKeys());
  TKey *key(0);
  if (fVerbose > 0) cout << "fitPsYield: looking for histogram ->" << hname << "<- in directory " << pDir->GetName() << endl;
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
    if (TString(key->GetName()).Contains(hname.c_str())) {
      string hname = key->GetName();
      TH2D *h2     = (TH2D*)pDir->Get(hname.c_str());
      cout << "found hname ->" << hname << "<- and h2 = " << h2 << endl;
      if (h2) {
	initFromHist(h2);
	break;
      }
    }
  }

}


// ----------------------------------------------------------------------
fitPsYield::fitPsYield(TH2D *h2, int verbose): fVerbose(verbose), fpIF(new initFunc), fCombined(0), fCombinedW8(0) {
  if (!h2) return;
  fData.clear();

  fBaseName = h2->GetName();
  initFromHist(h2);
}


// ----------------------------------------------------------------------
fitPsYield::~fitPsYield() {
  // for (vector<psd*>::iterator it = fData.begin(); it != fData.end(); ++it) {
  //   delete *it;
  // }
  fData.clear();
}

// ----------------------------------------------------------------------
void fitPsYield::initFromHist(TH2D *h2) {

  if (!h2) {
    cout << "histogram does not exist, returning!" << endl;
    fH2 = 0;
    return;
  }
  if (fVerbose > 0) cout << "fitPsYield: found histogram " << h2->GetName()
			 << " with integral = " << h2->Integral(1, h2->GetNbinsX(), 1, h2->GetNbinsY())
			 << endl;
  string hname = h2->GetName();
  fH2          = (TH2D*)h2->Clone(Form("fpy_%s", hname.c_str()));
  int chanB    = h2->GetYaxis()->FindBin(+0.1);
  int chan1    = h2->GetYaxis()->FindBin(1.1);
  int chan2    = h2->GetYaxis()->FindBin(h2->GetYaxis()->GetXmax()) - 1;
  fCombined    = h2->ProjectionX(Form("fpy_%s_comb", hname.c_str()), 2, 2);
  fCombinedW8  = h2->ProjectionX(Form("fpy_%s_combW8", hname.c_str()), 1, 1);
  // -- fW8Combined
  double ent   = fH2->Integral(1, fH2->GetNbinsX(), 1, 1);
  if (ent < 1) {
    cout << "no entries in W8Combined found, deleting and returning" << endl;
    delete fH2;
    delete fCombined;
    delete fCombinedW8;
    return;
  }
  fW8Combined             = new psd();
  fW8Combined->fPs        = -1;
  fW8Combined->fEntries   = ent;
  fW8Combined->fH1        = h2->ProjectionX(Form("fpy_%s_w8", hname.c_str()), 1, 1);
  fW8Combined->fH1->SetTitle(Form("fpy_%s_w8", hname.c_str()));
  // -- fUnW8Combined
  ent   = fH2->Integral(1, fH2->GetNbinsX(), 2, 2);
  if (ent < 1) {
    cout << "no entries in U8Combined found, deleting and returning" << endl;
    delete fH2;
    delete fW8Combined;
    delete fCombined;
    delete fCombinedW8;
    return;
  }
  fUnW8Combined           = new psd();
  fUnW8Combined->fPs      = 0;
  fUnW8Combined->fEntries = ent;
  fUnW8Combined->fH1      = h2->ProjectionX(Form("fpy_%s_u8", hname.c_str()), 2, 2);
  fUnW8Combined->fH1->SetTitle(Form("fpy_%s_u8", hname.c_str()));
  // -- per prescale value
  for (int ips = chan1; ips <= chan2; ++ips) {
    double ent   = fH2->Integral(1, fH2->GetNbinsX(), ips, ips);
    if (ent < 1) continue;
    psd *a       = new psd();
    a->fPs       = ips - chanB;
    a->fEntries  = ent;
    a->fH1       = h2->ProjectionX(Form("fpy_%s_ps%d", hname.c_str(), a->fPs), ips, ips);
    a->fH1->SetTitle(Form("fpy_%s_ps%d", hname.c_str(), a->fPs));
    fData.push_back(a);
  }

}



// ----------------------------------------------------------------------
void fitPsYield::printSummary() {
  for (unsigned int i = 0; i < fData.size(); ++i) {
    cout << Form("idx %2d   sg = %8.2f +/- %5.2f", i, fData[i]->fResults.fSg, fData[i]->fResults.fSgE)
	 << " ps = " << fData[i]->fPs
	 << Form(" (rel err = %3.2f)", fData[i]->fResults.fSgE/fData[i]->fResults.fSg)
	 << Form(" fit chi2/dof = %5.1f/%3d, prob = %4.3f", fData[i]->fChi2, fData[i]->fNdof, fData[i]->fProb)
	 << (fData[i]->fProb < 0.05? "**":"")
	 << endl;
  }

  cout << Form("  total = %8.1f +/- %3.2f", fSummary.fSg, fSummary.fSgE)
       << Form(" (rel err = %6.5f)", fSummary.fSgE/fSummary.fSg)
       << endl;
  cout << Form("w8c  sg = %8.1f +/- %3.2f", fW8Combined->fResults.fSg, fW8Combined->fResults.fSgE)
       << Form(" (rel err = %3.2f)", fW8Combined->fResults.fSgE/fW8Combined->fResults.fSg)
       << Form(" fit chi2/dof = %5.1f/%3d, prob = %4.3f", fW8Combined->fChi2, fW8Combined->fNdof, fW8Combined->fProb)
       << (fW8Combined->fProb < 0.05? "**":"")
       << endl;
  cout << Form("u8c  sg = %8.1f +/- %3.2f", fUnW8Combined->fResults.fSg, fUnW8Combined->fResults.fSgE)
       << Form(" (rel err = %3.2f)", fUnW8Combined->fResults.fSgE/fUnW8Combined->fResults.fSg)
       << Form(" fit chi2/dof = %5.1f/%3d, prob = %4.3f", fUnW8Combined->fChi2, fUnW8Combined->fNdof, fUnW8Combined->fProb)
       << (fUnW8Combined->fProb < 0.05? "**":"")
       << endl;

}



// ----------------------------------------------------------------------
void fitPsYield::fitBu2JpsiKp(int limitpars, string pdfprefix, int whichfit) {
  if (0 == fH2) {
    cout << "no histogram found/setup/defined, returning!" << endl;
    return;
  }

  void (fitPsYield::*pF)(psd *, int, string, bool);
  if (0 == whichfit) pF = &fitPsYield::fit0_Bu2JpsiKp;
  if (1 == whichfit) pF = &fitPsYield::fit1_Bu2JpsiKp;

  if (0 == fCombinedW8) {
    cout << "ERROR> no psd fCombinedW8 created (likely no entries in histogram), returning" << endl;
    fSummary.fSg      = 0.;
    fSummary.fSgE     = 0.;
    fSummary.fBg      = 0.;
    fSummary.fBgE     = 0.;
    fSummary.fSgSigma = 0.;
    fSummary.fSgPeak  = 0.;
    return;
  }
  // -- prefit weighted combination:
  (this->*pF)(fW8Combined, -1, pdfprefix, false);
  // -- prefit unweighted combination:
  (this->*pF)(fUnW8Combined, -1, pdfprefix, false);

  // -- fit all prescales
  for (unsigned int ips = 0; ips < fData.size(); ++ips) {
    (this->*pF)(fData[ips], limitpars, pdfprefix, false);
  }

  fSummary.clear();
  for (unsigned int i = 0; i < fData.size(); ++i) {
    // cout << " -> adding ps = " << fData[i]->fPs
    // 	 << " signal = " << fData[i]->fResults.fSg
    // 	 << " +/- " << fData[i]->fResults.fSgE;
    fSummary.fSg += fData[i]->fPs*fData[i]->fResults.fSg;
    fSummary.fSgE += fData[i]->fPs*fData[i]->fPs*fData[i]->fResults.fSgE*fData[i]->fResults.fSgE;
    // cout << " -> " << fSummary.fSg << " +/- " << TMath::Sqrt(fSummary.fSgE)
    //   	 << endl;

  }
  fSummary.fSgE = TMath::Sqrt(fSummary.fSgE);
  // cout << " => total: " << fSummary.fSg << " +/- " << fSummary.fSgE << endl;
  fSummary.fBg      = fUnW8Combined->fResults.fBg;
  fSummary.fBgE     = fUnW8Combined->fResults.fBgE;
  fSummary.fSgSigma = fUnW8Combined->fResults.fSgSigma;
  fSummary.fSgPeak  = fUnW8Combined->fResults.fSgPeak;
  printSummary();
}


// ----------------------------------------------------------------------
void fitPsYield::fit0_Bu2JpsiKp(TH1D *h1, int limitpars, string pdfprefix) {
  fUnW8Combined           = new psd();
  fUnW8Combined->fPs      = 0;
  fUnW8Combined->fEntries = h1->Integral(1, h1->GetNbinsX());
  fUnW8Combined->fH1      = h1;
  fData.push_back(fUnW8Combined);

  fit0_Bu2JpsiKp(fUnW8Combined, limitpars, pdfprefix, true);

  fSummary.fSg      = fUnW8Combined->fResults.fSg;
  fSummary.fSgE     = fUnW8Combined->fResults.fSgE;
  fSummary.fBg      = fUnW8Combined->fResults.fBg;
  fSummary.fBgE     = fUnW8Combined->fResults.fBgE;
  fSummary.fSgSigma = fUnW8Combined->fResults.fSgSigma;
  fSummary.fSgPeak  = fUnW8Combined->fResults.fSgPeak;

}


// ----------------------------------------------------------------------
// limitpars < 0: determine starting values for fit
//           = 0: only allow normalizations to float
//           > 0: constrain parameters within limitpars*sigma of prior (limitpars < 0) call
void fitPsYield::fit0_Bu2JpsiKp(psd *res, int limitpars, string pdfprefix, bool keepFunctions) {
  if (0 == res) {
    cout << "ERROR: calling for empty psd (prescale data), returning." << endl;
    return;
  }
  TH1D *h = res->fH1;
  setTitles(h, "#it{m}_{#it{#mu#mu K}} [GeV]", Form("Candidates/(%4.3f GeV)", h->GetBinWidth(1)), 0.05, 1.1, 1.9);
  if (0 == fData.size()) return;
  fSummary.clear();

  fpIF->fName = "fit";
  fpIF->fVerbose = true;
  TF1 *f1 = fpIF->bupsik(h);
  fpIF->fName = "comp";
  TF1 *fcnSig  = fpIF->gauss2c(h);
  fcnSig->SetLineColor(kBlue+1);
  TF1 *fcnExpo = fpIF->expo(0., 100.);
  fcnExpo->SetLineColor(kRed+1);
  fcnExpo->SetLineStyle(kSolid);
  TF1 *fcnErr2 = fpIF->err2(0., 100.);
  fcnErr2->SetLineColor(kRed+2);
  fcnErr2->SetLineStyle(kSolid);
  TF1 *fcnSat = fpIF->pisat(1.);
  fcnSat->SetLineColor(kRed+3);
  fcnSat->SetLineStyle(kSolid);

  TCanvas *c0(0);
  c0 = (TCanvas*)gROOT->FindObject("fpy_c0");
  if (!c0) c0 = new TCanvas("fpy_c0","--c0--",0,0,700,700);
  c0->cd();
  c0->Clear();
  shrinkPad(0.13, 0.19);

  fpIF->fLo = 5.0;
  fpIF->fHi = 6.0;
  cout << "==> fitPsYield::fit0_Bu2JpsiKp> FITTING " << h->GetName() << " with limitpars = " << limitpars << endl;
  double xmin(4.9), xmax(6.0), expoLo(5.15), expoHi(6.0);
  if (fVerbose) cout << "h->GetSumOfWeights() = " << h->GetSumOfWeights() << " h->GetEntries() = " << h->GetEntries() << endl;
  //  if (fVerbose) fpIF->dumpParameters(f1);
  string fitopt = "lr";
  if (0 == fVerbose) fitopt += "q";
  h->SetMinimum(0.);
  h->SetAxisRange(fpIF->fLo, fpIF->fHi);
  if (limitpars < 0) {
    fCombMax = h->GetMaximum();
    f1->SetParLimits(0, 0., 10.*h->GetMaximum());
    f1->SetParLimits(1, 5.25,  5.30); // it's a B+ fit!
    //    f1->SetParLimits(5, 0.,  1.e15);
    f1->SetParLimits(9, 0.,  10.*h->GetMaximum());
  } else {
    for (int ipar = 0; ipar < f1->GetNpar(); ++ipar) {
      f1->ReleaseParameter(ipar);
    }
    if (fCombMax < 1) {
      cout << "FIXME: fCombMax = " << fCombMax << endl;
      fCombMax = h->GetMaximum();
    }
    double scale = h->GetMaximum()/fCombMax;
    if (fVerbose) cout << "==> limitpars = " << limitpars << endl;
    if (limitpars > 0) {
      f1->SetParameter(0, scale*fPar[0]);  f1->SetParLimits(0, 0.,                           scale*(fPar[0] + fParE[0]));
      f1->SetParameter(1, fPar[1]);        f1->SetParLimits(1, fPar[1] - limitpars*fParE[1], fPar[1] + limitpars*fParE[1]);
      f1->SetParameter(2, fPar[2]);        f1->SetParLimits(2, fPar[2] - limitpars*fParE[2], fPar[2] + limitpars*fParE[2]);
      //      f1->SetParameter(3, fPar[3]);        f1->SetParLimits(3, fPar[3] - limitpars*fParE[3], fPar[3] + limitpars*fParE[3]);
      f1->FixParameter(3, fPar[3]);
      f1->SetParameter(4, fPar[4]);        f1->SetParLimits(4, fPar[4] - limitpars*fParE[4], fPar[4] + limitpars*fParE[4]);
      f1->SetParameter(5, scale*fPar[5]);  f1->SetParLimits(5, 0.,                           scale*(fPar[5] + fParE[5]));
      f1->SetParameter(6, fPar[6]);        f1->SetParLimits(6, fPar[6] - limitpars*fParE[6], fPar[6] + limitpars*fParE[6]);
      f1->SetParameter(7, fPar[7]);        f1->SetParLimits(7, fPar[7] - limitpars*fParE[7], fPar[7] + limitpars*fParE[7]);
      f1->SetParameter(8, fPar[8]);        f1->SetParLimits(8, fPar[8] - limitpars*fParE[8], fPar[8] + limitpars*fParE[8]);
      f1->SetParameter(9, scale*fPar[9]);  f1->SetParLimits(9, 0.,                           scale*(fPar[9] + fParE[9]));
    } else if (limitpars == 0) {
      f1->SetParameter(0, scale*fPar[0]);  f1->SetParLimits(0,  0.,                           1.e15);
      f1->FixParameter(1, fPar[1]);
      f1->FixParameter(2, fPar[2]);
      f1->FixParameter(3, fPar[3]);
      f1->FixParameter(4, fPar[4]);
      f1->SetParameter(5, scale*fPar[5]);   f1->SetParLimits(5,  0.,                           1.e15);
      f1->FixParameter(6, fPar[6]);
      f1->FixParameter(7, fPar[7]);
      f1->FixParameter(8, fPar[8]);
      f1->SetParameter(9, scale*fPar[9]);  f1->SetParLimits(9, 0.,                           1.e15);
    } else {
      // do nothing, boot strap!
    }
  }
  if (fVerbose > -1) fpIF->dumpParameters(f1);
  if (h->GetSumOfWeights() > 100) {
    if (h->GetSumOfWeights() > h->GetEntries()) {
      fitopt += "w";
      h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
    } else {
      h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
    }
  } else {
    h->Draw();
    if (limitpars < 0) {
      cout << "XX> you are screwed: fitting the (un)weighted combined data and less than 100 entries" << endl;
    } else {

    }
  }

  // -- store for possible later usage
  if (limitpars < 0) {
    fPar.clear();
    fParE.clear();
    for (int i = 0; i < f1->GetNpar(); ++i) {
      fPar.push_back(f1->GetParameter(i));
      fParE.push_back(f1->GetParError(i));
    }
  }

  // -- copy parameters into components and draw them
  for (int ipar = 0; ipar < fcnSig->GetNpar(); ++ipar) {
    fcnSig->SetParameter(ipar, f1->GetParameter(ipar));
    fcnSig->SetParError(ipar, f1->GetParError(ipar));
  }
  for (int ipar = 0; ipar < fcnExpo->GetNpar(); ++ipar) {
    fcnExpo->SetParameter(ipar, f1->GetParameter(ipar+5));
    fcnExpo->SetParError(ipar, f1->GetParError(ipar+5));
  }
  for (int ipar = 0; ipar < fcnErr2->GetNpar(); ++ipar) {
    fcnErr2->SetParameter(ipar, f1->GetParameter(ipar+7));
    fcnErr2->SetParError(ipar, f1->GetParameter(ipar+7));
  }


  double sqrt2pi = 2.506628275;
  double gintegral = sqrt2pi*f1->GetParameter(0)*(f1->GetParameter(2) + f1->GetParameter(3)*f1->GetParameter(4));
  gintegral = fcnSig->Integral(4.9, 6.0);
  double fracSat(0.04);
  double satintegral = initFunc::iF_int_pisat(1.);
  double norm        = fracSat*gintegral/satintegral;
  fcnSat->SetParameter(0, norm);

  cout << "f1->GetParameter(5): " << f1->GetParameter(5) << endl;
  cout << "NORM OF SIG: " << fcnSig->Integral(4.9, 6.0)
       << "NORM OF SAT: " << fcnSat->Integral(4.9, 6.0)
       << endl;

  fcnSig->Draw("same");
  fcnExpo->Draw("same");
  fcnErr2->Draw("same");
  fcnSat->Draw("same");

  // -- store for possible later usage
  if (limitpars < 0) {
    fPar.clear();
    fParE.clear();
    for (int i = 0; i < f1->GetNpar(); ++i) {
      fPar.push_back(f1->GetParameter(i));
      fParE.push_back(f1->GetParError(i));
    }
    double NSG   = fcnSig->Integral(5.1, 5.6)/h->GetBinWidth(1);
    double NTOT  = h->GetSumOfWeights();
    fCombS2All   = NSG/NTOT;
    fCombS2AllE  = TMath::Sqrt(1./NSG + 1./NTOT)*fCombS2All;
  }

  // -- other parameters of interest
  res->fChi2 = f1->GetChisquare();
  res->fNdof = f1->GetNDF();
  res->fProb = TMath::Prob(res->fChi2, res->fNdof);
  res->fResults.fSgPeak = f1->GetParameter(1);

  TH1D *hbla = new TH1D("hbla", "", 100, 4.5, 7.0);
  for (int i = 0; i < 10000; ++i) hbla->Fill(fcnSig->GetRandom());
  res->fResults.fSgSigma = hbla->GetRMS();
  delete hbla;

  // -- double gaussian integral over 3 sigma region
  double sig   = fcnSig->Integral(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
				  res->fResults.fSgPeak + 3.*res->fResults.fSgSigma);
  double sigE  = (fcnSig->GetParError(0)/fcnSig->GetParameter(0)) * sig;

  double bg   = f1->Integral(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
			     res->fResults.fSgPeak + 3.*res->fResults.fSgSigma)
    - sig ;
  double bgE  = fcnExpo->GetParError(0)/fcnExpo->GetParameter(0)*bg;

  cout << "XXXXX sig = " << sig << " +/- " << sigE
       << ", integrating from " << res->fResults.fSgPeak - 3.*res->fResults.fSgSigma << " to "
       << res->fResults.fSgPeak - 3.*res->fResults.fSgSigma
       << " XXXXX " << endl;

  // -- create 'sensible' errors
  sig  /= h->GetBinWidth(1);
  sigE /= h->GetBinWidth(1);
  bg   /= h->GetBinWidth(1);
  bgE  /= h->GetBinWidth(1);
  if (sigE > sig) {
    cout << "rescaled error: " << sigE << " (sig = " << sig << ") to ";
    sigE = TMath::Sqrt(sig);
    cout << sigE << endl;

  } else {
    cout << "integral error from " << res->fResults.fSgPeak - 3.*res->fResults.fSgSigma
	 << " .. " << res->fResults.fSgPeak + 3.*res->fResults.fSgSigma
	 << ": " << sigE << " (sig = " << sig << ")"
	 << endl;
  }
  res->fResults.fSg  = sig;
  res->fResults.fSgE = sigE;
  res->fResults.fBg  = bg;
  res->fResults.fBgE = bgE;

  TLatex tl;
  tl.SetTextSize(0.03);
  if (-1 == res->fPs) {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d, weighted)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  } else if (0 == res->fPs) {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d, unweighted)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  } else {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  }

  tl.SetTextAngle(90.);
  tl.DrawLatexNDC(0.93, 0.15, h->GetName());
  tl.SetTextAngle(0.);

  c0->Modified();
  c0->Update();
  c0->SaveAs(Form("%s%s.pdf", pdfprefix.c_str(), h->GetName()));
  delete c0;

  if (keepFunctions) {
    f1->SetName(Form("f1_%s", h->GetName()));        fFunctions.insert(make_pair(f1->GetName(), f1));
    fcnSig->SetName(Form("sig_%s", h->GetName()));   fFunctions.insert(make_pair(fcnSig->GetName(), fcnSig));
    fcnExpo->SetName(Form("expo_%s", h->GetName())); fFunctions.insert(make_pair(fcnExpo->GetName(), fcnExpo));
    fcnErr2->SetName(Form("err2_%s", h->GetName())); fFunctions.insert(make_pair(fcnErr2->GetName(), fcnErr2));
    fcnSat->SetName(Form("sat_%s", h->GetName()));   fFunctions.insert(make_pair(fcnSat->GetName(), fcnSat));
  } else {
    delete f1;
    delete fcnSig;
    delete fcnExpo;
    delete fcnErr2;
    delete fcnSat;
  }
}


// ----------------------------------------------------------------------
// limitpars < 0: determine starting values for fit
//           = 0: only allow normalizations to float
//           > 0: constrain parameters within limitpars*sigma of prior (limitpars < 0) call
void fitPsYield::fit1_Bu2JpsiKp(psd *res, int limitpars, string pdfprefix, bool keepFunctions) {
  TH1D *h = res->fH1;
  setTitles(h, "#it{m}_{#it{#mu#mu K}} [GeV]",
	    Form("Candidates/(%4.3f GeV)", h->GetBinWidth(1)), 0.05, 1.1, 1.9);
  if (0 == fData.size()) return;
  fSummary.clear();

  fpIF->fName = "fit";
  fpIF->fVerbose = false;
  TF1 *f1 = fpIF->bupsik(h);
  fpIF->fName = "comp";
  TF1 *fcnSig  = fpIF->crystalBallBiGauss(h);
  fcnSig->SetLineColor(kBlue+1);
  TF1 *fcnExpo = fpIF->expo(0., 100.);
  fcnExpo->SetLineColor(kRed+1);
  fcnExpo->SetLineStyle(kSolid);
  TF1 *fcnErr2 = fpIF->err2(0., 100.);
  fcnErr2->SetLineColor(kRed+2);
  fcnErr2->SetLineStyle(kSolid);
  TF1 *fcnSat = fpIF->pisat(1.);
  fcnSat->SetLineColor(kRed+3);
  fcnSat->SetLineStyle(kSolid);

  TCanvas *c0(0);
  c0 = (TCanvas*)gROOT->FindObject("c0");
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,700,700);
  c0->cd();
  c0->Clear();
  shrinkPad(0.13, 0.19);

  fpIF->fLo = 5.0;
  fpIF->fHi = 5.8;
  cout << "==> fitPsYield::fit0_Bu2JpsiKp> FITTING " << h->GetName()
       << " with limitpars = " << limitpars << endl;
  double xmin(5.0), xmax(5.8), expoLo(5.15), expoHi(5.85);
  if (fVerbose) cout << "h->GetSumOfWeights() = " << h->GetSumOfWeights()
		     << " h->GetEntries() = " << h->GetEntries() << endl;
  //  if (fVerbose) fpIF->dumpParameters(f1);
  string fitopt = "lr";
  if (0 == fVerbose) fitopt += "q";
  h->SetMinimum(0.);
  h->SetAxisRange(fpIF->fLo, fpIF->fHi);

  if (limitpars < 0) {
    fCombMax = h->GetMaximum();
  } else {
    for (int ipar = 0; ipar < f1->GetNpar(); ++ipar) {
      f1->ReleaseParameter(ipar);
    }
    if (fCombMax < 1) {
      cout << "FIXME: fCombMax = " << fCombMax << endl;
      fCombMax = h->GetMaximum();
    }
    double scale = h->GetMaximum()/fCombMax;
    if (fVerbose) cout << "==> limitpars = " << limitpars << endl;
    if (limitpars > 0) {
      f1->SetParameter(0, fPar[0]);
      f1->SetParLimits(0, 0., 1.);

      f1->SetParameter(1, fPar[1]);
      f1->SetParLimits(1, fPar[1] - limitpars*fParE[1], fPar[1] + limitpars*fParE[1]);

      f1->SetParameter(2, fPar[2]);
      f1->SetParLimits(2, fPar[2] - limitpars*fParE[2], fPar[2] + limitpars*fParE[2]);

      f1->SetParameter(3, fPar[3]);
      f1->SetParLimits(3, 0., 1.);

      f1->SetParameter(4, fPar[4]);
      f1->SetParLimits(4, fPar[4] + limitpars*fParE[4], fPar[4] + limitpars*fParE[4]);

      f1->SetParameter(5, fPar[5]);
      f1->SetParLimits(5, fPar[5] + limitpars*fParE[5], fPar[5] + limitpars*fParE[5]);

      f1->FixParameter(6, fPar[6]);

      f1->FixParameter(7, fPar[7]);

      f1->SetParameter(8, scale*fPar[8]);
      f1->SetParLimits(8, scale*fPar[8]*0.2, scale*fPar[8]*1.5);

      f1->FixParameter(9, fPar[9]);

      f1->SetParameter(10, scale*fPar[10]);
      f1->SetParLimits(10, scale*(fPar[10] - fParE[10]), scale*(fPar[10] + fParE[10]));

      f1->SetParameter(11, fPar[11]);
      f1->SetParLimits(11, fPar[11] + limitpars*fParE[11], fPar[11] + limitpars*fParE[11]);

      f1->SetParameter(12, fPar[12]);
      f1->SetParLimits(12, fPar[12] + limitpars*fParE[12], fPar[12] + limitpars*fParE[12]);

      f1->SetParameter(13, fPar[13]);
      f1->SetParLimits(13, fPar[13] + limitpars*fParE[13], fPar[13] + limitpars*fParE[13]);

      f1->SetParameter(14, scale*fPar[14]);f1->SetParLimits(14, 0., scale*(fPar[14] + fParE[14]));
    } else if (limitpars == 0) {
      f1->FixParameter(0, fPar[0]);
      f1->FixParameter(1, fPar[1]);
      f1->FixParameter(2, fPar[2]);
      f1->FixParameter(3, fPar[3]);
      f1->FixParameter(4, fPar[4]);
      f1->FixParameter(5, fPar[5]);
      f1->FixParameter(6, fPar[6]);
      f1->FixParameter(7, fPar[7]);
      f1->SetParameter(8, scale*fPar[8]); f1->SetParLimits(8, 0., 1.e7);
      f1->FixParameter(9, fPar[9]);
      f1->SetParameter(10, scale*fPar[10]); f1->SetParLimits(10, 0., 1.e7);
      f1->FixParameter(11, fPar[11]);
      f1->FixParameter(12, fPar[12]);
      f1->FixParameter(13, fPar[13]);
      f1->SetParameter(14, scale*fPar[14]); f1->SetParLimits(14, 0., 1.e7);
    } else {
      // do nothing, boot strap!
    }
    if (fVerbose > 1) fpIF->dumpParameters(f1);

  }

  if (h->GetSumOfWeights() > 100) {
    if (h->GetSumOfWeights() > h->GetEntries()) {
      fitopt += "w";
      h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
    } else {
      h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
    }
  } else {
    h->Draw();
    if (limitpars < 0) {
      cout << "XX> you are screwed: fitting the (un)weighted combined data and less than 100 entries" << endl;
    } else {

    }
  }

  // -- store for possible later usage
  if (limitpars < 0) {
    fPar.clear();
    fParE.clear();
    for (int i = 0; i < f1->GetNpar(); ++i) {
      fPar.push_back(f1->GetParameter(i));
      fParE.push_back(f1->GetParError(i));
    }
  }


  // -- copy parameters into components and draw them
  for (int ipar = 0; ipar < fcnSig->GetNpar(); ++ipar) {
    fcnSig->SetParameter(ipar, f1->GetParameter(ipar));
    fcnSig->SetParError(ipar, f1->GetParameter(ipar));
  }
  for (int ipar = 0; ipar < fcnExpo->GetNpar(); ++ipar) {
    fcnExpo->SetParameter(ipar, f1->GetParameter(ipar+10));
  }
  for (int ipar = 0; ipar < fcnErr2->GetNpar(); ++ipar) {
    fcnErr2->SetParameter(ipar, f1->GetParameter(ipar+12));
  }
  for (int ipar = 0; ipar < fcnSat->GetNpar(); ++ipar) {
    fcnSat->SetParameter(ipar, f1->GetParameter(ipar+9)*f1->GetParameter(ipar+8));
  }
  fcnSig->Draw("same");
  fcnExpo->Draw("same");
  fcnErr2->Draw("same");
  fcnSat->Draw("same");

  // -- store for possible later usage
  if (limitpars < 0) {
    fPar.clear();
    fParE.clear();
    for (int i = 0; i < f1->GetNpar(); ++i) {
      fPar.push_back(f1->GetParameter(i));
      fParE.push_back(f1->GetParError(i));
    }
    double NSG   = fcnSig->Integral(5.1, 5.6)/h->GetBinWidth(1);
    double NTOT  = h->GetSumOfWeights();
    fCombS2All   = NSG/NTOT;
    fCombS2AllE  = TMath::Sqrt(1./NSG + 1./NTOT)*fCombS2All;
  }

  // -- other parameters of interest
  res->fChi2 = f1->GetChisquare();
  res->fNdof = f1->GetNDF();
  res->fProb = TMath::Prob(res->fChi2, res->fNdof);
  res->fResults.fSgPeak = f1->GetParameter(4);

  TH1D *hbla = new TH1D("hbla", "", 100, 4.5, 7.0);
  for (int i = 0; i < 10000; ++i) hbla->Fill(fcnSig->GetRandom());
  res->fResults.fSgSigma = hbla->GetRMS();
  delete hbla;

  // -- double gaussian integral over 3 sigma region
  double sig   = fcnSig->Integral(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
				  res->fResults.fSgPeak + 3.*res->fResults.fSgSigma);
  double sigE  = TMath::Sqrt((fcnSig->GetParError(0)/fcnSig->GetParameter(0)) * (fcnSig->GetParError(0)/fcnSig->GetParameter(0))
			     + (fcnSig->GetParError(8)/fcnSig->GetParameter(8)) * (fcnSig->GetParError(8)/fcnSig->GetParameter(8))
			     )*sig;


  // -- create 'sensible' errors
  sig  /= h->GetBinWidth(1);
  sigE /= h->GetBinWidth(1);
  if (sigE > sig) {
    cout << "rescaled error: " << sigE << " (sig = " << sig << ") to ";
    sigE = TMath::Sqrt(sig);
    cout << sigE << endl;

  } else {
    cout << "integral error from " << res->fResults.fSgPeak - 3.*res->fResults.fSgSigma
	 << " .. " << res->fResults.fSgPeak + 3.*res->fResults.fSgSigma
	 << ": " << sigE << " (sig = " << sig << ")"
	 << endl;
  }
  res->fResults.fSg = sig;
  res->fResults.fSgE = sigE;

  TLatex tl;
  tl.SetTextSize(0.03);
  if (-1 == res->fPs) {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d, weighted)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  } else if (0 == res->fPs) {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d, unweighted)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  } else {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  }

  c0->Modified();
  c0->Update();
  cout << "fitPsYield::fit0_Bu2JpsiKp> printing "
       << Form("%s%s.pdf", pdfprefix.c_str(), h->GetName())
       << endl;
  c0->SaveAs(Form("%s%s.pdf", pdfprefix.c_str(), h->GetName()));

  delete f1;
  delete fcnSig;
  delete fcnExpo;
  delete fcnErr2;
  delete fcnSat;
}


// ----------------------------------------------------------------------
void fitPsYield::fitBs2JpsiPhi(int limitpars, string pdfprefix, int whichfit) {
  if (0 == fH2) {
    cout << "no histogram found/setup/defined, returning!" << endl;
    return;
  }

  void (fitPsYield::*pF)(psd *, int, string, bool);
  if (0 == whichfit) pF = &fitPsYield::fit0_Bs2JpsiPhi;
  //  if (1 == whichfit) pF = &fitPsYield::fit1_Bs2JpsiPhi;

  // -- prefit weighted combination:
  (this->*pF)(fW8Combined, -1, pdfprefix, false);
  // -- prefit unweighted combination:
  (this->*pF)(fUnW8Combined, -1, pdfprefix, false);

  // -- fit all prescales
  for (unsigned int ips = 0; ips < fData.size(); ++ips) {
    (this->*pF)(fData[ips], limitpars, pdfprefix, false);
  }

  fSummary.clear();
  for (unsigned int i = 0; i < fData.size(); ++i) {
    fSummary.fSg  += fData[i]->fPs*fData[i]->fResults.fSg;
    fSummary.fSgE += fData[i]->fPs*fData[i]->fPs*fData[i]->fResults.fSgE*fData[i]->fResults.fSgE;
  }
  fSummary.fSgE = TMath::Sqrt(fSummary.fSgE);
  fSummary.fBg      = fUnW8Combined->fResults.fBg;
  fSummary.fBgE     = fUnW8Combined->fResults.fBgE;
  fSummary.fSgSigma = fUnW8Combined->fResults.fSgSigma;
  fSummary.fSgPeak  = fUnW8Combined->fResults.fSgPeak;
}

// ----------------------------------------------------------------------
void fitPsYield::fit0_Bs2JpsiPhi(TH1D *h1, int limitpars, string pdfprefix) {
  fUnW8Combined           = new psd();
  fUnW8Combined->fPs      = 0;
  fUnW8Combined->fEntries = h1->Integral(1, h1->GetNbinsX());
  fUnW8Combined->fH1      = h1;
  fData.push_back(fUnW8Combined);

  fit0_Bs2JpsiPhi(fUnW8Combined, limitpars, pdfprefix, true);
}



// ----------------------------------------------------------------------
// limitpars < 0: determine starting values for fit
//           = 0: unconstrained fit, based on starting values of prior (limitpars < 0) call
//           > 0: constrain parameters within limitpars*sigma of prior (limitpars < 0) call
void fitPsYield::fit0_Bs2JpsiPhi(psd *res, int limitpars, string pdfprefix, bool keepFunctions) {
  TH1D *h = res->fH1;
  setTitles(h, "#it{m}_{#it{#mu#mu K K}} [GeV]",
	    Form("Candidates/(%4.3f GeV)", h->GetBinWidth(1)), 0.05, 1.1, 1.9);
  if (0 == fData.size()) return;
  fSummary.clear();

  fpIF->fName = "fit";
  TF1 *f1 = fpIF->expogauss2c(h, 5.37, 0.05, 1.1);
  fpIF->fName = "comp";
  TF1 *fg  = fpIF->gauss2c(0., 100.);
  fg->SetLineColor(kBlue+1);
  TF1 *fe = fpIF->expo(0., 100.);
  fe->SetLineColor(kRed+2);
  fe->SetLineStyle(kSolid);

  TCanvas *c0(0);
  c0 = (TCanvas*)gROOT->FindObject("fpy_c0");
  if (!c0) c0 = new TCanvas("fpy_c0","--c0--",0,0,700,700);
  c0->cd();
  c0->Clear();
  shrinkPad(0.13, 0.19);

  fpIF->fLo = 5.0;
  fpIF->fHi = 6.0;
  cout << "==> fitPsYield::fit0_Bs2JpsiPhi> FITTING " << h->GetName()
       << " with limitpars = " << limitpars << endl;
  double xmin(5.0), xmax(6.0), expoLo(5.15), expoHi(6.0);
  if (fVerbose) cout << "h->GetSumOfWeights() = " << h->GetSumOfWeights()
		     << " h->GetEntries() = " << h->GetEntries() << endl;
  //  if (fVerbose) fpIF->dumpParameters(f1);
  string fitopt = "lr";
  if (0 == fVerbose) fitopt += "q";
  h->SetMinimum(0.);
  h->SetAxisRange(fpIF->fLo, fpIF->fHi);

  if (limitpars < 0) {
    f1->SetParLimits(0, 0.,  1.e10);
    f1->SetParLimits(1, 5.3, 5.4); // this is a Bs fit, after all
    fCombMax = h->GetMaximum();
  } else {
    for (int ipar = 0; ipar < f1->GetNpar(); ++ipar) {
      f1->ReleaseParameter(ipar);
    }
    if (fCombMax < 1) {
      cout << "FIXME: fCombMax = " << fCombMax << endl;
      fCombMax = h->GetMaximum();
    }
    double scale = h->GetMaximum()/fCombMax;

    double s1lo = fPar[2] - limitpars*fParE[2];
    double s1hi = fPar[2] + limitpars*fParE[2];
    if (s1lo < 0.) s1lo = 0.001;
    double s2lo = fPar[4] - limitpars*fParE[4];
    if (s2lo < 0.) s2lo = 0.002;
    if (s2lo < s1hi) s2lo = 1.1*s1hi;
    if (fVerbose) cout << "==> limitpars = " << limitpars << endl;
    if (limitpars > 0) {
      f1->SetParameter(0, scale*fPar[0]);
      f1->SetParLimits(0, 0.,                           scale*(fPar[0] + fParE[0]));

      f1->SetParameter(1, fPar[1]);
      f1->SetParLimits(1, fPar[1] - limitpars*fParE[1], fPar[1] + limitpars*fParE[1]);

      f1->SetParameter(2, fPar[2]);
      f1->SetParLimits(2, s1lo,                         s1hi);

      f1->SetParameter(3, fPar[3]);
      f1->SetParLimits(3, fPar[3] - limitpars*fParE[3], fPar[3] + limitpars*fParE[3]);

      f1->SetParameter(4, fPar[4]);
      f1->SetParLimits(4, s2lo,                         fPar[4] + limitpars*fParE[4]);

      f1->SetParameter(5, scale*fPar[5]);
      f1->SetParLimits(5, 0.,                           scale*(fPar[5] + fParE[5]));

      f1->SetParameter(6, fPar[6]);
      f1->SetParLimits(6, fPar[6] - limitpars*fParE[6], fPar[6] + limitpars*fParE[6]);
    } else if (limitpars == 0) {
      f1->SetParameter(0, scale*fPar[0]);
      f1->SetParLimits(0, 0., 1.e10);
      f1->SetParameter(1, fPar[1]);
      f1->SetParameter(2, fPar[2]);
      f1->SetParameter(3, fPar[3]);
      f1->SetParameter(4, fPar[4]);
      f1->SetParameter(5, scale*fPar[5]);
      f1->SetParLimits(5, 0., 1.e15);
      f1->SetParameter(6, fPar[6]);
      f1->SetParLimits(6, -1.e7, 1.e7);
    } else {
      // do nothing, boot strap!
    }
    if (fVerbose > 1) fpIF->dumpParameters(f1);
  }

  if (h->GetSumOfWeights() > 100) {
    if (h->GetSumOfWeights() > h->GetEntries()) {
      fitopt += "w";
      h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
    } else {
      h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
    }
  } else {
    h->Draw();
    if (limitpars < 0) {
      cout << "XX> you are screwed: fitting the (un)weighted combined data and less than 100 entries" << endl;
    } else {

    }
  }

  // -- store for possible later usage
  if (limitpars < 0) {
    fPar.clear();
    fParE.clear();
    for (int i = 0; i < f1->GetNpar(); ++i) {
      fPar.push_back(f1->GetParameter(i));
      fParE.push_back(f1->GetParError(i));
    }
  }


  // -- copy parameters into components and draw them
  for (int ipar = 0; ipar < fg->GetNpar(); ++ipar) {
    fg->SetParameter(ipar, f1->GetParameter(ipar));
    fg->SetParError(ipar, f1->GetParError(ipar));
  }
  for (int ipar = 0; ipar < fe->GetNpar(); ++ipar) {
    fe->SetParameter(ipar, f1->GetParameter(ipar+5));
    fe->SetParError(ipar, f1->GetParError(ipar+5));
  }

  // -- other parameters of interest
  res->fChi2 = f1->GetChisquare();
  res->fNdof = f1->GetNDF();
  res->fProb = TMath::Prob(res->fChi2, res->fNdof);
  res->fResults.fSgPeak = f1->GetParameter(1);

  TH1D *hbla = new TH1D("hbla", "", 100, 4.5, 7.0);
  for (int i = 0; i < 10000; ++i) hbla->Fill(fg->GetRandom());
  res->fResults.fSgSigma = hbla->GetRMS();
  delete hbla;


  // -- double gaussian integral over 3 sigma region
  double sig   = fg->Integral(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
			      res->fResults.fSgPeak + 3.*res->fResults.fSgSigma);
  double sigE  = sig*(f1->GetParError(0)/f1->GetParameter(0));

  double bg   = f1->Integral(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
			     res->fResults.fSgPeak + 3.*res->fResults.fSgSigma)
    - sig ;
  double bgE  = fe->GetParError(0)/fe->GetParameter(0)*bg;

  // -- create 'sensible' errors
  sig  /= h->GetBinWidth(1);
  sigE /= h->GetBinWidth(1);
  bg   /= h->GetBinWidth(1);
  bgE  /= h->GetBinWidth(1);
  if (sigE > sig) {
    sigE = TMath::Sqrt(sig);
    cout << "rescaled error: " << sigE << " (sig = " << sig << ")"
	 << endl;
  } else {
    cout << "integral error from " << res->fResults.fSgPeak - 3.*res->fResults.fSgSigma
	 << " .. " << res->fResults.fSgPeak + 3.*res->fResults.fSgSigma
	 << ": " << sigE << " (sig = " << sig << ")"
	 << endl;
  }
  res->fResults.fSg  = sig;
  res->fResults.fSgE = sigE;
  res->fResults.fBg  = bg;
  res->fResults.fBgE = bgE;

  cout << "==> fitted " << h->GetTitle() << " signal = " << sig << " +/- " << sigE << endl;


  TLatex tl;
  tl.SetTextSize(0.03);
  if (-1 == res->fPs) {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d, weighted)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  } else if (0 == res->fPs) {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d, unweighted)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  } else {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  }

  c0->Modified();
  c0->Update();
  cout << "fitPsYield::fit0_Bs2JpsiPhi> printing "
       << Form("%s%s.pdf", pdfprefix.c_str(), h->GetName())
       << endl;
  c0->SaveAs(Form("%s%s.pdf", pdfprefix.c_str(), h->GetName()));

  if (keepFunctions) {
    f1->SetName(Form("f1_%s", h->GetName()));   fFunctions.insert(make_pair(f1->GetName(), f1));
    fg->SetName(Form("sig_%s", h->GetName()));  fFunctions.insert(make_pair(fg->GetName(), fg));
    fe->SetName(Form("expo_%s", h->GetName())); fFunctions.insert(make_pair(fe->GetName(), fe));
  } else {
    delete f1;
    delete fg;
    delete fe;
  }
}


// ----------------------------------------------------------------------
void fitPsYield::fitBd2JpsiKstar(int limitpars, string pdfprefix, int whichfit) {
  if (0 == fH2) {
    cout << "no histogram found/setup/defined, returning!" << endl;
    return;
  }

  void (fitPsYield::*pF)(psd *, int, string, bool);
  if (0 == whichfit) pF = &fitPsYield::fit0_Bd2JpsiKstar;

  // -- prefit weighted combination:
  (this->*pF)(fW8Combined, -1, pdfprefix, false);
  // -- prefit unweighted combination:
  (this->*pF)(fUnW8Combined, -1, pdfprefix, false);

  // -- fit all prescales
  for (unsigned int ips = 0; ips < fData.size(); ++ips) {
    (this->*pF)(fData[ips], limitpars, pdfprefix, false);
  }

  fSummary.clear();
  for (unsigned int i = 0; i < fData.size(); ++i) {
    fSummary.fSg  += fData[i]->fPs*fData[i]->fResults.fSg;
    fSummary.fSgE += fData[i]->fPs*fData[i]->fPs*fData[i]->fResults.fSgE*fData[i]->fResults.fSgE;
  }
  fSummary.fSgE = TMath::Sqrt(fSummary.fSgE);
}


// ----------------------------------------------------------------------
void fitPsYield::fit0_Bd2JpsiKstar(TH1D *h1, int limitpars, string pdfprefix) {
  fUnW8Combined           = new psd();
  fUnW8Combined->fPs      = 0;
  fUnW8Combined->fEntries = h1->Integral(1, h1->GetNbinsX());
  fUnW8Combined->fH1      = h1;
  fData.push_back(fUnW8Combined);

  fit0_Bd2JpsiKstar(fUnW8Combined, limitpars, pdfprefix, true);
}

// ----------------------------------------------------------------------
// limitpars < 0: determine starting values for fit
//           = 0: unconstrained fit, based on starting values of prior (limitpars < 0) call
//           > 0: constrain parameters within limitpars*sigma of prior (limitpars < 0) call
void fitPsYield::fit0_Bd2JpsiKstar(psd *res, int limitpars, string pdfprefix, bool keepFunctions) {
  TH1D *h = res->fH1;
  setTitles(h, "#it{m}_{#it{#mu#mu K#pi}} [GeV]", Form("Candidates/(%4.3f GeV)", h->GetBinWidth(1)), 0.05, 1.1, 1.9);
  if (0 == fData.size()) return;
  fSummary.clear();

  fpIF->fName = "fit";
  TF1 *f1 = fpIF->expogauss2c(h, 5.37, 0.04, 1.3);
  fpIF->fName = "comp";
  TF1 *fg  = fpIF->gauss2c(0., 100.);
  fg->SetLineColor(kBlue+1);
  TF1 *fe = fpIF->expo(0., 100.);
  fe->SetLineColor(kRed+2);
  fe->SetLineStyle(kSolid);

  TCanvas *c0(0);
  c0 = (TCanvas*)gROOT->FindObject("fpy_c0");
  if (!c0) c0 = new TCanvas("fpy_c0","--c0--",0,0,700,700);
  c0->cd();
  c0->Clear();
  shrinkPad(0.13, 0.19);

  fpIF->fLo = 5.0;
  fpIF->fHi = 6.0;
  cout << "==> fitPsYield::fit0_Bd2JpsiKstar> FITTING " << h->GetName() << " with limitpars = " << limitpars << endl;
  double xmin(4.9), xmax(6.0), expoLo(5.15), expoHi(6.0);
  if (fVerbose) cout << "h->GetSumOfWeights() = " << h->GetSumOfWeights() << " h->GetEntries() = " << h->GetEntries() << endl;
  //  if (fVerbose) fpIF->dumpParameters(f1);
  string fitopt = "lr";
  if (0 == fVerbose) fitopt += "q";
  h->SetMinimum(0.);
  h->SetAxisRange(fpIF->fLo, fpIF->fHi);

  if (limitpars <= 0) {
    f1->SetParLimits(0, 0.,  1.e10);
    f1->SetParLimits(1, 5.2, 5.3); // this is a Bd fit, after all
    fCombMax = h->GetMaximum();
  } else {
    for (int ipar = 0; ipar < f1->GetNpar(); ++ipar) {
      f1->ReleaseParameter(ipar);
    }
    if (fCombMax < 1) {
      cout << "FIXME: fCombMax = " << fCombMax << endl;
      fCombMax = h->GetMaximum();
    }
    double scale = h->GetMaximum()/fCombMax;

    double s1lo = fPar[2] - limitpars*fParE[2];
    double s1hi = fPar[2] + limitpars*fParE[2];
    if (s1lo < 0.) s1lo = 0.001;
    double s2lo = fPar[4] - limitpars*fParE[4];
    if (s2lo < 0.) s2lo = 0.002;
    if (s2lo < s1hi) s2lo = 1.1*s1hi;
    if (fVerbose) cout << "==> limitpars = " << limitpars << endl;
    if (limitpars > 0) {
      f1->SetParameter(0, scale*fPar[0]); if (limitpars) f1->SetParLimits(0, 0.,                           scale*(fPar[0] + fParE[0]));
      f1->SetParameter(1, fPar[1]);	  if (limitpars) f1->SetParLimits(1, fPar[1] - limitpars*fParE[1], fPar[1] + limitpars*fParE[1]);
      f1->SetParameter(2, fPar[2]);	  if (limitpars) f1->SetParLimits(2, s1lo,                         s1hi);
      f1->SetParameter(3, fPar[3]);       if (limitpars) f1->SetParLimits(3, fPar[3] - limitpars*fParE[3], fPar[3] + limitpars*fParE[3]);
      f1->SetParameter(4, fPar[4]);	  if (limitpars) f1->SetParLimits(4, s2lo,                         fPar[4] + limitpars*fParE[4]);
      f1->SetParameter(5, scale*fPar[5]); if (limitpars) f1->SetParLimits(5, 0.,                           scale*(fPar[5] + fParE[5]));
      f1->SetParameter(6, fPar[6]);       if (limitpars) f1->SetParLimits(6, fPar[6] - limitpars*fParE[6], fPar[6] + limitpars*fParE[6]);
    } else if (limitpars == 0) {
      f1->SetParameter(0, scale*fPar[0]); f1->SetParLimits(0, 0., 1.e10);
      f1->SetParameter(1, fPar[1]);
      f1->SetParameter(2, fPar[2]);
      f1->SetParameter(3, fPar[3]);
      f1->SetParameter(4, fPar[4]);
      f1->SetParameter(5, scale*fPar[5]); f1->SetParLimits(5, 0., 1.e15);
      f1->SetParameter(6, fPar[6]);       if (limitpars) f1->SetParLimits(6, -1.e7, 1.e7);
    } else {
      // do nothing, boot strap!
    }
    if (fVerbose > 1) fpIF->dumpParameters(f1);
  }

  if (h->GetSumOfWeights() > 100) {
    if (h->GetSumOfWeights() > h->GetEntries()) {
      fitopt += "w";
      h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
    } else {
      h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
    }
  } else {
    h->Draw();
    if (limitpars < 0) {
      cout << "XX> you are screwed: fitting the (un)weighted combined data and less than 100 entries" << endl;
    } else {

    }
  }

  // -- store for possible later usage
  if (limitpars < 0) {
    fPar.clear();
    fParE.clear();
    for (int i = 0; i < f1->GetNpar(); ++i) {
      fPar.push_back(f1->GetParameter(i));
      fParE.push_back(f1->GetParError(i));
    }
  }


  // -- copy parameters into components and draw them
  for (int ipar = 0; ipar < fg->GetNpar(); ++ipar) {
    fg->SetParameter(ipar, f1->GetParameter(ipar));
    fg->SetParError(ipar, f1->GetParError(ipar));
  }
  for (int ipar = 0; ipar < fe->GetNpar(); ++ipar) {
    fe->SetParameter(ipar, f1->GetParameter(ipar+5));
    fe->SetParError(ipar, f1->GetParError(ipar+5));
  }

  // -- other parameters of interest
  res->fChi2 = f1->GetChisquare();
  res->fNdof = f1->GetNDF();
  res->fProb = TMath::Prob(res->fChi2, res->fNdof);
  res->fResults.fSgPeak = f1->GetParameter(1);

  TH1D *hbla = new TH1D("hbla", "", 100, 4.5, 7.0);
  for (int i = 0; i < 10000; ++i) hbla->Fill(fg->GetRandom());
  res->fResults.fSgSigma = hbla->GetRMS();
  delete hbla;


  // -- double gaussian integral over 3 sigma region
  double sig   = fg->Integral(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
			      res->fResults.fSgPeak + 3.*res->fResults.fSgSigma);
  double sigE  = sig*(f1->GetParError(0)/f1->GetParameter(0));


  double bg   = f1->Integral(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
			     res->fResults.fSgPeak + 3.*res->fResults.fSgSigma)
    - sig ;
  double bgE  = fe->GetParError(0)/fe->GetParameter(0)*bg;

  // -- create 'sensible' errors
  sig  /= h->GetBinWidth(1);
  sigE /= h->GetBinWidth(1);
  bg   /= h->GetBinWidth(1);
  bgE  /= h->GetBinWidth(1);
  if (sigE > sig) {
    sigE = TMath::Sqrt(sig);
    cout << "rescaled error: " << sigE << " (sig = " << sig << ")"
	 << endl;
  } else {
    cout << "integral error from " << res->fResults.fSgPeak - 3.*res->fResults.fSgSigma
	 << " .. " << res->fResults.fSgPeak + 3.*res->fResults.fSgSigma
	 << ": " << sigE << " (sig = " << sig << ")"
	 << endl;
  }
  res->fResults.fBg  = bg;
  res->fResults.fBgE = bgE;
  res->fResults.fSg  = sig;
  res->fResults.fSgE = sigE;

  cout << "==> fitted " << h->GetTitle() << " signal = " << sig << " +/- " << sigE << endl;


  TLatex tl;
  tl.SetTextSize(0.03);
  if (-1 == res->fPs) {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d, weighted)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  } else if (0 == res->fPs) {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d, unweighted)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  } else {
    tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d)", res->fResults.fSg, res->fResults.fSgE, res->fPs));
  }

  c0->Modified();
  c0->Update();
  c0->SaveAs(Form("%s%s.pdf", pdfprefix.c_str(), h->GetName()));

  if (keepFunctions) {
    f1->SetName(Form("f1_%s", h->GetName()));   fFunctions.insert(make_pair(f1->GetName(), f1));
    fg->SetName(Form("sig_%s", h->GetName()));  fFunctions.insert(make_pair(fg->GetName(), fg));
    fe->SetName(Form("expo_%s", h->GetName())); fFunctions.insert(make_pair(fe->GetName(), fe));
  } else {
    delete f1;
    delete fg;
    delete fe;
  }
}


// ----------------------------------------------------------------------
TF1* fitPsYield::getFunction(string name) {
  for (map<string, TF1*>::iterator it = fFunctions.begin(); it != fFunctions.end(); ++it) {
    if (it->first == name) return it->second;
  }
  return 0;
}


// ----------------------------------------------------------------------
TF1* fitPsYield::listFunctions() {
  for (map<string, TF1*>::iterator it = fFunctions.begin(); it != fFunctions.end(); ++it) {
    cout << it->first << endl;
  }
}
