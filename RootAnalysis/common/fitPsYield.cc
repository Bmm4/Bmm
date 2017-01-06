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


using namespace std;

// ----------------------------------------------------------------------
fitPsYield::fitPsYield(string hname, TDirectory *pD, int verbose): fVerbose(verbose), fBaseName(hname), fpDir(pD), fpIF(new initFunc),
								   fCombined(0), fCombinedW8(0) {
  fData.clear();
  if (0 == fpDir) fpDir = gDirectory;
  TIter next(fpDir->GetListOfKeys());
  TKey *key(0);
  if (fVerbose > 0) cout << "fitPsYield: looking for histogram ->" << hname << "<- in directory " << fpDir->GetName() << endl;
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
    if (TString(key->GetName()).Contains(hname.c_str())) {
      string hname = key->GetName();
      TH2D *h2     = (TH2D*)fpDir->Get(hname.c_str());
      if (!h2) continue;
      if (fVerbose > 0) cout << "fitPsYield: found histogram " << h2->GetName()
			     << " with integral = " << h2->Integral(1, h2->GetNbinsX(), 1, h2->GetNbinsY())
			     << endl;
      fH2          = (TH2D*)h2->Clone(Form("fpy_%s", hname.c_str()));
      int chanB    = h2->GetYaxis()->FindBin(+0.1);
      int chan1    = h2->GetYaxis()->FindBin(1.1);
      int chan2    = h2->GetYaxis()->FindBin(h2->GetYaxis()->GetXmax()) - 1;
      fCombined    = h2->ProjectionX(Form("fpy_%s_comb", hname.c_str()), 2, 2);
      fCombinedW8  = h2->ProjectionX(Form("fpy_%s_combW8", hname.c_str()), 1, 1);
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
      // -- fW8Combined
      double ent   = fH2->Integral(1, fH2->GetNbinsX(), 1, 1);
      if (ent < 1) continue;
      fW8Combined             = new psd();
      fW8Combined->fPs        = -1;
      fW8Combined->fEntries   = ent;
      fW8Combined->fH1        = h2->ProjectionX(Form("fpy_%s_w8", hname.c_str()), 1, 1);
      fW8Combined->fH1->SetTitle(Form("fpy_%s_w8", hname.c_str()));
      // -- fUnW8Combined
      ent   = fH2->Integral(1, fH2->GetNbinsX(), 2, 2);
      if (ent < 1) continue;
      fUnW8Combined           = new psd();
      fUnW8Combined->fPs      = 0;
      fUnW8Combined->fEntries = ent;
      fUnW8Combined->fH1      = h2->ProjectionX(Form("fpy_%s_u8", hname.c_str()), 2, 2);
      fUnW8Combined->fH1->SetTitle(Form("fpy_%s_u8", hname.c_str()));
    }
  }

}


// ----------------------------------------------------------------------
fitPsYield::~fitPsYield() {
  // for (vector<psd*>::iterator it = fData.begin(); it != fData.end(); ++it) {
  //   delete *it;
  // }
  fData.clear();
}



// ----------------------------------------------------------------------
void fitPsYield::printSummary() {
  fSummary.clear();
  for (unsigned int i = 0; i < fData.size(); ++i) {
    cout << Form("%2d   sg = %8.1f +/- %3.2f", i, fData[i]->fResults.fSg, fData[i]->fResults.fSgE)
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
void fitPsYield::fitBu2JpsiKp(int limitpars, string pdfprefix) {
  // -- prefit weighted combination:
  fit0_Bu2JpsiKp(fW8Combined, -1, pdfprefix);
  // -- prefit unweighted combination:
  fit0_Bu2JpsiKp(fUnW8Combined, -1, pdfprefix);

  // -- fit all prescales
  for (unsigned int ips = 0; ips < fData.size(); ++ips) {
    fit0_Bu2JpsiKp(fData[ips], limitpars, pdfprefix);
  }

  fSummary.clear();
  for (unsigned int i = 0; i < fData.size(); ++i) {
    fSummary.fSg += fData[i]->fPs*fData[i]->fResults.fSg;
    fSummary.fSgE += fData[i]->fPs*fData[i]->fPs*fData[i]->fResults.fSgE*fData[i]->fResults.fSgE;
  }
  fSummary.fSgE = TMath::Sqrt(fSummary.fSgE);


}

// ----------------------------------------------------------------------
// limitpars < 0: determine starting values for fit
//           = 0: unconstrained fit, based on starting values of prior (limitpars < 0) call
//           > 0: constrain parameters within limitpars*sigma of prior (limitpars < 0) call
void fitPsYield::fit0_Bu2JpsiKp(psd *res, int limitpars, string pdfprefix) {
  TH1D *h = res->fH1;
  setTitles(h, "#it{m}_{#it{#mu #mu} K} #it{[GeV]}", Form("#it{Candidates/(%4.3f GeV)}", h->GetBinWidth(1)), 0.05, 1.1, 1.9);
  if (0 == fData.size()) return;
  fSummary.clear();

  fpIF->fName = "fit";
  TF1 *f1 = fpIF->pol1Err2gauss2c(0., 100.);
  fpIF->fName = "comp";
  TF1 *fg  = fpIF->gauss2c(0., 100.);
  fg->SetLineColor(kBlue+1);
  TF1 *fe = fpIF->err2(0., 100.);
  fe->SetLineColor(kRed+2);
  fe->SetLineStyle(kSolid);
  TF1 *fp = fpIF->pol1(0., 100.);
  fp->SetLineColor(kRed+2);
  fp->SetLineStyle(kSolid);


  TCanvas *c0(0);
  c0 = (TCanvas*)gROOT->FindObject("c0");
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);
  c0->cd();
  c0->Clear();
  shrinkPad(0.13, 0.19);

  fpIF->fLo = 4.9;
  fpIF->fHi = 5.9;

  double xmin(5.0), xmax(5.9), expoLo(5.16), expoHi(5.85);
  if (limitpars < 0) {
    fCombMax = h->GetMaximum();

    double mBp(5.28), sBp(0.015), stepBp(5.15);
    double p0, p1;
    fpIF->initPol1(p0, p1, h);
    fpIF->fLo = xmin;
    fpIF->fHi = xmax;
    double A   = 0.5*p1*(expoHi*expoHi - expoLo*expoLo) + p0*(expoHi - expoLo);
    double g0 = (h->Integral(h->FindBin(expoLo), h->FindBin(expoHi))*h->GetBinWidth(1) - A);
    double errN = 0.6*(h->GetBinContent(h->FindBin(expoLo - 0.1)) - h->GetBinContent(h->FindBin(expoLo)));
    f1->SetParameter(0, 0.8*g0); f1->SetParLimits(0, 0., h->GetMaximum());
    f1->SetParameter(1, mBp);    f1->SetParLimits(1, mBp - 1.*sBp, mBp + 1.*sBp);
    f1->SetParameter(2, sBp);    f1->SetParLimits(2, 0.010, 0.040);
    f1->SetParameter(3, 0.2);    f1->SetParLimits(3, 0.05, 0.60);
    f1->SetParameter(4, 5*sBp);  f1->SetParLimits(4, 0.050, 0.150);
    f1->SetParameter(5, p0);     //f1->SetParLimits(3, 0., 1.e10);
    f1->SetParameter(6, p1);     //f1->SetParLimits(4, -1.e10, 2.);
    f1->SetParameter(7, stepBp); f1->SetParLimits(7, stepBp - 2.*sBp, stepBp + 2.*sBp);
    f1->SetParameter(8, 2.*sBp); f1->SetParLimits(8, 2.*sBp - sBp, 2.*sBp + 1.5*sBp);
    f1->SetParameter(9, errN);   f1->SetParLimits(9, 0., h->GetMaximum());
  } else {
    for (int ipar = 0; ipar < f1->GetNpar(); ++ipar) {
      f1->ReleaseParameter(ipar);
    }
    if (fCombMax < 1) {
      cout << "FIXME: fCombMax = " << fCombMax << endl;
      fCombMax = h->GetMaximum();
    }
    double scale = h->GetMaximum()/fCombMax;

    f1->SetParameter(0, scale*fPar[0]); if (limitpars) f1->SetParLimits(0, 0., 1.e10);
    f1->SetParameter(1, fPar[1]);	if (limitpars) f1->SetParLimits(1, fPar[1] - limitpars*fParE[1], fPar[1] + limitpars*fParE[1]);
    f1->SetParameter(2, fPar[2]);	if (limitpars) f1->SetParLimits(2, fPar[2] - limitpars*fParE[2], fPar[2] + limitpars*fParE[2]);
    f1->SetParameter(3, fPar[3]);       if (limitpars) f1->SetParLimits(3, 0., 1.);
    f1->SetParameter(4, fPar[4]);	if (limitpars)  f1->SetParLimits(4, fPar[4] - limitpars*fParE[4], fPar[4] + limitpars*fParE[4]);
    f1->SetParameter(5, scale*fPar[5]); //f1->SetParLimits(5, pol0 - pol0E, pol0 + pol0E);
    f1->SetParameter(6, scale*fPar[6]); //f1->SetParLimits(6, pol1 - pol1E, pol1 + pol1E);
    f1->FixParameter(7, fPar[7]);
    f1->FixParameter(8, fPar[8]);
    double shift = (limitpars < 5?limitpars*0.1: 0.6);
    f1->SetParameter(9, scale*fPar[9]);	if (limitpars) f1->SetParLimits(9, scale*fPar[9]*(1.-shift), scale*fPar[9]*(1.+shift));
    if (fVerbose > 1) fpIF->dumpParameters(f1);
  }

  cout << "h->GetSumOfWeights() = " << h->GetSumOfWeights() << " h->GetEntries() = " << h->GetEntries() << endl;
  string fitopt = "lr";
  if (fVerbose) fitopt += "q";
  if (h->GetSumOfWeights() > h->GetEntries()) {
    fitopt += "w";
    h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
  } else {
    h->Fit(f1, fitopt.c_str(), "", xmin, xmax);
  }

  // -- copy parameters into components and draw them
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

  // -- store for possible later usage
  if (limitpars < 0) {
    fPar.clear();
    fParE.clear();
    for (int i = 0; i < f1->GetNpar(); ++i) {
      fPar.push_back(f1->GetParameter(i));
      fParE.push_back(f1->GetParError(i));
    }
    double NSIG  = fg->Integral(5.1, 5.6)/h->GetBinWidth(1);
    double NTOT  = h->GetSumOfWeights();
    fCombS2All   = NSIG/NTOT;
    fCombS2AllE  = TMath::Sqrt(1./NSIG + 1./NTOT)*fCombS2All;
  }

  // -- other parameters of interest
  res->fChi2 = f1->GetChisquare();
  res->fNdof = f1->GetNDF();
  res->fProb = TMath::Prob(res->fChi2, res->fNdof);
  res->fResults.fSgPeak = f1->GetParameter(1);
  res->fResults.fSgSigma = f1->GetParameter(2)*(1. - f1->GetParameter(3)) + f1->GetParameter(4)*f1->GetParameter(3);  // -- approximation!

  // -- double gaussian integral over 3 sigma region
  double sig   = fg->Integral(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
			    res->fResults.fSgPeak + 3.*res->fResults.fSgSigma);
  double sigE  = f1->IntegralError(res->fResults.fSgPeak - 3.*res->fResults.fSgSigma,
				   res->fResults.fSgPeak + 3.*res->fResults.fSgSigma);

  // -- create 'sensible' errors
  sig  /= h->GetBinWidth(1);
  sigE /= h->GetBinWidth(1);
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
  res->fResults.fSg = sig;
  res->fResults.fSgE = sigE;

  TLatex tl;
  tl.SetTextSize(0.03);
  tl.DrawLatexNDC(0.2, 0.91, Form("Sg: %.1f #pm %.1f (PS = %d)", res->fResults.fSg, res->fResults.fSgE, res->fPs));

  c0->Modified();
  c0->Update();
  c0->SaveAs(Form("%s%s.pdf", pdfprefix.c_str(), h->GetName()));

  delete f1;
  delete fg;
  delete fe;
  delete fp;
}

// ----------------------------------------------------------------------
void fitPsYield::fit0_Bs2JpsiPhi(int limitpars) {

}
