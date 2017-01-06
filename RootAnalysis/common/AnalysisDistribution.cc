#include "AnalysisDistribution.hh"

#include <iostream>
#include <iomanip>

#include "util.hh"
#include "PidTable.hh"
#include "initFunc.hh"

#include "TMath.h"
#include "TArrow.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TLine.h"

ClassImp(AnalysisDistribution)

using std::string;
using std::cout;
using std::endl;

// ----------------------------------------------------------------------
AnalysisDistribution::AnalysisDistribution(const char *name, const char *title, int nbins, double lo, double hi, double massLo, double massHi) {
  int NBINS(120);
  fSigLo = fSigHi = fBg1Lo = fBg1Hi = fBg2Lo = fBg2Hi = 0.0;
  fMassLo = massLo;
  fMassHi = massHi;
  fMassPeak = -1.;
  fMassSigma = -1.;

  fControlPlotsFileName = "controlPlot";

  fVerbose = 0;
  fDirectory = ".";

  fpIF = new initFunc();

  TH1::SetDefaultSumw2(kTRUE);

  string massbin[3];
  massbin[0] = "signal";
  massbin[1] = "sideband";
  massbin[2] = "all";
  int ndiv(504);
  for (int i = 0; i < 3; ++i) {
    hSi[i] = new TH1D(Form("%sSi%d", name, i), Form("%s, single cut, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hSi[i]->SetXTitle(title);
    hSi[i]->SetNdivisions(ndiv, "X");
  }
  for (int i = 0; i < 3; ++i) {
    hAo[i] = new TH1D(Form("%sAo%d", name, i), Form("%s, all other cuts, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hAo[i]->SetXTitle(title);
    hAo[i]->SetNdivisions(ndiv, "X");
  }
  for (int i = 0; i < 3; ++i) {
    hNm[i] = new TH1D(Form("%sNm%d", name, i), Form("%s, n-1 cuts, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hNm[i]->SetXTitle(title);
    hNm[i]->SetNdivisions(ndiv, "X");
  }
  for (int i = 0; i < 3; ++i) {
    hCu[i] = new TH1D(Form("%sCu%d", name, i), Form("%s, cumulative cuts, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hCu[i]->SetXTitle(title);
    hCu[i]->SetNdivisions(ndiv, "X");
  }
  for (int i = 0; i < 3; ++i) {
    hHLT[i] = new TH1D(Form("%sHLT%d", name, i), Form("%s, after HLT, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hHLT[i]->SetXTitle(title);
    hHLT[i]->SetNdivisions(ndiv, "X");
  }
  for (int i = 0; i < 3; ++i) {
    hPresel[i] = new TH1D(Form("%sPresel%d", name, i), Form("%s, after Presel, %s", name, massbin[i].c_str()), nbins, lo, hi);
    hPresel[i]->SetXTitle(title);
    hPresel[i]->SetNdivisions(ndiv, "X");
  }

  ndiv = 508;
  hMassSi    = new TH1D(Form("%sMassSi", name), Form("%sMassSi", name), NBINS, fMassLo, fMassHi);
  hMassSi->SetNdivisions(ndiv, "X");
  hMassAo    = new TH1D(Form("%sMassAo", name), Form("%sMassAo", name), NBINS, fMassLo, fMassHi);
  hMassAo->SetNdivisions(ndiv, "X");
  hMassNm    = new TH1D(Form("%sMassNm", name), Form("%sMassNm", name), NBINS, fMassLo, fMassHi);
  hMassNm->SetNdivisions(ndiv, "X");
  hMassCu    = new TH1D(Form("%sMassCu", name), Form("%sMassCu", name), NBINS, fMassLo, fMassHi);
  hMassCu->SetNdivisions(ndiv, "X");
  hMassHLT   = new TH1D(Form("%sMassHLT", name), Form("%sMassHLT", name), NBINS, fMassLo, fMassHi);
  hMassHLT->SetNdivisions(ndiv, "X");
  hMassPresel= new TH1D(Form("%sMassPresel", name), Form("%sMassPresel", name), NBINS, fMassLo, fMassHi);
  hMassPresel->SetNdivisions(ndiv, "X");

  hMassAll   = new TH1D(Form("%sMassAll", name), Form("%sMassALL", name), NBINS, fMassLo, fMassHi);
  hMassAll->SetNdivisions(ndiv, "X");
  hMassBGL   = new TH1D(Form("%sMassBGL", name), Form("%sMassBGL", name), NBINS, fMassLo, fMassHi);
  hMassBGL->SetNdivisions(ndiv, "X");
  hMassSG    = new TH1D(Form("%sMassSG", name), Form("%sMassSG", name), NBINS, fMassLo, fMassHi);
  hMassSG->SetNdivisions(ndiv, "X");
  hMassBGH   = new TH1D(Form("%sMassBGH", name), Form("%sMassGBH", name), NBINS, fMassLo, fMassHi);
  hMassBGH->SetNdivisions(ndiv, "X");

}



// ----------------------------------------------------------------------
AnalysisDistribution::AnalysisDistribution(const char *name, double SigLo, double SigHi, double Bg1Lo, double Bg1Hi, double Bg2Lo, double Bg2Hi) {

  fMassLo = 4.8;
  fMassHi = 6.0;

  fVerbose = 0;
  fDirectory = ".";

  fpIF = new initFunc();

  hMassSi    = (TH1D*)gDirectory->Get(Form("%sMassSi", name));
  hMassAo    = (TH1D*)gDirectory->Get(Form("%sMassAo", name));
  hMassNm    = (TH1D*)gDirectory->Get(Form("%sMassNm", name));
  hMassCu    = (TH1D*)gDirectory->Get(Form("%sMassCu", name));
  hMassHLT   = (TH1D*)gDirectory->Get(Form("%sMassHLT", name));
  hMassPresel= (TH1D*)gDirectory->Get(Form("%sMassPresel", name));

  hMassBGL    = (TH1D*)gDirectory->Get(Form("%sMassBGL", name));
  hMassBGH    = (TH1D*)gDirectory->Get(Form("%sMassBGH", name));
  hMassSG     = (TH1D*)gDirectory->Get(Form("%sMassSG", name));
  if (0 == hMassBGL) {
    cout << "%% Did NOT find " << Form("%sMassBGL", name) << " at " << hMassBGL << endl;
  }

  for (int i = 0; i < 3; ++i) {
    hSi[i]    = (TH1D*)gDirectory->Get(Form("%sSi%d", name, i));
    hAo[i]    = (TH1D*)gDirectory->Get(Form("%sAo%d", name, i));
    hNm[i]    = (TH1D*)gDirectory->Get(Form("%sNm%d", name, i));
    hCu[i]    = (TH1D*)gDirectory->Get(Form("%sCu%d", name, i));
    hHLT[i]   = (TH1D*)gDirectory->Get(Form("%sHLT%d", name, i));
    hPresel[i]= (TH1D*)gDirectory->Get(Form("%sPresel%d", name, i));
  }

  fSigLo = SigLo;
  fSigHi = SigHi;
  fBg1Lo = Bg1Lo;
  fBg1Hi = Bg1Hi;
  fBg2Lo = Bg2Lo;
  fBg2Hi = Bg2Hi;

}



// ----------------------------------------------------------------------
AnalysisDistribution::~AnalysisDistribution() {
//   if (fF0) delete fF0;
//   if (fF1) delete fF1;
//   if (fP1) delete fP1;
//   if (fPG1) delete fPG1;
//   if (fEG1) delete fEG1;
//   if (fEG2) delete fEG2;
  if (fpIF) delete fpIF;
}

// ----------------------------------------------------------------------
void AnalysisDistribution::setSigWindow(double lo, double hi) {
  fSigLo = lo;
  fSigHi = hi;
}

// ----------------------------------------------------------------------
void AnalysisDistribution::setBg1Window(double lo, double hi) {
  fBg1Lo = lo;
  fBg1Hi = hi;
}

// ----------------------------------------------------------------------
void AnalysisDistribution::setBg2Window(double lo, double hi) {
  fBg2Lo = lo;
  fBg2Hi = hi;
}

// ----------------------------------------------------------------------
void AnalysisDistribution::setAnalysisCuts(AnalysisCuts *p, const char *cutname) {
  fpAnaCuts = p;
  fCutName  = cutname;
  fCutIdx   = fpAnaCuts->getIndex(fCutName.c_str());
  if (fCutIdx < 0) {
    cout << "xx> AnalysisDistribution: ERROR "
	 << fCutName << " not found in AnalysisCuts"
	 << endl;
  }
  fHLTIdx   = fpAnaCuts->getIndex("fGoodHLT");
}


// ----------------------------------------------------------------------
double AnalysisDistribution::fitMass(TH1 *h1, double &error, int mode) {

  if (0 == h1) {
    cout << "no histogram provided. STOP." << endl;
    return -1.;
  }

  TF1 *f1;
  if (0 == mode) {
    // -- just count the number of events in the histogram
    double n = h1->GetSumOfWeights();
    error = TMath::Sqrt(n);
    h1->Draw();
    return n;
  } else if (1 == mode) {
    // -- blinded signal box
    f1 = fpIF->pol1(h1);
    h1->Fit(f1, "q", "", fMassLo, fMassHi);
    double p0  = f1->GetParameter(0);
    double p0E2= f1->GetParError(0)*f1->GetParError(0);
    double p1  = f1->GetParameter(1);
    double p1E2= f1->GetParError(1)*f1->GetParError(1);
    double dx  = fSigHi - fSigLo;
    double d2x = fSigHi*fSigHi - fSigLo*fSigLo;
    double yield = f1->Integral(fSigLo, fSigHi)/h1->GetBinWidth(1);
    double yieldAnalytical = (p0*dx + 0.5*d2x*p1)/h1->GetBinWidth(1);
    error = TMath::Sqrt((dx*dx*p0E2 + 0.25*d2x*d2x*p1E2)/h1->GetBinWidth(1));
    cout << "yield from TF1 integral: " << yield << " vs. " << yieldAnalytical << " +/- " << error << endl;
    return yield;
  } else if (10 == mode) {
    // -- One Gaussian plus pol1
    fpIF->fLo = fMassLo;
    fpIF->fHi = fMassHi;
    double peak = (fMassPeak>0.?fMassPeak:5.3);
    double sigma = (fMassSigma>0.?fMassSigma:0.04);
    f1 = fpIF->pol1gauss(h1, peak, sigma);
    h1->Fit(f1, "q", "", fMassLo, fMassHi);
    error = f1->GetParError(0)/h1->GetBinWidth(1);
    return f1->GetParameter(0)/h1->GetBinWidth(1);
  } else if (11 == mode) {
    // -- Double Gaussian with pol1
    fpIF->fLo = fMassLo;
    fpIF->fHi = fMassHi;
    double peak = (fMassPeak>0.?fMassPeak:5.3);
    double sigma = (fMassSigma>0.?fMassSigma:0.04);
    f1 = fpIF->pol1gauss2c(h1, peak, sigma);
    TFitResultPtr r;
    r = h1->Fit(f1, "lsq", "", fMassLo, fMassHi);
    double hlimit(0.), llimit(0.);
    f1->GetParLimits(3, llimit, hlimit);
    double relError = f1->GetParError(3)/f1->GetParameter(3);
    double limit = TMath::Abs(f1->GetParameter(3) - llimit);
    if (TMath::Abs(f1->GetParameter(3) - hlimit) < limit) limit = TMath::Abs(f1->GetParameter(3) - hlimit);
    if (limit < relError || relError > 0.1) {
      cout << "%%%%% REFITTING %%%%% " << endl;
      f1->SetParameters(f1->GetParameters());
      f1->SetParameter(0, h1->GetMaximum());
      f1->FixParameter(3, 0.);
      f1->FixParameter(4, 0.);
      r = h1->Fit(f1, "lsq", "", fMassLo, fMassHi);
    }
    f1->SetParameter(5, 0.);
    f1->SetParameter(6, 0.);
    peak   = f1->GetParameter(1);
    sigma  = TMath::Sqrt(f1->CentralMoment(2, fMassLo, fMassHi));
    double yield  = f1->Integral(peak-3.*sigma, peak+3.*sigma)/h1->GetBinWidth(1);
    cout << "integral:       " << yield << endl;
    double ierror = f1->IntegralError(peak-3.*sigma, peak+3.*sigma, r->GetParams(), r->GetCovarianceMatrix().GetMatrixArray())/h1->GetBinWidth(1);
    cout << "integral error: " << ierror << endl;
    double yieldE = f1->GetParError(0)/f1->GetParameter(0);
    error = yieldE*yield;
    cout << "yield error: " << error << endl;
    error = ierror;
    return yield;
  } else if (12 == mode) {
    // -- One Gaussian plus expo
    fpIF->fLo = fMassLo;
    fpIF->fHi = fMassHi;
    double peak = (fMassPeak>0.?fMassPeak:5.3);
    double sigma = (fMassSigma>0.?fMassSigma:0.04);
    f1 = fpIF->expoGauss(h1, peak, sigma);
    h1->Fit(f1, "q", "", fMassLo, fMassHi);
    error = f1->GetParError(0)/h1->GetBinWidth(1);
    return f1->GetParameter(0)/h1->GetBinWidth(1);
  } else {
    return -2.;
  }

  return -3.;
}


// ----------------------------------------------------------------------
TH1D* AnalysisDistribution::sbsDistributionPhiKK(const char *variable, const char *cut) {
  TCanvas *c0(0);
  if (fVerbose > 0) {
    gStyle->SetOptTitle(1);
    c0 = (TCanvas*)gROOT->FindObject("c1");

    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    c0->Divide(2,2);
  }

  // -- this is not really necessary, could use the class members instead
  TH1D *hm = (TH1D*)gDirectory->Get(Form("%sMass%s", variable, cut));
  TH1D *h0 = (TH1D*)gDirectory->Get(Form("%s%s0", variable, cut));
  TH1D *h1 = (TH1D*)gDirectory->Get(Form("%s%s1", variable, cut));
  if (0 == hm) {
    cout << "no histogram " << Form("%sMass%s", variable, cut) << " found in gDirectory = "; gDirectory->pwd();
    return 0;
  }

  if (fVerbose > 0) {
    c0->cd(1);
    h0->Draw("hist");
    c0->cd(2);
    h1->Draw("hist");
    c0->cd(3);
  }

  double l0   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  double l1   = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
  double s0   = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
  double s1   = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
  double u0   = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
  double u1   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);

  // -- fit mass distribution
  if (fMassLo > fMassHi) {
    fMassLo   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
    fMassHi   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);

    fMassPeak = 0.5*(s0+s1);
    fMassSigma= 0.2*(s1-s0);
  }
  fpIF->fLo = fMassLo;
  fpIF->fHi = fMassHi;

  cout << "fMass: " << fMassLo << " .. " << fMassHi << ", fMassPeak = " << fMassPeak << ", fMassSigma = " << fMassSigma << endl;

  TF1 *f1 = fpIF->phiKK(hm);
  hm->SetMinimum(0.);
  f1->SetParameter(0, 0.4*hm->GetMaximum());
  f1->SetParLimits(1, 1.017, 1.021);
  f1->SetParLimits(2, 0.001, 0.005);
  f1->SetParLimits(3, 0.05,  0.10);
  f1->SetParLimits(4, 0.006, 0.010);

  TFitResultPtr r;
  r = hm->Fit(f1, "lsq", "", fMassLo, fMassHi);
  if (fVerbose > 0) {
    hm->DrawCopy();
    cout << "Drawing histogram with name ->" << hm->GetName() << "<- and title ->" << hm->GetTitle() << "<-" << endl;
    c0->Update();
    TPaveStats *stats2 = (TPaveStats*)gPad->GetPrimitive("stats");
    if (stats2) {
      stats2->SetName("fitstatsHm");
      stats2->SetY1NDC(0.15);
      stats2->SetY2NDC(0.60);
      stats2->SetX1NDC(0.50);
      stats2->SetX2NDC(0.95);
    }
    hMassBGL->SetLineColor(kRed);
    hMassBGL->SetLineStyle(kSolid);
    hMassBGL->Draw("samehist");
    hMassBGH->SetLineColor(kRed);
    hMassBGH->SetLineStyle(kSolid);
    hMassBGH->Draw("samehist");
    hMassSG->SetLineColor(kBlue);
    hMassSG->SetLineStyle(kSolid);
    hMassSG->Draw("samehist");

  }

  // -- compute integrals
  TF1 *fs = fpIF->gauss2c(hm->GetBinLowEdge(1), hm->GetBinLowEdge(hm->GetNbinsX()+1));
  fs->SetLineColor(kBlue);
  fs->SetLineStyle(kDashed);
  fs->SetParameter(0, f1->GetParameter(0));
  fs->SetParameter(1, f1->GetParameter(1));
  fs->SetParameter(2, f1->GetParameter(2));
  fs->SetParameter(3, f1->GetParameter(3));
  fs->SetParameter(4, f1->GetParameter(4));
  fs->Draw("same");
  fpIF->dumpParameters(fs);



  // -- compute integrals
  TF1 *f2 = fpIF->argus(hm->GetBinLowEdge(1), hm->GetBinLowEdge(hm->GetNbinsX()+1));
  f2->SetLineColor(kRed);
  f2->SetLineStyle(kDashed);
  f2->SetParameter(0, f1->GetParameter(5));
  f2->SetParameter(1, f1->GetParameter(6));
  f2->SetParameter(2, f1->GetParameter(7));
  f2->Draw("same");
  fpIF->dumpParameters(f2);

  double bgl = f2->Integral(l0, l1);
  double sg  = f2->Integral(s0, s1);
  double bgh = f2->Integral(u0, u1);

  if (fVerbose > 0) {
    cout << "bgl (" << l0 << " .. " << l1 << ") = " << bgl << endl;
    cout << "sg  (" << s0 << " .. " << s1 << ") = " << sg << endl;
    cout << "bgh (" << u0 << " .. " << u1 << ") = " << bgh << endl;
  }

  TH1D *h = (TH1D*)h0->Clone(Form("sbs_%s%s", variable, cut));
  h->Sumw2();
  h->Add(h0, h1, 1., -sg/(bgl+bgh));

  if (fVerbose > 0) {
    c0->cd(4);
    h->Draw();
    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
    return h;

    c0->Clear();
    hm->SetTitle("");
    setTitles(hm, "m #it{[GeV]}", "#it{Entries/bin}", 0.05, 1.1, 1.6);
    hm->Draw();
    cout << "=========> "
	 << Form("%s/plotMass-%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/plotMass-%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  } else {
    gStyle->SetOptTitle(1);
    c0 = (TCanvas*)gROOT->FindObject("c1");
    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    hm->SetXTitle("mass [GeV]");
    hm->SetYTitle("candidates/bin");
    hm->Draw();
    TArrow aa;
    double x0 = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
    double y0 = 1.1*hm->GetBinContent(hMassBGL->FindFirstBinAbove(1.));
    double x1 = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
    double y1 = 1.1*hm->GetBinContent(hMassBGL->FindLastBinAbove(1.)+1);
    double x2 = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
    double y2 = 1.1*hm->GetBinContent(hMassBGH->FindFirstBinAbove(1.));
    double x3 = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
    double y3 = 1.1*hm->GetBinContent(hMassBGH->FindLastBinAbove(1.)+1);

    double x4 = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
    double y4 = 1.1*hm->GetBinContent(hMassSG->FindFirstBinAbove(1.));
    double x5 = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
    double y5 = 1.1*hm->GetBinContent(hMassSG->FindLastBinAbove(1.)+1);

    aa.SetLineColor(kRed);
    aa.DrawArrow(x0, y0, x0, 0.);
    aa.DrawArrow(x1, y1, x1, 0.);
    aa.SetLineColor(kBlue);
    aa.DrawArrow(x2, y2, x2, 0.);
    aa.DrawArrow(x3, y3, x3, 0.);

    aa.SetLineColor(kBlack);
    aa.DrawArrow(x4, y4, x4, 0.);
    aa.DrawArrow(x5, y5, x5, 0.);

    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  }

  return h;

}

// ----------------------------------------------------------------------
TH1D* AnalysisDistribution::sbsDistributionExpoErrGauss(const char *variable, const char *cut, double preco) {

  //  cout << "fVerbose: " << fVerbose << endl;

  TCanvas *c0(0);
  if (fVerbose > 0) {
    gStyle->SetOptTitle(1);
    c0 = (TCanvas*)gROOT->FindObject("c1");
    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    c0->Divide(2,2);
  }

  // -- this is not really necessary, could use the class members instead
  TH1D *hm = (TH1D*)gDirectory->Get(Form("%sMass%s", variable, cut));
  TH1D *h0 = (TH1D*)gDirectory->Get(Form("%s%s0", variable, cut));
  TH1D *h1 = (TH1D*)gDirectory->Get(Form("%s%s1", variable, cut));
  if (0 == hm) {
    cout << "no histogram " << Form("%sMass%s", variable, cut) << " found in gDirectory = "; gDirectory->pwd();
    return 0;
  }

  if (fVerbose > 0) {
    c0->cd(1);
    h0->Draw();
    c0->cd(2);
    h1->Draw();
  }

  double l0   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  double l1   = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
  double s0   = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
  double s1   = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
  double u0   = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
  double u1   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);

  // -- fit mass distribution
  fMassPeak = 0.5*(s0+s1);
  fMassSigma= 0.2*(s1-s0);
  fMassLo   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  fMassHi   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
  fpIF->fLo = fMassLo;
  fpIF->fHi = fMassHi;

  //  cout << "fMass: " << fMassLo << " .. " << fMassHi << ", fMassPeak = " << fMassPeak << ", fMassSigma = " << fMassSigma << endl;

  double peak  = (fMassPeak>0.?fMassPeak:5.3);
  double sigma = (fMassSigma>0.?fMassSigma:0.04);
  TF1 *f1 = fpIF->expoErrGauss(hm, peak, sigma, preco);
  hm->SetMinimum(0.);

  cout << "====> PRECO: " << preco << " in the function: " << f1->GetParameter(5) << endl;

  TFitResultPtr r;
  r = hm->Fit(f1, "lsq", "", fMassLo, fMassHi);
  if (fVerbose > 0) {
    c0->cd(3);
    hm->DrawCopy();
    //     hMassBGL->SetMinimum(0.);
    //     hMassBGL->SetMaximum(hm->GetMaximum());
    hMassBGL->Draw("same");
    hMassBGH->Draw("same");
    hMassSG->Draw("same");
  }

  TF1 *fpol1 = (TF1*)gROOT->FindObject("fpol1");
  if (fpol1) delete fpol1;
  fpol1 = fpIF->expoErr(fMassLo, fMassHi);
  fpol1->SetParameters(f1->GetParameters());

  double bgl = fpol1->Integral(l0, l1);
  double sg  = fpol1->Integral(s0, s1);
  double bgh = fpol1->Integral(u0, u1);

  if (fVerbose > 0) {
    cout << "bgl (" << l0 << " .. " << l1 << ") = " << bgl << endl;
    cout << "sg  (" << s0 << " .. " << s1 << ") = " << sg << endl;
    cout << "bgh (" << u0 << " .. " << u1 << ") = " << bgh << endl;
  }

  TH1D *h = (TH1D*)h0->Clone(Form("sbs_%s%s", variable, cut));
  h->Sumw2();
  h->Add(h0, h1, 1., -sg/(bgl+bgh));

  if (fVerbose > 0) {
    c0->cd(4);
    h->Draw();
    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  }

  if (0) {
    c0 = (TCanvas*)gROOT->FindObject("c1");
    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    hm->SetXTitle("mass [GeV]");
    hm->SetYTitle("candidates/bin");
    hm->Draw();
    TArrow aa;
    double x0 = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
    double y0 = 1.1*hm->GetBinContent(hMassBGL->FindFirstBinAbove(1.));
    double x1 = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
    double y1 = 1.1*hm->GetBinContent(hMassBGL->FindLastBinAbove(1.)+1);
    double x2 = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
    double y2 = 1.1*hm->GetBinContent(hMassBGH->FindFirstBinAbove(1.));
    double x3 = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
    double y3 = 1.1*hm->GetBinContent(hMassBGH->FindLastBinAbove(1.)+1);

    double x4 = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
    double y4 = 1.1*hm->GetBinContent(hMassSG->FindFirstBinAbove(1.));
    double x5 = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
    double y5 = 1.1*hm->GetBinContent(hMassSG->FindLastBinAbove(1.)+1);

    aa.SetLineColor(kRed);
    aa.DrawArrow(x0, y0, x0, 0.);
    aa.DrawArrow(x1, y1, x1, 0.);
    aa.SetLineColor(kBlue);
    aa.DrawArrow(x2, y2, x2, 0.);
    aa.DrawArrow(x3, y3, x3, 0.);

    aa.SetLineColor(kBlack);
    aa.DrawArrow(x4, y4, x4, 0.);
    aa.DrawArrow(x5, y5, x5, 0.);

    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  }

  return h;
}


// ----------------------------------------------------------------------
TH1D* AnalysisDistribution::sbsDistributionPol1ErrGauss(const char *variable, const char *cut, double preco) {

  cout << "fVerbose: " << fVerbose << endl;

  TCanvas *c0(0);
  if (fVerbose > 0) {
    gStyle->SetOptTitle(1);
    gStyle->SetOptFit(11111);
    c0 = (TCanvas*)gROOT->FindObject("c1");
    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    c0->Divide(2,2);
  }

  // -- this is not really necessary, could use the class members instead
  TH1D *hm = (TH1D*)gDirectory->Get(Form("%sMass%s", variable, cut));
  TH1D *h0 = (TH1D*)gDirectory->Get(Form("%s%s0", variable, cut));
  TH1D *h1 = (TH1D*)gDirectory->Get(Form("%s%s1", variable, cut));
  if (0 == hm) {
    cout << "no histogram " << Form("%sMass%s", variable, cut) << " found in gDirectory = "; gDirectory->pwd();
    return 0;
  }

  if (fVerbose > 0) {
    c0->cd(1);
    h0->Draw();
    c0->cd(2);
    h1->Draw();
    c0->cd(3);
  }

  double l0   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  double l1   = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
  double s0   = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
  double s1   = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
  double u0   = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
  double u1   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);

  // -- fit mass distribution
  fMassPeak = 0.5*(s0+s1);
  fMassSigma= 0.2*(s1-s0);
  fMassLo   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  fMassHi   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
  fpIF->fLo = fMassLo;
  fpIF->fHi = fMassHi;

  //  cout << "fMass: " << fMassLo << " .. " << fMassHi << ", fMassPeak = " << fMassPeak << ", fMassSigma = " << fMassSigma << endl;

  double peak  = (fMassPeak>0.?fMassPeak:5.3);
  double sigma = (fMassSigma>0.?fMassSigma:0.04);
  TF1 *f1 = fpIF->pol1ErrGauss(hm, peak, sigma, preco);
  hm->SetMinimum(0.);

  cout << "====> PRECO: " << preco << " in the function: " << f1->GetParameter(5) << endl;

  TFitResultPtr r;
  r = hm->Fit(f1, "ls", "", fMassLo, fMassHi);
  if (fVerbose > 0) {
    hm->DrawCopy();
//     hMassBGL->SetMinimum(0.);
//     hMassBGL->SetMaximum(hm->GetMaximum());
    hMassBGL->Draw("same");
    hMassBGH->Draw("same");
    hMassSG->Draw("same");
  }

  cout << " and after the fit: " << f1->GetParameter(5)
       << " .. " << f1->GetParameter(6)
       << " .. " << f1->GetParameter(7)
       << endl;

  TF1 *fpol1 = (TF1*)gROOT->FindObject("fpol1");
  if (fpol1) delete fpol1;
  fpol1 = fpIF->pol1Err(fMassLo, fMassHi);
  fpol1->SetParameters(f1->GetParameters());

  double bgl = fpol1->Integral(l0, l1);
  double sg  = fpol1->Integral(s0, s1);
  double bgh = fpol1->Integral(u0, u1);

  if (fVerbose > 0) {
    cout << "bgl (" << l0 << " .. " << l1 << ") = " << bgl << endl;
    cout << "sg  (" << s0 << " .. " << s1 << ") = " << sg << endl;
    cout << "bgh (" << u0 << " .. " << u1 << ") = " << bgh << endl;
  }

  TH1D *h = (TH1D*)h0->Clone(Form("sbs_%s%s", variable, cut));
  h->Sumw2();
  h->Add(h0, h1, 1., -sg/(bgl+bgh));

  if (fVerbose > 0) {
    c0->cd(4);
    h->Draw();
    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  }

  if (0) {
    c0 = (TCanvas*)gROOT->FindObject("c1");
    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    hm->SetXTitle("mass [GeV]");
    hm->SetYTitle("candidates/bin");
    hm->Draw();
    TArrow aa;
    double x0 = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
    double y0 = 1.1*hm->GetBinContent(hMassBGL->FindFirstBinAbove(1.));
    double x1 = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
    double y1 = 1.1*hm->GetBinContent(hMassBGL->FindLastBinAbove(1.)+1);
    double x2 = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
    double y2 = 1.1*hm->GetBinContent(hMassBGH->FindFirstBinAbove(1.));
    double x3 = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
    double y3 = 1.1*hm->GetBinContent(hMassBGH->FindLastBinAbove(1.)+1);

    double x4 = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
    double y4 = 1.1*hm->GetBinContent(hMassSG->FindFirstBinAbove(1.));
    double x5 = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
    double y5 = 1.1*hm->GetBinContent(hMassSG->FindLastBinAbove(1.)+1);

    aa.SetLineColor(kRed);
    aa.DrawArrow(x0, y0, x0, 0.);
    aa.DrawArrow(x1, y1, x1, 0.);
    aa.SetLineColor(kBlue);
    aa.DrawArrow(x2, y2, x2, 0.);
    aa.DrawArrow(x3, y3, x3, 0.);

    aa.SetLineColor(kBlack);
    aa.DrawArrow(x4, y4, x4, 0.);
    aa.DrawArrow(x5, y5, x5, 0.);


    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  }

  return h;
}


// ----------------------------------------------------------------------
TH1D* AnalysisDistribution::sbsDistributionExpoGauss(const char *variable, const char *cut) {

  TCanvas *c0(0);
  if (fVerbose > 0) {
    gStyle->SetOptTitle(1);
    c0 = (TCanvas*)gROOT->FindObject("c1");
    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    c0->Divide(2,2);
  }

  // -- this is not really necessary, could use the class members instead
  TH1D *hm = (TH1D*)gDirectory->Get(Form("%sMass%s", variable, cut));
  TH1D *h0 = (TH1D*)gDirectory->Get(Form("%s%s0", variable, cut));
  TH1D *h1 = (TH1D*)gDirectory->Get(Form("%s%s1", variable, cut));
  if (0 == hm) {
    cout << "no histogram " << Form("%sMass%s", variable, cut) << " found in gDirectory = "; gDirectory->pwd();
    return 0;
  }

  if (0 == h0) {
    cout << "no histogram " << Form("%sMass%s0", variable, cut) << " found in gDirectory = "; gDirectory->pwd();
    return 0;
  }

  if (0 == h1) {
    cout << "no histogram " << Form("%sMass%s1", variable, cut) << " found in gDirectory = "; gDirectory->pwd();
    return 0;
  }

  if (fVerbose > 0) {
    c0->cd(1);
    h0->SetMaximum(1.2*(h0->GetMaximum()>h1->GetMaximum()?h0->GetMaximum():h1->GetMaximum()));
    h0->Draw();
    h1->SetLineColor(kRed);
    h1->Draw("same");
  }

  double l0   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  double l1   = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
  double s0   = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
  double s1   = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
  double u0   = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
  double u1   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);

  // -- fit mass distribution
  fMassPeak = 0.5*(s0+s1);
  fMassSigma= 0.2*(s1-s0);
  fMassLo   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  fMassHi   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
  fpIF->fLo = fMassLo;
  fpIF->fHi = fMassHi;

  cout << "fMass: " << fMassLo << " .. " << fMassHi << ", fMassPeak = " << fMassPeak << ", fMassSigma = " << fMassSigma << endl;

  double peak  = (fMassPeak>0.?fMassPeak:5.3);
  double sigma = (fMassSigma>0.?fMassSigma:0.04);

  cout << " peak = " << peak << " sigma = " << sigma << endl;
  fpIF->resetLimits();
  TF1 *f1 = fpIF->expoGauss(hm, peak, sigma);
  hm->SetMinimum(0.);

  if (fVerbose > 0) {
    c0->cd(2);
    hm->Draw();
    c0->cd(3);
    f1->DrawCopy("");
  }

  TFitResultPtr r;
  string fitstring = (fVerbose>0?"ls":"lsq");
  hm->Fit(f1, fitstring.c_str(), "", fMassLo, fMassHi);
  if (fVerbose > 0) {
    hm->DrawCopy();
    hMassBGL->Draw("same");
    hMassBGH->Draw("same");
    hMassSG->Draw("same");
  }

  TF1 *fpol1 = (TF1*)gROOT->FindObject("fpol1");
  if (fpol1) delete fpol1;
  fpol1 = fpIF->expo(fMassLo, fMassHi);
  fpol1->SetParameters(f1->GetParameter(3), f1->GetParameter(4));

  double bgl = fpol1->Integral(l0, l1);
  double sg  = fpol1->Integral(s0, s1);
  double bgh = fpol1->Integral(u0, u1);

  if (fVerbose > 0) {
    cout << "bgl (" << l0 << " .. " << l1 << ") = " << bgl << endl;
    cout << "sg  (" << s0 << " .. " << s1 << ") = " << sg << endl;
    cout << "bgh (" << u0 << " .. " << u1 << ") = " << bgh << endl;
  }

  TH1D *h = (TH1D*)h0->Clone(Form("sbs_%s%s", variable, cut));
  h->Sumw2();
  h->Add(h0, h1, 1., -sg/(bgl+bgh));

  if (fVerbose > 0) {
    c0->cd(4);
    h->Draw();
    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  }

  if (0) {
    gStyle->SetOptTitle(1);
    c0 = (TCanvas*)gROOT->FindObject("c1");
    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    hm->SetXTitle("mass [GeV]");
    hm->SetYTitle("candidates/bin");
    hm->Draw();
    TArrow aa;
    double x0 = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
    double y0 = 1.1*hm->GetBinContent(hMassBGL->FindFirstBinAbove(1.));
    double x1 = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
    double y1 = 1.1*hm->GetBinContent(hMassBGL->FindLastBinAbove(1.)+1);
    double x2 = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
    double y2 = 1.1*hm->GetBinContent(hMassBGH->FindFirstBinAbove(1.));
    double x3 = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
    double y3 = 1.1*hm->GetBinContent(hMassBGH->FindLastBinAbove(1.)+1);

    double x4 = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
    double y4 = 1.1*hm->GetBinContent(hMassSG->FindFirstBinAbove(1.));
    double x5 = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
    double y5 = 1.1*hm->GetBinContent(hMassSG->FindLastBinAbove(1.)+1);

    aa.SetLineColor(kRed);
    aa.DrawArrow(x0, y0, x0, 0.);
    aa.DrawArrow(x1, y1, x1, 0.);
    aa.SetLineColor(kBlue);
    aa.DrawArrow(x2, y2, x2, 0.);
    aa.DrawArrow(x3, y3, x3, 0.);

    aa.SetLineColor(kBlack);
    aa.DrawArrow(x4, y4, x4, 0.);
    aa.DrawArrow(x5, y5, x5, 0.);

    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  }

  return h;
}


// ----------------------------------------------------------------------
TH1D* AnalysisDistribution::sbsDistribution(const char *variable, const char *cut) {
  TCanvas *c0(0);
  if (fVerbose > 0) {
    gStyle->SetOptTitle(1);
    c0 = (TCanvas*)gROOT->FindObject("c1");

    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    c0->Divide(2,2);
  }

  // -- this is not really necessary, could use the class members instead
  TH1D *hm = (TH1D*)gDirectory->Get(Form("%sMass%s", variable, cut));
  TH1D *h0 = (TH1D*)gDirectory->Get(Form("%s%s0", variable, cut));
  TH1D *h1 = (TH1D*)gDirectory->Get(Form("%s%s1", variable, cut));
  if (0 == hm) {
    cout << "no histogram " << Form("%sMass%s", variable, cut) << " found in gDirectory = "; gDirectory->pwd();
    return 0;
  }

  if (fVerbose > 0) {
    c0->cd(1);
    if (h1->GetMaximum() > h0->GetMaximum()) {
      h0->SetMaximum(1.2*h1->GetMaximum());
    }
    h0->SetLineColor(kBlue);
    h0->DrawCopy("hist");
    h1->SetLineColor(kRed);
    h1->DrawCopy("histsame");
    c0->cd(3);
  }

  double l0   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
  double l1   = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
  double s0   = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
  double s1   = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
  double u0   = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
  double u1   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);

  // -- fit mass distribution
  if (fMassLo > fMassHi) {
    fMassLo   = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
    fMassHi   = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);

    fMassPeak = 0.5*(s0+s1);
    fMassSigma= 0.2*(s1-s0);
  }
  fpIF->fLo = fMassLo;
  fpIF->fHi = fMassHi;
  cout << "fMass: " << fMassLo << " .. " << fMassHi << ", fMassPeak = " << fMassPeak << ", fMassSigma = " << fMassSigma << endl;
  double peak  = (fMassPeak>0.?fMassPeak:5.3);
  double sigma = (fMassSigma>0.?fMassSigma:0.04);
  TF1 *f1 = fpIF->pol1gauss(hm, peak, sigma);
  hm->SetMinimum(0.);
  fpIF->dumpParameters(f1);

  TFitResultPtr r;
  if (fVerbose > 0) cout << "now fitting" << endl;
  r = hm->Fit(f1, "ls", "", fMassLo, fMassHi);
  if (fVerbose > 0) {
    hm->DrawCopy();
//     hMassBGL->SetMinimum(0.);
//     hMassBGL->SetMaximum(hm->GetMaximum());
    hMassBGL->SetLineColor(kRed);
    hMassBGL->SetLineStyle(kSolid);
    hMassBGL->Draw("samehist");
    hMassBGH->SetLineColor(kRed);
    hMassBGH->SetLineStyle(kSolid);
    hMassBGH->Draw("samehist");
    hMassSG->SetLineColor(kBlue);
    hMassSG->SetLineStyle(kSolid);
    hMassSG->Draw("samehist");

  }

  // -- compute integrals
  TF1 *fpol1 = (TF1*)gROOT->FindObject("fpol1");
  if (fpol1) delete fpol1;
  fpol1 = fpIF->pol1(fMassLo, fMassHi);

  fpol1->SetParameters(f1->GetParameter(3), f1->GetParameter(4));

  double bgl = fpol1->Integral(l0, l1);
  double sg  = fpol1->Integral(s0, s1);
  double bgh = fpol1->Integral(u0, u1);

  if (fVerbose > 0) {
    cout << "bgl (" << l0 << " .. " << l1 << ") = " << bgl << endl;
    cout << "sg  (" << s0 << " .. " << s1 << ") = " << sg << endl;
    cout << "bgh (" << u0 << " .. " << u1 << ") = " << bgh << endl;
  }

  TH1D *h = (TH1D*)h0->Clone(Form("sbs_%s%s", variable, cut));
  h->Sumw2();
  h->Add(h0, h1, 1., -sg/(bgl+bgh));

  if (fVerbose > 0) {
    c0->cd(2);
    h1->Scale(sg/(bgl+bgh));
    if (h1->GetMaximum() > h0->GetMaximum()) {
      h0->SetMaximum(1.2*h1->GetMaximum());
    }
    h0->SetLineColor(kBlue);
    h0->DrawCopy("hist");
    h1->SetLineColor(kRed);
    h1->DrawCopy("histsame");
    c0->cd(3);
  }
  h->SetLineColor(kBlack);


  TLatex tl;
  tl.DrawLatexNDC(0.2, 0.70, Form("bgl: %5.3f", bgl));
  tl.DrawLatexNDC(0.2, 0.65, Form("sg:  %5.3f", sg));
  tl.DrawLatexNDC(0.2, 0.60, Form("bgh: %5.3f", bgh));
  tl.DrawLatexNDC(0.2, 0.55, Form("scl: %5.3f", -sg/(bgl+bgh)));

  tl.DrawLatexNDC(0.6, 0.70, Form("lo:  %5.3f-%5.3f", l0, l1));
  tl.DrawLatexNDC(0.6, 0.65, Form("sg:  %5.3f-%5.3f", s0, s1));
  tl.DrawLatexNDC(0.6, 0.60, Form("bgh: %5.3f-%5.3f", u0, u1));


  if (fVerbose > 0) {
    c0->cd(4);
    h->Draw();
    TLine pl;
    pl.SetNDC(kFALSE);
    pl.DrawLine(h->GetBinLowEdge(1), 0., h->GetBinLowEdge(h->GetNbinsX()+1), 0.);
    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
    return h;

    c0->Clear();
    hm->SetTitle("");
    setTitles(hm, "m #it{[GeV]}", "#it{Entries/bin}", 0.05, 1.1, 1.6);
    hm->Draw();
    cout << "=========> "
	 << Form("%s/plotMass-%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/plotMass-%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  } else {
    gStyle->SetOptTitle(1);
    c0 = (TCanvas*)gROOT->FindObject("c1");
    if (c0) {
      delete c0;
    }
    c0 = new TCanvas("c1");
    c0->Clear();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gStyle->SetOptTitle(0);
    hm->SetXTitle("mass #it{[GeV]}");
    hm->SetYTitle("candidates/bin");
    hm->Draw();
    TArrow aa;
    double ymax = hm->GetMaximum();
    double x0 = hMassBGL->GetBinLowEdge(hMassBGL->FindFirstBinAbove(1.));
    double y0 = 1.2*fpol1->Eval(x0);
    if (y0 < 0.2*ymax) y0 = 0.2*ymax;
    double x1 = hMassBGL->GetBinLowEdge(hMassBGL->FindLastBinAbove(1.)+1);
    double y1 = 1.2*fpol1->Eval(x1);
    if (y1 < 0.2*ymax) y1 = 0.2*ymax;

    double x2 = hMassBGH->GetBinLowEdge(hMassBGH->FindFirstBinAbove(1.));
    double y2 = 1.2*fpol1->Eval(x2);
    if (y2 < 0.2*ymax) y2 = 0.2*ymax;
    double x3 = hMassBGH->GetBinLowEdge(hMassBGH->FindLastBinAbove(1.)+1);
    double y3 = 1.2*fpol1->Eval(x3);
    if (y3 < 0.2*ymax) y3 = 0.2*ymax;

    double x4 = hMassSG->GetBinLowEdge(hMassSG->FindFirstBinAbove(1.));
    double y4 = 1.2*fpol1->Eval(x4);
    if (y4 < 0.2*ymax) y4 = 0.2*ymax;

    double x5 = hMassSG->GetBinLowEdge(hMassSG->FindLastBinAbove(1.)+1);
    double y5 = 1.2*fpol1->Eval(x5);
    if (y5 < 0.2*ymax) y5 = 0.2*ymax;

    aa.SetLineWidth(2);
    aa.SetLineColor(kRed);
    aa.DrawArrow(x0, y0, x0, 0.);
    aa.DrawArrow(x1, y1, x1, 0.);
    aa.SetLineColor(kBlack);
    aa.DrawArrow(x2, y2, x2, 0.);
    aa.DrawArrow(x3, y3, x3, 0.);

    aa.SetLineColor(kBlue);
    aa.DrawArrow(x4, y4, x4, 0.);
    aa.DrawArrow(x5, y5, x5, 0.);

    TString fname = fControlPlotsFileName;
    cout << "=========> "
	 << Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut)
	 << endl;
    c0->SaveAs(Form("%s/%s_%s-%s.pdf", fDirectory.c_str(), fname.Data(), variable, cut));
  }

  return h;
}




// ----------------------------------------------------------------------
void AnalysisDistribution::fill(double value, double mass) {
  int mBin(-1);

  // -- these histograms are just to KNOW afterwards what the mass windows were exactly
  //  cout <<   hMassAll << endl;
  hMassAll->Fill(mass);
  if ((fSigLo < mass) && (mass < fSigHi)) {
    hMassSG->Fill(mass);
    mBin = 0;
  }
  if ((fBg1Lo < mass) && (mass < fBg1Hi)) {
    hMassBGL->Fill(mass);
    mBin = 1;
  }
  if ((fBg2Lo < mass) && (mass < fBg2Hi)) {
    hMassBGH->Fill(mass);
    mBin = 1;
  }

  if (fVerbose > 0) {
    cout << "value: " << value
	 << " mass: " << mass
	 << " mBin: " << mBin
	 << " fCutIdx: " << fCutIdx
	 << " nm: " << fpAnaCuts->nMinus1CutsTrue(fCutIdx)
	 << " si: " << fpAnaCuts->singleCutTrue(fCutIdx)
	 << " presel: " << (*fpPreselCutTrue? 1 : 0)
	 << " this: " << this
	 << endl;
  }

  if (fpAnaCuts->singleCutTrue(fCutIdx)) {
    if (mBin > -1) hSi[mBin]->Fill(value);
    hSi[2]->Fill(value);
    hMassSi->Fill(mass);
  }
  if (fpAnaCuts->nMinus1CutsTrue(fCutIdx)) {
    if (mBin > -1) hNm[mBin]->Fill(value);
    hNm[2]->Fill(value);
    hMassNm->Fill(mass);
  }
  if (fpAnaCuts->cumulativeCutTrue(fCutIdx)) {
    if (mBin > -1) hCu[mBin]->Fill(value);
    hCu[2]->Fill(value);
    hMassCu->Fill(mass);
  }
  if (fpAnaCuts->allOtherCutsTrue(fCutIdx)) {
    if (mBin > -1) hAo[mBin]->Fill(value);
    hAo[2]->Fill(value);
    hMassAo->Fill(mass);
  }
  if (fpAnaCuts->singleCutTrue(fHLTIdx)) {
    if (mBin > -1) hHLT[mBin]->Fill(value);
    hHLT[2]->Fill(value);
    hMassHLT->Fill(mass);
  }

  if (true == *fpPreselCutTrue) {
    if (mBin > -1) hPresel[mBin]->Fill(value);
    hPresel[2]->Fill(value);
    hMassPresel->Fill(mass);
  }

}
