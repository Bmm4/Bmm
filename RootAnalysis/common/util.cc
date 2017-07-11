#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "util.hh"
#include "TMath.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"
#include "TH2.h"

#include<stdarg.h>

#include "dataset.hh"

using namespace std;

// ----------------------------------------------------------------------
void setMaximum(double scale, TH1 *h1, TH1 *h2) {
  double m(-99.), m1(-99.), m2(-99.);
  if (0 != h1) m1 = h1->GetMaximum();
  if (0 != h2) m2 = h2->GetMaximum();

  m = (m1 > m2? m1: m2);
  if (0 != h1) h1->SetMaximum(scale*m);
  if (0 != h2) h2->SetMaximum(scale*m);
}


// ----------------------------------------------------------------------
void setTitles(TH1 *h, const char *sx, const char *sy, float size,
	       float xoff, float yoff, float lsize, int font) {
  if (h == 0) {
    cout << " Histogram not defined" << endl;
  } else {
    h->SetXTitle(sx);                  h->SetYTitle(sy);
    h->SetTitleOffset(xoff, "x");      h->SetTitleOffset(yoff, "y");
    h->SetTitleSize(size, "x");        h->SetTitleSize(size, "y");
    h->SetLabelSize(lsize, "x");       h->SetLabelSize(lsize, "y");
    h->SetLabelFont(font, "x");        h->SetLabelFont(font, "y");
    h->GetXaxis()->SetTitleFont(font); h->GetYaxis()->SetTitleFont(font);
    h->SetNdivisions(508, "X");
  }
}

// ----------------------------------------------------------------------
void setTitles(TH1 *hold, TH1* hnew) {
  hnew->SetXTitle(hold->GetXaxis()->GetTitle());
  cout << "new x title: " << hnew->GetXaxis()->GetTitle() << endl;
  hnew->SetYTitle(hold->GetYaxis()->GetTitle());

  hnew->SetTitleOffset(hold->GetTitleOffset("X"), "x");
  hnew->SetTitleOffset(hold->GetTitleOffset("Y"), "y");

  hnew->SetTitleSize(hold->GetTitleSize("X"), "x");
  hnew->SetTitleSize(hold->GetTitleSize("Y"), "y");

  hnew->GetXaxis()->SetTitleFont(hold->GetTitleFont("X"));
  hnew->GetYaxis()->SetTitleFont(hold->GetTitleFont("Y"));

  hnew->SetLabelSize(hold->GetLabelSize("X"), "x");
  hnew->SetLabelSize(hold->GetLabelSize("Y"), "y");

  hnew->SetLabelFont(hold->GetLabelFont("X"), "x");
  hnew->SetLabelFont(hold->GetLabelFont("Y"), "y");

  hnew->SetNdivisions(508, "X");

}

// ----------------------------------------------------------------------
void addOverflow(TH1D *h) {
  int nbins = h->GetNbinsX();
  h->SetBinContent(nbins, h->GetBinContent(nbins) + h->GetBinContent(nbins+1));
  h->SetBinContent(nbins+1, 0.);
}


// ----------------------------------------------------------------------
void addUnderflow(TH1D *h) {
  h->SetBinContent(1, h->GetBinContent(1) + h->GetBinContent(0));
  h->SetBinContent(0, 0.);
}

// ----------------------------------------------------------------------
void setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size);
  h->SetStats(kFALSE);
  h->SetFillStyle(0); h->SetFillColor(color);
}

// ----------------------------------------------------------------------
void setHist(TH1 *h, dataset *ds) {
  if (0 == ds) {
    cout << "no valid dataset specified" << endl;
    return;
  }
  if (0 == h) {
    cout << "no valid histogram specified" << endl;
    return;
  }
  if (ds->fColor > -1) setHist(h, ds->fColor, ds->fSymbol, ds->fSize, ds->fWidth);
  if (ds->fFillStyle > -1) setFilledHist(h, ds->fColor, ds->fFcolor, ds->fFillStyle, ds->fWidth);
}

// ----------------------------------------------------------------------
void setFilledHist(TH1 *h, Int_t color, Int_t fillcolor, Int_t fillstyle, double width) {
  // Note: 3004, 3005 are crosshatches
  // ----- 1000       is solid
  //       kYellow    comes out gray on bw printers
  h->SetLineColor(color);     h->SetLineWidth(width);
  h->SetFillStyle(fillstyle); h->SetFillColor(fillcolor);
}


// ----------------------------------------------------------------------
void setGraph(TGraph *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(static_cast<Width_t>(width));
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size);
}


// ----------------------------------------------------------------------
void colors(int i) {
  TColor *a;

  if (i == 0) {
    a = new TColor(600, 0.66, 0.66, 0.66, "gray");
    a = new TColor(601, 0.90, 0.90, 0.90, "gray");
    a = new TColor(602, 0.26, 0.26, 0.26, "gray");
    a = new TColor(610, 0.00, 0.66, 0.00, "green");
    a = new TColor(612, 0.00, 0.33, 0.00, "darkgreen");
    a = new TColor(613, 0.00, 0.90, 0.00, "lightgreen");
    a = new TColor(624, 0.00, 0.00, 0.60, "blue");
    a = new TColor(625, 0.00, 0.00, 0.85, "lightblue");
    a = new TColor(626, 0.00, 0.00, 0.35, "darkblue");
    a = new TColor(637, 0.50, 0.00, 0.00, "red");
    a = new TColor(638, 0.99, 0.00, 0.00, "lightred");
    a = new TColor(639, 0.99, 0.00, 0.00, "darkred");
  }
  if (0) cout << a << endl;
}


// ----------------------------------------------------------------------
void printNonZero(TH1 *h) {
  if (0 == h) {
    cout << "histogram does not exist" << endl;
    return;
  }
  double con(0.), min(0.), max(0.);
  for (Int_t i = 0; i <= h->GetNbinsX()+1; ++i) {
    con = h->GetBinContent(i);
    if (con > 0.) {
      min = h->GetBinLowEdge(i);
      max = min + h->GetBinWidth(i);
      cout << Form("%3d ", i) << Form(" %7.3f ", min) << " .. " << Form(" %7.3f ", max) << ":"
	   << Form(" %12.3f", con) << " +/- " << Form("%12.3f ", h->GetBinError(i))
	   << h->GetXaxis()->GetBinLabel(i)
           << endl;
    }
  }
}

// ----------------------------------------------------------------------
void printAxesLabels(TH1 *h) {
  if (0 == h) {
    cout << "histogram does not exist" << endl;
    return;
  }
  double min(0.), max(0.);
  TString la;
  for (Int_t i = 0; i <= h->GetNbinsX()+1; ++i) {
    la = h->GetXaxis()->GetBinLabel(i);
    min = h->GetBinLowEdge(i);
    max = min + h->GetBinWidth(i);
    if (la != "") {
      cout << Form("%3d ", i) << Form(" %7.3f ", min) << " .. " << Form(" %7.3f ", max) << ":"
	   << h->GetXaxis()->GetBinLabel(i)
           << endl;
    }
  }
}


// ----------------------------------------------------------------------
void stampAndSave(TCanvas *fC, const char *s) {
  fC->cd();
  TLatex tl;
  TString filename(gFile->GetName()); Int_t index = filename.Index(".root"); filename.Remove(index);
  TString psname = filename + TString(s);
  double oldA = tl.GetTextAngle();
  double oldS= tl.GetTextSize();
  tl.SetTextAngle(90); tl.SetTextSize(0.02);
  tl.DrawTextNDC(0.001, 0.05, psname.Data());
  fC->SaveAs(psname.Data());
  tl.SetTextAngle(oldA);
  tl.SetTextSize(oldS);
}

// ----------------------------------------------------------------------
void shrinkPad(double b, double l, double r, double t) {
  gPad->SetBottomMargin(b);
  gPad->SetLeftMargin(l);
  gPad->SetRightMargin(r);
  gPad->SetTopMargin(t);
}

// ----------------------------------------------------------------------
void zone(int x, int y, TCanvas *c) {
  if (c == 0) {
    c = (TCanvas*)gROOT->FindObject("c0");
    if (c == 0) {
      cout << "TCanvas c0 not found. Creating my own version." << endl;
      c = new TCanvas("c0","--c0--",356,0,656,700);
    }
  }
  c->Clear();
  c->Divide(x, y);
  c->cd(1);
}


// ----------------------------------------------------------------------
int wait() {
  cout << " Continue [<RET>|q]?  ";
  char x;
  x = getchar();
  if ((x == 'q') || (x == 'Q')) return 1;
  return 0;
}


// ----------------------------------------------------------------------
double dEff(int in, int iN) {
  double n = (double)in;
  double N = (double)iN;
  return sqrt(((n+1)*(N-n+1))/((N+3)*(N+2)*(N+2)));
}

// ----------------------------------------------------------------------
double dEff(int iN, double eff) {
  double N = (double)iN;
  double n = N*eff;
  return sqrt(((n+1)*(N-n+1))/((N+3)*(N+2)*(N+2)));
}

// ----------------------------------------------------------------------
double dEff(double eff, int in) {
  double n = (double)in;
  double N = in/eff;
  return sqrt(((n+1)*(N-n+1))/((N+3)*(N+2)*(N+2)));
}

// ----------------------------------------------------------------------
double dEff(double n, double nE, double N, double NE) {

  if (n < N) {
    return sqrt((N*N - 2*n*N)*nE*nE + n*n*NE*NE)/(N*N);
  } else {
    return 0.;
  }
}

// ----------------------------------------------------------------------
double dBinomial(double n, double N) {
  double w = n/N;
  return TMath::Sqrt(TMath::Abs(w*(1-w)/N));
}

// ----------------------------------------------------------------------
double dRatio(double a, double b) {
  return TMath::Sqrt((a/(b*b)) + ((a*a)/(b*b*b)));
}

// ----------------------------------------------------------------------
double dRatio(double a, double ae, double b, double be) {
  return TMath::Sqrt(((ae*ae)/(b*b)) + ((a*a*be*be)/(b*b*b*b)));
}


// ----------------------------------------------------------------------
double dBF(double n, double nE, double N, double NE, double e, double eE) {
  double n2 = n*n;
  double e2 = e*e;
  double N2 = N*N;

  double dbdn2 = 1./(e2*N2);
  double dbde2 = (n2)/(e2*e2*N2);
  double dbdN2 = (n2)/(e2*N2*N2);

  return TMath::Sqrt(nE*nE*dbdn2 + NE*NE*dbdN2 + eE*eE*dbde2);
}

// ----------------------------------------------------------------------
double  getError(TH1* h) {
  double e(0.), t(0.);
  if (h->InheritsFrom(TH2::Class())) {
    for (int ix = 1; ix <= h->GetNbinsX(); ix++) {
      for (int iy = 1; iy <= h->GetNbinsY(); iy++) {
	t = h->GetCellError(ix, iy);
	e += t*t;
      }
    }
  } else {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      t = h->GetBinError(i);
      e += t*t;
    }
  }
  return TMath::Sqrt(e);
}

// ----------------------------------------------------------------------
double  getErrorRange(TH1* h, int min, int max) {
  double e(0.), t(0.);
  if (min < 0) min = 1;
  if (max < 0) max = h->GetNbinsX();
  for (int i = min; i <= max; ++i) {
    t = h->GetBinError(i);
    e += t*t;
  }
  return TMath::Sqrt(e);
}

// -------------------------------------------------------------------------------
void babar(double xpos, double ypos, double scale, int prel) {
  TText *Babar = new TText();
  Babar->SetNDC(kTRUE);
  Babar->SetTextFont(32);
  Babar->SetTextSize(0.10*scale);
  Babar->DrawText(xpos,ypos,"B");
  Babar->SetTextSize(0.075*scale);
  Babar->DrawText(xpos+0.042*scale,ypos,"A");
  Babar->SetTextSize(0.10*scale);
  Babar->DrawText(xpos+0.078*scale,ypos,"B");
  Babar->SetTextSize(0.075*scale);
  Babar->DrawText(xpos+0.120*scale,ypos,"AR");

  if (prel) {
    Babar->SetTextSize(0.10*0.8*scale);
    Babar->DrawText(xpos, ypos - 0.1*0.8*scale, "Preliminary");
  }
  delete Babar;
}

// ----------------------------------------------------------------------
TH1D *unmix(TH1D *rightSign, TH1D *wrongSign, double chid) {
  TH1D *prompt = new TH1D(*rightSign); prompt->SetName("prompt"); prompt->Reset();
  prompt->Add(rightSign, wrongSign, (1. - chid) / (1. - 2.*chid), -chid / (1.-2.*chid));
  return prompt;
}


// ----------------------------------------------------------------------
double chi2Prob(double chi2, double ndof) {
  return (1. - TMath::Gamma(0.5*ndof, 0.5*chi2));
}

// ----------------------------------------------------------------------
double chi2Test(TH1 *h1, TH1 *h2, double& chi2, double& ndof, int constrain) {
  int nbins = h1->GetNbinsX();
  if (nbins != h2->GetNbinsX()) {
    cout << "chi2Test: Number of bins not the same" << endl;
    return -99.;
  }
  double df = nbins - 1 - constrain;
  double chsq(0.), a1(0.), a2(0.);
  for (int i = 1; i <= nbins; ++i) {
    a1 = h1->GetBinContent(i);
    a2 = h2->GetBinContent(i);
    if ((TMath::Abs(a1) < 1.e-8) && (TMath::Abs(a2) < 1.e-8)) {
      df -= 1.;
    } else if ((a1 < 0.) || (a2 < 0.)) {
      df -= 1.;
    } else {
      cout << "Adding " << ((a1-a2)*(a1-a2))/(a1+a2)  << " from " << a1 << "  " << a2 << " for xmin =  " << h1->GetBinLowEdge(i) << endl;
      chsq += ((a1-a2)*(a1-a2))/(a1+a2);
    }
  }
  double gamma = 1. - TMath::Gamma(0.5*df, 0.5*chsq);
  chi2 = chsq;
  ndof = df;
  return gamma;
}

// ----------------------------------------------------------------------
double chi2TestErr(TH1 *h1, TH1 *h2, double& chi2, double& ndof, int constrain) {
  int nbins = h1->GetNbinsX();
  if (nbins != h2->GetNbinsX()) {
    cout << "chi2Test: Number of bins not the same" << endl;
    return -99.;
  }
  double df = nbins - 1 - constrain;
  double chsq(0.), a1(0.), a2(0.), e1(0.), e2(0.);
  for (int i = 1; i <= nbins; ++i) {
    a1 = h1->GetBinContent(i);
    e1 = h1->GetBinError(i) * h1->GetBinError(i);
    a2 = h2->GetBinContent(i);
    e2 = h2->GetBinError(i) *  h2->GetBinError(i);
    if ((TMath::Abs(a1) < 1.e-8) && (TMath::Abs(a2) < 1.e-8)) {
      df -= 1.;
    } else if ((a1 < 0.) || (a2 < 0.)) {
      df -= 1.;
    } else {
      chsq += ((a1 - a2) * (a1 - a2)) / (e1 + e2);
    }
  }
  double gamma = 1. - TMath::Gamma(0.5*df, 0.5*chsq);
  chi2 = chsq;
  ndof = df;
  return gamma;
}

// ----------------------------------------------------------------------
void average(double &av, double &error, int n, double *val, double *verr) {

  double e(0.), w8(0.), sumW8(0.), sumAve(0.);
  for (int i = 0; i < n; ++i) {
    //    cout << i << " " << val[i] << " +/- " << verr[i] << endl;

    // -- calculate mean and error
    e = verr[i];
    if (e > 0.) {
      w8 = 1./(e*e);
      sumW8  += w8;
      sumAve += w8*val[i];
    } else {
      cout << "average: Error = 0 for " << val[i] << endl;
      continue;
    }
  }
  if (sumW8 > 0.) {
    av = sumAve/sumW8;
    error = 1./TMath::Sqrt(sumW8);
  } else {
    av = -99.;
    error = -99.;
  }

}


// ----------------------------------------------------------------------
void replaceAll(string &str, const string &from, const string &to) {
  if (from.empty()) return;
  size_t start_pos = 0;
  while((start_pos = str.find(from, start_pos)) != string::npos) {
    str.replace(start_pos, from.length(), to);
    start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
  }
}

// ----------------------------------------------------------------------
bool isBeautyMesonWeak(int i) {
  return (5 == TMath::Abs(i)/100);
}

// ----------------------------------------------------------------------
bool isBeautyBaryonWeak(int i) {
  if (5122 == TMath::Abs(i)) return true;
  if (5132 == TMath::Abs(i)) return true;
  if (5232 == TMath::Abs(i)) return true;
  if (5332 == TMath::Abs(i)) return true;
  return false;
}


// ----------------------------------------------------------------------
bool isCharmMesonWeak(int i) {
  return (4 == TMath::Abs(i)/100);
}

// ----------------------------------------------------------------------
bool isCharmBaryonWeak(int i) {
  if (4122 == TMath::Abs(i)) return true;
  if (4132 == TMath::Abs(i)) return true;
  if (4232 == TMath::Abs(i)) return true;
  if (4332 == TMath::Abs(i)) return true;
  return false;
}


// ----------------------------------------------------------------------
bool isBeautyMeson(int i) {
  int rest = (TMath::Abs(i)%1000);
  return (rest > 499 && rest < 600);
}

// ----------------------------------------------------------------------
bool isBeautyBaryon(int i) {
  int rest = (TMath::Abs(i)%1000);
  return (rest > 4999 && rest < 6000);
}

// ----------------------------------------------------------------------
bool isCharmMeson(int i) {
  int rest = (TMath::Abs(i)%1000);
  return (rest > 399 && rest < 500);
}

// ----------------------------------------------------------------------
bool isCharmBaryon(int i) {
  int rest = (TMath::Abs(i)%1000);
  return (rest > 3999 && rest < 5000);
}

// ----------------------------------------------------------------------
bool isLightMeson(int i) {
  int rest = (TMath::Abs(i)%1000);
  return (rest > 99 && rest < 400);
}

// ----------------------------------------------------------------------
bool isStableCharged(int i) {
  if (211 == TMath::Abs(i)) return true;
  if (321 == TMath::Abs(i)) return true;
  if (11 == TMath::Abs(i)) return true;
  if (13 == TMath::Abs(i)) return true;
  if (2122 == TMath::Abs(i)) return true;
  return false;
}


// ----------------------------------------------------------------------
void showOverflow(TH1 *h) {
  h->SetBinContent(h->GetNbinsX(), h->GetBinContent(h->GetNbinsX()) + h->GetBinContent(h->GetNbinsX()+1));
  h->SetBinContent(h->GetNbinsX()+1, 0.);
  h->SetBinContent(1, h->GetBinContent(1) + h->GetBinContent(0));
  h->SetBinContent(0, 0.);
}


// ----------------------------------------------------------------------
vector<int> defVector(int n, ...) {
  va_list vl;
  va_start(vl, n);
  vector<int> vect;
  int a(0);
  for (int i = 0; i < n; ++i) {
    a = va_arg(vl, int);
    vect.push_back(a);
  }
  va_end(vl);
  return vect;
}


// ----------------------------------------------------------------------
string formatTex(double n, string name, int digits, int sgn) {

  char line[200];
  if ( TMath::IsNaN(n) ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str());
  } else if (0 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%i } } }", name.c_str(), static_cast<int>(n));
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.1f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.1f } } }", name.c_str(), n);
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.2f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.2f } } }", name.c_str(), n);
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.3f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.3f } } }", name.c_str(), n);
  } else if (4 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.4f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.4f } } }", name.c_str(), n);
  } else if (5 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%6.5f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+6.5f } } }", name.c_str(), n);
  } else if (6 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%7.6f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+7.6f } } }", name.c_str(), n);
  } else {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%f } } }", name.c_str(), n);
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+f } } }", name.c_str(), n);
  }

  return string(line);
}

// ----------------------------------------------------------------------
string formatTexErrSci(double n, double nE, string name, int digits, int sgn) {

  char line[200];

  if (TMath::Abs(n) < 1.e-12) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{0} }", name.c_str());
    return string(line);
  }

  bool negative(false);
  if (n < 0) {
    n *= -1.;
    negative = true;
  }
  double expoN  = TMath::Log10(n);
  //  cout << "original expoN = " << expoN << endl;
  //  if (negative) expoN -= 1;
  double expoNE = TMath::Log10(nE);
  //  cout << "original expoNE = " << expoNE << endl;
  if (expoNE < 0) expoNE -= 1;

  double ratio = expoNE - expoN;
  double mantN  = n / TMath::Power(10, static_cast<int>(expoN));
  if (negative) mantN *= -1.;
  double mantNE = nE / TMath::Power(10, static_cast<int>(expoNE-ratio));
  if (nE < 1.e-8) mantNE = 0.;

  //  cout << setw(2) << "expoN = " << expoN << " expoNE == " << expoNE << " mantN = " << mantN << " mantNE = " << mantNE << " ratio = " << ratio << endl;

  if (TMath::IsNaN(n) || TMath::IsNaN(nE)) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str());
  } else if (1 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.1f \\pm %5.1f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%+5.1f \\pm %5.1f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
  } else if (2 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.2f \\pm %5.2f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%+5.2f \\pm %5.2f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
  } else if (3 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.3f \\pm %5.3f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%+5.3f \\pm %5.3f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
  } else if (4 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.4f \\pm %5.4f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%+5.4f \\pm %5.4f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
  } else if (5 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.5f \\pm %5.5f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%+5.5f \\pm %5.5f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
  } else if (6 == digits ) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%5.6f \\pm %5.6f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{(%+5.6f \\pm %5.6f)\\times 10^{%d}} } }", name.c_str(), mantN, mantNE, static_cast<int>(expoN));
  } else if (digits <0) {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.1f \\times 10^{%d}} } }", name.c_str(), mantN, static_cast<int>(expoN));
    if (sgn) sprintf(line, "\\vdef{%s}   {\\ensuremath{{%+5.1f \\times 10^{%d}} } }", name.c_str(), mantN, static_cast<int>(expoN));
  } else {
    sprintf(line, "\\vdef{%s}   {\\ensuremath{{%5.0f \\times 10^{%d}} } }", name.c_str(), mantN, static_cast<int>(expoN));
  }

  return string(line);
}

// ----------------------------------------------------------------------
string formatTex(double n, string name, string tag) {

  char line[200];

  if (TMath::IsNaN(n) ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{\\mathrm{NaN} } } }", name.c_str(), tag.c_str());
    //   } else if ( n > 1.e10) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm2(n)).Data());
    //   } else if ( n > 1.e4) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm(n)).Data());
  } else if ( n > 100. ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%6.0f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1. ) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.2f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1.e-1) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name.c_str(), tag.c_str(), n);
  } else if ( n > 1.e-3) {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%4.3f } } }", name.c_str(), tag.c_str(), n);
    //   } else if ( n > 1.e-9 ) {
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm(n)).Data());
    //   } else if ( n > 1.e-19 ){
    //     sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{%s } } }", name.c_str(), tag.c_str(), (texForm2(n)).Data());
  } else {
    sprintf(line, "\\vdef{%s:%s}   {\\ensuremath{{0.0 } } }", name.c_str(), tag.c_str());
  }

  string result(line);
  return result;
}

// ----------------------------------------------------------------------
// -- based on the code downloaded on 2016/06/07 from the pubcomm twiki
//    (fails miserably when the top margin does not adhere to GHM's expectation)
void stamp(double x1, string text1, string text1a, double x2, string text2, int iPosX) {

  TLatex *tl = new TLatex();
  tl->SetNDC(kTRUE);

  float cmsTextFont   = 61;  // default is helvetic-bold
  float extraTextFont = 52;  // default is helvetica-italics

  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextSize     = 0.6;
  float lumiTextOffset   = 0.2;
  float cmsTextSize      = 0.75;

  // -- this seems better (at least for our working plots)
  lumiTextSize     = 0.5;
  lumiTextOffset   = 0.16;
  cmsTextSize      = 0.60;

  float relPosX    = 0.045;
  float relPosY    = 0.035;
  float relExtraDY = 1.2;

  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize  = 0.76;

  //  float H = gPad->GetWh();
  //  float W = gPad->GetWw();
  float l = gPad->GetLeftMargin();
  float t = gPad->GetTopMargin();
  float r = gPad->GetRightMargin();
  float b = gPad->GetBottomMargin();

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  int alignY_=3;
  int alignX_=2;
  if( iPosX/10==0 ) alignX_=1;
  if( iPosX==0    ) alignX_=1;
  if( iPosX==0    ) alignY_=1;
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  //if( iPosX == 0  ) relPosX = 0.12;
  int align_ = 10*alignX_ + alignY_;

  float posX_(0);
  if (iPosX%10 <= 1)  {
    posX_ = l + relPosX*(1-l-r);
  } else if (iPosX%10 == 2) {
    posX_ = l + 0.5*(1-l-r);
  } else if (iPosX%10 == 3) {
    posX_ = 1 - r - relPosX*(1-l-r);
  }
  float posY_ = 1 - t - relPosY*(1-t-b);

  // -- "CMS"
  tl->SetTextFont(cmsTextFont);
  tl->SetTextSize(cmsTextSize*t);
  tl->SetTextAlign(align_);
  tl->DrawLatex(posX_, posY_, text1.c_str());
  //  cout << "stamp() > " << posX_ << " " << posY_ << " " << text1 << endl;
  // -- text below "CMS"
  tl->SetTextFont(extraTextFont);
  tl->SetTextAlign(align_);
  tl->SetTextSize(extraTextSize*t);
  tl->DrawLatex(posX_, posY_- relExtraDY*cmsTextSize*t, text1a.c_str());
  //  cout << "stamp() > " << posX_ << " " << posY_- relExtraDY*cmsTextSize*t << " " << text1a << endl;
  // -- luminosity
  tl->SetTextFont(42);
  tl->SetTextAlign(31);
  tl->SetTextSize(lumiTextSize*t);
  tl->DrawLatex(1-r, 1-t+lumiTextOffset*t, text2.c_str());
  //  cout << "stamp() > " << 1-r << " " << 1-t+lumiTextOffset*t << " " << text2 << endl;



  // tl->SetNDC(kTRUE);
  // tl->SetTextSize(0.05);
  // tl->DrawLatex(x1, 0.91, text1.c_str());
  // tl->DrawLatex(x2, 0.91, text2.c_str());


}


// ----------------------------------------------------------------------=
void rmSubString(string &sInput, const string &sub) {
  string::size_type foundpos = sInput.find(sub);
  if (foundpos != string::npos)  sInput.erase(sInput.begin() + foundpos, sInput.begin() + foundpos + sub.length());
}

// ----------------------------------------------------------------------=
void rmPath(string &sInput) {
  while(string::size_type foundpos = sInput.find("/")) {
    if (foundpos != string::npos)  sInput.erase(sInput.begin(), sInput.begin() + foundpos);
  }
  sInput.erase(sInput.begin(), sInput.begin() + 1); // delete also leading /
}

// ----------------------------------------------------------------------
vector<string> glob(string basedir, string basename) {
  cout << "Looking in " << basedir << " for " << basename << endl;
  vector<string> lof;
  TString fname;
  const char *file;
  TSystem *lunix = gSystem; //new TUnixSystem();
  void *pDir = lunix->OpenDirectory(basedir.c_str());
  while ((file = lunix->GetDirEntry(pDir))) {
    fname = file;
    if (fname.Contains(basename.c_str())) {
      lof.push_back(string(fname));
    }
  }
  return lof;
}

// ----------------------------------------------------------------------
vector<string>& split(const string &s, char delim, vector<string> &elems) {
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}


// ----------------------------------------------------------------------
vector<string> split(const string &s, char delim) {
  vector<string> elems;
  split(s, delim, elems);
  return elems;
}

// ----------------------------------------------------------------------
void cleanupString(string &s) {
  replaceAll(s, "\t", " ");
  string::size_type s1 = s.find("#");
  if (string::npos != s1) s.erase(s1);
  if (0 == s.length()) return;
  string::iterator new_end = unique(s.begin(), s.end(), bothAreSpaces);
  s.erase(new_end, s.end());
  if (s.substr(0, 1) == string(" ")) s.erase(0, 1);
  if (s.substr(s.length()-1, 1) == string(" ")) s.erase(s.length()-1, 1);
}

// ----------------------------------------------------------------------
bool bothAreSpaces(char lhs, char rhs) {
  return (lhs == rhs) && (lhs == ' ');
}

// ----------------------------------------------------------------------
double quadraticSum(int n, ...) {
  va_list vl;
  va_start(vl, n);
  double a(0.), sum(0.);
  for (int i = 0; i < n; ++i) {
    a = va_arg(vl, double);
    sum += a*a;
  }

  va_end(vl);
  return TMath::Sqrt(sum);
}


// ----------------------------------------------------------------------
TLegend* newLegend(string title, double x1, double y1, double x2, double y2,
		   vector<TH1*> hists, vector<string> names, vector<string> options) {
  if (hists.size() != names.size()) {
    cout << "hists and names vectors do not match size" << endl;
    return 0;
  }
  if (hists.size() != options.size()) {
    cout << "hists and options vectors do not match size" << endl;
    return 0;
  }
  TLegend *legg = new TLegend(x1, y1, x2, y2, title.c_str());
  legg->SetFillStyle(0);
  legg->SetBorderSize(0);
  legg->SetTextSize(0.04);
  legg->SetFillColor(0);
  legg->SetTextFont(52);
  for (unsigned int i = 0; i < hists.size(); ++i) {
    legg->AddEntry(hists[i], names[i].c_str(), options[i].c_str());
  }
  return legg;
}
