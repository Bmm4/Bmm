#include "../common/util.hh"

// plot trigger efficiency for signal and normalization for 2011 and 2012 vs a variable

TFile *f2011(0);
TFile *f2012(0);

TFile *n2011(0);
TFile *n2012(0);

// ----------------------------------------------------------------------
TH1D* teff(string var, double max = -1., string prefix = "s", string cuts = "0==chan&&gmuid&&bdt>0.25", string tcuts = "hlt1") {
  TTree *t = (TTree*)gDirectory->Get("events");

  string allcuts = cuts + " && " + tcuts;

  TH1D *h0 = new TH1D(Form("%s0", prefix.c_str()), cuts.c_str(), 40, 0., max);
  TH1D *h1 = new TH1D(Form("%s1", prefix.c_str()), allcuts.c_str(), 40, 0., max);
  TH1D *he = new TH1D(Form("%se", prefix.c_str()), allcuts.c_str(), 40, 0., max); he->Sumw2(); setTitles(he, var.c_str(), "eff");

  t->Draw(Form("%s>>%s0", var.c_str(), prefix.c_str()), cuts.c_str(), "goff");
  t->Draw(Form("%s>>%s1", var.c_str(), prefix.c_str()), allcuts.c_str(), "goff");

  he->Divide(h1, h0, 1., 1., "b");

  return he;

  zone(2,2);
  h0->Draw();
  c0->cd(2);
  h1->Draw();
  c0->cd(3);
  he->Draw();
}

// ----------------------------------------------------------------------
void setSignal() {
  f2011 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BsToMuMu-s01.root");
  f2012 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BsToMuMu-s01.root");

}

// ----------------------------------------------------------------------
void setNorm() {
  n2011 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BuToJpsiKp-s01.root");
  n2012 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BuToJpsiKp-s01.root");

}


// ----------------------------------------------------------------------
void plotVar2(string var, double max = -1., string cuts = "0==chan&&gmuid&&bdt>0.25", string tcuts = "hlt1") {
  if (0 == n2011) {
    setNorm();
  }

  n2011->cd("candAnaBu2JpsiK");
  TH1D *hn2011 = teff(var, max, "n", cuts, tcuts);
  setHist(hn2011, kGreen+2, 20);
  n2012->cd("candAnaBu2JpsiK");
  TH1D *hn2012 = teff(var, max, "n", cuts, tcuts);
  setHist(hn2012, kGreen+3, 24);


  hn2011->SetMinimum(0.);  hn2011->SetMaximum(1.);
  hn2011->Draw("");
  hn2012->Draw("same");


  TLegend *tle = new TLegend(0.15, 0.77, 0.35, 0.85);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->SetTextSize(0.03);

  tle->AddEntry(hn2011, "no 2011", "p");
  tle->AddEntry(hn2012, "no 2012", "p");
  tle->Draw();

  replaceAll(var, "(", "");
  replaceAll(var, ")", "");
  replaceAll(var, ":", "");
  c0->SaveAs(Form("qnd/teff-%s.pdf", var.c_str()));

}

// ----------------------------------------------------------------------
void plotVar(string var, double max = -1., string cuts = "0==chan&&gmuid&&bdt>0.25", string tcuts = "hlt1") {
  if (0 == f2011) {
    setSignal();
  }

  f2011->cd("candAnaMuMu");
  TH1D *h2011 = teff(var, max, "s", cuts, tcuts);
  setHist(h2011, kBlue+1, 20);
  f2012->cd("candAnaMuMu");
  TH1D *h2012 = teff(var, max, "s", cuts, tcuts);
  setHist(h2012, kBlue+3, 24);

  h2011->SetMinimum(0.);
  h2011->SetMaximum(1.);
  h2011->Draw();
  h2012->Draw("same");

  if (0 == n2011) {
    setNorm();
  }

  n2011->cd("candAnaBu2JpsiK");
  TH1D *hn2011 = teff(var, max, "n", cuts, tcuts);
  setHist(hn2011, kGreen+2, 20);
  n2012->cd("candAnaBu2JpsiK");
  TH1D *hn2012 = teff(var, max, "n", cuts, tcuts);
  setHist(hn2012, kGreen+3, 24);


  hn2011->Draw("same");
  hn2012->Draw("same");

  TLegend *tle = new TLegend(0.15, 0.70, 0.35, 0.85);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->SetTextSize(0.03);

  tle->AddEntry(h2011, "sg 2011", "p");
  tle->AddEntry(h2012, "sg 2012", "p");
  tle->AddEntry(hn2011, "no 2011", "p");
  tle->AddEntry(hn2012, "no 2012", "p");
  tle->Draw();

  replaceAll(var, "(", "");
  replaceAll(var, ")", "");
  replaceAll(var, ":", "");
  c0->SaveAs(Form("qnd/teff-%s.pdf", var.c_str()));
}

// ----------------------------------------------------------------------
void plotAll(string cuts = "0==chan&&gmuid&&bdt>0.25", string tcuts = "hlt1") {
  plotVar("m1pt", 20., cuts, tcuts);
  plotVar("m2pt", 20., cuts, tcuts);
  plotVar("pt", 30., cuts, tcuts);
  plotVar("prob", 1., cuts, tcuts);
  plotVar("alpha", 0.1, cuts, tcuts);
  plotVar("flsxy", 30., cuts, tcuts);

  // -- no chan cuts
  plotVar("TMath::Abs(m1eta)", 2.2, "gmuid&&bdt>0.25");
  plotVar("TMath::Abs(m2eta)", 2.2, "gmuid&&bdt>0.25");

}
