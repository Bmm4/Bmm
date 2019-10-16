#include "../common/util.hh"

// plot trigger efficiency for signal and normalization for chan 0 and chan 1 vs a variable

TFile *f2012(0);
TFile *n2012(0);

// ----------------------------------------------------------------------
TH1D* teff(string var, double max = -1., string prefix = "s", string cuts = "gmuid&&bdt>0.34", string tcuts = "hlt1") {
  TTree *t = (TTree*)gDirectory->Get("events");

  string allcuts = cuts + " && " + tcuts;

  TH1D *h0 = new TH1D(Form("%s0", prefix.c_str()), cuts.c_str(), 40, 0., max);
  TH1D *h1 = new TH1D(Form("%s1", prefix.c_str()), allcuts.c_str(), 40, 0., max);
  TH1D *he = new TH1D(Form("%se", prefix.c_str()), allcuts.c_str(), 40, 0., max); he->Sumw2(); setTitles(he, var.c_str(), "eff");

  t->Draw(Form("%s>>%s0", var.c_str(), prefix.c_str()), cuts.c_str(), "goff");
  t->Draw(Form("%s>>%s1", var.c_str(), prefix.c_str()), allcuts.c_str(), "goff");

  he->Divide(h1, h0, 1., 1., "b");

  return he;
}

// ----------------------------------------------------------------------
void setSignal() {
  f2012 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BsToMuMu-s01.root");

}

// ----------------------------------------------------------------------
void setNorm() {
  n2012 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BuToJpsiKp-s01.root");
}

// ----------------------------------------------------------------------
void plotVar(string var, double max = -1., string cuts = "gmuid&&bdt>0.34", string tcuts = "hlt1") {
  if (0 == f2012) {
    setSignal();
  }

  f2012->cd("candAnaMuMu");
  TH1D *hchan0 = teff(var, max, "s", "chan==0&&" + cuts, tcuts);
  setHist(hchan0, kBlue+1, 20);
  f2012->cd("candAnaMuMu");
  TH1D *hchan1 = teff(var, max, "s", "chan==1&&" + cuts, tcuts);
  setHist(hchan1, kBlue+3, 24);

  hchan0->SetMinimum(0.);
  hchan0->SetMaximum(1.);
  hchan0->Draw();
  hchan1->Draw("same");

  if (0 == n2012) {
    setNorm();
  }

  n2012->cd("candAnaBu2JpsiK");
  TH1D *hnchan0 = teff(var, max, "n", "chan==0&&" + cuts, tcuts);
  setHist(hnchan0, kGreen+2, 20);
  n2012->cd("candAnaBu2JpsiK");
  TH1D *hnchan1 = teff(var, max, "n", "chan==1&&" + cuts, tcuts);
  setHist(hnchan1, kGreen+3, 24);


  hnchan0->Draw("same");
  hnchan1->Draw("same");

  TLegend *tle = new TLegend(0.20, 0.70, 0.35, 0.85);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->SetTextSize(0.03);

  tle->SetHeader("2012");
  tle->AddEntry(hchan0, "sg chan0", "p");
  tle->AddEntry(hchan1, "sg chan1", "p");
  tle->AddEntry(hnchan0, "no chan0", "p");
  tle->AddEntry(hnchan1, "no chan1", "p");
  tle->Draw();

  replaceAll(var, "(", "");
  replaceAll(var, ")", "");
  replaceAll(var, ":", "");
  c0->SaveAs(Form("qnd/teff2-%s.pdf", var.c_str()));
}

// ----------------------------------------------------------------------
void plotAll(string cuts = "gmuid&&bdt>0.34", string tcuts = "hlt1") {
  plotVar("m1pt", 20., cuts, tcuts);
  plotVar("m2pt", 20., cuts, tcuts);
  plotVar("pt", 30., cuts, tcuts);
  plotVar("prob", 1., cuts, tcuts);
  plotVar("alpha", 0.1, cuts, tcuts);
  plotVar("flsxy", 30., cuts, tcuts);

  // -- no chan cuts
  plotVar("TMath::Abs(m1eta)", 2.2, "gmuid&&bdt>0.34");
  plotVar("TMath::Abs(m2eta)", 2.2, "gmuid&&bdt>0.34");

}
