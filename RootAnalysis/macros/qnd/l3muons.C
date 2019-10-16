#include "../common/util.hh"

// ----------------------------------------------------------------------
TH1D *getEff(string filename, string suffix, double eps) {
  TFile *f = TFile::Open(filename.c_str());
  TTree *t = (TTree*)f->Get("candAnaMuMu/events");

  double xbins[] = {0., 2.e-12+eps, 4.e-12+eps, 6.e-12+eps, 8.e-12+eps, 10.e-12+eps, 20.e-12+eps};
  TH1D *h1 = new TH1D(Form("hall_%s", suffix.c_str()), "", 6, xbins);
  TH1D *h2 = new TH1D(Form("hcut_%s", suffix.c_str()), "", 6, xbins);
  TH1D *he = new TH1D(Form("heff_%s", suffix.c_str()), "", 6, xbins); he->Sumw2();

  t->Draw(Form("tau>>hall_%s", suffix.c_str()), "chan==0&&m2pt>4");
  t->Draw(Form("tau>>hcut_%s", suffix.c_str()), "chan==0&&m2pt>4&&hlt1");

  he->Divide(h2, h1, 1., 1., "b");
  return he;
}


// ----------------------------------------------------------------------
void l3muons() {

  TH1D *h02  = getEff("hlt-dbx28.root", "02", 0.); setHist(h02, kRed, 24); setTitles(h02, "t [sec]", "Efficiency(HLT_DoubleMu4_3_Bs_v7)", 0.04, 1.2, 1.4);
  TH1D *h05  = getEff("hlt-dbx30.root", "05", 0.1e-12); setHist(h05, kBlue, 24);
  TH1D *h10  = getEff("hlt-dbx31.root", "10", 0.2e-12); setHist(h10, kGreen+2, 24);
  TH1D *h20  = getEff("hlt-dbx33.root", "20", 0.4e-12); setHist(h20, kGreen+1, 24);
  TH1D *h100 = getEff("hlt-dbx34.root", "100", 0.5e-12); setHist(h100, kGreen+2, 24);
  TH1D *h200 = getEff("hlt-dbx29.root", "200", 0.6e-12); setHist(h200, kBlue, 24);

  shrinkPad(0.1, 0.13, 0.13);
  h02->SetMaximum(1.0);  h02->SetMinimum(0.0);
  h02->Draw("");
  h05->Draw("same");
  h10->Draw("same");
  //  h20->Draw("same");
  //  h100->Draw("same");
  //  h200->Draw("same");

  TLegend *tle = new TLegend(0.38, 0.70, 0.84, 0.85);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->SetTextSize(0.03);

  tle->SetHeader("L3Muons: tkTrajMaxDXYBeamSpot");
  tle->AddEntry(h02, "0.2cm", "p");
  tle->AddEntry(h05, "0.5cm", "p");
  tle->AddEntry(h10, "1.0cm", "p");
  tle->Draw();

  c0->SaveAs("qnd/l3muons-tkTrajMaxDXYBeamSpot.pdf");
}
