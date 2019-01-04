
// ----------------------------------------------------------------------
TH1D* fakeRateGM(string filename = "/scratch/ursl/bmm4/s01/fakeMC2012-s01.root", float etacut = 0.7, int id = 321, int q = +1, double eps = 0.) {
  TFile *f = TFile::Open(filename.c_str());
  TTree *t = (TTree*)f->Get("candAnaFakeMC/fakeTree");

  const Int_t nbins=7;
  Double_t xbins[nbins+1] = {0.+eps, 4.+eps, 6.+eps, 8.+eps, 10.+eps, 15.+eps, 20.+eps, 50.+eps};

  TH1D *ht = new TH1D("htgm", "ht", nbins, xbins); ht->Sumw2();
  TH1D *hm = new TH1D("hmgm", "hm", nbins, xbins); hm->Sumw2();

  t->Draw("pt>>htgm", Form("hp && TMath::Abs(eta)<%3.1f && id==%d && q==%d && !cowboy", etacut, id, q));
  t->Draw("pt>>hmgm", Form("hp && TMath::Abs(eta)<%3.1f && id==%d && q==%d && !cowboy && gm>0", etacut, id, q));

  TH1D  *hr = (TH1D*)ht->Clone("hrgm"); hr->Reset(); hr->SetTitle("fake rate");
  hr->SetMinimum(0.);
  hr->Divide(hm, ht, 1., 1., "b");

  return hr;
}

// ----------------------------------------------------------------------
TH1D* fakeRateMVA(string filename = "/scratch/ursl/bmm4/s01/fakeMC2012-s01.root", float etacut = 0.7, int id = 321, int q = +1) {
  TFile *f = TFile::Open(filename.c_str());
  TTree *t = (TTree*)f->Get("candAnaFakeMC/fakeTree");

  const Int_t nbins = 7;
  Double_t xbins[nbins+1] = {0., 4., 6., 8., 10., 15., 20., 50.};

  TH1D *ht = new TH1D("htmva", "ht", nbins, xbins); ht->Sumw2();
  TH1D *hm = new TH1D("hmmva", "hm", nbins, xbins); hm->Sumw2();

  float bdtcut(0.58);
  if (string::npos != filename.find("2012")) {
    bdtcut = 0.55;
  }
  t->Draw("pt>>htmva", Form("hp && TMath::Abs(eta)<%3.1f && id==%d && q==%d", etacut, id, q));
  t->Draw("pt>>hmmva", Form("hp && TMath::Abs(eta)<%3.1f && id==%d && q==%d && gm>0 && bdt > %4.3f", etacut, id, q, bdtcut));

  TH1D  *hr = (TH1D*)ht->Clone("hrmva"); hr->Reset(); hr->SetTitle("fake rate");
  hr->SetMinimum(0.);
  hr->Divide(hm, ht, 1., 1.);

  return hr;
}


// ----------------------------------------------------------------------
void plot(float etacut = 0.7, int id = 321, int q = +1) {
  TH1D *hr1 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-s01.root", etacut, id, q);
  TH1D *hr2 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2016-s01.root", etacut, id, q);

  TH1D *hm1 = fakeRateMVA("/scratch/ursl/bmm4/s01/fakeMC2012-s01.root", etacut, id, q);
  TH1D *hm2 = fakeRateMVA("/scratch/ursl/bmm4/s01/fakeMC2016-s01.root", etacut, id, q);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetLogy(1);

  hr1->SetLineColor(kRed);
  hr1->SetMarkerColor(kRed);
  hr1->SetMarkerStyle(24);
  hr1->SetMaximum(0.1);
  hr1->SetMinimum(0.0001);
  if (13 == id) {
    hr1->SetMaximum(1.0);
    hr1->SetMinimum(0.);
    gPad->SetLogy(0);
  }

  hr2->SetLineColor(kBlue);
  hr2->SetMarkerColor(kBlue);
  hr2->SetMarkerStyle(25);

  hr1->Draw("e");
  hr2->Draw("esame");

  hm1->SetLineColor(kRed);
  hm1->SetMarkerColor(kRed);
  hm2->SetLineColor(kBlue);
  hm2->SetMarkerColor(kBlue);

  hm1->Draw("esame");
  hm2->Draw("esame");

  TLegend *tle = new TLegend(0.5, 0.3, 0.8, 0.5);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->AddEntry(hr1, "Global muons, 2012", "p");
  tle->AddEntry(hr2, "Global muons, 2016", "p");

  tle->AddEntry(hm1, "MVA muons, 2012", "p");
  tle->AddEntry(hm2, "MVA muons, 2016", "p");
  tle->Draw();

  tl->DrawLatexNDC(0.2, 0.92, Form("ID = %d, q = %+d, |eta| < %3.2f", id, q, etacut));
  c0->SaveAs(Form("fakeRate-%3.2f-%d-%s.pdf", etacut, id, (q>0?"pos":"neg")));
}


// ----------------------------------------------------------------------
void plotAll(float etacut = 0.7) {
  plot(etacut, 321, +1);
  plot(etacut, 321, -1);

  plot(etacut, 211, +1);
  plot(etacut, 211, -1);

  plot(etacut, 13, +1);
  plot(etacut, 13, -1);
}


// ----------------------------------------------------------------------
void plotSplitFiles(float etacut = 0.7, int id = 321, int q = +1) {
  double eps(0.1);
  TH1D *hr0 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-split/cbdmm-s01.root", etacut, id, q, eps);
  TH1D *hr1 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-split/cbdkk-s01.root", etacut, id, q, 2*eps);
  TH1D *hr2 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-split/cbdkpi-s01.root", etacut, id, q, 3*eps);
  TH1D *hr3 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-split/cbdpipi-s01.root", etacut, id, q, 4*eps);
  TH1D *hr4 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-split/cbskpi-s01.root", etacut, id, q, 5*eps);
  TH1D *hr5 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-split/cbspipi-s01.root", etacut, id, q, 6*eps);
  TH1D *hr6 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-split/cbskk-s01.root", etacut, id, q, 7*eps);

  TH1D *hr10= fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-split/cbdkk-s01-a.root", etacut, id, q, 8*eps);

  //  TH1D *hr4 = fakeRateGM("/scratch/ursl/bmm4/s01/fakeMC2012-s01/cbdpipi-s01-00.root", etacut, id, q);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gPad->SetLogy(1);


  hr0->SetLineColor(kRed);
  hr0->SetMarkerColor(kRed);
  hr0->SetMarkerStyle(20);
  hr0->SetMaximum(1.0);
  hr0->SetMinimum(0.0001);
  if (13 == id) {
    hr0->SetMaximum(1.0);
    hr0->SetMinimum(0.);
    gPad->SetLogy(0);
  }

  hr1->SetLineColor(kBlack);
  hr1->SetMarkerColor(kBlack);
  hr1->SetMarkerStyle(20);

  hr2->SetLineColor(kBlue);
  hr2->SetMarkerColor(kBlue);
  hr2->SetMarkerStyle(21);

  hr3->SetLineColor(kMagenta);
  hr3->SetMarkerColor(kMagenta);
  hr3->SetMarkerStyle(22);

  hr4->SetLineColor(kCyan);
  hr4->SetMarkerColor(kCyan);
  hr4->SetMarkerStyle(24);

  hr5->SetLineColor(kGreen+1);
  hr5->SetMarkerColor(kGreen+1);
  hr5->SetMarkerStyle(25);

  hr6->SetLineColor(kBlue-2);
  hr6->SetMarkerColor(kBlue-2);
  hr6->SetMarkerStyle(26);

  hr10->SetLineColor(kBlack);
  hr10->SetMarkerColor(kBlack);
  hr10->SetMarkerStyle(27);

  hr0->Draw("e");
  hr1->Draw("esame");
  hr2->Draw("esame");
  hr3->Draw("esame");
  hr4->Draw("esame");
  hr5->Draw("esame");
  hr6->Draw("esame");
  hr10->Draw("esame");

  TLegend *tle = new TLegend(0.5, 0.2, 0.8, 0.4);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->SetHeader("Global muon fake rate");
  tle->AddEntry(hr0, "Bd2MM, 2012", "p");
  tle->AddEntry(hr1, "Bd2KK, 2012", "p");
  tle->AddEntry(hr10, "Bd2KK', 2012", "p");
  tle->AddEntry(hr2, "Bd2KPi, 2012", "p");
  tle->AddEntry(hr3, "Bd2PiPi, 2012", "p");
  tle->AddEntry(hr4, "Bs2KPi, 2012", "p");
  tle->AddEntry(hr5, "Bs2PiPi, 2012", "p");
  tle->AddEntry(hr6, "Bs2KK, 2012", "p");

  tle->Draw();

  tl->DrawLatexNDC(0.2, 0.92, Form("ID = %d, q = %+d, |eta| < %3.2f", id, q, etacut));
  c0->SaveAs(Form("fakeRate-%3.2f-%d-%s.pdf", etacut, id, (q>0?"pos":"neg")));
}

// ----------------------------------------------------------------------
void plotAllSplitFiles(float etacut = 0.7) {
  plotSplitFiles(etacut, 321, -1);
  plotSplitFiles(etacut, 321, +1);

  plotSplitFiles(etacut, 211, -1);
  plotSplitFiles(etacut, 211, +1);
}
