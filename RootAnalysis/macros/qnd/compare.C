void compare(string hname = "bdtScan_SSB_chan0") {
  TFile *f0 = TFile::Open("results-409/scanBDT-2016BF.root");
  TH1D *h0 = (TH1D*)f0->Get(hname.c_str());
  h0->SetMaximum(2.0);
  h0->Draw();

  TFile *f1 = TFile::Open("results-419/scanBDT-2016BF.root");
  TH1D *h1 = (TH1D*)f1->Get(hname.c_str());
  h1->SetLineColor(kBlue);
  h1->Draw("same");

  TFile *f2 = TFile::Open("results-429/scanBDT-2016BF.root");
  TH1D *h2 = (TH1D*)f2->Get(hname.c_str());
  h2->SetLineColor(kRed);
  h2->Draw("same");

}
