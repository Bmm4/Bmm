void jetht(string filename = "cbmm2-data-bmm2JetHT2016.bmmReader.2016.root") {
  TFile *file = TFile::Open(filename.c_str());
  int NBINS(7);
  double MAX(8.e-12);
  double eps(0.);
  double xbins[] = {0., 0.5e-12, 1.e-12, 1.5e-12, 2.e-12+eps, 4.e-12+eps, 6.e-12+eps, 10.e-12+eps};
  //  TH1D *heff = new TH1D("heff", "", NBINS, 0., MAX);
  TH1D *heff = new TH1D("heff", "", NBINS, xbins);
  heff->Sumw2();
  setTitles(heff, "decay time [sec]", "Eff((HLT_DoubleMu4_3_Jpsi_Displaced | m2pt>4)", 0.04, 1.2, 1.6, 0.04);
  // TH1D *h0 = new TH1D("h0", "", NBINS, 0., MAX);
  // TH1D *h1 = new TH1D("h1", "", NBINS, 0., MAX);
  TH1D *h0 = new TH1D("h0", "", NBINS, xbins);
  TH1D *h1 = new TH1D("h1", "", NBINS, xbins);
  TTree *events = (TTree*)file->Get("candAnaBu2JpsiK/events");
  events->Draw("tau>>h0", "m2pt>4.&&tau>0.");
  events->Draw("tau>>h1", "m2pt>4.&&tau>0.&&hlt1");
  heff->Divide(h1, h0, 1., 1., "b");
  shrinkPad(0.13, 0.15, 0.13);
  gStyle->SetErrorX(0.5);
  heff->Draw("b");
  heff->SetMinimum(0.01);
  gStyle->SetOptStat(0);
  tl->SetTextSize(0.04);
  tl->DrawLatexNDC(0.6, 0.85, Form("Entries: %4d", static_cast<int>(h1->GetSumOfWeights())));
  tl->DrawLatexNDC(0.6, 0.81, Form("Mean:   %.1e", h1->GetMean()));

  c0->SaveAs("qnd/eff-tau-jetht.pdf");
}
