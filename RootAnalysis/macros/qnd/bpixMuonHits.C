void bpixMuonHits(string era = "BF") {

  if ("all" == era) {
    bpixMuonHits("BF");
    bpixMuonHits("GH");
    return;
  }

  TFile *f0;
  if (era == "BF") {
    f0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016BF-s01.root");
  } else {
    f0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016GH-s01.root");
  }

  TTree *t0 = (TTree*)f0->Get("candAnaBu2JpsiK/events");

  TH1D *h0a = new TH1D("h0a", "BPIX hits (PU < 10)", 5, 0., 5.);
  TH1D *h0b = new TH1D("h0b", "BPIX hits (10 < PU < 20)", 5, 0., 5.);
  TH1D *h0c = new TH1D("h0c", "BPIX hits (30 < PU)", 5, 0., 5.);

  TH1D *hw0a = new TH1D("hw0a", "BPIX hits (PU < 10)", 5, 0., 5.);
  TH1D *hw0b = new TH1D("hw0b", "BPIX hits (10 < PU < 20)", 5, 0., 5.);
  TH1D *hw0c = new TH1D("hw0c", "BPIX hits (30 < PU)", 5, 0., 5.);

  t0->Draw("m1pix>>h0a", "pvn < 10");
  t0->Draw("m1pix>>h0b", "10 < pvn && pvn < 20");
  t0->Draw("m1pix>>h0c", "30 < pvn");

  t0->Draw("m1pix>>hw0a", "cw8*(pvn < 10)");
  t0->Draw("m1pix>>hw0b", "cw8*(10 < pvn && pvn < 20)");
  t0->Draw("m1pix>>hw0c", "cw8*(30 < pvn)");

  h0a->Scale(1./h0a->GetSumOfWeights());
  h0b->Scale(1./h0b->GetSumOfWeights());
  h0c->Scale(1./h0c->GetSumOfWeights());

  hw0a->Scale(1./hw0a->GetSumOfWeights());
  hw0b->Scale(1./hw0b->GetSumOfWeights());
  hw0c->Scale(1./hw0c->GetSumOfWeights());


  c0->Clear();
  c0->Divide(2,2);

  c0->cd(1);
  hw0a->Draw("hist");
  h0a->Draw("samee");

  c0->cd(2);
  hw0b->Draw("hist");
  h0b->Draw("samee");

  c0->cd(3);
  hw0c->Draw("hist");
  h0c->Draw("samee");

  TH1D *hu = new TH1D("hu", "Mean BPIX hits ", 3, 0., 3.); hu->SetMarkerColor(kBlack); hu->SetMarkerStyle(24);
  hu->GetXaxis()->SetBinLabel(1, "PU < 10");
  hu->GetXaxis()->SetBinLabel(2, "10 < PU < 20");
  hu->GetXaxis()->SetBinLabel(3, "30 < PU");
  TH1D *hw = new TH1D("hw", "", 3, 0., 3.); hw->SetMarkerColor(kBlue); hw->SetMarkerStyle(25);

  hu->SetBinContent(1, h0a->GetMean());
  hu->SetBinContent(2, h0b->GetMean());
  hu->SetBinContent(3, h0c->GetMean());

  hw->SetBinContent(1, hw0a->GetMean());
  hw->SetBinContent(2, hw0b->GetMean());
  hw->SetBinContent(3, hw0c->GetMean());

  c0->cd(4);
  hu->Draw("p");
  hw->Draw("psame");

  tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.3, 0.3, "weighted");
  tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.3, 0.35, "unweighted");
  c0->SaveAs(Form("s01/bpixmuonhits-%s.pdf", era.c_str()));

}
