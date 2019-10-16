#include "../common/util.hh"

// plot trigger efficiency vs a variable

// ----------------------------------------------------------------------
void bc(string mode = "BdToPiMuNu", string fname = "nada") {
  TFile *f0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer16_private-BcToJpsiMuNu-2016-s01.root");
  TFile *f1(0);
  if (fname == "nada") {
    f1 = TFile::Open(Form("/scratch/ursl/bmm4/s01/bmm-mc-combined-%s-2016BF-s01.root", mode.c_str()));
  } else {
    f1 = TFile::Open(Form("/scratch/ursl/bmm4/s01/%s", fname.c_str()));
  }

  TTree *t0 = (TTree*)f0->Get("candAnaMuMu/events");
  TTree *t1 = (TTree*)f1->Get("candAnaMuMu/events");

  TH1D *h0 = new TH1D("h0", "", 10, 5., 6.0); setHist(h0, kRed); setTitles(h0, "m_{#mu#mu} [GeV]", "a.u.");
  TH1D *h1 = new TH1D("h1", "", 10, 5., 6.0); setHist(h1, kBlue);
  TH1D *hr = new TH1D("hr", "", 10, 5., 6.0); setHist(hr, kBlack); hr->Sumw2();

  t0->Draw("m>>h0", "bdt>0.");
  t1->Draw("m>>h1", "bdt>0.");

  h1->Scale(h0->GetSumOfWeights()/h1->GetSumOfWeights());

  h0->Draw("e");
  h1->Draw("esame");
  //  h1e->Draw("same");
  hr->Divide(h0, h1, 1., 1.);
  hr->Scale(5.);
  hr->Fit("pol1", "", "same");
  hr->GetFunction("pol1")->SetLineColor(kBlack);

  tl->SetTextSize(0.04);
  tl->DrawLatexNDC(0.2, 0.92, mode.c_str());

  TLegend *tle = new TLegend(0.6, 0.70, 0.8, 0.85);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->SetTextSize(0.03);

  tle->SetHeader(Form("slope = %3.1f #pm %3.1f",
		      hr->GetFunction("pol1")->GetParameter(1),
		      hr->GetFunction("pol1")->GetParError(1)));
  tle->AddEntry(h0, "Bc2JpsiMuNu", "p");
  tle->AddEntry(h1, mode.c_str(), "p");
  tle->AddEntry(hr, "ratio #times 5", "p");
  tle->Draw();

  c0->SaveAs(Form("qnd/bc-%s.pdf", mode.c_str()));
}


// ----------------------------------------------------------------------
void plotAll() {
  bc("BdToPiMuNu");
  bc("BsToKMuNu");
  bc("LbToPMuNu");
  bc("BdToPiMuMu", "bmm-mc-Summer16_private-BdToPiMuMu-2016BF-s01.root");
  bc("BuToPiMuMu", "bmm-mc-Summer16_private-BuToPiMuMu-2016BF-s01.root");
}
