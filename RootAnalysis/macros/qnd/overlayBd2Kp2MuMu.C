void bla() {
  TFile *f0 = TFile::Open("cbmm-mc-Summer15_private-BdToKPi-2GM-v01-00.bmmReader.mix-Bd2KPi.root");
  TH1D *h0 = new TH1D("h0", "", 40, 4.8, 6.0);
  h0->SetLineColor(kRed);
  h0->SetLineWidth(2);
  TTree *t = (TTree*)f0->Get("candAnaMuMu/events");
  t->Draw("1.01*m>>h0", "fls3d>13 && alpha<0.05 && iso>0.8 && docatrk<0.015 && pvip<0.008 && pvips<2 && closetrk<3 && chi2/dof < 2.2");

  TFile *f1 = TFile::Open("/scratch/ursl/bmm4/v01/bmm-mc-Summer15_private-BdToMuMu-v01.root");
  TH1D *h1 = new TH1D("h1", "", 40, 4.8, 6.0);
  h1->GetXaxis()->SetTitle("m[GeV]");
  h1->SetLineWidth(2);
  t = (TTree*)f1->Get("candAnaMuMu/events");
  t->Draw("m>>h1", "fls3d>13 && alpha<0.05 && iso>0.8 && docatrk<0.015 && pvip<0.008 && pvips<2 && closetrk<3 && chi2/dof < 2.2");

  h1->Scale(h0->GetSumOfWeights()/h1->GetSumOfWeights());

  gStyle->SetOptStat(0);

  h1->Draw();

  h0->Draw("same");

  TLatex *tl = new TLatex();
  tl->SetTextSize(0.04);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.52, 0.82, "preselection and ");
  tl->DrawLatexNDC(0.52, 0.77, "2 global muons");
  tl->SetTextSize(0.04);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.52, 0.60, "B^{0} #rightarrow #mu#mu");
  tl->SetTextColor(kRed);
  tl->DrawLatexNDC(0.52, 0.55, "B^{0} #rightarrow K#pi (m #times= 1.01)");

  c0->SaveAs("lineshapeComparison-Bd2MuMu-Bd2KPi.pdf") ;

}
