void compareBDT(string test = "fls3d > 10", string base = "dcand && chan==0 && hlt1 && tos && bdt>-1", string fname = "/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BsToJpsiPhi-s01.root") {
  TFile *f0 = TFile::Open(fname.c_str());
  TTree *t = (TTree*)f0->Get("candAnaBs2JpsiPhi/events");
  TH1D *h0 = new TH1D("h0", "h0", 100, -1., 1.);
  TH1D *h1 = new TH1D("h1", "h1", 100, -1., 1.);

  t->Draw("bdt>>h0", base.c_str());
  t->Draw("bdt>>h1", Form("%s && %s", base.c_str(), test.c_str()));

  h0->Scale(1./h0->GetSumOfWeights());
  h1->Scale(1./h1->GetSumOfWeights());
  if (h1->GetMaximum() > h0->GetMaximum()) h0->SetMaximum(1.1*h1->GetMaximum());

  h0->SetLineColor(kBlack);
  h1->SetLineColor(kRed);


  h0->Draw("hist");
  h1->Draw("histsame");

  tl->SetTextSize(0.02);
  tl->SetTextColor(kBlack);
  tl->DrawLatexNDC(0.22, 0.8, base.c_str());
  tl->SetTextColor(kRed);
  tl->DrawLatexNDC(0.22, 0.75, ("&& " + test).c_str());

}
