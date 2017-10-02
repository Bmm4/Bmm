// ----------------------------------------------------------------------
void altCuts(string var, string cut1, string cut2, string cut3, double lo, double hi, string pdfname = "nada") {
  TFile *f1 = TFile::Open("/scratch/ursl/bmm4/v06/bmm-mc-Summer16_private-BuToJpsiKp-v06.root");
  f1->cd("candAnaBu2JpsiK");
  TH1D *h1 = new TH1D("h1", "", 40, lo, hi); h1->SetLineColor(kBlue);
  h1->SetXTitle(var.c_str());
  TH1D *h2 = new TH1D("h2", "", 40, lo, hi); h2->SetLineColor(kRed);

  TTree *t = (TTree*)gDirectory->Get("events");
  t->Draw(Form("%s >> h1", var.c_str()), Form("%s && %s", cut1.c_str(), cut2.c_str()));
  t->Draw(Form("%s >> h2", var.c_str()), Form("%s && %s", cut1.c_str(), cut3.c_str()));

  h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  double ymax = h1->GetMaximum();
  if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();
  h1->SetMaximum(1.3*ymax);

  gStyle->SetOptStat(0);
  h1->Draw();
  h2->Draw("samehist");

  TLatex tl;
  tl.SetTextSize(0.017);
  tl.SetTextColor(kBlack);
  tl.DrawLatexNDC(0.01, 0.97, Form("%s", cut1.c_str()));
  tl.SetTextColor(kBlue);
  tl.DrawLatexNDC(0.01, 0.94, Form("%s", cut2.c_str()));
  tl.SetTextColor(kRed);
  tl.DrawLatexNDC(0.01, 0.91, Form("%s", cut3.c_str()));

  if (pdfname != "nada") {
    c0->SaveAs(pdfname.c_str());
  }

}



// ----------------------------------------------------------------------
void twoCuts(string var, string cut1, string cut2, double lo, double hi, string pdfname = "nada") {
  TFile *f1 = TFile::Open("/scratch/ursl/bmm4/v06/bmm-mc-Summer16_private-BuToJpsiKp-v06.root");
  f1->cd("candAnaBu2JpsiK");
  TH1D *h1 = new TH1D("h1", "", 40, lo, hi); h1->SetLineColor(kBlue);
  TH1D *h2 = new TH1D("h2", "", 40, lo, hi); h2->SetLineColor(kRed);

  TTree *t = (TTree*)gDirectory->Get("events");
  t->Draw(Form("%s >> h1", var.c_str()), cut1.c_str());
  t->Draw(Form("%s >> h2", var.c_str()), Form("%s && %s", cut1.c_str(), cut2.c_str()));

  h1->Draw();
  h2->Scale(h1->GetSumOfWeights()/h2->GetSumOfWeights());
  h2->Draw("samehist");

  TLatex tl;
  tl.SetTextSize(0.017);
  tl.SetTextColor(kBlue);
  tl.DrawLatexNDC(0.01, 0.97, Form("%s", var.c_str()));
  tl.DrawLatexNDC(0.01, 0.94, Form("%s", cut1.c_str()));
  tl.SetTextColor(kRed);
  tl.DrawLatexNDC(0.01, 0.91, Form("&& %s", cut2.c_str()));

  if (pdfname != "nada") {
    c0->SaveAs(pdfname.c_str());
  }

}

// ----------------------------------------------------------------------
void allAltCuts() {
  altCuts("pvip", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 0.015, "pvip.pdf");
  altCuts("pvips", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 4, "pvips.pdf");
  altCuts("maxdoca", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 0.015, "maxdoca.pdf");
  altCuts("alpha", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 0.2, "alpha.pdf");
  altCuts("pt", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 40, "pt.pdf");
  altCuts("m1pt", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 25, "m1pt.pdf");
  altCuts("m2pt", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 20, "m2pt.pdf");
  altCuts("docatrk", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 0.3, "docatrk.pdf");
  altCuts("pvw8", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 1., "pvw8.pdf");

  altCuts("iso", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 1.02, "iso.pdf");
  altCuts("m1iso", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 1.02, "m1iso.pdf");
  altCuts("m2iso", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 1.02, "m2iso.pdf");
  altCuts("closetrk", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 40, "closetrk.pdf");
  altCuts("bdt", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", -1., 1., "bdt.pdf");

  altCuts("pvn", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 60., "pvn.pdf");
  altCuts("pvntrk", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 120., "pvntrk.pdf");
  altCuts("pv2ntrk", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 120., "pv2ntrk.pdf");
  altCuts("TMath::Abs(dzmin)", "chan==0&&m2pt>4&&psipt>7&&psiflsxy>3&&psicosa>0.9&&maxdoca<0.5&&psiprob>0.1&&hlt1&&l1t", "fls3d>50", "fls3d<12", 0., 2., "dzmin.pdf");
}
