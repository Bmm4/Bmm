// ----------------------------------------------------------------------
double adhocfunction(double *x, double *par) {
  // par[0] -> B0 const
  // par[1] -> B0 mean
  // par[2] -> B0 sigma

  // par[3] -> Bs const
  // par[4] -> Bs mean
  // par[5] -> Bs sigma

  // par[6] -> pol0

  double result(0.), arg(0.), fitval(0.);
  if (par[2] > 0.) {
    arg = (x[0] - par[1]) / par[2];
    fitval =  par[0]*TMath::Exp(-0.5*arg*arg);
    result += fitval;
  }

  if (par[5] > 0.) {
    arg = (x[0] - par[4]) / par[5];
    fitval =  par[3]*TMath::Exp(-0.5*arg*arg);
    result += fitval;
  }

  result += par[6];
  return result;

}


// ----------------------------------------------------------------------
void invertedMuonID2Comb(string era = "2016BF", float bdt = -99., TH1D *hBd = 0,  TH1D *hBs = 0) {
  if ("all" == era) {
    invertedMuonID2Comb("2016BF");
    invertedMuonID2Comb("2016GH");
    invertedMuonID2Comb("2011");
    invertedMuonID2Comb("2012");
    return;
  }

  TFile *fd(0), *fm0(0), *fm1(0);
  float bdtcut(0.);
  cout << "======================================================================" << endl;
  if (era == "2016BF") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016BF-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BdToMuMu-2016BF-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BsToMuMu-2016BF-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.30;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2016GH") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016GH-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BdToMuMu-2016GH-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BsToMuMu-2016GH-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.31;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2012") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2012-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BdToMuMu-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BsToMuMu-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.34;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2011") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BdToMuMu-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BsToMuMu-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.28;
    } else {
      bdtcut = bdt;
    }
  }

  TTree *td = (TTree*)fd->Get("candAnaMuMu/events");
  TTree *tm0 = (TTree*)fm0->Get("candAnaMuMu/events");
  TTree *tm1 = (TTree*)fm1->Get("candAnaMuMu/events");

  TH1D *hd0 = new TH1D("hd0", "mass", 60, 4.8, 6.0);
  TH1D *hm0 = new TH1D("hm0", "mass", 60, 4.8, 6.0);
  TH1D *hm1 = new TH1D("hm1", "mass", 60, 4.8, 6.0);
  int ichan(0);
  string cut(Form("chan==%d && !gmuid &&hlt1&&tos&&l1t &&bdt>%4.3f && m < 6 && m > 4.8", ichan, bdtcut));
  cout << cut << endl;
  td->Draw("m>>hd0", cut.c_str());
  tm0->Draw("m>>hm0", cut.c_str());
  tm1->Draw("m>>hm1", cut.c_str());

  c0->Clear();
  // c0->Divide(2,2);
  // c0->cd(1);
  hd0->Draw("hist");
  // B0 -> mu mu
  // c0->cd(2);
  hm0->Fit("gaus");
  double p0 = hm0->GetFunction("gaus")->GetParameter(1);
  double s0 = hm0->GetFunction("gaus")->GetParameter(2);
  cout << "p0: " << p0 << " s0: " << s0 << endl;
  // Bs -> mu mu
  // c0->cd(3);
  hm1->Fit("gaus");
  double p1 = hm1->GetFunction("gaus")->GetParameter(1);
  double s1 = hm1->GetFunction("gaus")->GetParameter(2);
  cout << "p1: " << p1 << " s1: " << s1 << endl;

  // c0->cd(4);
  TH1D *hFit = (TH1D*)hd0->Clone("hFit");
  hFit->SetTitle("");
  hFit->SetMinimum(0.);
  TF1 *f1 = new TF1("f1", adhocfunction, 4.8, 6.0, 7);
  f1->SetParameter(0, 1);  f1->SetParLimits(0, 0., 1000.);
  f1->SetParameter(1, p0); f1->FixParameter(1, p0);
  f1->SetParameter(2, s0); f1->FixParameter(2, s0);
  f1->SetParameter(3, 1);  f1->SetParLimits(3, 0., 1000.);
  f1->SetParameter(4, p1); f1->FixParameter(4, p1);
  f1->SetParameter(5, s1); f1->FixParameter(5, s1);
  f1->SetParameter(6, 1);

  hFit->Fit(f1, "RL", "", 5.1, 6.0);

  double nb0 = hFit->GetFunction("f1")->GetParameter(0);
  double eb0 = hFit->GetFunction("f1")->GetParError(0);
  double nbs = hFit->GetFunction("f1")->GetParameter(3);
  double ebs = hFit->GetFunction("f1")->GetParError(3);
  double nco = hFit->GetFunction("f1")->GetParameter(6);
  double eco = hFit->GetFunction("f1")->GetParError(6);

  TF1 *f0 = (TF1*)hm0->GetFunction("gaus")->Clone("f0");
  TF1 *fs = (TF1*)hm1->GetFunction("gaus")->Clone("fs");

  f0->SetParameter(0, nb0); f0->SetLineColor(kRed+2);
  fs->SetParameter(0, nbs); fs->SetLineColor(kBlue);
  f0->Draw("same");
  fs->Draw("same");
  tl->DrawLatexNDC(0.10, 0.92, era.c_str());
  tl->DrawLatexNDC(0.40, 0.92, Form("BDT > %3.2f", bdtcut));
  double iB0   = f0->Integral(5.1, 5.6);
  double iBs   = fs->Integral(5.1, 5.6);
  float  nB0   = f0->Integral(5.2, 5.3)/hFit->GetBinWidth(1);
  float  nB0E  = eb0/nb0*nB0;
  float  nBs   = fs->Integral(5.3, 5.45)/hFit->GetBinWidth(1);
  float  nBsE  = ebs/nbs*nBs;
  double iComb = nco*(6.0-5.6);
  double ratioVs(0.);
  double relErrs(0.);
  double ratioEs(0.01);
  double ratioVd(0.);
  double relErrd(0.);
  double ratioEd(0.01);
  if (iComb > 1.e-3) {
    ratioVs = iBs/iComb;
    relErrs = TMath::Sqrt((ebs/nbs)*(ebs/nbs) + (eco/nco)*(eco/nco));
    ratioEs = relErrs * ratioVs;

    ratioVd = iB0/iComb;
    relErrd = TMath::Sqrt((eb0/nb0)*(eb0/nb0) + (eco/nco)*(eco/nco));
    ratioEd = relErrd * ratioVd;
 }

  cout << "ratio of B0/comb: " << iB0 << "/" << iComb << " = " << ratioVd << " +/- " << ratioEd << " relErr: " << relErrd << endl;
  cout << "ratio of Bs/comb: " << iBs << "/" << iComb << " = " << ratioVs << " +/- " << ratioEs << " relErr: " << relErrs << endl;

  ofstream TEX(Form("invertedMuonIDStudies.tex"), ios::app);
  TEX << Form("\\vdef{%ss01:nb0:bdt%3.2f:chan%d}  {%3.1f}", era.c_str(), bdtcut, ichan, nB0) << endl;
  TEX << Form("\\vdef{%ss01:eb0:bdt%3.2f:chan%d}  {%3.1f}", era.c_str(), bdtcut, ichan, nB0E) << endl;
  TEX << Form("\\vdef{%ss01:nbs:bdt%3.2f:chan%d}  {%3.1f}", era.c_str(), bdtcut, ichan, nBs) << endl;
  TEX << Form("\\vdef{%ss01:ebs:bdt%3.2f:chan%d}  {%3.1f}", era.c_str(), bdtcut, ichan, nBsE) << endl;
  TEX.close();

  c0->SaveAs(Form("invertedMuonIDcomb-%s-%3.2f.pdf", era.c_str(), bdtcut));

  if (hBd) {
    hBd->SetBinContent(hBd->FindBin(bdtcut), ratioVd);
    hBd->SetBinError(hBd->FindBin(bdtcut), ratioEd);
  }
  if (hBs) {
    hBs->SetBinContent(hBs->FindBin(bdtcut), ratioVs);
    hBs->SetBinError(hBs->FindBin(bdtcut), ratioEs);
  }
}


// ----------------------------------------------------------------------
void invertedMuonIDRaw(string era = "2016BF", float bdt = -99., TH1D *hBd = 0,  TH1D *hBs = 0) {
  if ("all" == era) {
    invertedMuonIDRaw("2016BF");
    invertedMuonIDRaw("2016GH");
    invertedMuonIDRaw("2011");
    invertedMuonIDRaw("2012");
    return;
  }

  TFile *fd(0), *fm0(0), *fm1(0);
  float bdtcut(0.);
  cout << "======================================================================" << endl;
  if (era == "2016BF") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016BF-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BdToMuMu-2016BF-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BsToMuMu-2016BF-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.30;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2016GH") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016GH-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BdToMuMu-2016GH-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BsToMuMu-2016GH-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.31;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2012") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2012-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BdToMuMu-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BsToMuMu-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.34;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2011") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BdToMuMu-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BsToMuMu-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.28;
    } else {
      bdtcut = bdt;
    }
  }

  TTree *td = (TTree*)fd->Get("candAnaMuMu/events");
  TTree *tm0 = (TTree*)fm0->Get("candAnaMuMu/events");
  TTree *tm1 = (TTree*)fm1->Get("candAnaMuMu/events");

  TH1D *hd0 = new TH1D("hd0", "mass", 60, 4.8, 6.0);
  TH1D *hm0 = new TH1D("hm0", "mass", 60, 4.8, 6.0);
  TH1D *hm1 = new TH1D("hm1", "mass", 60, 4.8, 6.0);
  int ichan(0);
  string cut(Form("chan==%d && !gmuid &&hlt1&&tos&&l1t &&bdt>%4.3f && m < 6 && m > 4.8", ichan, bdtcut));
  cout << cut << endl;
  td->Draw("m>>hd0", cut.c_str());
  tm0->Draw("m>>hm0", cut.c_str());
  tm1->Draw("m>>hm1", cut.c_str());

  c0->Clear();
  // c0->Divide(2,2);
  // c0->cd(1);
  hd0->Draw("hist");
  // B0 -> mu mu
  // c0->cd(2);
  hm0->Fit("gaus");
  double p0 = hm0->GetFunction("gaus")->GetParameter(1);
  double s0 = hm0->GetFunction("gaus")->GetParameter(2);
  cout << "p0: " << p0 << " s0: " << s0 << endl;
  // Bs -> mu mu
  // c0->cd(3);
  hm1->Fit("gaus");
  double p1 = hm1->GetFunction("gaus")->GetParameter(1);
  double s1 = hm1->GetFunction("gaus")->GetParameter(2);
  cout << "p1: " << p1 << " s1: " << s1 << endl;

  // c0->cd(4);
  TH1D *hFit = (TH1D*)hd0->Clone("hFit");
  hFit->SetTitle("");
  hFit->SetMinimum(0.);
  TF1 *f1 = new TF1("f1", adhocfunction, 4.8, 6.0, 7);
  f1->SetParameter(0, 1);  f1->SetParLimits(0, 0., 1000.);
  f1->SetParameter(1, p0); f1->FixParameter(1, p0);
  f1->SetParameter(2, s0); f1->FixParameter(2, s0);
  f1->SetParameter(3, 1);  f1->SetParLimits(3, 0., 1000.);
  f1->SetParameter(4, p1); f1->FixParameter(4, p1);
  f1->SetParameter(5, s1); f1->FixParameter(5, s1);
  f1->SetParameter(6, 1);

  hFit->Fit(f1, "RL", "", 5.1, 6.0);

  double nb0 = hFit->GetFunction("f1")->GetParameter(0);
  double eb0 = hFit->GetFunction("f1")->GetParError(0);
  double nbs = hFit->GetFunction("f1")->GetParameter(3);
  double ebs = hFit->GetFunction("f1")->GetParError(3);

  TF1 *f0 = (TF1*)hm0->GetFunction("gaus")->Clone("f0");
  TF1 *fs = (TF1*)hm1->GetFunction("gaus")->Clone("fs");

  f0->SetParameter(0, nb0); f0->SetLineColor(kRed+2);
  fs->SetParameter(0, nbs); fs->SetLineColor(kBlue);
  f0->Draw("same");
  fs->Draw("same");
  tl->DrawLatexNDC(0.10, 0.92, era.c_str());
  tl->DrawLatexNDC(0.40, 0.92, Form("BDT > %3.2f", bdtcut));
  double iB0   = f0->Integral(5.1, 5.6);
  double iBs   = fs->Integral(5.1, 5.6);
  float  nB0   = f0->Integral(5.2, 5.3)/hFit->GetBinWidth(1);
  float  nB0E  = eb0/nb0*nB0;
  float  nBs   = fs->Integral(5.3, 5.45)/hFit->GetBinWidth(1);
  float  nBsE  = ebs/nbs*nBs;
  double ratioVs(0.);
  double relErrs(0.);
  double ratioEs(0.01);
  double ratioVd(0.);
  double relErrd(0.);
  double ratioEd(0.01);

  ofstream TEX(Form("invertedMuonIDStudies.tex"), ios::app);
  TEX << Form("\\vdef{%ss01:nb0:bdt%3.2f:chan%d}  {%3.1f}", era.c_str(), bdtcut, ichan, nB0) << endl;
  TEX << Form("\\vdef{%ss01:eb0:bdt%3.2f:chan%d}  {%3.1f}", era.c_str(), bdtcut, ichan, nB0E) << endl;
  TEX << Form("\\vdef{%ss01:nbs:bdt%3.2f:chan%d}  {%3.1f}", era.c_str(), bdtcut, ichan, nBs) << endl;
  TEX << Form("\\vdef{%ss01:ebs:bdt%3.2f:chan%d}  {%3.1f}", era.c_str(), bdtcut, ichan, nBsE) << endl;
  TEX.close();

  c0->SaveAs(Form("invertedMuonIDcomb-%s-%3.2f.pdf", era.c_str(), bdtcut));

  if (hBd) {
    hBd->SetBinContent(hBd->FindBin(bdtcut), iB0/hFit->GetBinWidth(1));
    hBd->SetBinError(hBd->FindBin(bdtcut), iB0*(nB0E/nB0)/hFit->GetBinWidth(1));
  }
  if (hBs) {
    hBs->SetBinContent(hBs->FindBin(bdtcut), iBs/hFit->GetBinWidth(1));
    hBs->SetBinError(hBs->FindBin(bdtcut), iBs*(nBsE/nBs)/hFit->GetBinWidth(1));
  }
}

// ----------------------------------------------------------------------
void invertedMuonID(string era = "2016BF", float bdt = -99., TH1D *hresult = 0) {
  if ("all" == era) {
    invertedMuonID("2016BF");
    invertedMuonID("2016GH");
    return;
  }

  TFile *fd(0), *fm0(0), *fm1(0);
  float bdtcut(0.);
  cout << "======================================================================" << endl;
  if (era == "2016BF") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016BF-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BdToMuMu-2016BF-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BsToMuMu-2016BF-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.30;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2016GH") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016GH-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BdToMuMu-2016GH-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-RunIISpring16DR80-BsToMuMu-2016GH-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.31;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2012") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2012-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BdToMuMu-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BsToMuMu-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.34;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2011") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BdToMuMu-s01.root");
    fm1 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BsToMuMu-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.28;
    } else {
      bdtcut = bdt;
    }
  }

  TTree *td = (TTree*)fd->Get("candAnaMuMu/events");
  TTree *tm0 = (TTree*)fm0->Get("candAnaMuMu/events");
  TTree *tm1 = (TTree*)fm1->Get("candAnaMuMu/events");

  TH1D *hd0 = new TH1D("hd0", "mass", 60, 4.8, 6.0);
  TH1D *hm0 = new TH1D("hm0", "mass", 60, 4.8, 6.0);
  TH1D *hm1 = new TH1D("hm1", "mass", 60, 4.8, 6.0);
  string cut(Form("chan==0 && !gmuid &&hlt1&&tos&&l1t &&bdt>%4.3f && m < 6 && m > 4.8", bdtcut));
  cout << cut << endl;
  td->Draw("m>>hd0", cut.c_str());
  tm0->Draw("m>>hm0", cut.c_str());
  tm1->Draw("m>>hm1", cut.c_str());

  c0->Clear();
  // c0->Divide(2,2);
  // c0->cd(1);
  hd0->Draw("hist");
  // B0 -> mu mu
  // c0->cd(2);
  hm0->Fit("gaus");
  double p0 = hm0->GetFunction("gaus")->GetParameter(1);
  double s0 = hm0->GetFunction("gaus")->GetParameter(2);
  cout << "p0: " << p0 << " s0: " << s0 << endl;
  // Bs -> mu mu
  // c0->cd(3);
  hm1->Fit("gaus");
  double p1 = hm1->GetFunction("gaus")->GetParameter(1);
  double s1 = hm1->GetFunction("gaus")->GetParameter(2);
  cout << "p1: " << p1 << " s1: " << s1 << endl;

  // c0->cd(4);
  TH1D *hFit = (TH1D*)hd0->Clone("hFit");
  hFit->SetTitle("");
  hFit->SetMinimum(0.);
  TF1 *f1 = new TF1("f1", adhocfunction, 4.8, 6.0, 7);
  f1->SetParameter(0, 1);  f1->SetParLimits(0, 0., 1000.);
  f1->SetParameter(1, p0); f1->FixParameter(1, p0);
  f1->SetParameter(2, s0); f1->FixParameter(2, s0);
  f1->SetParameter(3, 1);  f1->SetParLimits(3, 0., 1000.);
  f1->SetParameter(4, p1); f1->FixParameter(4, p1);
  f1->SetParameter(5, s1); f1->FixParameter(5, s1);
  f1->SetParameter(6, 1);

  hFit->Fit(f1, "RL", "", 5.1, 6.0);

  double nb0 = hFit->GetFunction("f1")->GetParameter(0);
  double eb0 = hFit->GetFunction("f1")->GetParError(0);
  double nbs = hFit->GetFunction("f1")->GetParameter(3);
  double ebs = hFit->GetFunction("f1")->GetParError(3);

  TF1 *f0 = (TF1*)hm0->GetFunction("gaus")->Clone("f0");
  TF1 *fs = (TF1*)hm1->GetFunction("gaus")->Clone("fs");

  f0->SetParameter(0, nb0); f0->SetLineColor(kRed+2);
  fs->SetParameter(0, nbs); fs->SetLineColor(kBlue);
  f0->Draw("same");
  fs->Draw("same");
  tl->DrawLatexNDC(0.10, 0.92, era.c_str());
  tl->DrawLatexNDC(0.40, 0.92, Form("BDT > %3.2f", bdtcut));
  double iB0 = f0->Integral(5.1, 5.6);
  double iBs = fs->Integral(5.1, 5.6);
  double ratioV(0.);
  double relErr(0.);
  double ratioE(0.01);
  if (iBs > 1.e-3) {
    ratioV = iB0/iBs;
    relErr = TMath::Sqrt((eb0/nb0)*(eb0/nb0) + (ebs/nbs)*(ebs/nbs));
    ratioE = relErr * ratioV;
  }

  cout << "ratio of B0 and Bs: " << iB0 << "/" << iBs << " = " << ratioV << " +/- " << ratioE << " relErr: " << relErr << endl;

  c0->SaveAs(Form("invertedMuonID-%s-%3.2f.pdf", era.c_str(), bdtcut));

  if (hresult) {
    hresult->SetBinContent(hresult->FindBin(bdtcut), ratioV);
    hresult->SetBinError(hresult->FindBin(bdtcut), ratioE);
  }
}


// ----------------------------------------------------------------------
void plotRatio(string era = "2016BF") {
  gStyle->SetOptStat(0);

  if (era == "all") {
    plotRatio("2011");
    plotRatio("2012");
    plotRatio("2016BF");
    plotRatio("2016GH");
    return;
  }


  TH1D *hresult = new TH1D("hresult", Form("n(gauss(5.28))/n(gauss(5.37)) %s", era.c_str()), 40, 0., 0.4);
  hresult->GetXaxis()->SetTitle("BDT > ");
  TH1D *hbdtcut = new TH1D("hbdtcut", Form("n(gauss(5.28))/n(gauss(5.37)) %s", era.c_str()), 40, 0., 0.4);
  hbdtcut->SetMarkerColor(kBlue);
  invertedMuonID(era, 0.15, hresult);
  invertedMuonID(era, 0.20, hresult);
  invertedMuonID(era, 0.25, hresult);
  invertedMuonID(era, 0.30, hresult);
  invertedMuonID(era, -99., hbdtcut);

  c0->Clear();
  double ymax(50.);
  if (hresult->GetMaximum() < 5.)  ymax = 10.;
  if (hresult->GetMaximum() < 0.5) ymax = 1.0;
  hresult->SetMinimum(-ymax);
  hresult->SetMaximum(ymax);
  hresult->Draw();
  hbdtcut->Draw("same");
  pl->DrawLine(0., 0., 0.4, 0.);
  c0->SaveAs(Form("invertedMuonID-%s-result.pdf", era.c_str()));
}


// ----------------------------------------------------------------------
void plotRatioRaw(string era = "2016BF") {
  gStyle->SetOptStat(0);

  if (era == "all") {
    plotRatioRaw("2011");
    plotRatioRaw("2012");
    plotRatioRaw("2016BF");
    plotRatioRaw("2016GH");
    return;
  }


  double eps(0.00001);
  TH1D *hresult = new TH1D("hresult", Form("%s", era.c_str()), 80, 0., 0.4);
  hresult->GetXaxis()->SetTitle("BDT > ");
  hresult->SetMarkerStyle(24); hresult->SetLineColor(kRed+2); hresult->SetMarkerColor(kRed+2);
  TH1D *hsresult = new TH1D("hsresult", Form("%s", era.c_str()), 80, 0.+eps, 0.4+eps);
  hsresult->SetMarkerStyle(25); hsresult->SetLineColor(kBlue);  hsresult->SetMarkerColor(kBlue);

  TH1D *hbdtcut = new TH1D("hbdtcut", Form("%s", era.c_str()), 80, 0., 0.4);
  hbdtcut->SetMarkerColor(kRed+2);  hbdtcut->SetMarkerStyle(20);
  TH1D *hsbdtcut = new TH1D("hsbdtcut", Form("%s", era.c_str()), 80, 0.+eps, 0.4+eps);
  hsbdtcut->SetMarkerColor(kBlue);  hsbdtcut->SetMarkerStyle(21);

  invertedMuonIDRaw(era, 0.15, hresult, hsresult);
  invertedMuonIDRaw(era, 0.20, hresult, hsresult);
  invertedMuonIDRaw(era, 0.25, hresult, hsresult);
  invertedMuonIDRaw(era, 0.30, hresult, hsresult);
  invertedMuonIDRaw(era, -99., hbdtcut, hsbdtcut);

  c0->Clear();
  double ymax(50.);
  double smax = hresult->GetMaximum();
  double dmax = hsresult->GetMaximum();
  double themax = (smax>dmax? smax: dmax);
  if (themax < 50.)  ymax = 50.;
  if (themax < 20.)  ymax = 30.;
  if (themax < 5.)  ymax = 10.;
  if (themax < 0.5) ymax = 1.0;
  hresult->SetMinimum(-ymax);
  hresult->SetMaximum(ymax);
  hresult->Draw();
  hsresult->Draw("same");
  hbdtcut->Draw("same");
  hsbdtcut->Draw("same");
  pl->DrawLine(0., 0., 0.4, 0.);

  TLegend *tle = new TLegend(0.25, 0.2, 0.50, 0.5);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->SetHeader(era.c_str());
  tle->AddEntry(hresult, "yield at B0", "p");
  tle->AddEntry(hsresult, "yield at Bs", "p");
  tle->Draw();

  c0->SaveAs(Form("invertedMuonIDRaw-%s-result.pdf", era.c_str()));
}


// ----------------------------------------------------------------------
void plotRatio2Comb(string era = "2016BF") {
  gStyle->SetOptStat(0);

  if (era == "all") {
    plotRatio2Comb("2011");
    plotRatio2Comb("2012");
    plotRatio2Comb("2016BF");
    plotRatio2Comb("2016GH");
    return;
  }

  double eps(0.00001);
  TH1D *hresultBd = new TH1D("hresultBd", era.c_str(), 80, 0., 0.4);
  TH1D *hresultBs = new TH1D("hresultBs", era.c_str(), 80, 0.+eps, 0.4+eps);
  hresultBd->GetYaxis()->SetTitle("a.u.");
  hresultBd->GetXaxis()->SetTitle("BDT > "); hresultBd->SetMarkerStyle(24); hresultBd->SetLineColor(kRed+2); hresultBd->SetMarkerColor(kRed+2);
  hresultBs->GetXaxis()->SetTitle("BDT > "); hresultBs->SetMarkerStyle(25); hresultBs->SetLineColor(kBlue);  hresultBs->SetMarkerColor(kBlue);
  TH1D *hbdtcutBd = new TH1D("hbdtcutbd", "", 80, 0., 0.4);
  TH1D *hbdtcutBs = new TH1D("hbdtcutBs", "", 80, 0.+eps, 0.4+eps);
  hbdtcutBs->SetMarkerColor(kBlue);  hbdtcutBd->SetMarkerStyle(24); hbdtcutBs->SetLineColor(kBlue);
  hbdtcutBd->SetMarkerColor(kRed+2); hbdtcutBd->SetMarkerStyle(25); hbdtcutBd->SetLineColor(kRed+2);

  invertedMuonID2Comb(era, 0.15, hresultBd, hresultBs);
  invertedMuonID2Comb(era, 0.20, hresultBd, hresultBs);
  invertedMuonID2Comb(era, 0.25, hresultBd, hresultBs);
  invertedMuonID2Comb(era, 0.30, hresultBd, hresultBs);
  invertedMuonID2Comb(era, -99., hbdtcutBd, hbdtcutBs);

  c0->Clear();
  double ymax(5.);
  hresultBd->SetMinimum(-ymax);
  hresultBd->SetMaximum(ymax);
  hresultBd->Draw();
  hbdtcutBd->Draw("same");
  pl->DrawLine(0., 0., 0.4, 0.);

  tl->SetTextSize(0.035);
  tl->SetTextColor(kRed+2); tl->DrawLatexNDC(0.2, 0.85, "n(gauss(5.28))/combinatorial");
  tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.2, 0.81, "n(gauss(5.37))/combinatorial");
  tl->SetTextColor(kBlack);

  hresultBs->Draw("same");
  hbdtcutBs->Draw("same");

  c0->SaveAs(Form("invertedMuonID2Comb-%s-result.pdf", era.c_str()));

}


// ----------------------------------------------------------------------
void fitCombinatorial(string era = "2016BF", float bdt = -99.) {
  if ("all" == era) {
    fitCombinatorial("2016BF");
    fitCombinatorial("2016GH");
    fitCombinatorial("2012");
    return;
  }

  if ("scan2012" == era) {
    fitCombinatorial("2012", 0.30);
    fitCombinatorial("2012", 0.25);
    fitCombinatorial("2012", 0.20);
    fitCombinatorial("2012", 0.15);
    return;
  }


  TFile *fd(0), *fm0(0), *fm1(0);
  float bdtcut(0.);
  cout << "======================================================================" << endl;
  if (era == "2016BF") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016BF-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.30;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2016GH") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016GH-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.31;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2012") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2012-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.34;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2011") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.28;
    } else {
      bdtcut = bdt;
    }
  }

  TTree *td = (TTree*)fd->Get("candAnaMuMu/events");

  TH1D *hd0 = new TH1D("hd0", "mass", 60, 4.8, 6.0);
  string cut(Form("chan==0 && !gmuid &&hlt1&&tos&&l1t &&bdt>%4.3f && m < 6 && m > 4.8", bdtcut));
  cout << cut << endl;
  td->Draw("m>>hd0", cut.c_str());

  TH1D *hd1 = (TH1D*)hd0->Clone("hd1");
  TH1D *hd2 = (TH1D*)hd0->Clone("hd2");

  TCanvas *c1 = new TCanvas("c1", "", 900, 300);
  c1->Modified();
  c1->Update();
  c1->Divide(3,1);

  c1->cd(1);
  hd0->Fit("pol0", "rl", "", 5.5, 6.0);
  float hd0chi2 = hd0->GetFunction("pol0")->GetChisquare();
  tl->SetTextSize(0.06);
  tl->DrawLatexNDC(0.4, 0.92, Form("#chi^{2} = %3.2f", hd0chi2));


  c1->cd(2);
  hd1->Fit("pol1", "rl", "", 5.5, 6.0);
  float hd1chi2 = hd1->GetFunction("pol1")->GetChisquare();
  tl->DrawLatexNDC(0.4, 0.92, Form("#chi^{2} = %3.2f", hd1chi2));

  c1->cd(3);
  hd2->Fit("expo", "rl", "", 5.5, 6.0);
  float hd2chi2 = hd2->GetFunction("expo")->GetChisquare();
  tl->DrawLatexNDC(0.4, 0.92, Form("#chi^{2} = %3.2f", hd2chi2));

  c1->SaveAs(Form("bg-fitting-%s-%3.2f.pdf", era.c_str(), bdtcut));
}
