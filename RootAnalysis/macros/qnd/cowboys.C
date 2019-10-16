#include "../common/util.hh"

// ----------------------------------------------------------------------
void cowboys(string era = "2016BF", int ichan = 0, string cuts = "gmuid &&hlt1&&tos&&l1t", double bdt = 0.1) {

  TFile *fd(0), *fm0(0), *fm1(0);
  float bdtcut(0.);
  cout << "======================================================================" << endl;
  if (era == "2016BF") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016BF-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-combined-BsToMuMu-2016BF-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.30;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2016GH") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmCharmonium2016GH-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-combined-BsToMuMu-2016GH-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.31;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2012") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2012-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BsToMuMu-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.34;
    } else {
      bdtcut = bdt;
    }
  } else if (era == "2011") {
    fd = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root");
    fm0 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BsToMuMu-s01.root");
    if (bdt < -1.) {
      bdtcut = 0.28;
    } else {
      bdtcut = bdt;
    }
  }

  TTree *td = (TTree*)fd->Get("candAnaMuMu/events");
  TTree *tm0 = (TTree*)fm0->Get("candAnaMuMu/events");
  //  TTree *tm1 = (TTree*)fm1->Get("candAnaMuMu/events");

  TH1D *hd0  = new TH1D("hd0 ", "", 20, 0., 2.0); hd0->Sumw2(); setHist(hd0, kBlack, 24); hd0->SetMinimum(0.);
  setTitles(hd0, "#Delta R", "cowboys/seagulls");
  TH1D *hd0c = new TH1D("hd0c", "dr", 20, 0., 2.0);
  TH1D *hd0s = new TH1D("hd0s", "dr", 20, 0., 2.0);

  double eps(0.01);
  TH1D *hm0  = new TH1D("hm0", "", 20, 0.+eps, 2.0+eps); hm0->Sumw2(); setHist(hm0, kRed, 24);
  TH1D *hm0c = new TH1D("hm0c", "dr", 20, 0.+eps, 2.0+eps);
  TH1D *hm0s = new TH1D("hm0s", "dr", 20, 0.+eps, 2.0+eps);
  //  TH1D *hm1 = new TH1D("hm1", "dr", 40, 0., 2.0);
  string cut(Form("chan==%d && %s &&bdt>%4.3f && m < 6 && m > 4.8", ichan, cuts.c_str(), bdtcut));
  cout << cut << endl;

  string ccut(Form("cb  && %s", cut.c_str()));
  string scut(Form("!cb && %s", cut.c_str()));
  td->Draw("dr>>hd0c", ccut.c_str());
  td->Draw("dr>>hd0s", scut.c_str());

  tm0->Draw("dr>>hm0c", ccut.c_str());
  tm0->Draw("dr>>hm0s", scut.c_str());

  hd0->Divide(hd0c, hd0s, 1., 1.); hd0->SetMinimum(0.);
  hm0->Divide(hm0c, hm0s, 1., 1.);



  hd0->Draw();
  hm0->Draw("same");

  tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.6, 0.92, Form("%s/chan%d", era.c_str(), ichan));
  tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.2, 0.92, "Data");
  tl->SetTextColor(kRed);   tl->DrawLatexNDC(0.4, 0.92, "MC");
  tl->SetTextSize(0.02);
  tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.2, 0.87, cut.c_str());
  tl->SetTextSize(0.05);
  c0->SaveAs(Form("qnd/cowboys-%s-%d.pdf", era.c_str(), ichan));

}


// ----------------------------------------------------------------------
void plotAll() {

  cowboys("2011", 0);
  cowboys("2011", 1);

  cowboys("2012", 0);
  cowboys("2012", 1);

  cowboys("2016BF", 0);
  cowboys("2016BF", 1);

  cowboys("2016GH", 0);
  cowboys("2016GH", 1);

}
