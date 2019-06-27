// ----------------------------------------------------------------------
// -- plots 'truth-identified cands' and applies rc or !rc (whether or not a reco candidate was present as well)
void plot3(string what = "pt", double xmin = 0., double xmax = 10., string cuts = "missedtracks==0",
	   string filename = "recoil-mc-2017-recoil12-v01-0000.BuToMuTauK:BuToMuTauKTruth.root") {
  TFile *f = TFile::Open(filename.c_str());
  TTree *t = (TTree*)f->Get("candAnaBuToMuTauKTruth/events");
  TTree *r = (TTree*)f->Get("candAnaBuToMuTauK/events");

  int NBINS(50);

  // -- require reconstructed signal cand, same PV for recoil and signal, and tracks in (or not in) recoil
  TH1D *h1(0), *h2(0), *h3(0), *hr(0);
  hr = new TH1D("hr", "", NBINS, xmin, xmax); setFilledHist(hr, kBlack, kYellow, 1000);
  setTitles(hr, what.c_str(), "Entries", 0.05, 1.0, 1.5);
  h1 = new TH1D("h1", "", NBINS, xmin, xmax); setFilledHist(h1, kBlack, kBlack, 3365);
  h2 = new TH1D("h2", "", NBINS, xmin, xmax); setFilledHist(h2, kBlue,  kBlue,  3365);
  h3 = new TH1D("h3", "", NBINS, xmin, xmax); setFilledHist(h3, kBlack,   kBlack,   3354);

  gStyle->SetOptStat(0);
  shrinkPad(0.11, 0.15);
  r->Draw(Form("%s>>hr", what.c_str()), Form("%s ", cuts.c_str()), "goff");
  t->Draw(Form("%s>>h1", what.c_str()), Form("%s ", cuts.c_str()), "goff");
  t->Draw(Form("%s>>h3", what.c_str()), Form("%s &&!rm", cuts.c_str()), "goff");
  t->Draw(Form("%s>>h2", what.c_str()), Form("%s &&rm", cuts.c_str()), "goff");

  double rmax = hr->GetMaximum();
  double tmax = h1->GetMaximum();
  if (tmax > rmax) {
    hr->SetMaximum(tmax*1.2);
  }
  hr->SetMinimum(0.5);
  hr->Draw();
  h1->Draw("esame");
  h2->Draw("same");
  h3->Draw("same");

  tl->SetTextSize(0.015);tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.05, 0.96, Form("cuts: %s", cuts.c_str()));
  tl->SetTextSize(0.03);tl->SetTextColor(kYellow+2); tl->DrawLatexNDC(0.20, 0.92, Form("reco: %.0f", hr->GetEntries()));
  tl->SetTextSize(0.03);tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.52, 0.92, Form("TRUTH:"));
  tl->SetTextSize(0.03);tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.65, 0.92, Form("rm: %.0f", h2->GetEntries()));
  tl->SetTextSize(0.03);tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.80, 0.92, Form("!rm: %.0f", h3->GetEntries()));

  c0->SaveAs(Form("qnd/plot3-%s.pdf", what.c_str()));
}

// ----------------------------------------------------------------------
void plot3All(string cuts = "missedtracks==0",
	      string filename = "recoil-mc-2017-recoil12-v01-0000.BuToMuTauK:BuToMuTauKTruth.root") {

  plot3("oahad", 0., 1.1, cuts, filename);
  plot3("abtau", 0., 3.3, cuts, filename);
  plot3("atauhad", 0., 3.3, cuts, filename);

  plot3("t0had", 0., 3.e-12, cuts, filename);
  plot3("t0muka", 0., 10.e-12, cuts, filename);

  plot3("mmukahad", 0., 7.0, cuts, filename);
  plot3("mmuka", 0., 4.0, cuts, filename);
  plot3("flsmuka", 0., 60.0, cuts, filename);
  plot3("flmuka", 0., 1.0, cuts, filename);
  plot3("abmuka", 0., 3.3, cuts, filename);

  plot3("ptmuka", 0., 30.0, cuts, filename);
  plot3("ptmu", 0., 20.0, cuts, filename);
  plot3("ptka", 0., 15.0, cuts, filename);
  plot3("pthad1", 0., 15.0, cuts, filename);
  plot3("pthad2", 0., 10.0, cuts, filename);
  plot3("pthad3", 0.,  8.0, cuts, filename);

  plot3("mhad", 0., 1.8, cuts, filename);
  plot3("liphad", 0., 0.03, cuts, filename);
  plot3("pthad", 0., 25., cuts, filename);
  plot3("flshad", 0., 25., cuts, filename);
  plot3("flhad", 0., 0.2, cuts, filename);

  plot3("mbr0pos", 0., 20, cuts, filename);
  plot3("mbr0neg", 0., 20, cuts, filename);

  plot3("doca3dmuka", 0., 0.025, cuts, filename);
  plot3("docamaxmuka", 0., 0.025, cuts, filename);
  plot3("lipmuka", 0., 0.05, cuts, filename);

}


// ----------------------------------------------------------------------
// -- plots 'data cands' and applies tm or !tm
void plot2(string what = "pt", double xmin = 0., double xmax = 10., string cuts = "ptb>0",
	   string filename = "bac-dbx-recoil-mc-2017-recoil12-v01-0000.BuToMuTauK.root") {
  TFile *f = TFile::Open(filename.c_str());
  TTree *t = (TTree*)f->Get("candAnaBuToMuTauK/events");

  int NBINS(50);

  // -- require reconstructed signal cand, same PV for recoil and signal, and tracks in (or not in) recoil
  TH1D *h1(0), *h2(0), *h3(0);
  h1 = new TH1D("h1", what.c_str(), NBINS, xmin, xmax);   setFilledHist(h1, kBlack, kBlack, 3365);
  h2 = new TH1D("h2", "tm", NBINS, xmin, xmax);    setFilledHist(h2, kBlue,  kBlue,  3365);
  h3 = new TH1D("h3", "!tmh3", NBINS, xmin, xmax); setFilledHist(h3, kRed,   kRed,   3354);

  gStyle->SetOptStat(0);
  t->Draw(Form("%s>>h1", what.c_str()), Form("%s ", cuts.c_str()), "e");
  t->Draw(Form("%s>>h3", what.c_str()), Form("%s &&!tm", cuts.c_str()), "same");
  t->Draw(Form("%s>>h2", what.c_str()), Form("%s &&tm", cuts.c_str()), "same");
  tl->SetTextSize(0.03);tl->SetTextColor(kBlack); tl->DrawLatexNDC(0.04, 0.92, cuts.c_str());
  tl->SetTextSize(0.03);tl->SetTextColor(kBlue);  tl->DrawLatexNDC(0.80, 0.96, Form("TM: %.0f", h2->GetEntries()));
  tl->SetTextSize(0.03);tl->SetTextColor(kRed);   tl->DrawLatexNDC(0.80, 0.92, Form("!TM: %.0f", h3->GetEntries()));

  c0->SaveAs(Form("plot2-%s.pdf", what.c_str()));
}

// ----------------------------------------------------------------------
void plot2All(string filename = "bac-dbx-recoil-mc-2017-recoil12-v01-0000.BuToMuTauK.root") {
  plot2("doca3dmuka", 0., 0.025, "1");
  plot2("docamaxmuka", 0., 0.025, "1");
  plot2("lipmuka", 0., 0.05, "1");
  plot2("mhad", 0., 1.6, "1");

}




// ----------------------------------------------------------------------
// -- plot cand tracks character
void plot1(string what = "pt", string filename = "dbx-recoil-mc-2017-recoil12-v01-0000.BuToMuTauK.root") {
  TFile *f = TFile::Open(filename.c_str());
  TTree *t = (TTree*)f->Get("candAnaBuToMuTauK/events");

  zone(2,3);
  // -- require reconstructed signal cand, same PV for recoil and signal, and tracks in (or not in) recoil
  TH1D *h1(0), *h2(0);
  if (what == "pt") {
    h1 = new TH1D("h1", "h1", 100, 0., 20.); setFilledHist(h1, kBlue, kBlue, 3056);
    h2 = new TH1D("h2", "h2", 100, 0., 20.); setFilledHist(h2, kRed, kRed, 3065);
  } else if (what == "doca") {
    h1 = new TH1D("h1", "h1", 40, 0., 0.2); setFilledHist(h1, kBlue, kBlue, 3056);
    h2 = new TH1D("h2", "h2", 40, 0., 0.2); setFilledHist(h2, kRed, kRed, 3065);
  }

  for (int is = 0; is < 5; ++is) {
    c0->cd(is+1);
    t->Draw(Form("%s[%d]>>h1", what.c_str(), is), Form("ptb>0 && pvr==pvc && inrecoil[%d]", is), "goff");
    t->Draw(Form("%s[%d]>>h2", what.c_str(), is), Form("ptb>0 && pvr==pvc && !inrecoil[%d]", is), "goff");
    h1->DrawCopy();
    h2->DrawCopy("same");
  }
  c0->cd(6);
  tl->SetTextSize(0.12);
  tl->DrawLatexNDC(0.5, 0.5, what.c_str());

  c0->SaveAs(Form("plot1-%s.pdf", what.c_str()));
}
