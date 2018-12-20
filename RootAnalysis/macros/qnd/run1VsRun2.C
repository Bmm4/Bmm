// ----------------------------------------------------------------------
struct data {
  int    run, evt, ls, chan, pvn;
  bool   muid, ok;
  double m1bdt, m2bdt;
  double m1eta, m2eta;
  double m1pt, m2pt;
  double m, bdt;

  void clear() {
    run = evt = ls = chan = pvn = -99;
    muid = ok = false;
    m1bdt = m2bdt = m1eta = m2eta = m1pt = m2pt = m = bdt = -99.;
  }
};

// ----------------------------------------------------------------------
void run1VsRun2(int year = 2012, double bdtmin = 0.36, double bdtmax = 99.) {

  TFile *f1(0), *f2(0);
  if (2012 == year) {
    f1 = TFile::Open("/t3home/ursl/oldhome/bmm3/ana/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/macros2/2012/small-SgData.root");
    f2 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2012-s01.root");
  } else if (2011 == year) {
    f1 = TFile::Open("/t3home/ursl/oldhome/bmm3/ana/CMSSW_5_3_6/src/HeavyFlavorAnalysis/Bs2MuMu/macros2/2011/small-SgData.root");
    f2 = TFile::Open("/scratch/ursl/bmm4/s01/bmm-data-bmmMuOnia2011-s01.root");
  }


  TTree *t1 = (TTree*)f1->Get("SgData_bdt");
  TTree *t2 = (TTree*)f2->Get("candAnaMuMu/events");

  int chan;
  int t1run, t1evt, t1ls, t2ls, t1pvn, t2pvn;
  Long64_t t2run, t2evt;
  double t1bdt, t2bdt;
  double t1m, t2m;
  double t1m1eta, t2m1eta;
  double t1m2eta, t2m2eta;
  double t1m1pt, t2m1pt;
  double t1m2pt, t2m2pt;
  double t1m1bdt, t2m1bdt;
  double t1m2bdt, t2m2bdt;
  bool   t2hlt1, t2tos, t2l1t;
  bool   t1muid, t2gmuid;

  t1->SetBranchAddress("run", &t1run);
  t2->SetBranchAddress("run", &t2run);

  t1->SetBranchAddress("evt", &t1evt);
  t2->SetBranchAddress("evt", &t2evt);

  t1->SetBranchAddress("ls", &t1ls);
  t2->SetBranchAddress("ls", &t2ls);

  t1->SetBranchAddress("pvn", &t1pvn);
  t2->SetBranchAddress("pvn", &t2pvn);

  t1->SetBranchAddress("bdt", &t1bdt);
  t2->SetBranchAddress("bdt", &t2bdt);

  t1->SetBranchAddress("m", &t1m);
  t2->SetBranchAddress("m", &t2m);

  t1->SetBranchAddress("muid", &t1muid);
  t2->SetBranchAddress("gmuid", &t2gmuid);

  t1->SetBranchAddress("m1eta", &t1m1eta);
  t2->SetBranchAddress("m1eta", &t2m1eta);

  t1->SetBranchAddress("m2eta", &t1m2eta);
  t2->SetBranchAddress("m2eta", &t2m2eta);

  t1->SetBranchAddress("m1pt", &t1m1pt);
  t2->SetBranchAddress("m1pt", &t2m1pt);

  t1->SetBranchAddress("m2pt", &t1m2pt);
  t2->SetBranchAddress("m2pt", &t2m2pt);

  t1->SetBranchAddress("m1bdt", &t1m1bdt);
  t2->SetBranchAddress("m1mvabdt", &t2m1bdt);

  t1->SetBranchAddress("m2bdt", &t1m2bdt);
  t2->SetBranchAddress("m2mvabdt", &t2m2bdt);

  t2->SetBranchAddress("hlt1", &t2hlt1);
  t2->SetBranchAddress("tos", &t2tos);
  t2->SetBranchAddress("l1t", &t2l1t);

  data a;
  vector<data> t1data;
  for (int i = 0; i < t1->GetEntries(); ++i) {
    a.clear();
    t1->GetEntry(i);
    // if (t1m1bdt < 0.2) continue;
    // if (t1m2bdt < 0.2) continue;
    if (t1m < 4.9) continue;
    if (t1m > 5.9) continue;
    if (t1m > 5.2 && t1m < 5.45) continue;
    chan = -1;
    if ((TMath::Abs(t1m1eta) < 1.4) && (TMath::Abs(t1m2eta) < 1.4)) {
      chan = 0;
    } else {
      if ((TMath::Abs(t1m1eta) < 2.4) && (TMath::Abs(t1m2eta) < 2.4)) {
	chan = 1;
      }
    }

    if (t1bdt < bdtmin) continue;
    if (t1bdt > bdtmax) continue;
    a.run   = t1run;
    a.evt   = t1evt;
    a.ls    = t1ls;
    a.pvn   = t1pvn;
    a.chan  = chan;
    a.bdt   = t1bdt;
    a.m1eta = t1m1eta;
    a.m2eta = t1m2eta;
    a.m1pt  = t1m1pt;
    a.m2pt  = t1m2pt;
    a.m1bdt = t1m1bdt;
    a.m2bdt = t1m2bdt;
    a.muid  = t1muid;
    a.ok    = true;
    a.m     = t1m;
    t1data.push_back(a);
  }

  cout << "BMM3: n = " << t1data.size() << endl;
  for (unsigned int i = 0; i < t1data.size(); ++i) {
    cout << t1data[i].run << " " << t1data[i].evt << " " << t1data[i].m << " m1bdt = " << t1data[i].m1bdt << " m2bdt = " << t1data[i].m2bdt << endl;
  }

  data b;
  vector<data> t2data;
  for (int i = 0; i < t2->GetEntries(); ++i) {
    b.clear();
    t2->GetEntry(i);
    // if (t2m1bdt < 0.5) continue;
    // if (t2m2bdt < 0.5) continue;
    if (t2m < 4.9) continue;
    if (t2m > 5.9) continue;
    if (t2m > 5.2 && t2m < 5.45) continue;
    chan = -1;
    if ((TMath::Abs(t2m1eta) < 1.4) && (TMath::Abs(t2m2eta) < 1.4)) {
      chan = 0;
    } else {
      if ((TMath::Abs(t2m1eta) < 2.4) && (TMath::Abs(t2m2eta) < 2.4)) {
	chan = 1;
      }
    }
    if (t2bdt < bdtmin) continue;
    if (t2bdt > bdtmax) continue;
    if (!t2hlt1) continue;
    if (!t2tos) continue;
    if (!t2l1t) continue;
    if (t2gmuid && t2bdt > 0.21) {
      cout << t2run << "/" << t2evt << " |m1eta| = " << TMath::Abs(t2m1eta) << " |m2eta| = " << TMath::Abs(t2m2eta) << " chan = " << chan << endl;
    }
    b.run   = t2run;
    b.evt   = t2evt;
    b.ls    = t2ls;
    b.pvn   = t2pvn;
    b.chan  = chan;
    b.bdt   = t2bdt;
    b.m1eta = t2m1eta;
    b.m2eta = t2m2eta;
    b.m1pt  = t2m1pt;
    b.m2pt  = t2m2pt;
    b.m1bdt = t2m1bdt;
    b.m2bdt = t2m2bdt;
    b.muid  = t2gmuid;
    b.ok    = true;
    b.m     = t2m;
    t2data.push_back(b);
  }


  cout << "BMM4: " << endl;
  for (unsigned int i = 0; i < t2data.size(); ++i) {
    cout << t2data[i].run << " " << t2data[i].evt << " " << t2data[i].m << " m1bdt = " << t2data[i].m1bdt << " m2bdt = " << t2data[i].m2bdt << endl;
  }



  // -- dump the data into a combined tree
  TFile *fnew = TFile::Open(Form("combined-%d.root", year), "RECREATE");
  TTree *tnew = new TTree("combined", "combined");
  tnew->Branch("run", &a.run, "run/I");
  tnew->Branch("evt", &a.evt, "evt/I");
  tnew->Branch("ls",  &a.ls,  "ls/I");
  tnew->Branch("bmm3chan",  &a.chan,  "bmm3chan/I");
  tnew->Branch("bmm3pvn",   &a.pvn,   "bmm3pvn/I");
  tnew->Branch("bmm3bdt",   &a.bdt,   "bmm3bdt/D");
  tnew->Branch("bmm3m1bdt", &a.m1bdt, "bmm3m1bdt/D");
  tnew->Branch("bmm3m2bdt", &a.m2bdt, "bmm3m2bdt/D");
  tnew->Branch("bmm3muid",  &a.muid,  "bmm3muid/O");
  tnew->Branch("bmm3m1eta", &a.m1eta, "bmm3m1eta/D");
  tnew->Branch("bmm3m2eta", &a.m2eta, "bmm3m2eta/D");
  tnew->Branch("bmm3m1pt",  &a.m1pt,  "bmm3m1pt/D");
  tnew->Branch("bmm3m2pt",  &a.m2pt,  "bmm3m2pt/D");
  tnew->Branch("bmm3m",     &a.m,     "bmm3m/D");
  tnew->Branch("bmm3ok",    &a.ok,     "bmm3m/O");

  tnew->Branch("bmm4chan",  &b.chan,  "bmm4chan/I");
  tnew->Branch("bmm4pvn",   &b.pvn,   "bmm4pvn/I");
  tnew->Branch("bmm4bdt",   &b.bdt,   "bmm4bdt/D");
  tnew->Branch("bmm4m1bdt", &b.m1bdt, "bmm4m1bdt/D");
  tnew->Branch("bmm4m2bdt", &b.m2bdt, "bmm4m2bdt/D");
  tnew->Branch("bmm4muid",  &b.muid,  "bmm4muid/O");
  tnew->Branch("bmm4m1eta", &b.m1eta, "bmm4m1eta/D");
  tnew->Branch("bmm4m2eta", &b.m2eta, "bmm4m2eta/D");
  tnew->Branch("bmm4m1pt",  &b.m1pt,  "bmm4m1pt/D");
  tnew->Branch("bmm4m2pt",  &b.m2pt,  "bmm4m2pt/D");
  tnew->Branch("bmm4m",     &b.m,     "bmm4m/D");
  tnew->Branch("bmm4ok",    &b.ok,    "bmm4m/O");

  // -- Fill BMM3 and overlap
  for (unsigned int i = 0; i < t1data.size(); ++i) {
    a = t1data[i];
    b.clear();
    for (unsigned int j = 0; j < t2data.size(); ++j) {
      if (t2data[j].run == a.run && t2data[j].evt == a.evt) {
	b = t2data[j];
	break;
      }
    }
    tnew->Fill();
  }

  // -- Fill non-overlapping BMM4
  a.clear();
  for (unsigned int j = 0; j < t2data.size(); ++j) {
    b = t2data[j];
    bool already(false);
    for (unsigned int i = 0; i < t1data.size(); ++i) {
      if (t1data[i].run == b.run && t1data[i].evt == b.evt) {
	already = true;
	break;
      }
    }
    if (already) {
      // skip, do nothing, already in tree
    } else {
      a.run = b.run;
      a.evt = b.evt;
      a.ls  = b.ls;
      tnew->Fill();
    }
  }


  fnew->Write();
  fnew->Save();
  fnew->Close();

  return;
}



// ----------------------------------------------------------------------
void plot(string var = "m", string sideband = "lo", int chan = 0, int year = 2012) {
  bool usebdt4(true);

  TFile *f1(0);
  if (2012 == year) {
    f1 = TFile::Open("combined-2012.root");
  } else if (2011 == year) {
    f1 = TFile::Open("combined-2011.root");
  }


  TTree *t1 = (TTree*)f1->Get("combined");


  TH1D *h1(0), *h2(0), *h0(0), *hd(0);
  TH2D *hc(0);
  if ("m" == var) {
    h1 = new TH1D("run1", Form("BMM3 %s sideband %d", sideband.c_str(), year), 40, 4.9, 5.9);
    h2 = new TH1D("run2", Form("BMM4 %s sideband %d", sideband.c_str(), year), 40, 4.9, 5.9);
    h0 = new TH1D("common", Form("overlap %s sideband %d", sideband.c_str(), year), 40, 4.9, 5.9);
    hd = new TH1D("diff", Form("overlap %s sideband %d", sideband.c_str(), year), 40, -0.1, 0.1);
    hc = new TH2D("corr", Form("overlap %s sideband %d", sideband.c_str(), year), 40, 4.9, 5.9, 40, 4.9, 5.9);
  }

  if ("bdt" == var) {
    h1 = new TH1D("run1", Form("BMM3 %s sideband %d", sideband.c_str(), year), 60, 0.2, 0.8);
    h2 = new TH1D("run2", Form("BMM4 %s sideband %d", sideband.c_str(), year), 60, 0.2, 0.8);
    h0 = new TH1D("common", Form("overlap %s sideband %d", sideband.c_str(), year), 60, 0.2, 0.8);
    hd = new TH1D("diff", Form("overlap %s sideband %d", sideband.c_str(), year), 40, -0.2, 0.2);
    hc = new TH2D("corr", Form("overlap %s sideband %d", sideband.c_str(), year), 60, 0.2, 0.8, 40, 0.2, 0.8);
  }

  if ("m1eta" == var) {
    h1 = new TH1D("run1", Form("BMM3 %s sideband %d", sideband.c_str(), year), 42, -2.1, 2.1);
    h2 = new TH1D("run2", Form("BMM4 %s sideband %d", sideband.c_str(), year), 42, -2.1, 2.1);
    h0 = new TH1D("common", Form("overlap %s sideband %d", sideband.c_str(), year), 42, -2.1, 2.1);
    hd = new TH1D("diff", Form("overlap %s sideband %d", sideband.c_str(), year), 40, -0.1, 0.1);
    hc = new TH2D("corr", Form("overlap %s sideband %d", sideband.c_str(), year), 42, -2.1, 2.1, 42, -2.1, 2.1);
  }
  h1->SetXTitle(var.c_str()); h1->SetYTitle("events/bin");
  h2->SetXTitle(var.c_str()); h2->SetYTitle("events/bin");
  h0->SetXTitle(var.c_str()); h0->SetYTitle("events/bin");
  hd->SetXTitle(Form("%s(bmm3) - %s(bmm4)", var.c_str(), var.c_str()));
  hc->SetXTitle(Form("%s (bmm4)", var.c_str())); hc->SetYTitle(Form("%s (bmm3)", var.c_str()));

  bool   ok3, ok4;
  bool   muid3, muid4;
  double m3, m4, bdt3, bdt4;
  int    chan3, chan4;
  double var3(0.), var4(0.);
  t1->SetBranchAddress("bmm3ok", &ok3);
  t1->SetBranchAddress("bmm4ok", &ok4);
  t1->SetBranchAddress("bmm3muid", &muid3);
  t1->SetBranchAddress("bmm4muid", &muid4);
  t1->SetBranchAddress("bmm3m", &m3);
  t1->SetBranchAddress("bmm4m", &m4);
  t1->SetBranchAddress("bmm3bdt", &bdt3);
  t1->SetBranchAddress("bmm4bdt", &bdt4);
  t1->SetBranchAddress("bmm3chan", &chan3);
  t1->SetBranchAddress("bmm4chan", &chan4);
  if (("m" != var) && ("bdt" != var)) {
    t1->SetBranchAddress(Form("bmm3%s", var.c_str()),  &var3);
    t1->SetBranchAddress(Form("bmm4%s", var.c_str()),  &var4);
  }

  double bdtcut(99.);
  double bdt4cut(99.);
  if (2012 == year) {
    if (0 == chan) {
      bdtcut = 0.36;
      if (usebdt4) {
	bdt4cut = 0.34;
      } else {
	bdt4cut = bdtcut;
      }
    } else if (1 == chan) {
      bdtcut = 0.38;
      if (usebdt4) {
	bdt4cut = 0.32;
      } else {
	bdt4cut = bdtcut;
      }
    } else {
      bdtcut = 0.38;
      if (usebdt4) {
	bdt4cut = 0.34;
      } else {
	bdt4cut = bdtcut;
      }
    }
  } else if (2011 == year) {
    if (0 == chan) {
      bdtcut = 0.29;
      if (usebdt4) {
	bdt4cut = 0.28;
      } else {
	bdt4cut = bdtcut;
      }
    } else if (1 == chan) {
      bdtcut = 0.29;
      if (usebdt4) {
	bdt4cut = 0.21;
      } else {
	bdt4cut = bdtcut;
      }
    } else {
      bdtcut = 0.29;
      if (usebdt4) {
	bdt4cut = 0.28;
      } else {
	bdt4cut = bdtcut;
      }
    }
  }


  for (int i = 0; i < t1->GetEntries(); ++i) {
    t1->GetEntry(i);
    // cout << "sideband = " << sideband << " chan = " << chan << " masses = " << m3 << " " << m4 << endl;
    if (("lo" == sideband) && (((m3 > -1.) && (m3 > 5.3)) || ((m4 > -1.) && (m4 > 5.3)))) continue;
    if (("hi" == sideband) && (((m3 > -1.) && (m3 < 5.3)) || ((m4 > -1.) && (m4 < 5.3)))) continue;
    if ((0 == chan) && (((chan3 > -1) && (chan3 != 0)) || ((chan4 > -1) && (chan4 != 0)))) continue;
    if ((1 == chan) && (((chan3 > -1) && (chan3 != 1)) || ((chan4 > -1) && (chan4 != 1)))) continue;
    if ("m" == var) {
      var3 = m3;
      var4 = m4;
    } else if ("bdt" == var) {
      var3 = bdt3;
      var4 = bdt4;
    }
    // cout << "  cont, ok = " << ok3 << " " << ok4 << " vars: " << var3 << " " << var4 << endl;

    //    if (ok3 && ok4 && muid3 && muid4 && (bdt3 > bdtcut) && (bdt4 > bdtcut)) {
    if (ok3 && ok4 && muid3 && muid4 && (bdt3 > bdtcut) && (bdt4 > bdt4cut)) {
      //      cout << " vars: " << var3 << " "  << var4 << " bdts: " << bdt3 << " " << bdt4 << endl;
      h0->Fill(var4);
      hc->Fill(var3, var4);
      hd->Fill(var3 - var4);
    } else if (ok3 && muid3 && (bdt3 > bdtcut)) {
	h1->Fill(var3);
	//    } else if (ok4 && muid4 && (bdt4 > bdtcut)) {
    } else if (ok4 && muid4 && (bdt4 > bdt4cut)) {
      //      cout << " var4: " << var4 << " bdt4: " << bdt4 << " " << endl;
      h2->Fill(var4);
    }
  }

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 250);
  c1->Clear();
  c1->Divide(4,1);

  c1->cd(1);
  h1->Draw();
  tl->DrawLatexNDC(0.7, 0.92, Form("N = %d", static_cast<int>(h1->Integral())));
  c1->cd(2);
  h2->Draw();
  tl->DrawLatexNDC(0.7, 0.92, Form("N = %d", static_cast<int>(h2->Integral())));

  c1->cd(3);
  h0->Draw();
  tl->DrawLatexNDC(0.7, 0.92, Form("N = %d", static_cast<int>(h0->Integral())));

  c1->cd(4);
  //  hc->Draw("colz");
  hd->Draw();
  tl->DrawLatexNDC(0.6, 0.92, Form("RMS = %4.3f", hd->GetRMS()));
  string fext("");
  if (usebdt4) fext = "Bmm4Cuts";
  if (var == "m") {
    ofstream TEX(Form("syncAnaVenn%s.tex", fext.c_str()), ios::app);
    TEX << Form("\\vdef{%sbmm3:sb%s:chan%d:year%d}  {%d}", fext.c_str(), sideband.c_str(), chan, year, static_cast<int>(h1->Integral())) << endl;
    TEX << Form("\\vdef{%sbmm4:sb%s:chan%d:year%d}  {%d}", fext.c_str(), sideband.c_str(), chan, year, static_cast<int>(h2->Integral())) << endl;
    TEX << Form("\\vdef{%sbmmX:sb%s:chan%d:year%d}  {%d}", fext.c_str(), sideband.c_str(), chan, year, static_cast<int>(h0->Integral())) << endl;
  }

  c1->SaveAs(Form("syncAna%s-%s-%s-chan%d-%d.pdf", fext.c_str(), var.c_str(), sideband.c_str(), chan, year));
}


// ----------------------------------------------------------------------
void plotAll(int year) {
  vector<string> sb; sb.push_back("lo"); sb.push_back("hi");
  vector<int>    chan; chan.push_back(0); chan.push_back(1);
  vector<string> vars; vars.push_back("m"); vars.push_back("bdt"); vars.push_back("m1eta"); //vars.push_back("m2bdt");

  for (unsigned isb = 0; isb < sb.size(); ++isb) {
    for (unsigned ichan = 0; ichan  < chan.size(); ++ichan) {
      for (unsigned ivar = 0; ivar  < vars.size(); ++ivar) {
	plot(vars[ivar], sb[isb], chan[ichan], year);
      }
    }
  }
}
