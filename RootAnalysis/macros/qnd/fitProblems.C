// ----------------------------------------------------------------------
void fitProblems(int bdtCut = 20) {

  TFile *f = TFile::Open(Form("s01/plotResults.2016BFs01.root"));
  f->cd("bspsiphiData");
  string hname(Form("hNorm_cnc2016BFs01_bspsiphiData_chan0"));
  fitPsYield *a = new fitPsYield(hname, 0);
  a->fitBs2JpsiPhi(-1, "bla");
}

// ----------------------------------------------------------------------
void fitProblems12() {
  int imode(0);
  TFile *f = TFile::Open(Form("s01/plotResults.2012s01.root"));
  f->cd("bupsikData");
  //  string hname(Form("hW8Norm_bdt_19_2016BFs01_bupsikData_chan1"));
  string hname(Form("hW8Norm_bdt_41_2012s01_bupsikData_chan0"));
  fitPsYield fpy(hname);
  fpy.fitBu2JpsiKp(5, "qnd/norm", imode, 5.0, 5.8, 0.025);

}

// ----------------------------------------------------------------------
void fitProblems11() {
  int imode(0);
  TFile *f = TFile::Open(Form("s01/plotResults.2016BFs01.root"));
  f->cd("bupsikData");
  //  string hname(Form("hW8Norm_bdt_19_2016BFs01_bupsikData_chan1"));
  string hname(Form("hW8Norm_bdt_0_2016BFs01_bupsikData_chan1"));
  fitPsYield fpy(hname);
  fpy.fitBu2JpsiKp(5, "qnd/norm", imode, 5.0, 5.8, 0.025);

}

// ----------------------------------------------------------------------
void fitProblems10() {
  int imode(0);
  TFile *f = TFile::Open(Form("s01/plotResults.2016GHs01.root"));
  f->cd("bupsikData");
  string hname(Form("hW8Norm_bdt_22_2016GHs01_bupsikData_chan1"));
  fitPsYield fpy(hname);
  fpy.fitBu2JpsiKp(5, "qnd/norm", imode, 5.0, 5.8, 0.025);

}

// ----------------------------------------------------------------------
void fitProblems9() {
  int imode(0);
  TFile *f = TFile::Open(Form("s01/plotResults.2011s01.root"));
  f->cd("bupsikData");
  string hname(Form("hNorm_bdt2011s01_bupsikData_chan0"));
  fitPsYield fpy(hname);
  fpy.fitBu2JpsiKp(5, "qnd/norm", imode, 5.0, 5.8, 0.025);

  string hname2(Form("hW8Norm_bdt2011s01_bupsikData_chan0"));
  fitPsYield fpy2(hname2);
  fpy2.fitBu2JpsiKp(5, "qnd/norm", imode, 5.0, 5.8, 0.025);

  cout << "XXXXXXX fpy.getSignalYield() = " << fpy.getSignalYield() << endl;
  cout << "XXXXXXX fpy.getSignalYield() = " << fpy2.getSignalYield() << endl;

}

// ----------------------------------------------------------------------
void fitProblems8() {

  vector<string> years;
  years.push_back("2011s01");
  years.push_back("2012s01");
  years.push_back("2016BFs01");
  years.push_back("2016GHs01");
  fitPsYield *a(0);
  for (unsigned int i = 0; i < years.size(); ++i) {
    TFile *f = TFile::Open(Form("s01/plotResults.%s.root", years[i].c_str()));
    f->cd("bupsikMcComb");
    for (int ichan = 0; ichan < 2; ++ichan) {
      string hname(Form("hNorm_bdt%s_bupsikMcComb_chan%d", years[i].c_str(), ichan));
      a = new fitPsYield(hname, 0);
      a->fitBu2JpsiKp(5, years[i]);
    }
  }
}


// ----------------------------------------------------------------------
void fitProblems7() {

  vector<string> years;
  years.push_back("2011s01");
  years.push_back("2012s01");
  years.push_back("2016BFs01");
  years.push_back("2016GHs01");
  fitPsYield *a(0);
  for (unsigned int i = 0; i < years.size(); ++i) {
    TFile *f = TFile::Open(Form("s01/plotResults.%s.root", years[i].c_str()));
    f->cd("bspsiphiMcComb");
    for (int ichan = 0; ichan < 2; ++ichan) {
      string hname(Form("hNorm_bdt%s_bspsiphiMcComb_chan%d", years[i].c_str(), ichan));
      a = new fitPsYield(hname, 0);
      a->fitBs2JpsiPhi(5, years[i]);
    }
  }
}

// ----------------------------------------------------------------------
void fitProblems6(int bdtCut = 20) {

  vector<string> years;
  years.push_back("2011");
  years.push_back("2012");
  years.push_back("2016BF");
  years.push_back("2016GH");
  fitPsYield *a(0);
  for (unsigned int i = 0; i < years.size(); ++i) {
    TFile *f = TFile::Open(Form("%s/plotResults.%s.root", years[i].c_str(), years[i].c_str()));
    f->cd("bspsiphiData");
    for (int ichan = 0; ichan < 2; ++ichan) {
      string hname(Form("hNorm_bdt_%d_%s_bspsiphiData_chan%d", bdtCut, years[i].c_str(), ichan));
      a = new fitPsYield(hname, 0);
      a->fitBs2JpsiPhi(5, years[i]);
    }
  }
}


// ----------------------------------------------------------------------
void fitProblems5(int bdtCut = 20) {

  vector<string> years;
  years.push_back("2011");
  years.push_back("2012");
  years.push_back("2016BF");
  years.push_back("2016GH");
  fitPsYield *a(0);
  for (unsigned int i = 0; i < years.size(); ++i) {
    TFile *f = TFile::Open(Form("%s/plotResults.%s.root", years[i].c_str(), years[i].c_str()));
    f->cd("bupsikData");
    for (int ichan = 0; ichan < 2; ++ichan) {
      string hname(Form("hNorm_bdt_%d_%s_bupsikData_chan%d", bdtCut, years[i].c_str(), ichan));
      a = new fitPsYield(hname, 0);
      a->fitBu2JpsiKp(5, years[i]);
    }
  }
}


// ----------------------------------------------------------------------
void fitProblems4() {
  TFile *f = TFile::Open("/shome/ursl/bmm4/CMSSW_9_2_7/src/Bmm/RootAnalysis/macros/results/plotReducedOverlays.2016BF.root");

  fitPsYield a(0, 1);

  string hname("ad0bdt_bupsikData_docatrkMassPresel");
  TH1D *hc     = (TH1D*)f->Get(hname.c_str());
  if (hc) {
    a.fit0_Bu2JpsiKp(hc, -5., "fam");
    cout << "bg     yield: " << a.getBackgroundYield() << " +/- " << a.getBackgroundError() << endl;
    cout << "signal yield: " << a.getSignalYield() << " +/- " << a.getSignalError() << endl;
    cout << "signal RMS: " << a.getSignalRMS() << " signal mass = " << a.getSignalPeak() << endl;
  } else {
    cout << hname << " not found" << endl;
  }
}


// ----------------------------------------------------------------------
void fitProblems3() {
  TFile *f = TFile::Open("/shome/ursl/bmm4/CMSSW_9_2_7/src/Bmm/RootAnalysis/macros/results/plotReducedOverlays.2012.root");

  fitPsYield a(0, 1);

  string hname("ad0bdt_bupsikData_fls3dMassPresel");
  TH1D *hc     = (TH1D*)f->Get(hname.c_str());
  if (hc) {
    a.fit0_Bu2JpsiKp(hc, -5., "fam");
    cout << "bg     yield: " << a.getBackgroundYield() << " +/- " << a.getBackgroundError() << endl;
    cout << "signal yield: " << a.getSignalYield() << " +/- " << a.getSignalError() << endl;
    cout << "signal RMS: " << a.getSignalRMS() << " signal mass = " << a.getSignalPeak() << endl;
  } else {
    cout << hname << " not found" << endl;
  }
}

// ----------------------------------------------------------------------
void fitProblems2() {
  TFile *f = TFile::Open("/shome/ursl/bmm4/CMSSW_9_2_7/src/Bmm/RootAnalysis/macrosresults/plotReducedOverlays.2012.root");

  fitPsYield a(0, 1);

  string hname("ad1cnc_bdpsikstarData_fls3dMassAo");
  TH1D *hc     = (TH1D*)f->Get(hname.c_str());
  if (hc) {
    a.fit0_Bd2JpsiKstar(hc, -5., "fam");
  } else {
    cout << hname << " not found" << endl;
  }
}


// ----------------------------------------------------------------------
void fitProblems1() {
  TFile *f = TFile::Open("/shome/ursl/bupsikData_badFits_1.root");

  fitPsYield a(0, 1);

  TDirectory *pDir = gDirectory;
  TIter next(pDir->GetListOfKeys());
  TKey *key(0);
  string hname("");
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
    hname = key->GetName();
    TH1D *h1     = (TH1D*)pDir->Get(hname.c_str());
    TH1D *hc = (TH1D*)h1->Clone(Form("hc_%s", h1->GetName()));
    hc->Reset();
    for (int i = 0; i <= h1->GetNbinsX(); ++i) {
      hc->SetBinContent(i, h1->GetBinContent(i));
    }
    cout << "found hname ->" << hname << "<- and h2 = " << hc << endl;
    a.fit0_Bu2JpsiKp(hc, -5., "fam");
  }
}
