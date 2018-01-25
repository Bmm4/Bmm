// ----------------------------------------------------------------------
void fitProblems(int bdtCut = 20) {

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
