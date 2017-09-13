// ----------------------------------------------------------------------
void fitProblems() {
  TFile *f = TFile::Open("results/plotReducedOverlays.2012.root");

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
