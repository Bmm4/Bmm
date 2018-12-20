void checkData(string era = "2011") {

  if ("all" == era) {
    checkData("2011");
    checkData("2012");
    checkData("2016BF");
    checkData("2016GH");
    return;
  }

  map<string, string> uFiles;
  uFiles.insert(make_pair("2011",   "/scratch/ursl/bmm4/s01/unblinded/bmm-data-bmmMuOnia2011-s01.root"));
  uFiles.insert(make_pair("2012",   "/scratch/ursl/bmm4/s01/unblinded/bmm-data-bmmMuOnia2012-s01.root"));
  uFiles.insert(make_pair("2016BF", "/scratch/ursl/bmm4/s01/unblinded/bmm-data-bmmCharmonium2016BF-s01.root"));
  uFiles.insert(make_pair("2016GH", "/scratch/ursl/bmm4/s01/unblinded/bmm-data-bmmCharmonium2016GH-s01.root"));

  map<string, string> bFiles;
  bFiles.insert(make_pair("2011", "/scratch/ursl/bmm4/s01/blinded/bmm-data-bmmMuOnia2011-s01.root"));
  bFiles.insert(make_pair("2012", "/scratch/ursl/bmm4/s01/blinded/bmm-data-bmmMuOnia2012-s01.root"));
  bFiles.insert(make_pair("2016BF", "/scratch/ursl/bmm4/s01/blinded/bmm-data-bmmCharmonium2016BF-s01.root"));
  bFiles.insert(make_pair("2016GH", "/scratch/ursl/bmm4/s01/blinded/bmm-data-bmmCharmonium2016GH-s01.root"));

  map<string, float> mmva;
  mmva.insert(make_pair("2011", 0.55));
  mmva.insert(make_pair("2012", 0.55));
  mmva.insert(make_pair("2016BF", 0.58));
  mmva.insert(make_pair("2016GH", 0.58));

  map<string, pair<float, float> > bdtCuts;
  bdtCuts.insert(make_pair("2011",   make_pair(0.28, 0.21)));
  bdtCuts.insert(make_pair("2012",   make_pair(0.34, 0.32)));
  bdtCuts.insert(make_pair("2016BF", make_pair(0.30, 0.30)));
  bdtCuts.insert(make_pair("2016GH", make_pair(0.31, 0.38)));

  map<string, pair<float, float> > expB0Window;
  map<string, pair<float, float> > expBsWindow;
  expB0Window.insert(make_pair("2011",   make_pair(1.24, 1.13)));
  expBsWindow.insert(make_pair("2011",   make_pair(4.47, 2.79)));

  expB0Window.insert(make_pair("2012",   make_pair(5.33, 3.89)));
  expBsWindow.insert(make_pair("2012",   make_pair(15.28, 7.82)));

  expB0Window.insert(make_pair("2016BF", make_pair(2.11, 8.37)));
  expBsWindow.insert(make_pair("2016BF", make_pair(6.38, 17.64)));

  expB0Window.insert(make_pair("2016GH", make_pair(1.43, 1.68)));
  expBsWindow.insert(make_pair("2016GH", make_pair(4.86, 4.87)));

  TFile *uf = TFile::Open(uFiles[era].c_str());
  TFile *bf = TFile::Open(bFiles[era].c_str());

  TTree *bt = (TTree*)bf->Get("candAnaMuMu/events");
  TTree *ut = (TTree*)uf->Get("candAnaMuMu/events");

  string masscut = "((m > 4.9 && m < 5.2) || (m > 5.45 && m < 5.9))";

  c0->Clear();
  c0->Divide(2, 2);

  TH1D *hb0 = new TH1D("hb0", "", 40, 4.9, 5.9);
  TH1D *hu0 = new TH1D("hu0", "", 40, 4.9, 5.9);

  TH1D *hb1 = new TH1D("hb1", "", 40, 4.9, 5.9);
  TH1D *hu1 = new TH1D("hu1", "", 40, 4.9, 5.9);
  int ipad(1);
  for (int ichan = 0; ichan < 2; ++ichan) {

    string cuts = Form("hlt1 && tos && l1t && (m1rmvabdt>%f && m2rmvabdt > %f) && %s && (%d == chan) && (bdt > %f)",
		       mmva[era], mmva[era], masscut.c_str(), ichan, (ichan==0?bdtCuts[era].first:bdtCuts[era].second));

    cout << "cuts: " << cuts << endl;


    c0->cd(ipad++);
    bt->Draw(Form("m>>hb%d", ichan), cuts.c_str());
    ut->Draw(Form("m>>hu%d", ichan), cuts.c_str(), "esame");
    //    tl->DrawLatexNDC(0.4, 0.2, "blind");

    if (0 == ichan) {
      tl->DrawLatexNDC(0.2, 0.92, Form("%.0f", hu0->Integral(1, 12)));
      tl->DrawLatexNDC(0.6, 0.92, Form("%.0f", hu0->Integral(23, 40)));
    } else {
      tl->DrawLatexNDC(0.2, 0.92, Form("%.0f", hu1->Integral(1, 12)));
      tl->DrawLatexNDC(0.6, 0.92, Form("%.0f", hu1->Integral(23, 40)));
    }

    c0->cd(ipad++);
    if (0 == ichan) {
      TH1D *hr0 = (TH1D*)hu0->Clone("hr0");
      hr0->Reset();
      hr0->Add(hb0, hu0, 1., -1.);
      hr0->Draw();

      tl->DrawLatexNDC(0.66, 0.36, "obs./exp.");
      tl->DrawLatexNDC(0.62, 0.29, Form("B0: %2.0f/%.2f", hu0->Integral(13, 16), expB0Window[era].first));
      tl->DrawLatexNDC(0.62, 0.22, Form("Bs: %2.0f/%.2f", hu0->Integral(17, 22), expBsWindow[era].first));
    } else {
      TH1D *hr1 = (TH1D*)hu1->Clone("hr1");
      hr1->Reset();
      hr1->Add(hb1, hu1, 1., -1.);
      hr1->Draw();
      tl->DrawLatexNDC(0.66, 0.36, "obs./exp.");
      tl->DrawLatexNDC(0.62, 0.29, Form("B0: %2.0f/%.2f", hu1->Integral(13, 16), expB0Window[era].second));
      tl->DrawLatexNDC(0.62, 0.22, Form("Bs: %2.0f/%.2f", hu1->Integral(17, 22), expBsWindow[era].second));
    }
    tl->DrawLatexNDC(0.4, 0.92, Form("%s/chan%d", era.c_str(), ichan));

  }
  c0->SaveAs(Form("s01/checkData-%s.pdf", era.c_str()));
}
