// ----------------------------------------------------------------------
double lhcbFunc(double *x, double *par) {
  // par[0] -> A
  // par[1] -> p1
  // par[2] -> p2
  // par[3] -> meanPt
  return (par[0]*(par[1] + par[2]*(x[0] - par[3])));

}

void fsfu(double xmin = 0., double xmax = 25.) {
  bool internal(false);
  double alpha(0.3);
  gPad->SetRightMargin(0.05);
  TH2D *frame = new TH2D("frame", "", 10, 0., 50., 10, 0., 0.30);
  frame->GetXaxis()->SetLimits(xmin, xmax);
  frame->GetYaxis()->SetLimits(0.15, 0.35);
  setTitles(frame, "#it{p}_{T} [GeV]", "#it{f_{#kern[-0.1]{s}} / f_{#kern[-0.1]{u}}}", 0.05, 1.1, 1.6);
  gStyle->SetOptStat(0);
  frame->Draw();

  TF1 *fLHCbOld = new TF1("fLHCbOld", lhcbFunc, 0., 30., 4);
  fLHCbOld->SetLineWidth(3);
  fLHCbOld->SetLineColorAlpha(kBlue+1, alpha);
  fLHCbOld->SetParameters(1., 0.256, -2.03e-3, 10.4);
  fLHCbOld->Draw("same");

  TF1 *fLHCbOldHi = new TF1("fLHCbOldHi", lhcbFunc, 0., 30., 4);
  fLHCbOldHi->SetLineColorAlpha(kBlue+1, alpha);
  fLHCbOldHi->SetFillColorAlpha(kBlue+1, alpha);

  fLHCbOldHi->SetParameters(1., 0.256+0.02, -2.03e-3, 10.4);
  fLHCbOldHi->SetFillStyle(3365);
  fLHCbOldHi->Draw("same");

  TF1 *fLHCbOldLo = new TF1("fLHCbOldLo", lhcbFunc, 0., 30., 4);
  fLHCbOldLo->SetLineColorAlpha(kBlue+1, alpha); fLHCbOldLo->SetFillColor(10);
  fLHCbOldLo->SetParameters(1., 0.256-0.02, -2.03e-3, 10.4);
  fLHCbOldLo->SetFillStyle(1000);
  fLHCbOldLo->Draw("same");

  TF1 *fLHCbNew = new TF1("fLHCbNew", lhcbFunc, 0., 30., 4);
  fLHCbNew->SetLineWidth(3);
  fLHCbNew->SetLineColorAlpha(kRed+1, alpha);
  fLHCbNew->SetParameters(2., 0.119, -0.91e-3, 10.1);
  fLHCbNew->Draw("samee");

  TF1 *fLHCbNewHi = new TF1("fLHCbNewHi", lhcbFunc, 0., 30., 4);
  fLHCbNewHi->SetParameters(2.+0.086, 0.119, -0.91e-3, 10.1);
  fLHCbNewHi->SetLineColorAlpha(kRed+1, alpha);  fLHCbNewHi->SetFillColorAlpha(kRed+1, alpha);
  fLHCbNewHi->SetFillStyle(3354);
  fLHCbNewHi->Draw("samee");
  TF1 *fLHCbNewLo = new TF1("fLHCbNewLo", lhcbFunc, 0., 30., 4);
  fLHCbNewLo->SetParameters(2.-0.086, 0.119, -0.91e-3, 10.1);
  fLHCbNewLo->SetLineColorAlpha(kRed+1, alpha);  fLHCbNewLo->SetFillColor(10);
  fLHCbNewLo->SetFillStyle(1000);
  fLHCbNewLo->Draw("samee");

  // -- PDG
  double mx[] = {18.5};
  double my[] = {0.252};
  double mxe[] = {0.0};
  double mye[] = {0.012};
  TGraphErrors *tge = new TGraphErrors(1, mx, my, mxe, mye);
  tge->SetMarkerStyle(20);
  tge->SetMarkerColor(kGray+3);
  tge->SetLineColor(kGray+3);
  tge->SetLineWidth(3);
  tge->SetMarkerSize(2);
  tge->Draw("psame");

  tl->SetTextSize(0.025);
  tl->SetTextColor(kBlack);
  tl->SetNDC(kFALSE);
  tl->DrawLatex(19.2, 0.250, "0.252");

  // -- ATLAS
  double ax[] = {18.};
  double ay[] = {0.240};
  double axe[] = {0.0};
  double aye[] = {0.020};
  TGraphErrors *tae = new TGraphErrors(1, ax, ay, axe, aye);
  tae->SetMarkerStyle(21);
  tae->SetMarkerColor(kGray+2);
  tae->SetLineColor(kGray+2);
  tae->SetLineWidth(3);
  tae->SetMarkerSize(2);
  tae->Draw("psame");
  tl->SetTextColor(kBlack);
  tl->DrawLatex(15.5, 0.238, "0.240");

  // -- CMS
  double cx[] = {19.};
  double cy[] = {0.231};
  double cxe[] = {0.0};
  double cye[] = {0.019};
  TGraphErrors *tce = new TGraphErrors(1, cx, cy, cxe, cye);
  tce->SetMarkerStyle(20);
  tce->SetMarkerColor(kGreen-2);
  tce->SetLineColor(kGreen-2);
  tce->SetLineWidth(3);
  tce->SetMarkerSize(2);
  if (internal) {
    tce->Draw("psame");
    tl->SetTextColor(kBlack);
    tl->DrawLatex(19.7, 0.229, "0.231");
  }

  // -- ARC
  double c1x[] = {21.5};
  double c1y[] = {0.252};
  double c1xe[] = {0.0};
  double c1ye[] = {0.019};
  TGraphErrors *tc1e = new TGraphErrors(1, c1x, c1y, c1xe, c1ye);
  tc1e->SetMarkerStyle(20);
  tc1e->SetMarkerColor(kCyan-2);
  tc1e->SetLineColor(kCyan-2);
  tc1e->SetLineWidth(3);
  tc1e->SetMarkerSize(2);
  if (internal) {
    tc1e->Draw("psame");
    tc1e->Draw("e1same");
  }


  // -- LHCb 13 TeV rescaled 1.068
  double lhcbx[] = {10.4};
  double lhcby[] = {0.259/1.068}; // NB: here it should be 0.259 which is the value for fs/fu in LHCb-CONF-2013-011!!!
  double lhcbxe[] = {0.0};
  double lhcbye[] = {0.019/1.068};
  TGraphErrors *tlhcb = new TGraphErrors(1, lhcbx, lhcby, lhcbxe, lhcbye);
  tlhcb->SetMarkerStyle(26);
  tlhcb->SetMarkerColor(kRed+1);
  tlhcb->SetLineColor(kRed+1);
  tlhcb->SetLineWidth(3);
  tlhcb->SetMarkerSize(2);
  tlhcb->Draw("psame");
  tl->SetTextColor(kRed+1);
  tl->DrawLatex(8.0, 0.233, "0.243");


  // -- LHCb 13 TeV as in paper
  double lhcb1x[] = {10.6};
  double lhcb1y[] = {0.244};
  double lhcb1xe[] = {0.0};
  double lhcb1ye[] = {0.012};
  TGraphErrors *tlhcb1 = new TGraphErrors(1, lhcb1x, lhcb1y, lhcb1xe, lhcb1ye);
  tlhcb1->SetMarkerStyle(20);
  tlhcb1->SetMarkerColor(kRed+1);
  tlhcb1->SetLineColor(kRed+1);
  tlhcb1->SetLineWidth(3);
  tlhcb1->SetMarkerSize(2);
  tlhcb1->Draw("psame");
  tl->SetTextColor(kRed+1);
  tl->DrawLatex(11.0, 0.244, "0.244");


  // -- LHCb 7 TeV
  double lhcb0x[] = {10.2};
  double lhcb0y[] = {0.259};
  double lhcb0xe[] = {0.0};
  double lhcb0ye[] = {0.015};
  TGraphErrors *tlhcb0 = new TGraphErrors(1, lhcb0x, lhcb0y, lhcb0xe, lhcb0ye);
  tlhcb0->SetMarkerStyle(22);
  tlhcb0->SetMarkerColor(kBlue+1);
  tlhcb0->SetLineColor(kBlue+1);
  tlhcb0->SetLineWidth(3);
  tlhcb0->SetMarkerSize(2);
  tlhcb0->Draw("psame");
  tl->SetTextColor(kBlue+1);
  tl->DrawLatex(8.0, 0.250, "0.259");


  double x = 15.;
  cout << "fs/(fu)[" << x << "] = " << fLHCbNew->Eval(x) << "(new) and " << fLHCbOld->Eval(x) << " (new)" << endl;
  x = 20.;
  cout << "fs/(fu)[" << x << "] = " << fLHCbNew->Eval(x) << "(new) and " << fLHCbOld->Eval(x) << " (new)" << endl;

  frame->Draw("axissame");

  TLegend *tle = new TLegend(0.17, 0.14, 0.45, 0.37);
  tle->SetFillStyle(0);
  tle->SetBorderSize(0);
  tle->SetTextSize(0.03);
  tle->SetTextFont(42);

  tle->AddEntry(tge, "PDG #sqrt{#it{s}} = 7TeV (PR,D98,030001; p728)", "p");
  tle->AddEntry(tae, "ATLAS #sqrt{#it{s}} = 7TeV (PRL,115,262001)", "p");
  if (internal) tle->AddEntry(tce, "CMS #sqrt{#it{s}} = 13TeV (BPH-16-004, AN-16-178), no p_{T} dep.", "p");
  tle->AddEntry(fLHCbOldHi, "LHCb #sqrt{#it{s}} = 7TeV (JHEP,1304,001)", "f");
  tle->AddEntry(fLHCbNewHi, "LHCb #sqrt{#it{s}} = 13TeV (PR, D100, 031102)", "f");
  tle->AddEntry(tlhcb, "LHCb scaled to 13 TeV (LHCb-CONF-2013-011)", "p");
  if (internal) tle->AddEntry(tc1e, "ARC BPH-16-004", "p");
  tle->Draw();

  TMarker *tma = new TMarker();
  tma->SetMarkerColor(kRed+1); tma->SetMarkerSize(2); tma->SetMarkerStyle(20);
  if (internal) {
    tma->DrawMarker(1.55, 0.174);
  } else {
    tma->DrawMarker(1.55, 0.1706);
  }
  tma->SetMarkerColor(kBlue+1); tma->SetMarkerSize(2); tma->SetMarkerStyle(22);
  if (internal) {
    tma->DrawMarker(1.55, 0.1829);
  } else {
    tma->DrawMarker(1.55, 0.1835);
  }
  double fsfu10 = fLHCbNew->Eval(10.1);
  double fsfu18 = fLHCbNew->Eval(17.2);

  if (internal) {
    tl->SetTextSize(0.03);
    tl->SetTextColor(kRed+1);
    tl->DrawLatexNDC(0.17, 0.85, Form("#it{f_{#kern[-0.1]{s}}/f_{#kern[-0.1]{u}}}(17.2 GeV) = %4.3f (= PDG #minus %2.1f%%)", fsfu18, 100.*(1.-fsfu18/0.252)));
    tl->DrawLatexNDC(0.17, 0.80, Form("#it{f_{#kern[-0.1]{s}}/f_{#kern[-0.1]{u}}}(10.1 GeV) = %4.3f (LHCb difference: #minus %2.1f%%)", fsfu10, 100.*(1.-fsfu18/fsfu10)));

    pl->SetLineColor(kGreen-2);
    pl->SetLineWidth(4);
    pl->DrawLine(15., 0.30, 21., 0.30);


    pl->SetLineColor(kGreen-2);
    pl->SetLineWidth(2);
    pl->DrawLine(17.2, 0.15, 17.2, 0.16);
    pl->DrawLine(17.2, 0.21, 17.2, 0.305);
    tl->SetTextColor(kGreen-2);
    tl->DrawLatexNDC(0.63, 0.75, "<#it{p}_{T}> = 17.2 GeV");
  }
  if (internal) {
    c0->SaveAs("qnd/fsfu.pdf");
  } else {
    c0->SaveAs("qnd/fsfu-notInternal.pdf");
  }
}


// ----------------------------------------------------------------------
void ptAverage() {
  vector<double> y, w;
  y.push_back(16.5);    w.push_back(4.1);
  y.push_back(15.3);    w.push_back(2.2);
  y.push_back(17.7);    w.push_back(11.1);
  y.push_back(15.3);    w.push_back(5.4);
  y.push_back(20.1);    w.push_back(4.6);
  y.push_back(17.9);    w.push_back(9.1);
  y.push_back(20.9);    w.push_back(4.1);
  y.push_back(18.2);    w.push_back(4.4);

  double result(0.), wtot(0.);
  for (unsigned int i = 0; i < y.size(); ++i) {
    wtot += w[i];
    result += w[i]*y[i];
  }

  result /= wtot;
  cout << "average pT = " << result << endl;

}


// ----------------------------------------------------------------------
void ptAverageFromTree(string era = "2011") {
  map<string, vector<double> > all;
  map<string, vector<double> > w8;
  vector<double> cuts, w;

  cuts.clear(); cuts.push_back(0.28); cuts.push_back(1.0);
  w.clear(); w.push_back(3.6);
  all.insert(make_pair("2011/chan0", cuts));
  w8.insert(make_pair("2011/chan0", w));
  cuts.clear(); cuts.push_back(0.21); cuts.push_back(1.0);
  w.clear(); w.push_back(2.0);
  all.insert(make_pair("2011/chan1", cuts));
  w8.insert(make_pair("2011/chan1", w));

  cuts.clear(); cuts.push_back(0.27); cuts.push_back(0.35); cuts.push_back(1.0);
  w.clear(); w.push_back(3.7); w.push_back(9.3);
  all.insert(make_pair("2012/chan0", cuts));
  w8.insert(make_pair("2012/chan0", w));
  cuts.clear(); cuts.push_back(0.23); cuts.push_back(0.32);  cuts.push_back(1.0);
  w.clear(); w.push_back(1.7); w.push_back(4.7);
  all.insert(make_pair("2012/chan1", cuts));
  w8.insert(make_pair("2012/chan1", w));

  cuts.clear(); cuts.push_back(0.19); cuts.push_back(0.30);  cuts.push_back(1.0);
  w.clear(); w.push_back(2.2); w.push_back(4.0);
  all.insert(make_pair("2016BF/chan0", cuts));
  w8.insert(make_pair("2016BF/chan0", w));
  cuts.clear(); cuts.push_back(0.19); cuts.push_back(0.30);  cuts.push_back(1.0);
  w.clear(); w.push_back(3.7); w.push_back(8.1);
  all.insert(make_pair("2016BF/chan1", cuts));
  w8.insert(make_pair("2016BF/chan1", w));

  cuts.clear(); cuts.push_back(0.18); cuts.push_back(0.31);  cuts.push_back(1.0);
  w.clear(); w.push_back(4.1); w.push_back(3.6);
  all.insert(make_pair("2016GH/chan0", cuts));
  w8.insert(make_pair("2016GH/chan0", w));
  cuts.clear(); cuts.push_back(0.23); cuts.push_back(0.38);  cuts.push_back(1.0);
  w.clear(); w.push_back(6.1); w.push_back(3.9);
  all.insert(make_pair("2016GH/chan1", cuts));
  w8.insert(make_pair("2016GH/chan1", w));

  TFile *f(0);
  TTree *t(0);
  double sum(0.), w8pt(0.);
  for (map<string, vector<double> >::iterator it = all.begin(); it != all.end(); ++it) {
    int chan = 0;
    if (string::npos != it->first.find("chan1")) chan = 1;
    if (string::npos != it->first.find("2011")) {
      f = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Summer17_private-BsToMuMu-s01.root");
    } else if (string::npos != it->first.find("2012")) {
      f = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-Winter17_private-BsToMuMu-s01.root");
    } else if (string::npos != it->first.find("2016BF")) {
      f = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-combined-BsToMuMu-2016BF-s01.root");
    } else if (string::npos != it->first.find("2016GH")) {
      f = TFile::Open("/scratch/ursl/bmm4/s01/bmm-mc-combined-BsToMuMu-2016GH-s01.root");
    }
    t = (TTree*)f->Get("candAnaMuMu/events");
    TH1D *h1 = new TH1D("h1", "", 200, 0., 50.);
    for (unsigned int ib = 0; ib  < it->second.size()-1; ++ib) {
      double bdtcut0 = it->second[ib];
      double bdtcut1 = it->second[ib+1];
      string cuts = Form("%d==chan && hlt1 && tos && l1t && gmuid && bdt > %f && bdt < %f", chan, bdtcut0, bdtcut1);
      t->Draw("pt>>h1", cuts.c_str());
      double weight = w8[it->first][ib];
      cout << it->first << ": <pT> = " << h1->GetMean() << " weight: " << weight << endl;
      w8pt += h1->GetMean()*weight;
      sum += weight;
    }
    delete h1;
  }
  cout << "average pt = " << w8pt/sum << endl;
}
