#include "util.hh"

// ----------------------------------------------------------------------
void scanBDT(string what = "CSBF", string fname = "results/scanBDT-2011.root") {

  if (what == "all") {
    scanBDT("CSBF", fname);
    c0->Modified();  c0->Update();
    scanBDT("BSMMBFS", fname);
    c0->Modified();  c0->Update();
    scanBDT("SSB", fname);
    c0->Modified();  c0->Update();
    scanBDT("SOB", fname);
    c0->Modified();  c0->Update();
    return;
  }

  c0->Clear();

  string bname("hBdt_bsmmMcComb");
  if (string::npos != what.find("CSBF")) {
    bname = "hBdt_bspsiphiMcComb";
  }

  TFile *f = TFile::Open(fname.c_str());
  string hname(Form("bdtScan_%s", what.c_str()));
  string pdfname(fname);
  replaceAll(pdfname, "results/", "");
  replaceAll(pdfname, "scanBDT-", "");
  replaceAll(pdfname, ".root", "");

  TH1D *h0(0), *h1(0);
  TH1D *b0(0), *b1(0);
  string s("");

  TH1D *hbdt = (TH1D*)f->Get("bdtCuts");
  double bdtCut0 = 100.*hbdt->GetBinContent(1);
  double bdtCut1 = 100.*hbdt->GetBinContent(2);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  s = hname + Form("_chan0");
  h0 = (TH1D*)f->Get(s.c_str());
  if (!h0) {
    cout << "Did not find histogram ->" << s << "<-" << endl;
    return;
  }
  double bdt0y   = h0->GetBinContent(h0->FindBin(bdtCut0));
  double bdt0x(0.);
  for (int ic = 0; ic < h0->GetNbinsX(); ++ic) {
    if (h0->GetBinContent(ic) > 0.) bdt0x = h0->GetBinCenter(ic);
  }
  bdt0x /= 100.;

  h0->SetMinimum(0.0);
  if (what == "CSBF") {
    h0->SetMaximum(5.e-5);
  } else if (what == "BSMMBFS") {
    h0->SetMaximum(7.e-9);
  } else {
    h0->SetMaximum(1.4*h0->GetMaximum());
  }
  h0->SetLineColor(kBlue);
  h0->Draw("hist");
  TMarker *pm = new TMarker(bdtCut0, bdt0y, 20);
  pm->SetMarkerColor(kBlue);
  pm->SetMarkerSize(2.);
  pm->Draw();
  h0->GetListOfFunctions()->Add(pm);

  tl->SetTextSize(0.05);
  tl->SetTextAngle(90.);  tl->SetTextColor(kBlue); tl->DrawLatex(0.93, 0.55, "chan 0"); tl->SetTextAngle(0.);

  s = hname + Form("_chan1");
  h1 = (TH1D*)f->Get(s.c_str());
  double bdt1y   = h1->GetBinContent(h1->FindBin(bdtCut1));
  h1->SetLineColor(kRed);
  h1->Draw("histsame");

  tl->SetTextAngle(90.); tl->SetTextColor(kRed); tl->DrawLatex(0.93, 0.75, "chan 1");  tl->SetTextAngle(0.);
  pm = new TMarker(bdtCut1, bdt1y, 21);
  pm->SetMarkerColor(kRed);
  pm->SetMarkerSize(2.);
  pm->Draw();
  h1->GetListOfFunctions()->Add(pm);

  tl->SetTextSize(0.05);
  tl->SetTextColor(kBlack); tl->DrawLatex(0.4, 0.92, Form("%s/%s", what.c_str(), pdfname.c_str()));


  s = bname + Form("_chan0");
  b0 = (TH1D*)f->Get(s.c_str());
  b0->Scale(1./b0->GetSumOfWeights());
  b0->SetMinimum(0.0);
  b0->SetLineColor(kBlue);
  b0->SetMaximum(1.4*b0->GetMaximum());

  s = bname + Form("_chan1");
  b1 = (TH1D*)f->Get(s.c_str());
  b1->Scale(1./b1->GetSumOfWeights());
  b1->SetMinimum(0.0);
  b1->SetMaximum(1.4*b1->GetMaximum());
  b1->SetLineColor(kRed);

  TPad *insetPad = new TPad("myname", "mytitle", 0.18, 0.75, 0.48, 0.88);
  insetPad->Draw();
  insetPad->cd();
  b0->Draw("hist");
  b1->Draw("samehist");

  pm = new TMarker(bdtCut0/100, b0->GetBinContent(b0->FindBin(bdtCut0/100)), 20);
  pm->SetMarkerColor(kBlue);
  pm->SetMarkerSize(1.);
  pm->Draw();
  b0->GetListOfFunctions()->Add(pm);

  pm = new TMarker(bdtCut1/100, b1->GetBinContent(b1->FindBin(bdtCut1/100)), 21);
  pm->SetMarkerColor(kRed);
  pm->SetMarkerSize(1.);
  pm->Draw();
  b1->GetListOfFunctions()->Add(pm);


  pl->DrawLine(0., b0->GetMaximum(), 0., 0.);
  cout << "bdt0x = " << bdt0x << endl;
  pl->DrawLine(bdt0x, b0->GetMaximum(), bdt0x, 0.);

  c0->SaveAs(Form("qnd/%s-%s.pdf", what.c_str(), pdfname.c_str()));

}



// ----------------------------------------------------------------------
void scanChan(string what = "CSBF", int mode = 0, int chan = 0) {

  c0->Clear();

  vector<string> inputFiles;
  vector<Color_t> colors;

  if (0 == mode) {
    inputFiles.push_back("results/scanBDT-2011.root");   colors.push_back(kRed);
    inputFiles.push_back("results/scanBDT-2012.root");   colors.push_back(kBlack);
    inputFiles.push_back("results/scanBDT-2016BF.root"); colors.push_back(kGreen+2);
    inputFiles.push_back("results/scanBDT-2016GH.root"); colors.push_back(kBlue);
  }

  string bname("hBdt_bsmmMcComb");
  if (string::npos != what.find("CSBF")) {
    bname = "hBdt_bspsiphiMcComb";
  }

  TFile *f(0);

  for (unsigned int ifile = 0; ifile < inputFiles.size(); ++ifile) {
    f = TFile::Open(inputFiles[ifile].c_str());
    string hname(Form("bdtScan_%s", what.c_str()));
    string pdfname(inputFiles[ifile].c_str());
    replaceAll(pdfname, "results/", "");
    replaceAll(pdfname, "scanBDT-", "");
    replaceAll(pdfname, ".root", "");

    TH1D *h0(0), *h1(0);
    TH1D *b0(0), *b1(0);
    string s("");

    TH1D *hbdt = (TH1D*)f->Get("bdtCuts");
    double bdtCut0 = 100.*hbdt->GetBinContent(1);
    double bdtCut1 = 100.*hbdt->GetBinContent(2);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    s = hname + Form("_chan%d", chan);
    h0 = (TH1D*)f->Get(s.c_str());
    if (!h0) {
      cout << "Did not find histogram ->" << s << "<-" << endl;
      return;
    }
    double bdt0y   = h0->GetBinContent(h0->FindBin(bdtCut0));
    double bdt0x(0.);
    for (int ic = 0; ic < h0->GetNbinsX(); ++ic) {
      if (h0->GetBinContent(ic) > 0.) bdt0x = h0->GetBinCenter(ic);
    }
    bdt0x /= 100.;

    h0->SetMinimum(0.0);
    if (what == "CSBF") {
      h0->SetMaximum(5.e-5);
    } else if (what == "BSMMBFS") {
      h0->SetMaximum(7.e-9);
    } else if (what == "SSB") {
      h0->SetMaximum(4.);
    } else {
      h0->SetMaximum(1.4*h0->GetMaximum());
    }
    h0->SetLineColor(colors[ifile]);
    h0->SetMarkerColor(colors[ifile]);
    h0->SetMarkerStyle(24);
    h0->SetMarkerSize(0.8);
    setTitles(h0, "100 #times BDT >", what.c_str(), 0.05, 1.2, 1.6);
    if (0 == ifile) {
      h0->Draw("p");
    } else {
      h0->Draw("psame");
    }
    TMarker *pm = new TMarker(bdtCut0, bdt0y, 20);
    pm->SetMarkerColor(colors[ifile]);
    pm->SetMarkerSize(2.);
    pm->Draw();
    h0->GetListOfFunctions()->Add(pm);

    tl->SetTextSize(0.03);
    tl->SetTextColor(colors[ifile]); tl->DrawLatex(0.22, 0.40 - ifile*0.035, inputFiles[ifile].c_str());
  }

  c0->SaveAs(Form("qnd/%s-%d-chan%d.pdf", what.c_str(), mode, chan));

}
