#include "/Users/ursl/macros/ana/Bmm/RootAnalysis/common/util.hh"

// ----------------------------------------------------------------------
TLegend* newLegend(string title, double x1, double y1, double x2, double y2,
		   vector<TH1D*> hists, vector<string> names, vector<string> options) {
  if (hists.size() != names.size()) {
    cout << "hists and names vectors do not match size" << endl;
    return 0;
  }
  if (hists.size() != options.size()) {
    cout << "hists and options vectors do not match size" << endl;
    return 0;
  }
  TLegend *legg = new TLegend(x1, y1, x2, y2, title.c_str());
  legg->SetFillStyle(0);
  legg->SetBorderSize(0);
  legg->SetTextSize(0.04);
  legg->SetFillColor(0);
  legg->SetTextFont(52);
  for (unsigned int i = 0; i < hists.size(); ++i) {
    legg->AddEntry(hists[i], names[i].c_str(), options[i].c_str());
  }
  return legg;
}


// ----------------------------------------------------------------------
void setFilledHist(TH1 *h, Int_t color, Int_t fillcolor, Int_t fillstyle, double width) {
  // Note: 3004, 3005 are crosshatches
  // ----- 1000       is solid
  //       kYellow    comes out gray on bw printers
  h->SetLineColor(color);     h->SetLineWidth(width);
  h->SetFillStyle(fillstyle); h->SetFillColor(fillcolor);
}

// ----------------------------------------------------------------------
void setHist(TH1 *h, Int_t color, Int_t symbol, Double_t size, Double_t width) {
  h->SetLineColor(color);   h->SetLineWidth(width);
  h->SetMarkerColor(color); h->SetMarkerStyle(symbol);  h->SetMarkerSize(size);
  h->SetStats(kFALSE);
  h->SetFillStyle(0); h->SetFillColor(color);
}

// ----------------------------------------------------------------------
void normalization(string fname, int chan = 0) {
  c0->Clear();
  c0->SetCanvasSize(700, 700);

  TFile *f = TFile::Open(fname.c_str());

  string suffix = "bdt2011s01";
  if (string::npos != fname.find("2012")) {
    suffix = "bdt2012s01";
  } else if (string::npos != fname.find("2016BF")) {
    suffix = "bdt2016BFs01";
  } else if (string::npos != fname.find("2016GH")) {
    suffix = "bdt2016GHs01";
  }
  string sample = "bupsikData";
  //  string  name = Form("hNorm_%s_%s_chan%d", modifier.c_str(), fSample.c_str(), chan);
  string  name = Form("hNorm_%s_%s_chan%d", suffix.c_str(), sample.c_str(), chan);
  bool ok = f->cd(sample.c_str());
  cout << "cd to " << sample << ": " << ok << ", normalization fitting: " << name << endl;

  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fitPsYield fpy(name, gDirectory, -1);
  fpy.fitBu2JpsiKp(5, "qnd/pub", 0, 5.01, 5.8, 0.025);

}

// ----------------------------------------------------------------------
void controlsample(string fname, int chan = 0) {
  c0->Clear();
  c0->SetCanvasSize(700, 700);

  TFile *f = TFile::Open(fname.c_str());

  string suffix = "bdt2011s01";
  if (string::npos != fname.find("2012")) {
    suffix = "bdt2012s01";
  } else if (string::npos != fname.find("2016BF")) {
    suffix = "bdt2016BFs01";
  } else if (string::npos != fname.find("2016GH")) {
    suffix = "bdt2016GHs01";
  }
  string sample = "bspsiphiData";
  //  string  name = Form("hNorm_%s_%s_chan%d", modifier.c_str(), fSample.c_str(), chan);
  string  name = Form("hNorm_%s_%s_chan%d", suffix.c_str(), sample.c_str(), chan);
  bool ok = f->cd(sample.c_str());
  cout << "cd to " << sample << ": " << ok << ", control sample fitting: " << name << endl;

  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  fitPsYield fpy(name, gDirectory, -1);
  fpy.fitBs2JpsiPhi(5, "qnd/pub", 1, 4.99, 5.8, 0.025);

}


// ----------------------------------------------------------------------
void overlayAndRatio(TH1D *h1, TH1D *h2, string year, int log = 0) {
  bool drawGrid(true), preliminary(false), cmsstyle(true);
  string sname = string(h1->GetName());
  bool isBmmData = (string::npos != sname.find("bmmData"));

  // -- Upper plot
  double splity(0.3);
  TCanvas *c = new TCanvas("c1", "c1", 0, 0, 700, 800);
  c->cd();
  TPad *pad1;
  if (!isBmmData) {
    pad1 = new TPad("pad1", "pad1", 0.0, splity, 1.0, 1.0);
    pad1->SetTopMargin(0.1);
    pad1->SetBottomMargin(0.);
    pad1->SetRightMargin(0.04);
    pad1->SetLeftMargin(0.24);
    //    if (drawGrid) pad1->SetGridx();
    pad1->Draw();
    pad1->cd();
  }

  if (1 == log) {
    c->SetLogy(1);
  }

  h1->SetTitle("");
  double ymax = h1->GetMaximum();
  if (h2->GetMaximum() > ymax) ymax = h2->GetMaximum();
  h1->SetStats(0);
  string hname = h1->GetName();
  if (string::npos != hname.find("Data")) {
    h1->Draw("e");
  } else {
    h1->Draw();
  }
  hname = h2->GetName();
  if (string::npos != hname.find("Data")) {
    h2->Draw("samee");
  } else {
    h2->Draw("samehist");
  }
  h1->Draw("esame");

  h1->SetMaximum(1.3*ymax);
  h1->SetMinimum(-0.01*ymax);
  h1->SetMinimum(0.001*ymax);
  if (1 == log) {
  h1->SetMaximum(10.*ymax);
  h1->SetMinimum(0.5);
  }
  // -- Lower plot
  TH1D *hr = (TH1D*)h1->Clone("hr");
  if (!isBmmData) {
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0., 1.0, splity);
    pad2->SetTopMargin(0.);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.04);
    pad2->SetLeftMargin(0.24);
    pad2->Draw();
    if (drawGrid) pad2->SetGridy();
    if (drawGrid) pad2->SetGridx();
    pad2->cd();

    hr->SetLineColor(kBlack);
    hr->SetMinimum(0.65);
    hr->SetMaximum(1.35);
    hr->SetStats(0);
    hr->Divide(h2);
    hr->SetMarkerStyle(24);

    hr->Draw("e0");
    TBox *lbox = new TBox();
    lbox->SetFillStyle(1000);
    lbox->SetFillColor(kGray);
    lbox->DrawBox(hr->GetBinLowEdge(1), 0.80, hr->GetBinLowEdge(hr->GetNbinsX()+1), 1.20);
    hr->Draw("e0same");
    hr->Draw("axissame");
    hr->Draw("axigsame");
    pl->DrawLine(hr->GetBinLowEdge(1), 1., hr->GetBinLowEdge(hr->GetNbinsX()+1), 1.0);
  } else {
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.22);

    TArrow *tar = new TArrow();
    tar->SetLineWidth(2);
    tar->SetAngle(15);
    tar->SetFillColor(kCyan+2);
    double ymin(0.57), asize(0.1);
    if (year == "2011") {
      tar->DrawArrow(0.28, 0.3*ymax, 0.28, ymin, asize, "|>");
    }
    if (year == "2012") {
      tar->DrawArrow(0.27, 0.3*ymax, 0.27, ymin, asize, "|>");
      tar->DrawArrow(0.35, 0.3*ymax, 0.35, ymin, asize, "|>");
    }
    if (year == "2016BF") {
      tar->DrawArrow(0.27, 0.3*ymax, 0.27, ymin, asize, "|>");
      tar->DrawArrow(0.35, 0.3*ymax, 0.35, ymin, asize, "|>");
    }
    if (year == "2016GH") {
      tar->DrawArrow(0.27, 0.3*ymax, 0.27, ymin, asize, "|>");
      tar->DrawArrow(0.35, 0.3*ymax, 0.35, ymin, asize, "|>");
    }

  }

  hr->GetYaxis()->SetNdivisions(504);
  h1->GetYaxis()->SetTitleOffset(1.5);
  string xtitle(h1->GetXaxis()->GetTitle()), ytitle("");
  bool doLegend(false);
  if (preliminary) doLegend = true;
  double shift(1.72);
  if (string::npos != hname.find("pt")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    xtitle = Form("#it{p}_{T} [GeV]");
    ytitle = Form("Candidates / %2.1f GeV", h1->GetBinWidth(1));
    doLegend = true;
  }  else if (string::npos != hname.find("m1iso")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    ytitle = Form("Candidates / %3.2f", h1->GetBinWidth(1));
    xtitle = Form("muon isolation");
  }  else if (string::npos != hname.find("m2iso")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    ytitle = Form("Candidates / %3.2f", h1->GetBinWidth(1));
    xtitle = Form("muon isolation");
  }  else if (string::npos != hname.find("iso")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    ytitle = Form("Candidates / %3.2f", h1->GetBinWidth(1));
    xtitle = Form("isolation");
  }  else if (string::npos != hname.find("alpha")) {
    h1->GetYaxis()->SetTitleOffset(1.85);
    ytitle = Form("Candidates / %4.3f", h1->GetBinWidth(1));
    xtitle = Form("#it{#alpha}_{3D}");
  }  else if (string::npos != hname.find("maxdoca")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    ytitle = Form("Candidates / %4.3f cm", h1->GetBinWidth(1));
  }  else if (string::npos != hname.find("chi2dof")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    ytitle = Form("Candidates / %2.1f", h1->GetBinWidth(1));
  }  else if (string::npos != hname.find("docatrk")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    ytitle = Form("Candidates / %4.3f cm", h1->GetBinWidth(1));
  }  else if (string::npos != hname.find("muhelicity")) {
    if (doLegend) h1->SetMaximum(1.2*h1->GetMaximum());
    h1->GetYaxis()->SetTitleOffset(shift);
    ytitle = Form("Candidates / %3.2f", h1->GetBinWidth(1));
    xtitle = Form("cos(#it{#theta_{#mu^{ #minus}}})");
  }  else if (string::npos != hname.find("ips")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    xtitle = Form("#delta_{3D} / #sigma(#delta_{3D})");
    ytitle = Form("Candidates/Bin");
    ytitle = Form("Candidates / %3.2f", h1->GetBinWidth(1));
  }  else if ((string::npos != hname.find("ip")) && (string::npos == hname.find("bspsiphi"))) {
    h1->GetYaxis()->SetTitleOffset(shift);
    xtitle = Form("#delta_{3D} [cm]");
    ytitle = Form("Candidates / %5.4f cm", h1->GetBinWidth(1));
    hr->GetXaxis()->SetNdivisions(504);
  }  else if (string::npos != hname.find("fls3d")) {
    doLegend = true;
    h1->GetYaxis()->SetTitleOffset(shift);
    xtitle = Form("#it{l}_{3D} / #it{#sigma}(#it{l}_{3D})");
    ytitle = Form("Candidates / Bin");
    //    ytitle = Form("Candidates / %1.0f", h1->GetBinWidth(1));
  }  else if (string::npos != hname.find("closetrk")) {
    h1->GetYaxis()->SetTitleOffset(1.85);
    ytitle = Form("Candidates / Bin");
    xtitle = Form("#it{N}_{ trk}^{ close}");
  }  else if (string::npos != hname.find("tau")) {
    h1->GetYaxis()->SetTitleOffset(shift);
    xtitle = Form("#it{t} [ps]");
    ytitle = Form("Candidates / %2.1f ps", h1->GetBinWidth(1));
  }  else if (string::npos != hname.find("bdtsel")) {
    xtitle = Form("BDT discriminator");
    ytitle = Form("Candidates / %3.2f", h1->GetBinWidth(1));
    h1->GetYaxis()->SetTitleOffset(shift);
    if ("2011" == year) {
      doLegend = true;
    }
  }
  double psize = 0.07;
  double pratio = (1-splity)/splity;
  h1->GetYaxis()->SetTitleSize(psize);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelSize(psize);
  h1->GetYaxis()->SetTitle(ytitle.c_str());
  h1->GetYaxis()->CenterTitle();
  h1->GetYaxis()->SetNoExponent();

  h1->GetXaxis()->SetLabelSize(psize);


  if (isBmmData) {
    h1->GetXaxis()->SetTitle(xtitle.c_str());
    h1->GetYaxis()->SetTitleSize(0.8*psize);
    h1->GetXaxis()->CenterTitle();

    h1->GetYaxis()->SetTitleFont(42);
    h1->GetYaxis()->SetTitleOffset(1.7);
    h1->GetYaxis()->SetLabelSize(0.65*psize);
    h1->GetXaxis()->SetLabelSize(0.65*psize);
    h1->GetYaxis()->SetTitle(ytitle.c_str());
    h1->GetYaxis()->CenterTitle();
  }
  // hr settings
  hr->SetTitle("");
  hr->GetXaxis()->SetTitle(xtitle.c_str());
  hr->GetXaxis()->SetTitleOffset(1.15);
  hr->GetXaxis()->CenterTitle();

  hr->GetYaxis()->SetTitle("Data / MC");
  hr->GetYaxis()->SetTitleOffset(0.5);
  hr->GetYaxis()->CenterTitle();

  hr->GetYaxis()->SetTitleFont(42);
  hr->GetXaxis()->SetTitleFont(42);
  hr->GetYaxis()->SetTitleSize(0.14);
  hr->GetXaxis()->SetTitleSize(pratio*psize);
  hr->GetYaxis()->SetLabelSize(pratio*psize);
  hr->GetXaxis()->SetLabelSize(pratio*psize);
  hr->GetXaxis()->SetLabelOffset(0.02);
  hr->GetXaxis()->CenterTitle();

  if (!isBmmData) pad1->cd();

  if (doLegend) {
    TLegend *legg(0);
    if (!isBmmData && (string::npos == hname.find("bdtsel"))) {
      legg = new TLegend(0.6, 0.65, 0.85, 0.85);
      legg->SetTextFont(42);
      if (string::npos != hname.find("bupsik")) {
	if (cmsstyle) {
	  legg->SetHeader("B^{+} #rightarrow J/#kern[-0.0]{#it{#psi}}K^{+}");
	} else {
	  legg->SetHeader("#it{B^{+}} #rightarrow #it{J}/#kern[-0.0]{#it{#psi}} #it{K^{+}}");
	}
      } else if (string::npos != hname.find("bspsiphi")) {
	if (cmsstyle) {
	  legg->SetHeader("B^{0}_{s} #rightarrow J/#kern[-0.0]{#it{#psi}}#it{#phi}");
	} else {
	  legg->SetHeader("#it{B^{0}_{s}} #rightarrow #it{J}/#kern[-0.0]{#it{#psi}} #it{#phi}");
	}
      }
      legg->AddEntry(h1, "Data", "pe");
      legg->AddEntry(h2, "Simulation", "f");
      legg->SetTextSize(0.06);
    }
    if (!isBmmData && (string::npos != hname.find("bdtsel"))) {
      legg = new TLegend(0.26, 0.70, 0.50, 0.87);
      legg->SetTextFont(42);
      if (string::npos != hname.find("bupsik")) {
      if (cmsstyle) {
	legg->SetHeader("B^{+} #rightarrow J/#kern[-0.0]{#it{#psi}}K^{+}");
      } else {
	legg->SetHeader("#it{B^{+}} #rightarrow #it{J}/#kern[-0.0]{#it{#psi}} #it{K^{+}}");
      }
      } else if (string::npos != hname.find("bspsiphi")) {
	if (cmsstyle) {
	  legg->SetHeader("B^{0}_{s} #rightarrow J/#kern[-0.0]{#it{#psi}}#it{#phi}");
	} else {
	  legg->SetHeader("#it{B^{0}_{s}} #rightarrow #it{J}/#kern[-0.0]{#it{#psi}}#it{#phi}");
	}
      }
      legg->AddEntry(h1, "Data", "pe");
      legg->AddEntry(h2, "Simulation", "f");
      legg->SetTextSize(0.064);
    }
    if (isBmmData && (string::npos != hname.find("bdtsel"))) {
      legg = new TLegend(0.26, 0.75, 0.60, 0.87);
      legg->SetTextFont(42);
      legg->AddEntry(h1, "Background (Data sideband)", "pe");
      if (cmsstyle) {
	legg->AddEntry(h2, "B^{0}_{s} #rightarrow #it{#mu^{+}#mu^{#minus}} (Simulation)", "f");
      } else {
	legg->AddEntry(h2, "#it{B^{0}_{s}} #rightarrow #it{#mu^{+}#mu^{#minus}} (Simulation)", "f");
      }
      legg->SetTextSize(0.05);
    }
    legg->SetFillStyle(0);
    legg->SetBorderSize(0);
    legg->SetFillColor(kWhite);
    legg->SetTextFont(42);
    legg->Draw();
  }

  TLatex tl;
  tl.SetTextAlign(11);
  if (isBmmData) {
    tl.SetTextSize(0.056);
  } else {
    tl.SetTextSize(0.07);
  }
  tl.SetTextFont(62);
  // if (string::npos != hname.find("chi2dof")) {
  //   tl.DrawLatexNDC(0.35, 0.91, "CMS");
  // } else {
  //   tl.DrawLatexNDC(0.20, 0.91, "CMS");
  // }
  tl.DrawLatexNDC(0.22, 0.92, "CMS");
  if (preliminary) {
    tl.SetTextFont(42);
    tl.DrawLatexNDC(0.33, 0.91, "#it{Preliminary}");
    tl.SetTextFont(62);
  }

  tl.SetTextFont(42);

  tl.SetTextAlign(31);
  if ("2011" == year) {
    tl.DrawLatexNDC(0.92, 0.91, "5 fb^{-1} (7 TeV)");
  } else if ("2012" == year) {
    tl.DrawLatexNDC(0.92, 0.91, "20 fb^{-1} (8 TeV)");
  } else if ("2016BF" == year) {
    tl.DrawLatexNDC(0.92, 0.91, "20 fb^{-1} (13 TeV)");
  } else if ("2016GH" == year) {
    tl.DrawLatexNDC(0.92, 0.91, "16 fb^{-1} (13 TeV)");
  }

  string slog("");
  if (1 == log) slog = "_log";
  c->SaveAs(Form("qnd/pub_overlay_%s_%s%s.pdf", year.c_str(), h1->GetName(), slog.c_str()));
  c->SaveAs(Form("qnd/pub_overlay_%s_%s%s.C", year.c_str(), h1->GetName(), slog.c_str()));

}


// ----------------------------------------------------------------------
void overlay(string var, string sample, string fname = "s01bis/plotSbsHistograms-2016BFs01.root", string year = "2016BF", int chan = 0, int log = 0) {
  TFile *f = TFile::Open(fname.c_str());
  vector<string> variables;
  variables.push_back(var);

  TH1D *hd(0), *hm(0);
  string sel("Presel");
  string mcsample = sample;
  if (sample == "bmm") mcsample = "bsmm";
  for (unsigned int iv = 0; iv < variables.size(); ++iv) {
    hm = (TH1D*)f->Get(Form("sbs_ad%dbdt_%sMcComb_%s%s", chan, mcsample.c_str(), variables[iv].c_str(), sel.c_str()));
    hd = (TH1D*)f->Get(Form("sbs_ad%dbdt_%sData_%s%s", chan, sample.c_str(), variables[iv].c_str(), sel.c_str()));
    if ("bsmm" == mcsample) {
      setFilledHist(hm, kRed, kRed, 3654, 2);
    } else if ("bupsik" == mcsample)  {
      setFilledHist(hm, kGreen-2, kGreen-2, 3654, 2);
    } else if ("bspsiphi" == mcsample)  {
      setFilledHist(hm, kRed, kRed, 3654, 2);
    }
    setHist(hd, kBlack, 20, 1.4);
    if (hm->GetSumOfWeights() > 0) {
      hm->Scale(hd->GetSumOfWeights()/hm->GetSumOfWeights());
    }

    overlayAndRatio(hd, hm, year, log);
  }

}

//----------------------------------------------------------------------
string htype(TH1D *h1) {
  string hname = h1->GetName();
  vector<string> elems = split(hname, '_');
  // cout << "name: " << elems[3] << endl;
  // h1->Print();
  return elems[5];
}

// ----------------------------------------------------------------------
void rareBgStack(string type = "Hh") {
  bool preliminary(false);
  vector<pair<string, TFile*> > files;
  files.push_back(make_pair("bdt_28_2011s01_chan0", TFile::Open("s01/stacks-bdt_28_2011s01.root")));
  files.push_back(make_pair("bdt_21_2011s01_chan1", TFile::Open("s01/stacks-bdt_21_2011s01.root")));
  files.push_back(make_pair("bdt_35_2012s01_chan0", TFile::Open("s01/stacks-bdt_35_2012s01.root")));
  files.push_back(make_pair("bdt_32_2012s01_chan1", TFile::Open("s01/stacks-bdt_32_2012s01.root")));
  files.push_back(make_pair("bdt_30_2016BFs01_chan0", TFile::Open("s01/stacks-bdt_30_2016BFs01.root")));
  files.push_back(make_pair("bdt_30_2016BFs01_chan1", TFile::Open("s01/stacks-bdt_30_2016BFs01.root")));
  files.push_back(make_pair("bdt_31_2016GHs01_chan0", TFile::Open("s01/stacks-bdt_31_2016GHs01.root")));
  files.push_back(make_pair("bdt_38_2016GHs01_chan1", TFile::Open("s01/stacks-bdt_38_2016GHs01.root")));

  vector<pair<string, TH1D*> > mhists;

  // -- get list of histograms associated with type
  THStack *hStack0 = (THStack*)(files[0].second->Get(Form("h%s_bdt_28_2011s01_chan0", type.c_str())));
  TList *hlist = hStack0->GetHists();
  TIter hnext(hlist);
  TObject *obj = hnext();
  while (obj) {
    TH1D *h1 = (TH1D*)obj;
    if (h1) {
      string hname = h1->GetName();
      vector<string> elems = split(hname, '_');
      TH1D *h2 = (TH1D*)h1->Clone(elems[5].c_str());
      h2->Reset();
      mhists.push_back(make_pair(elems[5], h2));
    }
    hnext.Next();
    obj = *hnext;
  }

  // -- fill all histograms from all files, applying an ad-hoc scaling factor determined by eye-balling Jack's sum
  for (int ifile = 0; ifile < files.size(); ++ifile) {
    TIter next(files[ifile].second->GetListOfKeys());
    TKey *key(0);
    while ((key = (TKey*)next())) {
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("THStack")) continue;
      hStack0 = (THStack*)key->ReadObj();
      if (string::npos == string(hStack0->GetName()).find(type)) continue;
      if (string::npos == string(hStack0->GetName()).find(files[ifile].first)) continue;
      cout << "file " << files[ifile].second->GetName() << " hStack0: " << hStack0->GetName() << endl;
      TList *list = hStack0->GetHists();

      vector<TH1D*> vhists;
      TIter next(list);
      TObject *obj = next();
      cout << "obj: " << obj << endl;
      while (obj) {
	TH1D *h1 = (TH1D*)obj;
	if (h1) {
	  vhists.push_back(h1);
	}
	next.Next();
	obj = *next;
      }

      cout << "vhists.size() = " << vhists.size() << endl;
      for (int i = 0; i < vhists.size(); ++i) {
	string shist = htype(vhists[i]);
	for (int im = 0; im < mhists.size(); ++im) {
	  if (mhists[im].first == shist) {
	    cout << "adding " << vhists[i]->GetName() << " to " << mhists[im].second->GetName() << endl;
	    mhists[im].second->Add(vhists[i], 0.6);
	  }
	}
      }
    }
  }

  THStack *hstack = new THStack(Form("combination_%s", type.c_str()), "");
  vector<string> vnames, voptions;
  vector<TH1D*> vh;
  for (int im = 0; im < mhists.size(); ++im) {
    cout << "-> adding to stack " << mhists[im].second->GetName() << " with w8s = " << mhists[im].second->GetSumOfWeights() << endl;
    //    mhists[im].second->SetMinimum(0.001);
    hstack->Add(mhists[im].second);
    if (string::npos != mhists[im].first.find("bskk")) vnames.insert(vnames.begin(), "#it{B^{0}_{s}} #rightarrow #it{K^{+}K^{#minus}}");
    if (string::npos != mhists[im].first.find("bskpi")) vnames.insert(vnames.begin(), "#it{B^{0}_{s}} #rightarrow #it{K^{#minus}}#it{#pi^{+}}");
    if (string::npos != mhists[im].first.find("bspipi")) vnames.insert(vnames.begin(), "#it{B^{0}_{s}} #rightarrow #it{#pi^{+}#pi^{#minus}}");
    if (string::npos != mhists[im].first.find("bdkk")) vnames.insert(vnames.begin(), "#it{B^{0}} #rightarrow #it{K^{+}K^{#minus}}");
    if (string::npos != mhists[im].first.find("bdkpi")) vnames.insert(vnames.begin(), "#it{B^{0}} #rightarrow #it{K^{+}}#it{#pi^{#minus}}");
    if (string::npos != mhists[im].first.find("bdpipi")) vnames.insert(vnames.begin(), "#it{B^{0}} #rightarrow #it{#pi^{+}#pi^{#minus}}");
    if (string::npos != mhists[im].first.find("bdpimunu")) vnames.insert(vnames.begin(), "#it{B^{0}} #rightarrow #it{#pi^{#minus}#mu^{+}#nu}");
    if (string::npos != mhists[im].first.find("bdpimumu")) vnames.insert(vnames.begin(), "#it{B^{0}} #rightarrow #it{#pi^{0}#mu^{+}#mu^{#minus}}");
    if (string::npos != mhists[im].first.find("bupimumu")) vnames.insert(vnames.begin(), "#it{B^{+}} #rightarrow #it{#pi^{+}#mu^{+}#mu^{#minus}}");
    if (string::npos != mhists[im].first.find("bskmunu")) vnames.insert(vnames.begin(), "#it{B^{0}_{s}} #rightarrow #it{K^{#minus}}#it{#mu^{+}}#it{#nu}");
    if (string::npos != mhists[im].first.find("lbppi")) vnames.insert(vnames.begin(), "#it{#Lambda_{b}} #rightarrow #it{p}#it{#pi^{#minus}}");
    if (string::npos != mhists[im].first.find("lbpk")) vnames.insert(vnames.begin(), "#it{#Lambda_{b}} #rightarrow #it{pK^{#minus}}");
    if (string::npos != mhists[im].first.find("lbpmunu")) vnames.insert(vnames.begin(), "#it{#Lambda_{b}} #rightarrow #it{p}#it{#mu^{#minus}}#it{#nu}");
    voptions.push_back("f");
    vh.insert(vh.begin(), mhists[im].second);
  }

  c0->Clear();
  shrinkPad(0.12, 0.20, 0.05);

  hstack->Draw("hist");
  hstack->GetXaxis()->SetLimits(4.9, 5.9);
  hstack->GetXaxis()->CenterTitle(kTRUE);
  hstack->GetYaxis()->CenterTitle(kTRUE);
  hstack->GetXaxis()->SetNdivisions(-405);
  hstack->GetXaxis()->SetTitle("#it{m}_{#it{#mu^{+}#mu^{#minus}}} [GeV]");
  //  hstack->GetXaxis()->SetTitle(" #bf{#it{#Beta}}  m_{#mu^{+}#mu^{-}} [GeV]");
  //  hstack->GetXaxis()->SetLabelSize(0.045);
  hstack->GetYaxis()->SetTitle(Form("Candidates / %4.3f GeV", hstack->GetXaxis()->GetBinWidth(1)));
  //  hstack->GetYaxis()->SetLabelSize(0.045);
  hstack->GetYaxis()->ChangeLabel(1, -1, -1, -1, -1, -1, " ");
  TLegend *lSl(0);
  if (type == "Sl") {
    lSl = newLegend("Semileptonic decays", 0.50, 0.55, 0.85, 0.85, vh, vnames, voptions);
    lSl->SetTextFont(42);
    lSl->Draw();

    tl->SetTextFont(62); tl->SetTextSize(0.05);
    tl->DrawLatexNDC(0.20, 0.91, "CMS simulation");
    if (preliminary) {
      tl->SetTextFont(42);
      tl->DrawLatexNDC(0.60, 0.91, "#it{Preliminary}");
    }
    tl->SetTextFont(42);

  }
  if (type == "Hh") {
    shrinkPad(0.12, 0.20, 0.05);
    lSl = newLegend("Hadronic decays", 0.56, 0.4, 0.85, 0.85, vh, vnames, voptions);
    lSl->SetTextFont(42);
    lSl->Draw();

    hstack->GetYaxis()->SetTitleOffset(1.7);
    tl->SetTextFont(62); tl->SetTextSize(0.05);
    tl->DrawLatexNDC(0.20, 0.91, "CMS simulation");
    if (preliminary) {
      tl->SetTextFont(42);
      tl->DrawLatexNDC(0.60, 0.91, "#it{Preliminary}");
    }
    tl->SetTextFont(42);
  }
  if (type == "Bg") {
    lSl = newLegend("Rare decays", 0.50, 0.2, 0.85, 0.85, vh, vnames, voptions);
    lSl->SetTextFont(42);
    lSl->Draw();

    tl->SetTextFont(62); tl->SetTextSize(0.05);
    tl->DrawLatexNDC(0.20, 0.91, "CMS simulation");
    if (preliminary) {
      tl->SetTextFont(42);
      tl->DrawLatexNDC(0.60, 0.91, "#it{Preliminary}");
    }
    tl->SetTextFont(42);
  }


  c0->Modified();
  c0->Update();

  c0->SaveAs(Form("qnd/pub_rareBgStack%s.pdf", type.c_str()));
  c0->SaveAs(Form("qnd/pub_rareBgStack%s.png", type.c_str()));

}


// ----------------------------------------------------------------------
void makeAll() {

  vector<string> vars;
  // -- kinematic plots
  vars.push_back("tau");
  vars.push_back("muhelicity");
  vars.push_back("muon2pt");
  //  vars.push_back("muon1pt");

  // -- BDT most discriminating variables
  vars.push_back("fls3d");
  vars.push_back("alpha");
  vars.push_back("closetrk");
  //  vars.push_back("chi2dof");
  // vars.push_back("ips");
  // vars.push_back("ip");
  // vars.push_back("maxdoca");

  // vars.push_back("iso");
  // vars.push_back("m1iso");
  // vars.push_back("m2iso");
  // vars.push_back("docatrk");

  vars.push_back("bdtsel2");

  for (unsigned int ivar = 0; ivar < vars.size(); ++ivar) {
    overlay(vars[ivar], "bupsik", "s01bis/plotSbsHistograms-2016BFs01.root", "2016BF", 0);
    overlay(vars[ivar], "bupsik", "s01bis/plotSbsHistograms-2016GHs01.root", "2016GH", 0);

    overlay(vars[ivar], "bspsiphi", "s01bis/plotSbsHistograms-2016BFs01.root", "2016BF", 0);
    overlay(vars[ivar], "bspsiphi", "s01bis/plotSbsHistograms-2016GHs01.root", "2016GH", 0);
  }
  overlay("bdtsel2", "bupsik",   "s01/plotSbsHistograms-2011s01.root", "2011", 0);
  overlay("bdtsel2", "bupsik",   "s01/plotSbsHistograms-2012s01.root", "2012", 0);
  overlay("bdtsel2", "bspsiphi", "s01/plotSbsHistograms-2011s01.root", "2011", 0);
  overlay("bdtsel2", "bspsiphi", "s01/plotSbsHistograms-2012s01.root", "2012", 0);

  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2011s01.root", "2011", 0);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2012s01.root", "2012", 0);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2016BFs01.root", "2016BF", 0);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2016GHs01.root", "2016GH", 0);


  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2011s01.root", "2011", 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2012s01.root", "2012", 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2016BFs01.root", "2016BF", 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2016GHs01.root", "2016GH", 1);

  // -- add version with log scale for dimuon distributions
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2011s01.root", "2011", 0, 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2012s01.root", "2012", 0, 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2016BFs01.root", "2016BF", 0, 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2016GHs01.root", "2016GH", 0, 1);

  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2011s01.root", "2011", 1, 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2012s01.root", "2012", 1, 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2016BFs01.root", "2016BF", 1, 1);
  overlay("bdtsel2", "bmm", "s01/plotSbsHistograms-2016GHs01.root", "2016GH", 1, 1);

  //  return;

  normalization("s01/plotResults.2016BFs01.root", 0);
  normalization("s01/plotResults.2016BFs01.root", 1);
  normalization("s01/plotResults.2016GHs01.root", 0);
  normalization("s01/plotResults.2016GHs01.root", 1);

  controlsample("s01/plotResults.2016BFs01.root", 0);
  controlsample("s01/plotResults.2016BFs01.root", 1);
  controlsample("s01/plotResults.2016GHs01.root", 0);
  controlsample("s01/plotResults.2016GHs01.root", 1);

  rareBgStack("Hh");
  rareBgStack("Sl");
  rareBgStack("Bg");
}
