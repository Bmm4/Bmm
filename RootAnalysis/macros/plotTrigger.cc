#include "plotTrigger.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>


#include "TROOT.h"
#include "TStyle.h"
#include "TKey.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TF1.h"
#include "TFitResult.h"

#include "common/dataset.hh"
#include "common/util.hh"
#include "common/Lumi.hh"

ClassImp(plotTrigger)

using namespace std;

// ----------------------------------------------------------------------
plotTrigger::plotTrigger(string dir, string files, string cuts, string setup, int year): plotClass(dir, files, cuts, setup, year) {
  plotClass::loadFiles(files);
  plotTrigger::loadFiles(files);

  changeSetup(dir, "plotTrigger", setup);
  init();

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

  if (2011 == fYear) {
    fLargeRuns.push_back(165993);
    fLargeRuns.push_back(166950);
    fLargeRuns.push_back(167898);
    fLargeRuns.push_back(172822);
    fLargeRuns.push_back(172868);
    fLargeRuns.push_back(177730);
    fLargeRuns.push_back(178100);
  } else if (2012 == fYear) {
    fLargeRuns.push_back(194050);
    fLargeRuns.push_back(194912);
    fLargeRuns.push_back(195552);
    fLargeRuns.push_back(199435);
    fLargeRuns.push_back(201191);
    fLargeRuns.push_back(201278);
    fLargeRuns.push_back(202178);
    fLargeRuns.push_back(203002);
    fLargeRuns.push_back(205344);
    fLargeRuns.push_back(206745);
    fLargeRuns.push_back(207099);
    fLargeRuns.push_back(207231);
    fLargeRuns.push_back(207454);
    fLargeRuns.push_back(207905);
  } else if (2016 == fYear) {
    fLargeRuns.push_back(273725); // B
    fLargeRuns.push_back(274241);
    fLargeRuns.push_back(274316);
    fLargeRuns.push_back(274335);
    fLargeRuns.push_back(274388);
    fLargeRuns.push_back(274422);
    fLargeRuns.push_back(274968);
    fLargeRuns.push_back(275001);
    fLargeRuns.push_back(275125);
    fLargeRuns.push_back(275310);
    fLargeRuns.push_back(275376);
    fLargeRuns.push_back(275782); // C
    fLargeRuns.push_back(275836);
    fLargeRuns.push_back(275847);         // very high
    fLargeRuns.push_back(275890);
    fLargeRuns.push_back(276242);         // very low
    fLargeRuns.push_back(276282);
    fLargeRuns.push_back(276363);
    fLargeRuns.push_back(276437); // D
    fLargeRuns.push_back(276501);
    fLargeRuns.push_back(276525);
    fLargeRuns.push_back(276542);
    fLargeRuns.push_back(276581);
    fLargeRuns.push_back(276655);
    fLargeRuns.push_back(276776);
    fLargeRuns.push_back(276811);
    fLargeRuns.push_back(276831);
    fLargeRuns.push_back(276870); // E
    fLargeRuns.push_back(276950);
    fLargeRuns.push_back(277096);
    fLargeRuns.push_back(277194);
    fLargeRuns.push_back(277305);        // catastrophic run
    fLargeRuns.push_back(278167); // F
    fLargeRuns.push_back(278308);
    fLargeRuns.push_back(278406);
    fLargeRuns.push_back(278509);
    fLargeRuns.push_back(278808);
    fLargeRuns.push_back(278820); // G
    fLargeRuns.push_back(278822);
    fLargeRuns.push_back(278969);
    fLargeRuns.push_back(279694);
    fLargeRuns.push_back(279931);
    fLargeRuns.push_back(279975);
    fLargeRuns.push_back(280249);
    fLargeRuns.push_back(280385);
    fLargeRuns.push_back(281693); // H
    fLargeRuns.push_back(281797);
    fLargeRuns.push_back(281976);
    fLargeRuns.push_back(282037);
    fLargeRuns.push_back(282092);
    fLargeRuns.push_back(282735);
    fLargeRuns.push_back(282814);
    fLargeRuns.push_back(283270);
    fLargeRuns.push_back(283408);
    fLargeRuns.push_back(283478);
    fLargeRuns.push_back(283865);
    fLargeRuns.push_back(283946);
  } else if (2016 == fYear) {
    fLargeRuns.push_back(297359);
    fLargeRuns.push_back(297411);
    fLargeRuns.push_back(297425);
    fLargeRuns.push_back(297429);
    fLargeRuns.push_back(297430);
    fLargeRuns.push_back(297432);
    fLargeRuns.push_back(297433);
    fLargeRuns.push_back(297469);
    fLargeRuns.push_back(297474);
    fLargeRuns.push_back(297483);
    fLargeRuns.push_back(297484);
    fLargeRuns.push_back(297485);
    fLargeRuns.push_back(297486);
    fLargeRuns.push_back(297487);
    fLargeRuns.push_back(297488);
    fLargeRuns.push_back(297503);
    fLargeRuns.push_back(297504);
    fLargeRuns.push_back(297505);
    fLargeRuns.push_back(297618);
    fLargeRuns.push_back(297670);
    fLargeRuns.push_back(297671);
    fLargeRuns.push_back(297674);
    fLargeRuns.push_back(297675);
  }


}


// ----------------------------------------------------------------------
plotTrigger::~plotTrigger() {

}


// ----------------------------------------------------------------------
void plotTrigger::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotTrigger::makeAll(string what) {

  if (what == "fill") {
    runStudy("bmmData", "fill");
    runStudy("bupsikData", "fill");
  }

  if (what == "ana") {
    runStudy("bmmData", "ana");
  }

}


// ----------------------------------------------------------------------
void plotTrigger::plotL1Seeds(std::string dsname) {
  if (string::npos != dsname.find("bupsik")) fMode = BU2JPSIKP;
  if (string::npos != dsname.find("bspsiphi")) fMode = BS2JPSIPHI;
  if (string::npos != dsname.find("bdpsikstar")) fMode = BD2JPSIKSTAR;

  fSample = dsname;
  string dir = "candAnaMuMu";
  if (string::npos != fSample.find("bspsiphi")) {
    fMode = BS2JPSIPHI;
    dir  = "candAnaBs2JpsiPhi";
  } else if (string::npos != fSample.find("bupsik")) {
    fMode = BU2JPSIKP;
    dir  = "candAnaBu2JpsiK";
  } else if (string::npos != fSample.find("bdpsikstar")) {
    fMode = BD2JPSIKSTAR;
    dir  = "candAnaBd2JpsiKstar";
  }

  TH1D *h1(0);
  fHistFile = TFile::Open(fHistFileName.c_str(), "");
  if (fHistFile) {
    h1 = (TH1D*)fHistFile->Get(Form("h_%s_%s", "s0", dsname.c_str()));
    fHistFile->Close();
  }
  if (!h1) {
    // -- create histograms
    cout << "fHistFile: " << fHistFileName;
    fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
    cout << " opened " << endl;


    TTree *t = getTree(fSample, dir);
    if (0 == t) {
      cout << "tree for sample = " << fSample << " not found" << endl;
      return;
    }
    //  bookHist(fSample);
    fpHL1s0 = new TH1D(Form("h_%s_%s", "s0", dsname.c_str()), Form("h_%s_%s", "s0", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s1 = new TH1D(Form("h_%s_%s", "s1", dsname.c_str()), Form("h_%s_%s", "s1", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s2 = new TH1D(Form("h_%s_%s", "s2", dsname.c_str()), Form("h_%s_%s", "s2", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s3 = new TH1D(Form("h_%s_%s", "s3", dsname.c_str()), Form("h_%s_%s", "s3", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s4 = new TH1D(Form("h_%s_%s", "s4", dsname.c_str()), Form("h_%s_%s", "s4", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s5 = new TH1D(Form("h_%s_%s", "s5", dsname.c_str()), Form("h_%s_%s", "s5", dsname.c_str()), 7236, 273150, 280385);
    fpHL1s6 = new TH1D(Form("h_%s_%s", "s6", dsname.c_str()), Form("h_%s_%s", "s6", dsname.c_str()), 7236, 273150, 280385);
    fpHL1All = new TH1D(Form("h_%s_%s", "All", dsname.c_str()), Form("h_%s_%s", "All", dsname.c_str()), 7236, 273150, 280385);

    setupTree(t, fSample);
    fCds = fDS[fSample];
    loopOverTree(t, 4);

    fpHL1s0->Write();
    fpHL1s1->Write();
    fpHL1s2->Write();
    fpHL1s3->Write();
    fpHL1s4->Write();
    fpHL1s5->Write();
    fpHL1s6->Write();
    fpHL1All->Write();

    fHistFile->Close();
  }

  fHistFile = TFile::Open(fHistFileName.c_str(), "");
  TH1D *hAll = (TH1D*)fHistFile->Get(Form("h_%s_%s", "All", dsname.c_str()));
  TH1D *hSeed(0), *hEff(0);
  gPad->SetGridx(1);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  vector<string> sseed;
  sseed.push_back("L1_DoubleMu0er1p6_dEtaMax1p8");
  sseed.push_back("L1_DoubleMu0er1p6_dEta_Max1p8_OS");
  sseed.push_back("L1_DoubleMu0er1p4_dEta_Max1p8_OS");
  sseed.push_back("L1_DoubleMu_10_0_dEta_Max1p8");
  sseed.push_back("L1_DoubleMu_11_4");
  sseed.push_back("L1_DoubleMu_12_5");

  for (int iseed = 0; iseed < 6; ++iseed) {
    hSeed = (TH1D*)fHistFile->Get(Form("h_%s_%s", Form("s%d", iseed), dsname.c_str()));
    hEff = (TH1D*)hSeed->Clone(Form("eff_%s", hSeed->GetName())); hEff->Sumw2(); hEff->Reset();
    hEff->SetTitle((sseed[iseed] + Form("(%s)", fDS[dsname]->fName.c_str())).c_str());
    cout << "hEff title = " << hEff->GetTitle() << endl;
    hEff->Divide(hSeed, hAll, 1., 1., "b");
    setTitles(hEff, "run", "fraction");
    hEff->Draw();
    savePad(Form("effL1Seed_%d_%s.pdf", iseed, dsname.c_str()));

  }
}



// ----------------------------------------------------------------------
void plotTrigger::plotTOSHistory(std::string dsname,unsigned int runMin, unsigned int runMax) {
  //if (string::npos != dsname.find("bupsik")) fMode = BU2JPSIKP;

  // fSample = dsname;
  // string dir = "candAnaMuMu";
  // if (string::npos != fSample.find("bupsik")) {
  //   fMode = BU2JPSIKP;
  //   dir  = "candAnaBu2JpsiK";
  // }

  if ( runMin>runMax ) {cout << "The maximum run is smaller than the minimum." << endl;return;}

  setup(dsname);
  fSample = dsname;
  zone();


  TTree *t = getTree(fSample, fTreeDir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  }
  setupTree(t, fSample);

  unsigned int histoBins = runMax - runMin;
  TH2D *h = new TH2D("h","htitle",histoBins,runMin,runMax,40,4.8,6.0);
  unsigned int tEntries = t->GetEntriesFast();
  cout << "Found " << tEntries << " entries to loop over." << endl;

  int json_tmp_run(0);
  for (unsigned int i=0;i<tEntries;i++)
    {
      t->GetEntry(i);
      //print out all the runs in the json file
      if ( fb.json )
	{
	  if ( json_tmp_run != fb.run )
	    {
	      json_tmp_run = fb.run;
	      cout << "passed json: run: " << fb.run << endl;
	    }
	}
      if ( fb.hlt1 && fb.tos && (fb.chan==0 || fb.chan==1) && fb.json && fb.pt>6 && fb.m1pt>4 && fb.m2pt>4 && fb.alpha < 0.05 && fb.fls3d>10 && fb.m1gmid && fb.m2gmid )
	{
	  // cout << "filled histo with: " << fb.run << "/" << fb.m << endl;
	  h->Fill(fb.run,fb.m);
	  // break;
	}
    }

  // c0->cd(1);
  // //h->Draw("colz");
  // t->Draw("m:run");
  // c0->cd(2);

  c0->cd();
  //used later for TGraph
  std::vector<double> vruns;
  std::vector<double> integrals;

  TAxis *x = h->GetXaxis();
  Double_t run[histoBins];
  x->GetLowEdge(run);
  int xbins = x->GetNbins();
  cout << "Generated " << xbins << " bins." << endl;
  int ybins = h->GetYaxis()->GetNbins();

  //get the lumi per run
  Lumi *pl = new Lumi("../common/json/json_DCSONLY.lumi");
  //cout << "lumi = " << pl->lumi(297723) << endl;

  //loop over the TH2 histo
  for (int i=0;i<=xbins;i++)
    {
      //take care of the binning...
      int currentRun = run[i]-1;
      // cout << "edge: " << run[i] << endl;
      double lumi = pl->lumi(currentRun);
      double integral = h->Integral(i,i,1,ybins);
      // if (run[i]>297500 && run[i]<297600)
      // 	{
      // 	  cout << "run/integral: " << run[i] << "/" << integral << " (" << h->GetBinContent(i,int(ybins/2)) << endl;
      // 	}

      //suppress entries == 0 or lumi == 0
      if (integral>0 && lumi >0)
	{
	  double result = integral/lumi;
	  //cout << "lumi: " << lumi << endl;
	  cout << "run: integral/lumi:  " << currentRun << "/" << integral << "/" << lumi << " == " << result << endl;
	  //Filled together. Must be of same size
	  vruns.push_back(currentRun);
	  integrals.push_back(result);
	}
      else if ( integral>0 )
	{cout << "Rejected run: " << currentRun << " with " << integral << " events because of 0 lumi." << endl;}
    }

  if ( vruns.size() != integrals.size() )
    {cout << "ERROR: Mismatch of the result vector sizes." << endl;return;}
  double *xData = &vruns[0];
  double *yData = &integrals[0];
  TGraph *gg = new TGraph(vruns.size(),xData,yData);
  gg->SetTitle(Form("%s;run number;#events/lumi [pb]",dsname.c_str()));
  gg->GetYaxis()->SetTitleOffset(0.7);
  gg->Draw("a*");


  return;
}


// ----------------------------------------------------------------------
void plotTrigger::plotTisEfficiency(string dsname) {

  // -- read histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str());
  cout << " opened " << endl;

  vector<string> vds;

  TIter next(fHistFile->GetListOfKeys());
  TKey *key(0);
  TH1D *hpass(0), *hnorm(0);
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;

    if (TString(key->GetName()).Contains("_norm_")) {
      string hname = key->GetName();
      replaceAll(hname, "h_norm_", "");
      vds.push_back(hname);
    }
  }

  string name;
  for (unsigned int i = 0; i < vds.size(); ++i) {
    if (dsname != "all" && vds[i] != dsname) continue;
    hpass = (TH1D*)(fHistFile->Get(Form("h_pass_%s", vds[i].c_str())));
    hnorm = (TH1D*)(fHistFile->Get(Form("h_norm_%s", vds[i].c_str())));

    hnorm->Draw();
    hpass->Draw("samee");

    name = vds[i];
    replaceAll(name, "Data", " Charmonium");
    replaceAll(name, "SingleMuon", " BLA");
    replaceAll(name, "BLA", " SingleMuon");

    tl->SetTextSize(0.03);
    tl->DrawLatexNDC(0.16, 0.92, name.c_str());
    tl->DrawLatexNDC(0.60, 0.92, Form("#varepsilon = %d/%d = %3.2f #pm %3.2f",
				     static_cast<int>(hpass->GetSumOfWeights()),
				     static_cast<int>(hnorm->GetSumOfWeights()),
				     hpass->GetSumOfWeights()/hnorm->GetSumOfWeights(),
				     dEff(static_cast<int>(hpass->GetSumOfWeights()),
					  static_cast<int>(hnorm->GetSumOfWeights())
					  ))
		     );

    savePad(Form("eff-%s.pdf", vds[i].c_str()));
  }



}


// ----------------------------------------------------------------------
void plotTrigger::runTisEfficiency(string dsname) {

  if (string::npos != dsname.find("bupsik")) fMode = BU2JPSIKP;
  if (string::npos != dsname.find("bspsiphi")) fMode = BS2JPSIPHI;
  if (string::npos != dsname.find("bdpsikstar")) fMode = BD2JPSIKSTAR;

  // -- dump histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;

  fSample = dsname;
  string dir = "candAnaBu2JpsiK";

  TTree *t = getTree(fSample, dir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  }
  //  bookHist(fSample);
  fpHnorm = new TH1D(Form("h_%s_%s", "norm", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);
  fpHpass = new TH1D(Form("h_%s_%s", "pass", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);


  setupTree(t, fSample);
  fCds = fDS[fSample];
  loopOverTree(t, 1);


  fpHnorm->Write();
  fpHpass->Write();
  fHistFile->Close();

}


// ----------------------------------------------------------------------
void plotTrigger::refTrgEfficiency(string selection, string dsname) {



  zone(2,2);

  string smode("bupsikMc");
  setup(smode);
  TTree *t = getTree(smode, fTreeDir);
  if (0 == t) {
    cout << "tree for sample = " << fSample << " not found" << endl;
    return;
  }
  setupTree(t, smode);

  fvHists.insert(make_pair(string("h_MC"), new TH1D("h_MC", "h_MC", 40, 4.8, 6.0)));
  fvHists.insert(make_pair(string("h_MC_hlt"), new TH1D("h_MC_hlt", "h_MC_hlt", 40, 4.8, 6.0)));
  fvHists.insert(make_pair(string("h_rt"), new TH1D("h_rt", "h_rt", 40, 4.8, 6.0)));
  fvHists.insert(make_pair(string("h_rt_hlt"), new TH1D("h_rt_hlt", "h_rt_hlt", 40, 4.8, 6.0)));
  fvHists.insert(make_pair(string("h_ps"), new TH1D("h_ps", "h_ps", 40, 0., 20.0)));

  gStyle->SetOptStat(0);

  double nNorm(0.), nPass(0.), effMc(0.), effRt(0.), effMcE(0.), effRtE(0.);
  string fitopt("lm");

  // -- basic MC HLT efficiency
  string tselection = selection;

  //  t->Draw("m >> h_0", tselection.c_str());
  //  tselection = selection + " && hlt1 && tos";
  //  t->Draw("m >> h_0_hlt", tselection.c_str());

  fSample = "h_MC";
  fDBX = false;
  loopOverTree(t, 1);


  TF1 *f1 = fIF->bupsik(fvHists["h_MC"]);
  c0->cd(1);
  fvHists["h_MC"]->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  f1->SetParameter(7, 0.);
  f1->SetParameter(8, 0.);
  f1->SetParameter(9, 0.);
  nNorm = f1->Integral(5.1, 5.4)/fvHists["h_MC"]->GetBinWidth(1);
  tl->SetTextSize(0.02);
  tl->DrawLatexNDC(0.25, 0.96, Form("sel0: %s", selection.c_str()));
  tl->SetTextSize(0.05);
  tl->DrawLatexNDC(0.2, 0.80, "MC");
  savePad("h_0.pdf");
  delete f1;

  c0->cd(2);
  fIF->limitPar(1, 5.1, 5.5);
  f1 = fIF->bupsik(fvHists["h_MC_hlt"]);
  fvHists["h_MC_hlt"]->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  f1->SetParameter(7, 0.);
  f1->SetParameter(8, 0.);
  f1->SetParameter(9, 0.);
  nPass = f1->Integral(5.1, 5.4)/fvHists["h_MC_hlt"]->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "sel. 0 && HLT");
  tl->DrawLatexNDC(0.2, 0.80, "MC");
  effMc = nPass/nNorm;
  effMcE = dEff(static_cast<int>(nPass), static_cast<int>(nNorm));
  tl->DrawLatexNDC(0.6, 0.50, Form("#varepsilon = %4.3f#pm %4.3f", effMc, effMcE));
  savePad("h_0_hlt.pdf");
  delete f1;
  cout << "==> efficiency = " << nPass << "/" << nNorm << " = " << nPass/nNorm << endl;


  smode = "bupsikData";
  setup(smode);
  t = getTree(smode, fTreeDir);
  fSample = "h_rt";
  fDBX = true;
  loopOverTree(t, 1);

  c0->cd(3);

  f1 = fIF->bupsik(fvHists["h_rt"]);
  fvHists["h_rt_hlt"]->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  f1->SetParameter(7, 0.);
  f1->SetParameter(8, 0.);
  f1->SetParameter(9, 0.);
  nNorm = f1->Integral(5.1, 5.4)/fvHists["h_rt_hlt"]->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "sel. 0 && reftrg");
  tl->DrawLatexNDC(0.2, 0.80, "2017B");
  savePad("h_rt.pdf");
  delete f1;

  c0->cd(4);
  f1 = fIF->bupsik(fvHists["h_rt_hlt"]);
  fvHists["h_rt_hlt"]->Fit(f1, fitopt.c_str());
  f1->SetParameter(5, 0.);
  f1->SetParameter(6, 0.);
  f1->SetParameter(7, 0.);
  f1->SetParameter(8, 0.);
  f1->SetParameter(9, 0.);
  nPass = f1->Integral(5.1, 5.4)/fvHists["h_rt_hlt"]->GetBinWidth(1);
  tl->DrawLatexNDC(0.2, 0.96, "sel. 0 && reftrg && HLT");
  tl->DrawLatexNDC(0.2, 0.80, "2017B");

  double ps = fvHists["h_ps"]->GetMean();
  double psE = fvHists["h_ps"]->GetRMS();
  effRt = ps*nPass/nNorm;
  effRtE = dEff(static_cast<int>(nPass), static_cast<int>(nNorm));
  effRtE = ps*effRtE;

  tl->DrawLatexNDC(0.6, 0.50, Form("#varepsilon = %4.3f#pm %4.3f", effRt, effRtE));
  tl->DrawLatexNDC(0.6, 0.44, Form("ps = %3.1f#pm %3.2f", ps, psE));
  savePad("h_rt_hlt.pdf");
  delete f1;

  savePad("reftrg-efficiencies.pdf");

  cout << "==> PRESCALE = " << ps << " +/- " << psE << endl;
  cout << "==> efficiency = " << nPass << "/" << nNorm << " = " << nPass/nNorm << endl;
  cout << "==> cuts: " << selection << endl;
  cout << "==> efficiencies: MC = " << Form("%4.2f +/- %4.2f", effMc, effMcE)
       << "; ref trigger = " << Form("%4.2f +/- %4.2f", effRt, effRtE)
       << endl;

}



// ----------------------------------------------------------------------
void plotTrigger::runStudy(string ds, string what) {
  int runStart(294927), runEnd(298000);
  int runBins = runEnd - runStart + 1;

  setup(ds);
  fSample = ds;
  fCds = fDS[fSample];
  if (what == "fill") {
    cout << "plotTrigger::runStudy(" << ds << "), fSample = " << fSample << ", fTreeDir = " << fTreeDir << endl;
    fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
    TDirectory *hDir = fHistFile->mkdir(fTreeDir.c_str());
    fHistFile->cd(fTreeDir.c_str());
    hDir = gDirectory;
    fProf.insert(make_pair(Form("npv_%s", ds.c_str()),
			   new TProfile(Form("npv_%s", ds.c_str()), Form("npv_%s", ds.c_str()), runBins, runStart, runEnd)));
    fProf.insert(make_pair(Form("pv2lip_%s", ds.c_str()),
			   new TProfile(Form("pv2lip_%s", ds.c_str()), Form("pv2lip_%s", ds.c_str()), runBins, runStart, runEnd)));
    fProf.insert(make_pair(Form("dzmin_%s", ds.c_str()),
			   new TProfile(Form("dzmin_%s", ds.c_str()), Form("dzmin_%s", ds.c_str()), runBins, runStart, runEnd)));

    TTree *t = getTree(fSample, fTreeDir);
    setupTree(t, fSample);
    loopOverTree(t, 2);

    cout << "writing output histograms: " << fYieldHists.size() << endl;
    fProf[Form("npv_%s", ds.c_str())]->SetDirectory(hDir);
    fProf[Form("npv_%s", ds.c_str())]->Write();
    delete fProf[Form("npv_%s", ds.c_str())];

    fProf[Form("pv2lip_%s", ds.c_str())]->SetDirectory(hDir);
    fProf[Form("pv2lip_%s", ds.c_str())]->Write();
    delete fProf[Form("pv2lip_%s", ds.c_str())];

    fProf[Form("dzmin_%s", ds.c_str())]->SetDirectory(hDir);
    fProf[Form("dzmin_%s", ds.c_str())]->Write();
    delete fProf[Form("dzmin_%s", ds.c_str())];

    fProf.clear();

    for (map<string, TH2D*>::iterator it = fHLs0.begin(); it != fHLs0.end(); ++it) {
      cout << "name = " << it->first << " hist " << it->second->GetName() << endl;
      it->second->SetDirectory(hDir);
      it->second->Write();
      delete it->second;
    }
    fHLs0.clear();

    for (map<string, TH2D*>::iterator it = fHLs1.begin(); it != fHLs1.end(); ++it) {
      cout << "name = " << it->first << " hist " << it->second->GetName() << endl;
      it->second->SetDirectory(hDir);
      it->second->Write();
      delete it->second;
    }
    fHLs1.clear();

    fHistFile->Close();
    return;
  } else {
    fHistFile = TFile::Open(fHistFileName.c_str());
  }

  Lumi *lumi(0);
  if (2016 == fYear) {
    if (string::npos != ds.find("psi")) {
      //	lumi = new Lumi("../common/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
      lumi = new Lumi("../common/json/2016-HLT_DoubleMu4_3_Jpsi_Displaced.lumi");
    } else {
      lumi = new Lumi("../common/json/2016-HLT_DoubleMu4_3_Bs.lumi");
    }
  } else if (2017 == fYear) {
    lumi = new Lumi("../common/json/json_DCSONLY.lumi");
  } else if (2012 == fYear) {
    lumi = new Lumi("../common/json/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.lumi");
  } else if (2011 == fYear) {
    lumi = new Lumi("../common/json/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_MuonPhys_v2.lumi");
  }


  shrinkPad(0.1, 0.15, 0.12, 0.1);
  TH1D *hl = new TH1D("lostLs", "lost lumi sections", runBins, runStart, runEnd);
  setTitles(hl, "run", Form("fractional L1 seed loss  (%s) / pb^{-1}", fDS[ds]->fName.c_str()), 0.05, 1.05, 1.5);
  TH1D *h0y = new TH1D("yields", "yields low pT", runBins, runStart, runEnd);
  setTitles(h0y, "run", Form("events with %s cand / pb^{-1}", fDS[ds]->fName.c_str()), 0.05, 1.05, 1.5);
  TH1D *h1y = new TH1D("yields", "yields high pT", runBins, runStart, runEnd);
  setTitles(h1y, "run", Form("events with %s cand / pb^{-1}", fDS[ds]->fName.c_str()), 0.05, 1.05, 1.5);
  double ymax0(1200.), ymax1(600.);
  if (string::npos != ds.find("bmm")) {
    ymax0 = 20.;
  }

  TH2D *hc0b = new TH2D("hc0b", "", 50, 0., 50., 50, 0., ymax0);
  TH2D *hc0c = new TH2D("hc0c", "", 50, 0., 50., 50, 0., ymax0);
  TH2D *hc0d = new TH2D("hc0d", "", 50, 0., 50., 50, 0., ymax0);
  TH2D *hc0e = new TH2D("hc0e", "", 50, 0., 50., 50, 0., ymax0);
  TH2D *hc0f = new TH2D("hc0f", "", 50, 0., 50., 50, 0., ymax0);
  TH2D *hc0g = new TH2D("hc0g", "", 50, 0., 50., 50, 0., ymax0);
  TH2D *hc0h = new TH2D("hc0h", "", 50, 0., 50., 50, 0., ymax0);

  TH2D *hc1b = new TH2D("hc1b", "", 50, 0., 50., 50, 0., ymax1);
  TH2D *hc1c = new TH2D("hc1c", "", 50, 0., 50., 50, 0., ymax1);
  TH2D *hc1d = new TH2D("hc1d", "", 50, 0., 50., 50, 0., ymax1);
  TH2D *hc1e = new TH2D("hc1e", "", 50, 0., 50., 50, 0., ymax1);
  TH2D *hc1f = new TH2D("hc1f", "", 50, 0., 50., 50, 0., ymax1);
  TH2D *hc1g = new TH2D("hc1g", "", 50, 0., 50., 50, 0., ymax1);
  TH2D *hc1h = new TH2D("hc1h", "", 50, 0., 50., 50, 0., ymax1);

  fHistFile->cd(fTreeDir.c_str());
  TProfile *hp = (TProfile*)gDirectory->Get(Form("npv_%s", ds.c_str()));
  TKey *key(0);
  TIter next = gDirectory->GetListOfKeys();
  while ((key = (TKey*)next())) {
    if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH2D")) continue;
    if (TString(key->GetName()).Contains(Form("ls0_%s_", ds.c_str()))) {
      string hname = key->GetName();
      TH2D *h0 = (TH2D*)((TH2D*)gDirectory->Get(hname.c_str()));
      replaceAll(hname, "ls0", "ls1");
      TH2D *h1 = (TH2D*)((TH2D*)gDirectory->Get(hname.c_str()));
      TH1D *h0l = h0->ProjectionX(Form("%s_x0", hname.c_str()), 1,1); h0l->SetLineColor(kBlue);
      TH1D *h1l = h1->ProjectionX(Form("%s_x1", hname.c_str()), 1,1); h1l->SetLineColor(kRed);
      h0l->Rebin(10);
      h1l->Rebin(10);
      replaceAll(hname, Form("ls1_%s_", ds.c_str()), "");
      int run = atoi(hname.c_str());
      string rr = era(run);
      double hltLumi = lumi->lumi(run);
      cout << "hp = " << hp  << " gDirectory = " << gDirectory->GetName() << endl;
      if (hltLumi > 10.) {
	if ("B" == rr) {
	  hc0b->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1b->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("C" == rr) {
	  hc0c->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1c->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("D" == rr) {
	  hc0d->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1d->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("E" == rr) {
	  hc0e->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1e->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("F" == rr) {
	  hc0f->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1f->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("G" == rr) {
	  hc0g->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1g->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	} else if ("H" == rr) {
	  hc0h->Fill(hp->GetBinContent(hp->FindBin(run)), h0l->Integral()/hltLumi);
	  hc1h->Fill(hp->GetBinContent(hp->FindBin(run)), h1l->Integral()/hltLumi);
	}
	h0y->SetBinContent(h1y->FindBin(run), h0l->Integral()/hltLumi);
	h1y->SetBinContent(h1y->FindBin(run), h1l->Integral()/hltLumi);

	cout << "run = " << run
	     << " lumi = " << hltLumi
	     << " nPV: " << hp->GetBinContent(hp->FindBin(run))
	     << " yield(0) = " <<  h0l->Integral()
	     << " yield(0)*pb = " <<  h0l->Integral()/hltLumi
	     << " yield(1) = " <<  h1l->Integral()
	     << " yield(1)*pb = " <<  h1l->Integral()/hltLumi
	     << endl;
      }
      // -- skip initial runs for the missed ls count
      if (run > 274000) {
	int lsCnt(0), lsTot(0);
	double ymax = h1l->GetMaximum();
	for (int ib = 0; ib < h0l->GetNbinsX(); ++ib) {
	  if (h1l->GetBinContent(ib) > 0.05*ymax) {
	    ++lsTot;
	    if (h0l->GetBinContent(ib) < h1l->GetBinContent(ib)) ++lsCnt;
	  }
	}
	if (lsCnt>1) {
	  cout << "XXXX Run " << run << " missed ls fraction: " << static_cast<double>(lsCnt)/lsTot << endl;
	  hl->SetBinContent(hl->FindBin(run), static_cast<double>(lsCnt)/lsTot);
	}
      }
      TH1D *hr = (TH1D*)h1l->Clone("h1r");
      hr->Divide(h0l, h1l, 1., 1.);
      hr->Scale(0.6*h0l->GetMaximum()/hr->GetMaximum());
      hr->SetLineColor(kBlack);
      h0l->Draw("hist");
      h1l->Draw("samehist");
      hr->Draw("samehist");
      tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.2, 0.93, h0l->GetName());
      tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.2, 0.90, h1l->GetName());
      //this produces yields vs LS for all runs:   savePad(Form("runStudy-%d-%s-l1seeds-%s.pdf", fYear, ds.c_str(), hname.c_str()));
    }
  }

  hl->SetMinimum(0.);
  hl->Draw();
  savePad(Form("runStudy-%d-%s-missedLS.pdf", fYear, ds.c_str()));
  h0y->SetMinimum(0.);
  h0y->Draw("e");
  savePad(Form("runStudy-%d-%s-yields-DoubleMu0.pdf", fYear, ds.c_str()));
  h1y->SetMinimum(0.);
  h1y->Draw("e");
  savePad(Form("runStudy-%d-%s-yields-highPt.pdf", fYear, ds.c_str()));

  cout << "MAXIMUM: " << hc0b->GetMaximum() << endl;
  //  hc0b->SetMaximum(1.2*hc0b->GetMaximum());
  hc0b->SetMarkerColor(kBlack);     hc0b->Draw();
  hc0c->SetMarkerColor(kGreen+1);   hc0c->Draw("same");
  hc0d->SetMarkerColor(kRed+1);     hc0d->Draw("same");
  hc0e->SetMarkerColor(kRed+2);     hc0e->Draw("same");
  hc0f->SetMarkerColor(kRed+3);     hc0f->Draw("same");
  hc0g->SetMarkerColor(kBlue+1);    hc0g->Draw("same");
  hc0h->SetMarkerColor(kBlue+3);    hc0h->Draw("same");

  tl->SetTextColor(kBlack);   tl->DrawLatexNDC(0.3, 0.92, "DoubleMu0 seeds");
  tl->SetTextColor(kBlack);   tl->DrawLatexNDC(0.7, 0.80, "2016B");
  tl->SetTextColor(kGreen+1); tl->DrawLatexNDC(0.7, 0.72, "2016C");
  tl->SetTextColor(kRed+1);   tl->DrawLatexNDC(0.7, 0.64, "2016D");
  tl->SetTextColor(kRed+2);   tl->DrawLatexNDC(0.7, 0.56, "2016E");
  tl->SetTextColor(kRed+3);   tl->DrawLatexNDC(0.7, 0.48, "2016F");
  tl->SetTextColor(kBlue+1);  tl->DrawLatexNDC(0.7, 0.40, "2016G");
  tl->SetTextColor(kBlue+3);  tl->DrawLatexNDC(0.7, 0.32, "2016H");

  savePad(Form("runStudy-%d-%s-yieldVsNpv-Seed0-allEras.pdf", fYear, ds.c_str()));

  hc1b->SetMarkerColor(kBlack);     hc1b->Draw();
  hc1c->SetMarkerColor(kGreen+1);   hc1c->Draw("same");
  hc1d->SetMarkerColor(kRed+1);     hc1d->Draw("same");
  hc1e->SetMarkerColor(kRed+2);     hc1e->Draw("same");
  hc1f->SetMarkerColor(kRed+3);     hc1f->Draw("same");
  hc1g->SetMarkerColor(kBlue+1);    hc1g->Draw("same");
  hc1h->SetMarkerColor(kBlue+3);    hc1h->Draw("same");

  tl->SetTextColor(kBlack);   tl->DrawLatexNDC(0.3, 0.92, "DoubleMu_1X_Y seeds");
  tl->SetTextColor(kBlack);   tl->DrawLatexNDC(0.7, 0.80, "2016B");
  tl->SetTextColor(kGreen+1); tl->DrawLatexNDC(0.7, 0.72, "2016C");
  tl->SetTextColor(kRed+1);   tl->DrawLatexNDC(0.7, 0.64, "2016D");
  tl->SetTextColor(kRed+2);   tl->DrawLatexNDC(0.7, 0.56, "2016E");
  tl->SetTextColor(kRed+3);   tl->DrawLatexNDC(0.7, 0.48, "2016F");
  tl->SetTextColor(kBlue+1);  tl->DrawLatexNDC(0.7, 0.40, "2016G");
  tl->SetTextColor(kBlue+3);  tl->DrawLatexNDC(0.7, 0.32, "2016H");

  savePad(Form("runStudy-%d-%s-yieldVsNpv-Seed1-allEras.pdf", fYear, ds.c_str()));

  delete hl;
  delete h0y;
  delete h1y;
  delete hc0b;
  delete hc0c;
  delete hc0d;
  delete hc0e;
  delete hc0f;
  delete hc0g;
  delete hc0h;

  delete hc1b;
  delete hc1c;
  delete hc1d;
  delete hc1e;
  delete hc1f;
  delete hc1g;
  delete hc1h;

  fHistFile->Close();

}


// ----------------------------------------------------------------------
// -- for runStudy: 2d chan vs lumisection
void plotTrigger::loopFunction2() {

  if (fLargeRuns.end() != find(fLargeRuns.begin(), fLargeRuns.end(), fb.run)) {
  } else {
    // -- comment the following line if ALL runs should be looked at
    //    return;
  }

  char hname0[200];
  char hname1[200];
  sprintf(hname0, "ls0_%s_%d", fSample.c_str(), static_cast<int>(fb.run));
  sprintf(hname1, "ls1_%s_%d", fSample.c_str(), static_cast<int>(fb.run));

  if (0 == fHLs0.count(hname0)) {
    fHLs0.insert(make_pair(hname0, new TH2D(hname0, hname0, 3500, 0., 3500., 6, -1., 5.)));
    fHLs1.insert(make_pair(hname1, new TH2D(hname1, hname1, 3500, 0., 3500., 6, -1., 5.)));
  }


  if (!fb.json) return;
  if (!fb.hlt1) return;
  if (!fb.tos) return;
  if (fb.m1pt < 4.0) return;
  if (fb.m2pt < 4.0) return;
  if (TMath::Abs(fb.m1eta) > 1.4) return;
  if (TMath::Abs(fb.m2eta) > 1.4) return;
  if (fb.m1q * fb.m2q > 0) return;

  if ((fMode == BU2JPSIKP) || (fMode == BD2JPSIKSTAR) || (fMode == BS2JPSIPHI)) {
    if (fb.mpsi < 3.04) return;
    if (fb.mpsi > 3.15) return;

    if (fb.psipt   < 7.0) return;
    if (fb.psiprob < 0.1) return;
    if (fMode == BU2JPSIKP) {
      if (fb.kpt < 0.60) return;
    }


    if (fMode == BS2JPSIPHI) {
      if (fb.mkk   < 1.01) return;
      if (fb.mkk   > 1.03) return;
      if (fb.phidr > 0.30) return;
      if (fb.k1pt  < 0.60) return;
      if (fb.k2pt  < 0.60) return;
    }

    if (fMode == BD2JPSIKSTAR) {
      if (!fb.kstarfail) return;
      if (fb.kpt < 0.60) return;
      if (fb.pipt < 0.60) return;
      if (fb.mkpi < 0.86) return;
      if (fb.mkpi > 0.94) return;
    }

  }


  if (((fb.l1s & 0x1) == 1) || ((fb.l1s & 0x2) == 2) || ((fb.l1s & 0x4) == 4)) {
    fHLs0[hname0]->Fill(fb.ls, -0.5);
    if (fb.chan > -1) fHLs0[hname0]->Fill(fb.ls, fb.chan);
  }

  if (((fb.l1s & 8) == 8) || ((fb.l1s & 16) == 16) || ((fb.l1s & 32) == 32) || ((fb.l1s & 64) == 64)) {
    fHLs1[hname1]->Fill(fb.ls, -0.5);
    if (fb.chan > -1) fHLs1[hname1]->Fill(fb.ls, fb.chan);
  }

  fProf[Form("npv_%s", fSample.c_str())]->Fill(fb.run, fb.pvn);
  fProf[Form("pv2lip_%s", fSample.c_str())]->Fill(fb.run, TMath::Abs(fb.pv2lip));
  fProf[Form("dzmin_%s", fSample.c_str())]->Fill(fb.run, TMath::Abs(fb.dzmin));


}



// ----------------------------------------------------------------------
void plotTrigger::loopFunction1() {

  // bool goodRun(false);

  // if (fLargeRuns.end() != find(fLargeRuns.begin(), fLargeRuns.end(), fb.run)) {
  //   goodRun = true;
  //   //    cout << "found " << fb.run << " in fLargeRuns" << endl;
  // } else {
  //   // -- comment the following if you want to look at large runs only
  //   //goodRun = true;
  // }


  if (fb.chan < 0 || fb.chan > 1) return;
  if (fb.m1pt < 8.) return;
  if (fb.m2pt < 4.) return;
  if (fb.flsxy < 5.) return;
  if (fb.maxdoca > 0.08) return;
  if (fb.pchi2dof < 0.10) return;
  if (!fb.gmugmid) return;
  if (fDBX) {
    if (!fb.reftrg) return;
  }
  fvHists[fSample]->Fill(fb.m);

  if (!fb.hlt1) return;
  if (!fb.tos) return;
  string hlt = Form("%s_hlt", fSample.c_str());
  fvHists[hlt]->Fill(fb.m);
  if (fDBX) {
    fvHists["h_ps"]->Fill(fb.ps);
  }
}


// ----------------------------------------------------------------------
void plotTrigger::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
  int nentries = Int_t(t->GetEntries());
  int nbegin(0), nend(nentries);
  if (nevts > 0 && nentries > nevts) {
    nentries = nevts;
    nbegin = 0;
    nend = nevts;
  }
  if (nevts > 0 && nstart > 0) {
    nentries = nstart + nevts;
    nbegin = nstart;
    if (nstart + nevts < t->GetEntries()) {
      nend = nstart + nevts;
    } else {
      nend = t->GetEntries();
    }
  }

  nentries = nend - nstart;

  int step(1000000);
  if (nentries < 5000000)  step = 500000;
  if (nentries < 1000000)  step = 100000;
  if (nentries < 100000)   step = 10000;
  if (nentries < 10000)    step = 1000;
  if (nentries < 1000)     step = 100;
  step = 500000;
  cout << "==> plotTrigger::loopOverTree> loop over dataset " << (fCds?fCds->fName:"undefined") << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotTrigger::*pF)(void);
  if (ifunc == 1) pF = &plotTrigger::loopFunction1;
  if (ifunc == 2) pF = &plotTrigger::loopFunction2;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotTrigger::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotTrigger::loadFile loading files listed in " << files << endl;

  char buffer[1000];
  ifstream is(files.c_str());
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}

    string sbuffer = string(buffer);
    replaceAll(sbuffer, " ", "");
    replaceAll(sbuffer, "\t", "");
    if (sbuffer.size() < 1) continue;

    string::size_type m1 = sbuffer.find("lumi=");
    string stype = sbuffer.substr(5, m1-5);

    string::size_type m2 = sbuffer.find("file=");
    string slumi = sbuffer.substr(m1+5, m2-m1-6);
    string sfile = sbuffer.substr(m2+5);
    string sname("nada"), sdecay("nada");

    TFile *pF(0);
    dataset *ds(0);

    if (string::npos != stype.find("SingleMuon")) {
      // -- SingleMuon
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("bmm")) {
        sname = "bmmSingleMuon";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bupsik")) {
        sname = "bupsikSingleMuon";
        sdecay = "bupsik";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bspsiphi")) {
        sname = "bspsiphiSingleMuon";
        sdecay = "bspsiphi";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bdpsikstar")) {
        sname = "bdpsikstarSingleMuon";
        sdecay = "bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

    } else if (string::npos != stype.find("Charmonium")) {
      // -- Charmonium
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("bmm")) {
        sname = "bmmCharmonium";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bupsik")) {
        sname = "bupsikCharmonium";
        sdecay = "bupsik";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bspsiphi")) {
        sname = "bspsiphiCharmonium";
        sdecay = "bspsiphi";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bdpsikstar")) {
        sname = "bdpsikstarCharmonium";
        sdecay = "bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

    } else if (string::npos != stype.find("mc")) {
      // -- MC
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2.;

      if (string::npos != stype.find("bctopsimunu,")) {
        sname = "bcpsimunuMc";
        sdecay = "bcpsimunu";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

    } else if (string::npos != stype.find("relval")) {
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2.;
      if (string::npos != stype.find("bsmm,relval")) {
	ds->fF      = pF;
	sname   = "bsmmrelval";
      }

      if (string::npos != stype.find("bdmm,relval")) {
	ds->fF      = pF;
	sname   = "bdmmrelval";
      }

      if (string::npos != stype.find("bupsik,relval")) {
	ds->fF      = pF;
	sname   = "bupsikrelval";
      }

      if (string::npos != stype.find("bspsiphi,relval")) {
	ds->fF      = pF;
	sname   = "bspsiphirelval";
      }
    }

    if (sname != "nada") {
      ds->fLcolor = ds->fColor;
      ds->fFcolor = ds->fColor;
      ds->fName   = sdecay;
      ds->fFullName = sname;
      insertDataset(sname, ds);
    } else {
      delete ds;
    }

  }

  is.close();

  int cnt(0);
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << Form("   %20s: %70s %15s %9s %8s %8s", "Dataset name", "Filename", "LaTeX", "Lumi", "Eff", "BF") << endl;
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    // cout << it->first << endl;
    // cout << it->second->fName << endl;
    // cout << it->second->fF->GetName() << endl;
    if (!it->second->fF) {
      cout << "missing " << it->first << endl;
    }
    cout << Form("%2d %20s: %70s %15s %8.1f %.2e %.2e",
		 cnt,
		 it->first.c_str(),
		 it->second->fF->GetName(),
		 it->second->fLatexName.c_str(),
		 it->second->fLumi,
		 it->second->fFilterEff,
		 it->second->fBf
		 )
	 << endl;
    ++cnt;
  }
  cout << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
}
