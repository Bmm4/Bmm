#include "plotStuff.hh"

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

ClassImp(plotStuff)

using namespace std;

// ----------------------------------------------------------------------
plotStuff::plotStuff(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotStuff::loadFiles(files);

  changeSetup(dir, "plotStuff", setup);
  init();

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;
}


// ----------------------------------------------------------------------
plotStuff::~plotStuff() {

}



// ----------------------------------------------------------------------
void plotStuff::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotStuff::makeAll(string what) {

  if (what == "dbx") {
    changeSetup("results", "yieldstability", "");
    yieldStability("bupsikData", "HLT");
    // yieldStability("bmmData", "HLT");
    // yieldStability("bspsiphiData", "HLT");
    // yieldStability("bdpsikstarData", "HLT");
  }

  if (what == "all" || what == "yieldstability") {
    yieldStability("bupsikData", "HLT");
    yieldStability("bmmData", "HLT");
    yieldStability("bspsiphiData", "HLT");
    yieldStability("bdpsikstarData", "HLT");
  }

  if (what == "all" || what == "pvstudy") {
    pvStudy("bdmmMcOff", "&& (fl1>0.01)", "fl1");
    pvStudy("bupsikMcOff", "&&(fl1>0.01)", "fl1");

    pvStudy("bdmmMcOff", "&& (fl1>0.01) && (idx1!=idx3)", "fl1_idx");
    pvStudy("bupsikMcOff", "&& (fl1>0.01) && (idx1!=idx3)", "fl1_idx");
  }

}



// ----------------------------------------------------------------------
void plotStuff::bookHist(string dsname) {


}


// ----------------------------------------------------------------------
void plotStuff::pvStudy(string dsname, string selection, string fmod) {
  fSample = dsname;
  string dir = "candAnaBu2JpsiK";
  if (string::npos != fSample.find("mm")) {
    dir = "candAnaMuMu";
  } else if (string::npos != fSample.find("bspsiphi")) {
    dir = "candAnaBs2JpsiPhi";
  } else if (string::npos != fSample.find("bdpsikstar")) {
    dir = "candAnaBd2JpsiKstar";
  }
  TTree *T = getTree(dsname, dir, "pvstudy");

  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");

  TH1D *h1(0);
  for (int i = 0; i < fNchan; ++i) {
    h1  = new TH1D(Form("dd1_ch%d", i), Form("dist(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "dist(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dd2_ch%d", i), Form("dist(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "dist(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dd3_ch%d", i), Form("dist(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "dist(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("dx1_ch%d", i), Form("dx(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta x(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dx2_ch%d", i), Form("dx(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta x(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dx3_ch%d", i), Form("dx(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta x(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("dy1_ch%d", i), Form("dy(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta y(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dy2_ch%d", i), Form("dy(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta y(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dy3_ch%d", i), Form("dy(PVreco, PVgen) ch%d", i), 100, 0., 0.25);  setTitles(h1, "#Delta y(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("dz1_ch%d", i), Form("dz(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "#Delta z(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dz2_ch%d", i), Form("dz(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "#Delta z(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dz3_ch%d", i), Form("dz(PVreco, PVgen) ch%d", i), 100, 0., 5.0);  setTitles(h1, "#Delta z(PVreco, PVgen) [cm]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("dt1_ch%d", i), Form("#delta t(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dt2_ch%d", i), Form("#delta t(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("dt3_ch%d", i), Form("#delta t(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);

    h1  = new TH1D(Form("ds1_ch%d", i), Form("#delta t^{2D}(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t^{2D}(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("ds2_ch%d", i), Form("#delta t^{2D}(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t^{2D}(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
    h1  = new TH1D(Form("ds3_ch%d", i), Form("#delta t^{2D}(reco, gen) ch%d", i), 100, -1.e-11, 1.e-11);  setTitles(h1, "#Delta t^{2D}(reco, gen) [sec]", "", 0.07, 0.9, 1.1, 0.06);
  }

  for (int i = 0; i < 4; ++i) {
    T->Draw(Form("d1 >> dd1_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("d2 >> dd2_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("d3 >> dd3_ch%d", i), Form("chan == %d %s", i, selection.c_str()));

    T->Draw(Form("p1z-gz >> dz1_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("p2z-gz >> dz2_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("p3z-gz >> dz3_ch%d", i), Form("chan == %d %s", i, selection.c_str()));

    T->Draw(Form("p1x-gx >> dx1_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("p2x-gx >> dx2_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("p3x-gx >> dx3_ch%d", i), Form("chan == %d %s", i, selection.c_str()));

    T->Draw(Form("p1y-gy >> dy1_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("p2y-gy >> dy2_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("p3y-gy >> dy3_ch%d", i), Form("chan == %d %s", i, selection.c_str()));

    T->Draw(Form("t1-gt >> dt1_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("t2-gt >> dt2_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("t3-gt >> dt3_ch%d", i), Form("chan == %d %s", i, selection.c_str()));

    T->Draw(Form("s1-gs >> ds1_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("s2-gs >> ds2_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
    T->Draw(Form("s3-gs >> ds3_ch%d", i), Form("chan == %d %s", i, selection.c_str()));
  }

  TH1D *h2(0);
  tl->SetTextSize(0.03);
  makeCanvas(1);
  zone(4, 1, c1);
  tl->SetTextSize(0.06);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dd1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dd3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.25, 0.80, Form("r1: %5.3f/%5.3f [cm]", h1->GetMean(), h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.25, 0.73, Form("r3: %5.3f/%5.3f [cm]", h2->GetMean(), h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  savePad(Form("pvstudy-dd-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dx1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dx3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.25, 0.80, Form("r1: %5.3f/%5.3f [cm]", h1->GetMean(), h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.25, 0.73, Form("r3: %5.3f/%5.3f [cm]", h2->GetMean(), h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dx.pdf"), c1);
  savePad(Form("pvstudy-dx-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dy1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dy3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.25, 0.80, Form("r1: %5.3f/%5.3f [cm]", h1->GetMean(), h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.25, 0.73, Form("r3: %5.3f/%5.3f [cm]", h2->GetMean(), h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dy.pdf"), c1);
  savePad(Form("pvstudy-dy-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dz1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dz3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.25, 0.80, Form("r1: %5.3f/%5.3f [cm]", h1->GetMean(), h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.25, 0.73, Form("r3: %5.3f/%5.3f [cm]", h2->GetMean(), h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dz.pdf"), c1);
  savePad(Form("pvstudy-dz-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);


  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("dt1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("dt3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->SetMaximum(13.*h1->GetMaximum());
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.2, 0.84, Form("r1: %5.3f/%5.3f [ps]", 1.e12*h1->GetMean(), 1.e12*h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.2, 0.79, Form("r3: %5.3f/%5.3f [ps]", 1.e12*h2->GetMean(), 1.e12*h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dt.pdf"), c1);
  savePad(Form("pvstudy-dt-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  zone(4, 1, c1);
  for (int i = 0; i < 4; ++i) {
    c1->cd(i+1); gPad->SetLogy(1);
    shrinkPad(0.15, 0.15, 0.20, 0.08);
    h1 = ((TH1D*)gDirectory->Get(Form("ds1_ch%d", i))); setFilledHist(h1, kBlue, kBlue, 3365);
    h2 = ((TH1D*)gDirectory->Get(Form("ds3_ch%d", i))); setFilledHist(h2, kRed, kRed, 3554);
    h1->SetMaximum(13.*h1->GetMaximum());
    h1->Draw();
    h2->Draw("same");
    tl->SetTextColor(kBlue); tl->DrawLatexNDC(0.2, 0.84, Form("r1: %5.3f/%5.3f [ps]", 1.e12*h1->GetMean(), 1.e12*h1->GetRMS()));
    tl->SetTextColor(kRed);  tl->DrawLatexNDC(0.2, 0.79, Form("r3: %5.3f/%5.3f [ps]", 1.e12*h2->GetMean(), 1.e12*h2->GetRMS()));
    tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.2, 0.96, Form("Channel %d ", i));
    tl->SetTextAngle(90.); tl->SetTextColor(kBlack);tl->DrawLatexNDC(0.88, 0.20, Form("Sel: %s", selection.c_str())); tl->SetTextAngle(0.);
  }
  //  savePad(Form("pvstudy-dt2d.pdf"), c1);
  savePad(Form("pvstudy-dt2d-%s-%s.pdf", fmod.c_str(), dsname.c_str()), c1);

  fHistFile->Close();

}


// ----------------------------------------------------------------------
void plotStuff::yieldStability(string dsname, string trg) {

  fIF->fVerbose = false;

  // -- dump histograms
  cout << "fHistFile: " << fHistFileName;
  fHistFile = TFile::Open(fHistFileName.c_str(), "UPDATE");
  cout << " opened " << endl;


  fSample = dsname;
  string dir = "candAnaBu2JpsiK";
  if (string::npos != fSample.find("bmm")) {
    dir = "candAnaMuMu";
  } else if (string::npos != fSample.find("bspsiphi")) {
    dir = "candAnaBs2JpsiPhi";
  } else if (string::npos != fSample.find("bdpsikstar")) {
    dir = "candAnaBd2JpsiKstar";
  }

  // -- check whether there are any histograms pre-produced already
  bool ok = fHistFile->cd(dir.c_str());
  if (ok) {
    cout << "histograms exist already, looping over them" << endl;
    TIter next(gDirectory->GetListOfKeys());
    TKey *key(0);
    vector<int> vruns;
    int run(-1), runMin(9999999), runMax(0);
    while ((key = (TKey*)next())) {
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
      if (TString(key->GetName()).Contains("_HLT_")) {
	string hname = key->GetName();
	replaceAll(hname, "h_HLT_", "");
	run = atoi(hname.c_str());
	if (run > runMax) runMax = run;
	if (run < runMin) runMin = run;
	if (find(vruns.begin(), vruns.end(), run) == vruns.end()) {
	  vruns.push_back(run);
	  // if (run == 274244) vruns.push_back(run);
	}
	//	if (run < 278800) continue;
	if (run > 274000) break;
      }
    }

    if (vruns.size() > 0) {
      double minLumi(200.);

      Lumi lumi("../common/json/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
      cout << "analyzing runs " << runMin << " .. " <<  runMax << endl;

      vector<TH1D *> vRunHLT;
      for (unsigned int ichan = 0; ichan < fNchan; ++ichan) {
	vRunHLT.push_back(new TH1D(Form("hRun%s_%s_chan%d", trg.c_str(), dsname.c_str(), ichan),
				   Form("hRun%s_%s_chan%d", trg.c_str(), dsname.c_str(), ichan),
				   runMax-runMin+1, runMin, runMax));
	vRunHLT[ichan]->Sumw2();
      }

      double mBp(5.28), sBp(0.04), stepBp(5.145);
      double xmin(5.0), xmax(5.9), ymax(0.);
      fIF->fLo = xmin;
      fIF->fHi = xmax;
      TH2D *h2 = (TH2D*)(gDirectory->Get(Form("h_%s_%d_1", trg.c_str(), vruns[0])));
      if (!h2) {
	cout << "histogram " << Form("h_%s_%d_1", trg.c_str(), vruns[0]) << " not found!?! returning" << endl;
	return;
      }

      // -- create map of chan histograms
      int MAXPS(10);
      TH2D *h2Sum = (TH2D*)h2->Clone("h2sum"); h2Sum->Reset();
      TH1D *h1 = h2->ProjectionX("chan_0", 1,1);
      h1->SetName("chan_0"); h1->Reset();
      map<string, TH1D*> allHists;
      for (unsigned int ichan = 0; ichan < fNchan; ++ichan) {
	for (unsigned int ips = 1; ips < MAXPS; ++ips) {
	  allHists.insert(make_pair(Form("chan%d_ps%d", ichan, ips), (TH1D*)h1->Clone(Form("chan%d_ps%d", ichan, ips))));
	  allHists[Form("chan%d_ps%d", ichan, ips)]->Reset();
	  allHists.insert(make_pair(Form("chan%d", ichan), (TH1D*)h1->Clone(Form("chan%d", ichan))));
	  allHists[Form("chan%d", ichan)]->Reset();
	}
      }
      double intLumi(0.);

      for (unsigned int irun = 0; irun < vruns.size(); ++irun) {
	for (int ips = 1; ips < MAXPS; ++ips) {
	  h2 = (TH2D*)(gDirectory->Get(Form("h_%s_%d_%d", trg.c_str(), vruns[irun], ips)));
	  if (h2) {
	    for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
	      h1 = h2->ProjectionX("bla", ichan+1, ichan+1);
	      allHists[Form("chan%d", ichan)]->Add(h1);
	      delete h1;
	    }
	  }
	}
      }

      // -- now loop over all runs, splitting into lumi blocks
      vector<int> prescales;
      for (unsigned int irun = 0; irun < vruns.size(); ++irun) {
	for (int ips = 1; ips < MAXPS; ++ips) {
	  h2 = (TH2D*)(gDirectory->Get(Form("h_%s_%d_%d", trg.c_str(), vruns[irun], ips)));
	  if (h2) {
	    if (find(prescales.begin(), prescales.end(), ips) == prescales.end()) {
	      prescales.push_back(ips);
	    }
	    for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
	      h1 = h2->ProjectionX("bla", ichan+1, ichan+1);
	      allHists[Form("chan%d_ps%d", ichan, ips)]->Add(h1);
	      //	      allHists[Form("chan%d", ichan)]->Add(h1);
	      if (0 == ichan) {
		if (h1->GetSumOfWeights() > 0) cout << "chan 0: ps " << ips << ", h1: " << h1->GetSumOfWeights();
	      }
	      delete h1;
	    }
	  }
	}
	cout << endl;
	intLumi += lumi.lumi(vruns[irun]);
	cout << "adding vruns[" << irun << "] = " << vruns[irun] << " with lumi = " << lumi.lumi(vruns[irun]) << " and intLumi = " << intLumi << endl;

	if (intLumi > minLumi) {
	  cout << "==== now fitting because intLumi = " << intLumi << endl;
	  for (int ips = 1; ips < MAXPS; ++ips) {
	    for (unsigned ichan = 0; ichan < 1; ++ichan) {
	      if (0 == ichan) {
		if (allHists[Form("chan%d_ps%d", ichan, ips)]->GetSumOfWeights() > 0) {
		  cout  << " chan 0, run" << vruns[irun] << ": ps " << ips << ", h1: " << allHists[Form("chan%d_ps%d", ichan, ips)]->GetSumOfWeights();
		}
	      }
	    }
	  }
	  if (allHists[Form("chan%d", 0)]->GetSumOfWeights() > 0) {
	    cout << " combined: " << allHists[Form("chan%d", 0)]->GetSumOfWeights();
	  }
	  cout << endl;

	  for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
	    // -- fit combined ps hist to determine signal parameters
	    h1 = allHists[Form("chan%d", ichan)];

	    TF1 *f1 = fIF->expoErrGauss(h1, mBp, sBp, stepBp);
	    f1->SetParLimits(8, 0., 1.e5);
	    h1->SetMinimum(0.);
	    cout << "========> Fitting combined ps for channel " << ichan << " h1->GetSumOfWeights() = " << h1->GetSumOfWeights() << endl;
	    if (h1->GetSumOfWeights() < 100) {
	      cout << "XXXXXXXXXXX not enough events for combined ps for channel = " << ichan << ", run hist = " << h1->GetName() << endl;
	      continue;
	    }
	    h1->Fit(f1, "lrq", "", xmin, xmax);
	    savePad(Form("yield-%s-%d-chan%d-allps.pdf", trg.c_str(), vruns[irun], ichan));
	    double A0(f1->GetParameter(0));
	    double peak(f1->GetParameter(1));
	    double peakE(f1->GetParError(1));
	    double sigma(f1->GetParameter(2));
	    double sigmaE(f1->GetParError(2));
	    double A1(f1->GetParameter(3));
	    double expslope(f1->GetParameter(4));
	    double expslopeE(f1->GetParError(4));
	    double err0(f1->GetParameter(5));
	    double err0E(f1->GetParError(5));
	    double err1(f1->GetParameter(6));
	    double err1E(f1->GetParError(6));
	    double err2(f1->GetParameter(7));
	    double err2E(f1->GetParError(7));
	    double A2(f1->GetParameter(8));
	    double nAll(h1->GetMaximum());
	    delete f1;
	    // -- now fit all histograms for the different prescales
	    double norm(0.), normE(0.);
	    cout << "prescales.size() = " << prescales.size() << endl;
	    for (unsigned ips = 0; ips < prescales.size(); ++ips) {
	      h1 = allHists[Form("chan%d_ps%d", ichan, prescales[ips])];
	      if (h1->GetSumOfWeights() < 100) {
		cout << "XXXXXXXXXXX not enough events in h1 for chan = " << ichan << ", run hist = " << h1->GetName() << endl;
		continue;
	      }
	      f1 = fIF->expoErrGauss(h1, peak, sigma, stepBp);
	      h1->SetMinimum(0.);
	      for (int ipar = 0; ipar < f1->GetNpar(); ++ipar) {
		f1->ReleaseParameter(ipar);
	      }
	      double scale = h1->GetMaximum()/nAll;
	      f1->SetParameter(0, A0*scale);
	      f1->SetParameter(1, peak);	      f1->SetParLimits(1, peak - peakE,         peak + peakE);
	      f1->SetParameter(2, sigma);	      f1->SetParLimits(2, sigma - sigmaE,       sigma + sigmaE);
	      f1->SetParameter(3, A1*scale);
	      f1->SetParameter(4, expslope);	      f1->SetParLimits(4, expslope - expslopeE, expslope + expslopeE);
	      f1->SetParameter(5, err0);	      f1->SetParLimits(5, err0 - err0E,         err0 + err0E);
	      f1->SetParameter(6, err1);	      f1->SetParLimits(6, err1 - err1E,         err1 + err1E);
	      f1->SetParameter(7, err2);	      f1->SetParLimits(7, err2 - err2E,         err2 + err2E);
	      f1->SetParameter(8, A2*scale); 	      f1->SetParLimits(8, 0., 1.e5);
	      cout << "========> Fitting ps: " << prescales[ips] << " for channel: " << ichan
		   << " histogram: " << h1->GetName() << "/" << h1->GetTitle()
		   << " histogram maximum: " << h1->GetMaximum()
		   << " sumofw8() = " << h1->GetSumOfWeights()
		   << endl;
	      //	      cout << "f1: " << endl;
	      //	      fIF->dumpParameters(f1);
	      h1->Fit(f1, "lrq", "", xmin, xmax);
	      double par0E = f1->GetParError(0);
	      double intError = f1->IntegralError(peak-3*sigma, peak+3*sigma);
	      double psNormE = intError;
	      TF1 *f2 = fIF->expoErr(xmin, xmax);
	      for (int ipar = 0; ipar < f2->GetNpar(); ++ipar) f2->SetParameter(ipar, f1->GetParameter(ipar+3));
	      f1->SetParameter(3, 0.);
	      f1->SetParameter(4, 0.);
	      f1->SetParameter(5, 0.);
	      f1->SetParameter(6, 0.);
	      f1->SetParameter(7, 0.);
	      f1->SetParameter(8, 0.);
	      f1->SetLineColor(kBlue);
	      f1->SetLineStyle(kDashed);
	      f1->Draw("same");
	      f2->SetLineColor(kRed);
	      f2->SetLineStyle(kDashed);
	      f2->Draw("same");
	      //	      cout << "f2: " << endl;
	      //	      fIF->dumpParameters(f2);
	      double psNorm = f1->Integral(peak-3*sigma, peak+3*sigma);
	      psNorm  /= h1->GetBinWidth(1);
	      psNormE /= h1->GetBinWidth(1);
	      if (psNormE < 0.01*psNorm) psNormE = par0E/h1->GetBinWidth(1);
	      if (psNormE > psNorm) psNormE = TMath::Sqrt(psNorm);
	      norm += prescales[ips]*psNorm;
	      normE += (prescales[ips]*psNormE)*(prescales[ips]*psNormE);
	      cout << "prescale: " << prescales[ips] << " Nsig = " << psNorm << " (area = " << f1->GetParameter(0)
		   << ") -> running sum: " << norm << " running error: " << TMath::Sqrt(normE) << endl;
	      if (TMath::Sqrt(normE) < 0.001) {
		cout << "XXXXXXXXX psNormE            = " << psNormE << endl;
		cout << "XXXXXXXXX f1->GetParError(0) = " << par0E << endl;
		cout << "XXXXXXXXX intError           = " << intError << endl;
		cout << "XXXXXXXXX sqrt(psNorm)       = " << TMath::Sqrt(psNorm)*h1->GetBinWidth(1) << endl;
	      }
	      tl->DrawLatex(0.2, 0.92, Form("yield-%s-%d-chan%d-ps%d.pdf", trg.c_str(), vruns[irun], ichan, prescales[ips]));
	      savePad(Form("yield-%s-%d-chan%d-ps%d.pdf", trg.c_str(), vruns[irun], ichan, prescales[ips]));
	      delete f1;
	      delete f2;
	    }

	    double result  = norm/intLumi;
	    double resultE = TMath::Sqrt(normE)/intLumi;
	    // -- normalize yield to 1/pb
	    cout << "==> Filling for chan = " << ichan << " into bin " << static_cast<double>(vruns[irun])
		 << " result = " << result << " +/- " << resultE
		 << " for intLumi = " << intLumi
		 << endl;
	    vRunHLT[ichan]->SetBinContent(vRunHLT[ichan]->FindBin(static_cast<double>(vruns[irun])), result);
	    vRunHLT[ichan]->SetBinError(vRunHLT[ichan]->FindBin(static_cast<double>(vruns[irun])), resultE);
	    vRunHLT[ichan]->Draw();
	  }
	  // -- reset lumi
	  intLumi = 0.;
	  for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
	    //allHists[Form("chan%d", ichan)]->Reset();
	    for (int ips = 1; ips < MAXPS; ++ips) {
	      allHists[Form("chan%d_ps%d", ichan, ips)]->Reset();
	    }
	  }
	  prescales.clear();
	}
      }

      for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
	vRunHLT[ichan]->SetMinimum(0.);
	vRunHLT[ichan]->Draw("e");
	savePad(Form("yield-vs-runs-%d-chan%d-%s-%s.pdf", fYear, ichan, fSample.c_str(), trg.c_str()));
      }
    }
    return;
  } else {
    TDirectory *hDir = fHistFile->mkdir(dir.c_str());
    fHistFile->cd(dir.c_str());
    hDir = gDirectory;
    cout << "created " << hDir->GetName() << endl;

    TTree *t = getTree(fSample, dir);
    if (0 == t) {
      cout << "tree for sample = " << fSample << " not found" << endl;
      return;
    }
    setupTree(t, fSample);
    fCds = fSample;
    loopOverTree(t, 1);

    cout << "writing output histograms: " << fYieldHLT.size() << endl;
    for (map<string, TH2D*>::iterator it = fYieldHLT.begin(); it != fYieldHLT.end(); ++it) {
      cout << "run " << it->first << endl;
      it->second->Draw("colz");
      it->second->SetDirectory(hDir);
      it->second->Write();
    }
    for (map<string, TH2D*>::iterator it = fYieldRTR.begin(); it != fYieldRTR.end(); ++it) {
      it->second->Draw("colz");
      it->second->SetDirectory(hDir);
      it->second->Write();
    }

    fYieldHLT.clear();
    fYieldRTR.clear();
  }

  //  fHistFile->Write();
  fHistFile->Close();
}



// ----------------------------------------------------------------------
void plotStuff::loopFunction1() {


  //  if (!fGoodMuonsID) return;

  if (!fGoodQ) return;
  if (!fGoodPvAveW8) return;
  if (!fGoodMaxDoca) return;

  if (!fb.json) return;

  if (fb.flsxy    < fCuts[fChan]->flsxy) return;
  if (fb.fls3d    < fCuts[fChan]->fls3d) return;

  if (fb.chi2dof  > fCuts[fChan]->chi2dof) return;
  if (fb.alpha    > fCuts[fChan]->alpha) return;
  if (fb.pvip     > fCuts[fChan]->pvip) return;
  if (fb.pvips    > fCuts[fChan]->pvips) return;

  if (fb.iso      < fCuts[fChan]->iso) return;
  if (fb.docatrk  < fCuts[fChan]->docatrk) return;
  if (fb.closetrk > fCuts[fChan]->closetrk) return;

  if (fb.m1pt < 4.0) return;
  if (fb.m2pt < 4.0) return;

  double m = fb.m;
  if ((fMode == BU2JPSIKP) || (fMode == BD2JPSIKSTAR) || (fMode == BS2JPSIPHI)) {
    if (TMath::Abs(fb.mpsi) < 2.9) return;
    if (TMath::Abs(fb.mpsi) > 3.3) return;

    if (TMath::Abs(fb.psipt) < 6.9) return;
    if (TMath::Abs(fb.psicosa) < 0.9) return;
    if (TMath::Abs(fb.psiprob) < 0.1) return;
    if (TMath::Abs(fb.psiflsxy) < 3) return;
    if (fMode == BS2JPSIPHI) {
      if (fb.mkk < 1.01) return;
      if (fb.mkk > 1.03) return;
      if (fb.dr  > 0.80) return;
      if (fb.kpt < 0.70) return;
    }

    if ((fMode == BU2JPSIKP) || (fMode == BD2JPSIKSTAR)) {
      if (fb.kpt < 0.70) return;
    }

    m = fb.cm;
  }


  if (0 == fYieldHLT.count(Form("%d_%d", static_cast<int>(fb.run), fb.ps))) {
    TH2D *h = new TH2D(Form("h_HLT_%d_%d", static_cast<int>(fb.run), fb.ps), Form("%d", fb.ps), 90, 5.0, 5.9, fNchan, 0., fNchan);
    fYieldHLT.insert(make_pair(Form("%d_%d", fb.run, fb.ps), h));

    h = new TH2D(Form("h_RTR_%d_%d", static_cast<int>(fb.run), fb.ps), Form("%d", fb.ps), 90, 5.0, 5.9, fNchan, 0., fNchan);
    fYieldRTR.insert(make_pair(Form("%d_%d", static_cast<int>(fb.run), fb.ps), h));
  }


  if (fb.hlt) {
    fYieldHLT[Form("%d_%d", static_cast<int>(fb.run), fb.ps)]->Fill(m, fb.chan);
  }

  if (fb.reftrg) {
    fYieldRTR[Form("%d_%d", static_cast<int>(fb.run), fb.ps)]->Fill(m, fb.chan);
  }

}


// ----------------------------------------------------------------------
void plotStuff::loopFunction2() { }
void plotStuff::loopFunction3() { }
void plotStuff::loopFunction4() { }

// ----------------------------------------------------------------------
void plotStuff::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotStuff::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotStuff::*pF)(void);
  if (ifunc == 1) pF = &plotStuff::loopFunction1;
  if (ifunc == 2) pF = &plotStuff::loopFunction2;
  if (ifunc == 3) pF = &plotStuff::loopFunction3;
  if (ifunc == 4) pF = &plotStuff::loopFunction4;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotStuff::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotStuff::loadFile loading files listed in " << files << endl;

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

      if (string::npos != stype.find("wrongreco,")) {
        sname = "wrongReco";
        sdecay = "wrongReco";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

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
  cout << "------------------------------------------------------------------------------------------" << endl;
  cout << Form("   %30s: %20s: ", "Dataset name", "Decay mode name") << "Filename:" << endl;
  cout << "------------------------------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    // cout << it->first << endl;
    // cout << it->second->fName << endl;
    // cout << it->second->fF->GetName() << endl;
    cout << Form("%2d %30s: %20s: ", cnt, it->first.c_str(), it->second->fName.c_str()) << it->second->fF->GetName() << endl;
    ++cnt;
  }
  cout << "------------------------------------------------------------------------------------------" << endl;
}
