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

  double minLumi(100.);

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
  cout << "OK? " << ok << endl;
  if (ok) {
    cout << "histograms exist already, looping over them" << endl;
    TIter next(gDirectory->GetListOfKeys());
    TKey *key(0);
    vector<int> vds;
    int run(-1), runMin(9999999), runMax(0);
    while ((key = (TKey*)next())) {
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
      if (TString(key->GetName()).Contains("_HLT_")) {
	string hname = key->GetName();
	replaceAll(hname, "h_HLT_", "");
	run = atoi(hname.c_str());
	if (run > runMax) runMax = run;
	if (run < runMin) runMin = run;
	vds.push_back(run);
      }
    }

    if (vds.size() > 0) {
      Lumi lumi("../common/json/Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_MuonPhys.lumi");
      cout << "runs " << runMin << " .. " <<  runMax << endl;

      vector<TH1D *> vRunHLT;
      for (unsigned int ichan = 0; ichan < fNchan; ++ichan) {
	vRunHLT.push_back(new TH1D(Form("hRun%s_chan%d", trg.c_str(), ichan), Form("hRun%s_chan%d", trg.c_str(), ichan), runMax-runMin+1, runMin, runMax));
	vRunHLT[ichan]->Sumw2();
      }

      double mBp(5.28), sBp(0.04), stepBp(5.145);
      double xmin(5.0), xmax(5.8), ymax(0.);
      fIF->fLo = xmin;
      fIF->fHi = xmax;
      TH2D *h2 = (TH2D*)(gDirectory->Get(Form("h_%s_%d", trg.c_str(), vds[0])));
      if (!h2) {
	cout << "histogram " << Form("h_%s_%d", trg.c_str(), vds[0]) << " not found!?! returning" << endl;
	return;
      }
      TH2D *h2Sum = (TH2D*)h2->Clone("h2sum"); h2Sum->Reset();
      TH1D *h1 = h2->ProjectionX("chan_0", 1,1);
      h1->SetName("base"); h1->Reset();
      double intLumi(0.);
      int oldPs(atoi(h2->GetTitle())), newPs(0);
      bool psChange(false);
      for (unsigned int i = 0; i < vds.size(); ++i) {
	h2 = (TH2D*)(gDirectory->Get(Form("h_%s_%d", trg.c_str(), vds[i])));
	newPs = atoi(h2->GetTitle());
	cout << "_____ get " << h2->GetName()
	     << " with lumi = " << lumi.lumi(vds[i])
	     << " and prescale = " << newPs
	     << endl;
	if (oldPs == newPs)  {
	  intLumi += lumi.lumi(vds[i]);
	  h2Sum->Add(h2);
	} else {
	  cout << "XXXXXXXX new PS = " << newPs << " while old PS = " << oldPs << ", intLumi = " << intLumi << endl;
	  psChange = true;
	}
	if (psChange || (intLumi > minLumi)) {
	  cout << "==== now fitting because intLumi = " << intLumi << " and/or psChange = " << (psChange? "TRUE": "FALSE") << endl;
	  for (unsigned ichan = 0; ichan < fNchan; ++ichan) {
	    // -- fit
	    h1 = h2Sum->ProjectionX("chan_0", ichan+1, ichan+1);
	    if (h1->GetSumOfWeights() < 100) {
	      cout << "XXXXXXXXXXX not enough events in h2Sum for chan = " << ichan << ", run hist = " << h2->GetName() << endl;
	      continue;
	    }
	    TF1 *f1 = fIF->expoErrGauss(h1, mBp, sBp, stepBp);
	    h1->SetMinimum(0.);
	    h1->Fit(f1, "lr", "", xmin, xmax);
	    double nNormE = f1->IntegralError(5.1, 5.4)/f1->Integral(5.1, 5.4);
	    if (nNormE < 0.1) nNormE = f1->GetParError(0)/f1->GetParameter(0);
	    f1->SetParameter(3, 0.);
	    f1->SetParameter(4, 0.);
	    f1->SetParameter(5, 0.);
	    f1->SetParameter(6, 0.);
	    double nNorm = f1->Integral(5.1, 5.4)/h1->GetBinWidth(1);
	    if (nNormE > nNorm) nNormE = TMath::Sqrt(nNorm)/nNorm;
	    double result  = oldPs*nNorm/intLumi;
	    double resultE = oldPs*nNormE/intLumi;
	    // -- normalize yield to 1/pb
	    cout << "==> Filling for chan = " << ichan << " into bin " << static_cast<double>(vds[i])
		 << " result = " << result << " +/- " << nNormE*nNorm/intLumi
		 << " for intLumi = " << intLumi
		 << " prescale = " << oldPs
		 << endl;
	    vRunHLT[ichan]->SetBinContent(vRunHLT[ichan]->FindBin(static_cast<double>(vds[i])), result);
	    vRunHLT[ichan]->SetBinError(vRunHLT[ichan]->FindBin(static_cast<double>(vds[i])), resultE);
	    savePad(Form("yield-%s-%d-chan%d.pdf", trg.c_str(), vds[i], ichan));
	  }
	  // -- reset
	  intLumi = 0.;
	  h2Sum->Reset();
	  if (psChange) {
	    intLumi += lumi.lumi(vds[i]);
	    h2Sum->Add(h2);
	    oldPs = newPs;
	    psChange = false;
	  }
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
