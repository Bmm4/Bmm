#include "umlLifetime.hh"

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

#include "RooPlot.h"
#include "RooDataSet.h"

#include "RooAbsPdf.h"
#include "RooRandom.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooKeysPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "RooMCStudy.h"
#include "RooMsgService.h"
#include "RooConstVar.h"
#include "RooMinuit.h"
#include "RooProfileLL.h"



#include "common/dataset.hh"
#include "common/util.hh"
#include "common/Lumi.hh"

ClassImp(umlLifetime)

using namespace std;
using namespace RooFit;
//using namespace RooStats;

// ----------------------------------------------------------------------
umlLifetime::umlLifetime(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  umlLifetime::loadFiles(files);

  changeSetup(dir, "umlLifetime", setup);
  init();

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

  ft = new RooRealVar("t", "t", TLO, THI, "ps");
  fm = new RooRealVar("m", "m", MLO, MHI, "GeV");

}


// ----------------------------------------------------------------------
umlLifetime::~umlLifetime() {

}


// ----------------------------------------------------------------------
void umlLifetime::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void umlLifetime::makeAll(string what) {

  if (what == "model1") {
    runToys("m1", 100, 1000, 900);
  }

}



// ----------------------------------------------------------------------
model* umlLifetime::createModel1(string name, int mode) {
  model1 *aModel = new model1(name);

  // -- parameters
  aModel->bsMassPeak  = new RooRealVar("m1_bsMassPeak", "Bs mass peak", 5.369, 5.359, 5.379);
  aModel->bsMassSigma = new RooRealVar("m1_bsMassSigma", "Bs mass width", 0.04, 0.035, 0.045);
  aModel->bdMassPeak  = new RooRealVar("m1_bdMassPeak", "Bd mass peak", 5.279, 5.269, 5.289);
  aModel->bdMassSigma = new RooRealVar("m1_bdMassSigma", "Bd mass width", 0.04, 0.035, 0.045);
  aModel->bgMassSlope = new RooRealVar("m1_bgMassSlope", "bg mass slope", -0.3, -10., 10.);

  // -- fit (fixed) parameters:
  aModel->bsTau  = new RooRealVar("m1_bsTau", "B signal lifetime", TAU0, 0., 10.);
  aModel->bdTau  = new RooRealVar("m1_bdTau", "B signal lifetime", 1.52, 0., 10.);
  aModel->bgTau  = new RooRealVar("m1_bgTau", "Background lifetime", 1.2, 0., 10.);

  aModel->bsN    = new RooRealVar("m1_bsN", "Bs signal yield", 1., 0., 1.e7);
  aModel->bdN    = new RooRealVar("m1_bdN", "Bd signal yield", 1., 0., 1.e7);
  aModel->bgN    = new RooRealVar("m1_bgN", "Background yield", 1., 0., 1.e7);


  // -- create PDFs
  aModel->tTruth = new RooTruthModel("m1_tTruth", "truth model", *ft); // Build a truth resolution model (delta function)
  aModel->bsPdfT = new RooDecay("m1_bsPdfT", "Bs t", *ft, *aModel->bsTau, *aModel->tTruth, RooDecay::SingleSided);
  aModel->bdPdfT = new RooDecay("m1_bdPdfT", "Bd t", *ft, *aModel->bdTau, *aModel->tTruth, RooDecay::SingleSided);
  aModel->bgPdfT = new RooDecay("m1_bgPdfT", "background t", *ft, *aModel->bgTau, *aModel->tTruth, RooDecay::SingleSided);

  aModel->bsPdfM = new RooGaussian("m1_bsPdfM", "Bs signal mass", *fm, *aModel->bsMassPeak, *aModel->bsMassSigma);
  aModel->bdPdfM = new RooGaussian("m1_bdPdfM", "Bd signal mass", *fm, *aModel->bdMassPeak, *aModel->bdMassSigma);
  aModel->bgPdfM = new RooExponential("m1_bgPdfM", "background mass", *fm, *aModel->bgMassSlope);


  aModel->bsPdf = new RooProdPdf("m1_bsPdf", "Bs pdf",         RooArgSet(*aModel->bsPdfM, *aModel->bsPdfT));
  aModel->bdPdf = new RooProdPdf("m1_bdPdf", "Bd pdf",         RooArgSet(*aModel->bdPdfM, *aModel->bdPdfT));
  aModel->bgPdf = new RooProdPdf("m1_bgPdf", "background pdf", RooArgSet(*aModel->bgPdfM, *aModel->bgPdfT));


  aModel->modelPdf = new RooAddPdf("m1_model", "model 1",
                                   RooArgList(*aModel->bsPdf, *aModel->bdPdf, *aModel->bgPdf),
                                   RooArgList(*aModel->bsN,   *aModel->bdN,   *aModel->bgN));


  return aModel;
}


// ----------------------------------------------------------------------
model2* umlLifetime::createModel2(string name, int mode) {
  model2 *m2 = new model2(name);

  const int nchan(2);
  double mres[nchan] = {0.040, 0.040};
  model *aModel(0);
  for (int ichan = 0; ichan < nchan; ++ichan) {

    aModel = new model(name);

    // -- parameters
    aModel->bsMassPeak  = new RooRealVar(Form("m2_chan%d_bsMassPeak", ichan), "Bs mass peak", 5.369, 5.359, 5.379);
    aModel->bsMassSigma = new RooRealVar(Form("m2_chan%d_bsMassSigma", ichan), "Bs mass width", mres[ichan], mres[ichan] - 0.010, mres[ichan] + 0.010);
    aModel->bdMassPeak  = new RooRealVar(Form("m2_chan%d_bdMassPeak", ichan), "Bd mass peak", 5.279, 5.269, 5.289);
    aModel->bdMassSigma = new RooRealVar(Form("m2_chan%d_bdMassSigma", ichan), "Bd mass width", mres[ichan], mres[ichan] - 0.010, mres[ichan] + 0.010);
    aModel->bgMassSlope = new RooRealVar(Form("m2_chan%d_bgMassSlope", ichan), "bg mass slope", -0.3, -10., 10.);

    // -- fit (fixed) parameters:
    aModel->bsTau  = new RooRealVar(Form("m2_chan%d_bsTau", ichan), "B signal lifetime", TAU0, 0., 10.);
    aModel->bdTau  = new RooRealVar(Form("m2_chan%d_bdTau", ichan), "B signal lifetime", 1.52, 0., 10.);
    aModel->bgTau  = new RooRealVar(Form("m2_chan%d_bgTau", ichan), "Background lifetime", 1.2, 0., 10.);

    aModel->bsN    = new RooRealVar(Form("m2_chan%d_bsN", ichan), "Bs signal yield", 1., 0., 1.e7);
    aModel->bdN    = new RooRealVar(Form("m2_chan%d_bdN", ichan), "Bd signal yield", 1., 0., 1.e7);
    aModel->bgN    = new RooRealVar(Form("m2_chan%d_bgN", ichan), "Background yield", 1., 0., 1.e7);


    // -- create PDFs
    aModel->tTruth = new RooTruthModel(Form("m2_chan%d_tTruth", ichan), "truth model", *ft); // Build a truth resolution model (delta function)
    aModel->bsPdfT = new RooDecay(Form("m2_chan%d_bsPdfT", ichan), "Bs t", *ft, *aModel->bsTau, *aModel->tTruth, RooDecay::SingleSided);
    aModel->bdPdfT = new RooDecay(Form("m2_chan%d_bdPdfT", ichan), "Bd t", *ft, *aModel->bdTau, *aModel->tTruth, RooDecay::SingleSided);
    aModel->bgPdfT = new RooDecay(Form("m2_chan%d_bgPdfT", ichan), "background t", *ft, *aModel->bgTau, *aModel->tTruth, RooDecay::SingleSided);

    aModel->bsPdfM = new RooGaussian(Form("m2_chan%d_bsPdfM", ichan), "Bs signal mass", *fm, *aModel->bsMassPeak, *aModel->bsMassSigma);
    aModel->bdPdfM = new RooGaussian(Form("m2_chan%d_bdPdfM", ichan), "Bd signal mass", *fm, *aModel->bdMassPeak, *aModel->bdMassSigma);
    aModel->bgPdfM = new RooExponential(Form("m2_chan%d_bgPdfM", ichan), "background mass", *fm, *aModel->bgMassSlope);


    aModel->bsPdf = new RooProdPdf(Form("m2_chan%d_bsPdf", ichan), "Bs pdf",         RooArgSet(*aModel->bsPdfM, *aModel->bsPdfT));
    aModel->bdPdf = new RooProdPdf(Form("m2_chan%d_bdPdf", ichan), "Bd pdf",         RooArgSet(*aModel->bdPdfM, *aModel->bdPdfT));
    aModel->bgPdf = new RooProdPdf(Form("m2_chan%d_bgPdf", ichan), "background pdf", RooArgSet(*aModel->bgPdfM, *aModel->bgPdfT));


    aModel->modelPdf = new RooAddPdf(Form("m2_chan%d_model", ichan), Form("m2_chan%d_model", ichan),
				     RooArgList(*aModel->bsPdf, *aModel->bdPdf, *aModel->bgPdf),
				     RooArgList(*aModel->bsN,   *aModel->bdN,   *aModel->bgN));

    m2->addChannel(aModel);
  }
  m2->createSimPdf();
  m2->fSimPdf->Print("t");

  return m2;
}


// ----------------------------------------------------------------------
void umlLifetime::runToys(string whichtoy, int ntoys, int nsg, int nbg) {

  bool doPlot(false); // setting to true will create a memory leak!

  double nbd = 0.1*nsg;

  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0");
  if (!c0) c0 = new TCanvas("c0","--c0--",0, 0, 656, 400);
  // -- summary plots
  c0->Clear();
  c0->Divide(2, 2);

  TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
  if (!c1) c1 = new TCanvas("c1","--c1--",0, 0, 800, 400);
  // -- example plots
  c1->Clear();
  c1->Divide(2, 1);



  TH1D *ht = new TH1D("ht", "", 100, TAU0-0.5, TAU0+0.5);
  TH1D *hs = new TH1D("hs", "", 100, 0.,      1.0);

  TH1D *hBs = new TH1D("hBs", "", 100, nsg-0.5*nsg, nsg+0.5*nsg);
  TH1D *hBd = new TH1D("hBd", "", 100, 0., 20.);

  cout << "======================================================================" << endl;

  model *pM(0);
  model2 *pM2(0);
  RooDataSet *d0(0);
  for (int i = 0; i < ntoys; ++i) {
    if ("m1" == whichtoy) {
      pM = createModel1("m1", 0);
      RooDataSet *bgData  = pM->bgPdf->generate(RooArgSet(*fm, *ft), nbg);
      RooDataSet *bsData  = pM->bsPdf->generate(RooArgSet(*fm, *ft), nsg);
      RooDataSet *bdData  = pM->bdPdf->generate(RooArgSet(*fm, *ft), nbd);

      d0 = new RooDataSet(*bgData);
      d0->append(*bsData);
      d0->append(*bdData);
      delete bgData;
      delete bsData;
      delete bdData;
      cout << " creating new toy run " << i << " for model " << whichtoy << ", sgTau = " << pM->bsTau->getVal() << endl;
    } else if ("m2" == whichtoy) {
      pM = createModel2("m2");
      pM2 = (model2*)pM;
      double nchan(pM2->fChan.size());
      RooDataSet *d2[2];
      for (unsigned int ichan = 0; ichan < nchan; ++ichan) {
	cout << "pM2->bgPdf = " << pM2->bgPdf << " nchan = " << nchan << endl;
	RooDataSet *bgData  = pM2->fChan[ichan]->bgPdf->generate(RooArgSet(*fm, *ft), nbg/nchan);
	RooDataSet *bsData  = pM2->fChan[ichan]->bsPdf->generate(RooArgSet(*fm, *ft), nsg/nchan);
	RooDataSet *bdData  = pM2->fChan[ichan]->bdPdf->generate(RooArgSet(*fm, *ft), nbd/nchan);
	d2[ichan] = new RooDataSet(*bgData);
	d2[ichan]->append(*bsData);
	d2[ichan]->append(*bdData);
	delete bgData;
	delete bsData;
	delete bdData;
      }

      d0 = new RooDataSet("combData","combined data", RooArgSet(*fm, *ft), Index(*(pM2->fChannels)),
			  Import("chan0", *d2[0]),
			  Import("chan1", *d2[1])
			  );
      cout << " creating new toy run " << i << " for model " << whichtoy << ", sgTau = " << pM2->fChan[0]->bsTau->getVal() << endl;
    }

    // Fit pdf. The normalization integral is calculated numerically.
    RooFitResult *r = pM->modelPdf->fitTo(*d0, Save()) ;

    ht->Fill(pM->bsTau->getVal());
    hs->Fill(pM->bsTau->getError());
    hBs->Fill(pM->bsN->getVal());
    hBd->Fill(pM->bdN->getVal());

    if (doPlot || (0 == i)) {
      tl->SetNDC(kTRUE);
      tl->SetTextSize(0.04);
      c1->cd(1);
      gPad->SetLogy(0);
      RooPlot *fd0m = fm->frame(Title("d0m"), Name("mass"), Range(MLO, MHI));
      d0->plotOn(fd0m);
      pM->modelPdf->plotOn(fd0m);
      pM->modelPdf->plotOn(fd0m, Components("m1_bgPdf"), LineStyle(kDashed), LineColor(kRed)) ;
      pM->modelPdf->plotOn(fd0m, Components("m1_bsPdf"), LineStyle(kDashed), LineColor(kBlue)) ;
      pM->modelPdf->plotOn(fd0m, Components("m1_bdPdf"), LineStyle(kDashed), LineColor(kGreen)) ;
      //    pM->modelPdf->paramOn(fd0m, Layout(0.5, 0.8, 0.45));
      fd0m->Draw();

      c1->cd(2);
      gPad->SetLogy(1);
      RooPlot *fd0t = ft->frame(Title("d0t"), Name("t"), Range(TLO, THI));
      d0->plotOn(fd0t);
      pM->modelPdf->plotOn(fd0t);
      pM->modelPdf->plotOn(fd0t, Components("m1_bgPdf"), LineStyle(kDashed), LineColor(kRed)) ;
      pM->modelPdf->plotOn(fd0t, Components("m1_bsPdf"), LineStyle(kDashed), LineColor(kBlue)) ;
      pM->modelPdf->plotOn(fd0t, Components("m1_bdPdf"), LineStyle(kDashed), LineColor(kGreen)) ;
      //    pM->modelPdf->paramOn(fd0t, Layout(0.5, 0.8, 0.45));
      fd0t->Draw();
      tl->SetTextSize(0.04);
      tl->DrawLatex(0.1+0.465, 0.85, Form("N^{0} = %d", nsg));
      tl->DrawLatex(0.1+0.70,  0.85, Form("#tau_{0} = %4.3f", TAU0));
      tl->DrawLatex(0.1+0.45,  0.80, Form("N_{Bs} = %3.1f #pm %3.1f", pM->bsN->getVal(), pM->bsN->getError()));
      tl->DrawLatex(0.1+0.45,  0.75, Form("N_{Bd} = %3.1f #pm %3.1f", pM->bdN->getVal(), pM->bdN->getError()));
      tl->DrawLatex(0.1+0.49,  0.70, Form("#tau = %4.3f #pm %4.3f",
					 pM->bsTau->getVal(), pM->bsTau->getError()));


      savePad(Form("runToys-example-%s-%d.pdf", whichtoy.c_str(), nsg), c1);

      if ("m2" == whichtoy) {

	c1->Divide(2, 2);
	c1->cd(1);
	gPad->SetLogy(0);
	d0->plotOn(fd0m, Cut("channel==channel::chan0")) ;
	pM2->fSimPdf->plotOn(fd0m, Slice(*pM2->fChannels, "chan0"), ProjWData(*pM2->fChannels, *d0)) ;

	c1->cd(2);
	gPad->SetLogy(1);
	d0->plotOn(fd0t, Cut("channel==channel::chan0")) ;

	c1->cd(3);
	gPad->SetLogy(0);
	d0->plotOn(fd0m, Cut("channel==channel::chan1")) ;

	c1->cd(2);
	gPad->SetLogy(1);
	d0->plotOn(fd0t, Cut("channel==channel::chan1")) ;

      }



      c0->cd(1);
      hBs->Draw();
      tl->SetTextSize(0.05);
      tl->DrawLatex(0.25, 0.87, Form("#mu = %4.3f #pm %4.3f", hBs->GetMean(), hBs->GetMeanError()));
      tl->SetTextSize(0.035);
      tl->DrawLatex(0.70, 0.96, Form("Nentries = %d", static_cast<int>(hBs->GetEntries())));
      tl->DrawLatex(0.25, 0.96, Form("N^{0} = %d", nsg));

      c0->cd(2);
      hBd->Draw();
      tl->SetTextSize(0.05);
      tl->DrawLatex(0.25, 0.87, Form("#mu = %4.3f #pm %4.3f", hBd->GetMean(), hBd->GetMeanError()));
      tl->SetTextSize(0.035);
      tl->DrawLatex(0.25, 0.96, Form("N^{0} = %4.1f", nbd));


      c0->cd(3);
      ht->Draw();
      tl->SetTextSize(0.05);
      tl->DrawLatex(0.25, 0.87, Form("#mu = %4.3f #pm %4.3f", ht->GetMean(), ht->GetMeanError()));
      tl->SetTextSize(0.035);
      tl->DrawLatex(0.25, 0.96, Form("tau^{0} = %4.3f", TAU0));

      c0->cd(4);
      hs->Draw();
      tl->SetTextSize(0.05);
      tl->DrawLatex(0.25, 0.87, Form("#mu = %4.3f #pm %4.3f", hs->GetMean(), hs->GetMeanError()));

      c0->Modified();
      c0->Update();

      c1->Modified();
      c1->Update();

    }
    delete r;
    delete pM;
    delete d0;
  }


  c0->cd(1);
  hBs->SetMaximum(1.3*hBs->GetMaximum());
  hBs->Draw();
  tl->SetTextSize(0.05);
  tl->DrawLatex(0.25, 0.87, Form("#mu = %4.1f #pm %4.1f", hBs->GetMean(), hBs->GetMeanError()));
  tl->SetTextSize(0.035);
  tl->DrawLatex(0.70, 0.96, Form("Nentries = %d", static_cast<int>(hBs->GetEntries())));
  tl->DrawLatex(0.25, 0.96, Form("N^{0} = %d", nsg));
  pa->DrawArrow(nsg, 0.5*hBs->GetMaximum(), nsg, 0.);

  c0->cd(2);
  hBd->SetMaximum(1.3*hBd->GetMaximum());
  hBd->Draw();
  tl->SetTextSize(0.05);
  tl->DrawLatex(0.25, 0.87, Form("#mu = %4.1f #pm %4.1f", hBd->GetMean(), hBd->GetMeanError()));
  tl->SetTextSize(0.035);
  pa->DrawArrow(nbd, 0.5*hBd->GetMaximum(), nbd, 0.);

  c0->cd(3);
  ht->SetMaximum(1.3*ht->GetMaximum());
  ht->Draw();
  tl->SetTextSize(0.05);
  tl->DrawLatex(0.25, 0.87, Form("#mu = %4.3f #pm %4.3f", ht->GetMean(), ht->GetMeanError()));
  tl->SetTextSize(0.035);
  tl->DrawLatex(0.70, 0.96, Form("Nentries = %d", static_cast<int>(ht->GetEntries())));
  pa->DrawArrow(TAU0, 0.5*ht->GetMaximum(), TAU0, 0.);

  c0->cd(4);
  hs->SetMaximum(1.3*hs->GetMaximum());
  hs->Draw();
  tl->SetTextSize(0.05);
  tl->DrawLatex(0.25, 0.87, Form("#mu = %4.3f #pm %4.3f", hs->GetMean(), hs->GetMeanError()));

  savePad(Form("runToys-summary-%s-%d.pdf", whichtoy.c_str(), nsg), c0);
}


// ----------------------------------------------------------------------
void umlLifetime::bookHist(string dsname) {


}


// ----------------------------------------------------------------------
void umlLifetime::loopFunction1() {

}


// ----------------------------------------------------------------------
void umlLifetime::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> umlLifetime::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (umlLifetime::*pF)(void);
  if (ifunc == 1) pF = &umlLifetime::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void umlLifetime::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> umlLifetime::loadFile loading files listed in " << files << endl;

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

    if (string::npos != stype.find("Data")) {
      // -- Charmonium
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("blablabla")) {
        sname = "bmmCharmonium";
        sdecay = "bmm";
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
