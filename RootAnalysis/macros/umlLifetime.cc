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
    model *pm = createModel1("m1", 0);
    runToy(pm, 1000, 900);
  }

}



// ----------------------------------------------------------------------
model* umlLifetime::createModel1(string name, int mode) {
  model *aModel = new model();

  aModel->name = "m1";

  // -- variables
  aModel->m  = new RooRealVar("m1_m", "m", MLO, MHI, "GeV");
  aModel->t  = new RooRealVar("m1_t", "t", TLO, THI, "ps");

  // -- parameters
  aModel->sgMassPeak  = new RooRealVar("m1_sgMassPeak", "B mass peak", 5.279, MLO, MHI);
  aModel->sgMassSigma = new RooRealVar("m1_sgMassSigma", "B mass width", 0.04, 0., 1.);
  aModel->bgMassSlope = new RooRealVar("m1_bgMassSlope", "bg mass slope", -0.3, -10., 10.);

  // -- fit (fixed) parameters:
  aModel->sgTau  = new RooRealVar("m1_sgTau", "B signal lifetime", 1.5, 0., 10.);
  aModel->bgTau  = new RooRealVar("m1_bgTau", "Background lifetime", 1.2, 0., 10.);

  aModel->sgN    = new RooRealVar("m1_sgN", "B signal yield", 1., 0., 1.e7);
  aModel->bgN    = new RooRealVar("m1_bgN", "Background yield", 1., 0., 1.e7);


  // -- create PDFs
  aModel->tTruth = new RooTruthModel("m1_tTruth", "truth model", *aModel->t); // Build a truth resolution model (delta function)
  aModel->sgPdfT = new RooDecay("m1_sgPdfT", "signal t", *aModel->t, *aModel->sgTau, *aModel->tTruth, RooDecay::SingleSided);
  aModel->bgPdfT = new RooDecay("m1_bgPdfT", "background t", *aModel->t, *aModel->bgTau, *aModel->tTruth, RooDecay::SingleSided);

  aModel->sgPdfM = new RooGaussian("m1_sgPdfM", "B signal mass", *aModel->m, *aModel->sgMassPeak, *aModel->sgMassSigma);
  aModel->bgPdfM = new RooExponential("m1_bgPdfM", "background mass", *aModel->m, *aModel->bgMassSlope);


  aModel->sgPdf = new RooProdPdf("m1_sgPdf", "signal pdf",     RooArgSet(*aModel->sgPdfM, *aModel->sgPdfT));
  aModel->bgPdf = new RooProdPdf("m1_bgPdf", "background pdf", RooArgSet(*aModel->bgPdfM, *aModel->bgPdfT));


  aModel->modelPdf = new RooAddPdf("m1_model", "model 1",
                                   RooArgList(*aModel->sgPdf, *aModel->bgPdf),
                                   RooArgList(*aModel->sgN, *aModel->bgN));


  return aModel;
}


// ----------------------------------------------------------------------
void umlLifetime::runToy(model *pM, int nsg, int nbg) {

  RooDataSet *bgData  = pM->bgPdf->generate(RooArgSet(*pM->m, *pM->t), nbg);
  RooDataSet *sgData  = pM->sgPdf->generate(RooArgSet(*pM->m, *pM->t), nsg);

  RooDataSet *d0 = new RooDataSet(*bgData);
  d0->append(*sgData);


  // Fit pdf. The normalization integral is calculated numerically.
  pM->modelPdf->fitTo(*d0) ;


  if (1) {
    tl->SetNDC(kTRUE);
    tl->SetTextSize(0.04);

    TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0");
    if (!c0) c0 = new TCanvas("c0","--c0--",0,0, 656, 400);

    // -- data 0
    c0->Clear();
    c0->Divide(2,2);

    c0->cd(1);
    gPad->SetLogy(0);
    RooPlot *fd0m = pM->m->frame(Title("d0m"), Name("mass"), Range(MLO, MHI));
    d0->plotOn(fd0m);
    pM->modelPdf->plotOn(fd0m);
    pM->modelPdf->plotOn(fd0m, Components("m1_bgPdf"), LineStyle(kDashed)) ;
    //    pM->modelPdf->paramOn(fd0m, Layout(0.5, 0.8, 0.45));
    fd0m->Draw();

    c0->cd(2);
    gPad->SetLogy(1);
    RooPlot *fd0t = pM->t->frame(Title("d0t"), Name("t"), Range(TLO, THI));
    d0->plotOn(fd0t);
    pM->modelPdf->plotOn(fd0t);
    pM->modelPdf->plotOn(fd0t, Components("m1_bgPdf"), LineStyle(kDashed)) ;
    //    pM->modelPdf->paramOn(fd0t, Layout(0.5, 0.8, 0.45));
    fd0t->Draw();

  }



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
