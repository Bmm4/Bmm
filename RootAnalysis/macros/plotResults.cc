#include "plotResults.hh"

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

ClassImp(plotResults)

using namespace std;

// ----------------------------------------------------------------------
plotResults::plotResults(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotResults::loadFiles(files);

  changeSetup(dir, "plotResults", setup);
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
plotResults::~plotResults() {

}


// ----------------------------------------------------------------------
void plotResults::init() {
  // cout << Form("/bin/rm -f %s", fHistFileName.c_str()) << endl;
  // system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotResults::makeAll(string what) {
  if (what == "all" || what == "dumpdatasets") {
    dumpDatasets();
  }

  if (what == "all" || what == "genvalidation") {
    genSummary("bdmmMcOffAcc", "candAnaMuMu");
    genSummary("bdmmMcOff", "candAnaMuMu");
    genSummary("bdmmMc", "candAnaMuMu");
    genSummary("bsmmMcOffAcc", "candAnaMuMu");
    genSummary("bsmmMcOff", "candAnaMuMu");
    genSummary("bsmmMc", "candAnaMuMu");

    genSummary("bupsikMcOffAcc", "candAnaBu2JpsiK");
    genSummary("bupsikMcOff", "candAnaBu2JpsiK");
    genSummary("bupsikMc", "candAnaBu2JpsiK");
    genSummary("bspsiphiMcOffAcc", "candAnaBs2JpsiPhi");
    genSummary("bspsiphiMcOff", "candAnaBs2JpsiPhi");
    genSummary("bspsiphiMc", "candAnaBs2JpsiPhi");
    genSummary("bdpsikstarMcAcc", "candAnaBd2JpsiKstar");
    genSummary("bdpsikstarMc", "candAnaBd2JpsiKstar");


    // -- loop over all (effective) two-body backgrounds
    for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
      if (string::npos == it->first.find("Bg")) continue;
      cout << "===> Create genSummary for " << it->first << endl;
      genSummary(it->first, "candAnaMuMu");
    }

  }


}



// ----------------------------------------------------------------------
void plotResults::bookHist(string dsname) {


}


// ----------------------------------------------------------------------
void plotResults::dumpDatasets() {
  fTEX << "% ----------------------------------------------------------------------" << endl;
  fTEX << formatTex(fBfPsiMuMu, Form("%s:BfPsiMuMu:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(fBfPsiMuMuE, Form("%s:BfPsiMuMu:err", fSuffix.c_str()), 5) << endl;

  fTEX << formatTex(fBfPhiKpKm, Form("%s:BfPhiKpKm:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(fBfPhiKpKmE, Form("%s:BfPhiKpKm:err", fSuffix.c_str()), 5) << endl;

  fTEX << formatTex(fBfKstarKpPim, Form("%s:BfKstarKpPim:val", fSuffix.c_str()), 5) << endl;
  fTEX << formatTex(fBfKstarKpPim, Form("%s:BfKstarKpPim:err", fSuffix.c_str()), 5) << endl;

  fTEX << formatTexErrSci(fCrossSection, 0., Form("%s:PythiaCrossSection:val", fSuffix.c_str()), -1) << endl;

  fTEX << "% ----------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    TH1D *h1 = it->second->getHist("monEvents", false);
    int nEvtFile = static_cast<int>(h1->GetBinContent(1));
    int nCands   = static_cast<int>(h1->GetBinContent(2));
    double epscand  = static_cast<double>(nCands)/nEvtFile;
    double epscandE = dEff(nCands, nEvtFile);
    fTEX << Form("\\vdef{%s:%s:name} {%s}", fSuffix.c_str(), it->first.c_str(), it->first.c_str()) << endl;
    fTEX << Form("\\vdef{%s:%s:decay} {%s}", fSuffix.c_str(), it->first.c_str(), it->second->fLatexName.c_str()) << endl;
    fTEX << formatTex(nEvtFile, Form("%s:%s:nEvtFile", fSuffix.c_str(), it->first.c_str()), 0) << endl;
    fTEX << formatTex(nCands, Form("%s:%s:nCands", fSuffix.c_str(), it->first.c_str()), 0) << endl;
    fTEX << formatTex(epscand, Form("%s:%s:epsCand:val", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    fTEX << formatTex(epscandE, Form("%s:%s:epsCand:err", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    fTEX << formatTex(it->second->fFilterEff, Form("%s:%s:filterEff:val", fSuffix.c_str(), it->first.c_str()), 6) << endl;
    fTEX << formatTex(0., Form("%s:%s:filterEff:err", fSuffix.c_str(), it->first.c_str()), 4) << endl;
    if (it->second->fBf > 0.) {
      fTEX << formatTexErrSci(it->second->fBf, it->second->fBfE, Form("%s:%s:bf", fSuffix.c_str(), it->first.c_str()), 2) << endl;
      double eqLumi(0.);
      if (it->second->fFilterEff > 0) {
	eqLumi = nEvtFile/fCrossSection/it->second->fBf/it->second->fFilterEff;
	fTEX << formatTex(eqLumi, Form("%s:%s:eqLumi:val", fSuffix.c_str(), it->first.c_str()), 1) << endl;
	cout << it->first << ": eqLumi = " << eqLumi << " filterEff: " << it->second->fFilterEff << endl;
      }
    }
    if (it->second->fLumi > 0.) {
      fTEX << formatTex(it->second->fLumi, Form("%s:%s:lumi", fSuffix.c_str(), it->first.c_str()), 1) << endl;
    }
  }

}


// ----------------------------------------------------------------------
void plotResults::genSummary(std::string dsname, std::string dir) {
  TH1D *hpt    = new TH1D("pt", "pt", 50, 0., 50.0);
  TH1D *heta   = new TH1D("eta", "eta", 40, -4., 4.0);
  TH1D *tpt    = new TH1D("tpt", "pt (HLT)", 50, 0., 50.0); setFilledHist(tpt, kBlue, kYellow, 1000);
  TH1D *teta   = new TH1D("teta", "eta (HLT)", 40, -4., 4.0); setFilledHist(teta, kBlue, kYellow, 1000);
  TH1D *hm1eta = new TH1D("m1eta", "m1 eta", 40, -4., 4.0);
  TH1D *hm2eta = new TH1D("m2eta", "m2 eta", 40, -4., 4.0);
  TH1D *tm1eta = new TH1D("tm1eta", "m1 eta (HLT)", 40, -4., 4.0); setFilledHist(tm1eta, kBlue, kYellow, 1000);
  TH1D *tm2eta = new TH1D("tm2eta", "m2 eta (HLT)", 40, -4., 4.0); setFilledHist(tm2eta, kBlue, kYellow, 1000);
  TH1D *hketa  = new TH1D("keta", "kaon eta", 40, -4., 4.0);
  TH1D *hm1pt  = new TH1D("m1pt", "m1 pt", 100, 0., 10.0);
  TH1D *hm2pt  = new TH1D("m2pt", "m2 pt", 100, 0., 10.0);
  TH1D *tm1pt  = new TH1D("tm1pt", "m1 pt (HLT)", 100, 0., 10.0); setFilledHist(tm1pt, kBlue, kYellow, 1000);
  TH1D *tm2pt  = new TH1D("tm2pt", "m2 pt (HLT)", 100, 0., 10.0); setFilledHist(tm2pt, kBlue, kYellow, 1000);
  TH1D *hkpt   = new TH1D("kpt", "kaon pt", 50, 0., 10.0);
  TH1D *htau   = new TH1D("tau", "tau", 100, 0., 15.e-12);

  TTree *T = getTree(dsname, dir, "effTree");
  T->Draw("gtau>>tau");
  T->Draw("gpt>>pt");
  T->Draw("geta>>eta");

  T->Draw("g1pt>>m1pt");
  T->Draw("g2pt>>m2pt");
  T->Draw("g1eta>>m1eta");
  T->Draw("g2eta>>m2eta");

  T->Draw("gpt>>tpt", "hlt");
  T->Draw("g1pt>>tm1pt", "hlt");
  T->Draw("g2pt>>tm2pt", "hlt");

  T->Draw("geta>>teta", "hlt");
  T->Draw("g1eta>>tm1eta", "hlt");
  T->Draw("g2eta>>tm2eta", "hlt");


  bool addKaon(false);
  bool addHLT(true);
  if (string::npos != dsname.find("Bg")) {
    addHLT = false;
  }

  if (string::npos != dsname.find("bupsik")) {
    T->Draw("g3eta>>keta");
    T->Draw("g3pt>>kpt");
    addKaon = true;
  }
  if (string::npos != dsname.find("bspsiphi")) {
    T->Draw("g3eta>>keta");
    T->Draw("g4eta>>keta");
    T->Draw("g3pt>>kpt");
    T->Draw("g4pt>>kpt");
    addKaon = true;
  }

  tl->SetTextSize(0.05);
  makeCanvas(1);
  int ncol(4);
  if (addKaon) {
    c1->Divide(5,2);
    ncol = 5;
  } else {
    c1->Divide(4,2);
  }

  c1->cd(1);
  setTitles(hpt, "p_{T} [GeV]", "entries/bin");
  hpt->Draw();
  if (addHLT) {
    tpt->Draw("same");
    hpt->Draw("sameaxis");
  }

  c1->cd(2);
  hm1pt->Draw();
  setTitles(hm1pt, "p_{T} [GeV]", "entries/bin");
  if (addHLT) {
    tm1pt->Draw("same");
    hm1pt->Draw("sameaxis");
  }

  c1->cd(3);
  hm2pt->Draw();
  setTitles(hm2pt, "p_{T} [GeV]", "entries/bin");
  if (addHLT) {
    tm2pt->Draw("same");
    hm2pt->Draw("sameaxis");
  }

  if (addKaon) {
    c1->cd(ncol-1);
    setTitles(hkpt, "p_{T} [GeV]", "entries/bin");
    hkpt->Draw();
  }

  c1->cd(ncol);
  gPad->SetLogy(1);
  htau->Fit("expo", "l");
  TF1 *f = (TF1*)htau->GetFunction("expo");
  double chi2 = f->GetChisquare();
  int    ndf  = f->GetNDF();
  double t    = -1./f->GetParameter(1);
  double tE   = -t*f->GetParError(1)/f->GetParameter(1);
  t  *= 1.e12;
  tE *= 1.e12;

  c1->cd(ncol+1);
  setTitles(heta, "#eta", "entries/bin");
  heta->Draw();
  if (addHLT) {
    teta->Draw("same");
    heta->Draw("axissame");
  }
  tl->DrawLatexNDC(0.55, 0.3, "B");

  c1->cd(ncol+2);
  setTitles(hm1eta, "#eta", "entries/bin");
  hm1eta->Draw();
  if (addHLT) {
    tm1eta->Draw("same");
    hm1eta->Draw("sameaxis");
  }
  tl->DrawLatexNDC(0.40, 0.3, "leading muon");

  c1->cd(ncol+3);
  setTitles(hm2eta, "#eta", "entries/bin");
  hm2eta->Draw();
  if (addHLT) {
    tm2eta->Draw("same");
    hm2eta->Draw("sameaxis");
  }
  tl->DrawLatexNDC(0.35, 0.3, "subleading muon");

  if (addKaon) {
    c1->cd(2*ncol-1);
    setTitles(hketa, "#eta", "entries/bin");
    hketa->Draw();
    tl->DrawLatexNDC(0.35, 0.3, "kaon(s)");
  }

  c1->cd(2*ncol);
  tl->SetTextSize(0.07);
  if (addHLT) {
    tl->DrawLatexNDC(0.2, 0.9, Form("%s (HLT)", dsname.c_str()));
  } else {
    tl->DrawLatexNDC(0.2, 0.9, Form("%s ", dsname.c_str()));
  }
  tl->DrawLatexNDC(0.2, 0.8, Form("Events: %d", T->GetEntries()));

  tl->DrawLatexNDC(0.2, 0.35, Form("#tau"));
  tl->DrawLatexNDC(0.25, 0.35, Form("= (%3.2f #pm %5.2f)ps", t, tE));

  tl->DrawLatexNDC(0.2, 0.25, Form("#varepsilon"));
  tl->DrawLatexNDC(0.25, 0.25, Form("= %4.3f ", teta->GetSumOfWeights()/heta->GetSumOfWeights()));

  c1->SaveAs(Form("%s/genSummary-%s.pdf", fDirectory.c_str(), dsname.c_str()));

}




// ----------------------------------------------------------------------
void plotResults::loopFunction1() {

}

// ----------------------------------------------------------------------
void plotResults::loopFunction2() {

}


// ----------------------------------------------------------------------
void plotResults::loopFunction3() {


}



// ----------------------------------------------------------------------
void plotResults::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotResults::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotResults::*pF)(void);
  if (ifunc == 1) pF = &plotResults::loopFunction1;
  if (ifunc == 2) pF = &plotResults::loopFunction2;
  if (ifunc == 3) pF = &plotResults::loopFunction3;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotResults::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotResults::loadFile loading files listed in " << files << endl;

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
