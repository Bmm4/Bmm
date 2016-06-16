#include "plotWork.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
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

ClassImp(plotWork)

using namespace std;

// ----------------------------------------------------------------------
plotWork::plotWork(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotWork::loadFiles(files);

  changeSetup(dir, "plotWork", setup);

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

}


// ----------------------------------------------------------------------
plotWork::~plotWork() {

}



// ----------------------------------------------------------------------
void plotWork::init() {
  system(Form("/bin/rm -f %s", fHistFileName.c_str()));
  system(Form("/bin/rm -f %s", fNumbersFileName.c_str()));
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
}


// ----------------------------------------------------------------------
void plotWork::makeAll(int bitmask) {

  if (bitmask & 0x1) {
  }

}


// ----------------------------------------------------------------------
void plotWork::bookHist(string dsname) {

  fpHnorm = new TH1D(Form("h_%s_%s", "norm", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);
  fpHpass = new TH1D(Form("h_%s_%s", "pass", dsname.c_str()), Form("h_%s_%s", "all", dsname.c_str()), 40, 4.8, 6.0);

}




// ----------------------------------------------------------------------
void plotWork::tisEfficiency(string dsname) {

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
  bookHist(fSample);
  setupTree(t, fSample);
  fCds = fSample;
  loopOverTree(t, 1);


  fpHnorm->Write();
  fpHpass->Write();
  fHistFile->Close();

}



// ----------------------------------------------------------------------
void plotWork::loopFunction1() {

  if (!fb.tis) return;

  if (!fGoodQ) return;
  if (!fGoodPvAveW8) return;
  if (!fGoodMaxDoca) return;
  if (!fGoodIp) return;
  if (!fGoodIpS) return;

  if (!fGoodLip) return;
  if (!fGoodLipS) return;

  if (!fGoodAlpha) return;
  if (!fGoodChi2) return;
  if (!fGoodFLS) return;

  if (!fGoodCloseTrack) return;
  if (!fGoodIso) return;
  if (!fGoodDocaTrk) return;

  if (TMath::Abs(fb.flsxy) < 3.0) return;

  if (TMath::Abs(fb.m1eta) > 1.6) return;
  if (TMath::Abs(fb.m2eta) > 1.6) return;


  fpHnorm->Fill(fb.m);

  if (!fb.hlt) return;
  fpHpass->Fill(fb.m);

}



// ----------------------------------------------------------------------
void plotWork::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotWork::loopOverTree> loop over dataset " << fCds << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotWork::*pF)(void);
  if (ifunc == 1) pF = &plotWork::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
// NOTE: This works with the output of genAnalysis!!!!
void plotWork::prodSummary(string ds1, int year) {
  static double M511(0.), M521(0.), M531(0.), M541(0.), M5122(0.);
  static double L511(0.), L521(0.), L531(0.), L541(0.), L5122(0.);

  if (2014 == year) {
    M511 = 5.27958;
    M521 = 5.27926;
    M531 = 5.36677;
    M541 = 6.2756;
    M5122= 5.6195;

    L511 = 455.4;
    L521 = 491.1;
    L531 = 453.3;
    L541 = 135.5;
    L5122= 435;
  }

  static const double aparticles[] = {511, 521, 531, 5122};
  vector<int> particles(aparticles, aparticles + sizeof(aparticles)/sizeof(aparticles[0]));

  if (fDS[ds1]) {
    fDS[ds1]->cd("");
  } else {
    cout << "fDS[" << ds1 << "] not found" << endl;
    return;
  }
  // -- Masses
  cout << "Masses" << endl;
  cout << Form("B0: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
	       ((TH1D*)gFile->Get("m511"))->GetMean(),
	       M511,
	       ((TH1D*)gFile->Get("m511"))->GetMean() - M511
	       )
       << endl;

  cout << Form("B+: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
	       ((TH1D*)gFile->Get("m521"))->GetMean(),
	       M521,
	       ((TH1D*)gFile->Get("m521"))->GetMean() - M521
	       )
       << endl;

  cout << Form("Bs: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
	       ((TH1D*)gFile->Get("m531"))->GetMean(),
	       M531,
	       ((TH1D*)gFile->Get("m531"))->GetMean() - M531
	       )
       << endl;

  cout << Form("Lb: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
	       ((TH1D*)gFile->Get("m5122"))->GetMean(),
	       M5122,
	       ((TH1D*)gFile->Get("m5122"))->GetMean() - M5122
	       )
       << endl;

  if (((TH1D*)gFile->Get("m541"))->GetEntries() > 1)
    cout << Form("Bc: %6.5f (PDG = %6.5f, diff = %+6.5f GeV)",
		 ((TH1D*)gFile->Get("m541"))->GetMean(),
		 M541,
		 ((TH1D*)gFile->Get("m541"))->GetMean() - M541
		 )
	 << endl;



  // -- Lifetimes
  cout << "Lifetime" << endl;
  double t(1.), tE(1.), chi2(0);
  int ndf(0);
  TH1D *h = (TH1D*)gFile->Get("t521");
  TF1 *f(0);
  gStyle->SetOptFit(1);
  if (h->GetEntries() > 100) {
    h->Fit("expo", "lq");
    //      gPad->SaveAs("t521.pdf");
    f = (TF1*)h->GetFunction("expo");
    chi2 = f->GetChisquare()/f->GetNDF();
    t    = -1./f->GetParameter(1);
    tE   = -t*f->GetParError(1)/f->GetParameter(1);
    cout << Form("B+: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L521, (t-L521)/tE, chi2) << endl;
  }
  h = (TH1D*)gFile->Get("t511");
  if (h->GetEntries() > 100) {
    h->Fit("expo", "ql");
    //      gPad->SaveAs("t511.pdf");
    f = (TF1*)h->GetFunction("expo");
    chi2 = f->GetChisquare()/f->GetNDF();
    t    = -1./f->GetParameter(1);
    tE   = -t*f->GetParError(1)/f->GetParameter(1);
    cout << Form("B0: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L511, (t-L511)/tE, chi2) << endl;
  }

  h = (TH1D*)gFile->Get("t531");
  if (h->GetEntries() > 100) {
    h->Fit("expo", "ql");
    //      gPad->SaveAs("t531.pdf");
    f = (TF1*)h->GetFunction("expo");
    chi2 = f->GetChisquare()/f->GetNDF();
    t    = -1./f->GetParameter(1);
    tE   = -t*f->GetParError(1)/f->GetParameter(1);
    cout << Form("Bs: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L531, (t-L531)/tE, chi2) << endl;
  }

  h = (TH1D*)gFile->Get("t5122");
  if (h->GetEntries() > 100) {
    h->Fit("expo", "ql");
    //      gPad->SaveAs("t5122.pdf");
    f = (TF1*)h->GetFunction("expo");
    chi2 = f->GetChisquare()/f->GetNDF();
    t    = -1./f->GetParameter(1);
    tE   = -t*f->GetParError(1)/f->GetParameter(1);
    cout << Form("Lb: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L5122, (t-L5122)/tE, chi2) << endl;
  }

  h = (TH1D*)gFile->Get("t541");
  if (h->GetEntries() > 100) {
    h->Fit("expo", "ql");
    //      gPad->SaveAs("t541.pdf");
    f = (TF1*)h->GetFunction("expo");
    chi2 = f->GetChisquare()/f->GetNDF();
    t    = -1./f->GetParameter(1);
    tE   = -t*f->GetParError(1)/f->GetParameter(1);
    cout << Form("Bc: %4.2f+/-%4.2f (PDG = %4.2f, pull = %+4.2f, chi2 = %4.1f)", t, tE, L541, (t-L541)/tE, chi2) << endl;
  }

  //  -- save at the end to remove the intermittent root printout
  tl->SetTextSize(0.03);
  gPad->SetLogy(1);
  h = ((TH1D*)gFile->Get("t521"));
  f = (TF1*)h->GetFunction("expo");
  chi2 = f->GetChisquare();
  ndf  = f->GetNDF();
  t    = -1./f->GetParameter(1);
  tE   = -t*f->GetParError(1)/f->GetParameter(1);
  h->SetTitle("");
  setTitles(h, "c#tau [#mum]", "Entries/bin");
  h->Draw();
  setItalic(); tl->DrawLatexNDC(0.2, 0.45, "B^{+}"); setRoman();
  tl->DrawLatexNDC(0.2, 0.40, Form("#chi^{2}/ndf")); tl->DrawLatexNDC(0.3, 0.40, Form("= %3.1f/%d", chi2, ndf));
  tl->DrawLatexNDC(0.2, 0.35, Form("#tau(fit)"));    tl->DrawLatexNDC(0.3, 0.35, Form("= %5.2f #pm %5.2f", t, tE));
  tl->DrawLatexNDC(0.2, 0.30, Form("#tau(pdg)"));    tl->DrawLatexNDC(0.3, 0.30, Form("= %5.2f ", L521));
  tl->DrawLatexNDC(0.2, 0.25, Form("pull"));         tl->DrawLatexNDC(0.3, 0.25, Form("= %+5.2f ", (t-L521)/tE));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  savePad("t521.pdf");

  h = ((TH1D*)gFile->Get("t511"));
  f = (TF1*)h->GetFunction("expo");
  chi2 = f->GetChisquare();
  ndf  = f->GetNDF();
  t    = -1./f->GetParameter(1);
  tE   = -t*f->GetParError(1)/f->GetParameter(1);
  h->SetTitle("");
  setTitles(h, "c#tau [#mum]", "Entries/bin");
  h->Draw();
  setItalic(); tl->DrawLatexNDC(0.2, 0.45, "B^{0}"); setRoman();
  tl->DrawLatexNDC(0.2, 0.40, Form("#chi^{2}/ndf")); tl->DrawLatexNDC(0.3, 0.40, Form("= %3.1f/%d", chi2, ndf));
  tl->DrawLatexNDC(0.2, 0.35, Form("#tau(fit)"));    tl->DrawLatexNDC(0.3, 0.35, Form("= %5.2f #pm %5.2f", t, tE));
  tl->DrawLatexNDC(0.2, 0.30, Form("#tau(pdg)"));    tl->DrawLatexNDC(0.3, 0.30, Form("= %5.2f ", L511));
  tl->DrawLatexNDC(0.2, 0.25, Form("pull"));         tl->DrawLatexNDC(0.3, 0.25, Form("= %+5.2f ", (t-L511)/tE));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  savePad("t511.pdf");

  h = ((TH1D*)gFile->Get("t531"));
  f = (TF1*)h->GetFunction("expo");
  chi2 = f->GetChisquare();
  ndf  = f->GetNDF();
  t    = -1./f->GetParameter(1);
  tE   = -t*f->GetParError(1)/f->GetParameter(1);
  h->SetTitle("");
  setTitles(h, "c#tau [#mum]", "Entries/bin");
  h->Draw();
  setItalic(); tl->DrawLatexNDC(0.2, 0.45, "B_{s}"); setRoman();
  tl->DrawLatexNDC(0.2, 0.40, Form("#chi^{2}/ndf")); tl->DrawLatexNDC(0.3, 0.40, Form("= %3.1f/%d", chi2, ndf));
  tl->DrawLatexNDC(0.2, 0.35, Form("#tau(fit)"));    tl->DrawLatexNDC(0.3, 0.35, Form("= %5.2f #pm %5.2f", t, tE));
  tl->DrawLatexNDC(0.2, 0.30, Form("#tau(pdg)"));    tl->DrawLatexNDC(0.3, 0.30, Form("= %5.2f ", L531));
  tl->DrawLatexNDC(0.2, 0.25, Form("pull"));         tl->DrawLatexNDC(0.3, 0.25, Form("= %+5.2f ", (t-L531)/tE));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  savePad("t531.pdf");

  h = ((TH1D*)gFile->Get("t5122"));
  f = (TF1*)h->GetFunction("expo");
  chi2 = f->GetChisquare();
  ndf  = f->GetNDF();
  t    = -1./f->GetParameter(1);
  tE   = -t*f->GetParError(1)/f->GetParameter(1);
  h->SetTitle("");
  setTitles(h, "c#tau [#mum]", "Entries/bin");
  h->Draw();
  setItalic(); tl->DrawLatexNDC(0.2, 0.45, "#Lambda_{b}"); setRoman();
  tl->DrawLatexNDC(0.2, 0.40, Form("#chi^{2}/ndf")); tl->DrawLatexNDC(0.3, 0.40, Form("= %3.1f/%d", chi2, ndf));
  tl->DrawLatexNDC(0.2, 0.35, Form("#tau(fit)"));    tl->DrawLatexNDC(0.3, 0.35, Form("= %5.2f #pm %5.2f", t, tE));
  tl->DrawLatexNDC(0.2, 0.30, Form("#tau(pdg)"));    tl->DrawLatexNDC(0.3, 0.30, Form("= %5.2f ", L5122));
  tl->DrawLatexNDC(0.2, 0.25, Form("pull"));         tl->DrawLatexNDC(0.3, 0.25, Form("= %+5.2f ", (t-L5122)/tE));
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  savePad("t5122.pdf");
}


// ----------------------------------------------------------------------
void plotWork::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotWork::loadFile loading files listed in " << files << endl;

  bool readMoreFiles(false);
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

    if (string::npos != stype.find("data")) {
      // -- DATA
      pF = loadFile(sfile);

      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("XXX")) {
        sname = "bdpsikstarData";
        sdecay = "bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
	fDS.insert(make_pair(sname, ds));
      }

    } else {
      // -- MC
      pF = loadFile(sfile);

      ds = new dataset();
      ds->fSize = 0.1;
      ds->fWidth = 2.;

      if (string::npos != stype.find("YYY")) {
        sname = "bupsikMc";
        sdecay = "bupsik";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
	fDS.insert(make_pair(sname, ds));
      }

    }
    if (sname != "nada") {
      readMoreFiles = true;
      ds->fLcolor = ds->fColor;
      ds->fFcolor = ds->fColor;
      ds->fName   = sdecay;
      ds->fFullName = sname;
      fDS.insert(make_pair(sname, ds));
    }


  }

  is.close();
  if (readMoreFiles) {
    for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
      cout << it->first << ": " << it->second->fName << ", " << it->second->fF->GetName() << endl;
    }
  }
}
