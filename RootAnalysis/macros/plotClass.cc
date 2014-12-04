#include "plotClass.hh"

#include <fstream>
#include <iostream>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"

#include "util.hh"

ClassImp(plotClass)

using namespace std; 

// ----------------------------------------------------------------------
plotClass::plotClass(string dir,  string files, string setup) {

  fDBX = true; 
  fVerbose = true;

  fDirectory = dir; 
  fSetup = setup;

  delete gRandom;
  gRandom = (TRandom*) new TRandom3;

  fEpsilon = 0.00001; 
  fLumi = 20.; 

  legg = 0;
  c0 = c1 = c2 = c3 = c4 = c5 =0;
  tl = new TLatex();
  box = new TBox();
  pa = new TArrow();
  pl = new TLine(); 
  legge = 0;

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);

  fHistFile = 0; // this must be opened in a derived class!
}

// ----------------------------------------------------------------------
plotClass::~plotClass() {
}

// ----------------------------------------------------------------------
// see http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=15054
void plotClass::closeHistFile() {
  fHistFile->Write(); 
}

// ----------------------------------------------------------------------
void plotClass::cd(std::string dataset, std::string dir) {
  if (0 == fDS.count(dataset)) {
    cout << "unknown dataset: " << dataset << endl;
  } else {
    fDS[dataset]->cd(dir.c_str());
  }
}

// ----------------------------------------------------------------------
void plotClass::bookHist(string name) {
  cout << "==> plotClass: bookHist " << name << endl;
}

// ----------------------------------------------------------------------
void plotClass::makeAll(int bitmask) {
  cout << "==> plotClass: makeAll " << bitmask << endl;
}

// ----------------------------------------------------------------------
void plotClass::treeAnalysis() {
  cout << "==> plotClass: treeAnalysis " << endl;
}

// ----------------------------------------------------------------------
void plotClass::normHist(TH1 *h, string ds, int method) {
  double scale(1.); 
  string smethod("");
  // -- normalize to 1
  if (method == UNITY) {
    smethod = "unity"; 
    scale = (h->Integral() > 0 ? 1./h->Integral() : 1.); 
    setTitles(h, h->GetXaxis()->GetTitle(), "normalized to 1");
  } else if (method == SOMETHING) {
    smethod = "something"; 
    scale = fNorm * (h->Integral() > 0 ? fNorm/h->Integral() : 1.); 
    setTitles(h, h->GetXaxis()->GetTitle(), "weighted events");
  } else if (method == XSECTION) {
    smethod = "xsection"; 
    // -- normalize to EFFECTIVE xsec*bf (EFFECTIVE to account for cuts)
    //    the cross section is known for ds
    //    ds corresponds to know lumi
    //    
    //    n = xsec * L
    //    "integral" over histogram should be EFFECTIVE xsec
    scale = (h->Integral() > 0 ? fDS[ds]->fXsec*fDS[ds]->fBf/h->Integral() : 1.); 
    setTitles(h, h->GetXaxis()->GetTitle(), "pb");
  } else if (method == LUMI) {
    smethod = "lumi"; 
    // -- normalize to xsec*bf
    //    n = xsec * L
    //    "integral" over histogram should be events expected in fLumi
    scale = (h->Integral() > 0 ? fLumi/fDS[ds]->fLumi : 1.); 
    setTitles(h, h->GetXaxis()->GetTitle(), Form("events in %4.0f/fb", fLumi));
  } else if (method == NONORM) {
    smethod = "nonorm"; 
    scale = 1.;
  } else {
    scale = 1.;
  }

  cout << "==>plotClass:  normHist scaling by " << scale << ", based on method " << smethod << endl;

  double c(0.), e(0.); 
  for (int i = 0; i <= h->GetNbinsX(); ++i) {
    c = h->GetBinContent(i); 
    e = h->GetBinError(i); 
    h->SetBinContent(i, c*scale);
    h->SetBinError(i, e*scale);
  }
}


// ----------------------------------------------------------------------
void plotClass::overlayAll() {
}

// ----------------------------------------------------------------------
void plotClass::overlay(TH1* h1, string f1, TH1* h2, string f2, int method, bool log, bool legend, double xleg, double yleg) {

  normHist(h1, f1, method); 
  normHist(h2, f2, method); 

  double hmax(1.2*h1->GetMaximum()); 
  if (h2->GetMaximum() > hmax) hmax = 1.2*h2->GetMaximum(); 
  if (log) {
    gPad->SetLogy(1); 
    hmax *= 2.;
    h1->SetMinimum(0.1*h1->GetMinimum(1.e-6));
  } else {
    h1->SetMinimum(0.); 
  }
  h1->SetMaximum(hmax); 

  h1->DrawCopy("hist"); 
  h2->DrawCopy("histsame");
  cout << "==> plotClass: overlay(" << f1 << ", " << h1->GetName() << " integral= " << h1->Integral()
       << ", " << f2 << ", " << h2->GetName() << " integral= " << h2->Integral()
       << ") legend = " << legend << " log: " << log 
       << endl;
  
  if (legend) {
    newLegend(xleg, yleg, xleg+0.25, yleg+0.10); 
    legg->SetTextSize(0.03);
    legg->AddEntry(h1, fDS[f1]->fName.c_str(), "f"); 
    legg->AddEntry(h2, fDS[f2]->fName.c_str(), "f"); 
    legg->Draw();
    if (fDBX) {
      tl->SetNDC(kTRUE);
      tl->SetTextSize(0.05);
      tl->SetTextColor(fDS[f1]->fColor); 
      tl->DrawLatex(0.15, 0.92, Form("%.2e", h1->Integral())); 
      tl->SetTextColor(fDS[f2]->fColor); 
      tl->DrawLatex(0.40, 0.92, Form("%.2e", h2->Integral())); 
    }
  }
}

// ----------------------------------------------------------------------
void plotClass::overlay(string h1name, string f1, string h2name, string f2, int method, bool log, bool legend, double xleg, double yleg) {
  TH1D *h1 = fDS[f1]->getHist(Form("%s", h1name.c_str())); 
  TH1D *h2 = fDS[f2]->getHist(Form("%s", h2name.c_str())); 
  overlay(h1, f1, h2, f2, method, log, legend); 
}

// ----------------------------------------------------------------------
void plotClass::overlay(TH1* h1, string f1, TH1* h2, string f2, TH1* h3, string f3, int method, bool log, bool legend, double xleg, double yleg) {

  normHist(h1, f1, method); 
  normHist(h2, f2, method); 
  normHist(h3, f3, method); 
  double ymin(0.0001);
  double hmax(1.2*h1->GetMaximum()); 
  if (h2->GetMaximum() > hmax) hmax = 1.2*h2->GetMaximum(); 
  if (h3->GetMaximum() > hmax) hmax = 1.2*h3->GetMaximum(); 
  if (log) {
    gPad->SetLogy(1); 
    hmax *= 2.;
    double hmin(h1->GetMinimum(ymin)); 
    cout << "hmin1 = " << hmin << endl;
    if (h2->GetMinimum(ymin) < hmin) {
      hmin = h2->GetMinimum(ymin);
      cout << "hmin2 = " << hmin << endl;
    }
    if (h3->GetMinimum(ymin) < hmin) {
      hmin = h3->GetMinimum(ymin);
      cout << "hmin3 = " << hmin << endl;
    }
    h1->SetMinimum(0.1*hmin); 
    cout << "hmin = " << hmin << endl;
  } else {
    h1->SetMinimum(0.); 
  }
  h1->SetMaximum(hmax); 

  h1->DrawCopy("hist"); 
  h2->DrawCopy("histsame");
  h3->DrawCopy("histsame");
  cout << "==> plotClass: overlay(" << f1 << ", " << h1->GetName() << " integral= " << h1->Integral()
       << ", " << f2 << ", " << h2->GetName() << " integral= " << h2->Integral()
       << ", " << f3 << ", " << h3->GetName() << " integral= " << h3->Integral()
       << ") legend = " << legend << " log: " << log 
       << endl;
  
  if (legend) {
    newLegend(xleg, yleg, xleg+0.25, yleg+0.15); 
    legg->SetTextSize(0.03);
    legg->AddEntry(h1, fDS[f1]->fName.c_str(), "f"); 
    legg->AddEntry(h2, fDS[f2]->fName.c_str(), "f"); 
    legg->AddEntry(h3, fDS[f3]->fName.c_str(), "f"); 
    legg->Draw();
  }
  if (fDBX) {
    tl->SetNDC(kTRUE);
    tl->SetTextSize(0.05);
    tl->SetTextColor(fDS[f1]->fColor); 
    tl->DrawLatex(0.15, 0.92, Form("%.2e", h1->Integral())); 
    tl->SetTextColor(fDS[f2]->fColor); 
    tl->DrawLatex(0.40, 0.92, Form("%.2e", h2->Integral())); 
    tl->SetTextColor(fDS[f3]->fColor); 
    tl->DrawLatex(0.65, 0.92, Form("%.2e", h3->Integral())); 
  }
}

// ----------------------------------------------------------------------
void plotClass::overlay(string h1name, string f1, string h2name, string f2, string h3name, string f3, int method, bool log, 
			bool legend, double xleg, double yleg) {
  TH1D *h1 = fDS[f1]->getHist(Form("%s", h1name.c_str())); 
  TH1D *h2 = fDS[f2]->getHist(Form("%s", h2name.c_str())); 
  TH1D *h3 = fDS[f3]->getHist(Form("%s", h3name.c_str())); 
  overlay(h1, f1, h2, f2, h3, f3, method, log, legend); 
}


// ----------------------------------------------------------------------
void plotClass::loopFunction1() {
}


// ----------------------------------------------------------------------
void plotClass::loopFunction2() {
}

// ----------------------------------------------------------------------
void plotClass::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  if (2 == ifunc)          step = 10000; 
  cout << "==> plotClass::loopOverTree> loop over dataset " << fCds << " in file " 
       << t->GetDirectory()->GetName() 
       << " with " << nentries << " entries"  << " looping from  " << nbegin << " .. " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  void (plotClass::*pF)(void);
  if (ifunc == 1) pF = &plotClass::loopFunction1;
  if (ifunc == 2) pF = &plotClass::loopFunction2;

  cout << "pF: " << pF << endl;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotClass::setupTree(TTree *t) {
}

// ----------------------------------------------------------------------
TTree* plotClass::getTree(string ds, string dir) {
  TTree *t(0);
  if (!dir.compare("")) {
    t = (TTree*)fDS[ds]->fF->Get("events"); 
  } else {
    t = (TTree*)fDS[ds]->fF->Get(Form("%s/events", dir.c_str())); 
  }
  return t; 
}

// ----------------------------------------------------------------------
TFile* plotClass::loadFile(string file) {
  TFile *f = TFile::Open(file.c_str());
  return f; 
}



// ----------------------------------------------------------------------
void plotClass::replaceAll(string &sInput, const string &oldString, const string &newString) {
  string::size_type foundpos = sInput.find(oldString);
  while (foundpos != string::npos)  {
    sInput.replace(sInput.begin() + foundpos, sInput.begin() + foundpos + oldString.length(), newString);
    foundpos = sInput.find(oldString);
  }
}

// ----------------------------------------------------------------------
void plotClass::newLegend(double x1, double y1, double x2, double y2, string title) {
  if (legg) delete legg;
  legg = new TLegend(x1, y1, x2, y2, title.c_str());
  legg->SetFillStyle(0); 
  legg->SetBorderSize(0); 
  legg->SetTextSize(0.04);  
  legg->SetFillColor(0); 
  legg->SetTextFont(42); 
}

// ----------------------------------------------------------------------
void plotClass::makeCanvas(int i) {
  if (i & 16) { 
    c5 = new TCanvas("c5", "c5", 210,   0, 800, 900);
    c5->ToggleEventStatus();
  }
  if (i & 8) { 
    c4 = new TCanvas("c4", "c4", 210,   0, 800, 600);
    c4->ToggleEventStatus();
  }
  if (i & 4) {
    c3 = new TCanvas("c3", "c3", 200,  20, 800, 800);
    c3->ToggleEventStatus();
  }
  if (i & 1) {
    //    c1 = new TCanvas("c1", "c1", 20,  60, 1200, 400);
    c1 = new TCanvas("c1", "c1", 20,  60, 1000, 400);
    c1->ToggleEventStatus();
  }
  if (i & 2) { 
    c2 = new TCanvas("c2", "c2", 300, 200, 400, 800);
    c2->ToggleEventStatus();
  }
}

#include "plotClass.icc"
