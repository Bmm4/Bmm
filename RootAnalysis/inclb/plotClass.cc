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

#include "common/util.hh"

#include "preselection.hh"

ClassImp(plotClass)

using namespace std; 

// ----------------------------------------------------------------------
plotClass::plotClass(string dir,  string files, string setup) {

  gStyle->SetHatchesSpacing(2);

  fDBX = true; 
  fDoUseBDT = false;
  fVerbose = true;

  fDirectory = dir; 
  fSetup = "A";
  fSuffix = setup; 
  

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

  fAccPt = 3.5;
  fAccEtaGen = 2.5;
  fAccEtaRec = 2.4;

  c0 = (TCanvas*)gROOT->FindObject("c0"); 
  if (!c0) c0 = new TCanvas("c0","--c0--",0,0,656,700);

  fHistFile = 0; // this must be opened in a derived class!

  fStampString = "preliminary";
  if (fDoUseBDT) {
    fStampString = "BDT preliminary"; 
  } else {
    fStampString = "CNC preliminary"; 
  }
  fStampCms = "BMM4";

  string sfiles(files);
  if (string::npos != sfiles.find("2011")) {
    fYear = 2011; 
    fStampCms = "L = 5 fb^{-1} (#sqrt{s} = 7 TeV)";
  } 
  if (string::npos != sfiles.find("2012")) {
    fYear = 2012; 
    fStampCms = "L = 20 fb^{-1} (#sqrt{s} = 8 TeV)";
  } 

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
    setTitles(h, h->GetXaxis()->GetTitle(), "normalized to 1", 1.1, 1.5);
  } else if (method == SOMETHING) {
    smethod = "something"; 
    scale = fNorm * (h->Integral() > 0 ? fNorm/h->Integral() : 1.); 
    setTitles(h, h->GetXaxis()->GetTitle(), "weighted events", 1.1, 1.5);
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
void plotClass::overlay(TH1* h1, string f1, TH1* h2, string f2, TH1* h3, string f3, int method, bool loga, bool legend, double xleg, double yleg) {
  const bool verbose(false);

  showOverflow(h1);
  showOverflow(h2);
  if (h3) showOverflow(h3);

  normHist(h1, f1, method); 
  normHist(h2, f2, method); 
  if (h3) normHist(h3, f3, method); 
  double ymin(0.0001);
  double h1max(h1->GetBinContent(h1->GetMaximumBin())); 
  double h2max(h2->GetBinContent(h2->GetMaximumBin())); 
  double h3max(h3->GetBinContent(h3->GetMaximumBin())); 
  double hmax(h1max);
  int imax(1);
  if (h2max > h1max) {
    hmax = h2max; 
    imax = 2;
  }
  if (h3max > h2max) {
    hmax = h3max; 
    imax = 3;
  }
  hmax *= 1.2; 
  if (verbose)  {
    cout << "hmax = " << hmax << " from imax = " << imax;
    if (1 == imax) cout << " bin " << h1->GetMaximumBin() << " with maximum " << h1->GetBinContent(h1->GetMaximumBin()) << endl;
    if (2 == imax) cout << " bin " << h2->GetMaximumBin() << " with maximum " << h2->GetBinContent(h2->GetMaximumBin()) << endl;
    if (3 == imax) cout << " bin " << h3->GetMaximumBin() << " with maximum " << h3->GetBinContent(h3->GetMaximumBin()) << endl;
  }
  if (loga) {
    gPad->SetLogy(1); 
    hmax *= 2.;
    double hmin(h1->GetMinimum(ymin)); 
    if (verbose) cout << "hmin1 = " << hmin << endl;
    if (h2->GetMinimum(ymin) < hmin) {
      hmin = h2->GetMinimum(ymin);
      if (verbose) cout << "hmin2 = " << hmin << endl;
    }
    if (h3 && h3->GetMinimum(ymin) < hmin) {
      hmin = h3->GetMinimum(ymin);
      if (verbose) cout << "hmin3 = " << hmin << endl;
    }
    h1->SetMinimum(0.1*hmin); 
    if (verbose) cout << "hmin = " << hmin << endl;
  } else {
    gPad->SetLogy(0); 
    h1->SetMinimum(0.); 
  }
  
  h1->SetMaximum(hmax); 

  h1->DrawCopy("hist"); 
  h2->DrawCopy("histsame");
  if (h3) h3->DrawCopy("histsame");

  if (legend) {
    newLegend(xleg, yleg, xleg+0.25, yleg+0.15); 
    legg->SetTextSize(0.03);
    string text; 
    text = fDS[f1]->fName.c_str(); 
    if (fDBX) {
      text = Form("%s: %4.3f#pm%4.3f, %4.3f", fDS[f1]->fName.c_str(), h1->GetMean(), h1->GetMeanError(), h1->GetRMS()); 
    }
    legg->AddEntry(h1, text.c_str(), "f"); 

    text = fDS[f2]->fName.c_str(); 
    if (fDBX) {
      text = Form("%s: %4.3f#pm%4.3f, %4.3f", fDS[f2]->fName.c_str(), h2->GetMean(), h2->GetMeanError(), h2->GetRMS()); 
    }
    legg->AddEntry(h2, text.c_str(), "f"); 

    if (h3) {
      text = fDS[f3]->fName.c_str(); 
      if (fDBX) {
	text = Form("%s: %4.3f#pm%4.3f, %4.3f", fDS[f3]->fName.c_str(), h3->GetMean(), h3->GetMeanError(), h3->GetRMS()); 
      }
      legg->AddEntry(h3, text.c_str(), "f"); 
    }

    legg->Draw();
  }


  if (verbose) cout << "==>plotClass: overlay(" << f1 << ", " << h1->GetName() << " integral= " << h1->Integral()
		    << ", " << f2 << ", " << h2->GetName() << " integral= " << h2->Integral()
		    << (h3? Form(", %s, %s, %f", f3.c_str(), h3->GetName(), h3->Integral()) : "")
		    << ")  log: " << loga << " legend = " << legend
		    << endl;
}

// ----------------------------------------------------------------------
void plotClass::overlay(string h1name, string f1, string h2name, string f2, string h3name, string f3, int method, bool loga, 
			bool legend, double xleg, double yleg) {

  cout << h1name << " from " << f1 << " vs. " << h2name << " from " << f2 << " vs. " << h3name << " from " << f3 << endl;

  TH1D *h1 = fDS[f1]->getHist(Form("%s", h1name.c_str()), true); 
  TH1D *h2 = fDS[f2]->getHist(Form("%s", h2name.c_str()), true); 
  TH1D *h3(0); 
  if (h3name != "") {
    h3 = fDS[f3]->getHist(Form("%s", h3name.c_str()), true); 
  } else {
    h3= 0; 
  }

  overlay(h1, f1, h2, f2, h3, f3, method, loga, legend, xleg, yleg); 
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
void plotClass::setupTree(TTree *t, string mode) {

}





// ----------------------------------------------------------------------
void plotClass::candAnalysis(int mode) {

}


// ----------------------------------------------------------------------
int plotClass::detChan(double m1eta, double m2eta) {
  // -- simple two channel analysis: channel 0 if both muons in barrel, channel 1 else
  if (TMath::Abs(m1eta) < fCuts[0]->etaMax && TMath::Abs(m2eta) < fCuts[0]->etaMax) return 0; 
  if (TMath::Abs(m1eta) < 2.4 && TMath::Abs(m2eta) < 2.4) return 1; 
  return -1; 
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
  //  if (legg) delete legg;
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


// ----------------------------------------------------------------------
void plotClass::calcBDT() {
}


// ----------------------------------------------------------------------
TMVA::Reader* plotClass::setupReader(string xmlFile, readerData &rd) {
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );
  
  TString dir    = "weights/";
  TString methodNameprefix = "BDT";
  //  TString methodName = TString(fBdt) + TString(" method");
  //  TString weightfile = dir + fBdt + "_" + methodNameprefix + TString(".weights.xml");
  TString weightfile = xmlFile;

  // -- read in variables from weight file
  vector<string> allLines; 
  char  buffer[2000];
  cout << "setupReader, open file " << weightfile << endl;
  ifstream is(weightfile); 
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  int nvars(-1); 
  string::size_type m1, m2;
  string stype; 
  cout << "  read " << allLines.size() << " lines " << endl;
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add variables
    if (string::npos != allLines[i].find("Variables NVar")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      cout << "  " << stype << " variables" << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
        m1 = allLines[j].find("Expression=\"")+10; 
        m2 = allLines[j].find("\" Label=\"");
        stype = allLines[j].substr(m1+2, m2-m1-2); 
        //      cout << "ivar " << j-i << " variable string: ->" << stype << "<-" << endl;
        if (stype == "m1pt") {
          cout << "  adding m1pt" << endl;
          reader->AddVariable( "m1pt", &rd.m1pt);
        }
        if (stype == "m2pt") {
          cout << "  adding m2pt" << endl;
          reader->AddVariable( "m2pt", &rd.m2pt);
        }
        if (stype == "m1eta") {
          cout << "  adding m1eta" << endl;
          reader->AddVariable( "m1eta", &rd.m1eta);
        }
        if (stype == "m2eta") {
          reader->AddVariable( "m2eta", &rd.m2eta);
          cout << "  adding m2eta" << endl;
        }
        if (stype == "pt") {
          cout << "  adding pt" << endl;
          reader->AddVariable( "pt", &rd.pt);
        }
        if (stype == "eta") {
          cout << "  adding eta" << endl;
          reader->AddVariable( "eta", &rd.eta);
        }
        if (stype == "fls3d") {
          cout << "  adding fls3d" << endl;
          reader->AddVariable( "fls3d", &rd.fls3d);
        }
        if (stype == "alpha") {
          cout << "  adding alpha" << endl;
          reader->AddVariable( "alpha", &rd.alpha);
        }
        if (stype == "maxdoca") {
          cout << "  adding maxdoca" << endl;
          reader->AddVariable( "maxdoca", &rd.maxdoca);
        }
        if (stype == "pvip") {
          cout << "  adding pvip" << endl;
          reader->AddVariable( "pvip", &rd.pvip);
        }
        if (stype == "pvips") {
          cout << "  adding pvips" << endl;
          reader->AddVariable( "pvips", &rd.pvips);
        }
        if (stype == "iso") {
          cout << "  adding iso" << endl;
          reader->AddVariable( "iso", &rd.iso);
        }
        if (stype == "docatrk") {
          cout << "  adding docatrk" << endl;
          reader->AddVariable( "docatrk", &rd.docatrk);
        }
        if (stype == "closetrk") {
          cout << "  adding closetrk" << endl;
          reader->AddVariable( "closetrk", &rd.closetrk);
        }
        if (stype == "chi2dof") {
          cout << "  adding chi2dof" << endl;
          reader->AddVariable( "chi2dof", &rd.chi2dof);
        }
        if (stype == "closetrks1") {
          cout << "  adding closetrks1" << endl;
          reader->AddVariable( "closetrks1", &rd.closetrks1);
        }
        if (stype == "closetrks2") {
          cout << "  adding closetrks2" << endl;
          reader->AddVariable( "closetrks2", &rd.closetrks2);
        }
        if (stype == "closetrks3") {
          cout << "  adding closetrks3" << endl;
          reader->AddVariable( "closetrks3", &rd.closetrks3);
        }
        if (stype == "m1iso") {
          cout << "  adding m1iso" << endl;
          reader->AddVariable( "m1iso", &rd.m1iso);
        }
        if (stype == "m2iso") {
          cout << "  adding m2iso" << endl;
          reader->AddVariable( "m2iso", &rd.m2iso);
        }
        if (stype == "othervtx") {
          cout << "  adding othervtx" << endl;
          reader->AddVariable( "othervtx", &rd.othervtx);
        }
        if (stype == "pvdchi2") {
          cout << "  adding pvdchi2" << endl;
          reader->AddVariable( "pvdchi2", &rd.pvdchi2);
        }
        if (stype == "pvlip2") {
          cout << "  adding pvlip2" << endl;
          reader->AddVariable( "pvlip2", &rd.pvlip2);
        }
        if (stype == "pvlips2") {
          cout << "  adding pvlips2" << endl;
          reader->AddVariable( "pvlips2", &rd.pvlips2);
        }
      }
      break;
    }
  }
  
  nvars = -1; 
  for (unsigned int i = 0; i < allLines.size(); ++i) {
    // -- parse and add spectators
    if (string::npos != allLines[i].find("Spectators NSpec")) {
      m1 = allLines[i].find("=\""); 
      stype = allLines[i].substr(m1+2, allLines[i].size()-m1-2-2); 
      //      cout << "==> " << stype << endl;
      nvars = atoi(stype.c_str());
      if (-1 == nvars) continue;
      for (unsigned int j = i+1; j < i+nvars+1; ++j) {
        m1 = allLines[j].find("Expression=\"")+10; 
        m2 = allLines[j].find("\" Label=\"");
        stype = allLines[j].substr(m1+2, m2-m1-2); 
        cout << "ivar " << j-i << " spectator string: ->" << stype << "<-" << endl;
        if (stype == "m") {
          cout << "  adding m as spectator" << endl;
          reader->AddSpectator( "m", &rd.m);  
        }
      }
      break;
    }
  }

  // --- Book the MVA methods
  reader->BookMVA("BDT", weightfile); 
  return reader; 
}


// ----------------------------------------------------------------------
void plotClass::readCuts(string filename) {
  cout << "==> plotClass: Reading " << filename << " for cut settings" << endl;
  vector<string> cutLines; 
  char  buffer[200];
  ifstream is(filename);
  while (is.getline(buffer, 200, '\n')) {
    cutLines.push_back(string(buffer));
  }

  char CutName[100], XmlName[1000];
  float CutValue;
  int dump(1), ok(0);

  cuts *a = 0;
  
  fCuts.clear();

  for (unsigned int i = 0; i < cutLines.size(); ++i) {
    sprintf(buffer, "%s", cutLines[i].c_str()); 
    
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "index")) {
      ok = 1;
      if (dump) cout << "index:            " << CutValue << endl;
      if (a) fCuts.push_back(a); 
      a = new cuts; 
      a->index = static_cast<int>(CutValue); 
    }
    
    if (!strcmp(CutName, "mBdLo")) {
      a->mBdLo = CutValue; ok = 1;
      if (dump) cout << "mBdLo:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "bdtpt")) {
      a->bdtpt = CutValue; ok = 1;
      if (dump) cout << "bdtpt:              " << CutValue << endl;
    }

    if (!strcmp(CutName, "bdt")) {
      a->bdt = CutValue; ok = 1;
      if (dump) cout << "bdt:              " << CutValue << endl;
    }

    if (!strcmp(CutName, "bdtMax")) {
      a->bdtMax = CutValue; ok = 1;
      if (dump) cout << "bdtMax:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBdHi")) {
      a->mBdHi = CutValue; ok = 1;
      if (dump) cout << "mBdHi:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBsLo")) {
      a->mBsLo = CutValue; ok = 1;
      if (dump) cout << "mBsLo:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "mBsHi")) {
      a->mBsHi = CutValue; ok = 1;
      if (dump) cout << "mBsHi:            " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMin")) {
      a->etaMin = CutValue; ok = 1;
      if (dump) cout << "etaMin:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMax")) {
      a->etaMax = CutValue; ok = 1;
      if (dump) cout << "etaMax:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "pt")) {
      a->pt = CutValue; ok = 1;
      if (dump) cout << "pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1pt")) {
      a->m1pt = CutValue; ok = 1;
      if (dump) cout << "m1pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2pt")) {
      a->m2pt = CutValue; ok = 1;
      if (dump) cout << "m2pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1eta")) {
      a->m1eta = CutValue; ok = 1;
      if (dump) cout << "m1eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2eta")) {
      a->m2eta = CutValue; ok = 1;
      if (dump) cout << "m2eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "iso")) {
      a->iso = CutValue; ok = 1;
      if (dump) cout << "iso:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "chi2dof")) {
      a->chi2dof = CutValue; ok = 1;
      if (dump) cout << "chi2dof:             " << CutValue << endl;
    }

    if (!strcmp(CutName, "alpha")) {
      a->alpha = CutValue; ok = 1;
      if (dump) cout << "alpha:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "fls3d")) {
      a->fls3d = CutValue; ok = 1;
      if (dump) cout << "fls3d:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "docatrk")) {
      a->docatrk = CutValue; ok = 1;
      if (dump) cout << "docatrk:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "closetrk")) {
      a->closetrk = static_cast<int>(CutValue); ok = 1;
      if (dump) cout << "closetrk:              " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlip")) {
      a->pvlip = CutValue; ok = 1;
      if (dump) cout << "pvlip:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "etaMax")) {
      a->etaMax = CutValue; ok = 1;
      if (dump) cout << "etaMax:           " << CutValue << endl;
    }

    if (!strcmp(CutName, "pt")) {
      a->pt = CutValue; ok = 1;
      if (dump) cout << "pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1pt")) {
      a->m1pt = CutValue; ok = 1;
      if (dump) cout << "m1pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2pt")) {
      a->m2pt = CutValue; ok = 1;
      if (dump) cout << "m2pt:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m1eta")) {
      a->m1eta = CutValue; ok = 1;
      if (dump) cout << "m1eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "m2eta")) {
      a->m2eta = CutValue; ok = 1;
      if (dump) cout << "m2eta:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "iso")) {
      a->iso = CutValue; ok = 1;
      if (dump) cout << "iso:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "chi2dof")) {
      a->chi2dof = CutValue; ok = 1;
      if (dump) cout << "chi2dof:             " << CutValue << endl;
    }

    if (!strcmp(CutName, "alpha")) {
      a->alpha = CutValue; ok = 1;
      if (dump) cout << "alpha:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "fls3d")) {
      a->fls3d = CutValue; ok = 1;
      if (dump) cout << "fls3d:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "docatrk")) {
      a->docatrk = CutValue; ok = 1;
      if (dump) cout << "docatrk:               " << CutValue << endl;
    }

    if (!strcmp(CutName, "closetrk")) {
      a->closetrk = static_cast<int>(CutValue); ok = 1;
      if (dump) cout << "closetrk:              " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlip")) {
      a->pvlip = CutValue; ok = 1;
      if (dump) cout << "pvlip:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlips")) {
      a->pvlips = CutValue; ok = 1;
      if (dump) cout << "pvlips:                " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlip2")) {
      a->pvlip2 = CutValue; ok = 1;
      if (dump) cout << "pvlip2:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvlips2")) {
      a->pvlips2 = CutValue; ok = 1;
      if (dump) cout << "pvlips2:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "maxdoca")) {
      a->maxdoca = CutValue; ok = 1;
      if (dump) cout << "maxdoca:                 " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvip")) {
      a->pvip = CutValue; ok = 1;
      if (dump) cout << "pvip:                    " << CutValue << endl;
    }

    if (!strcmp(CutName, "pvips")) {
      a->pvips = CutValue; ok = 1;
      if (dump) cout << "pvips:                   " << CutValue << endl;
    }

    sscanf(buffer, "%s %s", CutName, XmlName);
    string ctmp = CutName; 
    string sXmlName;
    replaceAll(ctmp, " ", ""); 
    if (!strcmp(ctmp.c_str(), "xml")) {
      a->xmlFile = XmlName; ok = 1; 
      sXmlName = "weights/" + a->xmlFile + "-Events0_BDT.weights.xml"; 
      //      fReaderEvents0.push_back(setupReader(sXmlName, frd)); 
      TMVA::Reader *ar = setupReader(sXmlName, frd); 
      fReaderEvents0[a->index] = ar;
      if (dump) cout << "xml:                   " << sXmlName << endl;
      sXmlName = "weights/" + a->xmlFile + "-Events1_BDT.weights.xml"; 
      //      fReaderEvents1.push_back(setupReader(sXmlName, frd)); 
      ar = setupReader(sXmlName, frd); 
      fReaderEvents1[a->index] = ar;
      if (dump) cout << "xml:                   " << sXmlName << endl;
      sXmlName = "weights/" + a->xmlFile + "-Events2_BDT.weights.xml"; 
      //      fReaderEvents2.push_back(setupReader(sXmlName, frd)); 
      ar = setupReader(sXmlName, frd); 
      fReaderEvents2[a->index] = ar;
      if (dump) cout << "xml:                   " << sXmlName << endl;
    }

    if (!ok) cout << "==> what about " << CutName << endl;
  }

  if (a) fCuts.push_back(a); 

  cout << "==> finished reading cut setting, fCuts.size() =  " << fCuts.size() << endl;
  
}



// ----------------------------------------------------------------------
void plotClass::printCuts(ostream &OUT) {

  OUT << "----------------------------------------------------------------------" << endl;
  cout << "printCuts ... fCuts.size() = " << fCuts.size() << endl;
  for (unsigned int i = 0; i < fCuts.size(); ++i) {
    cuts *a = fCuts[i]; 
    OUT << "# -- channel " << a->index << endl;
    OUT << "index    " << a->index << endl;
    fTEX << "% ----------------------------------------------------------------------" << endl;
    fTEX << "% -- Cuts for channel " << a->index << endl;

    OUT << "xml      " << Form("%s", a->xmlFile.c_str()) << endl;
    OUT << "bdt      " << Form("%4.3f", a->bdt) << endl;
    OUT << "bdtMax   " << Form("%4.3f", a->bdtMax) << endl;
    fTEX <<  Form("\\vdef{%s:bdt:%d}     {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->bdt) << endl;
    fTEX <<  Form("\\vdef{%s:bdtMax:%d}  {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->bdtMax) << endl;

    OUT << "mBdLo    " << Form("%4.3f", a->mBdLo) << endl;
    OUT << "mBdHi    " << Form("%4.3f", a->mBdHi) << endl;
    fTEX <<  Form("\\vdef{%s:mBdLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBdLo) << endl;
    fTEX <<  Form("\\vdef{%s:mBdHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBdHi) << endl;

    OUT << "mBsLo    " << Form("%4.3f", a->mBsLo) << endl;
    OUT << "mBsHi    " << Form("%4.3f", a->mBsHi) << endl;
    fTEX <<  Form("\\vdef{%s:mBsLo:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBsLo) << endl;
    fTEX <<  Form("\\vdef{%s:mBsHi:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->mBsHi) << endl;
 
    OUT << "etaMin   " << Form("%3.1f", a->etaMin) << endl;
    OUT << "etaMax   " << Form("%3.1f", a->etaMax) << endl;
    fTEX <<  Form("\\vdef{%s:etaMin:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->etaMin) << endl;
    fTEX <<  Form("\\vdef{%s:etaMax:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->etaMax) << endl;

    OUT << "pt       " << Form("%3.1f", a->pt) << endl;
    fTEX <<  Form("\\vdef{%s:pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->pt) << endl;
    OUT << "m1pt     " << a->m1pt << endl;
    fTEX <<  Form("\\vdef{%s:m1pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m1pt) << endl;
    OUT << "m2pt     " << a->m2pt << endl;
    fTEX <<  Form("\\vdef{%s:m2pt:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m2pt) << endl;
    OUT << "m1eta    " << a->m1eta << endl;
    fTEX <<  Form("\\vdef{%s:m1eta:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m1eta) << endl;
    OUT << "m2eta    " << a->m2eta << endl;
    fTEX <<  Form("\\vdef{%s:m2eta:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->m2eta) << endl;

    OUT << "iso      " << a->iso << endl;
    fTEX <<  Form("\\vdef{%s:iso:%d}   {\\ensuremath{{%3.2f } } }", fSuffix.c_str(), a->index, a->iso) << endl;
    OUT << "chi2dof  " << a->chi2dof << endl;
    fTEX <<  Form("\\vdef{%s:chi2dof:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->chi2dof) << endl;
    OUT << "alpha    " << a->alpha << endl;
    fTEX <<  Form("\\vdef{%s:alpha:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->alpha) << endl;
    OUT << "fls3d    " << a->fls3d << endl;
    fTEX <<  Form("\\vdef{%s:fls3d:%d}   {\\ensuremath{{%3.1f } } }", fSuffix.c_str(), a->index, a->fls3d) << endl;
    OUT << "docatrk  " << a->docatrk << endl;
    fTEX <<  Form("\\vdef{%s:docatrk:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->docatrk) << endl;

    OUT << "closetrk " << a->closetrk << endl;
    fTEX <<  Form("\\vdef{%s:closetrk:%d}   {\\ensuremath{{%d } } }", fSuffix.c_str(), a->index, static_cast<int>(a->closetrk)) << endl;
    OUT << "pvlip    " << a->pvlip << endl;
    fTEX <<  Form("\\vdef{%s:pvlip:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvlip) << endl;
    OUT << "pvlips   " << a->pvlips << endl;
    fTEX <<  Form("\\vdef{%s:pvlips:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvlips) << endl;
    OUT << "pvlip2   " << a->pvlip2 << endl;
    fTEX <<  Form("\\vdef{%s:pvlip2:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvlip2) << endl;
    OUT << "pvlips2  " << a->pvlips2 << endl;
    fTEX <<  Form("\\vdef{%s:pvlips2:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvlips2) << endl;
    OUT << "maxdoca  " << a->maxdoca << endl;
    fTEX <<  Form("\\vdef{%s:maxdoca:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->maxdoca) << endl;
    OUT << "pvip     " << a->pvip << endl;
    fTEX <<  Form("\\vdef{%s:pvip:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvip) << endl;
    OUT << "pvips    " << a->pvips << endl;
    fTEX <<  Form("\\vdef{%s:pvips:%d}   {\\ensuremath{{%4.3f } } }", fSuffix.c_str(), a->index, a->pvips) << endl;

  }
  OUT.flush();
}

#include "plotClass.icc"


