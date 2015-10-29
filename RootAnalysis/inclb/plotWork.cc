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
#include "TPad.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFractionFitter.h"

#include "common/dataset.hh"
#include "common/util.hh"

ClassImp(plotWork)

using namespace std; 

// ----------------------------------------------------------------------
plotWork::plotWork(string dir,  string files, string setup): plotClass(dir, files, setup) {
  loadFiles(files);

  if (setup == "") {
    fHistFileName = Form("%s/plotWork.root", dir.c_str()); 
  } else {
    fHistFileName = Form("%s/plotWork-%s.root", dir.c_str(), setup.c_str()); 
  }

  fTexFileName = fHistFileName; 
  replaceAll(fTexFileName, ".root", ".tex"); 
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));

}


// ----------------------------------------------------------------------
plotWork::~plotWork() {

}


// ----------------------------------------------------------------------
void plotWork::makeAll(int bitmask) {

  if (bitmask & 0x1) {
    zone(2,2); 
    validation("ptrel_0_6", "candAnaMuHIL2Mu3", "data_HIL2Mu3", "bSignalLowPU", "cSignalLowPU"); 
    c0->cd(2);
    validation("ptrel_6_8", "candAnaMuHIL2Mu3", "data_HIL2Mu3", "bSignalLowPU", "cSignalLowPU"); 
    c0->cd(3);
    validation("ptrel_8_10", "candAnaMuHIL2Mu3", "data_HIL2Mu3", "bSignalLowPU", "cSignalLowPU"); 
    c0->cd(4);
    validation("ptrel_10_15", "candAnaMuHIL2Mu3", "data_HIL2Mu3", "bSignalLowPU", "cSignalLowPU"); 
  } else if (bitmask & 0x2) {
    zone(2,2); 
    validation("ptrel_8_10", "candAnaMu8", "data_Mu8", "bSignalMu8", "cSignalMu8"); 
    c0->cd(2);
    validation("ptrel_10_15", "candAnaMu8", "data_Mu8", "bSignalMu8", "cSignalMu8"); 
    c0->cd(3);
    validation("ptrel_15_20", "candAnaMu8", "data_Mu8", "bSignalMu8", "cSignalMu8"); 
    c0->cd(4);
    validation("ptrel_20_30", "candAnaMu8", "data_Mu8", "bSignalMu8", "cSignalMu8"); 
  } else if (bitmask & 0x4) {
    zone(2,2); 
    validation("ptrel_8_10", "candAnaMu24", "data_Mu24", "bSignalMu8", "cSignalMu8"); 
    c0->cd(4);
    validation("ptrel_10_15", "candAnaMu24", "data_Mu24", "bSignalMu8", "cSignalMu8"); 
    c0->cd(5);
    validation("ptrel_15_20", "candAnaMu24", "data_Mu24", "bSignalMu8", "cSignalMu8"); 
    c0->cd(6);
    validation("ptrel_20_30", "candAnaMu24", "data_Mu24", "bSignalMu8", "cSignalMu8"); 
  }

}



// ----------------------------------------------------------------------
void plotWork::validation(string hist, string dir, string dname, string bname, string cname) {

  TH1D *hd = fDS[dname]->getHist(Form("%s/%s", dir.c_str(), hist.c_str()));
  TH1D *hb = fDS[bname]->getHist(Form("%s/%s", dir.c_str(), hist.c_str()));
  TH1D *hc = fDS[cname]->getHist(Form("%s/%s", dir.c_str(), hist.c_str()));
  
  TObjArray *mc = new TObjArray(2);        // MC histograms are put in this array
  mc->Add(hb);
  mc->Add(hc);
  TFractionFitter* fit = new TFractionFitter(hd, mc); // initialise
  //  fit->Constrain(1, 0.0, 1.0);               // constrain fraction 1 to be between 0 and 1
  fit->SetRangeX(1, 50);                    // use only the first 15 bins in the fit
  Int_t status = fit->Fit();               // perform the fit
  cout << "fit status: " << status << endl;
  if (status == 0) {                       // check on fit status
    double fracB, fracC, err;
    hd->Draw("Ep");
    fit->GetResult(0, fracB, err);
    fit->GetResult(1, fracC, err);
    hb->Scale(fracB*hd->GetSumOfWeights()/hb->GetSumOfWeights());
    hc->Scale(fracC*hd->GetSumOfWeights()/hc->GetSumOfWeights());
    hb->Draw("samehist");
    hc->Draw("samehist");
    TH1D *result = (TH1D*)hb->Clone("sum"); 
    result->SetLineColor(kBlack); 
    result->SetFillStyle(0);
    result->Add(hc); 
    result->Draw("samehist");

  }




}


// ----------------------------------------------------------------------
void plotWork::loopFunction1() {

}

// ----------------------------------------------------------------------
void plotWork::bookHist(int mode) {

}



// ----------------------------------------------------------------------
void plotWork::candAnalysis() {
  fGoodCand = true; 
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
void plotWork::setupTree(TTree *t) {

}


// ----------------------------------------------------------------------
void plotWork::setCuts(string cuts) {
  cout << "==> plotWork::setCuts: " << cuts << endl;

  istringstream ss(cuts);
  string token, name, sval;

  while (getline(ss, token, ',')) {
    
    string::size_type m1 = token.find("="); 
    name = token.substr(0, m1);
    sval = token.substr(m1+1);

    if (string::npos != name.find("PTLO")) {
      float val; 
      val = atof(sval.c_str()); 
      PTLO = val;
    }

  }
}


// ----------------------------------------------------------------------
void plotWork::loadFiles(string afiles) {
  
  string files = fDirectory + string("/") + afiles;
  cout << "==> Loading files listed in " << files << endl;

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
    string sname, sdecay; 
    
    cout << "stype: ->" << stype << "<-" << endl;
    
    TFile *pF(0); 
    if (string::npos != stype.find("data")) {
      // -- DATA
      pF = loadFile(sfile); 
      
      dataset *ds = new dataset(); 
      ds->fSize = 1; 
      ds->fWidth = 2; 
      
      if (string::npos != stype.find("HIL2Mu3")) {
        sname = "data_HIL2Mu3"; 
        sdecay = "HIL2Mu3"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("Mu8")) {
        sname = "data_Mu8"; 
        sdecay = "Mu8"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("Mu24")) {
        sname = "data_Mu24"; 
        sdecay = "Mu24"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }

      if (string::npos != stype.find("Mu50")) {
        sname = "data_Mu24"; 
        sdecay = "Mu24"; 
	ds->fColor = kBlack; 
	ds->fSymbol = 20; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3350; 
      }
      
      ds->fLcolor = ds->fColor; 
      ds->fFcolor = ds->fColor; 
      ds->fName   = sdecay; 
      ds->fFullName = sname; 
      fDS.insert(make_pair(sname, ds)); 
      
      
    } else {
      // -- MC
      pF = loadFile(sfile); 
      cout << "  " << sfile << ": " << pF << endl;
      
      dataset *ds = new dataset(); 
      ds->fSize = 1; 
      ds->fWidth = 2; 
 
      if (string::npos != stype.find("bsignal,lowpu")) {
        sname = "bSignalLowPU"; 
        sdecay = "bSignalLowPU"; 
	ds->fColor = kBlue-7; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      }

      if (string::npos != stype.find("csignal,lowpu")) {
        sname = "cSignalLowPU"; 
        sdecay = "cSignalLowPU"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      }


      if (string::npos != stype.find("bsignal,Mu8")) {
        sname = "bSignalMu8"; 
        sdecay = "bSignalMu8"; 
	ds->fColor = kBlue-7; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365; 
      }

      if (string::npos != stype.find("csignal,Mu8")) {
        sname = "cSignalMu8"; 
        sdecay = "cSignalMu8"; 
	ds->fColor = kGreen+3; 
	ds->fSymbol = 24; 
	ds->fF      = pF; 
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3356; 
      }


      ds->fLcolor = ds->fColor; 
      ds->fFcolor = ds->fColor; 
      ds->fName   = sdecay; 
      ds->fFullName = sname; 
      fDS.insert(make_pair(sname, ds)); 
      
    }
  }

  is.close();
  cout << "Summary: " << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    cout << "===> " << it->first << endl;
    cout << "       " << it->second->fName << endl;
    cout << "       " << it->second->fF->GetName() << endl;
    cout << "       " << it->first << ": " << it->second->fName << ", " << it->second->fF->GetName() << endl;
  }
}

  
