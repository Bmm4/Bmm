#include "plotRecoil.hh"

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

#include "common/dataset.hh"
#include "common/util.hh"

ClassImp(plotRecoil)

using namespace std;

// ----------------------------------------------------------------------
plotRecoil::plotRecoil(string dir, string files, string cuts, string setup, int year): plotClass(dir, files, cuts, setup, year) {
  loadFiles(files);

  if (setup == "") {
    fHistFileName = Form("%s/plotRecoil.root", dir.c_str());
  } else {
    fHistFileName = Form("%s/plotRecoil-%s.root", dir.c_str(), setup.c_str());
  }

  fTexFileName = fHistFileName;
  replaceAll(fTexFileName, ".root", ".tex");
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));

}


// ----------------------------------------------------------------------
plotRecoil::~plotRecoil() {

}


// ----------------------------------------------------------------------
void plotRecoil::makeAll(int bitmask) {

}




// ----------------------------------------------------------------------
void plotRecoil::recoil0() {


}



// ----------------------------------------------------------------------
void plotRecoil::loopFunction1() {

}

// ----------------------------------------------------------------------
void plotRecoil::bookHist(int mode) {

}



// ----------------------------------------------------------------------
void plotRecoil::candAnalysis() {
  fGoodCand = true;
}


// ----------------------------------------------------------------------
void plotRecoil::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
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
  cout << "==> plotRecoil::loopOverTree> loop over dataset " << fCds->fName << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotRecoil::*pF)(void);
  if (ifunc == 1) pF = &plotRecoil::loopFunction1;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotRecoil::setupTree(TTree *t) {

}


// ----------------------------------------------------------------------
void plotRecoil::setCuts(string cuts) {
  cout << "==> plotRecoil::setCuts: " << cuts << endl;

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
void plotRecoil::loadFiles(string afiles) {

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

      if (string::npos != stype.find("recoil0")) {
        sname = "data_bmm";
        sdecay = "bmm";
	ds->fColor = kBlack;
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


    } else {
      // -- MC
      pF = loadFile(sfile);
      cout << "  " << sfile << ": " << pF << endl;

      dataset *ds = new dataset();
      ds->fSize = 1;
      ds->fWidth = 2;

      if (string::npos != stype.find("recoil0")) {
        sname = "recoil0";
        sdecay = "recoil0";
	ds->fColor = kBlue-7;
	ds->fSymbol = 24;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }


      cout << "  inserting as " << sname << " and " << sdecay << endl;
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
