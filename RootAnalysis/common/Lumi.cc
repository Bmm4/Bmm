#include "Lumi.hh"
#include "util.hh"

#include <sstream>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

// ----------------------------------------------------------------------
Lumi::Lumi(string fname, int verbose) {
  fVerbose = verbose;

  cout << "Lumi initializing from " << fname << endl;
  parse(fname);

  if (fVerbose > 0) print();
}


// ----------------------------------------------------------------------
Lumi::~Lumi() {
}


// ----------------------------------------------------------------------
void Lumi::parse(std::string fname) {
  vector<string> lFile;
  string line;
  ifstream file(fname.c_str());
  if (!file.is_open()) {
    cout << "ERROR: could not open " << fname << endl;
    return;
  }
  while (getline(file, line)) {
    lFile.push_back(line);
  }

  string::size_type p1;
  string::size_type p2;
  string::size_type p3;

  bool doNormTag(true);
  for (unsigned int i = 0; i < lFile.size(); ++i) {
    if (doNormTag && string::npos != lFile[i].find("Norm tag:")) {
      fNormtag  = lFile[i];
      doNormTag = false;
    }
    if (string::npos != lFile[i].find("#Summary")) break;
    vector<string> elems = split(lFile[i], '|');
    if (elems.size() > 4) {
      string run = elems[1].substr(0, elems[1].find(':'));
      // cout << lFile[i] << ", run = " << run << " lumi = " << elems[6] << endl;
      int irun = atoi(run.c_str());
      double dlumi = atof(elems[6].c_str()) * 1.e-6;
      if (irun > 0) fMapRunLumi.insert(make_pair(irun, dlumi));
    }
  }
}


// ----------------------------------------------------------------------
double Lumi::lumi(int run) {
  return fMapRunLumi[run];
}


// ----------------------------------------------------------------------
double Lumi::totalLumi(int run1, int run2) {

  map<int, double>::iterator it = fMapRunLumi.begin();
  double total = 0.;
  for (; it != fMapRunLumi.end(); ++it) {
    if (run1 < run2) {
      if ((it->first >= run1) && (it->first <= run2)) {
	total += it->second;
      }
    } else {
      total += it->second;
    }
  }
  return total;
}


// ----------------------------------------------------------------------
bool Lumi::contains(int run) {
  return (fMapRunLumi.count(run) > 0);
}


// ----------------------------------------------------------------------
void Lumi::print() {
  map<int, double>::iterator it = fMapRunLumi.begin();
  cout << "Lumi for normtag " << fNormtag << endl;
  for (; it != fMapRunLumi.end(); ++it) {
    cout << "run: " << it->first << " lumi: " << it->second << endl;
  }

}
