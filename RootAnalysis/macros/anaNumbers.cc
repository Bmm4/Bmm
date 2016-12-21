#include "anaNumbers.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "TMath.h"

using namespace std;

// ----------------------------------------------------------------------
void number::calcEtot() {
  etot = TMath::Sqrt(estat*estat + esyst*esyst);
}


// ----------------------------------------------------------------------
void number::setErrors(double sta, double sys) {
  estat = sta;
  esyst = sys;
  calcEtot();
}


// ----------------------------------------------------------------------
void number::setErrors(double sta, double sys, double tot) {
  estat = sta;
  esyst = sys;
  etot  = tot;
}



// ----------------------------------------------------------------------
anaNumbers::anaNumbers(string name, int ichan) : fName(name), fChan(ichan)
					       , fMcYield(NWIN), fObsYield(NWIN), fFitYield(NWIN)
{
  clear();
};

// ----------------------------------------------------------------------
void anaNumbers::clear() {
  //  cout << "clearing numbers for name = " << fName << " chan = " << fChan << endl;
  for (int iwin = 0; iwin < NWIN; ++iwin) {
    fMcYield[iwin].clear(fName);
    fObsYield[iwin].clear(fName);
    fFitYield[iwin].clear(fName);
  }

  vector<string> names;
  names.push_back("lo");
  names.push_back("bd");
  names.push_back("bs");
  names.push_back("hi");
  vector<vector<number> *> p;
  p.push_back(&fMcYield);
  p.push_back(&fObsYield);
  p.push_back(&fFitYield);
  for (unsigned int iv = 0; iv < p.size(); ++iv) {
    for (unsigned int i = 0; i < names.size(); ++i) {
      p[iv]->at(i).name = names[i];
    }
  }
}
