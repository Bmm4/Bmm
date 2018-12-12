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
void number::add2Errors(number &n) {
  estat = TMath::Sqrt(estat*estat + n.estat*n.estat);
  esyst = TMath::Sqrt(esyst*esyst + n.esyst*n.esyst);
  etot  = TMath::Sqrt(etot*etot + n.etot*n.etot);
}

// ----------------------------------------------------------------------
void number::scaleErrors(double sf) {
  estat = sf*estat;
  esyst = sf*esyst;
  etot  = sf*etot;
}



// ----------------------------------------------------------------------
anaNumbers::anaNumbers(string name, int ichan) : fName(name), fChan(ichan)
					       , fMcYield(NWIN), fObsYield(NWIN), fObsYieldAm(NWIN), fFitYield(NWIN), fFitYieldAm(NWIN), fAmYield(NWIN)
{
  clear();
};

// ----------------------------------------------------------------------
void anaNumbers::clear() {
  //  cout << "clearing numbers for name = " << fName << " chan = " << fChan << endl;
  for (int iwin = 0; iwin < NWIN; ++iwin) {
    fMcYield[iwin].clear(fName);
    fAmYield[iwin].clear(fName);
    fObsYield[iwin].clear(fName);
    fObsYieldAm[iwin].clear(fName);
    fFitYield[iwin].clear(fName);
    fFitYieldAm[iwin].clear(fName);
  }

  fSignalFit.clear(fName);
  fW8SignalFit.clear(fName);
  fScaledYield.clear(fName);
  fScaledYieldAm.clear(fName);

  vector<string> names;
  names.push_back("lo");
  names.push_back("bd");
  names.push_back("bs");
  names.push_back("hi");
  names.push_back("all");
  vector<vector<number> *> p;
  p.push_back(&fMcYield);
  p.push_back(&fObsYield);
  p.push_back(&fObsYieldAm);
  p.push_back(&fFitYield);
  p.push_back(&fFitYieldAm);
  for (unsigned int iv = 0; iv < p.size(); ++iv) {
    for (unsigned int i = 0; i < names.size(); ++i) {
      p[iv]->at(i).name = names[i];
    }
  }
}
