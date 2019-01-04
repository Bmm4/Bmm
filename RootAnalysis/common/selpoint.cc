#include <iostream>
#include "selpoint.hh"

using namespace std;

// ----------------------------------------------------------------------
selpoint::selpoint() {
  fSmallerThan = vector<pair<double *, double> >();
  fLargerThan = vector<pair<double *, double> >();

  for (int icat = 0; icat < 10; ++icat) {
    fCnt[icat] = 0;
  }

}

// ----------------------------------------------------------------------
selpoint::~selpoint() {
  fSmallerThan.clear();
  fLargerThan.clear();

}

// ----------------------------------------------------------------------
void selpoint::eval(int cat, double w8) {

  for (unsigned int i = 0; i < fSmallerThan.size(); ++i) {
    if ((*(fSmallerThan[i].first)) > fSmallerThan[i].second) {
      return;
    }
  }

  for (unsigned int i = 0; i < fLargerThan.size(); ++i) {
    if ((*(fLargerThan[i].first)) < fLargerThan[i].second) {
      return;
    }
  }

  if (cat > -1 && cat < 10) {
    fCnt[cat] += w8;
  } else {
    cout << "selpoint: invalid cat " << cat << endl;
  }

}
