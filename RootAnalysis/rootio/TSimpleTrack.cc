#include "TSimpleTrack.hh"
#include <iostream>

ClassImp(TSimpleTrack)

using namespace std;

TSimpleTrack::TSimpleTrack(int index) { 
  clear();
}

void TSimpleTrack::clear() {
  fIndices = fBits = 0; 
  fPx = fPy = fPz = 0;
}


void TSimpleTrack::dump() {
  cout << Form("SimpleTrack idx = %d q = %+d, p = (%+6.4f, %+6.4f, %+6.3f), HP = %d, muid = %d, PV = %d, genIdx = %d", 
	       getIndex(), getCharge(), fPx, fPy, fPz, getHighPurity(), getMuonID(), getPvIndex(), getGenIndex())
       << endl;
}
  
