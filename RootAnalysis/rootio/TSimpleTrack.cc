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
  TVector3 p(fPx, fPy, fPz); 
  cout << Form("SimpleTrack idx = %d q = %+d, p/t = (%+6.4f, %+6.4f, %+6.3f), HP = %d, muid = %d, PV = %d, genIdx = %d", 
	       getIndex(), getCharge(), p.Pt(), p.Eta(), p.Phi(), getHighPurity(), getMuonID(), getPvIndex(), getGenIndex())
       << endl;
}
  
