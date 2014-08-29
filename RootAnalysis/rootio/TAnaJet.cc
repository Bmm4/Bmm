#include "TAnaJet.hh"
#include <iostream>

ClassImp(TAnaJet)

using namespace std;

TAnaJet::TAnaJet(int index) { 
  fIndex = index;
  fNtracks = 0; 

}

void TAnaJet::clear() {
  fIndex    = -1;
  
}


void TAnaJet::dump() {
  //   cout << " q = " << fQ
  //        << " p = " << fPlab.Mag() 
  //        << " f = " << fPlab.Phi() 
  //        << " t = " << fPlab.Theta();

  cout << Form("eT=%6.2f f=%+4.3f eta=%+4.3f ", fEt, fPlab.Phi(), fPlab.Eta());
 
  if (fIndex > -1) {
    cout << " idx = " << fIndex;
  }

  cout << " trks: (" << fNtracks << ")";
  for (int i = 0; i < fNtracks; ++i) cout << getTrack(i) << ", ";
    
  cout << endl;
}
  
void TAnaJet::addTrack(int index) {
  if (fNtracks == TANAJET_MAXTRK-1) {
    cout << "TAnaJet: Too many tracks in jet. Not adding track!" << endl;
    return;
  }

  //  cout << "TAnaVertex: Adding track index " << index << " at array position " << fNtracks << endl;
  fTracksIndex[fNtracks] = index;
  ++fNtracks;
}
