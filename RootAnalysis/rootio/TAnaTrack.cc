#include "TAnaTrack.hh"
#include <iostream>

ClassImp(TAnaTrack)

using namespace std;

TAnaTrack::TAnaTrack(int index) { 
  fIndex = index;
}

void TAnaTrack::clear() {
  fMCID     = -99999;
  fIndex    = -1;
  fGenIndex = -1;
  fMuIndex  = -1;
  fMuID     = -1;
  fAlgorithm= -1; 

  fQ = -99; 

  fDof = fValidHits = -99; 
  fChi2 = -99.;
  fTrackQuality = -99; 
  fAlgorithm = -99;

  fdz = fdzE = fd0 = fd0E = -99.;
  fLip = fLipE = fTip = fTipE = -99.; 

  fMuID = 0; 
  fMuIndex = -99; 

  fDeDx = fDeDx2 = -99.; 
  fDeDxNmeas = fDeDxNsat = -99; 

  fInt1 = fInt2 = fInt3 -99; 
  fDouble1 = fDouble2 = fDouble3 = -99.;
  
}


void TAnaTrack::dump() {
  //   cout << " q = " << fQ
  //        << " p = " << fPlab.Mag() 
  //        << " f = " << fPlab.Phi() 
  //        << " t = " << fPlab.Theta();

  cout << Form("q=%+2d pT=%6.2f f=%+4.3f eta=%+4.3f a=%2d", 
	       fQ, fPlab.Perp(), fPlab.Phi(), fPlab.Eta(), fAlgorithm
	       );

  if (fMCID != -99999) {
    cout << Form(" mcid=%+6d", fMCID);
  }

  if (fGenIndex > -1) {
    cout << Form(" mcidx=%3d", fGenIndex);
  }

  if (fIndex > -1) {
    cout << Form(" idx=%3d", fIndex);
  }

  if (fMuIndex > -1) {
    cout << Form(" *** midx=%3d", fMuIndex);
  }

  if (fMuID > -1) {
    cout << Form(" mid=%3x", fMuID);
  }
    
  cout << endl;
}
  
