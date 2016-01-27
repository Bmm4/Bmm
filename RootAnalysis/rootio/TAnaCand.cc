#include "TAnaCand.hh"
#include <iostream>

ClassImp(TAnaCand)

using namespace std;

TAnaCand::TAnaCand(int index) {
  clear();
  fIndex = index;
}

TAnaCand::~TAnaCand() {
}


void TAnaCand::clear() {
  fMass  = -99.;
  fMom   = -99;
  fDau1  = -99;
  fDau2  = -99;
  fSig1  = -99;
  fSig2  = -99;
  fQ     = -99; 
  fIndex = -99; 
  fPvIdx = -99;
  fDeltaChi2 = -99;
  fPvLip = fPvLipE = fPvTip = fPvTipE = -99.;
  fPvLip2 = fPvLipE2 = -99.;
  fTauxy = fTauxyE = -99;
  fTau3d = fTau3dE = -99.;
  fMassE = -99.;

  fInt1 = fInt2 = fInt3 -99; 
  fDouble1 = fDouble2 = fDouble3 = -99.;

  fVtx.clear();
  fNstTracks.clear();
}


void TAnaCand::dump() {
  cout << Form("Cand: idx = %3d m = %5.3f, pT=%6.2f f=%+4.3f eta=%+4.3f ", fIndex, fMass, fPlab.Perp(), fPlab.Phi(), fPlab.Eta()) << endl;
  cout << "daughter cands: " << fDau1 << " .. " << fDau2  << " sig tracks: " << fSig1 << " .. " << fSig2 << endl;
  fVtx.dump();
}

  
