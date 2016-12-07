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

  fPvIdx = fPv2Idx = -99;
  fDeltaChi2 = -99;
  fPvLip = fPvLipE = fPvTip = fPvTipE = -99.;
  fPv2Lip = fPv2LipE = fPv2Tip = fPv2TipE = -99.;
  fPvIP3d = fPvIP3dE = fPv2IP3d = fPv2IP3dE = -99.;
  fTauxy = fTauxyE = -99;
  fTau3d = fTau3dE = -99.;
  fMassE = -99.;

  fInt1 = fInt2 = fInt3 -99;
  fDouble1 = fDouble2 = fDouble3 = -99.;

  fVtx.clear();
  fNstTracks.clear();
}


void TAnaCand::dump() {
  cout << Form("Cand: idx=%3d type=%d m=%5.3f cm=%5.3f, pT=%6.2f f=%+4.3f eta=%+4.3f maxdoca=%4.3f pvips=%4.3f",
	       fIndex, fType, fMass, fMassC, fPlab.Perp(), fPlab.Phi(), fPlab.Eta(), fMaxDoca, (fPvIP3dE>0.?fPvIP3d/fPvIP3dE:-1.))
       << endl;

  cout << "mother cand: " << fMom << " daughter cands: " << fDau1 << " .. " << fDau2  << " sig tracks: " << fSig1 << " .. " << fSig2 << endl;
  fVtx.dump();
}
