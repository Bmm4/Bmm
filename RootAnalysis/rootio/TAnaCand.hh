#ifndef TANACAND
#define TANACAND

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "TAnaVertex.hh"

#include <vector>
#include <utility>

class TAnaCand: public TObject {

public:
  TAnaCand(int index = -99);
  virtual ~TAnaCand();
  void clear();
  void dump();


  // ----------------------------------------------------------------------
  int         fSig1, fSig2;  // first and last (second) index in Signal Track block
  int         fMom;          // mother in Cand block
  int         fDau1, fDau2;  // daughters in Cand block

  int         fType, fQ;
  int         fIndex;
  double      fMass;
  double      fMassE;
  double      fMassC;
  double      fMassCE;

  TAnaVertex  fVtx;

  TVector3    fPlab;
  double      fPtE, fPhiE, fEtaE;
  double      fAlpha;

  // -- point of closest approach and related quantities
  double      fMinDoca, fMaxDoca;
  TVector3    fPoca;

  // -- information on the associated primary vertex
  int         fPvIdx, fPv2Idx;
  double      fPvLip, fPvLipE, fPvTip, fPvTipE;
  double      fPv2Lip, fPv2LipE, fPv2Tip, fPv2TipE; // Second best PV for pile-up detection
  double      fPvIP3d, fPvIP3dE, fPv2IP3d, fPv2IP3dE; // True 3D impact parameter

  // -- proper lifetime and its error
  double      fTauxy, fTauxyE; // from 2d (r-phi) measurement
  double      fTau3d, fTau3dE; // from 3d measurement

  // -- Delta chi2 PV fitting
  double      fDeltaChi2;

  // information on tracks near this candidate
  std::vector<std::pair<int,std::pair<float,float> > > fNstTracks; // usage: (trackIx,doca(value,error))

  // -- reserve
  int            fInt1,    fInt2,    fInt3;
  double         fDouble1, fDouble2, fDouble3;

private:
  ClassDef(TAnaCand,1)

};


#endif
