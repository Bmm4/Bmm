#ifndef TANAVERTEX
#define TANAVERTEX

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"

class TAnaVertex: public TObject {
public:
  TAnaVertex();
  virtual ~TAnaVertex();
  void clear();
  void dump();
  void setInfo(double chi2, double ndof, double prob, int isFake, int ntracks) {
    fChi2 = chi2;
    fNdof = ndof;
    fProb = prob; 
    fStatus = isFake;
    fNtracks = ntracks;
  }

  void setCovXX(Double_t *x);
  void getCovXX(TMatrixD& x);

  int getNtracks() {return fNtracks;}

  float fChi2;
  float fNdof;
  float fProb;
  int   fStatus;
  int   fNtracks; 

  // -- Distance to primary vertex
  float fDxy, fDxyE;
  float fD3d, fD3dE;
  
  TVector3 fPoint;

private:
  float fCovXX[9];
  ClassDef(TAnaVertex,1)
};

#endif
