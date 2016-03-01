#include "TAnaVertex.hh"
#include <iostream>

ClassImp(TAnaVertex)

using namespace std;

// ----------------------------------------------------------------------
TAnaVertex::TAnaVertex() { 
  //  cout << ".X.X.X.X.X: In the ctor of TAnaVertex" << endl;
  fPoint.SetXYZ(0., 0., 0.);
}


// ----------------------------------------------------------------------
TAnaVertex::~TAnaVertex() {
  //  cout << ".Y.Y.Y.Y.Y: In the dtor of TAnaVertex" << endl;
}


// ----------------------------------------------------------------------
void TAnaVertex::clear() {
  //  cout << ".X.X.X.X.X: Clearing TAnaVertex" << endl;
  setInfo(-99.,-99,-99.,-99,-99); 
  fPoint.SetXYZ(-99.,-99.,-99.);
  fDxy = fDxyE = fD3d = fD3dE = -99.;
  fCovXX[0] = -99.;
  fCovXX[1] = -99.;
  fCovXX[2] = -99.;
  fCovXX[3] = -99.;
  fCovXX[4] = -99.;
  fCovXX[5] = -99.;
  fCovXX[6] = -99.;
  fCovXX[7] = -99.;
  fCovXX[8] = -99.; 
  //  cout << ".X.X.X.X.X: Clearing TAnaVertex END" << endl;
}

// ----------------------------------------------------------------------
void TAnaVertex::setCovXX(Double_t *xx) {
  fCovXX[0] = xx[0];
  fCovXX[1] = xx[1];
  fCovXX[2] = xx[2];
  fCovXX[3] = xx[3];
  fCovXX[4] = xx[4];
  fCovXX[5] = xx[5];
  fCovXX[6] = xx[6];
  fCovXX[7] = xx[7];
  fCovXX[8] = xx[8];
}

// ----------------------------------------------------------------------
void TAnaVertex::getCovXX(TMatrixD &bla) {
  if ( (bla.GetNcols() != 3) || (bla.GetNrows() != 3)) bla.ResizeTo(3,3);
  bla(0,0) = fCovXX[0];  bla(0,1) = fCovXX[1];  bla(0,2) = fCovXX[2];
  bla(1,0) = fCovXX[3];  bla(1,1) = fCovXX[4];  bla(1,2) = fCovXX[5];
  bla(2,0) = fCovXX[6];  bla(2,1) = fCovXX[7];  bla(2,2) = fCovXX[8];
}
 
// ----------------------------------------------------------------------
void TAnaVertex::dump() {
  char line[200];
  sprintf(line, "chi2/dof=%5.3f/%3.0f prob=%5.3f ntrk=%2d vtx=(%5.3f,%5.3f,%5.3f), flsxy = %5.3f",
	  fChi2, fNdof, fProb, fNtracks, fPoint.X(), fPoint.Y(), fPoint.Z(), fDxy/fDxyE);
  cout << line;
  cout << endl;
}

