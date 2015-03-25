#include "TGenCand.hh"
#include <iostream>

ClassImp(TGenCand)

using namespace std;

TGenCand::TGenCand() { }

TGenCand::TGenCand(Int_t Option) { }

TGenCand::TGenCand(const  TGenCand &other) { 
  fID     = other.fID; 
  fNumber = other.fNumber;
  fStatus = other.fStatus;
  fMom1   = other.fMom1; 
  fMom2   = other.fMom2;
  fDau1   = other.fDau1; 
  fDau2   = other.fDau2;
  fTag    = other.fTag;
  fQ      = other.fQ;
  fP      = other.fP; 
  fMass   = other.fMass;
  fV      = other.fV;
  fTime   = other.fTime;
}


void TGenCand::dump(int printPt) {
  char line[200];
  if (1 == printPt) {
    sprintf(line, "%4d %+6d S%2d mom(%4d,%4d) dau(%5d,%5d) P=%8.3f(%+9.3f,%+9.3f,%+9.3f) v=(%+8.6f,%+8.6f,%+8.6f)", 
	    fNumber, fID, fStatus, fMom1, fMom2, fDau1, fDau2,
	    fP.Rho(), 
	    fP.Perp(), (fP.Perp()>0.01?fP.Eta():99.), fP.Phi(),
	    fV.X(), fV.Y(), fV.Z());
  } else {
    sprintf(line, "%4d %+6d S%2d mom(%4d,%4d) dau(%5d,%5d) p=%8.3f(%+9.3f,%+9.3f,%+9.3f) v=(%+8.6f,%+8.6f,%+8.6f)", 
	    fNumber, fID, fStatus, fMom1, fMom2, fDau1, fDau2,
	    fP.Rho(), 
	    fP.X(), fP.Y(), fP.Z(),
	    fV.X(), fV.Y(), fV.Z());

  }
  cout << line << endl;

}

void TGenCand::dump(ofstream &OUT) {
  char line[200];
  sprintf(line, "%4d %+6d S%2d mom(%4d,%4d) dau(%5d,%5d) p=%8.3f(%+9.3f,%+9.3f,%+9.3f) v=(%+8.6f,%+8.6f,%+8.6f)", 
	  fNumber, fID, fStatus, fMom1, fMom2, fDau1, fDau2,
	  fP.Rho(), 
	  fP.X(), fP.Y(), fP.Z(),
	  fV.X(), fV.Y(), fV.Z());

  OUT << line << endl;

}
  
