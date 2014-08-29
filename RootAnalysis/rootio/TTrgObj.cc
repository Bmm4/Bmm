#include "TTrgObj.hh"
#include <iostream>

ClassImp(TTrgObj)

using namespace std;

TTrgObj::TTrgObj() { }

TTrgObj::TTrgObj(Int_t Option) { }


void TTrgObj::dump() {
  char line[200];
  sprintf(line, "%4d %+4d p=(%+9.3f,%+9.3f,%+9.3f)=(%+9.3f,%+9.3f,%+9.3f) %s", 
	  fNumber, fID, fP.X(), fP.Y(), fP.Z(), fP.Pt(), fP.Eta(), fP.Phi(),
	  fLabel.Data());
  cout << line << endl;

}

void TTrgObj::dump(ofstream &OUT) {
  char line[200];
  sprintf(line, "%4d %+4d p=(%+9.3f,%+9.3f,%+9.3f)=(%+9.3f,%+9.3f,%+9.3f) %s", 
	  fNumber, fID, fP.X(), fP.Y(), fP.Z(), fP.Pt(), fP.Eta(), fP.Phi(),
	  fLabel.Data());
  OUT << line << endl;

}
  
