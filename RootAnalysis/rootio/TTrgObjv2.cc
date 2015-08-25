#include "TTrgObjv2.hh"
#include <iostream>

ClassImp(TTrgObjv2)

using namespace std;

TTrgObjv2::TTrgObjv2() { }

TTrgObjv2::TTrgObjv2(Int_t Option) { }


void TTrgObjv2::dump() {
  char line[200];
  sprintf(line, "HLT= %s %+4d Filter= %s %s %+4d ", 
	  fHltPath.Data(), fHltIndex, fLabel.Data(), fType.Data(), fNumber); 
  cout << line << endl;

  int num = fIndex.size();
  vector<int>::iterator i2=fID.begin();
  vector<int>::iterator i1=fIndex.begin();
  vector<TLorentzVector>::iterator i3=fP.begin();
  for(int n=0;n<num;++n) {
    sprintf(line, "index= %+4d id= %+4d p=(%+9.3f,%+9.3f,%+9.3f)", 
	    *i1, *i2, i3->Pt(),i3->Eta(),i3->Phi());  
    cout << line << endl;
    i1++; i2++; i3++;
  }

}

void TTrgObjv2::dump(ofstream &OUT) {
  char line[200];

  sprintf(line, "HLT= %s %+4d Filter= %s %s %+4d ", 
	  fHltPath.Data(), fHltIndex, fLabel.Data(), fType.Data(), fNumber); 
  OUT << line << endl;

  int num = fIndex.size();
  vector<int>::iterator i2=fID.begin();
  vector<int>::iterator i1=fIndex.begin();
  vector<TLorentzVector>::iterator i3=fP.begin();
  for(int n=0;n<num;++n) {
    sprintf(line, "index= %+4d id= %+4d p=(%+9.3f,%+9.3f,%+9.3f)", 
	    *i1, *i2, i3->Pt(),i3->Eta(),i3->Phi());  
    OUT << line << endl;
    i1++; i2++; i3++;
  }




}
  
