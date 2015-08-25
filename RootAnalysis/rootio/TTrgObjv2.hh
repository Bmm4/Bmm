#ifndef TTRGOBJV2
#define TTRGOBJV2


#include <fstream>
#include <vector>

#include "TObject.h"
#include "TLorentzVector.h"

class TTrgObjv2: public TObject {

public:

  TTrgObjv2();
  TTrgObjv2(int Option);
  ~TTrgObjv2() { };
  void     clear() {fNumber = -1; fHltIndex = -1; fP.clear(); fID.clear();
    fIndex.clear(); fHltPath=""; fLabel=""; fType="";}

  // ----------------------------------------------------------------------
  void dump();
  void dump(std::ofstream &);

  // ----------------------------------------------------------------------
  int            fNumber, fHltIndex;
  std::vector<TLorentzVector> fP;
  TString        fLabel, fType, fHltPath;
  std::vector<int>    fID, fIndex;

private:

  ClassDef(TTrgObjv2,1)

};

#endif
