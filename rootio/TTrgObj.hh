#ifndef TTRGOBJ
#define TTRGOBJ


#include <fstream>

#include "TObject.h"
#include "TLorentzVector.h"

class TTrgObj: public TObject {

public:

  TTrgObj();
  TTrgObj(int Option);
  ~TTrgObj() { };
  void     clear() {fNumber = -1; fID = -123;}

  // ----------------------------------------------------------------------
  void dump();
  void dump(ofstream &);

  // ----------------------------------------------------------------------
  int            fNumber, fID;
  TLorentzVector fP;
  TString        fLabel;

private:

  ClassDef(TTrgObj,1)

};

#endif
