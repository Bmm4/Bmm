#include <iostream>
#include <fstream>

void lsDir(string oname) {
  ofstream OUT(oname);
  TKey *key(0);
  string sname;
  TIter next(gDirectory->GetListOfKeys());
  while ((key = (TKey*)next())) {
    sname = key->GetName();
    OUT << sname << endl;
  }
  OUT.close();
}
