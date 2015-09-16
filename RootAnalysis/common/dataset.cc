#include <iostream>
#include <fstream>
#include <iomanip>

#include "dataset.hh"
#include "util.hh"

using namespace std;

// ----------------------------------------------------------------------
dataset::dataset() {
  fF = 0; 
  fXsec = fFilterEff = fBf = fLumi = fMass = fLambda = 0.; 
  fName = "";
  fColor = fLcolor = fFcolor = fSymbol = fFillStyle = fSize = fWidth = -1; 
}


// ----------------------------------------------------------------------
TH1D*  dataset::getHist(string name, bool clone) {
  if (!fF) return 0;
  TH1D *h(0); 
  if (clone) 
    h = (TH1D*)((TH1D*)(fF->Get(Form("%s", name.c_str()))))->Clone();
  else 
    h = (TH1D*)(fF->Get(Form("%s", name.c_str())));

  if (!h) return 0; 
  setHist(h); 
  if (fColor > -1) setHist(h, fColor, fSymbol, fSize, fWidth); 
  if (fFillStyle > -1) setFilledHist(h, fColor, fFcolor, fFillStyle, fWidth); 
  return h; 
}

