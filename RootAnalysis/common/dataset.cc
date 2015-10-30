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
  fFullName = "";
  fColor = fLcolor = fFcolor = fSymbol = fFillStyle = fSize = fWidth = -1; 

  fNclone = 0; 
}


// ----------------------------------------------------------------------
TH1D*  dataset::getHist(string name, bool clone) {
  if (!fF) return 0;
  TH1D *h(0); 
  if (clone) 
    h = (TH1D*)((TH1D*)(fF->Get(Form("%s", name.c_str()))))->Clone(Form("%s_%s_%d", fFullName.c_str(), name.c_str(), fNclone));
  else 
    h = (TH1D*)(fF->Get(Form("%s", name.c_str())));

  cout << fName << ": hist " << h->GetName() << " with entries = " << h->GetEntries() << " from " << fF->GetName() << endl;

  if (!h) return 0; 
  setHist(h); 
  if (fColor > -1) setHist(h, fColor, fSymbol, fSize, fWidth); 
  if (fFillStyle > -1) setFilledHist(h, fColor, fFcolor, fFillStyle, fWidth); 
  return h; 
}

// ----------------------------------------------------------------------
TH2D*  dataset::getHist2(string name, bool clone) {
  if (!fF) return 0;
  TH2D *h(0); 
  if (clone) 
    h = (TH2D*)((TH2D*)(fF->Get(Form("%s", name.c_str()))))->Clone(Form("%s_%s_%d", fFullName.c_str(), name.c_str(), fNclone));
  else 
    h = (TH2D*)(fF->Get(Form("%s", name.c_str())));

  cout << fName << ": hist " << h->GetName() << " with entries = " << h->GetEntries() << " from " << fF->GetName() << endl;

  if (!h) return 0; 
  setHist(h); 
  if (fColor > -1) setHist(h, fColor, fSymbol, fSize, fWidth); 
  if (fFillStyle > -1) setFilledHist(h, fColor, fFcolor, fFillStyle, fWidth); 
  return h; 
}

