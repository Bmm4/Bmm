#include "TH1.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"

//extern TROOT  *gROOT;
//extern TStyle *gStyle;
// extern TFile  *gFile;
// extern TDirectory  *gDirectory;

void  hpl(TH1 *h, const char *options="");
void  hpl(const char *h, const char *options="");
