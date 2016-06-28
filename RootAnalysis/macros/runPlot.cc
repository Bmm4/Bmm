#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>

#include "TROOT.h"
#include "TRint.h"
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TRandom.h"
#include "TUnixSystem.h"

#include "plotReducedOverlays.hh"
#include "plotWork.hh"

using namespace std;

// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  string progName  = argv[0];

  string dir("nada"), cuts("nada"), files("nada"), mode("nada"), setup("nada");
  int year(2016), plot(0);
  bool remove(false);

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-x"))  {remove= true;}
    if (!strcmp(argv[i], "-y"))  {year  = atoi(argv[++i]);}
    if (!strcmp(argv[i], "-p"))  {plot  = atoi(argv[++i]);}
    if (!strcmp(argv[i], "-m"))  {mode  = argv[++i];}
    if (!strcmp(argv[i], "-s"))  {setup = argv[++i];}
  }

  if (2016 == year) {
    if ("nada" == files) files = "plotReducedOverlays.2016.files";
    if ("nada" == cuts)  cuts  = "plotClass.2016.cuts";
    if ("nada" == dir)   dir   = "results";
    if ("nada" == setup) setup = "";
  }

  // -- overlays
  if (plot & 1) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotReducedOverlays a(dir, files, cuts, setup);
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }


  return 0;
}
