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
#include "plotResults.hh"
#include "plotStuff.hh"
#include "plotFake.hh"
#include "plotBDT.hh"
#include "plotTrigger.hh"

#include "tmva1.hh"

#include "umlLifetime.hh"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

using namespace std;

int TMVAClassification( TString myMethodList);

// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  string progName  = argv[0];

  string dir("nada"), cuts("nada"), files("nada"), plot("nada"), mode("nada"), setup("nada"), rootfilename("nada");
  int year(2016);
  bool remove(false);

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-d"))  {dir   = argv[++i];}
    if (!strcmp(argv[i], "-f"))  {files = argv[++i];} // for tmva1: offset
    if (!strcmp(argv[i], "-m"))  {mode  = argv[++i];} // for tmva1: BDT parameters
    if (!strcmp(argv[i], "-p"))  {plot  = argv[++i];}
    if (!strcmp(argv[i], "-r"))  {rootfilename  = argv[++i];}
    if (!strcmp(argv[i], "-s"))  {setup = argv[++i];} // for tmva1: vars
    if (!strcmp(argv[i], "-w"))  {mode  = argv[++i];}
    if (!strcmp(argv[i], "-x"))  {remove= true;}
    if (!strcmp(argv[i], "-y"))  {year  = atoi(argv[++i]);}
  }

  if (2017 == year) {
    if ("nada" == files) files = "plotResults.2017.files";
    if ("nada" == cuts)  cuts  = "baseCuts.2017.cuts";
    if ("nada" == dir)   dir   = "results";
    if ("nada" == setup) setup = "";
  }

  if (2016 == year) {
    if ("nada" == files) files = "plotResults.2016.files";
    if ("nada" == cuts)  cuts  = "baseCuts.2016.cuts";
    if ("nada" == dir)   dir   = "results";
    if ("nada" == setup) setup = "";
  }

  if (2012 == year) {
    if ("nada" == files) files = "plotResults.2012.files";
    if ("nada" == cuts)  cuts  = "baseCuts.2012.cuts";
    if ("nada" == dir)   dir   = "results";
    if ("nada" == setup) setup = "";
  }

  if (2011 == year) {
    if ("nada" == files) files = "plotResults.2011.files";
    if ("nada" == cuts)  cuts  = "baseCuts.2011.cuts";
    if ("nada" == dir)   dir   = "results";
    if ("nada" == setup) setup = "";
  }


  // -- run everything
  if ("nada" == plot) {
    cout << "The FULL show" << endl;
    {
      gROOT->Clear();  gROOT->DeleteAll();
      plotResults a(dir, files, cuts, setup);
      a.makeAll();
      return 0;
    }

    {
      gROOT->Clear();  gROOT->DeleteAll();
      plotReducedOverlays a(dir, files, cuts, setup);
      a.makeAll();
    }

    {
      gROOT->Clear();  gROOT->DeleteAll();
      plotStuff a(dir, files, cuts, setup);
      a.makeAll();
    }

  }



  // -- results
  if (string::npos != plot.find("results")) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotResults a(dir, files, cuts, setup);
    if (rootfilename != "nada") a.changeSetup(dir, rootfilename, setup);
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }

  // -- overlays
  if (string::npos != plot.find("overlays")) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotReducedOverlays a(dir, files, cuts, setup);
    if (rootfilename != "nada") a.changeSetup(dir, rootfilename, setup);
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }


  // -- stuff
  if (string::npos != plot.find("stuff")) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotStuff a(dir, files, cuts, setup);
    if (rootfilename != "nada") a.changeSetup(dir, rootfilename, setup);
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }

  // -- fake
  if (string::npos != plot.find("fake")) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotFake a(dir, files, cuts, setup);
    if (rootfilename != "nada") a.changeSetup(dir, rootfilename, setup);
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }

  // -- work
  if (string::npos != plot.find("work")) {
    gROOT->Clear();  gROOT->DeleteAll();
    plotWork a(dir, files, cuts, setup);
    if (rootfilename != "nada") a.changeSetup(dir, rootfilename, setup);
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }

  // -- TMVA training
  if (string::npos != plot.find("tmva1")) {
    gROOT->Clear();  gROOT->DeleteAll();
    tmva1 a(year, setup, mode);
    int ifiles(100);
    if ("nada" == files) {
      ifiles = 100;
    } else {
      ifiles = atoi(files.c_str());
    }
    int chan = ifiles%10;
    if (9 == chan) chan = -1;
    cout << "calling tmva1::makeAll(" << ifiles << ", \"\", " << chan << ")" << endl;
    a.makeAll(ifiles, "", chan);
  }


  // -- BDT
  if (string::npos != plot.find("bdt")) {
    cout << "files: " << files << " cuts: " << cuts << " setup: " << setup << endl;
    gROOT->Clear();  gROOT->DeleteAll();
    plotBDT a(dir, files, cuts, setup);
    if (rootfilename != "nada") a.changeSetup(dir, rootfilename, setup);
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }

  // -- trigger
  if (string::npos != plot.find("trigger")) {
    cout << "files: " << files << " cuts: " << cuts << " setup: " << setup << endl;
    gROOT->Clear();  gROOT->DeleteAll();
    plotBDT a(dir, files, cuts, setup);
    if (rootfilename != "nada") a.changeSetup(dir, rootfilename, setup);
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }


  // -- umlLifetime
  if (string::npos != plot.find("umllifetime")) {
    gROOT->Clear();  gROOT->DeleteAll();
    files = "plotResults.2016.files";
    cuts  = "baseCuts.2016.cuts";
    setup = "";
    umlLifetime a(dir, files, cuts, setup);
    if (rootfilename == "nada") {
      a.fHistFileName = Form("%s.root", mode.c_str());
    } else {
      a.changeSetup(dir, rootfilename, setup);
    }
    if (mode != "nada") {
      a.makeAll(mode);
    } else {
      a.makeAll();
    }
  }

  return 0;
}
