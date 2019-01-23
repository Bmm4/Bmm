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
#include "TTimeStamp.h"

#include "tmva1.hh"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/TMVARegGui.h"

using namespace std;

int TMVAClassification( TString myMethodList);

// ----------------------------------------------------------------------
// -- Usage:
// --------
//  bin/runPlot -p tmva1 -f 200 -s gmtauhad:gm2hadmax:gm2hadmin:gmbr0pos:gmbr0neg -r testRecoil.recoilReader.jpsix.root
//
//
// ----------------------------------------------------------------------
int main(int argc, char *argv[]) {

  TTimeStamp ts0;
  cout << "start time: " << ts0.AsString("lc") << endl;

  string progName  = argv[0];

  string dir("nada"), cuts("nada"), files("nada"), plot("nada"), mode("nada"), setup("nada"), rootfilename("nada"), syear("0");
  int year(0);

  // -- command line arguments
  for (int i = 0; i < argc; i++){                             // tmva1:               trainingfiles:
    if (!strcmp(argv[i], "-c"))  {cuts  = argv[++i];}         //
    if (!strcmp(argv[i], "-d"))  {dir   = argv[++i];}         //                      data input filename
    if (!strcmp(argv[i], "-f"))  {files = argv[++i];}         // offset
    if (!strcmp(argv[i], "-m"))  {mode  = argv[++i];}         // BDT parameters
    if (!strcmp(argv[i], "-p"))  {plot  = argv[++i];}         //
    if (!strcmp(argv[i], "-r"))  {rootfilename  = argv[++i];} // input rootfilename   output rootfilename
    if (!strcmp(argv[i], "-s"))  {setup = argv[++i];}         // vars                 signal input filename
    if (!strcmp(argv[i], "-w"))  {mode  = argv[++i];}         // BDT parameters
    if (!strcmp(argv[i], "-y"))  {syear = argv[++i];}         //
  }


  // -- TMVA training
  if (string::npos != plot.find("tmva1")) {
    int ioffset(100);
    if ("nada" == files) {
      ioffset = 100;
    } else {
      ioffset = atoi(files.c_str());
    }

    gROOT->Clear();  gROOT->DeleteAll();
    tmva1 a(year, setup, mode);
    cout << "calling tmva1::makeAll(" << ioffset << ", " << rootfilename << ")" << endl;
    a.makeAll(ioffset, rootfilename);
  }

  // -- TMVA create training files
  // if (string::npos != plot.find("trainingfiles")) {
  //   gROOT->Clear();  gROOT->DeleteAll();
  //   tmva1 a(0, "", "");
  //   cout << "calling tmva1::createInputFile(\"" << rootfilename << "\", \""
  // 	 << setup << "\", \""
  // 	 << dir << "\")" << endl;
  //   a.createInputFile(rootfilename, setup, dir);

  // }



  TTimeStamp ts1;
  cout << "end time: " << ts1.AsString("lc") << ", this is the end, my friend." << endl;

  return 0;
}
