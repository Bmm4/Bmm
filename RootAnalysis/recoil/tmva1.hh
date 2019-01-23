#ifndef TMVA1
#define TMVA1

#include <iostream>
#include <fstream>
#include <vector>
#include <map>


#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TLine.h"
#include "TArrow.h"
#include "TBox.h"

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1.h"

class tmva1: public TObject {

  public:

  tmva1(int year = 2016, std::string vars = "gmtauhad:gm2hadmax:gm2hadmin:gpttauhad:gatauhad:gperptauhad:gparatauhad:gpttaur0pos:gpttaur0neg",
	std::string pars = "NTrees=800");
  ~tmva1();

  void train(std::string oname = "TMVA-0", std::string filename = "testRecoil12.recoilReader.jpsix.root", std::string cut = "");
  void makeAll(int offset, std::string filename = "");
  void make(int offset, std::string filename);

  void setBDTParameters(std::string pars) {fBDTParameters = pars;}

  std::string fVariables, fBDTParameters;

  ClassDef(tmva1,1) //Testing tmva1
};


#endif
