#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TText.h"
#include "TH2.h"
#include "TGraph.h"
#include "TPaveStats.h"
#include "TMinuit.h"
#include "TVirtualFitter.h"

#include "tmva1.hh"
#include "setupReader.hh"
#include "common/initFunc.hh"

#include "TMVA/Config.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/tmvaglob.h"

//2011
//#define LUMISCALE 2.01e-4

//2012  12/7405 = 0.00162
//#define LUMISCALE 0.00162

#define MINKS 0.05

ClassImp(tmva1)

using namespace std;

// ----------------------------------------------------------------------
// --
// -- USAGE: a.makeAll(0, 1); > TMVA-0.log
// --
// ----------------------------------------------------------------------

tmva1::tmva1(int year, string vars, string pars) {

  fVariables     = vars;
  fBDTParameters = pars;

  cout << "tmva1 hello: setup for year = " << year << endl
       << " with variables:  " << vars << endl
       << " with parameters: " << fBDTParameters
       << endl;

  fVariables = vars;

  legg = 0;
  legge = 0;
  tl = new TLatex();
  tl->SetTextFont(42);
  tl->SetTextSize(0.03);
  tl->SetNDC(kTRUE);

  fYear = year;
  int nbins(40);
  fH1s = new TH1D("h1s", "signal", nbins, -1., 1.);
  fH1b = new TH1D("h1b", "background", nbins, -1., 1.);
  fH1r = new TH1D("h1r", "ratio", nbins, -1., 1.);

  fLumiScale = 1.86e-3; // 37/19845
  if (year == 2011) {
    fLumiScale = 3.1e-4; // 4.9/16000
  } else if (year == 2012) {
    fLumiScale = 2.8e-4; // 20/714000
  } else if (year == 2016) {
    fLumiScale = 1.86e-3; // 37/19845 I think this was wrong, 19845 refers to the bsmmMcComb lumi, not the official one!
    fLumiScale = 8.8e-4; // 37/42185 (combined bsmm und bdmm signal MC)
  }

  fApplyOn0 = false;
  fApplyOn1 = false;
  fApplyOn2 = false;
  fTrainAntiMuon = false;
  fChannel = 0;

}


// ----------------------------------------------------------------------
tmva1::~tmva1() {
  cout << "tmva1 good bye " << endl;
}


// ----------------------------------------------------------------------
void tmva1::makeAll(int offset, string filename, int chan) {
  fChannel = chan;

  if (filename == "") {
    if (2011 == fYear) {
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-0-2011.root";
    }
    if (2012 == fYear) {
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-0-2012.root";
    }
    if (2016 == fYear) {
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-0-2016.root"; // first shot
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-1-2016.root"; // preselection on m?iso, "All" = Chan0 + Chan1 also contained
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-3-2016.root"; // corrected m1iso and m2iso calculation that also covers rare backgrounds
                                                                  //                                                 Nevt(chan0Events0): 10225
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-4-2016.root"; // harder preselection: fls3d>12, added tos&&l1t.  Nevt(chan0Events0):  1422
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-5-2016.root"; // harder preselection: fls3d>7                    Nevt(chan0Events0):  5644
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-6-2016.root"; // harder preselection: fls3d>5                    Nevt(chan0Events0): 10139
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-7-2016.root"; // more signal MC
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-8-2016.root"; // fls3d>10
    }
    if (0 == fYear) {
      filename = "/scratch/ursl/bmm4/bdt/tmva-trees-0-0.root";    // harder preselection: fls3d>5 && tos && hlt1
    }
  }

  string logname = Form("TMVA-%d", offset);
  string oname("");
  oname = Form("TMVA-%d", offset);
  fKS.clear();
  make(offset, filename, 0);
  make(offset, filename, 1);
  make(offset, filename, 2);
  cout << "----------------------------------------------------------------------" << endl;
  cout << "KS probabilities: " << endl;
  bool good(true), badSSB(false), badBDT(false);
  double minKS(99.);
  for (unsigned int i = 0; i < fKS.size(); i += 2) {
    cout << "sb/bg = " << fKS[i] << "/" << fKS[i+1] << endl;
    if (fKS[i] < MINKS) good = false;
    if (fKS[i+1] < MINKS) good = false;
    if (fKS[i] < minKS)   minKS = fKS[i];
    if (fKS[i+1] < minKS) minKS = fKS[i+1];
  }
  double  avbdt =  (fMaxBdt[0]+fMaxBdt[1]+fMaxBdt[2])/3;
  if (TMath::Abs(avbdt - fMaxBdt[0]) > 0.1*avbdt) {
    good = false;
    badBDT = true;
  }
  if (TMath::Abs(avbdt - fMaxBdt[1]) > 0.1*avbdt) {
    good = false;
    badBDT = true;
  }
  if (TMath::Abs(avbdt - fMaxBdt[2]) > 0.1*avbdt) {
    good = false;
    badBDT = true;
  }
  double  avssb =  (fMaxSSB[0]+fMaxSSB[1]+fMaxSSB[2])/3;
  if (TMath::Abs(avssb - fMaxSSB[0]) > 0.1*avssb) {
    good = false;
    badSSB = true;
  }
  if (TMath::Abs(avssb - fMaxSSB[1]) > 0.1*avssb) {
    good = false;
    badSSB = true;
  }
  if (TMath::Abs(avssb - fMaxSSB[2]) > 0.1*avssb) {
    good = false;
    badSSB = true;
  }
  cout << "ssb: " << fMaxSSB[0] << "/" << fMaxSSB[1] << "/" << fMaxSSB[2] << " -> avssb = " << Form("%3.2f", avssb) << endl;
  cout << "bdt: " << fMaxBdt[0] << "/" << fMaxBdt[1] << "/" << fMaxBdt[2] << " -> avbdt = " << Form("%3.2f", avbdt) << endl;
  cout << "offset = " << offset << " is a " << (good? "good": "bad") << " BDT, minKS = " << minKS << ", ssb0 = " << fMaxSSB[0]
       << " spreads:  " << badSSB << " " << badBDT
       << endl;
  cout << "----------------------------------------------------------------------" << endl;

  double minssb(99.);
  if (fMaxSSB[0] < minssb)  minssb = fMaxSSB[0];
  if (fMaxSSB[1] < minssb)  minssb = fMaxSSB[1];
  if (fMaxSSB[2] < minssb)  minssb = fMaxSSB[2];
  ofstream OUT;
  OUT.open("/shome/ursl/abdt.log", ios::app);
  OUT << "offset = " << offset << "/" << (good? "good": "bad") << " BDT, minKS = " << minKS << ", avssb = " << avssb << ", minssb = " << minssb << ", avbdt = " << avbdt
      << "/" << fBDTParameters << "/" << fVariables << "/" << filename
      << endl;
  OUT.close();

  OUT.open("append-basecuts.txt");
  OUT << Form("bdtxml   TMVA-%d   TMVA-%d    TMVA-%d   TMVA-%d", offset, offset, offset, offset) << endl;
  OUT << Form("bdtcut         %3.2f          %3.2f          %3.2f          %3.2f", avbdt, avbdt, avbdt, avbdt) << endl;
  OUT.close();

}

// ----------------------------------------------------------------------
void tmva1::make(int offset, string filename, int evt) {

  if (0 == evt)  setApply0();
  if (1 == evt)  setApply1();
  if (2 == evt)  setApply2();

  string type;
  switch (evt) {
  case 0:
    type = "Events0";
    break;
  case 1:
    type = "Events1";
    break;
  case 2:
    type = "Events2";
    break;
  default:
    cout << "All hell break loose" << endl;
  }

  string oname = Form("TMVA-%d-%s", offset, type.c_str());
  cout << "======================================================================" << endl;
  cout << "==> tmva1(" << oname << ", " << filename << ") " << endl;
  cout << "======================================================================" << endl;
  train(oname, filename);
}


// ----------------------------------------------------------------------
TCanvas* tmva1::getC0() {
  TCanvas *c0 = (TCanvas*)gROOT->FindObject("c0");
  if (0 == c0) c0 = new TCanvas("c0","--c0--",2303,0,656,700);
  return c0;
}


// ----------------------------------------------------------------------
void tmva1::train(string oname, string filename, int nsg, int nbg, string cut) {
  // This loads the library
  TMVA::Tools::Instance();

  (TMVA::gConfig().GetVariablePlotting()).fNbins1D = 40;
  (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput = 40;

  int offset;
  string sint = oname.substr(oname.find("TMVA-")+5, oname.rfind("-") - oname.find("TMVA-")-5);
  offset = atoi(sint.c_str());
  cout << "offset= " << offset << " from ->" << sint << "<-" << endl;

  // -- Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName(Form("%s.root", oname.c_str()));
  TFile* outputFile = TFile::Open(outfileName, "RECREATE" );

  TH1D *hSetup = new TH1D("hSetup", "hSetup", 100, 0., 100.);
  int i(0);
  i =  1; hSetup->SetBinContent(i, fTrainAntiMuon?1:0); hSetup->GetXaxis()->SetBinLabel(i, "antimuon");
  i =  3; hSetup->SetBinContent(i, fRsigma); hSetup->GetXaxis()->SetBinLabel(i, "rsigma");
  i =  5; hSetup->SetBinContent(i, fApplyOn0?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn0");
  i =  6; hSetup->SetBinContent(i, fApplyOn1?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn1");
  i =  7; hSetup->SetBinContent(i, fApplyOn2?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn2");
  i = 14; hSetup->SetBinContent(i, offset); hSetup->GetXaxis()->SetBinLabel(i, "bdtname");

  cout << "----------------------------------------------------------------------" << endl;
  cout << "==> oname: " << oname << " antimuon: " << fTrainAntiMuon <<  endl;
  string optstring = "V:!Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
  optstring        = "V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification";
  cout << "==> Factory: " << optstring << endl;
  TMVA::Factory *factory = new TMVA::Factory(Form("%s", oname.c_str()), outputFile,  optstring.c_str());
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

  // -- parse string with all variables into a vector
  vector<string> vVar;
  string allvar = fVariables;
  string::size_type m0(allvar.size());
  cout << "--> allvar = " << allvar << endl;
  string var;
  while (string::npos != m0) {
    m0 = allvar.find(":");
    var = allvar.substr(0, m0);
    vVar.push_back(var);
    allvar = allvar.substr(m0+1);
    m0 = allvar.find(":");
  }
  var = allvar;
  vVar.push_back(var);

  for (unsigned int i = 0; i < vVar.size(); ++i) {
    cout << " addVariable:  " << vVar[i] << endl;
    //      if (string::npos != vVar[i].find("closetrk")) {
    //        factory->AddVariable(vVar[i].c_str(), 'I');
    //      } else {
    dataloader->AddVariable(vVar[i].c_str(), 'F');
    //      }
  }

  dataloader->AddSpectator("m",  "mass", "GeV", 'F' );

  TFile* inFile;
  TTree *applySg(0), *trainSg(0), *testSg(0), *applyBg(0), *trainBg(0), *testBg(0);

  inFile = TFile::Open(filename.c_str());
  cout << "opening " << filename << " inFile = " << inFile << endl;
  string sChannel("");
  if (fChannel > -1) {
    sChannel = Form("Chan%d", fChannel);
  } else {
    sChannel = Form("All");
  }

  if (fApplyOn0) {
    cout << "==============> Apply on events0, train on events1, test on events2" << endl;
    cout << Form("signal%sEvents0/events", sChannel.c_str()) << " "
	 << Form("sideband%sEvents0/events", sChannel.c_str()) << endl;
    applySg = (TTree*)inFile->Get(Form("signal%sEvents0/events", sChannel.c_str()));
    trainSg = (TTree*)inFile->Get(Form("signal%sEvents1/events", sChannel.c_str()));
    testSg  = (TTree*)inFile->Get(Form("signal%sEvents2/events", sChannel.c_str()));
    applyBg = (TTree*)inFile->Get(Form("sideband%sEvents0/events", sChannel.c_str()));
    trainBg = (TTree*)inFile->Get(Form("sideband%sEvents1/events", sChannel.c_str()));
    testBg  = (TTree*)inFile->Get(Form("sideband%sEvents2/events", sChannel.c_str()));
    cout << "==============> trainSg =  " << trainSg->GetDirectory()->GetName() << " entries: " << trainSg->GetEntries() << endl;
    cout << "==============> testSg  =  " << testSg->GetDirectory()->GetName()  << " entries: " << testSg->GetEntries() << endl;
    cout << "==============> applySg =  " << applySg->GetDirectory()->GetName()  << " entries: " << applySg->GetEntries() << endl;
    cout << "==============> trainBg =  " << trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
    cout << "==============> testBg  =  " << testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
    cout << "==============> applyBg =  " << applyBg->GetDirectory()->GetName()  << " entries: " << applyBg->GetEntries() << endl;
  }


  if (fApplyOn1) {
    cout << "==============> Apply on events1, train on events2, test on events0" << endl;
    cout << Form("signal%sEvents1/events", sChannel.c_str()) << " "
	 << Form("sideband%sEvents1/events", sChannel.c_str()) << endl;
    applySg = (TTree*)inFile->Get(Form("signal%sEvents1/events", sChannel.c_str()));
    trainSg = (TTree*)inFile->Get(Form("signal%sEvents2/events", sChannel.c_str()));
    testSg  = (TTree*)inFile->Get(Form("signal%sEvents0/events", sChannel.c_str()));
    applyBg = (TTree*)inFile->Get(Form("sideband%sEvents1/events", sChannel.c_str()));
    trainBg = (TTree*)inFile->Get(Form("sideband%sEvents2/events", sChannel.c_str()));
    testBg  = (TTree*)inFile->Get(Form("sideband%sEvents0/events", sChannel.c_str()));
    cout << "==============> trainSg =  " << trainSg->GetDirectory()->GetName() << " entries: " << trainSg->GetEntries() << endl;
    cout << "==============> testSg  =  " << testSg->GetDirectory()->GetName()  << " entries: " << testSg->GetEntries() << endl;
    cout << "==============> applySg =  " << applySg->GetDirectory()->GetName()  << " entries: " << applySg->GetEntries() << endl;
    cout << "==============> trainBg =  "<< trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
    cout << "==============> testBg  =  "<< testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
    cout << "==============> applyBg =  "<< applyBg->GetDirectory()->GetName()  << " entries: " << applyBg->GetEntries() << endl;
  }

  if (fApplyOn2) {
    cout << "==============> Apply on events2, train on events0, test on events1" << endl;
    cout << Form("signal%sEvents2/events", sChannel.c_str()) << " "
	 << Form("sideband%sEvents2/events", sChannel.c_str()) << endl;
    applySg = (TTree*)inFile->Get(Form("signal%sEvents2/events", sChannel.c_str()));
    trainSg = (TTree*)inFile->Get(Form("signal%sEvents0/events", sChannel.c_str()));
    testSg  = (TTree*)inFile->Get(Form("signal%sEvents1/events", sChannel.c_str()));
    applyBg = (TTree*)inFile->Get(Form("sideband%sEvents2/events", sChannel.c_str()));
    trainBg = (TTree*)inFile->Get(Form("sideband%sEvents0/events", sChannel.c_str()));
    testBg  = (TTree*)inFile->Get(Form("sideband%sEvents1/events", sChannel.c_str()));
    cout << "==============> trainSg =  " << trainSg->GetDirectory()->GetName() << " entries: " << trainSg->GetEntries() << endl;
    cout << "==============> testSg  =  " << testSg->GetDirectory()->GetName()  << " entries: " << testSg->GetEntries() << endl;
    cout << "==============> applySg =  " << applySg->GetDirectory()->GetName()  << " entries: " << applySg->GetEntries() << endl;
    cout << "==============> trainBg =  "<< trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
    cout << "==============> testBg  =  "<< testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
    cout << "==============> applyBg =  "<< applyBg->GetDirectory()->GetName()  << " entries: " << applyBg->GetEntries() << endl;
  }

  i = 30; hSetup->SetBinContent(i, applySg->GetEntries()); hSetup->GetXaxis()->SetBinLabel(i, "sgcnt");
  i = 31; hSetup->SetBinContent(i, applyBg->GetEntries()); hSetup->GetXaxis()->SetBinLabel(i, "bgcnt");
  writeOut(outputFile, hSetup);

  Double_t signalWeight      = 1.; //= LUMISCALE; // 0.000388
  Double_t rbackgroundWeight = 1.;
  Double_t cbackgroundWeight = 1.;
  Double_t tbackgroundWeight = cbackgroundWeight;

  cout << "--> signal weight:     " << signalWeight << endl;
  cout << "--> cbackground weight: " << cbackgroundWeight << endl;
  cout << "--> rbackground weight: " << rbackgroundWeight << endl;

  dataloader->AddSignalTree(trainSg,     signalWeight,      TMVA::Types::kTraining);
  dataloader->AddSignalTree(testSg,      signalWeight,      TMVA::Types::kTesting);
  dataloader->AddBackgroundTree(trainBg, cbackgroundWeight, TMVA::Types::kTraining);
  dataloader->AddBackgroundTree(testBg,  tbackgroundWeight, TMVA::Types::kTesting);

  int nSgTrain = trainSg->GetEntries();
  int nSgTest  = testSg->GetEntries();

  int nBgTrain = trainBg->GetEntries();
  int nBgTest  = testBg->GetEntries();

  if (nsg > -1) {
    nSgTrain = nsg;
    nSgTest = nsg;
  }

  if (nbg > -1) {
    nBgTrain = nbg;
    nBgTest = nbg;
  }

  if (nSgTrain > nBgTrain) {
    nSgTrain =  nBgTrain;
  } else {
    nBgTrain =  nSgTrain;
  }

  int seed = static_cast<int>(100*gRandom->Rndm());

  // optstring=Form("nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None:V");
  // optstring=Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=random:SplitSeed=%d:NormMode=None:V",
  // nSgTrain, nSgTest, nBgTrain, nBgTest, seed);

  optstring = Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=Block:NormMode=None:V",
		   nSgTrain, nSgTest, nBgTrain, nBgTest);
  cout << "==> PrepareTrainingAndTestTree: " << optstring << " cuts ->" << cut << "<-" << endl;
  TCut tCut(cut.c_str());
  dataloader->PrepareTrainingAndTestTree(tCut, optstring.c_str());

  optstring = "!H:V" + fBDTParameters;

  cout << "==> BookMethod: " << optstring << endl;
  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT", optstring);

  cout << "==> TrainAllMethods " << endl;
  factory->TrainAllMethods();

  // ---- Evaluate all MVAs using the set of test events
  cout << "==> TestAllMethods " << endl;
  factory->TestAllMethods();

  // ----- Evaluate and compare performance of all configured MVAs
  cout << "==> EvaluateAllMethods " << endl;
  factory->EvaluateAllMethods();

  // Save the output
  outputFile->Close();
  delete factory;
  gROOT->Clear();
  gROOT->DeleteAll();
  cout << "==> TMVAClassification is done!" << endl;

  // -- dump preselection information into XML files
  //  </GeneralInfo>
  TH1D *hpresel = (TH1D*)inFile->Get("hpresel");
  if (hpresel) {
    string xmlName = Form("dataset/weights/%s_BDT.weights.xml", oname.c_str());
    cout << "xmlName = " << xmlName  << endl;
    vector<string> lines;
    char  buffer[1000];
    ifstream is(xmlName);
    if (!is) {
      cout << "file ->" << xmlName << "<- not found, exit(1)" << endl;
      exit(1);
    }
    char input[1000];
    int icnt(0);
    while (is.getline(buffer, 1000, '\n')) {
      lines.push_back(string(buffer));
      ++icnt;
    }
    is.close();
    cout << "read a total of " << lines.size() << " lines, icnt = " << icnt << ", now adding preselection info" << endl;

    system(Form("/bin/mv %s %s.bac", xmlName.c_str(), xmlName.c_str()));

    ofstream os;
    os.open(xmlName, ios::out);

    for (unsigned int i = 0; i < lines.size(); ++i) {
      if (string::npos != lines[i].find("</GeneralInfo>")) {
	for (unsigned int ibin = 1; ibin <= hpresel->GetNbinsX(); ++ibin) {
	  string label = hpresel->GetXaxis()->GetBinLabel(ibin);
	  if (string::npos != label.find("cut:")) {
	    replaceAll(label, "cut:", "");
	    os <<  Form("    <Info name=\"Preselection:%s\" value=\"%f\"/>", label.c_str(), hpresel->GetBinContent(ibin)) << endl;
	  }
	  if (string::npos != label.find("file")) {
	    os <<  Form("    <Info name=\"Files:%s\" value=\"%f\"/>", label.c_str(), hpresel->GetBinContent(ibin)) << endl;
	  }
	}
      }
      os << lines[i] << endl;
    }
    os.close();
  }

  outputFile = TFile::Open(outfileName);
  TString hname = "dataset/Method_BDT/BDT/MVA_BDT";
  TH1 *sig = dynamic_cast<TH1*>(outputFile->Get(hname + "_S" ));
  TH1 *bgd = dynamic_cast<TH1*>(outputFile->Get(hname + "_B" ));
  cout << "==> Looking for Kolmogorov-Smirnov input: gDirectory = " << gDirectory->GetName()
       << " sig = " << sig << " bgd = " << bgd << endl;

  TH1D *hs = (TH1D*)sig->Clone("hs");
  hs->Scale(1./hs->Integral());
  TH1D *hb = (TH1D*)bgd->Clone("hb");
  hb->Scale(1./hb->Integral());

  sig = dynamic_cast<TH1*>(outputFile->Get(hname + "_Train_S" ));
  bgd = dynamic_cast<TH1*>(outputFile->Get(hname + "_Train_B" ));
  TH1D *hS = (TH1D*)sig->Clone("hS");
  hS->Scale(1./hS->Integral());
  TH1D *hB = (TH1D*)bgd->Clone("hB");
  hB->Scale(1./hB->Integral());

  double kolS = hs->KolmogorovTest(hS);
  double kolB = hb->KolmogorovTest(hB);

  fKS.push_back(kolS);
  fKS.push_back(kolB);
  cout << "KS-probability /" << gDirectory->GetName() << "/ ks-sg = " << kolS << " ks-bg = " << kolB << endl;

  outfileName.ReplaceAll(".root", ".pdf");
  TTree *t = dynamic_cast<TTree*>(gDirectory->Get("dataset/TestTree"));
  int id;
  float bdt;
  t->SetBranchAddress("classID", &id);
  t->SetBranchAddress("BDT", &bdt);
  fH1s->Reset();
  fH1b->Reset();
  fH1r->Reset();
  for (int jentry = 0; jentry < t->GetEntries(); jentry++) {
    t->GetEntry(jentry);
    if (0 == id) {
      fH1s->Fill(bdt);
    } else {
      fH1b->Fill(bdt);
    }
  }
  fH1s->Scale(fLumiScale);

  fH1s->Draw();
  gPad->SaveAs(Form("h1s-%s", outfileName.Data()));
  fH1b->Draw();
  gPad->SaveAs(Form("h1b-%s", outfileName.Data()));

  double rMax(-1.), rBdt(99);
  int nbins(fH1s->GetNbinsX()+1);
  for (int ibin = 1; ibin < nbins; ++ibin) {
    double s = fH1s->Integral(ibin, nbins);
    double b = fH1b->Integral(ibin, nbins);
    double r = 0.;
    if (s+b > 0.) {
      r = s/TMath::Sqrt(s+b);
      fH1r->SetBinContent(ibin, r);
      if (r > rMax) {
	rMax = r;
	rBdt = fH1s->GetBinCenter(ibin);
      }
    }
    cout << "s(" << ibin << "," << nbins << ") = " << s << " b = " << b << " r = " << r << " (bin center = " << fH1s->GetBinCenter(ibin) << ")" << endl;
  }
  fMaxSSB.push_back(rMax);
  fMaxBdt.push_back(rBdt);
  cout << "hello" << endl;
  // h1r->Write();
  cout << "hello2" << endl;
  outputFile->Close();
  cout << "hello3" << endl;
  fH1r->Draw();
  gPad->SaveAs(Form("h1r-%s", outfileName.Data()));

  if ((kolS < MINKS) || (kolB < MINKS)) {
    cout << "bad BDT" << endl;
    exit(0);
  }

}



// ----------------------------------------------------------------------
void tmva1::createInputFile(string filename, string sfile, string dfile, int randomSeed) {

  fInputFiles.sname = sfile;
  fInputFiles.dname = dfile;

  cout << "signal: " << fInputFiles.sname << endl;
  cout << "data:   " << fInputFiles.dname << endl;
  TFile *sinput = TFile::Open(fInputFiles.sname.c_str());
  TFile *dinput = TFile::Open(fInputFiles.dname.c_str());

  fPresel.setCut("ALPHAMAX", 0.2);
  fPresel.setCut("PVIPSMAX", 4.);
  fPresel.setCut("PVIPMAX", 0.02);
  TCut sgcut = fPresel.preselection().c_str();

  cout << "new: " << endl;
  cout << sgcut << endl;

  TCut masscut = "m>4.9&&m<5.9";
  TCut massbg  = "!(5.2<m&&m<5.45)";
  TCut muonid  = "gmugmid && hlt1 && l1t && tos";

  cout << "==> signal input file:     " << sinput->GetName() << std::endl;
  cout << "==> background input file: " << dinput->GetName() << std::endl;

  TTree *signal      = (TTree*)sinput->Get("candAnaMuMu/events");
  TTree *cbackground = (TTree*)dinput->Get("candAnaMuMu/events");

  TFile *outFile = TFile::Open(filename.c_str(),"RECREATE");

  // -- channel selection/definition
  string chanDef[] = {"chan == 0", "chan == 1"};

  int nchan = 2;
  string sdir, type;
  TTree *copyTree(0);
  TCut copyCuts;
  TCut chanCut, typeCut;
  if (randomSeed > -1) gRandom->SetSeed(randomSeed);
  for (int j = 0; j < 3; ++j) {
    if (0 == j) {
      type = "Events0";
      typeCut = "TMath::Abs(evt%3)==0";
      if (randomSeed > -1) typeCut = "3*rndm%3==0";
    } else if (1 == j) {
      type = "Events1";
      typeCut = "TMath::Abs(evt%3)==1";
      if (randomSeed > -1) typeCut = "3*rndm%3==1";
    } else if (2 == j) {
      type = "Events2";
      typeCut = "TMath::Abs(evt%3)==2";
      if (randomSeed > -1) typeCut = "3*rndm%3==2";
    }

    if (0 /*no statistics*/) {
      for (int i = 0; i < nchan; ++i) {
	// -- signal
	sdir = Form("signalChan%d%s", i, type.c_str());
	chanCut = chanDef[i].c_str();
	outFile->mkdir(sdir.c_str());
	outFile->cd(sdir.c_str());
	copyCuts = sgcut + muonid + chanCut + typeCut;
	cout << "sg copyCuts: " << copyCuts << endl;
	copyTree = signal->CopyTree(copyCuts);
	cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;

	// -- background
	sdir = Form("sidebandChan%d%s", i, type.c_str());
	chanCut = chanDef[i].c_str();
	outFile->mkdir(sdir.c_str());
	outFile->cd(sdir.c_str());
	copyCuts = sgcut + muonid + massbg + masscut + chanCut + typeCut;
	cout << "bg copyCuts: " << copyCuts << endl;
	copyTree = cbackground->CopyTree(copyCuts);
	cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;
      }
    }

    // -- combined version:
    chanCut = "(0 == chan) || (1 == chan)";
    // -- signal
    sdir = Form("signalAll%s", type.c_str());
    outFile->mkdir(sdir.c_str());
    outFile->cd(sdir.c_str());
    copyCuts = sgcut + muonid + chanCut + typeCut;
    cout << "sg copyCuts: " << copyCuts << endl;
    copyTree = signal->CopyTree(copyCuts);
    cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;

    // -- background
    sdir = Form("sidebandAll%s", type.c_str());
    outFile->mkdir(sdir.c_str());
    outFile->cd(sdir.c_str());
    copyCuts = sgcut + muonid + massbg + masscut + chanCut + typeCut;
    cout << "bg copyCuts: " << copyCuts << endl;
    copyTree = cbackground->CopyTree(copyCuts);
    cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;


  }

  TH1D *h1 = fPresel.getPreselectionNumbers();
  h1->GetXaxis()->SetBinLabel(98, Form("file-sg:%s", fInputFiles.sname.c_str()));
  h1->GetXaxis()->SetBinLabel(99, Form("file-bg:%s", fInputFiles.dname.c_str()));

  h1->SetDirectory(outFile);
  h1->Write();

  outFile->Write();
  outFile->Close();

  sinput->Close();
  dinput->Close();
}


// ----------------------------------------------------------------------
void tmva1::writeOut(TFile *f, TH1 *h) {
  TDirectory *pD = gDirectory;
  f->cd();
  h->SetDirectory(f);
  h->Write();
  pD->cd();
}
