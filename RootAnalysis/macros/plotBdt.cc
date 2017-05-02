#include "plotBdt.hh"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>


#include "TROOT.h"
#include "TStyle.h"
#include "TKey.h"
#include "TMath.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TPaveStats.h"

#include "TMVA/tmvaglob.h"
#include "TMVA/Config.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Config.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"

#include "preselection.hh"
#include "common/dataset.hh"
#include "common/util.hh"
#include "common/Lumi.hh"

ClassImp(plotBdt)

using namespace std;

// ----------------------------------------------------------------------
plotBdt::plotBdt(string dir, string files, string cuts, string setup): plotClass(dir, files, cuts, setup) {
  plotClass::loadFiles(files);
  plotBdt::loadFiles(files);

  changeSetup(dir, "plotBdt", setup);
  init();

  // -- initialize cuts
  string cutfile = Form("%s/%s", dir.c_str(), cuts.c_str());
  cout << "===> Reading cuts from " << cutfile << endl;
  readCuts(cutfile);
  fNchan = fCuts.size();

  printCuts(cout);

  fChan = 0;

  fLumiScale = 37.5/10879.6;
  if (2012 == fYear) fLumiScale = 20./10000.;

  fVariables = "pt:eta:alpha:fls3d:maxdoca:pvip:pvips:iso:m1iso:m2iso";

}


// ----------------------------------------------------------------------
plotBdt::~plotBdt() {

}

// ----------------------------------------------------------------------
void plotBdt::init() {
  fTEX.close();
  cout << Form("/bin/rm -f %s", fTexFileName.c_str()) << endl;
  system(Form("/bin/rm -f %s", fTexFileName.c_str()));
  cout << Form("open for TeX output: %s", fTexFileName.c_str()) << endl;
  fTEX.open(fTexFileName.c_str(), ios::app);

}


// ----------------------------------------------------------------------
void plotBdt::makeAll(string what) {
  string filename("empty.root");
  if (string::npos != what.find("create")) {
    filename = Form("tmva-trees-0-%d.root", fYear);
    createInputFiles(filename);
  } else {
    if (2011 == fYear) {
      filename = "/scratch/ursl/bdt/tmva-trees-0-2011.root";
    }
    if (2012 == fYear) {
      filename = "/scratch/ursl/bdt/tmva-trees-0-2012.root";
    }
    if (2016 == fYear) {
      filename = "/scratch/ursl/bdt/tmva-trees-0-2016.root";
    }
  }

  int offset(0), clean(0);
  cout << "make(" << offset << ", " << filename << " ..." << endl;
  make(offset, filename, 0, clean);
  make(offset, filename, 1, clean);
  make(offset, filename, 2, clean);

  return;
  string oname = Form("TMVA-%d", offset);
  cout << "-->apply(" << oname.c_str() << ")" << endl;
  apply(oname.c_str());
  cout << "-->analyze(\"" << oname.c_str() << "\")" << endl;
  analyze(oname.c_str());

  cout << "-->mvas(...)" << endl;
  string sEvents = oname + "-Events0";
  mvas(sEvents.c_str());
  sEvents = oname + "-Events1";
  mvas(sEvents.c_str());
  sEvents = oname + "-Events2";
  mvas(sEvents.c_str());

  if (clean) {
    cout << "-->cleanup(...)" << endl;
    cleanup(oname);
  }

}




// ----------------------------------------------------------------------
void plotBdt::bookHist(string dsname) {

}



// ----------------------------------------------------------------------
void plotBdt::loopFunction1() {

}


// ----------------------------------------------------------------------
void plotBdt::loopFunction2() {

}



// ----------------------------------------------------------------------
void plotBdt::loopFunction3() {

}


// ----------------------------------------------------------------------
void plotBdt::loopFunction4() {

}

// ----------------------------------------------------------------------
void plotBdt::loopOverTree(TTree *t, int ifunc, int nevts, int nstart) {
  int nentries = Int_t(t->GetEntries());
  int nbegin(0), nend(nentries);
  if (nevts > 0 && nentries > nevts) {
    nentries = nevts;
    nbegin = 0;
    nend = nevts;
  }
  if (nevts > 0 && nstart > 0) {
    nentries = nstart + nevts;
    nbegin = nstart;
    if (nstart + nevts < t->GetEntries()) {
      nend = nstart + nevts;
    } else {
      nend = t->GetEntries();
    }
  }

  nentries = nend - nstart;

  int step(1000000);
  if (nentries < 5000000)  step = 500000;
  if (nentries < 1000000)  step = 100000;
  if (nentries < 100000)   step = 10000;
  if (nentries < 10000)    step = 1000;
  if (nentries < 1000)     step = 100;
  step = 500000;
  cout << "==> plotBdt::loopOverTree> loop over dataset " << (fCds?fCds->fName:"undefined") << " in file "
       << t->GetDirectory()->GetName()
       << " with " << nentries << " entries"
       << " nbegin = " << nbegin << " nend = " << nend
       << endl;

  // -- setup loopfunction through pointer to member functions
  //    (this is the reason why this function is NOT in plotClass!)
  void (plotBdt::*pF)(void);
  if (ifunc == 1) pF = &plotBdt::loopFunction1;
  if (ifunc == 2) pF = &plotBdt::loopFunction2;
  if (ifunc == 3) pF = &plotBdt::loopFunction3;
  if (ifunc == 4) pF = &plotBdt::loopFunction4;

  // -- the real loop starts here
  for (int jentry = nbegin; jentry < nend; jentry++) {
    t->GetEntry(jentry);
    if (jentry%step == 0) cout << Form(" .. evt = %d", jentry) << endl;

    candAnalysis();
    (this->*pF)();
  }

}


// ----------------------------------------------------------------------
void plotBdt::make(int offset, string filename, int evt, int clean) {

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

  cout << "-->train(...) with oname = " << oname << " and filename = " << filename << endl;
  train(oname, filename);
}


// ----------------------------------------------------------------------
void plotBdt::createInputFiles(string filename, int randomSeed) {
  TFile *sinput = fDS["bsmmMcComb"]->getFile();
  TFile *dinput = fDS["bmmData"]->getFile();

  TCut sgcut = preselection().c_str();

  cout << "new: " << endl;
  cout << sgcut << endl;

  TCut masscut = "m>4.9&&m<5.9";
  TCut massbg  = "!(5.2<m&&m<5.45)";

  cout << "==> signal input file:     " << sinput->GetName() << endl;
  cout << "==> background input file: " << dinput->GetName() << endl;

  TTree *signal      = (TTree*)sinput->Get("candAnaMuMu/events");
  TTree *cbackground = (TTree*)dinput->Get("candAnaMuMu/events");

  TFile *outFile = TFile::Open(filename.c_str(),"RECREATE");

  // -- channel selection/definition
  float etaCut(1.6);
  if (fYear == 2012) etaCut = 2.1;
  string chanDef[] = {"TMath::Abs(m1eta) < 0.7 && TMath::Abs(m2eta) < 0.7",
                      Form("(TMath::Abs(m1eta)>0.7 || TMath::Abs(m2eta)>0.7) && TMath::Abs(m1eta)<%3.1f && TMath::Abs(m2eta)<%3.1f",
			   etaCut, etaCut)
  };

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

    for (int i = 0; i < nchan; ++i) {
      // -- signal
      sdir = Form("signalChan%d%s", i, type.c_str());
      chanCut = chanDef[i].c_str();
      outFile->mkdir(sdir.c_str());
      outFile->cd(sdir.c_str());
      copyCuts = sgcut + chanCut + typeCut;
      cout << "sg copyCuts: " << copyCuts << endl;
      copyTree = signal->CopyTree(copyCuts);
      cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;

      // -- background
      sdir = Form("sidebandChan%d%s", i, type.c_str());
      chanCut = chanDef[i].c_str();
      outFile->mkdir(sdir.c_str());
      outFile->cd(sdir.c_str());
      copyCuts = sgcut + massbg + masscut + chanCut + typeCut;
      cout << "bg copyCuts: " << copyCuts << endl;
      copyTree = cbackground->CopyTree(copyCuts);
      cout << "--> " << copyTree->GetEntries() << " events in tree" << endl;
    }

  }

  outFile->Write();
  outFile->Close();

  sinput->Close();
  dinput->Close();
}


// ----------------------------------------------------------------------
void plotBdt::train(string oname, string infilename, int nsg, int nbg) {
   // This loads the library
   TMVA::Tools::Instance();

   // (TMVA::gConfig().GetVariablePlotting()).fNbins1D = 40;
   // (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput = 40;

   // -- Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName(Form("%s.root", oname.c_str()));
   TFile* outputFile = TFile::Open(outfileName, "RECREATE" );

   TH1D *hSetup = new TH1D("hSetup", "hSetup", 100, 0., 100.);
   int i(0);
   i =  1; hSetup->SetBinContent(i, fTrainAntiMuon?1:0); hSetup->GetXaxis()->SetBinLabel(i, "antimuon");
   i =  5; hSetup->SetBinContent(i, fApplyOn0?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn0");
   i =  6; hSetup->SetBinContent(i, fApplyOn1?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn1");
   i =  7; hSetup->SetBinContent(i, fApplyOn2?1:0); hSetup->GetXaxis()->SetBinLabel(i, "applyOn2");
   i = 10; hSetup->SetBinContent(i, fBdtSetup.NTrees); hSetup->GetXaxis()->SetBinLabel(i, "NTrees");
   i = 11; hSetup->SetBinContent(i, fBdtSetup.nEventsMin); hSetup->GetXaxis()->SetBinLabel(i, "nEventsMin");
   i = 12; hSetup->SetBinContent(i, fBdtSetup.nCuts); hSetup->GetXaxis()->SetBinLabel(i, "nCuts");
   i = 13; hSetup->SetBinContent(i, fBdtSetup.AdaBoostBeta); hSetup->GetXaxis()->SetBinLabel(i, "AdaBoostBeta");

   i = 20; hSetup->SetBinContent(i, fBdtSetup.MaxDepth); hSetup->GetXaxis()->SetBinLabel(i, "MaxDepth");
   i = 21; hSetup->SetBinContent(i, fBdtSetup.NNodesMax); hSetup->GetXaxis()->SetBinLabel(i, "NNodesMax");

   cout << "----------------------------------------------------------------------" << endl;

   string optstring = "bla"; //V:!Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
   optstring        = "V:!Silent:Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification";
   cout << "==> oname: " << oname << " outputFile: " << outputFile << " antimuon: " << fTrainAntiMuon <<  endl;
   cout << "==> Factory: " << optstring << endl;
   TMVA::Factory *factory = new TMVA::Factory(oname, outputFile,  optstring.c_str());

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
     cout << "  adding to factory: " << vVar[i] << endl;
     factory->AddVariable(vVar[i].c_str(), 'F');
   }

   factory->AddSpectator("m",  "mass", "GeV", 'F' );

   TTree *applySg(0), *trainSg(0), *testSg(0), *applyBg(0), *trainBg(0), *testBg(0);

   TFile *inFile = TFile::Open(infilename.c_str());
   if (fApplyOn0) {
     cout << "==============> Apply on events0, train on events1, test on events2" << endl;
     applySg = (TTree*)inFile->Get(Form("signalChan%dEvents0/events", fChannel));
     trainSg = (TTree*)inFile->Get(Form("signalChan%dEvents1/events", fChannel));
     testSg  = (TTree*)inFile->Get(Form("signalChan%dEvents2/events", fChannel));
     applyBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents0/events", fChannel));
     trainBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents1/events", fChannel));
     testBg  = (TTree*)inFile->Get(Form("sidebandChan%dEvents2/events", fChannel));
     cout << "==============> trainBg =  " << trainBg << ": " << trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
     cout << "==============> testBg  =  " << testBg  << ": " << testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
     cout << "==============> applyBg =  " << applyBg << ": " << applyBg->GetDirectory()->GetName()  << " entries: " << applyBg->GetEntries() << endl;
   }


   if (fApplyOn1) {
     cout << "==============> Apply on events1, train on events2, test on events0" << endl;
     applySg = (TTree*)inFile->Get(Form("signalChan%dEvents1/events", fChannel));
     trainSg = (TTree*)inFile->Get(Form("signalChan%dEvents2/events", fChannel));
     testSg  = (TTree*)inFile->Get(Form("signalChan%dEvents0/events", fChannel));
     applyBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents1/events", fChannel));
     trainBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents2/events", fChannel));
     testBg  = (TTree*)inFile->Get(Form("sidebandChan%dEvents0/events", fChannel));
     cout << "==============> trainBg =  "<< trainBg->GetDirectory()->GetName() << " entries: " << trainBg->GetEntries() << endl;
     cout << "==============> testBg  =  "<< testBg->GetDirectory()->GetName()  << " entries: " << testBg->GetEntries() << endl;
     cout << "==============> applyBg =  "<< applyBg->GetDirectory()->GetName()  << " entries: " << applyBg->GetEntries() << endl;
   }

   if (fApplyOn2) {
     cout << "==============> Apply on events2, train on events0, test on events1" << endl;
     applySg = (TTree*)inFile->Get(Form("signalChan%dEvents2/events", fChannel));
     trainSg = (TTree*)inFile->Get(Form("signalChan%dEvents0/events", fChannel));
     testSg  = (TTree*)inFile->Get(Form("signalChan%dEvents1/events", fChannel));
     applyBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents2/events", fChannel));
     trainBg = (TTree*)inFile->Get(Form("sidebandChan%dEvents0/events", fChannel));
     testBg  = (TTree*)inFile->Get(Form("sidebandChan%dEvents1/events", fChannel));
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

   factory->AddTree(trainSg,     "Signal",     signalWeight,  "", "train");
   factory->AddTree(testSg,      "Signal",     signalWeight,  "", "test");
   factory->AddTree(trainBg, "Background", cbackgroundWeight, "", "train");
   factory->AddTree(testBg,  "Background", tbackgroundWeight, "", "test");

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

   int seed = static_cast<int>(100*gRandom->Rndm());

   optstring = Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=Block:NormMode=None:V",
                    nSgTrain, nSgTest, nBgTrain, nBgTest);
   cout << "==> PrepareTrainingAndTestTree: " << optstring << endl;
   factory->PrepareTrainingAndTestTree("", "", optstring.c_str());

   if (0) {
     optstring = Form("!H:V:NTrees=%d", fBdtSetup.NTrees);
     optstring += Form(":nCuts=%d:PruneMethod=NoPruning", fBdtSetup.nCuts);
     optstring += Form(":PruneMethod=NoPruning");
     optstring += Form(":BoostType=AdaBoost:AdaBoostBeta=%f:SeparationType=GiniIndex", fBdtSetup.AdaBoostBeta);
     optstring += Form(":MaxDepth=%d", fBdtSetup.MaxDepth);
     optstring += Form(":NNodesMax=%d", fBdtSetup.NNodesMax);
     //     optstring += Form(":nEventsMin=%d", fBdtSetup.nEventsMin);
   } else {
     optstring = "VerbosityLevel=Verbose:" + fBDTParameters;
   }

   cout << "==> BookMethod: " << optstring << endl;
   factory->BookMethod(TMVA::Types::kBDT, "Bmm4BDT", optstring);

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

   cout << "==> Wrote root file: " << outputFile->GetName() << endl;
   cout << "==> TMVAClassification is done!" << endl;

   delete factory;
   gROOT->Clear();  gROOT->DeleteAll();

}



// ----------------------------------------------------------------------
void plotBdt::apply(const char *fname) {

  // --- Book the MVA methods
  string methodName("BDT");
  string dir("weights");
  string XmlName = Form("%s/%s-Events0_%s.weights.xml", dir.c_str(), fname, methodName.c_str());
  fReader.push_back(setupReader(XmlName, frd));
  XmlName = Form("%s/%s-Events1_%s.weights.xml", dir.c_str(), fname, methodName.c_str());
  fReader.push_back(setupReader(XmlName, frd));
  XmlName = Form("%s/%s-Events2_%s.weights.xml", dir.c_str(), fname, methodName.c_str());
  fReader.push_back(setupReader(XmlName, frd));

  // -- open files
  TFile *dfile = fDS["bmmData"]->getFile();
  if (!dfile) {
    cout << "ERROR: could not open data file" << endl;
    exit(1);
  }

  TFile *sfile = fDS["bsmmMcComb"]->getFile();
  if (!sfile) {
    cout << "ERROR: could not open signal file" << endl;
    exit(1);
  }

  TFile f(Form("%s-combined.root", fname), "RECREATE");

  TTree *tt = new TTree("bdtTree", "bdtTree");
  int classID;
  double w8;
  tt->Branch("run",     &fb.run,   "run/I");
  tt->Branch("evt",     &fb.evt,   "evt/I");
  tt->Branch("m",       &fb.m,     "m/D");
  tt->Branch("me",      &fb.me,    "me/D");
  tt->Branch("weight",  &w8,        "weight/D");
  tt->Branch("bdt",     &fBDT,      "bdt/D");
  tt->Branch("bdt0",    &fBDT0,     "bdt0/D");
  tt->Branch("bdt1",    &fBDT1,     "bdt1/D");
  tt->Branch("bdt2",    &fBDT2,     "bdt2/D");
  tt->Branch("classID", &classID,   "classID/I");
  tt->Branch("hlt1",    &fb.hlt1,   "hlt1/O");
  tt->Branch("tos",     &fb.tos,    "tos/O");
  tt->Branch("gmuid",   &fb.gmuid, "gmuid/O");
  tt->Branch("gtqual",  &fb.gtqual,"gtqual/O");
  tt->Branch("json",    &fb.json,  "json/O");
  tt->Branch("pvw8",    &fb.pvw8,  "pvw8/D");

    // -- add the usual variables
  tt->Branch("m1pt",    &fb.m1pt,     "m1pt/D");
  tt->Branch("m2pt",    &fb.m2pt,     "m2pt/D");
  tt->Branch("m1eta",   &fb.m1eta,    "m1eta/D");
  tt->Branch("m2eta",   &fb.m2eta,    "m2etat/D");
  tt->Branch("m1q",     &fb.m1q,      "m1q/I");
  tt->Branch("m2q",     &fb.m2q,      "m2q/I");
  tt->Branch("pt",      &fb.pt,       "pt/D");
  tt->Branch("eta",     &fb.eta,      "eta/D");

  tt->Branch("chi2dof", &fb.chi2dof,  "chi2dof/D");
  tt->Branch("maxdoca", &fb.maxdoca,  "maxdoca/D");
  tt->Branch("fls3d",   &fb.fls3d,  "fls3d/D");
  tt->Branch("fl3d",    &fb.fl3d,  "fl3d/D");
  tt->Branch("flsxy",   &fb.flsxy,  "flsxy/D");
  tt->Branch("alpha",   &fb.alpha,  "alpha/D");
  tt->Branch("pvip",    &fb.pvip,  "pvip/D");
  tt->Branch("pvips",   &fb.pvips,  "pvips/D");
  tt->Branch("pvlip",   &fb.pvlip,  "pvlip/D");
  tt->Branch("pvlips",  &fb.pvlips,  "pvlips/D");

  tt->Branch("iso",     &fb.iso,  "iso/D");
  tt->Branch("m1iso",   &fb.m1iso,  "m1iso/D");
  tt->Branch("m2iso",   &fb.m2iso,  "m2iso/D");
  tt->Branch("docatrk", &fb.docatrk,  "docatrk/D");
  tt->Branch("closetrk",&fb.closetrk,  "closetrk/D");


  TH1D *hd = new TH1D("bdtTreeData", "", 20, 0., 20.);
  TH1D *hs = new TH1D("bdtTreeSignal", "", 20, 0., 20.);


  // -- data processing
  TTree* t = (TTree*)dfile->Get("candAnaMuMu/events");
  setupTree(t, fb);
  w8 = 1.;
  classID = 1;
  cout << "--- Processing data: " << t->GetEntries() << " events" << endl;
  int nEvent = t->GetEntries();
  double lostEvents(0);
  double totalEvents(0);
  bool otherCuts(false);
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    if (ievt%1000000 == 0) std::cout << "--- ... Processing data event: " << ievt << std::endl;
    t->GetEntry(ievt);
    fBDT = fBDT0 = fBDT1 =  fBDT2 = -99.;
    totalEvents += w8;
    otherCuts = (fb.m > 4.9 && fb.m < 5.9) && (fb.m1q*fb.m2q < 0) && (fb.pvw8 > 0.7) && (fb.gtqual);

    if (fb.json && otherCuts && fb.gmuid && fb.hlt1 && fb.tos && preselection(fb)) {
      calcBDT();
      tt->Fill();
    } else {
      lostEvents += w8;
      hd->Fill(TMath::Abs(fb.evt%3), w8);
      hd->Fill(9, w8);
    }
  }
  hd->SetBinContent(10, lostEvents);
  cout << "lost events: " << lostEvents << " out of " << totalEvents << " events in total" << endl;


  // -- signal MC processing
  t = (TTree*)sfile->Get("candAnaMuMu/events");
  classID = 0;
  setupTree(t, fb);
  w8 = fLumiScale;
  cout << "--- Processing signal: " << t->GetEntries() << " events and weight " << w8 << endl;
  nEvent = t->GetEntries();
  lostEvents = 0;
  totalEvents = 0;
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    if (ievt%1000000 == 0) std::cout << "--- ... Processing signal event: " << ievt << std::endl;
    t->GetEntry(ievt);

    fBDT = fBDT0 = fBDT1 =  fBDT2 = -99.;
    totalEvents += w8;
    otherCuts = (fb.m > 4.9 && fb.m < 5.9) && (fb.m1q*fb.m2q < 0) && (fb.pvw8 > 0.7) && (fb.gtqual);
    if (otherCuts && fb.gmuid && fb.hlt1 && fb.tos && preselection(fb)) {
      calcBDT();
      tt->Fill();
    } else {
      lostEvents += w8;
      hs->Fill(TMath::Abs(fb.evt%3), w8);
      hs->Fill(9, w8);
    }
  }
  hs->SetBinContent(10, lostEvents);
  cout << "lost events: " << lostEvents << " out of " << totalEvents << " events in total" << endl;

  hd->Write();
  hs->Write();
  f.Write();
  f.Close();

}



// ----------------------------------------------------------------------
void plotBdt::setupTree(TTree *t, redTreeData &b) {
  t->SetBranchAddress("evt", &b.evt);
  t->SetBranchAddress("run", &b.run);
  t->SetBranchAddress("json", &b.json);
  t->SetBranchAddress("pvw8", &b.pvw8);
  t->SetBranchAddress("gmuid", &b.gmuid);
  t->SetBranchAddress("gmugmid", &b.gmugmid);
  t->SetBranchAddress("gtqual", &b.gtqual);
  t->SetBranchAddress("hlt1", &b.hlt1);
  t->SetBranchAddress("tos", &b.tos);
  t->SetBranchAddress("m1pt", &b.m1pt);
  t->SetBranchAddress("m1eta", &b.m1eta);
  t->SetBranchAddress("m1q", &b.m1q);
  t->SetBranchAddress("m2pt", &b.m2pt);
  t->SetBranchAddress("m2eta", &b.m2eta);
  t->SetBranchAddress("m2q", &b.m2q);
  t->SetBranchAddress("pt", &b.pt);
  t->SetBranchAddress("eta", &b.eta);
  t->SetBranchAddress("pvlip", &b.pvlip);
  t->SetBranchAddress("pvlips", &b.pvlips);
  t->SetBranchAddress("fl3d", &b.fl3d);
  t->SetBranchAddress("fls3d", &b.fls3d);
  t->SetBranchAddress("flsxy", &b.flsxy);
  t->SetBranchAddress("alpha", &b.alpha);
  t->SetBranchAddress("maxdoca", &b.maxdoca);
  t->SetBranchAddress("pvip", &b.pvip);
  t->SetBranchAddress("pvips", &b.pvips);
  t->SetBranchAddress("iso", &b.iso);
  t->SetBranchAddress("docatrk", &b.docatrk);
  t->SetBranchAddress("chi2dof", &b.chi2dof);
  t->SetBranchAddress("closetrk", &b.closetrk);
  t->SetBranchAddress("m", &b.m);

  t->SetBranchAddress("m1iso",&b.m1iso);
  t->SetBranchAddress("m2iso",&b.m2iso);
  t->SetBranchAddress("closetrks1", &b.closetrks1);
  t->SetBranchAddress("closetrks2", &b.closetrks2);
  t->SetBranchAddress("closetrks3", &b.closetrks3);
  t->SetBranchAddress("pvdchi2",&b.pvdchi2);
  t->SetBranchAddress("othervtx",&b.othervtx);
  t->SetBranchAddress("m1xpdist",&b.m1xpdist);
  t->SetBranchAddress("m2xpdist",&b.m2xpdist);

  t->SetBranchAddress("pv2lips", &b.pv2lips);
  t->SetBranchAddress("pv2lip", &b.pv2lip);
}


// ----------------------------------------------------------------------
void plotBdt::loadFiles(string afiles) {

  string files = fDirectory + string("/") + afiles;
  cout << "==> plotBdt::loadFile loading files listed in " << files << endl;

  char buffer[1000];
  ifstream is(files.c_str());
  while (is.getline(buffer, 1000, '\n')) {
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}

    string sbuffer = string(buffer);
    replaceAll(sbuffer, " ", "");
    replaceAll(sbuffer, "\t", "");
    if (sbuffer.size() < 1) continue;

    string::size_type m1 = sbuffer.find("lumi=");
    string stype = sbuffer.substr(5, m1-5);

    string::size_type m2 = sbuffer.find("file=");
    string slumi = sbuffer.substr(m1+5, m2-m1-6);
    string sfile = sbuffer.substr(m2+5);
    string sname("nada"), sdecay("nada");

    TFile *pF(0);
    dataset *ds(0);

    if (string::npos != stype.find("SingleMuon")) {
      // -- SingleMuon
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("bmm")) {
        sname = "bmmSingleMuon";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bupsik")) {
        sname = "bupsikSingleMuon";
        sdecay = "bupsik";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bspsiphi")) {
        sname = "bspsiphiSingleMuon";
        sdecay = "bspsiphi";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bdpsikstar")) {
        sname = "bdpsikstarSingleMuon";
        sdecay = "bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

    } else if (string::npos != stype.find("Charmonium")) {
      // -- Charmonium
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2;
      if (string::npos != stype.find("bmm")) {
        sname = "bmmCharmonium";
        sdecay = "bmm";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bupsik")) {
        sname = "bupsikCharmonium";
        sdecay = "bupsik";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bspsiphi")) {
        sname = "bspsiphiCharmonium";
        sdecay = "bspsiphi";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

      if (string::npos != stype.find("bdpsikstar")) {
        sname = "bdpsikstarCharmonium";
        sdecay = "bdpsikstar";
	ds->fColor = kBlack;
	ds->fSymbol = 20;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3365;
      }

    } else if (string::npos != stype.find("mc")) {
      // -- MC
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2.;

      if (string::npos != stype.find("bctopsimunu,")) {
        sname = "bcpsimunuMc";
        sdecay = "bcpsimunu";
	ds->fColor = kGreen-2;
	ds->fSymbol = 24;
	ds->fWidth  = 2.;
	ds->fF      = pF;
	ds->fBf     = 1.;
	ds->fMass   = 1.;
	ds->fFillStyle = 3354;
      }

    } else if (string::npos != stype.find("relval")) {
      pF = loadFile(sfile);
      ds = new dataset();
      ds->fSize = 1.2;
      ds->fWidth = 2.;
      if (string::npos != stype.find("bsmm,relval")) {
	ds->fF      = pF;
	sname   = "bsmmrelval";
      }

      if (string::npos != stype.find("bdmm,relval")) {
	ds->fF      = pF;
	sname   = "bdmmrelval";
      }

      if (string::npos != stype.find("bupsik,relval")) {
	ds->fF      = pF;
	sname   = "bupsikrelval";
      }

      if (string::npos != stype.find("bspsiphi,relval")) {
	ds->fF      = pF;
	sname   = "bspsiphirelval";
      }
    }

    if (sname != "nada") {
      ds->fLcolor = ds->fColor;
      ds->fFcolor = ds->fColor;
      ds->fName   = sdecay;
      ds->fFullName = sname;
      insertDataset(sname, ds);
    } else {
      delete ds;
    }

  }

  is.close();

  int cnt(0);
  cout << "------------------------------------------------------------------------------------------" << endl;
  cout << Form("   %30s: %20s: ", "Dataset name", "Decay mode name") << "Filename:" << endl;
  cout << "------------------------------------------------------------------------------------------" << endl;
  for (map<string, dataset*>::iterator it = fDS.begin(); it != fDS.end(); ++it) {
    // cout << it->first << endl;
    // cout << it->second->fName << endl;
    // cout << it->second->fF->GetName() << endl;
    cout << Form("%2d %30s: %20s: ", cnt, it->first.c_str(), it->second->fName.c_str()) << it->second->fF->GetName() << endl;
    ++cnt;
  }
  cout << "------------------------------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void plotBdt::writeOut(TFile *f, TH1 *h) {
  TDirectory *pD = gDirectory;
  f->cd();
  h->SetDirectory(f);
  h->Write();
  pD->cd();
}


// ----------------------------------------------------------------------
void plotBdt::redrawStats(double x, double y, const char *newname, int color) {

  TPaveStats *st = (TPaveStats*)gPad->GetPrimitive("stats");
  st->SetName(newname);
  st->SetX1NDC(x);
  st->SetY1NDC(y);
  st->SetTextColor(color);
  //  st->Draw("sames");
}


// ----------------------------------------------------------------------
void plotBdt::cleanup(string fname) {
  cout << "     cleanup of " << Form("%s-Events0.root, %s-Events1.root, %s-Events2.root", fname.c_str(), fname.c_str(), fname.c_str()) << endl;
  TFile *f = TFile::Open(Form("%s-Events0.root", fname.c_str()));
  if (!f) {
    cout << "      file not found " << endl;
    return;
  }
  TH1D* h0 = (TH1D*)f->Get("res_ks");
  if (0 == h0) {
    cout << "      hist res_ks not found " << endl;
    return;
  }

  f = TFile::Open(Form("%s-Events1.root", fname.c_str()));
  if (!f) {
    cout << "      file not found " << endl;
    return;
  }
  TH1D* h1 = (TH1D*)f->Get("res_ks");
  if (0 == h1) {
    cout << "      hist res_ks not found " << endl;
    return;
  }

  f = TFile::Open(Form("%s-Events2.root", fname.c_str()));
  if (!f) {
    cout << "      file not found " << endl;
    return;
  }
  TH1D* h2 = (TH1D*)f->Get("res_ks");
  if (0 == h2) {
    cout << "      hist res_ks not found " << endl;
    return;
  }

  // -- change!
  f = TFile::Open(Form("%s-combined.root", fname.c_str()));
  TH1D *H = (TH1D*)f->Get("res_ssb");
  if (0 == H) {
    cout << "      hist res_ssb not found " << endl;
    return;
  }

  double kssg0 = h0->GetBinContent(1);
  double ksbg0 = h0->GetBinContent(2);
  double kssg1 = h1->GetBinContent(1);
  double ksbg1 = h1->GetBinContent(2);
  double kssg2 = h2->GetBinContent(1);
  double ksbg2 = h2->GetBinContent(2);
  double ssb  = H->GetBinContent(1);
  cout << fname << " performance: kssg0 = " << kssg0 << " ksbg0 = " << ksbg0
       << " kssg1 = " << kssg1 << " ksbg1 = " << ksbg1
       << " kssg2 = " << kssg2 << " ksbg2 = " << ksbg2
       << " ssb = "  << ssb << endl;
  double ssbCut = 1.3;
  if (1 == fChannel) ssbCut = 0.9;
  if (kssg0 < 0.05 || ksbg0 < 0.05 || kssg1 < 0.05 || ksbg1 < 0.05 || kssg2 < 0.05 || ksbg2 < 0.05 || ssb < ssbCut) {
    system(Form("/bin/rm -f %s-Events0.root", fname.c_str()));
    system(Form("/bin/rm -f %s-Events1.root", fname.c_str()));
    system(Form("/bin/rm -f %s-Events2.root", fname.c_str()));
    cout << Form("===> REMOVED %s-[Events0,Events1,Events2].root, kssg0 = %f ksbg0 = %f kssg1 = %f ksbg1 = %f kssg2 = %f ksbg2 = %f SSB = %f",
                 fname.c_str(), kssg0, ksbg0, kssg1, ksbg1, kssg2, ksbg2, ssb) << endl;
    system(Form("/bin/rm -f weights/%s-Events0_BDT.weights.xml", fname.c_str()));
    system(Form("/bin/rm -f weights/%s-Events1_BDT.weights.xml", fname.c_str()));
    system(Form("/bin/rm -f weights/%s-Events2_BDT.weights.xml", fname.c_str()));
    cout << Form("     REMOVED weights/%s-[Events0,Events1,Events2]_BDT.weights.xml", fname.c_str()) << endl;
    system(Form("/bin/rm -f weights/%s-Events0_BDT.class.C", fname.c_str()));
    system(Form("/bin/rm -f weights/%s-Events1_BDT.class.C", fname.c_str()));
    system(Form("/bin/rm -f weights/%s-Events2_BDT.class.C", fname.c_str()));
    cout << Form("     REMOVED weights/%s-[Events0,Events1,Events2]_BDT.class.C", fname.c_str()) << endl;
  } else {
    cout << Form("===> KEEP %s.root, kssg0 = %f ksbg0 = %f kssg1 = %f ksbg1 = %f kssg0 = %f ksbg0 = %f SSB = %f",
                 fname.c_str(), kssg0, ksbg0, kssg1, ksbg1, kssg2, ksbg2, ssb) << endl;
  }

}


// ----------------------------------------------------------------------
void plotBdt::analyze(string fname) {

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  int bdtBins(200);
  double bdtMin(-1.0), bdtMax(1.0);

  TH1D *sm = new TH1D("sm", "m(signal)", 100, 4.9, 5.9); sm->Sumw2();
  TH1D *dm = new TH1D("dm", "m(data)", 100, 4.9, 5.9); dm->Sumw2();

  TH1D *sm0 = new TH1D("sm0", "m(signal)", 100, 4.9, 5.9); sm0->Sumw2();
  TH1D *dm0 = new TH1D("dm0", "m(data)", 100, 4.9, 5.9); dm0->Sumw2();

  TH1D *sm1 = new TH1D("sm1", "m(signal)", 100, 4.9, 5.9); sm1->Sumw2();
  TH1D *dm1 = new TH1D("dm1", "m(data)", 100, 4.9, 5.9); dm1->Sumw2();

  TH1D *sm2 = new TH1D("sm2", "m(signal)", 100, 4.9, 5.9); sm2->Sumw2();
  TH1D *dm2 = new TH1D("dm2", "m(data)", 100, 4.9, 5.9); dm2->Sumw2();

  cout << "Open " << Form("%s-combined.root", fname.c_str()) << endl;
  TFile *f = TFile::Open(Form("%s-combined.root", fname.c_str()), "UPDATE");
  TH1D *hd = (TH1D*)f->Get("bdtTreeData");
  double dLostEvents  = hd->GetBinContent(10);
  TH1D *hs = (TH1D*)f->Get("bdtTreeSignal");
  double sLostEvents  = hs->GetBinContent(10);

  TH1D *h = new TH1D("s1", "S/Sqrt(S+B)", bdtBins, bdtMin, bdtMax); h->Sumw2();
  setTitles(h, "b > ", "S/#sqrt{S+B}");
  TH1D *hs2 = new TH1D("s2", "S/Sqrt(S+B)", bdtBins, bdtMin, bdtMax); hs2->Sumw2();
  setTitles(hs2, "b > ", "S/#sqrt{S+B'}");
  TGraph *grocop = new TGraph(1); grocop->SetMarkerStyle(34); grocop->SetMarkerSize(3.); grocop->SetMarkerColor(kBlack);
  TGraph *groc   = new TGraph(bdtBins);   groc->SetMarkerStyle(24);  groc->SetMarkerSize(1.5);    groc->SetMarkerColor(kBlack);
  TGraph *groc0  = new TGraph(bdtBins);   groc0->SetMarkerStyle(20); groc0->SetMarkerSize(0.8);   groc0->SetMarkerColor(kMagenta);
  TGraph *groc1  = new TGraph(bdtBins);   groc1->SetMarkerStyle(20); groc1->SetMarkerSize(0.8);   groc1->SetMarkerColor(kRed);
  TGraph *groc2  = new TGraph(bdtBins);   groc2->SetMarkerStyle(20); groc2->SetMarkerSize(0.8);   groc2->SetMarkerColor(kBlue);

  TH1D *hroc = new TH1D("hroc", "", 200, 0., 1.);
  TH1D *hroc0 = new TH1D("hroc0", "", 200, 0., 1.);
  TH1D *hroc1 = new TH1D("hroc1", "", 200, 0., 1.);
  TH1D *hroc2 = new TH1D("hroc2", "", 200, 0., 1.);
  TFile *f0 = TFile::Open(Form("%s-Events0.root", fname.c_str()));
  TFile *f1 = TFile::Open(Form("%s-Events1.root", fname.c_str()));
  TFile *f2 = TFile::Open(Form("%s-Events2.root", fname.c_str()));

  if (!f0 || !f1 || !f2) {
    cout << "ERROR: could not open a file" << endl;
    exit(1);
  }

  TH1D *hr01 = getRanking(fname.c_str(), "IdTransformation", "events0");
  hr01->SetDirectory(f);
  TH1D *hr02 =  getRanking(fname.c_str(), "BDT", "events0");
  hr02->SetDirectory(f);

  TH1D *hr11 = getRanking(fname.c_str(), "IdTransformation", "events1");
  hr11->SetDirectory(f);
  TH1D *hr12 =  getRanking(fname.c_str(), "BDT", "events1");
  hr12->SetDirectory(f);

  TH1D *hr21 = getRanking(fname.c_str(), "IdTransformation", "events2");
  hr21->SetDirectory(f);
  TH1D *hr22 =  getRanking(fname.c_str(), "BDT", "events2");
  hr22->SetDirectory(f);

  TH1F *trainBDT0 = (TH1F*)f0->Get("Method_BDT/BDT/MVA_BDT_Train_B"); trainBDT0->SetLineColor(kBlack);
  TH1F *testBDT0 = (TH1F*)f0->Get("Method_BDT/BDT/MVA_BDT_B"); testBDT0->SetLineColor(kBlack);

  TH1F *trainBDT1 = (TH1F*)f1->Get("Method_BDT/BDT/MVA_BDT_Train_B"); trainBDT1->SetLineColor(kRed);
  TH1F *testBDT1 = (TH1F*)f1->Get("Method_BDT/BDT/MVA_BDT_B"); testBDT1->SetLineColor(kRed);

  TH1F *trainBDT2 = (TH1F*)f2->Get("Method_BDT/BDT/MVA_BDT_Train_B"); trainBDT2->SetLineColor(kBlue);
  TH1F *testBDT2 = (TH1F*)f2->Get("Method_BDT/BDT/MVA_BDT_B"); testBDT2->SetLineColor(kBlue);

  TH1F *trainSgBDT0 = (TH1F*)f0->Get("Method_BDT/BDT/MVA_BDT_Train_S"); trainSgBDT0->SetLineColor(kBlack);
  TH1F *testSgBDT0 = (TH1F*)f0->Get("Method_BDT/BDT/MVA_BDT_S"); testSgBDT0->SetLineColor(kBlack);

  TH1F *trainSgBDT1 = (TH1F*)f1->Get("Method_BDT/BDT/MVA_BDT_Train_S"); trainSgBDT1->SetLineColor(kRed);
  TH1F *testSgBDT1 = (TH1F*)f1->Get("Method_BDT/BDT/MVA_BDT_S"); testSgBDT1->SetLineColor(kRed);

  TH1F *trainSgBDT2 = (TH1F*)f2->Get("Method_BDT/BDT/MVA_BDT_Train_S"); trainSgBDT2->SetLineColor(kBlue);
  TH1F *testSgBDT2 = (TH1F*)f2->Get("Method_BDT/BDT/MVA_BDT_S"); testSgBDT2->SetLineColor(kBlue);

  TH1D *bBDT = new TH1D("bBDT", "", trainBDT0->GetNbinsX(), trainBDT0->GetBinLowEdge(1), trainBDT0->GetBinLowEdge(trainBDT0->GetNbinsX()+1));
  TH1D *cBDT = new TH1D("cBDT", "", trainBDT0->GetNbinsX(), trainBDT0->GetBinLowEdge(1), trainBDT0->GetBinLowEdge(trainBDT0->GetNbinsX()+1));
  TH1D *h2 = new TH1D("res_ssb", "max(S/Sqrt(S+B))", 10, 0., 10.); h2->Sumw2();  h2->SetDirectory(f);

  // -- new rebinned versions of the BDT output plots. These ARE identical to the TMVA output IF the binning IS the same!
  const int BINS(100);
  TH1D *ap0bgBDT = new TH1D("ap0bgBDT", "",    BINS, -1., 1.);
  TH1D *tr0bgBDT = new TH1D("tr0bgBDT", "",    BINS, -1., 1.);
  TH1D *te0bgBDT = new TH1D("te0bgBDT", "",    BINS, -1., 1.);
  ap0bgBDT->SetLineColor(kMagenta);   te0bgBDT->SetLineColor(kMagenta);   tr0bgBDT->SetLineColor(kMagenta);
  ap0bgBDT->SetMarkerColor(kMagenta); te0bgBDT->SetMarkerColor(kMagenta); tr0bgBDT->SetMarkerColor(kMagenta);

  TH1D *ap1bgBDT = new TH1D("ap1bgBDT", "",    BINS, -1., 1.);
  TH1D *tr1bgBDT = new TH1D("tr1bgBDT", "",    BINS, -1., 1.);
  TH1D *te1bgBDT = new TH1D("te1bgBDT", "",    BINS, -1., 1.);
  ap1bgBDT->SetLineColor(kRed);   te1bgBDT->SetLineColor(kRed);   tr1bgBDT->SetLineColor(kRed);
  ap1bgBDT->SetMarkerColor(kRed); te1bgBDT->SetMarkerColor(kRed); tr1bgBDT->SetMarkerColor(kRed);

  TH1D *ap2bgBDT = new TH1D("ap2bgBDT", "",    BINS, -1., 1.);
  TH1D *tr2bgBDT = new TH1D("tr2bgBDT", "",    BINS, -1., 1.);
  TH1D *te2bgBDT = new TH1D("te2bgBDT", "",    BINS, -1., 1.);
  ap2bgBDT->SetLineColor(kBlue);   te2bgBDT->SetLineColor(kBlue);   tr2bgBDT->SetLineColor(kBlue);
  ap2bgBDT->SetMarkerColor(kBlue); te2bgBDT->SetMarkerColor(kBlue); tr2bgBDT->SetMarkerColor(kBlue);

  TH1D *ap3bgBDT = new TH1D("ap3bgBDT", "",    BINS, -1., 1.);
  ap3bgBDT->SetLineColor(kBlack);
  ap3bgBDT->SetMarkerColor(kBlack);

  TH1D *ap0sgBDT = new TH1D("ap0sgBDT", "",    BINS, -1., 1.);
  TH1D *tr0sgBDT = new TH1D("tr0sgBDT", "",    BINS, -1., 1.);
  TH1D *te0sgBDT = new TH1D("te0sgBDT", "",    BINS, -1., 1.);
  ap0sgBDT->SetLineColor(kMagenta);   te0sgBDT->SetLineColor(kMagenta);   tr0sgBDT->SetLineColor(kMagenta);
  ap0sgBDT->SetMarkerColor(kMagenta); te0sgBDT->SetMarkerColor(kMagenta); tr0sgBDT->SetMarkerColor(kMagenta);

  TH1D *ap1sgBDT = new TH1D("ap1sgBDT", "",    BINS, -1., 1.);
  TH1D *tr1sgBDT = new TH1D("tr1sgBDT", "",    BINS, -1., 1.);
  TH1D *te1sgBDT = new TH1D("te1sgBDT", "",    BINS, -1., 1.);
  ap1sgBDT->SetLineColor(kRed);   te1sgBDT->SetLineColor(kRed);   tr1sgBDT->SetLineColor(kRed);
  ap1sgBDT->SetMarkerColor(kRed); te1sgBDT->SetMarkerColor(kRed); tr1sgBDT->SetMarkerColor(kRed);

  TH1D *ap2sgBDT = new TH1D("ap2sgBDT", "",    BINS, -1., 1.);
  TH1D *tr2sgBDT = new TH1D("tr2sgBDT", "",    BINS, -1., 1.);
  TH1D *te2sgBDT = new TH1D("te2sgBDT", "",    BINS, -1., 1.);
  ap2sgBDT->SetLineColor(kBlue);   te2sgBDT->SetLineColor(kBlue);   tr2sgBDT->SetLineColor(kBlue);
  ap2sgBDT->SetMarkerColor(kBlue); te2sgBDT->SetMarkerColor(kBlue); tr2sgBDT->SetMarkerColor(kBlue);

  TH1D *ap3sgBDT = new TH1D("ap3sgBDT", "",    BINS, -1., 1.);
  ap3sgBDT->SetLineColor(kBlack);
  ap3sgBDT->SetMarkerColor(kBlack);

  TTree *t = (TTree*)f->Get("bdtTree");
  cout << "bdtTree with entries = " << t->GetEntries() << endl;
  double bdt, m, w8;
  double bdt0, bdt1, bdt2;
  int classID, evt;
  bool gmuid, hlt1, tos;
  t->SetBranchAddress("bdt", &bdt);
  t->SetBranchAddress("bdt0", &bdt0);
  t->SetBranchAddress("bdt1", &bdt1);
  t->SetBranchAddress("bdt2", &bdt2);
  t->SetBranchAddress("classID", &classID);
  t->SetBranchAddress("m", &m);
  t->SetBranchAddress("weight", &w8);
  t->SetBranchAddress("hlt1", &hlt1);
  t->SetBranchAddress("tos", &tos);
  t->SetBranchAddress("gmuid", &gmuid);
  t->SetBranchAddress("evt", &evt);

  // -- data (overall) distribution
  int nEvent(0);
  nEvent = t->GetEntries();
  for (Long64_t ievt=0; ievt<nEvent; ievt++) {
    t->GetEntry(ievt);
    if (1 == classID) {
      ap0bgBDT->Fill(bdt0);
      ap1bgBDT->Fill(bdt1);
      ap2bgBDT->Fill(bdt2);
      ap3bgBDT->Fill(bdt);

      if (evt%3==0) {
	te1bgBDT->Fill(bdt1);
	tr2bgBDT->Fill(bdt2);
      }
      if (evt%3==1) {
	te2bgBDT->Fill(bdt2);
	tr0bgBDT->Fill(bdt0);
      }
      if (evt%3==2) {
	te0bgBDT->Fill(bdt0);
	tr1bgBDT->Fill(bdt1);
      }
    }

    if (0 == classID) {
      ap0sgBDT->Fill(bdt0);
      ap1sgBDT->Fill(bdt1);
      ap2sgBDT->Fill(bdt2);
      ap3sgBDT->Fill(bdt);

      if (evt%3==0) {
	te1sgBDT->Fill(bdt1);
	tr2sgBDT->Fill(bdt2);
      }
      if (evt%3==1) {
	te2sgBDT->Fill(bdt2);
	tr0sgBDT->Fill(bdt0);
      }
      if (evt%3==2) {
	te0sgBDT->Fill(bdt0);
	tr1sgBDT->Fill(bdt1);
      }
    }
  }

  trainBDT0->Scale(1./trainBDT0->GetSumOfWeights());
  testBDT0->Scale(1./testBDT0->GetSumOfWeights());
  trainBDT1->Scale(1./trainBDT1->GetSumOfWeights());
  testBDT1->Scale(1./testBDT1->GetSumOfWeights());
  trainBDT2->Scale(1./trainBDT2->GetSumOfWeights());
  testBDT2->Scale(1./testBDT2->GetSumOfWeights());

  ap0bgBDT->Scale(1./ap0bgBDT->GetSumOfWeights());
  tr0bgBDT->Scale(1./tr0bgBDT->GetSumOfWeights());
  te0bgBDT->Scale(1./te0bgBDT->GetSumOfWeights());

  ap1bgBDT->Scale(1./ap1bgBDT->GetSumOfWeights());
  tr1bgBDT->Scale(1./tr1bgBDT->GetSumOfWeights());
  te1bgBDT->Scale(1./te1bgBDT->GetSumOfWeights());

  ap2bgBDT->Scale(1./ap2bgBDT->GetSumOfWeights());
  tr2bgBDT->Scale(1./tr2bgBDT->GetSumOfWeights());
  te2bgBDT->Scale(1./te2bgBDT->GetSumOfWeights());

  ap0sgBDT->Scale(1./ap0sgBDT->GetSumOfWeights());
  tr0sgBDT->Scale(1./tr0sgBDT->GetSumOfWeights());
  te0sgBDT->Scale(1./te0sgBDT->GetSumOfWeights());

  ap1sgBDT->Scale(1./ap1sgBDT->GetSumOfWeights());
  tr1sgBDT->Scale(1./tr1sgBDT->GetSumOfWeights());
  te1sgBDT->Scale(1./te1sgBDT->GetSumOfWeights());

  ap2sgBDT->Scale(1./ap2sgBDT->GetSumOfWeights());
  tr2sgBDT->Scale(1./tr2sgBDT->GetSumOfWeights());
  te2sgBDT->Scale(1./te2sgBDT->GetSumOfWeights());

  ap3sgBDT->Scale(1./ap3sgBDT->GetSumOfWeights());
  ap3bgBDT->Scale(1./ap3bgBDT->GetSumOfWeights());

  double hmax = tr0bgBDT->GetMaximum();
  if (tr1bgBDT->GetMaximum() > hmax) hmax = tr1bgBDT->GetMaximum();
  if (tr2bgBDT->GetMaximum() > hmax) hmax = tr2bgBDT->GetMaximum();
  tr0bgBDT->SetMaximum(1.3*hmax);

  /*
  tr0bgBDT->SetMarkerStyle(20);
  tr0bgBDT->SetMarkerSize(0.7);
  tr0bgBDT->Draw("p");
  te0bgBDT->Draw("same");

  tr1bgBDT->SetMarkerStyle(20);
  tr1bgBDT->SetMarkerSize(0.7);
  tr1bgBDT->Draw("psame");
  te1bgBDT->Draw("same");

  tr2bgBDT->SetMarkerStyle(20);
  tr2bgBDT->SetMarkerSize(0.7);
  tr2bgBDT->Draw("psame");
  te2bgBDT->Draw("same");
  */

  c0->Clear();
  shrinkPad(0.15, 0.15);
  setTitles(ap0bgBDT, "b", "a.u.", 0.05, 1.1, 1.5);
  ap0bgBDT->Draw();
  ap1bgBDT->Draw("same");
  ap2bgBDT->Draw("same");
  ap3bgBDT->Draw("same");

  newLegend(0.60, 0.67, 0.93, 0.87);
  legg->SetTextSize(0.035);
  legg->SetHeader("Background events");
  legg->AddEntry(ap3sgBDT, "combined", "l");
  legg->AddEntry(ap0sgBDT, "BDT 0", "l");
  legg->AddEntry(ap1sgBDT, "BDT 1", "l");
  legg->AddEntry(ap2sgBDT, "BDT 2", "l");
  legg->Draw();

  c0->SaveAs(Form("plots/%s-rebinned-bg-overlays.pdf", fname.c_str()));

  hmax = tr0sgBDT->GetMaximum();
  if (tr1sgBDT->GetMaximum() > hmax) hmax = tr1sgBDT->GetMaximum();
  if (tr2sgBDT->GetMaximum() > hmax) hmax = tr2sgBDT->GetMaximum();
  tr0sgBDT->SetMaximum(1.3*hmax);

  /*
  tr0sgBDT->SetMarkerStyle(20);
  tr0sgBDT->SetMarkerSize(0.7);
  tr0sgBDT->Draw("p");
  te0sgBDT->Draw("same");

  tr1sgBDT->SetMarkerStyle(20);
  tr1sgBDT->SetMarkerSize(0.7);
  tr1sgBDT->Draw("psame");
  te1sgBDT->Draw("same");

  tr2sgBDT->SetMarkerStyle(20);
  tr2sgBDT->SetMarkerSize(0.7);
  tr2sgBDT->Draw("psame");
  te2sgBDT->Draw("same");
  */

  shrinkPad(0.15, 0.15);
  setTitles(ap0sgBDT, "b", "a.u.", 0.05, 1.1, 1.5);
  ap0sgBDT->Draw();
  ap1sgBDT->Draw("same");
  ap2sgBDT->Draw("same");
  ap3sgBDT->Draw("same");

  newLegend(0.20, 0.67, 0.5, 0.87);

  legg->SetTextSize(0.035);
  legg->SetHeader("Signal events");
  legg->AddEntry(ap3sgBDT, "combined", "l");
  legg->AddEntry(ap0sgBDT, "BDT 0", "l");
  legg->AddEntry(ap1sgBDT, "BDT 1", "l");
  legg->AddEntry(ap2sgBDT, "BDT 2", "l");
  legg->Draw();

  c0->SaveAs(Form("plots/%s-rebinned-sg-overlays.pdf", fname.c_str()));


  writeOut(f, ap0bgBDT);
  writeOut(f, tr0bgBDT);
  writeOut(f, te0bgBDT);

  writeOut(f, ap1bgBDT);
  writeOut(f, tr1bgBDT);
  writeOut(f, te1bgBDT);

  writeOut(f, ap2bgBDT);
  writeOut(f, tr2bgBDT);
  writeOut(f, te2bgBDT);

  writeOut(f, ap0sgBDT);
  writeOut(f, tr0sgBDT);
  writeOut(f, te0sgBDT);

  writeOut(f, ap1sgBDT);
  writeOut(f, tr1sgBDT);
  writeOut(f, te1sgBDT);

  writeOut(f, ap2sgBDT);
  writeOut(f, tr2sgBDT);
  writeOut(f, te2sgBDT);


  // -- new world: compute S, B, and SSB for each event category separately
  gStyle->SetOptTitle(1);
  int vDmax[4], vDhiMax[4];
  double vSeffTot[4], vSeffMax[4], vSmax[4];
  double vDeffMax[4];
  double vSeffBinD99[4], vBmax[4], vMaxSSB[4], vMaxBDT[4];
  double vMaxSSBsimple[4], vMaxBDTsimple[4];
  double vMaxSSBfit[4], vMaxBDTfit[4];

  TH1D *hvsg[4], *hvbg[4];
  TH1D *hvsgAll[4], *hvbgAll[4];

  TH1D *hrocs[4];
  TH1D *hssb[4], *h2ssb[4];
  TGraph *grocs[4];
  TGraph *grocsOp[4];
  for (int ie = 0; ie < 4; ++ie) {
    vDhiMax[ie] = vDmax[ie] = 0;

    vSeffTot[ie] = vSeffMax[ie] = vSmax[ie] =
      vDeffMax[ie] =
      vSeffBinD99[ie] = vBmax[ie] = vMaxSSB[ie] = vMaxBDT[ie] =
      vMaxSSBsimple[ie] = vMaxBDTsimple[ie] =
      vMaxSSBfit[ie] = vMaxBDTfit[ie] = 0;

    hrocs[ie] = new TH1D(Form("hroc%d", ie), Form("hroc%d", ie), 200, 0., 1.);
    hssb[ie]  = new TH1D(Form("hssb%d", ie), Form("hssb%d", ie), bdtBins, bdtMin, bdtMax); hssb[ie]->Sumw2();
    setTitles(hssb[ie], "b > ", "S/#sqrt{S+B}");
    h2ssb[ie] = new TH1D(Form("h2ssb%d", ie), "S/Sqrt(S+B)", bdtBins, bdtMin, bdtMax); h2ssb[ie]->Sumw2();
    setTitles(h2ssb[ie], "b > ", "S/#sqrt{S+B'}");

    grocs[ie]   = new TGraph(bdtBins);
    grocsOp[ie] = new TGraph(1); grocsOp[ie]->SetMarkerStyle(kOpenCross); grocsOp[ie]->SetMarkerSize(3.); grocsOp[ie]->SetMarkerColor(kBlack);

    grocs[ie]->SetMarkerStyle(20); grocs[ie]->SetMarkerSize(0.8);
    if (0 == ie) grocs[ie]->SetMarkerColor(kMagenta);
    if (1 == ie) grocs[ie]->SetMarkerColor(kRed);
    if (2 == ie) grocs[ie]->SetMarkerColor(kBlue);
    if (3 == ie) grocs[ie]->SetMarkerColor(kBlack);

    hvsg[ie] = new TH1D(Form("hvsg%d", ie), Form("hvsg%d", ie), 100, 4.9, 5.9); hvsg[ie]->Sumw2();
    hvbg[ie] = new TH1D(Form("hvbg%d", ie), Form("hvbg%d", ie), 100, 4.9, 5.9); hvbg[ie]->Sumw2();

    hvsgAll[ie] = new TH1D(Form("hvsgAll%d", ie), Form("hvsgAll%d", ie), 100, 4.9, 5.9); hvsgAll[ie]->Sumw2();
    hvbgAll[ie] = new TH1D(Form("hvbgAll%d", ie), Form("hvbgAll%d", ie), 100, 4.9, 5.9); hvbgAll[ie]->Sumw2();
  }

  int ibin(0);
  double bdtCut;
  for (ibin = bdtBins; ibin >=0; --ibin) {
    bdtCut = bdtMin + ibin*(bdtMax-bdtMin)/bdtBins;
    cout << "+-+-+-+-+-+-+- bin " << ibin << " cutting at bdt > " << bdtCut << endl;
    for (int ie = 0; ie < 4; ++ie) {
      hvsg[ie]->Reset();
      hvbg[ie]->Reset();
      hvsgAll[ie]->Reset();
      hvbgAll[ie]->Reset();
    }

    // -- loop over tree for BDT cut
    nEvent = t->GetEntries();
    for (Long64_t ievt=0; ievt<nEvent; ievt++) {
      t->GetEntry(ievt);
      if (0 == classID) {
	hvsgAll[3]->Fill(m, w8);
	if (0 == TMath::Abs(evt%3)) hvsgAll[0]->Fill(m, w8);
	if (1 == TMath::Abs(evt%3)) hvsgAll[1]->Fill(m, w8);
	if (2 == TMath::Abs(evt%3)) hvsgAll[2]->Fill(m, w8);
      } else {
	hvbgAll[3]->Fill(m, w8);
	if (0 == TMath::Abs(evt%3)) hvbgAll[0]->Fill(m, w8);
	if (1 == TMath::Abs(evt%3)) hvbgAll[1]->Fill(m, w8);
	if (2 == TMath::Abs(evt%3)) hvbgAll[2]->Fill(m, w8);
      }
      if (bdt <= bdtCut) continue;
      if (0 == classID) {
	hvsg[3]->Fill(m, w8);
	if (0 == TMath::Abs(evt%3)) hvsg[0]->Fill(m, w8);
	if (1 == TMath::Abs(evt%3)) hvsg[1]->Fill(m, w8);
	if (2 == TMath::Abs(evt%3)) hvsg[2]->Fill(m, w8);
      } else {
	hvbg[3]->Fill(m, w8);
	if (0 == TMath::Abs(evt%3)) hvbg[0]->Fill(m, w8);
	if (1 == TMath::Abs(evt%3)) hvbg[1]->Fill(m, w8);
	if (2 == TMath::Abs(evt%3)) hvbg[2]->Fill(m, w8);
      }
    }

    // -- determine per event type the characteristics
    for (int ie = 0; ie < 4; ++ie) {
      double sCnt = hvsgAll[ie]->Integral(1, hvsgAll[ie]->GetNbinsX());
      double dCnt = hvbgAll[ie]->Integral(1, hvbgAll[ie]->GetNbinsX());

      h = hvsg[ie];
      double s = h->Integral(h->FindBin(5.3), h->FindBin(5.45));

      h = hvbg[ie];
      double d = h->Integral(1, h->GetNbinsX());
      double dhi = h->Integral(h->FindBin(5.45), h->GetNbinsX());
      double bsimple = d*(5.45-5.30)/(5.9-4.9-0.25);

      double seff = s/sCnt;
      double deff = bsimple/(dLostEvents + dCnt);
      double pbg  = 0.07*s;

      double b = 0.1; // FIXME!!!!!!! bgBlind(h, 3, 4.9, 5.9);

      if ((1.-deff) > 0.999998) vSeffBinD99[ie] = seff;
      grocs[ie]->SetPoint(ibin, seff, 1.-deff);
      hrocs[ie]->SetBinContent(hrocs[ie]->FindBin(seff), 1.-deff);

      double ssb(0.);
      cout << "** " << ie << "** bdt> " << bdtCut << " d = " << d << " s = " << s << " b = " << b
	   << " seff = " << seff << " deff = " << deff << endl;
      if (s+b+pbg >0) {
	ssb = s/TMath::Sqrt(s+b+pbg);
	cout << "** ** ssb " << ssb << endl;
	if (ssb > vMaxSSB[ie]) {
	  vSmax[ie] = s;
	  vDmax[ie] = static_cast<int>(d);
	  vDhiMax[ie] = static_cast<int>(dhi);
	  vBmax[ie] = b+pbg;
	  vMaxSSB[ie] = ssb;
	  vMaxBDT[ie] = bdtCut;
	  vSeffMax[ie] = seff;
	  vDeffMax[ie] = deff;
	  vSeffTot[ie] = s/(sLostEvents + sCnt);
	}
	hssb[ie]->SetBinContent(ibin, ssb);
      } else {
	hssb[ie]->SetBinContent(ibin, 0);
      }

      double ssbs(0.);
      if (s+bsimple >0) {
	ssbs = s/TMath::Sqrt(s+bsimple);
	cout << "** ** ssbs " << ssbs << endl;
	if (ssbs > vMaxSSBsimple[ie]) {
	  vMaxSSBsimple[ie] = ssbs;
	  vMaxBDTsimple[ie] = bdtCut;
	}
	h2ssb[ie]->SetBinContent(ibin, ssbs);
      } else {
	h2ssb[ie]->SetBinContent(ibin, 0);
      }


      //    cout << "S = " << s << " B = " << b << " => S/sqrt(S+B) = " << s/TMath::Sqrt(s+b) << endl;
      c0->Clear();
      if (3 == ie) {
	hvbg[ie]->SetTitle(Form("evt type %d, bdt >%3.2f, S/B/D = %5.2f/%5.2f/%5.0f ssb = %4.3f/%4.3f",
				ie, bdtCut, s, b, d, ssb, ssbs));
	hvbg[ie]->Draw("e");
	hvsg[ie]->Draw("samehist");
	c0->Modified();
	c0->Update();
      }
    }
  }

  for (int ie = 0; ie < 4; ++ie) {
    grocsOp[ie]->SetPoint(0, vSeffMax[ie], 1.-vDeffMax[ie]);
  }

  // -- patch empty bins
  double okVal(1.);
  for (int ie = 0; ie < 4; ++ie) {
    for (int i = 1; i< hrocs[ie]->GetNbinsX(); ++i) {
      if (hrocs[ie]->GetBinContent(i) < 1.e-5) hrocs[ie]->SetBinContent(i, okVal);
      okVal = hrocs[ie]->GetBinContent(i);
    }

    for (int i = hrocs[ie]->GetNbinsX()-1; i > 0; --i) {
      if (!(hrocs[ie]->GetBinContent(i) < hrocs[ie]->GetBinContent(i-1))) {
	hrocs[ie]->SetBinContent(i, 0);
      } else {
	break;
      }
    }
  }


  // -- fit for the maximum ssb and bdt cut
  initFunc *pFunc  = new initFunc();
  gStyle->SetOptFit(0);
  for (int ie = 0; ie < 4; ++ie) {
    double xmax(0.), xmin(0.);
    int nbins(0);
    double maxVal = -1.;
    for (int i = 1; i < hssb[ie]->GetNbinsX(); ++i)
      if (hssb[ie]->GetBinContent(i) > maxVal) maxVal = hssb[ie]->GetBinContent(i);

    for (int i = 1; i < hssb[ie]->GetNbinsX(); ++i) {
      if (hssb[ie]->GetBinContent(i) > 0.7*maxVal) {
	xmax = hssb[ie]->GetBinCenter(i);
	nbins = i - hssb[ie]->GetMaximumBin();
      }
      hssb[ie]->SetBinError(i, 0.03*hssb[ie]->GetBinContent(i));
    }
    xmin = hssb[ie]->GetBinCenter(hssb[ie]->GetMaximumBin() - TMath::Abs(nbins));

    cout << "maxval: " << hssb[ie]->GetMaximum() << endl;
    cout << "maxbin: " << hssb[ie]->GetMaximumBin() << endl;
    cout << "xmax: " << xmax << endl;
    cout << "xmin: " << xmin << endl;
    cout << "nbins: " << nbins << endl;

    TF1 *f1 = pFunc->pol2local(hssb[ie], 0.05);
    hssb[ie]->Fit(f1, "r", "", xmin, xmax);
    double maxfitssbX = hssb[ie]->GetFunction("iF_pol2local")->GetParameter(2);
    double maxfitssbY = hssb[ie]->GetFunction("iF_pol2local")->GetParameter(0);
    vMaxSSBfit[ie] = maxfitssbX;
    vMaxBDTfit[ie] = maxfitssbY;
  }
  delete pFunc;

  // -- BG overlays
  // --------------
  gStyle->SetOptTitle(0);
  c0->Clear();
  TH2F* frame(0);
  frame = new TH2F("frame", "BDT output distributions", 100, 0., 0.65, 100, 0.99999, 1.000001);
  frame->GetXaxis()->SetTitle(" #epsilon_{S}");
  frame->GetYaxis()->SetTitle(" 1 - #epsilon_{B}");

  string texname = fname + ".tex";
  system(Form("/bin/rm -f %s", texname.c_str()));


  int dMaxSum(0);
  double bMaxSum(0.), sMaxSum(0.);
  double a(0.), dMaxBDT(0.), dMaxBDT3(0.);

  dMaxBDT = TMath::Abs(vMaxBDT[0] - vMaxBDT[1]);
  a = TMath::Abs(vMaxBDT[0] - vMaxBDT[2]);
  if (a > dMaxBDT) dMaxBDT = a;
  a = TMath::Abs(vMaxBDT[1] - vMaxBDT[2]);
  if (a > dMaxBDT) dMaxBDT = a;

  a = 0.;
  for (int ie = 0; ie < 4; ++ie) {

    if (ie < 3) {
      a = TMath::Abs(vMaxBDT[3] - vMaxBDT[ie]);
      if (a > dMaxBDT3) dMaxBDT3 = a;
    }

    c0->Clear();
    frame->Draw();
    gPad->SetLogy(0);
    gPad->SetLeftMargin(0.15);
    grocs[ie]->Draw("p");
    grocs[ie]->GetXaxis()->SetTitle("#epsilon_{ S}");
    grocs[ie]->GetYaxis()->SetTitle("1 - #epsilon_{ B}");
    double rocInt = hrocs[ie]->Integral(1, hrocs[ie]->GetNbinsX())*hrocs[ie]->GetBinWidth(1);
    double rocInt2= hrocs[ie]->Integral(1, hrocs[ie]->FindBin(vSeffBinD99[ie]))*hrocs[ie]->GetBinWidth(1);
    cout << " ==> seffBinD99 " << vSeffBinD99[ie] << " in bin " << hrocs[ie]->FindBin(vSeffBinD99[ie]) << endl;
    grocs[ie]->SetName("groc");
    grocs[ie]->SetTitle(Form("integral = %5.3f", rocInt));

    grocsOp[ie]->SetName(Form("grocsOp%d", ie));
    grocs[ie]->SetName(Form("groc%d", ie));

    grocs[ie]->Draw("p");
    grocsOp[ie]->Draw("p");

    int dMax = vDmax[ie];
    int dhiMax = vDhiMax[ie];
    double bMax = vBmax[ie];
    double sMax = vSmax[ie];
    double maxSSB = vMaxSSB[ie];
    double maxBDT = vMaxBDT[ie];
    double maxBDTsimple = vMaxBDTsimple[ie];
    double maxSSBsimple = vMaxSSBsimple[ie];
    double maxBDTfit = vMaxSSBfit[ie];
    double maxSSBfit = vMaxBDTfit[ie];
    double seffMax = vSeffMax[ie];
    double seffTot = vSeffTot[ie];

    if (ie < 3) {
      dMaxSum += dMax;
      sMaxSum += sMax;
      bMaxSum += bMax;
    }

    tl->DrawLatex(0.25, 0.44, Form("D/B/S = %d/%2.1f/%2.1f", dMax, bMax, sMax));
    tl->DrawLatex(0.25, 0.40, Form("S/#sqrt{S+B}(MC) = %4.3f (%4.3f)", maxSSB, maxSSBsimple));
    tl->DrawLatex(0.25, 0.36, Form("b_{max}(MC) = %4.3f (%4.3f)", maxBDT, maxBDTsimple));
    tl->DrawLatex(0.25, 0.32, Form("#epsilon_{BDT} = %4.3f", seffMax));
    tl->DrawLatex(0.25, 0.28, Form("#epsilon_{tot} = %6.5f", seffTot));
    tl->DrawLatex(0.25, 0.24, Form("I_{tot} = %6.5f", rocInt));
    tl->DrawLatex(0.25, 0.20, Form("I_{part} = %6.5f", rocInt2));

    ofstream TEX(texname.c_str(), ios::app);

    TEX << Form("\\vdef{s%s:ie%d:string}       {%s}", fname.c_str(), ie, fname.c_str()) << endl;
    TEX << Form("\\vdef{s%s:ie%d:ssb}          {%4.3f}", fname.c_str(), ie, maxSSB) << endl;
    TEX << Form("\\vdef{s%s:ie%d:maxbdt}       {%4.3f}", fname.c_str(), ie, maxBDT) << endl;
    TEX << Form("\\vdef{s%s:ie%d:ssbsimple}    {%4.3f}", fname.c_str(), ie, maxSSBsimple) << endl;
    TEX << Form("\\vdef{s%s:ie%d:maxbdtsimple} {%4.3f}", fname.c_str(), ie, maxBDTsimple) << endl;
    TEX << Form("\\vdef{s%s:ie%d:ssbfit}       {%4.3f}", fname.c_str(), ie, maxSSBfit) << endl;
    TEX << Form("\\vdef{s%s:ie%d:maxbdtfit}    {%4.3f}", fname.c_str(), ie, maxBDTfit) << endl;
    TEX << Form("\\vdef{s%s:ie%d:Smc}          {%4.3f}", fname.c_str(), ie, sMax) << endl;
    TEX << Form("\\vdef{s%s:ie%d:D}            {%d}", fname.c_str(), ie, dMax) << endl;
    TEX << Form("\\vdef{s%s:ie%d:Dhi}          {%d}", fname.c_str(), ie, dhiMax) << endl;
    TEX << Form("\\vdef{s%s:ie%d:B}            {%4.3f}", fname.c_str(), ie, bMax) << endl;
    TEX << Form("\\vdef{s%s:ie%d:ipart}        {%6.5f}", fname.c_str(), ie, rocInt2) << endl;
    TEX << Form("\\vdef{s%s:ie%d:itot}         {%6.5f}", fname.c_str(), ie, rocInt) << endl;
    TEX << Form("\\vdef{s%s:ie%d:epstot}       {%6.5f}", fname.c_str(), ie, seffTot) << endl;
    TEX << Form("\\vdef{s%s:ie%d:epsbdt}       {%6.5f}", fname.c_str(), ie, seffMax) << endl;
    if (ie == 3) {
      TEX << Form("\\vdef{s%s:BDTparameters}  {%s}", fname.c_str(), fBDTParameters.c_str()) << endl;
      TEX << Form("\\vdef{s%s:BDTvariables}  {%s}", fname.c_str(), fVariables.c_str()) << endl;
      TEX << Form("\\vdef{s%s:sum:Smc}     {%4.3f}", fname.c_str(), sMaxSum) << endl;
      TEX << Form("\\vdef{s%s:sum:D}       {%d}", fname.c_str(), dMaxSum) << endl;
      TEX << Form("\\vdef{s%s:sum:B}       {%4.3f}", fname.c_str(), bMaxSum) << endl;
      TEX << Form("\\vdef{s%s:sum:ssb}     {%4.3f}", fname.c_str(), sMaxSum/sqrt(sMaxSum+bMaxSum)) << endl;

      double ks0 = trainBDT0->KolmogorovTest(testBDT0);
      double ks1 = trainBDT1->KolmogorovTest(testBDT1);
      double ks2 = trainBDT2->KolmogorovTest(testBDT2);
      TEX << Form("\\vdef{s%s:ie0:ksBg}    {%4.3f}", fname.c_str(), ks0) << endl;
      TEX << Form("\\vdef{s%s:ie1:ksBg}    {%4.3f}", fname.c_str(), ks1) << endl;
      TEX << Form("\\vdef{s%s:ie2:ksBg}    {%4.3f}", fname.c_str(), ks2) << endl;

      ks0 = trainSgBDT0->KolmogorovTest(testSgBDT0);
      ks1 = trainSgBDT1->KolmogorovTest(testSgBDT1);
      ks2 = trainSgBDT2->KolmogorovTest(testSgBDT2);
      TEX << Form("\\vdef{s%s:ie0:ksSg}    {%4.3f}", fname.c_str(), ks0) << endl;
      TEX << Form("\\vdef{s%s:ie1:ksSg}    {%4.3f}", fname.c_str(), ks1) << endl;
      TEX << Form("\\vdef{s%s:ie2:ksSg}    {%4.3f}", fname.c_str(), ks2) << endl;
      TEX << Form("\\vdef{s%s:dMaxBDT}     {%4.3f}", fname.c_str(), dMaxBDT) << endl;
      TEX << Form("\\vdef{s%s:dMaxBDT3}    {%4.3f}", fname.c_str(), dMaxBDT3) << endl;
    }

    TEX.close();
    system(Form("/bin/cp %s plots", texname.c_str()));

    newLegend(0.22, 0.47, 0.55, 0.67);

    legg->SetTextSize(0.035);
    legg->AddEntry(grocsOp[ie], Form("operating point b > %4.3f", maxBDT), "p");
    if (ie < 3) {
      legg->AddEntry(grocs[ie],  Form("BDT %d", ie), "p");
    } else {
      legg->AddEntry(grocs[ie],  "combined", "p");
    }
    legg->Draw();

    c0->SaveAs(Form("plots/%s-roc-ie%d.pdf", fname.c_str(), ie));

    c0->Clear();
    hssb[ie]->Draw();
    h2ssb[ie]->Draw("same");
    tl->DrawLatex(0.2, 0.90, Form("event type %d", ie));
    tl->DrawLatex(0.2, 0.85, Form("SSB_{max} = %4.3f (%4.3f/%4.3f)", maxSSB, maxSSBsimple, maxSSBfit));
    tl->DrawLatex(0.2, 0.80, Form("BDT_{max} > %4.3f (%4.3f/%4.3f)", maxBDT, maxBDTsimple, maxBDTfit));
    tl->DrawLatex(0.2, 0.75, Form("ROC_{int} = %4.3f", rocInt));
    c0->SaveAs(Form("plots/%s-ssb-ie%d.pdf", fname.c_str(), ie));

    cout << "Write out SSB histograms" << endl;
    cout << "  maxSSB: " << maxSSB << " at BDT > " << maxBDT << endl;
    hssb[ie]->SetDirectory(f);

    f->cd();
    grocsOp[ie]->Write();
    grocs[ie]->Write();
  }


  f->Write();
  f->Close();

}



// ----------------------------------------------------------------------
void plotBdt::mvas(string fname) { //, HistType htype, Bool_t useTMVAStyle ) {
   // set style and remove existing canvas'
  TMVA::TMVAGlob::Initialize( kTRUE );

   // checks if file with name "fin" is already open, and if not opens one
   TString fin = Form("%s.root", fname.c_str());
   TFile* file = TFile::Open(fin.Data(), "UPDATE");
   TH1D *h2 = new TH1D("res_ks", "KS probabilities", 10, 0., 10.); h2->Sumw2();

   // define Canvas layout here!
   const Int_t width = 600;   // size of canvas

   // this defines how many canvases we need
   TCanvas *c = 0;

   // counter variables
   Int_t countCanvas = 0;

   // search for the right histograms in full list of keys
   TIter next(file->GetListOfKeys());
   TKey *key(0);
   while ((key = (TKey*)next())) {

      if (!TString(key->GetName()).BeginsWith("Method_")) continue;
      if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TDirectory")) continue;

      TString methodName;
      TMVA::TMVAGlob::GetMethodName(methodName,key);

      TDirectory* mDir = (TDirectory*)key->ReadObj();

      TIter keyIt(mDir->GetListOfKeys());
      TKey *titkey;
      while ((titkey = (TKey*)keyIt())) {

         if (!gROOT->GetClass(titkey->GetClassName())->InheritsFrom("TDirectory")) continue;

         TDirectory *titDir = (TDirectory *)titkey->ReadObj();
         TString methodTitle;
         TMVA::TMVAGlob::GetMethodTitle(methodTitle,titDir);

         cout << "--- Found directory for method: " << methodName << "::" << methodTitle << flush;
         TString hname = "MVA_" + methodTitle;
         TH1* sig = dynamic_cast<TH1*>(titDir->Get( hname + "_S" ));
         TH1* bgd = dynamic_cast<TH1*>(titDir->Get( hname + "_B" ));


         cout << " containing " << hname << "_S/_B" << endl;
         // chop off useless stuff
	 sig->SetTitle( Form("TMVA overtraining check for classifier: %s", methodTitle.Data()) );

         // create new canvas
         TString ctitle = Form("TMVA comparison %s",methodTitle.Data()) ;

         c = new TCanvas( Form("canvas%d", countCanvas+1), ctitle,
                          countCanvas*50+200, countCanvas*20, width, static_cast<int>(width*0.78) );

         // set the histogram style
         TMVA::TMVAGlob::SetSignalAndBackgroundStyle( sig, bgd );

         // normalise both signal and background
         TMVA::TMVAGlob::NormalizeHists( sig, bgd );

         // frame limits (choose judicuous x range)
         Float_t nrms = 10;
         cout << "--- Mean and RMS (S): " << sig->GetMean() << ", " << sig->GetRMS() << endl;
         cout << "--- Mean and RMS (B): " << bgd->GetMean() << ", " << bgd->GetRMS() << endl;
         Float_t xmin = TMath::Max( TMath::Min(sig->GetMean() - nrms*sig->GetRMS(),
                                               bgd->GetMean() - nrms*bgd->GetRMS() ),
                                    sig->GetXaxis()->GetXmin() );
         Float_t xmax = TMath::Min( TMath::Max(sig->GetMean() + nrms*sig->GetRMS(),
                                               bgd->GetMean() + nrms*bgd->GetRMS() ),
                                    sig->GetXaxis()->GetXmax() );
         Float_t ymin = 0;
         Float_t maxMult = 2.0;
         Float_t ymax = TMath::Max( sig->GetMaximum(), bgd->GetMaximum() )*maxMult;

	 xmin = -1.;
	 xmax = 1.;

         // build a frame
         Int_t nb = 500;
         TString hFrameName(TString("frame") + methodTitle);
         TObject *o = gROOT->FindObject(hFrameName);
         if(o) delete o;
         TH2F* frame = new TH2F( hFrameName, sig->GetTitle(),
                                 nb, xmin, xmax, nb, ymin, ymax );
         frame->GetXaxis()->SetTitle( methodTitle + " response"  );
         frame->GetYaxis()->SetTitle("(1/N) dN^{ }/^{ }dx");
         TMVA::TMVAGlob::SetFrameStyle( frame );

         // eventually: draw the frame
         frame->Draw();

         c->GetPad(0)->SetLeftMargin( 0.105 );
         frame->GetYaxis()->SetTitleOffset( 1.2 );


         // Draw legend
         TLegend *legend= new TLegend( c->GetLeftMargin(), 1 - c->GetTopMargin() - 0.12,
                                       c->GetLeftMargin() +  0.40, 1 - c->GetTopMargin() );
         legend->SetFillStyle( 1 );
	 cout << "legend->AddEntry(sig...)" << endl;
         legend->AddEntry(sig, TString("Signal")     + " (test sample)", "F");
         legend->AddEntry(bgd, TString("Background") + " (test sample)", "F");
         legend->SetBorderSize(1);
         legend->SetMargin(0.2);
         legend->Draw();

         // overlay signal and background histograms
         sig->Draw("samehist");
         bgd->Draw("samehist");

	 TH1* sigOv = 0;
	 TH1* bgdOv = 0;

	 TString ovname = hname += "_Train";
	 sigOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_S" ));
	 bgdOv = dynamic_cast<TH1*>(titDir->Get( ovname + "_B" ));

	 if (sigOv == 0 || bgdOv == 0) {
	   cout << "+++ Problem in \"mvas.C\": overtraining check histograms do not exist" << endl;
	 }
	 else {
	   cout << "--- Found comparison histograms for overtraining check" << endl;

	   TLegend *legend2= new TLegend( 1 - c->GetRightMargin() - 0.42, 1 - c->GetTopMargin() - 0.12,
					  1 - c->GetRightMargin(), 1 - c->GetTopMargin() );
	   legend2->SetFillStyle( 1 );
	   legend2->SetBorderSize(1);
	   legend2->AddEntry(sigOv,"Signal (training sample)","P");
	   legend2->AddEntry(bgdOv,"Background (training sample)","P");
	   legend2->SetMargin( 0.1 );
	   legend2->Draw();
	 }
	 // normalise both signal and background
	 TMVA::TMVAGlob::NormalizeHists( sigOv, bgdOv );

	 Int_t col = sig->GetLineColor();
	 sigOv->SetMarkerColor( col );
	 sigOv->SetMarkerSize( 0.7 );
	 sigOv->SetMarkerStyle( 20 );
	 sigOv->SetLineWidth( 1 );
	 sigOv->SetLineColor( col );
	 sigOv->Draw("e1same");

	 col = bgd->GetLineColor();
	 bgdOv->SetMarkerColor( col );
	 bgdOv->SetMarkerSize( 0.7 );
	 bgdOv->SetMarkerStyle( 20 );
	 bgdOv->SetLineWidth( 1 );
	 bgdOv->SetLineColor( col );
	 bgdOv->Draw("e1same");

	 ymax = TMath::Max(ymax,
			   TMath::Max( static_cast<Float_t>(sigOv->GetMaximum()), static_cast<Float_t>(bgdOv->GetMaximum()) )*maxMult);
	 frame->GetYaxis()->SetLimits( 0, ymax );

	 // for better visibility, plot thinner lines
	 sig->SetLineWidth( 1 );
	 bgd->SetLineWidth( 1 );

	 // perform K-S test
	 cout << "--- Perform Kolmogorov-Smirnov tests" << endl;
	 Double_t kolS = sig->KolmogorovTest( sigOv );
	 Double_t kolB = bgd->KolmogorovTest( bgdOv );
	 cout << "--- Goodness of signal (background) consistency: " << kolS << " (" << kolB << ")" << endl;

	 TString probatext = Form( "Kolmogorov-Smirnov test: signal (background) probability = %5.3g (%5.3g)", kolS, kolB );
	 TText* tt = new TText( 0.12, 0.91, probatext );
	 tt->SetNDC(); tt->SetTextSize( 0.032 ); tt->AppendPad();

	 h2->SetBinContent(1, kolS); h2->GetXaxis()->SetBinLabel(1, "KS(sg)");
	 h2->SetBinContent(2, kolB); h2->GetXaxis()->SetBinLabel(2, "KS(bg)");
	 writeOut(file, h2);
	 // h2->Write();



         // redraw axes
         frame->Draw("sameaxis");

	 //          // text for overflows
	 //          Int_t    nbin = sig->GetNbinsX();
	 //          Double_t dxu  = sig->GetBinWidth(0);
	 //          Double_t dxo  = sig->GetBinWidth(nbin+1);
	 //          TString uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%",
	 //                                 sig->GetBinContent(0)*dxu*100, bgd->GetBinContent(0)*dxu*100,
	 //                                 sig->GetBinContent(nbin+1)*dxo*100, bgd->GetBinContent(nbin+1)*dxo*100 );
	 TText* t = new TText( 0.975, 0.115, fname.c_str());
	 t->SetNDC();
	 t->SetTextSize( 0.030 );
	 t->SetTextAngle( 90 );
	 t->AppendPad();

         // update canvas
         c->Update();

         // save canvas to file
	 c->SaveAs(Form("plots/%s-overtrain0.pdf", fname.c_str()));
	 frame->GetYaxis()->SetLimits(1.e-5, 5.e1);
	 c->SetLogy(1);
	 c->SaveAs(Form("plots/%s-overtrain1.pdf", fname.c_str()));

         countCanvas++;
      }
      cout << "";
   }
   cout << endl;
   file->Write();
   file->Close();
}


// ----------------------------------------------------------------------
TH1D* plotBdt::getRanking(string fname, string prefix, string type) {
  TH1D *h1 = new TH1D(Form("rank_%s_%s", type.c_str(), prefix.c_str()), Form("rank_%s", prefix.c_str()), 100, 0., 100.);
  // -- read in variable ranking from logfile
  vector<string> allLines;
  char  buffer[2000];
  cout << "getRanking: open file " << Form("%s.log", fname.c_str()) << endl;
  ifstream is(Form("%s.log", fname.c_str()));
  while (is.getline(buffer, 2000, '\n')) allLines.push_back(string(buffer));
  string::size_type m1, m2;
  string varn, vars, varp;
  int bail(0), istart(0);
  string after;

  if (type == "events0") {
    after = "Apply on events0, train on events1, test on events2";
  } else if (type == "events1") {
    after = "Apply on events1, train on events2, test on events0";
  } else if (type == "events2") {
    after = "Apply on events2, train on events0, test on events1";
  }

  for (unsigned int i = 0; i < allLines.size(); ++i) {
    if (string::npos != allLines[i].find(after)) {
      istart = i;
      break;
    }
  }

  cout << "start after: " << after << " at line " << istart << endl;
  for (unsigned int i = istart; i < allLines.size(); ++i) {
    // -- method unspecific classification
    if ((string::npos != allLines[i].find(Form("--- %s", prefix.c_str())))
	&& (string::npos != allLines[i].find(": Rank : Variable "))
	) {
      bail = 0;
      for (unsigned int j = i+2; j < i+100; ++j) {
	if (string::npos != allLines[j].find(": ---------------------------------")) {
	  bail = 1;
	  cout << "  -> breaking out " << endl;
	  break;
	}

	m1 = allLines[j].find(":");
	m2 = allLines[j].find(":", m1+1);
	varn = allLines[j].substr(m1+2, m2-m1-2);
	m1 = m2;
	m2 = allLines[j].find(":", m1+1);
	vars = allLines[j].substr(m1+2, m2-m1-2);
	m1 = m2;
	m2 = allLines[j].find(":", m1+1);
	varp = allLines[j].substr(m1+2, m2-m1-2);
	cout << varn << "-> " << vars << " -> " << varp << endl;
	//	cout << allLines[j] << endl;
	int ibin = atoi(varn.c_str());
	if (0 == ibin) cout << "===> bin == 0: " << varn << endl;
	h1->GetXaxis()->SetBinLabel(ibin, vars.c_str());
	h1->SetBinContent(ibin, atof(varp.c_str()));
      }
      if (1 == bail) break;
    }
  }

  return h1;
}
