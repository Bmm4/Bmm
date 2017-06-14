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
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/tmvaglob.h"

//2011
//#define LUMISCALE 2.01e-4

//2012  12/7405 = 0.00162
//#define LUMISCALE 0.00162

ClassImp(tmva1)

using namespace std;

// ----------------------------------------------------------------------
// --
// -- USAGE: a.makeAll(0, 1); > TMVA-0.log
// --
// ----------------------------------------------------------------------

tmva1::tmva1(int year, string vars, string pars) {

  fYear          = 2016;
  fVariables     = vars;
  fBDTParameters = pars;

  cout << "tmva1 hello: setup for year = " << year << " with variables: " << vars << endl;

  fVariables = vars;
  fBDTParameters = "";

  legg = 0;
  legge = 0;
  tl = new TLatex();
  tl->SetTextFont(42);
  tl->SetTextSize(0.03);
  tl->SetNDC(kTRUE);

  fYear = year;

  if (year == 2011) {
    fLumiScale = 3.1e-4; // 4.9/16000
    fInputFiles.sname = "/scratch/ursl/bmm4/v05/";
    fInputFiles.dname = "/scratch/ursl/bmm4/v05/";
  } else if (year == 2012) {
    fLumiScale = 2.8e-4; // 20/714000
    fInputFiles.sname = "/scratch/ursl/bmm4/v05/";
    fInputFiles.dname = "/scratch/ursl/bmm4/v05/";
  } else if (year == 2016) {
    fLumiScale = 2.8e-4; // 20/714000
    fInputFiles.sname = "/scratch/ursl/bmm4/v06/bmm-mc-RunIISpring16DR80-BsToMuMu_BMuonFilter-v06.root";
    fInputFiles.dname = "/scratch/ursl/bmm4/v06/bmm-data-bmmCharmonium2016-v06.root";
  }

//   // -- BDT setup 108/109
//   fBdtSetup.NTrees = 800;
//   fBdtSetup.nEventsMin = 50;
//   fBdtSetup.MaxDepth = 2;
//   //  fBdtSetup.MaxDepth = 3;
//   fBdtSetup.nCuts = 20;
//   fBdtSetup.AdaBoostBeta = 1.0;
//   fBdtSetup.NNodesMax = 5;
//   //  fBdtSetup.NNodesMax = 20;

  // -- TMVA default
  fBdtSetup.NTrees = 800;
  fBdtSetup.nEventsMin = 100;
  fBdtSetup.MaxDepth = 3;
  fBdtSetup.nCuts = 20;
  fBdtSetup.AdaBoostBeta = 1.0;
  fBdtSetup.NNodesMax = 100000;

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
  // createInputFile(filename);

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
    }
  }

  string logname = Form("TMVA-%d", offset);
  string oname("");
  oname = Form("TMVA-%d", offset);
  make(offset, filename, 0);
  make(offset, filename, 1);
  make(offset, filename, 2);
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
void tmva1::train(string oname, string filename, int nsg, int nbg) {
   // This loads the library
   TMVA::Tools::Instance();

   (TMVA::gConfig().GetVariablePlotting()).fNbins1D = 40;
   (TMVA::gConfig().GetVariablePlotting()).fNbinsMVAoutput = 40;

   // -- Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName(Form("%s.root", oname.c_str()));
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TH1D *hSetup = new TH1D("hSetup", "hSetup", 100, 0., 100.);
   int i(0);
   i =  1; hSetup->SetBinContent(i, fTrainAntiMuon?1:0); hSetup->GetXaxis()->SetBinLabel(i, "antimuon");
   i =  3; hSetup->SetBinContent(i, fRsigma); hSetup->GetXaxis()->SetBinLabel(i, "rsigma");
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
   cout << "==> oname: " << oname << " antimuon: " << fTrainAntiMuon <<  endl;

   string optstring = "V:!Silent:!Color:!DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification";
   optstring        = "V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification";
   cout << "==> Factory: " << optstring << endl;
   TMVA::Factory *factory = new TMVA::Factory(Form("%s", oname.c_str()), outputFile,  optstring.c_str());

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
     factory->AddVariable(vVar[i].c_str(), 'F');
     //      }
   }

   factory->AddSpectator("m",  "mass", "GeV", 'F' );

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
     cout << "trainBg = " << trainBg << endl;
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

   // optstring=Form("nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=None:V");
   // optstring=Form("nTrain_Signal=%d:nTest_Signal=%d:nTrain_Background=%d:nTest_Background=%d:SplitMode=random:SplitSeed=%d:NormMode=None:V",
   // nSgTrain, nSgTest, nBgTrain, nBgTest, seed);

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
     optstring = "!H:V" + fBDTParameters;
   }

   // -- Josh's (modified) proposal
   //    optstring = "!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1";
   //    optstring += ":UseBaggedGrad=F:nCuts=200:MaxDepth=3:NNodesMax=100000:UseYesNoLeaf=F:nEventsMin=1000:";

   cout << "==> BookMethod: " << optstring << endl;
   factory->BookMethod( TMVA::Types::kBDT, "BDT", optstring);

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
void tmva1::createInputFile(string filename, int randomSeed) {
  TFile *sinput = TFile::Open(fInputFiles.sname.c_str());
  TFile *dinput = TFile::Open(fInputFiles.dname.c_str());

  TCut sgcut = preselection().c_str();

  cout << "new: " << endl;
  cout << sgcut << endl;

  TCut masscut = "m>4.9&&m<5.9";
  TCut massbg  = "!(5.2<m&&m<5.45)";
  TCut muonid  = "gmugmid";

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

    // -- combined version:
    // -- signal
    sdir = Form("signalAll%s", type.c_str());
    chanCut = "(0 == chan) || (1 == chan)";
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
