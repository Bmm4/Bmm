#include "FWCore/Framework/interface/MakerMacros.h"
#include "Bmm/CmsswAnalysis/interface/HFTree.hh"

#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"
#include "Bmm/RootAnalysis/rootio/TAnaVertex.hh"
#include "Bmm/RootAnalysis/rootio/TAnaMuon.hh"

// -- Yikes!
TAna01Event  *gHFEvent;
TFile        *gHFFile;

using namespace::std;

// ----------------------------------------------------------------------
HFTree::HFTree(const edm::ParameterSet& iConfig) :
  fRequireCand(iConfig.getUntrackedParameter<bool>("requireCand", true)),
  fFullGenBlock(iConfig.getUntrackedParameter<bool>("fullGenBlock", false)),
  fFileName(iConfig.getUntrackedParameter<string>("fileName", string("hfa.root"))),
  fTreeName(iConfig.getUntrackedParameter<string>("treeName", string("T1"))),
  fVerbose(iConfig.getUntrackedParameter<int>("verbose", 1)),
  fPrintFrequency(iConfig.getUntrackedParameter<int>("printFrequency", 1000)) {
  using namespace std;
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFTree constructor" << endl;
  cout << "---  verbose:                         " << fVerbose << endl;
  cout << "---  printFrequency:                  " << fPrintFrequency << endl;
  cout << "---  fileName:                        " << fFileName << endl;
  cout << "---  treeName:                        " << fTreeName << endl;
  cout << "---  requireCand:                     " << (fRequireCand?"true":"false") << endl;
  cout << "---  fullGenBlock:                    " << (fFullGenBlock?"true":"false") << endl;
  cout << "----------------------------------------------------------------------" << endl;
  fFile = TFile::Open(fFileName.c_str(), "RECREATE");
  fTree = new TTree(fTreeName.c_str(), "CMSSW HF tree");
  fEvent = new TAna01Event(0);
  TAna01Event::Class()->SetCanSplit(1);
  fTree->Branch("TAna01Event", "TAna01Event", &fEvent, 64000, 99);

  //   fTree->Branch("TAna01Event", "TAna01Event", &fEvent, 64000, 1);
  //  fTree->Branch("TAna01Event", "TAna01Event", &fEvent, 256000/8, 0);
  //  fTree->Bronch("TAna01Event", "TAna01Event", &fEvent, 256000/8, 1);

  fH1 = new TH1D("h1", "h1", 20, 0., 20.);
  fH1->SetDirectory(fFile);

  gHFEvent = fEvent;
  gHFFile  = fFile;

  fEventCounter = -1;
}


// ----------------------------------------------------------------------
HFTree::~HFTree() {

  // -- Save output
  fFile->cd();
  //  fTree->Write();
  fFile->Write();
  fFile->Close();
  delete fFile;
}


// ----------------------------------------------------------------------
void HFTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  ++fEventCounter;
  fH1->Fill(10.);
  if (fVerbose > 0) {
    if (fEventCounter%fPrintFrequency == 0) {
      pid_t pid = getpid();
      char line[100];
      sprintf(line, "ps -F %i", pid);
      cout << "==>HFTree: analyze() in event #" << fEventCounter
	   << "  run: " << gHFEvent->fRunNumber
	   << " event: "  << gHFEvent->fEventNumber
	   << endl;
      system(line);
    }
  }

  if (false == fFullGenBlock) {
    gHFEvent->clearGenBlock();
  }


  if (fRequireCand){
    if (gHFEvent->nCands() > 0) {
      if (fVerbose > 1) {
	cout << "HFTree> filling tree for run: " << gHFEvent->fRunNumber
	     << " event: "  << gHFEvent->fEventNumber ;

	if (fVerbose > 2) {
	  cout << " Cand: " << gHFEvent->nCands() << endl;

	  for (int i = 0; i < gHFEvent->nCands(); ++i) {
	    cout << gHFEvent->getCand(i)->fType << " ";
	  }
	}
	cout << endl;
      }
      fTree->Fill();
      fH1->Fill(11.);
    }
  } else {
    if (fVerbose > 1) {
      cout << "HFTree> filling tree for run: " << gHFEvent->fRunNumber
	   << " event: "  << gHFEvent->fEventNumber ;

      if (fVerbose > 2) {
	cout << " GENT Cand: " << gHFEvent->nGenT() << endl;
	for (int i = 0; i < gHFEvent->nGenT(); ++i) {
	  cout << gHFEvent->getGenCand(i)->fID << " ";
	}
	cout << endl;

	cout << " GEN Cand: " << gHFEvent->nGenCands() << endl;
	for (int i = 0; i < gHFEvent->nGenCands(); ++i) {
	  cout << gHFEvent->getGenCand(i)->fID << " ";
	}
	cout << endl;
      }
    }
    fTree->Fill();
    fH1->Fill(2.);
  }

  gHFEvent->Clear();
}

// ------------ method called once each job just before starting event loop  ------------
void  HFTree::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void  HFTree::endJob() {

  pid_t pid = getpid();
  char line[100];
  sprintf(line, "ps -F %i", pid);
  cout << "==>HFTree: endJob():" << endl;
  system(line);

}

//define this as a plug-in
DEFINE_FWK_MODULE(HFTree);
