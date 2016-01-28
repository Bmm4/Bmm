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
#include "TSystem.h"
#include "TKey.h"

#include "bmmReader.hh"
#include "genAnalysis.hh"

#include "common/util.hh"


using namespace std;

void skimEvents(TChain *);


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: bin/runBmm [-r genAnalysis] -f test.root
// %%        bin/runBmm -c chains/bg-test -D root
// %%        bin/runBmm -c chains/bg-test -s skim
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "Running under process ID  " << processID << endl;

  string progName  = argv[0]; 
  string writeName, fileName, jsonName;
  int file(0), json(0);
  int dirspec(0);
  int nevents(-1), start(-1);
  int randomSeed(processID);
  int verbose(-99); 
  int blind(1); 
  int isMC(0);
  int year(0); 

  // Change the MaxTreeSize to 100 GB (default since root v5.26)
  TTree::SetMaxTreeSize(100000000000ll); // 100 GB

  // -- Some defaults
  string dirBase("./");               // this could point to "/home/ursl/data/root/."
  string dirName("."); dirspec = 0;   // and this to, e.g. "bmm", "bee", "bem", ...
  string cutFile("tree.defaults.cuts");

  string treeName("T1");
  string evtClassName("TAna01Event");

  string readerName("bmmReader"); 
  TString histfile("");

  // -- command line arguments
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i],"-h")) {
	cout << "List of arguments:" << endl;
	cout << "-b {0,1}      run blind?" << endl; 
	cout << "-c filename   chain definition file" << endl;
	cout << "-C filename   file with cuts" << endl;
	cout << "-D path       where to put the output" << endl;
	cout << "-f filename   single file instead of chain" << endl;
	cout << "-m            use MC information" << endl;
	cout << "-n integer    number of events to run on" << endl;
	cout << "-r class name which tree reader class to run" << endl;
	cout << "-s number     seed for random number generator" << endl;
	cout << "-S start     starting event number" << endl;
	cout << "-o filename   set output file" << endl;
	cout << "-v level      set verbosity level" << endl;
	cout << "-y year       set year" << endl;
	cout << "-h            prints this message and exits" << endl;
	return 0;
    }
    if (!strcmp(argv[i],"-b"))  {blind      = atoi(argv[++i]); }                 // run blind?
    if (!strcmp(argv[i],"-c"))  {fileName   = string(argv[++i]); file = 0; }     // file with chain definition
    if (!strcmp(argv[i],"-C"))  {cutFile    = string(argv[++i]);           }     // file with cuts
    if (!strcmp(argv[i],"-D"))  {dirName    = string(argv[++i]);  dirspec = 1; } // where to put the output
    if (!strcmp(argv[i],"-f"))  {fileName   = string(argv[++i]); file = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-j"))  {jsonName   = string(argv[++i]); json = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-m"))  {isMC       = 1; }                               // use MC information?
    if (!strcmp(argv[i],"-n"))  {nevents    = atoi(argv[++i]); }                 // number of events to run 
    if (!strcmp(argv[i],"-r"))  {readerName = string(argv[++i]); }               // which tree reader class to run
    if (!strcmp(argv[i],"-s"))  {randomSeed = atoi(argv[++i]); }                 // set seed for random gen.
    if (!strcmp(argv[i],"-S"))  {start = atoi(argv[++i]); }                      // set start event number
    if (!strcmp(argv[i],"-o"))  {histfile   = TString(argv[++i]); }              // set output file
    if (!strcmp(argv[i],"-v"))  {verbose    = atoi(argv[++i]); }                 // set verbosity level
    if (!strcmp(argv[i],"-y"))  {year       = atoi(argv[++i]); }                 // set year
  }


  // -- Prepare histfilename variation with (part of) cut file name
  TString fn(cutFile);
  fn.ReplaceAll("cuts/", "");
  fn.ReplaceAll(".cuts", "");
  fn.ReplaceAll("tree", "");

  // -- Determine filename for output histograms and 'final' small/reduced tree
  TString meta = fileName;
  if(histfile == "") {
    TString  barefile(fileName), chainFile, meta;
    if (file == 0) {
      // -- input from chain
      if (barefile.Contains("chains/")) {
	barefile.ReplaceAll("chains/", "");
	histfile = barefile + "." + fn + ".root";
	if (dirspec) {
	  if (dirName[0] == '/') {
	    histfile = dirName + "/" + histfile;
	  } else {
	    histfile = dirBase + "/" + dirName + "/" + histfile;
	  }
	}
      } else {
	histfile =  barefile + "." + fn + ".root";
	if (dirspec) {
	  if (dirName[0] == '/') {
	    histfile = dirName + "/" + histfile;
	  } else {
	    histfile = dirBase + "/" + dirName + "/" + histfile;
	  }
	}
      }
      // -- The following lines strip everything from the string up to and including the last '/'
      int fl = barefile.Last('/');
      TString bla(barefile);
      bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
      histfile =  bla + "." + fn + ".root";
      if (dirspec) {
	histfile = dirBase + "/" + dirName + "/" + histfile;
      }
    }  else if (file == 1) {
      // -- single file input
      // -- The following lines strip everything from the string up to and including the last '/'
      int fl = barefile.Last('/');
      TString bla(barefile);
      bla.Replace(0, fl+1, ' '); bla.Strip(TString::kLeading, ' ');  bla.Remove(0,1);
      histfile =  bla;
      histfile.ReplaceAll(".root", "");
      histfile +=  "." + fn + ".root";
      if (dirspec) {
	if (dirName[0] == '/') {
	  histfile = dirName + "/" + histfile;
	} else {
	  histfile = dirBase + "/" + dirName + "/" + histfile;
	}
      }
    }
  }
  string shistfile = histfile.Data();
  replaceAll(shistfile, "..", ".");
  histfile = shistfile.c_str(); 
  cout << "Opening " << histfile.Data() << " for output histograms" << endl;
  cout << "Opening " << fileName.c_str() << " for input" << endl;


  // -- Set up chain
  TChain *chain = new TChain(TString(treeName));
  cout << "Chaining ... " << treeName << endl;
  char pName[2000]; 
  int nentries; 
  if (file == 0) {
    // -- non-trivial chain input
    ifstream is(meta);  
    while(meta.ReadLine(is) && (!meta.IsNull())){ 
      nentries = -1;
      if (meta.Data()[0] == '#') continue; 
      sscanf(meta.Data(), "%s %d", pName, &nentries); 
      if (nentries > -1) {
        cout << pName << " -> " << nentries << " entries" << endl; 
        chain->Add(pName, nentries); 
      } else {
        cout << meta << endl;
        chain->Add(meta); 
      }
    }
    is.close();
  }
  else if (file == 1) {
    // -- single file input
    cout << fileName << endl;
    chain->Add(TString(fileName));
  }

  // -- Now instantiate the tree-analysis class object, initialize, and run it ...
  //treeReader01 *a = new bmm2Reader(chain, TString(evtClassName));  
  treeReader01 *a = NULL;
  if ("bmmReader" == readerName) {
    cout << "instantiating bmmReader" << endl;
    a = new bmmReader(chain, TString(evtClassName));
  } else if ("genAnalysis" == readerName) {
    cout << "instantiating genAnalysis" << endl;
    a = new genAnalysis(chain, TString(evtClassName));
  } else if ("skim" == readerName) {
    skimEvents(chain);
  }
  
  
  if (a) {
    a->setYear(year); 
    if (verbose > -99) a->setVerbosity(verbose); 
    a->openHistFile(histfile); 

    if (isMC) {
      a->setMC(1);
      blind = 0; 
    } else {
      a->setMC(0); 
    }

    cout << "blind? " << blind << endl;
    if (1 == blind) a->runBlind();

    a->readCuts(cutFile.c_str(), 1);
    a->bookHist();
    if (json) {
      a->setJSONFile(jsonName.c_str()); 
      a->forceJSON();
    }


    a->startAnalysis(); 
    a->loop(nevents, start);
    a->endAnalysis();
    a->closeHistFile(); 
  } 

  delete a; // so we can dump some information in the destructor
  
  return 0;
}



// ----------------------------------------------------------------------
void skimEvents(TChain *chain) {
  int oldRun(-1), run(-1), evt(-1), ls(-1); 
  vector<pair<int, int> > events; 
  ifstream INS; 
  string sline, bla; 
  INS.open("skim.events");
  while (1) {
    INS >> bla >> run >> bla >> evt; 
    if (INS.eof()) break;
    events.push_back(make_pair(run, evt)); 
  }
      
  INS.close();

  cout << "skim event list" << endl;
  cout << "----------------------------------------------------------------------" << endl;
  if (events.size() < 1) {
    cout << "no file ./skim.events found! exit(0)" << endl;
    exit(0); 
  } else {
    for (unsigned int i = 0; i < events.size(); ++i) {
      cout << "run: " << events[i].first << " event: " << events[i].second << endl;
    }
  }
  cout << "----------------------------------------------------------------------" << endl;

  TAna01Event *pEvt = new TAna01Event(0);
  chain->SetBranchAddress("TAna01Event", &pEvt);

  TFile *newfile = new TFile("skimevents.root", "recreate");
  TTree *newtree = chain->CloneTree(0);
  
  int nb(0); 
  cout << "chain entries: " << chain->GetEntries() << endl;
  for (int jEvent = 0; jEvent < chain->GetEntries(); ++jEvent) {
    pEvt->Clear();                        
    nb += chain->GetEvent(jEvent);          
    
    evt = static_cast<long int>(pEvt->fEventNumber);		
    run = static_cast<long int>(pEvt->fRunNumber);
    ls = pEvt->fLumiSection;
    
    if (run != oldRun) {
      cout << "new run: " << run << " (event " << jEvent << " in chain)" << endl;
      oldRun = run;

      // -- dump all pd histograms
      TKey *key(0);   
      TH1D *h(0);
      string sname; 
      string::size_type m1, m2;
      string pd, lpd, hk; 
      
      TIter next(chain->GetFile()->GetListOfKeys());
      vector<string> triggers;
      while ((key = (TKey*)next())) {
	sname = key->GetName();
	
	if (string::npos == sname.find("triggers_")
	    && string::npos == sname.find("pd_run")
	    && string::npos == sname.find("_run")) continue;
	
	h = (TH1D*)chain->GetFile()->Get(sname.c_str()); 
	cout << sname << endl;
	h->SetDirectory(newfile);
	h->Write();
	delete h; 
      }

    }

    for (unsigned int i = 0; i < events.size(); ++i) {
      if ((run == events[i].first) && (evt == events[i].second)) {
	cout << "skimming run = " << run <<  " ls = " << ls << " evt = " << evt
	     << " file: " << chain->GetFile()->GetName()
	     << endl;
	newtree->Fill();
	pEvt->Clear();
      }
    }
  }

  newtree->Print();
  newtree->AutoSave();

  delete newfile;

}
