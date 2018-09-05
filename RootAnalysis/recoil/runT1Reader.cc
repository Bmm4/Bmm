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

#include "t1Reader.hh"
#include "genAnalysis.hh"

#include "common/util.hh"


using namespace std;

void skimEvents(TChain *);
void dumpTriggers(string);


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% Usage: bin/runT1Reader [-r genAnalysis] -f test.root
// %%        bin/runT1Reader -c chains/bg-test -D root
// %%        bin/runT1Reader -c chains/bg-test -s skim
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char *argv[]) {

  int processID = gSystem->GetPid();
  cout << "Running under process ID  " << processID << endl;

  string progName  = argv[0];
  string writeName, fileName;
  int file(0);
  int dirspec(0);
  int nevents(-1), start(-1);
  int randomSeed(processID);
  int verbose(-99);
  int blind(1);
  int isMC(0);
  int year(0);
  int json(0);
  bool testType(0);
  string era("");

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
	cout << "-e era        BF, GH, ...?" << endl;
	cout << "-f filename   single file instead of chain" << endl;
	cout << "-m            use MC information" << endl;
	cout << "-n integer    number of events to run on" << endl;
	cout << "-r class name which tree reader class to run" << endl;
	cout << "-s number     seed for random number generator" << endl;
	cout << "-S start      starting event number" << endl;
	cout << "-t            test for cand type before calling corresponing candAna* (careful: screws mon*!)" << endl;
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
    if (!strcmp(argv[i],"-e"))  {era        = string(argv[++i]);  }              // where to put the output
    if (!strcmp(argv[i],"-f"))  {fileName   = string(argv[++i]); file = 1; }     // single file instead of chain
    if (!strcmp(argv[i],"-j"))  {json       = 1; }                               // ignore JSON status
    if (!strcmp(argv[i],"-m"))  {isMC       = 1; }                               // use MC information?
    if (!strcmp(argv[i],"-n"))  {nevents    = atoi(argv[++i]); }                 // number of events to run
    if (!strcmp(argv[i],"-o"))  {histfile   = TString(argv[++i]); }              // set output file
    if (!strcmp(argv[i],"-r"))  {readerName = string(argv[++i]); }               // which tree reader class to run
    if (!strcmp(argv[i],"-s"))  {randomSeed = atoi(argv[++i]); }                 // set seed for random gen.
    if (!strcmp(argv[i],"-S"))  {start = atoi(argv[++i]); }                      // set start event number
    if (!strcmp(argv[i],"-t"))  {testType   = true; }                            // check for type before running candAna* class
    if (!strcmp(argv[i],"-v"))  {verbose    = atoi(argv[++i]); }                 // set verbosity level
    if (!strcmp(argv[i],"-y"))  {year       = atoi(argv[++i]); }                 // set year
  }


  // -- Prepare histfilename variation with (part of) cut file name
  TString fn(cutFile);
  fn.ReplaceAll("cuts/", "");
  fn.ReplaceAll(".cuts", "");
  fn.ReplaceAll("tree", "");

  if ("trigger" == readerName) {
    dumpTriggers(fileName);
    return 0;
  }

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
  t1Reader *a = NULL;
  if ("t1Reader" == readerName) {
    cout << "instantiating t1Reader" << endl;
    a = new t1Reader(chain, TString(evtClassName));
  } else if ("genAnalysis" == readerName) {
    cout << "instantiating genAnalysis" << endl;
    a = new genAnalysis(chain, TString(evtClassName));
  } else if ("skim" == readerName) {
    skimEvents(chain);
  }


  if (a) {
    a->setYear(year);
    if (era != "") a->setEra(era);
    a->setCheckCandTypes(testType);
    if (verbose > -99) a->setVerbosity(verbose);
    a->openHistFile(histfile);

    if (json) {
      a->ignoreJSON();
    }
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


// ----------------------------------------------------------------------
void dumpTriggers(string fileName) {
  vector<string> tnames;
  tnames.push_back("Bs");
  tnames.push_back("Jpsi_Displaced");

  vector<string> files;
  ifstream INS;
  string sline;
  INS.open(fileName.c_str());
  while (getline(INS, sline)) {
    string::size_type m1 = sline.rfind(".root ");
    sline = sline.substr(0, m1+6);
    files.push_back(sline);
  }


  map<string, pair<int, int> > tranges;

  cout << "----------------------------------------------------------------------" << endl;
  cout << "dumpTriggers: " << endl;
  for (unsigned int i = 0; i < tnames.size(); ++i) {
    cout << "trigger name part: " << tnames[i] << endl;
  }
  cout << "----------------------------------------------------------------------" << endl;

  TFile *f(0);
  TH1D *h1(0);
  TKey *key(0);
  string st("");
  int run(0);
  for (unsigned int i = 0; i < files.size(); ++i) {
    //  for (unsigned int i = 0; i < 20; ++i) {
    cout << files[i] << endl;
    f = TFile::Open(files[i].c_str());
    if (f && f->IsOpen()) {
      TIter next(f->GetListOfKeys());
      while ((key = (TKey*)next())) {
	if (!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1")) continue;
	if (1 != key->GetCycle()) {
	  continue;
	}
	TString skey = TString(key->GetName());
	if (skey.Contains("triggers_MuOnia") || skey.Contains("triggers_Charmonium")) {
	  h1 = (TH1D*)((TH1D*)f->Get(skey));
	  if (0 == h1) continue;
	  st = h1->GetName();
	  st = st.substr(st.rfind("_run")+4);
	  run = atoi(st.c_str());
	  for (int ib = 0; ib <= h1->GetNbinsX(); ++ib) {
	    st = h1->GetXaxis()->GetBinLabel(ib);
	    for (unsigned int ip = 0; ip < tnames.size(); ++ip) {
	      if (string::npos != st.find(tnames[ip])) {
		if (tranges.count(st) != 0) {
		  int rmin = tranges[st].first;
		  int rmax = tranges[st].second;
		  if (run < rmin) tranges[st].first = run;
		  if (run > rmax) tranges[st].second = run;
		} else {
		  tranges.insert(make_pair(st, make_pair(run, run)));
		}
	      }
	    }
	  }
	}
      }
      f->Close();
    }
  }

  for (map<string, pair<int, int> >::iterator it = tranges.begin(); it != tranges.end(); ++it) {
    cout << it->first << ": " << it->second.first << " .. " << it->second.second << endl;
  }
}
