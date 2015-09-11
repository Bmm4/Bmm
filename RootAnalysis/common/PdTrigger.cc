#include "PdTrigger.hh"

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include <TString.h>
#include <TSystem.h>

#include "util.hh"

using namespace std;

// ----------------------------------------------------------------------
PdTrigger::PdTrigger(string basedir, int verbose) {
 
  vector<string> pdlist = glob(basedir, ".triggers"); 
  vector<string> trigger; 
  for (unsigned int i = 0; i < pdlist.size(); ++i) {
    trigger.clear(); 
    readFile(basedir + string("/") + pdlist[i], trigger); 
    string pdname = pdlist[i];
    replaceAll(pdname, ".triggers", ""); 
    fPdTrigger.insert(make_pair(pdname, trigger)); 
  }

  if (verbose) print();
  fVerbose = verbose;
}


// ----------------------------------------------------------------------
PdTrigger::~PdTrigger() {
}


// ----------------------------------------------------------------------
void PdTrigger::print() {
  map<string, vector<string> >::iterator it = fPdTrigger.begin(); 
  for (; it != fPdTrigger.end(); ++it) {
    cout << it->first << endl;
    for (unsigned int i = 0; i < it->second.size(); ++i) {
      cout << "  " << it->second[i] << endl;
    }
  }

}

// ----------------------------------------------------------------------
bool PdTrigger::triggerInPd(std::string pdname, std::string triggername) {
  size_t pos = triggername.rfind("_v");
  if (pos != std::string::npos) {
    triggername = triggername.substr(0, pos);
    if(fVerbose>1) cout << "using truncated triggername: " << triggername << endl;
  }

  vector<string> triggers = fPdTrigger[pdname];
  for (unsigned int i = 0; i < triggers.size(); ++i) {
    if (triggername == triggers[i]) return true;
  }
  return false;
}


// ----------------------------------------------------------------------
void PdTrigger::readFile(string filename, vector<string> &lines) {
  cout << "    readFile " << filename << endl;
  char  buffer[200];
  ifstream is(filename.c_str());
  if (!is) {
    exit(1);
  }
  while (is.getline(buffer, 200, '\n')) {
    if (buffer[0] == '#') continue;
    lines.push_back(string(buffer));
  }

  is.close();
}


// ----------------------------------------------------------------------
vector<string> PdTrigger::glob(string basedir, string basename) {
  cout << "Looking in " << basedir << " for " << basename << endl;
  vector<string> lof; 
  TString fname;
  const char *file;
  TSystem *lunix = gSystem; //new TUnixSystem();
  void *pDir = lunix->OpenDirectory(basedir.c_str());
  while ((file = lunix->GetDirEntry(pDir))) {
    fname = file;
    if (fname.Contains(basename.c_str())) {
      lof.push_back(string(fname));
    }
  }  
  return lof; 
}
