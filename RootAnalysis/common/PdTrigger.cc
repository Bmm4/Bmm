#include "PdTrigger.hh"

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>

#include "util.hh"

using namespace std;

// ----------------------------------------------------------------------
PdTrigger::PdTrigger(string directory, int verbose): fHLTKey("nada"), fVerbose(verbose) {
  //  readPdTriggers(Form("%s/pdtriggers.txt", directory.c_str()));
  mkPdTriggersNoV();
  if (verbose) print();
}


// ----------------------------------------------------------------------
PdTrigger::~PdTrigger() {
}


// ----------------------------------------------------------------------
void PdTrigger::print() {
  for (map<string, vector<string> >::iterator it = fPdTriggers.begin(); it != fPdTriggers.end(); ++it) {
    cout << it->first << endl;
    for (unsigned int i = 0; i < it->second.size(); ++i) {
      cout << "  " << it->second[i] << endl;
    }
  }

}


// ----------------------------------------------------------------------
string PdTrigger::getHLTKey(int run, TFile *f) {
  TH1D *h(0);
  h = (TH1D*)f->Get(Form("pd_run%d", run)); 
  string key("nada");
  if (h) {
    key = h->GetTitle(); 
  } else {
    cout << "did not find " << Form("pd_run%d", run) << " in file " << f->GetName() << endl;
  }
  return key;
}


// ----------------------------------------------------------------------
void PdTrigger::setHLTKey(int run, TFile *f) {
  if (run == fRun) return;

  fHLTKey = getHLTKey(run, f); 
  fRun = run; 
  cout << "setHLTKey: " << fHLTKey << " for run " << fRun << endl;

  bool notFound(true);
  string pdfound(""); 
  if (0 == fPdTriggers.size()) {
    allPdTriggersFromFile(f);
    mkPdTriggersNoV();
    notFound = false; 
  } else {
    for (map<string, vector<string> >::iterator it = fPdTriggers.begin(); it != fPdTriggers.end(); ++it) {
      if (string::npos == it->first.find(fHLTKey)) {
	pdfound = it->first; 
      } else {
	notFound = false; 
	break;
      }
    }
  }

  if (notFound) {
    cout << "PdTrigger::setHLTKey> did not find " << fHLTKey << " for run " << run << "; inserting into map" << endl;
    allPdTriggersFromFile(f);
    mkPdTriggersNoV();
  } else {
    cout << "PdTrigger::setHLTKey> did find " << fHLTKey << " for run " << run << " in " << pdfound << endl;
  }
}

// ----------------------------------------------------------------------
bool PdTrigger::triggerInPd(string pdname, string triggername) {
  size_t pos = triggername.rfind("_v");
  map<string , vector<string> > *pMap(0);
  if (string::npos != pos) {
    pMap = &fPdTriggers;
    if (fVerbose > 1) cout << "using v triggername: " << triggername << endl;
  } else {
    pMap = &fPdTriggersNoV;
  }

  // -- add HLT key to pdname in case it is not contained yet
  string key(pdname); 
  if (string::npos == key.find(":")) {
    if (string::npos == fHLTKey.find("nada")) {
      key = fHLTKey + string(":") + pdname; 
    } else {
      cout << "error: HLTKey not set, and argument to triggerInPd does not contain HLT key; returning false" << endl;
      return false;
    }
  }
  if (0 == pMap->count(key)) {
    cout << "error: HLTKey:PD = " << key << " not found in map; returning false (but you should add this to PdTriggers!)" << endl;
    
    return false; 
  }

  vector<string> triggers = pMap->at(key);
  for (unsigned int i = 0; i < triggers.size(); ++i) {
    if (triggername == triggers[i]) return true;
  }
  return false;
}


// ----------------------------------------------------------------------
bool PdTrigger::triggerInPd(string hltkey, string pdname, string triggername) {
  size_t pos = triggername.rfind("_v");
  map<string , vector<string> > *pMap(0);
  if (string::npos != pos) {
    pMap = &fPdTriggers;
    if (fVerbose > 1) cout << "using v triggername: " << triggername << endl;
  } else {
    pMap = &fPdTriggersNoV;
  }

  string key = hltkey+string(":")+pdname;

  if (0 == pMap->count(key)) {
    cout << "error: HLTKey:PD = " << key << ", from " << hltkey << ":" << pdname 
	 << " not found in map; returning false (but you should add this to PdTriggers!)" 
	 << endl;
    return false; 
  }

  vector<string> triggers = pMap->at(key);
  for (unsigned int i = 0; i < triggers.size(); ++i) {
    if (triggername == triggers[i]) return true;
  }
  return false;
}


// ----------------------------------------------------------------------
void PdTrigger::addPdTriggersFromChain(string chain) {
  
  vector<string> vChain;
  cout << "reading all files form chain " << chain << endl;
  char pName[2000]; 
  int nentries; 
  int verbose(0); 
  char  buffer[1000];

  if ("all" == chain) {
    addPdTriggersFromChain("chains/v01/jobs/cbmm-v01-rereco-Run2012A__MuOnia__22Jan2013-v1");
    addPdTriggersFromChain("chains/v01/jobs/cbmm-v01-rereco-Run2012B__MuOnia__22Jan2013-v1");
    addPdTriggersFromChain("chains/v01/jobs/cbmm-v01-rereco-Run2012C__MuOnia__22Jan2013-v1");
    addPdTriggersFromChain("chains/v01/jobs/cbmm-v01-rereco-Run2012D__MuOnia__22Jan2013-v1");

    addPdTriggersFromChain("chains/v01/jobs/cbmm-v01-prompt-Run2015B__Charmonium__PromptReco-v1");
    addPdTriggersFromChain("chains/v01/jobs/cbmm-v01-prompt-Run2015D__Charmonium__PromptReco-v3");
    return;    
  }

  ifstream is(chain);
  while (is.getline(buffer, 1000, '\n')) {
    sscanf(buffer, "%s %d", pName, &nentries); 
    if (nentries > -1) {
      if (verbose) cout << pName << " -> " << nentries << " entries" << endl; 
      vChain.push_back(pName); 
    }
  }
  is.close();

  TFile *f(0); 
  for (unsigned int i = 0; i < vChain.size(); ++i) {
    //  for (unsigned int i = 0; i < 1; ++i) {
    cout << vChain[i] << endl;
    f = TFile::Open(vChain[i].c_str()); 

    allPdTriggersFromFile(f);
  }
  
  writePdTriggers("../common/pd/pdtriggers.txt");
}


// ----------------------------------------------------------------------
void PdTrigger::allPdTriggersFromFile(TFile *f) {
  TKey *key(0);   
  TH1D *h(0);
  string sname; 
  string::size_type m1, m2;
  string pd, lpd, hk; 

  TIter next(f->GetListOfKeys());
  vector<string> triggers;
  while ((key = (TKey*)next())) {
    sname = key->GetName();
    
    if (string::npos == sname.find("triggers_") || string::npos == sname.find("_run")) continue;
    m1 = sname.find("triggers_") + string("triggers_").size(); 
    m2 = sname.rfind("_run");       
    pd = sname.substr(m1, m2-m1); 
    
    lpd = pd;
    std::transform(lpd.begin(), lpd.end(), lpd.begin(), ::tolower);
    if (string::npos != lpd.find("alca")) {
      continue;
    }
    if (string::npos != lpd.find("test")) continue;
    if (string::npos != lpd.find("commissioning")) continue;
    if (string::npos != lpd.find("cosmics")) continue;
    if (string::npos != lpd.find("hcal")) continue;
    if (string::npos != lpd.find("monitor")) continue;
    if (string::npos != lpd.find("zerobias")) continue;
    if (string::npos != lpd.find("laser")) continue;
    if (string::npos != lpd.find("scouting")) continue;
    if (string::npos != lpd.find("l1accept")) continue;
    if (string::npos != lpd.find("express")) continue;
    if (string::npos != lpd.find("onlinehltresults")) continue;
    
    h = (TH1D*)f->Get(sname.c_str()); 
    sname = h->GetTitle(); 
    m1 = sname.find("("); 
    m2 = sname.find(")"); 
    hk = sname.substr(m1+1, m2-m1-1);
    
    for (int ibin = 1; ibin < h->GetNbinsX(); ++ibin) {
      triggers.push_back(h->GetXaxis()->GetBinLabel(ibin));
      //	cout << "   " << h->GetXaxis()->GetBinLabel(ibin) << endl;
    }
    
    if (0 == fPdTriggers.count(hk+string(":")+pd)) {
      cout << "inserting new hlt key: " << hk << " pd: " << pd << endl;
      fPdTriggers.insert(make_pair(hk+string(":")+pd, triggers));
    }
    triggers.clear();
    
    
  }
  f->Close();
}


// ----------------------------------------------------------------------
void PdTrigger::readPdTriggers(string file) {

  char  buffer[1000];
  ifstream is(file);
  string sbuffer, hk, pd;
  string::size_type m1, m2;
  vector<string> triggers;
  while (is.getline(buffer, 1000, '\n')) {
    sbuffer = buffer;
    if (string::npos != sbuffer.find(":")) {
      // -- if filled, put away old PD and its triggers
      if (triggers.size() > 0) {
	fPdTriggers.insert(make_pair(hk+string(":")+pd, triggers));
	triggers.clear();
      }
      m1 = sbuffer.find(":");
      hk = sbuffer.substr(0, m1);
      pd = sbuffer.substr(m1+1);
    } else{
      replaceAll(sbuffer, " ", "");
      triggers.push_back(sbuffer);
    }
  }
  is.close();
}




// ----------------------------------------------------------------------
void PdTrigger::writePdTriggers(string file) {

  ofstream os(file);
  for (map<string, vector<string> >::iterator it = fPdTriggers.begin(); it != fPdTriggers.end(); ++it) {
    os << (*it).first << endl;
    for (unsigned int i = 0; i < (*it).second.size(); ++i) {
      os << (*it).second[i] << endl;
    }
    
  }

  os.close();
}


// ----------------------------------------------------------------------
void PdTrigger::mkPdTriggersNoV() {

  fPdTriggersNoV.clear();
  
  string t; 
  for (map<string, vector<string> >::iterator it = fPdTriggers.begin(); it != fPdTriggers.end(); ++it) {
    vector<string> v; 
    for (unsigned int i = 0; i < (*it).second.size(); ++i) {
      t = (*it).second[i];
      size_t pos = t.rfind("_v");
      if (pos != string::npos) {
	t = t.substr(0, pos);
	v.push_back(t);
      }
    }
    fPdTriggersNoV.insert(make_pair((*it).first, v)); 
  }
  
}
