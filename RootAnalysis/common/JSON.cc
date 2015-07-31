#include "JSON.hh"

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

// ----------------------------------------------------------------------
JSON::JSON(const char *fname, int verbose) {
  fVerbose = verbose; 
  fCountGood = fCountBad = 0; 

  cout << "JSON initializing from " << fname << endl;

  vector<string> jsonFile; 
  string line; 
  ifstream file(fname);  
  if (!file.is_open()) {
    cout << "ERROR: could not open " << fname << endl;
    return;
  }
  while (getline(file, line)) { 
    jsonFile.push_back(line); 
    //    cout<<line<<endl;
  }

  string::size_type p1;
  string::size_type p2;
  string::size_type p3;

  //cout<<" size "<<jsonFile.size()<<endl;

  for (unsigned int i = 0; i < jsonFile.size(); ++i) {
    p1 = jsonFile[i].find("\""); 
    p2 = jsonFile[i].find("\"", p1+1); 
    p3 = jsonFile[i].find("]],", p2+1); 
    if(p3 == string::npos ) p3 = jsonFile[i].find("]]}", p2+1); // take care of the last line  
    while (p1 != string::npos && p2 != string::npos && p3 != string::npos ) {
      string sRun = jsonFile[i].substr(p1+1, p2-p1-1); 
      string sLs  = jsonFile[i].substr(p2+3, p3-p2-1); 
      //cout << "run: " << sRun << " -> LS: " << sLs << " "<<i<<endl;
      int run = atoi(sRun.c_str());
      vector<string> sections = parseLS(sLs); 

      pair<int, int> lumisection;
      vector<pair<int, int> > runLS; 
      for (unsigned int j = 0; j < sections.size(); ++j) {
	lumisection = ls(sections[j]); 
	runLS.push_back(lumisection); 
	//cout << "  " << sections[j] << " -> " << lumisection.first << " ... " << lumisection.second << endl;
      }
      fRunLsList.insert(make_pair(run, runLS)); 
     
      p1 = jsonFile[i].find("\"", p3+1); 
      p2 = jsonFile[i].find("\"", p1+1); 
      p3 = jsonFile[i].find("]]", p2+1); 
    }
  }

  fBegin = fRunLsList.begin(); 
  fEnd = fRunLsList.end(); 

  if (fVerbose > 0) print();

}


// ----------------------------------------------------------------------
JSON::~JSON() {
}


// ----------------------------------------------------------------------
void JSON::print() {
  map<int, vector<pair<int, int> > >::iterator it = fRunLsList.begin(); 
  for (; it != fRunLsList.end(); ++it) {
    cout << it->first << " -> ";
    for (unsigned int i = 0; i < it->second.size(); ++i) {
      cout << it->second[i].first << ".." << it->second[i].second << "  ";
    }
    cout << endl;
  }

}


// ----------------------------------------------------------------------
bool JSON::good(int run, int lumisection) {
  map<int, vector<pair<int, int> > >::iterator it;
  // FIXME: Is it faster to directly access the pair with the key in the map?
  for (it = fBegin; it != fEnd; ++it) {
    if (run == it->first) {
      for (unsigned int i = 0; i < it->second.size(); ++i) {
	if (it->second[i].first <= lumisection && lumisection <= it->second[i].second) {
	  return true;
	}
      }
      return false;
    } 
  }

  return false;
}


// ----------------------------------------------------------------------
bool JSON::goodRun(int run) {
  map<int, vector<pair<int, int> > >::iterator it;
  // FIXME: Is it faster to directly access the pair with the key in the map?
  for (it = fBegin; it != fEnd; ++it) {
    if (run == it->first) {
      return true;
    }
  }

  return false;
}


// ----------------------------------------------------------------------
pair<int, int> JSON::ls(string ls) {
  int ls1(-1), ls2(-1); 
  sscanf(ls.c_str(), "[%d, %d]", &ls1, &ls2);
  return make_pair(ls1, ls2); 
}


// ----------------------------------------------------------------------
vector<string> JSON::parseLS(string &sLs) {
  vector<string> a;
  string sls = sLs;
  string::size_type s1 = sls.find("]");
  string::size_type s2 = sls.find("]", s1+1);
  if (s2 == s1+1) {
    sls = sls.replace(0, 1, ""); 
    sls = sls.replace(sls.size()-1, sls.size(), ""); 
    a.push_back(sls); 
  } else {
    while (s1 < sls.size()) {
      string ls = sls.substr(1, s1); 
      sls = sls.replace(0, s1+2, ""); 
      s1 = sls.find("]"); 
      a.push_back(ls); 
    }
  }
  return a;
}
