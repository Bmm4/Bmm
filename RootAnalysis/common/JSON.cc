#include "JSON.hh"

#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>

#include "util.hh"

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

  p1 = jsonFile[0].find("[[");
  if (string::npos == p1) {
    jsonFile = transform(jsonFile);
  }

  //cout<<" size "<<jsonFile.size()<<endl;

  int nrun(0), nls(0);
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
	nls += lumisection.second - lumisection.first + 1;
	//cout << "  " << sections[j] << " -> " << lumisection.first << " ... " << lumisection.second << endl;
      }
      fRunLsList.insert(make_pair(run, runLS));
      nrun++;

      p1 = jsonFile[i].find("\"", p3+1);
      p2 = jsonFile[i].find("\"", p1+1);
      p3 = jsonFile[i].find("]]", p2+1);
    }
  }

  fBegin = fRunLsList.begin();
  fEnd = fRunLsList.end();

  if (fVerbose > 0) print();
  cout << "     read " << nrun << " runs with " << nls << " lumisections" << endl;

}


// ----------------------------------------------------------------------
JSON::~JSON() {
}


// ----------------------------------------------------------------------
vector<string> JSON::transform(vector<string> dcs) {
  vector<string> runList, newList;
  string::size_type p1;
  string::size_type p2;
  string::size_type p3;
  string run("{");
  for (unsigned int i = 0; i < dcs.size(); ++i) {
    p1 = dcs[i].find("\"");
    p2 = dcs[i].find("\"", p1+1);
    if (string::npos != p1 && string::npos != p2) {
      for (unsigned ils = i; ils < dcs.size(); ++ils) {
	// -- break out if next run found
	if ((ils > i) && (string::npos != dcs[ils].find("\"")) && (string::npos != dcs[ils].find("\"", p1+1))) break;
	cleanupString(dcs[ils]);
	replaceAll(dcs[ils], ",", ", ");
	run += dcs[ils];
      }
      // cout <<  run << endl;
      newList.push_back(run);
      run = "";
    }
  }

  return newList;
}


// ----------------------------------------------------------------------
void JSON::print() {
  cout << "JSON::print() start" << endl;
  map<int, vector<pair<int, int> > >::iterator it = fRunLsList.begin();
  for (; it != fRunLsList.end(); ++it) {
    cout << it->first << " -> ";
    for (unsigned int i = 0; i < it->second.size(); ++i) {
      cout << it->second[i].first << ".." << it->second[i].second << "  ";
    }
    cout << endl;
  }
  cout << "JSON::print() end" << endl;
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
