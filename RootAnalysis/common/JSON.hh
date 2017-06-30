#ifndef JSON_HH
#define JSON_HH

#include <utility>
#include <vector>
#include <string>
#include <map>

// ----------------------------------------------------------------------
// -- JSON
// -------
//
// class to parse a JSON file and provide an answer whether a run and lumisection are
// considered good (by being in that JSON file)
//
// Usage:
//          JSON a("/shome/ursl/json/1.json");
//          if (a.good(fRun, fLS)) cout << "This is a good run and lumisection" << endl;
//
// ----------------------------------------------------------------------

class JSON {
public:
  JSON(const char *fname, int verbose = 1);
  ~JSON();

  void print();
  bool good(int run, int lumisection);
  bool goodRun(int run);

private:

  std::vector<std::string> transform(std::vector<std::string> dcsList);
  std::vector<std::string> parseLS(std::string &sLs);
  std::pair<int, int> ls(std::string ls);

  std::map<int, std::vector<std::pair<int, int> > > fRunLsList;
  std::map<int, std::vector<std::pair<int, int> > >::iterator fBegin;
  std::map<int, std::vector<std::pair<int, int> > >::iterator fEnd;

  int fVerbose;
  int fCountGood, fCountBad;
};


#endif
