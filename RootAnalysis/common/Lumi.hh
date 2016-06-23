#ifndef LUMI_HH
#define LUMI_HH

#include <utility>
#include <vector>
#include <string>
#include <map>

// ----------------------------------------------------------------------
// -- Lumi
// -------
//
// class to parse the output of bril_calc and provide a map run:lumi
//
// Usage:
//          Lumi a("/shome/ursl/json/1.lumi");
//          double lumi = a.lumi(run);
//
// ----------------------------------------------------------------------

class Lumi {
public:
  Lumi(std::string fname, int verbose = 1);
  ~Lumi();

  void parse(std::string fname);

  void print();
  double lumi(int run);
  bool contains(int run);

private:

  std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems);
  std::vector<std::string> split(const std::string &s, char delim);


  std::map<int, double> fMapRunLumi;
  std::string           fNormtag;
  int fVerbose;
};


#endif
