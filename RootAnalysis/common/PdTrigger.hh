#ifndef PDTRIGGER_HH
#define PDTRIGGER_HH

#include <utility>
#include <vector>
#include <string>
#include <map>

// ----------------------------------------------------------------------
// -- PdTrigger
// ------------
//
// class to provide information about trigger/primary datasets
//
// ----------------------------------------------------------------------

class PdTrigger {
public:
  PdTrigger(std::string basedir, int verbose = 1);
  ~PdTrigger();

  void print(); 
  bool triggerInPd(std::string pdname, std::string triggername); 

  std::vector<std::string> glob(std::string basedir, std::string basename);
  void readFile(std::string filename, std::vector<std::string> &lines);

private: 

  std::map<std::string /*pd*/, std::vector<std::string> /*triggers*/> fPdTrigger;

  int fVerbose; 
};


#endif
