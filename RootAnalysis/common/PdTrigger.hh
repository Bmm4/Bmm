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

class TFile; 

class PdTrigger {
public:
  PdTrigger(std::string dir, int verbose = 1);
  ~PdTrigger();

  void        print(); 

  bool        triggerInPd(std::string pdname, std::string triggername); 
  bool        triggerInPd(std::string hltkey, std::string pdname, std::string triggername); 
  std::string getHLTKey(int run, TFile *f);
  void        setHLTKey(int run, TFile *f);
  void        setHLTKey(std::string hltkey) {fHLTKey = hltkey;}

  // -- create information for PdTriggers
  void        addPdTriggers(std::string chain); 
  void        readPdTriggers(std::string file); 
  void        writePdTriggers(std::string file); 
  void        mkPdTriggersNoV();
  
private: 

  std::string fHLTKey;

  std::map<std::string /*hltKey:PD*/, std::vector<std::string> /*triggers*/> fPdTriggers;
  std::map<std::string /*hltKey:PD*/, std::vector<std::string> /*triggers*/> fPdTriggersNoV;

  int fVerbose; 
};


#endif
