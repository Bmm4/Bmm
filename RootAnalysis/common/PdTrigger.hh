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
  // -- query map whether a triggername is an ingredient if a primary dataset
  bool        triggerInPd(std::string hltkey, std::string pdname, std::string triggername); 
  // -- query map whether a triggername is an ingredient if a primary dataset
  //    This version assume that fHLTKey has been set beforehand
  bool        triggerInPd(std::string pdname, std::string triggername); 
  // -- set the HLT key 
  void        setHLTKey(int run, TFile *f);
  std::string getHLTKey(int run, TFile *f);

  void        mkPdTriggersNoV();
  // -- read in all pd/trigger information from filename and fill into fPdTriggers/fPdTriggersNoV
  void        allPdTriggersFromFile(TFile *f);
  // -- create information for PdTriggers from chains, stored in summary file (deprecated)
  void        addPdTriggersFromChain(std::string chain); 
  void        writePdTriggers(std::string file); 
  void        readPdTriggers(std::string file); 
  
private: 
  int         fRun; 
  std::string fHLTKey;

  std::map<std::string /*hltKey:PD*/, std::vector<std::string> /*triggers*/> fPdTriggers;
  std::map<std::string /*hltKey:PD*/, std::vector<std::string> /*triggers*/> fPdTriggersNoV;

  int fVerbose; 
};


#endif
