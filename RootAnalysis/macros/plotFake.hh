#ifndef PLOTFAKE_h
#define PLOTFAKE_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotFake: public plotClass {

public :
  plotFake(std::string dir = "results", std::string files = "plotFake.files", std::string setup = "default");
  virtual        ~plotFake();

  void   setCuts(std::string cuts);

  // -- Main analysis methods 
  void   loadFiles(std::string afiles);
  void   makeAll(int bitmask = 0);
  void   bookHist(int mode);

  // -- Fake rate from light hadrons
  void   fake1(int id = 211); 
  void   loopFunction1(); 
 

  void   setupTree(TTree *t); 

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  void   candAnalysis(); 

private: 

  TTree* fTree;

  struct redTreeData fb; 


  bool fGoodCand;
  double PTLO;

  
  // ----------------------------------------------------------------------
  ClassDef(plotFake,1) 

};


#endif
