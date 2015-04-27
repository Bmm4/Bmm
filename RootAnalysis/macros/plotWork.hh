#ifndef PLOTWORK_h
#define PLOTWORK_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotWork: public plotClass {

public :
  plotWork(std::string dir = "hpt0", std::string files = "plotWork.files", std::string setup = "default");
  virtual        ~plotWork();

  void   setCuts(std::string cuts);

  // -- Main analysis methods 
  void   makeAll(int bitmask = 0);
  void   bookHist(int mode);

  // -- trigger efficiency vs various variables
  void   triggerEff(std::string ds = "bssg", std::string dir = "candAnaMuMu"); 
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
  ClassDef(plotWork,1) 

};


#endif
