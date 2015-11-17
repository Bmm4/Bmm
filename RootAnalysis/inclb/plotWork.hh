#ifndef PLOTWORK_h
#define PLOTWORK_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "RedTreeData.hh"

// ----------------------------------------------------------------------
class plotWork: public plotClass {

public :
  plotWork(std::string dir = "results", std::string files = "plotWork.files", std::string setup = "default");
  virtual        ~plotWork();

  void   setCuts(std::string cuts);

  // -- Main analysis methods 
  void   loadFiles(std::string afiles);
  void   makeAll(int bitmask = 0, std::string what = "run251721");
  void   bookHist(int mode);

  // -- validate production
  void   validation(std::string hist, std::string dir, std::string dname, std::string bname, std::string cname); 
  void   validation(std::string hist1, std::string hist2, double xmin, double xmax,
		    std::string dname, std::string bname, std::string cname);
  TH1D*  getPtRel(std::string histname, std::string dir, std::string ds, double xmin, double xmax);
  void   loopFunction1(); 


  void   setupTree(TTree *t); 

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  void   candAnalysis(); 

private: 

  TTree* fTree;

  struct RedTreeData fb; 


  bool fGoodCand;
  double PTLO;

  
  // ----------------------------------------------------------------------
  ClassDef(plotWork,1) 

};


#endif
