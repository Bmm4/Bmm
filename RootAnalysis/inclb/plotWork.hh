#ifndef PLOTWORK_h
#define PLOTWORK_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "RedTreeData.hh"

struct numbers {
  void init() {hD = hSum = hB = hC = hL = 0; status = -1;}
  int status; 
  double minPt, maxPt, meanPt, minEta, maxEta, meanEta;
  double nData2, nB2, nC2, nL2;     // Integral()
  double nData, nB, nC, nL;         // GetSumOfWeights()
  double nDataE0, nBE0, nCE0, nLE0; // central value
  double fracB, fracC, fracL;       // central value
  double fracBE0, fracCE0, fracLE0; // statistical error
  double effB; 
  std::string sname;         // pd name: dataMu24, dataMuHIL2Mu3, ...
  std::string dname;         // directory name: candAnaMu8, candAnaMuHIL2Mu3, ...
  TH1D *hD, *hSum, *hB, *hC, *hL, *hPlot;   // histograms with scale factors from fit applied
};


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
  void   validation(std::string hist1, std::string hist2, double xmin, double xmax,
		    std::string dname, std::string bname, std::string cname);
  void   loopFunction1(); 

  // -- data analysis
  void   hinValidation(); 
  void   dSigmadPt(); 
  
  TH1D*  getPtRel(std::string histname, std::string dir, std::string ds, double xmin, double xmax);
  void   fitPtRel(numbers *, TH1D* hd, TH1D* hb, TH1D* hc, TH1D* hl = 0); 
  void   efficiency(numbers *, std::string bString); 

  void   setupTree(TTree *t); 

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0); 
  void   candAnalysis(); 

private: 

  TTree* fTree;

  struct RedTreeData fb; 


  bool fGoodCand;
  double PTLO;

  std::vector<numbers*> fVectorResult; 
  // ----------------------------------------------------------------------
  ClassDef(plotWork,1) 

};


#endif
