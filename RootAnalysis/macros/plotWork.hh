#ifndef PLOTWORK_h
#define PLOTWORK_h

#include "plotClass.hh"
#include "common/dataset.hh"
#include "common/AnalysisDistribution.hh"
#include "redTreeData.hh"

// ----------------------------------------------------------------------
class plotWork: public plotClass {

public :
  plotWork(std::string dir = "results",
	   std::string files = "plotWork.files",
	   std::string cuts = "plotClass.cuts",
	   std::string setup = "default");
  virtual        ~plotWork();

  void   setCuts(std::string cuts);

  // -- Main analysis methods
  void   loadFiles(std::string afiles);
  void   makeAll(int bitmask = 0);

  // -- validate gen/mc production
  void   prodSummary(string ds1, int year = 2014);
  void   privateVsOfficial(string mode = "bdmm");
  void   bookDistributions(std::string sample);
  AnalysisDistribution* bookDistribution(std::string hn, std::string ht, std::string hc, int nbins, double lo, double hi);

  // -- code for loops
  void   loopFunction1();

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);
  //  void   candAnalysis();

private:

  TTree* fTree;

#define NAD 2
  AnalysisDistribution
  *fpMuonsEta[NAD], *fpMuon1Pt[NAD], *fpMuon2Pt[NAD]
    , *fpPt[NAD], *fpP[NAD], *fpPz[NAD], *fpEta[NAD]
    , *fpAlpha[NAD]
    , *fpIso[NAD], *fpCloseTrk[NAD], *fpDocaTrk[NAD]
    , *fpChi2Dof[NAD], *fpPChi2Dof[NAD]
    , *fpFLS3d[NAD], *fpFL3d[NAD], *fpFL3dE[NAD]
    , *fpMaxDoca[NAD], *fpIp[NAD], *fpIpS[NAD]
    , *fpPvZ[NAD], *fpPvN[NAD], *fpPvAveW8[NAD]
    ;

  int fOffset;
  int fChanMode;



  bool fGoodCand;
  double PTLO;


  // ----------------------------------------------------------------------
  ClassDef(plotWork,1)

};


#endif
