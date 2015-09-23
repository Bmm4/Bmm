#ifndef PLOTREDUCEDOVERLAYS
#define PLOTREDUCEDOVERLAYS

#include "plotClass.hh"

#include <TDirectory.h>
#include <TH1.h>
#include <TTree.h>

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "common/AnalysisDistribution.hh"


class plotReducedOverlays: public plotClass {

public:

  plotReducedOverlays(std::string dir = "results", std::string files = "plotReducedOverlays.files", std::string setup = "2015");
  ~plotReducedOverlays();

  virtual void   loadFiles(std::string afiles);
  virtual void   loopFunction1(); 
  virtual void   loopFunction2(); 


  void makeAll(std::string selection = "Presel"); 
  void makeSampleOverlay(std::string sample1, std::string sample2, std::string channel, std::string selection);
  void makeSample(std::string sample1, std::string selection, std::string = "B", int nevents = -1, int nstart = 0);

  void makeOverlay(std::string sample1, std::string sample2, std::string channel, std::string selection);
  void makeOverlay2Channels(string sample, string channel1, string channel2, string selection);
  void allSystematics();
  void systematics(std::string sample1, std::string selection, int chan = 0);
  void compareTheYears(std::string sample, std::string channel, std::string file1, std::string file2);   
  void compareBsAndBp(std::string file = "2012/anaBmm.plotReducedOverlaysSystematics.2012.root"); 

  // -- 
  void bookDistributions(std::string sample);
  void fillDistributions();
  void sbsDistributions(std::string sample, std::string selection = "Presel", std::string what = "");
  void overlay(std::string sample1, std::string sample2, std::string selection = "Presel", std::string what = ""); 
  void sbsSingleFile(string file1, string sample1 = "NoData", string channel = "B", string selection = "Presel");
  void overlay2Files(std::string file1, std::string sample1, 
		     std::string file2, std::string sample2, 
		     std::string chan1 = "BLoPU", std::string chan2 = "BHiPU", 
		     std::string selection = "Presel", std::string what = ""); 
  
  AnalysisDistribution* bookDistribution(std::string hn, std::string ht, std::string hc, int nbins, double lo, double hi); 
  AnalysisDistribution* bookSpecialDistribution(string hn, string ht, string hc, int nbins, double lo, double hi, bool *presel);


  TDirectory *fHistDir; 

#define NAD 2
  AnalysisDistribution   
  *fpMuonsEta[NAD], *fpMuon1Pt[NAD], *fpMuon2Pt[NAD]
    , *fpPt[NAD], *fpP[NAD], *fpPz[NAD], *fpEta[NAD] 
    , *fpAlpha[NAD]
    , *fpIso[NAD], *fpCloseTrk[NAD], *fpDocaTrk[NAD]   
    , *fpChi2Dof[NAD], *fpPChi2Dof[NAD] 
    , *fpFLS3d[NAD], *fpFL3d[NAD], *fpFL3dE[NAD] 
    , *fpMaxDoca[NAD], *fpIp[NAD], *fpIpS[NAD]
    , *fpBDT[NAD]   
    , *fpBDTSel0[NAD], *fpBDTSel1[NAD], *fpBDTSel2[NAD]
    , *fpPvZ[NAD], *fpPvN[NAD], *fpPvAveW8[NAD]
    , *fpLip[NAD], *fpLipS[NAD], *fpLip2[NAD], *fpLipS2[NAD]
    , *fpCloseTrkS1[NAD], *fpCloseTrkS2[NAD], *fpCloseTrkS3[NAD]
    , *fpM1Iso[NAD], *fpM2Iso[NAD]
    , *fpPvDchi2[NAD], *fpOtherVtx[NAD]
    ;
  int fOffset;
  int fChanMode;

  vector<std::string> fDoList; 

  bool fSel0, fSel1, fSel2; 

  ClassDef(plotReducedOverlays,1) //Testing plotReducedOverlays

};



#endif

