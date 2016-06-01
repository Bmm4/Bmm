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


struct adset {
  AnalysisDistribution
  *fpMuonsEta, *fpMuon1Pt, *fpMuon2Pt
    , *fpPt, *fpP, *fpPz, *fpEta
    , *fpAlpha
    , *fpIso, *fpCloseTrk, *fpDocaTrk
    , *fpChi2Dof, *fpPChi2Dof
    , *fpFLS3d, *fpFL3d, *fpFL3dE
    , *fpMaxDoca, *fpIp, *fpIpS
    , *fpBDT
    , *fpBDTSel0, *fpBDTSel1, *fpBDTSel2
    , *fpPvZ, *fpPvN, *fpPvAveW8
    , *fpLip, *fpLipS, *fpLip2, *fpLipS2
    , *fpCloseTrkS1, *fpCloseTrkS2, *fpCloseTrkS3
    , *fpM1Iso, *fpM2Iso
    , *fpPvDchi2, *fpOtherVtx
    , *fpLastCut
    ;

};



class plotReducedOverlays: public plotClass {
public:
  plotReducedOverlays(std::string dir = "results",
		      std::string files = "plotReducedOverlays.files",
		      std::string cuts = "plotClass.2016.cuts",
		      std::string setup = "");
  ~plotReducedOverlays();

  virtual void   init();
  virtual void   loadFiles(std::string afiles);
  virtual void   loopFunction1();
  virtual void   loopFunction2();
  virtual void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);


  void makeAll(std::string selection = "Presel");

  /// run over sample with specified channel and selection and dump AD into rootfile
  void makeSample(std::string sample1,
		  std::string selection = "Presel",
		  int nevents = -1, int nstart = 0);

  /// overlay two sample previously makeSample'd (read AD in from rootfile)
  void makeOverlay(std::string sample1,
		   std::string sample2,
		   std::string selection);

  /// combine the previous two methods
  void makeSampleOverlay(std::string sample1,
			 std::string sample2,
			 std::string selection = "Presel");

  /// generate mass plots for illustration of the sbs-subtracted plots
  void plotMass(std::string sample = "bupsikData", std::string selection = "Presel");


  /// I forgot what this is for
  void makeOverlay2Channels(std::string sample,
			    std::string channel1,
			    std::string channel2,
			    std::string selection);

  void allSystematics();
  void systematics(std::string sample1,
		   std::string selection,
		   int chan = 0);
  void compareBsAndBp(std::string file = "2012/anaBmm.plotReducedOverlaysSystematics.2012.root");

  // --
  void bookDistributions();  void bookDistributions0();
  void fillDistributions();  void fillDistributions0();
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

  std::string fChannel;

  std::vector<std::string> fDoList, fChannelList;

  std::map<string, adset*> fAdMap;

  bool fSel0, fSel1, fSel2;

  ClassDef(plotReducedOverlays,1) //Testing plotReducedOverlays

};



#endif
