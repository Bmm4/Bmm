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
    , *fpTau
    , *fpLastCut
    ;

};



class plotReducedOverlays: public plotClass {
public:
  plotReducedOverlays(std::string dir = "results",
		      std::string files = "plotResults.2016.files",
		      std::string cuts = "baseCuts.cuts",
		      std::string setup = "");
  ~plotReducedOverlays();

  virtual void   init();
  virtual void   loadFiles(std::string afiles);
  virtual void   loopFunction1();
  virtual void   loopFunction2();
  virtual void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);


  void makeAll(std::string what = "all");

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


  /// Validation of MC samples
  void comparePrivateAndOfficial();

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
  void bookDistributions();
  void fillDistributions();
  void sbsDistributions(std::string sample, std::string selection = "Presel", std::string what = "");
  void overlay(std::string sample1, std::string sample2, std::string selection = "Presel", std::string what = "");
  void overlayAndRatio(TCanvas *c, TH1D *h1, TH1D *h2);
  void overlayOld(std::string sample1, std::string sample2, std::string selection = "Presel", std::string what = "");
  void sbsSingleFile(string file1, string sample1 = "NoData", string channel = "B", string selection = "Presel");
  void overlay2Files(std::string file1, std::string sample1,
		     std::string file2, std::string sample2,
		     std::string chan1 = "BLoPU", std::string chan2 = "BHiPU",
		     std::string selection = "Presel", std::string what = "");

  AnalysisDistribution* bookDistribution(std::string hn, std::string ht, std::string hc, int nbins, double lo, double hi);
  AnalysisDistribution* bookSpecialDistribution(string hn, string ht, string hc, int nbins, double lo, double hi, bool *presel);


  TDirectory *fHistDir;

  std::string fChannel;

  std::vector<std::string> fDoList, fChannelList;

  std::map<string, adset*> fAdMap;

  bool fIncludeOverflowInLastBin;
  bool fSel0, fSel1, fSel2;

  ClassDef(plotReducedOverlays,1) //Testing plotReducedOverlays

};



#endif
