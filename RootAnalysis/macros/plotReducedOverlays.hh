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
  *fpMuonsEta, *fpMuonsPhi, *fpMuon1Pt, *fpMuon2Pt
    , *fpKaonsEta, *fpKaonsPhi, *fpKaonsPt
    , *fpPt, *fpP, *fpPz, *fpEta, *fpPhi
    , *fpAlpha
    , *fpIso, *fpCloseTrk, *fpDocaTrk
    , *fpChi2Dof, *fpPChi2Dof
    , *fpFLS3d, *fpFL3d, *fpFL3dE
    , *fpFLSxy, *fpFLxy, *fpFLxyE
    , *fpMaxDoca, *fpIp, *fpIpS
    , *fpBDT
    , *fpBDTSel0, *fpBDTSel1, *fpBDTSel2
    , *fpPvZ, *fpPvN, *fpPvNtrk, *fpPv2Ntrk, *fpPvAveW8
    , *fpLip, *fpLipS, *fpLip2, *fpLipS2, *fpOsmDr, *fpDzmin,  *fpDz12
    , *fpCloseTrkS1, *fpCloseTrkS2, *fpCloseTrkS3
    , *fpM1Iso, *fpM2Iso
    , *fpPvDchi2, *fpOtherVtx
    , *fpTau
    , *fpLastCut
    // -- BEGIN for tracking studies
    , *fpMupixhits, *fpMutrkhits, *fpMutkqual, *fpMualg, *fpMuvalhits, *fpMuvalhitfraction, *fpMulayerswithhits
    , *fpMuchi2, *fpMudz, *fpMudzE, *fpMud0, *fpMud0E, *fpMudsz, *fpMudszE, *fpMudxy, *fpMudxyE, *fpMuptE, *fpMuptEpt, *fpMuetaE, *fpMuphiE
    , *fpKapixhits, *fpKatrkhits, *fpKatkqual, *fpKaalg, *fpKavalhits, *fpKavalhitfraction, *fpKalayerswithhits
    , *fpKachi2, *fpKadz, *fpKadzE, *fpKad0, *fpKad0E, *fpKadsz, *fpKadszE, *fpKadxy, *fpKadxyE, *fpKaptE, *fpKaptEpt, *fpKaetaE, *fpKaphiE
    // -- END for tracking studies
    ;
};



class plotReducedOverlays: public plotClass {
public:
  plotReducedOverlays(std::string dir = "results",
		      std::string files = "plotResults.2016.files",
		      std::string cuts = "baseCuts.2016.cuts",
		      std::string setup = "",
		      int year = 0);
  ~plotReducedOverlays();

  virtual void   init();
  virtual void   loadFiles(std::string afiles);
  virtual void   loopFunction1();
  virtual void   loopFunction2();
  virtual void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);


  void makeAll(std::string what = "all");

  // -- run over sample with specified channel and selection and dump AD into rootfile
  void makeSample(std::string sample1,
		  int nevents = -1, int nstart = 0);

  // -- overlay two sample previously makeSample'd (read AD in from rootfile)
  void makeOverlay(std::string sample1, std::string sample2, std::string selection = "");

  // -- combine the previous two methods
  void makeSampleOverlay(std::string sample1, std::string sample2, int nevents = -1);

  // -- generate mass plots for illustration of the sbs-subtracted plots
  void plotMass(std::string sample = "bupsikData", std::string selection = "Presel");


  // -- Validation of MC samples
  void comparePrivateAndOfficial();

  // -- Derivation of scale factors for 2016 MC
  void scaleFactors2016(std::string filename);

  // -- I forgot what this is for
  void makeOverlay2Channels(std::string sample,
			    std::string channel1,
			    std::string channel2,
			    std::string selection);

  void allSystematics();
  void sysBdtCut(std::string sample1, std::string sample2, std::string selection, std::string file2 = "nada", bool bdtscan = false);
  void sysComparison(std::string sample1, std::string sample2, std::string selection, std::string file2 = "nada");
  void compareBsAndBp(std::string file = "2012/anaBmm.plotReducedOverlaysSystematics.2012.root");

  // --
  void bookDistributions(std::string selmode = "");
  void fillDistributions(std::string selmode = "");
  void sbsDistributions(std::string sample, std::string selection = "Presel", std::string what = "");
  void overlay(std::string sample1, std::string sample2, std::string selection = "Presel", std::string what = "");
  void overlayAndRatio(TCanvas *c, TH1D *h1, TH1D *h2);
  void overlay2Files(std::string file1, std::string sample1, std::string chan1, std::string selection1,
		     std::string file2, std::string sample2, std::string chan2, std::string selection2,
		     std::string what = "");
  void overlay3Samples(std::string sample1, std::string file1,
		       std::string sample2, std::string file2,
		       std::string sample3, std::string file3,
		       std::string selection);

  AnalysisDistribution* bookDistribution(string hn, string ht, string hc, AnalysisCuts *pCuts,
					 int nbins, double lo, double hi, bool *presel = 0);


  TDirectory *fHistDir;

  std::string fChannel;

  std::vector<std::string> fDoList, fChannelList;

  std::map<string, adset*> fAdMap;

  bool fIncludeOverflowInLastBin, fDoCNC;
  bool fSel0, fSel1, fSel2;

  ClassDef(plotReducedOverlays,1) //Testing plotReducedOverlays

};



#endif
