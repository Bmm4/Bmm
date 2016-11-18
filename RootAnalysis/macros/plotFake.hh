#ifndef PLOTFAKE_h
#define PLOTFAKE_h

#include "plotClass.hh"
#include "common/AnalysisDistribution.hh"
#include "common/dataset.hh"
#include "redTreeData.hh"

struct adsetFake {
  AnalysisDistribution
  *fpAllEta,
    *fpAllPt,

    *fpFakeEta,
    *fpFakePt,
    *fpFakeBdt,
    *fpFakeTip,
    *fpFakeLip,
    *fpFakeInnerChi2,
    *fpFakeOuterChi2,
    *fpFakeChi2LocalPosition,
    *fpFakeChi2LocalMomentum,
    *fpFakeStaTrkMult,
    *fpFakeTmTrkMult,
    *fpFakeDeltaR,
    *fpFakeItrkValidFraction,
    *fpFakeSegmentComp,
    *fpFakeGtrkNormChi2,
    *fpFakeDzRef,
    *fpFakeDxyRef,
    *fpFakeGtrkProb,
    *fpFakeMuonChi2,
    *fpFakeGlbKinkFinder,
    *fpFakeStaRelChi2,
    *fpFakeTrkRelChi2,
    *fpFakeGlbDeltaEtaPhi,
    *fpFakeTimeInOut,
    *fpFakeTimeInOutE,
    *fpFakeTimeInOutS,

    *fpFakeNvalidMuonHits,
    *fpFakeNmatchedStations,
    *fpFakeLayersWithHits,
    *fpFakeNumberOfValidTrkHits,
    *fpFakeNumberOfLostTrkHits,
    *fpFakeNumberOfValidPixHits,
    *fpFakeRPChits1,
    *fpFakeRPChits2,
    *fpFakeRPChits3,
    *fpFakeRPChits4,
    *fpFakeCombHits,

  // require only TIS:
    *fpFakeTisAllEta,
    *fpFakeTisAllPt,
    *fpFakeTisFakeEta,
    *fpFakeTisFakePt,

  // require TIS && separation to trigger
    *fpFakeTisDtFakePt,
    *fpFakeTisDtFakeEta,
    *fpFakeTisDtAllPt,
    *fpFakeTisDtAllEta,

  // require TIS && separation to trigger && distance to muon
    *fpFakeTisDtDmFakePt,
    *fpFakeTisDtDmFakeEta,
    *fpFakeTisDtDmAllPt,
    *fpFakeTisDtDmAllEta;

};


// ----------------------------------------------------------------------
class plotFake: public plotClass {

public :
  plotFake(std::string dir = "results",
	   std::string files = "plotResults.2016.files",
	   std::string cuts = "baseCuts.cuts",
	   std::string setup = "");
  virtual        ~plotFake();

  // -- Main analysis methods
  virtual void bookHist(int mode);
  virtual void init();
  virtual void loadFiles(std::string afiles);
  virtual void makeAll(std::string what = "all");

  // -- validation of GM distributions for hadrons and muons from control samples
  void makeSample(std::string dataset, std::string cand, int nevents = -1, int nstart = 0);
  void makeOverlay(std::string what1, std::string what2, std::string selection = "Cu");
  void fillDistributions();
  void bookDistributions();
  AnalysisDistribution* bookDistribution(std::string hn, std::string ht, std::string hc, int nbins, double lo, double hi);
  void sbsDistributions(std::string sample, std::string selection = "Cu", std::string what = "");
  void overlay(std::string sample1, std::string sample2, std::string what = "");
  void fakeRate(std::string dataset1 = "fakeData_ks", std::string dataset2 = "fakeMc_ks",
		std::string var = "FakeTisFakePt", std::string varA = "FakeTisAllPt",
		double ymax = 0.02);

  // -- utilities for control sample optimization
  void   playKs(std::string cuts = "flxy<4&&flsxy>15&&pvips<3&&prob>0.01&&abs(m1ips)>5&&abs(m2ips)>5&&m1q*m2q<0", std::string name = "");
  void   playPhi(std::string cuts = "m1q*m2q<0&&prob>0.01&&maxdoca<0.002&&m2pt>4&&m1pv==m2pv&&m1bpixl1&&m2bpixl1", std::string name = "");
  void   playLambda(std::string cuts = "flxy<4&&flsxy>15&&fls3d>15&&chi2<3&&pvips<3", std::string name = "");
  void   fitKs(TH1D*);
  void   fitPhi(TH1D*);
  void   fitLambda(TH1D*);


  void   loopFunction1();

  void   setupTree(TTree *t);

  void   loopOverTree(TTree *t, int ifunc, int nevts = -1, int nstart = 0);
  void   analysis();

private:

  TTree* fTree;

  struct redTreeData fb;

  std::map<string, adsetFake*> fAdMap;

  std::string fChannel;
  std::vector<std::string> fDoList, fChannelList, fChannelSample;

  double PTLO;

  static const int NTRKMAX = 10;
  double fCandM, fCandPvIp, fCandPvIpS, fCandChi2Dof, fCandFLS3d, fCandFL3d, fCandFLxy, fCandFL3dE, fCandFLSxy, fCandDoca;
  bool   fTIS, fCowboy;
  int    fFakeNtrk, fFakeId[NTRKMAX], fFakeQ[NTRKMAX], fFakeGm[NTRKMAX];
  float  fFakePt[NTRKMAX], fFakeEta[NTRKMAX], fFakePhi[NTRKMAX], fFakeBdt[NTRKMAX], fFakeDtrig[NTRKMAX], fFakeDmuon[NTRKMAX];
  float  fFakeTip[NTRKMAX], fFakeLip[NTRKMAX], fFakeInnerChi2[NTRKMAX]
    , fFakeOuterChi2[NTRKMAX]
    , fFakeChi2LocalPosition[NTRKMAX]
    , fFakeChi2LocalMomentum[NTRKMAX]
    , fFakeStaTrkMult[NTRKMAX]
    , fFakeTmTrkMult[NTRKMAX]
    , fFakeDeltaR[NTRKMAX]
    , fFakeItrkValidFraction[NTRKMAX]
    , fFakeSegmentComp[NTRKMAX]
    , fFakeGtrkNormChi2[NTRKMAX]
    , fFakeDzRef[NTRKMAX]
    , fFakeDxyRef[NTRKMAX]
    , fFakeGtrkProb[NTRKMAX]
    , fFakeMuonChi2[NTRKMAX]
    , fFakeGlbKinkFinder[NTRKMAX]
    , fFakeStaRelChi2[NTRKMAX]
    , fFakeTrkRelChi2[NTRKMAX]
    , fFakeGlbDeltaEtaPhi[NTRKMAX]
    , fFakeTimeInOut[NTRKMAX]
    , fFakeTimeInOutE[NTRKMAX]
    ;

  int    fFakeNvalidMuonHits[NTRKMAX]
    , fFakeNmatchedStations[NTRKMAX]
    , fFakeLayersWithHits[NTRKMAX]
    , fFakeNumberOfValidTrkHits[NTRKMAX]
    , fFakeNumberOfLostTrkHits[NTRKMAX]
    , fFakeNumberOfValidPixHits[NTRKMAX]
    , fFakeRPChits1[NTRKMAX]
    , fFakeRPChits2[NTRKMAX]
    , fFakeRPChits3[NTRKMAX]
    , fFakeRPChits4[NTRKMAX]
    , fFakeCombHits[NTRKMAX]
    ;


  bool fGoodCand, fGlobalMuon, fGoodPt, fGoodDtrig, fGoodDmuon, fGoodCowboy;
  bool fGood, fGoodFake;
  bool fGoodTIS, fGoodTISFake;
  bool fGoodTISDT, fGoodTISDTDM, fGoodTISDTFake, fGoodTISDTDMFake;

  double fYield, fYieldE;


  // ----------------------------------------------------------------------
  ClassDef(plotFake,1)

};


#endif
