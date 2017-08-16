#ifndef CANDANA_H
#define CANDANA_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>

#include <TH1.h>
#include <TH2.h>

#include "bmmReader.hh"
#include "common/AnalysisCuts.hh"

#include "redTreeData.hh"
#include "preselection.hh"
#include "ReaderData.hh"
#include "cuts.hh"


#define NTRGMAX 40

struct isoNumbers {
  double iso;
  int    pvTracks;
  int    clTracks;
  int    Tracks;
};

// -- TMVA related
#include "TMVA/Reader.h"

struct muonData {
  float mbdt;
  // -- tight muon variables
  float validMuonHits, glbNChi2;
  int   nMatchedStations, validPixelHits, validPixelHits2, trkLayerWithHits;
  // -- Luca's additional variables
  float trkValidFract;
  float pt, eta;
  float segComp, chi2LocMom, chi2LocPos, glbTrackProb;
  float NTrkVHits, NTrkEHitsOut;
  // -- Mario's and other additional variables
  float kink;
  float dpt, dptrel, deta, dphi, dr;
};

struct mvaMuonIDData {
  // -- Luca's original setup
  float trkValidFract, glbNChi2;
  float pt, eta;
  float segComp, chi2LocMom, chi2LocPos, glbTrackProb;
  float NTrkVHits, NTrkEHitsOut;
  // -- Stephan's other variables
  float glbTrackTailProb, glbDeltaEtaPhi, iValFrac;
  float LWH; //int
  float dxyRef, dzRef;
  float kinkFinder, glbKinkFinder, glbKinkFinderLOG, timeAtIpInOutErr, outerChi2;
  float valPixHits, TMTrkMult100; //int
  float innerChi2, trkRelChi2;
  float vMuonHitComb;
  float Qprod; //int
  float spectatorDummy;
};

class TTree;
class TDirectory;

// ----------------------------------------------------------------------
class candAna {

public:
  candAna(bmmReader *pReader, std::string name, std::string cutsFile);
  virtual ~candAna();

  virtual void        evtAnalysis(TAna01Event *evt);
  virtual void        candAnalysis();
  virtual void        candEvaluation();
  virtual void        endAnalysis();
  virtual void        efficiencyCalculation();
  virtual void        setupReducedTree(TTree *);
  virtual void        setupMuonIdTree(TTree *);

  virtual bool        tis(TAnaCand *pC);
  virtual bool        tos(TAnaCand *pC);
  virtual bool        refTrigger(TAnaCand *pC, std::string refTriggerPath);
  virtual int         matchTrgObj2Trk(TVector3 t, double &dr);
  //  virtual double      dist2PdTrigger(const TVector3 &t);
  virtual void        dist2PdTrigger(TSimpleTrack *pS, double &dr, double &dm1, double &dm2);
  virtual bool        triggerFired(std::string triggerPath);
  virtual void        nTriggers();
  virtual void        triggerL1T();
  virtual void        triggerHLT();
  virtual int         nearestPV(int pvIdx, double maxDist = 99.);
  virtual void        getSigTracks(std::vector<int> &v, TAnaCand *pC);
  virtual double      constrainedMass();
  virtual void        muScaleCorrectedMasses();
  virtual void        pvStudy(bool bookHist = false);
  virtual void        runRange();
  virtual void        genMatch();
  virtual void        genMatchOld();
  virtual void        recoMatch();
  virtual void        candMatch();
  virtual bool        checkDataCand(TAnaCand *);
  virtual void        triggerSelection();
  virtual void        fillCandidateHistograms(int offset);
  virtual void        fillRedTreeData();
  virtual void        fillMuonData(muonData &a, TAnaMuon *pM);
  virtual void        replaceAll(std::string &s, std::string a, std::string b);

  // virtual TMVA::Reader* setupReader(std::string xmlFile, readerData &rd);
  virtual TMVA::Reader* setupMuonMvaReader(std::string xmlFile, mvaMuonIDData &rd);
  virtual void        calcBDT();
  virtual int         detChan(double m1eta, double m2eta);

  virtual void        bookHist();

  virtual void        basicCuts();
  virtual void        moreBasicCuts();
  virtual void        candidateCuts();
  virtual void        moreCandidateCuts();

  virtual void        readCuts(std::string fileName, int dump = 1);
  virtual void        readFile(std::string fileName, std::vector<std::string> &lines);

  virtual bool        highPurity(TAnaTrack *pt);
  virtual bool        tightMuon(TAnaTrack *pt, bool hadronsPass = true);
  virtual bool        tightMuon(TSimpleTrack *pt, bool hadronsPass = true);

  virtual bool        mvaMuon(TAnaMuon *pt, double &result, bool hadronsPass = true);
  virtual bool        mvaMuon(TSimpleTrack *pt, double &result, bool hadronsPass = true);
  virtual bool        mvaMuon(TAnaTrack *pt, double &result, bool hadronsPass = true);
  double              getDetVarComb(TAnaMuon *mu);
  bool                mvaMuonPassedPreselection(mvaMuonIDData mu);

  virtual std::string splitTrigRange(std::string tl, int &r1, int &r2);

  virtual double      isoClassicWithDOCA(TAnaCand*, double dca, double r = 0.7, double ptmin = 0.9);
  virtual std::pair<int, int> nCloseTracks(TAnaCand*, double dca, double dcaS, double pt = 0.5);
  virtual double      isoMuon(TAnaCand *, TAnaMuon *);
  virtual double      isoTrack(TAnaCand *pC, TAnaTrack *pTrack, double docaCut = 0.1, double r = 0.5, double pmin = 0.5);
  virtual void        xpDistMuons();
  virtual void        findAllTrackIndices(TAnaCand* pCand, std::map<int,int> *indices);

  virtual TAnaCand*   osCand(TAnaCand *pC);
  virtual double      osIsolation(TAnaCand *pC, double r = 1.0, double ptmin = 0.9);
  virtual int         osMuon(TAnaCand *pC, double r = 1.0);
  virtual void        boostGames();
  virtual double      matchToMuon(TAnaTrack *pt, bool skipSame = false); // match a single track to ALL muons
  virtual double      distToMuon(TSimpleTrack *pt); // calculate DR separation of a track to closest (other) muon
  virtual void        play();
  virtual void        play2();
  virtual void        play3();
  // match the 2 muons from the dimuon to HLT
  virtual bool        doTriggerMatching(TAnaTrack *pt1, TAnaTrack *pt2);
  // match a single track to HLT
  virtual bool        doTriggerMatching(TAnaTrack *pt, bool anyTrig = false,
					bool muonsOnly=true, bool anyModule=false);
  // To return the full deltaR not just a bool
  virtual double      doTriggerMatchingR(TAnaTrack *pt, bool anyTrig = false,
                                        bool muonsOnly=true, bool anyModule=false);
  // Veto track too close to a trigger object
  virtual bool        doTriggerVeto(TAnaTrack *pt, bool singleMatch,
				    bool matchPt, bool anyModule, float drCut, int histoOffset=0);
  virtual void        print1();
  virtual void        printCuts(ostream &OUT);



  std::string fName;
  std::string fCutFile;
  TDirectory *fHistDir;
  bmmReader *fpReader;
  TTree *fTree;
  TAna01Event *fpEvt;
  TAnaCand *fpCand, *fpOsCand;
  int fCandIdx;

  int fVerbose, fDbx;
  int fIsMC;

  Long64_t fRun, fEvt;
  int fLS;
  int fEvent;
  int fRunRange;
  int fYear;
  double fLumi;

  double       MASSMIN,   MASSMAX;
  double       SIGBOXMIN, SIGBOXMAX;
  double       BGLBOXMIN, BGLBOXMAX;
  double       BGHBOXMIN, BGHBOXMAX;

  double
  CANDPTLO, CANDETALO, CANDETAHI
    , CANDCOSALPHA, CANDALPHA
    , CANDFLS3D, CANDFLSXY, CANDVTXCHI2
    , CANDISOLATION, CANDDOCATRK, CANDCLOSETRK
    , PVAVEW8, CANDLIP, CANDLIPS, CANDLIP2, CANDLIPS2
    , CANDDOCA, CANDIP, CANDIPS
    , TRACKPTLO, TRACKPTHI, TRACKETALO, TRACKETAHI
    , TRACKTIP, TRACKLIP
    , MUPTLO, MUPTHI
    , MUETALO, MUETAHI, MUIP
    , MUBDT, MUBDT1
    ;

  int BLIND, TYPE, SELMODE, MUIDMASK, MUIDRESULT, TRACKQUALITY, TRUTHCAND, DATACAND, IGNORETRIGGER, NOPRESELECTION;

  std::map<std::string, pair<int, int> > HLTRANGE, NTRIGGERS;
  string DSNAME; // data set name

  bool fBarrel, fWideMass;
  AnalysisCuts fAnaCuts;
  std::vector<cuts*> fCuts;

  // -- TM
  int     fGenBTmi;
  int     fGenM1Tmi, fGenM2Tmi, fNGenPhotons;
  int     fRecM1Tmi, fRecM2Tmi;
  int     fCandTmi;
  int     fGenBpartial;
  int     fProcessType;
  double  fGenLifeTime, fGenMass, fGenFl3d;

  // -- BDT
  std::vector<TMVA::Reader*> fReaderEvents0;
  std::vector<TMVA::Reader*> fReaderEvents1;
  std::vector<TMVA::Reader*> fReaderEvents2;
  double  fBDT;
  ReaderData frd;

  TMVA::Reader *fMvaMuonID, *fMvaMuonID1;
  mvaMuonIDData mrd, mrd1;
  double  fMu1BDT, fMu2BDT, fMu1rBDT, fMu2rBDT;

  muonData fMuonData;
  muonData fMu1Data, fMu2Data;

  // -- variables for reduced tree, they are from fpCand
  bool    fJSON, fCowboy, fCandDcand;
  int     fCandTM, fCandType;
  int     fMu1TkQuality, fMu2TkQuality, fMu1Q, fMu2Q, fCandQ, fMu1PV, fMu2PV;
  bool    fMu1Id, fMu2Id, fHLTmatch, fMu1MvaId, fMu2MvaId, fMu1GmId, fMu2GmId;
  bool    fMu1rMvaId, fMu2rMvaId;
  double  fMu1TrigM, fMu2TrigM;

  double  fDeltaR, fMuDist, fMuDeltaR, fMu1Chi2, fMu2Chi2;
  double  fHltMu1Pt, fHltMu1Eta, fHltMu1Phi, fHltMu2Pt, fHltMu2Eta, fHltMu2Phi;
  double  fMu1Pt, fMu1Eta, fMu1Phi, fMu2Pt, fMu2Eta, fMu2Phi;
  double  fMu1PtGen, fMu2PtGen, fMu1EtaGen, fMu2EtaGen, fMu1PhiGen, fMu2PhiGen;
  int     fMu1GenID, fMu2GenID;
  double  fMu1PtNrf, fMu2PtNrf, fMu1EtaNrf, fMu2EtaNrf, fMu1PhiNrf, fMu2PhiNrf; // "not refitted"
  int     fMu1TrkLayer, fMu1Pix, fMu1BPix, fMu1BPixL1, fMu2TrkLayer, fMu2Pix, fMu2BPix, fMu2BPixL1;
  double  fMu1W8Mu, fMu1W8Tr, fMu2W8Mu, fMu2W8Tr;
  double  fMu1Iso, fMu2Iso, fMu1Isold, fMu2Isold;
  double  fMu1VtxProb, fMu2VtxProb;
  bool    fMu1OtherVtx, fMu2OtherVtx;
  double  fMu1XpDist, fMu2XpDist;
  double  fPvX, fPvY, fPvZ, fPvNdof, fPvAveW8, fPv2Ndof, fPv2AveW8, fPvDzmin, fPvDz12;
  int     fPvN, fPvIdx, fPv2Idx, fPvNtrk, fPv2Ntrk;
  double  fCandPt, fCandP, fCandTau, fCandTauE, fCandTauxy, fCandTauxyE, fCandEta, fCandPhi, fCandM, fCandME, fCandM2, fCandM3, fCandM4, fCandW8Tr, fCandW8Mu;
  double  fCandCosA, fCandA;
  double  fCandChi2, fCandDof, fCandChi2Dof, fCandProb, fCandFL3d, fCandFL3dE, fCandFLS3d, fCandFLxy, fCandFLSxy, fCandDoca;
  double  fCandIso;
  int     fCandIsoTrk, fCandPvTrk, fCandI0trk, fCandI1trk, fCandI2trk;
  int     fCandCloseTrk, fCandCloseTrkS1, fCandCloseTrkS2, fCandCloseTrkS3, fCandCloseTrkS4, fCandCloseTrkS5;
  double  fCandDocaTrk, fCandDocaTrkBdt, fMu1IPS, fMu1IP, fMu1IPE, fMu2IPS, fMu2IP, fMu2IPE;
  double  fCandPvTip, fCandPvTipE, fCandPvTipS, fCandPvLip, fCandPvLipE, fCandPvLipS, fCandPvIp, fCandPvIpE, fCandPvIpS;
  double  fCandPvIp3D, fCandPvIpE3D, fCandPvIpS3D;
  double  fCandPv2Lip, fCandPv2LipS, fCandPv12Lip, fCandPv12LipE, fCandPv12LipS;
  double  fCandPvDeltaChi2;
  double  fCandOtherVtx;
  // NO: They are filled into the "Jpsi" variables in the derived classes!

  double  fOsMuonPt, fOsMuonPtRel, fOsIso, fOsRelIso, fOsMuonDeltaR;

  int     fChan, fNchan;

  // -- another reduced tree
  TTree       *fEffTree;
  bool fETm1gt, fETm2gt, fETm1id, fETm1mvaid, fETm2id, fETm2mvaid;
  int fETm1q, fETm2q, fETchan;
  float fETgm, fETgpt, fETgeta;
  float fETcandMass;
  float fETm1pt, fETm1eta, fETg1pt, fETg1eta;
  float fETm2pt, fETm2eta, fETg2pt, fETg2eta;
  float fETgtau, fETtau;
  bool  fETdcand;

  TTree       *fMuonIdTree, *fPvStudyTree;

  TAnaTrack *fpMuon1, *fpMuon2;

  // -- isolation study
  isoNumbers fIsoR03Pt03, fIsoR03Pt05, fIsoR03Pt07, fIsoR03Pt09, fIsoR03Pt11;
  isoNumbers fIsoR05Pt03, fIsoR05Pt05, fIsoR05Pt07, fIsoR05Pt09, fIsoR05Pt11;
  isoNumbers fIsoR07Pt03, fIsoR07Pt05, fIsoR07Pt07, fIsoR07Pt09, fIsoR07Pt11;
  isoNumbers fIsoR09Pt03, fIsoR09Pt05, fIsoR09Pt07, fIsoR09Pt09, fIsoR09Pt11;
  isoNumbers fIsoR10Pt03, fIsoR10Pt05, fIsoR10Pt07, fIsoR10Pt09, fIsoR10Pt11;
  isoNumbers fIsoR11Pt03, fIsoR11Pt05, fIsoR11Pt07, fIsoR11Pt09, fIsoR11Pt11;

  string  fHLTPath;
  int     fL1Seeds, fHltPrescale;
  bool    fGoodCNC;
  bool    fGoodHLT, fGoodMuonsID, fGoodMuonsGmID, fGoodMuonsMvaID, fGoodMuonsPt, fGoodMuonsEta, fGoodTracks, fGoodTracksPt, fGoodTracksEta;
  bool    fGoodPvAveW8, fGoodPvLip, fGoodPvLipS, fGoodPv2Lip, fGoodPv2LipS, fGoodMaxDoca, fGoodIp, fGoodIpS, fGoodLip, fGoodLipS;
  bool    fGoodQ, fGoodPt, fGoodEta, fGoodCosA, fGoodAlpha, fGoodIso, fGoodM1Iso, fGoodM2Iso, fGoodChi2, fGoodFLS;
  bool    fGoodCloseTrack, fGoodCloseTrackS1, fGoodCloseTrackS2, fGoodCloseTrackS3;
  bool    fGoodDocaTrk, fGoodLastCut, fTIS, fTOS, fL1T, fRefTrigger, fGoodAcceptance, fGoodJpsiCuts;
  bool    fPreselAlpha, fPreselFLS, fPreselOther;
  // -- the following are for UL's trigger selection and matching (as a cross check for DK's original)
  bool    fGoodHLT1;
  string  fHLT1Path, fL1SeedString;
  // --
  double  fHltD1, fHltD2;

  bool    fPreselection;
  bool    fBadEvent;
  int     fhltType; // to hold the HLT information d.k.
  double  fTrigMatchDeltaPt;
  map<unsigned int, unsigned int, less<unsigned int> > hltObjMap;

  int fNtrg, fNtrgPs[NTRGMAX], fNtrgTos[NTRGMAX];

  struct redTreeData fRTD;
};

#endif //  CANDANA_H
