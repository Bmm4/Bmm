#ifndef HFBDT_H
#define HFBDT_H

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TString.h"
#include "TMVA/Reader.h"
#include <string>
#include <map>
#include <exception>
#include <vector> //debug

#define _USE_MATH_DEFINES

class BDTmuon {

public:
  BDTmuon();
  ~BDTmuon();

  void createVariableMap();
  void fillBDTmuon(const reco::Muon& recoMu, const reco::VertexCollection* vc, reco::BeamSpot* bs, int STATrkMult, int TMTrkMult);
  void getMuonHitsPerStation(const reco::TrackRef gTrack);
  float getVar(std::string);
  float* getPtr(std::string);
  float* getSpecDummy() {return &spectatorDummy;}
  bool varsSet() {return varsAreSet;}
  void unsetVars() {varsAreSet = false;}

private:
  float pt, eta;
  float deltaR,gNchi2,vMuHits,mMuStations,dxyRef,dzRef,LWH,valPixHits;
  float innerChi2,outerChi2,iValFrac,segComp,chi2LocMom,chi2LocPos;
  float glbTrackTailProb,NTrkVHits,kinkFinder,vRPChits,vDThits,vCSChits;
  float glbKinkFinder,glbKinkFinderLog,staRelChi2,trkRelChi2,glbDeltaEtaPhi;
  float inverseBeta,inverseBetaErr,timeAtIpInOut,timeAtIpInOutErr;
  //per station (1-4)
  float vDThits_1,vDThits_2,vDThits_3,vDThits_4;
  float vRPChits_1,vRPChits_2,vRPChits_3,vRPChits_4;
  float vCSChits_1,vCSChits_2,vCSChits_3,vCSChits_4;
  float vMuonHitComb,Qprod;
  //stand alone muon && tracker muon multiplicity in a radius given in cm
  float STATrkMult_150,TMTrkMult_100; 
  float spectatorDummy;
  
  std::map<std::string,float*> variableMap; //default constructor fine?
  bool mapIsSet;
  bool varsAreSet;
};

//By default this class does not create any output in addition 
//to the TMVA output. 
class HFBDT {

public:
  HFBDT();
  HFBDT(edm::FileInPath wf, bool verbose);
  ~HFBDT();

  void setupReader();
  double evaluate();
  void setWeightFile(edm::FileInPath weightFile) {weightFile_ = weightFile;}
  void setVerbose(bool verbose) {verbose_ = verbose;}
  BDTmuon* getMuon() {return muon;}

private:
  void addVarToReader(TMVA::Reader* reader,std::string str);
  void addSpecToReader(TMVA::Reader* reader,std::string str);

  bool verbose_;
  edm::FileInPath weightFile_;
  TMVA::Reader *BDTreader;
  bool isSetup;

  BDTmuon *muon;
  
};

class MapError: public std::exception {

public:
  MapError(const std::string& message):message_(message){};

  virtual const char* what() const throw()
  {
    return message_.c_str();
  }

private:
  std::string message_;

};

double getDeltaR(reco::Track track1,reco::Track track2);

#endif
