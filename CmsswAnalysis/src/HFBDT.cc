// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFBDT
// ---------------------
//
// 2017/05/08 Stephan Wiederkehr      added variables glbKinkFinderLog & Qprod
// 2016/09/27 Stephan Wiederkehr      updated variables
// 2016/08/25 Stephan Wiederkehr      first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Bmm/CmsswAnalysis/interface/HFBDT.hh"
#include <map>
#include <string>
#include <fstream>

#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

using namespace std;

HFBDT::HFBDT() {
  verbose_ = false;
  weightFile_ = edm::FileInPath();
  BDTreader = 0;
  isSetup = false;
  muon = new BDTmuon();
}

HFBDT::HFBDT(edm::FileInPath weightFile, bool verbose = false): verbose_(verbose),weightFile_(weightFile) {
  BDTreader = 0;
  isSetup = false;
  muon = new BDTmuon();
  ifstream wfile(weightFile_.fullPath());
  if (!wfile.is_open())
    {cout << "Could not open weight file: " << weightFile_.fullPath() << endl;}
}

HFBDT::~HFBDT(){
  //delete can lead to segfault: is the muon being deleted before this 
  //destructor is called?

  //delete muon;
}

void HFBDT::addVarToReader(TMVA::Reader* reader,string str) {

  if (isSetup)
    {cout << "ERROR: The reader is already set up." << endl;return;}

  cout << " Adding variable " << str << " to the reader." << endl;
  // cout << "value: " << muon->getVar(str) 
  //      << "(" << *(muon->getPtr(str)) << ")" << endl;
  //float *variable = muon->getPtr(str);
  reader->AddVariable(str,muon->getPtr(str));
}

void HFBDT::addSpecToReader(TMVA::Reader* reader,string str) {

  if (isSetup)
    {cout << "ERROR: The reader is already set up." << endl;return;}

  cout << " Adding spectator " << str << " to the reader (dummy value)." << endl;
  // cout << "value: " << muon->getVar(str) 
  //      << "(" << *(muon->getPtr(str)) << ")" << endl;
  //float *variable = muon->getPtr(str);
  reader->AddSpectator(str,muon->getSpecDummy());
}

void HFBDT::setupReader() {
  //reads the first MAX_CHARS characters of a line until line MAX_LINE of the
  //weight file and adds the found variables to the reader.
 
  if (verbose_) {cout << "Setting up the TMVAreader." << endl;}
  //Reader needs to be replaced. Overwriting can lead to bad states.
  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  muon->createVariableMap();

  ifstream fileContent;
  int NumberOfVariables = -1;
  if (verbose_) {cout << "The weight file: " << weightFile_.fullPath() << endl;}
  fileContent.open(weightFile_.fullPath());
  if (!fileContent.good())
    {
      cout << "Could not open file: " << weightFile_.fullPath() << endl;
      return;
    }
  const unsigned int MAX_CHARS = 100;
  //the lineCounter is a backup/debugging tool
  const unsigned int MAX_LINES = 20000;
  unsigned int lineCounter=0; 
  char buffer[MAX_CHARS];

  while (!fileContent.getline(buffer,MAX_CHARS,'\n').eof())
    {
      //if getline does not find the delimiter it sets a fail bit stopping 
      //all input. -->Clear the failbit and continue to parse
      fileContent.clear();
      string line = string(buffer);
      //cout << "line " << lineCounter << ": " << line << endl; //debugging
      //find number of variables
      if (line.find("Variables NVar") != string::npos)
	{
	  NumberOfVariables = atoi((line.substr(line.find("\"")+1,line.rfind("\"")-line.find("\"")-1)).c_str());
	  if (NumberOfVariables == -1) {break;}
	  cout << "Found "
	       <<  NumberOfVariables
	       << " variables." << endl;
	}
      //find variable name
      if (line.find("Variable VarIndex") != string::npos)
	{
	  string variableName = line.substr(line.find("Expression=\"")+12,line.rfind("Label")-line.find("Expression=\"")-14);
	  cout << "Found variable: '" << variableName << "'" << endl;
	  addVarToReader(reader, variableName);
	}
     if (line.find("Spectator SpecIndex") != string::npos)
	{
	  string spectatorName = line.substr(line.find("Expression=\"")+12,line.rfind("Label")-line.find("Expression=\"")-14);
	  cout << "Found spectator: '" << spectatorName << "'" << endl;
	  addSpecToReader(reader, spectatorName);
	}
      //stop search
      if (line.find("</Variables") != string::npos)
	{break;}

      //debug
      if (lineCounter==MAX_LINES) {break;}
      else {lineCounter++;}
    }
  if (verbose_) {cout << "Booking the Method." << endl;}
  BDTreader = reader;
  BDTreader->BookMVA("BDT",weightFile_.fullPath());
  isSetup = true;
}

double HFBDT::evaluate() {

  if (BDTreader == 0) return -42;
  if (isSetup == false)
    {cout << "WARNING: The BDT was not set up correctly beforehand." << endl;}
  if (!muon->varsSet())
    {cout << "WARNING: The values are not up to date, e.g. old or not set." << endl;}
  if ( !muon->areVarsValid() )
    {return -1;}

  muon->unsetVars();
  return BDTreader->EvaluateMVA("BDT");
}


////////////////////////////////////////////////////////////////////////////

BDTmuon::BDTmuon() {
  mapIsSet = false;
  varsAreSet_ = false;
  varsValid_ = false;
  pt=-42;
  eta=-42;
  deltaR=-42;
  gNchi2=-42;
  vMuHits=-42;
  mMuStations=-42;
  dxyRef=-42;
  dzRef=-42;
  LWH=-42;
  valPixHits=-42;
  innerChi2=-42;
  outerChi2=-42;
  iValFrac=-42;
  segComp=-42;
  chi2LocMom=-42;
  chi2LocPos=-42;
  glbTrackTailProb=-42;
  NTrkVHits=-42;
  kinkFinder=-42;
  vRPChits=-42;
  vDThits=-42;
  vCSChits=-42;
  glbKinkFinder=-42;
  staRelChi2=-42;
  trkRelChi2=-42;
  glbDeltaEtaPhi=-42;
  timeAtIpInOut=-42;
  timeAtIpInOutErr=-42;
  vDThits_1=-42;
  vDThits_2=-42;
  vDThits_3=-42;
  vDThits_4=-42;
  vRPChits_1=-42;
  vRPChits_2=-42;
  vRPChits_3=-42;
  vRPChits_4=-42;
  vCSChits_1=-42;
  vCSChits_2=-42;
  vCSChits_3=-42;
  vCSChits_4=-42;
  vMuonHitComb=-42;
  Qprod=-42;
  spectatorDummy=-42;
}

BDTmuon::~BDTmuon(){
}

void BDTmuon::createVariableMap() {
  //creates a map to link the strings from the weight file to the variables
  //in the BDTmuon struct to be inserted in the TMVA reader.
  map<string,float*> tmpMap;
  tmpMap["pt"] = &pt;
  tmpMap["eta"] = &eta;
  tmpMap["deltaR"] = &deltaR;
  tmpMap["gNchi2"] = &gNchi2;
  tmpMap["vMuHits"] = &vMuHits;
  tmpMap["mMuStations"] = &mMuStations;
  tmpMap["dxyRef"] = &dxyRef;
  tmpMap["dzRef"] = &dzRef;
  tmpMap["LWH"] = &LWH;
  tmpMap["valPixHits"] = &valPixHits;
  tmpMap["innerChi2"] = &innerChi2;
  tmpMap["outerChi2"] = &outerChi2;
  tmpMap["iValFrac"] = &iValFrac;
  tmpMap["segComp"] = &segComp;
  tmpMap["chi2LocMom"] = &chi2LocMom;
  tmpMap["chi2LocPos"] = &chi2LocPos;
  tmpMap["glbTrackTailProb"] = &glbTrackTailProb;
  tmpMap["NTrkVHits"] = &NTrkVHits;
  tmpMap["kinkFinder"] = &kinkFinder;
  tmpMap["vRPChits"] = &vRPChits;
  tmpMap["vDThits"] = &vDThits;
  tmpMap["vCSChits"] = &vCSChits;
  tmpMap["glbKinkFinder"] = &glbKinkFinder;
  tmpMap["TMath::Log(2+glbKinkFinder)"] = &glbKinkFinderLog;
  tmpMap["staRelChi2"] = &staRelChi2;
  tmpMap["trkRelChi2"] = &trkRelChi2;
  tmpMap["glbDeltaEtaPhi"] = &glbDeltaEtaPhi;
  tmpMap["timeAtIpInOut"] = &timeAtIpInOut;
  tmpMap["timeAtIpInOutErr"] = &timeAtIpInOutErr;
  tmpMap["vDT_1"] = &vDThits_1;
  tmpMap["vDT_2"] = &vDThits_2;
  tmpMap["vDT_3"] = &vDThits_3;
  tmpMap["vDT_4"] = &vDThits_4;
  tmpMap["vRPC_1"] = &vRPChits_1;
  tmpMap["vRPC_2"] = &vRPChits_2;
  tmpMap["vRPC_3"] = &vRPChits_3;
  tmpMap["vRPC_4"] = &vRPChits_4;
  tmpMap["vCSC_1"] = &vCSChits_1;
  tmpMap["vCSC_2"] = &vCSChits_2;
  tmpMap["vCSC_3"] = &vCSChits_3;
  tmpMap["vCSC_4"] = &vCSChits_4;
  tmpMap["vMuonHitComb"] = &vMuonHitComb;
  tmpMap["STATrkMult150"] = &STATrkMult_150;
  tmpMap["TMTrkMult100"] = &TMTrkMult_100;
  tmpMap["Qprod"] = &Qprod;

  variableMap = tmpMap;
  mapIsSet=true;
}

float BDTmuon::getVar(string str) {
  if (!mapIsSet)
    {cout << "ERROR: Variable map is not set." << endl;return -42;}
  else
    {
      if (variableMap.find(str) != variableMap.end())
	{return *(variableMap[str]);}
      else
	{throw MapError("The variable '"+str+"' is not in the variable map.");}
    }
}

float* BDTmuon::getPtr(string str) {
  if (!mapIsSet)
    {cout << "ERROR: Variable map is not set." << endl;return 0;}
  else
    {
      if (variableMap.find(str) != variableMap.end())
	{return variableMap[str];}
      else
	{throw MapError("The variable '"+str+"' is not in the variable map.");}
    }
}

void BDTmuon::fillBDTmuon(const reco::Muon& recoMuon, const reco::VertexCollection* vc, reco::BeamSpot* bs, int staTrkMult_150, int tmTrkMult_100) {

  reco::TrackRef gTrack = recoMuon.globalTrack();
  reco::TrackRef iTrack = recoMuon.innerTrack();
  reco::TrackRef oTrack = recoMuon.outerTrack();
  const reco::HitPattern gHits = gTrack->hitPattern();
  const reco::HitPattern iHits = iTrack->hitPattern();
  const reco::MuonQuality muonQuality = recoMuon.combinedQuality();
  if ( iTrack.isNonnull() && oTrack.isNonnull() && gTrack.isNonnull() )
    {varsValid_ = true;}
  else {return;}

  int pvIndex = getPv(iTrack.index(),vc); //HFDumpUtitilies
  math::XYZPoint refPoint;
  if (pvIndex > -1)
    {refPoint = vc->at(pvIndex).position();}
  else
    {
      if (bs)
	{refPoint = bs->position();}
      else
	{cout << "ERROR: No beam sport found!" << endl;}
    }

  pt = iTrack->pt();
  eta = iTrack->eta();
  deltaR = getDeltaR(*iTrack,*oTrack);
  gNchi2 = gTrack->normalizedChi2();
  vMuHits = gHits.numberOfValidMuonHits();
  mMuStations = recoMuon.numberOfMatchedStations();
  dxyRef = iTrack->dxy(refPoint);
  dzRef = iTrack->dz(refPoint);
  LWH = iHits.trackerLayersWithMeasurement();
  valPixHits = iHits.numberOfValidPixelHits();
  innerChi2 = iTrack->normalizedChi2();
  outerChi2 = oTrack->normalizedChi2();
  iValFrac = iTrack->validFraction();
  segComp = muon::segmentCompatibility(recoMuon);
  chi2LocMom = muonQuality.chi2LocalMomentum;
  chi2LocPos = muonQuality.chi2LocalPosition;
  glbTrackTailProb = muonQuality.glbTrackProbability;
  NTrkVHits = iHits.numberOfValidTrackerHits();
  kinkFinder = muonQuality.trkKink;
  vRPChits = gHits.numberOfValidMuonRPCHits();
  glbKinkFinder = muonQuality.glbKink;
  glbKinkFinderLog = TMath::Log(2+muonQuality.glbKink);
  staRelChi2 = muonQuality.staRelChi2;
  glbDeltaEtaPhi = muonQuality.globalDeltaEtaPhi;
  trkRelChi2 = muonQuality.trkRelChi2;
  vDThits = gHits.numberOfValidMuonDTHits();
  vCSChits = gHits.numberOfValidMuonCSCHits();
  timeAtIpInOut = recoMuon.time().timeAtIpInOut;
  timeAtIpInOutErr = recoMuon.time().timeAtIpInOutErr;
  getMuonHitsPerStation(gTrack); //also fills vMuonHitComb
  STATrkMult_150 = staTrkMult_150;
  TMTrkMult_100 = tmTrkMult_100;
  Qprod = (iTrack->charge() * oTrack->charge());

  varsAreSet_ = true;
}

void BDTmuon::getMuonHitsPerStation(const reco::TrackRef gTrack) {

  unsigned int dt1(0),dt2(0),dt3(0),dt4(0);
  unsigned int rpc1(0),rpc2(0),rpc3(0),rpc4(0);
  unsigned int csc1(0),csc2(0),csc3(0),csc4(0);
  float comb(0);
  const reco::HitPattern &pattern = gTrack->hitPattern();
  for (int i=0;i<pattern.numberOfAllHits(reco::HitPattern::TRACK_HITS);i++)
    {
      uint32_t hit = pattern.getHitPattern(reco::HitPattern::TRACK_HITS,i);
      if (pattern.validHitFilter(hit) != 1) {continue;}
      if (pattern.getMuonStation(hit) == 1)
	{
	  if (pattern.muonDTHitFilter(hit)) {dt1++;}
	  if (pattern.muonRPCHitFilter(hit)) {rpc1++;}
	  if (pattern.muonCSCHitFilter(hit)) {csc1++;}
	}
      else if (pattern.getMuonStation(hit) == 2)
	{
	  if (pattern.muonDTHitFilter(hit)) {dt2++;}
	  if (pattern.muonRPCHitFilter(hit)) {rpc2++;}
	  if (pattern.muonCSCHitFilter(hit)) {csc2++;}
	}
      else if (pattern.getMuonStation(hit) == 3)
	{
	  if (pattern.muonDTHitFilter(hit)) {dt3++;}
	  if (pattern.muonRPCHitFilter(hit)) {rpc3++;}
	  if (pattern.muonCSCHitFilter(hit)) {csc3++;}
	}
      else if (pattern.getMuonStation(hit) == 4)
	{
	  if (pattern.muonDTHitFilter(hit)) {dt4++;}
	  if (pattern.muonRPCHitFilter(hit)) {rpc4++;}
	  if (pattern.muonCSCHitFilter(hit)) {csc4++;}
	}      
    }//for
  comb = (dt1+dt2+dt3+dt4)/2. + (rpc1+rpc2+rpc3+rpc4);
  csc1>6 ? comb+=6 : comb+=csc1;
  csc2>6 ? comb+=6 : comb+=csc2;
  csc3>6 ? comb+=6 : comb+=csc3;
  csc4>6 ? comb+=6 : comb+=csc4;
  //assignments
  vDThits_1 = dt1;
  vDThits_2 = dt2;
  vDThits_3 = dt3;
  vDThits_4 = dt4;
  vRPChits_1 = rpc1;
  vRPChits_2 = rpc2;
  vRPChits_3 = rpc3;
  vRPChits_4 = rpc4;
  vCSChits_1 = csc1;
  vCSChits_2 = csc2;
  vCSChits_3 = csc3;
  vCSChits_4 = csc4;
  vMuonHitComb = comb;
}

double getDeltaR(reco::Track track1,reco::Track track2) {

  DeltaR<reco::Track,reco::Track> DR_init;
  double DR = DR_init(track1,track2);
  return DR;
}



