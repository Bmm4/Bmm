#include "TAna01Event.hh"
#include <iostream>
#include <stdexcept>

ClassImp(TAna01Event)
 
using namespace std;

#define NGEN 1000
#define NGENCAND 100000
#define NRECTRACK 50000
#define NSIMPLETRACK 50000
#define NMUON 1000
#define NTRGOBJ 1000
#define NSIGTRACK 500000
#define NCALOJET 1000
#define NPFJET 1000
#define NGENJET 1000
#define NTRACKJET 1000
#define NCANDIDATES 200000
#define NVERTEX 1000

// ----------------------------------------------------------------------
TAna01Event::TAna01Event(Int_t Option) {
  fGenCands        = new TClonesArray("TGenCand", NGENCAND);
  fnGenCands       = 0;
 
  fGenT            = new TClonesArray("TGenCand", NGEN);
  fnGenT           = 0;

  fRecTracks       = new TClonesArray("TAnaTrack", NRECTRACK);
  fnRecTracks      = 0;

  fSimpleTracks    = new TClonesArray("TSimpleTrack", NSIMPLETRACK);
  fnSimpleTracks   = 0;

  fMuons           = new TClonesArray("TAnaMuon", NMUON);
  fnMuons          = 0;

  fTrgObj          = new TClonesArray("TTrgObj", NTRGOBJ);
  fnTrgObj         = 0;

  fSigTracks       = new TClonesArray("TAnaTrack", NSIGTRACK);
  fnSigTracks      = 0;

  fCaloJets        = new TClonesArray("TAnaJet", NCALOJET);
  fnCaloJets       = 0;

  fGenJets         = new TClonesArray("TAnaJet", NGENJET);
  fnGenJets        = 0;

  fTrackJets       = new TClonesArray("TAnaJet", NTRACKJET);
  fnTrackJets      = 0;

  fPFJets          = new TClonesArray("TAnaJet", NPFJET);
  fnPFJets         = 0;

  fCandidates      = new TClonesArray("TAnaCand", NCANDIDATES);
  fnCandidates     = 0;

  fPV              = new TClonesArray("TAnaVertex", NVERTEX);
  fnPV             = 0;

  Clear();

}

// ----------------------------------------------------------------------
void TAna01Event::Clear(Option_t *option) {

  fLumi = fLumiInt = fLumiSection = fBx = fOrbit = -9999; 
  fTimeLo = fTimeHi = 0; 

  for (int i = 0; i < NL1T; ++i) {
    fL1TPrescale[i] = 0;
    fL1TResult[i] = fL1TMask[i] = fL1TError[i] = false; 
  }

  for (int i = 0; i < NHLT; ++i) {
    fHLTPrescale[i] = 0; 
    fHLTResult[i] = fHLTWasRun[i] = fHLTError[i] = false; 
  }


  TAnaTrack *pTrack;
  for (int i = 0; i < fnRecTracks; i++) {
    pTrack = getRecTrack(i);
    pTrack->clear();
  }
  fRecTracks->Clear(option);
  fnRecTracks = 0;


  TSimpleTrack *sTrack;
  for (int i = 0; i < fnSimpleTracks; i++) {
    sTrack = getSimpleTrack(i);
    sTrack->clear();
  }
  fSimpleTracks->Clear(option);
  fnSimpleTracks = 0;

  TAnaMuon *pMuon;
  for (int i = 0; i < fnMuons; i++) {
    pMuon = getMuon(i);
    pMuon->clear();
    pMuon->~TAnaMuon();
  }
  fMuons->Clear(option);
  fnMuons = 0;

  TTrgObj *pTrgObj;
  for (int i = 0; i < fnTrgObj; i++) {
    pTrgObj = getTrgObj(i);
    pTrgObj->clear();
  }
  fTrgObj->Clear(option);
  fnTrgObj = 0;

  for (int i = 0; i < NL1T; ++i) {
    fL1TPrescale[i] = 0; 
    fL1TResult[i] = fL1TMask[i] = fL1TError[i] = false; 
  }

  for (int i = 0; i < NLTT; ++i) {
    fLTTPrescale[i] = 0; 
    fLTTResult[i] = fLTTMask[i] = fLTTError[i] = false; 
  }

  for (int i = 0; i < NHLT; ++i) {
    fHLTPrescale[i] = 0; 
    fHLTResult[i] = fHLTWasRun[i] = fHLTError[i] = false; 
  }


  TAnaTrack *pSigTrack;
  for (int i = 0; i < fnSigTracks; i++) {
    pSigTrack = getSigTrack(i);
    pSigTrack->clear();
  }
  fSigTracks->Clear(option);
  fnSigTracks = 0;

  TAnaJet *pCaloJet;
  for (int i = 0; i < fnCaloJets; i++) {
    pCaloJet = getCaloJet(i);
    pCaloJet->clear();
  }
  fCaloJets->Clear(option);
  fnCaloJets = 0;

  TAnaJet *pGenJet;
  for (int i = 0; i < fnGenJets; i++) {
    pGenJet = getGenJet(i);
    pGenJet->clear();
  }
  fGenJets->Clear(option);
  fnGenJets = 0;

  TAnaJet *pTrackJet;
  for (int i = 0; i < fnTrackJets; i++) {
    pTrackJet = getTrackJet(i);
    pTrackJet->clear();
  }
  fTrackJets->Clear(option);
  fnTrackJets = 0;

  TAnaJet *pPFJet;
  for (int i = 0; i < fnPFJets; i++) {
    pPFJet = getPFJet(i);
    pPFJet->clear();
  }
  fPFJets->Clear(option);
  fnPFJets = 0;

  TAnaCand *pCand;
  for (int i = 0; i < fnCandidates; i++) {
    pCand = getCand(i);
    pCand->~TAnaCand(); // call destructor, not only clear(). Prevents memleak due to std::vector
  }
  fCandidates->Clear(option);
  fnCandidates = 0;

  clearGenBlock();

  TGenCand *pGenT;
  for (int i = 0; i < fnGenT; i++) {
    pGenT = getGenT(i);
    pGenT->clear();
  }
  fGenT->Clear(option);
  fnGenT = 0;

  TAnaVertex *pPV;
  for (int i = 0; i < fnPV; i++) {
    pPV = getPV(i);
    pPV->clear();
  }
  fPV->Clear(option);
  fnPV = 0;

}


// ----------------------------------------------------------------------
void TAna01Event::clearGenBlock() {
  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->clear();
  }
  fGenCands->Clear();
  fnGenCands = 0;

}


// ----------------------------------------------------------------------
TGenCand* TAna01Event::getGenCand(Int_t n) { 
  return (TGenCand*)fGenCands->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TGenCand* TAna01Event::addGenCand() {
  TClonesArray& d = *fGenCands; 
  new(d[d.GetLast()+1]) TGenCand();
  ++fnGenCands;
  return (TGenCand*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
TGenCand* TAna01Event::getGenT(Int_t n) { 
  return (TGenCand*)fGenT->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TGenCand *TAna01Event::getGenTWithIndex(Int_t n) {
  TGenCand *pG; 
  for (int i = 0; i < fnGenT; ++i) {
    pG = (TGenCand*)fGenT->UncheckedAt(i); 
    if (n == pG->fNumber) return pG; 
  }
  return 0; 
}

// ----------------------------------------------------------------------
TGenCand* TAna01Event::addGenT() {
  TClonesArray& d = *fGenT; 
  new(d[d.GetLast()+1]) TGenCand();
  ++fnGenT;
  return (TGenCand*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
void TAna01Event::fillGenT(int theParticleIdx) {
  
  const int verbose(0);
  
  set<int> daughters, mothers;
  TGenCand *pDau, *pTmp; 
  int iMom; 
  
  // -- daughters of theParticleIdx
  for (int id = theParticleIdx+1; id < nGenCands(); ++id) {
    pDau = getGenCand(id);
    iMom = pDau->fMom1;
    while (iMom > theParticleIdx) {
      pTmp = getGenCand(iMom);
      iMom = pTmp->fMom1;
    }
    if (iMom == theParticleIdx) {
      daughters.insert(id); 
    }
  }
  
  mothers.insert(theParticleIdx); 
  int cnt(0); 
  while (1) {
    ++cnt;
    for (set<int>::iterator i = mothers.begin(); i != mothers.end(); ++i) {
      pTmp = getGenCand(*i); 
      for (int id = pTmp->fMom1; id <= pTmp->fMom2; ++id) {
	mothers.insert(id);
      }
    }
    if (mothers.end() != find(mothers.begin(), mothers.end(), 0)) break;
    if (mothers.end() != find(mothers.begin(), mothers.end(), 1)) break;
    if (mothers.end() != find(mothers.begin(), mothers.end(), -1)) break;
  }

  // -- add documentation block
  for (int i = 0; i < nGenCands(); ++i) {
    pTmp = getGenCand(i); 
    if (pTmp->fStatus == 3) mothers.insert(i); 
  }

  // -- now combine all sets into one set
  set<int> allParticles;
  for (set<int>::iterator i = mothers.begin(); i != mothers.end(); ++i) {
    allParticles.insert(*i); 
  }

  for (set<int>::iterator i = daughters.begin(); i != daughters.end(); ++i) {
    allParticles.insert(*i); 
  }

  if (verbose) {
    cout << "all: " << allParticles.size() << endl;
    for (set<int>::iterator i = allParticles.begin(); i != allParticles.end(); ++i) {
      cout << "a: "; getGenCand(*i)->dump(); 
    }
  }

  // -- fill genT
  TGenCand *pG, *pOld; 
  bool alreadyFilled(false); 
  for (set<int>::iterator i = allParticles.begin(); i != allParticles.end(); ++i) {
    alreadyFilled = false; 
    for (int it = 0; it < nGenT(); ++it) {
      if (*i == getGenT(it)->fNumber) {
	alreadyFilled = true; 
	break;
      }
    }
    if (alreadyFilled) continue;

    pOld = getGenCand(*i); 
    pG = addGenT(); 
    pG->fID     = pOld->fID; 
    pG->fNumber = pOld->fNumber;
    pG->fStatus = pOld->fStatus;
    pG->fMom1   = pOld->fMom1; 
    pG->fMom2   = pOld->fMom2;
    pG->fDau1   = pOld->fDau1; 
    pG->fDau2   = pOld->fDau2;
    pG->fTag    = pOld->fTag;
    pG->fQ      = pOld->fQ;
    pG->fP      = pOld->fP; 
    pG->fMass   = pOld->fMass;
    pG->fV      = pOld->fV;
    pG->fTime   = pOld->fTime;
  }

}

// ----------------------------------------------------------------------
void TAna01Event::dumpGenBlock() {
  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->dump();
  }
}

// ----------------------------------------------------------------------
int TAna01Event::getGenIndex(double px, double py, double pz, int charge, double precision) {
  TGenCand *pGenCand;
  int index(-1);

  //   cout << "Searching among " << fnGenCands << " cands for match to p = (" 
  //        << px << ", " << py << ", " << pz << ") with ID = " << id << endl;

  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    if (pGenCand->fQ * charge < 0) continue;
    if ((TMath::Abs((pGenCand->fP.X() - px))/px) > precision) continue;
    if ((TMath::Abs((pGenCand->fP.Y() - py))/py) > precision) continue;
    if ((TMath::Abs((pGenCand->fP.Z() - pz))/pz) > precision) continue;
    
    index = i; 
    //     cout << "I think this is a match at index= " << index 
    // 	 << "  " << px << " " << py << " " << pz << " " << endl;
    //     pGenCand->dump();
    break;
  }
  return index;
}

// ----------------------------------------------------------------------
// This truth matcher gives the particle with minimum deltaR below dRthrsh
// If there are more than one track below dRthrsh, the best one is returned
int TAna01Event::getGenIndexWithDeltaR(const TLorentzVector &tlv, double charge, double dRthrsh)
{
  TGenCand *pGenCand;
  int ret(-1);
  double dRmin(99999.9);
  int dRminIdx(-1);

  for (int i = 0; i < fnGenCands; i++) {
      pGenCand = getGenCand(i);
      if ((pGenCand->fStatus == 1 || pGenCand->fStatus == 8) && charge * pGenCand->fQ > 0) {
	  double dR = tlv.DeltaR(pGenCand->fP);
	  if(dR < dRthrsh && dR < dRmin)
	  {
	      dRmin = dR;
	      dRminIdx = i;
	  }
      }
  }
  ret = dRminIdx;
  return ret;
}

// ----------------------------------------------------------------------
// This truth matcher gives the particle with minimum deltaR below dRthrsh
// If there are more than one track below dRthrsh, the best one is returned
int TAna01Event::getGenIndexWithDeltaR(const TLorentzVector &tlv, const TVector3 &vtx, double charge, double dRthrsh, double dVtxRatioThrsh)
{
  if(0.0==dVtxRatioThrsh) throw std::invalid_argument("TAna01Event::getGenIndexWithDeltaR: double dVtxRatioThrsh == 0 leads to division by zero");
  TGenCand *pGenCand;
  double dRmin(99999.9);
  int ret(-1);
  int dRminIdx(-1);

  for (int i = 0; i < fnGenCands; i++) {
      pGenCand = getGenCand(i);
      if ((pGenCand->fStatus == 1 || pGenCand->fStatus == 8) && charge * pGenCand->fQ > 0) {
	  const double dR = tlv.DeltaR(pGenCand->fP);
	  double dVtxRatio;
	  try
	  {
	    dVtxRatio = vtx.Mag()/pGenCand->fV.Mag();
	  }
	  catch (std::runtime_error)
	  {
	      dVtxRatio = 0; // in case we have div by zero we make sure the next check fails
	  }
	  if(dR < dRthrsh && dR < dRmin && dVtxRatio < dVtxRatioThrsh && dVtxRatio > 1./dVtxRatioThrsh)
	  {
	      dRmin = dR;
	      dRminIdx = i;
	  }
      }
  }
  ret = dRminIdx;
  return ret;
}

// ----------------------------------------------------------------------
int TAna01Event::getGenIndexWithDeltaR(double pt, double eta, double phi, double charge) {
  TGenCand *pGenCand;
  int index(-1);
  int tmp_index[30];
  TVector3 track, gencand;
  int count(0);
  track.SetPtEtaPhi(pt, eta, phi);
  
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    if ((pGenCand->fStatus == 1 || pGenCand->fStatus == 8) && (charge*pGenCand->fQ > 0)) {
      double dR = track.DeltaR(pGenCand->fP.Vect());
      double dpt = TMath::Abs((track.Perp() - pGenCand->fP.Perp())/track.Perp());
      if (dR < 0.12 && dpt < 0.3  && count < 30) {
	tmp_index[count] = i;
	count++;
      }
    }
  }
  
  if (0 == count) return index;
  if (1 == count) return tmp_index[0];	

  double deltaR_min(0.12);
  int  i_min(-1); 	
  double dr(1.); 
  for (int i = 0; i < count; i++) {
    pGenCand = getGenCand(tmp_index[i]);
    dr = track.DeltaR(pGenCand->fP.Vect()); 
    if (dr < deltaR_min) {
      deltaR_min = dr;
      i_min = tmp_index[i];
    }
  }
  index = i_min;
  return index;

}


// ----------------------------------------------------------------------
TAnaTrack* TAna01Event::getRecTrack(Int_t n) { 
  return (TAnaTrack*)fRecTracks->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaTrack* TAna01Event::addRecTrack() {
  TClonesArray& d = *fRecTracks; 
  new(d[d.GetLast()+1]) TAnaTrack(fnRecTracks);
  ++fnRecTracks;
  return (TAnaTrack*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TSimpleTrack* TAna01Event::getSimpleTrack(Int_t n) { 
  return (TSimpleTrack*)fSimpleTracks->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TSimpleTrack* TAna01Event::addSimpleTrack() {
  TClonesArray& d = *fSimpleTracks; 
  new(d[d.GetLast()+1]) TSimpleTrack(fnSimpleTracks);
  ++fnSimpleTracks;
  TSimpleTrack *st = (TSimpleTrack*)d[d.GetLast()];
  st->clear(); // to flush the bits and indices
  return st; 
}

// ----------------------------------------------------------------------
int TAna01Event::getSimpleTrackMCID(int n) {
  TSimpleTrack *st = getSimpleTrack(n); 
  int gidx = st->getGenIndex(); 
  if (gidx > -1 && gidx < nGenCands()) {
    return getGenCand(gidx)->fID; 
  } 
  return -1; 
}


// ----------------------------------------------------------------------
int TAna01Event::getSimpleTrackMuonIdx(int n) {
  TSimpleTrack *st = getSimpleTrack(n); 
  if (st->getMuonID()) {
    TAnaMuon *pm; 
    for (int i = 0; i < fnMuons; ++i) {
      pm = getMuon(i); 
      if (n == pm->fIndex) return i;
    }
  } 
  return -1; 
}

// ----------------------------------------------------------------------
TAnaMuon* TAna01Event::getSimpleTrackMuon(int n) {
  TSimpleTrack *st = getSimpleTrack(n); 
  if (st->getMuonID()) {
    TAnaMuon *pm; 
    for (int i = 0; i < fnMuons; ++i) {
      pm = getMuon(i); 
      if (n == pm->fIndex) return pm;
    }
  } 
  return 0; 
}


// ----------------------------------------------------------------------
TAnaMuon* TAna01Event::getMuon(Int_t n) { 
  return (TAnaMuon*)fMuons->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaMuon* TAna01Event::addMuon() {
  TClonesArray& d = *fMuons; 
  new(d[d.GetLast()+1]) TAnaMuon(fnMuons);
  ++fnMuons;
  return (TAnaMuon*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
TTrgObj* TAna01Event::getTrgObj(Int_t n) { 
  return (TTrgObj*)fTrgObj->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TTrgObj* TAna01Event::addTrgObj() {
  TClonesArray& d = *fTrgObj; 
  new(d[d.GetLast()+1]) TTrgObj(fnTrgObj);
  ++fnTrgObj;
  return (TTrgObj*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
TAnaTrack* TAna01Event::getSigTrack(Int_t n) { 
  return (TAnaTrack*)fSigTracks->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaTrack* TAna01Event::addSigTrack() {
  TClonesArray& d = *fSigTracks; 
  new(d[d.GetLast()+1]) TAnaTrack(fnSigTracks);
  ++fnSigTracks;
  return (TAnaTrack*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::getCaloJet(Int_t n) { 
  return (TAnaJet*)fCaloJets->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::addCaloJet() {
  TClonesArray& d = *fCaloJets; 
  new(d[d.GetLast()+1]) TAnaJet(fnCaloJets);
  ++fnCaloJets;
  return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::getTrackJet(Int_t n) {
 return (TAnaJet*)fTrackJets->UncheckedAt(n);
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::getPFJet(Int_t n) {
  return (TAnaJet*)fPFJets->UncheckedAt(n);
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::addTrackJet() {
 TClonesArray& d = *fTrackJets;
 new(d[d.GetLast()+1]) TAnaJet(fnTrackJets);
 ++fnTrackJets;
 return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::addPFJet() {
 TClonesArray& d = *fPFJets;
 new(d[d.GetLast()+1]) TAnaJet(fnPFJets);
 ++fnPFJets;
 return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::getGenJet(Int_t n) { 
  return (TAnaJet*)fGenJets->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaJet* TAna01Event::addGenJet() {
  TClonesArray& d = *fGenJets; 
  new(d[d.GetLast()+1]) TAnaJet(fnGenJets);
  ++fnGenJets;
  return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaCand* TAna01Event::getCand(Int_t n) { 
  return (TAnaCand*)fCandidates->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaCand* TAna01Event::addCand() {
  TClonesArray& d = *fCandidates; 
  new(d[d.GetLast()+1]) TAnaCand(fnCandidates);
  ++fnCandidates;
  return (TAnaCand*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaVertex* TAna01Event::getPV(Int_t n) { 
  return (TAnaVertex*)fPV->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaVertex* TAna01Event::addPV() {
  TClonesArray& d = *fPV; 
  new(d[d.GetLast()+1]) TAnaVertex();
  ++fnPV;
  return (TAnaVertex*)d[d.GetLast()];
}



// ----------------------------------------------------------------------
void TAna01Event::dump() {
  cout << "Run: " << fRunNumber << ", Event: " << fEventNumber 
       << " nGenCand = " << nGenCands()
       << " nRecTrk = " << nRecTracks()
       << endl;

  cout << "Primary Vtx: "; 
  TAnaVertex *pPV;
  for (int i = 0; i < fnPV; i++) {
    pPV = getPV(i);
    if (i > 0) cout << "             ";
    pPV->dump();
  }
  
  for (int j = 0; j < NL1T; ++j) {
    cout << Form("%15s %2d ", fL1TNames[j].Data(), (fL1TResult[j]? 1: 0)) << endl; 
  }

  for (int i = 0; i < NLTT; ++i) {
    cout << Form("%15s %2d ", fLTTNames[i].Data(), (fLTTResult[i]? 1: 0)) << endl; 
  }

  for (int i = 0; i < NHLT; ++i) {
    cout << Form("%15s %2d ", fHLTNames[i].Data(), (fHLTResult[i]? 1: 0)); 
  }

  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    cout << "genT: ";
    pGenCand->dump();
  }
  
  TSimpleTrack *sTrk;
  for (int i = 0; i < fnSimpleTracks; ++i) {
    sTrk = getSimpleTrack(i);
    cout << "simpleT: ";
    sTrk->dump();
  }

  TAnaTrack *pTrk;
  TAnaCand  *pCand;
  for (int i = 0; i < fnRecTracks; ++i) {
    pTrk = getRecTrack(i);
    cout << "recT: ";
    pTrk->dump();
  }

  TAnaMuon  *pMuon;
  for (int i = 0; i < fnMuons; ++i) {
    pMuon = getMuon(i);
    cout << "muon: ";
    pMuon->dump();
  }
  
  TTrgObj  *pTrgObj;
  for (int i = 0; i < fnTrgObj; ++i) {
    pTrgObj = getTrgObj(i);
    cout << "trgObj: ";
    pTrgObj->dump();
  }
  
  for (int i = 0; i < fnSigTracks; ++i) {
    pTrk = getSigTrack(i);
    cout << "sigT: ";
    pTrk->dump();
  }

  for (int i = 0; i < fnCandidates; ++i) {
    pCand = getCand(i);
    cout << "cand: ";
    pCand->dump();
  }

  cout << endl;
}


// ----------------------------------------------------------------------
bool TAna01Event::isAncestor(TGenCand *pMother, TGenCand *pParticle) {
  if (pParticle == pMother) return true;
  
  TGenCand *pGen = pParticle; 
  if (pGen->fMom1 < 0) return false; 
  if (pGen->fMom1 > nGenCands()) return false; 

  while ((pGen = getGenCand(pGen->fMom1))) {
    if (pGen == pMother) return true;
    if (pGen->fMom1 < 0) break;
    if (pGen->fMom1 > nGenCands()) break;
  }
  return false; 
}
