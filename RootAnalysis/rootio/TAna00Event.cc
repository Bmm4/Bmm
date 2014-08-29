#include "TAna00Event.hh"
#include <iostream>

ClassImp(TAna00Event)

using namespace std;

// ----------------------------------------------------------------------
TAna00Event::TAna00Event(Int_t Option) {
  fRecTracks       = new TClonesArray("TAnaTrack", 1000);
  fnRecTracks      = 0;

  fSigTracks       = new TClonesArray("TAnaTrack", 1000);
  fnSigTracks      = 0;

  fCaloJets        = new TClonesArray("TAnaJet", 1000);
  fnCaloJets       = 0;

  fGenJets         = new TClonesArray("TAnaJet", 1000);
  fnGenJets        = 0;

  fCandidates      = new TClonesArray("TAnaCand", 1000);
  fnCandidates     = 0;

  fGenCands        = new TClonesArray("TGenCand", 1000);
  fnGenCands       = 0;
}

// ----------------------------------------------------------------------
void TAna00Event::Clear(Option_t *option) {
  int i;

  fPrimaryVertex.clear();
  fPrimaryVertex2.clear();

  TAnaTrack *pTrack;
  //  cout << "Starting to clear tracks ";
  for (i = 0; i < fnRecTracks; i++) {
    pTrack = getRecTrack(i);
    pTrack->clear();
  }
  fRecTracks->Clear(option);
  fnRecTracks = 0;
  //  cout << " ... done " << endl;

  TAnaTrack *pSigTrack;
  //  cout << "Starting to clear signal tracks ";
  for (i = 0; i < fnSigTracks; i++) {
    pSigTrack = getSigTrack(i);
    pSigTrack->clear();
  }
  fSigTracks->Clear(option);
  fnSigTracks = 0;
  //  cout << " ... done " << endl;

  TAnaJet *pCaloJet;
  //  cout << "Starting to clear jets ";
  for (i = 0; i < fnCaloJets; i++) {
    pCaloJet = getCaloJet(i);
    pCaloJet->clear();
  }
  fCaloJets->Clear(option);
  fnCaloJets = 0;
  //  cout << " ... done " << endl;

  TAnaJet *pGenJet;
  //  cout << "Starting to clear jets ";
  for (i = 0; i < fnGenJets; i++) {
    pGenJet = getGenJet(i);
    pGenJet->clear();
  }
  fGenJets->Clear(option);
  fnGenJets = 0;
  //  cout << " ... done " << endl;

  TAnaCand *pCand;
  //  cout << "Starting to clear candidates ";
  for (i = 0; i < fnCandidates; i++) {
    pCand = getCand(i);
    pCand->clear();
  }
  fCandidates->Clear(option);
  fnCandidates = 0;
  //  cout << " ... done " << endl;

  TGenCand *pGenCand;
  //  cout << "Starting to clear gencands ";
  for (i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->clear();
  }
  fGenCands->Clear(option);
  fnGenCands = 0;
  //  cout << " ... done " << endl;

}


// ----------------------------------------------------------------------
TGenCand* TAna00Event::getGenCand(Int_t n) { 
  return (TGenCand*)fGenCands->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TGenCand* TAna00Event::addGenCand() {
  TClonesArray& d = *fGenCands; 
  new(d[d.GetLast()+1]) TGenCand();
  ++fnGenCands;
  return (TGenCand*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
void TAna00Event::dumpGenBlock() {
  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->dump();
  }
}

// ----------------------------------------------------------------------
int TAna00Event::getGenIndex(double px, double py, double pz, int id, double precision) {
  TGenCand *pGenCand;
  int index(-1);

  //   cout << "Searching among " << fnGenCands << " cands for match to preco = (" 
  //         << px << ", " << py << ", " << pz << ") with ID = " << id << " and precision = " << precision << endl;

  int aid(0), bestMatch(-1); 
  double bestDeltaR(99.); 
  TVector3 vr(0., 0., 0.);
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    
    // -- Look only at stable particles!
    aid = TMath::Abs(pGenCand->fID); 
    if ((aid != 211)
	&& (aid != 321)
	&& (aid != 2212)
	&& (aid != 11)
	&& (aid != 13)
	) {
      continue;
    }

    if (id != -1) {
      if (pGenCand->fID != id) {
	continue;
      }
    }
    if (TMath::Abs((pGenCand->fP.X() - px)/px) > precision) {
      continue;
    }

    if (TMath::Abs((pGenCand->fP.Y() - py)/py) > precision) {
      continue;
    }

    if (TMath::Abs((pGenCand->fP.Z() - pz)/pz) > precision) {
      continue;
    }
    
    vr = TVector3(px, py, pz); 
    double deltaR = vr.DeltaR(pGenCand->fP.Vect()); 
    //    cout << deltaR << endl;
    if (deltaR < bestDeltaR) {
      if (bestMatch > -1) {
	cout << "--> TAna00Event::getGenIndex() Replacing matched object, delta(r) = " << bestDeltaR << ",  with new delta(r) = " << deltaR << endl;
      }
      bestDeltaR = deltaR; 
      bestMatch = i;
      index = i; 
    }

    //     cout << "I think there is a match at index= " << index << " i = " << i 
    //      	 << " for reco p = " << px << " " << py << " " << pz << " " 
    // 	 << ", gen p = " << pGenCand->fP.X() << " " << pGenCand->fP.Y() << " " << pGenCand->fP.Z() << " " << endl;
    //    pGenCand->dump();
    break;
  }
  return index;
}


// ----------------------------------------------------------------------
TAnaTrack* TAna00Event::getRecTrack(Int_t n) { 
  return (TAnaTrack*)fRecTracks->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaTrack* TAna00Event::addRecTrack() {
  TClonesArray& d = *fRecTracks; 
  new(d[d.GetLast()+1]) TAnaTrack(fnRecTracks);
  ++fnRecTracks;
  return (TAnaTrack*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
TAnaTrack* TAna00Event::getSigTrack(Int_t n) { 
  return (TAnaTrack*)fSigTracks->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaTrack* TAna00Event::addSigTrack() {
  TClonesArray& d = *fSigTracks; 
  new(d[d.GetLast()+1]) TAnaTrack(fnSigTracks);
  ++fnSigTracks;
  return (TAnaTrack*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna00Event::getCaloJet(Int_t n) { 
  return (TAnaJet*)fCaloJets->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaJet* TAna00Event::addCaloJet() {
  TClonesArray& d = *fCaloJets; 
  new(d[d.GetLast()+1]) TAnaJet(fnCaloJets);
  ++fnCaloJets;
  return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaJet* TAna00Event::getGenJet(Int_t n) { 
  return (TAnaJet*)fGenJets->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaJet* TAna00Event::addGenJet() {
  TClonesArray& d = *fGenJets; 
  new(d[d.GetLast()+1]) TAnaJet(fnGenJets);
  ++fnGenJets;
  return (TAnaJet*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TAnaCand* TAna00Event::getCand(Int_t n) { 
  return (TAnaCand*)fCandidates->UncheckedAt(n); 
}

// ----------------------------------------------------------------------
TAnaCand* TAna00Event::addCand() {
  TClonesArray& d = *fCandidates; 
  new(d[d.GetLast()+1]) TAnaCand();
  ++fnCandidates;
  return (TAnaCand*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
void TAna00Event::dump() {
  cout << "Run: " << fRunNumber << ", Event: " << fEventNumber 
       << " nGenCand = " << nGenCands()
       << " nRecTrk = " << nRecTracks()
       << endl;

  cout << "Primary Vtx: "; 
  fPrimaryVertex.dump();
  
  Int_t i;
  
  TGenCand *pGenCand;
  for (i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    cout << "genT: ";
    pGenCand->dump();
  }
  
  TAnaTrack *pTrk;
  TAnaCand  *pCand;
  for (i = 0; i < nRecTracks(); ++i) {
    pTrk = getRecTrack(i);
    cout << "recT: ";
    pTrk->dump();
  }
  
  for (i = 0; i < nSigTracks(); ++i) {
    pTrk = getSigTrack(i);
    cout << "sigT: ";
    pTrk->dump();
  }

  for (i = 0; i < nCands(); ++i) {
    pCand = getCand(i);
    cout << "cand: ";
    pCand->dump();
  }

  cout << endl;
}


