#include "TAnaMuon.hh"
#include <iostream>
#include <bitset>

ClassImp(TAnaMuon)

using namespace std;

TAnaMuon::TAnaMuon(int index) { 
  fIndex = index;
}

TAnaMuon::~TAnaMuon() {  
}

// ----------------------------------------------------------------------
void TAnaMuon::clear() {
  fMuonChi2 = fMuonZ =  fMuonR = fCaloComp = fSegmComp = -999.;
  fTimeInOut = fTimeInOutE = fTimeOutIn = fTimeOutInE = -999.;
  fTimeNdof = -999;
  
  fNstTracks.clear();
  
  fMuonTrackPosAtM1.SetXYZ(-1e30,-1e30,-1e30);
  fMuonTrackPlabAtM1.SetXYZ(-1e30,-1e30,-1e30);

  fInnerPlab.SetXYZ(-1e30,-1e30,-1e30);
  fOuterPlab.SetXYZ(-1e30,-1e30,-1e30);
  fGlobalPlab.SetXYZ(-1e30,-1e30,-1e30);

  fPositionAtM1.SetXYZ(-1e30,-1e30,-1e30); 
  fPositionAtM2.SetXYZ(-1e30,-1e30,-1e30);
  
  fItrkValidFraction = fGtrkNormChi2 = fSegmentComp = fChi2LocalPosition = fChi2LocalMomentum = -999.;
  fNumberOfLostTrkHits = fNumberOfValidTrkHits = -99;
  fNvalidMuonHits = -99; 

  fGtrkProb = -99.;

  fInt1 = fInt2 = fInt3 -99; 
  fDouble1 = fDouble2 = fDouble3 = -99.;

  for (unsigned int i = 0; i < NXPTRACKS; ++i) {
    fXpTracks[i].dist = 99.;
  }
  
  fVtxProb = -99.0;
  fVtxTracks.clear();
}


// ----------------------------------------------------------------------
void TAnaMuon::dump() {
  //   cout << " q = " << fQ
  //        << " p = " << fPlab.Mag() 
  //        << " f = " << fPlab.Phi() 
  //        << " t = " << fPlab.Theta();

  cout << "Muon: idx = " << fIndex;
  cout << Form(" midx=%3d", fMuIndex) << " ID=" << static_cast<std::bitset<32> >(fMuID) << endl;

  if (fPlab.Perp() > 0) {
    cout << Form(" P:      q=%+2d pT=%6.2f eta=%+4.3f f=%+4.3f ", fQ, fPlab.Perp(), fPlab.Eta(), fPlab.Phi()) << endl;
  }
  if (fInnerPlab.Perp() > 0) {
    cout << Form(" inner:  q=%+2d pT=%6.2f eta=%+4.3f f=%+4.3f ", fQ, fInnerPlab.Perp(), fInnerPlab.Eta(), fInnerPlab.Phi()) << endl;
  } 
  if (fOuterPlab.Perp() > 0) {
    cout << Form(" outer:  q=%+2d pT=%6.2f eta=%+4.3f f=%+4.3f ", fQ, fOuterPlab.Perp(), fOuterPlab.Eta(), fOuterPlab.Phi()) << endl;
  } 
  if (fGlobalPlab.Perp() > 0) {
    cout << Form(" global: q=%+2d pT=%6.2f eta=%+4.3f f=%+4.3f ", fQ, fGlobalPlab.Perp(), fGlobalPlab.Eta(), fGlobalPlab.Phi()) << endl;
  } 
  cout << Form(" X=%+4.1f z=%4.1f r=%+4.1f cl=%4.1f c2=%4.1f", fMuonChi2, fMuonZ, fMuonR, fCaloComp, fSegmComp);
  cout << Form(" DT=%3d CSC=%3d RPC=%3d t=%4.1f+/-%4.1f (%3d)", fNhitsDT, fNhitsCSC, fNhitsRPC, fTimeInOut, fTimeInOutE, fTimeNdof);
  cout << endl;

  if (fPositionAtM1.Perp() > 0) {
    cout << " Position at M1: x1=" << fPositionAtM1.X() << " y1=" << fPositionAtM1.Y() << " z1=" << fPositionAtM1.Z() << endl;
  }
  if (fPositionAtM2.Perp() > 0) {
    cout << " Position at M2: x2=" << fPositionAtM2.X() << " y2=" << fPositionAtM2.Y() << " z2=" << fPositionAtM2.Z() << endl;
  }
  
  if (fMuonTrackPosAtM1.X() > -1e30 && fMuonTrackPosAtM1.Y() > -1e30 && fMuonTrackPosAtM1.Z() > -1e30) {
	  cout << " Position at M1 of muon track: x=" << fMuonTrackPosAtM1.X() << " y=" << fMuonTrackPosAtM1.Y() << " z=" << fMuonTrackPosAtM1.Z() << endl;
  }
  
  if (fMuonTrackPlabAtM1.X() > -1e30 && fMuonTrackPlabAtM1.Y() > -1e30 && fMuonTrackPlabAtM1.Z() > -1e30) {
	  cout << " Momentum at M1 of muon track: x=" << fMuonTrackPlabAtM1.X() << " y=" << fMuonTrackPlabAtM1.Y() << " z=" << fMuonTrackPlabAtM1.Z() << endl;
  }
  
  cout << " closeby tracks (" << fNstTracks.size() << ") are" << endl;
  for (std::map<int,float>::const_iterator it = fNstTracks.begin(); it != fNstTracks.end(); ++it)
	  cout << "   (ix=" << it->first << ", doca=" << it->second << ")" << endl;

  cout << " vertex probability vtxProb=" << fVtxProb << " with tracks ( ";
  for (std::set<unsigned>::const_iterator it = fVtxTracks.begin(); it != fVtxTracks.end(); ++it) cout << *it << " ";
  cout << ")" << endl;

}
  
