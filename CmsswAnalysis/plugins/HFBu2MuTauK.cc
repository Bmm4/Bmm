#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBu2MuTauK.h"

#include <iostream>
#include <utility>

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"

#include "Bmm/CmsswAnalysis/interface/HFSequentialVertexFit.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"
#include "Bmm/CmsswAnalysis/interface/HFThreeParticleCombinatorics.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFDecayTree.hh"
#include "Bmm/CmsswAnalysis/interface/HFTrackListBuilder.hh"
#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// -- Yikes!
extern TAna01Event *gHFEvent;

using namespace std;
using namespace reco;
using namespace edm;

// 32000 = mu K tau(hhh)
// 33000 = mu+ pi+ D-(K+ pi- pi-)

// ----------------------------------------------------------------------
HFBu2MuTauK::HFBu2MuTauK(const ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fMcType(iConfig.getUntrackedParameter<int>("mcType", 3002000)),
  fMuKaVtxProb(iConfig.getUntrackedParameter<double>("mukaVtxProb", 0.01)),
  fTauVtxProb(iConfig.getUntrackedParameter<double>("tauVtxProb", 0.01)) {
  dumpConfiguration();
}

// ----------------------------------------------------------------------
void HFBu2MuTauK::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBu2MuTauK configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  mukaVtexProb                " << fMuKaVtxProb << endl;
  cout << "---  tauVtxProb                  " << fTauVtxProb << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFBu2MuTauK::analyze(const Event& iEvent, const EventSetup& iSetup) {
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  } catch(HFSetupException e) {
    cout << "==>HFBu2MuTauK> " << e.fMsg << endl;
    return;
  }

  // -----------------------------------
  // -- build truth-identified candidate
  // -----------------------------------
  if (fMcType > 0) {
    genMatch();
    vector<int> truthTracks = truthIndices();
    bool foundTruthCand(false);
    if (truthTracks.size() > 0) {
      foundTruthCand = true;
      for (unsigned int i = 0; i < truthTracks.size(); ++i) {
	if (truthTracks[i] < 0) {
	  foundTruthCand = false;
	  break;
	}
      }
    }
    if (foundTruthCand) {
      bool filled = buildCandidate(truthTracks, fMcType, false);
      if (!filled) {
	if (fVerbose > 0) {
	  cout << "no MC truth cand built?!" << endl;
	}
      }
    } else {
      if (fVerbose > 0) {
	cout << "no MC truth cand found?! ";
	for (unsigned int i = 0; i < truthTracks.size(); ++i) {
	  cout << truthTracks[i] << " ";
	}
	cout << endl;
      }
    }
  }

  // ------------------------
  // -- build reco candidates
  // ------------------------
  fListBuilder->setMinPt(fMuonPt);
  vector<int> muonList = fListBuilder->getMuonList();

  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setMaxDocaToTracks(0.1);
  if (fVerbose > 0) cout << "==>HFBu2MuTauK> muonList  size: " << muonList.size() << endl;
  fListBuilder->setCloseTracks(&muonList);
  vector<int> trkList = fListBuilder->getTrackList();
  if (fVerbose > 0) cout << "==>HFBu2MuTauK> trkList  size: " << trkList.size() << endl;
  if (trkList.size() < 1) return;

  TLorentzVector tlv;
  vector<pair<int, TLorentzVector> > piList, kaList;
  for (unsigned int i = 0; i < trkList.size(); ++i) {
    TrackBaseRef trackView(fTracksHandle, trkList[i]);
    Track trk(*trackView);
    tlv.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), MPION);
    piList.push_back(make_pair(trkList[i], tlv));
    tlv.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), MKAON);
    kaList.push_back(make_pair(trkList[i], tlv));
  }

  HFTwoParticleCombinatoricsNew a(fTracksHandle, fVerbose);
  HFTwoParticleCombinatoricsSet promptList;
  if (32000 == fType) {
    promptList = a.combine(muonList, MMUON, trkList, MKAON, 0.6, 3.6, 1);
  } else if (33000 == fType) {
    promptList = a.combine(muonList, MMUON, trkList, MPION, 0.2, 5.2, 1);
  }
  if (fVerbose > 0) cout << "==>HFBu2MuTauK> promptList  size: " << promptList.size() << endl;
  if (promptList.size() < 1) return;

  // -- build 3-prong
  vector<triplet> tauList;
  HFThreeParticleCombinatorics tau(fVerbose);
  if (32000 == fType) {
    tau.combine(tauList, piList, 0.3, 2.0, 0.6, 1.1); // FIXME: Are the loresmass and hiresmass cuts good?!
  } else if (33000 == fType) {
    tau.combine(tauList, piList, kaList, 0.3, 2.0); // FIXME: Are the mass cuts good?
  }
  vector<triplet> filledTriplets;
  if (fVerbose > 0) cout << "==>HFBu2MuTauK> tripletList  size: " << tauList.size() << endl;
  if (tauList.size() < 1) return;

  // -- Loop over (mu,K) pairs and all associated triplets
  for (HFTwoParticleCombinatoricsNew::iterator promptIt = promptList.begin(); promptIt != promptList.end(); ++promptIt) {
    int iMuon = promptIt->first;
    int iKaon = promptIt->second;
    for (unsigned int it = 0; it < tauList.size(); ++it) {
      // -- check against overlap with muons and/or other 3-hadron combinations already written to the tree
      if (tauList[it].pi1() == iMuon || tauList[it].pi1() == iKaon) continue;
      if (tauList[it].pi2() == iMuon || tauList[it].pi2() == iKaon) continue;
      if (tauList[it].pi3() == iMuon || tauList[it].pi3() == iKaon) continue;
      bool skip(false);
      for (unsigned ik = 0; ik < filledTriplets.size(); ++ik) {
	if (tauList[it].isPermutation(filledTriplets[ik])) {
	  if (fVerbose > 0) {
	    cout << "not filling " << tauList[it].pi1() << "/" << tauList[it].pi2() << "/" << tauList[it].pi3()
		 << " because have already " << filledTriplets[ik].pi1() << "/" << filledTriplets[ik].pi2() << "/" << filledTriplets[ik].pi3()
		 << endl;
	  }
	  skip = true;
	  break;
	}
      }
      if (skip) continue;

      vector<int> allTrks;
      allTrks.push_back(iMuon);
      allTrks.push_back(iKaon);
      allTrks.push_back(tauList[it].pi1());
      allTrks.push_back(tauList[it].pi2());
      allTrks.push_back(tauList[it].pi3());

      bool filled = buildCandidate(allTrks, fType, true);
      if (filled) {
	filledTriplets.push_back(tauList[it]);
      }
    }
  }
}


// ----------------------------------------------------------------------
bool HFBu2MuTauK::buildCandidate(std::vector<int> trkIdx, int type, bool cut) {
  if (fVerbose > 0) cout << "==HFBu2MuTauK::buildCandidate> " << type << endl;
  int iMuon = trkIdx[0];
  int iKaon = trkIdx[1];
  int iPion1 = trkIdx[2];
  int iPion3 = trkIdx[3];
  int iPion2 = trkIdx[4];

  vector<TransientTrack> vtt;
  vector<int> vid;

  reco::TrackBaseRef muTrackView(fTracksHandle, iMuon);
  reco::Track tMuon(*muTrackView);
  if (cut && tMuon.pt() < fMuonPt) return false;

  TrackBaseRef kaTrackView(fTracksHandle, iKaon);
  Track tKaon(*kaTrackView);
  if (cut && tKaon.pt() < fTrackPt) return false;

  if (cut && (fMaxD0 < 90.) && (tMuon.d0() > fMaxD0)) return false;
  if (cut && (fMaxDz < 90.) && (tMuon.dz() > fMaxDz)) return false;
  if (cut && (fMaxD0 < 90.) && (tKaon.d0() > fMaxD0)) return false;
  if (cut && (fMaxDz < 90.) && (tKaon.dz() > fMaxDz)) return false;

  TransientTrack mTT = fTTB->build(*muTrackView);
  TransientTrack kTT = fTTB->build(*kaTrackView);
  vtt.clear(); vid.clear();
  vtt.push_back(mTT); vid.push_back(13);
  vtt.push_back(kTT);
  if (32000 == fType) {
    vid.push_back(321);
  } else if (33000 == fType) {
    vid.push_back(211);
  }
  RefCountedKinematicTree mukaTree = fitTree(vtt, vid);
  if (mukaTree->isEmpty()) return false;
  RefCountedKinematicVertex mukaVtx = mukaTree->currentDecayVertex();
  double mkmass = mukaTree->currentParticle()->currentState().mass();
  double mkprob = 0.0;

  TVector3 mkplab = TVector3(mukaTree->currentParticle()->currentState().globalMomentum().x(),
			     mukaTree->currentParticle()->currentState().globalMomentum().y(),
			     mukaTree->currentParticle()->currentState().globalMomentum().z());
  if ((mukaVtx->chiSquared() >= 0.0) && (mukaVtx->degreesOfFreedom() > 0)) {
    mkprob = TMath::Prob(mukaVtx->chiSquared(), mukaVtx->degreesOfFreedom());
  }

  if (fVerbose > 0)
    cout << "muka mass(" << iMuon << "/" << iKaon << "):  mvtx = " << mkmass
	 << " and pT(mu) = " << tMuon.pt() << " pT(K) = " << tKaon.pt()
	 << " vtx prob = " << mkprob
	 << endl;

  if (cut && mkprob < fMuKaVtxProb) {
    if (fVerbose > 0) cout << "muKa vtx probability too low, skipping" << endl;
    return false;
  }
  pair<int, double> mukaPV = findBestPV(vtt, mukaTree, fVertexCollection, fMagneticField);

  Vertex mypVertex = fVertexCollection[mukaPV.first];
  Vertex mysVertex = mkVertex(mukaTree);
  pair<double, double> psSep = vtxSeparation(mypVertex, mysVertex);
  if (fVerbose > 0) cout << "vtxsep(PV, muka) = " << psSep.first << " +/- " << psSep.second << endl;

  TrackBaseRef pi1TrackView(fTracksHandle, iPion1);
  reco::Track pi1T(*pi1TrackView);
  TransientTrack pi1TT = fTTB->build(*pi1TrackView);

  TrackBaseRef pi2TrackView(fTracksHandle, iPion2);
  reco::Track pi2T(*pi2TrackView);
  TransientTrack pi2TT = fTTB->build(*pi2TrackView);

  TrackBaseRef pi3TrackView(fTracksHandle, iPion3);
  reco::Track pi3T(*pi3TrackView);
  TransientTrack pi3TT = fTTB->build(*pi3TrackView);

  if (cut && pi1T.pt() < fTrackPt)  return false;
  if (cut && pi2T.pt() < fTrackPt)  return false;
  //NO!      if (pi3T.pt() < fTrackPt)  return;

  int qtot = pi1T.charge() + pi2T.charge() + pi3T.charge();
  if (TMath::Abs(qtot) != 1) {
    if (fVerbose > 0) cout << "qtot(3 hadrons) = " << qtot << ", skipping" << endl;
    return false;
  }

  vtt.clear(); vid.clear();
  if (32000 == fType) {
    // -- For tau+ -> pi- pi+ pi+ decays, all hadrons are pions. Ignore those decays with kaons included!?
    vtt.push_back(pi1TT); vid.push_back(211);
    vtt.push_back(pi2TT); vid.push_back(211);
    vtt.push_back(pi3TT); vid.push_back(211);
  } else if (33000 == fType) {
    // -- For D+ -> K- pi+ pi+ decays, the same-sign particles are the pions. DCS decays are ignored!
    if (pi1T.charge() == pi2T.charge()) {
      vtt.push_back(pi1TT); vid.push_back(211);
      vtt.push_back(pi2TT); vid.push_back(211);
      vtt.push_back(pi3TT); vid.push_back(321);
    } else if (pi1T.charge() == pi3T.charge()) {
      vtt.push_back(pi1TT); vid.push_back(211);
      vtt.push_back(pi2TT); vid.push_back(321);
      vtt.push_back(pi3TT); vid.push_back(321);
    } else if (pi2T.charge() == pi3T.charge()) {
      vtt.push_back(pi1TT); vid.push_back(321);
      vtt.push_back(pi2TT); vid.push_back(211);
      vtt.push_back(pi3TT); vid.push_back(321);
    }
  }

  RefCountedKinematicTree tauTree = fitTree(vtt, vid);
  if (tauTree->isEmpty()) return false;
  RefCountedKinematicVertex tauVtx = tauTree->currentDecayVertex();
  TVector3 tauplab = TVector3(tauTree->currentParticle()->currentState().globalMomentum().x(),
			      tauTree->currentParticle()->currentState().globalMomentum().y(),
			      tauTree->currentParticle()->currentState().globalMomentum().z());
  double tauprob = 0.0;
  if ((tauVtx->chiSquared() >= 0.0) && (tauVtx->degreesOfFreedom() > 0)) {
    tauprob = TMath::Prob(tauVtx->chiSquared(), tauVtx->degreesOfFreedom());
  }

  if (cut && tauprob < fTauVtxProb) {
    if (fVerbose > 0) cout << "tau vtx probability too low (" << tauprob << ") for pion idxs = "
			   << iPion1 << "/" << iPion2 << "/" << iPion3
			   << " skipping" << endl;
    return false;
  }

  Vertex mytVertex = mkVertex(tauTree);
  pair<double, double> stSep = vtxSeparation(mytVertex, mysVertex);

  double taumass = tauTree->currentParticle()->currentState().mass();
  if (fVerbose > 0) {
    TLorentzVector tlv,  muka, mu, ka, bu, pi1, pi2, pi3;
    mu.SetPtEtaPhiM(tMuon.pt(), tMuon.eta(), tMuon.phi(), MMUON);
    if (32000 == fType) {
      ka.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MKAON);
    } else if (33000 == fType) {
      ka.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MPION);
    }
    muka = mu + ka;
    pi1.SetPtEtaPhiM(pi1T.pt(), pi1T.eta(), pi1T.phi(), MPION);
    pi2.SetPtEtaPhiM(pi2T.pt(), pi2T.eta(), pi2T.phi(), MPION);
    pi3.SetPtEtaPhiM(pi3T.pt(), pi3T.eta(), pi3T.phi(), MPION);
    tlv = pi1 + pi2 + pi3;
    cout << "indices(" << iMuon << "/" << iKaon << "/"
	 << iPion1 << "/" << iPion2 << "/" << iPion3
	 << ": m4v(h3) = " << tlv.M() << " mvtx(h3) = " << taumass
	 << " vtx prob = " << tauprob
	 << " separation = " << stSep.first << " +/- " << stSep.second << endl
	 << " FILLING CANDIDATE " << type
	 << endl;
  }

  vector<int> mkidx;
  mkidx.push_back(iMuon);
  mkidx.push_back(iKaon);

  TVector3 plab = mkplab + tauplab;
  TAnaCand *pCand = fillCand(type, mukaTree, mkidx, mypVertex);
  pCand->fPvIdx = mukaPV.first;

  vector<int> h3idx;
  h3idx.push_back(iPion1);
  h3idx.push_back(iPion2);
  h3idx.push_back(iPion3);
  TAnaCand *tauCand = fillCand(type+15, tauTree, h3idx, mysVertex);
  tauCand->fPvIdx = -1;

  // -- now fix the pointers to daughter and mother
  pCand->fDau1  = tauCand->fIndex;
  pCand->fDau2  = tauCand->fIndex;
  tauCand->fMom = pCand->fIndex;

  // -- calculate tau and B effective decay tiume
  double massOverC = tauCand->fMass/TMath::Ccgs();

  const TVector3 p2(mytVertex.position().x(), mytVertex.position().y(), mytVertex.position().z());
  const TVector3 p1(mysVertex.position().x(), mysVertex.position().y(), mysVertex.position().z());
  TVector3 pDiff = p2-p1;

  double taucosa   = tauCand->fPlab.Dot(pDiff) / (tauCand->fPlab.Mag() * pDiff.Mag());
  tauCand->fTau3d = stSep.first / tauCand->fPlab.Mag() * taucosa * massOverC;


  massOverC = MBPLUS/TMath::Ccgs();
  const TVector3 p0(mypVertex.position().x(), mypVertex.position().y(), mypVertex.position().z());
  pDiff = p1 - p0;
  double bcosa  = plab.Dot(pDiff) / (plab.Mag() * pDiff.Mag());
  pCand->fTau3d = psSep.first / pCand->fPlab.Mag() * bcosa * massOverC;

  if (fVerbose > 0) {
    cout << "HFListCand dump: " << endl;
    pCand->dump();
    for (int is = pCand->fDau1; is <= pCand->fDau2; ++is) {
      TAnaCand *tc = gHFEvent->getCand(is);
      tc->dump();
      for (int is = tc->fSig1; is <= tc->fSig2; ++is) {
	gHFEvent->getSigTrack(is)->dump();
      }
    }
    for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
      gHFEvent->getSigTrack(is)->dump();
    }
  }
  return true;
}


// ----------------------------------------------------------------------
TAnaCand* HFBu2MuTauK::fillCand(int type, RefCountedKinematicTree &tree, std::vector<int> trkIdx, reco::Vertex &otherVtx) {
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fType  = type;
  pCand->fSig1  = gHFEvent->nSigTracks();
  pCand->fSig2  = pCand->fSig1 + trkIdx.size() - 1;

  pCand->fPlab  = TVector3(tree->currentParticle()->currentState().globalMomentum().x(),
			   tree->currentParticle()->currentState().globalMomentum().y(),
			   tree->currentParticle()->currentState().globalMomentum().z());
  pCand->fMass  = tree->currentParticle()->currentState().mass();

  const reco::BeamSpot *pBeamSpot = &fBeamSpot;
  TAnaTrack *pTrack(0);
  int q(0);
  for (unsigned int i = 0; i < trkIdx.size(); ++i) {
    pTrack = fillSigTrack(trkIdx[i], fTracksHandle, &fVertexCollection, fMuonCollection, pBeamSpot);
    q += pTrack->fQ;
  }
  pCand->fQ = q;


  // -- decay vertex
  RefCountedKinematicVertex vtx = tree->currentDecayVertex();
  double prob = 0.0;
  if ((vtx->chiSquared() >= 0.0) && (vtx->degreesOfFreedom() > 0)) {
    prob = TMath::Prob(vtx->chiSquared(), vtx->degreesOfFreedom());
  }
  pCand->fVtx.setInfo(vtx->chiSquared(), vtx->degreesOfFreedom(), prob, 0, -1);

  Vertex myVertex = mkVertex(tree);
  const TVector3 p1(myVertex.position().x(), myVertex.position().y(), myVertex.position().z());
  pCand->fVtx.fPoint   = p1;

  pair<double, double> vtxSep = vtxSeparation(myVertex, otherVtx);
  pCand->fVtx.fD3d     = vtxSep.first;
  pCand->fVtx.fD3dE    = vtxSep.second;

  return pCand;
}

// ----------------------------------------------------------------------
vector<int>  HFBu2MuTauK::truthIndices() {
  vector<int> v;
  if (0 == fpGenKa) return v;
  if (0 == fpGenMu) return v;
  if (0 == fpGenHad1) return v;
  if (0 == fpGenHad2) return v;
  if (0 == fpGenHad3) return v;
  int iMuon(-1), iKaon(-1), iPion1(-1), iPion2(-1), iPion3(-1);
  TSimpleTrack *ps(0);
  for (int i = 0; i < gHFEvent->nSimpleTracks(); ++i) {
    ps = gHFEvent->getSimpleTrack(i);
    if (ps->getGenIndex() == fpGenMu->fNumber) {
      iMuon = i;
    }
    if (ps->getGenIndex() == fpGenKa->fNumber) {
      iKaon = i;
    }
    if (ps->getGenIndex() == fpGenHad1->fNumber) {
      iPion1 = i;
    }
    if (ps->getGenIndex() == fpGenHad2->fNumber) {
      iPion2 = i;
    }
    if (ps->getGenIndex() == fpGenHad3->fNumber) {
      iPion3 = i;
    }
  }

  v.push_back(iMuon);
  v.push_back(iKaon);
  v.push_back(iPion1);
  v.push_back(iPion2);
  v.push_back(iPion3);
  return v;
}

// ----------------------------------------------------------------------
void HFBu2MuTauK::genMatch() {
  // -- now try to find signal decay
  bool goodMatch(false);
  TGenCand *pCand, *pD, *pT;
  for (int iC = 0; iC < gHFEvent->nGenCands(); ++iC) {
    pCand = gHFEvent->getGenCand(iC);
    if (521 == TMath::Abs(pCand->fID)) {
      fpGenB = fpGenMu = fpGenKa = fpGenTau = fpGenHad1 = fpGenHad2 = fpGenHad3 = 0;
      fpGenB = pCand;
      for (int iD = pCand->fDau1; iD <= pCand->fDau2; ++iD) {
	pD = gHFEvent->getGenCand(iD);
	if (321 == TMath::Abs(pD->fID)) fpGenKa = pD;
	if (13 == TMath::Abs(pD->fID)) fpGenMu = pD;
	if (15 == TMath::Abs(pD->fID)) {
	  fpGenTau = pD;
	  for (int iT = pD->fDau1; iT <= pD->fDau2; ++iT) {
	    pT = gHFEvent->getGenCand(iT);
	    if (0 != pT->fQ) {
	      if (0 == fpGenHad1) {
		fpGenHad1 = pT;
	      } else if (0 == fpGenHad2) {
		fpGenHad2 = pT;
	      } else if (0 == fpGenHad3) {
		fpGenHad3 = pT;
	      } else {
		cout << "XXXXXXXXXXXXXXXX what should I do here?????" << endl;
	      }
	    }
	  }
	}

      }
      if (fpGenB && fpGenMu && fpGenKa && fpGenTau && fpGenHad1 && fpGenHad2 && fpGenHad3) {
	goodMatch = true;
	break;
      }
    }
  }

  if (goodMatch) {
    // -- sort according to pT of hadronic tracks
    double pt1 = fpGenHad1->fP.Perp();
    double pt2 = fpGenHad2->fP.Perp();
    double pt3 = fpGenHad3->fP.Perp();
    TGenCand *tgc1(0), *tgc2(0), *tgc3(0);
    if (pt1 > pt2 && pt1 > pt3) {
      tgc1 = fpGenHad1;
      if (pt2 > pt3) {
	tgc2 = fpGenHad2;
	tgc3 = fpGenHad3;
      } else {
	tgc2 = fpGenHad3;
	tgc3 = fpGenHad2;
      }
    } else if (pt2 > pt1 && pt2 > pt3) {
      tgc1 = fpGenHad2;
      if (pt1 > pt3) {
	tgc2 = fpGenHad1;
	tgc3 = fpGenHad3;
      } else {
	tgc2 = fpGenHad3;
	tgc3 = fpGenHad1;
      }
    } else if (pt3 > pt1 && pt3 > pt2) {
      tgc1 = fpGenHad3;
      if (pt1 > pt2) {
	tgc2 = fpGenHad1;
	tgc3 = fpGenHad2;
      } else {
	tgc2 = fpGenHad2;
	tgc3 = fpGenHad1;
      }
    } else {
      cout << "THIS SHOULD NOT HAPPEN!!!" << endl;
    }
    fpGenHad1 = tgc1;
    fpGenHad2 = tgc2;
    fpGenHad3 = tgc3;
  }
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBu2MuTauK);
