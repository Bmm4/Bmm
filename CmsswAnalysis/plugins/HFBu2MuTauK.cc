#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBu2MuTauK.h"

#include <iostream>
#include <utility>

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"

#include "Bmm/CmsswAnalysis/interface/HFSequentialVertexFit.hh"
#include "Bmm/CmsswAnalysis/interface/HFCombinatorics.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFDecayTree.hh"
#include "Bmm/CmsswAnalysis/interface/HFTrackListBuilder.hh"
#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"


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

// recoil + 32 = mu K tau(hhh)
// recoil + 33 = mu+ pi+ D-(K+ pi- pi-)
// in both cases, recoil + X + 100 is the three prong candidate

// ----------------------------------------------------------------------
HFBu2MuTauK::HFBu2MuTauK(const ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fMcType(iConfig.getUntrackedParameter<int>("mcType", 4000032)),
  fMuKaVtxProb(iConfig.getUntrackedParameter<double>("mukaVtxProb", 0.01)),
  fMuKaMaxDoca(iConfig.getUntrackedParameter<double>("mukaMaxDoca", 0.02)),
  fTauVtxProb(iConfig.getUntrackedParameter<double>("tauVtxProb", 0.01)),
  fTauMaxDoca(iConfig.getUntrackedParameter<double>("tauMaxDoca", 0.02)),
  fTauMaxIP(iConfig.getUntrackedParameter<double>("tauMaxIP", 0.02)) {

  if (fType%100 == 32) {
    fMode = B2MUKTAU;
  } else if (fType%100 == 41) {
    fMode = B2PIPID;
  } else if (fType%100 == 51) {
    fMode = B2MUPID;
  }
  dumpConfiguration();
}

// ----------------------------------------------------------------------
void HFBu2MuTauK::dumpConfiguration() {
  string cMode("");
  if (B2MUKTAU == fMode) cMode = "B2MUKTAU";
  if (B2PIPID == fMode) cMode = "B2PIPID";
  if (B2MUPID == fMode) cMode = "B2MUPID";

  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFBu2MuTauK configuration" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  mukaVtexProb                " << fMuKaVtxProb << endl;
  cout << "---  mukaMaxDoca                 " << fMuKaMaxDoca << endl;
  cout << "---  tauVtxProb                  " << fTauVtxProb  << endl;
  cout << "---  tauMaxDoca                  " << fTauMaxDoca  << endl;
  cout << "---  tauMaxIP                    " << fTauMaxIP  << endl;
  cout << "---  MODE                        " << cMode  << endl;
  cout << "----------------------------------------------------------------------" << endl;
}


// ----------------------------------------------------------------------
void HFBu2MuTauK::analyze(const Event& iEvent, const EventSetup& iSetup) {
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
    if (B2MUKTAU == fMode) {
      genMatchB2MUKTAU();
    }
    if (B2MUPID == fMode) {
      genMatchB2MUPID();
    }
    if (B2PIPID == fMode) {
      genMatchB2PIPID();
    }
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
      if (fVerbose > 0) {
	cout << "MC truth cand found! ";
	for (unsigned int i = 0; i < truthTracks.size(); ++i) {
	  cout << truthTracks[i] << " ";
	}
	cout << endl;
      }
      bool filled = buildCandidate(truthTracks, fMcType, false);
      if (!filled) {
	if (fVerbose > 0) {
	  cout << "no MC truth cand built?!?!?!?!?!" << endl;
	}
      }
    } else {
      if (fVerbose > 0) {
	cout << "no MC truth cand found! ";
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
  vector<int> muonList;

  if (B2MUKTAU == fMode) {
    muonList = fListBuilder->getMuonList();
  }

  // -- for this mode we use the muonList variable but use getTrackList with the higher pT threshold
  if (B2PIPID == fMode) {
    muonList = fListBuilder->getTrackList();
  }
  if (fVerbose > 0) {
    cout << "==>HFBu2MuTauK> muonList  size: " << muonList.size() << endl;

    for (vector<int>::iterator it = muonList.begin(); it != muonList.end(); ++it) {
      cout << *it << " ";
    }
    cout << endl;
  }

  if (muonList.size() < 1) return;
  fListBuilder->setCloseTracks(&muonList);
  fListBuilder->setMinPt(0.5);

  vector<int> trkList = fListBuilder->getTrackList();
  if (fVerbose > 0) {
    cout << "==>HFBu2MuTauK> trkList  size: " << trkList.size() << endl;
    for (vector<int>::iterator it = trkList.begin(); it != trkList.end(); ++it) {
      cout << *it << " ";
    }
    cout << endl;
  }
  if (trkList.size() < 1) return;

  HFCombinatorics tau2(fTracksHandle, fTTB, fVerbose);
  vector<HFCombinatorics::doublet> promptList;
  double docacut(fMuKaMaxDoca); // any security margin??
  if (B2MUKTAU == fMode) {
    //    promptList = a.combine(muonList, MMUON, trkList, MKAON, 0.6, 3.6, 1);
    tau2.combine2Tracks(promptList, muonList, MMUON, trkList, MKAON, 0.6, 3.6, docacut, -1);
  } else if (B2PIPID == fMode) {
    //    promptList = a.combine(muonList, MMUON, trkList, MPION, 0.2, 5.2, 1);
    tau2.combine2Tracks(promptList, muonList, MPION, trkList, MPION, 0.2, 3.6, docacut, +1);
  }
  if (fVerbose > 0) cout << "==>HFBu2MuTauK> promptList  size: " << promptList.size() << endl;
  if (promptList.size() < 1) return;

  // -- build 3-prong
  vector<HFCombinatorics::triplet> tauList2;
  docacut = fTauMaxDoca; // any security margin??
  if (B2MUKTAU == fMode) {
    tau2.combine3Tracks(tauList2, trkList, MPION, 0.4, 1.8, docacut);
  } else if (B2PIPID == fMode) {
    tau2.combineDp2Km2Pip(tauList2, trkList, 1.75, 2.05, docacut);
  }

  if (fVerbose > 0) cout << "==>HFBu2MuTauK> tripletList  size: " << tauList2.size() << endl;
  if (tauList2.size() < 1) return;

  // -- Loop over (mu,K) pairs and all associated triplets
  vector<TransientTrack> vtt;
  vector<int> mkid;
  for (vector<HFCombinatorics::doublet>::iterator promptIt = promptList.begin(); promptIt != promptList.end(); ++promptIt) {
    int iMuon = promptIt->p1();
    int iKaon = promptIt->p2();

    if (fVerbose > 0) cout << Form("%d/%d", iMuon, iKaon) << " muka indices " << endl;
    // -- check basic pT requirements for muon and kaon
    reco::TrackBaseRef muTrackView(fTracksHandle, iMuon);
    reco::Track tMuon(*muTrackView);
    if (tMuon.pt() < fMuonPt) {
      if (fVerbose > 0) cout << Form("%d/%d", iMuon, iKaon) << " muon pT = " << tMuon.pt() << endl;
      continue;
    }

    TrackBaseRef kaTrackView(fTracksHandle, iKaon);
    Track tKaon(*kaTrackView);
    if (tKaon.pt() < fTrackPt) {
      if (fVerbose > 0) cout << Form("%d/%d", iMuon, iKaon) << " kaon pT = " << tKaon.pt() << endl;
      continue;
    }

    TransientTrack mTT = fTTB->build(*muTrackView);
    TransientTrack kTT = fTTB->build(*kaTrackView);
    vtt.clear(); mkid.clear();
    if (B2MUKTAU == fMode) {
      vtt.push_back(mTT); mkid.push_back(13);
      vtt.push_back(kTT); mkid.push_back(321);
    }
    if (B2PIPID == fMode) {
      vtt.push_back(mTT); mkid.push_back(211);
      vtt.push_back(kTT); mkid.push_back(211);
    }

    // -- check that the (mu,Ka) vertex is (1) meaningful and (2) well separated from best PV
    RefCountedKinematicTree mukaTree = fitTree(vtt, mkid);
    if (mukaTree->isEmpty()) {
      if (fVerbose > 0) cout << Form("%d/%d", iMuon, iKaon) << " muka tree empty" << endl;
      continue;
    }
    RefCountedKinematicVertex mukaVtx = mukaTree->currentDecayVertex();
    double mkprob = 0.0;

    if ((mukaVtx->chiSquared() >= 0.0) && (mukaVtx->degreesOfFreedom() > 0)) {
      mkprob = TMath::Prob(mukaVtx->chiSquared(), mukaVtx->degreesOfFreedom());
    }

    if (mkprob < fMuKaVtxProb) {
      if (fVerbose > 0) cout << Form("%d/%d", iMuon, iKaon) << " muka vtx probability = " << mkprob << endl;
      continue;
    }
    pair<int, double> mukaPV = findBestPV(vtt, mukaTree, fVertexCollection, fMagneticField);

    Vertex mypVertex = fVertexCollection[mukaPV.first];
    Vertex mysVertex = mkVertex(mukaTree);
    pair<double, double> psSep = vtxSeparation(mypVertex, mysVertex);
    if (psSep.first/psSep.second < fFlsxy) {
      if (fVerbose > 0) cout << Form("%d/%d", iMuon, iKaon) << " vtxsep(PV, muka) = " << psSep.first << " +/- " << psSep.second << endl;
      continue;
    }
    TVector3 mkplab = TVector3(mukaTree->currentParticle()->currentState().globalMomentum().x(),
			       mukaTree->currentParticle()->currentState().globalMomentum().y(),
			       mukaTree->currentParticle()->currentState().globalMomentum().z());
    const TVector3 p0(mypVertex.position().x(), mypVertex.position().y(), mypVertex.position().z());
    const TVector3 p1(mysVertex.position().x(), mysVertex.position().y(), mysVertex.position().z());
    TVector3 pDiff10 = p1 - p0;
    double cosa   = TMath::Cos(pDiff10.Angle(mkplab));
    if (cosa < 0.) {
      if (fVerbose > 0) cout << Form("%d/%d", iMuon, iKaon) <<  " p(muka) opposite to B flight direction, cosa = " << cosa << endl;
      continue;
    }


    double mkmass = mukaTree->currentParticle()->currentState().mass();
    if (fVerbose > 0) cout << Form("%d/%d", iMuon, iKaon) << " muka mass(" << iMuon << "/" << iKaon << "):  mvtx = " << mkmass
			   << " and pT(mu) = " << tMuon.pt() << " pT(K) = " << tKaon.pt()
			   << " vtx prob = " << mkprob
			   << endl;

    // -- now add 3-prongs!
    for (unsigned int it = 0; it < tauList2.size(); ++it) {
      // -- check against overlap with muons and/or other 3-hadron combinations already written to the tree
      if (tauList2[it].p1() == iMuon || tauList2[it].p1() == iKaon) continue;
      if (tauList2[it].p2() == iMuon || tauList2[it].p2() == iKaon) continue;
      if (tauList2[it].p3() == iMuon || tauList2[it].p3() == iKaon) continue;

      vector<int> allTrks;
      allTrks.push_back(iMuon);
      allTrks.push_back(iKaon);
      allTrks.push_back(tauList2[it].p1());
      allTrks.push_back(tauList2[it].p2());
      allTrks.push_back(tauList2[it].p3());

      TAnaCand *filledCand = buildCandidate(allTrks, fType, true);
      if (filledCand) {
	if (fVerbose > 0) {
	  cout << " .. filled!" << endl;
	  filledCand->dump();
	}
      }
    }
  }
}


// ----------------------------------------------------------------------
TAnaCand* HFBu2MuTauK::buildCandidate(std::vector<int> trkIdx, int type, bool cut) {
  if (fVerbose > 0) cout << "==HFBu2MuTauK::buildCandidate> " << type
			 << " " << trkIdx[0] << "/" << trkIdx[1] << "/" << trkIdx[2] << "/" << trkIdx[3] << "/" << trkIdx[4]
			 << endl;
  int iMuon = trkIdx[0];
  int iKaon = trkIdx[1];
  int iPion1 = trkIdx[2];
  int iPion3 = trkIdx[3];
  int iPion2 = trkIdx[4];

  vector<TransientTrack> vtt;
  vector<int> mkid, h3id;

  reco::TrackBaseRef muTrackView(fTracksHandle, iMuon);
  reco::Track tMuon(*muTrackView);
  if (cut && tMuon.pt() < fMuonPt) return 0;

  TrackBaseRef kaTrackView(fTracksHandle, iKaon);
  Track tKaon(*kaTrackView);
  if (cut && tKaon.pt() < fTrackPt) return 0;

  if (cut && (fMaxD0 < 90.) && (tMuon.d0() > fMaxD0)) return 0;
  if (cut && (fMaxDz < 90.) && (tMuon.dz() > fMaxDz)) return 0;
  if (cut && (fMaxD0 < 90.) && (tKaon.d0() > fMaxD0)) return 0;
  if (cut && (fMaxDz < 90.) && (tKaon.dz() > fMaxDz)) return 0;

  TransientTrack mTT = fTTB->build(*muTrackView);
  TransientTrack kTT = fTTB->build(*kaTrackView);
  vtt.clear(); mkid.clear();
  if (B2MUKTAU == fMode) {
    vtt.push_back(mTT); mkid.push_back(13);
    vtt.push_back(kTT); mkid.push_back(321);
  }
  if (B2PIPID == fMode) {
    vtt.push_back(mTT); mkid.push_back(211);
    vtt.push_back(kTT); mkid.push_back(211);
  }
  RefCountedKinematicTree mukaTree = fitTree(vtt, mkid);
  if (mukaTree->isEmpty()) return 0;
  RefCountedKinematicVertex mukaVtx = mukaTree->currentDecayVertex();
  double mkmass = mukaTree->currentParticle()->currentState().mass();
  double mkprob = 0.0;

  TVector3 mkplab = TVector3(mukaTree->currentParticle()->currentState().globalMomentum().x(),
			     mukaTree->currentParticle()->currentState().globalMomentum().y(),
			     mukaTree->currentParticle()->currentState().globalMomentum().z());
  if ((mukaVtx->chiSquared() >= 0.0) && (mukaVtx->degreesOfFreedom() > 0)) {
    mkprob = TMath::Prob(mukaVtx->chiSquared(), mukaVtx->degreesOfFreedom());
  }

  if (cut && mkprob < fMuKaVtxProb) {
    if (fVerbose > 0) cout << "muKa vtx probability too low, skipping" << endl;
    return 0;
  }
  pair<int, double> mukaPV = findBestPV(vtt, mukaTree, fVertexCollection, fMagneticField);

  Vertex mypVertex = fVertexCollection[mukaPV.first];
  Vertex mysVertex = mkVertex(mukaTree);
  pair<double, double> psSep = vtxSeparation(mypVertex, mysVertex);
  if (cut) {
    if (psSep.first/psSep.second < fFlsxy) {
      if (fVerbose > 0) cout << "vtxsep(PV, muka) = " << psSep.first << " +/- " << psSep.second << endl;
      return 0;
    }
  }

  TrackBaseRef pi1TrackView(fTracksHandle, iPion1);
  reco::Track pi1T(*pi1TrackView);
  TransientTrack pi1TT = fTTB->build(*pi1TrackView);

  TrackBaseRef pi2TrackView(fTracksHandle, iPion2);
  reco::Track pi2T(*pi2TrackView);
  TransientTrack pi2TT = fTTB->build(*pi2TrackView);

  TrackBaseRef pi3TrackView(fTracksHandle, iPion3);
  reco::Track pi3T(*pi3TrackView);
  TransientTrack pi3TT = fTTB->build(*pi3TrackView);

  if (cut && pi1T.pt() < fTrackPt)  return 0;
  if (cut && pi2T.pt() < fTrackPt)  return 0;
  //NO!      if (pi3T.pt() < fTrackPt)  return;

  // -- calculate maxDoca and potentially cut on it.
  double maxDocaHad(-99.), maxDocaMuKa(-99.);
  TwoTrackMinimumDistance md;
  FreeTrajectoryState freestates[3];
  freestates[0]=pi1TT.initialFreeState();
  freestates[1]=pi2TT.initialFreeState();
  freestates[2]=pi3TT.initialFreeState();
  double doca(99.);
  for (int i = 0; i < 2; ++i) {
    md.calculate(freestates[i], freestates[i+1]);
    doca = md.distance();
    if (doca > maxDocaHad) {
      maxDocaHad = doca;
    }
  }
  if (cut && (maxDocaHad > fTauMaxDoca)) {
    if (fVerbose > 0) cout << "tau cand with bad doca: " << doca << endl;
    return 0;
  }

  md.calculate(mTT.initialFreeState(), kTT.initialFreeState());
  maxDocaMuKa = md.distance();
  if (cut && (maxDocaMuKa > fMuKaMaxDoca)) {
    if (fVerbose > 0) cout << "muka cand with bad doca: " << maxDocaMuKa << endl;
    return 0;
  }

  int qtot = pi1T.charge() + pi2T.charge() + pi3T.charge();
  if (cut && TMath::Abs(qtot) != 1) {
    if (fVerbose > 0) cout << "qtot(3 hadrons) = " << qtot << ", skipping" << endl;
    return 0;
  }

  vtt.clear();
  h3id.clear();
  if (B2MUKTAU == fMode) {
    // -- For tau+ -> pi- pi+ pi+ decays, all hadrons are pions. Ignore those decays with kaons included!?
    vtt.push_back(pi1TT); h3id.push_back(211);
    vtt.push_back(pi2TT); h3id.push_back(211);
    vtt.push_back(pi3TT); h3id.push_back(211);
  }
  if (B2PIPID == fMode) {
    // -- For D+ -> K- pi+ pi+ decays, the same-sign particles are the pions. DCS decays are ignored!
    vtt.push_back(pi1TT); h3id.push_back(321);
    vtt.push_back(pi2TT); h3id.push_back(211);
    vtt.push_back(pi3TT); h3id.push_back(211);
  }

  RefCountedKinematicTree tauTree = fitTree(vtt, h3id);
  if (tauTree->isEmpty()) return 0;
  double doca3dhad =  doca2Vtx(tauTree, mysVertex).second.value();
  if (cut && (doca3dhad > fTauMaxIP)) {
    if (fVerbose > 0) cout << "tau 3D IP = " << doca3dhad << " larger than " << fTauMaxIP << endl;
    return 0;
  }
  RefCountedKinematicVertex tauVtx = tauTree->currentDecayVertex();
  TVector3 tauplab = TVector3(tauTree->currentParticle()->currentState().globalMomentum().x(),
			      tauTree->currentParticle()->currentState().globalMomentum().y(),
			      tauTree->currentParticle()->currentState().globalMomentum().z());
  double tauprob = 0.0;
  if ((tauVtx->chiSquared() >= 0.0) && (tauVtx->degreesOfFreedom() > 0)) {
    tauprob = TMath::Prob(tauVtx->chiSquared(), tauVtx->degreesOfFreedom());
  }

  if (cut && (tauprob < fTauVtxProb)) {
    if (fVerbose > 0) cout << "tau vtx probability too low (" << tauprob << ") for pion idxs = "
			   << iPion1 << "/" << iPion2 << "/" << iPion3
			   << " skipping" << endl;
    return 0;
  }

  Vertex mytVertex = mkVertex(tauTree);
  pair<double, double> stSep = vtxSeparation(mytVertex, mysVertex);

  if (cut) {
    if (stSep.first/stSep.second < fFlsxy) {
      if (fVerbose > 0) cout << "vtxsep(muka, had) = " << stSep.first << " +/- " << stSep.second << endl;
      return 0;
    }
  }

  double taumass = tauTree->currentParticle()->currentState().mass();
  double massCut(1.5);
  if (B2PIPID == fMode) {
    massCut = 2.1;
  }
  if (cut && (taumass > massCut)) {
    if (fVerbose > 0) cout << "had mass = " << taumass << " larger than " << massCut << endl;
    return 0;
  }

  // -- flight directions from vertices
  const TVector3 p0(mypVertex.position().x(), mypVertex.position().y(), mypVertex.position().z());
  const TVector3 p1(mysVertex.position().x(), mysVertex.position().y(), mysVertex.position().z());
  TVector3 pDiff10 = p1 - p0;
  const TVector3 p2(mytVertex.position().x(), mytVertex.position().y(), mytVertex.position().z());
  TVector3 pDiff21 = p2-p1;
  double cosa1     = TMath::Cos(pDiff10.Angle(pDiff21));
  double cosa2     = TMath::Cos(tauplab.Angle(pDiff21));
  if (cut && (cosa1 < 0.)) {
    if (fVerbose > 0) cout << "cosa1 = " << cosa1 << " smaller than 0" << endl;
    return 0;
  }
  if (cut && (cosa2 < 0.)) {
    if (fVerbose > 0) cout << "cosa2 = " << cosa2 << " smaller than 0" << endl;
    return 0;
  }

  if (fVerbose > 0) {
    TLorentzVector tlv,  muka, mu, ka, bu, pi1, pi2, pi3;
    if (B2MUKTAU == fMode) {
      mu.SetPtEtaPhiM(tMuon.pt(), tMuon.eta(), tMuon.phi(), MMUON);
      ka.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MKAON);
      pi1.SetPtEtaPhiM(pi1T.pt(), pi1T.eta(), pi1T.phi(), MPION);
      pi2.SetPtEtaPhiM(pi2T.pt(), pi2T.eta(), pi2T.phi(), MPION);
      pi3.SetPtEtaPhiM(pi3T.pt(), pi3T.eta(), pi3T.phi(), MPION);
    } else if (B2PIPID == fMode) {
      mu.SetPtEtaPhiM(tMuon.pt(), tMuon.eta(), tMuon.phi(), MPION);
      ka.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MPION);
      pi1.SetPtEtaPhiM(pi1T.pt(), pi1T.eta(), pi1T.phi(), MKAON);
      pi2.SetPtEtaPhiM(pi2T.pt(), pi2T.eta(), pi2T.phi(), MPION);
      pi3.SetPtEtaPhiM(pi3T.pt(), pi3T.eta(), pi3T.phi(), MPION);
    }
    muka = mu + ka;
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
  TAnaCand *pCand = fillCand(type, mukaTree, mkidx, mkid, mypVertex);
  pCand->fPvIdx = mukaPV.first;
  pCand->fMaxDoca = maxDocaMuKa;

  // -- fill doca to PV
  gHFEvent->getSigTrack(pCand->fSig1)->fDouble1 = doca2Vtx(mTT, mypVertex).second.value();
  gHFEvent->getSigTrack(pCand->fSig2)->fDouble1 = doca2Vtx(kTT, mypVertex).second.value();

  vector<int> h3idx;
  h3idx.push_back(iPion1);
  h3idx.push_back(iPion2);
  h3idx.push_back(iPion3);
  TAnaCand *tauCand = fillCand(type+100, tauTree, h3idx, h3id, mysVertex);
  tauCand->fPvIdx = -1;
  tauCand->fMaxDoca = maxDocaHad;

  // -- now fix the pointers to daughter and mother
  pCand->fDau1  = tauCand->fIndex;
  pCand->fDau2  = tauCand->fIndex;
  tauCand->fMom = pCand->fIndex;

  gHFEvent->getSigTrack(tauCand->fSig1)->fDouble1   =  doca2Vtx(pi1TT, mypVertex).second.value();
  gHFEvent->getSigTrack(tauCand->fSig1+1)->fDouble1 =  doca2Vtx(pi2TT, mypVertex).second.value();
  gHFEvent->getSigTrack(tauCand->fSig2)->fDouble1   =  doca2Vtx(pi3TT, mypVertex).second.value();

  // -- calculate tau and B effective decay time
  double taucosa   = TMath::Cos(tauCand->fPlab.Angle(pDiff21));
  double massOverC = tauCand->fMass/TMath::Ccgs();
  tauCand->fTau3d = stSep.first / tauCand->fPlab.Mag() * taucosa * massOverC;
  massOverC = MBPLUS/TMath::Ccgs();

  double bcosa  = TMath::Cos(plab.Angle(pDiff10));
  pCand->fTau3d = psSep.first / pCand->fPlab.Mag() * bcosa * massOverC;

  if (fVerbose > 100) {
    cout << "HFListCand dump: " << endl;
    pCand->dump();
    for (int is = pCand->fDau1; is <= pCand->fDau2; ++is) {
      TAnaCand *tc = gHFEvent->getCand(is);
      tc->dump();
      for (int is = tc->fSig1; is <= tc->fSig2; ++is) {
	cout << "sigIdx = " << is << ": "; gHFEvent->getSigTrack(is)->dump();
      }
    }
    for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
      cout << "sigIdx = " << is << ": ";   gHFEvent->getSigTrack(is)->dump();
    }
  }
  return pCand;
}


// ----------------------------------------------------------------------
TAnaCand* HFBu2MuTauK::fillCand(int type, RefCountedKinematicTree &tree, std::vector<int> trkIdx, std::vector<int> trkId, reco::Vertex &refVtx) {
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
    pTrack->fMCID = trkId[i];
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

  pair<double, double> vtxSep = vtxSeparation(myVertex, refVtx);
  pCand->fVtx.fD3d     = vtxSep.first;
  pCand->fVtx.fD3dE    = vtxSep.second;

  pCand->fPvIP3d  = doca2Vtx(tree, refVtx).second.value();
  pCand->fPvIP3dE = doca2Vtx(tree, refVtx).second.error();

  // -- calculate flight direction:
  //    o B:    PV          -> Vtx(Mu, Ka)
  //    o tau:  Vtx(Mu, Ka) -> Vtx(3-prong)

  const TVector3 p0(refVtx.position().x(), refVtx.position().y(), refVtx.position().z());
  const TVector3 d = p1 - p0;

  RefCountedKinematicParticle kinParticle = tree->currentParticle();
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TrajectoryStateOnSurface tsos = extrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),
							   RecoVertex::convertPos(refVtx.position()));
  GlobalVector direction(d.X(), d.Y(), d.Z());
  pair<bool, Measurement1D> ip = IPTools::signedDecayLength3D(tsos, direction, refVtx);
  pCand->fPvLip  = ip.second.value();
  pCand->fPvLipE = ip.second.error();

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
void HFBu2MuTauK::genMatchB2PIPID() {
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
	if ((211 == TMath::Abs(pD->fID)) && (iC == pD->fMom1)) {
	  if (0 == fpGenMu) {
	    fpGenMu = pD;
	  } else if (0 == fpGenKa) {
	    fpGenKa = pD;
	  }
	}
	if (411 == TMath::Abs(pD->fID)) {
	  fpGenTau = pD;
	  for (int iT = pD->fDau1; iT <= pD->fDau2; ++iT) {
	    pT = gHFEvent->getGenCand(iT);
	    if (0 == pT->fQ) continue;
	    if (321 == TMath::Abs(pT->fID)) {
	      fpGenHad1 = pT;
	    } else if (211 == TMath::Abs(pT->fID)) {
	      if (0 == fpGenHad2) {
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
    if (fVerbose > 0) {
      cout << "==> HFBu2MuTauK::genMatchB2PIPID> successfully matched gen decay" << endl;
    }
    // -- sort pions ONLY according to pT of hadronic tracks
    double pt2 = fpGenHad2->fP.Perp();
    double pt3 = fpGenHad3->fP.Perp();
    TGenCand *tgc1(0), *tgc2(0), *tgc3(0);
    if (pt2 > pt3) {
      tgc2 = fpGenHad2;
      tgc3 = fpGenHad3;
    } else {
      tgc2 = fpGenHad3;
      tgc3 = fpGenHad2;
    }
    fpGenHad2 = tgc2;
    fpGenHad3 = tgc3;
  }
}


// ----------------------------------------------------------------------
void HFBu2MuTauK::genMatchB2MUPID() {
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
	if (211 == TMath::Abs(pD->fID)) fpGenKa = pD;
	if (13 == TMath::Abs(pD->fID)) fpGenMu = pD;
	if (411 == TMath::Abs(pD->fID)) {
	  fpGenTau = pD;
	  for (int iT = pD->fDau1; iT <= pD->fDau2; ++iT) {
	    pT = gHFEvent->getGenCand(iT);
	    if (0 != pT->fQ) {
	      if (321 == TMath::Abs(pT->fID)) {
		fpGenHad1 = pD;
	      } else if (211 == TMath::Abs(pT->fID)) {
		if (0 == fpGenHad2) {
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
      }
      if (fpGenB && fpGenMu && fpGenKa && fpGenTau && fpGenHad1 && fpGenHad2 && fpGenHad3) {
	goodMatch = true;
	break;
      }
    }
  }
  if (goodMatch) {
    // -- sort pions ONLY according to pT of hadronic tracks
    double pt2 = fpGenHad2->fP.Perp();
    double pt3 = fpGenHad3->fP.Perp();
    TGenCand *tgc1(0), *tgc2(0), *tgc3(0);
    if (pt2 > pt3) {
      tgc2 = fpGenHad2;
      tgc3 = fpGenHad3;
    } else {
      tgc2 = fpGenHad3;
      tgc3 = fpGenHad2;
    }
    fpGenHad2 = tgc2;
    fpGenHad3 = tgc3;
  }
}


// ----------------------------------------------------------------------
void HFBu2MuTauK::genMatchB2MUKTAU() {
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

// ----------------------------------------------------------------------
pair<bool, Measurement1D> HFBu2MuTauK::doca2Vtx(TrajectoryStateOnSurface &tsos, Vertex &vtx) {
  VertexDistance3D a3d;
  return IPTools::absoluteImpactParameter(tsos, vtx, a3d);

}

// ----------------------------------------------------------------------
pair<bool, Measurement1D> HFBu2MuTauK::doca2Vtx(TransientTrack &tt, Vertex &vtx) {
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TrajectoryStateOnSurface tsos = extrapolator.extrapolate(tt.initialFreeState(), RecoVertex::convertPos(vtx.position()));
  //  TrajectoryStateOnSurface tsos = extrapolator.extrapolate(tt.initialFreeState(), vtx);
  return doca2Vtx(tsos, vtx);
}

// ----------------------------------------------------------------------
pair<bool, Measurement1D> HFBu2MuTauK::doca2Vtx(RefCountedKinematicTree &tree, Vertex &vtx) {
  RefCountedKinematicParticle kinParticle = tree->currentParticle();
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TrajectoryStateOnSurface tsos = extrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),
							   RecoVertex::convertPos(vtx.position()));
  return doca2Vtx(tsos, vtx);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HFBu2MuTauK);
