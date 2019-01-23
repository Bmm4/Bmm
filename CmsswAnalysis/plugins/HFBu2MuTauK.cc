#include "FWCore/Framework/interface/MakerMacros.h"
#include "HFBu2MuTauK.h"

#include <iostream>
#include <utility>

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TAnaCand.hh"
#include "Bmm/RootAnalysis/rootio/TGenCand.hh"

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


// ----------------------------------------------------------------------
HFBu2MuTauK::HFBu2MuTauK(const ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
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

  fListBuilder->setMinPt(fMuonPt);
  vector<int> muonList = fListBuilder->getMuonList();

  fListBuilder->setMinPt(fTrackPt);
  fListBuilder->setMaxDocaToTracks(0.1);
  if (fVerbose > 0) cout << "==>HFBu2MuTauK> muonList  size: " << muonList.size() << endl;
  fListBuilder->setCloseTracks(&muonList);
  vector<int> trkList = fListBuilder->getTrackList();

  vector<pair<int, TLorentzVector> > piList;
  TLorentzVector tlv;
  for (unsigned int i = 0; i < trkList.size(); ++i) {
    TrackBaseRef trackView(fTracksHandle, trkList[i]);
    Track trk(*trackView);
    tlv.SetPtEtaPhiM(trk.pt(), trk.eta(), trk.phi(), MPION);
    piList.push_back(make_pair(trkList[i], tlv));
  }

  if (fVerbose > 0) cout << "==>HFBu2MuTauK> trkList  size: " << trkList.size() << endl;
  HFTwoParticleCombinatoricsNew a(fTracksHandle, 1);
  HFTwoParticleCombinatoricsSet promptList = a.combine(muonList, MMUON, trkList, MKAON, 0.6, 3.6, 1 );

  if (fVerbose > 0) cout << "==>HFBu2MuTauK> promptList  size: " << promptList.size() << endl;

  // -- Build mu + ka
  TAnaCand *pCand(0);
  double mass(0.1);
  TLorentzVector muka, mu, ka, bu, pi1, pi2, pi3;
  //  HFDecayTree theTree(300521, true, MBPLUS, false, -1.0, true);
  vector<TransientTrack> vtt;
  vector<int> vid;
  for (HFTwoParticleCombinatoricsNew::iterator promptIt = promptList.begin(); promptIt != promptList.end(); ++promptIt) {
    int iMuon = promptIt->first;
    int iKaon = promptIt->second;

    reco::TrackBaseRef muTrackView(fTracksHandle, iMuon);
    reco::Track tMuon(*muTrackView);
    if (tMuon.pt() < fMuonPt)  continue;
    mu.SetPtEtaPhiM(tMuon.pt(), tMuon.eta(), tMuon.phi(), MMUON);

    TrackBaseRef kaTrackView(fTracksHandle, iKaon);
    Track tKaon(*kaTrackView);
    if (tKaon.pt() < fTrackPt)  continue;
    ka.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MKAON);

    if ((fMaxD0 < 90.) && (tMuon.d0() > fMaxD0)) continue;
    if ((fMaxDz < 90.) && (tMuon.dz() > fMaxDz)) continue;
    if ((fMaxD0 < 90.) && (tKaon.d0() > fMaxD0)) continue;
    if ((fMaxDz < 90.) && (tKaon.dz() > fMaxDz)) continue;

    muka = mu + ka;

    TransientTrack mTT = fTTB->build(*muTrackView);
    TransientTrack kTT = fTTB->build(*kaTrackView);
    vtt.clear(); vid.clear();
    vtt.push_back(mTT); vid.push_back(13);
    vtt.push_back(kTT); vid.push_back(321);
    RefCountedKinematicTree mukaTree = fitTree(vtt, vid);
    RefCountedKinematicVertex mukaVtx = mukaTree->currentDecayVertex();
    double mkmass = mukaTree->currentParticle()->currentState().mass();
    double prob = 0.0;
    if ((mukaVtx->chiSquared() >= 0.0) && (mukaVtx->degreesOfFreedom() > 0)) {
      prob = TMath::Prob(mukaVtx->chiSquared(), mukaVtx->degreesOfFreedom());
    }

    cout << "muka mass(" << iMuon << "/" << iKaon << "): m4v = " << muka.M() << " mvtx = " << mkmass
	 << " and pT(mu) = " << tMuon.pt() << " pT(K) = " << tKaon.pt()
	 << " vtx prob = " << prob
	 << endl;

    if (prob < fMuKaVtxProb) {
      if (fVerbose > 0) cout << "muKa vtx probability too low, skipping" << endl;
      continue;
    }
    pair<int, double> mukaPV = findBestPV(vtt, mukaTree, fVertexCollection, fMagneticField);

    Vertex mypVertex = fVertexCollection[mukaPV.first];
    Vertex mysVertex = mkVertex(mukaTree);
    pair<double, double> psSep = vtxSeparation(mypVertex, mysVertex);
    cout << "vtxsep(PV, muka) = " << psSep.first << " +/- " << psSep.second << endl;

    // -- build 3-prong
    vector<triplet> tauList;
    HFThreeParticleCombinatorics tau(10);
    tau.combine(tauList, piList, 0.3, 2.0, 0.6, 1.1);
    vector<triplet> filledTriplets;
    for (unsigned int it = 0; it < tauList.size(); ++it) {
      // -- check against overlap with muons and/or other 3-hadron combinations already written to the tree
      if (tauList[it].pi1() == iMuon || tauList[it].pi1() == iKaon) continue;
      if (tauList[it].pi2() == iMuon || tauList[it].pi2() == iKaon) continue;
      if (tauList[it].pi3() == iMuon || tauList[it].pi3() == iKaon) continue;
      for (unsigned ik = 0; ik < filledTriplets.size(); ++ik) {
	if (false == tauList[it].different(filledTriplets[ik])) {
	  if (fVerbose > 0) {
	    cout << "not filling " << tauList[it].pi1() << "/" << tauList[it].pi2() << "/" << tauList[it].pi3()
		 << " because have already " << filledTriplets[ik].pi1() << "/" << filledTriplets[ik].pi2() << "/" << filledTriplets[ik].pi3()
		 << endl;
	  }
	  continue;
	}
      }
      TrackBaseRef pi1TrackView(fTracksHandle, tauList[it].pi1());
      reco::Track pi1T(*pi1TrackView);
      TransientTrack pi1TT = fTTB->build(*pi1TrackView);

      TrackBaseRef pi2TrackView(fTracksHandle, tauList[it].pi2());
      reco::Track pi2T(*pi2TrackView);
      TransientTrack pi2TT = fTTB->build(*pi2TrackView);

      TrackBaseRef pi3TrackView(fTracksHandle, tauList[it].pi3());
      reco::Track pi3T(*pi3TrackView);
      TransientTrack pi3TT = fTTB->build(*pi3TrackView);

      if (pi1T.pt() < fTrackPt)  continue;
      if (pi2T.pt() < fTrackPt)  continue;
      //NO!      if (pi3T.pt() < fTrackPt)  continue;

      pi1.SetPtEtaPhiM(pi1T.pt(), pi1T.eta(), pi1T.phi(), MPION);
      pi2.SetPtEtaPhiM(pi2T.pt(), pi2T.eta(), pi2T.phi(), MPION);
      pi3.SetPtEtaPhiM(pi3T.pt(), pi3T.eta(), pi3T.phi(), MPION);
      tlv = pi1 + pi2 + pi3;

      vtt.clear(); vid.clear();
      vtt.push_back(pi1TT); vid.push_back(211);
      vtt.push_back(pi2TT); vid.push_back(211);
      vtt.push_back(pi3TT); vid.push_back(211);
      RefCountedKinematicTree tauTree = fitTree(vtt, vid);
      RefCountedKinematicVertex tauVtx = tauTree->currentDecayVertex();
      double tauprob = 0.0;
      if ((tauVtx->chiSquared() >= 0.0) && (tauVtx->degreesOfFreedom() > 0)) {
	tauprob = TMath::Prob(tauVtx->chiSquared(), tauVtx->degreesOfFreedom());
      }

      if (tauprob < fTauVtxProb) {
	if (fVerbose > 0) cout << "tau vtx probability too low, skipping" << endl;
	continue;
      }

      Vertex mytVertex = mkVertex(tauTree);
      pair<double, double> stSep = vtxSeparation(mysVertex,  mytVertex);

      double taumass = tauTree->currentParticle()->currentState().mass();
      cout << "h3mass(" << tauList[it].pi1() << "/" << tauList[it].pi2() << "/" << tauList[it].pi3()
	   << "): m4v = " << tlv.M() << " mvtx = " << taumass
	   << " vtx prob = " << tauprob
	   << " separation = " << stSep.first << " +/- " << stSep.second << endl
	   << " FILLING CANDIDATE"
	   << endl;

      filledTriplets.push_back(tauList[it]);

      TAnaCand *pCand  = gHFEvent->addCand();
      pCand->fType  = fType;
      pCand->fSig1  = gHFEvent->nSigTracks();
      pCand->fSig2  = pCand->fSig1 + 2 - 1;
      // -- only the visible part!(?)
      pCand->fPlab  = TVector3(mukaTree->currentParticle()->currentState().globalMomentum().x(),
			       mukaTree->currentParticle()->currentState().globalMomentum().y(),
			       mukaTree->currentParticle()->currentState().globalMomentum().z());
      pCand->fMass  = mukaTree->currentParticle()->currentState().mass();
      pCand->fMassC = (muka + tlv).M();
      pCand->fQ     = 0;
      // pCand->fDau1 set below

      int gidx(-1);
      const reco::BeamSpot *pBeamSpot = &fBeamSpot;
      TAnaTrack *pTrack(0);
      TSimpleTrack *sTrack(0);

      // -- muon
      pTrack = gHFEvent->addSigTrack();
      fillSigTrack(pTrack, iMuon, fTracksHandle, &fVertexCollection, fMuonCollection, pBeamSpot);
      // -- kaon
      pTrack = gHFEvent->addSigTrack();
      fillSigTrack(pTrack, iKaon, fTracksHandle, &fVertexCollection, fMuonCollection, pBeamSpot);



      TAnaCand *tauCand = gHFEvent->addCand();
      pCand->fDau1  = tauCand->fIndex;
      pCand->fDau2  = tauCand->fIndex;
      tauCand->fMom = pCand->fIndex;

      tauCand->fPlab  = TVector3(tauTree->currentParticle()->currentState().globalMomentum().x(),
			       tauTree->currentParticle()->currentState().globalMomentum().y(),
			       tauTree->currentParticle()->currentState().globalMomentum().z());

      tauCand->fMass = taumass;

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
      TVector3 plab = pCand->fPlab + tauCand->fPlab;
      double bcosa  = plab.Dot(pDiff) / (plab.Mag() * pDiff.Mag());
      pCand->fTau3d = psSep.first / pCand->fPlab.Mag() * bcosa * massOverC;

      tauCand->fSig1 = gHFEvent->nSigTracks();
      tauCand->fSig2 = pCand->fSig1 + 3 - 1;

      // -- pions
      pTrack = gHFEvent->addSigTrack();
      fillSigTrack(pTrack, tauList[it].pi1(), fTracksHandle, &fVertexCollection, fMuonCollection, pBeamSpot);
      pTrack = gHFEvent->addSigTrack();
      fillSigTrack(pTrack, tauList[it].pi2(), fTracksHandle, &fVertexCollection, fMuonCollection, pBeamSpot);
      pTrack = gHFEvent->addSigTrack();
      fillSigTrack(pTrack, tauList[it].pi3(), fTracksHandle, &fVertexCollection, fMuonCollection, pBeamSpot);

      // -- vertices
      pCand->fVtx.setInfo(mukaVtx->chiSquared(), mukaVtx->degreesOfFreedom(), prob, 0, -1);
      pCand->fVtx.fPoint   = p1;
      pCand->fVtx.fD3d     = psSep.first;
      pCand->fVtx.fD3dE    = psSep.second;

      tauCand->fVtx.setInfo(tauVtx->chiSquared(), tauVtx->degreesOfFreedom(), tauprob, 0, -1);
      tauCand->fVtx.fPoint = p2;
      tauCand->fVtx.fD3d   = stSep.first;
      tauCand->fVtx.fD3dE  = stSep.second;


      cout << "HFListCand dump: " << endl;
      pCand->dump();
      for (int is = pCand->fSig1; is <= pCand->fSig2; ++is) {
	gHFEvent->getSigTrack(is)->dump();
      }
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFBu2MuTauK);
