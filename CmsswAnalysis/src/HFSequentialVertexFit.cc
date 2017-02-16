// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFSequentialVertexFit
// ---------------------
//
// 2016/02/28 Urs Langenegger        modularization
// 2010/04/10 Christoph Naegeli      first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include "Bmm/CmsswAnalysis/interface/HFSequentialVertexFit.hh"
#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "CommonTools/Statistics/interface/ChiSquared.h"

#include <algorithm>

#define DBX 0

const static unsigned kNBR_CLOSED_TRACKS = 20;

extern TAna01Event *gHFEvent;

struct EmptyTreeError {
  EmptyTreeError() {}
};
struct MassNotFoundException {
  MassNotFoundException() {}
};
struct PVRefitException {
  PVRefitException() {}
};

using namespace std;
using namespace edm;
using namespace reco;


// ----------------------------------------------------------------------
HFSequentialVertexFit::HFSequentialVertexFit(Handle<View<Track> > hTracks,
					     const MuonCollection* muons,
					     const TransientTrackBuilder *TTB,
					     Handle<VertexCollection> pvCollection,
					     const MagneticField *field,
					     BeamSpot beamSpot,
					     int verbose,
					     bool removeCandTracksFromVtx) :
  fVerbose(verbose),
  fpTTB(TTB),
  fhTracks(hTracks),
  fPVCollection(pvCollection),
  fMuons(muons),
  fMagneticField(field),
  fBeamSpot(beamSpot),
  fRemoveCandTracksFromVtx(removeCandTracksFromVtx),
  fPvW8(0.6)
{}


// ----------------------------------------------------------------------
HFSequentialVertexFit::~HFSequentialVertexFit() {}



// ----------------------------------------------------------------------
void HFSequentialVertexFit::doFit(HFDecayTree *tree) {
  if (fVerbose > 5) cout << "==>HFSequentialVertexFit> doFit() for particleID = " << tree->particleID() << endl;

  try {
    tree->resetKinematicTree(1);
    if (0) tree->dump(1);
    // -- fit it
    fitTree(tree);
    // -- calculate all node variables
    calculateAll(tree);
    // -- calculate top tree variables wrt PV
    calculateStuff(tree, 0);
    if (0) tree->dump(1);
    // -- check cuts
    if (passAllCuts(tree)) {
      saveTree(tree);
    }
  } catch (cms::Exception &ex) {
    if (fVerbose > 0) cout << "==> HFSequentialVertexFit: cms exception caught: " << ex.what() << endl;
  } catch (VertexException &ex) {
    if (fVerbose > 0) cout << "==> HFSequentialVertexFit: vertex exception caught: " << ex.what() << endl;
  } catch (EmptyTreeError& ex) {
    if (fVerbose > 0) cout << "==> HFSequentialVertexFit: empty tree." << endl;
  } catch (PVRefitException& ex) {
    if (fVerbose > 0) cout << "==> HFSequentialVertexFit: unable to refit PV." << endl;
  }
}


// ----------------------------------------------------------------------
void HFSequentialVertexFit::calculateAll(HFDecayTree *tree) {
  if (DBX) cout << "now in HFSequentialVertexFit::calculateAll for tree  " << tree->particleID() << " at " << tree << endl;

  RefCountedKinematicTree kinTree = *(tree->getKinematicTree());
  kinTree->movePointerToTheTop();
  RefCountedKinematicVertex kinVertex = kinTree->currentDecayVertex();
  VertexState vState = kinVertex->vertexState();

  if ((tree->getVerticesEndIterator() - tree->getVerticesBeginIterator()) > 0) {
    for (HFDecayTreeIterator treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
      if (DBX) cout << "call HFSequentialVertexFit::calculateStuff for tree  " << treeIt->particleID() << " at " << &(*treeIt) << endl;
      calculateStuff(&(*treeIt), &vState);
    }
  }
}


// ----------------------------------------------------------------------
bool HFSequentialVertexFit::passAllCuts(HFDecayTree *tree) {
  if (DBX) cout << "now in HFSequentialVertexFit::passAllCuts for tree  " << tree->particleID()
		<< " at " << tree << " with tree->fTV = " << &(tree->fTV)
		<< endl;
  if ((tree->getVerticesEndIterator() - tree->getVerticesBeginIterator()) > 0) {
    for (HFDecayTreeIterator treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
      //      treeIt->dump();
      bool result = passAllCuts(&(*treeIt));
      if (!result) {
	if (DBX) cout << "HFSequentialVertexFit::passAllCuts returned false for tree  " << treeIt->particleID() << " at " << &(*treeIt)
		      << " with tree->fTV = " << &(tree->fTV)
		      << endl;
	return false;
      }
    }
  }
  //  tree->dump();
  if (tree->passAllCuts()) {
    if (DBX) cout << "HFSequentialVertexFit::passAllCuts returning true for tree  " << tree->particleID()
		  << " at " << tree  << " with tree->fTV = " << &(tree->fTV) << endl;
    return true;
  } else {
    if (DBX) cout << "HFSequentialVertexFit::passAllCuts returning false for tree  " << tree->particleID()
		  << " at " << tree  << " with tree->fTV = " << &(tree->fTV) << endl;
    return false;
  }
  if (DBX) cout << "HFSequentialVertexFit::passAllCuts reached end of function; returning false for tree  " << tree->particleID()
		<< " at " << tree  << " with tree->fTV = " << &(tree->fTV) << endl;
  return false;
}


// ----------------------------------------------------------------------
// if this node does not survive the nodeCut, then it returns false and the fitting sequence stops
bool HFSequentialVertexFit::fitTree(HFDecayTree *tree) {
  // -- propagate the fit message down the tree
  if (DBX) cout << "now in HFSequentialVertexFit::fitTree for tree  " << tree->particleID() << " at " << tree
		<< " &fTV.mass = " << &(tree->fTV.mass)
		<< endl;
  for (HFDecayTreeIterator treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
    fitTree(&(*treeIt));
  }

  // -- set up the kinParticleMap for 'tree'
  map<int,int> *kinParticleMap = tree->getKinParticleMap();
  kinParticleMap->clear();

  // -- add the particles from the tracks, we have to add them first for the KinematicConstrainedVertexFitter
  vector<track_entry_t> allTreeTracks = tree->getAllTracks(1);
  sort(allTreeTracks.begin(), allTreeTracks.end()); // sort such that all mass fit tracks are first
  KinematicParticleFactoryFromTransientTrack pFactory;
  vector<RefCountedKinematicParticle> kinParticles;
  int mass_constrained_tracks(0);
  for (vector<track_entry_t>::const_iterator trackIt = allTreeTracks.begin(); trackIt != allTreeTracks.end(); ++trackIt) {
    float sigma;
    float mass = getParticleMass(trackIt->particleID, &sigma);
    TrackBaseRef baseRef(fhTracks, trackIt->trackIx);
    (*kinParticleMap)[trackIt->trackIx] = kinParticles.size();
    kinParticles.push_back(pFactory.particle(fpTTB->build(*baseRef), mass, 0.0f, 0.0f, sigma));
    if (trackIt->massFit) mass_constrained_tracks++;
  }

  // -- add the particles from the sub-trees with vertexing
  for (HFDecayTreeIterator treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
    if (DBX) cout << "Start assembling kinematic particles for tree->particleID() = " << (int)tree->particleID() << endl;
    // get all the already fit particles
    addFittedParticles(&kinParticles, &(*treeIt));
  }

  // -- do the actual fit of this vertex
  if (DBX) cout << "do the actual fit for tree  " << tree->particleID() << " at " << tree << endl;
  RefCountedKinematicTree kinTree;
  if (tree->mass_tracks() > 0 && mass_constrained_tracks > 0) {
    KinematicConstrainedVertexFitter kcvFitter;
    auto_ptr<MultiTrackKinematicConstraint> tr_c(new MultiTrackMassKinematicConstraint(tree->mass_tracks(), mass_constrained_tracks));
    kinTree = kcvFitter.fit(kinParticles, &(*tr_c));
  } else {
    KinematicParticleVertexFitter kpvFitter;
    kinTree = kpvFitter.fit(kinParticles);
  }

  if (!kinTree->isEmpty() && tree->massConstraint()) {
    KinematicParticleFitter csFitter;
    auto_ptr<KinematicConstraint> con(new MassKinematicConstraint(tree->mass(), tree->massSigma()));
    kinTree = csFitter.fit(&(*con), kinTree);
  }

  tree->setKinematicTree(kinTree);

  kinTree->movePointerToTheTop();
  RefCountedKinematicParticle kinPart = kinTree->currentParticle();
  RefCountedKinematicVertex kinVertex = kinTree->currentDecayVertex();

  tree->fTV.setMaxDoca(getMaxDoca(kinParticles));
  tree->fTV.setMinDoca(getMinDoca(kinParticles));
  tree->fTV.setChi2(kinPart->chiSquared());

  return true;

}


// ----------------------------------------------------------------------
void HFSequentialVertexFit::addFittedParticles(vector<RefCountedKinematicParticle> *kinParticles, HFDecayTree *decayTree) {
  if (DBX) cout << "now in addFittedParticles for tree " << decayTree->particleID() << " at " << decayTree
		<< " &fTV.mass = " << &(decayTree->fTV.mass)
		<< endl;
  if (decayTree->vertexing()) {
    fitTree(decayTree); // fit again (no more dependencies on mothers)
    RefCountedKinematicTree kinTree = *(decayTree->getKinematicTree());
    if (kinTree->isEmpty()) throw EmptyTreeError();
    kinTree->movePointerToTheTop();
    kinParticles->push_back(kinTree->currentParticle());
    if (DBX) {
      RefCountedKinematicParticle p = kinTree->currentParticle();
      KinematicState state = p->currentState();
      TLorentzVector plab;
      plab.SetXYZM(state.globalMomentum().x(),state.globalMomentum().y(),state.globalMomentum().z(),state.mass());
      cout << "  vertexing::addFittedParticles(): Adding particle with\n\t x = " << plab.X() << ", y = " << plab.Y() << ", z = " << plab.Z()
	   << ", mass = " << plab.M() << endl;
    }
  } else {
    // recursively step down
    for (HFDecayTreeIterator treeIt = decayTree->getVerticesBeginIterator(); treeIt != decayTree->getVerticesEndIterator(); ++treeIt) {
      if (DBX) {
	cout << "adding particles for tree  " << treeIt->particleID() << " at " << &(*treeIt)
	     << " &fTV.mass = " << &(treeIt->fTV.mass)
	     << endl;
	RefCountedKinematicTree kinTree = *(treeIt->getKinematicTree());
	if (kinTree->isEmpty()) throw EmptyTreeError();
	kinTree->movePointerToTheTop();
	RefCountedKinematicParticle p = kinTree->currentParticle();
	KinematicState state = p->currentState();
	TLorentzVector plab;
	plab.SetXYZM(state.globalMomentum().x(),state.globalMomentum().y(),state.globalMomentum().z(),state.mass());
	cout << "  novertexing::addFittedParticles(): Adding particle with\n\t x = " << plab.X() << ", y = " << plab.Y() << ", z = " << plab.Z()
	     << ", mass = " << plab.M() << endl;
      }
      addFittedParticles(kinParticles, &(*treeIt));
    }
  }
}


// ----------------------------------------------------------------------
void HFSequentialVertexFit::saveTree(HFDecayTree *tree) {
  int dau1 = -1,dau2 = -1;
  TAnaCand *pCand,*pMomCand;
  RefCountedKinematicTree subTree;

  // -- create the Ana Candidate of the node if requested  (top candidate w.r.t. primary vertex)
  if (tree->particleID() && !tree->getAnaCand()) {
    TAnaCand *a = addCandidate(tree);
    if (a) tree->setAnaCand(a);
  }

  // get the current vertex state
  subTree = *(tree->getKinematicTree());
  subTree->movePointerToTheTop();
  VertexState vState = subTree->currentDecayVertex()->vertexState();

  // -- create all the requested candidates of the daughters
  HFDecayTreeIterator treeIt;
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
    if (treeIt->particleID() && !treeIt->getAnaCand())
      treeIt->setAnaCand(addCandidate(&(*treeIt), &vState));
  }

  // link the candidates
  pMomCand = tree->getAnaCand();
  if (pMomCand) {
    // now link the daughters
    for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
      pCand = treeIt->getAnaCand();
      if (pCand) {
    	// set mother
	pCand->fMom = pMomCand->fIndex;

     	if (dau1 == -1) dau1 = pCand->fIndex;
	else            dau1 = (pCand->fIndex < dau1) ? pCand->fIndex : dau1;

	if (dau2 == -1) dau2 = pCand->fIndex;
	else            dau2 = (pCand->fIndex > dau2) ? pCand->fIndex : dau2;
      }
    }

    pMomCand->fDau1 = dau1;
    pMomCand->fDau2 = dau2;
  }

  // -- override the dxy of the daughters if requested
  if (tree->daughtersToPV()) computeDaughterDistance(tree);

  // -- recursively continue
  for (treeIt = tree->getVerticesBeginIterator(); treeIt != tree->getVerticesEndIterator(); ++treeIt) {
    saveTree(&(*treeIt));
  }

  if (fVerbose > 5) {
    TAnaTrack *pTrk(0);
    pMomCand->dump();
    for (int i = pMomCand->fSig1; i <= pMomCand->fSig2; ++i) {
      pTrk = gHFEvent->getSigTrack(i);
      pTrk->dump();
    }
  }


}


// ----------------------------------------------------------------------
void HFSequentialVertexFit::computeDaughterDistance(HFDecayTree *tree) {
  VertexDistanceXY axy;
  VertexDistance3D a3d;
  RefCountedKinematicTree dauTree;
  RefCountedKinematicVertex dauVertex;

  // Load the Primary Vertex Collection
  TAnaCand *mom = tree->getAnaCand();
  for(HFDecayTreeIterator it = tree->getVerticesBeginIterator(); it != tree->getVerticesEndIterator(); ++it) {
    TAnaCand *dau = it->getAnaCand();
    if (!dau) continue;

    dauTree = *(it->getKinematicTree());
    dauTree->movePointerToTheTop();
    dauVertex = dauTree->currentDecayVertex();

    // Vertex Distanz neu berechnen zum PV: xy
    dau->fVtx.fDxy  = axy.distance((*fPVCollection)[mom->fPvIdx], dauVertex->vertexState() ).value();
    dau->fVtx.fDxyE = axy.distance((*fPVCollection)[mom->fPvIdx], dauVertex->vertexState() ).error();

    // Vertex Distanz neu berechnen zum PV: 3d
    dau->fVtx.fD3d  = a3d.distance((*fPVCollection)[mom->fPvIdx], dauVertex->vertexState() ).value();
    dau->fVtx.fD3dE = a3d.distance((*fPVCollection)[mom->fPvIdx], dauVertex->vertexState() ).error();
  }
}


// ----------------------------------------------------------------------
// Utility routine to sort the mindoca array of the candidate...
static bool doca_less(pair<int,pair<double,double> > x,pair<int,pair<double,double> > y) {
  return x.second.first < y.second.first;
}


// ----------------------------------------------------------------------
TAnaCand *HFSequentialVertexFit::addCandidate(HFDecayTree *tree, VertexState *wrtVertexState) {
  if (DBX) cout << "now in HFSequentialVertexFit::addCandidate for tree  " << tree->particleID() << " at " << tree << endl;
  unsigned int j;
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TransverseImpactPointExtrapolator transverseExtrapolator(fMagneticField);
  TrajectoryStateOnSurface tsos;

  if (tree->particleID() == 0) return 0;

  RefCountedKinematicTree kinTree = *(tree->getKinematicTree());
  if (kinTree->isEmpty()) return 0;
  kinTree->movePointerToTheTop();

  RefCountedKinematicParticle  kinParticle = kinTree->currentParticle();
  vector<RefCountedKinematicParticle> daughterParticles = kinTree->daughterParticles();

  RefCountedKinematicVertex  kinVertex = kinTree->currentDecayVertex();
  map<int,int> *kinParticleMap = tree->getKinParticleMap();

  if (!kinVertex->vertexIsValid()) return 0;

  TVector3 plab = TVector3(kinParticle->currentState().globalMomentum().x(),
			   kinParticle->currentState().globalMomentum().y(),
			   kinParticle->currentState().globalMomentum().z());
  double mass = kinParticle->currentState().mass();

  if (kinParticle->chiSquared() < 0) return 0;

  ChiSquared chi(kinParticle->chiSquared(),kinParticle->degreesOfFreedom());

  if (fVerbose > 2) {
    cout << "-----------------------------------------" << endl;
    cout << "==> HFSequentialVertexFit: Filling candidate " << tree->particleID() << " with mass = " << mass << endl;
    cout << "-----------------------------------------" << endl;
  }

  tree->fAnaVertex.setInfo(kinParticle->chiSquared(), kinParticle->degreesOfFreedom(), chi.probability(), 0, -1);
  tree->fAnaVertex.fPoint.SetXYZ(kinVertex->position().x(), kinVertex->position().y(), kinVertex->position().z());

  vector<track_entry_t> completeTrackList = tree->getAllTracks(0);
  vector<track_entry_t> allTreeTracks = tree->getAllTracks(1);

  if (DBX) {
    cout << "DBX addCand:  pid = " << tree->particleID() << " mass = " << mass << " flsxy = " << tree->fAnaVertex.fDxy/tree->fAnaVertex.fDxyE
	 << " flxy = " << tree->fAnaVertex.fDxy
	 << " pvips = " << tree->fTV.pvImpParams.ip3d.value()/tree->fTV.pvImpParams.ip3d.error()
	 << " pt/eta/phi = " << plab.Perp() << "/" << plab.Eta() << "/" << plab.Phi()
	 << " x/y/z = " << tree->fAnaVertex.fPoint.X() << "/" << tree->fAnaVertex.fPoint.Y() << "/" << tree->fAnaVertex.fPoint.Z()
	 << endl;
    cout << " tracks: " ;
    for (unsigned int i = 0; i < allTreeTracks.size(); ++i) cout << allTreeTracks[i].trackIx << " ";
    cout << endl;
    cout << "  going on to filling cand into tree" << endl;
  }


  if (TMath::Abs(mass - tree->fTV.mass) > 0.001*mass) {
    cout << "error in HFSequentialVertexFit: mass = " << mass << " and tree->fTV.mass = " << tree->fTV.mass << " are not equal" << endl;
  }

  // -- create and fill candidate
  TAnaCand *pCand = gHFEvent->addCand();
  pCand->fPlab  = plab;
  pCand->fMass  = tree->fTV.mass;
  pCand->fMassE = tree->fTV.masserr;
  pCand->fVtx = tree->fAnaVertex;
  pCand->fType = tree->particleID();

  pCand->fMom = -1; // Mom gets linked later.
  pCand->fDau1 = -1; // Daughters get linked later
  pCand->fDau2 = -1;

  pCand->fSig1 = gHFEvent->nSigTracks();
  pCand->fSig2 = pCand->fSig1 + allTreeTracks.size() - 1;

  pCand->fMaxDoca = tree->fTV.maxDoca;
  pCand->fMinDoca = tree->fTV.minDoca;

  pCand->fPvIdx = tree->fTV.pvIx;
  pCand->fPv2Idx = tree->fTV.pvIx2;
  pCand->fPvLip = tree->fTV.pvImpParams.lip.value();
  pCand->fPvLipE = tree->fTV.pvImpParams.lip.error();
  pCand->fPvTip = tree->fTV.pvImpParams.tip.value();
  pCand->fPvTipE = tree->fTV.pvImpParams.tip.error();
  pCand->fPv2Lip = tree->fTV.pvImpParams2nd.lip.value();
  pCand->fPv2LipE = tree->fTV.pvImpParams2nd.lip.error();
  pCand->fPv2Tip = tree->fTV.pvImpParams2nd.tip.value();
  pCand->fPv2TipE = tree->fTV.pvImpParams2nd.tip.error();

  pCand->fPvIP3d  = tree->fTV.pvImpParams.ip3d.value();
  pCand->fPvIP3dE = tree->fTV.pvImpParams.ip3d.error();

  pCand->fPv2IP3d  = tree->fTV.pvImpParams2nd.ip3d.value();
  pCand->fPv2IP3dE = tree->fTV.pvImpParams2nd.ip3d.error();

  pCand->fDeltaChi2 = tree->fTV.diffChi2;

  pCand->fAlpha = TMath::ACos(tree->fTV.vtxDistanceCosAlphaPlab);


  pCand->fTau3d  = tree->fTV.tau3d;
  pCand->fTau3dE = tree->fTV.tau3dE;

  pCand->fTauxy  = tree->fTV.tauxy;
  pCand->fTauxyE = tree->fTV.tauxyE;

  TAnaTrack *pTrack;
  for (j = 0; j < allTreeTracks.size(); j++) {
    TransientTrack fitTrack = daughterParticles[(*kinParticleMap)[allTreeTracks[j].trackIx]]->refittedTransientTrack();
    TSimpleTrack *sTrack = gHFEvent->getSimpleTrack(allTreeTracks[j].trackIx);
    pTrack = gHFEvent->addSigTrack();
    TrackBaseRef baseRef(fhTracks, allTreeTracks[j].trackIx);
    Track trackView(*baseRef);

    int gidx = sTrack->getGenIndex();
    const reco::BeamSpot *pBeamSpot = &fBeamSpot;
    fillAnaTrack(pTrack, trackView, allTreeTracks[j].trackIx, gidx, fPVCollection.product(), fMuons, pBeamSpot);

    pTrack->fIndex = allTreeTracks[j].trackIx;
    pTrack->fMCID = allTreeTracks[j].particleID; // use the sigTrack MCID to store the assumed particle ID for the mass hypothesis
    pTrack->fRefPlab = TVector3(fitTrack.track().px(),fitTrack.track().py(),fitTrack.track().pz());
    pTrack->fRefDof = fitTrack.ndof();
    pTrack->fRefValidHits = fitTrack.numberOfValidHits();
    pTrack->fRefChi2 = fitTrack.chi2();
    //    pTrack->fQ = fitTrack.charge();
    //    pTrack->fMuID = recTrack->fMuID;
    TAnaMuon *pM = gHFEvent->getSimpleTrackMuon(allTreeTracks[j].trackIx);
    if (pM) {
      pTrack->fMuID = pM->fMuID;
    } else {
      pTrack->fMuID = 0;
    }
    // get the reco muon if available
    if (fMuons) {
      for (MuonCollection::const_iterator muonIt = fMuons->begin(); muonIt != fMuons->end(); ++muonIt) {
	if ((int)muonIt->track().index() == allTreeTracks[j].trackIx) {
	  Vertex secVertex(RecoVertex::convertPos(kinVertex->vertexState().position()),
			   kinVertex->vertexState().error().matrix_new(),
			   kinVertex->chiSquared(),
			   kinVertex->degreesOfFreedom(),
			   daughterParticles.size());

	  if (tree->fTV.pvIx >= 0 && muon::isTightMuon(*muonIt, (*fPVCollection)[tree->fTV.pvIx]))
	    pTrack->fMuID |= 0x1<<15;

	  if (muon::isTightMuon(*muonIt, secVertex))
	    pTrack->fMuID |= 0x1<<16;

	  break;
	}
      }
    }
  }

  // fill the closest approaching tracks -- only if this is supposed to be a SV
  VertexDistance3D a3d;
  set<int> allUsedTrackIndices = tree->getAllTracksIndices();
  if (!wrtVertexState) {
    for (j = 0; j < fhTracks->size(); j++) {
      if (allUsedTrackIndices.count(j)>0) continue; // this tracks belongs to the candidate
      TrackBaseRef baseRef(fhTracks,j);
      TransientTrack transTrack = fpTTB->build(*baseRef);
      tsos = extrapolator.extrapolate(transTrack.initialFreeState(), kinVertex->position());
      if (!tsos.isValid()) {
	if (fVerbose > 0)  cout << "==>HFSequentialVertexFit> tsos not valid" << endl;
	continue;
      }
      Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), kinVertex->vertexState());
      pCand->fNstTracks.push_back(make_pair(j, make_pair(doca.value(), doca.error())));
    }

    sort(pCand->fNstTracks.begin(), pCand->fNstTracks.end(), doca_less);
    if (pCand->fNstTracks.size() > kNBR_CLOSED_TRACKS) {
      pCand->fNstTracks.erase(pCand->fNstTracks.begin() + kNBR_CLOSED_TRACKS, pCand->fNstTracks.end());
    }
  }

  return pCand;
}



// ----------------------------------------------------------------------
float HFSequentialVertexFit::getParticleMass(int particleID, float *mass_sigma) {
  float mass;
  float sigma = 0.0;
  particleID = abs(particleID);

  // sigma corresponds to standard uncertainty as can be found in the PDG
  switch(particleID) {
  case 11: // electron
    mass = MELECTRON;
    sigma = 0.013E-9f;
    break;
  case 13: // muon
    mass = MMUON;
    sigma = 4E-9f;
    break;
  case 211: // pion
    mass = MPION;
    sigma = 3.5E-7f;
    break;
  case 321: // kaon
    mass = MKAON;
    sigma = 1.6E-5f;
    break;
  case 2212: // proton
    mass = MPROTON;
    sigma = 8E-8f;
    break;
  default:
    throw MassNotFoundException();
    break;
  }

  if (mass_sigma) *mass_sigma = sigma;

  return mass;
}


// ----------------------------------------------------------------------
double HFSequentialVertexFit::getMaxDoca(vector<RefCountedKinematicParticle> &kinParticles) {
  double maxDoca = -1.0;
  TwoTrackMinimumDistance md;
  vector<RefCountedKinematicParticle>::iterator in_it, out_it;

  for (out_it = kinParticles.begin(); out_it != kinParticles.end(); ++out_it) {
    for (in_it = out_it + 1; in_it != kinParticles.end(); ++in_it) {
      md.calculate((*out_it)->currentState().freeTrajectoryState(),(*in_it)->currentState().freeTrajectoryState());
      if (md.distance() > maxDoca)
	maxDoca = md.distance();
    }
  }

  return maxDoca;
}


// ----------------------------------------------------------------------
double HFSequentialVertexFit::getMinDoca(vector<RefCountedKinematicParticle> &kinParticles) {
  double minDoca = 99999.9;
  TwoTrackMinimumDistance md;
  unsigned j,k,n;

  n = kinParticles.size();
  for (j = 0; j < n; j++) {
    for (k = j+1; k < n; k++) {
      md.calculate(kinParticles[j]->currentState().freeTrajectoryState(),kinParticles[k]->currentState().freeTrajectoryState());
      if (md.distance() < minDoca)
	minDoca = md.distance();
    }
  }

  return minDoca;
}


// ----------------------------------------------------------------------
cov33_t HFSequentialVertexFit::GlobalError2SMatrix_33(GlobalError m_in) {
  cov33_t m_out;
  for (int i=0; i<3; i++) {
    for (int j=i; j<3; j++)  {
      m_out(i,j) = m_in.matrix()(i,j);
    }
  }
  return m_out;
}


// ----------------------------------------------------------------------
cov99_t HFSequentialVertexFit::makeCovarianceMatrix(const cov33_t cov_vtx1,
						    const cov77_t cov_vtx2) {
  cov99_t cov;
  cov.Place_at(cov_vtx1,0,0);
  cov.Place_at(cov_vtx2.Sub<cov66_t>(0,0),3,3);
  return cov;
}


// ----------------------------------------------------------------------
jac9_t HFSequentialVertexFit::makeJacobianVector3d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2, const AlgebraicVector3 &momentum) {
  jac9_t jac;
  const AlgebraicVector3 dist = vtx2 - vtx1;
  const double factor2 = 1. / ROOT::Math::Mag2(momentum);
  const double lifetime = ROOT::Math::Dot(dist, momentum) * factor2;
  jac.Place_at(-momentum*factor2,0);
  jac.Place_at( momentum*factor2,3);
  jac.Place_at( factor2*(dist-2*lifetime*momentum*factor2),6);
  return jac;
}


// ----------------------------------------------------------------------
jac9_t HFSequentialVertexFit::makeJacobianVector3d(const GlobalPoint &vtx1, const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
  // TODO: Update 2d calculation to projected version as in 3d
  return makeJacobianVector3d(AlgebraicVector3(vtx1.x(),vtx1.y(),vtx1.z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}


// ----------------------------------------------------------------------
jac9_t HFSequentialVertexFit::makeJacobianVector3d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
									  ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
									  const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
  return makeJacobianVector3d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}


// -------------------------------------------------------------
jac9_t HFSequentialVertexFit::makeJacobianVector2d(const GlobalPoint &vtx1, const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
  return makeJacobianVector2d(AlgebraicVector3(vtx1.x(),vtx1.y(),vtx1.z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}


// ----------------------------------------------------------------------
jac9_t HFSequentialVertexFit::makeJacobianVector2d(const ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,
									  ROOT::Math::DefaultCoordinateSystemTag> &vtx1,
									  const GlobalPoint &vtx2, const TVector3 &tv3momentum) {
  return makeJacobianVector2d(AlgebraicVector3(vtx1.X(),vtx1.Y(),vtx1.Z()),
			      AlgebraicVector3(vtx2.x(),vtx2.y(),vtx2.z()),
			      AlgebraicVector3(tv3momentum.x(),tv3momentum.y(),tv3momentum.z()));
}


// ----------------------------------------------------------------------
jac9_t HFSequentialVertexFit::makeJacobianVector2d(const AlgebraicVector3 &vtx1, const AlgebraicVector3 &vtx2,
									  const AlgebraicVector3 &momentum) {
  jac9_t jac;
  const double momentumMag = ROOT::Math::Mag(momentum);
  const AlgebraicVector3 dist = vtx2 - vtx1;
  const double distMag = ROOT::Math::Mag(dist);
  const double factorPositionComponent = 1./(distMag*momentumMag);
  const double factorMomentumComponent = 1./pow(momentumMag,3);
  jac(0)=-dist(0)*factorPositionComponent;
  jac(1)=-dist(1)*factorPositionComponent;
  jac(3)= dist(0)*factorPositionComponent;
  jac(4)= dist(1)*factorPositionComponent;
  jac(6)= momentum(0)*factorMomentumComponent;
  jac(7)= momentum(1)*factorMomentumComponent;
  return jac;
}


// ----------------------------------------------------------------------
void HFSequentialVertexFit::calculateStuff(HFDecayTree *tree, VertexState *wrtVertexState) {

  if (tree->fTV.valid) {
    if (DBX) cout << "HFSequentialVertexFit::calculateStuff> already calculated stuff for tree " << tree->particleID()
		  << " at " << tree << " with VertexState at " << wrtVertexState
		  << endl;
    return;
  }

  RefCountedKinematicTree kinTree = *(tree->getKinematicTree());
  kinTree->movePointerToTheTop();
  RefCountedKinematicParticle kinParticle = kinTree->currentParticle();
  RefCountedKinematicVertex   kinVertex  = kinTree->currentDecayVertex();

  TVector3 plab = TVector3(kinParticle->currentState().globalMomentum().x(),
			   kinParticle->currentState().globalMomentum().y(),
			   kinParticle->currentState().globalMomentum().z());

  if (DBX) {
    cout << "now in HFSequentialVertexFit::calculateStuff for tree " << tree->particleID() << " at " << tree
	 << " tree->fTV = " << &(tree->fTV)
	 << " ("  << kinVertex->position().x() << "/" << kinVertex->position().y() << "/" << kinVertex->position().z() << ")"
	 << " and vertexstate " << wrtVertexState;
    if (wrtVertexState) {
      cout << " (" << wrtVertexState->position().x() << "/" << wrtVertexState->position().y() << "/" << wrtVertexState->position().z() << ")";
    } else {
      cout << " (0)";
    }
    cout << endl;
  }

  tree->fTV.setPt(plab.Perp());
  tree->fTV.setMass(kinParticle->currentState().mass());

  tree->fTV.setMassE(TMath::Sqrt(kinParticle->currentState().kinematicParametersError().matrix()(6,6)));   // mass error from Keith

  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  TransverseImpactPointExtrapolator transverseExtrapolator(fMagneticField);

  VertexDistanceXY axy;
  VertexDistance3D a3d;

  vector<track_entry_t> completeTrackList = tree->getAllTracks(0);
  if (!wrtVertexState) {
    if (fVerbose > 2)  cout << "==> HFSequentialVertexFit: Number of PV vertices to compare is " << fPVCollection->size() << endl;
    unsigned int nGoodVtx;
    unsigned int j;
    VertexCollection::const_iterator vertexIt;
    for (vertexIt = fPVCollection->begin(), j = 0, nGoodVtx = 0; vertexIt != fPVCollection->end(); ++vertexIt,++j) {
      std::pair<bool,Measurement1D> currentIp;
      std::pair<bool,Measurement1D> cur3DIP;

      AdaptiveVertexFitter avf;
      vector<TransientTrack> vrtxRefit;

      for (vector<TrackBaseRef>::const_iterator itTBR = vertexIt->tracks_begin(); itTBR != vertexIt->tracks_end(); ++itTBR) {
	TrackRef tref = itTBR->castTo<TrackRef>();
	vector<track_entry_t>::const_iterator trackIt;
	for (trackIt = completeTrackList.begin(); fRemoveCandTracksFromVtx && trackIt != completeTrackList.end(); ++trackIt) {
	  TrackBaseRef curTr(fhTracks, trackIt->trackIx);
	  TrackRef curTref = curTr.castTo<TrackRef>();
	  if (tref == curTref)
	    break;
	}

	if (!fRemoveCandTracksFromVtx || trackIt == completeTrackList.end()) {
	  TransientTrack tTrk = fpTTB->build(*(*itTBR));
	  vrtxRefit.push_back(tTrk);
	}
      }

      if (vrtxRefit.size() < 2) continue; // next one
      TransientVertex newVtx = avf.vertex(vrtxRefit, fBeamSpot);
      if (!newVtx.isValid()) continue; // no valid refit, take next one
      Vertex currentPV = reco::Vertex(newVtx);

      // extrapolate to PCA
      TrajectoryStateOnSurface tsos = extrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),
							       RecoVertex::convertPos(currentPV.position()));

      // compute with iptools
      currentIp = IPTools::signedDecayLength3D(tsos, GlobalVector(0,0,1), currentPV);
      cur3DIP = IPTools::absoluteImpactParameter(tsos, currentPV, a3d);

      // Compute the PV weight
      double weight = (currentPV.ndof()+2.)/(2. * currentPV.tracksSize());

      if( weight < fPvW8) { // Check the PV weight and skip it if lower than the cut
	if (fVerbose > 2) cout<<"==>HFSequentialVertexFit: PV "<<j<<" rejected because of too low weight "<<weight<<endl;
	continue;
      }

      if (!currentIp.first) {
	if (fVerbose > 2) cout << "==>HFSequentialVertexFit: Unable to compute lip to vertex at index " << j << endl;
	continue;
      }

      // store?
      // the first PV, this is currently the best one ;)
      if (nGoodVtx == 0 ) {
	tree->fTV.pvIx = j;
	tree->fTV.pvImpParams.lip = currentIp.second;
	tree->fTV.pvImpParams.ip3d = cur3DIP.second;
	tree->fTV.setPvips(tree->fTV.pvImpParams.ip3d.value()/tree->fTV.pvImpParams.ip3d.error());
      } else if (nGoodVtx == 1) { // now the second PV
	if (fabs(currentIp.second.value()) >= fabs(tree->fTV.pvImpParams.lip.value())) {// not the best but the second best
	  tree->fTV.pvImpParams2nd.lip = currentIp.second;
	  tree->fTV.pvImpParams2nd.ip3d = cur3DIP.second;
	  tree->fTV.pvIx2 = j;
	}
	else {   // the best, the previous one is the current 2nd best
	  tree->fTV.pvIx2 = tree->fTV.pvIx;
	  tree->fTV.pvIx = j;
	  tree->fTV.pvImpParams2nd.lip = tree->fTV.pvImpParams.lip;
	  tree->fTV.pvImpParams2nd.ip3d = tree->fTV.pvImpParams.ip3d;
	  tree->fTV.pvImpParams.lip = currentIp.second;
	  tree->fTV.pvImpParams.ip3d = cur3DIP.second;
	  tree->fTV.setPvips(tree->fTV.pvImpParams.ip3d.value()/tree->fTV.pvImpParams.ip3d.error());
	}
      } else { // we have more than 2 PV
	if (fabs(currentIp.second.value()) >= fabs(tree->fTV.pvImpParams.lip.value())) {
	  // not the best
	  if (fabs(currentIp.second.value()) < fabs(tree->fTV.pvImpParams2nd.lip.value())) {
	    // but the second best
	    tree->fTV.pvImpParams2nd.lip = currentIp.second;
	    tree->fTV.pvImpParams2nd.ip3d = cur3DIP.second;
	    tree->fTV.pvIx2 = j;
	  }
	} else {
	  // this is currently the best one, keep it and put the old best one to 2nd best
	  tree->fTV.pvIx2 = tree->fTV.pvIx;
	  tree->fTV.pvIx = j;
	  tree->fTV.pvImpParams2nd.lip = tree->fTV.pvImpParams.lip;
	  tree->fTV.pvImpParams2nd.ip3d = tree->fTV.pvImpParams.ip3d;
	  tree->fTV.pvImpParams.lip = currentIp.second;
	  tree->fTV.pvImpParams.ip3d = cur3DIP.second;
	  tree->fTV.setPvips(tree->fTV.pvImpParams.ip3d.value()/tree->fTV.pvImpParams.ip3d.error());
	}
      }
      nGoodVtx++; // Count the no. of good vertices
    }
  }

  if (wrtVertexState) {
    // -- Distance to mother vertex
    tree->fAnaVertex.fDxy = axy.distance(*wrtVertexState, kinVertex->vertexState()).value();
    tree->fAnaVertex.fDxyE = axy.distance(*wrtVertexState, kinVertex->vertexState()).error();
    tree->fTV.setFlxy(tree->fAnaVertex.fDxy);
    tree->fTV.setFlsxy(tree->fAnaVertex.fDxy/tree->fAnaVertex.fDxyE);

    tree->fAnaVertex.fD3d = a3d.distance(*wrtVertexState, kinVertex->vertexState()).value();
    tree->fAnaVertex.fD3dE = a3d.distance(*wrtVertexState, kinVertex->vertexState()).error();
    // -- get covariance matrix for error propagation in lifetime calculation
    tree->fTV.vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(wrtVertexState->error()),
						kinParticle->currentState().kinematicParametersError().matrix());
    tree->fTV.vtxDistanceJac3d = makeJacobianVector3d(wrtVertexState->position(), kinVertex->vertexState().position(), plab);
    tree->fTV.vtxDistanceJac2d = makeJacobianVector2d(wrtVertexState->position(), kinVertex->vertexState().position(), plab);
    // -- get sign of distance
    const GlobalVector diff = kinVertex->vertexState().position() - wrtVertexState->position() ;
    const TVector3 tv3diff = TVector3(diff.x(),diff.y(),diff.z());
    tree->fTV.vtxDistanceCosAlphaPlab = plab.Dot(tv3diff) / (plab.Mag() * tv3diff.Mag());

    TrajectoryStateOnSurface tsos = transverseExtrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),
								       wrtVertexState->position());
    if (!tsos.isValid()) {
      if (fVerbose > 0)  cout << "==>HFSequentialVertexFit> tsos not valid" << endl;
      tree->fTV.pvImpParams.tip = (Measurement1D(-9999.,-9999.), Measurement1D(-9999.,-9999.), Measurement1D(-9999.,-9999.));
    } else {
      tree->fTV.pvImpParams.tip = axy.distance(VertexState(tsos.globalPosition(),tsos.cartesianError().position()),
					       *wrtVertexState);
    }
  } else if (tree->fTV.pvIx >= 0) {
    // -- Distance w.r.t primary vertex
    AdaptiveVertexFitter avf; // PV refit
    KalmanVertexFitter kvf; // delta chi2
    Vertex currentPV = (*fPVCollection)[tree->fTV.pvIx];
    Vertex currentPVWithSignal;
    vector<TransientTrack> vrtxRefit;

    if (fVerbose > 5) cout << "HFSequentialVertexFit::addCandidate(): currentPV.tracksSize() = " << currentPV.tracksSize() << endl;
    vector<track_entry_t>::const_iterator trackIt;
    for (vector<TrackBaseRef>::const_iterator itTBR = currentPV.tracks_begin(); itTBR != currentPV.tracks_end(); ++itTBR) {
      TrackRef tref = itTBR->castTo<TrackRef>();
      for (trackIt = completeTrackList.begin(); fRemoveCandTracksFromVtx && trackIt != completeTrackList.end(); ++trackIt) {
	TrackBaseRef curTr(fhTracks, trackIt->trackIx);
	TrackRef curTref = curTr.castTo<TrackRef>();
	if (tref == curTref)
	  break;
      }

      if (!fRemoveCandTracksFromVtx || trackIt == completeTrackList.end()) {
	TransientTrack tTrk = fpTTB->build(*(*itTBR));
	vrtxRefit.push_back(tTrk);
      }
    }
    if (vrtxRefit.size() < 5) throw PVRefitException(); // do not try to fit with less than five tracks, else problems in CrossingPtBasedLinPtFinder
    if (fVerbose > 5) cout << "==> HFSequentialVertexFit::addCandidate(): refitting with vrtxRefit.size() = " << vrtxRefit.size() << endl;

    TransientVertex newVtx = avf.vertex(vrtxRefit, fBeamSpot);
    if (newVtx.isValid()) {
      currentPV = reco::Vertex(newVtx);
    } else {
      throw PVRefitException();
    }

    tree->fAnaVertex.fDxy = axy.distance(currentPV,  kinVertex->vertexState()).value();
    tree->fAnaVertex.fDxyE = axy.distance(currentPV, kinVertex->vertexState()).error();
    tree->fTV.setFlxy(tree->fAnaVertex.fDxy);
    tree->fTV.setFlsxy(tree->fAnaVertex.fDxy/tree->fAnaVertex.fDxyE);

    tree->fAnaVertex.fD3d = a3d.distance(currentPV,  kinVertex->vertexState()).value();
    tree->fAnaVertex.fD3dE = a3d.distance(currentPV, kinVertex->vertexState()).error();
    tree->fTV.setFls3d(tree->fAnaVertex.fD3d/tree->fAnaVertex.fD3dE);

    // -- recompute impact parameter: lip & ip3d
    TrajectoryStateOnSurface tsos = extrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(),
							     RecoVertex::convertPos(currentPV.position()));
    tree->fTV.pvImpParams.lip = IPTools::signedDecayLength3D(tsos, GlobalVector(0,0,1), currentPV).second;
    tree->fTV.pvImpParams.ip3d = IPTools::absoluteImpactParameter(tsos, currentPV, a3d).second;
    tree->fTV.setPvips(tree->fTV.pvImpParams.ip3d.value()/tree->fTV.pvImpParams.ip3d.error());

    // -- recompute impact parameter: tip
    tsos = transverseExtrapolator.extrapolate(kinParticle->currentState().freeTrajectoryState(), RecoVertex::convertPos(currentPV.position()));
    if (!tsos.isValid()) {
      if (fVerbose > 0)  cout << "==>HFSequentialVertexFit> tsos not valid" << endl;
      tree->fTV.pvImpParams.tip = (Measurement1D(-9999.,-9999.), Measurement1D(-9999.,-9999.), Measurement1D(-9999.,-9999.));
    } else {
      tree->fTV.pvImpParams.tip = axy.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()),
					       VertexState(RecoVertex::convertPos(currentPV.position()),
							   RecoVertex::convertError(currentPV.error())));
    }

    // -- get covariance matrix for error propagation in lifetime calculation
    tree->fTV.vtxDistanceCov = makeCovarianceMatrix(GlobalError2SMatrix_33(currentPV.error()),
						    kinParticle->currentState().kinematicParametersError().matrix());
    tree->fTV.vtxDistanceJac3d = makeJacobianVector3d(currentPV.position(), kinVertex->vertexState().position(), plab);
    tree->fTV.vtxDistanceJac2d = makeJacobianVector2d(currentPV.position(), kinVertex->vertexState().position(), plab);
    // -- get sign of distance
    const TVector3 p1(currentPV.position().x(), currentPV.position().y(), currentPV.position().z());
    const TVector3 p2(kinVertex->vertexState().position().x(), kinVertex->vertexState().position().y(), kinVertex->vertexState().position().z());
    const TVector3 pDiff = p2-p1;
    tree->fTV.vtxDistanceCosAlphaPlab = plab.Dot(pDiff) / (plab.Mag() * pDiff.Mag());

    // compute the delta chi2 using the Kalman vertex fitter and only tracks with weight > 0.5
    vrtxRefit.clear();
    currentPV = (*fPVCollection)[tree->fTV.pvIx];
    for (vector<TrackBaseRef>::const_iterator itTBR = currentPV.tracks_begin(); itTBR != currentPV.tracks_end(); ++itTBR) {
      if (currentPV.trackWeight(*itTBR) <= 0.5)
	continue;

      TrackRef tref = itTBR->castTo<TrackRef>();
      vector<track_entry_t>::const_iterator trackIt;
      for (trackIt = completeTrackList.begin(); trackIt != completeTrackList.end(); ++trackIt) {
	TrackBaseRef curTr(fhTracks, trackIt->trackIx);
	TrackRef curTref = curTr.castTo<TrackRef>();
	if (tref == curTref) {
	  break;
	}
      }
      if (trackIt == completeTrackList.end()) {
	TransientTrack tTrk = fpTTB->build( *(*itTBR) );
	vrtxRefit.push_back(tTrk);
      }
    }
    newVtx = kvf.vertex(vrtxRefit);
    if (newVtx.isValid()) {
      currentPV = reco::Vertex(newVtx);
    } else	{
      throw PVRefitException();
    }
    for (vector<track_entry_t>::const_iterator trackIt = completeTrackList.begin(); trackIt != completeTrackList.end(); ++trackIt) {
      TrackBaseRef curTr(fhTracks, trackIt->trackIx);
      vrtxRefit.push_back(fpTTB->build(*curTr));
    }
    newVtx = kvf.vertex(vrtxRefit);
    if (newVtx.isValid()) {
      currentPVWithSignal = reco::Vertex(newVtx);
    } else {
      throw PVRefitException();
    }

    tree->fTV.diffChi2 = currentPVWithSignal.chi2() - currentPV.chi2();
  } else if (fVerbose > 2) {
    cout << "==> HFSequentialVertexFit: No idea what distance to compute in TAnaVertex.fDxy and TAnaVertex.fD3d" << endl;
  }

  // -- set covariance matrix
  double cov[9];
  cov[0] = kinVertex->error().cxx();
  cov[1] = kinVertex->error().cyx();
  cov[2] = kinVertex->error().czx();
  cov[3] = kinVertex->error().cyx();
  cov[4] = kinVertex->error().cyy();
  cov[5] = kinVertex->error().czy();
  cov[6] = kinVertex->error().czx();
  cov[7] = kinVertex->error().czy();
  cov[8] = kinVertex->error().czz();
  tree->fAnaVertex.setCovXX(cov);


  // -- calculate lifetime
  // TMath::Ccgs() is to convert from cm to s (speed of light in cgs system, CMS uses cm)
  if (plab.Mag() > 0) {
    const double massOverC = tree->fTV.mass/TMath::Ccgs();
    // from 3d vertexing
    tree->fTV.tau3d = tree->fAnaVertex.fD3d / plab.Mag() * tree->fTV.vtxDistanceCosAlphaPlab * massOverC;
    tree->fTV.tau3dE = TMath::Sqrt(ROOT::Math::Similarity(tree->fTV.vtxDistanceCov, tree->fTV.vtxDistanceJac3d)) * massOverC;
    // from 2d vertexing
    const double sinTheta = TMath::Sin(plab.Theta());
    const double flightlength2d = (sinTheta != 0 ? tree->fAnaVertex.fDxy / sinTheta : 0);
    tree->fTV.tauxy = flightlength2d / plab.Mag() * tree->fTV.vtxDistanceCosAlphaPlab * massOverC;
    tree->fTV.tauxyE = TMath::Sqrt(ROOT::Math::Similarity(tree->fTV.vtxDistanceCov, tree->fTV.vtxDistanceJac2d)) * massOverC;
  } else {
    tree->fTV.tau3d = tree->fTV.tauxy = -99.;
  }

  tree->fTV.valid = true;
  //  cout << "set  tree->fTV.valid " <<  (tree->fTV.valid?"true":"false") << endl;
}
