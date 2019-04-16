#include "Bmm/CmsswAnalysis/plugins/HfrGenRecoilProducer.h"
#include "Bmm/RootAnalysis/common/HFMasses.hh"

#include <memory>
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "Bmm/CmsswAnalysis/interface/HFDumpUtilities.hh"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"
#include "Bmm/RootAnalysis/common/util.hh"

#include <TLorentzVector.h>

using namespace std;
using namespace edm;
using namespace reco;


// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HfrGenRecoilProducer::HfrGenRecoilProducer(const ParameterSet& iConfig) :
  HfrBaseProducer(iConfig),
  fVtxProb(iConfig.getUntrackedParameter<double>("vtxProb", 0.01)),
  fMotherID(iConfig.getUntrackedParameter("motherID", 0)),
  fPartialDecayMatching(iConfig.getUntrackedParameter<bool>("partialDecayMatching", false)),
  fStableDaughters(0) {

  vector<int> defaultIDs;
  defaultIDs.push_back(0);
  fDaughtersID = iConfig.getUntrackedParameter<vector<int> >("daughtersID", defaultIDs);

  dumpConfiguration();
  produces<vector<int> >(); // vector of track indices
  produces<int>();          // possibly a PV index
}


// ----------------------------------------------------------------------
void HfrGenRecoilProducer::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HfrGenRecoilProducer::dumpConfiguration()" << endl;
  HfrBaseProducer::dumpConfiguration();
  fDaughtersSet.clear();
  cout << "--- vtxProb:               " << fVtxProb << endl;
  cout << "--- motherID:              " << fMotherID << endl;
  cout << "--- partialDecayMatching:  " << fPartialDecayMatching << endl;
  for (unsigned int i = 0; i < fDaughtersID.size(); ++i) {
    cout << "---   daughterID:              " << fDaughtersID[i] << endl;
    if (TMath::Abs(fDaughtersID[i]) == 11)   ++fStableDaughters;
    if (TMath::Abs(fDaughtersID[i]) == 13)   ++fStableDaughters;
    if (TMath::Abs(fDaughtersID[i]) == 211)  ++fStableDaughters;
    if (TMath::Abs(fDaughtersID[i]) == 321)  ++fStableDaughters;
    if (TMath::Abs(fDaughtersID[i]) == 2212) ++fStableDaughters;
    fDaughtersSet.insert(TMath::Abs(fDaughtersID[i]));
  }
  cout << "---    total stable particles: " << fStableDaughters << endl;
  cout << "----------------------------------------------------------------------" << endl;
}

// ----------------------------------------------------------------------
void HfrGenRecoilProducer::put(edm::Event &iEvent,
				  std::unique_ptr<std::vector<int> > &trkIdxColl,
				  std::unique_ptr<int> &pvIdx) {
  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
}

// ----------------------------------------------------------------------
void HfrGenRecoilProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  auto trkIdxColl = std::make_unique<vector<int> >();
  auto pvIdx = std::make_unique<int>(1234);

  if (fVerbose > 0)  cout << " ======================================================================"
			  << endl
			  << "=== HfrGenRecoilProducer run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event() << endl
			  << " ----------------------------------------------------------------------"
			  << endl;

  HfrBaseProducer::analyze(iEvent, iSetup);

  // -- (1) In generator block, find mother with declared decay channel
  multiset<int> genDaughters;
  multiset<int> genIndices;
  multiset<pair<int, int> > genMap;
  TGenCand *pGen(0), *pDau(0), *pTmp(0), *pBsignal(0), *pBtag(0);
  int matchedDecay(0);
  int iMom(-1);

  vector<int> bla(100);
  vector<int>::iterator blaIt;

  if (fVerbose > 0) cout << " ngenCands: " << gHFEvent->nGenCands() << endl;
  for (int ig = 0; ig < gHFEvent->nGenCands(); ++ig) {
    pGen = gHFEvent->getGenCand(ig);
    if (TMath::Abs(pGen->fID) == fMotherID) {
      pBsignal = pGen;
      if (fVerbose > 1) {
	cout << "mother ";
	pGen->dump();
      }
      genDaughters.clear();
      genIndices.clear();
      genMap.clear();

      // -- version with descendants
      for (int id = ig+1; id < gHFEvent->nGenCands(); ++id) {
	pDau = gHFEvent->getGenCand(id);
	iMom = pDau->fMom1;
	while (iMom > ig) {
	  pTmp = gHFEvent->getGenCand(iMom);
	  iMom = pTmp->fMom1;
	}
	if (iMom == ig) {
	  if (fVerbose > 1) {
	    cout << "  daug: ";
	    pDau->dump();
	  }
	  if (fPartialDecayMatching) {
	    if (fDaughtersSet.find(TMath::Abs(pDau->fID)) != fDaughtersSet.end()) {
	      genDaughters.insert(TMath::Abs(pDau->fID));
	      genIndices.insert(id);
	      genMap.insert(make_pair(id, TMath::Abs(pDau->fID)));
	    }
	  } else {
	    genDaughters.insert(TMath::Abs(pDau->fID));
	    genIndices.insert(id);
	    genMap.insert(make_pair(id, TMath::Abs(pDau->fID)));
	  }
	}
      }

      // -- now check whether this is PARTIALLY the decay channel in question
      if (fPartialDecayMatching) {
	blaIt = set_intersection(genDaughters.begin(), genDaughters.end(), fDaughtersSet.begin(), fDaughtersSet.end(), bla.begin());
	if (static_cast<unsigned int>(blaIt - bla.begin()) == fDaughtersSet.size()) {
	  matchedDecay = 1;
	  if (fVerbose > 0) {
	    cout << "matched partial decay: ";
	    for (vector<int>::iterator it = bla.begin(); it != blaIt; ++it) cout << *it << " ";
	    cout << endl;
	  }
	  break;
	}
      }

      // -- allow for arbitrary number of photons
      multiset<int> genDaughters2 = genDaughters;
      bool missingParticles(false);
      for (multiset<int>::iterator it = fDaughtersSet.begin(); it != fDaughtersSet.end(); ++it) {
	std::multiset<int>::iterator hit(genDaughters2.find(*it));
	if (hit!= genDaughters2.end()) {
	  genDaughters2.erase(hit);
	} else {
	  missingParticles = true;
	}
      }
      if (!missingParticles) {
	bool onlyPhotonsLeft(true);
	int nphotons(0);
	for (multiset<int>::iterator it = genDaughters2.begin(); it != genDaughters2.end(); ++it) {
	  if (*it != 22) {
	    onlyPhotonsLeft = false;
	  } else {
	    ++nphotons;
	  }
	}
	if (onlyPhotonsLeft) {
	  matchedDecay = 1;
	  if (fVerbose > 0) cout << "matched decay with " << nphotons << " bremsstrahlung photons" << endl;
	  break;
	}
      }
    }
  }


  // -- If nothing matched, put 'nothing' into event
  if (0 == matchedDecay) {
    put(iEvent, trkIdxColl, pvIdx);
    return;
  }

  // -- (2) find other gen-level B with 'same' PV as pBsignal
  int hq(quarkFlavor(pBsignal->fID));
  cout << "quark flavor: " << hq << endl;
  for (int ig = 0; ig < gHFEvent->nGenCands(); ++ig) {
    pGen = gHFEvent->getGenCand(ig);
    if (pGen == pBsignal) continue;
    if (quarkFlavor(pGen->fID) != hq) continue;
    if (5 == hq && !(isBeautyMesonWeak(pGen->fID) || isBeautyBaryonWeak(pGen->fID))) continue;
    if (4 == hq && !(isCharmMesonWeak(pGen->fID) || isCharmBaryonWeak(pGen->fID))) continue;
    if ((pGen->fV - pBsignal->fV).Mag() < 1.e-3) {
      pBtag = pGen;
      if (0) {
	pBsignal->dump();
	pGen->dump();
      }
      for (int im = pBtag->fMom1; im <= pBtag->fMom2; ++im) {
	if (pBsignal->fNumber == im) {
	  cout << "FIIIIIIIIXXXXXXXMMMMMMMMEEEEEEEE! pBsignal->fNumber == im" << endl;
	}
      }
      for (int im = pBsignal->fMom1; im <= pBsignal->fMom2; ++im) {
	if (pBtag->fNumber == im) {
	  cout << "FIIIIIIIIXXXXXXXMMMMMMMMEEEEEEEE! pBtag->fNumber == im " << endl;
	}
      }

      break;
    }
  }

  if (0 == pBtag) {
    cout << "HfrGenRecoilProducer: ERRRRRRROOOOORRRRR no Btag at gen-level found!!!!!!!!" << endl;
    put(iEvent, trkIdxColl, pvIdx);
    return;
  }

  // -- (3) create Breco (Btag) index with all simple tracks gen-matched to descendents of Btag
  vector<int> brecoIdx;
  TSimpleTrack *ps(0);
  //  pBtag->dump();
  for (int ig = 0; ig < gHFEvent->nGenCands(); ++ig) {
    pGen = gHFEvent->getGenCand(ig);
    if (gHFEvent->isAncestor(pBtag, pGen)) {
      //      pGen->dump();
      for (int it = 0; it < gHFEvent->nSimpleTracks(); ++it) {
	ps = gHFEvent->getSimpleTrack(it);
	if (ps->getGenIndex() == pGen->fNumber) {
	  //	  cout << "matched to track idx = " << it << endl;
	  brecoIdx.push_back(it);
	}
      }
    }
  }

  if (fVerbose > 0) {
    cout << " ======= brecoIdx.size() = " << brecoIdx.size()  << endl;
  }

  // -- find reco PV
  TVector3 pvGen(pBtag->fV);
  double bestDist(999.);
  int bestPV(999);
  for (int ipv = 0; ipv < gHFEvent->nPV(); ++ipv) {
    TAnaVertex *pV = gHFEvent->getPV(ipv);
    double  dist = (pV->fPoint - pvGen).Mag();
    if (dist < bestDist) {
      bestDist = dist;
      bestPV   = ipv;
    }
  }

  // -- get candidate and PV from gHFEvent
  pvIdx = std::make_unique<int>(bestPV);

  Vertex myVertex = fVertexCollection[bestPV];
  GlobalPoint pVertex(myVertex.position().x(), myVertex.position().y(), myVertex.position().z());
  TVector3    p3Vertex(myVertex.position().x(), myVertex.position().y(), myVertex.position().z());

  if (fVerbose > 0) {
    cout << " ======= genPV = " << formatTVector3(pvGen)
	 << " reco PV = " << formatTVector3(gHFEvent->getPV(bestPV)->fPoint)
	 << " cmssw PV = " <<formatTVector3(p3Vertex)
	 << endl;
  }

  // -- create and fill the track collections
  TrajectoryStateOnSurface tsos;
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  VertexDistance3D a3d;
  int ntracks(0), cnt(0);
  for (unsigned int ix = 0; ix < fTracksHandle->size(); ix++) {
    // -- skip tracks from BRECO
    for (unsigned int ib = 0; ib < brecoIdx.size(); ++ib) {
      if (static_cast<int>(ix) == brecoIdx[ib]) continue;
    }
    reco::TrackBaseRef rTrackView(fTracksHandle, ix);
    const reco::Track track(*rTrackView);
    TransientTrack transTrack = fTTB->build(*rTrackView);
    ++ntracks;
    if (!track.quality(reco::TrackBase::qualityByName(fTrackQualityString))) continue;
    if (track.pt() < fTrackMinPt) continue;

    tsos = extrapolator.extrapolate(transTrack.initialFreeState(), pVertex);
    Measurement1D doca = a3d.distance(VertexState(tsos.globalPosition(), tsos.cartesianError().position()), myVertex);
    if (doca.value() > 0.1) continue;

    trkIdxColl->push_back(ix);
    ++cnt;
  }
  if (fVerbose > 0) {
    cout << "==>HfrGenRecoilProducer> put into event tracklist with " << trkIdxColl->size()
	 << " tracks, removed "  << ntracks - trkIdxColl->size()
	 << endl;
  }
  // -- put into event
  put(iEvent, trkIdxColl, pvIdx);
}

// ----------------------------------------------------------------------
void HfrGenRecoilProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HfrGenRecoilProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}

DEFINE_FWK_MODULE(HfrGenRecoilProducer);
