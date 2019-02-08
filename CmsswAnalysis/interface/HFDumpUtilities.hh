#ifndef HFDUMPUTILITIES_H
#define HFDUMPUTILITIES_H

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "Bmm/RootAnalysis/rootio/TAnaTrack.hh"
#include "Bmm/RootAnalysis/rootio/TSimpleTrack.hh"


void fillSimpleTrack(TSimpleTrack *pTrack, const reco::Track &trackView,
		     int tidx, int mid, int gidx, const reco::VertexCollection *vc);

void fillAnaTrack(TAnaTrack *pTrack, const reco::Track &trackView, int tidx, int gidx,
		  const reco::VertexCollection *vc, const reco::MuonCollection *mc, const reco::BeamSpot *bs);

TAnaTrack* fillSigTrack(int tidx,
			edm::Handle<edm::View<reco::Track> > &,
			const reco::VertexCollection *vc, const reco::MuonCollection *mc, const reco::BeamSpot *bs);


int getPv(int tidx, const reco::VertexCollection *vc);

int muonID(const reco::Muon &rm);

//void cleanupTruthMatching();
void cleanupTruthMatching(edm::Handle<edm::View<reco::Track> > &hTracks, const MagneticField *magfield);

std::pair<float, float> getParticleMassAndSigma(int particleID);
RefCountedKinematicTree fitTree(std::vector<reco::TransientTrack> &, std::vector<int> &);
std::pair<int, double> findBestPV(std::vector<reco::TransientTrack> &, RefCountedKinematicTree &, reco::VertexCollection&, const MagneticField *);

reco::Vertex mkVertex(RefCountedKinematicTree &);
std::pair<double, double> vtxSeparation(const reco::Vertex &, const reco::Vertex &);

#endif
