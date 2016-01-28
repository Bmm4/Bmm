// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// HFTrackListBuilder
// ------------------
//
// 2012/11/10 Christoph Naegeli    first shot
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#include <algorithm>

#include "Bmm/CmsswAnalysis/interface/HFTrackListBuilder.hh"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"

using namespace std;

// ----------------------------------------------------------------------
HFTrackListBuilder::HFTrackListBuilder(edm::Handle<edm::View<reco::Track> > &hTracks, const reco::MuonCollection *muons, const TransientTrackBuilder *ttb, int verbose) :
  fhTracks(hTracks),
  fMuons(muons),
  fTTB(ttb),
  fVerbose(verbose),
  fCallerName("HFTrackListBuilder"),
  fMaxD0(999.),
  fMaxDz(999.),
  fMinPt(-1.),
  fMaxDocaToTrks(999.),
  fTrackQuality("highPurity"),
  fMuonQuality(muon::AllGlobalMuons), 
  fCloseTracks(NULL)
{}


// ----------------------------------------------------------------------
std::vector<int> HFTrackListBuilder::getMuonList() {
  reco::MuonCollection::const_iterator muonIt;
  vector<int> trackList;
  vector<int>::iterator trackIt;
	
  trackList.reserve(20); // 20 muons should be enough

  int ix(0); 
  for (muonIt = fMuons->begin(); muonIt != fMuons->end(); ++muonIt) {
    if (!muon::isGoodMuon(*muonIt, fMuonQuality)) continue;
    int ixMu = muonIt->track().index();
    if (ixMu >= 0) trackList.push_back(ixMu);
    if (fVerbose > 0) cout << "==>" << fCallerName << "> muon " << ix << " with track index " << ixMu << endl;
    ++ix; 
  }
	
  if (fVerbose > 0) {
    cout << "==>" << fCallerName << "> nMuons = " << fMuons->size() << endl;
    cout << "==>" << fCallerName << "> nMuonIndices = " << trackList.size() << endl;
  }
	
  trackIt = std::remove_if(trackList.begin(), trackList.end(), *this);
  trackList.erase(trackIt,trackList.end());
	
  return trackList;
} 


// ----------------------------------------------------------------------
std::vector<int> HFTrackListBuilder::getTrackList() {
  vector<int> trackList; // allocate capacity
  int ix;
	
  trackList.reserve(300);
  for (ix = 0; (unsigned)ix < fhTracks->size(); ix++) {
    reco::TrackBaseRef rTrackView(fhTracks, ix);
    const reco::Track trackView(*rTrackView);
    if (!trackView.quality(reco::TrackBase::qualityByName(fTrackQuality))) continue;
    if (!(*this)(ix)) trackList.push_back(ix);
  }
	
  return trackList;
}


// ----------------------------------------------------------------------
bool HFTrackListBuilder::operator()(int ix) {
  reco::TrackBaseRef rTrackView(fhTracks,ix);
  reco::Track tTrack(*rTrackView);
  bool result = tTrack.d0() > fMaxD0 || tTrack.dz() > fMaxDz || tTrack.pt() < fMinPt;
	
  if (!result && fCloseTracks) {
    // check wether this track is nearby anyone in the fCloseTracks vector
    reco::TransientTrack tTrkCur = fTTB->build(tTrack);
    TwoTrackMinimumDistance md;
    double minDoca = DBL_MAX;
    size_t j;
		
    for (j = 0; j < fCloseTracks->size(); j++) {
      if ((*fCloseTracks)[j] == ix) // compare to itself
	continue;
			
      reco::TransientTrack tTrkCompare = fTTB->build((*fhTracks)[(*fCloseTracks)[j]]);
      md.calculate(tTrkCur.initialFreeState(),tTrkCompare.initialFreeState());
      if (md.distance() < minDoca)
	minDoca = md.distance();
    }
		
    result = minDoca > fMaxDocaToTrks;
  }
	
  return result;
}
