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
  fCloseTracks(NULL),
  fDoFilter(false)
{}


// ----------------------------------------------------------------------
std::vector<int> HFTrackListBuilder::getMuonList() {
  reco::MuonCollection::const_iterator muonIt;
  vector<int> trackList;
  vector<int>::iterator trackIt;

  if (fDoFilter && fVerbose > 0) {
    cout << fCallerName
	 << "::getMuonList() filter track list has " << fRecoilTrkIdx.size()
	 << " entries" << endl;
  }

  trackList.reserve(50); // 50 muons should be enough
  int ix(0);

  for (muonIt = fMuons->begin(); muonIt != fMuons->end(); ++muonIt) {
    if (!muon::isGoodMuon(*muonIt, fMuonQuality)) continue;
    int ixMu = muonIt->track().index();
    if (ixMu < 0) continue;
    // -- recoil filtering
    if (fDoFilter) {
      if (fRecoilTrkIdx.end() !=  find(fRecoilTrkIdx.begin(), fRecoilTrkIdx.end(), ixMu)) {
	if (fVerbose > 0) {
	  cout << "found ixMu = " << ixMu << " in fRecoilTrkIdx" << endl;
	}
      } else {
	continue;
      }
    }

    if (!(*this)(ixMu)) {
      trackList.push_back(ixMu);
    } else {
      reco::TrackBaseRef rTrackView(fhTracks, ixMu);
      reco::Track tTrack(*rTrackView);
    }
    ++ix;
  }

  if (fVerbose > 0) {
    cout << "==>" << fCallerName << "> nMuons = " << fMuons->size() << endl;
    cout << "==>" << fCallerName << "> nMuonIndices = " << trackList.size() << ": ";
    for (unsigned int im = 0; im < trackList.size(); ++im) {
      cout << trackList[im] << " ";
    }
    cout << endl;
  }

  return trackList;
}


// ----------------------------------------------------------------------
std::vector<int> HFTrackListBuilder::getTrackList() {
  vector<int> trackList; // allocate capacity
  int ix;

  if (fDoFilter && (fVerbose > 0)) {
    cout << fCallerName
	 << "::getTrackList() filter track list has " << fRecoilTrkIdx.size()
	 << " entries, fCloseTracks->size() = " << (fCloseTracks? fCloseTracks->size(): 0)
	 << endl;
  }

  trackList.reserve(300);
  if (fVerbose > 0) {
    cout << "==>" << fCallerName << "> nTracks = " << fhTracks->size() << " filter: " << fDoFilter << endl;
  }

  for (ix = 0; (unsigned)ix < fhTracks->size(); ix++) {
    reco::TrackBaseRef rTrackView(fhTracks, ix);
    const reco::Track trackView(*rTrackView);
    if (!trackView.quality(reco::TrackBase::qualityByName(fTrackQuality))) continue;
    // -- recoil filtering
    if (fDoFilter) {
      if (fRecoilTrkIdx.end() !=  find(fRecoilTrkIdx.begin(), fRecoilTrkIdx.end(), ix)) {
	//cout << "found ix = " << ix << " in fRecoilTrkIdx" << endl;
      } else {
	continue;
      }
    }
    if (!(*this)(ix)) {
      // cout << "added ix = " << ix << " to trackList" << endl;
      trackList.push_back(ix);
    } else {
      // cout << "did NOT add ix = " << ix << " to trackList" << endl;
    }
  }

  return trackList;
}


// ----------------------------------------------------------------------
bool HFTrackListBuilder::operator()(int ix) {
  reco::TrackBaseRef rTrackView(fhTracks,ix);
  reco::Track tTrack(*rTrackView);
  bool skip = tTrack.d0() > fMaxD0 || tTrack.dz() > fMaxDz || tTrack.pt() < fMinPt;
  if (!skip && fCloseTracks) {
    // check wether this track is nearby anyone in the fCloseTracks vector
    reco::TransientTrack tTrkCur = fTTB->build(tTrack);
    TwoTrackMinimumDistance md;
    double minDoca = DBL_MAX;
    size_t j;

    for (j = 0; j < fCloseTracks->size(); j++) {
      if ((*fCloseTracks)[j] == ix) // compare to itself
	continue;

      reco::TransientTrack tTrkCompare = fTTB->build((*fhTracks)[(*fCloseTracks)[j]]);
      md.calculate(tTrkCur.initialFreeState(), tTrkCompare.initialFreeState());
      if (md.distance() < minDoca)
	minDoca = md.distance();
    }
    skip = minDoca > fMaxDocaToTrks;
  }

  return skip;
}
