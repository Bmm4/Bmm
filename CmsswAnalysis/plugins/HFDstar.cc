#include "HFDstar.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include "Bmm/RootAnalysis/common/HFMasses.hh"
#include "Bmm/CmsswAnalysis/interface/HFTwoParticleCombinatoricsNew.hh"

using std::cout;
using std::endl;

// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HFDstar::HFDstar(const edm::ParameterSet& iConfig) :
  HFVirtualDecay(iConfig),
  fSlowPionPt(iConfig.getUntrackedParameter<double>("slowPionPt", 0.4)),
  fD0Window(iConfig.getUntrackedParameter<double>("D0Window", 0.1)),
  fDeltaM(iConfig.getUntrackedParameter<double>("deltaM", 0.03))
{
  dumpConfiguration();
} // HFDstar()


// ----------------------------------------------------------------------
void HFDstar::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HFDstar constructor" << endl;
  HFVirtualDecay::dumpConfiguration();
  cout << "---  slowPionPt:               " << fSlowPionPt << endl;
  cout << "---  D0Window:                 " << fD0Window << endl;
  cout << "---  deltaM:                   " << fDeltaM << endl;
  cout << "----------------------------------------------------------------------" << endl;
} // dumpConfiguration()

// ----------------------------------------------------------------------
void HFDstar::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  const double MDSTAR0(1.9), MDSTAR1(2.2), MD00(1.8), MD01(2.0);
  typedef HFTwoParticleCombinatoricsNew::HFTwoParticleCombinatoricsSet HFTwoParticleCombinatoricsSet;
	
  try {
    HFVirtualDecay::analyze(iEvent,iSetup);
  }
  catch (HFSetupException e) {
    cout << "==>HFDstar> " << e.fMsg << endl;
    return; // problem with setup
  }
	
  std::vector<int> trkFastList = fListBuilder->getTrackList();
  fListBuilder->setMinPt(fSlowPionPt);
  fListBuilder->setMaxDocaToTracks(fMaxDoca);
  fListBuilder->setCloseTracks(&trkFastList);
  std::vector<int> trkSlowList = fListBuilder->getTrackList();
	
  HFTwoParticleCombinatoricsNew a(fTracksHandle,fVerbose);
  HFTwoParticleCombinatoricsSet kapiSet = a.combine(trkFastList, MKAON, trkFastList, MPION, MD0-fD0Window,MD0+fD0Window, 0);
	
  TLorentzVector ka,pi,d0,pis,dstar;
  for (HFTwoParticleCombinatoricsNew::iterator d0It = kapiSet.begin(); d0It != kapiSet.end(); ++d0It) {
		
    unsigned int iKaon = d0It->first;
    unsigned int iPion = d0It->second;
		
    reco::TrackBaseRef kaonTrackView(fTracksHandle, iKaon);
    reco::Track tKaon(*kaonTrackView);
    ka.SetPtEtaPhiM(tKaon.pt(), tKaon.eta(), tKaon.phi(), MKAON);
		
    reco::TrackBaseRef pionTrackView(fTracksHandle, iPion);
    reco::Track tPion(*pionTrackView);
    pi.SetPtEtaPhiM(tPion.pt(), tPion.eta(), tPion.phi(), MPION);
		
    d0 = ka + pi;
		
    for (unsigned iTrack = 0; iTrack < trkSlowList.size(); iTrack++) {
		
      if (iTrack == iKaon || iTrack == iPion) continue;
      reco::TrackBaseRef rTrackView(fTracksHandle, iTrack);
      reco::Track tSlowPion(*rTrackView);
			
      // -- the slow pion has the same charge like the fast pion
      if (tSlowPion.charge()*tPion.charge() < 0) continue;
			
      pis.SetXYZM(tSlowPion.px(), tSlowPion.py(), tSlowPion.pz(), MPION);
      if (d0.DeltaR(pis) > fDeltaR) continue;
			
      dstar = d0 + pis;
      if (TMath::Abs(dstar.M() - d0.M()) > 0.16) continue;
			
      // -- sequential fit: D0 pi_slow
      if (fVerbose > 5) {
	cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
	cout << "==>HFDstar> going to sequential fit with track indices: " 
	     << iKaon << " " << iPion << " " << iTrack
	     << endl;
      }
			
      HFDecayTree theTree(fType, true, MDSTARPLUS, false, -1.0, true);
      HFDecayTreeIterator iterator = theTree.addDecayTree(300050, true, MD0, false);
      iterator->addTrack(iKaon,321);
      iterator->addTrack(iPion,211);
      //      iterator->setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.maxDoca), &(iterator->fTV.maxDocaV), -1., fMaxDoca, "300050 maxdoca"));
      iterator->addSimpleCut(HFSimpleCut(&(iterator->fTV.mass), &(iterator->fTV.massV), MD00, MD01, "300050 mass"));
      theTree.addTrack(iTrack,211);
      //      theTree.setNodeCut(RefCountedHFNodeCut(new HFMaxDocaCut(fMaxDoca)));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.maxDoca), &(theTree.fTV.maxDocaV), -1., fMaxDoca, Form("%d maxdoca", fType)));
      theTree.addSimpleCut(HFSimpleCut(&(theTree.fTV.mass), &(theTree.fTV.massV), MDSTAR0, MDSTAR1, Form("%d mass", fType)));
			
      fSequentialFitter->doFit(&theTree);
    }
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HFDstar);
