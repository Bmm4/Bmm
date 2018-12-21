#include "Bmm/CmsswAnalysis/plugins/HfrSimpleRecoilProducer.h"
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

#include "Bmm/RootAnalysis/rootio/TAna01Event.hh"

#include <TLorentzVector.h>

using namespace std;
using namespace edm;
using namespace reco;


// -- Yikes!
extern TAna01Event *gHFEvent;

// ----------------------------------------------------------------------
HfrSimpleRecoilProducer::HfrSimpleRecoilProducer(const ParameterSet& iConfig) :
  HfrBaseProducer(iConfig),
  fVtxProb(iConfig.getUntrackedParameter<double>("vtxProb", 0.01)),
  fCandType(iConfig.getUntrackedParameter<int>("candtype", 3000068)) {
  produces<vector<int> >(); // vector of track indices
  produces<int>();          // possibly a PV index
}


// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::dumpConfiguration() {
  cout << "----------------------------------------------------------------------" << endl;
  cout << "--- HfrBaseProducer::dumpConfiguration()" << endl;
  HfrBaseProducer::dumpConfiguration();
}

// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::put(edm::Event &iEvent,
				  std::unique_ptr<std::vector<int> > &trkIdxColl,
				  std::unique_ptr<int> &pvIdx) {
  iEvent.put(std::move(trkIdxColl));
  iEvent.put(std::move(pvIdx));
}

// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::produce(Event& iEvent, const EventSetup& iSetup) {
  auto trkIdxColl = std::make_unique<vector<int> >();
  auto pvIdx = std::make_unique<int>(1234);

  if (fVerbose > 0)  cout << " ======================================================================"
			  << endl
			  << "=== HfrSimpleRecoilProducer run = " << iEvent.id().run()
			  << " evt = " << iEvent.id().event() << endl
			  << " ----------------------------------------------------------------------"
			  << endl;

  HfrBaseProducer::analyze(iEvent, iSetup);
  // -- get candidate and PV from gHFEvent
  TAnaCand *pCand(0);
  for (int ic = 0; ic < gHFEvent->nCands(); ++ic) {
    pCand = gHFEvent->getCand(ic);
    if (fCandType == pCand->fType) {
      cout << "found one" << endl;
    }
  }

  // -- put into event
  put(iEvent, trkIdxColl, pvIdx);
}

// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::beginRun(const edm::Run& run, const EventSetup& iSetup) {
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
}

// ----------------------------------------------------------------------
void HfrSimpleRecoilProducer::endRun(const edm::Run& run, const EventSetup& iSetup) {
}

DEFINE_FWK_MODULE(HfrSimpleRecoilProducer);
