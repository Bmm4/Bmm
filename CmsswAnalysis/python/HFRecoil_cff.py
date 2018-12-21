import FWCore.ParameterSet.Config as cms

hfrtracks = cms.EDProducer(
    "HfrTrackListProducer",
    verbose    = cms.untracked.int32(1),
    trackMinPt = cms.untracked.double(1.0)
)

hfrjpsix = cms.EDProducer(
    "HfrJpsiXListProducer",
    verbose    = cms.untracked.int32(1),
    trackMinPt = cms.untracked.double(1.0),
    muonMinPt = cms.untracked.double(4.123)
)

# ----------------------------------------------------------------------
hfrTest = cms.EDAnalyzer(
    "HfrTest",
    verbose            = cms.untracked.int32(1),
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
)
