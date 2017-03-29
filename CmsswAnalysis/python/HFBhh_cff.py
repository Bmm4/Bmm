import FWCore.ParameterSet.Config as cms


# ----------------------------------------------------------------------
hhDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"),
    muonPt             = cms.untracked.double(4.0),
    pt                 = cms.untracked.double(5.),
    chi2               = cms.untracked.double(10.),
    pvips              = cms.untracked.double(5.),
    flsxy              = cms.untracked.double(4.),
    flxy               = cms.untracked.double(9999.),
    type               = cms.untracked.int32(1313),
    candlo             = cms.untracked.double(4.9),
    candhi             = cms.untracked.double(5.9),
    extra              = cms.untracked.double(0.3),
    maxDoca            = cms.untracked.double(0.08),
    pvWeight           = cms.untracked.double(0.6),
    nbrMuons           = cms.untracked.int32(0),
    closeToMuons       = cms.untracked.bool(False)
    )

# ######################################################################
# Sequences
# ######################################################################
bhhSequence     = cms.Sequence(hhDump)
