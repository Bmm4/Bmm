import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
psiDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    unlikeCharge       = cms.untracked.bool(True),
    muonPt             = cms.untracked.double(4.0),
    chi2               = cms.untracked.double(5.),
    type               = cms.untracked.int32(1313),
    maxDoca            = cms.untracked.double(0.08),
    candlo             = cms.untracked.double(2.7),
    candhi             = cms.untracked.double(4.5),
    pvWeight           = cms.untracked.double(0.6),
    nbrMuons           = cms.untracked.int32(2),
    closeToMuons       = cms.untracked.bool(False)
    )

# ----------------------------------------------------------------------
upsDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    unlikeCharge       = cms.untracked.bool(True),
    chi2               = cms.untracked.double(5.),
    muonPt             = cms.untracked.double(4.0),
    type               = cms.untracked.int32(1313),
    maxDoca            = cms.untracked.double(0.08),
    candlo             = cms.untracked.double(8.0),
    candhi             = cms.untracked.double(12.0),
    pvWeight           = cms.untracked.double(0.6),
    nbrMuons           = cms.untracked.int32(2),
    closeToMuons       = cms.untracked.bool(False)
    )

# ######################################################################
# Sequences
# ######################################################################
oniaSequence     = cms.Sequence(psiDump*upsDump)
