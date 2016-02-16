import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
dstarpiDump = cms.EDAnalyzer(
    "HFBd2DstarPi",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    trackPt            = cms.untracked.double(4.0),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    slowPionPt         = cms.untracked.double(0.2),
    deltaM             = cms.untracked.double(0.01),
    D0Window           = cms.untracked.double(0.1),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.1),
    maxD0              = cms.untracked.double(2.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(300030)
    )

# ----------------------------------------------------------------------
dstarDump = cms.EDAnalyzer(
    "HFDstar",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    trackPt            = cms.untracked.double(3.5),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    slowPionPt         = cms.untracked.double(0.1),
    D0Window           = cms.untracked.double(0.1),
    deltaM             = cms.untracked.double(0.01),
    deltaR             = cms.untracked.double(0.3),
    maxDoca            = cms.untracked.double(0.1),
    maxD0              = cms.untracked.double(2.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(300054)
    )

# ----------------------------------------------------------------------
ksDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    trackPt            = cms.untracked.double(2.5),
    track1Mass         = cms.untracked.double(0.1396),
    track2Mass         = cms.untracked.double(0.1396),
    massLow            = cms.untracked.double(0.44),
    massHigh           = cms.untracked.double(0.56),
    maxDoca            = cms.untracked.double(0.1),
    pvWeight           = cms.untracked.double(0.60),
    type               = cms.untracked.int32(11310),
    nbrMuons           = cms.untracked.int32(0),
    closeToMuons       = cms.untracked.bool(False)
    )

# ----------------------------------------------------------------------
phiDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    trackPt            = cms.untracked.double(2.5),
    track1Mass         = cms.untracked.double(0.4937),
    track2Mass         = cms.untracked.double(0.4937),
    massLow            = cms.untracked.double(0.98),
    massHigh           = cms.untracked.double(1.10),
    maxDoca            = cms.untracked.double(0.1),
    pvWeight           = cms.untracked.double(0.60),
    type               = cms.untracked.int32(11333),
    nbrMuons           = cms.untracked.int32(0),
    closeToMuons       = cms.untracked.bool(False)
    )

# ----------------------------------------------------------------------
lambdaDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    trackPt            = cms.untracked.double(2.5),
    track1Mass         = cms.untracked.double(0.9383),
    track2Mass         = cms.untracked.double(0.1396),
    massLow            = cms.untracked.double(1.0),
    massHigh           = cms.untracked.double(1.2),
    maxDoca            = cms.untracked.double(0.1),
    pvWeight           = cms.untracked.double(0.60),
    type               = cms.untracked.int32(113122),
    nbrMuons           = cms.untracked.int32(0),
    closeToMuons       = cms.untracked.bool(False)
    )


# ----------------------------------------------------------------------
hhDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0), 
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    trackPt            = cms.untracked.double(4.0),
    track1Mass         = cms.untracked.double(0.1057),
    track2Mass         = cms.untracked.double(0.1057),
    massLow            = cms.untracked.double(4.5),
    massHigh           = cms.untracked.double(6.0),
    maxDoca            = cms.untracked.double(0.025),
    pvWeight           = cms.untracked.double(0.60),
    type               = cms.untracked.int32(211211),
    nbrMuons           = cms.untracked.int32(0),
    closeToMuons       = cms.untracked.bool(False)
    )


# ######################################################################
# Sequences
# ######################################################################
hadronicSequence     = cms.Sequence(dstarpiDump
                                    *dstarDump
                                    *hhDump)
