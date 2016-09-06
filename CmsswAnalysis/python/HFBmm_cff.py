import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
bmmDump = cms.EDAnalyzer(
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
    nbrMuons           = cms.untracked.int32(2),
    closeToMuons       = cms.untracked.bool(False)
    )

# ----------------------------------------------------------------------
bupsikpDump = cms.EDAnalyzer(
    "HFBu2JpsiKp",
    verbose            = cms.untracked.int32(0),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"),
    muonPt             = cms.untracked.double(4.0),
    chi2               = cms.untracked.double(10.),
    pvips              = cms.untracked.double(5.),
    flsxy              = cms.untracked.double(4.0),
    flxy               = cms.untracked.double(9999.),
    psiMuons           = cms.untracked.int32(2),
    psiLo              = cms.untracked.double(2.8),
    psiHi              = cms.untracked.double(3.3),
    dimuonpt           = cms.untracked.double(7.0),
    candpt             = cms.untracked.double(-1.0),
    candlo             = cms.untracked.double(4.9),
    candhi             = cms.untracked.double(5.9),
    trackPt            = cms.untracked.double(0.6),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.08),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(521)
    )

# ----------------------------------------------------------------------
bspsiphiDump = cms.EDAnalyzer(
    "HFBs2JpsiPhi",
    verbose            = cms.untracked.int32(0),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"),
    muonPt             = cms.untracked.double(4.0),
    chi2               = cms.untracked.double(10.),
    pvips              = cms.untracked.double(5.),
    flsxy              = cms.untracked.double(4.0),
    flxy               = cms.untracked.double(9999.),
    psiMuons           = cms.untracked.int32(2),
    psiLo              = cms.untracked.double(2.8),
    psiHi              = cms.untracked.double(3.3),
    phiLo              = cms.untracked.double(0.99),
    phiHi              = cms.untracked.double(1.04),
    dimuonpt           = cms.untracked.double(7.0),
    candpt             = cms.untracked.double(-1.0),
    candlo             = cms.untracked.double(4.9),
    candhi             = cms.untracked.double(5.9),
    trackPt            = cms.untracked.double(0.6),
    deltaR             = cms.untracked.double(0.8),
    maxDoca            = cms.untracked.double(0.08),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(531)

    )

# ----------------------------------------------------------------------
bdpsikstarDump = cms.EDAnalyzer(
    "HFBd2JpsiKstar",
    verbose            = cms.untracked.int32(0),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"),
    muonPt             = cms.untracked.double(4.0),
    chi2               = cms.untracked.double(10.),
    pvips              = cms.untracked.double(5.),
    flsxy              = cms.untracked.double(4.0),
    flxy               = cms.untracked.double(9999.),
    psiMuons           = cms.untracked.int32(2),
    psiLo              = cms.untracked.double(2.8),
    psiHi              = cms.untracked.double(3.3),
    kstarLo            = cms.untracked.double(0.75),
    kstarHi            = cms.untracked.double(1.05),
    dimuonpt           = cms.untracked.double(7.0),
    candpt             = cms.untracked.double(-1.0),
    candlo             = cms.untracked.double(4.9),
    candhi             = cms.untracked.double(5.9),
    trackPt            = cms.untracked.double(0.8),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.08),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(511)
    )

# ######################################################################
# Sequences
# ######################################################################
bmmSequence     = cms.Sequence(bmmDump*bupsikpDump*bspsiphiDump*bdpsikstarDump)
