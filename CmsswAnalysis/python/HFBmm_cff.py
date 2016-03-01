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
    flsxy              = cms.untracked.double(2.5),
    flxy               = cms.untracked.double(9999.),
    type               = cms.untracked.int32(1313), 
    massLow            = cms.untracked.double(4.2), 
    massHigh           = cms.untracked.double(6.7),
    maxDoca            = cms.untracked.double(0.06),
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
    flsxy              = cms.untracked.double(3.0),
    flxy               = cms.untracked.double(9999.),
    psiMuons           = cms.untracked.int32(2),
    psiWindow          = cms.untracked.double(0.2),
    dimuonpt           = cms.untracked.double(6.9),
    BuWindow           = cms.untracked.double(1.0),
    candpt             = cms.untracked.double(-1.0),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.06),
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
    flsxy              = cms.untracked.double(2.5),
    flxy               = cms.untracked.double(9999.),
    psiMuons           = cms.untracked.int32(2),
    psiWindow          = cms.untracked.double(0.2),
    phiWindow          = cms.untracked.double(0.1),
    dimuonpt           = cms.untracked.double(6.9),
    BsWindow           = cms.untracked.double(0.7),
    candpt             = cms.untracked.double(-1.0),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.06),
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
    psiMuons           = cms.untracked.int32(2),
    psiWindow          = cms.untracked.double(0.2),
    kstarWindow        = cms.untracked.double(0.1),
    BdWindow           = cms.untracked.double(0.7),
    trackPt            = cms.untracked.double(0.5),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.06),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(511)
    )

# ######################################################################
# Sequences
# ######################################################################
bmmSequence     = cms.Sequence(bmmDump*bupsikpDump*bspsiphiDump*bdpsikstarDump)
