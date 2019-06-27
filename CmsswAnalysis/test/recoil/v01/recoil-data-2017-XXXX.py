import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# ----------------------------------------------------------------------
# -- Database configuration
process.load("CondCore.CondDB.CondDB_cfi")
# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11', '')

# ----------------------------------------------------------------------
# POOLSOURCE



process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "recoil-data-2017-XXXX.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose        = cms.untracked.int32(0),
    printFrequency = cms.untracked.int32(100),
    requireCand    =  cms.untracked.bool(False),
    fullGenBlock   = cms.untracked.bool(True),
    fileName       = cms.untracked.string(rootFileName)
    )


# ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Bmm.CmsswAnalysis.HFRecoStuff_cff")


# ----------------------------------------------------------------------
process.hfrjpsixrecoil = cms.EDProducer(
    "HfrJpsiXRecoilProducer",
    verbose    = cms.untracked.int32(0),
    type       = cms.untracked.int32(10000),
    trackMinPt = cms.untracked.double(1.0),
    muonMinPt  = cms.untracked.double(3.0),
    docaMaxPv  = cms.untracked.double(0.1),
    docaMaxTk  = cms.untracked.double(0.05)
)


# ----------------------------------------------------------------------
process.butomutaukDump = cms.EDAnalyzer(
    "HFBu2MuTauK",
    verbose            = cms.untracked.int32(0),
    filterLabel        = cms.untracked.InputTag("hfrjpsixrecoil"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"),
    muonPt             = cms.untracked.double(4.0),
    chi2               = cms.untracked.double(5.),
    pvips              = cms.untracked.double(5.),
    flsxy              = cms.untracked.double(1.0),
    flxy               = cms.untracked.double(9999.),
    candpt             = cms.untracked.double(-1.0),
    candlo             = cms.untracked.double(0.0),
    candhi             = cms.untracked.double(6.0),
    trackPt            = cms.untracked.double(0.4),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.08),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    mukaMaxDoca        = cms.untracked.double(0.02),
    tauMaxDoca         = cms.untracked.double(0.02),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(10032),
    mcType             = cms.untracked.int32(-1)
)


# ----------------------------------------------------------------------
process.butopipidDump = cms.EDAnalyzer(
    "HFBu2MuTauK",
    verbose            = cms.untracked.int32(0),
    filterLabel        = cms.untracked.InputTag("hfrjpsixrecoil"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"),
    muonPt             = cms.untracked.double(4.0),
    chi2               = cms.untracked.double(5.),
    pvips              = cms.untracked.double(5.),
    flsxy              = cms.untracked.double(1.0),
    flxy               = cms.untracked.double(9999.),
    candpt             = cms.untracked.double(-1.0),
    candlo             = cms.untracked.double(0.0),
    candhi             = cms.untracked.double(6.0),
    trackPt            = cms.untracked.double(0.4),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.08),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    mukaMaxDoca        = cms.untracked.double(0.02),
    tauMaxDoca         = cms.untracked.double(0.02),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(10041),
    mcType             = cms.untracked.int32(-1)
)


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.recoStuffSequence*
    process.hfrjpsixrecoil*
    process.butomutaukDump*
    process.butopipidDump*
    process.tree
)
