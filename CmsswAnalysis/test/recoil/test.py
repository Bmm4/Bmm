import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


# ----------------------------------------------------------------------
# -- Database configuration
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

#CMSSW_5_3_12_patch3:
#CMSSW_7_0_7_patch1:  process.GlobalTag.globaltag = "START70_V7A::All"

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
    fileNames = cms.untracked.vstring("/store/user/ursl/recoil/mcprod2017/aodsim/recoil0/recoil0_step2-40000.root"

    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "test.root"

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
process.load("Bmm.CmsswAnalysis.HFRecoil_cff")

process.load("Bmm.CmsswAnalysis.HFMCTruth_cff")

process.bmmDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(1),
    filterLabel        = cms.untracked.InputTag("hfrtracks"),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"),
    muonPt             = cms.untracked.double(4.123),
    pt                 = cms.untracked.double(5.),
    chi2               = cms.untracked.double(5.),
    pvips              = cms.untracked.double(5.),
    flsxy              = cms.untracked.double(4.),
    flxy               = cms.untracked.double(9999.),
    type               = cms.untracked.int32(1313),
    candlo             = cms.untracked.double(4.8),
    candhi             = cms.untracked.double(6.0),
    extra              = cms.untracked.double(0.3),
    maxDoca            = cms.untracked.double(0.08),
    pvWeight           = cms.untracked.double(0.6),
    nbrMuons           = cms.untracked.int32(2),
    closeToMuons       = cms.untracked.bool(False)
    )

# ----------------------------------------------------------------------
process.bupsikpDump = cms.EDAnalyzer(
    "HFBu2JpsiKp",
    verbose            = cms.untracked.int32(2),
    filterLabel        = cms.untracked.InputTag("hfrtracks"),
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
    psiMuons           = cms.untracked.int32(2),
    psiLo              = cms.untracked.double(2.8),
    psiHi              = cms.untracked.double(3.3),
    dimuonpt           = cms.untracked.double(7.0),
    candpt             = cms.untracked.double(-1.0),
    candlo             = cms.untracked.double(4.8),
    candhi             = cms.untracked.double(6.0),
    trackPt            = cms.untracked.double(0.0),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.08),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(521)
    )


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.recoStuffSequence*
#    process.hfrtracks*
    process.hfrjpsix*
    ###    process.hfrTest*
#    process.bmmDump*
#    process.bupsikpDump*
    process.tree
)
