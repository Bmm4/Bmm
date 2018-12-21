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
#import CondCore.CondDB.CondDB_cfi
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
    fileNames = cms.untracked.vstring(
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil2/recoil2_step2-40004.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil2/recoil2_step2-40006.root"
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "testRecoil2.root"

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
#process.load("Bmm.CmsswAnalysis.HFRecoil_cff")

process.load("Bmm.CmsswAnalysis.HFMCTruth_cff")
process.load("Bmm.CmsswAnalysis.HFTruthCandidates_cff")
# ----------------------------------------------------------------------
process.genDump = cms.EDAnalyzer(
    "HFDumpGenerator",
    generatorCandidates = cms.untracked.string('genParticles'),
    generatorEvent = cms.untracked.string('generator')
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
process.hfrjpsix = cms.EDProducer(
    "HfrJpsiXListProducer",
    verbose    = cms.untracked.int32(0),
    trackMinPt = cms.untracked.double(0.1),
    muonMinPt = cms.untracked.double(4.123)
)

# ----------------------------------------------------------------------
process.hfrTest = cms.EDAnalyzer(
    "HfrTest",
    verbose            = cms.untracked.int32(0),
    trkIdxLabel        = cms.untracked.InputTag("hfrjpsix"),
    pvIdxLabel         = cms.untracked.InputTag("hfrjpsix"),
    vertexLabel         = cms.untracked.InputTag("hfrjpsix")
)

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genDump*
    process.recoStuffSequence*
    process.truthBsToMuMuDump*
#    process.hfrTest*
    process.hfrjpsix*
    process.bupsikpDump*
    process.tree
)
