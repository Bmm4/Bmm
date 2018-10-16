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
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40001.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40002.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40003.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40004.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40005.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40007.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40008.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40009.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40010.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40011.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40012.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40013.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40014.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40015.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40016.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40018.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40019.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40020.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40021.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40022.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40023.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40024.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40025.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40026.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40027.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40028.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40029.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40030.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40031.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40032.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40033.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40034.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40035.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40036.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40037.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40038.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40039.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40040.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40041.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40042.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40043.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40044.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40045.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40046.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40047.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40048.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40049.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40050.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40051.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40052.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40053.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40054.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40055.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40056.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40057.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40058.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40059.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40060.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40061.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40062.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40063.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40064.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40065.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40066.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40067.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40068.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40069.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40070.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40071.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40072.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40073.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40074.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40075.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40077.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40078.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40079.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40080.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40081.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40082.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40083.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40084.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40085.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40086.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40087.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40088.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40089.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40090.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40091.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40092.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40093.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40094.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40095.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40096.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40097.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40098.root",
        "/store/user/ursl/recoil/mcprod2017/aodsim/recoil10/recoil10_step2-40099.root"
    )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "testRecoil10.root"

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

process.load("Bmm.CmsswAnalysis.HFMCTruth_cff")
process.load("Bmm.CmsswAnalysis.HFTruthCandidates_cff")
# ----------------------------------------------------------------------
process.genDump = cms.EDAnalyzer(
    "HFDumpGenerator",
    generatorCandidates = cms.untracked.string('genParticles'),
    generatorEvent = cms.untracked.string('generator')
    )


# ----------------------------------------------------------------------
process.bmmDump = cms.EDAnalyzer(
    "HFDiTracks",
    verbose            = cms.untracked.int32(0),
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag("generalTracks"),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"),
    muonPt             = cms.untracked.double(4.0),
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
process.hfrsimple = cms.EDProducer(
    "HfrSimpleRecoilProducer",
    verbose    = cms.untracked.int32(0),
)

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genDump*
    process.recoStuffSequence*
    process.truthBsToMuMuDump*
    process.truthBu2JpsiKpDump*
    process.bupsikpDump*
    process.bmmDump*
    process.hfrsimple*
    process.tree
)
