import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('HLTrigReport')
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

# requires >= CMSSW_8_0_19
process.GlobalTag.globaltag = "92X_dataRun2_Prompt_v4"

# ----------------------------------------------------------------------
# POOLSOURCE
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())



# ----------------------------------------------------------------------
rootFileName = "bmm-prompt-Run2017B-XXXX.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose        = cms.untracked.int32(0),
    printFrequency = cms.untracked.int32(100),
    requireCand    =  cms.untracked.bool(True),
    fullGenBlock   = cms.untracked.bool(False),
    fileName       = cms.untracked.string(rootFileName)
    )


# ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Bmm.CmsswAnalysis.HFRecoStuff_cff")
process.load("Bmm.CmsswAnalysis.HFBmm_cff")
process.load("Bmm.CmsswAnalysis.HFOnia_cff")
process.load("Bmm.CmsswAnalysis.HFHadronic_cff")
process.load("Bmm.CmsswAnalysis.HFPhysicsDeclared_cff")
process.load("Bmm.CmsswAnalysis.HFSkipEvents_cff")

process.skipEvents.filterOnJson = cms.untracked.int32(1)
process.skipEvents.JSONFile     = cms.untracked.FileInPath("Bmm/RootAnalysis/common/json/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt")

# ----------------------------------------------------------------------
process.p = cms.Path(
#    process.skipEvents*
    process.recoStuffSequence*
    process.bmmSequence*
    process.psiDump*
    process.upsDump*
    process.dstarDump*
    process.ksDump*
    process.phiDump*
    process.lambdaDump*
    process.tree
)
