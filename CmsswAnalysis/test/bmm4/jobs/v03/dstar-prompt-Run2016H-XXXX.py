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

# requires >= CMSSW_8_0_19
process.GlobalTag.globaltag = "80X_dataRun2_Prompt_v12"

# ----------------------------------------------------------------------
# POOLSOURCE
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())



# ----------------------------------------------------------------------
rootFileName = "dstar-prompt-Run2016G-XXXX.root"

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
process.load("Bmm.CmsswAnalysis.HFHadronic_cff")
process.load("Bmm.CmsswAnalysis.HFOnia_cff")
process.load("Bmm.CmsswAnalysis.HFPhysicsDeclared_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.recoStuffSequence*
    process.ksDump*
    process.phiDump*
    process.lambdaDump*
    process.psiDump*
    process.tree
)
