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
process.load("CondCore.DBCommon.CondDBCommon_cfi")
process.load("CondCore.DBCommon.CondDBSetup_cfi")
# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

# GT for 76X (X>=3) 2015 End of the Yerar data-reprocessing
process.GlobalTag.globaltag = "76X_dataRun2_v15"

# ----------------------------------------------------------------------
# POOLSOURCE



# ----------------------------------------------------------------------
rootFileName = "dstar-rereco-Run2015-XXXX.root"

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
process.load("Bmm.CmsswAnalysis.HFHadronic_cff")
process.load("Bmm.CmsswAnalysis.HFPhysicsDeclared_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.recoStuffSequence*
    process.dstarDump*
    process.ksDump*
    process.phiDump*
    process.lambdaDump*
    process.tree
)
