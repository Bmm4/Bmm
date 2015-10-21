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
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag.globaltag = "74X_dataRun2_Prompt_v1"

# ----------------------------------------------------------------------
# POOLSOURCE


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "inclb-prompt-Run2015B-XXXX.root"

# load all processes
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Bmm.CmsswAnalysis.HFRecoStuff_cff")
process.load("Bmm.CmsswAnalysis.HFMCTruth_cff")

process.load("Bmm.CmsswAnalysis.HFInclBSignal_cff")
process.load("Bmm.CmsswAnalysis.HFGenJets_cff")
#process.load("Bmm.CmsswAnalysis.HFPFJets_cff")
process.load("Bmm.CmsswAnalysis.HFTrackJets_cff")
process.load("Bmm.CmsswAnalysis.HFCandidate_cff")
process.load("Bmm.CmsswAnalysis.HFCandidateNew_cff")


# ----------------------------------------------------------------------
process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose        = cms.untracked.int32(0),
    printFrequency = cms.untracked.int32(100),
    requireCand    =  cms.untracked.bool(False),
    fileName       = cms.untracked.string(rootFileName)
    )


#process.HFInclBMuonTrackJets.verbose = cms.untracked.int32(10)

# ----------------------------------------------------------------------
process.p = cms.Path(
#    process.genDump
    process.recoStuffSequence
#    *process.PFJetDumpAOD
    *process.TrackJetDumpAOD
    *process.HFInclBSignalSequence
    *process.tree
    )
