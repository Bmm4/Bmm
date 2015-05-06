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

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_mc', '')

# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
         "/store/mc/Summer12_DR53X/BuToJPsiK_K2MuPtEtaEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v2/0000/521047B4-E0DD-E111-8146-02215E94E7DF.root"
#        "file:./7CA64471-C16B-E411-BB89-02163E00F463.root"
#        "file:./F49EE52E-8D6C-E411-A720-02163E010DC8.root"
 )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "test.root"

# load all processes
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Bmm.CmsswAnalysis.HFRecoStuff_cff")
process.load("Bmm.CmsswAnalysis.HFMCTruth_cff")

process.load("Bmm.CmsswAnalysis.HFSignal_cff")
process.load("Bmm.CmsswAnalysis.HFGenJets_cff")
process.load("Bmm.CmsswAnalysis.HFPFJets_cff")
process.load("Bmm.CmsswAnalysis.HFCaloJets_cff")
process.load("Bmm.CmsswAnalysis.HFTrackJets_cff")
process.load("Bmm.CmsswAnalysis.HFCandidate_cff")
process.load("Bmm.CmsswAnalysis.HFCandidateNew_cff")


# ----------------------------------------------------------------------
rootFileName = "test.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose        = cms.untracked.int32(1),
    printFrequency = cms.untracked.int32(100),
    requireCand    =  cms.untracked.bool(False),
    fileName       = cms.untracked.string(rootFileName)
    )


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genDump
    *process.recoStuffSequence
    *process.PFJetDumpAOD
    *process.CaloJetDumpAOD
    *process.SignalDump
    *process.tree
    )
