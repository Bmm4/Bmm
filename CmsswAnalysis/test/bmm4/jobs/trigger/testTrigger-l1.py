# ######################################################################
# t3ui01
# /mnt/t3nfs01/data01/shome/ursl/bmm4/CMSSW_8_0_22/src/Bmm/CmsswAnalysis/test/bmm4/jobs/v05/RunIISpring16DR80
# file list contains 39355 events
# mkPyFiles -t ../bmm-mc-RunIISpring16DR80-XXXX.py -f ../../../catalogs/RunIISpring16DR80/BdToPiPi_BMuonFilter -s v05 -e 30000
# ./bmm-mc-RunIISpring16DR80-BdToPiPi_BMuonFilter-v05-0024.py with 39355 events
# ######################################################################
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

FILE1 = os.environ.get('FILE1')

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

process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')

# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
         "file:/scratch/ursl/bmm5-trigger/$FILE1"
 )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )



process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "T1-$FILE1"

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
process.load("Bmm.CmsswAnalysis.HFBmm_cff")
process.load("Bmm.CmsswAnalysis.HFPhysicsDeclared_cff")

process.load("Bmm.CmsswAnalysis.HFMCTruth_cff")
process.load("Bmm.CmsswAnalysis.HFTruthCandidates_cff")

# ----------------------------------------------------------------------
process.genDump = cms.EDAnalyzer(
    "HFDumpGenerator",
    generatorCandidates = cms.untracked.string('genParticles'),
    generatorEvent = cms.untracked.string('generator')
    )

process.triggerDump.HLTProcessName = cms.untracked.string('REHLT')
process.triggerDump.HLTResultsLabel = cms.untracked.InputTag('TriggerResults::REHLT')
process.triggerDump.TriggerEventLabel = cms.untracked.InputTag('hltTriggerSummaryAOD::REHLT')


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genDump*
    process.recoStuffSequence*
    process.bmmSequence*
    process.truthBmmSequence*
    process.tree
)
