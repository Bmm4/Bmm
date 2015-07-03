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

#CMSSW_5_3_12_patch3: 
#CMSSW_7_0_7_patch1:  process.GlobalTag.globaltag = "START70_V7A::All"

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100000.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100001.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100002.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100003.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100004.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100005.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100006.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100007.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100008.root",
         "file:/shome/ursl/inclb/CMSSW_7_1_16_patch2/src/PYTHIA8_InclB2Mu_CUEP8M1_8TeV_hardQCD_GEN-100009.root"
 )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "test-pythia8-hard.root"

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


# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genDump*
    process.tree
)
