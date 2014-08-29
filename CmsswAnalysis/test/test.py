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
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#CMSSW_5_3_12_patch3: 
#CMSSW_7_0_7_patch1:  process.GlobalTag.globaltag = "START70_V7A::All"

version = os.environ.get('CMSSW_VERSION')
if version.find('CMSSW_5_3') != -1:
    globaltag = "START53_V7G::All"
else:
    globaltag = "START70_V7A::All"

#globaltag = "START70_V7A::All"
print(globaltag)
process.GlobalTag.globaltag = globaltag

# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
         "/store/mc/Summer12_DR53X/BuToJPsiK_K2MuPtEtaEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v2/0000/521047B4-E0DD-E111-8146-02215E94E7DF.root"
 )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "test.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose        = cms.untracked.int32(1),
    printFrequency = cms.untracked.int32(100),
    requireCand    =  cms.untracked.bool(True),
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
    process.recoStuffSequence*
    process.bmmSequence*
    process.truthBmmSequence*
    process.tree
)
