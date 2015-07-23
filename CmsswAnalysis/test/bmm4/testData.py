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

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_mc', '')


# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
#         "/store/user/ursl/files/Summer12_DR53X/BuToJPsiK_K2MuPtEtaEtaFilter_8TeV-pythia6-evtgen-0000.root"
#         "/store/mc/Summer12_DR53X/BuToJPsiK_K2MuPtEtaEtaFilter_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7A-v2/0000/521047B4-E0DD-E111-8146-02215E94E7DF.root"
#	 "/store/relval/CMSSW_7_3_0_pre2/RelValBsToMuMu_13/GEN-SIM-RECO/MCRUN2_73_V1-v1/00000/7CA64471-C16B-E411-BB89-02163E00F463.root"
#	"/store/data/Run2012A/MuOnia/AOD/22Jan2013-v1/30000/000D2FF5-EE82-E211-BEBA-0026189438A5.root"
#    "/store/data/Run2012A/MuOnia/RECO/22Jan2013-v1/30000/0002CEE5-4183-E211-84FB-00304867915A.root"
        "/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/562/00000/066732D9-752A-E511-B293-02163E0121CC.root",
        "/store/data/Run2015B/Charmonium/AOD/PromptReco-v1/000/251/562/00000/2C938991-692A-E511-9F4A-02163E014761.root"
 ),
 skipEvents = cms.untracked.uint32(704)
)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )

process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "output.root"

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
process.load("Bmm.CmsswAnalysis.HFLambdas_cff")

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
#    process.genDump*
    process.recoStuffSequence*
    process.bmmSequence*
#    process.HFLambdasSequence*
    process.truthBmmSequence*
    process.tree
)
