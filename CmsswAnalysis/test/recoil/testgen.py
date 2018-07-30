import os
import sys
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
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


# ----------------------------------------------------------------------
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)
cmsFile  = "file:%s" % sys.argv[2]
rootFile = "file:%s" % sys.argv[3]
print cmsFile
print rootFile


process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(cmsFile)
    )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose        = cms.untracked.int32(0),
    printFrequency = cms.untracked.int32(1),
    requireCand    = cms.untracked.bool(False),
    fullGenBlock   = cms.untracked.bool(True),
    fileName       = cms.untracked.string(rootFile)
    )

# ----------------------------------------------------------------------
process.stuffDump = cms.EDAnalyzer(
    "HFDumpStuff",
    verbose                  = cms.untracked.int32(0),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    PrimaryVertexTracksLabel = cms.untracked.InputTag("generalTracks")
    )

# ----------------------------------------------------------------------
process.load("Bmm.CmsswAnalysis.HFMCTruth_cff")
process.genDump.verbose = cms.untracked.int32(0)

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.stuffDump*
    process.MCTruthSequence*
    process.tree
)
