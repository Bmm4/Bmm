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

version = os.environ.get('CMSSW_VERSION')
if version.find('CMSSW_5_3') != -1:
    globaltag = "START53_V7G::All"
else:
    globaltag = "START72_V1::All"

#globaltag = "START70_V7A::All"
print('==> GT = %s ' % globaltag)
process.GlobalTag.globaltag = globaltag

# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
         "file:./PYTHIA8_Bs2MuMu_NoFilter_CUEP8M1_13TeV_cff_py_GEN.root"
 )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

# ----------------------------------------------------------------------
rootFileName = "gen.root"
process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose        = cms.untracked.int32(10),
    printFrequency = cms.untracked.int32(1),
    requireCand    = cms.untracked.bool(False),
    fullGenBlock   = cms.untracked.bool(True),
    fileName       = cms.untracked.string(rootFileName)
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
process.genDump.verbose = cms.untracked.int32(3)

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.stuffDump*
    process.MCTruthSequence*
    process.tree
)
