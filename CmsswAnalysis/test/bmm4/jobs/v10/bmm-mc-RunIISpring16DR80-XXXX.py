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

process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v8', '')

# ----------------------------------------------------------------------
# POOLSOURCE
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())



# ----------------------------------------------------------------------
rootFileName = "bmm-mc-RunIISpring16DR80-XXXX.root"

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


process.bmmDump.useBeamspotConstraint = cms.untracked.bool(False)
process.bupsikpDump.useBeamspotConstraint = cms.untracked.bool(False)
process.bspsiphiDump.useBeamspotConstraint = cms.untracked.bool(False)
process.bdpsikstarDump.useBeamspotConstraint = cms.untracked.bool(False)
process.bspsif0Dump.useBeamspotConstraint = cms.untracked.bool(False)

process.truthBsToMuMuDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBdToMuMuDump.useBeamspotConstraint = cms.untracked.bool(False)

process.truthBsToMuMuGaDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBsToKKDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBsToKPiDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBsToPiPiDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBsToPiMuNuDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBsToKMuNuDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBsToPhiMuMuDump.useBeamspotConstraint = cms.untracked.bool(False)

process.truthBdToPiPiDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBdToKPiDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBdToKKDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBdToMuMuPi0Dump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBdToPiMuNuDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBdToMuMuGaDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBdToMuMuK0Dump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBdToRhoPiDump.useBeamspotConstraint = cms.untracked.bool(False)


process.truthBcToJpsiMuNuDump.useBeamspotConstraint = cms.untracked.bool(False)

process.truthBuTo3MuNuDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBuToMuMuKDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBuToMuMuPiDump.useBeamspotConstraint = cms.untracked.bool(False)


process.truthLambdaBToPPiDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthLambdaBToPKDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthLambdaBToPMuNuDump.useBeamspotConstraint = cms.untracked.bool(False)

process.truthBd2DstarPiDump.useBeamspotConstraint = cms.untracked.bool(False)

process.truthBs2JpsiPhiDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBs2Jpsif0Dump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBs2JpsiPhiAsBpDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBs2JpsiPhiAsBdDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBs2JpsiPiPiAsBsDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBd2JpsiKsDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBd2JpsiKstarDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBd2JpsiKstarAsBpDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBd2JpsiKstarAsBsDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBu2JpsiKpDump.useBeamspotConstraint = cms.untracked.bool(False)
process.truthBu2JpsiPiDump.useBeamspotConstraint = cms.untracked.bool(False)

process.truthPsiToMuMu.useBeamspotConstraint = cms.untracked.bool(False)
process.truthPsi2SToMuMu.useBeamspotConstraint = cms.untracked.bool(False)
process.truthUps1SToMuMu.useBeamspotConstraint = cms.untracked.bool(False)
process.truthUps2SToMuMu.useBeamspotConstraint = cms.untracked.bool(False)
process.truthUps3SToMuMu.useBeamspotConstraint = cms.untracked.bool(False)

process.truthD0ToKPi.useBeamspotConstraint = cms.untracked.bool(False)
process.truthDpToKPiPi.useBeamspotConstraint = cms.untracked.bool(False)
process.truthDpToKstarPi.useBeamspotConstraint = cms.untracked.bool(False)
process.truthDsToPhiPi.useBeamspotConstraint = cms.untracked.bool(False)
process.truthDstarToD0PiToKPiPi.useBeamspotConstraint = cms.untracked.bool(False)
process.truthDpToKKPi.useBeamspotConstraint = cms.untracked.bool(False)
process.truthLambdaCToPrKPi.useBeamspotConstraint = cms.untracked.bool(False)

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.genDump*
    process.recoStuffSequence*
    process.bmmSequence*
    process.truthBmmSequence*
    process.tree
)
