# ######################################################################
# t3ui01
# /mnt/t3nfs01/data01/shome/ursl/bmm4/CMSSW_8_0_29/src/Bmm/CmsswAnalysis/test/bmm4/jobs/v11/bmmCharmonium2016H
# file list contains 12881 events
# mkPyFiles -t ../bmm-legacy-Run2016-XXXX.py -f ../../../catalogs/Run2016__Charmonium__07Aug17/Charmonium__Run2016H-07Aug17 -s v11 -n 1
# ./bmm-legacy-Run2016-Charmonium__Run2016H-07Aug17-v11-0000.py with 12881 events
# ######################################################################
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("HFA")

# ----------------------------------------------------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.categories.append('HLTrigReport')
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

# requires >= CMSSW_8_0_29
process.GlobalTag.globaltag = "80X_dataRun2_2016LegacyRepro_v4"

# ----------------------------------------------------------------------
process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
         "/store/data/Run2016H/Charmonium/AOD/07Aug17-v1/10000/0003E4A7-F19A-E711-963D-1866DAEA8178.root"
  )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )



# ----------------------------------------------------------------------
rootFileName = "bmm-legacy-Run2016-Charmonium__Run2016H-07Aug17-v11-0000.root"

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
process.load("Bmm.CmsswAnalysis.HFBmm_cff")
process.load("Bmm.CmsswAnalysis.HFOnia_cff")
process.load("Bmm.CmsswAnalysis.HFHadronic_cff")
process.load("Bmm.CmsswAnalysis.HFPhysicsDeclared_cff")

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.recoStuffSequence*
    process.bmmSequence*
    process.psiDump*
    process.upsDump*
    process.dstarDump*
    process.ksDump*
    process.phiDump*
    process.lambdaDump*
    process.tree
)
