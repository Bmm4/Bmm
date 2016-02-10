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

# requires >= CMSSW_7_4_6
process.GlobalTag.globaltag = "76X_dataRun2_v15"

# ----------------------------------------------------------------------
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        #        "file:pickEvent-258694-57654257.root"
        # "/store/user/ursl/files/Summer12_DR53X/BsToMuMu_EtaPtFilter_8TeV-pythia6-evtgen-0000.root"
        #"/store/data/Run2015D/Charmonium/AOD/PromptReco-v4/000/258/694/00000/0A4E5B2B-8D70-E511-BDA7-02163E01233E.root"
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-256676-1062887410.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-256675-55481373.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-256843-433406866.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-256843-1150323527.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-257400-2027967968.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-256868-726929360.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-256867-144900350.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-257400-841821568.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-257613-403905837.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-257816-399246850.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-257645-1091133716.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-257816-422098776.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-257969-235736000.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-257822-1385587789.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-258694-57654257.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-258403-76137916.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-258158-1475541228.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-258712-9242142.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-258706-37127791.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-258706-1058311978.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-258703-202382186.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-259686-586255673.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-259683-151558946.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-259626-74695669.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-260424-1200254056.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-260431-667345014.root",
        "/store/user/ursl/files/daneks-D0-skim0/pickEvent-260426-15692016.root"
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )


# ----------------------------------------------------------------------
rootFileName = "T1-test.root"

process.tree = cms.EDAnalyzer(
    "HFTree",
    verbose        = cms.untracked.int32(0),
    printFrequency = cms.untracked.int32(100),
    requireCand    =  cms.untracked.bool(False),
    fullGenBlock   = cms.untracked.bool(False),
    fileName       = cms.untracked.string(rootFileName)
    )


process.debug = cms.EDAnalyzer(
    "HFDebug",
    verbose                 = cms.untracked.int32(100),
    tracksLabel             = cms.untracked.InputTag('generalTracks'),
    HLTProcessName          = cms.untracked.string('HLT'), 
    L1GTReadoutRecordLabel  = cms.untracked.InputTag("gtDigis"), 
    hltL1GtObjectMap        = cms.untracked.InputTag("hltL1GtObjectMap"), 
    L1MuonsLabel            = cms.untracked.InputTag("hltL1extraParticles"), 
    HLTResultsLabel         = cms.untracked.InputTag("TriggerResults::HLT"), 
    TriggerEventLabel       = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"), 
    )


# ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Bmm.CmsswAnalysis.HFRecoStuff_cff")
process.load("Bmm.CmsswAnalysis.HFBmm_cff")
process.load("Bmm.CmsswAnalysis.HFOnia_cff")
process.load("Bmm.CmsswAnalysis.HFPhysicsDeclared_cff")


process.psiDump.verbose   = cms.untracked.int32(100)
process.psiDump.massLow   = cms.untracked.double(3.02)
process.psiDump.massHigh  = cms.untracked.double(3.17)

# ----------------------------------------------------------------------
process.p = cms.Path(
    process.recoStuffSequence*
#    process.psiDump*
    process.debug*
    process.tree
)
