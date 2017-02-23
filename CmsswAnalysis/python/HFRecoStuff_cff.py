import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
stuffDump = cms.EDAnalyzer(
    "HFDumpStuff",
    verbose                  = cms.untracked.int32(0),
    PrimaryVertexLabel       = cms.untracked.InputTag("offlinePrimaryVertices"),
    PrimaryVertexTracksLabel = cms.untracked.InputTag("generalTracks"),
    type                     = cms.untracked.int32(0)
    )


# ----------------------------------------------------------------------
from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import *
#from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *

trkDump = cms.EDAnalyzer(
    "HFDumpTracks",
    verbose                        = cms.untracked.int32(0),
    tracksLabel                    = cms.untracked.InputTag('generalTracks'),
    primaryVertexCollectionLabel   = cms.untracked.InputTag('offlinePrimaryVertices'),
    muonsLabel                     = cms.untracked.InputTag("muons"),
    doTruthMatching                = cms.untracked.int32(3),
    dumpSimpleTracks               = cms.untracked.bool(True),
    dumpRecTracks                  = cms.untracked.bool(False),
    type                           = cms.untracked.int32(0)
)

# ----------------------------------------------------------------------
muonDump = cms.EDAnalyzer(
    "HFDumpMuons",
    verbose                        = cms.untracked.int32(0),
    primaryVertexCollectionLabel   = cms.untracked.InputTag('offlinePrimaryVertices'),
    tracksLabel                    = cms.untracked.InputTag("generalTracks"),
    muonsLabel                     = cms.untracked.InputTag("muons"),
    calomuonsLabel                 = cms.untracked.InputTag("calomuons"),
    maxTrackDist                   = cms.untracked.double(0.1),
    docaVertex                     = cms.untracked.double(0.05),
    keepBest                       = cms.untracked.int32(3),
    maxCandTracks                  = cms.untracked.int32(3),
    type                           = cms.untracked.int32(0),
    # Configuration for the extrapolation at the muon system
    propM1 = cms.PSet(
        useStation2 = cms.bool(False),
        useTrack = cms.string("tracker"),
        useState = cms.string("atVertex"),  # in AOD; this is for matching to other tracks?
        useSimpleGeometry = cms.bool(True), # use just one cylinder and two planes, not all the fancy chambers
    ),
    propM2 = cms.PSet(
        useStation2 = cms.bool(True),
        useTrack = cms.string("tracker"),
        useState = cms.string("outermost"), # in AOD; this is for matching to the L1 muon!
        useSimpleGeometry = cms.bool(True), # use just one cylinder and two planes, not all the fancy chambers
    ),
    #propagate the global track from vertex outwards
    OutwardPropM1 = cms.PSet(
        useStation2 = cms.bool(False),
        useTrack = cms.string("global"),
        useState = cms.string("atVertex"),
        useSimpleGeometry = cms.bool(True)
        ),
    #propagate the stand alone muon track inwards
    InwardPropM1 = cms.PSet(
        useStation2 = cms.bool(False),
        useTrack = cms.string("muon"),
        useState = cms.string("outermost"),
        useSimpleGeometry = cms.bool(True)
        ),
    weightFileBarrel               = cms.untracked.FileInPath("Bmm/RootAnalysis/macros/weights/TMVA-muonid-bmm4-B-19.weights.xml"),
    weightFileEndcap               = cms.untracked.FileInPath("Bmm/RootAnalysis/macros/weights/TMVA-muonid-bmm4-E-19.weights.xml"),
)

# ----------------------------------------------------------------------
triggerDump = cms.EDAnalyzer(
    "HFDumpTrigger",
    verbose                 = cms.untracked.int32(0),
    HLTProcessName          = cms.untracked.string('HLT'),
    L1GTReadoutRecordLabel  = cms.untracked.InputTag("gtDigis"),
    hltL1GtObjectMap        = cms.untracked.InputTag("hltL1GtObjectMap"),
    L1MuonsLabel            = cms.untracked.InputTag("hltL1extraParticles"),
    HLTResultsLabel         = cms.untracked.InputTag("TriggerResults::HLT"),
    TriggerEventLabel       = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
    )


# ----------------------------------------------------------------------
hltrep = cms.EDAnalyzer(
    "HLTrigReport",
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT")
    )


l1trep = cms.EDAnalyzer(
    "L1GtTrigReport",
    UseL1GlobalTriggerRecord = cms.bool( False ),
    L1GtRecordInputTag = cms.InputTag( "gtDigis" )
    )


# ######################################################################
# Sequences
# ######################################################################
#recoStuffSequence     = cms.Sequence(stuffDump*trkDump*muonDump*triggerDump*hltrep*l1trep)
#recoStuffSequence     = cms.Sequence(stuffDump*trkDump*muonDump*triggerDump)
recoStuffSequence     = cms.Sequence(stuffDump*trkDump*muonDump*triggerDump)
