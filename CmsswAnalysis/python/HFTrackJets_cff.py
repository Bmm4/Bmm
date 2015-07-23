import FWCore.ParameterSet.Config as cms


#from RecoJets.JetProducers.TrackJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from  SimGeneral.HepPDTESSource.pythiapdt_cfi import *


trkjetflavourByRef = cms.EDProducer("JetPartonMatcher",
    #jets = cms.InputTag("sis4TrackJets"),
    jets = cms.InputTag("myak4TrackJets"),
    coneSizeToAssociate = cms.double(0.5),
    partons = cms.InputTag("myPartons") #myPartons defined in HFJets_cff.py
)

selectTracks = cms.EDFilter("TrackSelector",
    src = cms.InputTag("generalTracks"),
    cut = cms.string('pt > 0.3 & pt<500 & hitPattern.pixelLayersWithMeasurement > 1 & hitPattern.trackerLayersWithMeasurement > 4')
 )


trackCandidates = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("selectTracks"),
    particleType = cms.string('pi+')
)

alltrackCandidates = cms.EDProducer("ConcreteChargedCandidateProducer",
    src = cms.InputTag("generalTracks"),
    particleType = cms.string('pi+')
)

trackjetDump = cms.EDAnalyzer("HFDumpTrackJets",
    doflavortagging = cms.untracked.int32(0),
    verbose = cms.untracked.int32(0),
    jetsLabel = cms.untracked.string('myak4TrackJets'),
    jetsTagLabel  = cms.untracked.string('simpleSecondaryVertexHighPurBJetTags'),
    tracksLabel = cms.untracked.string('alltrackCandidates'),#for indices (need all tracks in the event, not only the ones used in the jet algorithm
    sourceByRefer = cms.InputTag("trkjetflavourByRef"),
    genparticlesLabel = cms.untracked.string('genParticles')
)



TrackJetParameters = cms.PSet(
    src            = cms.InputTag("trackCandidates"),
    srcPVs         = cms.InputTag('offlinePrimaryVerticesWithBS'),
    jetType        = cms.string('BasicJet'),
    jetPtMin       = cms.double(1.0),
    inputEMin      = cms.double(0.0),
    inputEtMin     = cms.double(0.0),
    doPVCorrection = cms.bool(False),
    # pileup with offset correction
    doPUOffsetCorr = cms.bool(False),
       # if pileup is false, these are not read:
       nSigmaPU = cms.double(1.0),
       radiusPU = cms.double(0.5),
    # fastjet-style pileup
    doPUFastjet    = cms.bool(False),
      # if doPU is false, these are not read:
      Active_Area_Repeats = cms.int32(5),
      GhostArea = cms.double(0.01),
      Ghost_EtaMax = cms.double(6.0),
      doAreaFastjet = cms.bool( False ),
      doRhoFastjet = cms.bool( False ),
    # only use the tracks that were used to fit the vertex
    UseOnlyVertexTracks = cms.bool(False),
    # only consider the highest-sum-pT PV for clustering
    UseOnlyOnePV        = cms.bool(True),
    # maximum z-distance between track and vertex for association (in cm)
    DzTrVtxMax          = cms.double(2),
    # maximum xy-distance between track and vertex for association (in cm)
    DxyTrVtxMax         = cms.double(10000),
    # new in 37x
    MinVtxNdof=cms.int32(4),  # default is 5,  358,36 had no cut
    MaxVtxZ=cms.double(15.)   # default is 15
)

myak4TrackJets = cms.EDProducer(
    "FastjetJetProducer",
    TrackJetParameters,
    AnomalousCellParameters,
    jetAlgorithm = cms.string("AntiKt"),
    rParam       = cms.double(0.4)
)




TrackJetDumpAOD = cms.Sequence(selectTracks*trackCandidates*alltrackCandidates*myak4TrackJets*trackjetDump)


