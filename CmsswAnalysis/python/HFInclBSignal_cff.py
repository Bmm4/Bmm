import FWCore.ParameterSet.Config as cms



HFInclBMuonTrackJets = cms.EDAnalyzer("HFInclBMuonTrackJets",
    usejetforip = cms.untracked.int32(1),
    muonLabel = cms.untracked.string('muons'),
#    jetsLabel = cms.untracked.string('ak4TrackJets'),
    jetsLabel = cms.untracked.string('myak5TrackJets'),
    tracksLabel = cms.untracked.string('alltrackCandidates'),
    vertexLabel = cms.untracked.string('offlinePrimaryVerticesWithBS'),
    simvertexLabel = cms.untracked.string('simG4'),
    jetmatch =       cms.untracked.double(0.5),
    jetetmin =        cms.untracked.double(1.0),
    verbose = cms.untracked.int32(0)
)

#HFInclBSignalSequence = cms.Sequence(HFInclBMuonTrackJets)

HFInclB = cms.EDAnalyzer("HFInclB",
    verbose            = cms.untracked.int32(0), 
    muonsLabel         = cms.untracked.InputTag("muons"),
    tracksLabel        = cms.untracked.InputTag('generalTracks'),
    PrimaryVertexLabel = cms.untracked.InputTag("offlinePrimaryVertices"),
    BeamSpotLabel      = cms.untracked.InputTag("offlineBeamSpot"),
    muonQualityString  = cms.untracked.string("AllGlobalMuons"), 
    muonPt             = cms.untracked.double(4.0),
    nMuons             = cms.untracked.int32(1),
    trackPt            = cms.untracked.double(1.0),
    deltaR             = cms.untracked.double(99.0),
    maxDoca            = cms.untracked.double(0.01),
    maxD0              = cms.untracked.double(99.0),
    maxDz              = cms.untracked.double(99.0),
    pvWeight           = cms.untracked.double(0.6),
    type               = cms.untracked.int32(1300)
)

HFInclBSignalSequence = cms.Sequence(HFInclB*HFInclBMuonTrackJets)

 
