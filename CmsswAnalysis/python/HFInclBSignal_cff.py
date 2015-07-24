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

HFInclBSignalSequence = cms.Sequence(HFInclBMuonTrackJets)

 
