import FWCore.ParameterSet.Config as cms


candidateDump = cms.EDAnalyzer("HFDumpCandidate",
    muonLabel = cms.untracked.string('muons'),
    trackLabel = cms.untracked.string('generalTracks'),
    trackptmin = cms.untracked.double(1.0),
    deltaR  = cms.untracked.double(1.5),
    verbose = cms.untracked.int32(0)
)

CandidateDump = cms.Sequence(candidateDump)

 
