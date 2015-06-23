import FWCore.ParameterSet.Config as cms


candidateDumpNew = cms.EDAnalyzer("HFDumpCandidateNew",
    muonLabel = cms.untracked.string('muons'),
    trackLabel = cms.untracked.string('generalTracks'),
    trackptmin = cms.untracked.double(1.0),
    muonptmin = cms.untracked.double(3.0),
    muoneta = cms.untracked.double(2.1),
    deltaR  = cms.untracked.double(1),
    verbose = cms.untracked.int32(0)
)

CandidateDumpNew = cms.Sequence(candidateDumpNew)

 
