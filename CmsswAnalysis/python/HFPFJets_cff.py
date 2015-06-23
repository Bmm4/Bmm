import FWCore.ParameterSet.Config as cms




pfjflavourByRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5PFJets"),
    coneSizeToAssociate = cms.double(0.5),
    partons = cms.InputTag("myPartons")
)

pfjetDump = cms.EDAnalyzer("HFDumpPFJets",
    jetsLabel = cms.untracked.string('ak5PFJets'),
#    jetsTagLabel  = cms.untracked.string('simpleSecondaryVertexBJetTags'),
#    jetsTagLabel  = cms.untracked.string('trackProbabilityJetTags'),
    jetsTagLabel  = cms.untracked.string('simpleSecondaryVertexHighPurBJetTags'),
    sourceByRefer = cms.InputTag("pfjflavourByRef"),
    verbose = cms.untracked.int32(0),
    doflavortagging = cms.untracked.int32(0),
    genparticlesLabel = cms.untracked.string('genParticles')
)


PFJetDump = cms.Sequence(pfjflavourByRef*pfjetDump)
PFJetDumpAOD = cms.Sequence(pfjetDump)
