import FWCore.ParameterSet.Config as cms


from PhysicsTools.JetMCAlgos.SelectPartons_cff import *

myPartons = cms.EDProducer("PartonSelector",
    withLeptons = cms.bool(False)
)


calojflavourByRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5CaloJets"),
    coneSizeToAssociate = cms.double(0.5),
    partons = cms.InputTag("myPartons")
)

calojetDump = cms.EDAnalyzer("HFDumpJets",
    jetsLabel = cms.untracked.string('ak5CaloJets'),
#    jetsTagLabel  = cms.untracked.string('simpleSecondaryVertexBJetTags'),
#    jetsTagLabel  = cms.untracked.string('trackProbabilityJetTags'),
    jetsTagLabel  = cms.untracked.string('simpleSecondaryVertexHighPurBJetTags'),
    sourceByRefer = cms.InputTag("calojflavourByRef"),
    verbose = cms.untracked.int32(0),
    doflavortagging = cms.untracked.int32(0),
    genparticlesLabel = cms.untracked.string('genParticles')
)


#CaloJetDump = cms.Sequence(myPartons*calojflavourByRef*calojetDump)

CaloJetDump = cms.Sequence(calojflavourByRef*calojetDump)
CaloJetDumpAOD = cms.Sequence(calojetDump)
