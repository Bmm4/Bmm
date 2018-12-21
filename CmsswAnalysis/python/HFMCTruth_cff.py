import FWCore.ParameterSet.Config as cms

# ----------------------------------------------------------------------
HepPDTESSource = cms.ESSource(
    "HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/particle.tbl')
    )

# ----------------------------------------------------------------------
genParticles = cms.EDProducer(
    "GenParticleProducer",
    saveBarCodes          = cms.untracked.bool(True),
    src                   = cms.InputTag("generatorSmeared"),
#    src                   = cms.InputTag("generator"),
#    src                   = cms.InputTag("source"),
    abortOnUnknownPDGCode = cms.untracked.bool(False)
    )


# ----------------------------------------------------------------------
genDump = cms.EDAnalyzer(
    "HFDumpGenerator",
    generatorCandidates = cms.untracked.string('genParticles'),
    generatorEvent      = cms.untracked.string('generator'),
    verbose             = cms.untracked.int32(0)
    )

# ######################################################################
# Sequences
# ######################################################################
MCTruthSequence     = cms.Sequence(genParticles*genDump)
