import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from GeneratorInterface.EvtGenInterface.EvtGenSetting_cff import *

generator = cms.EDFilter(
    "Pythia8GeneratorFilter", 
    pythiaPylistVerbosity = cms.untracked.int32(0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.double(13000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    crossSection = cms.untracked.double(540000000.),
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            operates_on_particles = cms.vint32( 0 ), # 0 (zero) means default list (hardcoded)
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2010.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt.pdl'),
        ),
        parameterSets = cms.vstring('EvtGen130')
    ),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring("SoftQCD:nonDiffractive = on"),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CUEP8M1Settings',
            'processParameters',
        )
    )
)

generator.PythiaParameters.processParameters.extend(EvtGenExtraParticles)

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: Configuration/Generator/python/PYTHIA8_InclMu_nonDiffractive_CUEP8M1_13TeV_cff.py $'),
    annotation = cms.untracked.string('Pythia8+EvtGen130 generation, non-diffractive, pT(mu) > 3GeV, 13TeV, Tune CUEP8M1')
)

# -- Filters
muFilter = cms.EDFilter(
    "PythiaFilter",
    MinPt = cms.untracked.double(3.),
    MaxPt = cms.untracked.double(9999.),
    MaxEta = cms.untracked.double(3.),
    MinEta = cms.untracked.double(-3.),
    ParticleID = cms.untracked.int32(13)
)

ProductionFilterSequence = cms.Sequence(generator*muFilter)



