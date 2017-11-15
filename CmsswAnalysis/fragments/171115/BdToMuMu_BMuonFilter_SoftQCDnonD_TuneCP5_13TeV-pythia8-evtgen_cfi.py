import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         maxEventsToPrint = cms.untracked.int32(0),
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.38e-3),
                         crossSection = cms.untracked.double(540000000.),
                         comEnergy = cms.double(13000.0),
                         ExternalDecays = cms.PSet(
                             EvtGen130 = cms.untracked.PSet(
                                 decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
                                 particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt_2014.pdl'),
                                 user_decay_file = cms.vstring('GeneratorInterface/ExternalDecays/data/Bd_mumu.dec'),
                                 list_forced_decays = cms.vstring('MyB0', 'Myanti-B0'),
                                 operates_on_particles = cms.vint32(),
                                 convertPythiaCodes = cms.untracked.bool(False),
                             ),
                             parameterSets = cms.vstring('EvtGen130')
                         ),

                         PythiaParameters = cms.PSet(pythia8CommonSettingsBlock,
                                                     pythia8CP5SettingsBlock,
                                                     processParameters = cms.vstring("SoftQCD:nonDiffractive = on"),
                                                     parameterSets = cms.vstring('pythia8CommonSettings',
                                                                                 'pythia8CP5Settings',
                                                                                 'processParameters',
                                                     )
                         )
)

bfilter = cms.EDFilter(
    "PythiaFilter",
    MaxEta = cms.untracked.double(9999.),
    MinEta = cms.untracked.double(-9999.),
    ParticleID = cms.untracked.int32(511)
)

decayfilter = cms.EDFilter(
    "PythiaDauVFilter",
    verbose         = cms.untracked.int32(1),
    NumberDaughters = cms.untracked.int32(2),
    ParticleID      = cms.untracked.int32(511),
    DaughterIDs     = cms.untracked.vint32(13, -13),
    MinPt           = cms.untracked.vdouble(3.5, 3.5),
    MinEta          = cms.untracked.vdouble(-2.5, -2.5),
    MaxEta          = cms.untracked.vdouble( 2.5,  2.5)
)

ProductionFilterSequence = cms.Sequence(generator*bfilter*decayfilter)
