import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *

generator = cms.EDFilter(
    "Pythia8GeneratorFilter",
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(2.040e-01),
    crossSection = cms.untracked.double(540000000.),
    comEnergy = cms.double(13000.0),
    ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2014_NOLONGLIFE.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_2014.pdl'),
            user_decay_file = cms.vstring('GeneratorInterface/EvtGenInterface/data/Bs_Jpsif0_mumupipi.dec'),
            list_forced_decays = cms.vstring('MyB_s0', 'Myanti-B_s0'),
            operates_on_particles = cms.vint32(),
            convertPythiaCodes = cms.untracked.bool(False),
        ),
        parameterSets = cms.vstring('EvtGen130')
    ),
    PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CP5SettingsBlock,
        processParameters = cms.vstring(
            "SoftQCD:nonDiffractive = on",
            'PTFilter:filter = on', # this turn on the filter
            'PTFilter:quarkToFilter = 5', # PDG id of q quark
            'PTFilter:scaleToFilter = 1.0'),
        parameterSets = cms.vstring(
            'pythia8CommonSettings',
            'pythia8CP5Settings',
            'processParameters',
        )
    )
)

bfilter = cms.EDFilter(
    "PythiaFilter",
    MaxEta = cms.untracked.double(9999.),
    MinEta = cms.untracked.double(-9999.),
    ParticleID = cms.untracked.int32(531)
)

jpsifilter = cms.EDFilter(
    "PythiaDauVFilter",
    MotherID = cms.untracked.int32(531),
    ParticleID = cms.untracked.int32(443),
    NumberDaughters = cms.untracked.int32(2),
    DaughterIDs = cms.untracked.vint32(13, -13),
    MinPt = cms.untracked.vdouble(-99., -99.),
    MinEta = cms.untracked.vdouble(-9999., -9999.),
    MaxEta = cms.untracked.vdouble(9999., 9999.),
    verbose = cms.untracked.int32(1)
)

f0filter = cms.EDFilter(
    "PythiaDauVFilter",
    MotherID = cms.untracked.int32(531),
    ParticleID = cms.untracked.int32(9010221),
    NumberDaughters = cms.untracked.int32(2),
    DaughterIDs = cms.untracked.vint32(211, -211),
    MinPt = cms.untracked.vdouble(-99., -99.),
    MaxEta = cms.untracked.vdouble(9999., 9999.),
    MinEta = cms.untracked.vdouble(-9999., -9999.),
    verbose = cms.untracked.int32(1)
)


ProductionFilterSequence = cms.Sequence(generator*bfilter*jpsifilter*f0filter)
