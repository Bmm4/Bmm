import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from GeneratorInterface.EvtGenInterface.EvtGenSetting_cff import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         maxEventsToPrint = cms.untracked.int32(0),
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.38e-3),
                         crossSection = cms.untracked.double(567900000.),
                         comEnergy = cms.double(13000.0),

                         ExternalDecays = cms.PSet(
        EvtGen130 = cms.untracked.PSet(
            decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2010.DEC'),
            particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt_Bmm.pdl'),
            user_decay_file = cms.vstring('GeneratorInterface/ExternalDecays/data/Bd_pimunu.dec'),
            list_forced_decays = cms.vstring('MyB0',
                                             'Myanti-B0'),
            operates_on_particles = cms.vint32(),
            ),
        parameterSets = cms.vstring('EvtGen130')
        ),
                         
                         PythiaParameters = cms.PSet(pythia8CommonSettingsBlock,
                                                     pythia8CUEP8M1SettingsBlock,
                                                     processParameters = cms.vstring("SoftQCD:nonDiffractive = on"),
                                                     parameterSets = cms.vstring('pythia8CommonSettings',
                                                                                 'pythia8CUEP8M1Settings',
                                                                                 'processParameters',
                                                                                 )
                                                     )
                         )

generator.PythiaParameters.processParameters.extend(EvtGenExtraParticles)

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: Configuration/Generator/python/PYTHIA8_Bd2PiMuNu_EtaPtFilter_CUEP8M1_13TeV_cff.py $'),
    annotation = cms.untracked.string('Spring 2015: Pythia8+EvtGen130 generation of Bd --> Pi- Mu+ Nu, 13TeV, Tune CUETP8M1')
    )



bfilter = cms.EDFilter(
        "PythiaFilter",
        MaxEta = cms.untracked.double(9999.),
        MinEta = cms.untracked.double(-9999.),
        ParticleID = cms.untracked.int32(511)
        )

decayfilter = cms.EDFilter(
        "PythiaDauVFilter",
	verbose         = cms.untracked.int32(0), 
	NumberDaughters = cms.untracked.int32(3), 
	ParticleID      = cms.untracked.int32(511),  
        DaughterIDs     = cms.untracked.vint32(-211, -13, 14),
	MinPt           = cms.untracked.vdouble(-1., -1., -1.), 
	MinEta          = cms.untracked.vdouble(-9999., -9999., -9999.), 
	MaxEta          = cms.untracked.vdouble( 9999.,  9999., 9999.)
        )


ProductionFilterSequence = cms.Sequence(generator*bfilter*decayfilter)

