import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from GeneratorInterface.EvtGenInterface.EvtGenSetting_cff import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         maxEventsToPrint = cms.untracked.int32(0),
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.38e-3),
                         crossSection = cms.untracked.double(568000000.),
                         comEnergy = cms.double(13000.0),

                         ExternalDecays = cms.PSet(
                         EvtGen130 = cms.untracked.PSet(
                         operates_on_particles = cms.vint32(), # 0 (zero) means default list (hardcoded)
                                                                  # you can put here the list of particles (PDG IDs)
                                                                  # that you want decayed by EvtGen
                         decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2010.DEC'),
                         particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt_Bmm.pdl'),
                         user_decay_file = cms.vstring('GeneratorInterface/ExternalDecays/data/Bd_mumu.dec'),
                         list_forced_decays = cms.vstring('MyB0',
                                                          'Myanti-B0'),
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

#new line
generator.PythiaParameters.processParameters.extend(EvtGenExtraParticles)

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: Configuration/Generator/python/PYTHIA8_Bd2MuMu_PtEtaFilter_CUEP8M1_13TeV_cff.py $'),
    annotation = cms.untracked.string('Spring 2015: Pythia8+EvtGen130 generation of Bd --> Mu+Mu-, 13TeV, Tune CUETP8M1')
    )


# Filters

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
	MinPt           = cms.untracked.vdouble(2.5, 2.5), 
	MinEta          = cms.untracked.vdouble(-2.5, -2.5), 
	MaxEta          = cms.untracked.vdouble( 2.5,  2.5)
        )


ProductionFilterSequence = cms.Sequence(generator*bfilter*decayfilter)
   


