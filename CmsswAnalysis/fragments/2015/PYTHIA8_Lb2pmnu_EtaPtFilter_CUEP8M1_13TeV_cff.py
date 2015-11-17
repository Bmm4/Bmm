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
                         user_decay_file = cms.vstring('GeneratorInterface/ExternalDecays/data/LambdaB_pmnuLCSR.dec'),
                         list_forced_decays = cms.vstring('MyLambda_b0',
                                                          'Myanti-Lambda_b0'),
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
    name = cms.untracked.string('$Source: Configuration/Generator/python/PYTHIA8_Lb2pmnu_NoFilter_CUEP8M1_13TeV_cff.py $'),
    annotation = cms.untracked.string('Spring 2015: Pythia8+EvtGen130 generation of Lb --> p+ mu- nu, 13TeV, Tune CUETP8M1')
    )
# Filters

bfilter = cms.EDFilter(
        "PythiaFilter",
        MaxEta = cms.untracked.double(9999.),
        MinEta = cms.untracked.double(-9999.),
        ParticleID = cms.untracked.int32(5122)
        )

decayfilter = cms.EDFilter(
        "PythiaDauVFilter",
        verbose         = cms.untracked.int32(0), 
        NumberDaughters = cms.untracked.int32(3), 
        ParticleID      = cms.untracked.int32(5122),  
        DaughterIDs     = cms.untracked.vint32(2212, 13, -14),
        MinPt           = cms.untracked.vdouble(3.5, 3.5, -1.), 
        MinEta          = cms.untracked.vdouble(-2.5, -2.5, -9999.), 
        MaxEta          = cms.untracked.vdouble( 2.5,  2.5, 9999.)
        )


ProductionFilterSequence = cms.Sequence(generator*bfilter*decayfilter)

