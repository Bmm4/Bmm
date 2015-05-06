import FWCore.ParameterSet.Config as cms

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from GeneratorInterface.EvtGenInterface.EvtGenSetting_cff import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         maxEventsToPrint = cms.untracked.int32(0),
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.38e-3),
                         crossSection = cms.untracked.double(540000000.),
                         comEnergy = cms.double(13000.0),

                         ExternalDecays = cms.PSet(
                         EvtGen130 = cms.untracked.PSet(
#                         EvtGen = cms.untracked.PSet(
                         operates_on_particles = cms.vint32( 0 ), # 0 (zero) means default list (hardcoded)
                                                                  # you can put here the list of particles (PDG IDs)
                                                                  # that you want decayed by EvtGen
#                         use_default_decay = cms.untracked.bool(False),
#                         decay_table = cms.FileInPath('GeneratorInterface/ExternalDecays/data/DECAY_NOLONGLIFE.DEC'),
                         decay_table = cms.string('GeneratorInterface/EvtGenInterface/data/DECAY_2010.DEC'),
                         particle_property_file = cms.FileInPath('GeneratorInterface/EvtGenInterface/data/evt.pdl'),
#                         particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt.pdl'),
#                         user_decay_file = cms.vstring('GeneratorInterface/EvtGenInterface/data/Bs_mumu.dec'),
#                         user_decay_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/Bs_mumu.dec'),
#                         list_forced_decays = cms.vstring('MyB_s0',
#                                                          'Myanti-B_s0'),
                         ),
#                         parameterSets = cms.vstring('EvtGen')
                         parameterSets = cms.vstring('EvtGen130')
                         ),

                         PythiaParameters = cms.PSet(pythia8CommonSettingsBlock,
                                                     pythia8CUEP8M1SettingsBlock,
#                                                     processParameters = cms.vstring('HardQCD:all = on'),
                                                     processParameters = cms.vstring("SoftQCD:nonDiffractive = on"),
#                                                     processParameters = cms.vstring(            
#                                                                                     'HardQCD:gg2bbbar = on ',
#                                                                                     'HardQCD:qqbar2bbbar = on ',
#                                                                                     'HardQCD:hardbbbar = on',
#                                                                                     'PhaseSpace:pTHatMin = 1.5',
#                                                                                     ),
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
    name = cms.untracked.string('$Source: Configuration/Generator/python/PYTHIA8_InclC2Mu_CUEP8M1_13TeV_cff.py $'),
    annotation = cms.untracked.string('Spring 2015: Pythia8+EvtGen130 generation enriched with c -> mu, 13TeV, Tune CUETP8M1')
    )

# -- Filters
cFilter = cms.EDFilter(
        "PythiaFilter",
        MaxEta = cms.untracked.double(9999.),
        MinEta = cms.untracked.double(-9999.),
        ParticleID = cms.untracked.int32(4)
        )

muFilter = cms.EDFilter(
        "PythiaFilter",
        MaxEta = cms.untracked.double(3.0),
        MinEta = cms.untracked.double(-3.0),
        MaxPt = cms.untracked.double(9999.),
        MinPt = cms.untracked.double(2.5),
        ParticleID = cms.untracked.int32(13)
        )

ProductionFilterSequence = cms.Sequence(generator*cFilter*muFilter)
   


