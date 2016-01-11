import FWCore.ParameterSet.Config as cms

from Configuration.Generator.PythiaUEZ2starSettings_cfi import *

generator = cms.EDFilter(
    "Pythia6GeneratorFilter",
    comEnergy = cms.double(13000.0),
    crossSection = cms.untracked.double(540000000.),
    filterEfficiency = cms.untracked.double(1.38e-3),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    maxEventsToPrint = cms.untracked.int32(0),
    pythiaPylistVerbosity = cms.untracked.int32(0),
    ExternalDecays = cms.PSet(
        EvtGen = cms.untracked.PSet(
             operates_on_particles = cms.vint32( 0 ), # 0 (zero) means default list (hardcoded)
                                                      # you can put here the list of particles (PDG IDs)
                                                      # that you want decayed by EvtGen
             use_default_decay = cms.untracked.bool(False),
             decay_table = cms.FileInPath('GeneratorInterface/ExternalDecays/data/DECAY_NOLONGLIFE.DEC'),
             particle_property_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/evt.pdl'),
             user_decay_file = cms.FileInPath('GeneratorInterface/ExternalDecays/data/Bs_mumu.dec'),
             list_forced_decays = cms.vstring('MyB_s0',
                                              'Myanti-B_s0'),
        ),
        parameterSets = cms.vstring('EvtGen')
    ),

    
    PythiaParameters = cms.PSet(
    pythiaUESettingsBlock,
         bbbarSettings = cms.vstring('MSEL = 1'), 
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring(
             'pythiaUESettings',
             'bbbarSettings')
       
    )
    )
bfilter = cms.EDFilter(
        "PythiaFilter",
        MaxEta = cms.untracked.double(9999.),
        MinEta = cms.untracked.double(-9999.),
        ParticleID = cms.untracked.int32(531)
        )

decayfilter = cms.EDFilter(
        "PythiaDauVFilter",
	verbose         = cms.untracked.int32(1), 
	NumberDaughters = cms.untracked.int32(2), 
	ParticleID      = cms.untracked.int32(531),  
        DaughterIDs     = cms.untracked.vint32(13, -13),
	MinPt           = cms.untracked.vdouble(-99., -99.), 
	MinEta          = cms.untracked.vdouble(-9999., -9999.), 
	MaxEta          = cms.untracked.vdouble( 9999.,  9999.)
        )

#configurationMetadata = cms.untracked.PSet(
#    version = cms.untracked.string('$Revision: 1.1 $'),
#    name = cms.untracked.string
#    ('$Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/GenProduction/python/EightTeV/PYTHIA6_Bs2MuMu_EtaPtFilter_TuneZ2star_13TeV_cff.py,v $'),
#    annotation = cms.untracked.string('Bs -> mu mu at 13TeV')
#    )

ProductionFilterSequence = cms.Sequence(generator*bfilter*decayfilter)


