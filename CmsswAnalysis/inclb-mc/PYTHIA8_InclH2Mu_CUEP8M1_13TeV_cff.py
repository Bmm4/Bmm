import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter",
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         maxEventsToPrint = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.0),
                         comEnergy = cms.double(13000.0),
                         crossSection = cms.untracked.double(540000000.),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring(
            'HardQCD:all = on',
            'ParticleDecays:limitCylinder   = true',
            'ParticleDecays:xyMax = 1500',
            'ParticleDecays:zMax = 3000',
            '321:onMode = on', #kaons
            '211:onMode = on', #pion
            '130:onMode = on' #KL0
            ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'processParameters',
                                    )
        )
                         )


configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: Configuration/Generator/python/PYTHIA8_InclH2Mu_CUEP8M1_13TeV_cff.py $'),
    annotation = cms.untracked.string('Spring 2015: Pythia8 generation enriched with h -> mu, 13TeV, Tune CUEP8M1')
    )

mugenfilter = cms.EDFilter("MCSmartSingleParticleFilter",
                           MinPt = cms.untracked.vdouble(3.0,3.0),
                           MinEta = cms.untracked.vdouble(-3.0,-3.0),
                           MaxEta = cms.untracked.vdouble(3.0,3.0),
                           ParticleID = cms.untracked.vint32(13,-13),
                           Status = cms.untracked.vint32(1,1),
                           # Decay cuts are in mm
                           MaxDecayRadius = cms.untracked.vdouble(1500.,1500.),
                           MinDecayZ = cms.untracked.vdouble(-3000.,-3000.),
                           MaxDecayZ = cms.untracked.vdouble(3000.,3000.)
                           )

ProductionFilterSequence = cms.Sequence(generator*mugenfilter)
