import FWCore.ParameterSet.Config as cms
from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *

generator = cms.EDFilter("Pythia8GeneratorFilter", 
                         pythiaPylistVerbosity = cms.untracked.int32(0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13000.0),
                         maxEventsToPrint = cms.untracked.int32(0),
                         filterEfficiency = cms.untracked.double(1.0),
                         crossSection = cms.untracked.double(540000000.),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        processParameters = cms.vstring("SoftQCD:nonDiffractive = on"),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'processParameters',
                                    )
        )
                         )

configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    name = cms.untracked.string('$Source: Configuration/Generator/python/PYTHIA8_InclB2Mu_CUEP8M1_13TeV_cff.py $'),
    annotation = cms.untracked.string('Spring 2015: Pythia8 generation enriched with b -> mu, 13TeV, Tune CUEP8M1')
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
    MinPt = cms.untracked.double(3.0),
    ParticleID = cms.untracked.int32(13)
    )

ProductionFilterSequence = cms.Sequence(generator*cFilter*muFilter)



