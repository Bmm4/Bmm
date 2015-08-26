# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/PYTHIA8_InclB2Mu_NOEVTGEN_CUEP8M1_13TeV_cfi.py --conditions auto:run2_mc -n 10000 --eventcontent FEVTDEBUG -s GEN,SIM --datatier GEN-SIM --beamspot NominalCollision2015 --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --magField 38T_PostLS1 --fileout file:step1.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
	input = cms.untracked.int32(500000)
	)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

	)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('$Revision: 1.1 $'),
	annotation = cms.untracked.string('Spring 2015: Pythia8 generation enriched with b -> mu, 13TeV, Tune CUEP8M1'),
	name = cms.untracked.string('$Source: Configuration/Generator/python/PYTHIA8_InclB2Mu_CUEP8M1_13TeV_cff.py $')
	)

# -- Reduced verbosity
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Random number
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
						   generator = cms.PSet(initialSeed = cms.untracked.uint32(10000)),
						   VtxSmeared = cms.PSet(initialSeed = cms.untracked.uint32(10000)),
						   g4SimHits =  cms.PSet(initialSeed = cms.untracked.uint32(10000))
						   )

# Output definition

process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
					   splitLevel = cms.untracked.int32(0),
					   eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
					   outputCommands = process.FEVTDEBUGEventContent.outputCommands,
					   fileName = cms.untracked.string('file:PYTHIA8_InclC2Mu_NOEVTGEN_CUEP8M1_13TeV_GEN_SIM-20000.root'),
					   dataset = cms.untracked.PSet(
		filterName = cms.untracked.string(''),
		dataTier = cms.untracked.string('GEN-SIM')
		),
					   SelectEvents = cms.untracked.PSet(
		SelectEvents = cms.vstring('generation_step')
		)
					   )

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.muFilter = cms.EDFilter("PythiaFilter",
				MaxEta = cms.untracked.double(3.0),
				MinEta = cms.untracked.double(-3.0),
				MaxPt = cms.untracked.double(9999.0),
				MinPt = cms.untracked.double(3.0),
				ParticleID = cms.untracked.int32(13)
				)


process.cFilter = cms.EDFilter("PythiaFilter",
			       MaxEta = cms.untracked.double(9999.0),
			       MinEta = cms.untracked.double(-9999.0),
			       ParticleID = cms.untracked.int32(4)
			       )


process.generator = cms.EDFilter("Pythia8GeneratorFilter",
				 pythiaPylistVerbosity = cms.untracked.int32(0),
				 filterEfficiency = cms.untracked.double(1.0),
				 pythiaHepMCVerbosity = cms.untracked.bool(False),
				 comEnergy = cms.double(13000.0),
				 crossSection = cms.untracked.double(540000000.0),
				 maxEventsToPrint = cms.untracked.int32(0),
				 PythiaParameters = cms.PSet(
		pythia8CommonSettings = cms.vstring('Tune:preferLHAPDF = 2', 
						    'Main:timesAllowErrors = 10000', 
						    'Check:epTolErr = 0.01', 
						    'Beams:setProductionScalesFromLHEF = off', 
						    'SLHA:keepSM = on', 
						    'SLHA:minMassSM = 1000.', 
						    'ParticleDecays:limitTau0 = on', 
						    'ParticleDecays:tau0Max = 10', 
						    'ParticleDecays:allowPhotonRadiation = on'),
		pythia8CUEP8M1Settings = cms.vstring('Tune:pp 14', 
						     'Tune:ee 7', 
						     'MultipartonInteractions:pT0Ref=2.4024', 
						     'MultipartonInteractions:ecmPow=0.25208', 
						     'MultipartonInteractions:expPow=1.6'),
		processParameters = cms.vstring('SoftQCD:nonDiffractive = on'),
		parameterSets = cms.vstring('pythia8CommonSettings', 
					    'pythia8CUEP8M1Settings', 
					    'processParameters')
		)
				 )


process.ProductionFilterSequence = cms.Sequence(process.generator+process.cFilter+process.muFilter)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.FEVTDEBUGoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
