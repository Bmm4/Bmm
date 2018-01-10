# Auto generated configuration file
# using: 
# Revision: 1.381.2.28 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: step1 --filein file:BdToMuMu_step1-40000.root --fileout file:BdToMuMu_step2-40000.root --pileup_input dbs:/MinBias_TuneZ2_7TeV-pythia6/Summer11Leg-START53_LV4-v1/GEN-SIM --mc --eventcontent RAWSIM --pileup 2011_FinalDist_OOTPU --datatier GEN-RAW --conditions START53_LV6::All --step DIGI,L1,DIGI2RAW,HLT:2011 --python_filename step2.py --no_exec
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_2011_FinalDist_OOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_2011_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
FILE1 = os.environ.get('FILE1')
FILE2 = os.environ.get('FILE2')
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(FILE1)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.28 $'),
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string(FILE2),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-RAW')
    )
)

# Additional output definition

# Other statements
process.mix.input.fileNames = cms.untracked.vstring(['/store/user/ursl/files/mix/Summer11/00064CCC-A218-E311-A2E9-D485646A4E1A.root', '/store/user/ursl/files/mix/Summer11/000C8549-8C18-E311-85EB-80000048FE80.root', '/store/user/ursl/files/mix/Summer11/0015EB9A-AD18-E311-8F14-003048D47662.root', '/store/user/ursl/files/mix/Summer11/003139D6-B218-E311-B020-0022194FFF5A.root', '/store/user/ursl/files/mix/Summer11/0061B6CC-A518-E311-A147-0023AEFDEEA8.root', '/store/user/ursl/files/mix/Summer11/0079A0C1-A418-E311-A93A-D4AE527EE013.root', '/store/user/ursl/files/mix/Summer11/0094049B-B018-E311-BC4A-009C02AAB258.root', '/store/user/ursl/files/mix/Summer11/00979D6B-A818-E311-B37D-00145E5521AD.root', '/store/user/ursl/files/mix/Summer11/00A8B0F3-B018-E311-A735-009C02AAB484.root', '/store/user/ursl/files/mix/Summer11/00D6482E-B018-E311-87E6-00266CF89130.root'])
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_LV6::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWSIMoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
