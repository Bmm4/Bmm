# Auto generated configuration file
# using:
# Revision: 1.381.2.28
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v
# with command line options: step2 --filein file:BPH-Summer12DR53X-00214_step1.root --fileout file:BPH-Summer12DR53X-00214.root --mc --eventcontent AODSIM,DQM --datatier AODSIM,DQM --conditions START53_V19F::All --step RAW2DIGI,L1Reco,RECO,VALIDATION:validation_prod,DQM:DQMOfflinePOGMC --python_filename BPH-Summer12DR53X-00214_2_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 144
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
FILE2 = os.environ.get('FILE2')
FILE3 = os.environ.get('FILE3')
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring(FILE2)
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.28 $'),
    annotation = cms.untracked.string('step2 nevts:144'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = process.AODSIMEventContent.outputCommands,
    fileName = cms.untracked.string(FILE3),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('AODSIM')
    )
)

#process.DQMoutput = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    outputCommands = process.DQMEventContent.outputCommands,
#    fileName = cms.untracked.string('file:BPH-Summer12DR53X-00214_inDQM.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('DQM')
#    )
#)

# Additional output definition

# Other statements
process.mix.playback = True
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V19F::All', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
#process.dqmoffline_step = cms.Path(process.DQMOfflinePOGMC)
process.validation_step = cms.EndPath(process.validation_prod)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)
#process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.validation_step,process.dqmoffline_step,process.endjob_step,process.AODSIMoutput_step,process.DQMoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.validation_step,process.endjob_step,process.AODSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
#from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
#process = addMonitoring(process)

# End of customisation functions
