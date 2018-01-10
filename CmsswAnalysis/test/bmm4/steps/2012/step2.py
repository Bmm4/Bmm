# Auto generated configuration file
# using:
# Revision: 1.381.2.28
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v
# with command line options: step1 --filein dbs:/LambdaBToJpsiPK_BMuonFilter_MSEL5_TuneZ2star_8TeV-pythia6-evtgen/Summer12-START53_V7C-v1/GEN-SIM --fileout file:BPH-Summer12DR53X-00214_step1.root --pileup_input dbs:/MinBias_TuneZ2star_8TeV-pythia6/Summer12-START50_V13-v3/GEN-SIM --mc --eventcontent RAWSIM --runsAndWeightsForMC [(190482,0.924) , (194270,4.811), (200466,7.21), (207214,7.631)] --pileup fromDB --datatier GEN-SIM-RAW --conditions START53_V19F::All --step DIGI,L1,DIGI2RAW,HLT:cached:fromSource --python_filename BPH-Summer12DR53X-00214_1_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 144
import os
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_fromDB_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
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
    annotation = cms.untracked.string('step1 nevts:144'),
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
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    )
)

# Additional output definition

# Other statements
import SimGeneral.Configuration.ThrowAndSetRandomRun as ThrowAndSetRandomRun
ThrowAndSetRandomRun.throwAndSetRandomRun(process.source,[(190482, 0.92400000000000004), (194270, 4.8109999999999999), (200466, 7.21), (207214, 7.6310000000000002)])
process.mix.input.fileNames = cms.untracked.vstring(['/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/0005E496-3661-E111-B31E-003048F0E426.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/003EEBD4-8061-E111-9A23-003048D437F2.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/005825F1-F260-E111-BD97-003048C692DA.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/0065594C-B35E-E111-8B8C-003048C693EA.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/0091FFD0-6B5E-E111-92FE-003048C693DA.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/00C69AE3-FE60-E111-BC48-0030487D8633.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/02079692-AC61-E111-97BB-0025901D4D54.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/021BD915-2D61-E111-8BAD-002481E76052.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/02386E4D-DC5E-E111-9413-00266CF1074C.root', '/store/user/ursl/files/mix/Summer12/MinBias_TuneZ2star_8TeV-pythia6/02446A08-515E-E111-82C7-00266CF330B8.root'])
import HLTrigger.Configuration.Utilities
import Configuration.HLT.cachedHLT
process.loadCachedHltConfiguration( process.source.setRunNumber.value(), False )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V19F::All', '')

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

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
