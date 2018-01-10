# Auto generated configuration file
# using:
# Revision: 1.19
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v
# with command line options: step1 --filein file:bla.root --fileout file:BPH-RunIISummer16DR80Premix-00009_step1.root --pileup_input file:bla.root --mc --eventcontent PREMIXRAW --datatier GEN-SIM-RAW --conditions 80X_mcRun2_asymptotic_2016_TrancheIV_v6 --step DIGIPREMIX_S2,DATAMIX,L1,DIGI2RAW,HLT:@frozen2016 --nThreads 1 --datamix PreMix --era Run2_2016 --no_exec
import os
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HLT',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.DigiDMPreMix_cff')
process.load('SimGeneral.MixingModule.digi_MixPreMix_cfi')
process.load('Configuration.StandardSequences.DataMixerPreMix_cff')
process.load('Configuration.StandardSequences.SimL1EmulatorDM_cff')
process.load('Configuration.StandardSequences.DigiToRawDM_cff')
process.load('HLTrigger.Configuration.HLT_25ns15e33_v4_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
FILE1 = os.environ.get('FILE1')
FILE2 = os.environ.get('FILE2')
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(FILE1),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.PREMIXRAWoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RAW'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string(FILE2),
    outputCommands = process.PREMIXRAWEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.mix.digitizers = cms.PSet(process.theDigitizersMixPreMix)
process.mixData.input.fileNames = cms.untracked.vstring([
    '/store/user/ursl/files/mix/2016quater/B227E846-2981-E611-8308-549F3525B220.root',
    '/store/user/ursl/files/mix/2016quater/5CB96541-2A81-E611-B55E-D4AE526DF45D.root',
    '/store/user/ursl/files/mix/2016quater/A6214B5E-2A81-E611-ABC6-782BCB5364D6.root',
    '/store/user/ursl/files/mix/2016quater/3617F0C4-2A81-E611-925F-B083FED04D68.root',
    '/store/user/ursl/files/mix/2016quater/EA8633EA-2A81-E611-956F-90B11C0BCE26.root',
    '/store/user/ursl/files/mix/2016quater/2C594B10-2981-E611-B2E8-549F3525A64C.root',
    '/store/user/ursl/files/mix/2016quater/C6EF5827-2981-E611-BE11-141877411367.root',
    '/store/user/ursl/files/mix/2016quater/48C74D24-2981-E611-92B8-141877410522.root',
    '/store/user/ursl/files/mix/2016quater/78498A3A-2981-E611-BD9E-1418774121A1.root',
    '/store/user/ursl/files/mix/2016quater/7CB78739-2981-E611-8D4B-549F3525DD6C.root',
    '/store/user/ursl/files/mix/2016quater/78B2B1C3-2881-E611-913B-0025905A60B2.root',
    '/store/user/ursl/files/mix/2016quater/960DEA65-2A81-E611-896D-0025905B85B2.root',
    '/store/user/ursl/files/mix/2016quater/A2263AF5-2A81-E611-931E-0CC47A4D762A.root',
    '/store/user/ursl/files/mix/2016quater/32C9391D-2B81-E611-886C-0CC47A78A468.root',
    '/store/user/ursl/files/mix/2016quater/8287BD2C-2A81-E611-8067-0025905A606A.root',
    '/store/user/ursl/files/mix/2016quater/C036F36A-2D81-E611-B1A3-0CC47A1DF7E4.root',
    '/store/user/ursl/files/mix/2016quater/0E7D9F22-3081-E611-B630-0025905A612E.root',
    '/store/user/ursl/files/mix/2016quater/3603F4A6-3081-E611-AE83-0025905A60E0.root',
    '/store/user/ursl/files/mix/2016quater/20963A66-3081-E611-90EE-0025905B85DA.root',
    '/store/user/ursl/files/mix/2016quater/E41C8765-3081-E611-A23D-0CC47A4C8EC6.root',
    '/store/user/ursl/files/mix/2016quater/B014FCA4-3081-E611-9F02-0025905A60B4.root',
    '/store/user/ursl/files/mix/2016quater/FCD95A23-2E81-E611-A7FD-0CC47A4DEE70.root',
    '/store/user/ursl/files/mix/2016quater/7C2C3FB7-2E81-E611-9C99-0CC47A4DEE70.root',
    '/store/user/ursl/files/mix/2016quater/5AFA1EC2-2E81-E611-97CF-0090FAA57FA4.root',
    '/store/user/ursl/files/mix/2016quater/16022BC0-2E81-E611-9DB1-0090FAA57E24.root',
    '/store/user/ursl/files/mix/2016quater/10A78FCA-2F81-E611-913E-0025904CDDF8.root',
    '/store/user/ursl/files/mix/2016quater/B62029D0-2F81-E611-998D-0025905C42FE.root',
    '/store/user/ursl/files/mix/2016quater/4C8A01D0-2F81-E611-A948-0025905C9724.root',
    '/store/user/ursl/files/mix/2016quater/34529AE6-2F81-E611-8D92-0025905C5474.root',
    '/store/user/ursl/files/mix/2016quater/6EB28006-3081-E611-9AD7-0025905C54F6.root',
    '/store/user/ursl/files/mix/2016quater/2A3CB2EF-2981-E611-A073-002590A37116.root',
    '/store/user/ursl/files/mix/2016quater/46E3D874-2981-E611-8B81-001E67792558.root',
    '/store/user/ursl/files/mix/2016quater/D4BDD5A4-2F81-E611-AD66-001E673973EB.root',
    '/store/user/ursl/files/mix/2016quater/A0F620AF-2F81-E611-ADE1-001E67E6F5EE.root',
    '/store/user/ursl/files/mix/2016quater/A6248481-1B81-E611-AACC-A4BADB1E65B1.root',
    '/store/user/ursl/files/mix/2016quater/64DD48AC-1E81-E611-B120-0025905A60DE.root',
    '/store/user/ursl/files/mix/2016quater/B6E57CF1-1D81-E611-B27C-90B11C0BCE80.root',
    '/store/user/ursl/files/mix/2016quater/8EDE02F5-1D81-E611-AC76-001C23C0B763.root',
    '/store/user/ursl/files/mix/2016quater/22C18622-1E81-E611-B9FD-90B11C0BCE80.root',
    '/store/user/ursl/files/mix/2016quater/5039F00C-1F81-E611-AFAC-A4BADB1CF89C.root',
    '/store/user/ursl/files/mix/2016quater/AA2271A7-1F81-E611-8C1E-B083FED045ED.root',
    '/store/user/ursl/files/mix/2016quater/94AD1AB4-1C81-E611-8042-001E677924BA.root',
    '/store/user/ursl/files/mix/2016quater/D6B647D2-1C81-E611-9A97-001E677927A2.root',
    '/store/user/ursl/files/mix/2016quater/EEF6F911-1E81-E611-AAD7-001E67E713AE.root',
    '/store/user/ursl/files/mix/2016quater/18FDF69A-2181-E611-AE10-001E67E6F5AD.root',
    '/store/user/ursl/files/mix/2016quater/64B74EAE-2381-E611-911F-002590FC5ACC.root',
    '/store/user/ursl/files/mix/2016quater/BA36014C-2B81-E611-B32B-0CC47A4D7638.root',
    '/store/user/ursl/files/mix/2016quater/3A2B3652-2B81-E611-AA86-0025905A60E4.root',
    '/store/user/ursl/files/mix/2016quater/A0B3BD5B-2E81-E611-B095-0CC47A4C8E22.root',
    '/store/user/ursl/files/mix/2016quater/9246D75D-2E81-E611-A989-0CC47A78A3D8.root',
    '/store/user/ursl/files/mix/2016quater/620AEC78-2E81-E611-914B-0CC47A4D7674.root',
    '/store/user/ursl/files/mix/2016quater/7AA30FD5-2581-E611-8102-0025905B8580.root',
    '/store/user/ursl/files/mix/2016quater/FACBF30A-2C81-E611-9197-0CC47A4C8E56.root',
    '/store/user/ursl/files/mix/2016quater/1224C00F-2C81-E611-86B4-0CC47A7C345C.root',
    '/store/user/ursl/files/mix/2016quater/E675E50F-2C81-E611-BE92-0CC47A4D7664.root',
    '/store/user/ursl/files/mix/2016quater/12B86E0E-2C81-E611-897D-0CC47A4D7632.root',
    '/store/user/ursl/files/mix/2016quater/FE91B6DC-2C81-E611-93B1-0025905A611C.root',
    '/store/user/ursl/files/mix/2016quater/0283F1CF-2C81-E611-9B9D-0CC47A4D7654.root',
    '/store/user/ursl/files/mix/2016quater/A0619ADC-2C81-E611-A26E-0CC47A4D7644.root',
    '/store/user/ursl/files/mix/2016quater/BAEC41D1-2C81-E611-AD7D-0CC47A4D7606.root',
    '/store/user/ursl/files/mix/2016quater/7660C2DF-2C81-E611-8610-0025905B8576.root',
    '/store/user/ursl/files/mix/2016quater/98105F8C-2F81-E611-AD8B-14187741121F.root',
    '/store/user/ursl/files/mix/2016quater/344FFA74-3081-E611-877A-001C23C0F203.root',
    '/store/user/ursl/files/mix/2016quater/B2970274-3081-E611-AC07-141877411970.root',
    '/store/user/ursl/files/mix/2016quater/74E06989-3081-E611-AF4F-782BCB20E8A5.root',
    '/store/user/ursl/files/mix/2016quater/BEAF7997-3081-E611-BB35-842B2B172901.root',
    '/store/user/ursl/files/mix/2016quater/E8CDD275-2981-E611-9F24-0CC47A78A2F6.root',
    '/store/user/ursl/files/mix/2016quater/9AA05B4F-2E81-E611-9057-0CC47A4D76D2.root',
    '/store/user/ursl/files/mix/2016quater/3E5B0D6C-2E81-E611-AF64-0025905B85C6.root',
    '/store/user/ursl/files/mix/2016quater/BC8BA016-2E81-E611-BB61-0025905A60D6.root',
    '/store/user/ursl/files/mix/2016quater/3A3B8414-2E81-E611-BDCF-0CC47A4C8E14.root',
    '/store/user/ursl/files/mix/2016quater/B49A50BE-2E81-E611-9DB4-0CC47A4D76D6.root',
    '/store/user/ursl/files/mix/2016quater/BAD8CAC3-2E81-E611-9C2F-0025905B8564.root',
    '/store/user/ursl/files/mix/2016quater/C26E1CC0-2E81-E611-B501-0025905B857A.root',
    '/store/user/ursl/files/mix/2016quater/FCBE0999-2E81-E611-B2AB-0025905B85FE.root',
    '/store/user/ursl/files/mix/2016quater/182B13C1-2E81-E611-BCF9-0CC47A7C3424.root',
    '/store/user/ursl/files/mix/2016quater/E2E52386-3381-E611-9937-0090FAA57770.root',
    '/store/user/ursl/files/mix/2016quater/BAD9C24C-3481-E611-BD55-0090FAA58864.root',
    '/store/user/ursl/files/mix/2016quater/3E255F86-3481-E611-A980-0025907B4F48.root',
    '/store/user/ursl/files/mix/2016quater/F6B31A4C-3581-E611-973A-0090FAA57560.root',
    '/store/user/ursl/files/mix/2016quater/12FE9747-3381-E611-AD02-0025905C54BA.root',
    '/store/user/ursl/files/mix/2016quater/12AB1611-3481-E611-B896-0025905D1CAE.root',
    '/store/user/ursl/files/mix/2016quater/1ED7C826-3481-E611-A298-0025904C6224.root',
    '/store/user/ursl/files/mix/2016quater/8A3F1A52-3381-E611-92D2-0025905C54FC.root',
    '/store/user/ursl/files/mix/2016quater/A49B2046-1381-E611-A5FF-549F3525A184.root',
    '/store/user/ursl/files/mix/2016quater/E68A9134-1481-E611-8F36-141877411FCD.root',
    '/store/user/ursl/files/mix/2016quater/CE4C6E3F-1481-E611-BCE7-782BCB78621D.root',
    '/store/user/ursl/files/mix/2016quater/562CA136-1481-E611-9CA4-A4BADB1E6031.root',
    '/store/user/ursl/files/mix/2016quater/0E37B2B7-2881-E611-8246-001E67E71CA4.root',
    '/store/user/ursl/files/mix/2016quater/4C52B757-2A81-E611-8429-001E677923F0.root',
    '/store/user/ursl/files/mix/2016quater/4C63A2FA-2B81-E611-85C2-002590A3C978.root',
    '/store/user/ursl/files/mix/2016quater/AE3E07A9-2B81-E611-B192-001E67E71DEE.root',
    '/store/user/ursl/files/mix/2016quater/2AB134B6-2A81-E611-8CBB-001E67E6F67F.root',
    '/store/user/ursl/files/mix/2016quater/287947AC-3881-E611-B7A2-0CC47A4D76C6.root',
    '/store/user/ursl/files/mix/2016quater/642F5CF7-3681-E611-8D99-0025905A60DE.root',
    '/store/user/ursl/files/mix/2016quater/529892F9-3881-E611-BA4D-0025905B85A0.root',
    '/store/user/ursl/files/mix/2016quater/222A735C-3981-E611-A74C-0025905A60BC.root',
    '/store/user/ursl/files/mix/2016quater/8E9110A0-3981-E611-BA9A-0025905A60DE.root',
    '/store/user/ursl/files/mix/2016quater/1A1328F5-1481-E611-9C31-782BCB2058AD.root'
])

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_TrancheIV_v6', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.datamixing_step = cms.Path(process.pdatamix)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.PREMIXRAWoutput_step = cms.EndPath(process.PREMIXRAWoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.datamixing_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.PREMIXRAWoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforFullSim

#call to customisation function customizeHLTforFullSim imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforFullSim(process)

# End of customisation functions
