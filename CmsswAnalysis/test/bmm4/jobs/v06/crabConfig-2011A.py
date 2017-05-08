from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'BMM4_Run2011A'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'bmm-rereco-Run2011-XXXX.py'
config.JobType.outputFiles = ['bmm-rereco-Run2011-XXXX.root']

config.Data.inputDataset = '/MuOnia/Run2011A-12Oct2013-v1/AOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.outLFNDirBase = '/store/user/ursl/bmm4/cmsRun/v06/bmmMuOnia2011A'
config.Data.publication = False

config.Site.storageSite = 'T3_CH_PSI'
