from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'DSTAR_Run2012C'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'dstar-rereco-Run2012-XXXX.py'
config.JobType.outputFiles = ['dstar-rereco-Run2012-XXXX.root']

config.Data.inputDataset = '/MuOnia/Run2012C-22Jan2013-v1/AOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.outLFNDirBase = '/store/user/ursl/bmm4/cmsRun/v01/dstar-2012'
config.Data.publication = False

config.Site.storageSite = 'T3_CH_PSI'
