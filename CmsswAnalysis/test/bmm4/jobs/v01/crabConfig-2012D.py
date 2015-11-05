from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'BMM4_Run2012D'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'bmm-rereco-Run2012-XXXX.py'
config.JobType.outputFiles = ['bmm-rereco-Run2012-XXXX.root']

config.Data.inputDataset = '/MuOnia/Run2012D-22Jan2013-v1/AOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.outLFNDirBase = '/store/user/ursl/bmm4/cmsRun/v01/Run2012D__MuOnia__22Jan2013-v1'
config.Data.publication = False

config.Site.storageSite = 'T3_CH_PSI'
