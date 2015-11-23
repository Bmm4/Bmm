from CRABClient.UserUtilities import config
config = config()

config.General.requestName     = 'singleMuon-2015Dv4'
config.General.workArea        = 'crab_data'
config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'inclb-prompt-Run2015Dv4-crab.py'
config.JobType.outputFiles = ['inclb-prompt-Run2015Dv4-XXXX.root']

config.Data.inputDataset  = '/SingleMuon/Run2015D-PromptReco-v4/AOD'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'LumiBased'
config.Data.unitsPerJob   = 50
config.Data.outLFNDirBase = '/store/user/ursl/inclb/cmsRun/v01'
config.Data.publication   = False

config.Site.storageSite = 'T3_CH_PSI'
