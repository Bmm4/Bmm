from CRABClient.UserUtilities import config
config = config()

config.General.requestName     = 'doubleMuon-2015B'
config.General.workArea        = 'crab_data'
config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = 'inclb-prompt-Run2015B-crab.py'
config.JobType.outputFiles = ['inclb-prompt-Run2015B-XXXX.root']

config.Data.inputDataset  = '/DoubleMuon/Run2015B-PromptReco-v1/AOD'
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'LumiBased'
config.Data.unitsPerJob   = 200
#config.Data.runRange      = '251168-251493,251496-251548,251559-251562,251638-252126'
config.Data.outLFNDirBase = '/store/user/ursl/inclb/cmsRun/v01'
config.Data.publication   = False

config.Site.storageSite = 'T3_CH_PSI'
