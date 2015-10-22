from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'b2mu-noevtgen-25ns'
config.General.workArea = 'crab_mc'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'inclb-mc-RunIISpring15DR74-XXXX.py'
config.JobType.outputFiles = ['inclb-mc-RunIISpring15DR74-XXXX.root']

config.Data.inputDataset = '/InclusivebtoMu_bMuonFilter_TuneCUEP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.outLFNDirBase = '/store/user/ursl/inclb/cmsRun/v01/b2mu-noevtgen-25ns'
config.Data.publication = False

config.Site.storageSite = 'T3_CH_PSI'
