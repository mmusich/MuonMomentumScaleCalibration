from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'Wplus_calibration'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.numCores = 1
config.JobType.maxMemoryMB = 2000
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../configs/config.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['output_0.root']

config.Data.inputDataset = '/WplusJetsToMuNu_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.outLFNDirBase = '/store/group/cmst3/group/wmass/w-mass-13TeV/NanoAOD' 
#config.Data.publication = True
config.Data.outputDatasetTag = 'NanoV9MCPostVFP_TrackFitV722_NanoProdv3'
config.Data.inputDBS = 'global'
config.Data.useParent = False

config.Site.storageSite = 'T2_IT_Pisa'
