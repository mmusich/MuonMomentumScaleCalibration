from CRABClient.UserUtilities import config

config = config()

label = "Charmonium2017D_v722_layerbylayer"

config.General.requestName = label
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'pset.py'
config.JobType.psetName = 'RunGlobalCorRecJpsi_2017Data.py'

ncores = 1
outfiles = []
for i in range(ncores):
    outfiles.append("globalcor_data_%i.root" % i)
    #outfiles.append("output_%i.txt" % i)
print(outfiles)

#config.JobType.inputFiles = ["FrameworkJobReport.xml", "RunGlobalCorRecJpsi_2017Data.py", "command.sh"]
#config.JobType.inputFiles = ["RunGlobalCorRecJpsi_2017Data.py", "command.sh"]
config.JobType.outputFiles = outfiles
config.JobType.numCores = ncores
config.JobType.allowUndistributedCMSSW = True
#config.JobType.scriptExe = 'myscript.sh'

config.Data.inputDataset = "/Charmonium/Run2017D-TkAlJpsiMuMu-09Aug2019_UL2017-v1/ALCARECO"
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 10*1000
config.Data.publication = False
config.Data.outputDatasetTag = label


config.Site.storageSite = "T2_IT_Pisa"
