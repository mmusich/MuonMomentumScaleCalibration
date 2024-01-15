# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step1 --mc --runUnscheduled --conditions auto:run2_design --step RECO --nThreads 32 --geometry DB:Extended --era Run2_2016 --dasquery file dataset=/MuonGunDesignFwd/bendavid-MuonGunDesignFwd-449be899dc2f7c04b70495711fe9dfd6/USER instance=prod/phys03 --processName RECO2 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('RECO2',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.XMLFromDBSource.label = cms.string("Extended")

process.load("TrackPropagation.Geant4e.geantRefit_cff")
    
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/mc/RunIISummer20UL18RECO/JPsiToMuMu_Pt0to8-pythia8/ALCARECO/TkAlJpsiMuMu-106X_upgrade2018_realistic_v11_L1v1-v2/2820000/0126B990-069F-0541-B6F5-203814BCB153.root',
        '/store/mc/RunIISummer20UL18RECO/JPsiToMuMu_Pt8toInf-pythia8/ALCARECO/TkAlJpsiMuMu-106X_upgrade2018_realistic_v11_L1v1-v2/2830000/0069C60B-EAD7-4243-9951-82D1C097E922.root',
    ),
    secondaryFileNames = cms.untracked.vstring(),
)

process.options = cms.untracked.PSet(

)

process.SiteLocalConfigService = cms.Service("SiteLocalConfigService",
  overrideSourceCacheHintDir = cms.untracked.string("lazy-download"),
)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('step1_RECO.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')

process.offlineBeamSpot = cms.EDProducer("BeamSpotProducer")

process.globalCor = cms.EDProducer('ResidualGlobalCorrectionMakerTwoTrackG4e',
                                   src = cms.InputTag("ALCARECOTkAlJpsiMuMu"),
                                   fitFromGenParms = cms.bool(False),
                                   fitFromSimParms = cms.bool(False),
                                   fillTrackTree = cms.bool(True),
                                   
                                   fillGrads = cms.bool(True),
                                   fillJac = cms.bool(True),
                                   
                                   fillRunTree = cms.bool(True),
                                   
                                   doGen = cms.bool(True),
                                   genParticles = cms.InputTag("genParticles"),
                                   pileupInfo = cms.InputTag("addPileupInfo"),
                                   doSim = cms.bool(False),
                                   requireGen = cms.bool(False),
                                   doMuons = cms.bool(False),
                                   doMuonAssoc = cms.bool(False),
                                   doTrigger = cms.bool(True),
                                   doRes = cms.bool(False),
                                   useIdealGeometry = cms.bool(True),
                                   # useIdealGeometry = cms.bool(False),
                                   bsConstraint = cms.bool(False),
                                   applyHitQuality = cms.bool(True),

                                   doVtxConstraint = cms.bool(True),

                                   doMassConstraint = cms.bool(True),
                                   #pdg values
                                   #massConstraint = cms.double(3.092666073886687),
                                   #massConstraintWidth = cms.double(0.0195),
                                   massConstraint = cms.double(3.092403924504445),
                                   massConstraintWidth = cms.double(0.0203),

                                   corFiles = cms.vstring(),
                                   
                                   #triggers = cms.vstring(["HLT_Dimuon0_Jpsi_Muon",
                                   #                          "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing",
                                   #                          "HLT_Dimuon0er16_Jpsi_NoVertexing",
                                   #                          "HLT_Dimuon10_Jpsi_Barrel",
                                   #                          "HLT_Dimuon13_PsiPrime",
                                   #                          "HLT_Dimuon16_Jpsi",
                                   #                          "HLT_Dimuon20_Jpsi",
                                   #                          "HLT_Dimuon8_PsiPrime_Barrel",
                                   #                          "HLT_DoubleMu4_3_Bs",
                                   #                          "HLT_DoubleMu4_3_Jpsi_Displaced",
                                   #                          "HLT_DoubleMu4_JpsiTrk_Displaced",
                                   #                          "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
                                   #                          "HLT_Mu7p5_Track2_Jpsi",
                                   #                          "HLT_Mu7p5_Track3p5_Jpsi",
                                   #                          "HLT_Mu7p5_Track7_Jpsi"]),

                                   triggers = cms.vstring(["HLT_Dimuon0_Jpsi",
                                                             "HLT_Dimuon0_Jpsi3p5_Muon2",
                                                             "HLT_Dimuon0_Jpsi_L1_4R_0er1p5R",
                                                             "HLT_Dimuon0_Jpsi_L1_NoOS",
                                                             "HLT_Dimuon0_Jpsi_NoVertexing",
                                                             "HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R",
                                                             "HLT_Dimuon0_Jpsi_NoVertexing_NoOS",
                                                             "HLT_Dimuon10_PsiPrime_Barrel_Seagulls",
                                                             "HLT_Dimuon18_PsiPrime",
                                                             "HLT_Dimuon18_PsiPrime_noCorrL1",
                                                             "HLT_Dimuon20_Jpsi_Barrel_Seagulls",
                                                             "HLT_Dimuon25_Jpsi",
                                                             "HLT_Dimuon25_Jpsi_noCorrL1",
                                                             "HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi",
                                                             "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",
                                                             "HLT_DoubleMu4_3_Bs",
                                                             "HLT_DoubleMu4_3_Jpsi",
                                                             "HLT_DoubleMu4_JpsiTrkTrk_Displaced",
                                                             "HLT_DoubleMu4_JpsiTrk_Displaced",
                                                             "HLT_DoubleMu4_Jpsi_Displaced",
                                                             "HLT_DoubleMu4_Jpsi_NoVertexing",
                                                             "HLT_DoubleMu4_PsiPrimeTrk_Displaced",
                                                             "HLT_Mu30_TkMu0_Psi",
                                                             "HLT_Mu7p5_Track2_Jpsi",
                                                             "HLT_Mu7p5_Track3p5_Jpsi",
                                                             "HLT_Mu7p5_Track7_Jpsi",
                                                             "HLT_Mu7p5_L2Mu2_Jpsi",
                                                             "HLT_Dimuon0_LowMass_L1_0er1p5R",
                                                             "HLT_Dimuon0_LowMass_L1_4R",
                                                             "HLT_Dimuon0_LowMass_L1_0er1p5",
                                                             "HLT_Dimuon0_LowMass",
                                                             "HLT_Dimuon0_LowMass_L1_4"]),

                                   ## upsilon triggers (MuOnia PD)
                                   # triggers = cms.vstring(['HLT_Dimuon0_Phi_Barrel',
                                   #     'HLT_Dimuon0_Upsilon_Muon',
                                   #     'HLT_Dimuon13_Upsilon',
                                   #     'HLT_Dimuon8_Upsilon_Barrel',
                                   #     'HLT_Mu16_TkMu0_dEta18_Onia',
                                   #     'HLT_Mu16_TkMu0_dEta18_Phi',
                                   #     'HLT_Mu25_TkMu0_dEta18_Onia',
                                   #     'HLT_Mu7p5_L2Mu2_Upsilon',
                                   #     'HLT_Mu7p5_Track2_Upsilon',
                                   #     'HLT_Mu7p5_Track3p5_Upsilon',
                                   #     'HLT_Mu7p5_Track7_Upsilon',
                                   #     'HLT_QuadMuon0_Dimuon0_Upsilon']),

                                   MagneticFieldLabel = cms.string(""),

)


# Path and EndPath definitions
#process.reconstruction_step = cms.Path(process.reconstruction_fromRECO)
#process.reconstruction_step = cms.Path(process.globalCor*process.globalCorOut)
#process.reconstruction_step = cms.Path(process.trackExtraRekeyer*process.globalCor)
#process.reconstruction_step = cms.Path(process.globalCor)
#process.reconstruction_step = cms.Path(process.geopro*process.globalCor)
process.reconstruction_step = cms.Path(process.geopro*process.offlineBeamSpot*process.globalCor)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

