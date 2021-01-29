from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')


config.General.requestName = 'QCD_Validate_11_2_0'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.section_('JobType')

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'pfanalyzerOnCrab_cfg.py'
#config.JobType.psetName = 'scoutingPF_HLTFullTrackingPF_forCrab.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.numCores = 1
config.section_('Data')

#config.Data.inputDataset = '/QCD_Pt-15to3000_TuneCP5_Flat_14TeV_pythia8/Run3Winter20DRMiniAOD-DRFlatPU30to80_110X_mcRun3_2021_realistic_v6-v2/GEN-SIM-RAW'
#config.Data.inputDataset = '/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/Run3Winter20DRPremixMiniAOD-110X_mcRun3_2021_realistic_v6-v1/GEN-SIM-RAW'
#config.Data.userInputFiles = open('/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/list_QCD_Patatrack_TrackAlgoCut5_ZetaCut5.txt').readlines()
config.Data.userInputFiles = open('/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/Producer_11_2_0/CMSSW_11_2_0_Patatrack/src/pfpatatrackvalidation/ProducerTest/test/list_QCD_11_2_0.txt').readlines()
#config.Data.userInputFiles = open('/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/list_QCD_FullTracking.txt').readlines()
config.Data.outputPrimaryDataset = 'QCD_Validate_11_2_0'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
config.Data.outLFNDirBase = '/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/IncludingMet/Validate_11_2_0/'
config.Data.outputDatasetTag = 'test'
config.section_('Site')

config.Site.storageSite = 'T2_CH_CERN'
