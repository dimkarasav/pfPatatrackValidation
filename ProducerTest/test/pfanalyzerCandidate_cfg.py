import FWCore.ParameterSet.Config as cms

process = cms.Process('LL')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')



## ----------------- Global Tag -----------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


#--------------------- Report and output ---------------------------   

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService=cms.Service("TFileService",
#                                 fileName=cms.string("full_hist_zmmPU_Dec19th.root")
                                 fileName=cms.string("Patatrack_QCD_PU.root")
#                                 fileName=cms.string("delme.root")
#                                 fileName=cms.string("FullTracking_Candidates_pigun60GeV_noPU.root")
                                 )

process.options = cms.untracked.PSet(
        allowUnscheduled = cms.untracked.bool(True),
        wantSummary = cms.untracked.bool(False),
)



##-------------------- Define the source  ----------------------------

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(

####'file:outputScoutingPF_pixelTracking_gpu_DevMC_ZJetsPUcleanedClassified.root'
#        'file:outputScoutingPF_fullTrackingDevMCFullZmmNoPU.root',
 #######       'file:outputScoutingPF_fullTrackingDevMCFullZmmPU.root',
#        'file:outputScoutingPF_fullTrackingDevMCFullttbarNoPU.root',
#        'file:outputScoutingPF_fullTrackingDevMCFullttbarPU.root',
#        'file:outputScoutingPF_pixelTracking_gpu_DevMC_ZmmNoPU.root',
#        'file:outputScoutingPF_pixelTracking_gpu_DevMC_ZmmPU.root',
#        'file:outputScoutingPF_pixelTracking_gpu_DevMC_ZmmPUZeta01.root',
#        'file:outputScoutingPF_pixelTracking_gpu_DevMC_ttbarNoPU.root',
#        'file:outputScoutingPF_pixelTracking_gpu_DevMC_ttbarPU.root',
#        'file:outputScoutingPF_pixelTracking_gpu_DevMC_ttbarPUZeta01.root'
#        'file:outputScoutingPF_pixelTracking_gpu_DevMC_zmmPUcleanedClassified.root'
#        'file:outputScoutingPF_pixelTracking_gpu_DevMC_ttbarPUcleanedClassified.root'


#        'file:outputScoutingPF_fullTrackingDevMCFullttbarPU.root'
#         'file:outputScoutingPF_pixelTracking_gpu_DevMC_ZJetsPUcleanedClassified.root'

#         'file:outputScoutingPF_pixelTracking_gpu_DevMC_ttbarPUcleanedClassified.root'

       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_0.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_1.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_2.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_3.root',
 #      '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_4.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_5.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_6.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_7.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_8.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_9.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_10.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_11.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_12.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_13.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_14.root',
       '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_15.root'
 #   'file:/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre7/src/Jakobs_producer/ProducerTest/jakobs_root_files/outputScoutingPF_DijetFullTracking.root'
#    '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/outputScoutingPF_PatatrackPixelTracks_ttbar_noPU_LooseCuts_linkFix.root'
 #   '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/QCD_sample/outputScoutingPF_PatatrackPixelTracks_QCD_FlatPU30to80_linkFix_optimized_0.root'
#    '/store/user/dkarasav/ScoutingFiles/edmOutputs/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/outputScoutingPF_PatatrackPixelTracks_ttbar_PU_linkFix_optimized.root'
#    'file:/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_1_0_pre2_Patatrack/src/pfPatatrackValidation/ProducerTest/delme.root'
#    'file:/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre13/outputScoutingPF_delme.root'
 #       'file:outputScoutingPF_pixelTracking_gpu_DevMC_ZJetsPU_CleanedZeta03_3vtx_Dec19.root'
    )
)

##--- l1 stage2 digis ---
process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" )


##-------------------- User analyzer  --------------------------------
#import trigger conf


process.dijetscouting = cms.EDAnalyzer(
    'pfTreeProducer_AddedCandidates',
    ## JETS/MET ########################################
    jetsAK4    = cms.InputTag('hltScoutingPFPacker'),
    #jetsAK4    = cms.InputTag('hltScoutingCaloPacker'),
    muons      = cms.InputTag('hltScoutingMuonPacker'),
    pfcands    = cms.InputTag("hltScoutingPFPacker"),
    metpt      = cms.InputTag('hltScoutingPFPacker:pfMetPt'),
    metphi     = cms.InputTag('hltScoutingPFPacker:pfMetPhi'),
#    metpt      = cms.InputTag('hltScoutingCaloPacker:caloMetPt'),
#    metphi     = cms.InputTag('hltScoutingCaloPacker:caloMetPhi'),
    genJet     = cms.InputTag('ak4GenJets'),
    primVtx     = cms.InputTag('hltScoutingPrimaryVertexPacker:primaryVtx'),
    genpart    = cms.InputTag('genParticles'),
    ptHat            = cms.untracked.InputTag('generator'),
    pu               = cms.untracked.InputTag('addPileupInfo') # 
    
)


# ------------------ path --------------------------

process.p = cms.Path(process.dijetscouting)
