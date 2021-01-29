import FWCore.ParameterSet.Config as cms

process = cms.Process('LL')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')



## ----------------- Global Tag -----------------
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


#--------------------- Report and output ---------------------------   

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.TFileService=cms.Service("TFileService",
#                                 fileName=cms.string("full_hist_zmmPU_Dec19th.root")
#                                 fileName=cms.string("Patatrack_QCD_reclusteredAK8_looseCuts.root")
                                 fileName=cms.string("Patatrack_QCD_11_2_0.root")
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
    fileNames = cms.untracked.vstring()
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
#    jetsAK4    = cms.InputTag('hltScoutingPFPackerAK8'),
    #jetsAK4    = cms.InputTag('hltScoutingCaloPacker'),
    muons      = cms.InputTag('hltScoutingMuonPacker'),
    pfcands    = cms.InputTag("hltScoutingPFPacker"),
    metpt      = cms.InputTag('hltScoutingPFPacker:pfMetPt'),
    metphi     = cms.InputTag('hltScoutingPFPacker:pfMetPhi'),
#    metpt      = cms.InputTag('hltScoutingCaloPacker:caloMetPt'),
#    metphi     = cms.InputTag('hltScoutingCaloPacker:caloMetPhi'),
#    genJet     = cms.InputTag('ak8GenJets'),
    genJet     = cms.InputTag('ak4GenJets'),
    primVtx     = cms.InputTag('hltScoutingPrimaryVertexPacker:primaryVtx'),
    genpart    = cms.InputTag('genParticles'),
    ptHat            = cms.untracked.InputTag('generator'),
    pu               = cms.untracked.InputTag('addPileupInfo') # 
    
)


# ------------------ path --------------------------

process.p = cms.Path(process.dijetscouting)
