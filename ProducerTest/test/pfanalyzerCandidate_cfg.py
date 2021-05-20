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
                                 fileName=cms.string("FullTracking_Zprime_M50_pT300_reclusteredAK8.root")
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



#	'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/201217_125524/0000/outputScoutingPFRun3_v17newPath_Zprime_1.root', #tighht
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/201217_125524/0000/outputScoutingPFRun3_v17newPath_Zprime_3.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/201217_125524/0000/outputScoutingPFRun3_v17newPath_Zprime_4.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/201217_125524/0000/outputScoutingPFRun3_v17newPath_Zprime_5.root'

  #     '/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/201217_130130/0000/outputScoutingPFRun3_v17newPath_Zprime_1.root',# loose
  #      '/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/201217_130130/0000/outputScoutingPFRun3_v17newPath_Zprime_3.root',
# '/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/201217_130130/0000/outputScoutingPFRun3_v17newPath_Zprime_4.root',
# '/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/201217_130130/0000/outputScoutingPFRun3_v17newPath_Zprime_5.root'

#	'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/210325_134915/0000/outputScoutingPFRun3_1.root',
#	'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/210325_134915/0000/outputScoutingPFRun3_3.root',
#	'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test/210325_134915/0000/outputScoutingPFRun3_2.root'

#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M200_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackLoose_ZprimeM200/210422_153648/0000/outputScoutingPFRun3_v17newPath_Zprime_1.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M200_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackLoose_ZprimeM200/210422_153648/0000/outputScoutingPFRun3_v17newPath_Zprime_2.root'

# '/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M200_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackTight_ZprimeM200/210422_150755/0000/outputScoutingPFRun3_v17newPath_Zprime_1.root'

#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M200_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM200/210422_154401/0000/outputScoutingPFRun3_1.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/VectorZPrimeToQQ_M200_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM200/210422_154401/0000/outputScoutingPFRun3_2.root'

#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M100_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM100/210504_112150/0000/outputScoutingPFRun3_1.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M100_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM100/210504_112150/0000/outputScoutingPFRun3_2.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M100_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM100/210504_112150/0000/outputScoutingPFRun3_4.root'


#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M100_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackLoose_ZprimeM100/210504_111814/0000/outputScoutingPFRun3_v17newPath_Zprime_1.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M100_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackLoose_ZprimeM100/210504_111814/0000/outputScoutingPFRun3_v17newPath_Zprime_3.root'

#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M100_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackTight_ZprimeM100/210504_105309/0000/outputScoutingPFRun3_v17newPath_Zprime_1.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M100_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackTight_ZprimeM100/210504_105309/0000/outputScoutingPFRun3_v17newPath_Zprime_3.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M100_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackTight_ZprimeM100/210504_105309/0000/outputScoutingPFRun3_v17newPath_Zprime_4.root'

#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackTight_ZprimeM50/210507_170938/0000/outputScoutingPFRun3_v17newPath_Zprime_3.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackTight_ZprimeM50/210507_170938/0000/outputScoutingPFRun3_v17newPath_Zprime_4.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackTight_ZprimeM50/210507_170938/0000/outputScoutingPFRun3_v17newPath_Zprime_5.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackTight_ZprimeM50/210507_170938/0000/outputScoutingPFRun3_v17newPath_Zprime_6.root'


#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackLoose_ZprimeM50/210507_170822/0000/outputScoutingPFRun3_v17newPath_Zprime_3.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackLoose_ZprimeM50/210507_170822/0000/outputScoutingPFRun3_v17newPath_Zprime_4.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackLoose_ZprimeM50/210507_170822/0000/outputScoutingPFRun3_v17newPath_Zprime_5.root',
#'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_PatatrackLoose_ZprimeM50/210507_170822/0000/outputScoutingPFRun3_v17newPath_Zprime_6.root'

'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM50/210507_170711/0000/outputScoutingPFRun3_1.root',
'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM50/210507_170711/0000/outputScoutingPFRun3_2.root',
'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM50/210507_170711/0000/outputScoutingPFRun3_3.root',
'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM50/210507_170711/0000/outputScoutingPFRun3_5.root',
'/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/QCD_forSubstructure/VectorZPrimeToQQ_M50_pT300_TuneCP5_14TeV_madgraph_pythia8/test_FullTracking_ZprimeM50/210507_170711/0000/outputScoutingPFRun3_6.root'
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
    genJet     = cms.InputTag('ak8GenJets'),
    primVtx     = cms.InputTag('hltScoutingPrimaryVertexPacker:primaryVtx'),
    genpart    = cms.InputTag('genParticles'),
    ptHat            = cms.untracked.InputTag('generator'),
    pu               = cms.untracked.InputTag('addPileupInfo') # 
    
)


# ------------------ path --------------------------

process.p = cms.Path(process.dijetscouting)
