static char input_files[][800] ={

//"/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/IncludingMet/FullTracking/FullTracking_QCD_reclusteredAK8.root",
//"/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/IncludingMet/Patatrack_tightCuts/Patatrack_QCD_reclusteredAK8_tightCuts.root",
//"/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/IncludingMet/Patatrack_looseCuts/Patatrack_QCD_reclusteredAK8_looseCuts.root"

//"/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/IncludingMet/Validate_11_1_6/QCD_Validate_11_1_6_withL1/Patatrack_QCD_11_1_6_withL1.root",
//"/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/IncludingMet/Patatrack_looseCuts/Patatrack_QCD_reclusteredAK8_looseCuts.root"

"/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/IncludingMet/Validate_11_2_0/QCD_Validate_11_2_0/test/210122_134030/0000/Patatrack_QCD_11_2_0.root",
"/eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/dimitris/ScoutingQCDtoBeTransfered/IncludingMet/Patatrack_looseCuts/Patatrack_QCD_reclusteredAK8_looseCuts.root"

//"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8_Patatrack/src/pfpatatrackvalidation/ProducerTest/FullTracking_Zprime_M50_pT300_reclusteredAK8.root",
//"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8_Patatrack/src/pfpatatrackvalidation/ProducerTest/PatatrackTight_Zprime_M50_pT300_reclusteredAK8.root",
//"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8_Patatrack/src/pfpatatrackvalidation/ProducerTest/PatatrackLoose_Zprime_M50_pT300_reclusteredAK8.root"

//"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/Producer_11_2_0/CMSSW_11_2_0_Patatrack/src/pfpatatrackvalidation/ProducerTest/test/Fulltracking_QCD_ttbar_11_2_0.root",
//"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/Producer_11_2_0/CMSSW_11_2_0_Patatrack/src/pfpatatrackvalidation/ProducerTest/test/PatatrackLoose_QCD_ttbar_11_2_0.root"


};

static char analyzer_path[500] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/"}; 
//static char output_directory[200] = {"Zprime_plots/"}; 
static char output_directory[200] = {"Validate_11_2_0/1120_vs_111X/pT_inclusive//"}; 
//static char output_directory[200] = {"QCD_sample_plots/Full_statistics/matced/All_pT_regime/FT_vs_PatLoose_vs_PatTight/"}; 
//static char image_name[200] = {"btb"}; 
static char image_name[200] = {"matched_ak4_GausFits"}; 


//static char legend_array[][500] = { "FullTracking" ,  "PatatrackTight", "PatatrackLoose", "8_25GeVthreshold", "8_45GeVthreshold", "15_45GeVthreshold"  };
//static char legend_array[][500] = { "FullTracking" ,  "PatatrackLoose"  };
//static char legend_array[][500] = { "FullTracking" ,  "PatatrackTight", "PatatrackLoose"  };
static char legend_array[][500] = { "Patatrack_11_2_0", "Patatrack_11_1_0"  };

static double yBnd[]={0.0, 1.3, 2.4, 2.7, 3.0};  
//static double yBnd[]={0.0, 1.3, 2.4};  
//static double yBnd[]={0.0, 0.3, 0.6, 0.9, 1.1,  1.3};  

//this is only used for the jet pT resolution calculations 
//static double ptBnd[] = {30.,   70., 110., 150. , 190. , 230., 270., 310., 350.}; //low pT regime //this is only used for the jet pT resolution calculations 
//static double ptBnd[] = {350.,   410., 470., 530. , 590. , 650., 710., 770., 830.}; //middle pT regime  
//static double ptBnd[] = {830.,   900., 1000., 1150. , 1300. , 1550., 1850., 2100., 3000.}; //high pT regime 


//static double ptBnd[] = {30.,  110.,  190. , 230., 270.,  350., 500., 1000.}; //ttbar relval  //this is only used for the jet pT resolution calculations 


static double ptBnd[] = {30.,   70., 110., 150. , 190. , 230., 270., 310., 350.,  410., 470., 530. , 590. , 650., 710., 770., 830.,   900., 1000., 1150. , 1300. , 1550., 1850., 2100., 3000.}; //All pT regime //this is only used for the jet pT resolution calculations 
//static double ptBnd[] = {30., 110.,  190. , 250.,  350., 450., 600.}; //Zprime pT range 

static double pTlowCut  = 30; //cuts on jet pT's to be considered
static double pThighCut = 3000;

static bool useWeights  = true;
static bool useGausFits = true; //only used for jet pT resolution calculations
static bool useRecoBins = true; //calc resolution in reco or gen pT bins


static int EventsToProcess = 6000000; //negative value will process all events in each tree

static double DR_threshold = 0.2; //used for gen-to-reco matching


// ===================================================== Plotting options ==============================================================
static bool plot_matched = true; // choose to plot btb or matched 
static bool scale_histos = true ; //scale all histos to the number of entries of the histo of the 1st input
static bool Save_Plots = true; 


static int Colors[] = { 1, 4, 2 , 6, 3, 7 , 28, 46} ; // black, blue, red , magenta, green, light blue ,  brown, redish. // define colors for each input
static int MarkerStyle[] = { 8, 2, 5 , 4, 22, 21, 27, 28 } ; // 


static int PadColumnsEtaBins = 3;  // define number of pad columns for the canvas
static int PadRowsEtaBins = 2;     // define number of pad rows for the canvas

//for nicely looking plots it should be : eta_bins <= PadColumnsEtaBins * PadRowsEtaBins  

static int Canvas_XpixelsEtaBins = PadColumnsEtaBins*333; //define canvas X pixels
static int Canvas_YpixelsEtaBins = PadRowsEtaBins*500;    //define canvas Y pixels

//static int Canvas_XpixelsEtaBins = PadColumnsEtaBins*450; //define canvas X pixels
//static int Canvas_XpixelsEtaBins = PadColumnsEtaBins*450; //define canvas X pixels
//static int Canvas_YpixelsEtaBins = PadRowsEtaBins*450;    //define canvas Y pixels

static double YaxisLowEndMultiplier = 0.00005; // define lower bound of Y axis for plot: Lower Bound = YaxisLowEndMultiplier* histo_maximum
static double YaxisHighEndMultiplier = 10.;    // define upper bound of Y axis for plot: Upper Bound = YaxisHighEndMultiplier* histo_maximum



//==============================additional options only used for the jet pT resolution & response plots ===============================

static int PadColumnsPtBins = 5; 
static int PadRowsPtBins = 5;

//static int Canvas_XpixelsPtBins = PadColumnsPtBins*333;
//static int Canvas_YpixelsPtBins = PadRowsPtBins*500;

static int Canvas_XpixelsPtBins = PadColumnsPtBins*450;
static int Canvas_YpixelsPtBins = PadRowsPtBins*450;













