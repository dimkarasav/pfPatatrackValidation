cmsrel CMSSW_11_2_0_Patatrack
cd CMSSW_11_2_0_Patatrack/src/
git clone https://gitlab.cern.ch/dkarasav/pfpatatrackvalidation.git
cmsenv
scram b

Instructions on how to perform comparisons of jet characteristics between two samples reconstructed with different cmssw configurations. 

(The different configurations usually are reconstructing jets with patatrack pixel tracks or with full tracks.)


1) Producing ntuples

	This is a two-step procedure at the moment. 

	a) Running the cmssw reconstruction configuration

	Follow instructions on scouting twiki to get the latest patatrack release and get the configuration:
	https://twiki.cern.ch/twiki/bin/viewauth/CMS/CMSScouting

	In the module "hltOutputScoutingPF" in the configuration, some needed collections are missing. Make sure the following collections are saved (add by hand the ones missing) :

		  'keep *_hltFEDSelectorL1_*_*',
		  'keep *_hltScoutingEgammaPacker_*_*',
		  'keep *_hltScoutingMuonPackerNoVtx_*_*',
		  'keep *_hltScoutingPFPacker_*_*',
		  'keep *_hltScoutingPrimaryVertexPacker_*_*',
		  'keep *_hltScoutingTrackPacker_*_*',
		  'keep ScoutingPFJets_hltScoutingPFPackerAK8_*_*',
		  'keep edmTriggerResults_*_*_*',
		  'keep recoGenJets_ak4GenJetsNoNu_*_*',
		  'keep recoGenJets_ak4GenJets_*_*',
		  'keep recoGenMETs_genMetTrue_*_*',
		  'keep GenEventInfoProduct_*_*_*',
		  'keep GenRunInfoProduct_*_*_*',
		  'keep *_addPileupInfo__*',
		  'keep recoGenParticles_genParticles_*_*' 


	To run the configuration do either cmsRun or submit the jobs to crab (suggested option for large input samples)

	This gives you one (or many in case of job splitting) output.root file containing all the collections saved above for the processed events.

	b) Produce ntuples with pfTreeProducer


	Producer is located at pfpatatrackvalidation/ProducerTest/plugins/pfTreeProducer.cc

	The configuration for the producer is : 

	pfpatatrackvalidation/ProducerTest/test/pfanalyzerCandidate_cfg.py to run locally via cmsRun
	pfpatatrackvalidation/ProducerTest/test/pfanalyzerOnCrab_cfg.py + crab_test.py   to send the jobs to crab

	The inputs here are the outputs of step (1a) and can be given either directly on "pfanalyzerCandidate_cfg.py"  to run locally or 
	via a .txt on crab_test.py to run on crab.



2) Making comparison plots. 

Instructions on how to run the scripts to produce comparisons between two samples that were recontructed via (1) :

i) Comparisons of jet pT resolution and response
ii) Comparisons of PF jet characteristics
iii) Comparisons of PF candidate characteristics


The inputs for these scripts are the outputs of step (1b). These scripts are used to compare jet characteristics of two (or more) samples (e.g QCD jets reconstructed with pixel tracks vs QCD jets reconstructed with full tracking).

 Each sample should have just one input ntuple file (if crab was used and the splitted jobs resulted in more than one output.root file per sample you should hadd them into a single file for each sample ).



!!!All the inputs and arguments are defined in the file "C_scripts/include/Input_definition.h"

In the directory "C_scripts" there are all the scripts to make and plot histograms for the comparisons.


	a) Edit the "C_scripts/include/Input_definition.h" file.
	
	Make sure that all the inputs are given correctly (path to input ntuple files, output directory, legend_array, pTcuts, etc..)

	b) Make .root files with histograms and create the plots.

	
		i) Jet pT resolution/response comparison plots

		cd C_scripts
		root -l
		.L Make_pTreco_ov_pTgen_histos.C++
		Make_pTreco_ov_pTgen_histos()
		.q
		root -l
		.L Make_pTreco_ov_pTgen_JEScorrected_histos.C++
		Make_pTreco_ov_pTgen_JEScorrected_histos()
		.q
		root -l
		.L Plot_resolution_and_response_histos.C++
		Plot_resolution_and_response_histos()

		Also option to plot the pTreco/pTgen distributions with the Gaussian fits -- helps in debugging
		
		root -l
		.L Plot_pTreco_ov_pTgen.C++
		Plot_pTreco_ov_pTgen(false)   ## "false" for raw pTreco/pTgen and "true" for the response-corrected ones.



		ii) Jet characteristic comparison plots

		(For comparison of matched gen-reco jets)
		root -l
		.L Make_Tracking_Comparisons_matched_histos.C++ 
		Make_Tracking_Comparisons_matched_histos()
		.q
		root -l
		.L Plot_Tracking_Comparisons_matched_histos.C++
		Plot_Tracking_Comparisons_matched_histos()


		(For comparison of back-to-back dijets)  

		root -l 
		.L Make_Tracking_Comparisons_histos_btb.C++
		Make_Tracking_Comparisons_histos_btb()
		.q
		.L Plot_Tracking_Comparisons_matched_histos.C++
		Plot_Tracking_Comparisons_matched_histos()
		.q


		iii) PF candidate characteristic comparisons

		root -l 
		.L Make_PFCandidate_comparison_matched_histos.C++
		Make_PFCandidate_comparison_matched_histos()
		.q

		root -l 
		.L Plot_PFCandidate_comparison_matched_histos.C++
		Plot_PFCandidate_comparison_matched_histos()		
		.q		





