#include "include_functions.h"




void PFCandidate_comparison_btb()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
//	gStyle->SetOptStat(0);

	


	char analyzer_path[200] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre7/src/Jakobs_producer/ProducerTest/"}; 
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/noPU/scanning_ptError_PFmodule/"};         
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/PUvsNoPU//"};      
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/noPU/scanning_DptOverPtCut_10_3/"};      
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/mkFit/FullVsmkFit/"};        
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/testing_PFClusterProducer_iteration1/"};
   	char output_directory[200] = {"Realistic_RunIII_14TeV/QCD_sample/Weighted_btb/"};                      
//	char output_directory[200] = {"Realistic_RunIII_14TeV/pigun_60GeV/Relaxing_ZetaCuts/"};             
	char image_name[150] = {"Patatrack"}; 

//	double yBnd[]={0.0, 1.3, 2.4, 2.7, 3.0, 5.0}; 
	double yBnd[]={0.0, 1.3, 2.4, 2.7, 3.0}; 
//	double yBnd[]={0.0, 1.3, 2.4}; 




//	int Colors[NoFiles] = { 1, 4, 2 , 6, 3 , 68} ; // black, blue, red , magenta
	int Colors[] = { 1, 4, 2 , 6, 3, 7 , 28, 46} ; // black, blue, red , magenta, green, light blue ,  brown, redish
	bool scale_on = true ;
	bool Save_Plots = true; 
	bool useWeights = true;

	double pTlowCut_leading = 60;
	double pTlowCut_subleading = 30;
	double pThighCut = 5000;
	double pThistoMax = 450;
	double deltaPhi = 2.7;


	char input_files[][500] ={

//mkFit comparisons
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_LooseTrackParams_noPU_onlyZeta.root" ,
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/mkFit_Candidates_ttbar_noPU.root"
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/PatatrackPixels_Candidates_ttbar_14TeV_noPU.root"

//optimized vs full tracking noPU
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut2_Zeta0p5.root"

//optimized vs full tracking
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/QCD_samples/FullTracking_QCD_PU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/QCD_samples/Patatrack_QCD_PU.root"

//TrackAlgoCut scan
/*
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_ttbar_noPU_TightCuts_linkFix.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut1.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut2.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut3.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut4.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1.root"
*/
//ZetaCut scan
/*
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_ttbar_noPU_TightCuts_linkFix.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut2_Zeta0p2.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut2_Zeta0p5.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut2_Zeta1p0.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut2.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_ttbar_noPU_LooseCuts_linkFix.root"
*/
//optimized vs full PU
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_PU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_PU_linkFix_optimized.root"
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_PU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_ttbar_iteration1.root"
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_PU_linkFix_optimized.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_ttbar_iteration1.root"

//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/PatatrackPixels_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_LooseTrackParams_noPU_onlyZeta.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_ttbar_noPU_TightCuts_linkFix.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_ttbar_noPU_LooseCuts_linkFix.root"
//PU vs no PU
/*
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_PU_nominal.root", 
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/PatatrackPixels_Candidates_ttbar_14TeV_noPU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_PU_LooseCuts.root", 
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_LooseTrackParams_noPU_onlyZeta.root"
*/

//pigun samples
/*
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/pigun_60GeV/FullTracking_Candidates_pigun60GeV_noPU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/pigun_60GeV/Patatrack_Candidates_pigun60GeV_noPU_TightCuts.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/pigun_60GeV/Patatrack_Candidates_pigun60GeV_noPU_LooseCuts.root"
*/
};


	int PadColumnsEtaBins = 3; 
	int PadRowsEtaBins = 2;
	int Canvas_XpixelsEtaBins = PadColumnsEtaBins*333;
	int Canvas_YpixelsEtaBins = PadRowsEtaBins*500;


	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);
	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;

//	char legend_array[NoFiles][500] = { "FullTracking" , "TightCuts", "ZetaCut0p2", "ZetaCut0p5", "ZetaCut1p0", "ZetaCut2p0", "LooseCuts" };
	char legend_array[NoFiles][500] = { "FullTracking" ,  "PatatrackPixelTracks"  };
//	char legend_array[NoFiles][500] = { "MaxIteration_50" ,  "MaxIteration_1"  };
//	char legend_array[NoFiles][500] = { "FullTracking" ,  "mkFitTracking"  };
//	char legend_array[NoFiles][500] = { "FullTracking" , "LooseCuts", "TrackAlgoCut2", "TrackAlgoCut1" ,"TrackCut1-ptError10", "TightCuts"  };
//	char legend_array[NoFiles][500] = { "FullTracking" , "TrackAlgoCut5", "TrackAlgoCut2", "TrackAlgoCut1"  };

//	TFile *f1 = TFile::Open("PatatrackPixels_Candidates_ttbar_14TeV_noPU.root","READ"); //full tracking file

	TFile *f[NoFiles];
	TTree *tree[NoFiles];

//	TTree *tree1 = (TTree*)f1->Get("dijetscouting/tree");


	double pTmax = 300;
	double DR_threshold = 0.2;

	char eta_bins_legend[eta_bins][25];
	char name[256]; 

	int treeEntries[NoFiles];
	int maxEntries;
	bool flag_for_max_entries = true;

	for (int iy=0; iy< eta_bins; iy++)
	{
		if (iy==0) sprintf( eta_bins_legend[iy], "|#eta|<%3.1f" , yBnd[iy+1] );  
		else
		{
			sprintf( eta_bins_legend[iy], "%3.1f<|#eta|<%3.1f" , yBnd[iy] ,yBnd[iy+1] );
		}
	}

	TH1D *h_ChargedHadron_pT[NoFiles][eta_bins], *h_NeutralHadron_pT[NoFiles][eta_bins], *h_ChargedHadron_SumpT[NoFiles][eta_bins], *h_NeutralHadron_SumpT[NoFiles][eta_bins], *h_NumberOfNeutralHadrons[NoFiles][eta_bins], *h_NumberOfChargedHadrons[NoFiles][eta_bins], *h_ChargedHadron_DR[NoFiles][eta_bins], *h_NeutralHadron_DR[NoFiles][eta_bins], *h_girth[NoFiles][eta_bins], *h_quadratic_moment[NoFiles][eta_bins];

	TH2D *h_ChargedHadron_eta_phi[NoFiles][eta_bins], *h_NeutralHadron_eta_phi[NoFiles][eta_bins];


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int iy=0; iy< eta_bins; iy++)
		{
	
			sprintf(name,"h_ChargedHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_pT[NoFile][iy] = new TH1D(name, "", 50,0,0.66*pTmax);
			if (useWeights) 	h_ChargedHadron_pT[NoFile][iy]->Sumw2();
				
			sprintf(name,"h_NeutralHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_pT[NoFile][iy] = new TH1D(name, "",  50,0,0.66*pTmax);
			if (useWeights) 	h_NeutralHadron_pT[NoFile][iy]->Sumw2();
	
			sprintf(name,"h_ChargedHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_SumpT[NoFile][iy] = new TH1D(name, "", 50,0,pTmax);
			if (useWeights) 	h_ChargedHadron_SumpT[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_SumpT[NoFile][iy] = new TH1D(name, "", 50,0,pTmax);
			if (useWeights) 	h_NeutralHadron_SumpT[NoFile][iy]->Sumw2();

			sprintf(name,"h_NumberOfChargedHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfChargedHadrons[NoFile][iy] = new TH1D(name, "", 50,0,50);
			if (useWeights) 	h_NumberOfChargedHadrons[NoFile][iy]->Sumw2();

			sprintf(name,"h_NumberOfNeutralHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfNeutralHadrons[NoFile][iy] = new TH1D(name, "",  50,0,50);
			if (useWeights) 	h_NumberOfNeutralHadrons[NoFile][iy]->Sumw2();

			sprintf(name,"h_ChargedHadron_DR_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_DR[NoFile][iy] = new TH1D(name, "",  50,0, 0.5);
			if (useWeights) 	h_ChargedHadron_DR[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_DR_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_DR[NoFile][iy] = new TH1D(name, "",  50,0, 0.5);
			if (useWeights) 	h_NeutralHadron_DR[NoFile][iy]->Sumw2();

			sprintf(name,"h_ChargedHadron_eta_phi_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_eta_phi[NoFile][iy] = new TH2D(name, "",  30,-0.5, 0.5, 30, -0.5, 0.5);
			if (useWeights) 	h_ChargedHadron_eta_phi[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_eta_phi_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_eta_phi[NoFile][iy] = new TH2D(name, "",  30,-0.5, 0.5, 30, -0.5, 0.5);
			if (useWeights) 	h_NeutralHadron_eta_phi[NoFile][iy]->Sumw2();

			sprintf(name,"h_girth_%s_bin%i",legend_array[NoFile],iy);
			h_girth[NoFile][iy] = new TH1D(name, "",  50, 0., 1.0 );	
			if (useWeights) 	h_girth[NoFile][iy]->Sumw2();

			sprintf(name,"h_quadratic_moment_%s_bin%i",legend_array[NoFile],iy);
			h_quadratic_moment[NoFile][iy] = new TH1D(name, "",  50, 0., 1.0 );	
			if (useWeights) 	h_quadratic_moment[NoFile][iy]->Sumw2();

		}
	}




	std::vector<float> *jpt = 0;
	std::vector<float> *eta = 0;
	std::vector<float> *phi = 0;
	std::vector<float> *nhf = 0;
	std::vector<float> *chf = 0;
	std::vector<float> *cemf = 0;
	std::vector<float> *nemf = 0;
	std::vector<float> *npr = 0;
	std::vector<float> *muf = 0;
	std::vector<float> *chMult = 0;
	std::vector<float> *neMult = 0;

	
	std::vector < std::vector<int> >  	 *jetconstituents;
	std::vector<int> *cand_vector = 0;

	std::vector<float> *npu = 0;
	std::vector<int> *PileupInteractions = 0;
	std::vector<int> *PileupOriginBX = 0;

	float weight;


	const int 	max_pfcand = 10000;
	int n_pfcand;
	float	pfcandpt[max_pfcand];
	float	pfcandeta[max_pfcand];
	float	pfcandphi[max_pfcand];
	float	pfcandm[max_pfcand];
	float	pfcandpdgid[max_pfcand];
	float	pfcandvertex[max_pfcand];

 	double SumEt,gen_SumEt;
 	int nVtx;


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		f[NoFile] = TFile::Open(input_files[NoFile],"READ"); //read the .root files
		tree[NoFile] = (TTree*)f[NoFile]->Get("dijetscouting/tree"); // get the trees from the files

		tree[NoFile]->SetBranchAddress("jpt",&jpt);
		tree[NoFile]->SetBranchAddress("eta",&eta);
		tree[NoFile]->SetBranchAddress("phi",&phi);
		tree[NoFile]->SetBranchAddress("nhf",&nhf);
		tree[NoFile]->SetBranchAddress("chf",&chf);
		tree[NoFile]->SetBranchAddress("cemf",&cemf);
		tree[NoFile]->SetBranchAddress("nemf",&nemf);
		tree[NoFile]->SetBranchAddress("npr",&npr);
		tree[NoFile]->SetBranchAddress("muf",&muf);
		tree[NoFile]->SetBranchAddress("chMult",&chMult);
		tree[NoFile]->SetBranchAddress("neMult",&neMult);
		tree[NoFile]->SetBranchAddress("jetconstituents",&jetconstituents);
		tree[NoFile]->SetBranchAddress("pfcandpt",&pfcandpt);
		tree[NoFile]->SetBranchAddress("pfcandeta",&pfcandeta);
		tree[NoFile]->SetBranchAddress("pfcandphi",&pfcandphi);
		tree[NoFile]->SetBranchAddress("pdcandm",&pfcandm);
		tree[NoFile]->SetBranchAddress("pfcandpdgid",&pfcandpdgid);
		tree[NoFile]->SetBranchAddress("pfcandvertex",&pfcandvertex);
		tree[NoFile]->SetBranchAddress("SumEt",&SumEt);
		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);

		tree[NoFile]->SetBranchAddress("SumEt",&SumEt);
		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);
		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);


		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		////////////////find the tree with the maximum entries to loop over them.

		if (flag_for_max_entries)  // get the entry of the first tree
		{
			maxEntries = treeEntries[NoFile];
			flag_for_max_entries = false;
		}
		if ( treeEntries[NoFile] > maxEntries ) maxEntries = treeEntries[NoFile];

		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
		cout << " Max entries = " << maxEntries << endl;


		int reco_size, gen_size;

//		for (int i=0; i<maxEntries; i++) //event loop
		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
	//	for (int i=0; i<5; i++) //event loop
		{

			if ( i >=  treeEntries[NoFile] ) continue; 

			if(i>0 && i%100000==0) cout<<"Event reached: "<<i<<endl;

			tree[NoFile]->GetEntry(i);

			reco_size = jpt->size();

			if(reco_size < 2) continue;	 									// skip event with less than 2 reco jets
			if ( fabs( phi->at(0) - phi->at(1) ) < deltaPhi ) continue; 	// skip event if the two leading jets are not back-to-back
			if ( jpt->at(0) < pTlowCut_leading || jpt->at(0) > pThighCut ) continue;// skip event if leading jet does not pass pT cuts
			if ( jpt->at(1) < pTlowCut_subleading || jpt->at(1) > pThighCut ) continue;// skip event if sub-leading jet does not pass pT cuts



//			if (npu->at(0)>80 || npu->at(0)<55) continue;





//			cout << "\n Processing event " << i << " number of reco_jets : "<< reco_size <<endl; 
			for ( int j=0; j<2; j++) // loop over two leading jets
			{

				int jet_ybin = getBin(fabs(eta->at(j)),yBnd, eta_bins); //get the jet eta bin

				int cand_size = ( jetconstituents->at( j ) ).size();
//				cout << "jet "<< reco_jets_matched_sequence[j] << " cand size = " << cand_size << endl;

				float chargedCand_SumpT = 0;
				float neutralCand_SumpT = 0;
				int neutralPFcand_number = 0;
				int chargedPFcand_number = 0;
				double jet_girth = 0;
				double jet_quadr_moment = 0;

				for (int l=0; l< cand_size; l++)
				{
					int jet_candidate = (jetconstituents->at(j)).at(l) ; //get the jet candidate number
					//cout << "\njet candidate number: " << jet_candidate << endl;

					int 	cand_pdgID = pfcandpdgid[jet_candidate];
					float 	cand_eta = pfcandeta[jet_candidate];
					float 	cand_phi = pfcandphi[jet_candidate];
					float 	cand_pT = pfcandpt[jet_candidate];
					double 	cand_DR = sqrt( pow( eta->at(j)-cand_eta ,2) +pow( phi->at(j)-cand_phi ,2) );

					jet_girth = jet_girth + (cand_pT/jpt->at(j))*cand_DR;
					jet_quadr_moment = jet_quadr_moment + (cand_pT/jpt->at(j))*cand_DR*cand_DR;
	//				cout << "\n candidate pdgID = " << cand_pdgID << endl;


					int cand_ybin = getBin( fabs(pfcandeta[ jet_candidate ]),yBnd, eta_bins );  //get the pf candidate eta bin

					if (jet_ybin > -1)//fill hist's in the corresponding eta bin
					{

						if ( cand_pdgID==-211 || cand_pdgID==211 || cand_pdgID==11 || cand_pdgID==-11 || cand_pdgID==13 || cand_pdgID==-13) chargedPFcand_number = chargedPFcand_number +1;
						if ( cand_pdgID==-211 || cand_pdgID==211 )
						{
							chargedCand_SumpT = chargedCand_SumpT + cand_pT; // get the pT sum of charged candidates
							if(!useWeights)
							{

								h_ChargedHadron_pT[NoFile][jet_ybin]->Fill(cand_pT); // fill charged hadron pT's
								h_ChargedHadron_DR[NoFile][jet_ybin]->Fill( cand_DR );
								h_ChargedHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(j), cand_phi-phi->at(j)  );
							}
							else
							{	
								h_ChargedHadron_pT[NoFile][jet_ybin]->Fill(cand_pT, weight); // fill charged hadron pT's with weight
								h_ChargedHadron_DR[NoFile][jet_ybin]->Fill( cand_DR, weight );
								h_ChargedHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(j), cand_phi-phi->at(j), weight );
							}
						}
						if ( cand_pdgID==130 || cand_pdgID==22 || cand_pdgID==0 ) 	neutralPFcand_number = neutralPFcand_number +1;
						if ( cand_pdgID==130 )
						{
							neutralCand_SumpT = neutralCand_SumpT + cand_pT; // get the pT sum of charged candidates
							if(!useWeights)
							{
								h_NeutralHadron_pT[NoFile][jet_ybin]->Fill(cand_pT); // fill neutral hadron pT's
								h_NeutralHadron_DR[NoFile][jet_ybin]->Fill( cand_DR );
								h_NeutralHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(j), cand_phi-phi->at(j)  );
							}
							else
							{
								h_NeutralHadron_pT[NoFile][jet_ybin]->Fill(cand_pT, weight); // fill neutral hadron pT's with weight
								h_NeutralHadron_DR[NoFile][jet_ybin]->Fill( cand_DR, weight );
								h_NeutralHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(j), cand_phi-phi->at(j), weight  );
							}
						}
					}
			//		cout << " cand " << cand_vector->at(l) << " , pt = " << pfcandpt[cand_vector->at(l)] << endl;  
				} // end of candidate loop

				if(!useWeights)
				{
					h_NumberOfNeutralHadrons[NoFile][jet_ybin]->Fill(neutralPFcand_number);
					h_NumberOfChargedHadrons[NoFile][jet_ybin]->Fill(chargedPFcand_number);
					h_ChargedHadron_SumpT[NoFile][jet_ybin]->Fill(chargedCand_SumpT);
					h_NeutralHadron_SumpT[NoFile][jet_ybin]->Fill(neutralCand_SumpT);
					h_girth[NoFile][jet_ybin]->Fill(jet_girth);
					h_quadratic_moment[NoFile][jet_ybin]->Fill(jet_quadr_moment);
				}
				else 
				{
					h_NumberOfNeutralHadrons[NoFile][jet_ybin]->Fill(neutralPFcand_number, weight);
					h_NumberOfChargedHadrons[NoFile][jet_ybin]->Fill(chargedPFcand_number, weight);
					h_ChargedHadron_SumpT[NoFile][jet_ybin]->Fill(chargedCand_SumpT, weight);
					h_NeutralHadron_SumpT[NoFile][jet_ybin]->Fill(neutralCand_SumpT, weight);
					h_girth[NoFile][jet_ybin]->Fill(jet_girth, weight);
					h_quadratic_moment[NoFile][jet_ybin]->Fill(jet_quadr_moment, weight);
				}
			} //  end of jet loop


			reco_size = 0;


		} // end of event loop


		delete tree[NoFile];
		delete f[NoFile];
	} // end of file loop



// =========================================================== Plotting staff ===================================================================

	TCanvas *pad1 = new TCanvas("pad1", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad1->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad2 = new TCanvas("pad2", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad2->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad3 = new TCanvas("pad3", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad3->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad4 = new TCanvas("pad4", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad4->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad5 = new TCanvas("pad5", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad5->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad6 = new TCanvas("pad6", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad6->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_ChargedCand_DR = new TCanvas("pad_ChargedCand_DR", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_ChargedCand_DR->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_NeutralCand_DR = new TCanvas("pad_NeutralCand_DR", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_NeutralCand_DR->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_girth = new TCanvas("pad_girth", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_girth->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_quadr_moment = new TCanvas("pad_quadr_moment", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_quadr_moment->Divide(PadColumnsEtaBins, PadRowsEtaBins);


	TCanvas *pad_ChargedHadron_eta_phi[NoFiles], *pad_NeutralHadron_eta_phi[NoFiles];

	TPaveText *paveCMS = new TPaveText(0.45,0.95,0.5,1.0,"NDC");
	// paveCMS->AddText("CMS Preliminary L=9.2 fb^{-1} #sqrt{s} = 13 TeV");
//	paveCMS->AddText("CMS Preliminary#sqrt{s} = 13 TeV");
	paveCMS->AddText("CMS Simulation #sqrt{s} = 14 TeV");
	paveCMS->SetFillColor(0);
	paveCMS->SetBorderSize(0);
	paveCMS->SetTextSize(0.04);

	TLegend *leg1 =new TLegend(0.2,0.4,0.9,0.9); //7899  //4899
	leg1->SetTextSize(0.055);
	leg1->SetFillColor(0); 
	leg1->SetBorderSize(0);  


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		sprintf(name,"pad_ChargedHadron_eta_phi_%s",legend_array[NoFile]);
		pad_ChargedHadron_eta_phi[NoFile] = new TCanvas(name, name , Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
		pad_ChargedHadron_eta_phi[NoFile]->Divide(PadColumnsEtaBins, PadRowsEtaBins);

		sprintf(name,"pad_NeutralHadron_eta_phi_%s",legend_array[NoFile]);
		pad_NeutralHadron_eta_phi[NoFile] = new TCanvas(name, name , Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
		pad_NeutralHadron_eta_phi[NoFile]->Divide(PadColumnsEtaBins, PadRowsEtaBins);

 
		leg1->AddEntry(h_ChargedHadron_SumpT[NoFile][0], legend_array[NoFile],  "L");

		for(int iy=0; iy<eta_bins; iy++)
		{

			double etamin = yBnd[iy];
			double etamax = yBnd[iy+1];
			const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
			TLatex *teta = new TLatex(0.38,0.86,seta); //cout<<seta<<endl;
			teta->SetNDC();
			teta->SetTextSize(0.045);

			if(eta_bins > 1) pad1->cd(iy+1);
			else pad1->cd();
			if(eta_bins > 1) pad1->cd(iy+1)->SetLogy(1);
			else pad1->cd()->SetLogy(1);
			h_ChargedHadron_SumpT[NoFile][iy]->GetXaxis()->SetTitle("Sum pT of charged hadrons per jet (GeV) ");
			h_ChargedHadron_SumpT[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_ChargedHadron_SumpT[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_ChargedHadron_SumpT[NoFile][iy]->SetLineWidth(2);
		//	h_ChargedHadron_SumpT[NoFile][iy]->SetMinimum(0.1);
		//	h_ChargedHadron_SumpT[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_ChargedHadron_SumpT[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 1000000);
			else			h_ChargedHadron_SumpT[NoFile][iy]->GetYaxis()->SetRangeUser(0.001*h_ChargedHadron_SumpT[NoFile][iy]->GetMaximum(), 10*h_ChargedHadron_SumpT[NoFile][iy]->GetMaximum());
			h_ChargedHadron_SumpT[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_ChargedHadron_SumpT[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_ChargedHadron_SumpT[NoFile][iy]->Scale(h_ChargedHadron_SumpT[0][iy]->Integral()/h_ChargedHadron_SumpT[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_ChargedHadron_SumpT[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_ChargedHadron_SumpT[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");



			if(eta_bins > 1) pad2->cd(iy+1);
			else pad2->cd();
			if(eta_bins > 1) pad2->cd(iy+1)->SetLogy(1);
			else pad2->cd()->SetLogy(1);
			h_NeutralHadron_SumpT[NoFile][iy]->GetXaxis()->SetTitle("Sum pT of neutral hadrons per jet (GeV) ");
			h_NeutralHadron_SumpT[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_NeutralHadron_SumpT[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NeutralHadron_SumpT[NoFile][iy]->SetLineWidth(2);
		//	h_NeutralHadron_SumpT[NoFile][iy]->SetMinimum(0.1);
		//	h_NeutralHadron_SumpT[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_NeutralHadron_SumpT[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 1000000);
			else			h_NeutralHadron_SumpT[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_NeutralHadron_SumpT[NoFile][iy]->GetMaximum(), 10*h_NeutralHadron_SumpT[NoFile][iy]->GetMaximum());
			h_NeutralHadron_SumpT[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NeutralHadron_SumpT[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NeutralHadron_SumpT[NoFile][iy]->Scale(h_NeutralHadron_SumpT[0][iy]->Integral()/h_NeutralHadron_SumpT[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NeutralHadron_SumpT[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NeutralHadron_SumpT[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			if(eta_bins > 1) pad3->cd(iy+1);
			else pad3->cd();
			if(eta_bins > 1) pad3->cd(iy+1)->SetLogy(1);
			else pad3->cd()->SetLogy(1);
			h_ChargedHadron_pT[NoFile][iy]->GetXaxis()->SetTitle("pT of charged hadrons (GeV) ");
			h_ChargedHadron_pT[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_ChargedHadron_pT[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_ChargedHadron_pT[NoFile][iy]->SetLineWidth(2);
		//	h_ChargedHadron_pT[NoFile][iy]->SetMinimum(0.1);
		//	h_ChargedHadron_pT[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_ChargedHadron_pT[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 1000000);
			else			h_ChargedHadron_pT[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_ChargedHadron_pT[NoFile][iy]->GetMaximum(), 10*h_ChargedHadron_pT[NoFile][iy]->GetMaximum());
			h_ChargedHadron_pT[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_ChargedHadron_pT[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_ChargedHadron_pT[NoFile][iy]->Scale(h_ChargedHadron_pT[0][iy]->Integral()/h_ChargedHadron_pT[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_ChargedHadron_pT[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_ChargedHadron_pT[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");



			if(eta_bins > 1) pad4->cd(iy+1);
			else pad4->cd();
			if(eta_bins > 1) pad4->cd(iy+1)->SetLogy(1);
			else pad4->cd()->SetLogy(1);
			h_NeutralHadron_pT[NoFile][iy]->GetXaxis()->SetTitle("pT of neutral hadrons (GeV) ");
			h_NeutralHadron_pT[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_NeutralHadron_pT[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NeutralHadron_pT[NoFile][iy]->SetLineWidth(2);
		//	h_NeutralHadron_pT[NoFile][iy]->SetMinimum(0.1);
		//	h_NeutralHadron_pT[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_NeutralHadron_pT[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 1000000);
			else			h_NeutralHadron_pT[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_NeutralHadron_pT[NoFile][iy]->GetMaximum(), 10*h_NeutralHadron_pT[NoFile][iy]->GetMaximum());
			h_NeutralHadron_pT[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NeutralHadron_pT[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NeutralHadron_pT[NoFile][iy]->Scale(h_NeutralHadron_pT[0][iy]->Integral()/h_NeutralHadron_pT[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NeutralHadron_pT[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NeutralHadron_pT[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			


			if(eta_bins > 1) pad5->cd(iy+1);
			else pad5->cd();
			if(eta_bins > 1) pad5->cd(iy+1)->SetLogy(1);
			else pad5->cd()->SetLogy(1);
			h_NumberOfChargedHadrons[NoFile][iy]->GetXaxis()->SetTitle("# of charged candidates per jet ");
			h_NumberOfChargedHadrons[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_NumberOfChargedHadrons[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NumberOfChargedHadrons[NoFile][iy]->SetLineWidth(2);
		//	h_NumberOfChargedHadrons[NoFile][iy]->SetMinimum(0.1);
		//	h_NumberOfChargedHadrons[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_NumberOfChargedHadrons[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 1000000);
			else			h_NumberOfChargedHadrons[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_NumberOfChargedHadrons[NoFile][iy]->GetMaximum(), 10*h_NumberOfChargedHadrons[NoFile][iy]->GetMaximum());
			h_NumberOfChargedHadrons[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NumberOfChargedHadrons[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NumberOfChargedHadrons[NoFile][iy]->Scale(h_NumberOfChargedHadrons[0][iy]->Integral()/h_NumberOfChargedHadrons[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NumberOfChargedHadrons[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NumberOfChargedHadrons[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");




			if(eta_bins > 1) pad6->cd(iy+1);
			else pad6->cd();
			if(eta_bins > 1) pad6->cd(iy+1)->SetLogy(1);
			else pad6->cd()->SetLogy(1);
			h_NumberOfNeutralHadrons[NoFile][iy]->GetXaxis()->SetTitle("# of neutral candidates per jet");
			h_NumberOfNeutralHadrons[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_NumberOfNeutralHadrons[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NumberOfNeutralHadrons[NoFile][iy]->SetLineWidth(2);
		//	h_NumberOfNeutralHadrons[NoFile][iy]->SetMinimum(0.1);
		//	h_NumberOfNeutralHadrons[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_NumberOfNeutralHadrons[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 1000000);
			else			h_NumberOfNeutralHadrons[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_NumberOfNeutralHadrons[NoFile][iy]->GetMaximum(), 10*h_NumberOfNeutralHadrons[NoFile][iy]->GetMaximum());
			h_NumberOfNeutralHadrons[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NumberOfNeutralHadrons[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NumberOfNeutralHadrons[NoFile][iy]->Scale(h_NumberOfNeutralHadrons[0][iy]->Integral()/h_NumberOfNeutralHadrons[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NumberOfNeutralHadrons[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NumberOfNeutralHadrons[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			if(eta_bins > 1) pad_ChargedCand_DR->cd(iy+1);
			else pad_ChargedCand_DR->cd();
			if(eta_bins > 1) pad_ChargedCand_DR->cd(iy+1)->SetLogy(1);
			else pad_ChargedCand_DR->cd()->SetLogy(1);
			h_ChargedHadron_DR[NoFile][iy]->GetXaxis()->SetTitle("Jet-ChargedCandidate DR");
			h_ChargedHadron_DR[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_ChargedHadron_DR[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_ChargedHadron_DR[NoFile][iy]->SetLineWidth(2);
		//	h_ChargedHadron_DR[NoFile][iy]->SetMinimum(0.1);
		//	h_ChargedHadron_DR[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_ChargedHadron_DR[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 100000);
			else			h_ChargedHadron_DR[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_ChargedHadron_DR[NoFile][iy]->GetMaximum(), 10*h_ChargedHadron_DR[NoFile][iy]->GetMaximum());
			h_ChargedHadron_DR[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_ChargedHadron_DR[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_ChargedHadron_DR[NoFile][iy]->Scale(h_ChargedHadron_DR[0][iy]->Integral()/h_ChargedHadron_DR[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_ChargedHadron_DR[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_ChargedHadron_DR[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			if(eta_bins > 1) pad_NeutralCand_DR->cd(iy+1);
			else pad_NeutralCand_DR->cd();
			if(eta_bins > 1) pad_NeutralCand_DR->cd(iy+1)->SetLogy(1);
			else pad_NeutralCand_DR->cd()->SetLogy(1);
			h_NeutralHadron_DR[NoFile][iy]->GetXaxis()->SetTitle("Jet-NeutralCandidate DR");
			h_NeutralHadron_DR[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_NeutralHadron_DR[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NeutralHadron_DR[NoFile][iy]->SetLineWidth(2);
		//	h_ChargedHadron_DR[NoFile][iy]->SetMinimum(0.1);
		//	h_ChargedHadron_DR[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_NeutralHadron_DR[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 100000);
			else			h_NeutralHadron_DR[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_NeutralHadron_DR[NoFile][iy]->GetMaximum(), 10*h_NeutralHadron_DR[NoFile][iy]->GetMaximum());
			h_NeutralHadron_DR[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NeutralHadron_DR[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NeutralHadron_DR[NoFile][iy]->Scale(h_NeutralHadron_DR[0][iy]->Integral()/h_NeutralHadron_DR[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NeutralHadron_DR[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NeutralHadron_DR[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			if(eta_bins > 1) pad_girth->cd(iy+1);
			else pad_girth->cd();
			if(eta_bins > 1) pad_girth->cd(iy+1)->SetLogy(1);
			else pad_girth->cd()->SetLogy(1);
			h_girth[NoFile][iy]->GetXaxis()->SetTitle("Jet-girth");
			h_girth[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_girth[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_girth[NoFile][iy]->SetLineWidth(2);
		//	h_ChargedHadron_DR[NoFile][iy]->SetMinimum(0.1);
		//	h_ChargedHadron_DR[NoFile][iy]->SetMaximum(100000);
			if(!useWeights)	h_girth[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 100000);
			else			h_girth[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_girth[NoFile][iy]->GetMaximum(), 10*h_girth[NoFile][iy]->GetMaximum());
			h_girth[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_girth[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_girth[NoFile][iy]->Scale(h_girth[0][iy]->Integral()/h_girth[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_girth[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_girth[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");



			if(eta_bins > 1) pad_quadr_moment->cd(iy+1);
			else pad_quadr_moment->cd();
			if(eta_bins > 1) pad_quadr_moment->cd(iy+1)->SetLogy(1);
			else pad_quadr_moment->cd()->SetLogy(1);
			h_quadratic_moment[NoFile][iy]->GetXaxis()->SetTitle("Quadratic Radial Moment");
			h_quadratic_moment[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_quadratic_moment[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_quadratic_moment[NoFile][iy]->SetLineWidth(2);
			if(!useWeights)	h_quadratic_moment[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 100000);
			else			h_quadratic_moment[NoFile][iy]->GetYaxis()->SetRangeUser(0.01*h_quadratic_moment[NoFile][iy]->GetMaximum(), 10*h_quadratic_moment[NoFile][iy]->GetMaximum());
			h_quadratic_moment[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_quadratic_moment[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_quadratic_moment[NoFile][iy]->Scale(h_quadratic_moment[0][iy]->Integral()/h_quadratic_moment[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_quadratic_moment[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_quadratic_moment[NoFile][iy]->Draw("same hist");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");



//============================================== 2D plots ==============================================

			if(eta_bins > 1) pad_ChargedHadron_eta_phi[NoFile]->cd(iy+1);
			else pad_ChargedHadron_eta_phi[NoFile]->cd();
			if(eta_bins > 1)
			{
				pad_ChargedHadron_eta_phi[NoFile]->cd(iy+1)->SetLogz(1);
				pad_ChargedHadron_eta_phi[NoFile]->cd(iy+1)->SetRightMargin(0.18);
				pad_ChargedHadron_eta_phi[NoFile]->cd(iy+1)->SetTopMargin(0.07);
			}
			else
			{
				pad_ChargedHadron_eta_phi[NoFile]->cd()->SetLogz(1);
				pad_ChargedHadron_eta_phi[NoFile]->cd()->SetRightMargin(0.18);
				pad_ChargedHadron_eta_phi[NoFile]->cd()->SetTopMargin(0.07);
			}
			
			h_ChargedHadron_eta_phi[NoFile][iy]->SetStats(true);
		    

			h_ChargedHadron_eta_phi[NoFile][iy]->SetTitle("ChargedCandidate eta-phi");
			h_ChargedHadron_eta_phi[NoFile][iy]->GetXaxis()->SetTitle("#eta_{cand} - #eta_{jet}");
			h_ChargedHadron_eta_phi[NoFile][iy]->GetYaxis()->SetTitle("#phi_{cand} - #phi_{jet}");
			h_ChargedHadron_eta_phi[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_ChargedHadron_eta_phi[NoFile][iy]->SetMinimum(0.001*h_ChargedHadron_eta_phi[NoFile][iy]->GetMaximum() );
			h_ChargedHadron_eta_phi[NoFile][iy]->SetMaximum(3*h_ChargedHadron_eta_phi[NoFile][iy]->GetMaximum());

			h_ChargedHadron_eta_phi[NoFile][iy]->Draw("colz");
			paveCMS ->Draw("same");
			teta->Draw();



			if(eta_bins > 1) pad_NeutralHadron_eta_phi[NoFile]->cd(iy+1);
			else pad_NeutralHadron_eta_phi[NoFile]->cd();
			if(eta_bins > 1)
			{
				pad_NeutralHadron_eta_phi[NoFile]->cd(iy+1)->SetLogz(1);
				pad_NeutralHadron_eta_phi[NoFile]->cd(iy+1)->SetRightMargin(0.18);
				pad_NeutralHadron_eta_phi[NoFile]->cd(iy+1)->SetTopMargin(0.07);
			}
			else
			{
				pad_NeutralHadron_eta_phi[NoFile]->cd()->SetLogz(1);
				pad_NeutralHadron_eta_phi[NoFile]->cd()->SetRightMargin(0.18);
				pad_NeutralHadron_eta_phi[NoFile]->cd()->SetTopMargin(0.07);
			}
			
			h_NeutralHadron_eta_phi[NoFile][iy]->SetStats(true);
		    

			h_NeutralHadron_eta_phi[NoFile][iy]->SetTitle("NeutralCandidate eta-phi");
			h_NeutralHadron_eta_phi[NoFile][iy]->GetXaxis()->SetTitle("#eta_{cand} - #eta_{jet}");
			h_NeutralHadron_eta_phi[NoFile][iy]->GetYaxis()->SetTitle("#phi_{cand} - #phi_{jet}");
			h_NeutralHadron_eta_phi[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NeutralHadron_eta_phi[NoFile][iy]->SetMinimum(0.001*h_NeutralHadron_eta_phi[NoFile][iy]->GetMaximum() );
			h_NeutralHadron_eta_phi[NoFile][iy]->SetMaximum(3*h_NeutralHadron_eta_phi[NoFile][iy]->GetMaximum());

			h_NeutralHadron_eta_phi[NoFile][iy]->Draw("colz");
			paveCMS ->Draw("same");
			teta->Draw();






			if(eta_bins > 1 &&  NoFile == NoFiles-1 )
			{
				pad1->cd(eta_bins+1);    			leg1->Draw();
				pad2->cd(eta_bins+1);    			leg1->Draw();
				pad3->cd(eta_bins+1);    			leg1->Draw();
				pad4->cd(eta_bins+1);   			leg1->Draw();
				pad5->cd(eta_bins+1);   			leg1->Draw();
				pad6->cd(eta_bins+1);   			leg1->Draw();
				pad_girth->cd(eta_bins+1);   		leg1->Draw();
				pad_quadr_moment->cd(eta_bins+1);	leg1->Draw();

				pad_ChargedCand_DR->cd(eta_bins+1);	leg1->Draw();
				pad_NeutralCand_DR->cd(eta_bins+1);	leg1->Draw();
			}

		} //end of loop on eta bins
	} // end of loop on files

/*
		pad1->SaveAs("Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/Comparisons_PatatrackVsLegacyVsFull_ChargedPFcand_SumpT.png");
		pad2->SaveAs("Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/Comparisons_PatatrackVsLegacyVsFull_NeutralPFcand_SumpT.png");
		pad3->SaveAs("Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/Comparisons_PatatrackVsLegacyVsFull_ChargedPFcand_pT.png");
		pad4->SaveAs("Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/Comparisons_PatatrackVsLegacyVsFull_NeutralPFcand_pT.png");
		pad5->SaveAs("Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/Comparisons_PatatrackVsLegacyVsFull_ChargedPFcand_number.png");
		pad6->SaveAs("Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/Comparisons_PatatrackVsLegacyVsFull_NeutralPFcand_number.png");
*/

	char filename[1500]; 
	if(Save_Plots)
	{ 
		sprintf(filename,"%s/%s/%s_ChargedPFcand_SumpT.png",analyzer_path,output_directory,image_name);
		pad1->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_SumpT.png",analyzer_path,output_directory,image_name);
		pad2->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_ChargedPFcand_pT.png",analyzer_path,output_directory,image_name);
		pad3->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_pT.png",analyzer_path,output_directory,image_name);
		pad4->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_ChargedPFcand_number.png",analyzer_path,output_directory,image_name);
		pad5->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_number.png",analyzer_path,output_directory,image_name);
		pad6->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_ChargedPFcand_DR.png",analyzer_path,output_directory,image_name);
		pad_ChargedCand_DR->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_DR.png",analyzer_path,output_directory,image_name);
		pad_NeutralCand_DR->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_jet_girth.png",analyzer_path,output_directory,image_name);
		pad_girth->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_jet_quadratic_moment.png",analyzer_path,output_directory,image_name);
		pad_quadr_moment->SaveAs(filename);

		for (int NoFile=0; NoFile<NoFiles; NoFile++)
		{
			sprintf(filename,"%s/%s/%s_ChargedPFcand_eta_phi_%s.png",analyzer_path,output_directory,image_name,legend_array[NoFile]);
			pad_ChargedHadron_eta_phi[NoFile]->SaveAs(filename);
		
			sprintf(filename,"%s/%s/%s_NeutralPFcand_eta_phi_%s.png",analyzer_path,output_directory,image_name,legend_array[NoFile]);
			pad_NeutralHadron_eta_phi[NoFile]->SaveAs(filename);
		}

	}

}
	
