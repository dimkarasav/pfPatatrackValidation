#include "include_functions.h"


void Tracking_Comparisons_btb()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);



//	double yBnd[]={0.0, 1.3, 2.4};
	double yBnd[]={0.0, 1.3, 2.4, 2.7, 3.0};  




	char input_files[][800] ={


//mkFit comparisons
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_LooseTrackParams_noPU_onlyZeta.root" ,
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/mkFit_Candidates_ttbar_noPU.root"
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/PatatrackPixels_Candidates_ttbar_14TeV_noPU.root"

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
//optimized vs full noPU
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_noPU_linkFix_pTError1_TrackAlgoCut2_Zeta0p5.root"

//optimized vs full PU
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_PU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_ttbar_iteration1.root"
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FixTrackCaloLink/Patatrack_ttbar_PU_linkFix_optimized.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_ttbar_iteration1.root"
//"delme_ttbar.root"


//optimized vs full tracking QCD
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/QCD_samples/FullTracking_QCD_PU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/QCD_samples/Patatrack_QCD_PU.root"

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
//no PU comparisons
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

	double pTlowCut_leading = 60;
	double pTlowCut_subleading = 30;
	double pThighCut = 5000;
	double pThistoMax = 450;

	int PadColumns = 3;
	int PadRows = 2;
	int Canvas_Xpixels = PadColumns*333;
	int Canvas_Ypixels = PadRows*500;

	char analyzer_path[500] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre7/src/Jakobs_producer/ProducerTest/"}; 
//	char *output_directory = "Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/JetpT0_70"; 
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/TrackCaloFix/ScanZetaCut/"};   
   	char output_directory[200] = {"Realistic_RunIII_14TeV/QCD_sample/Weighted_btb/"};          
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/PUvsNoPU//"};          
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/mkFit/FullVsmkFit"};              
//	char output_directory[200] = {"deleteme/"};             
	char image_name[200] = {"btb_balanced"}; 


	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);
	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;
	char legend_array[NoFiles][500] = { "FullTracking" ,  "PatatrackPixelTracks"  };
//	char legend_array[NoFiles][500] = { "MaxIteration_50" ,  "MaxIteration_1"  };
//	char legend_array[NoFiles][500] = { "FullTracking" , "TightCuts", "ZetaCut0p2", "ZetaCut0p5", "ZetaCut1p0", "ZetaCut2p0", "LooseCuts" };
//	char legend_array[NoFiles][500] = { "FullTracking" ,  "mkFitTracking"  };
//	char legend_array[NoFiles][500] = { "nominal_PU" , "nominal_noPU","LooseCuts_PU" , "LooseCuts_noPU"   };
//	char legend_array[NoFiles][500] = { "FullTracking" , "ptError90", "ptError45", "ptError10"  };
//	char legend_array[NoFiles][500] = { "FullTracking" , "LooseCuts", "TrackAlgoCut2", "TrackAlgoCut1" ,"TrackCut1-ptError10", "TightCuts"  };


	int Colors[] = { 1, 4, 2 , 6, 3, 7 , 28, 46} ; // black, blue, red , magenta, green, light blue ,  brown, redish
//	int Colors[] = { 1, 4, 2 , 6 } ; // black, blue, red , magenta
	bool scale_histos = false;
	bool Save_Plots = true ;
	bool useWeights = true ;
	bool plot_minDR = true;
	double deltaPhi = 2.7;
	double DR_threshold = 0.2;


	char eta_bins_legend[eta_bins][25];



	for (int iy=0; iy< eta_bins; iy++)
	{
		if (iy==0) sprintf( eta_bins_legend[iy], "|#eta|<%3.1f" , yBnd[iy+1] );  
		else
		{
			sprintf( eta_bins_legend[iy], "%3.1f<|#eta|<%3.1f" , yBnd[iy] ,yBnd[iy+1] );
		}
	}


	TH1D *h_METovSUMET[NoFiles][eta_bins], *h_CHFJet[NoFiles][eta_bins], *h_NHFJet[NoFiles][eta_bins], *h_CEMFJet[NoFiles][eta_bins], *h_NEMFJet[NoFiles][eta_bins], *h_MUFJet[NoFiles][eta_bins], *h_CMJet[NoFiles][eta_bins], *h_NMJet[NoFiles][eta_bins], *h_ptJet[NoFiles][eta_bins], *h_PHIJet[NoFiles][eta_bins], *h_pT_res[NoFiles][eta_bins], *h_gen_pT[NoFiles][eta_bins],*h_gen_pT_all[NoFiles][eta_bins], *h_reco_pT_unmatched[NoFiles][eta_bins], *h_reco_pT_all[NoFiles][eta_bins];

	

	TH1D *h_ETAJet[NoFiles], *h_SumEt[NoFiles], *h_minDR[NoFiles];

	TGraphAsymmErrors *Reco_eff[NoFiles][eta_bins], *Fake_rate[NoFiles][eta_bins]; 



	char filename[256]; 
	char name[256]; 

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko();


for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		sprintf(name,"h_ETAJet_%s",legend_array[NoFile]);
		h_ETAJet[NoFile] = new TH1D(name, "", 60,-5,5);

	
		sprintf(name,"h_SumEt_%s",legend_array[NoFile]);
		h_SumEt[NoFile] = new TH1D(name, "", 100, 0, 7000); // 40,0,1.0
	
		sprintf(name,"h_minDR_%s",legend_array[NoFile]);
		h_minDR[NoFile] = new TH1D(name, "", 100, 0.0, 5.0); // 40,0,1.0

		if (useWeights)
		{
			h_ETAJet[NoFile]->Sumw2();
			h_SumEt[NoFile]->Sumw2();
		}
	

		for(Int_t h=0; h<eta_bins;h++)
		{ 
		//========== reco jets===========
			sprintf(name,"h_METovSUMET_%s_bin%i",legend_array[NoFile],h);
			h_METovSUMET[NoFile][h] = new TH1D(name, "", 50, 0, 1); // 40,0,1.0
			if (useWeights) h_METovSUMET[NoFile][h]->Sumw2();

			sprintf(name,"h_CHFJet_%s_bin%i",legend_array[NoFile],h);
			h_CHFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2); //40, 0,1.2
			if (useWeights) h_CHFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_NHFJet_%s_bin%i",legend_array[NoFile],h);
			h_NHFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2);
			if (useWeights) h_NHFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_CEMFJet_%s_bin%i",legend_array[NoFile],h);
			h_CEMFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2);
			if (useWeights) h_CEMFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_NEMFJet_%s_bin%i",legend_array[NoFile],h);
			h_NEMFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2);
			if (useWeights) h_NEMFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_MUFJet_%s_bin%i",legend_array[NoFile],h);
			h_MUFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2);
			if (useWeights) h_MUFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_CMJet_%s_bin%i",legend_array[NoFile],h);
			h_CMJet[NoFile][h] = new TH1D(name, "", 80, 0.0-0.5, 80.0 - 0.5);
		//	if (useWeights) h_CMJet[NoFile][h]->Sumw2();
		
			sprintf(name,"h_PHIJet_%s_bin%i",legend_array[NoFile],h);
			h_PHIJet[NoFile][h] = new TH1D(name, "", 40,-4,4);
			if (useWeights) h_PHIJet[NoFile][h]->Sumw2();

			sprintf(name,"h_NMJet_%s_bin%i",legend_array[NoFile],h);
			h_NMJet[NoFile][h] = new TH1D(name, "", 80, 0.-0.5, 80.-0.5);
			if (useWeights) h_NMJet[NoFile][h]->Sumw2();

			sprintf(name,"h_ptJet_%s_bin%i",legend_array[NoFile],h);
			h_ptJet[NoFile][h] = new TH1D(name, "", 50,0,pThistoMax);
			if (useWeights) h_ptJet[NoFile][h]->Sumw2();

			sprintf(name,"h_pT_res_%s_bin%i",legend_array[NoFile],h);
			h_pT_res[NoFile][h] = new TH1D(name, "", 50, 0., 3.5); // 40,0,1.0
			if (useWeights) h_pT_res[NoFile][h]->Sumw2();

			sprintf(name,"h_gen_pT_all_%s_bin%i",legend_array[NoFile],h);
			h_gen_pT_all[NoFile][h] = new TH1D(name, "", 50, 0., pThistoMax); // 40,0,1.0
			if (useWeights) h_gen_pT_all[NoFile][h]->Sumw2();

			sprintf(name,"h_gen_pT_%s_bin%i",legend_array[NoFile],h);
			h_gen_pT[NoFile][h] = new TH1D(name, "", 50, 0., pThistoMax); // 40,0,1.0
			if (useWeights) h_gen_pT[NoFile][h]->Sumw2();

			sprintf(name,"h_reco_pT_unmatched_%s_bin%i",legend_array[NoFile],h);
			h_reco_pT_unmatched[NoFile][h] = new TH1D(name, "", 50, 0., pThistoMax); // 40,0,1.0
			if (useWeights) h_reco_pT_unmatched[NoFile][h]->Sumw2();

			sprintf(name,"h_reco_pT_all_%s_bin%i",legend_array[NoFile],h);
			h_reco_pT_all[NoFile][h] = new TH1D(name, "", 50, 0., pThistoMax); // 40,0,1.0
			if (useWeights) h_reco_pT_all[NoFile][h]->Sumw2();

		}// end of etabin loop
	} // end of file loop

	int treeEntries[NoFiles];
	TFile *f[NoFiles];
	TTree *tree[NoFiles];


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


	std::vector<float> *npu = 0;
	std::vector<int> *PileupInteractions = 0;
	std::vector<int> *PileupOriginBX = 0;

	float weight;
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

		tree[NoFile]->SetBranchAddress("SumEt",&SumEt);
		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);
		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);



		int size, reco_size;

		//======================= f1 jets ===================================

		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
//		for (int i=0; i<200000; i++) //event loop
		{
			tree[NoFile]->GetEntry(i);
			reco_size = jpt->size();

			if(reco_size < 2) continue;	 									// skip event with less than 2 jets
			if ( fabs( phi->at(0) - phi->at(1) ) < deltaPhi ) continue; 	// skip event if the two leading jets are not back-to-back
			if ( jpt->at(0) < pTlowCut_leading || jpt->at(0) > pThighCut ) continue;// skip event if leading jet does not pass pT cuts
			if ( jpt->at(1) < pTlowCut_subleading || jpt->at(1) > pThighCut ) continue;// skip event if sub-leading jet does not pass pT cuts
			if ( jpt->at(0)/jpt->at(1)>1.3 ) continue; // skip events with 2 leading jets more than 30% unbalanced in terms of pt 
			
			for ( int j=0; j<2; j++) 
			{

				if (j==0 && !useWeights) h_SumEt[NoFile]->Fill(SumEt);
				if (j==0 &&  useWeights) h_SumEt[NoFile]->Fill(SumEt, weight);

				if (!useWeights)	h_ETAJet[NoFile]->Fill( eta->at(j) );
				else 				h_ETAJet[NoFile]->Fill( eta->at(j), weight );
				int ybin = getBin(fabs(eta->at(j)), yBnd, eta_bins);
				if (ybin > -1)//fill hist's in the corresponding eta bin
				{
					if (!useWeights)
					{
						h_CHFJet[NoFile][ybin] ->Fill( chf->at(j)  ); 
						h_NHFJet[NoFile][ybin] ->Fill( nhf->at(j)  );
						h_CEMFJet[NoFile][ybin]->Fill( cemf->at(j) );
						h_NEMFJet[NoFile][ybin]->Fill( nemf->at(j) );
						h_MUFJet[NoFile][ybin] ->Fill( muf->at(j)  );
						h_PHIJet[NoFile][ybin] ->Fill( phi->at(j)  );
						h_ptJet[NoFile][ybin]  ->Fill( jpt->at(j)  );
						h_CMJet[NoFile][ybin]  ->Fill( chMult->at(j) );
						h_NMJet[NoFile][ybin]  ->Fill( neMult->at(j) );
					}
					else 
					{
						h_CHFJet[NoFile][ybin] ->Fill( chf->at(j),weight ); 
						h_NHFJet[NoFile][ybin] ->Fill( nhf->at(j),weight  );
						h_CEMFJet[NoFile][ybin]->Fill( cemf->at(j),weight );
						h_NEMFJet[NoFile][ybin]->Fill( nemf->at(j),weight );
						h_MUFJet[NoFile][ybin] ->Fill( muf->at(j),weight  );
						h_PHIJet[NoFile][ybin] ->Fill( phi->at(j),weight  );
						h_ptJet[NoFile][ybin]  ->Fill( jpt->at(j),weight  );
						h_CMJet[NoFile][ybin]  ->Fill( chMult->at(j),weight );
						h_NMJet[NoFile][ybin]  ->Fill( neMult->at(j),weight );
					}
				}// valid ybin if 
			} // end of jet loop
		} // end of event loop
	} // end of file loop



//======================================= create efficiency & fake rate assymetric error graphs ===============================


	TCanvas *pad_chf = new TCanvas("pad_chf", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_chf->Divide(PadColumns,PadRows);
	TCanvas *pad_nhf = new TCanvas("pad_nhf", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_nhf->Divide(PadColumns,PadRows);
	TCanvas *pad_cemf = new TCanvas("pad_cemf", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_cemf->Divide(PadColumns,PadRows);
	TCanvas *pad_nemf = new TCanvas("pad_nemf", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_nemf->Divide(PadColumns,PadRows);
	TCanvas *pad_muf = new TCanvas("pad_muf", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_muf->Divide(PadColumns,PadRows);
	TCanvas *pad_phi = new TCanvas("pad_phi", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_phi->Divide(PadColumns,PadRows);
	TCanvas *pad_CM = new TCanvas("pad_CM", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_CM->Divide(PadColumns,PadRows);
	TCanvas *pad_NM = new TCanvas("pad_NM", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_NM->Divide(PadColumns,PadRows);
	TCanvas *pad_pt = new TCanvas("pad_pt", "",Canvas_Xpixels,Canvas_Ypixels);
	pad_pt->Divide(PadColumns,PadRows);


	TCanvas *c_eta = new TCanvas("c_eta", "",1);
	TCanvas *c_SumEt = new TCanvas("c_SumEt", "",1);


	TPaveText *paveCMS = new TPaveText(0.45,0.95,0.5,1.0,"NDC");
	// paveCMS->AddText("CMS Preliminary L=9.2 fb^{-1} #sqrt{s} = 13 TeV");
//	paveCMS->AddText("CMS Preliminary#sqrt{s} = 13 TeV");
	paveCMS->AddText("CMS Simulation #sqrt{s} = 14 TeV");
	paveCMS->SetFillColor(0);
	paveCMS->SetBorderSize(0);
	paveCMS->SetTextSize(0.04);

	TLegend *leg1 =new TLegend(.1, .6, .9, .9);//7899//4899
	leg1->SetTextSize(0.055);
	leg1->SetFillColor(0); 
	leg1->SetBorderSize(0);


	TLegend *leg3[eta_bins];

//	TLegend *leg3 =new TLegend(.2, .75, .8, .85);//7899//4899
//	leg3->SetTextSize(0.04);
//	leg3->SetFillColor(0); 
//	leg3->SetBorderSize(0);
//	char res_text[NoFiles][200];

	char res_text[200];
	TLegend *leg2 =new TLegend(.5, .7, .7, .9);//7899//4899
	leg2->SetTextSize(0.03);
	leg2->SetFillColor(0); 
	leg2->SetBorderSize(0);


	 //dummy histograms to be used as frames in the plots
	TH1D *frameCHF[eta_bins], *frameNHF[eta_bins], *frameCEMF[eta_bins], *frameNEMF[eta_bins], *frameMUF[eta_bins], *frameCM[eta_bins], *frameNM[eta_bins], *framePhi[eta_bins], *framePt[eta_bins], *frameEta, *frameSumEt, *frameMinDR;


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		leg1->AddEntry(h_CHFJet[NoFile][0], legend_array[NoFile], "L");
		leg2->AddEntry(h_CHFJet[NoFile][0], legend_array[NoFile], "L");

		if (NoFile==0 )  frameEta = InitiateFrameOnCanvasPad(c_eta, 0 , "frameEta", "Jet eta", "Entries", -5., 5., 0.1, 100000., true, paveCMS);
		if (scale_histos && h_ETAJet[NoFile]->Integral()>0 ) h_ETAJet[NoFile]->Scale(h_ETAJet[0]->Integral() / h_ETAJet[NoFile]->Integral());
		DrawHistoToCanvasPad(c_eta, 0, h_ETAJet[NoFile], Colors[NoFile], 1);
		if( NoFile== 0)  leg2->Draw("same"); 


		if (NoFile==0 )  frameSumEt = InitiateFrameOnCanvasPad(c_SumEt, 0 , "frameSumEt", "SumEt", "Entries", 0., 7000., 0.1, 100000., true, paveCMS);
		if (scale_histos && h_SumEt[NoFile]->Integral()>0 ) h_SumEt[NoFile]->Scale(h_SumEt[0]->Integral() / h_SumEt[NoFile]->Integral());
		DrawHistoToCanvasPad(c_SumEt, 0, h_SumEt[NoFile], Colors[NoFile], 1);
		if( NoFile== 0)  leg2->Draw("same"); 

		
		for(int iy=0; iy<eta_bins; iy++)
		{
			
			double etamin = yBnd[iy];
			double etamax = yBnd[iy+1];
			const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
			TLatex *teta = new TLatex(0.38,0.86,seta); //cout<<seta<<endl;
			teta->SetNDC();
			teta->SetTextSize(0.06);


			if (NoFile==0 )  frameCHF[iy] = InitiateFrameOnCanvasPad(pad_chf, iy+1, "frameCHF", "Charged Hadron Fraction", "Entries", 0., 1.1, 0.001, 10., true, paveCMS);
			if (scale_histos && h_CHFJet[NoFile][iy]->Integral()>0 ) h_CHFJet[NoFile][iy]->Scale(h_CHFJet[0][iy]->Integral()/h_CHFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_chf, iy+1, h_CHFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad_chf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_chf->cd(iy+1); teta->Draw(); }

			
			if (NoFile==0 )  frameNHF[iy] = InitiateFrameOnCanvasPad(pad_nhf,iy+1, "frameNHF", "Neutral Hadron Fraction", "Entries", 0., 1.1, 0.001, 10., true, paveCMS);
			if (scale_histos && h_NHFJet[NoFile][iy]->Integral()>0 ) h_NHFJet[NoFile][iy]->Scale(h_NHFJet[0][iy]->Integral()/h_NHFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_nhf, iy+1, h_NHFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_nhf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_nhf->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameCEMF[iy] = InitiateFrameOnCanvasPad(pad_cemf,iy+1, "frameCEMF", "Charged E/M Fraction", "Entries", 0., 1.1, 0.001, 10., true, paveCMS);
			if (scale_histos && h_CEMFJet[NoFile][iy]->Integral()>0 ) h_CEMFJet[NoFile][iy]->Scale(h_CEMFJet[0][iy]->Integral()/h_CEMFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_cemf, iy+1, h_CEMFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_cemf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_cemf->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNEMF[iy] = InitiateFrameOnCanvasPad(pad_nemf,iy+1, "frameNEMF", "Neutral E/M Fraction", "Entries", 0., 1.1, 0.001, 10., true, paveCMS);
			if (scale_histos && h_NEMFJet[NoFile][iy]->Integral()>0 ) h_NEMFJet[NoFile][iy]->Scale(h_NEMFJet[0][iy]->Integral()/h_NEMFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_nemf, iy+1, h_NEMFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_nemf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_nemf->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameMUF[iy] = InitiateFrameOnCanvasPad(pad_muf,iy+1, "frameMUF", "Muon Fraction", "Entries", 0., 1.1, 0.001, 10., true, paveCMS);
			if (scale_histos && h_MUFJet[NoFile][iy]->Integral()>0 ) h_MUFJet[NoFile][iy]->Scale(h_MUFJet[0][iy]->Integral()/h_MUFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_muf, iy+1, h_MUFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_muf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_muf->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  framePhi[iy] = InitiateFrameOnCanvasPad(pad_phi,iy+1, "framePhi", "Jet Phi", "Entries", -3.14, 3.14, 0.001, 10., true, paveCMS);
			if (scale_histos && h_PHIJet[NoFile][iy]->Integral()>0 ) h_PHIJet[NoFile][iy]->Scale(h_PHIJet[0][iy]->Integral()/h_PHIJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_phi, iy+1, h_PHIJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_phi->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_phi->cd(iy+1); teta->Draw(); }


			if (NoFile==0 ) frameCM[iy] = InitiateFrameOnCanvasPad(pad_CM,iy+1, "frameCM", "Charged Multiplicity", "Entries", 0., 80., 0.001, 10., true, paveCMS);
			if (scale_histos && h_CMJet[NoFile][iy]->Integral()>0 ) h_CMJet[NoFile][iy]->Scale(h_CMJet[0][iy]->Integral()/h_CMJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_CM, iy+1, h_CMJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_CM->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_CM->cd(iy+1); teta->Draw(); }


			if (NoFile==0 ) frameNM[iy] = InitiateFrameOnCanvasPad(pad_NM,iy+1, "frameNM", "Neutral Multiplicity", "Entries", 0., 80., 0.001, 10., true, paveCMS);
			if (scale_histos && h_NMJet[NoFile][iy]->Integral()>0 ) h_NMJet[NoFile][iy]->Scale(h_NMJet[0][iy]->Integral()/h_NMJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_NM, iy+1, h_NMJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_NM->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_NM->cd(iy+1); teta->Draw(); }


			if (NoFile==0 ) framePt[iy] = InitiateFrameOnCanvasPad(pad_pt,iy+1, "framepT", "Jet pT (GeV)", "Entries", 0., pThistoMax, 0.001, 10., true, paveCMS);
			if (scale_histos && h_ptJet[NoFile][iy]->Integral()>0 ) h_ptJet[NoFile][iy]->Scale(h_ptJet[0][iy]->Integral()/h_ptJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_pt, iy+1, h_ptJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_pt->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_pt->cd(iy+1); teta->Draw(); }

		}
	 

		if(eta_bins > 1 &&  NoFile == NoFiles-1 )
		{
			pad_chf->cd(eta_bins+1);	leg1->Draw();
			pad_nhf->cd(eta_bins+1);	leg1->Draw();
			pad_cemf->cd(eta_bins+1);	leg1->Draw();
			pad_nemf->cd(eta_bins+1);	leg1->Draw();
			pad_muf->cd(eta_bins+1);	leg1->Draw();
			pad_phi->cd(eta_bins+1);	leg1->Draw();
			pad_CM->cd(eta_bins+1);		leg1->Draw();
			pad_NM->cd(eta_bins+1);		leg1->Draw();
			pad_pt->cd(eta_bins+1);		leg1->Draw();
		}

	} // end of loop on files

	if (Save_Plots)
	{ 
		sprintf(filename,"%s/%s/%s_chf.png",analyzer_path,output_directory,image_name);
		pad_chf->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_nhf.png",analyzer_path,output_directory,image_name);
		pad_nhf->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_cemf.png",analyzer_path,output_directory,image_name);
		pad_cemf->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_nemf.png",analyzer_path,output_directory,image_name);
		pad_nemf->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_muf.png",analyzer_path,output_directory,image_name);
		pad_muf->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_phi.png",analyzer_path,output_directory,image_name);
		pad_phi->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_CM.png",analyzer_path,output_directory,image_name);
		pad_CM->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NM.png",analyzer_path,output_directory,image_name);
		pad_NM->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_pt.png",analyzer_path,output_directory,image_name);
		pad_pt->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_eta.png",analyzer_path,output_directory,image_name);
		c_eta->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_SumEt.png",analyzer_path,output_directory,image_name);
		c_SumEt->SaveAs(filename);
   }
}

