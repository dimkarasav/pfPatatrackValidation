#include "include/include_functions.h"
#include "include/Input_definition.h"


void Make_Tracking_Comparisons_histos_btb()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);


	double pTlowCut_leading = 60;
	double pTlowCut_subleading = 30;

	double deltaPhi = 2.7;


	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);
	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;

	char eta_bins_legend[eta_bins][25];


	for (int iy=0; iy< eta_bins; iy++)
	{
		if (iy==0) sprintf( eta_bins_legend[iy], "|#eta|<%3.1f" , yBnd[iy+1] );  
		else
		{
			sprintf( eta_bins_legend[iy], "%3.1f<|#eta|<%3.1f" , yBnd[iy] ,yBnd[iy+1] );
		}
	}


	TH1D *h_CHFJet[NoFiles][eta_bins], *h_NHFJet[NoFiles][eta_bins], *h_CEMFJet[NoFiles][eta_bins], *h_NEMFJet[NoFiles][eta_bins], *h_MUFJet[NoFiles][eta_bins], *h_CMJet[NoFiles][eta_bins], *h_NMJet[NoFiles][eta_bins], *h_ptJet[NoFiles][eta_bins], *h_PHIJet[NoFiles][eta_bins];

	TH1D *h_ETAJet[NoFiles], *h_METovSumEt[NoFiles], *h_minDR[NoFiles], *h_MET[NoFiles] ;


	char name[256]; 
	

	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		sprintf(name,"h_ETAJet_%s",legend_array[NoFile]);
		h_ETAJet[NoFile] = new TH1D(name, "", 60,-5,5);

	
		sprintf(name,"h_METovSumEt_%s",legend_array[NoFile]);
		h_METovSumEt[NoFile] = new TH1D(name, "", 60, 0, 1.2); // 40,0,1.0

		sprintf(name,"h_MET_%s",legend_array[NoFile]);
		h_MET[NoFile] = new TH1D(name, "", 120, 0, pThighCut); // 40,0,1.0

		if (useWeights)
		{
			h_ETAJet[NoFile]->Sumw2();
			h_MET[NoFile]->Sumw2();
			h_METovSumEt[NoFile]->Sumw2();
		}
	

		for(Int_t h=0; h<eta_bins;h++)
		{ 
		//========== reco jets===========
			

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
			if (useWeights) h_CMJet[NoFile][h]->Sumw2();
		
			sprintf(name,"h_PHIJet_%s_bin%i",legend_array[NoFile],h);
			h_PHIJet[NoFile][h] = new TH1D(name, "", 40,-4,4);
			if (useWeights) h_PHIJet[NoFile][h]->Sumw2();

			sprintf(name,"h_NMJet_%s_bin%i",legend_array[NoFile],h);
			h_NMJet[NoFile][h] = new TH1D(name, "", 80, 0.-0.5, 80.-0.5);
			if (useWeights) h_NMJet[NoFile][h]->Sumw2();

			sprintf(name,"h_ptJet_%s_bin%i",legend_array[NoFile],h);
			h_ptJet[NoFile][h] = new TH1D(name, "", 50,0,pThighCut);
			if (useWeights) h_ptJet[NoFile][h]->Sumw2();

			
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
 	double SumEt,gen_SumEt,met, gen_met;
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

		tree[NoFile]->SetBranchAddress("met",&met);
		tree[NoFile]->SetBranchAddress("SumEt",&SumEt);
		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);
		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);



		int size, reco_size;

		//======================= f1 jets ===================================

		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		int LoopEntries;
		if (EventsToProcess< 0 || EventsToProcess > treeEntries[NoFile]) LoopEntries = treeEntries[NoFile];
		else LoopEntries = EventsToProcess;
		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
		int LoopSize20percent = 0.2*LoopEntries; 

		for (int i=0; i<LoopEntries; i++) //event loop

		{
			if (i<5 || i%(LoopSize20percent)==0 || i==(LoopEntries-1)) cout << " Reached entry: "<< i<< "   Loop complete: "<< 100*(i+1)/(LoopEntries) << "\%" <<"\n";
			tree[NoFile]->GetEntry(i);
			reco_size = jpt->size();

			if ( reco_size < 2) continue;	 									// skip event with less than 2 jets
			if ( fabs( phi->at(0) - phi->at(1) ) < deltaPhi ) continue; 	// skip event if the two leading jets are not back-to-back
			if (jpt->at(0) < pTlowCut_leading ) continue; //leading jet pT default cut
			if (jpt->at(1) < pTlowCut_subleading ) continue; //subleading jet pT default cut
			if ( !(jpt->at(0)> pTlowCut && jpt->at(0) < pThighCut) ) continue;// skip event if leading jet does not pass pT cuts
			if ( !(jpt->at(1)> pTlowCut && jpt->at(1) < pThighCut) ) continue;// skip event if sub-leading jet does not pass pT cuts
			//if ( jpt->at(0)/jpt->at(1)>1.3 ) continue; // skip events with 2 leading jets more than 30% unbalanced in terms of pt 
			
			for ( int j=0; j<2; j++) 
			{

				if (j==0 && !useWeights && SumEt>0) { h_METovSumEt[NoFile]->Fill(met/SumEt); h_MET[NoFile]->Fill(met); }
				if (j==0 &&  useWeights && SumEt>0) { h_METovSumEt[NoFile]->Fill(met/SumEt, weight); h_MET[NoFile]->Fill(met, weight);}

				if (!useWeights)	h_ETAJet[NoFile]->Fill( eta->at(j) );
				else 				h_ETAJet[NoFile]->Fill( eta->at(j), weight );
				int ybin = getBin(fabs(eta->at(j)), yBnd, eta_bins);
				if (ybin > -1)//fill hist's in the corresponding eta bin
				{
					if (!useWeights)
					{
						h_CHFJet[NoFile][ybin] ->Fill( chf->at(j)    ); 
						h_NHFJet[NoFile][ybin] ->Fill( nhf->at(j)    );
						h_CEMFJet[NoFile][ybin]->Fill( cemf->at(j)   );
						h_NEMFJet[NoFile][ybin]->Fill( nemf->at(j)   );
						h_MUFJet[NoFile][ybin] ->Fill( muf->at(j)    );
						h_PHIJet[NoFile][ybin] ->Fill( phi->at(j)    );
						h_ptJet[NoFile][ybin]  ->Fill( jpt->at(j)    );
						h_CMJet[NoFile][ybin]  ->Fill( chMult->at(j) );
						h_NMJet[NoFile][ybin]  ->Fill( neMult->at(j) );
					}
					else 
					{
						h_CHFJet[NoFile][ybin] ->Fill( chf->at(j)   , weight  ); 
						h_NHFJet[NoFile][ybin] ->Fill( nhf->at(j)   , weight  );
						h_CEMFJet[NoFile][ybin]->Fill( cemf->at(j)  , weight  );
						h_NEMFJet[NoFile][ybin]->Fill( nemf->at(j)  , weight  );
						h_MUFJet[NoFile][ybin] ->Fill( muf->at(j)   , weight  );
						h_PHIJet[NoFile][ybin] ->Fill( phi->at(j)   , weight  );
						h_ptJet[NoFile][ybin]  ->Fill( jpt->at(j)   , weight  );
						h_CMJet[NoFile][ybin]  ->Fill( chMult->at(j), weight  );
						h_NMJet[NoFile][ybin]  ->Fill( neMult->at(j), weight  );
					}
				}// valid ybin if 
			} // end of jet loop
		} // end of event loop
	} // end of file loop

	char filename[1500];
	sprintf(filename, "%s/%s/Jet_characteristics_histos_%s.root",analyzer_path,output_directory,image_name );

	TFile *fout = new TFile(filename, "recreate");
	fout->cd();

	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		h_ETAJet[NoFile]->Write();
		h_METovSumEt[NoFile]->Write();
		h_MET[NoFile]->Write();

		for(int ybin=0; ybin<eta_bins; ybin++)
		{
			h_CHFJet[NoFile][ybin]           ->Write(); 
			h_NHFJet[NoFile][ybin]           ->Write();
			h_CEMFJet[NoFile][ybin]          ->Write();
			h_NEMFJet[NoFile][ybin]          ->Write();
			h_MUFJet[NoFile][ybin]           ->Write();
			h_PHIJet[NoFile][ybin]           ->Write();
			h_ptJet[NoFile][ybin]            ->Write();
			h_CMJet[NoFile][ybin]  		     ->Write();
			h_NMJet[NoFile][ybin]            ->Write();
		}
	}


}

