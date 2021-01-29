//==========================================
// Initial Author: Dimitrios Karasavvas
// Date:   29 Jan 2021
//==========================================

#include "include/include_functions.h"
#include "include/Input_definition.h"


void Make_Tracking_Comparisons_matched_histos()
{

//	gROOT->LoadMacro("setTDRStyle_teliko.C");
//	setTDRStyle_teliko(); 
//	gStyle->SetOptStat(0);

	double pThistoMax = 1000;
	if(pThistoMax< pThighCut) pThistoMax = pThighCut;

	cout << pThistoMax <<endl;

	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);
	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;

	char dir[700];
	sprintf(dir,"%s/%s", analyzer_path,output_directory );
	createDirectory(dir);

	bool plot_minDR = true;
	char eta_bins_legend[eta_bins][25];

	for (int iy=0; iy< eta_bins; iy++)
	{
		if (iy==0) sprintf( eta_bins_legend[iy], "|#eta|<%3.1f" , yBnd[iy+1] );  
		else
		{
			sprintf( eta_bins_legend[iy], "%3.1f<|#eta|<%3.1f" , yBnd[iy] ,yBnd[iy+1] );
		}
	}


	TH1D *h_CHFJet[NoFiles][eta_bins], *h_NHFJet[NoFiles][eta_bins], *h_CEMFJet[NoFiles][eta_bins], *h_NEMFJet[NoFiles][eta_bins], *h_MUFJet[NoFiles][eta_bins], *h_CMJet[NoFiles][eta_bins], *h_NMJet[NoFiles][eta_bins], *h_ptJet[NoFiles][eta_bins], *h_PHIJet[NoFiles][eta_bins], *h_gen_pT[NoFiles][eta_bins],*h_gen_pT_all[NoFiles][eta_bins], *h_reco_pT_unmatched[NoFiles][eta_bins], *h_reco_pT_all[NoFiles][eta_bins];

	

	TH1D *h_ETAJet[NoFiles], *h_METovSumEt[NoFiles], *h_minDR[NoFiles], *h_MET[NoFiles] ;

	TGraphAsymmErrors *Reco_eff[NoFiles][eta_bins], *Fake_rate[NoFiles][eta_bins]; 




	char name[1024]; 

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko();


for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		sprintf(name,"h_ETAJet_%s",legend_array[NoFile]);
		h_ETAJet[NoFile] = new TH1D(name, "", 60,-5,5);

	
		sprintf(name,"h_METovSumEt_%s",legend_array[NoFile]);
		h_METovSumEt[NoFile] = new TH1D(name, "", 60, 0, 1.2); // 40,0,1.0

		sprintf(name,"h_MET_%s",legend_array[NoFile]);
		h_MET[NoFile] = new TH1D(name, "", 120, 0, pThistoMax); // 40,0,1.0
	
		sprintf(name,"h_minDR_%s",legend_array[NoFile]);
		h_minDR[NoFile] = new TH1D(name, "", 100, 0.0, 5.0); // 40,0,1.0

		if (useWeights)
		{
			h_ETAJet[NoFile]->Sumw2();
			h_METovSumEt[NoFile]->Sumw2();
			h_MET[NoFile]->Sumw2();
		}
	

		for(Int_t h=0; h<eta_bins;h++)
		{ 
		//========== reco jets===========
			

			sprintf(name,"h_CHFJet_%s_bin%i",legend_array[NoFile],h);
			h_CHFJet[NoFile][h] = new TH1D(name, "", 30, 0, 1.2); //40, 0,1.2
			if (useWeights) h_CHFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_NHFJet_%s_bin%i",legend_array[NoFile],h);
			h_NHFJet[NoFile][h] = new TH1D(name, "", 30, 0, 1.2);
			if (useWeights) h_NHFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_CEMFJet_%s_bin%i",legend_array[NoFile],h);
			h_CEMFJet[NoFile][h] = new TH1D(name, "", 30, 0, 1.2);
			if (useWeights) h_CEMFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_NEMFJet_%s_bin%i",legend_array[NoFile],h);
			h_NEMFJet[NoFile][h] = new TH1D(name, "", 30, 0, 1.2);
			if (useWeights) h_NEMFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_MUFJet_%s_bin%i",legend_array[NoFile],h);
			h_MUFJet[NoFile][h] = new TH1D(name, "", 30, 0, 1.2);
			if (useWeights) h_MUFJet[NoFile][h]->Sumw2();

			sprintf(name,"h_CMJet_%s_bin%i",legend_array[NoFile],h);
			h_CMJet[NoFile][h] = new TH1D(name, "", 120, 0.0-0.5, 120.0 - 0.5);
			if (useWeights) h_CMJet[NoFile][h]->Sumw2();
		
			sprintf(name,"h_PHIJet_%s_bin%i",legend_array[NoFile],h);
			h_PHIJet[NoFile][h] = new TH1D(name, "", 40,-4,4);
			if (useWeights) h_PHIJet[NoFile][h]->Sumw2();

			sprintf(name,"h_NMJet_%s_bin%i",legend_array[NoFile],h);
			h_NMJet[NoFile][h] = new TH1D(name, "", 120, 0.-0.5, 120.-0.5);
			if (useWeights) h_NMJet[NoFile][h]->Sumw2();

			sprintf(name,"h_ptJet_%s_bin%i",legend_array[NoFile],h);
			h_ptJet[NoFile][h] = new TH1D(name, "", 50,0,pThistoMax);
			if (useWeights) h_ptJet[NoFile][h]->Sumw2();

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

	std::vector<float> *gen_jpt = 0;
	std::vector<float> *gen_eta = 0;
	std::vector<float> *gen_phi = 0;
	std::vector<float> *gen_nhf = 0;
	std::vector<float> *gen_chf = 0;
	std::vector<float> *gen_cemf = 0;
	std::vector<float> *gen_nemf = 0;
	std::vector<float> *gen_npr = 0;
	std::vector<float> *gen_muf = 0;
	std::vector<float> *gen_chMult = 0;
	std::vector<float> *gen_neMult = 0;


	std::vector<float> *npu = 0;
	std::vector<int> *PileupInteractions = 0;
	std::vector<int> *PileupOriginBX = 0;

	float weight;
 	double SumEt, gen_SumEt, met, gen_met;
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
//		tree[NoFile]->SetBranchAddress("chm",&chMult);
		tree[NoFile]->SetBranchAddress("neMult",&neMult);

		tree[NoFile]->SetBranchAddress("SumEt",&SumEt);
		tree[NoFile]->SetBranchAddress("met",&met);
		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);
		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);


		tree[NoFile]->SetBranchAddress("gen_jpt",&gen_jpt);
		tree[NoFile]->SetBranchAddress("gen_eta",&gen_eta);
		tree[NoFile]->SetBranchAddress("gen_phi",&gen_phi);
		tree[NoFile]->SetBranchAddress("gen_nhf",&gen_nhf);
		tree[NoFile]->SetBranchAddress("gen_chf",&gen_chf);
		tree[NoFile]->SetBranchAddress("gen_cemf",&gen_cemf);
		tree[NoFile]->SetBranchAddress("gen_nemf",&gen_nemf);
		tree[NoFile]->SetBranchAddress("gen_npr",&gen_npr);
		tree[NoFile]->SetBranchAddress("gen_muf",&gen_muf);
		tree[NoFile]->SetBranchAddress("gen_chMult",&gen_chMult);
//		tree[NoFile]->SetBranchAddress("gen_chm",&gen_chMult);
		tree[NoFile]->SetBranchAddress("gen_neMult",&gen_neMult);
		tree[NoFile]->SetBranchAddress("gen_SumEt",&gen_SumEt);
		tree[NoFile]->SetBranchAddress("gen_met",&gen_met);


		int size, reco_size;

		//======================= f1 jets ===================================

		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		int LoopEntries;
		if (EventsToProcess< 0 || EventsToProcess > treeEntries[NoFile]) LoopEntries = treeEntries[NoFile];
		else LoopEntries = EventsToProcess;
		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
		int LoopSize20percent = 0.2*LoopEntries; 
//		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
		for (int i=0; i<LoopEntries; i++) //event loop
		{

			if (i<5 || i%(LoopSize20percent)==0 || i==(LoopEntries-1)) cout << " Reached entry: "<< i<< "   Loop complete: "<< 100*(i+1)/(LoopEntries) << "\%" <<"\n";
			tree[NoFile]->GetEntry(i);
			size = gen_jpt->size();

			int *reco_jets_matched_sequence;
			int *reco_jets_matched_sequence_forPlot;


//			if (npu->at(0)>80 || npu->at(0)<55) continue;
		
			reco_jets_matched_sequence = GetRecoToGenMatchSequence(gen_eta, gen_phi, eta, phi, DR_threshold);
			if (plot_minDR)	reco_jets_matched_sequence_forPlot = GetRecoToGenMatchSequence(gen_eta, gen_phi, eta, phi, 5.0);


			if ( size > 0 )
			{
				for ( int j=0; j<size; j++) 
				{

				//	cout << " reco matched jet = " << reco_jets_matched_sequence[j] << endl; ;
					if (reco_jets_matched_sequence[j]>-0.1 && ( jpt->at(reco_jets_matched_sequence[j]) < pTlowCut || jpt->at(reco_jets_matched_sequence[j]) > pThighCut ) ) continue;
					if (reco_jets_matched_sequence_forPlot[j]>-0.1 && ( jpt->at(reco_jets_matched_sequence_forPlot[j]) < pTlowCut || jpt->at(reco_jets_matched_sequence_forPlot[j]) > pThighCut ) ) continue;
					if(plot_minDR && reco_jets_matched_sequence_forPlot[j]>=0 )
					{
						double deltaR = sqrt( pow( gen_eta->at(j) - eta->at(reco_jets_matched_sequence_forPlot[j]), 2) + pow( gen_phi->at(j) - phi->at(reco_jets_matched_sequence_forPlot[j]), 2) );
						h_minDR[NoFile]->Fill(deltaR);
					}
					int gen_ybin = getBin(fabs(gen_eta->at(j)),yBnd, eta_bins);
					if (gen_ybin > -1) h_gen_pT_all[NoFile][gen_ybin]->Fill( gen_jpt->at(j) ); // this is filled only with every gen jet.



					if (j==0 && !useWeights && SumEt>0) { h_METovSumEt[NoFile]->Fill(met/SumEt); h_MET[NoFile]->Fill(met);}
					if (j==0 &&  useWeights && SumEt>0) { h_METovSumEt[NoFile]->Fill(met/SumEt, weight); h_MET[NoFile]->Fill(met, weight);}
					if (reco_jets_matched_sequence[j]<0 ) continue ; // if no match was made, skip this gen jet

					if (gen_ybin > -1) h_gen_pT[NoFile][gen_ybin]->Fill( gen_jpt->at(j) ); // this is filled only with matched gen jets.



					if (!useWeights)	h_ETAJet[NoFile]->Fill( eta->at(reco_jets_matched_sequence[j]) );
					else 			h_ETAJet[NoFile]->Fill( eta->at(reco_jets_matched_sequence[j]),weight );
					int ybin = getBin(fabs(eta->at(reco_jets_matched_sequence[j])),yBnd, eta_bins);
					if (ybin > -1)//fill hist's in the corresponding eta bin
					{
						if (!useWeights)
						{
							h_CHFJet[NoFile][ybin] ->Fill( chf->at(reco_jets_matched_sequence[j])  ); 
							h_NHFJet[NoFile][ybin] ->Fill( nhf->at(reco_jets_matched_sequence[j])  );
							h_CEMFJet[NoFile][ybin]->Fill( cemf->at(reco_jets_matched_sequence[j]) );
							h_NEMFJet[NoFile][ybin]->Fill( nemf->at(reco_jets_matched_sequence[j]) );
							h_MUFJet[NoFile][ybin] ->Fill( muf->at(reco_jets_matched_sequence[j])  );
							h_PHIJet[NoFile][ybin] ->Fill( phi->at(reco_jets_matched_sequence[j])  );
							h_ptJet[NoFile][ybin]  ->Fill( jpt->at(reco_jets_matched_sequence[j])  );
							h_CMJet[NoFile][ybin]  ->Fill( chMult->at(reco_jets_matched_sequence[j]) );
							h_NMJet[NoFile][ybin]  ->Fill( neMult->at(reco_jets_matched_sequence[j]) );
						}
						else 
						{
							h_CHFJet[NoFile][ybin] ->Fill( chf->at(reco_jets_matched_sequence[j]),weight ); 
							h_NHFJet[NoFile][ybin] ->Fill( nhf->at(reco_jets_matched_sequence[j]),weight  );
							h_CEMFJet[NoFile][ybin]->Fill( cemf->at(reco_jets_matched_sequence[j]),weight );
							h_NEMFJet[NoFile][ybin]->Fill( nemf->at(reco_jets_matched_sequence[j]),weight );
							h_MUFJet[NoFile][ybin] ->Fill( muf->at(reco_jets_matched_sequence[j]),weight  );
							h_PHIJet[NoFile][ybin] ->Fill( phi->at(reco_jets_matched_sequence[j]),weight  );
							h_ptJet[NoFile][ybin]  ->Fill( jpt->at(reco_jets_matched_sequence[j]),weight  );
							h_CMJet[NoFile][ybin]  ->Fill( chMult->at(reco_jets_matched_sequence[j]),weight );
							h_NMJet[NoFile][ybin]  ->Fill( neMult->at(reco_jets_matched_sequence[j]),weight );
						}
					} 
				} // jet loop
			} // gen size >0 if


			reco_size = jpt->size();
			if (reco_size > 0)
			{
				std::vector<int> unmatched_sequence = GetRecoUnmatchedSequence( reco_size, size, reco_jets_matched_sequence ) ;

				int unmatched_size = unmatched_sequence.size();
				for (int k=0; k<unmatched_size; k++)
				{
					if (  jpt->at( unmatched_sequence.at(k) ) < pTlowCut || jpt->at(  unmatched_sequence.at(k) ) > pThighCut  ) continue;
					int ybin = getBin(fabs( eta->at(  unmatched_sequence.at(k) ) ),yBnd, eta_bins);
					h_reco_pT_unmatched[NoFile][ybin]->Fill( jpt->at(  unmatched_sequence.at(k) ) );

				}

				for (int k=0; k<reco_size; k++)
				{
					if (  jpt->at(k) < pTlowCut || jpt->at(k)  > pThighCut  ) continue;
					int ybin = getBin(fabs( eta->at(k) ),yBnd, eta_bins);
					h_reco_pT_all[NoFile][ybin]->Fill( jpt->at(k) );
				}

			}

			size = 0.;
			reco_size = 0;
			delete reco_jets_matched_sequence;
			delete reco_jets_matched_sequence_forPlot;
		} // end of entries loop
		f[NoFile]->Close();
	} // end of file loop


/*
//======================================= create efficiency & fake rate assymetric error graphs ===============================
	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for(int iy=0; iy<eta_bins; iy++)
		{
			Reco_eff[NoFile][iy] = GetEfficiencyGraph(h_gen_pT[NoFile][iy] , h_gen_pT_all[NoFile][iy] );
			Fake_rate[NoFile][iy] = GetEfficiencyGraph(h_reco_pT_unmatched[NoFile][iy],h_reco_pT_all[NoFile][iy]  ) ; 
		} 
	}
*/
	char filename[1500];
	sprintf(filename, "%s/%s/Jet_characteristics_histos_%s.root",analyzer_path,output_directory,image_name );

	TFile *fout = new TFile(filename, "recreate");
	fout->cd();

	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		h_ETAJet[NoFile]->Write();
		h_METovSumEt[NoFile]->Write();
		h_MET[NoFile]->Write();
		h_minDR[NoFile]->Write();

		for(int ybin=0; ybin<eta_bins; ybin++)
		{
			h_CHFJet[NoFile][ybin]           ->Write(); 
			h_NHFJet[NoFile][ybin]           ->Write();
			h_CEMFJet[NoFile][ybin]          ->Write();
			h_NEMFJet[NoFile][ybin]          ->Write();
			h_MUFJet[NoFile][ybin]           ->Write();
			h_PHIJet[NoFile][ybin]           ->Write();
			h_ptJet[NoFile][ybin]            ->Write();
			h_CMJet[NoFile][ybin]  		 ->Write();
			h_NMJet[NoFile][ybin]            ->Write();
			h_gen_pT_all[NoFile][ybin]       ->Write();
			h_gen_pT[NoFile][ybin]           ->Write();
			h_reco_pT_unmatched[NoFile][ybin]->Write();
			h_reco_pT_all[NoFile][ybin]      ->Write();
		}
	}

	
}

