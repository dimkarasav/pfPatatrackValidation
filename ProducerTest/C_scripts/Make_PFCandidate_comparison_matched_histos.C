//==========================================
// Initial Author: Dimitrios Karasavvas
// Date:   29 Jan 2021
//==========================================




#include "include/include_functions.h"
#include "include/Input_definition.h"



void Make_PFCandidate_comparison_matched_histos()
{

//	gROOT->LoadMacro("setTDRStyle_teliko.C");
//	setTDRStyle_teliko(); 
//	gStyle->SetOptStat(0);

	


	double pTmax = 1000;

	char dir[700];
	sprintf(dir,"%s/%s", analyzer_path,output_directory );
	createDirectory(dir);

	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);
	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;


	TFile *f[NoFiles];
	TTree *tree[NoFiles];

	char eta_bins_legend[eta_bins][25];
	char name[256]; 

	int treeEntries[NoFiles];


	for (int iy=0; iy< eta_bins; iy++)
	{
		if (iy==0) sprintf( eta_bins_legend[iy], "|#eta|<%3.1f" , yBnd[iy+1] );  
		else
		{
			sprintf( eta_bins_legend[iy], "%3.1f<|#eta|<%3.1f" , yBnd[iy] ,yBnd[iy+1] );
		}
	}

	TH1D *h_ChargedHadron_pT[NoFiles][eta_bins], *h_NeutralHadron_pT[NoFiles][eta_bins], *h_ChargedHadron_SumpT[NoFiles][eta_bins], *h_NeutralHadron_SumpT[NoFiles][eta_bins], *h_NumberOfNeutralHadrons[NoFiles][eta_bins], *h_NumberOfChargedHadrons[NoFiles][eta_bins], *h_ChargedHadron_DR[NoFiles][eta_bins], *h_NeutralHadron_DR[NoFiles][eta_bins], *h_girth[NoFiles][eta_bins], *h_quadratic_moment[NoFiles][eta_bins], *h_ChargedHadron_AveragepT[NoFiles][eta_bins], *h_NeutralHadron_AveragepT[NoFiles][eta_bins];

	TH2D *h_ChargedHadron_eta_phi[NoFiles][eta_bins], *h_NeutralHadron_eta_phi[NoFiles][eta_bins];



	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int iy=0; iy< eta_bins; iy++)
		{
	
			sprintf(name,"h_ChargedHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_pT[NoFile][iy] = new TH1D(name, "", 64,0,pTmax);
			if (useWeights) 	h_ChargedHadron_pT[NoFile][iy]->Sumw2();
				
			sprintf(name,"h_NeutralHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_pT[NoFile][iy] = new TH1D(name, "",  64,0,pTmax);
			if (useWeights) 	h_NeutralHadron_pT[NoFile][iy]->Sumw2();
	
			sprintf(name,"h_ChargedHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_SumpT[NoFile][iy] = new TH1D(name, "", 64,0,pTmax);
			if (useWeights) 	h_ChargedHadron_SumpT[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_SumpT[NoFile][iy] = new TH1D(name, "", 64,0,pTmax);
			if (useWeights) 	h_NeutralHadron_SumpT[NoFile][iy]->Sumw2();

			sprintf(name,"h_ChargedHadron_AveragepT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_AveragepT[NoFile][iy] = new TH1D(name, "", 500,0,pTmax);
			if (useWeights) 	h_ChargedHadron_AveragepT[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_AveragepT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_AveragepT[NoFile][iy] = new TH1D(name, "", 500,0,pTmax);
			if (useWeights) 	h_NeutralHadron_AveragepT[NoFile][iy]->Sumw2();

			sprintf(name,"h_NumberOfChargedHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfChargedHadrons[NoFile][iy] = new TH1D(name, "", 50,0,50);
			if (useWeights) 	h_NumberOfChargedHadrons[NoFile][iy]->Sumw2();

			sprintf(name,"h_NumberOfNeutralHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfNeutralHadrons[NoFile][iy] = new TH1D(name, "",  50,0,50);
			if (useWeights) 	h_NumberOfNeutralHadrons[NoFile][iy]->Sumw2();

			sprintf(name,"h_ChargedHadron_DR_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_DR[NoFile][iy] = new TH1D(name, "",  50,0, 0.9);
			if (useWeights) 	h_ChargedHadron_DR[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_DR_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_DR[NoFile][iy] = new TH1D(name, "",  50,0, 0.9);
			if (useWeights) 	h_NeutralHadron_DR[NoFile][iy]->Sumw2();

			sprintf(name,"h_ChargedHadron_eta_phi_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_eta_phi[NoFile][iy] = new TH2D(name, "",  50,-1.0, 1.0, 50, -1.0, 1.0);
			if (useWeights) 	h_ChargedHadron_eta_phi[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_eta_phi_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_eta_phi[NoFile][iy] = new TH2D(name, "",  50,-1.0, 1.0, 50, -1.0, 1.0);
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
		tree[NoFile]->SetBranchAddress("gen_neMult",&gen_neMult);
		tree[NoFile]->SetBranchAddress("gen_SumEt",&gen_SumEt);

		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);


		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		int LoopEntries;
		if (EventsToProcess< 0 || EventsToProcess > treeEntries[NoFile]) LoopEntries = treeEntries[NoFile];
		else LoopEntries = EventsToProcess;

		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
		if ( LoopEntries != treeEntries[NoFile] )  cout<<"\nNumber of entries to be processed: "<< LoopEntries << endl;

		int reco_size, gen_size;
		int LoopSize20percent = 0.2*LoopEntries; 
		for (int i=0; i<LoopEntries; i++) //event loop
		{

//			if (i<5 || i%500000==0 || i==(LoopEntries-1)) cout << " Reached entry: "<< i<< "   Loop complete: "<< 100*i/(LoopEntries-1) << "\%" <<"\n";
			
			if (i<5 || i%(LoopSize20percent)==0 || i==(LoopEntries-1)) cout << " Reached entry: "<< i<< "   Loop complete: "<< 100*i/(LoopEntries-1) << "\%" <<"\n";
			
			tree[NoFile]->GetEntry(i);

			reco_size = jpt->size();
			gen_size  = gen_jpt->size();

			int *reco_jets_matched_sequence;	
			reco_jets_matched_sequence = GetRecoToGenMatchSequence(gen_eta, gen_phi, eta, phi, DR_threshold);

//			cout << "\n Processing event " << i << " number of reco_jets : "<< reco_size <<endl; 
			for ( int j=0; j<gen_size; j++) // loop over jets
			{

				if (reco_jets_matched_sequence[j]>-0.1 && ( jpt->at(reco_jets_matched_sequence[j]) < pTlowCut || jpt->at(reco_jets_matched_sequence[j]) > pThighCut ) ) continue;
				if (reco_jets_matched_sequence[j]<0 ) continue ; // if no match was made, skip this gen jet
	//			int cand_size = cand_vector->size();
				int jet_ybin = getBin(fabs(eta->at(reco_jets_matched_sequence[j])),yBnd, eta_bins); //get the jet eta bin

				int cand_size = ( jetconstituents->at( reco_jets_matched_sequence[j] ) ).size();
//				cout << "jet "<< reco_jets_matched_sequence[j] << " cand size = " << cand_size << endl;

				float chargedCand_SumpT = 0;
				float neutralCand_SumpT = 0;
				int neutralPFcand_number = 0;
				int chargedPFcand_number = 0;
				double jet_girth = 0;
				double jet_quadr_moment = 0;

				for (int l=0; l< cand_size; l++)
				{
					int jet_candidate = (jetconstituents->at(reco_jets_matched_sequence[j])).at(l) ; //get the jet candidate number
					//cout << "\njet candidate number: " << jet_candidate << endl;

					int cand_pdgID = pfcandpdgid[jet_candidate];
					float cand_eta = pfcandeta[jet_candidate];
					float cand_phi = pfcandphi[jet_candidate];
					float cand_pT = pfcandpt[jet_candidate];
					double 	cand_DR = sqrt( pow( eta->at(reco_jets_matched_sequence[j])-cand_eta ,2) + pow( phi->at(reco_jets_matched_sequence[j])-cand_phi ,2) );

					jet_girth = jet_girth + (cand_pT/jpt->at(reco_jets_matched_sequence[j]))*cand_DR;
					jet_quadr_moment = jet_quadr_moment + (cand_pT/jpt->at(reco_jets_matched_sequence[j]))*cand_DR*cand_DR;

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
								h_ChargedHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(reco_jets_matched_sequence[j]), cand_phi-phi->at(reco_jets_matched_sequence[j])  );
							}
							else
							{
								h_ChargedHadron_pT[NoFile][jet_ybin]->Fill(cand_pT, weight); // fill charged hadron pT's
								h_ChargedHadron_DR[NoFile][jet_ybin]->Fill( cand_DR, weight );
								h_ChargedHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(reco_jets_matched_sequence[j]), cand_phi-phi->at(reco_jets_matched_sequence[j]), weight  );
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
								h_NeutralHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(reco_jets_matched_sequence[j]), cand_phi-phi->at(reco_jets_matched_sequence[j])  );
							}
							else
							{
								h_NeutralHadron_pT[NoFile][jet_ybin]->Fill(cand_pT, weight); // fill neutral hadron pT's with weight
								h_NeutralHadron_DR[NoFile][jet_ybin]->Fill( cand_DR, weight );
								h_NeutralHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(reco_jets_matched_sequence[j]), cand_phi-phi->at(reco_jets_matched_sequence[j]), weight  );
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
					if (chargedPFcand_number>0) h_ChargedHadron_AveragepT[NoFile][jet_ybin]->Fill(chargedCand_SumpT/chargedPFcand_number);
					if (neutralPFcand_number>0) h_NeutralHadron_AveragepT[NoFile][jet_ybin]->Fill(neutralCand_SumpT/neutralPFcand_number);
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
					if (chargedPFcand_number>0) h_ChargedHadron_AveragepT[NoFile][jet_ybin]->Fill(chargedCand_SumpT/chargedPFcand_number, weight);
					if (neutralPFcand_number>0) h_NeutralHadron_AveragepT[NoFile][jet_ybin]->Fill(neutralCand_SumpT/neutralPFcand_number, weight);
					h_girth[NoFile][jet_ybin]->Fill(jet_girth, weight);
					h_quadratic_moment[NoFile][jet_ybin]->Fill(jet_quadr_moment, weight);
				}

			} //  end of jet loop

			gen_size = 0.;
			reco_size = 0;
			delete reco_jets_matched_sequence;

		} // end of event loop


		delete tree[NoFile];
		delete f[NoFile];
	} // end of file loop


	char filename[1500];
	sprintf(filename, "%s/%s/PF_candidates_matched_histos_%s.root",analyzer_path,output_directory,image_name );

	TFile *fout = new TFile(filename, "recreate");
	fout->cd();



	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for(int eta_bin=0; eta_bin<eta_bins; eta_bin++)
		{
			h_ChargedHadron_pT[NoFile][eta_bin]       ->Write();
			h_NeutralHadron_pT[NoFile][eta_bin]       ->Write();
			h_ChargedHadron_SumpT[NoFile][eta_bin]    ->Write();
			h_NeutralHadron_SumpT[NoFile][eta_bin]    ->Write();
			h_ChargedHadron_AveragepT[NoFile][eta_bin]->Write();
			h_NeutralHadron_AveragepT[NoFile][eta_bin]->Write();
			h_NumberOfNeutralHadrons[NoFile][eta_bin] ->Write();
			h_NumberOfChargedHadrons[NoFile][eta_bin] ->Write();
			h_ChargedHadron_DR[NoFile][eta_bin]       ->Write();
			h_NeutralHadron_DR[NoFile][eta_bin]       ->Write();
			h_girth[NoFile][eta_bin]                  ->Write();
			h_quadratic_moment[NoFile][eta_bin]       ->Write();
			h_ChargedHadron_eta_phi[NoFile][eta_bin]  ->Write();
			h_NeutralHadron_eta_phi[NoFile][eta_bin]  ->Write();
		}	
	}

}

