#include "include/include_functions.h"
#include "include/Input_definition.h"



void Make_pTreco_ov_pTgen_JEScorrected_histos()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);
	gROOT->ForceStyle(kTRUE);
	gStyle->SetOptFit(0);


	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;
	const int pT_bins = sizeof(ptBnd)/sizeof(ptBnd[0])-1;
	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);



	char dir[700];
	sprintf(dir,"%s/%s", analyzer_path,output_directory );
	createDirectory(dir);
//	double DR_threshold = 0.2; //reco_gen matching threshold 
	if (pThighCut > ptBnd[pT_bins-1]) pThighCut = ptBnd[pT_bins];
	if (pTlowCut  < ptBnd[0])         pTlowCut  = ptBnd[0];  

	double pTbinCenter[pT_bins], pTbinCenter_error[pT_bins];

	double pTraw_responses[NoFiles][eta_bins][pT_bins], pTraw_resolution[NoFiles][eta_bins][pT_bins];
	double pTraw_responses_error[NoFiles][eta_bins][pT_bins], pTraw_resolution_error[NoFiles][eta_bins][pT_bins];

	double pTraw_resolution_Gaus[NoFiles][eta_bins][pT_bins], pTraw_response_Gaus[NoFiles][eta_bins][pT_bins];
	double pTraw_resolution_Gaus_error[NoFiles][eta_bins][pT_bins], pTraw_response_Gaus_error[NoFiles][eta_bins][pT_bins];

	TH1D *h_pTreco_ov_pTgen[NoFiles][eta_bins][pT_bins], *h_pTreco_ov_pTgen_JEScorrected[NoFiles][eta_bins][pT_bins], *h_pTreco_ov_pTgen_JEScorrected_noWeights[NoFiles][eta_bins][pT_bins];

//	cout << "eta bins = " << eta_bins;
//	cout << "pT bins = " << pT_bins;
	char name[1500]; 
	sprintf(name, "%s/%s/pTreco_ov_pTgen_%s_histos.root",analyzer_path,output_directory,image_name );
	TFile *f_input = new TFile (name,"READ");


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int  eta_bin=0; eta_bin<eta_bins; eta_bin++)
		{
			for (int  pT_bin=0; pT_bin<pT_bins; pT_bin++)
			{
				sprintf(name,"h_pTreco_ov_pTgen_JEScorrected_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin] = new TH1D(name, "", 600, 0., 6.0); // 40,0,1.0
				if (useWeights) h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Sumw2();

				sprintf(name,"h_pTreco_ov_pTgen_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin] = (TH1D*)(f_input->Get(name));

				sprintf(name,"h_pTreco_ov_pTgen_JEScorrected_noWeights_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pTreco_ov_pTgen_JEScorrected_noWeights[NoFile][eta_bin][pT_bin] = new TH1D(name, "", 600, 0., 6.0); // 40,0,1.0
			}
		}
	}


	for (int  pT_bin=0; pT_bin<pT_bins; pT_bin++)
	{
		pTbinCenter[pT_bin] = 0.5*(ptBnd[pT_bin+1] + ptBnd[pT_bin] );
		pTbinCenter_error[pT_bin] = 0.5* ( ptBnd[pT_bin+1] - ptBnd[pT_bin] );
	}


	char eta_bins_legend[eta_bins][25];
	for (int iy=0; iy< eta_bins; iy++)
	{
		if (iy==0) sprintf( eta_bins_legend[iy], "|#eta|<%3.1f" , yBnd[iy+1] );  
		else 
		{
			sprintf( eta_bins_legend[iy], "%3.1f<|#eta|<%3.1f" , yBnd[iy] ,yBnd[iy+1] );
		}
	}

	char pT_bins_legend[pT_bins][25];
	for (int iy=0; iy< pT_bins; iy++)
	{
		if (iy==0) sprintf( pT_bins_legend[iy], "pT<%f" , ptBnd[iy+1] );  
		else if (iy < pT_bins -1 )
		{
			sprintf( pT_bins_legend[iy], "%f<pT<%f" , ptBnd[iy] ,ptBnd[iy+1] );
		}
		else sprintf( pT_bins_legend[iy], "p_T>%f" , ptBnd[iy] );  
	}


	int treeEntries[NoFiles];
	TFile *f[NoFiles];
	TTree *tree[NoFiles];


	std::vector<float> *jpt = 0;
	std::vector<float> *eta = 0;
	std::vector<float> *phi = 0;


	std::vector<float> *gen_jpt = 0;
	std::vector<float> *gen_eta = 0;
	std::vector<float> *gen_phi = 0;


	std::vector<float> *npu = 0;
	std::vector<int> *PileupInteractions = 0;
	std::vector<int> *PileupOriginBX = 0;

	float weight;


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		f[NoFile] = TFile::Open(input_files[NoFile],"READ"); //read the .root files
		tree[NoFile] = (TTree*)f[NoFile]->Get("dijetscouting/tree"); // get the trees from the files

		tree[NoFile]->SetBranchAddress("jpt",&jpt);
		tree[NoFile]->SetBranchAddress("eta",&eta);
		tree[NoFile]->SetBranchAddress("phi",&phi);
		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);
		tree[NoFile]->SetBranchAddress("gen_jpt",&gen_jpt);
		tree[NoFile]->SetBranchAddress("gen_eta",&gen_eta);
		tree[NoFile]->SetBranchAddress("gen_phi",&gen_phi);



		int gen_size, reco_size;

		treeEntries[NoFile] = tree[NoFile]->GetEntries();


		cout << " Calculating raw responses and resolutions for file" << input_files[NoFile] << " ...  "<< endl;
		TF1 *gausTF1;

		for ( int eta_bin=0; eta_bin< eta_bins; eta_bin++ )
		{

			for (int pT_bin=0; pT_bin< pT_bins; pT_bin++ )
			{
				double Nsize = h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->Integral();
				if ( !(Nsize>0) ) 
				{
					cout << "\nWarning!!! Histogram empty!!  eta bin : " << eta_bins_legend[eta_bin] << "  pT : " << pT_bins_legend[pT_bin] << endl;
					continue;
				}

				double res_mean = h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->GetMean();
				double res_rms = h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->GetRMS(); //
				double res_mean_error = res_rms / Nsize; // error of mean value.
				double res_rms_error = 	CalculateRMSerror( h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin] );


//error of standard deviation found here : https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation


				pTraw_responses[NoFile][eta_bin][pT_bin] = res_mean;
				pTraw_resolution[NoFile][eta_bin][pT_bin] = res_rms;
				pTraw_responses_error[NoFile][eta_bin][pT_bin] = res_mean_error;
				pTraw_resolution_error[NoFile][eta_bin][pT_bin] = res_rms_error;

				gausTF1 = FindBestGaussianCoreFit(h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]);
				pTraw_response_Gaus[NoFile][eta_bin][pT_bin] = gausTF1->GetParameter(1);
				pTraw_response_Gaus_error[NoFile][eta_bin][pT_bin] = gausTF1->GetParError(1);
			//	cout << "\n\n\nHERE"<<pTraw_response_Gaus[NoFile][eta_bin][pT_bin] << endl;
				pTraw_resolution_Gaus[NoFile][eta_bin][pT_bin] = fabs(gausTF1->GetParameter(2));
				pTraw_resolution_Gaus_error[NoFile][eta_bin][pT_bin] = gausTF1->GetParError(2);

			} // end of loop of pT bins

		}// end of loop of eta bins

		cout << " Initializing loop to calculate the JES corrected resolutions ..." << endl;
		int LoopEntries;
		if (EventsToProcess< 0 || EventsToProcess > treeEntries[NoFile]) LoopEntries = treeEntries[NoFile];
		else LoopEntries = EventsToProcess;

		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
//		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
		int LoopSize20percent = 0.2*LoopEntries;
		int StartMeasureTimeEntry =0 ;
		int TimeMeasurementEntries = 0.1* LoopSize20percent;
    		struct timeval begin, end;

		for (int i=0; i<LoopEntries; i++) //event loop
		{

			if (i<5 || i%(LoopSize20percent)==0 || i==(LoopEntries-1)) cout << " Reached entry: "<< i<< "   Loop complete :"<< 100*i/LoopEntries << "\% \n" ;
			if ( i%(LoopSize20percent)==0) //start measuring time
			{
				StartMeasureTimeEntry =i;
				gettimeofday(&begin, 0);
			}
			if (i==StartMeasureTimeEntry+TimeMeasurementEntries) //after 1000 events, stop time measurement and print 
			{
				gettimeofday(&end, 0);
				long seconds = end.tv_sec - begin.tv_sec;
				long microseconds = end.tv_usec - begin.tv_usec;
				double elapsed = seconds + microseconds*1e-6;
				cout << ",  Average time spent per event =  " << elapsed/TimeMeasurementEntries << " seconds\n\n";

			}
			tree[NoFile]->GetEntry(i);
			gen_size = gen_jpt->size();

		//	cout << "hi 1" << endl;
			int *reco_jets_matched_sequence;
			reco_jets_matched_sequence = GetRecoToGenMatchSequence(gen_eta, gen_phi, eta, phi, DR_threshold);

			for ( int j=0; j<gen_size; j++) 
			{
				if (reco_jets_matched_sequence[j]< 0 )  continue;
				if (  jpt->at(reco_jets_matched_sequence[j]) < pTlowCut || jpt->at(reco_jets_matched_sequence[j]) > pThighCut  ) continue;

				int reco_ybin = getBin(fabs(eta->at(reco_jets_matched_sequence[j])),yBnd, eta_bins);
				int reco_pTbin = getBin(jpt->at(reco_jets_matched_sequence[j]),ptBnd, pT_bins);
				double Response_correction;
				if (useGausFits) Response_correction = pTraw_response_Gaus[NoFile][reco_ybin][reco_pTbin];
				else Response_correction = pTraw_responses[NoFile][reco_ybin][reco_pTbin];

				if (reco_ybin>-1 && reco_pTbin>-1 && Response_correction>0 )
				{

					if(!useWeights ) h_pTreco_ov_pTgen_JEScorrected[NoFile][reco_ybin][reco_pTbin]->Fill( (jpt->at( reco_jets_matched_sequence[j] ) / gen_jpt->at(j)) / Response_correction );
					else 
					{
						h_pTreco_ov_pTgen_JEScorrected[NoFile][reco_ybin][reco_pTbin]->Fill( (jpt->at( reco_jets_matched_sequence[j] ) / gen_jpt->at(j)) / Response_correction , weight);
						h_pTreco_ov_pTgen_JEScorrected_noWeights[NoFile][reco_ybin][reco_pTbin]->Fill( (jpt->at( reco_jets_matched_sequence[j] ) / gen_jpt->at(j)) / Response_correction);
					}
				}
			} //end of loop on jets
			delete reco_jets_matched_sequence;
		} // end loop on events  
		f[NoFile]->Close();
	} // end loop on files





	char filename[1500];
	sprintf(filename, "%s/%s/pTreco_ov_pTgen_%s_JEScorrected_histos.root",analyzer_path,output_directory,image_name );

	TFile *fout = new TFile(filename, "recreate");
	fout->cd();

	
	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int eta_bin = 0; eta_bin<eta_bins; eta_bin++) 
		{
			for (int pT_bin = 0; pT_bin<pT_bins; pT_bin++)
			{
				h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Write();
				if(useWeights) 	h_pTreco_ov_pTgen_JEScorrected_noWeights[NoFile][eta_bin][pT_bin]->Write();
			}
		}
	}
	




}

