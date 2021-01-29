#include "include/include_functions.h"
#include "include/Input_definition.h"


void Make_pTreco_ov_pTgen_histos()
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


	if (pThighCut > ptBnd[pT_bins]) pThighCut = ptBnd[pT_bins];
	if (pTlowCut  < ptBnd[0])         pTlowCut  = ptBnd[0];  

	cout << ptBnd[pT_bins] << endl;

	TH1D *h_pTreco_ov_pTgen[NoFiles][eta_bins][pT_bins];

	char name[256]; 
	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int  eta_bin=0; eta_bin<eta_bins; eta_bin++)
		{
			for (int  pT_bin=0; pT_bin<pT_bins; pT_bin++)
			{
				sprintf(name,"h_pTreco_ov_pTgen_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin] = new TH1D(name, "", 600, 0., 6.0); // 40,0,1.0
				if (useWeights) h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->Sumw2();

			}
		}
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

 	int nVtx;


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		f[NoFile] = TFile::Open(input_files[NoFile],"READ"); //read the .root files
		tree[NoFile] = (TTree*)f[NoFile]->Get("dijetscouting/tree"); // get the trees from the files

		tree[NoFile]->SetBranchAddress("jpt",&jpt);
		tree[NoFile]->SetBranchAddress("eta",&eta);
		tree[NoFile]->SetBranchAddress("phi",&phi);

		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);
		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);

		tree[NoFile]->SetBranchAddress("gen_jpt",&gen_jpt);
		tree[NoFile]->SetBranchAddress("gen_eta",&gen_eta);
		tree[NoFile]->SetBranchAddress("gen_phi",&gen_phi);



		int gen_size, reco_size;

		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		int LoopEntries;
		if (EventsToProcess< 0 || EventsToProcess > treeEntries[NoFile] ) LoopEntries = treeEntries[NoFile];
		else LoopEntries = EventsToProcess;
		cout << " Initializing first loop over events to calculate pTreco/pTgen histos..." << endl;
		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
//		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
		int LoopSize20percent = 0.2*LoopEntries; 
		for (int i=0; i<LoopEntries; i++) //event loop
		{

			if (i<5 || i%(LoopSize20percent)==0 || i==(LoopEntries-1)) cout << " Reached entry: "<< i<< "   Loop complete: "<< 100*i/(LoopEntries-1) << "\%" <<"\n";
			tree[NoFile]->GetEntry(i);
			gen_size = gen_jpt->size();

			int *reco_jets_matched_sequence;
			reco_jets_matched_sequence = GetRecoToGenMatchSequence(gen_eta, gen_phi, eta, phi, DR_threshold);

			for ( int j=0; j<gen_size; j++) 
			{
				if (reco_jets_matched_sequence[j]< 0 )  continue; // skip unmatched reco jets
				if ( jpt->at(reco_jets_matched_sequence[j]) < pTlowCut || jpt->at(reco_jets_matched_sequence[j]) > pThighCut ) continue;

				int reco_ybin = getBin(fabs(eta->at(reco_jets_matched_sequence[j])),yBnd, eta_bins);
				int reco_pTbin = getBin(jpt->at(reco_jets_matched_sequence[j]),ptBnd, pT_bins);

				int gen_pTbin = getBin(gen_jpt->at(j), ptBnd, eta_bins);
				int gen_ybin = getBin(fabs(gen_eta->at(j)), yBnd, pT_bins);
				int ybin, pTbin;
				if(useRecoBins) 
				{
					ybin  = reco_ybin;
					pTbin = reco_pTbin;
				}
				else
				{
					ybin  = gen_ybin;
					pTbin = gen_pTbin;
				}
				

				if (ybin>-1 && pTbin>-1)
				{
					if(!useWeights)		h_pTreco_ov_pTgen[NoFile][ybin][pTbin] ->Fill( jpt->at( reco_jets_matched_sequence[j] ) / gen_jpt->at(j) );
					else			h_pTreco_ov_pTgen[NoFile][ybin][pTbin] ->Fill( jpt->at( reco_jets_matched_sequence[j] ) / gen_jpt->at(j), weight );
				}

			} //end of loop on jets
			delete reco_jets_matched_sequence;
		} // end of first loop on events  

	f[NoFile]->Close();
	} // end of loop on files



	char filename[1500];
	sprintf(filename, "%s/%s/pTreco_ov_pTgen_%s_histos.root",analyzer_path,output_directory,image_name );

	TFile *fout = new TFile(filename, "recreate");
	fout->cd();

	
	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int eta_bin = 0; eta_bin<eta_bins; eta_bin++) 
		{
			for (int pT_bin = 0; pT_bin<pT_bins; pT_bin++)
			{
				h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->Write();
			}
		}
	}
	




}


