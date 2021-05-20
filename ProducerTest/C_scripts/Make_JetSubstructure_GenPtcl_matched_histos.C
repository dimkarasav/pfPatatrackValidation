#include "include/include_functions.h"
#include "include/Input_definition.h"

void Make_JetSubstructure_GenPtcl_matched_histos()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);





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


	TH1D *h_ak8_pT[NoFiles][eta_bins], *h_ak8_phi[NoFiles][eta_bins], *h_ak8_m[NoFiles][eta_bins], *h_ak8_softDrop_m[NoFiles][eta_bins], *h_ak8_trimmed_m[NoFiles][eta_bins], *h_ak8_tau1[NoFiles][eta_bins], *h_ak8_tau2[NoFiles][eta_bins], *h_ak8_tau3[NoFiles][eta_bins], *h_ak8_tau4[NoFiles][eta_bins], *h_ak8_tau5[NoFiles][eta_bins], *h_ak8_tau21[NoFiles][eta_bins], *h_ak8_tau32[NoFiles][eta_bins], *h_ak8_tau43[NoFiles][eta_bins], *h_ak8_tau54[NoFiles][eta_bins], *h_genJet_pt[NoFiles][eta_bins], *h_genJet_mass[NoFiles][eta_bins];

	TH1D *h_ak8_eta[NoFiles], *h_minDR[NoFiles];

	TH2D *h_softDrop_m_qqDR[NoFiles][eta_bins], *h_tau21_qqDR[NoFiles][eta_bins];


	char filename[1500]; 
	char name[256]; 

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko();

	bool plot_minDR=true;


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		sprintf(name,"h_ak8_eta_%s",legend_array[NoFile]);
		h_ak8_eta[NoFile] = new TH1D(name, "", 60,-5,5);
		if (useWeights)		h_ak8_eta[NoFile]->Sumw2();

		sprintf(name,"h_minDR_%s",legend_array[NoFile]);
		h_minDR[NoFile] = new TH1D(name, "", 100, 0.0, 5.0); // 40,0,1.0
		if (useWeights)		h_minDR[NoFile]->Sumw2();

		for(Int_t h=0; h<eta_bins; h++)
		{ 
			//========== reco jets===========
			sprintf(name,"h_ak8_pT_%s_bin%i",legend_array[NoFile],h);
			h_ak8_pT[NoFile][h] = new TH1D(name, "", 500, 0., 3000); // 40,0,1.0
			if (useWeights) h_ak8_pT[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_phi_%s_bin%i",legend_array[NoFile],h);
			h_ak8_phi[NoFile][h] = new TH1D(name, "", 50, -3.14, 3.14); // 40,0,1.0
			if (useWeights) h_ak8_phi[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_m_%s_bin%i",legend_array[NoFile],h);
			h_ak8_m[NoFile][h] = new TH1D(name, "", 500, 0, 3000); // 40,0,1.0
			if (useWeights) h_ak8_m[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_softDrop_m_%s_bin%i",legend_array[NoFile],h);
			h_ak8_softDrop_m[NoFile][h] = new TH1D(name, "", 500, 0, 3000); // 40,0,1.0
			if (useWeights) h_ak8_softDrop_m[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_trimmed_m_%s_bin%i",legend_array[NoFile],h);
			h_ak8_trimmed_m[NoFile][h] = new TH1D(name, "", 125, 0, 3000); // 40,0,1.0
			if (useWeights) h_ak8_trimmed_m[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau1_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau1[NoFile][h] = new TH1D(name, "", 400, 0., 400.); // 40,0,1.0
			if (useWeights) h_ak8_tau1[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau2_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau2[NoFile][h] = new TH1D(name, "", 400, 0., 400.); // 40,0,1.0
			if (useWeights) h_ak8_tau2[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau3_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau3[NoFile][h] = new TH1D(name, "", 400, 0., 400.); // 40,0,1.0
			if (useWeights) h_ak8_tau3[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau4_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau4[NoFile][h] = new TH1D(name, "", 400, 0., 400.); // 40,0,1.0
			if (useWeights) h_ak8_tau4[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau5_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau5[NoFile][h] = new TH1D(name, "", 400, 0., 400.); // 40,0,1.0
			if (useWeights) h_ak8_tau5[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau21_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau21[NoFile][h] = new TH1D(name, "", 50, 0., 1.2); // 40,0,1.0
			if (useWeights) h_ak8_tau21[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau32_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau32[NoFile][h] = new TH1D(name, "", 50, 0., 1.2); // 40,0,1.0
			if (useWeights) h_ak8_tau32[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau43_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau43[NoFile][h] = new TH1D(name, "", 50, 0., 1.2); // 40,0,1.0
			if (useWeights) h_ak8_tau43[NoFile][h]->Sumw2();

			sprintf(name,"h_ak8_tau54_%s_bin%i",legend_array[NoFile],h);
			h_ak8_tau54[NoFile][h] = new TH1D(name, "", 50, 0., 1.2); // 40,0,1.0
			if (useWeights) h_ak8_tau54[NoFile][h]->Sumw2();

			sprintf(name,"h_genJet_pt_%s_bin%i",legend_array[NoFile],h);
			h_genJet_pt[NoFile][h] = new TH1D(name, "", 50, 0., pThighCut); // 40,0,1.0
			if (useWeights) h_genJet_pt[NoFile][h]->Sumw2();

			sprintf(name,"h_genJet_mass_%s_bin%i",legend_array[NoFile],h);
			h_genJet_mass[NoFile][h] = new TH1D(name, "", 125, 0., pThighCut); // 40,0,1.0
			if (useWeights) h_genJet_mass[NoFile][h]->Sumw2();

			sprintf(name,"h_softDrop_m_qqDR_%s_bin%i",legend_array[NoFile],h);
			h_softDrop_m_qqDR[NoFile][h] = new TH2D(name, "",  60,0, 300, 50, 0, 5.0);
			if (useWeights) 	h_softDrop_m_qqDR[NoFile][h]->Sumw2();


			sprintf(name,"h_tau21_qqDR_%s_bin%i",legend_array[NoFile],h);
			h_tau21_qqDR[NoFile][h] = new TH2D(name, "",  25,0, 1.2, 50, 0, 5.0);
			if (useWeights) 	h_tau21_qqDR[NoFile][h]->Sumw2();

		}// end of etabin loop
	} // end of file loop

	int treeEntries[NoFiles];
	TFile *f[NoFiles];
	TTree *tree[NoFiles];

	std::vector<float> *gen_mass = 0;
	std::vector<float> *gen_jpt = 0;
	std::vector<float> *gen_eta = 0;
	std::vector<float> *gen_phi = 0;

	std::vector<float> *ak8_jpt = 0;
	std::vector<float> *ak8_eta = 0;
	std::vector<float> *ak8_phi = 0;
	std::vector<float> *ak8_m = 0;
	std::vector<float> *ak8_softDrop_m = 0;
	std::vector<float> *ak8_trimmed_m = 0;
	std::vector<float> *ak8_tau1 = 0;
	std::vector<float> *ak8_tau2 = 0;
	std::vector<float> *ak8_tau3 = 0;
	std::vector<float> *ak8_tau4 = 0;
	std::vector<float> *ak8_tau5 = 0;

	std::vector<float>           *genpart_pt =0;
	std::vector<float>           *genpart_eta=0;
	std::vector<float>           *genpart_phi=0;
	std::vector<float>           *genpart_m=0;
	std::vector<int>	     *genpart_pdgid=0;
	std::vector<float>           *genpart_px=0;
	std::vector<float>           *genpart_py=0;
	std::vector<float>           *genpart_pz=0;
	std::vector<float>           *genpart_p=0;
	std::vector<float>           *genpart_E=0;
	std::vector<int>             *genpart_numDaught=0;
	std::vector<bool>            *genpart_isQfromZ=0;

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

		tree[NoFile]->SetBranchAddress("gen_mass",&gen_mass);
		tree[NoFile]->SetBranchAddress("gen_jpt",&gen_jpt);
		tree[NoFile]->SetBranchAddress("gen_eta",&gen_eta);
		tree[NoFile]->SetBranchAddress("gen_phi",&gen_phi);

		tree[NoFile]->SetBranchAddress("ak8_pt",&ak8_jpt);
		tree[NoFile]->SetBranchAddress("ak8_eta",&ak8_eta);
		tree[NoFile]->SetBranchAddress("ak8_phi",&ak8_phi);
		tree[NoFile]->SetBranchAddress("ak8_m",&ak8_m);
		tree[NoFile]->SetBranchAddress("ak8_softDrop_m",&ak8_softDrop_m);
		tree[NoFile]->SetBranchAddress("ak8_trimmed_m",&ak8_trimmed_m);
		tree[NoFile]->SetBranchAddress("ak8_tau1",&ak8_tau1);
		tree[NoFile]->SetBranchAddress("ak8_tau2",&ak8_tau2);
		tree[NoFile]->SetBranchAddress("ak8_tau3",&ak8_tau3);
		tree[NoFile]->SetBranchAddress("ak8_tau4",&ak8_tau4);
		tree[NoFile]->SetBranchAddress("ak8_tau5",&ak8_tau5);

		tree[NoFile]->SetBranchAddress("genpart_pt" ,&genpart_pt);
		tree[NoFile]->SetBranchAddress("genpart_eta",&genpart_eta);
		tree[NoFile]->SetBranchAddress("genpart_phi",&genpart_phi);
		tree[NoFile]->SetBranchAddress("genpart_m",&genpart_m);
		tree[NoFile]->SetBranchAddress("genpart_pdgid",&genpart_pdgid);
		tree[NoFile]->SetBranchAddress("genpart_px" ,&genpart_px);
		tree[NoFile]->SetBranchAddress("genpart_py" ,&genpart_py);
		tree[NoFile]->SetBranchAddress("genpart_pz" ,&genpart_pz);
		tree[NoFile]->SetBranchAddress("genpart_p" ,&genpart_p);
		tree[NoFile]->SetBranchAddress("genpart_E" ,&genpart_E);
		tree[NoFile]->SetBranchAddress("genpart_numDaught" ,&genpart_numDaught);
		tree[NoFile]->SetBranchAddress("genpart_isQfromZ" ,&genpart_isQfromZ);

		tree[NoFile]->SetBranchAddress("SumEt",&SumEt);
		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);
		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);



		int gen_size, reco_size;

		//======================= f1 jets ===================================

		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		int LoopEntries;
		if (EventsToProcess< 0 || EventsToProcess > treeEntries[NoFile]) LoopEntries = treeEntries[NoFile];
		else LoopEntries = EventsToProcess;

		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
		if ( LoopEntries != treeEntries[NoFile] )  cout<<"\nNumber of entries to be processed: "<< LoopEntries << endl;


		int LoopSize20percent = 0.2*LoopEntries; 

		for (int i=0; i<LoopEntries; i++) //event loop
//		for (int i=0; i<100000; i++) //event loop
		{

			if (i<5 || i%(LoopSize20percent)==0 || i==(LoopEntries-1)) cout << " Reached entry: "<< i<< "   Loop complete: "<< 100*i/(LoopEntries-1) << "\%" <<"\n";

			tree[NoFile]->GetEntry(i);
			reco_size = ak8_jpt->size();
			gen_size  = genpart_pt->size();

			int *reco_jets_matched_sequence;
			int *reco_jets_matched_sequence_forPlot;
			double deltaPhi;
			double QQ_DR=0;
//			reco_jets_matched_sequence = GetRecoToGenMatchSequence(genpart_eta, genpart_phi, ak8_eta, ak8_phi, DR_threshold);
//			if (plot_minDR)	reco_jets_matched_sequence_forPlot = GetRecoToGenMatchSequence(genpart_eta, genpart_phi, ak8_eta, ak8_phi, 5.0);
				
			reco_jets_matched_sequence = GetRecoToGenMatchSequenceAdvanced(genpart_eta, genpart_phi, ak8_eta, ak8_phi, DR_threshold);
			if (plot_minDR)	reco_jets_matched_sequence_forPlot = GetRecoToGenMatchSequenceAdvanced(genpart_eta, genpart_phi, ak8_eta, ak8_phi, 5.0);

			int gen_ybin = getBin(fabs(genpart_eta->at(1)),yBnd, eta_bins);
			if (fabs(genpart_pdgid->at(1)) <= 6.1 && fabs(genpart_pdgid->at(2)) <= 6.1) //make sure they are quarks
			{
				if(!useWeights) weight = 1.0;
				deltaPhi = fabs(genpart_phi->at(1)-genpart_phi->at(2));
				if(deltaPhi>TMath::Pi()) deltaPhi=2.*TMath::Pi()-deltaPhi;

				QQ_DR = sqrt( pow( genpart_eta->at(1) - genpart_eta->at(2) ,2) + pow( deltaPhi ,2)  ); //calculate DR between quarks
			}

//			if(QQ_DR>0.8) continue; //this is a test

			for ( int j=1; j<gen_size; j++) //starting from j=1, cause j=0 corresponds to Z' in these ntuples.
			{

				if(!useWeights) weight = 1.0;

				if (reco_jets_matched_sequence[j]>-0.1 && ( ak8_jpt->at(reco_jets_matched_sequence[j]) < pTlowCut || ak8_jpt->at(reco_jets_matched_sequence[j]) > pThighCut ) ) continue;
				if (reco_jets_matched_sequence_forPlot[j]>-0.1 && ( ak8_jpt->at(reco_jets_matched_sequence_forPlot[j]) < pTlowCut || ak8_jpt->at(reco_jets_matched_sequence_forPlot[j]) > pThighCut ) ) continue;

				if(plot_minDR && reco_jets_matched_sequence_forPlot[j]>=0 )
				{
					deltaPhi = fabs(genpart_phi->at(j)-ak8_phi->at(reco_jets_matched_sequence_forPlot[j]));
					if(deltaPhi>TMath::Pi()) deltaPhi=2.*TMath::Pi()-deltaPhi;					
					double deltaR = sqrt( pow( genpart_eta->at(j) - ak8_eta->at(reco_jets_matched_sequence_forPlot[j]), 2) + pow( deltaPhi , 2) );
					h_minDR[NoFile]->Fill(deltaR);
				}

				int gen_ybin = getBin(fabs(genpart_eta->at(j)),yBnd, eta_bins);
				if (reco_jets_matched_sequence[j]<0 ) continue ; // if no match was made, skip this gen jet/Ptcle

				h_ak8_eta[NoFile]->Fill( ak8_eta->at(reco_jets_matched_sequence[j]), weight );
				int ybin = getBin(fabs(ak8_eta->at(reco_jets_matched_sequence[j])),yBnd, eta_bins);


				if (ybin > -1)//fill hist's in the corresponding eta bin
				{

					h_ak8_phi[NoFile][ybin] ->Fill( ak8_phi->at(reco_jets_matched_sequence[j]), weight  );
					h_ak8_pT[NoFile][ybin]  ->Fill( ak8_jpt->at(reco_jets_matched_sequence[j]), weight  );
					h_ak8_m[NoFile][ybin]  ->Fill(  ak8_m->at(reco_jets_matched_sequence[j]), weight );
					h_ak8_softDrop_m[NoFile][ybin]->Fill( ak8_softDrop_m->at(reco_jets_matched_sequence[j]), weight );
					h_ak8_trimmed_m[NoFile][ybin]->Fill( ak8_trimmed_m->at(reco_jets_matched_sequence[j]), weight );
					h_ak8_tau1[NoFile][ybin]->Fill( ak8_tau1->at(reco_jets_matched_sequence[j]), weight );
					h_ak8_tau2[NoFile][ybin]->Fill( ak8_tau2->at(reco_jets_matched_sequence[j]), weight );
					h_ak8_tau3[NoFile][ybin]->Fill( ak8_tau3->at(reco_jets_matched_sequence[j]), weight );
					h_ak8_tau4[NoFile][ybin]->Fill( ak8_tau4->at(reco_jets_matched_sequence[j]), weight );
					h_ak8_tau5[NoFile][ybin]->Fill( ak8_tau5->at(reco_jets_matched_sequence[j]), weight );

					if (ak8_tau1->at(reco_jets_matched_sequence[j])>0) h_ak8_tau21[NoFile][ybin]->Fill( ak8_tau2->at(reco_jets_matched_sequence[j]) / ak8_tau1->at(reco_jets_matched_sequence[j]), weight  );
 					if (ak8_tau2->at(reco_jets_matched_sequence[j])>0) h_ak8_tau32[NoFile][ybin]->Fill( ak8_tau3->at(reco_jets_matched_sequence[j]) / ak8_tau2->at(reco_jets_matched_sequence[j]), weight  );
					if (ak8_tau3->at(reco_jets_matched_sequence[j])>0) h_ak8_tau43[NoFile][ybin]->Fill( ak8_tau4->at(reco_jets_matched_sequence[j]) / ak8_tau3->at(reco_jets_matched_sequence[j]), weight  );
					if (ak8_tau4->at(reco_jets_matched_sequence[j])>0) h_ak8_tau54[NoFile][ybin]->Fill( ak8_tau5->at(reco_jets_matched_sequence[j]) / ak8_tau4->at(reco_jets_matched_sequence[j]), weight  );

					if (ak8_tau1->at(reco_jets_matched_sequence[j])>0) h_tau21_qqDR[NoFile][ybin]->Fill(ak8_tau2->at(reco_jets_matched_sequence[j]) / ak8_tau1->at(reco_jets_matched_sequence[j]), QQ_DR,weight);
					
					h_softDrop_m_qqDR[NoFile][ybin]->Fill(ak8_softDrop_m->at(reco_jets_matched_sequence[j]), QQ_DR,weight);

				}//end of valid ybin if 
			} // end of Ptcl loop
			delete reco_jets_matched_sequence;
		} // end of event loop
	} // end of file loop

	sprintf(filename, "%s/%s/JetSubstructure_matched_histos_%s.root",analyzer_path,output_directory,image_name);
	TFile *fout = new TFile(filename, "recreate");
	fout->cd();

	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		h_ak8_eta[NoFile]->Write();
		h_minDR[NoFile]->Write();

		for(int eta_bin=0; eta_bin<eta_bins; eta_bin++)
		{
			h_ak8_phi[NoFile][eta_bin]->Write();
			h_ak8_pT[NoFile][eta_bin]->Write();
			h_ak8_m[NoFile][eta_bin]->Write();
			h_ak8_softDrop_m[NoFile][eta_bin]->Write();
			h_ak8_trimmed_m[NoFile][eta_bin]->Write();
			h_ak8_tau1[NoFile][eta_bin]->Write();
			h_ak8_tau2[NoFile][eta_bin]->Write();
			h_ak8_tau3[NoFile][eta_bin]->Write();
			h_ak8_tau4[NoFile][eta_bin]->Write();
			h_ak8_tau5[NoFile][eta_bin]->Write();
			h_ak8_tau21[NoFile][eta_bin]->Write();
			h_ak8_tau32[NoFile][eta_bin]->Write();
			h_ak8_tau43[NoFile][eta_bin]->Write();
			h_ak8_tau54[NoFile][eta_bin]->Write();
			h_genJet_pt[NoFile][eta_bin]->Write();
			h_genJet_mass[NoFile][eta_bin]->Write();
			h_tau21_qqDR[NoFile][eta_bin]->Write();
			h_softDrop_m_qqDR[NoFile][eta_bin]->Write();
		}
	}

}

