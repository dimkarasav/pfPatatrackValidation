#include "include_functions.h"




void PFCandidate_comparison_matched()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
//	gStyle->SetOptStat(0);

	


	char analyzer_path[200] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre7/src/Jakobs_producer/ProducerTest/"}; 
//	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/noPU/scanning_ptError_PFmodule/"};         
	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/PUvsNoPU//"};                  
	char image_name[150] = {"LooseCuts"}; 



	const int NoFiles = 4;
	const int eta_bins = 4;

//	int Colors[NoFiles] = { 1, 4, 2 , 6, 3 , 68} ; // black, blue, red , magenta
	int Colors[NoFiles] = { 1, 4, 2 , 6 } ; // black, blue, red , magenta
	bool scale_on = false ;
	bool Save_Plots = true; 

//	char input_files[NoFiles][500] ={"FullTracking_Candidates_ttbar_14TeV_PU.root","LegacyPixelsTracks_Candidates_ttbar_14TeV_PU.root","PatatrackPixels_Candidates_ttbar_14TeV_PU.root" };
	char input_files[NoFiles][500] ={
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_LooseTrackParams_noPU_onlyZeta.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_noPU_LooseCuts_DptOverPtCut2.root" , 
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_noPU_LooseCuts_DptOverPtCut1.root" , 
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_noPU_LooseCuts_DptOverPtCut0p8.root" , 
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_noPU_LooseCuts_DptOverPtCut1_ptError10.root"  

//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_LooseTrackParams_noPU_onlyZeta.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/PatatrackPixels_Candidates_ttbar_14TeV_noPU.root",
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_noPU_LooseCuts_ptError45.root" , 
//"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_noPU_LooseCuts_ptError10.root" 
 
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_PU_nominal.root", 
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/PatatrackPixels_Candidates_ttbar_14TeV_noPU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_PU_LooseCuts.root", 
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_LooseTrackParams_noPU_onlyZeta.root"
};

	char legend_array[NoFiles][500] = { "nominal_PU" , "nominal_noPU","LooseCuts_PU" , "LooseCuts_noPU"   };
//	char legend_array[NoFiles][500] = { "FullTracking" , "TrackAlgoCut5", "TrackAlgoCut2", "TrackAlgoCut1", "TrackAlgoCut0p8" ,"TrackCut1-ptError10"  };
//	char legend_array[NoFiles][500] = { "FullTracking" , "TrackAlgoCut5", "TrackAlgoCut2", "TrackAlgoCut1"  };

//	TFile *f1 = TFile::Open("PatatrackPixels_Candidates_ttbar_14TeV_noPU.root","READ"); //full tracking file

	TFile *f[NoFiles];
	TTree *tree[NoFiles];

//	TTree *tree1 = (TTree*)f1->Get("dijetscouting/tree");


	double pTmax = 300;
	double DR_threshold = 0.2;
//	double yBnd[eta_bins+1]={0.0, 1.3, 2.4, 2.7, 3.0, 5.0}; 
	double yBnd[eta_bins+1]={0.0, 1.3, 2.4, 2.7, 3.0}; 
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

	TH1D *h_ChargedHadron_pT[NoFiles][eta_bins], *h_NeutralHadron_pT[NoFiles][eta_bins], *h_ChargedHadron_SumpT[NoFiles][eta_bins], *h_NeutralHadron_SumpT[NoFiles][eta_bins], *h_NumberOfNeutralHadrons[NoFiles][eta_bins], *h_NumberOfChargedHadrons[NoFiles][eta_bins];

	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int iy=0; iy< eta_bins; iy++)
		{
	
			sprintf(name,"h_ChargedHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_pT[NoFile][iy] = new TH1D(name, "", 50,0,0.66*pTmax);

				
			sprintf(name,"h_NeutralHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_pT[NoFile][iy] = new TH1D(name, "",  50,0,0.66*pTmax);

	
			sprintf(name,"h_ChargedHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_SumpT[NoFile][iy] = new TH1D(name, "", 50,0,pTmax);

			sprintf(name,"h_NeutralHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_SumpT[NoFile][iy] = new TH1D(name, "", 50,0,pTmax);

			sprintf(name,"h_NumberOfChargedHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfChargedHadrons[NoFile][iy] = new TH1D(name, "", 50,0,50);

			sprintf(name,"h_NumberOfNeutralHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfNeutralHadrons[NoFile][iy] = new TH1D(name, "",  50,0,50);
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

			tree[NoFile]->GetEntry(i);

			reco_size = jpt->size();
			gen_size  = gen_jpt->size();

			int *reco_jets_matched_sequence;	
			reco_jets_matched_sequence = GetRecoToGenMatchSequence(gen_eta, gen_phi, eta, phi, DR_threshold);

//			cout << "\n Processing event " << i << " number of reco_jets : "<< reco_size <<endl; 
			for ( int j=0; j<gen_size; j++) // loop over jets
			{

				if (reco_jets_matched_sequence[j]<0 ) continue ; // if no match was made, skip this gen jet
	//			int cand_size = cand_vector->size();
				int jet_ybin = getBin(fabs(eta->at(reco_jets_matched_sequence[j])),yBnd, eta_bins); //get the jet eta bin

				int cand_size = ( jetconstituents->at( reco_jets_matched_sequence[j] ) ).size();
//				cout << "jet "<< reco_jets_matched_sequence[j] << " cand size = " << cand_size << endl;

				float chargedCand_SumpT = 0;
				float neutralCand_SumpT = 0;
				int neutralPFcand_number = 0;
				int chargedPFcand_number = 0;

				for (int l=0; l< cand_size; l++)
				{
					int jet_candidate = (jetconstituents->at(reco_jets_matched_sequence[j])).at(l) ; //get the jet candidate number
					//cout << "\njet candidate number: " << jet_candidate << endl;

					int cand_pdgID = pfcandpdgid[jet_candidate];
					float cand_eta = pfcandeta[jet_candidate];
					float cand_phi = pfcandphi[jet_candidate];
					float cand_pT = pfcandpt[jet_candidate];

	//				cout << "\n candidate pdgID = " << cand_pdgID << endl;


					int cand_ybin = getBin( fabs(pfcandeta[ jet_candidate ]),yBnd, eta_bins );  //get the pf candidate eta bin

					if (jet_ybin > -1)//fill hist's in the corresponding eta bin
					{

						if ( cand_pdgID==-211 || cand_pdgID==211 || cand_pdgID==11 || cand_pdgID==-11 || cand_pdgID==13 || cand_pdgID==-13) chargedPFcand_number = chargedPFcand_number +1;
						if ( cand_pdgID==-211 || cand_pdgID==211 )
						{
							chargedCand_SumpT = chargedCand_SumpT + cand_pT; // get the pT sum of charged candidates
							h_ChargedHadron_pT[NoFile][jet_ybin]->Fill(cand_pT); // fill charged hadron pT's
						}
						if ( cand_pdgID==130 || cand_pdgID==22 || cand_pdgID==0 ) 	neutralPFcand_number = neutralPFcand_number +1;
						if ( cand_pdgID==130 )
						{
							neutralCand_SumpT = neutralCand_SumpT + cand_pT; // get the pT sum of charged candidates
							h_NeutralHadron_pT[NoFile][jet_ybin]->Fill(cand_pT); // fill neutral hadron pT's
						}
					}
			//		cout << " cand " << cand_vector->at(l) << " , pt = " << pfcandpt[cand_vector->at(l)] << endl;  
				} // end of candidate loop


				h_NumberOfNeutralHadrons[NoFile][jet_ybin]->Fill(neutralPFcand_number);
				h_NumberOfChargedHadrons[NoFile][jet_ybin]->Fill(chargedPFcand_number);
				h_ChargedHadron_SumpT[NoFile][jet_ybin]->Fill(chargedCand_SumpT);
				h_NeutralHadron_SumpT[NoFile][jet_ybin]->Fill(neutralCand_SumpT);

			} //  end of jet loop

			gen_size = 0.;
			reco_size = 0;
			delete reco_jets_matched_sequence;

		} // end of event loop


		delete tree[NoFile];
		delete f[NoFile];
	} // end of file loop



// =========================================================== Plotting staff ===================================================================

	TCanvas *pad1 = new TCanvas("pad1", "",1);
	pad1->Divide(3,2);
	TCanvas *pad2 = new TCanvas("pad2", "",1);
	pad2->Divide(3,2);
	TCanvas *pad3 = new TCanvas("pad3", "",1);
	pad3->Divide(3,2);
	TCanvas *pad4 = new TCanvas("pad4", "",1);
	pad4->Divide(3,2);
	TCanvas *pad5 = new TCanvas("pad5", "",1);
	pad5->Divide(3,2);
	TCanvas *pad6 = new TCanvas("pad6", "",1);
	pad6->Divide(3,2);

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
			h_ChargedHadron_SumpT[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 1000000);
			h_ChargedHadron_SumpT[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_ChargedHadron_SumpT[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_ChargedHadron_SumpT[NoFile][iy]->Scale(h_ChargedHadron_SumpT[0][iy]->Integral()/h_ChargedHadron_SumpT[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_ChargedHadron_SumpT[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_ChargedHadron_SumpT[NoFile][iy]->Draw("same");
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
			h_NeutralHadron_SumpT[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 100000);
			h_NeutralHadron_SumpT[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NeutralHadron_SumpT[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NeutralHadron_SumpT[NoFile][iy]->Scale(h_NeutralHadron_SumpT[0][iy]->Integral()/h_NeutralHadron_SumpT[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NeutralHadron_SumpT[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NeutralHadron_SumpT[NoFile][iy]->Draw("same");
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
			h_ChargedHadron_pT[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 1000000);
			h_ChargedHadron_pT[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_ChargedHadron_pT[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_ChargedHadron_pT[NoFile][iy]->Scale(h_ChargedHadron_pT[0][iy]->Integral()/h_ChargedHadron_pT[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_ChargedHadron_pT[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_ChargedHadron_pT[NoFile][iy]->Draw("same");
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
			h_NeutralHadron_pT[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 100000);
			h_NeutralHadron_pT[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NeutralHadron_pT[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NeutralHadron_pT[NoFile][iy]->Scale(h_NeutralHadron_pT[0][iy]->Integral()/h_NeutralHadron_pT[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NeutralHadron_pT[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NeutralHadron_pT[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			


			if(eta_bins > 1) pad5->cd(iy+1);
			else pad5->cd();
			if(eta_bins > 1) pad5->cd(iy+1)->SetLogy(1);
			else pad5->cd()->SetLogy(1);
			h_NumberOfChargedHadrons[NoFile][iy]->GetXaxis()->SetTitle("# of charged candidates per jet (GeV) ");
			h_NumberOfChargedHadrons[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_NumberOfChargedHadrons[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NumberOfChargedHadrons[NoFile][iy]->SetLineWidth(2);
		//	h_NumberOfChargedHadrons[NoFile][iy]->SetMinimum(0.1);
		//	h_NumberOfChargedHadrons[NoFile][iy]->SetMaximum(100000);
			h_NumberOfChargedHadrons[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 100000);
			h_NumberOfChargedHadrons[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NumberOfChargedHadrons[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NumberOfChargedHadrons[NoFile][iy]->Scale(h_NumberOfChargedHadrons[0][iy]->Integral()/h_NumberOfChargedHadrons[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NumberOfChargedHadrons[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NumberOfChargedHadrons[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");




			if(eta_bins > 1) pad6->cd(iy+1);
			else pad6->cd();
			if(eta_bins > 1) pad6->cd(iy+1)->SetLogy(1);
			else pad6->cd()->SetLogy(1);
			h_NumberOfNeutralHadrons[NoFile][iy]->GetXaxis()->SetTitle("# of neutral candidates per jet (GeV) ");
			h_NumberOfNeutralHadrons[NoFile][iy]->GetYaxis()->SetTitle("# jets");
			h_NumberOfNeutralHadrons[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NumberOfNeutralHadrons[NoFile][iy]->SetLineWidth(2);
		//	h_NumberOfNeutralHadrons[NoFile][iy]->SetMinimum(0.1);
		//	h_NumberOfNeutralHadrons[NoFile][iy]->SetMaximum(100000);
			h_NumberOfNeutralHadrons[NoFile][iy]->GetYaxis()->SetRangeUser(0.1, 100000);
			h_NumberOfNeutralHadrons[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NumberOfNeutralHadrons[NoFile][iy]->SetMarkerSize(0.04);
			if(scale_on) h_NumberOfNeutralHadrons[NoFile][iy]->Scale(h_NumberOfNeutralHadrons[0][iy]->Integral()/h_NumberOfNeutralHadrons[NoFile][iy]->Integral());
			if ( NoFile==0 ) 
			{
				h_NumberOfNeutralHadrons[NoFile][iy]->Draw("hist");
				paveCMS ->Draw("same");
			}
			else h_NumberOfNeutralHadrons[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");




			if(eta_bins > 1 &&  NoFile == NoFiles-1 )
			{
				pad1->cd(eta_bins+1);    leg1->Draw();
				pad2->cd(eta_bins+1);    leg1->Draw();
				pad3->cd(eta_bins+1);    leg1->Draw();
				pad4->cd(eta_bins+1);    leg1->Draw();
				pad5->cd(eta_bins+1);    leg1->Draw();
				pad6->cd(eta_bins+1);    leg1->Draw();
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

	}

}

