#include "include_functions.h"



void Plot_Resolution_and_response()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);
	gROOT->ForceStyle(kTRUE);
	gStyle->SetOptFit(0);


	char input_files[][800] ={
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_PU_nominal.root", 
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/PatatrackPixels_Candidates_ttbar_14TeV_noPU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_PU_LooseCuts.root", 
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/Patatrack_Candidates_ttbar_14TeV_LooseTrackParams_noPU_onlyZeta.root"
};

	double yBnd[] = {0.0, 1.3, 2.4, 2.7, 3.0}; 
	double ptBnd[] = {10.,  60.,   100.,   150.,   220.,   450.}; 
//	double ptBnd[] = {10.,   150.,    1500.}; 

	double pTlowCut = 0;
	double pThighCut = 5000;

	bool ResfromGaus = false;
	bool Save_Plots = true ;

	int PadColumns = 3;
	int PadRows = 2;
	int Canvas_Xpixels = 1000;
	int Canvas_Ypixels = 1000;

	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;
	const int pT_bins = sizeof(ptBnd)/sizeof(ptBnd[0])-1;
	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);

	char analyzer_path[500] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre7/src/Jakobs_producer/ProducerTest/"}; 
	int Colors[] = { 1, 4, 2 , 6, 3, 7 , 28, 46} ; // black, blue, red , magenta, green, light blue ,  brown, redish
	int MarkerStyle[] = { 8, 2, 5 , 4, 22, 21, 27, 28 } ; // 
	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/PUvsNoPU/"};                  
	char image_name[200] = {"LooseCuts"}; 
	char legend_array[NoFiles][500] = { "nominal_PU" , "nominal_noPU","LooseCuts_PU" , "LooseCuts_noPU"   };
	double DR_threshold = 0.2; //reco_gen matching threshold 



	double pTraw_responses[NoFiles][eta_bins][pT_bins], pTraw_resolution[NoFiles][eta_bins][pT_bins], pTraw_resolution_Gaus[NoFiles][eta_bins][pT_bins];
	double pT_resolution_corrected[NoFiles][eta_bins][pT_bins], pT_resolution_corrected_Gaus[NoFiles][eta_bins][pT_bins];

	double pTbinCenter[pT_bins];



	
	TGraph *gr_pT_response[NoFiles][eta_bins];
	TH1D *h_pT_res[NoFiles][eta_bins][pT_bins], *h_pT_res_JEScorrected[NoFiles][eta_bins][pT_bins];

//	cout << "eta bins = " << eta_bins;
//	cout << "pT bins = " << pT_bins;

	TF1 *gauss = new TF1("gauss","[0]*exp(-0.5*((x-[1])/[2])**2)",-10,10);     // gaussian  distribution
	gauss->SetParameter(0, 1/sqrt(2*3.14159));
	gauss->SetParameter(1,0);
	gauss->SetParameter(2,1);


	char name[256]; 
	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{


		for (int  eta_bin=0; eta_bin<eta_bins; eta_bin++)
		{
			for (int  pT_bin=0; pT_bin<pT_bins; pT_bin++)
			{
				sprintf(name,"h_pT_res_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pT_res[NoFile][eta_bin][pT_bin] = new TH1D(name, "", 50, 0., 3.5); // 40,0,1.0

				sprintf(name,"h_pT_res_JEScorrected_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin] = new TH1D(name, "", 50, 0., 3.5); // 40,0,1.0
			}
		}
	}

	for (int  pT_bin=0; pT_bin<pT_bins; pT_bin++)
	{
		pTbinCenter[pT_bin] = 0.5*(ptBnd[pT_bin+1] + ptBnd[pT_bin] );
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

 	double SumEt,gen_SumEt;
 	int nVtx;
 	const Int_t kNotDraw = 1<<9;
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


		int gen_size, reco_size;

		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		cout << " Initializing first loop over events to calculate pTreco/pTgen histos..." << endl;
		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
//		for (int i=0; i<8000; i++) //event loop
		{
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

			//	int gen_pTbin = getBin(gen_pt->at(j),ptBnd);
			//	int gen_ybin = getBin(fabs(gen_eta->at(j)),yBnd);

				if (reco_ybin>-1 && reco_pTbin>-1)
				{
					h_pT_res[NoFile][reco_ybin][reco_pTbin] ->Fill( jpt->at( reco_jets_matched_sequence[j] ) / gen_jpt->at(j) );

				}

			} //end of loop on jets
			delete reco_jets_matched_sequence;
		} // end of first loop on events  

		cout << " Calculating raw responses and resolutions for file" << input_files[NoFile] << " ...  "<< endl;

		for ( int eta_bin=0; eta_bin< eta_bins; eta_bin++ )
		{

			for (int pT_bin=0; pT_bin< pT_bins; pT_bin++ )
			{
				double res_mean = h_pT_res[NoFile][eta_bin][pT_bin]->GetMean();
				pTraw_responses[NoFile][eta_bin][pT_bin] = res_mean;

				double res_rms = h_pT_res[NoFile][eta_bin][pT_bin]->GetRMS(); //
				pTraw_resolution[NoFile][eta_bin][pT_bin] = res_rms;

				if ( !(h_pT_res[NoFile][eta_bin][pT_bin]->Integral()>0) ) 
				{
					cout << "Warning!!! Histogram empty!!  eta bin : " << eta_bins_legend[eta_bin] << "  pT : " << pT_bins_legend[pT_bin] << endl;
					continue;
				}

				gauss->SetParameter(1, res_mean);
				gauss->SetParameter(2, res_rms);
				h_pT_res[NoFile][eta_bin][pT_bin]->Fit("gauss","0","0",res_mean-1.5*res_rms,res_mean+1.5*res_rms);
				pTraw_resolution_Gaus[NoFile][eta_bin][pT_bin] = gauss->GetParameter(2);

			} // end of loop of pT bins

			gr_pT_response[NoFile][eta_bin] = new TGraph(pT_bins, pTbinCenter, pTraw_responses[NoFile][eta_bin] ); // creating response graphs

		}// end of loop of eta bins

		cout << " Initializing 2nd loop to calculate resolutions corrected by response..." << endl;


		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
//		for (int i=0; i<8000; i++) //event loop
		{
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
			//	int gen_pTbin = getBin(gen_pt->at(j),ptBnd);
			//	int gen_ybin = getBin(fabs(gen_eta->at(j)),yBnd);
			//cout << "hi 2" << endl;
				if (reco_ybin>-1 && reco_pTbin>-1 && pTraw_responses[NoFile][reco_ybin][reco_pTbin]>0 )
				{
				//	cout << "hi 3" << "raw response = "<<pTraw_responses[NoFile][reco_ybin][reco_pTbin] <<  "   " <<  gen_jpt->at(j) <<endl;
				//	cout << NoFile<<  "  " << reco_ybin << "  " << reco_pTbin << endl;
					h_pT_res_JEScorrected[NoFile][reco_ybin][reco_pTbin]->Fill( (jpt->at( reco_jets_matched_sequence[j] ) / gen_jpt->at(j)) / pTraw_responses[NoFile][reco_ybin][reco_pTbin] );
				}
			//cout << "hi 4" << endl;
			} //end of loop on jets
			delete reco_jets_matched_sequence;
		} // end of first loop on events  

		cout << " Calculating (JES corrected) resolutions for file" << input_files[NoFile] << " ...  "<< endl;

		for ( int eta_bin=0; eta_bin< eta_bins; eta_bin++ )
		{
			for (int pT_bin=0; pT_bin< pT_bins; pT_bin++ )
			{
				double res_rms = h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->GetRMS(); //
				pT_resolution_corrected[NoFile][eta_bin][pT_bin] = res_rms;

				if ( !(h_pT_res[NoFile][eta_bin][pT_bin]->Integral()>0) ) 
				{
					cout << "Warning!!! Histogram empty!!  eta bin : " << eta_bins_legend[eta_bin] << "  pT : " << pT_bins_legend[pT_bin] << endl;
					continue;
				}


				gauss->SetParameter(1, 1.);
				gauss->SetParameter(2, res_rms);
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->Fit("gauss","","",1.-1.5*res_rms,1.+1.5*res_rms);
				pT_resolution_corrected_Gaus[NoFile][eta_bin][pT_bin] = fabs(gauss->GetParameter(2));

			}
		}
	} // end of loop on files



// ======================================================plotting stuff======================================//
	char res_text[200];
	

	TCanvas *ptResponsePad = new TCanvas("ptResponsePad", "",Canvas_Xpixels,Canvas_Ypixels);
	ptResponsePad->Divide(PadColumns,PadRows);


	TCanvas *ptResolutionPad[eta_bins];
	
	for (int eta_bin = 0; eta_bin<eta_bins; eta_bin++)
	{
		sprintf(name,"ptResolutionPad_etaBin%i",eta_bin);
		ptResolutionPad[eta_bin] = new TCanvas(name, eta_bins_legend[eta_bin],Canvas_Xpixels,Canvas_Ypixels);
		ptResolutionPad[eta_bin]->Divide(PadColumns,PadRows);
	}


	TLegend *leg3[eta_bins][pT_bins];
	TPaveText *paveCMS = new TPaveText(0.45,0.95,0.5,1.0,"NDC");
	// paveCMS->AddText("CMS Preliminary L=9.2 fb^{-1} #sqrt{s} = 13 TeV");
//	paveCMS->AddText("CMS Preliminary#sqrt{s} = 13 TeV");
	paveCMS->AddText("CMS Simulation #sqrt{s} = 14 TeV");
	paveCMS->SetFillColor(0);
	paveCMS->SetBorderSize(0);
	paveCMS->SetTextSize(0.04);

	TLegend *leg1[eta_bins];

	TLegend *leg2 =new TLegend(.1, .6, .9, .9);//7899//4899
	leg2->SetTextSize(0.055);
	leg2->SetFillColor(0); 
	leg2->SetBorderSize(0);



	for (int eta_bin = 0; eta_bin<eta_bins; eta_bin++) 
	{

		leg1[eta_bin]=new TLegend(.1, .6, .9, .9);//7899//4899
		leg1[eta_bin]->SetTextSize(0.06);
		leg1[eta_bin]->SetHeader(eta_bins_legend[eta_bin],"C"); 
		leg1[eta_bin]->SetFillColor(0); 
		leg1[eta_bin]->SetBorderSize(0);


		for (int NoFile=0; NoFile<NoFiles; NoFile++)
		{

			if(eta_bin==0) leg2->AddEntry(gr_pT_response[NoFile][0], legend_array[NoFile], "L");

			const char *etaText = (yBnd[eta_bin]==0 ? Form("|#eta| < %1.2g",yBnd[eta_bin+1]) :
			Form("%1.2g #leq |#eta| < %1.2g",yBnd[eta_bin],yBnd[eta_bin+1]));
			TLatex *etaTextLatx = new TLatex(0.28,0.86,etaText); //cout<<seta<<endl;
			etaTextLatx->SetNDC();
			etaTextLatx->SetTextSize(0.06);

			ptResponsePad->cd(eta_bin+1);
			gr_pT_response[NoFile][eta_bin]->SetLineColor(Colors[NoFile]);
			gr_pT_response[NoFile][eta_bin]->SetLineWidth(2);
			gr_pT_response[NoFile][eta_bin]->SetLineStyle(1);
			gr_pT_response[NoFile][eta_bin]->SetMarkerColor(Colors[NoFile]);
			gr_pT_response[NoFile][eta_bin]->SetMarkerStyle(MarkerStyle[NoFile]);
			gr_pT_response[NoFile][eta_bin]->GetYaxis()->SetRangeUser(0.5, 1.5);
			gr_pT_response[NoFile][eta_bin]->GetYaxis()->SetTitle("Response");
			gr_pT_response[NoFile][eta_bin]->GetXaxis()->SetTitle("jet pT (GeV)");
			if ( NoFile==0 ) 	{	gr_pT_response[NoFile][eta_bin]->Draw("");	paveCMS ->Draw("same");  }
			else gr_pT_response[NoFile][eta_bin]->Draw("same LP");
			paveCMS->Draw("same");
			etaTextLatx->Draw();


			leg1[eta_bin]->AddEntry(h_pT_res_JEScorrected[NoFile][0][0], legend_array[NoFile], "L");

			for (int pT_bin = 0; pT_bin<pT_bins; pT_bin++)
			{


				const char *seta = (ptBnd[pT_bin]==0 ? Form("|y| < %3.0f",ptBnd[pT_bin+1]) :
				Form("%3.0f #leq pT < %3.0f",ptBnd[pT_bin],ptBnd[pT_bin+1]));
				TLatex *teta = new TLatex(0.28,0.86,seta); //cout<<seta<<endl;
				teta->SetNDC();
				teta->SetTextSize(0.06);
				
				if( NoFile == 0 )
				{
					leg3[eta_bin][pT_bin] =new TLegend(.2, .7, .8, .85);//7899//4899
					leg3[eta_bin][pT_bin]->SetTextSize(0.04);
					leg3[eta_bin][pT_bin]->SetFillColor(0); 
					leg3[eta_bin][pT_bin]->SetBorderSize(0);
		 		}

				
				if (!ResfromGaus) sprintf(res_text, " RES = %3.1f %%", 100*pT_resolution_corrected[NoFile][eta_bin][pT_bin]);
				else sprintf(res_text, " RES = %3.1f %%", 100*pT_resolution_corrected_Gaus[NoFile][eta_bin][pT_bin]);
				leg3[eta_bin][pT_bin]->AddEntry(h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin], res_text , "L");


				ptResolutionPad[eta_bin]->cd(pT_bin+1);
				ptResolutionPad[eta_bin]->cd(pT_bin+1)->SetLogy(1);
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->GetXaxis()->SetTitle("pT_reco / pT_gen");
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->GetYaxis()->SetTitle("Entries");
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->GetYaxis()->SetTitleOffset(1.3);

				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->SetLineColor(Colors[NoFile]); 
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->SetLineStyle(1);
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->SetMinimum(0.9);
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->SetMaximum(100000);
				h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->GetFunction("gauss")->SetBit(TF1::kNotDraw); // do not draw fitted functions
				if ( NoFile==0 ) 	{	h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->Draw("hist");	paveCMS ->Draw("same");  }
				else h_pT_res_JEScorrected[NoFile][eta_bin][pT_bin]->Draw("same");
				if ( NoFile == NoFiles-1 )	leg3[eta_bin][pT_bin]->Draw("same");
				teta->Draw();


			}  //end of pT bin loops
		} // end of file loops

		ptResolutionPad[eta_bin]->cd(pT_bins+1);	leg1[eta_bin]->Draw();

	} // end of eta bin loops

	ptResponsePad->cd(eta_bins+1); 
	leg2->Draw();

	char filename[500];
	if (Save_Plots)
	{
		sprintf(filename,"%s/%s/pTresponse_%s.png",analyzer_path,output_directory,image_name);
		ptResponsePad->SaveAs(filename);

		for (int eta_bin = 0; eta_bin<eta_bins; eta_bin++) 
		{
			sprintf(filename,"%s/%s/pTresolution_pTbinned_etaBin%i_%s.png",analyzer_path,output_directory,eta_bin,image_name);
			ptResolutionPad[eta_bin]->SaveAs(filename);
		}
	}



}


