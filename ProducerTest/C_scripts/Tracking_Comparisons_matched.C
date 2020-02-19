#include "include_functions.h"


void Tracking_Comparisons_matched()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);

	const int NoFiles = 3;
	const int eta_bins = 4;

	char input_files[NoFiles][500] ={"FullTracking_Candidates_ttbar_14TeV_PU.root","Patatrack_Candidates_ttbar_14TeV_ParamCheck1_onlyZeta.root","PatatrackPixels_Candidates_ttbar_14TeV_PU.root" };

	double pTlowCut = 0;
	double pThighCut = 5000;
	double pThistoMax = 350;

	int PadColumns = 3;
	int PadRows = 2;
	int Canvas_Xpixels = 1000;
	int Canvas_Ypixels = 1000;

	char analyzer_path[500] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre7/src/Jakobs_producer/ProducerTest/"}; 
//	char *output_directory = "Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/JetpT0_70"; 
	char output_directory[200] = {"Realistic_RunIII_14TeV/ttbar/Relaxing_ZetaCuts/"};         
	char image_name[200] = {"Comparisons_PatatrackVsLegacyVsFull"}; 

	char legend_array[NoFiles][500] = { "FullTracking" , "LooseZetaCuts", "PatatrackPixelTracking"  };


	int Colors[NoFiles] = { 1, 4, 2 } ; // black, blue, red , magenta
	bool scale_histos = false;
	bool Save_Plots = true ;
	double DR_threshold = 0.2;

//	double yBnd[eta_bins+1]={0.0, 1.3, 2.4, 2.7, 3.0, 5.0}; 
	double yBnd[eta_bins+1]={0.0, 1.3, 2.4, 2.7, 3.0}; 

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

	TH1D *h_ETAJet[NoFiles], *h_SumEt[NoFiles];

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

	

		for(Int_t h=0; h<eta_bins;h++)
		{ 
		//========== reco jets===========
			sprintf(name,"h_METovSUMET_%s_bin%i",legend_array[NoFile],h);
			h_METovSUMET[NoFile][h] = new TH1D(name, "", 50, 0, 1); // 40,0,1.0

			sprintf(name,"h_CHFJet_%s_bin%i",legend_array[NoFile],h);
			h_CHFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2); //40, 0,1.2

			sprintf(name,"h_NHFJet_%s_bin%i",legend_array[NoFile],h);
			h_NHFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2);

			sprintf(name,"h_CEMFJet_%s_bin%i",legend_array[NoFile],h);
			h_CEMFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2);

			sprintf(name,"h_NEMFJet_%s_bin%i",legend_array[NoFile],h);
			h_NEMFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2);

			sprintf(name,"h_MUFJet_%s_bin%i",legend_array[NoFile],h);
			h_MUFJet[NoFile][h] = new TH1D(name, "", 60, 0, 1.2);


			sprintf(name,"h_CMJet_%s_bin%i",legend_array[NoFile],h);
			h_CMJet[NoFile][h] = new TH1D(name, "", 80, 0, 80);
		
			sprintf(name,"h_PHIJet_%s_bin%i",legend_array[NoFile],h);
			h_PHIJet[NoFile][h] = new TH1D(name, "", 40,-4,4);

			sprintf(name,"h_NMJet_%s_bin%i",legend_array[NoFile],h);
			h_NMJet[NoFile][h] = new TH1D(name, "", 80, 0, 80);

			sprintf(name,"h_ptJet_%s_bin%i",legend_array[NoFile],h);
			h_ptJet[NoFile][h] = new TH1D(name, "", 50,0,pThistoMax);

			sprintf(name,"h_pT_res_%s_bin%i",legend_array[NoFile],h);
			h_pT_res[NoFile][h] = new TH1D(name, "", 50, 0., 3.5); // 40,0,1.0

			sprintf(name,"h_gen_pT_all_%s_bin%i",legend_array[NoFile],h);
			h_gen_pT_all[NoFile][h] = new TH1D(name, "", 50, 0., pThistoMax); // 40,0,1.0

			sprintf(name,"h_gen_pT_%s_bin%i",legend_array[NoFile],h);
			h_gen_pT[NoFile][h] = new TH1D(name, "", 50, 0., pThistoMax); // 40,0,1.0

			sprintf(name,"h_reco_pT_unmatched_%s_bin%i",legend_array[NoFile],h);
			h_reco_pT_unmatched[NoFile][h] = new TH1D(name, "", 50, 0., pThistoMax); // 40,0,1.0

			sprintf(name,"h_reco_pT_all_%s_bin%i",legend_array[NoFile],h);
			h_reco_pT_all[NoFile][h] = new TH1D(name, "", 50, 0., pThistoMax); // 40,0,1.0

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


		int size, reco_size;

		//======================= f1 jets ===================================

		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;
//		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
		for (int i=0; i<8000; i++) //event loop
		{
			tree[NoFile]->GetEntry(i);
			size = gen_jpt->size();

			int *reco_jets_matched_sequence;

		
			reco_jets_matched_sequence = GetRecoToGenMatchSequence(gen_eta, gen_phi, eta, phi, DR_threshold);


			if ( size > 0 )
			{
				for ( int j=0; j<size; j++) 
				{


				//	cout << " reco matched jet = " << reco_jets_matched_sequence[j] << endl; ;
					if (reco_jets_matched_sequence[j]>-0.1 && ( jpt->at(reco_jets_matched_sequence[j]) < pTlowCut || jpt->at(reco_jets_matched_sequence[j]) > pThighCut ) ) continue;
					int gen_ybin = getBin(fabs(gen_eta->at(j)),yBnd, eta_bins);
					if (gen_ybin > -1) h_gen_pT_all[NoFile][gen_ybin]->Fill( gen_jpt->at(j) ); // this is filled only with every gen jet.


					if (j==0) h_SumEt[NoFile]->Fill(SumEt);
					if (reco_jets_matched_sequence[j]<0 ) continue ; // if no match was made, skip this gen jet

					if (gen_ybin > -1) h_gen_pT[NoFile][gen_ybin]->Fill( gen_jpt->at(j) ); // this is filled only with matched gen jets.


					h_ETAJet[NoFile]->Fill( eta->at(reco_jets_matched_sequence[j]) );
					int ybin = getBin(fabs(eta->at(reco_jets_matched_sequence[j])),yBnd, eta_bins);
					if (ybin > -1)//fill hist's in the corresponding eta bin
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
						h_pT_res[NoFile][ybin] ->Fill( jpt->at( reco_jets_matched_sequence[j] ) / gen_jpt->at(j) );

					} 
				}
			}


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
		}

	} // end of file loop



//======================================= create efficiency & fake rate assymetric error graphs ===============================
	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for(int iy=0; iy<eta_bins; iy++)
		{
			Reco_eff[NoFile][iy] = GetEfficiencyGraph(h_gen_pT[NoFile][iy] , h_gen_pT_all[NoFile][iy] );
			Fake_rate[NoFile][iy] = GetEfficiencyGraph(h_reco_pT_unmatched[NoFile][iy],h_reco_pT_all[NoFile][iy]  ) ; 
		} 
	}

	TCanvas *pad1 = new TCanvas("pad1", "",Canvas_Xpixels,Canvas_Ypixels);
	pad1->Divide(PadColumns,PadRows);
	TCanvas *pad2 = new TCanvas("pad2", "",Canvas_Xpixels,Canvas_Ypixels);
	pad2->Divide(PadColumns,PadRows);
	TCanvas *pad3 = new TCanvas("pad3", "",Canvas_Xpixels,Canvas_Ypixels);
	pad3->Divide(PadColumns,PadRows);
	TCanvas *pad4 = new TCanvas("pad4", "",Canvas_Xpixels,Canvas_Ypixels);
	pad4->Divide(PadColumns,PadRows);
	TCanvas *pad5 = new TCanvas("pad5", "",Canvas_Xpixels,Canvas_Ypixels);
	pad5->Divide(PadColumns,PadRows);
	TCanvas *pad6 = new TCanvas("pad6", "",Canvas_Xpixels,Canvas_Ypixels);
	pad6->Divide(PadColumns,PadRows);
	TCanvas *pad7 = new TCanvas("pad7", "",Canvas_Xpixels,Canvas_Ypixels);
	pad7->Divide(PadColumns,PadRows);
	TCanvas *pad8 = new TCanvas("pad8", "",Canvas_Xpixels,Canvas_Ypixels);
	pad8->Divide(PadColumns,PadRows);
	TCanvas *pad9 = new TCanvas("pad9", "",Canvas_Xpixels,Canvas_Ypixels);
	pad9->Divide(PadColumns,PadRows);
	TCanvas *pad10 = new TCanvas("pad10", "",Canvas_Xpixels,Canvas_Ypixels);
	pad10->Divide(PadColumns,PadRows);
	TCanvas *pad11 = new TCanvas("pad11", "",Canvas_Xpixels,Canvas_Ypixels);
	pad11->Divide(PadColumns,PadRows);
	TCanvas *pad12 = new TCanvas("pad12", "",Canvas_Xpixels,Canvas_Ypixels);
	pad12->Divide(PadColumns,PadRows);

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


	


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		leg1->AddEntry(h_CHFJet[NoFile][0], legend_array[NoFile], "L");
		leg2->AddEntry(h_CHFJet[NoFile][0], legend_array[NoFile], "L");

		c_eta->cd();
		c_eta->SetLogy(1);
		h_ETAJet[NoFile]->GetXaxis()->SetTitle("Jet eta");
		h_ETAJet[NoFile]->GetYaxis()->SetTitle("Entries");
		h_ETAJet[NoFile]->GetYaxis()->SetTitleOffset(1.3);
		h_ETAJet[NoFile]->SetLineColor(Colors[NoFile]); 
		h_ETAJet[NoFile]->SetLineStyle(1);
		h_ETAJet[NoFile]->SetMinimum(0.1);
		h_ETAJet[NoFile]->SetMaximum(100000);
		if (scale_histos && h_ETAJet[NoFile]->Integral()>0 ) h_ETAJet[NoFile]->Scale(h_ETAJet[0]->Integral() /h_ETAJet[NoFile]->Integral());

		if (NoFile==0)	
		{
			h_ETAJet[NoFile]->Draw("");
			paveCMS ->Draw("same");
			leg2->Draw("same");
		}
		else h_ETAJet[NoFile]->Draw("same hist");

		


		c_SumEt->cd();
		c_SumEt->SetLogy(1);
		h_SumEt[NoFile]->GetXaxis()->SetTitle("SumEt");
		h_SumEt[NoFile]->GetYaxis()->SetTitle("Entries");
		h_SumEt[NoFile]->GetYaxis()->SetTitleOffset(1.3);
		h_SumEt[NoFile]->SetLineColor(Colors[NoFile]); 
		h_SumEt[NoFile]->SetLineStyle(1);
		h_SumEt[NoFile]->SetMinimum(0.1);
		h_SumEt[NoFile]->SetMaximum(100000);
		if (scale_histos && h_SumEt[NoFile]->Integral()>0 ) h_SumEt[NoFile]->Scale(h_SumEt[0]->Integral() /h_SumEt[NoFile]->Integral());
		if (NoFile==0)	
		{
			h_SumEt[NoFile]->Draw("");
			paveCMS ->Draw("same");
			leg2->Draw("same");
		}
		else h_SumEt[NoFile]->Draw("same hist");

		for(int iy=0; iy<eta_bins; iy++)
		{
			
			
			double etamin = yBnd[iy];
			double etamax = yBnd[iy+1];
			const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
			TLatex *teta = new TLatex(0.38,0.86,seta); //cout<<seta<<endl;
			teta->SetNDC();
			teta->SetTextSize(0.06);

			pad1->cd(iy+1);
			pad1->cd(iy+1)->SetLogy(1);
			h_CHFJet[NoFile][iy]->GetXaxis()->SetTitle("Charged Hadron Fraction");
			h_CHFJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_CHFJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_CHFJet[NoFile][iy]->SetLineColor(Colors[NoFile]); 
			h_CHFJet[NoFile][iy]->SetLineStyle(1);
			h_CHFJet[NoFile][iy]->SetMinimum(0.1);
			h_CHFJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_CHFJet[NoFile][iy]->Integral()>0 ) h_CHFJet[NoFile][iy]->Scale(h_CHFJet[0][iy]->Integral()/h_CHFJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) 	{	h_CHFJet[NoFile][iy]->Draw("hist");		paveCMS ->Draw("same");	}
			else h_CHFJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");



			pad2->cd(iy+1);
			pad2->cd(iy+1)->SetLogy(1);
			h_NHFJet[NoFile][iy]->GetXaxis()->SetTitle("Neutral Hadron Fraction");
			h_NHFJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_NHFJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NHFJet[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_NHFJet[NoFile][iy]->SetLineStyle(1);
			h_NHFJet[NoFile][iy]->SetMinimum(0.1);
			h_NHFJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_NHFJet[NoFile][iy]->Integral()>0 ) h_NHFJet[NoFile][iy]->Scale(h_NHFJet[0][iy]->Integral()/h_NHFJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) {	h_NHFJet[NoFile][iy]->Draw("hist");				paveCMS ->Draw("same");	}
			else h_NHFJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			pad3->cd(iy+1);
			pad3->cd(iy+1)->SetLogy(1);
			h_CEMFJet[NoFile][iy]->GetXaxis()->SetTitle("Charged E/M Fraction");
			h_CEMFJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_CEMFJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_CEMFJet[NoFile][iy]->SetLineColor(Colors[NoFile]); 

			h_CEMFJet[NoFile][iy]->SetLineStyle(1);
			h_CEMFJet[NoFile][iy]->SetMinimum(0.1);
			h_CEMFJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_CEMFJet[NoFile][iy]->Integral()>0 ) h_CEMFJet[NoFile][iy]->Scale(h_CEMFJet[0][iy]->Integral()/h_CEMFJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) 	{	h_CEMFJet[NoFile][iy]->Draw("hist");	paveCMS ->Draw("same");  }
			else h_CEMFJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			pad4->cd(iy+1);
			pad4->cd(iy+1)->SetLogy(1);
			h_NEMFJet[NoFile][iy]->GetXaxis()->SetTitle("Neutral E/M Fraction");
			h_NEMFJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_NEMFJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NEMFJet[NoFile][iy]->SetLineColor(Colors[NoFile]); 
			h_NEMFJet[NoFile][iy]->SetLineStyle(1);
			h_NEMFJet[NoFile][iy]->SetMinimum(0.1);
			h_NEMFJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_NEMFJet[NoFile][iy]->Integral()>0 ) h_NEMFJet[NoFile][iy]->Scale(h_NEMFJet[0][iy]->Integral()/h_NEMFJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) 	{	h_NEMFJet[NoFile][iy]->Draw("hist");	paveCMS ->Draw("same");  }
			else h_NEMFJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");

			pad5->cd(iy+1);
			pad5->cd(iy+1)->SetLogy(1);
			h_MUFJet[NoFile][iy]->GetXaxis()->SetTitle("Muon Fraction");
			h_MUFJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_MUFJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_MUFJet[NoFile][iy]->SetLineColor(Colors[NoFile]); 
			h_MUFJet[NoFile][iy]->SetLineStyle(1);
			h_MUFJet[NoFile][iy]->SetMinimum(0.1);
			h_MUFJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_MUFJet[NoFile][iy]->Integral()>0 ) h_MUFJet[NoFile][iy]->Scale(h_MUFJet[0][iy]->Integral()/h_MUFJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) 	{	h_MUFJet[NoFile][iy]->Draw("hist");	paveCMS ->Draw("same");  }
			else h_MUFJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");

			pad6->cd(iy+1);
			pad6->cd(iy+1)->SetLogy(1);
			h_PHIJet[NoFile][iy]->GetXaxis()->SetTitle("Jet phi");
			h_PHIJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_PHIJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_PHIJet[NoFile][iy]->SetLineColor(Colors[NoFile]);
			h_PHIJet[NoFile][iy]->SetLineStyle(1);
			h_PHIJet[NoFile][iy]->SetMinimum(0.1);
			h_PHIJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_PHIJet[NoFile][iy]->Integral()>0 ) h_PHIJet[NoFile][iy]->Scale(h_PHIJet[0][iy]->Integral()/h_PHIJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) 	{	h_PHIJet[NoFile][iy]->Draw("hist");	paveCMS ->Draw("same");  }
			else h_PHIJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");

			pad7->cd(iy+1);
			pad7->cd(iy+1)->SetLogy(1);
			h_CMJet[NoFile][iy]->GetXaxis()->SetTitle("Charged Multiplicity");
			h_CMJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_CMJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_CMJet[NoFile][iy]->SetLineColor(Colors[NoFile]); 
			h_CMJet[NoFile][iy]->SetLineStyle(1);
			h_CMJet[NoFile][iy]->SetMinimum(0.1);
			h_CMJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_CMJet[NoFile][iy]->Integral()>0 ) h_CMJet[NoFile][iy]->Scale(h_CMJet[0][iy]->Integral()/h_CMJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) 	{	h_CMJet[NoFile][iy]->Draw("hist");	paveCMS ->Draw("same");  }
			else h_CMJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");

			pad8->cd(iy+1);
			pad8->cd(iy+1)->SetLogy(1);
			h_NMJet[NoFile][iy]->GetXaxis()->SetTitle("Neutral Multiplicity");
			h_NMJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_NMJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NMJet[NoFile][iy]->SetLineColor(Colors[NoFile]); 
			h_NMJet[NoFile][iy]->SetLineStyle(1);
			h_NMJet[NoFile][iy]->SetMinimum(0.1);
			h_NMJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_NMJet[NoFile][iy]->Integral()>0 ) h_NMJet[NoFile][iy]->Scale(h_NMJet[0][iy]->Integral()/h_NMJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) 	{	h_NMJet[NoFile][iy]->Draw("hist");	paveCMS ->Draw("same");  }
			else h_NMJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");

			pad9->cd(iy+1);
			pad9->cd(iy+1)->SetLogy(1);
			h_ptJet[NoFile][iy]->GetXaxis()->SetTitle("Jet pT (GeV)");
			h_ptJet[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_ptJet[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_ptJet[NoFile][iy]->SetLineColor(Colors[NoFile]); 
			h_ptJet[NoFile][iy]->SetLineStyle(1);
			h_ptJet[NoFile][iy]->SetMinimum(0.1);
			h_ptJet[NoFile][iy]->SetMaximum(100000);
			if (scale_histos && h_ptJet[NoFile][iy]->Integral()>0 ) h_ptJet[NoFile][iy]->Scale(h_ptJet[0][iy]->Integral()/h_ptJet[NoFile][iy]->Integral());
			if ( NoFile==0 ) 	{	h_ptJet[NoFile][iy]->Draw("hist");	paveCMS ->Draw("same");  }
			else h_ptJet[NoFile][iy]->Draw("same");
			teta->Draw();
			//leg1->Draw("same");
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) leg1->Draw("same");


			double res_mean = h_pT_res[NoFile][iy]->GetMean();
			double res_rms = h_pT_res[NoFile][iy]->GetRMS();

			if( NoFile == 0)
			{
				leg3[iy] =new TLegend(.2, .75, .8, .85);//7899//4899
				leg3[iy]->SetTextSize(0.04);
				leg3[iy]->SetFillColor(0); 
				leg3[iy]->SetBorderSize(0);
	 		}


			
//			sprintf(res_text[NoFile], "Offset = %3.2f ,  RES = %3.1f %%",1.0-res_mean,100*res_rms );
			sprintf(res_text, "Offset = %3.2f ,  RES = %3.1f %%",res_mean-1.0,100*res_rms );
			leg3[iy]->AddEntry(h_pT_res[NoFile][iy], res_text , "L");


			pad10->cd(iy+1);
			pad10->cd(iy+1)->SetLogy(1);
			h_pT_res[NoFile][iy]->GetXaxis()->SetTitle("pT_reco / pT_gen");
			h_pT_res[NoFile][iy]->GetYaxis()->SetTitle("Entries");
			h_pT_res[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);

			h_pT_res[NoFile][iy]->SetLineColor(Colors[NoFile]); 
			h_pT_res[NoFile][iy]->SetLineStyle(1);
			h_pT_res[NoFile][iy]->SetMinimum(0.1);
			h_pT_res[NoFile][iy]->SetMaximum(100000);
			if ( NoFile==0 ) 	{	h_pT_res[NoFile][iy]->Draw("hist");	paveCMS ->Draw("same");  }
			else h_pT_res[NoFile][iy]->Draw("same");
			if ( NoFile == NoFiles-1 )	leg3[iy]->Draw("same");
			teta->Draw();



			pad11->cd(iy+1);
			if(NoFile == 0 )
			{	
				TH1F *hr = pad11->cd(iy+1)->DrawFrame(0,0.6,pThistoMax,1.2);
				//if(eta_bin_counter==5) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.8,3000,1.2);
				//if(eta_bin_counter==6) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.0,3000,1.5);
				hr->SetXTitle("p_{T} gen jet (GeV)");
				hr->SetYTitle("Reconstruction Efficiency");
				hr->GetYaxis()->SetTitleOffset(1.3);
				hr->GetXaxis()->SetTitleOffset(1.3);
				hr->SetTitle(eta_bins_legend[iy]);
			}
			Reco_eff[NoFile][iy]->SetMinimum(0.4);
			Reco_eff[NoFile][iy]->SetMaximum(1.2);
			Reco_eff[NoFile][iy]->SetMarkerStyle(24);
			Reco_eff[NoFile][iy]->SetMarkerColor(Colors[NoFile]);
			Reco_eff[NoFile][iy]->SetMarkerSize(0.3);
			Reco_eff[NoFile][iy]->SetLineColor(Colors[NoFile]);
			if ( NoFile==0 ) 	{	Reco_eff[NoFile][iy]->Draw("p");	paveCMS ->Draw("same");  }
			else Reco_eff[NoFile][iy]->Draw("same p");
			teta->Draw();



			pad12->cd(iy+1);
			//pad12->cd(iy+1)->SetLogy(1);
			if(NoFile == 0 )
			{			
				TH1F *hr12 = pad12->cd(iy+1)->DrawFrame(0,-0.1,pThistoMax,1.0);
			//if(eta_bin_counter==5) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.8,3000,1.2);
			//if(eta_bin_counter==6) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.0,3000,1.5);
				hr12->SetXTitle("reco jet p_{T} (GeV)");
				hr12->SetYTitle("Reconstruction Fake rate");
				hr12->GetYaxis()->SetTitleOffset(1.3);
				hr12->GetXaxis()->SetTitleOffset(1.3);
				hr12->SetTitle(eta_bins_legend[iy]);
			}
			Fake_rate[NoFile][iy]->SetMinimum(0.4);
			Fake_rate[NoFile][iy]->SetMaximum(1.2);
			Fake_rate[NoFile][iy]->SetMarkerStyle(24);
			Fake_rate[NoFile][iy]->SetMarkerColor(Colors[NoFile]);
			Fake_rate[NoFile][iy]->SetMarkerSize(0.3);
			Fake_rate[NoFile][iy]->SetLineColor(Colors[NoFile]);
			if ( NoFile==0 ) 	{	Fake_rate[NoFile][iy]->Draw("p");	paveCMS ->Draw("same");  }
			else Fake_rate[NoFile][iy]->Draw("same p");
			teta->Draw();



		}
	 

		if(eta_bins > 1 &&  NoFile == NoFiles-1 )
		{
			pad1->cd(eta_bins+1);	leg1->Draw();
			pad2->cd(eta_bins+1);	leg1->Draw();
			pad3->cd(eta_bins+1);	leg1->Draw();
			pad4->cd(eta_bins+1);	leg1->Draw();
			pad5->cd(eta_bins+1);	leg1->Draw();
			pad6->cd(eta_bins+1);	leg1->Draw();
			pad7->cd(eta_bins+1);	leg1->Draw();
			pad8->cd(eta_bins+1);	leg1->Draw();
			pad9->cd(eta_bins+1);	leg1->Draw();
			pad10->cd(eta_bins+1);	leg1->Draw();
			pad11->cd(eta_bins+1);	leg1->Draw();
			pad12->cd(eta_bins+1);	leg1->Draw();
		}

	} // end of loop on files

	if (Save_Plots)
	{ 
		sprintf(filename,"%s/%s/%s_chf.png",analyzer_path,output_directory,image_name);
		pad1->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_nhf.png",analyzer_path,output_directory,image_name);
		pad2->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_cemf.png",analyzer_path,output_directory,image_name);
		pad3->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_nemf.png",analyzer_path,output_directory,image_name);
		pad4->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_muf.png",analyzer_path,output_directory,image_name);
		pad5->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_phi.png",analyzer_path,output_directory,image_name);
		pad6->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_CM.png",analyzer_path,output_directory,image_name);
		pad7->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NM.png",analyzer_path,output_directory,image_name);
		pad8->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_pt.png",analyzer_path,output_directory,image_name);
		pad9->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_pT_resolution.png",analyzer_path,output_directory,image_name);
		pad10->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reco_efficiency.png",analyzer_path,output_directory,image_name);
		pad11->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reco_Fake_rate.png",analyzer_path,output_directory,image_name);
		pad12->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_eta.png",analyzer_path,output_directory,image_name);
		c_eta->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_SumEt.png",analyzer_path,output_directory,image_name);
		c_SumEt->SaveAs(filename);

   }
}

