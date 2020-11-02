#include "include/include_functions.h"
#include "include/Input_definition.h"

void Plot_Tracking_Comparisons_matched_histos()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);

	double pThistoMax = 1000;


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

	char name[1500]; 
	sprintf(name, "%s/%s/Jet_characteristics_histos_%s.root",analyzer_path,output_directory,image_name );
	TFile *f_input = new TFile (name,"READ");

	f_input->ls();
	char filename[1024]; 


	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko();


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		sprintf(name,"h_ETAJet_%s",legend_array[NoFile]);
		h_ETAJet[NoFile] = (TH1D*)(f_input->Get(name));

	
		sprintf(name,"h_METovSumEt_%s",legend_array[NoFile]);
		h_METovSumEt[NoFile] = (TH1D*)(f_input->Get(name));

		sprintf(name,"h_MET_%s",legend_array[NoFile]);
		h_MET[NoFile] = (TH1D*)(f_input->Get(name)); // 40,0,1.0
	
		sprintf(name,"h_minDR_%s",legend_array[NoFile]);
		h_minDR[NoFile] = (TH1D*)(f_input->Get(name));

		for(Int_t h=0; h<eta_bins;h++)
		{ 
		//========== reco jets===========
			

			sprintf(name,"h_CHFJet_%s_bin%i",legend_array[NoFile],h);
			h_CHFJet[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NHFJet_%s_bin%i",legend_array[NoFile],h);
			h_NHFJet[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_CEMFJet_%s_bin%i",legend_array[NoFile],h);
			h_CEMFJet[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NEMFJet_%s_bin%i",legend_array[NoFile],h);
			h_NEMFJet[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_MUFJet_%s_bin%i",legend_array[NoFile],h);
			h_MUFJet[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_CMJet_%s_bin%i",legend_array[NoFile],h);
			h_CMJet[NoFile][h] = (TH1D*)(f_input->Get(name));
		
			sprintf(name,"h_PHIJet_%s_bin%i",legend_array[NoFile],h);
			h_PHIJet[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NMJet_%s_bin%i",legend_array[NoFile],h);
			h_NMJet[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_ptJet_%s_bin%i",legend_array[NoFile],h);
			h_ptJet[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_gen_pT_all_%s_bin%i",legend_array[NoFile],h);
			h_gen_pT_all[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_gen_pT_%s_bin%i",legend_array[NoFile],h);
			h_gen_pT[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_reco_pT_unmatched_%s_bin%i",legend_array[NoFile],h);
			h_reco_pT_unmatched[NoFile][h] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_reco_pT_all_%s_bin%i",legend_array[NoFile],h);
			h_reco_pT_all[NoFile][h] = (TH1D*)(f_input->Get(name));


		}// end of etabin loop
	} // end of file loop


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for(Int_t h=0; h<eta_bins;h++)
		{
				Reco_eff[NoFile][h] = GetEfficiencyGraph(h_gen_pT[NoFile][h] , h_gen_pT_all[NoFile][h] );
				Fake_rate[NoFile][h] = GetEfficiencyGraph(h_reco_pT_unmatched[NoFile][h],h_reco_pT_all[NoFile][h]  ) ; 
		}
	}

//========================================== plotting stuff ==========================

	TCanvas *pad_chf = new TCanvas("pad_chf", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_chf->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_nhf = new TCanvas("pad_nhf", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_nhf->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_cemf = new TCanvas("pad_cemf", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_cemf->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_nemf = new TCanvas("pad_nemf", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_nemf->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_muf = new TCanvas("pad_muf", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_muf->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_phi = new TCanvas("pad_phi", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_phi->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_CM = new TCanvas("pad_CM", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_CM->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_NM = new TCanvas("pad_NM", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_NM->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_pt = new TCanvas("pad_pt", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_pt->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad11 = new TCanvas("pad11", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad11->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad12 = new TCanvas("pad12", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad12->Divide(PadColumnsEtaBins,PadRowsEtaBins);

	TCanvas *c_eta = new TCanvas("c_eta", "",1);
	TCanvas *c_minDR = new TCanvas("c_minDR", "",1);
	TCanvas *c_METovSumEt = new TCanvas("c_METovSumEt", "",1);
	TCanvas *c_MET = new TCanvas("c_MET", "",1);


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
	TH1D *frameCHF[eta_bins], *frameNHF[eta_bins], *frameCEMF[eta_bins], *frameNEMF[eta_bins], *frameMUF[eta_bins], *frameCM[eta_bins], *frameNM[eta_bins], *framePhi[eta_bins], *framePt[eta_bins],*frameRecoEff[eta_bins] , *frameFakeRate[eta_bins], *frameEta, *frameMETovSumEt, *frameMET, *frameMinDR;


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		leg1->AddEntry(h_CHFJet[NoFile][0], legend_array[NoFile], "L");
		leg2->AddEntry(h_CHFJet[NoFile][0], legend_array[NoFile], "L");



		if (NoFile==0 )  frameEta = InitiateFrameOnCanvasPad(c_eta, 0 , "frameEta", "Jet eta", "Entries", -5., 5., YaxisLowEndMultiplier*h_ETAJet[NoFile]->GetMaximum(), YaxisHighEndMultiplier*h_ETAJet[NoFile]->GetMaximum(), true, paveCMS);
		if (scale_histos && h_ETAJet[NoFile]->Integral()>0 ) h_ETAJet[NoFile]->Scale(h_ETAJet[0]->Integral() / h_ETAJet[NoFile]->Integral());
		DrawHistoToCanvasPad(c_eta, 0, h_ETAJet[NoFile], Colors[NoFile], 1);
		if( NoFile== 0)  leg2->Draw("same"); 


		if (NoFile==0 )  frameMETovSumEt = InitiateFrameOnCanvasPad(c_METovSumEt, 0 , "frameMETovSumEt", "MET / SumEt", "Entries", 0., 1.2, YaxisLowEndMultiplier*h_METovSumEt[NoFile]->GetMaximum(), YaxisHighEndMultiplier*h_METovSumEt[NoFile]->GetMaximum(), true, paveCMS);
		if (scale_histos && h_METovSumEt[NoFile]->Integral()>0 ) { h_METovSumEt[NoFile]->Scale(h_METovSumEt[0]->Integral() / h_METovSumEt[NoFile]->Integral()); h_METovSumEt[NoFile]->Rebin(2); }
		DrawHistoToCanvasPad(c_METovSumEt, 0, h_METovSumEt[NoFile], Colors[NoFile], 1);
		if( NoFile== 0)  leg2->Draw("same"); 

		if (NoFile==0 )  frameMET = InitiateFrameOnCanvasPad(c_MET, 0 , "frameMET", "MET (GeV)", "Entries", 0., 600, 0.0001*YaxisLowEndMultiplier*h_MET[NoFile]->GetMaximum(), YaxisHighEndMultiplier*h_MET[NoFile]->GetMaximum(), true, paveCMS);
		if (scale_histos && h_MET[NoFile]->Integral()>0 ) h_MET[NoFile]->Scale(h_MET[0]->Integral() / h_MET[NoFile]->Integral());
		DrawHistoToCanvasPad(c_MET, 0, h_MET[NoFile], Colors[NoFile], 1);
		if( NoFile== 0)  leg2->Draw("same"); 

		if(plot_minDR)
		{
			if (NoFile==0 )  frameMinDR = InitiateFrameOnCanvasPad(c_minDR, 0 , "frameMinDR", "minDR", "Entries", 0., 5., YaxisLowEndMultiplier*h_minDR[NoFile]->GetMaximum(), YaxisHighEndMultiplier*h_minDR[NoFile]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_minDR[NoFile]->Integral()>0 ) h_minDR[NoFile]->Scale(h_minDR[0]->Integral() /h_minDR[NoFile]->Integral());
			DrawHistoToCanvasPad(c_minDR, 0,  h_minDR[NoFile], Colors[NoFile], 1);
			if( NoFile== 0)  leg2->Draw("same"); 
		}

		for(int iy=0; iy<eta_bins; iy++)
		{
			

			double etamin = yBnd[iy];
			double etamax = yBnd[iy+1];
			const char *seta = (etamin==0 ? Form("|y| < %1.2g",etamax) :
			Form("%1.2g #leq |y| < %1.2g",etamin,etamax));
			TLatex *teta = new TLatex(0.38,0.86,seta); //cout<<seta<<endl;
			teta->SetNDC();
			teta->SetTextSize(0.06);


			if (NoFile==0 )  frameCHF[iy] = InitiateFrameOnCanvasPad(pad_chf, iy+1, "frameCHF", "Charged Hadron Fraction", "Entries", 0., 1.1, YaxisLowEndMultiplier*h_CHFJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_CHFJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_CHFJet[NoFile][iy]->Integral()>0 ) h_CHFJet[NoFile][iy]->Scale(h_CHFJet[0][iy]->Integral()/h_CHFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_chf, iy+1, h_CHFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad_chf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_chf->cd(iy+1); teta->Draw(); }

			
			if (NoFile==0 )  frameNHF[iy] = InitiateFrameOnCanvasPad(pad_nhf,iy+1, "frameNHF", "Neutral Hadron Fraction", "Entries", 0., 1.1,YaxisLowEndMultiplier*h_NHFJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NHFJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_NHFJet[NoFile][iy]->Integral()>0 ) h_NHFJet[NoFile][iy]->Scale(h_NHFJet[0][iy]->Integral()/h_NHFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_nhf, iy+1, h_NHFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_nhf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_nhf->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameCEMF[iy] = InitiateFrameOnCanvasPad(pad_cemf,iy+1, "frameCEMF", "Charged E/M Fraction", "Entries", 0., 1.1, YaxisLowEndMultiplier*h_CEMFJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_CEMFJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_CEMFJet[NoFile][iy]->Integral()>0 ) h_CEMFJet[NoFile][iy]->Scale(h_CEMFJet[0][iy]->Integral()/h_CEMFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_cemf, iy+1, h_CEMFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_cemf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_cemf->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNEMF[iy] = InitiateFrameOnCanvasPad(pad_nemf,iy+1, "frameNEMF", "Neutral E/M Fraction", "Entries", 0., 1.1, YaxisLowEndMultiplier*h_NEMFJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NEMFJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_NEMFJet[NoFile][iy]->Integral()>0 ) h_NEMFJet[NoFile][iy]->Scale(h_NEMFJet[0][iy]->Integral()/h_NEMFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_nemf, iy+1, h_NEMFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_nemf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_nemf->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameMUF[iy] = InitiateFrameOnCanvasPad(pad_muf,iy+1, "frameMUF", "Muon Fraction", "Entries", 0., 1.1, YaxisLowEndMultiplier*h_MUFJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_MUFJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_MUFJet[NoFile][iy]->Integral()>0 ) h_MUFJet[NoFile][iy]->Scale(h_MUFJet[0][iy]->Integral()/h_MUFJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_muf, iy+1, h_MUFJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_muf->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_muf->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  framePhi[iy] = InitiateFrameOnCanvasPad(pad_phi,iy+1, "framePhi", "Jet Phi", "Entries", -3.14, 3.14,YaxisLowEndMultiplier*h_PHIJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_PHIJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_PHIJet[NoFile][iy]->Integral()>0 ) h_PHIJet[NoFile][iy]->Scale(h_PHIJet[0][iy]->Integral()/h_PHIJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_phi, iy+1, h_PHIJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_phi->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_phi->cd(iy+1); teta->Draw(); }


			if (NoFile==0 ) frameCM[iy] = InitiateFrameOnCanvasPad(pad_CM,iy+1, "frameCM", "Charged Multiplicity", "Entries", 0., 80., YaxisLowEndMultiplier*h_CMJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_CMJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_CMJet[NoFile][iy]->Integral()>0 ) h_CMJet[NoFile][iy]->Scale(h_CMJet[0][iy]->Integral()/h_CMJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_CM, iy+1, h_CMJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_CM->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_CM->cd(iy+1); teta->Draw(); }


			if (NoFile==0 ) frameNM[iy] = InitiateFrameOnCanvasPad(pad_NM,iy+1, "frameNM", "Neutral Multiplicity", "Entries", 0., 80., YaxisLowEndMultiplier*h_NMJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NMJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_NMJet[NoFile][iy]->Integral()>0 ) h_NMJet[NoFile][iy]->Scale(h_NMJet[0][iy]->Integral()/h_NMJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_NM, iy+1, h_NMJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_NM->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_NM->cd(iy+1); teta->Draw(); }


			if (NoFile==0 ) framePt[iy] = InitiateFrameOnCanvasPad(pad_pt,iy+1, "framepT", "Jet pT (GeV)", "Entries", 0., pThistoMax,0.01*YaxisLowEndMultiplier*h_ptJet[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ptJet[NoFile][iy]->GetMaximum(), true, paveCMS);
			if (scale_histos && h_ptJet[NoFile][iy]->Integral()>0 ) h_ptJet[NoFile][iy]->Scale(h_ptJet[0][iy]->Integral()/h_ptJet[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_pt, iy+1, h_ptJet[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_pt->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_pt->cd(iy+1); teta->Draw(); }


			if (NoFile==0 ) frameRecoEff[iy] = InitiateFrameOnCanvasPad(pad11,iy+1, "frameRecoEff", "p_{T} gen jet (GeV)", "Reconstruction Efficiency", 0., pThistoMax, 0.6, 1.2, false, paveCMS);
			DrawGraphToCanvasPad(pad11, iy+1, Reco_eff[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad11->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad11->cd(iy+1); teta->Draw(); }


			if (NoFile==0 ) frameFakeRate[iy] = InitiateFrameOnCanvasPad(pad12,iy+1,"frameFakeRate","p_{T} reco jet (GeV)", "Reconstruction Fake rate", 0., pThistoMax, -0.1, 1.0, false, paveCMS);
			DrawGraphToCanvasPad(pad12, iy+1, Fake_rate[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad12->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad12->cd(iy+1); teta->Draw(); }


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
			pad11->cd(eta_bins+1);	leg1->Draw();
			pad12->cd(eta_bins+1);	leg1->Draw();
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


		sprintf(filename,"%s/%s/%s_Reco_efficiency.png",analyzer_path,output_directory,image_name);
		pad11->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reco_Fake_rate.png",analyzer_path,output_directory,image_name);
		pad12->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_eta.png",analyzer_path,output_directory,image_name);
		c_eta->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_METovSumEt.png",analyzer_path,output_directory,image_name);
		c_METovSumEt->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_MET.png",analyzer_path,output_directory,image_name);
		c_MET->SaveAs(filename);

		if (plot_minDR) 
		{
			sprintf(filename,"%s/%s/%s_minDR.png",analyzer_path,output_directory,image_name);
			c_minDR->SaveAs(filename);
		}

   }
}

