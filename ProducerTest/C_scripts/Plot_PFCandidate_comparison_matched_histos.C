#include "include/include_functions.h"
#include "include/Input_definition.h"



void Plot_PFCandidate_comparison_matched_histos()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
//	gStyle->SetOptStat(0);

	double pTmax = 1000;




	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);
	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;




	TFile *f[NoFiles];
	TTree *tree[NoFiles];

	char eta_bins_legend[eta_bins][25];
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

	char name[1500]; 
	sprintf(name, "%s/%s/PF_candidates_matched_histos_%s.root",analyzer_path,output_directory,image_name );
	TFile *f_input = new TFile (name,"READ");

	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int iy=0; iy< eta_bins; iy++)
		{
	
			sprintf(name,"h_ChargedHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_pT[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NeutralHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_pT[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_ChargedHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_SumpT[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NeutralHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_SumpT[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_ChargedHadron_AveragepT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_AveragepT[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NeutralHadron_AveragepT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_AveragepT[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NumberOfChargedHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfChargedHadrons[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NumberOfNeutralHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfNeutralHadrons[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_ChargedHadron_DR_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_DR[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_NeutralHadron_DR_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_DR[NoFile][iy] = (TH1D*)(f_input->Get(name));

			sprintf(name,"h_ChargedHadron_eta_phi_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_eta_phi[NoFile][iy] = (TH2D*)(f_input->Get(name));

			sprintf(name,"h_NeutralHadron_eta_phi_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_eta_phi[NoFile][iy] = (TH2D*)(f_input->Get(name));

			sprintf(name,"h_girth_%s_bin%i",legend_array[NoFile],iy);
			h_girth[NoFile][iy] = (TH1D*)(f_input->Get(name)); 

			sprintf(name,"h_quadratic_moment_%s_bin%i",legend_array[NoFile],iy);
			h_quadratic_moment[NoFile][iy] = (TH1D*)(f_input->Get(name));
		}
	}



// =========================================================== Plotting staff ===================================================================

	TCanvas *pad_ChargedHadron_SumpT = new TCanvas("pad_ChargedHadron_SumpT", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_ChargedHadron_SumpT->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_NeutralHadron_SumpT = new TCanvas("pad_NeutralHadron_SumpT", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_NeutralHadron_SumpT->Divide(PadColumnsEtaBins, PadRowsEtaBins);

	TCanvas *pad_ChargedHadron_AveragepT = new TCanvas("pad_ChargedHadron_AveragepT", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_ChargedHadron_AveragepT->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	
	TCanvas *pad_NeutralHadron_AveragepT = new TCanvas("pad_NeutralHadron_AveragepT", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_NeutralHadron_AveragepT->Divide(PadColumnsEtaBins, PadRowsEtaBins);

	TCanvas *pad_ChargedHadron_pT = new TCanvas("pad_ChargedHadron_pT", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_ChargedHadron_pT->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_NeutralHadron_pT = new TCanvas("pad_NeutralHadron_pT", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_NeutralHadron_pT->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_NumberOfChargedHadrons = new TCanvas("pad_NumberOfChargedHadrons", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_NumberOfChargedHadrons->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_NumberOfNeutralHadrons = new TCanvas("pad_NumberOfNeutralHadrons", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_NumberOfNeutralHadrons->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_ChargedCand_DR = new TCanvas("pad_ChargedCand_DR", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_ChargedCand_DR->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_NeutralCand_DR = new TCanvas("pad_NeutralCand_DR", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_NeutralCand_DR->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_girth = new TCanvas("pad_girth", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_girth->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad_quadr_moment = new TCanvas("pad_quadr_moment", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_quadr_moment->Divide(PadColumnsEtaBins, PadRowsEtaBins);


	TCanvas *pad_ChargedHadron_eta_phi[NoFiles], *pad_NeutralHadron_eta_phi[NoFiles];



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

	TH1D *frameChargedHadron_SumpT[eta_bins], *frameNeutralHadron_SumpT[eta_bins], *frameChargedHadron_AveragepT[eta_bins], *frameNeutralHadron_AveragepT[eta_bins], *frameChargedHadron_pT[eta_bins], *frameNeutralHadron_pT[eta_bins], *frameNumberOfChargedHadrons[eta_bins], *frameNumberOfNeutralHadrons[eta_bins], *frameChargedHadron_DR[eta_bins], *frameNeutralHadron_DR[eta_bins], *framegirth[eta_bins], *framequadratic_moment[eta_bins];



	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		sprintf(name,"pad_ChargedHadron_eta_phi_%s",legend_array[NoFile]);
		pad_ChargedHadron_eta_phi[NoFile] = new TCanvas(name, name , Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
		pad_ChargedHadron_eta_phi[NoFile]->Divide(PadColumnsEtaBins, PadRowsEtaBins);

		sprintf(name,"pad_NeutralHadron_eta_phi_%s",legend_array[NoFile]);
		pad_NeutralHadron_eta_phi[NoFile] = new TCanvas(name, name , Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
		pad_NeutralHadron_eta_phi[NoFile]->Divide(PadColumnsEtaBins, PadRowsEtaBins);


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

			
			if (NoFile==0 )  frameChargedHadron_SumpT[iy] = InitiateFrameOnCanvasPad(pad_ChargedHadron_SumpT, iy+1, "frameChargedHadron_SumpT", "Sum pT of Charged hadrons per jet (GeV)", "#jets", 0., pTmax, 10*YaxisLowEndMultiplier*h_ChargedHadron_SumpT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ChargedHadron_SumpT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ChargedHadron_SumpT[NoFile][iy]->Integral()>0 ) h_ChargedHadron_SumpT[NoFile][iy]->Scale(h_ChargedHadron_SumpT[0][iy]->Integral()/h_ChargedHadron_SumpT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_ChargedHadron_SumpT, iy+1, h_ChargedHadron_SumpT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad_ChargedHadron_SumpT->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_ChargedHadron_SumpT->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNeutralHadron_SumpT[iy] = InitiateFrameOnCanvasPad(pad_NeutralHadron_SumpT, iy+1, "frameNeutralHadron_SumpT", "Sum pT of Neutral hadrons per jet (GeV)", "#jets", 0., pTmax, 10*YaxisLowEndMultiplier*h_NeutralHadron_SumpT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NeutralHadron_SumpT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NeutralHadron_SumpT[NoFile][iy]->Integral()>0 ) h_NeutralHadron_SumpT[NoFile][iy]->Scale(h_NeutralHadron_SumpT[0][iy]->Integral()/h_NeutralHadron_SumpT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_NeutralHadron_SumpT, iy+1, h_NeutralHadron_SumpT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad_NeutralHadron_SumpT->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_NeutralHadron_SumpT->cd(iy+1); teta->Draw(); }



			if (NoFile==0 )  frameChargedHadron_AveragepT[iy] = InitiateFrameOnCanvasPad(pad_ChargedHadron_AveragepT, iy+1, "frameChargedHadron_AveragepT", "Average pT of Charged hadrons per jet (GeV)", "#jets", 0., 80, 0.1*YaxisLowEndMultiplier*h_ChargedHadron_AveragepT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ChargedHadron_AveragepT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ChargedHadron_AveragepT[NoFile][iy]->Integral()>0 ) h_ChargedHadron_AveragepT[NoFile][iy]->Scale(h_ChargedHadron_AveragepT[0][iy]->Integral()/h_ChargedHadron_AveragepT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_ChargedHadron_AveragepT, iy+1, h_ChargedHadron_AveragepT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad_ChargedHadron_AveragepT->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_ChargedHadron_AveragepT->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNeutralHadron_AveragepT[iy] = InitiateFrameOnCanvasPad(pad_NeutralHadron_AveragepT, iy+1, "frameNeutralHadron_AveragepT", "Average pT of Neutral hadrons per jet (GeV)", "#jets", 0., 300, 0.1*YaxisLowEndMultiplier*h_NeutralHadron_AveragepT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NeutralHadron_AveragepT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NeutralHadron_AveragepT[NoFile][iy]->Integral()>0 ) h_NeutralHadron_AveragepT[NoFile][iy]->Scale(h_NeutralHadron_AveragepT[0][iy]->Integral()/h_NeutralHadron_AveragepT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_NeutralHadron_AveragepT, iy+1, h_NeutralHadron_AveragepT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad_NeutralHadron_AveragepT->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_NeutralHadron_AveragepT->cd(iy+1); teta->Draw(); }



			if (NoFile==0 )  frameChargedHadron_pT[iy] = InitiateFrameOnCanvasPad(pad_ChargedHadron_pT, iy+1, "frameChargedHadron_pT", "Charged Hadron pT (GeV)", "# Charged Hadrons", 0., pTmax, YaxisLowEndMultiplier*h_ChargedHadron_pT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ChargedHadron_pT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ChargedHadron_pT[NoFile][iy]->Integral()>0 ) h_ChargedHadron_pT[NoFile][iy]->Scale(h_ChargedHadron_pT[0][iy]->Integral()/h_ChargedHadron_pT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_ChargedHadron_pT, iy+1, h_ChargedHadron_pT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad_ChargedHadron_pT->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_ChargedHadron_pT->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNeutralHadron_pT[iy] = InitiateFrameOnCanvasPad(pad_NeutralHadron_pT, iy+1, "frameNeutralHadron_pT", "Neutral Hadron pT (GeV)", "# Neutral Hadrons", 0., pTmax, YaxisLowEndMultiplier*h_NeutralHadron_pT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NeutralHadron_pT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NeutralHadron_pT[NoFile][iy]->Integral()>0 ) h_NeutralHadron_pT[NoFile][iy]->Scale(h_NeutralHadron_pT[0][iy]->Integral()/h_NeutralHadron_pT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_NeutralHadron_pT, iy+1, h_NeutralHadron_pT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_NeutralHadron_pT->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_NeutralHadron_pT->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNumberOfChargedHadrons[iy] = InitiateFrameOnCanvasPad(pad_NumberOfChargedHadrons, iy+1, "frameNumberOfChargedHadrons", "# of charged candidates per jet", "# jets", 0., 50., YaxisLowEndMultiplier*h_NumberOfChargedHadrons[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NumberOfChargedHadrons[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NumberOfChargedHadrons[NoFile][iy]->Integral()>0 ) h_NumberOfChargedHadrons[NoFile][iy]->Scale(h_NumberOfChargedHadrons[0][iy]->Integral()/h_NumberOfChargedHadrons[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_NumberOfChargedHadrons, iy+1, h_NumberOfChargedHadrons[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_NumberOfChargedHadrons->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_NumberOfChargedHadrons->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNumberOfNeutralHadrons[iy] = InitiateFrameOnCanvasPad(pad_NumberOfNeutralHadrons, iy+1, "frameNumberOfNeutralHadrons", "# of Neutral candidates per jet", "# jets", 0., 50., YaxisLowEndMultiplier*h_NumberOfNeutralHadrons[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NumberOfNeutralHadrons[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NumberOfNeutralHadrons[NoFile][iy]->Integral()>0 ) h_NumberOfNeutralHadrons[NoFile][iy]->Scale(h_NumberOfNeutralHadrons[0][iy]->Integral()/h_NumberOfNeutralHadrons[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_NumberOfNeutralHadrons, iy+1, h_NumberOfNeutralHadrons[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_NumberOfNeutralHadrons->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_NumberOfNeutralHadrons->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameChargedHadron_DR[iy] = InitiateFrameOnCanvasPad(pad_ChargedCand_DR, iy+1, "frameChargedHadron_DR", "Jet-ChargedCandidate DR", "# jets", 0., 0.9, YaxisLowEndMultiplier*h_ChargedHadron_DR[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ChargedHadron_DR[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ChargedHadron_DR[NoFile][iy]->Integral()>0 ) h_ChargedHadron_DR[NoFile][iy]->Scale(h_ChargedHadron_DR[0][iy]->Integral()/h_ChargedHadron_DR[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_ChargedCand_DR, iy+1, h_ChargedHadron_DR[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_ChargedCand_DR->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_ChargedCand_DR->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNeutralHadron_DR[iy] = InitiateFrameOnCanvasPad(pad_NeutralCand_DR, iy+1, "frameNeutralHadron_DR", "Jet-NeutralCandidate DR", "# jets", 0., 0.9, YaxisLowEndMultiplier*h_NeutralHadron_DR[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NeutralHadron_DR[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NeutralHadron_DR[NoFile][iy]->Integral()>0 ) h_NeutralHadron_DR[NoFile][iy]->Scale(h_NeutralHadron_DR[0][iy]->Integral()/h_NeutralHadron_DR[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_NeutralCand_DR, iy+1, h_NeutralHadron_DR[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_NeutralCand_DR->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_NeutralCand_DR->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  framegirth[iy] = InitiateFrameOnCanvasPad(pad_girth, iy+1, "framegirth", "Jet girth", "# jets", 0., 1.0, YaxisLowEndMultiplier*h_girth[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_girth[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_girth[NoFile][iy]->Integral()>0 ) h_girth[NoFile][iy]->Scale(h_girth[0][iy]->Integral()/h_girth[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_girth, iy+1, h_girth[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_girth->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_girth->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  framequadratic_moment[iy] = InitiateFrameOnCanvasPad(pad_quadr_moment, iy+1, "framequadratic_moment", "Quadratic Radial Moment", "# jets", 0., 1.0, YaxisLowEndMultiplier*h_quadratic_moment[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_quadratic_moment[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_quadratic_moment[NoFile][iy]->Integral()>0 ) h_quadratic_moment[NoFile][iy]->Scale(h_quadratic_moment[0][iy]->Integral()/h_quadratic_moment[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_quadr_moment, iy+1, h_quadratic_moment[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_quadr_moment->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_quadr_moment->cd(iy+1); teta->Draw(); }




//============================================== 2D plots ==============================================

			if(eta_bins > 1) pad_ChargedHadron_eta_phi[NoFile]->cd(iy+1);
			else pad_ChargedHadron_eta_phi[NoFile]->cd();
			if(eta_bins > 1)
			{
				pad_ChargedHadron_eta_phi[NoFile]->cd(iy+1)->SetLogz(1);
				pad_ChargedHadron_eta_phi[NoFile]->cd(iy+1)->SetRightMargin(0.18);
				pad_ChargedHadron_eta_phi[NoFile]->cd(iy+1)->SetTopMargin(0.07);
			}
			else
			{
				pad_ChargedHadron_eta_phi[NoFile]->cd()->SetLogz(1);
				pad_ChargedHadron_eta_phi[NoFile]->cd()->SetRightMargin(0.18);
				pad_ChargedHadron_eta_phi[NoFile]->cd()->SetTopMargin(0.07);
			}
			
			h_ChargedHadron_eta_phi[NoFile][iy]->SetStats(true);
		    

			h_ChargedHadron_eta_phi[NoFile][iy]->SetTitle("ChargedCandidate eta-phi");
			h_ChargedHadron_eta_phi[NoFile][iy]->GetXaxis()->SetTitle("#eta_{cand} - #eta_{jet}");
			h_ChargedHadron_eta_phi[NoFile][iy]->GetYaxis()->SetTitle("#phi_{cand} - #phi_{jet}");
			h_ChargedHadron_eta_phi[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_ChargedHadron_eta_phi[NoFile][iy]->SetMinimum(0.001*h_ChargedHadron_eta_phi[NoFile][iy]->GetMaximum() );
			h_ChargedHadron_eta_phi[NoFile][iy]->SetMaximum(3*h_ChargedHadron_eta_phi[NoFile][iy]->GetMaximum());

			h_ChargedHadron_eta_phi[NoFile][iy]->Draw("colz");
			paveCMS ->Draw("same");
			teta->Draw();



			if(eta_bins > 1) pad_NeutralHadron_eta_phi[NoFile]->cd(iy+1);
			else pad_NeutralHadron_eta_phi[NoFile]->cd();
			if(eta_bins > 1)
			{
				pad_NeutralHadron_eta_phi[NoFile]->cd(iy+1)->SetLogz(1);
				pad_NeutralHadron_eta_phi[NoFile]->cd(iy+1)->SetRightMargin(0.18);
				pad_NeutralHadron_eta_phi[NoFile]->cd(iy+1)->SetTopMargin(0.07);
			}
			else
			{
				pad_NeutralHadron_eta_phi[NoFile]->cd()->SetLogz(1);
				pad_NeutralHadron_eta_phi[NoFile]->cd()->SetRightMargin(0.18);
				pad_NeutralHadron_eta_phi[NoFile]->cd()->SetTopMargin(0.07);
			}
			
			h_NeutralHadron_eta_phi[NoFile][iy]->SetStats(true);
		    

			h_NeutralHadron_eta_phi[NoFile][iy]->SetTitle("NeutralCandidate eta-phi");
			h_NeutralHadron_eta_phi[NoFile][iy]->GetXaxis()->SetTitle("#eta_{cand} - #eta_{jet}");
			h_NeutralHadron_eta_phi[NoFile][iy]->GetYaxis()->SetTitle("#phi_{cand} - #phi_{jet}");
			h_NeutralHadron_eta_phi[NoFile][iy]->GetYaxis()->SetTitleOffset(1.3);
			h_NeutralHadron_eta_phi[NoFile][iy]->SetMinimum(0.001*h_NeutralHadron_eta_phi[NoFile][iy]->GetMaximum() );
			h_NeutralHadron_eta_phi[NoFile][iy]->SetMaximum(3*h_NeutralHadron_eta_phi[NoFile][iy]->GetMaximum());

			h_NeutralHadron_eta_phi[NoFile][iy]->Draw("colz");
			paveCMS ->Draw("same");
			teta->Draw();




			if(eta_bins > 1 &&  NoFile == NoFiles-1 )
			{
				pad_ChargedHadron_SumpT->cd(eta_bins+1);        leg1->Draw();
				pad_NeutralHadron_SumpT->cd(eta_bins+1);        leg1->Draw();
				pad_ChargedHadron_AveragepT->cd(eta_bins+1);    leg1->Draw();
				pad_NeutralHadron_AveragepT->cd(eta_bins+1);    leg1->Draw();
				pad_ChargedHadron_pT->cd(eta_bins+1);           leg1->Draw();
				pad_NeutralHadron_pT->cd(eta_bins+1);           leg1->Draw();
				pad_NumberOfChargedHadrons->cd(eta_bins+1);     leg1->Draw();
				pad_NumberOfNeutralHadrons->cd(eta_bins+1);     leg1->Draw();
				pad_girth->cd(eta_bins+1);   		        leg1->Draw();
				pad_quadr_moment->cd(eta_bins+1);	        leg1->Draw();

				pad_ChargedCand_DR->cd(eta_bins+1);	        leg1->Draw();
				pad_NeutralCand_DR->cd(eta_bins+1);	        leg1->Draw();
			}

		} //end of loop on eta bins
	} // end of loop on files



	char filename[1500]; 
	if(Save_Plots)
	{ 
		sprintf(filename,"%s/%s/%s_ChargedPFcand_SumpT.png",analyzer_path,output_directory,image_name);
		pad_ChargedHadron_SumpT->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_SumpT.png",analyzer_path,output_directory,image_name);
		pad_NeutralHadron_SumpT->SaveAs(filename);


		sprintf(filename,"%s/%s/%s_ChargedPFcand_AveragepT.png",analyzer_path,output_directory,image_name);
		pad_ChargedHadron_AveragepT->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_AveragepT.png",analyzer_path,output_directory,image_name);
		pad_NeutralHadron_AveragepT->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_ChargedPFcand_pT.png",analyzer_path,output_directory,image_name);
		pad_ChargedHadron_pT->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_pT.png",analyzer_path,output_directory,image_name);
		pad_NeutralHadron_pT->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_ChargedPFcand_number.png",analyzer_path,output_directory,image_name);
		pad_NumberOfChargedHadrons->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_number.png",analyzer_path,output_directory,image_name);
		pad_NumberOfNeutralHadrons->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_ChargedPFcand_DR.png",analyzer_path,output_directory,image_name);
		pad_ChargedCand_DR->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_DR.png",analyzer_path,output_directory,image_name);
		pad_NeutralCand_DR->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_jet_girth.png",analyzer_path,output_directory,image_name);
		pad_girth->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_jet_quadratic_moment.png",analyzer_path,output_directory,image_name);
		pad_quadr_moment->SaveAs(filename);

		for (int NoFile=0; NoFile<NoFiles; NoFile++)
		{
			sprintf(filename,"%s/%s/%s_ChargedPFcand_eta_phi_%s.png",analyzer_path,output_directory,image_name,legend_array[NoFile]);
			pad_ChargedHadron_eta_phi[NoFile]->SaveAs(filename);
		
			sprintf(filename,"%s/%s/%s_NeutralPFcand_eta_phi_%s.png",analyzer_path,output_directory,image_name,legend_array[NoFile]);
			pad_NeutralHadron_eta_phi[NoFile]->SaveAs(filename);
		}

	}

}

