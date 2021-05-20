#include "include/include_functions.h"
//#include "include/Input_definition.h"

void Plot_JetSubstructure_QCD_Z()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);

	double yBnd[]={0.0, 1.3, 2.4, 2.7, 3.0};  
	const int NoFiles = 4;
	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;
	bool plot_minDR =false;
	bool scale_histos = true ; //scale all histos to the number of entries of the histo of the 1st input
	bool Save_Plots = true; 
	
	char analyzer_path[500] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/"}; 
	char output_directory[200] = {"Zprime_plots/AK8_matched_to_genPtcles/M50/"};
	char image_name[200] = {"linear_scale"}; 


	TFile *f_input[NoFiles];
	
	f_input[0]=new TFile("/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/Zprime_plots/AK8_matched_to_genPtcles/M50/JetSubstructure_matched_histos_test.root","READ");

	f_input[1]=new TFile("/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/Zprime_plots/AK8_matched_to_genPtcles/M50/JetSubstructure_matched_histos_test.root","READ");

	f_input[2]=new TFile("/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/Zprime_plots/AK8_matched_to_genPtcles/M50/JetSubstructure_matched_histos_test.root","READ");

/*		
	f_input[0]=new TFile("/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/Zprime_plots/AK8_matched_to_genPtcles/qqDR_lessThan_0p8/M50/JetSubstructure_matched_histos_test.root","READ");

	f_input[1]=new TFile("/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/Zprime_plots/AK8_matched_to_genPtcles/qqDR_lessThan_0p8/M50/JetSubstructure_matched_histos_test.root","READ");

	f_input[2]=new TFile("/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/Zprime_plots/AK8_matched_to_genPtcles/qqDR_lessThan_0p8/M50/JetSubstructure_matched_histos_test.root","READ");
*/

	f_input[3]=new TFile("/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/test2_Abijith/CMSSW_11_1_0_pre8/src/pfpatatrackvalidation/ProducerTest/QCD_sample_plots/Full_statistics/matced/New_pT_split/pT_30_3000/JetSubstructure_matched_histos_matched_ak4_GausFits_v2.root","READ");


//	char legend_array[][500] = {"Z' (50 GeV)", "Z' (100 GeV)", "Z' (200 GeV)", "QCD"  };

	char legend_array[][500] = {"FullTracking", "PatatrackTight", "PatatrackLoose", "QCD"  };
	
	char eta_bins_legend[eta_bins][25];

	for (int iy=0; iy< eta_bins; iy++)
	{
		if (iy==0) sprintf( eta_bins_legend[iy], "|#eta|<%3.1f" , yBnd[iy+1] );  
		else
		{
			sprintf( eta_bins_legend[iy], "%3.1f<|#eta|<%3.1f" , yBnd[iy] ,yBnd[iy+1] );
		}
	}

	double pTmax = 1000 ;

	TH1D *h_ak8_pT[NoFiles][eta_bins], *h_ak8_phi[NoFiles][eta_bins], *h_ak8_m[NoFiles][eta_bins], *h_ak8_softDrop_m[NoFiles][eta_bins], *h_ak8_trimmed_m[NoFiles][eta_bins], *h_ak8_tau1[NoFiles][eta_bins], *h_ak8_tau2[NoFiles][eta_bins], *h_ak8_tau3[NoFiles][eta_bins], *h_ak8_tau4[NoFiles][eta_bins], *h_ak8_tau5[NoFiles][eta_bins], *h_ak8_tau21[NoFiles][eta_bins], *h_ak8_tau32[NoFiles][eta_bins], *h_ak8_tau43[NoFiles][eta_bins], *h_ak8_tau54[NoFiles][eta_bins], *h_genJet_pt[NoFiles][eta_bins], *h_genJet_mass[NoFiles][eta_bins];

	TH1D *h_ak8_eta[NoFiles], *h_minDR[NoFiles];


	char name[256]; 

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko();

//	char histo_type[NoFiles][100] = {"FullTracking", "FullTracking", "FullTracking", "FullTracking"};

	char histo_type[NoFiles][100] = {"FullTracking", "PatatrackTight", "PatatrackLoose", "FullTracking"};

	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		sprintf(name,"h_ak8_eta_%s",histo_type[NoFile]);
		h_ak8_eta[NoFile] = (TH1D*)(f_input[NoFile]->Get(name));

		sprintf(name,"h_minDR_%s",histo_type[NoFile]);
		h_minDR[NoFile] = (TH1D*)(f_input[NoFile]->Get(name));

		for(Int_t h=0; h<eta_bins; h++)
		{ 
			//========== reco jets===========
			sprintf(name,"h_ak8_pT_%s_bin%i",histo_type[NoFile],h);
			h_ak8_pT[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_phi_%s_bin%i",histo_type[NoFile],h);
			h_ak8_phi[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_m_%s_bin%i",histo_type[NoFile],h);
			h_ak8_m[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_softDrop_m_%s_bin%i",histo_type[NoFile],h);
			h_ak8_softDrop_m[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_trimmed_m_%s_bin%i",histo_type[NoFile],h);
			h_ak8_trimmed_m[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau1_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau1[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau2_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau2[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau3_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau3[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau4_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau4[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau5_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau5[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau21_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau21[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau32_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau32[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau43_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau43[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_ak8_tau54_%s_bin%i",histo_type[NoFile],h);
			h_ak8_tau54[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));


			sprintf(name,"h_genJet_pt_%s_bin%i",histo_type[NoFile],h);
			h_genJet_pt[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

			sprintf(name,"h_genJet_mass_%s_bin%i",histo_type[NoFile],h);
			h_genJet_mass[NoFile][h] = (TH1D*)(f_input[NoFile]->Get(name));

		}// end of etabin loop
	} // end of file loop
//======================================= Create ROC curves for tau21 ratio =========================================================

	const int ROCpoints = 20;
	
	double sigEfficiency[NoFiles-1][eta_bins][ROCpoints];
	double bkgRejection[eta_bins][ROCpoints];

	
	for (int NoFile=0; NoFile<NoFiles-1; NoFile++)
	{
		for(int eta_bin=0; eta_bin<eta_bins; eta_bin++)
		{
			for (int p=0; p<ROCpoints; p++)
			{
				double workingPoint = 0.5*(1./ROCpoints) + p*(1./ROCpoints);		
				int workingPointBin = h_ak8_tau21[NoFile][eta_bin]->FindBin(workingPoint);
				sigEfficiency[NoFile][eta_bin][p] = h_ak8_tau21[NoFile][eta_bin]->Integral(0,workingPointBin) / h_ak8_tau21[NoFile][eta_bin]->Integral();
				if(NoFile==0)	bkgRejection[eta_bin][p] = 1. - h_ak8_tau21[3][eta_bin]->Integral(0,workingPointBin) / h_ak8_tau21[3][eta_bin]->Integral();

			} 
		}
	}

	TGraph *gr_ROC_curve[NoFiles-1][eta_bins];

	for (int NoFile=0; NoFile<NoFiles-1; NoFile++)
	{
		for(int eta_bin=0; eta_bin<eta_bins; eta_bin++)
		{
			gr_ROC_curve[NoFile][eta_bin] = new TGraph(ROCpoints,sigEfficiency[NoFile][eta_bin], bkgRejection[eta_bin]);
		}
	}



//=======================================	    Create plots	    =========================================================

	int Colors[] = { 1, 4, 2 , 6, 3, 7 , 28, 46} ; // black, blue, red , magenta, green, light blue ,  brown, redish. // define colors for each input
	int MarkerStyle[] = { 8, 2, 5 , 4, 22, 21, 27, 28 } ; // 

	int PadColumnsEtaBins = 3;  // define number of pad columns for the canvas
	int PadRowsEtaBins = 2;     // define number of pad rows for the canvas

//for nicely looking plots it should be : eta_bins <= PadColumnsEtaBins * PadRowsEtaBins  

	int Canvas_XpixelsEtaBins = PadColumnsEtaBins*333; //define canvas X pixels
	int Canvas_YpixelsEtaBins = PadRowsEtaBins*500;    //define canvas Y pixels

	double YaxisLowEndMultiplier = 0.00005; // define lower bound of Y axis for plot: Lower Bound = YaxisLowEndMultiplier* histo_maximum
	double YaxisHighEndMultiplier = 2.2;    // define upper bound of Y axis for plot: Upper Bound = YaxisHighEndMultiplier* histo_maximum

	TCanvas *pad_pt = new TCanvas("pad_pt", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_pt->Divide(PadColumnsEtaBins,PadRowsEtaBins);

	TCanvas *pad_gen_pt = new TCanvas("pad_gen_pt", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_gen_pt->Divide(PadColumnsEtaBins,PadRowsEtaBins);

	TCanvas *pad_gen_mass = new TCanvas("pad_gen_mass", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_gen_mass->Divide(PadColumnsEtaBins,PadRowsEtaBins);


	TCanvas *pad_phi = new TCanvas("pad_phi", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_phi->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_m = new TCanvas("pad_m", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_m->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_softDrop_m = new TCanvas("pad_softDrop_m", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_softDrop_m->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_trimmed_m = new TCanvas("pad_trimmed_m", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_trimmed_m->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_tau1 = new TCanvas("pad_tau1", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau1->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_tau2 = new TCanvas("pad_tau2", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau2->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_tau3 = new TCanvas("pad_tau3", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau3->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_tau4 = new TCanvas("pad_tau4", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau4->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_tau5 = new TCanvas("pad_tau5", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau5->Divide(PadColumnsEtaBins,PadRowsEtaBins);

	TCanvas *pad_tau21 = new TCanvas("pad_tau21", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau21->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_tau32 = new TCanvas("pad_tau32", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau32->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_tau43 = new TCanvas("pad_tau43", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau43->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	TCanvas *pad_tau54 = new TCanvas("pad_tau54", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_tau54->Divide(PadColumnsEtaBins,PadRowsEtaBins);
	
	TCanvas *pad_ROC_curve = new TCanvas("pad_ROC_curve", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad_ROC_curve->Divide(PadColumnsEtaBins,PadRowsEtaBins);

	TCanvas *c_eta = new TCanvas("c_eta", "",1);
	TCanvas *c_minDR = new TCanvas("c_minDR", "",1);



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

	char res_text[200];
	TLegend *leg2 =new TLegend(.5, .7, .7, .9);//7899//4899
	leg2->SetTextSize(0.03);
	leg2->SetFillColor(0); 
	leg2->SetBorderSize(0);

	TLegend *leg4 =new TLegend(.1, .6, .9, .9);//7899//4899
	leg4->SetTextSize(0.055);
	leg4->SetFillColor(0); 
	leg4->SetBorderSize(0);

	 //dummy histograms to be used as frames in the plots
	TH1D *framept[eta_bins], *frameGenpt[eta_bins], *frameGenMass[eta_bins], *framephi[eta_bins], *framem[eta_bins], *framesoftDrop_m[eta_bins], *frametrimmed_m[eta_bins], *frametau1[eta_bins], *frametau2[eta_bins], *frametau3[eta_bins], *frametau4[eta_bins], *frametau5[eta_bins], *frametau21[eta_bins], *frametau32[eta_bins], *frametau43[eta_bins], *frametau54[eta_bins], *frameROCurve[eta_bins], *frameEta, *frameMinDR;

	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		leg1->AddEntry(h_ak8_pT[NoFile][0], legend_array[NoFile], "L");
		leg2->AddEntry(h_ak8_pT[NoFile][0], legend_array[NoFile], "L");
		if(NoFile<NoFiles-1) leg4->AddEntry(h_ak8_pT[NoFile][0], legend_array[NoFile], "L");

		if (NoFile==0 )  frameEta = InitiateFrameOnCanvasPad(c_eta, 0 , "frameEta", "Jet eta", "Entries", -5., 5., YaxisLowEndMultiplier*h_ak8_eta[NoFile]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_eta[NoFile]->GetMaximum(), true, paveCMS);
		if (scale_histos && h_ak8_eta[NoFile]->Integral()>0 ) h_ak8_eta[NoFile]->Scale(h_ak8_eta[0]->Integral() / h_ak8_eta[NoFile]->Integral());
		DrawHistoToCanvasPad(c_eta, 0, h_ak8_eta[NoFile], Colors[NoFile], 1);
		if( NoFile== 0)  leg2->Draw("same"); 

		if(plot_minDR)
		{
			if (NoFile==0 )  frameMinDR = InitiateFrameOnCanvasPad(c_minDR, 0 , "frameMinDR", "minDR", "Entries", 0., 5., 0.1*YaxisLowEndMultiplier*h_minDR[NoFile]->GetMaximum(), YaxisHighEndMultiplier*h_minDR[NoFile]->GetMaximum(), true, paveCMS);
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

			if (NoFile==0 )  framept[iy] = InitiateFrameOnCanvasPad(pad_pt, iy+1, "framept", "Jet pT (GeV)", "Entries", 0., pTmax, YaxisLowEndMultiplier*h_ak8_pT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_pT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_pT[NoFile][iy]->Integral()>0 ) h_ak8_pT[NoFile][iy]->Scale(h_ak8_pT[0][iy]->Integral()/h_ak8_pT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_pt, iy+1, h_ak8_pT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_pt->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_pt->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frameGenpt[iy] = InitiateFrameOnCanvasPad(pad_gen_pt, iy+1, "frameGenpt", "gen Jet pT (GeV)", "Entries", 0., pTmax, YaxisLowEndMultiplier*h_genJet_pt[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_genJet_pt[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_genJet_pt[NoFile][iy]->Integral()>0 ) h_genJet_pt[NoFile][iy]->Scale(h_genJet_pt[0][iy]->Integral()/h_genJet_pt[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_gen_pt, iy+1, h_genJet_pt[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_gen_pt->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_gen_pt->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frameGenMass[iy] = InitiateFrameOnCanvasPad(pad_gen_mass, iy+1, "frameGenMass", "gen jet mass (GeV)", "Entries", 0., 200, 0.8, 1.3*h_genJet_mass[NoFile][iy]->GetMaximum() , false, paveCMS);
			if (scale_histos && h_genJet_mass[NoFile][iy]->Integral()>0 ) h_genJet_mass[NoFile][iy]->Scale(h_genJet_mass[0][iy]->Integral()/h_genJet_mass[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_gen_mass, iy+1, h_genJet_mass[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_gen_mass->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_gen_mass->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  framephi[iy] = InitiateFrameOnCanvasPad(pad_phi, iy+1, "framephi", "Jet phi", "Entries", -3.14, 3.14 , YaxisLowEndMultiplier*h_ak8_phi[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_phi[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_phi[NoFile][iy]->Integral()>0 ) h_ak8_phi[NoFile][iy]->Scale(h_ak8_phi[0][iy]->Integral()/h_ak8_phi[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_phi, iy+1, h_ak8_phi[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_phi->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_phi->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  framem[iy] = InitiateFrameOnCanvasPad(pad_m, iy+1, "framem", "Jet mass (GeV) ", "Entries", 0., 200 , 0.8, 1.3*h_ak8_m[NoFile][iy]->GetMaximum() , false, paveCMS);
			if (scale_histos && h_ak8_m[NoFile][iy]->Integral()>0 ) h_ak8_m[NoFile][iy]->Scale(h_ak8_m[0][iy]->Integral()/h_ak8_m[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_m, iy+1, h_ak8_m[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_m->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_m->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  framesoftDrop_m[iy] = InitiateFrameOnCanvasPad(pad_softDrop_m, iy+1, "framesoftDrop_m", "Jet SoftDrop mass (GeV) ", "Entries", 0., 300 , 0.5, YaxisHighEndMultiplier*h_ak8_softDrop_m[NoFile][iy]->GetMaximum() , false, paveCMS);
			if (scale_histos && h_ak8_softDrop_m[NoFile][iy]->Integral()>0 ) h_ak8_softDrop_m[NoFile][iy]->Scale(h_ak8_softDrop_m[0][iy]->Integral()/h_ak8_softDrop_m[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_softDrop_m, iy+1, h_ak8_softDrop_m[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_softDrop_m->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_softDrop_m->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametrimmed_m[iy] = InitiateFrameOnCanvasPad(pad_trimmed_m, iy+1, "frametrimmed_m", "Jet trimmed mass (GeV) ", "Entries", 0., 200 , 0.8, 1.3*h_ak8_trimmed_m[NoFile][iy]->GetMaximum() , false, paveCMS);
			if (scale_histos && h_ak8_trimmed_m[NoFile][iy]->Integral()>0 ) h_ak8_trimmed_m[NoFile][iy]->Scale(h_ak8_trimmed_m[0][iy]->Integral()/h_ak8_trimmed_m[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_trimmed_m, iy+1, h_ak8_trimmed_m[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_trimmed_m->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_trimmed_m->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau1[iy] = InitiateFrameOnCanvasPad(pad_tau1, iy+1, "frametau1", "#tau_{1} ", "Entries", 0., 400.0 , YaxisLowEndMultiplier*h_ak8_tau1[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_tau1[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_tau1[NoFile][iy]->Integral()>0 ) h_ak8_tau1[NoFile][iy]->Scale(h_ak8_tau1[0][iy]->Integral()/h_ak8_tau1[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau1, iy+1, h_ak8_tau1[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau1->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau1->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau2[iy] = InitiateFrameOnCanvasPad(pad_tau2, iy+1, "frametau2", "#tau_{2} ", "Entries", 0., 400.0 , YaxisLowEndMultiplier*h_ak8_tau2[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_tau2[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_tau2[NoFile][iy]->Integral()>0 ) h_ak8_tau2[NoFile][iy]->Scale(h_ak8_tau2[0][iy]->Integral()/h_ak8_tau2[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau2, iy+1, h_ak8_tau2[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau2->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau2->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau3[iy] = InitiateFrameOnCanvasPad(pad_tau3, iy+1, "frametau3", "#tau_{3} ", "Entries", 0., 400.0 , YaxisLowEndMultiplier*h_ak8_tau3[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_tau3[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_tau3[NoFile][iy]->Integral()>0 ) h_ak8_tau3[NoFile][iy]->Scale(h_ak8_tau3[0][iy]->Integral()/h_ak8_tau3[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau3, iy+1, h_ak8_tau3[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau3->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau3->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau4[iy] = InitiateFrameOnCanvasPad(pad_tau4, iy+1, "frametau4", "#tau_{4} ", "Entries", 0., 400.0 , YaxisLowEndMultiplier*h_ak8_tau4[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_tau4[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_tau4[NoFile][iy]->Integral()>0 ) h_ak8_tau4[NoFile][iy]->Scale(h_ak8_tau4[0][iy]->Integral()/h_ak8_tau4[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau4, iy+1, h_ak8_tau4[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau4->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau4->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau5[iy] = InitiateFrameOnCanvasPad(pad_tau5, iy+1, "frametau5", "#tau_{5} ", "Entries", 0., 400.0 , YaxisLowEndMultiplier*h_ak8_tau5[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_tau5[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_tau5[NoFile][iy]->Integral()>0 ) h_ak8_tau5[NoFile][iy]->Scale(h_ak8_tau5[0][iy]->Integral()/h_ak8_tau5[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau5, iy+1, h_ak8_tau5[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau5->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau5->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau21[iy] = InitiateFrameOnCanvasPad(pad_tau21, iy+1, "frametau21", "#tau_{2} / #tau_{1} ", "Entries", 0., 1.2 , 0.5, YaxisHighEndMultiplier*h_ak8_tau21[NoFile][iy]->GetMaximum() , false, paveCMS);
			if (scale_histos && h_ak8_tau21[NoFile][iy]->Integral()>0 ) h_ak8_tau21[NoFile][iy]->Scale(h_ak8_tau21[0][iy]->Integral()/h_ak8_tau21[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau21, iy+1, h_ak8_tau21[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau21->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau21->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau32[iy] = InitiateFrameOnCanvasPad(pad_tau32, iy+1, "frametau32", "#tau_{3} / #tau_{2} ", "Entries", 0., 1.2 , YaxisLowEndMultiplier*h_ak8_tau32[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_tau32[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_tau32[NoFile][iy]->Integral()>0 ) h_ak8_tau32[NoFile][iy]->Scale(h_ak8_tau32[0][iy]->Integral()/h_ak8_tau32[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau32, iy+1, h_ak8_tau32[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau32->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau32->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau43[iy] = InitiateFrameOnCanvasPad(pad_tau43, iy+1, "frametau43", "#tau_{4} / #tau_{3} ", "Entries", 0., 1.2 , YaxisLowEndMultiplier*h_ak8_tau43[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_tau43[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_tau43[NoFile][iy]->Integral()>0 ) h_ak8_tau43[NoFile][iy]->Scale(h_ak8_tau43[0][iy]->Integral()/h_ak8_tau43[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau43, iy+1, h_ak8_tau43[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau43->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau43->cd(iy+1); teta->Draw(); }

			if (NoFile==0 )  frametau54[iy] = InitiateFrameOnCanvasPad(pad_tau54, iy+1, "frametau54", "#tau_{5} / #tau_{4} ", "Entries", 0., 1.2 , YaxisLowEndMultiplier*h_ak8_tau54[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ak8_tau54[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ak8_tau54[NoFile][iy]->Integral()>0 ) h_ak8_tau54[NoFile][iy]->Scale(h_ak8_tau54[0][iy]->Integral()/h_ak8_tau54[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_tau54, iy+1, h_ak8_tau54[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_tau54->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_tau54->cd(iy+1); teta->Draw(); }


			if(NoFile<NoFiles-1)
			{
				if (NoFile==0 ) frameROCurve[iy] = InitiateFrameOnCanvasPad(pad_ROC_curve, iy+1, "frameROCurve", "Signal efficiency", "Background rejection", 0., 1.0, 0., 1.2, false, paveCMS);
				DrawGraphToCanvasPad(pad_ROC_curve, iy+1, gr_ROC_curve[NoFile][iy], Colors[NoFile], 1, false);
				if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_ROC_curve->cd(iy+1); leg4->Draw("same");  }
				else if( eta_bins > 1 && NoFile== 0) { pad_ROC_curve->cd(iy+1); teta->Draw(); }
			}

		}
	 

		if(eta_bins > 1 &&  NoFile == NoFiles-1 )
		{
			pad_pt->cd(eta_bins+1);	leg1->Draw();
			pad_gen_pt->cd(eta_bins+1);	leg1->Draw();
			pad_gen_mass->cd(eta_bins+1);	leg1->Draw();
			pad_phi->cd(eta_bins+1);	leg1->Draw();
			pad_m->cd(eta_bins+1);	leg1->Draw();
			pad_softDrop_m->cd(eta_bins+1);	leg1->Draw();
			pad_trimmed_m->cd(eta_bins+1);	leg1->Draw();
			pad_tau1->cd(eta_bins+1);	leg1->Draw();
			pad_tau2->cd(eta_bins+1);	leg1->Draw();
			pad_tau3->cd(eta_bins+1);	leg1->Draw();
			pad_tau4->cd(eta_bins+1);	leg1->Draw();
			pad_tau5->cd(eta_bins+1);	leg1->Draw();
			pad_tau21->cd(eta_bins+1);	leg1->Draw();
			pad_tau32->cd(eta_bins+1);	leg1->Draw();
			pad_tau43->cd(eta_bins+1);	leg1->Draw();
			pad_tau54->cd(eta_bins+1);	leg1->Draw();

			pad_ROC_curve->cd(eta_bins+1); leg4->Draw();
			
		}

	} // end of loop on files

	char filename[1500];
	if (Save_Plots)
	{ 
		sprintf(filename,"%s/%s/%s_Reclustered_pT.png",analyzer_path,output_directory,image_name);
		pad_pt->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_genpT.png",analyzer_path,output_directory,image_name);
		pad_gen_pt->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_genMass.png",analyzer_path,output_directory,image_name);
		pad_gen_mass->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_phi.png",analyzer_path,output_directory,image_name);
		pad_phi->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_m.png",analyzer_path,output_directory,image_name);
		pad_m->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_softDrop_m.png",analyzer_path,output_directory,image_name);
		pad_softDrop_m->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_trimmed_m.png",analyzer_path,output_directory,image_name);
		pad_trimmed_m->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau1.png",analyzer_path,output_directory,image_name);
		pad_tau1->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau2.png",analyzer_path,output_directory,image_name);
		pad_tau2->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau3.png",analyzer_path,output_directory,image_name);
		pad_tau3->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau4.png",analyzer_path,output_directory,image_name);
		pad_tau4->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau5.png",analyzer_path,output_directory,image_name);
		pad_tau5->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau21.png",analyzer_path,output_directory,image_name);
		pad_tau21->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau32.png",analyzer_path,output_directory,image_name);
		pad_tau32->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau43.png",analyzer_path,output_directory,image_name);
		pad_tau43->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_tau54.png",analyzer_path,output_directory,image_name);
		pad_tau54->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_tau21_ROC_curve.png",analyzer_path,output_directory,image_name);
		pad_ROC_curve->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_eta.png",analyzer_path,output_directory,image_name);
		c_eta->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_Reclustered_minDR.png",analyzer_path,output_directory,image_name);
		c_minDR->SaveAs(filename);

		
   }
}

