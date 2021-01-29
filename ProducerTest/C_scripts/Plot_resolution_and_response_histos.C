#include "include/include_functions.h"
#include "include/Input_definition.h"


void Plot_resolution_and_response_histos()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);
	gROOT->ForceStyle(kTRUE);
	gStyle->SetOptFit(0);


	double NsigmaTail = 3.0;


	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;
	const int pT_bins = sizeof(ptBnd)/sizeof(ptBnd[0])-1;
	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);


	double pTraw_responses[NoFiles][eta_bins][pT_bins], pTraw_resolution[NoFiles][eta_bins][pT_bins];
	double pTraw_responses_error[NoFiles][eta_bins][pT_bins], pTraw_resolution_error[NoFiles][eta_bins][pT_bins];

	double pTraw_resolution_Gaus[NoFiles][eta_bins][pT_bins], pTraw_response_Gaus[NoFiles][eta_bins][pT_bins];
	double pTraw_resolution_Gaus_error[NoFiles][eta_bins][pT_bins], pTraw_response_Gaus_error[NoFiles][eta_bins][pT_bins];

	double pT_resolution_corrected[NoFiles][eta_bins][pT_bins], pT_resolution_corrected_Gaus[NoFiles][eta_bins][pT_bins];
	double pTreco_ov_pTgen_JEScorrected_error[NoFiles][eta_bins][pT_bins], pTreco_ov_pTgen_JEScorrected_Gaus_error[NoFiles][eta_bins][pT_bins];
	double tailFrac_pTreco_ov_pTgen[NoFiles][eta_bins][pT_bins], tailFrac_pTreco_ov_pTgen_highErr[NoFiles][eta_bins][pT_bins], tailFrac_pTreco_ov_pTgen_lowErr[NoFiles][eta_bins][pT_bins];


	double pTbinCenter[pT_bins], pTbinCenter_error[pT_bins];

	
	TGraphAsymmErrors *gr_pT_response[NoFiles][eta_bins], *gr_pT_resolution[NoFiles][eta_bins], *gr_tailFrac_pTreco_ov_pTgen[NoFiles][eta_bins];

	TH1D *h_pTreco_ov_pTgen[NoFiles][eta_bins][pT_bins], *h_pTreco_ov_pTgen_JEScorrected[NoFiles][eta_bins][pT_bins], *h_pTreco_ov_pTgen_JEScorrected_noWeights[NoFiles][eta_bins][pT_bins];

//	cout << "eta bins = " << eta_bins;
//	cout << "pT bins = " << pT_bins;


	char name[1500]; 
	sprintf(name, "%s/%s/pTreco_ov_pTgen_%s_histos.root",analyzer_path,output_directory,image_name );
	TFile *f_input_raw = new TFile (name,"READ");

	sprintf(name, "%s/%s/pTreco_ov_pTgen_%s_JEScorrected_histos.root", analyzer_path,output_directory,image_name );
	TFile *f_input_JEScorrected = new TFile (name,"READ");


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int  eta_bin=0; eta_bin<eta_bins; eta_bin++)
		{
			for (int  pT_bin=0; pT_bin<pT_bins; pT_bin++)
			{
				sprintf(name,"h_pTreco_ov_pTgen_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin] = (TH1D*)(f_input_raw->Get(name));

				sprintf(name,"h_pTreco_ov_pTgen_JEScorrected_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin] = (TH1D*)(f_input_JEScorrected->Get(name));

				sprintf(name,"h_pTreco_ov_pTgen_JEScorrected_noWeights_%s_EtaBin%i_ptBin%i",legend_array[NoFile],eta_bin,pT_bin);
				h_pTreco_ov_pTgen_JEScorrected_noWeights[NoFile][eta_bin][pT_bin] = (TH1D*)(f_input_JEScorrected->Get(name));
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

	char pT_bins_legend[pT_bins][50];
	for (int iy=0; iy< pT_bins; iy++)
	{

		sprintf( pT_bins_legend[iy], "%3.0f<pT<%3.0f" , ptBnd[iy] ,ptBnd[iy+1] );
/*
		if (iy==0) sprintf( pT_bins_legend[iy], "pT<%f" , ptBnd[iy+1] );  
		else if (iy < pT_bins -1 )
		{
			sprintf( pT_bins_legend[iy], "%f<pT<%f" , ptBnd[iy] ,ptBnd[iy+1] );
		}
		else sprintf( pT_bins_legend[iy], "p_T>%f" , ptBnd[iy] );  
*/
	}


	TF1 *gausTF1; 
 	const Int_t kNotDraw = 1<<9;
	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{

		cout << " Calculating raw responses and resolutions for file" << input_files[NoFile] << " ...  "<< endl;

		for ( int eta_bin=0; eta_bin< eta_bins; eta_bin++ )
		{

			for (int pT_bin=0; pT_bin< pT_bins; pT_bin++ )
			{
				cout << "Now processing:   " << eta_bins_legend[eta_bin] << "   &   " << pT_bins_legend[pT_bin] << endl;
				double Nsize = h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->Integral();
				if ( !(Nsize>0) ) 
				{
					cout << "Warning!!! Histogram empty!!  eta bin : " << eta_bins_legend[eta_bin] << "  pT : " << pT_bins_legend[pT_bin] << endl;
					continue;
				}

				double res_mean = h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->GetMean();
				double res_rms = h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->GetRMS(); //
				double res_mean_error = res_rms / Nsize; // error of mean value.
				double res_rms_error = 	CalculateRMSerror( h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin] ); //error of standard deviation found here : https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation


				pTraw_responses[NoFile][eta_bin][pT_bin] = res_mean;
				pTraw_resolution[NoFile][eta_bin][pT_bin] = res_rms;
				pTraw_responses_error[NoFile][eta_bin][pT_bin] = res_mean_error;
				pTraw_resolution_error[NoFile][eta_bin][pT_bin] = res_rms_error;

				if (res_rms>0)
				{
				
					//if( pTbinCenter[pT_bin] < maxPtForGausReducedRangeFit) h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->Fit("gaus","0","0",res_mean-2.2*res_rms,res_mean+0.8*res_rms);
					//else h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->Fit("gaus","0","0",res_mean-2.2*res_rms,res_mean+1.8*res_rms);
					//TF1 *gausTF1 = (TF1*)h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->GetListOfFunctions()->FindObject("gaus");
					gausTF1 = FindBestGaussianCoreFit(h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]);

					pTraw_response_Gaus[NoFile][eta_bin][pT_bin] = gausTF1->GetParameter(1);
					pTraw_response_Gaus_error[NoFile][eta_bin][pT_bin] = gausTF1->GetParError(1);

					pTraw_resolution_Gaus[NoFile][eta_bin][pT_bin] = fabs(gausTF1->GetParameter(2));
					pTraw_resolution_Gaus_error[NoFile][eta_bin][pT_bin] = gausTF1->GetParError(2);
				}
				else continue;
				

			} // end of loop of pT bins

			

			if (useGausFits) gr_pT_response[NoFile][eta_bin] = new TGraphAsymmErrors( pT_bins, pTbinCenter, pTraw_response_Gaus[NoFile][eta_bin], pTbinCenter_error, pTbinCenter_error, pTraw_response_Gaus_error[NoFile][eta_bin], pTraw_response_Gaus_error[NoFile][eta_bin] );	
			else gr_pT_response[NoFile][eta_bin] = new TGraphAsymmErrors( pT_bins, pTbinCenter, pTraw_responses[NoFile][eta_bin], pTbinCenter_error, pTbinCenter_error, pTraw_responses_error[NoFile][eta_bin], pTraw_responses_error[NoFile][eta_bin] );
			//gr_pT_response[NoFile][eta_bin] = new TGraphAsymmErrors( pT_bins, pTbinCenter, pTraw_responses[NoFile][eta_bin], pTbinCenter_error, pTbinCenter_error, pTraw_responses_error[NoFile][eta_bin], pTraw_responses_error[NoFile][eta_bin] );	

		}// end of loop of eta bins



		cout << " Calculating (JES corrected) resolutions for file" << input_files[NoFile] << " ...  "<< endl;

		for ( int eta_bin=0; eta_bin< eta_bins; eta_bin++ )
		{
			for (int pT_bin=0; pT_bin< pT_bins; pT_bin++ )
			{

				cout << "Now processing:   " << eta_bins_legend[eta_bin] << "   &   " << pT_bins_legend[pT_bin] << endl;
				if ( !(h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin]->Integral()>0) ) 
				{
					cout << "Warning!!! Histogram empty!!  eta bin : " << eta_bins_legend[eta_bin] << "  pT : " << pT_bins_legend[pT_bin] << endl;
					continue;
				}
				double res_rms = h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->GetRMS(); //
				double res_rms_error = 	CalculateRMSerror( h_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin] );
				pT_resolution_corrected[NoFile][eta_bin][pT_bin] = res_rms;
				pTreco_ov_pTgen_JEScorrected_error[NoFile][eta_bin][pT_bin] = res_rms_error;

				if (res_rms>0)
				{
				
					//if( pTbinCenter[pT_bin] < maxPtForGausReducedRangeFit) h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Fit("gaus","0","0",1.-2.2*res_rms,1.+0.8*res_rms);
					//else h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Fit("gaus","0","0",1.-2.2*res_rms,1.+1.8*res_rms);
					//TF1 *gausTF1 = (TF1*)h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->GetListOfFunctions()->FindObject("gaus");
					gausTF1 = FindBestGaussianCoreFit(h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]);

					pT_resolution_corrected_Gaus[NoFile][eta_bin][pT_bin] = fabs(gausTF1->GetParameter(2));
					pTreco_ov_pTgen_JEScorrected_Gaus_error[NoFile][eta_bin][pT_bin] = gausTF1->GetParError(2);
				
					double scaleForTail = h_pTreco_ov_pTgen_JEScorrected_noWeights[NoFile][eta_bin][pT_bin]->Integral() / h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Integral();
					h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Scale(scaleForTail); // scale to unweighted statistics to calc Wilson intervals
					
					double tail_Njet_low  = h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Integral( h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->FindBin(0.), h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->FindBin( gausTF1->GetParameter(1) - NsigmaTail*fabs(gausTF1->GetParameter(2)) ) );
					double tail_Njet_high = h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Integral( h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->FindBin( gausTF1->GetParameter(1) + NsigmaTail*fabs(gausTF1->GetParameter(2)) ), h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->FindBin( 6.0 ) );

					tailFrac_pTreco_ov_pTgen[NoFile][eta_bin][pT_bin] = (tail_Njet_low + tail_Njet_high) / h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Integral();
					tailFrac_pTreco_ov_pTgen_highErr[NoFile][eta_bin][pT_bin] = error1( h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Integral(), tail_Njet_low + tail_Njet_high );
					tailFrac_pTreco_ov_pTgen_lowErr[NoFile][eta_bin][pT_bin]  = error2( h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Integral(), tail_Njet_low + tail_Njet_high );

				//	h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Scale(1./scaleForTail); // scale back to weighted statistics (currently useless)					
				}
				else continue;

			} // end of pT bin loop

			gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin] = new TGraphAsymmErrors(pT_bins, pTbinCenter,tailFrac_pTreco_ov_pTgen[NoFile][eta_bin], pTbinCenter_error, pTbinCenter_error, tailFrac_pTreco_ov_pTgen_lowErr[NoFile][eta_bin], tailFrac_pTreco_ov_pTgen_highErr[NoFile][eta_bin] );

			if (useGausFits) gr_pT_resolution[NoFile][eta_bin] = new TGraphAsymmErrors(pT_bins, pTbinCenter, pT_resolution_corrected_Gaus[NoFile][eta_bin], pTbinCenter_error, pTbinCenter_error, pTreco_ov_pTgen_JEScorrected_Gaus_error[NoFile][eta_bin], pTreco_ov_pTgen_JEScorrected_Gaus_error[NoFile][eta_bin] );
			else gr_pT_resolution[NoFile][eta_bin] = new TGraphAsymmErrors(pT_bins, pTbinCenter, pT_resolution_corrected[NoFile][eta_bin], pTbinCenter_error, pTbinCenter_error, pTreco_ov_pTgen_JEScorrected_error[NoFile][eta_bin], pTreco_ov_pTgen_JEScorrected_error[NoFile][eta_bin] ); // creating resolution graphs

		} // end of eta bin loop
	} // end of loop on files




// ======================================================plotting stuff======================================//
	char res_text[400];
	TH1D *frame_h_pTreco_ov_pTgen[eta_bins][pT_bins];


	TCanvas *ptResponsePad = new TCanvas("ptResponsePad", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	ptResponsePad->Divide(PadColumnsEtaBins,PadRowsEtaBins);

	TCanvas *ptResolutionPad = new TCanvas("ptResolutionPad", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	ptResolutionPad->Divide(PadColumnsEtaBins,PadRowsEtaBins);

	TCanvas *pTresTailPad = new TCanvas("pTresTailPad", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pTresTailPad->Divide(PadColumnsEtaBins,PadRowsEtaBins);

	TCanvas *EtaPtResolutionPad[eta_bins];
	
	for (int eta_bin = 0; eta_bin<eta_bins; eta_bin++)
	{
		sprintf(name,"EtaPtResolutionPad_etaBin%i",eta_bin);
		EtaPtResolutionPad[eta_bin] = new TCanvas(name, eta_bins_legend[eta_bin],Canvas_XpixelsPtBins,Canvas_YpixelsPtBins);
		EtaPtResolutionPad[eta_bin]->Divide(PadColumnsPtBins,PadRowsPtBins);
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




			pTresTailPad->cd(eta_bin+1);
			pTresTailPad->cd(eta_bin+1)->SetLogx(1);
			if(NoFile == 0 )
			{	
				TH1F *hr = pTresTailPad->cd(eta_bin+1)->DrawFrame(ptBnd[0],0.,ptBnd[pT_bins],0.3);
				//if(eta_bin_counter==5) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.8,3000,1.2);
				//if(eta_bin_counter==6) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.0,3000,1.5);
				hr->SetXTitle("jet pT (GeV)");
				hr->SetYTitle("pT RES tail frac");
				hr->GetYaxis()->SetTitleOffset(1.3);
				hr->GetXaxis()->SetTitleOffset(1.3);
				hr->SetTitle(eta_bins_legend[eta_bin]);
			}

			gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin]->SetLineColor(Colors[NoFile]);
			gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin]->SetLineWidth(2);
			gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin]->SetLineStyle(1);
			gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin]->SetMarkerColor(Colors[NoFile]);
//			gr_pT_response[NoFile][eta_bin]->SetMarkerStyle(MarkerStyle[NoFile]);
			gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin]->SetMarkerStyle(24);
			gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin]->SetMarkerSize(0.3);
	//		gr_pT_response[NoFile][eta_bin]->GetYaxis()->SetRangeUser(0.5, 1.5);
//			gr_pT_response[NoFile][eta_bin]->GetYaxis()->SetTitle("Response");
//			gr_pT_response[NoFile][eta_bin]->GetXaxis()->SetTitle("jet pT (GeV)");
			if ( NoFile==0 ) 	{	gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin]->Draw("p");	paveCMS ->Draw("same");  }
			else gr_tailFrac_pTreco_ov_pTgen[NoFile][eta_bin]->Draw("same p");
			paveCMS->Draw("same");
			etaTextLatx->Draw();




			ptResponsePad->cd(eta_bin+1);
			ptResponsePad->cd(eta_bin+1)->SetLogx(1);
			if(NoFile == 0 )
			{	
				TH1F *hr = ptResponsePad->cd(eta_bin+1)->DrawFrame(ptBnd[0],0.5,ptBnd[pT_bins],1.5);
				//if(eta_bin_counter==5) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.8,3000,1.2);
				//if(eta_bin_counter==6) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.0,3000,1.5);
				hr->SetXTitle("jet pT (GeV)");
				hr->SetYTitle("jet pT Response");
				hr->GetYaxis()->SetTitleOffset(1.3);
				hr->GetXaxis()->SetTitleOffset(1.3);
				hr->SetTitle(eta_bins_legend[eta_bin]);
			}



			gr_pT_response[NoFile][eta_bin]->SetLineColor(Colors[NoFile]);
			gr_pT_response[NoFile][eta_bin]->SetLineWidth(2);
			gr_pT_response[NoFile][eta_bin]->SetLineStyle(1);
			gr_pT_response[NoFile][eta_bin]->SetMarkerColor(Colors[NoFile]);
//			gr_pT_response[NoFile][eta_bin]->SetMarkerStyle(MarkerStyle[NoFile]);
			gr_pT_response[NoFile][eta_bin]->SetMarkerStyle(24);
			gr_pT_response[NoFile][eta_bin]->SetMarkerSize(0.3);
	//		gr_pT_response[NoFile][eta_bin]->GetYaxis()->SetRangeUser(0.5, 1.5);
//			gr_pT_response[NoFile][eta_bin]->GetYaxis()->SetTitle("Response");
//			gr_pT_response[NoFile][eta_bin]->GetXaxis()->SetTitle("jet pT (GeV)");
			if ( NoFile==0 ) 	{	gr_pT_response[NoFile][eta_bin]->Draw("p");	paveCMS ->Draw("same");  }
			else gr_pT_response[NoFile][eta_bin]->Draw("same p");
			paveCMS->Draw("same");
			etaTextLatx->Draw();


			
			ptResolutionPad->cd(eta_bin+1);
			ptResolutionPad->cd(eta_bin+1)->SetLogx(1);
			if(NoFile == 0 )
			{	
				TH1F *hr = ptResolutionPad->cd(eta_bin+1)->DrawFrame(ptBnd[0],0.001,ptBnd[pT_bins],0.4);
				//if(eta_bin_counter==5) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.8,3000,1.2);
				//if(eta_bin_counter==6) TH1F *hr = pad1->cd(eta_bin_counter+1)->DrawFrame(0,0.0,3000,1.5);
				hr->SetXTitle("jet pT (GeV)");
				hr->SetYTitle("jet pT Resolution");
				hr->GetYaxis()->SetTitleOffset(1.3);
				hr->GetXaxis()->SetTitleOffset(1.3);
				hr->SetTitle(eta_bins_legend[eta_bin]);
			}

			gr_pT_resolution[NoFile][eta_bin]->SetLineColor(Colors[NoFile]);
			gr_pT_resolution[NoFile][eta_bin]->SetLineWidth(2);
			gr_pT_resolution[NoFile][eta_bin]->SetLineStyle(1);
			gr_pT_resolution[NoFile][eta_bin]->SetMarkerColor(Colors[NoFile]);
//			gr_pT_resolution[NoFile][eta_bin]->SetMarkerStyle(MarkerStyle[NoFile]);
			gr_pT_resolution[NoFile][eta_bin]->SetMarkerStyle(24);
			gr_pT_resolution[NoFile][eta_bin]->SetMarkerSize(0.3);
//			gr_pT_resolution[NoFile][eta_bin]->GetYaxis()->SetRangeUser(0.02, 0.4);
//			gr_pT_resolution[NoFile][eta_bin]->GetYaxis()->SetTitle("jet pT resolution");
//			gr_pT_resolution[NoFile][eta_bin]->GetXaxis()->SetTitle("jet pT (GeV)");
			if ( NoFile==0 ) 	{	gr_pT_resolution[NoFile][eta_bin]->Draw("p");	paveCMS ->Draw("same");  }
			else gr_pT_resolution[NoFile][eta_bin]->Draw("same p");
			paveCMS->Draw("same");
			etaTextLatx->Draw();



			leg1[eta_bin]->AddEntry(h_pTreco_ov_pTgen_JEScorrected[NoFile][0][0], legend_array[NoFile], "L");

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

				
				if (!useGausFits) sprintf(res_text, " RES = %3.1f %%", 100*pT_resolution_corrected[NoFile][eta_bin][pT_bin]);
				else sprintf(res_text, " RES = %3.1f %%", 100*pT_resolution_corrected_Gaus[NoFile][eta_bin][pT_bin]);
				leg3[eta_bin][pT_bin]->AddEntry(h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin], res_text , "L");


				h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Rebin(5);

				if(NoFile==0)
				{
					sprintf(name,"frame_h_pTreco_ov_pTgen_etaBin%i_ptBin%i",eta_bin, pT_bin);
					frame_h_pTreco_ov_pTgen[eta_bin][pT_bin] = InitiateFrameOnCanvasPad(EtaPtResolutionPad[eta_bin], pT_bin+1 , name, "pT_reco / pT_gen", "Entries", 0., 2.1, 0.9, 800*h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->GetMaximum(), true, paveCMS);
					
				}

				
				if(h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->Integral()>0 && h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->GetRMS()>0 )
				{
					gausTF1 = (TF1*)h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->GetListOfFunctions()->FindObject("gaus");
					h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin]->GetFunction("gaus")->SetBit(TF1::kNotDraw); // do not
				}
				DrawHistoToCanvasPad(EtaPtResolutionPad[eta_bin], pT_bin+1, h_pTreco_ov_pTgen_JEScorrected[NoFile][eta_bin][pT_bin], Colors[NoFile], 1);

				if ( NoFile == NoFiles-1 )	leg3[eta_bin][pT_bin]->Draw("same");
				teta->Draw();



			}  //end of pT bin loops
		} // end of file loops

		EtaPtResolutionPad[eta_bin]->cd(pT_bins+1);	leg1[eta_bin]->Draw();

	} // end of eta bin loops

	ptResponsePad->cd(eta_bins+1); 
	leg2->Draw();

	ptResolutionPad->cd(eta_bins+1); 
	leg2->Draw();

	pTresTailPad->cd(eta_bins+1); 
	leg2->Draw();



	char filename[1500];
	if (Save_Plots)
	{
		sprintf(filename,"%s/%s/pTresolution_TailFrac_3p0_%s.png",analyzer_path,output_directory,image_name);
		pTresTailPad->SaveAs(filename);

		sprintf(filename,"%s/%s/pTresponse_%s.png",analyzer_path,output_directory,image_name);
		ptResponsePad->SaveAs(filename);

		sprintf(filename,"%s/%s/pTresolution_%s.png",analyzer_path,output_directory,image_name);
		ptResolutionPad->SaveAs(filename);

		for (int eta_bin = 0; eta_bin<eta_bins; eta_bin++) 
		{
			sprintf(filename,"%s/%s/pTresolution_pTbinned_etaBin%i_%s.png",analyzer_path,output_directory,eta_bin,image_name);
			EtaPtResolutionPad[eta_bin]->SaveAs(filename);
		}
	}



}


