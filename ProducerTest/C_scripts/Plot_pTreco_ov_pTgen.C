//==========================================
// Initial Author: Dimitrios Karasavvas
// Date:   29 Jan 2021
//==========================================

#include "include/include_functions.h"
#include "include/Input_definition.h"



void Plot_pTreco_ov_pTgen(bool JEScorrected = false)
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
	gStyle->SetOptStat(0);
	gROOT->ForceStyle(kTRUE);
	gStyle->SetOptFit(0);


	char input_file[500];
	if(!JEScorrected) sprintf(input_file, "%s/%s/pTreco_ov_pTgen_%s_histos.root",analyzer_path,output_directory,image_name );
	else sprintf(input_file, "%s/%s/pTreco_ov_pTgen_%s_JEScorrected_histos.root", analyzer_path,output_directory,image_name );
	TFile *f_input = new TFile (input_file,"READ");
	


/*
	bool JEScorrected;
	char *output = NULL; 
	char JEScorrChar[] = "JEScorrected";
	output = strstr (input_file, JEScorrChar);
	if(output)	JEScorrected = true; // search for "JEScorrected" in input file's name
	else		JEScorrected = false;
*/

	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;
	const int pT_bins = sizeof(ptBnd)/sizeof(ptBnd[0])-1;
    	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);    
	
	cout << "Availabe legends in current file: " << endl;   
	for(int NoFile=0; NoFile<NoFiles; NoFile++)		cout << "For " << legend_array[NoFile] << ", type " << NoFile<< endl;

	cout << "Pick an integer: " << endl;
	int chosenLegend;
//	cin >> chosenLegend;


	while (!(cin >> chosenLegend) || chosenLegend < 0 || chosenLegend >NoFiles-1)
	{
		cin.clear();//to clear the buffer memory
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cout << " Invalid input.! What you typed does not correspond to any of the available legends. Pick an integer between 0 and " << NoFiles-1 << ". "<< endl;
		cout << "Availabe legends in current file: " << endl;   
		for(int NoFile=0; NoFile<NoFiles; NoFile++)		cout << "For " << legend_array[NoFile] << ", type " << NoFile<< endl;
	//	cin >> chosenLegend;
	}
	
	cout << "Will plot the pTreco/pTgen distributions for " << legend_array[chosenLegend]<< endl;


	char legend[500];
	char imageName[100];
	sprintf(legend, "%s", legend_array[chosenLegend]);
	sprintf(imageName, "%s", legend);	

	if (JEScorrected) strcat(imageName, "_JEScorrected");
	else strcat(imageName,"_raw"); 

	double pTraw_responses[eta_bins][pT_bins], pTraw_resolution[eta_bins][pT_bins];
	double Xmax[eta_bins][pT_bins], WidthHalfMax[eta_bins][pT_bins];

	double pTraw_responses_error[eta_bins][pT_bins], pTraw_resolution_error[eta_bins][pT_bins];

	double pTraw_resolution_Gaus[eta_bins][pT_bins], pTraw_response_Gaus[eta_bins][pT_bins];
	double pTraw_resolution_Gaus_error[eta_bins][pT_bins], pTraw_response_Gaus_error[eta_bins][pT_bins];

	double GausFits_ChiSq[eta_bins][pT_bins], GausFits_prob[eta_bins][pT_bins];
	int    GausFits_NDF[eta_bins][pT_bins];


	double pT_resolution_corrected[eta_bins][pT_bins], pT_resolution_corrected_Gaus[eta_bins][pT_bins];
	double pTreco_ov_pTgen_JEScorrected_error[eta_bins][pT_bins], pTreco_ov_pTgen_JEScorrected_Gaus_error[eta_bins][pT_bins];



	double pTbinCenter[pT_bins], pTbinCenter_error[pT_bins];

	

	TH1D *h_pTreco_ov_pTgen[eta_bins][pT_bins];

//	cout << "eta bins = " << eta_bins;
//	cout << "pT bins = " << pT_bins;


	char name[1500]; 
	



	for (int  eta_bin=0; eta_bin<eta_bins; eta_bin++)
	{
		for (int  pT_bin=0; pT_bin<pT_bins; pT_bin++)
		{
			if(JEScorrected)
			{
				sprintf(name,"h_pTreco_ov_pTgen_JEScorrected_%s_EtaBin%i_ptBin%i",legend,eta_bin,pT_bin);
				h_pTreco_ov_pTgen[eta_bin][pT_bin] = (TH1D*)(f_input->Get(name));
				//h_pTreco_ov_pTgen[eta_bin][pT_bin]->Rebin(2);
			}
			else
			{
				sprintf(name,"h_pTreco_ov_pTgen_%s_EtaBin%i_ptBin%i",legend,eta_bin,pT_bin);
				h_pTreco_ov_pTgen[eta_bin][pT_bin] = (TH1D*)(f_input->Get(name));
				//h_pTreco_ov_pTgen[eta_bin][pT_bin]->Rebin(2);
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


	TF1 *gausTF1;

	cout << " fitting pTreco / pTgen histos for file " << input_file << " ...  "<< endl;

	for ( int eta_bin=0; eta_bin< eta_bins; eta_bin++ )
	{
		for (int pT_bin=0; pT_bin< pT_bins; pT_bin++ )
		{
			cout << " Now processing: Eta bin: " <<  eta_bins_legend[eta_bin] << ",  pT bin: " << pT_bins_legend[pT_bin] << endl;
			
			double Nsize = h_pTreco_ov_pTgen[eta_bin][pT_bin]->Integral();
			if ( !(Nsize>0) ) 
			{
				cout << "Warning!!! Histogram empty!!  eta bin : " << eta_bins_legend[eta_bin] << "  pT : " << pT_bins_legend[pT_bin] << endl;
				continue;
			}
			double res_mean = h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetMean();
			double res_rms = h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetRMS(); //
			double res_mean_error = res_rms / Nsize; // error of mean value.
			double res_rms_error = 	CalculateRMSerror( h_pTreco_ov_pTgen[eta_bin][pT_bin] ); //error of standard deviation found here : https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation
			int HalfMaxBinLow = h_pTreco_ov_pTgen[eta_bin][pT_bin]->FindFirstBinAbove(h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetMaximum()/2);
			int HalfMaxBinHigh = h_pTreco_ov_pTgen[eta_bin][pT_bin]->FindLastBinAbove(h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetMaximum()/2);
			double WidthAtHalfMaximum = 0.5*(h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetBinCenter(HalfMaxBinHigh) - h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetBinCenter(HalfMaxBinLow));
    		double xmax  = h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetXaxis()->GetBinCenter(h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetMaximumBin());


			pTraw_responses[eta_bin][pT_bin] = res_mean;
			pTraw_resolution[eta_bin][pT_bin] = res_rms;
			pTraw_responses_error[eta_bin][pT_bin] = res_mean_error;
			pTraw_resolution_error[eta_bin][pT_bin] = res_rms_error;
			Xmax[eta_bin][pT_bin] = xmax;
			WidthHalfMax[eta_bin][pT_bin] = WidthAtHalfMaximum;

			if (res_rms>0 && h_pTreco_ov_pTgen[eta_bin][pT_bin]->Integral()>0 )
			{
				//if( pTbinCenter[pT_bin] < maxPtForGausReducedRangeFit) h_pTreco_ov_pTgen[eta_bin][pT_bin]->Fit("gaus","0","0",res_mean-2.2*res_rms,res_mean+0.8*res_rms);
				//else h_pTreco_ov_pTgen[eta_bin][pT_bin]->Fit("gaus","0","0",res_mean-2.2*res_rms,res_mean+1.8*res_rms);
				//TF1 *gausTF1 = (TF1*)	h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetListOfFunctions()->FindObject("gaus");
				gausTF1 = FindBestGaussianCoreFit(h_pTreco_ov_pTgen[eta_bin][pT_bin]);
				GausFits_ChiSq[eta_bin][pT_bin] = gausTF1->GetChisquare();
				GausFits_NDF[eta_bin][pT_bin]   = gausTF1->GetNDF();
				GausFits_prob[eta_bin][pT_bin]  = TMath::Prob( gausTF1->GetChisquare() , gausTF1->GetNDF() );


				pTraw_response_Gaus[eta_bin][pT_bin] = gausTF1->GetParameter(1);
				pTraw_response_Gaus_error[eta_bin][pT_bin] = gausTF1->GetParError(1);

				pTraw_resolution_Gaus[eta_bin][pT_bin] = fabs(gausTF1->GetParameter(2));
				pTraw_resolution_Gaus_error[eta_bin][pT_bin] = gausTF1->GetParError(2);
			}
			else continue;
				

		} // end of loop of pT bins
	}// end of loop of eta bins


// ======================================================plotting stuff======================================//
	char res_text[400];
	



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


		

			const char *etaText = (yBnd[eta_bin]==0 ? Form("|#eta| < %1.2g",yBnd[eta_bin+1]) :
			Form("%1.2g #leq |#eta| < %1.2g",yBnd[eta_bin],yBnd[eta_bin+1]));
			TLatex *etaTextLatx = new TLatex(0.28,0.86,etaText); //cout<<seta<<endl;
			etaTextLatx->SetNDC();
			etaTextLatx->SetTextSize(0.06);


			for (int pT_bin = 0; pT_bin<pT_bins; pT_bin++)
			{

				const char *seta = (ptBnd[pT_bin]==0 ? Form("|y| < %3.0f",ptBnd[pT_bin+1]) :
				Form("%3.0f #leq pT < %3.0f",ptBnd[pT_bin],ptBnd[pT_bin+1]));
				TLatex *teta = new TLatex(0.28,0.86,seta); //cout<<seta<<endl;
				teta->SetNDC();
				teta->SetTextSize(0.06);

				
				leg3[eta_bin][pT_bin] =new TLegend(.2, .7, .8, .85);//7899//4899
				leg3[eta_bin][pT_bin]->SetTextSize(0.04);
				leg3[eta_bin][pT_bin]->SetFillColor(0); 
				leg3[eta_bin][pT_bin]->SetBorderSize(0);
		 		
			
				gausTF1 = (TF1*)h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetListOfFunctions()->FindObject("gaus");
				
				sprintf(res_text, "Mean = %3.3f , RMS = %3.3f ", pTraw_responses[eta_bin][pT_bin] ,pTraw_resolution[eta_bin][pT_bin]);
				leg3[eta_bin][pT_bin]->AddEntry(h_pTreco_ov_pTgen[eta_bin][pT_bin], res_text , "p");
				sprintf(res_text, "Xmax = %3.3f , WHF = %3.3f ", Xmax[eta_bin][pT_bin] ,WidthHalfMax[eta_bin][pT_bin]);
				leg3[eta_bin][pT_bin]->AddEntry((TObject*)0, res_text , "");
				sprintf(res_text, "Mean = %3.3f , Sigma = %3.3f ", pTraw_response_Gaus[eta_bin][pT_bin] , pTraw_resolution_Gaus[eta_bin][pT_bin]);
				leg3[eta_bin][pT_bin]->AddEntry(gausTF1, res_text , "L");
				sprintf(res_text, "#chi^{2}= %3.1f, ndf= %i, Prob= %3.3f",GausFits_ChiSq[eta_bin][pT_bin], GausFits_NDF[eta_bin][pT_bin] ,GausFits_prob[eta_bin][pT_bin] );
				leg3[eta_bin][pT_bin]->AddEntry( (TObject*)0, res_text, "");

				EtaPtResolutionPad[eta_bin]->cd(pT_bin+1);
				//EtaPtResolutionPad[eta_bin]->cd(pT_bin+1)->SetLogy(1);
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetXaxis()->SetTitle("pT_reco / pT_gen");
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetXaxis()->SetRangeUser( 0., 2.);
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetYaxis()->SetTitle("Entries");
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetYaxis()->SetTitleOffset(1.3);

				h_pTreco_ov_pTgen[eta_bin][pT_bin]->SetLineColor(1); 
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->SetMarkerSize(0.45); 
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->SetLineWidth(1.05); 
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->SetLineStyle(1);

//				h_pTreco_ov_pTgen[eta_bin][pT_bin]->SetMinimum(0.01*h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetMaximum());
//				h_pTreco_ov_pTgen[eta_bin][pT_bin]->SetMaximum(20*h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetMaximum());
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->SetMinimum(0.05*h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetMaximum());
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->SetMaximum(1.5*h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetMaximum());
							
				h_pTreco_ov_pTgen[eta_bin][pT_bin]->Draw();
				if( h_pTreco_ov_pTgen[eta_bin][pT_bin]->Integral()>0 && h_pTreco_ov_pTgen[eta_bin][pT_bin]->GetRMS()>0 ) gausTF1->Draw("same");
				paveCMS ->Draw("same");  
				leg3[eta_bin][pT_bin]->Draw("same");
				teta->Draw();

			}  //end of pT bin loops
		EtaPtResolutionPad[eta_bin]->cd(pT_bins+1);	leg1[eta_bin]->Draw();
	} // end of eta bin loops

	char filename[1500];
	if (Save_Plots)
	{
		char GausFitDir[500];
		sprintf(GausFitDir,"%s/%s/GausFits",analyzer_path,output_directory );
		createDirectory(GausFitDir);
		sprintf(GausFitDir,"%s/%s/GausFits/%s",analyzer_path,output_directory,legend );
		createDirectory(GausFitDir);


		for (int eta_bin = 0; eta_bin<eta_bins; eta_bin++) 
		{
			sprintf(filename,"%s/%s/GausFits/%s/pTreco_ov_pTgen_pTbinned_etaBin%i_%s.png",analyzer_path,output_directory,legend,eta_bin,imageName);
			EtaPtResolutionPad[eta_bin]->SaveAs(filename);
		}
	}



}


