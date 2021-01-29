//template <typename T,unsigned S>
//inline unsigned arraysize(const T (&v)[S]) { return S; }

#include <vector>
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<int> >+;
#endif
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
//#include "TH1F.h"
#include "TH2F.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TF1.h"
#include "TROOT.h"
#include "setTDRStyle_teliko.C"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <bits/stdc++.h> 
#include <sys/stat.h> 
#include <sys/types.h> 
#include <sys/time.h>

using namespace std;

// C++ program to create a directory in Linux 

  
void createDirectory(char dir_name[500]) 
  
{ 
  
    // Creating a directory 
    if (mkdir(dir_name, 0777) == -1) 
        cerr << "Warning: Tried to create directory :  " <<  dir_name << ".   "<< strerror(errno) << endl; 
  
    else
        cout << "Directory created: " << dir_name << endl; 
} 



/*
int getBin(double x, double boundaries[] ) 
{
	int i;

	//int n = sizeof(boundaries)/sizeof(boundaries[0])-1;
	int n = arraysize(boundaries) -1;
	cout << n << endl;
	cout <<    << endl; 

	if (x<boundaries[0] || x>=boundaries[n]) return -1;
	for(i=0; i<n; i++)
	{
		if (x>=boundaries[i] && x<boundaries[i+1]) return i;
	}
	return 0;
}
*/



int getBin(double x, double boundaries[], int bins) 
{
	int i;
	int n = bins; //sizeof(boundaries)/sizeof(boundaries[0])-1;
//	int n = sizeof(boundaries)/sizeof(boundaries[0])-1;
//	cout << n << endl;
	if (x<boundaries[0] || x>=boundaries[n]) return -1;
	for(i=0;i<n;i++)
	{
		if (x>=boundaries[i] && x<boundaries[i+1]) return i;
	}
	return 0;
}


double CalculateMu4(TH1D* histo)  //calculates E( (X-μ)^4 )
{
	double mean = histo->GetMean();
	double Ntot = histo->Integral();
	double Nbins = histo->GetNbinsX();

	double Mu4=0.; 

	for( int i=0; i<Nbins; i++ )
	{
		double BinWeight = histo->GetBinContent(i);
		double BinCenter = histo->GetBinCenter(i);

		Mu4 = Mu4 + BinWeight*pow( BinCenter-mean ,4);		
	}

	Mu4 = Mu4 / Ntot;

	return Mu4;
}

//error of standard deviation found here : https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation
double CalculateRMSerror(TH1D* histo)
{
	double RMS = histo->GetRMS();
	double Ntot = histo->Integral();
	double mu4 = CalculateMu4(histo);
	double RMSerror =  0.5*(1/RMS)*sqrt( (1/Ntot)*( mu4-(Ntot-3)*pow(RMS,4)/(Ntot-1) )  );
//	cout << "\n\n RMS = " <<  RMS << " + - " <<  RMSerror << endl;
	return RMSerror;
}


//=====================================================================
//The error bars are determined from Wilson interval of binomial errors. 
//define upper error
double error1 (double n1, double n2)
{
	if (n1 == 0)
	{
		return (0);
	}
	else
	{
	return ( (n2/n1 + 0.5/n1 + sqrt(n2/pow(n1,2)*(1-n2/n1) + 0.25/pow(n1,2)))/(1+1.0/n1) - n2/n1 );
	}
}
//define lower error
double error2 (double n1, double n2)
{
	if (n1 == 0)
	{
		return (0);
	}
	else
	{
		return ( n2/n1 - (n2/n1 + 0.5/n1 - sqrt(n2/pow(n1,2)*(1-n2/n1) + 0.25/pow(n1,2)))/(1+1.0/n1) );
	}
}





int * GetRecoToGenMatchSequence (std::vector<float> *gen_eta, std::vector<float> *gen_phi, std::vector<float> *reco_eta, std::vector<float> *reco_phi, double DR_threshold )
{

	const int gen_size = gen_eta->size();
	const int reco_size = reco_eta->size();
//	if (gen_size == 0 || reco_size == 0 ) return -1;

//	static int reco_jets_matched_sequence[gen_size];
	//const static int reco_jets_matched_sequence[gen_size];

	int * reco_jets_matched_sequence = NULL;    //define a pointer and initialize it
	reco_jets_matched_sequence = new int[gen_size];  //pointer to an int array ( define it like this to be of variable size)

//	double matched_minDR[gen_size];

	for ( int j=0; j<gen_size; j++) 
	{
		double minDR = 100;
		int matched_index = -1;
		reco_jets_matched_sequence[j] = -1; //initilize the array. 
//		matched_minDR[j] = -1;

		for(int k=0; k<reco_size; k++)
		{
			double deltaR = sqrt( pow( gen_eta->at(j) - reco_eta->at(k), 2) + pow( gen_phi->at(j) - reco_phi->at(k), 2) );
		//	cout << "gen jet "<< j << "   reco jet " << k <<"    deltaR = " << deltaR << endl;
			if ( minDR > deltaR )
			{
				minDR = deltaR;
				matched_index = k;  
		//		cout << "minDR = " << minDR << ",   matched indexes: " << j << " - " << matched_index << endl;
			}		
		}

//			matched_minDR[j] = minDR;
			if (minDR < DR_threshold)
			{
				reco_jets_matched_sequence[j] = matched_index;
	//			matched_minDR[j] = minDR;
			}
//		cout << "At function, gen jet " << j  <<"matched with reco jet = " << reco_jets_matched_sequence[j] <<endl;
		}

	return reco_jets_matched_sequence;
}


std::vector<int>  GetRecoUnmatchedSequence ( int reco_size, int gen_size, int *reco_jets_matched_sequence )
{
	std::vector<int> unmatched_sequence ;

	for (int l=0; l< reco_size; l++)
	{
		bool matched = false;
		for (int k=0; k<gen_size; k++) if (l == reco_jets_matched_sequence[k] ) { matched = true; break; }
		if (!matched) unmatched_sequence.push_back(l);
	}

	return unmatched_sequence;
}


double * GetMatchedMinDR (std::vector<float> *gen_eta, std::vector<float> *gen_phi, std::vector<float> *reco_eta, std::vector<float> *reco_phi, int *reco_jets_matched_sequence )
{

	const int gen_size = gen_eta->size();
	double * matched_minDR = NULL;    //define a pointer and initialize it
	matched_minDR = new double[gen_size];  //pointer is an int array ( define it like this to be of variable size)




	for (int l=0; l<gen_size; l++) // get the
	{
		double minDR;
		if (reco_jets_matched_sequence[l]<0) matched_minDR[l] = -1.; // if no match was found set value to -1
		else 
		{
			minDR = sqrt( pow( gen_eta->at(l) - reco_eta->at(reco_jets_matched_sequence[l]), 2) + pow( gen_phi->at(l) - reco_phi->at(reco_jets_matched_sequence[l]),2) );
			matched_minDR[l] = minDR;
		}
	}

	return matched_minDR;
}

TF1* FindBestGaussianCoreFit(TH1D* histo)
{

	double mean = histo->GetMean();
	double rms = histo->GetRMS();
	bool quiet_mode = true; //minimum prints on terminal
	
	int HalfMaxBinLow = histo->FindFirstBinAbove(histo->GetMaximum()/2);
	int HalfMaxBinHigh = histo->FindLastBinAbove(histo->GetMaximum()/2);
	double WidthAtHalfMaximum = 0.5*(histo->GetBinCenter(HalfMaxBinHigh) - histo->GetBinCenter(HalfMaxBinLow));
    	double Xmax  = histo->GetXaxis()->GetBinCenter(histo->GetMaximumBin());

	TF1 *gausTF1;

	double Pvalue = 0;
	double ChiSquare;
	double ndf;

	double rms_step_plus;	
	double RangeLow = 0.;
	double RangeUp = 5.;
	double meanForRange;
	double spreadForRange;

	double PvalueBest = 0;
	double RangeLowBest = 0;
	double RangeUpBest = 0;
	double ChiSquareBest;
	double ndfBest;
	double StepMinusBest;
	double StepPlusBest;

	if(Xmax < mean) meanForRange = Xmax;
	else meanForRange = mean; //because some entries with LARGE errors sometimes falsely become the maximum

	if(WidthAtHalfMaximum < rms && WidthAtHalfMaximum>0) spreadForRange = WidthAtHalfMaximum;
	else spreadForRange = rms; //because WHF does not take into account weights and sometimes it turns LARGE

	double rms_step_minus = 2.2;
	while (rms_step_minus>1.1)
	{ 
		RangeLow = meanForRange - rms_step_minus*spreadForRange;
//		rms_step_plus = 2.2;
		rms_step_plus = rms_step_minus;  

		while ( rms_step_plus>0.7  )
		{	
			RangeUp = meanForRange + rms_step_plus*spreadForRange;
			if(quiet_mode)	histo->Fit("gaus","0Q","0",RangeLow, RangeUp);
			else 		histo->Fit("gaus","0","0",RangeLow, RangeUp);
			gausTF1 = (TF1*)	histo->GetListOfFunctions()->FindObject("gaus");
			ChiSquare = gausTF1->GetChisquare();
			ndf       = gausTF1->GetNDF();
			Pvalue = TMath::Prob(ChiSquare, ndf);			
	
			if (Pvalue > PvalueBest)
			{
				PvalueBest = Pvalue;
				RangeLowBest = RangeLow;
				RangeUpBest = RangeUp;
				ndfBest = ndf;
				ChiSquareBest = ChiSquare;
				StepMinusBest = rms_step_minus; 
				StepPlusBest = rms_step_plus; 
				meanForRange = gausTF1->GetParameter(1);
			}

			if(!quiet_mode)
			{
				cout << "\n\nFitting range used: [Mean - " << rms_step_minus  << " sigma,  Mean + " << rms_step_plus << " sigma ] "<< endl;
				cout << "ChiSquare = " << ChiSquare << " NDF = " << ndf << " Prob =  " << Pvalue << "  Best Prob so far = " << PvalueBest << endl;
			}
			rms_step_plus = rms_step_plus - 0.1;
		}
		rms_step_minus = rms_step_minus - 0.1;
	}

	if (quiet_mode) histo->Fit("gaus","0Q","0",RangeLowBest, RangeUpBest);
	else 
	{
		histo->Fit("gaus","0","0",RangeLowBest, RangeUpBest);
		cout << "\n\n\nMean =  " << mean << "  Xmax = " << Xmax << "  RMS = " << rms << "  WidthAtHalfMax = " << WidthAtHalfMaximum <<  endl;
		cout << "Fit found!" << endl;
		cout << "Final fitting range used: [Mean(Xmax) - " << StepMinusBest << " rms(WHF), Mean(Xmax) + " << StepPlusBest << " rms(WHF) ] "<< endl;
		cout << "ChiSquare = " << ChiSquareBest << " NDF = " << ndfBest << " Prob =  " << PvalueBest << "\n\n" << endl;
	}	
	gausTF1 = (TF1*)	histo->GetListOfFunctions()->FindObject("gaus");

	return gausTF1;
}




TGraphAsymmErrors* GetEfficiencyGraph(TH1D* numerator, TH1D* denominator )
{
	const int nbins = denominator->GetNbinsX(); 
	double eff_x[nbins], eff_y[nbins], eff_eyu[nbins], eff_eyl[nbins], eff_ex[nbins];

	for (int i=0; i<nbins; i++)
	{
		eff_x[i] = 0.0; eff_y[i] = 0.0; eff_eyu[i] = 0.0; eff_eyl[i] = 0.0; eff_ex[i] = 0.0;

		if( denominator->GetBinContent(i)>0 )
		{
			eff_x[i]   = denominator->GetBinCenter(i);
			eff_y[i]   = numerator->GetBinContent(i) / denominator->GetBinContent(i);
			eff_eyu[i] = error1(denominator->GetBinContent(i), numerator->GetBinContent(i) );
			eff_eyl[i] = error2(denominator->GetBinContent(i), numerator->GetBinContent(i) );
			eff_ex[i]  = denominator->GetBinWidth(i)/2;
		}
	}
	TGraphAsymmErrors *Eff = new TGraphAsymmErrors( nbins, eff_x, eff_y,eff_ex, eff_ex, eff_eyl, eff_eyu );	

	return Eff;
}

TH1D* InitiateFrameOnCanvasPad(TCanvas* c,int padNo, const char frameName[50], const char TitleX[50], const char TitleY[50], double Xmin, double Xmax, double Ymin, double Ymax, bool logYOn, TPaveText *pave)
{

	char name[100];
	sprintf(name,"%s_%d",frameName,padNo);
	TH1D *frame = new TH1D(name,name,100, Xmin, Xmax);

	frame->SetMinimum(Ymin);
	frame->SetMaximum(Ymax);
	frame->GetXaxis()->SetRangeUser(Xmin,Xmax);
	frame->GetXaxis()->SetTitle(TitleX);
	frame->GetYaxis()->SetTitle(TitleY);
	frame->SetLineColor(0);
	frame->SetMarkerColor(0);
	frame->GetXaxis()->SetTitleOffset(1.3);
	frame->GetYaxis()->SetTitleOffset(1.4);

	c->cd(padNo);
	if(logYOn) c->cd(padNo)->SetLogy(1);
	frame->Draw("");
	//pave->Draw("same");

	return frame;
}

void DrawHistoToCanvasPad(TCanvas* c, int padNo, TH1D *hist, int ColorNo, int LineStyleNo)
{
	hist->SetLineColor(ColorNo);	
	hist->SetLineStyle(LineStyleNo);	
	c->cd(padNo);
	hist->Draw("same hist");
}

void DrawTH1FToCanvasPad(TCanvas* c, int padNo, TH1F *hist, int ColorNo, int LineStyleNo)
{
	hist->SetLineColor(ColorNo);	
	hist->SetLineStyle(LineStyleNo);	
	c->cd(padNo);
	hist->Draw("same hist");
}

void DrawTH2FToCanvasPad(TCanvas* c, int padNo, TH2F *hist, const char TitleX[50], const char TitleY[50], const char TitlePad[50], double Xmin, double Xmax, double Ymin, double Ymax, bool logXOn, bool logYOn,bool logZOn, TPaveText *pave)
{
		c->cd(padNo);
		if(logZOn)c->cd(padNo)->SetLogz(1);
		if(logYOn)c->cd(padNo)->SetLogy(1);
		if(logXOn)c->cd(padNo)->SetLogx(1);
		c->cd(padNo)->SetRightMargin(0.18);
		c->cd(padNo)->SetTopMargin(0.07);					
		c->SetTitle(TitlePad);

		hist->SetStats(true);
		hist->GetXaxis()->SetRangeUser(Xmin,Xmax);
		hist->GetYaxis()->SetRangeUser(Ymin,Ymax);
		hist->GetXaxis()->SetTitle(TitleX);
		hist->GetYaxis()->SetTitle(TitleY);
		hist->GetYaxis()->SetTitleOffset(1.3);
		hist->SetMinimum(0.9 );
		hist->SetMaximum(5*hist->GetMaximum());

		hist->Draw("colz");
		//pave ->Draw("same");

}

void DrawGraphToCanvasPad(TCanvas* c, int padNo, TGraphAsymmErrors *graph, int ColorNo, int LineStyleNo)
{
	graph->SetLineColor(ColorNo);
	graph->SetMarkerStyle(24);	
	graph->SetMarkerColor(ColorNo);
	graph->SetLineStyle(LineStyleNo);
	graph->SetMarkerSize(0.3);	
	c->cd(padNo);
	graph->Draw("same p");
}


