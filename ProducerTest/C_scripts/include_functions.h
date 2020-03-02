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
#include "TH1F.h"
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

//#define eta_bins 4

using namespace std;



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

