#include "include_functions.h"


void Reco_gen_matching_Comparisons_function()
{
	bool Save_Plots = true;
//	TFile *f1 = TFile::Open("FullTracking_Scouting_dijet_pre13.root","READ"); 
//	TFile *f1 = TFile::Open("PixelTracking_PatatrackPixelTracks_Cleaned.root","READ"); //full tracking file
//	TFile *f1 = TFile::Open("PixelTracking_Scouting_dijet_pre13.root","READ"); 
	TFile *f1 = TFile::Open("/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/RelVal_ttbar/FullTracking_ttbar_14TeV_noPU.root","READ"); 

	char analyzer_path[600] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre7/src/Jakobs_producer/ProducerTest/"}; 
//	char *output_directory = "deleteme/";     
	char output_directory[300] = {"Realistic_RunIII_14TeV/ttbar/Reco_vs_gen/FullTracking/"};     
	char image_name[200] = {"Comparison_reco_vs_gen_ttbar_DR0p2_noPU"}; 


	const int eta_bins = 5;

	TH1D *h_METovSUMET[eta_bins], *h_CHFJet[eta_bins], *h_NHFJet[eta_bins], *h_CEMFJet[eta_bins], *h_NEMFJet[eta_bins], *h_MUFJet[eta_bins], *h_CMJet[eta_bins], *h_NMJet[eta_bins], *h_ptJet[eta_bins], *h_PHIJet[eta_bins];

	TH1D *h_METovSUMET_gen[eta_bins], *h_CHFJet_gen[eta_bins], *h_NHFJet_gen[eta_bins], *h_CEMFJet_gen[eta_bins], *h_NEMFJet_gen[eta_bins], *h_MUFJet_gen[eta_bins], *h_CMJet_gen[eta_bins], *h_NMJet_gen[eta_bins], *h_ptJet_gen[eta_bins], *h_PHIJet_gen[eta_bins];

	TH1D *h_ETAJet, *h_ETAJet_gen, *h_minDR, *h_SumEt, *h_SumEt_gen;

	TH1D *h_pT_resolution[eta_bins], *h_reco_efficiency[eta_bins], *h_gen_pT_all[eta_bins];



	double DR_threshold = 0.2;

	char filename[256]; 
	char name[256]; 

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko();


	sprintf(name,"h_ETAJet");
	h_ETAJet = new TH1D(name, "", 60,-5,5);

	sprintf(name,"h_ETAJet_gen");
	h_ETAJet_gen = new TH1D(name, "", 60,-5,5);

	sprintf(name,"h_minDR");
	h_minDR = new TH1D(name, "", 120,0,0.5);

	sprintf(name,"h_SumEt");
	h_SumEt = new TH1D(name, "", 100, 0, 3000); // 40,0,1.0

	sprintf(name,"h_SumEt_gen");
	h_SumEt_gen = new TH1D(name, "", 100, 0, 3000); // 40,0,1.0

	for(Int_t h=0; h <eta_bins; h++)
	{ 
	//========== reco jets===========
		sprintf(name,"h_METovSUMET%i",h);
		h_METovSUMET[h] = new TH1D(name, "", 50, 0, 1); // 40,0,1.0

		sprintf(name,"h_CHFJet%i",h);
		h_CHFJet[h] = new TH1D(name, "", 60, 0, 1.2); //40, 0,1.2

		sprintf(name,"h_NHFJet%i",h);
		h_NHFJet[h] = new TH1D(name, "", 60, 0, 1.2);

		sprintf(name,"h_CEMFJet%i",h);
		h_CEMFJet[h] = new TH1D(name, "", 60, 0, 1.2);

		sprintf(name,"h_NEMFJet%i",h);
		h_NEMFJet[h] = new TH1D(name, "", 60, 0, 1.2);

		sprintf(name,"h_MUFJet%i",h);
		h_MUFJet[h] = new TH1D(name, "", 60, 0, 1.2);

		sprintf(name,"h_CMJet%i",h);
		h_CMJet[h] = new TH1D(name, "", 80, 0, 80);

		sprintf(name,"h_PHIJet%i",h);
		h_PHIJet[h] = new TH1D(name, "", 40,-4,4);

		sprintf(name,"h_NMJet%i",h);
		h_NMJet[h] = new TH1D(name, "", 80, 0, 80);

		sprintf(name,"h_ptJet%i",h);
		h_ptJet[h] = new TH1D(name, "", 30,0,1000);

		//========== gen jets===========
		sprintf(name,"h_METovSUMET_gen%i",h);
		h_METovSUMET_gen[h] = new TH1D(name, "", 50, 0, 1); // 40,0,1.0

		sprintf(name,"h_CHFJet_gen%i",h);
		h_CHFJet_gen[h] = new TH1D(name, "", 60, 0, 1.2); //40, 0,1.2
		sprintf(name,"h_NHFJet_gen%i",h);
		h_NHFJet_gen[h] = new TH1D(name, "", 60, 0, 1.2);
		sprintf(name,"h_CEMFJet_gen%i",h);
		h_CEMFJet_gen[h] = new TH1D(name, "", 60, 0, 1.2);
		sprintf(name,"h_NEMFJet_gen%i",h);
		h_NEMFJet_gen[h] = new TH1D(name, "", 60, 0, 1.2);
		sprintf(name,"h_MUFJet_gen%i",h);
		h_MUFJet_gen[h] = new TH1D(name, "", 60, 0, 1.2);

		sprintf(name,"h_CMJet_gen%i",h);
		h_CMJet_gen[h] = new TH1D(name, "", 80, 0, 80);
		
		sprintf(name,"h_PHIJet_gen%i",h);
		h_PHIJet_gen[h] = new TH1D(name, "", 40,-4,4);

		sprintf(name,"h_NMJet_gen%i",h);
		h_NMJet_gen[h] = new TH1D(name, "", 80, 0, 80);

		sprintf(name,"h_ptJet_gen%i",h);
		h_ptJet_gen[h] = new TH1D(name, "", 30,0,1000);

		sprintf(name,"h_pT_resolution%i",h);
		h_pT_resolution[h] = new TH1D(name, "", 50, 0., 2.); // 40,0,1.0

		sprintf(name,"h_reco_efficiency%i",h);
		h_reco_efficiency[h] = new TH1D(name, "", 30, 0., 1000); // 40,0,1.0
		h_reco_efficiency[h]->Sumw2();

		sprintf(name,"h_gen_pT_all%i",h);
		h_gen_pT_all[h] = new TH1D(name, "", 30, 0., 1000); // 40,0,1.0

	}


	double yBnd[eta_bins+1]={0.0,1.3,2.4,2.7, 3.0, 5.0}; 

 
	TTree *tree = (TTree*)f1->Get("dijetscouting/tree");


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


 	double SumEt, gen_SumEt;
   
 
 

    tree->SetBranchAddress("jpt",&jpt);
    tree->SetBranchAddress("eta",&eta);
    tree->SetBranchAddress("phi",&phi);
    tree->SetBranchAddress("nhf",&nhf);
    tree->SetBranchAddress("chf",&chf);
    tree->SetBranchAddress("cemf",&cemf);
    tree->SetBranchAddress("nemf",&nemf);
    tree->SetBranchAddress("npr",&npr);
    tree->SetBranchAddress("muf",&muf);
    tree->SetBranchAddress("chMult",&chMult);
    tree->SetBranchAddress("neMult",&neMult);
    tree->SetBranchAddress("SumEt",&SumEt);


    tree->SetBranchAddress("gen_jpt",&gen_jpt);
    tree->SetBranchAddress("gen_eta",&gen_eta);
    tree->SetBranchAddress("gen_phi",&gen_phi);
    tree->SetBranchAddress("gen_nhf",&gen_nhf);
    tree->SetBranchAddress("gen_chf",&gen_chf);
    tree->SetBranchAddress("gen_cemf",&gen_cemf);
    tree->SetBranchAddress("gen_nemf",&gen_nemf);
    tree->SetBranchAddress("gen_npr",&gen_npr);
    tree->SetBranchAddress("gen_muf",&gen_muf);
    tree->SetBranchAddress("gen_chMult",&gen_chMult);
    tree->SetBranchAddress("gen_neMult",&gen_neMult);
    tree->SetBranchAddress("gen_SumEt",&gen_SumEt);


	int nentries = tree->GetEntries(); 

	cout << " Number of Entries:  " << nentries << endl;

	for (int i=0; i<nentries; i++) //event loop
//	for (int i=0; i<60; i++) //event loop
	{

		tree->GetEntry(i);

	//======================= gen jets ===================================
		const int gen_size = gen_jpt->size();
		const int reco_size = jpt->size();
		if (gen_size ==0 || reco_size == 0 ) continue;


		int *reco_jets_matched_sequence;
		double *matched_minDR;

		reco_jets_matched_sequence = GetRecoToGenMatchSequence(gen_eta, gen_phi, eta, phi, DR_threshold);

		matched_minDR = GetMatchedMinDR(gen_eta, gen_phi, eta, phi, reco_jets_matched_sequence);



		for ( int j=0; j<gen_size; j++) 
		{

			int gen_ybin = getBin(fabs(gen_eta->at(j)),yBnd, eta_bins);
			if (gen_ybin > -1) h_gen_pT_all[gen_ybin]->Fill( gen_jpt->at(j) );


			if (reco_jets_matched_sequence[j] < 0 ) continue; // if gen jet has no match continue;
			int reco_ybin = getBin( fabs( eta->at( reco_jets_matched_sequence[j] ) ),yBnd, eta_bins );
//			cout << "gen_ybin = " << gen_ybin <<"   reco_ybin = " << reco_ybin << endl;

			h_minDR->Fill(matched_minDR[j]);
			h_ETAJet_gen->Fill( gen_eta->at(j) );


			if (j==0 ) h_SumEt_gen->Fill(gen_SumEt);

			if (gen_ybin > -1)//fill hist's in the corresponding eta bin
			{
				h_CHFJet_gen[gen_ybin]->Fill( gen_chf->at(j) ); 
				h_NHFJet_gen[gen_ybin]->Fill( gen_nhf->at(j) );
				h_CEMFJet_gen[gen_ybin]->Fill( gen_cemf->at(j) );
				h_NEMFJet_gen[gen_ybin]->Fill( gen_nemf->at(j) );
				h_MUFJet_gen[gen_ybin]->Fill( gen_muf->at(j) );
				h_PHIJet_gen[gen_ybin]->Fill( gen_phi->at(j) );
				h_CMJet_gen[gen_ybin]->Fill( gen_chMult->at(j) );
				h_NMJet_gen[gen_ybin]->Fill( gen_neMult->at(j) );
				h_ptJet_gen[gen_ybin]->Fill( gen_jpt->at(j) );

				h_pT_resolution[gen_ybin]->Fill(  jpt->at( reco_jets_matched_sequence[j] )  / gen_jpt->at(j) );
			} 

			h_ETAJet->Fill( eta->at(reco_jets_matched_sequence[j]) );
			if (j==0 ) h_SumEt->Fill(SumEt);

		//	if(reco_jets_matched_sequence[j]>-0.1)cout << " pT gen = " << gen_jpt->at(j) <<  " pT reco matched = " << jpt->at(reco_jets_matched_sequence[j]) << endl; 

			if (reco_ybin > -1)//fill hist's in the corresponding eta bin
			{
				h_CHFJet[reco_ybin]->Fill( chf->at(reco_jets_matched_sequence[j]) ); 
				h_NHFJet[reco_ybin]->Fill( nhf->at(reco_jets_matched_sequence[j]) );
				h_CEMFJet[reco_ybin]->Fill( cemf->at(reco_jets_matched_sequence[j]) );
				h_NEMFJet[reco_ybin]->Fill( nemf->at(reco_jets_matched_sequence[j]) );
				h_MUFJet[reco_ybin]->Fill( muf->at(reco_jets_matched_sequence[j]) );
				h_PHIJet[reco_ybin]->Fill( phi->at(reco_jets_matched_sequence[j]) );
				h_CMJet[reco_ybin]->Fill( chMult->at(reco_jets_matched_sequence[j]) );
				h_NMJet[reco_ybin]->Fill( neMult->at(reco_jets_matched_sequence[j]) );
				h_ptJet[reco_ybin]->Fill( jpt->at(reco_jets_matched_sequence[j]) );
			}
			
		}	



	delete reco_jets_matched_sequence;
	delete matched_minDR;
	}

	for(Int_t h=0; h <eta_bins; h++)
	{ 

		cout << " reco pt = " << h_ptJet[h]->Integral() << "  gen pT all = " << h_gen_pT_all[h]->Integral() << "  matched gen pT = "<< h_ptJet_gen[h]->Integral() << endl;
		h_reco_efficiency[h]->Divide(h_ptJet_gen[h] , h_gen_pT_all[h] );
	}

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
	TCanvas *pad7 = new TCanvas("pad7", "",1);
	pad7->Divide(3,2);
	TCanvas *pad8 = new TCanvas("pad8", "",1);
	pad8->Divide(3,2);
	TCanvas *pad9 = new TCanvas("pad9", "",1);
	pad9->Divide(3,2);
	TCanvas *pad10 = new TCanvas("pad10", "",1);
	pad10->Divide(3,2);
	TCanvas *pad11 = new TCanvas("pad11", "",1);
	pad11->Divide(3,2);

	TCanvas *c_eta = new TCanvas("c_eta", "",1);
	TCanvas *c_SumEt = new TCanvas("c_SumEt", "",1);
	TCanvas *c_minDR = new TCanvas("c_minDR", "",1);


	TPaveText *paveCMS = new TPaveText(0.45,0.95,0.5,1.0,"NDC");
	// paveCMS->AddText("CMS Preliminary L=9.2 fb^{-1} #sqrt{s} = 13 TeV");
	paveCMS->AddText("CMS Preliminary#sqrt{s} = 13 TeV");
	paveCMS->SetFillColor(0);
	paveCMS->SetBorderSize(0);
	paveCMS->SetTextSize(0.04);

	TLegend *leg1 =new TLegend(.2, .6, .9, .9);//7899//4899
	leg1->SetTextSize(0.06);
	leg1->SetFillColor(0); 
	leg1->SetBorderSize(0);
	leg1->AddEntry(h_CHFJet[0], "Reco jets", "L");
	leg1->AddEntry(h_CHFJet_gen[0], "Gen jets", "L");


	TLegend *leg2 =new TLegend(.7, .7, .9, .9);//7899//4899
	leg2->SetTextSize(0.03);
	leg2->SetFillColor(0); 
	leg2->SetBorderSize(0);
	leg2->AddEntry(h_CHFJet[0], "Reco jets", "L");
	leg2->AddEntry(h_CHFJet_gen[0], "Gen jets", "L");


	c_eta->cd();
	c_eta->SetLogy(1);
	h_ETAJet->GetXaxis()->SetTitle("Jet eta");
	h_ETAJet->GetYaxis()->SetTitle("Entries");
	h_ETAJet->GetYaxis()->SetTitleOffset(1.3);
	h_ETAJet_gen->SetLineColor(4); 
	h_ETAJet->SetLineStyle(1);
	h_ETAJet->SetMinimum(0.1);
	h_ETAJet->SetMaximum(10000);
	h_ETAJet->Draw("");
	h_ETAJet_gen->Draw("same hist");
	paveCMS ->Draw("same");
	leg2->Draw("same");


	c_SumEt->cd();
	c_SumEt->SetLogy(1);
	h_SumEt->GetXaxis()->SetTitle("SumEt");
	h_SumEt->GetYaxis()->SetTitle("Entries");
	h_SumEt->GetYaxis()->SetTitleOffset(1.3);
	h_SumEt_gen->SetLineColor(4); 
	h_SumEt->SetLineStyle(1);
	h_SumEt->SetMinimum(0.1);
	h_SumEt->SetMaximum(10000);
	h_SumEt->Draw("");
	h_SumEt_gen->Draw("same hist");
	paveCMS ->Draw("same");
	leg2->Draw("same");

	c_minDR->cd();
	c_minDR->SetLogy(1);
	h_minDR->GetXaxis()->SetTitle("minDR");
	h_minDR->GetYaxis()->SetTitle("Entries");
	h_minDR->GetYaxis()->SetTitleOffset(1.3);
	h_minDR->SetLineStyle(1);
	h_minDR->SetMinimum(0.1);
	h_minDR->SetMaximum(10000);
	h_minDR->Draw("");
	paveCMS ->Draw("same");


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
		h_CHFJet[iy]->GetXaxis()->SetTitle("Charged Hadron Fraction");
		h_CHFJet[iy]->GetYaxis()->SetTitle("Entries");
		h_CHFJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_CHFJet_gen[iy]->SetLineColor(4); 
		h_CHFJet[iy]->SetLineStyle(1);
		h_CHFJet[iy]->SetMinimum(0.1);
		h_CHFJet[iy]->SetMaximum(10000);
		h_CHFJet[iy]->Draw("");
		h_CHFJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();

		pad2->cd(iy+1);
		pad2->cd(iy+1)->SetLogy(1);
		h_NHFJet[iy]->GetXaxis()->SetTitle("Neutral Hadron Fraction");
		h_NHFJet[iy]->GetYaxis()->SetTitle("Entries");
		h_NHFJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_NHFJet[iy]->SetLineColor(1);
		h_NHFJet_gen[iy]->SetLineColor(4);
		h_NHFJet[iy]->SetLineStyle(1);
		if (iy<4) h_NHFJet[iy]->SetMinimum(0.1);
		h_NHFJet[iy]->SetMaximum(10000);
		h_NHFJet[iy]->Draw();
		h_NHFJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();


		pad3->cd(iy+1);
		pad3->cd(iy+1)->SetLogy(1);
		h_CEMFJet[iy]->GetXaxis()->SetTitle("Charged E/M Fraction");
		h_CEMFJet[iy]->GetYaxis()->SetTitle("Entries");
		h_CEMFJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_CEMFJet_gen[iy]->SetLineColor(4); 
		h_CEMFJet[iy]->SetLineStyle(1);
		h_CEMFJet[iy]->SetMinimum(0.1);
		h_CEMFJet[iy]->SetMaximum(10000);
		h_CEMFJet[iy]->Draw("");
		h_CEMFJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();


		pad4->cd(iy+1);
		pad4->cd(iy+1)->SetLogy(1);
		h_NEMFJet[iy]->GetXaxis()->SetTitle("Neutral E/M Fraction");
		h_NEMFJet[iy]->GetYaxis()->SetTitle("Entries");
		h_NEMFJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_NEMFJet_gen[iy]->SetLineColor(4); 
		h_NEMFJet[iy]->SetLineStyle(1);
		h_NEMFJet[iy]->SetMinimum(0.1);
		h_NEMFJet[iy]->SetMaximum(10000);
		h_NEMFJet[iy]->Draw("");
		h_NEMFJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();

		pad5->cd(iy+1);
		pad5->cd(iy+1)->SetLogy(1);
		h_MUFJet[iy]->GetXaxis()->SetTitle("Muon Fraction");
		h_MUFJet[iy]->GetYaxis()->SetTitle("Entries");
		h_MUFJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_MUFJet_gen[iy]->SetLineColor(4); 
		h_MUFJet[iy]->SetLineStyle(1);
		h_MUFJet[iy]->SetMinimum(0.1);
		h_MUFJet[iy]->SetMaximum(10000);
		h_MUFJet[iy]->Draw("");
		h_MUFJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();

		pad6->cd(iy+1);
		pad6->cd(iy+1)->SetLogy(1);
		h_PHIJet[iy]->GetXaxis()->SetTitle("Jet phi");
		h_PHIJet[iy]->GetYaxis()->SetTitle("Entries");
		h_PHIJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_PHIJet_gen[iy]->SetLineColor(4); 
		h_PHIJet[iy]->SetLineStyle(1);
		h_PHIJet[iy]->SetMinimum(0.1);
		h_PHIJet[iy]->SetMaximum(10000);
		h_PHIJet[iy]->Draw("");
		h_PHIJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();

		pad7->cd(iy+1);
		pad7->cd(iy+1)->SetLogy(1);
		h_CMJet[iy]->GetXaxis()->SetTitle("Charged Multiplicity");
		h_CMJet[iy]->GetYaxis()->SetTitle("Entries");
		h_CMJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_CMJet_gen[iy]->SetLineColor(4); 
		h_CMJet[iy]->SetLineStyle(1);
		h_CMJet[iy]->SetMinimum(0.1);
		h_CMJet[iy]->SetMaximum(10000);
		h_CMJet[iy]->Draw("");
		h_CMJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();

		pad8->cd(iy+1);
		pad8->cd(iy+1)->SetLogy(1);
		h_NMJet[iy]->GetXaxis()->SetTitle("Neutral Multiplicity");
		h_NMJet[iy]->GetYaxis()->SetTitle("Entries");
		h_NMJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_NMJet_gen[iy]->SetLineColor(4); 
		h_NMJet[iy]->SetLineStyle(1);
		h_NMJet[iy]->SetMinimum(0.1);
		h_NMJet[iy]->SetMaximum(10000);
		h_NMJet[iy]->Draw("");
		h_NMJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();

		pad9->cd(iy+1);
		pad9->cd(iy+1)->SetLogy(1);
		h_ptJet[iy]->GetXaxis()->SetTitle("Jet pT (GeV)");
		h_ptJet[iy]->GetYaxis()->SetTitle("Entries");
		h_ptJet[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_ptJet_gen[iy]->SetLineColor(4); 
		h_ptJet[iy]->SetLineStyle(1);
		h_ptJet[iy]->SetMinimum(0.1);
		h_ptJet[iy]->SetMaximum(10000);
		h_ptJet[iy]->Draw("");
		h_ptJet_gen[iy]->Draw("same hist");
		paveCMS ->Draw("same");
		teta->Draw();


		double res_mean = h_pT_resolution[iy]->GetMean();
		double res_rms = h_pT_resolution[iy]->GetRMS();

		char res_text[200];
		sprintf(res_text, "Offset = %3.2f ,  RES = %3.1f %%",1.0-res_mean,100*res_rms );

		TLegend *leg3 =new TLegend(.2, .75, .8, .85);//7899//4899
		leg3->SetTextSize(0.045);
		leg3->SetFillColor(0); 
		leg3->SetBorderSize(0);
		leg3->AddEntry(h_pT_resolution[iy], res_text , "p");


		pad10->cd(iy+1);
		pad10->cd(iy+1)->SetLogy(1);
		h_pT_resolution[iy]->GetXaxis()->SetTitle("pT_reco / pT_gen");
		h_pT_resolution[iy]->GetYaxis()->SetTitle("Entries");
		h_pT_resolution[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_pT_resolution[iy]->SetLineStyle(1);
		h_pT_resolution[iy]->SetMinimum(0.1);
		h_pT_resolution[iy]->SetMaximum(10000);
		h_pT_resolution[iy]->Draw("");
		paveCMS ->Draw("same");
		leg3->Draw("same");
		teta->Draw();

		pad11->cd(iy+1);
		h_reco_efficiency[iy]->GetXaxis()->SetTitle("Generated pT (GeV)");
		h_reco_efficiency[iy]->GetYaxis()->SetTitle("Reconstruction Efficiency");
		h_reco_efficiency[iy]->GetYaxis()->SetTitleOffset(1.3);
		h_reco_efficiency[iy]->SetLineStyle(1);
		h_reco_efficiency[iy]->SetMinimum(0.4);
		h_reco_efficiency[iy]->SetMaximum(1.2);
		h_reco_efficiency[iy]->Draw("hist");
		paveCMS ->Draw("same");
		teta->Draw();

	}
 
	pad1->cd(eta_bins+1);
	leg1->Draw();
	pad2->cd(eta_bins+1);
	leg1->Draw();
	pad3->cd(eta_bins+1);
	leg1->Draw();
	pad4->cd(eta_bins+1);
	leg1->Draw();
	pad5->cd(eta_bins+1);
	leg1->Draw();
	pad6->cd(eta_bins+1);
	leg1->Draw();
	pad7->cd(eta_bins+1);
	leg1->Draw();
	pad8->cd(eta_bins+1);
	leg1->Draw();
	pad9->cd(eta_bins+1);
	leg1->Draw();
	

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

		sprintf(filename,"%s/%s/%s_reco_efficiency.png",analyzer_path,output_directory,image_name);
		pad11->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_eta.png",analyzer_path,output_directory,image_name);
		c_eta->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_SumEt.png",analyzer_path,output_directory,image_name);
		c_SumEt->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_minDR.png",analyzer_path,output_directory,image_name);
		c_minDR->SaveAs(filename);
	}


	sprintf(filename,"%s/%s/%s_minDR.png",analyzer_path,output_directory,image_name);
		c_minDR->SaveAs(filename);

}

