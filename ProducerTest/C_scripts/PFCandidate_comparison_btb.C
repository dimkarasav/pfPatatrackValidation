//==========================================
// Initial Author: Dimitrios Karasavvas
// Date:   29 Jan 2021
//==========================================

#include "include_functions.h"




void PFCandidate_comparison_btb()
{

	gROOT->LoadMacro("setTDRStyle_teliko.C");
	setTDRStyle_teliko(); 
//	gStyle->SetOptStat(0);

	


	char analyzer_path[200] = {"/afs/cern.ch/work/d/dkarasav/public/ScoutingTriggers/CMSSW_11_0_0_pre7/src/Jakobs_producer/ProducerTest/"}; 
                  
	char output_directory[200] = {"deleteme/"};             
	char image_name[150] = {"Patatrack"}; 

//	double yBnd[]={0.0, 1.3, 2.4, 2.7, 3.0, 5.0}; 
	double yBnd[]={0.0, 1.3, 2.4, 2.7, 3.0}; 
//	double yBnd[]={0.0, 1.3, 2.4}; 




//	int Colors[NoFiles] = { 1, 4, 2 , 6, 3 , 68} ; // black, blue, red , magenta
	int Colors[] = { 1, 4, 2 , 6, 3, 7 , 28, 46} ; // black, blue, red , magenta, green, light blue ,  brown, redish
	bool scale_histos = true ;
	bool Save_Plots = true; 
	bool useWeights = true;

	double pTlowCut_leading = 60;
	double pTlowCut_subleading = 30;
	double pTmax = 450;
	double pThighCut = 5000;
	double deltaPhi = 2.7;

	double YaxisLowEndMultiplier = 0.00001;
	double YaxisHighEndMultiplier = 10.;


	char input_files[][500] ={

"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/QCD_samples/FullTracking_QCD_PU.root",
"/eos/cms/store/user/dkarasav/ScoutingFiles/Producer_files/110X_mcRun3_2021_realistic_v6-v1/QCD_samples/Patatrack_QCD_PU.root"

};


	int PadColumnsEtaBins = 3; 
	int PadRowsEtaBins = 2;
	int Canvas_XpixelsEtaBins = PadColumnsEtaBins*333;
	int Canvas_YpixelsEtaBins = PadRowsEtaBins*500;


	const int NoFiles = sizeof(input_files)/sizeof(input_files[0]);
	const int eta_bins = sizeof(yBnd)/sizeof(yBnd[0])-1;


	char legend_array[NoFiles][500] = { "FullTracking" ,  "PatatrackPixelTracks"  };


	TFile *f[NoFiles];
	TTree *tree[NoFiles];


	char eta_bins_legend[eta_bins][25];
	char name[256]; 

	int treeEntries[NoFiles];

	for (int iy=0; iy< eta_bins; iy++)
	{
		if (iy==0) sprintf( eta_bins_legend[iy], "|#eta|<%3.1f" , yBnd[iy+1] );  
		else
		{
			sprintf( eta_bins_legend[iy], "%3.1f<|#eta|<%3.1f" , yBnd[iy] ,yBnd[iy+1] );
		}
	}

	TH1D *h_ChargedHadron_pT[NoFiles][eta_bins], *h_NeutralHadron_pT[NoFiles][eta_bins], *h_ChargedHadron_SumpT[NoFiles][eta_bins], *h_NeutralHadron_SumpT[NoFiles][eta_bins], *h_NumberOfNeutralHadrons[NoFiles][eta_bins], *h_NumberOfChargedHadrons[NoFiles][eta_bins], *h_ChargedHadron_DR[NoFiles][eta_bins], *h_NeutralHadron_DR[NoFiles][eta_bins], *h_girth[NoFiles][eta_bins], *h_quadratic_moment[NoFiles][eta_bins];

	TH2D *h_ChargedHadron_eta_phi[NoFiles][eta_bins], *h_NeutralHadron_eta_phi[NoFiles][eta_bins];


	for (int NoFile=0; NoFile<NoFiles; NoFile++)
	{
		for (int iy=0; iy< eta_bins; iy++)
		{
	
			sprintf(name,"h_ChargedHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_pT[NoFile][iy] = new TH1D(name, "", 50,0,0.66*pTmax);
			if (useWeights) 	h_ChargedHadron_pT[NoFile][iy]->Sumw2();
				
			sprintf(name,"h_NeutralHadron_pT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_pT[NoFile][iy] = new TH1D(name, "",  50,0,0.66*pTmax);
			if (useWeights) 	h_NeutralHadron_pT[NoFile][iy]->Sumw2();
	
			sprintf(name,"h_ChargedHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_SumpT[NoFile][iy] = new TH1D(name, "", 50,0,pTmax);
			if (useWeights) 	h_ChargedHadron_SumpT[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_SumpT_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_SumpT[NoFile][iy] = new TH1D(name, "", 50,0,pTmax);
			if (useWeights) 	h_NeutralHadron_SumpT[NoFile][iy]->Sumw2();

			sprintf(name,"h_NumberOfChargedHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfChargedHadrons[NoFile][iy] = new TH1D(name, "", 50,0,50);
			if (useWeights) 	h_NumberOfChargedHadrons[NoFile][iy]->Sumw2();

			sprintf(name,"h_NumberOfNeutralHadrons_%s_bin%i",legend_array[NoFile],iy);
			h_NumberOfNeutralHadrons[NoFile][iy] = new TH1D(name, "",  50,0,50);
			if (useWeights) 	h_NumberOfNeutralHadrons[NoFile][iy]->Sumw2();

			sprintf(name,"h_ChargedHadron_DR_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_DR[NoFile][iy] = new TH1D(name, "",  50,0, 0.5);
			if (useWeights) 	h_ChargedHadron_DR[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_DR_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_DR[NoFile][iy] = new TH1D(name, "",  50,0, 0.5);
			if (useWeights) 	h_NeutralHadron_DR[NoFile][iy]->Sumw2();

			sprintf(name,"h_ChargedHadron_eta_phi_%s_bin%i",legend_array[NoFile],iy);
			h_ChargedHadron_eta_phi[NoFile][iy] = new TH2D(name, "",  30,-0.5, 0.5, 30, -0.5, 0.5);
			if (useWeights) 	h_ChargedHadron_eta_phi[NoFile][iy]->Sumw2();

			sprintf(name,"h_NeutralHadron_eta_phi_%s_bin%i",legend_array[NoFile],iy);
			h_NeutralHadron_eta_phi[NoFile][iy] = new TH2D(name, "",  30,-0.5, 0.5, 30, -0.5, 0.5);
			if (useWeights) 	h_NeutralHadron_eta_phi[NoFile][iy]->Sumw2();

			sprintf(name,"h_girth_%s_bin%i",legend_array[NoFile],iy);
			h_girth[NoFile][iy] = new TH1D(name, "",  50, 0., 1.0 );	
			if (useWeights) 	h_girth[NoFile][iy]->Sumw2();

			sprintf(name,"h_quadratic_moment_%s_bin%i",legend_array[NoFile],iy);
			h_quadratic_moment[NoFile][iy] = new TH1D(name, "",  50, 0., 1.0 );	
			if (useWeights) 	h_quadratic_moment[NoFile][iy]->Sumw2();

		}
	}




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

	
	std::vector < std::vector<int> >  	 *jetconstituents;
	std::vector<int> *cand_vector = 0;

	std::vector<float> *npu = 0;
	std::vector<int> *PileupInteractions = 0;
	std::vector<int> *PileupOriginBX = 0;

	float weight;


	const int 	max_pfcand = 10000;
	int n_pfcand;
	float	pfcandpt[max_pfcand];
	float	pfcandeta[max_pfcand];
	float	pfcandphi[max_pfcand];
	float	pfcandm[max_pfcand];
	float	pfcandpdgid[max_pfcand];
	float	pfcandvertex[max_pfcand];

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
		tree[NoFile]->SetBranchAddress("jetconstituents",&jetconstituents);
		tree[NoFile]->SetBranchAddress("pfcandpt",&pfcandpt);
		tree[NoFile]->SetBranchAddress("pfcandeta",&pfcandeta);
		tree[NoFile]->SetBranchAddress("pfcandphi",&pfcandphi);
		tree[NoFile]->SetBranchAddress("pdcandm",&pfcandm);
		tree[NoFile]->SetBranchAddress("pfcandpdgid",&pfcandpdgid);
		tree[NoFile]->SetBranchAddress("pfcandvertex",&pfcandvertex);
		tree[NoFile]->SetBranchAddress("SumEt",&SumEt);
		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);

		tree[NoFile]->SetBranchAddress("SumEt",&SumEt);
		tree[NoFile]->SetBranchAddress("nVtx",&nVtx);
		tree[NoFile]->SetBranchAddress("weight",&weight);
		tree[NoFile]->SetBranchAddress("npu",&npu);
		tree[NoFile]->SetBranchAddress("PileupInteractions",&PileupInteractions);
		tree[NoFile]->SetBranchAddress("PileupOriginBX",&PileupOriginBX);


		treeEntries[NoFile] = tree[NoFile]->GetEntries();
		////////////////find the tree with the maximum entries to loop over them.


		cout<<"\nNumber of entries for tree "<< input_files[NoFile] << "\n  =  " << treeEntries[NoFile] <<endl;



		int reco_size, gen_size;


//		for (int i=0; i<treeEntries[NoFile]; i++) //event loop
		for (int i=0; i<100000; i++) //event loop
		{

			if ( i >=  treeEntries[NoFile] ) continue; 

			if(i>0 && i%100000==0) cout<<"Event reached: "<<i<<endl;

			tree[NoFile]->GetEntry(i);

			reco_size = jpt->size();

			if(reco_size < 2) continue;	 									// skip event with less than 2 reco jets
			if ( fabs( phi->at(0) - phi->at(1) ) < deltaPhi ) continue; 	// skip event if the two leading jets are not back-to-back
			if ( jpt->at(0) < pTlowCut_leading || jpt->at(0) > pThighCut ) continue;// skip event if leading jet does not pass pT cuts
			if ( jpt->at(1) < pTlowCut_subleading || jpt->at(1) > pThighCut ) continue;// skip event if sub-leading jet does not pass pT cuts


//			cout << "\n Processing event " << i << " number of reco_jets : "<< reco_size <<endl; 
			for ( int j=0; j<2; j++) // loop over two leading jets
			{

				int jet_ybin = getBin(fabs(eta->at(j)),yBnd, eta_bins); //get the jet eta bin

				int cand_size = ( jetconstituents->at( j ) ).size();
//				cout << "jet "<< reco_jets_matched_sequence[j] << " cand size = " << cand_size << endl;

				float chargedCand_SumpT = 0;
				float neutralCand_SumpT = 0;
				int neutralPFcand_number = 0;
				int chargedPFcand_number = 0;
				double jet_girth = 0;
				double jet_quadr_moment = 0;

				for (int l=0; l< cand_size; l++)
				{
					int jet_candidate = (jetconstituents->at(j)).at(l) ; //get the jet candidate number
					//cout << "\njet candidate number: " << jet_candidate << endl;

					int 	cand_pdgID = pfcandpdgid[jet_candidate];
					float 	cand_eta = pfcandeta[jet_candidate];
					float 	cand_phi = pfcandphi[jet_candidate];
					float 	cand_pT = pfcandpt[jet_candidate];
					double 	cand_DR = sqrt( pow( eta->at(j)-cand_eta ,2) +pow( phi->at(j)-cand_phi ,2) );

					jet_girth = jet_girth + (cand_pT/jpt->at(j))*cand_DR;
					jet_quadr_moment = jet_quadr_moment + (cand_pT/jpt->at(j))*cand_DR*cand_DR;
	//				cout << "\n candidate pdgID = " << cand_pdgID << endl;


					int cand_ybin = getBin( fabs(pfcandeta[ jet_candidate ]),yBnd, eta_bins );  //get the pf candidate eta bin

					if (jet_ybin > -1)//fill hist's in the corresponding eta bin
					{

						if ( cand_pdgID==-211 || cand_pdgID==211 || cand_pdgID==11 || cand_pdgID==-11 || cand_pdgID==13 || cand_pdgID==-13) chargedPFcand_number = chargedPFcand_number +1;
						if ( cand_pdgID==-211 || cand_pdgID==211 )
						{
							chargedCand_SumpT = chargedCand_SumpT + cand_pT; // get the pT sum of charged candidates
							if(!useWeights)
							{

								h_ChargedHadron_pT[NoFile][jet_ybin]->Fill(cand_pT); // fill charged hadron pT's
								h_ChargedHadron_DR[NoFile][jet_ybin]->Fill( cand_DR );
								h_ChargedHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(j), cand_phi-phi->at(j)  );
							}
							else
							{	
								h_ChargedHadron_pT[NoFile][jet_ybin]->Fill(cand_pT, weight); // fill charged hadron pT's with weight
								h_ChargedHadron_DR[NoFile][jet_ybin]->Fill( cand_DR, weight );
								h_ChargedHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(j), cand_phi-phi->at(j), weight );
							}
						}
						if ( cand_pdgID==130 || cand_pdgID==22 || cand_pdgID==0 ) 	neutralPFcand_number = neutralPFcand_number +1;
						if ( cand_pdgID==130 )
						{
							neutralCand_SumpT = neutralCand_SumpT + cand_pT; // get the pT sum of charged candidates
							if(!useWeights)
							{
								h_NeutralHadron_pT[NoFile][jet_ybin]->Fill(cand_pT); // fill neutral hadron pT's
								h_NeutralHadron_DR[NoFile][jet_ybin]->Fill( cand_DR );
								h_NeutralHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(j), cand_phi-phi->at(j)  );
							}
							else
							{
								h_NeutralHadron_pT[NoFile][jet_ybin]->Fill(cand_pT, weight); // fill neutral hadron pT's with weight
								h_NeutralHadron_DR[NoFile][jet_ybin]->Fill( cand_DR, weight );
								h_NeutralHadron_eta_phi[NoFile][jet_ybin]->Fill( cand_eta-eta->at(j), cand_phi-phi->at(j), weight  );
							}
						}
					}
			//		cout << " cand " << cand_vector->at(l) << " , pt = " << pfcandpt[cand_vector->at(l)] << endl;  
				} // end of candidate loop

				if(!useWeights)
				{
					h_NumberOfNeutralHadrons[NoFile][jet_ybin]->Fill(neutralPFcand_number);
					h_NumberOfChargedHadrons[NoFile][jet_ybin]->Fill(chargedPFcand_number);
					h_ChargedHadron_SumpT[NoFile][jet_ybin]->Fill(chargedCand_SumpT);
					h_NeutralHadron_SumpT[NoFile][jet_ybin]->Fill(neutralCand_SumpT);
					h_girth[NoFile][jet_ybin]->Fill(jet_girth);
					h_quadratic_moment[NoFile][jet_ybin]->Fill(jet_quadr_moment);
				}
				else 
				{
					h_NumberOfNeutralHadrons[NoFile][jet_ybin]->Fill(neutralPFcand_number, weight);
					h_NumberOfChargedHadrons[NoFile][jet_ybin]->Fill(chargedPFcand_number, weight);
					h_ChargedHadron_SumpT[NoFile][jet_ybin]->Fill(chargedCand_SumpT, weight);
					h_NeutralHadron_SumpT[NoFile][jet_ybin]->Fill(neutralCand_SumpT, weight);
					h_girth[NoFile][jet_ybin]->Fill(jet_girth, weight);
					h_quadratic_moment[NoFile][jet_ybin]->Fill(jet_quadr_moment, weight);
				}
			} //  end of jet loop


			reco_size = 0;


		} // end of event loop


		delete tree[NoFile];
		delete f[NoFile];
	} // end of file loop



// =========================================================== Plotting staff ===================================================================

	TCanvas *pad1 = new TCanvas("pad1", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad1->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad2 = new TCanvas("pad2", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad2->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad3 = new TCanvas("pad3", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad3->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad4 = new TCanvas("pad4", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad4->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad5 = new TCanvas("pad5", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad5->Divide(PadColumnsEtaBins, PadRowsEtaBins);
	TCanvas *pad6 = new TCanvas("pad6", "",Canvas_XpixelsEtaBins,Canvas_YpixelsEtaBins);
	pad6->Divide(PadColumnsEtaBins, PadRowsEtaBins);
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


 	//dummy histograms to be used as frames in the plots
	TH1D *frameChargedHadron_SumpT[eta_bins], *frameNeutralHadron_SumpT[eta_bins], *frameChargedHadron_pT[eta_bins], *frameNeutralHadron_pT[eta_bins], *frameNumberOfChargedHadrons[eta_bins], *frameNumberOfNeutralHadrons[eta_bins], *frameChargedHadron_DR[eta_bins], *frameNeutralHadron_DR[eta_bins], *framegirth[eta_bins], *framequadratic_moment[eta_bins];

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

			
			if (NoFile==0 )  frameChargedHadron_SumpT[iy] = InitiateFrameOnCanvasPad(pad1, iy+1, "frameChargedHadron_SumpT", "Sum pT of Charged hadrons per jet (GeV)", "#jets", 0., pTmax, YaxisLowEndMultiplier*h_ChargedHadron_SumpT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ChargedHadron_SumpT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ChargedHadron_SumpT[NoFile][iy]->Integral()>0 ) h_ChargedHadron_SumpT[NoFile][iy]->Scale(h_ChargedHadron_SumpT[0][iy]->Integral()/h_ChargedHadron_SumpT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad1, iy+1, h_ChargedHadron_SumpT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad1->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad1->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNeutralHadron_SumpT[iy] = InitiateFrameOnCanvasPad(pad2, iy+1, "frameNeutralHadron_SumpT", "Sum pT of Neutral hadrons per jet (GeV)", "#jets", 0., pTmax, YaxisLowEndMultiplier*h_NeutralHadron_SumpT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NeutralHadron_SumpT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NeutralHadron_SumpT[NoFile][iy]->Integral()>0 ) h_NeutralHadron_SumpT[NoFile][iy]->Scale(h_NeutralHadron_SumpT[0][iy]->Integral()/h_NeutralHadron_SumpT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad2, iy+1, h_NeutralHadron_SumpT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad2->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad2->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameChargedHadron_pT[iy] = InitiateFrameOnCanvasPad(pad3, iy+1, "frameChargedHadron_pT", "Charged Hadron pT (GeV)", "# Charged Hadrons", 0., pTmax, YaxisLowEndMultiplier*h_ChargedHadron_pT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ChargedHadron_pT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ChargedHadron_pT[NoFile][iy]->Integral()>0 ) h_ChargedHadron_pT[NoFile][iy]->Scale(h_ChargedHadron_pT[0][iy]->Integral()/h_ChargedHadron_pT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad3, iy+1, h_ChargedHadron_pT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) {pad3->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad3->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNeutralHadron_pT[iy] = InitiateFrameOnCanvasPad(pad4, iy+1, "frameNeutralHadron_pT", "Neutral Hadron pT (GeV)", "# Neutral Hadrons", 0., pTmax, YaxisLowEndMultiplier*h_NeutralHadron_pT[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NeutralHadron_pT[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NeutralHadron_pT[NoFile][iy]->Integral()>0 ) h_NeutralHadron_pT[NoFile][iy]->Scale(h_NeutralHadron_pT[0][iy]->Integral()/h_NeutralHadron_pT[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad4, iy+1, h_NeutralHadron_pT[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad4->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad4->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNumberOfChargedHadrons[iy] = InitiateFrameOnCanvasPad(pad5, iy+1, "frameNumberOfChargedHadrons", "# of charged candidates per jet", "# jets", 0., 50., YaxisLowEndMultiplier*h_NumberOfChargedHadrons[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NumberOfChargedHadrons[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NumberOfChargedHadrons[NoFile][iy]->Integral()>0 ) h_NumberOfChargedHadrons[NoFile][iy]->Scale(h_NumberOfChargedHadrons[0][iy]->Integral()/h_NumberOfChargedHadrons[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad5, iy+1, h_NumberOfChargedHadrons[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad5->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad5->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNumberOfNeutralHadrons[iy] = InitiateFrameOnCanvasPad(pad6, iy+1, "frameNumberOfNeutralHadrons", "# of Neutral candidates per jet", "# jets", 0., 50., YaxisLowEndMultiplier*h_NumberOfNeutralHadrons[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NumberOfNeutralHadrons[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_NumberOfNeutralHadrons[NoFile][iy]->Integral()>0 ) h_NumberOfNeutralHadrons[NoFile][iy]->Scale(h_NumberOfNeutralHadrons[0][iy]->Integral()/h_NumberOfNeutralHadrons[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad6, iy+1, h_NumberOfNeutralHadrons[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad6->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad6->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameChargedHadron_DR[iy] = InitiateFrameOnCanvasPad(pad_ChargedCand_DR, iy+1, "frameChargedHadron_DR", "Jet-ChargedCandidate DR", "# jets", 0., 0.5, YaxisLowEndMultiplier*h_ChargedHadron_DR[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_ChargedHadron_DR[NoFile][iy]->GetMaximum() , true, paveCMS);
			if (scale_histos && h_ChargedHadron_DR[NoFile][iy]->Integral()>0 ) h_ChargedHadron_DR[NoFile][iy]->Scale(h_ChargedHadron_DR[0][iy]->Integral()/h_ChargedHadron_DR[NoFile][iy]->Integral());
			DrawHistoToCanvasPad(pad_ChargedCand_DR, iy+1, h_ChargedHadron_DR[NoFile][iy], Colors[NoFile], 1);
			if(eta_bins <= 1 && NoFile == NoFiles-1 ) { pad_ChargedCand_DR->cd(iy+1); leg1->Draw("same");  }
			else if( eta_bins > 1 && NoFile== 0) { pad_ChargedCand_DR->cd(iy+1); teta->Draw(); }


			if (NoFile==0 )  frameNeutralHadron_DR[iy] = InitiateFrameOnCanvasPad(pad_NeutralCand_DR, iy+1, "frameNeutralHadron_DR", "Jet-NeutralCandidate DR", "# jets", 0., 0.5, YaxisLowEndMultiplier*h_NeutralHadron_DR[NoFile][iy]->GetMaximum(), YaxisHighEndMultiplier*h_NeutralHadron_DR[NoFile][iy]->GetMaximum() , true, paveCMS);
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
				pad1->cd(eta_bins+1);    			leg1->Draw();
				pad2->cd(eta_bins+1);    			leg1->Draw();
				pad3->cd(eta_bins+1);    			leg1->Draw();
				pad4->cd(eta_bins+1);   			leg1->Draw();
				pad5->cd(eta_bins+1);   			leg1->Draw();
				pad6->cd(eta_bins+1);   			leg1->Draw();
				pad_girth->cd(eta_bins+1);   		leg1->Draw();
				pad_quadr_moment->cd(eta_bins+1);	leg1->Draw();

				pad_ChargedCand_DR->cd(eta_bins+1);	leg1->Draw();
				pad_NeutralCand_DR->cd(eta_bins+1);	leg1->Draw();
			}

		} //end of loop on eta bins
	} // end of loop on files


	char filename[1500]; 
	if(Save_Plots)
	{ 
		sprintf(filename,"%s/%s/%s_ChargedPFcand_SumpT.png",analyzer_path,output_directory,image_name);
		pad1->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_SumpT.png",analyzer_path,output_directory,image_name);
		pad2->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_ChargedPFcand_pT.png",analyzer_path,output_directory,image_name);
		pad3->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_pT.png",analyzer_path,output_directory,image_name);
		pad4->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_ChargedPFcand_number.png",analyzer_path,output_directory,image_name);
		pad5->SaveAs(filename);

		sprintf(filename,"%s/%s/%s_NeutralPFcand_number.png",analyzer_path,output_directory,image_name);
		pad6->SaveAs(filename);

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
	
