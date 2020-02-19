#ifndef PFTREEPRODUCER2_H
#define PFTREEPRODUCER2_H

// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>

// ROOT includes
#include "TTree.h"
#include "TLorentzVector.h"
#include "TPRegexp.h"

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
//#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingParticle.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Scouting/interface/ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/ScoutingMuon.h"
#include "DataFormats/Scouting/interface/ScoutingVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h" //edw
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h" //edw
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include <DataFormats/TrackReco/interface/TrackBase.h> //edw
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"



class pfTreeProducer_AddedCandidates : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
	public:
		explicit pfTreeProducer_AddedCandidates(const edm::ParameterSet&);
		~pfTreeProducer_AddedCandidates();
		
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	
	
	private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        
        const edm::EDGetTokenT<std::vector<ScoutingPFJet> >     jetsToken;
        const edm::EDGetTokenT<std::vector<reco::GenJet> >      genjetToken;
        const edm::EDGetTokenT<std::vector<ScoutingMuon> >      muonsToken;
	const edm::EDGetTokenT<std::vector<ScoutingParticle> >  	pfcandsToken;

        
		edm::EDGetTokenT<double> metpttoken, metphitoken;

        const edm::EDGetTokenT<std::vector<ScoutingVertex> >      VtxToken;
        const edm::EDGetTokenT<GenEventInfoProduct>             genEvtInfoToken;
        const edm::EDGetTokenT<std::vector<reco::GenParticle> > gensToken;
       
        TTree* tree;

  	//Run and lumisection
  	int run;
  	int lumSec;
	double recx;
	double recy;
	double my;
	double mx;
	double u1;
	double u2;
	double met;

	double dijetMass;

	std::vector<float>           gen_jpt;
	std::vector<float>           gen_chf;
	std::vector<float>           gen_nhf;
	std::vector<float>           gen_nemf;
	std::vector<float>           gen_cemf;
	std::vector<float>           gen_muf;
//	std::vector<float>           gen_hf_hf;
//	std::vector<float>           gen_hf_emf;
//	std::vector<float>           gen_hof;
	std::vector<float>           gen_chm;
	std::vector<float>           gen_chMult;
	std::vector<float>           gen_neMult;
	std::vector<float>           gen_npr;
	std::vector<float>           gen_eta;
	std::vector<float>           gen_phi;
	std::vector<float>           gen_jetSumEt;



	std::vector<float>           jpt;
	std::vector<float>           chf;
	std::vector<float>           nhf;
	std::vector<float>           nemf;
	std::vector<float>           cemf;
	std::vector<float>           muf;
	std::vector<float>           hf_hf;
	std::vector<float>           hf_emf;
	std::vector<float>           hof;
	std::vector<float>           chm;
	std::vector<float>           chMult;
	std::vector<float>           neMult;
	std::vector<float>           npr;
	std::vector<float>           eta;
	std::vector<float>           phi;
	std::vector<float>           jetSumEt;
	std::vector < std::vector<int> >  	 jetconstituents;

	int nVtx;
	double SumEt, gen_SumEt;

	int nch_;
	int nnh_;
	int nph_;


	  //PFCand
	const static int 	max_pfcand = 10000;
	UInt_t n_pfcand;
	Float_t 	    pfcandpt[max_pfcand];
	Float_t           pfcandeta[max_pfcand];
	Float_t           pfcandphi[max_pfcand];
	Float_t	    pdcandm[max_pfcand];
	Float_t	    pfcandpdgid[max_pfcand];
	Float_t	    pfcandvertex[max_pfcand];
  

};

//Constructor
pfTreeProducer_AddedCandidates::pfTreeProducer_AddedCandidates(const edm::ParameterSet& iConfig): 
  jetsToken            (consumes<std::vector<ScoutingPFJet> >           (iConfig.getParameter<edm::InputTag>("jetsAK4"))),
//  jetsToken            (consumes<std::vector<ScoutingPFJet> >           (iConfig.getParameter<edm::InputTag>("hltAK4PFJets"))),
  genjetToken            (consumes<std::vector<reco::GenJet> >           (iConfig.getParameter<edm::InputTag>("genJet"))),
  muonsToken            (consumes<std::vector<ScoutingMuon> >           (iConfig.getParameter<edm::InputTag>("muons"))),
  pfcandsToken             (consumes<std::vector<ScoutingParticle> >         (iConfig.getParameter<edm::InputTag>("pfcands"))), 
  metpttoken            (consumes<double>           (iConfig.getParameter<edm::InputTag>("metpt"))),
  metphitoken            (consumes<double>           (iConfig.getParameter<edm::InputTag>("metphi"))),
  VtxToken            (consumes<std::vector<ScoutingVertex> >           (iConfig.getParameter<edm::InputTag>("primVtx"))),
  gensToken               (consumes<std::vector<reco::GenParticle> >        (iConfig.getParameter<edm::InputTag>("genpart")))

   
{
usesResource("TFileService");	
}

//destructor
pfTreeProducer_AddedCandidates::~pfTreeProducer_AddedCandidates() {
}

void pfTreeProducer_AddedCandidates::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    
}

void pfTreeProducer_AddedCandidates::beginJob() {
    // Access the TFileService
    edm::Service<TFileService> fs; //edw

	std::cout << " HELLO \n" ;
    //book Histos

    
    // Create the TTree
    tree = fs->make<TTree>("tree"      , "tree");
 //   tree->Branch("recx"                , &recx                         , "recx/D"        );
//    tree->Branch("recy"                , &recy                         , "recy/D"        );
    tree->Branch("mx"                  , &mx                           , "mx/D"        );
    tree->Branch("my"                  , &my                           , "my/D"        );
//    tree->Branch("u1"                  , &u1                           , "u1/D"        );
//    tree->Branch("u2"                  , &u2                           , "u2/D"        );



    tree->Branch("jpt"                 , "std::vector<float>"          , &jpt      , 32000, 0);
    tree->Branch("eta"                 , "std::vector<float>"          , &eta      , 32000, 0);
    tree->Branch("phi"                 , "std::vector<float>"          , &phi      , 32000, 0);
    tree->Branch("chf"                 , "std::vector<float>"          , &chf      , 32000, 0);
    tree->Branch("nhf"                 , "std::vector<float>"          , &nhf      , 32000, 0);
    tree->Branch("nemf"                , "std::vector<float>"          , &nemf     , 32000, 0);
    tree->Branch("cemf"                , "std::vector<float>"          , &cemf     , 32000, 0);
    tree->Branch("muf"                 , "std::vector<float>"          , &muf      , 32000, 0);
    tree->Branch("hf_hf"               , "std::vector<float>"          , &hf_hf    , 32000, 0);
    tree->Branch("hf_emf"              , "std::vector<float>"          , &hf_emf   , 32000, 0);
    tree->Branch("hof"                 , "std::vector<float>"          , &hof      , 32000, 0);
    tree->Branch("chm"                 , "std::vector<float>"          , &chm      , 32000, 0);
    tree->Branch("chMult"              , "std::vector<float>"          , &chMult   , 32000, 0);
    tree->Branch("neMult"              , "std::vector<float>"          , &neMult   , 32000, 0);
    tree->Branch("npr"                 , "std::vector<float>"          , &npr      , 32000, 0);
    tree->Branch("jetSumEt"            , "std::vector<float>"          , &jetSumEt , 32000, 0);
    tree->Branch("jetconstituents"            	, "std::vector< vector<int> >"   , &jetconstituents 		, 32000, 0);
    tree->Branch("nVtx"                , &nVtx                         , "nVtx/I"         );
    tree->Branch("SumEt"               , &SumEt            		       , "SumEt/D"        );




    tree->Branch("gen_jpt"                 , "std::vector<float>"          , &gen_jpt      , 32000, 0);
    tree->Branch("gen_eta"                 , "std::vector<float>"          , &gen_eta      , 32000, 0);
    tree->Branch("gen_phi"                 , "std::vector<float>"          , &gen_phi      , 32000, 0);
    tree->Branch("gen_chf"                 , "std::vector<float>"          , &gen_chf      , 32000, 0);
    tree->Branch("gen_nhf"                 , "std::vector<float>"          , &gen_nhf      , 32000, 0);
    tree->Branch("gen_nemf"                , "std::vector<float>"          , &gen_nemf     , 32000, 0);
    tree->Branch("gen_cemf"                , "std::vector<float>"          , &gen_cemf     , 32000, 0);
    tree->Branch("gen_muf"                 , "std::vector<float>"          , &gen_muf      , 32000, 0);
    tree->Branch("gen_chm"                 , "std::vector<float>"          , &gen_chm      , 32000, 0);
    tree->Branch("gen_chMult"              , "std::vector<float>"          , &gen_chMult   , 32000, 0);
    tree->Branch("gen_neMult"              , "std::vector<float>"          , &gen_neMult   , 32000, 0);
    tree->Branch("gen_npr"                 , "std::vector<float>"          , &gen_npr      , 32000, 0);
    tree->Branch("gen_jetSumEt"            , "std::vector<float>"          , &gen_jetSumEt , 32000, 0);
    tree->Branch("gen_SumEt"                  , &gen_SumEt                 , "gen_SumEt/D"        );

	tree->Branch("n_pfcand"            	   ,&n_pfcand 		,"n_pfcand/i"		);	
	tree->Branch("pfcandpt"        	   ,pfcandpt 		,"pfcandpt[n_pfcand]/f" );
	tree->Branch("pfcandeta"            	   ,pfcandeta 		,"pfcandeta[n_pfcand]/f" );
	tree->Branch("pfcandphi"            	   ,pfcandphi		,"pfcandphi[n_pfcand]/f" );
	tree->Branch("pdcandm"            	   ,pdcandm 		,"pfcandm[n_pfcand]/f" );
	tree->Branch("pfcandpdgid"               ,pfcandpdgid		,"pfcandpdgid[n_pfcand]/f" );
	tree->Branch("pfcandvertex"              ,pfcandvertex 	,"pfcandvertex[n_pfcand]/f" );




    // Event weights
    return ;
}


void pfTreeProducer_AddedCandidates::endJob() {
}



void pfTreeProducer_AddedCandidates::endRun(edm::Run const&, edm::EventSetup const&) {
}


void pfTreeProducer_AddedCandidates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;
    using namespace reco;
    
	//std::cout << "entered here \n" ;


    jpt.clear();
	eta.clear();
	phi.clear();
	chf.clear();
	nhf.clear();
	nemf.clear();
	cemf.clear();
	muf.clear();
	jetSumEt.clear();
	hf_hf.clear();
	hf_emf.clear();
	hof.clear();
	chm.clear();
	chMult.clear();
	neMult.clear();
	npr.clear();
	nVtx = 0;
	SumEt = 0;

	gen_jpt.clear();
	gen_eta.clear();
	gen_phi.clear();
	gen_chf.clear();
	gen_nhf.clear();
	gen_nemf.clear();    
	gen_cemf.clear();
	gen_muf.clear();
	gen_jetSumEt.clear();
	gen_chm.clear();
	gen_chMult.clear();
	gen_neMult.clear();
	gen_npr.clear();
	jetconstituents.clear();
	gen_SumEt = 0;

//	std::cout << "1 \n" ;
    
    //recx.clear();
    //recy.clear();
    // Handles to the EDM content
    
    Handle<vector<ScoutingPFJet>> jetsH;
    iEvent.getByToken(jetsToken, jetsH);

    Handle<vector<GenJet>> genjetsH;
    iEvent.getByToken(genjetToken, genjetsH);

    Handle<vector<ScoutingMuon>> muonsH;
    iEvent.getByToken(muonsToken, muonsH);
    Handle<double> metptH;
    iEvent.getByToken(metpttoken, metptH);
    Handle<double> metphiH;
    iEvent.getByToken(metphitoken, metphiH);
    
    Handle<vector<GenParticle> > gensH;
    iEvent.getByToken(gensToken, gensH);
    
	Handle<vector<ScoutingVertex>> VtxH;
    iEvent.getByToken(VtxToken, VtxH);

  Handle<vector<ScoutingParticle> > pfcandsH;
  iEvent.getByToken(pfcandsToken, pfcandsH);
    
    TLorentzVector mass;
    int jetcount=0;
	double gen_tot_energy = 0; 
	double tot_energy = 0;
	int gjet_count = 0;

	for (auto gjet = genjetsH->begin(); gjet != genjetsH->end(); ++gjet) 
	{	
		//std::cout << "2 \n" ;
	//	std::cout <<" gjet count = " << gjet_count <<" \n" ;
		TLorentzVector gtemp;
		gtemp.SetPtEtaPhiM(gjet->pt(),gjet->eta(),gjet->phi(),gjet->mass());	    

		if(gtemp.Pt()<30.) continue;

		double jet_energy = gjet->neutralEmEnergy() + gjet->chargedHadronEnergy() + gjet->neutralHadronEnergy() + gjet->chargedEmEnergy() + gjet->muonEnergy();

		gen_jetSumEt.push_back( jet_energy );
		gen_jpt.push_back( gjet->pt() );
		gen_eta.push_back( gjet->eta() );
		gen_phi.push_back( gjet->phi() );
		//gen_SumEt.push_back( jet_energy );
		gen_nhf.push_back(    gjet->neutralHadronEnergy() /  jet_energy  );
		gen_chf.push_back(    gjet->chargedHadronEnergy() /  jet_energy  );
		gen_nemf.push_back(   gjet->neutralEmEnergy()     /  jet_energy  );
		gen_cemf.push_back(   gjet->chargedEmEnergy()     /  jet_energy  );
		gen_muf.push_back(    gjet->muonEnergy()          /  jet_energy  );
	//	std::cout << "3 \n" ;
	//	gen_hf_hf.push_back(  gjet->HFHadronEnergy()      /  jet_energy  );
	//      gen_hf_emf.push_back( gjet->HFEMEnergy()          /  jet_energy  );
	//      gen_hof.push_back(    gjet->HOEnergy()            /  jet_energy  );
		gen_chm.push_back(    gjet->chargedHadronMultiplicity() );
		gen_chMult.push_back( gjet->chargedHadronMultiplicity() + gjet->chargedEmMultiplicity() + gjet->muonMultiplicity() );
		gen_neMult.push_back( gjet->neutralHadronMultiplicity() + gjet->neutralEmMultiplicity() );  
		gen_npr.push_back(    gjet->chargedHadronMultiplicity() + gjet->chargedEmMultiplicity() + gjet->muonMultiplicity() + gjet->neutralHadronMultiplicity() + gjet->neutralEmMultiplicity() );

		gen_tot_energy = gen_tot_energy + jet_energy; 
//		std::cout << "pT = " << gjet->pt() << ",  eta = " << gjet->eta() << ",  phi = " << gjet->phi() << ",  nhf = " << gjet->neutralHadronEnergy() <<  ",  chf = " << gjet->chargedHadronEnergy() <<   ",  nemf = " << gjet->neutralEmEnergy() <<  ",  cemf = " << gjet->chargedEmEnergy() <<  ",  muf = " << gjet->muonEnergy() <<  ",  CM = " << gjet->chargedHadronMultiplicity() + gjet->chargedEmMultiplicity() + gjet->muonMultiplicity() <<  ",  NM = " << gjet->neutralHadronMultiplicity() + gjet->neutralEmMultiplicity() << "\n" ;
	gjet_count++;
	} //test




		TLorentzVector rtemp;
		for (auto ijet = jetsH->begin(); ijet != jetsH->end(); ++ijet) 
		{
		//	std::cout << "entered here \n" ;
		//	std::cout << ijet->pt() <<"\n" ;
		//	std::cout << "4 \n" ;
			if(ijet->pt()<1.) continue;

			double jet_energy = ijet->photonEnergy() + ijet->chargedHadronEnergy() + ijet->neutralHadronEnergy() + ijet->electronEnergy() + ijet->muonEnergy();

			jetSumEt.push_back( jet_energy );
			jpt.push_back(    ijet->pt() );
			eta.push_back( ijet->eta() );
			phi.push_back( ijet->phi() );
			//SumEt.push_back( jet_energy );
			nhf.push_back(    ijet->neutralHadronEnergy() /  jet_energy  );
		  	chf.push_back(    ijet->chargedHadronEnergy() /  jet_energy  );
			nemf.push_back(   ijet->photonEnergy()     /  jet_energy  );
			cemf.push_back(   ijet->electronEnergy()     /  jet_energy  );
			muf.push_back(    ijet->muonEnergy()          /  jet_energy  );
			hf_hf.push_back(  ijet->HFHadronEnergy()      /  jet_energy  );
		  	hf_emf.push_back( ijet->HFEMEnergy()          /  jet_energy  );
			hof.push_back(    ijet->HOEnergy()            /  jet_energy  );
			jetconstituents.push_back( ( ijet->constituents() ) ); 
			chm.push_back(    ijet->chargedHadronMultiplicity() );
			chMult.push_back( ijet->chargedHadronMultiplicity() + ijet->electronMultiplicity() + ijet->muonMultiplicity() );
			neMult.push_back( ijet->neutralHadronMultiplicity() + ijet->photonMultiplicity() );  
			npr.push_back(    ijet->chargedHadronMultiplicity() + ijet->electronMultiplicity() + ijet->muonMultiplicity() + ijet->neutralHadronMultiplicity() + ijet->photonMultiplicity() );
			tot_energy = tot_energy + jet_energy;


	//		std::cout << "pT = " << ijet->pt() << ",  eta = " << ijet->eta() << ",  phi = " << ijet->phi() << ",  nhf = " << ijet->neutralHadronEnergy() <<  ",  chf = " << ijet->chargedHadronEnergy() <<   ",  nemf = " << ijet->photonEnergy() <<  ",  cemf = " << ijet->electronEnergy() <<  ",  muf = " << ijet->muonEnergy() <<  ",  CM = " << ijet->chargedHadronMultiplicity() + ijet->electronMultiplicity() + ijet->muonMultiplicity() <<  ",  NM = " << ijet->neutralHadronMultiplicity() + ijet->photonMultiplicity() << "\n" ;

		    float eta  = ijet->eta();
		    float pt   = ijet->pt();
		}

		
	
	n_pfcand = 0;
	for (auto pfcands_iter = pfcandsH->begin(); pfcands_iter != pfcandsH->end(); ++pfcands_iter)
	{
	//	pfcandpt[n_pfcand]=MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->pt()));
	//	pfcandeta[n_pfcand]=MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->eta()));
	//	pfcandphi[n_pfcand]=MiniFloatConverter::float16to32(MiniFloatConverter::float32to16(pfcands_iter->phi()));

		pfcandpt[n_pfcand]=pfcands_iter->pt();
		pfcandeta[n_pfcand]=pfcands_iter->eta();
		pfcandphi[n_pfcand]= pfcands_iter->phi();
		pdcandm[n_pfcand]=pfcands_iter->m();
		pfcandpdgid[n_pfcand]=pfcands_iter->pdgId();
		pfcandvertex[n_pfcand]=pfcands_iter->vertex();
		n_pfcand++;
	} 


 


	
//	std::cout << VtxH->size() << "\n" ;
	nVtx = VtxH->size();
	SumEt = tot_energy; 
	gen_SumEt = gen_tot_energy;
	//		std::cout << "6 \n" ;


    //cout<"made it to Z"<<endl;
    
    //cout<"made it to after Z"<<endl;
    mx=TMath::Sin((*metphiH))*(*metptH);
    my=TMath::Cos((*metphiH))*(*metptH);

	
    
/*
    TVector2 vMet(my, mx);
    TVector2 vZPt((muonsV[0]+muonsV[1]).Px(),(muonsV[0]+muonsV[1]).Py());
    TVector2 vU = -1.0*(vMet+vZPt);
    
    u1 = (((muonsV[0]+muonsV[1]).Px())*(vU.Px()) + ((muonsV[0]+muonsV[1]).Py())*(vU.Py()))/((muonsV[0]+muonsV[1]).Pt());  // u1 = (pT . u)/|pT|
    u2 = (((muonsV[0]+muonsV[1]).Px())*(vU.Py()) - ((muonsV[0]+muonsV[1]).Py())*(vU.Px()))/((muonsV[0]+muonsV[1]).Pt());  // u2 = (pT x u)/|pT|
   // uperp->Fill(u1);
  //  upara->Fill(u2);
*/
	//		std::cout << "7 \n" ;
    tree->Fill();	
	//		std::cout << "8 \n" ;
  //  }	
}




void pfTreeProducer_AddedCandidates::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void pfTreeProducer_AddedCandidates::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void pfTreeProducer_AddedCandidates::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(pfTreeProducer_AddedCandidates);

#endif 
