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

//fastjet includes
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/contrib/RecursiveSoftDrop.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/JadePlugin.hh"
#include "fastjet/contrib/SoftKiller.hh"



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
        const edm::EDGetTokenT<std::vector<reco::GenParticle> > gensToken;
        const edm::EDGetTokenT<GenEventInfoProduct>             genEvtInfoToken;
	const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PUInfoToken;

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
	double met, met_phi;

	double dijetMass;

	std::vector<float>           gen_jpt;
	std::vector<float>           gen_mass;
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

	std::vector<float>           npu_;
	std::vector<int>           Number_interactions;
	std::vector<int>           OriginBX;

	int nVtx;
	double SumEt, gen_SumEt;

	int nch_;
	int nnh_;
	int nph_;
	float weight_;



	std::vector<float>           genpart_pt;
	std::vector<float>           genpart_eta;
	std::vector<float>           genpart_phi;
	std::vector<float>           genpart_m;
	std::vector<int>	     genpart_pdgid;
	std::vector<float>           genpart_px;
	std::vector<float>           genpart_py;
	std::vector<float>           genpart_pz;
	std::vector<float>           genpart_p;
	std::vector<float>           genpart_E;
	std::vector<int>             genpart_numDaught;
	std::vector<bool>            genpart_isQfromZ;



	//AK8 stuff
	std::vector<float>           ak8_pt;
	std::vector<float>           ak8_eta;
	std::vector<float>           ak8_phi;
	std::vector<float>           ak8_m;
	std::vector<float>           ak8_softDrop_m;
	std::vector<float>           ak8_trimmed_m;
	std::vector<float>           ak8_tau1;
	std::vector<float>           ak8_tau2;
	std::vector<float>           ak8_tau3;
	std::vector<float>           ak8_tau4;
	std::vector<float>           ak8_tau5;


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
  gensToken               (consumes<std::vector<reco::GenParticle> >        (iConfig.getParameter<edm::InputTag>("genpart"))),

  genEvtInfoToken           ( consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>  ("ptHat"))),

  PUInfoToken              ( consumes<std::vector<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>    ("pu")))
//  genEvtInfoToken       (consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("PdfInfoTag", edm::InputTag("generator"))))

//  genEvtInfoToken       (consumes<GenEventInfoProduct>(iConfig.getUntrackedParameter<edm::InputTag>("generator")))
 // genEvtInfoToken  (mayConsume<GenEventInfoProduct>(edm::InputTag("generator")))
   
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
    tree->Branch("met"               , &met            		       , "met/D"        );
    tree->Branch("met_phi"               , &met_phi            		       , "met_phi/D"        );

	tree->Branch("npu"                  ,"vector<float>"       , &npu_, 32000, 0 );
	tree->Branch("PileupInteractions"   ,"vector<int>"       , &Number_interactions, 32000, 0 );
	tree->Branch("PileupOriginBX"       ,"vector<int>"       , &OriginBX , 32000, 0);
	tree->Branch("weight"               ,&weight_            ,"weight/F");


    tree->Branch("genpart_pt"				,"std::vector<float>"           , &genpart_pt		, 32000, 0 );
    tree->Branch("genpart_eta"				,"std::vector<float>"           , &genpart_eta	 	, 32000, 0 );
    tree->Branch("genpart_phi"				,"std::vector<float>"           , &genpart_phi	        , 32000, 0 );
    tree->Branch("genpart_m"				,"std::vector<float>"           , &genpart_m		, 32000, 0 );
    tree->Branch("genpart_pdgid"			,"std::vector<int>"             , &genpart_pdgid	, 32000, 0 );
    tree->Branch("genpart_px"				,"std::vector<float>"           , &genpart_px		, 32000, 0 );
    tree->Branch("genpart_py"				,"std::vector<float>"           , &genpart_py		, 32000, 0 );
    tree->Branch("genpart_pz"				,"std::vector<float>"           , &genpart_pz		, 32000, 0 );
    tree->Branch("genpart_p"				,"std::vector<float>"           , &genpart_p		, 32000, 0 );
    tree->Branch("genpart_E"				,"std::vector<float>"           , &genpart_E		, 32000, 0 );
    tree->Branch("genpart_numDaught"			,"std::vector<int>"             , &genpart_numDaught	, 32000, 0 );
    tree->Branch("genpart_isQfromZ"			,"std::vector<bool>"            , &genpart_isQfromZ	, 32000, 0 );


	//AK8 stuff
	tree->Branch("ak8_pt"				,"std::vector<float>"           , &ak8_pt		 , 32000, 0 );
	tree->Branch("ak8_eta"				,"std::vector<float>"           , &ak8_eta	 , 32000, 0 );
	tree->Branch("ak8_phi"				,"std::vector<float>"           , &ak8_phi	 , 32000, 0 );
	tree->Branch("ak8_m"				,"std::vector<float>"           , &ak8_m		 , 32000, 0 );
	tree->Branch("ak8_softDrop_m"		,"std::vector<float>"           , &ak8_softDrop_m		 , 32000, 0 );
	tree->Branch("ak8_trimmed_m"		,"std::vector<float>"           , &ak8_trimmed_m		 , 32000, 0 );
	tree->Branch("ak8_tau1"				,"std::vector<float>"           , &ak8_tau1	 , 32000, 0 );
	tree->Branch("ak8_tau2"				,"std::vector<float>"           , &ak8_tau2	 , 32000, 0 );
	tree->Branch("ak8_tau3"				,"std::vector<float>"           , &ak8_tau3	 , 32000, 0 );
	tree->Branch("ak8_tau4"				,"std::vector<float>"           , &ak8_tau4	 , 32000, 0 );
	tree->Branch("ak8_tau5"				,"std::vector<float>"           , &ak8_tau5	 , 32000, 0 );


    tree->Branch("gen_mass"                , "std::vector<float>"          , &gen_mass     , 32000, 0);
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
	using namespace fastjet;
	using namespace fastjet::contrib;
    
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
	met = 0;
	met_phi = 0;

	weight_ = 0;
	Number_interactions.clear();
	OriginBX.clear();
	npu_.clear();

	genpart_pt.clear();
	genpart_eta.clear();
	genpart_phi.clear();
	genpart_m.clear();
	genpart_pdgid.clear();
	genpart_px.clear();
	genpart_py.clear();
	genpart_pz.clear();
	genpart_E.clear();
	genpart_p.clear();
	genpart_numDaught.clear();
	genpart_isQfromZ.clear();


	ak8_pt.clear();
	ak8_eta.clear();
	ak8_phi.clear();
	ak8_m.clear();
	ak8_softDrop_m.clear();
	ak8_trimmed_m.clear();
	ak8_tau1.clear();
	ak8_tau2.clear();
	ak8_tau3.clear();
	ak8_tau4.clear();
	ak8_tau5.clear();

	gen_mass.clear();
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

	Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(genEvtInfoToken,genEvtInfo);

  edm::Handle<std::vector<PileupSummaryInfo> > PUInfo;
    iEvent.getByToken(PUInfoToken,PUInfo);

//    iEvent.getByLabel( "generator","PdfInfoTag",genEvtInfo);
    

 //-------------- Gen Event Info -----------------------------------
//	if (!iEvent.isRealData())
//	{
//	cout << genEvtInfo->weight() << endl;
	weight_ = genEvtInfo->weight(); 

	for( std::vector<PileupSummaryInfo>::const_iterator it = PUInfo->begin(); it != PUInfo->end(); ++it )
	{
		npu_.push_back ( it -> getTrueNumInteractions() );
		Number_interactions.push_back ( it->getPU_NumInteractions() ); 
		OriginBX.push_back ( it -> getBunchCrossing());                
	
      }

/*           
		if( !genEvtInfo.isValid() ) edm::LogInfo("GenEvtInfo") << "ERROR: genEvtInfo not valid! " << genEvtInfo;
    
		if( genEvtInfo.isValid() )
		{
			edm::LogInfo("GenEvtInfo") << "Successfully obtained " << genEvtInfo;
	//      ptHat_ = (genEvtInfo->hasBinningValues() ? genEvtInfo->binningValues()[0] : -999.);
	//      processID_ = genEvtInfo->signalProcessID();
	weight_ = genEvtInfo->weight();            
	    }
*/
//	}    


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
		gen_mass.push_back ( gjet->mass() );
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
			float nergy_test = ijet->chargedHadronEnergy() /  jet_energy;
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


			if(nergy_test>1 || chf.back()>1)std::cout << "pT = " << ijet->pt() << ",  eta = " << ijet->eta() << ",  phi = " << ijet->phi() << ",  nhf = " << ijet->neutralHadronEnergy() <<  ",  chf = " << ijet->chargedHadronEnergy() <<   ",  nemf = " << ijet->photonEnergy() <<  ",  cemf = " << ijet->electronEnergy() <<  ",  muf = " << ijet->muonEnergy() <<  ",  CM = " << ijet->chargedHadronMultiplicity() + ijet->electronMultiplicity() + ijet->muonMultiplicity() <<  ",  NM = " << ijet->neutralHadronMultiplicity() + ijet->photonMultiplicity() << " chf =  "<<  chf.back()<< "energy test"<< nergy_test <<"\n" ;

		    float eta  = ijet->eta();
		    float pt   = ijet->pt();
		}

		

	//Save info about all Gen Particles
	for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) 
	{
		if (gens_iter->pdgId()==55 && gens_iter->numberOfDaughters()>1) 
		{
			genpart_pt.push_back(gens_iter->pt());
			genpart_eta.push_back(gens_iter->eta());
			genpart_phi.push_back(gens_iter->phi());
			genpart_m.push_back(gens_iter->mass());
			genpart_pdgid.push_back(gens_iter->pdgId());
			genpart_px.push_back(gens_iter->px());
			genpart_py.push_back(gens_iter->py());
			genpart_pz.push_back(gens_iter->pz());
			genpart_E.push_back(gens_iter->energy());
			genpart_p.push_back(gens_iter->p());
			genpart_numDaught.push_back(gens_iter->numberOfDaughters());		
			genpart_isQfromZ.push_back(false);

			
				
			int Ndaughters = gens_iter->numberOfDaughters();
		//	std::cout << "Found a Z'. It has  "<<Ndaughters<<" daughters. \n" ;
			for(int j = 0; j < Ndaughters; ++ j)
			{
	      			const Candidate * d = gens_iter->daughter( j );
	      			int dauId = d->pdgId();
				genpart_pt.push_back(d->pt());
				genpart_eta.push_back(d->eta());
				genpart_phi.push_back(d->phi());
				genpart_m.push_back(d->mass());
				genpart_pdgid.push_back(d->pdgId());
				genpart_px.push_back(d->px());
				genpart_py.push_back(d->py());
				genpart_pz.push_back(d->pz());
				genpart_E.push_back(d->energy());
				genpart_p.push_back(d->p());
				genpart_numDaught.push_back(d->numberOfDaughters());
				if(fabs(d->pdgId())<=6) genpart_isQfromZ.push_back(true);
	
			//	std::cout << "Daughter number "<<j<<"   ID = " << dauId<<"\n" ;
	       
	     		}	
		}
	}


	
	n_pfcand = 0;
	vector<PseudoJet> fj_part;
	for (auto pfcands_iter = pfcandsH->begin(); pfcands_iter != pfcandsH->end(); ++pfcands_iter)
	{

		pfcandpt[n_pfcand]=pfcands_iter->pt();
		pfcandeta[n_pfcand]=pfcands_iter->eta();
		pfcandphi[n_pfcand]= pfcands_iter->phi();
		pdcandm[n_pfcand]=pfcands_iter->m();
		pfcandpdgid[n_pfcand]=pfcands_iter->pdgId();
		pfcandvertex[n_pfcand]=pfcands_iter->vertex();
		n_pfcand++;
		
		PseudoJet temp_jet = PseudoJet(0, 0, 0, 0);
		temp_jet.reset_PtYPhiM(pfcands_iter->pt(), pfcands_iter->eta(), pfcands_iter->phi(), pfcands_iter->m());
		temp_jet.set_user_index(pfcands_iter->pdgId());
		fj_part.push_back(temp_jet);
	} 

	//AK8 stuff
	ClusterSequence ak8_cs(fj_part,JetDefinition(kt_algorithm, 0.8));
//	ClusterSequence ak8_cs(fj_part, ak8_def);
	vector<PseudoJet> ak8_jets = sorted_by_pt(ak8_cs.inclusive_jets(30.0)); // This is a Pt cut of 100. Change as required.

	//Groomer parameters from JET/MET , can change them as we like
	double sd_z_cut = 0.10;
    	double sd_beta = 0;
	SoftDrop sd_groomer = SoftDrop(sd_z_cut, sd_beta, 1.0);
	Filter trimmer = Filter(JetDefinition(kt_algorithm, 0.2), SelectorPtFractionMin(0.03));

	double beta = 1.0;
    	Nsubjettiness nSub1 = Nsubjettiness(1, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    	Nsubjettiness nSub2 = Nsubjettiness(2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    	Nsubjettiness nSub3 = Nsubjettiness(3, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    	Nsubjettiness nSub4 = Nsubjettiness(4, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
    	Nsubjettiness nSub5 = Nsubjettiness(5, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));

	for(auto &j: ak8_jets) {
		ak8_pt.push_back(j.pt());
        ak8_eta.push_back(j.pseudorapidity());
        ak8_phi.push_back(j.phi_std());
        ak8_m.push_back(j.m());

		PseudoJet trimmed_ak8 = trimmer(j);
        ak8_trimmed_m.push_back(trimmed_ak8.m());

		PseudoJet sd_ak8 = sd_groomer(j);
        ak8_softDrop_m.push_back(sd_ak8.m());

		ak8_tau1.push_back(nSub1.result(j));
        ak8_tau2.push_back(nSub2.result(j));
        ak8_tau3.push_back(nSub3.result(j));
        ak8_tau4.push_back(nSub4.result(j));
        ak8_tau5.push_back(nSub5.result(j));
	}

	//done with AK8 stuff

	




	
//	std::cout << VtxH->size() << "\n" ;
	nVtx = VtxH->size();
	SumEt = tot_energy; 
	gen_SumEt = gen_tot_energy;
	//		std::cout << "6 \n" ;


    //cout<"made it to Z"<<endl;
    
    //cout<"made it to after Z"<<endl;
    met =*metptH;
    met_phi = *metphiH;

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

