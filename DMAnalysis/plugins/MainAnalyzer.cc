// -*- C++ -*-
//
// Package:    DMAnalysis/MainAnalyzer
// Class:      MainAnalyzer
//
/**\class MainAnalyzer MainAnalyzer.cc llvvAnalysis/DMAnalysis/plugins/MainAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Renjie Wang
//         Created:  Wed, 28 Jan 2015 23:26:14 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


//
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"


//
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"


#include "llvvAnalysis/DMAnalysis/interface/SelectionMonitor.h"
#include "llvvAnalysis/DMAnalysis/interface/TSelectionMonitor.h"
#include "llvvAnalysis/DMAnalysis/interface/DataEvtSummaryHandler.h"


#include "DataFormats/Common/interface/ValueMap.h"
//#include "llvvAnalysis/DMAnalysis/interface/EGammaMvaEleEstimatorCSA14.h"


//
// class declaration
//
using namespace edm;
using namespace reco;
using namespace pat;
using namespace std;

class MainAnalyzer : public edm::EDAnalyzer {
public:
    explicit MainAnalyzer(const edm::ParameterSet&);
    ~MainAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);

    enum ElectronMatchType {UNMATCHED = 0,
                            TRUE_PROMPT_ELECTRON,
                            TRUE_ELECTRON_FROM_TAU,
                            TRUE_NON_PROMPT_ELECTRON
                           }; // The last does not include tau parents


private:

    edm::EDGetTokenT<reco::VertexCollection> vtxTag_;
    edm::EDGetTokenT<reco::BeamSpot> beamSpotTag_;
    edm::EDGetTokenT<pat::MuonCollection> muonTag_;
//    edm::EDGetTokenT<pat::ElectronCollection> electronTag_;
    edm::EDGetTokenT<edm::View<pat::Electron> > electronTag_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdTag_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdTag_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdTag_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdTag_;
    edm::EDGetTokenT<edm::ValueMap<bool> > electronHEEPIdTag_;
    edm::EDGetTokenT<pat::TauCollection> tauTag_;
    edm::EDGetTokenT<pat::PhotonCollection> photonTag_;
    edm::EDGetTokenT<pat::JetCollection> jetTag_;
    edm::EDGetTokenT<pat::JetCollection> jetPuppiTag_;
//    edm::EDGetTokenT<pat::JetCollection> fatjetTag_;
    edm::EDGetTokenT<pat::METCollection> metTag_;
    edm::EDGetTokenT<pat::METCollection> metNoHFTag_;
    edm::EDGetTokenT<pat::METCollection> metPuppiTag_;

    edm::EDGetTokenT<edm::TriggerResults> metFilterBitsTag_;

    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenTag_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoTag_;
    edm::EDGetTokenT<GenEventInfoProduct> genInfoTag_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenTag_;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genjetTag_;



    edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;


    std::vector<std::string> DoubleMuTrigs_, DoubleEleTrigs_, SingleMuTrigs_, SingleEleTrigs_, MuEGTrigs_;// DoubleTauTrigs_;
    DataEvtSummaryHandler summaryHandler_;
    TSelectionMonitor controlHistos_;
    bool isPythia8_;
    bool isMC_;
    bool verbose_;


    //edm::EDGetTokenT<double> rhoAllTag_;
    edm::EDGetTokenT<double> rhoFastjetAllTag_;
    //edm::EDGetTokenT<double> rhoFastjetAllCaloTag_;
    //edm::EDGetTokenT<double> rhoFastjetCentralCaloTag_;
    //edm::EDGetTokenT<double> rhoFastjetCentralChargedPileUpTag_;
    //edm::EDGetTokenT<double> rhoFastjetCentralNeutralTag_;


//    std::map<std::string, edm::ParameterSet> objConfig_;

    //MVAs for triggering and non-triggering electron ID
    //EGammaMvaEleEstimatorCSA14* myMVATrig;
    //EGammaMvaEleEstimatorCSA14* myMVANonTrig;



    // ID decisions objects
    edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapTokenTrig_;
    edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapTokenTrig_;

    // MVA values and categories (optional)
    edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapTokenTrig_;
    edm::EDGetTokenT<edm::ValueMap<int> >mvaCategoriesMapTokenTrig_;


    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    void getMCtruth(const edm::Event&, const edm::EventSetup&);
    bool checkIfTriggerFired(edm::Handle<edm::TriggerResults> &allTriggerBits, const edm::TriggerNames &triggerNames, std::string triggerPath);
    // MC truth matching utilities
    // The function that uses algorith from Josh Bendavid with
    // an explicit loop over gen particles.
    int matchToTruth(const pat::Electron &el, const edm::Handle<edm::View<reco::GenParticle>>  &prunedGenParticles);
    void findFirstNonElectronMother(const reco::Candidate *particle, int &ancestorPID, int &ancestorStatus);

    //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
    //virtual void beginLuminosityBlock(edm::LuminosityBlock &iLumi, edm::EventSetup &iSetup);
    //virtual void endLuminosityBlock(edm::LuminosityBlock &iLumi, edm::EventSetup &iSetup);

    // ----------member data ---------------------------
    float curAvgInstLumi_;
    float curIntegLumi_;

};


//
// constants, enums and typedefs
//
// Effective areas for electrons from Giovanni P. and Cristina
// distributed as private slides in Jan 2015, derived for PHYS14
namespace EffectiveAreas {
const int nEtaBins = 5;
const float etaBinLimits[nEtaBins+1] = {
    0.0, 0.8, 1.3, 2.0, 2.2, 2.5
};
const float effectiveAreaValues[nEtaBins] = {
    0.1013, 0.0988, 0.0572, 0.0842, 0.1530
};
}

//
// static data member definitions
//


//
// constructors and destructor
//
MainAnalyzer::MainAnalyzer(const edm::ParameterSet& iConfig):
    vtxTag_(		consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("verticesTag"))		),
    beamSpotTag_(	consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))			),
    muonTag_(		consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonsTag"))			),
//    electronTag_(	consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronsTag"))		),
    electronTag_(	consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronsTag"))	),
    electronVetoIdTag_(	consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronVetoIdTag"))	),
    electronLooseIdTag_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronLooseIdTag"))       ),
    electronMediumIdTag_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronMediumIdTag"))     ),
    electronTightIdTag_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronTightIdTag"))	),
    electronHEEPIdTag_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronHEEPIdTag"))         ),
    tauTag_(		consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tausTag"))			),
    photonTag_(		consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photonsTag"))		),
    jetTag_(		consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsTag"))			),
    jetPuppiTag_(       consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsPuppiTag"))               ),
//    fatjetTag_(		consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjetsTag"))			),
    metTag_(		consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsTag"))			),
    metNoHFTag_(        consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsNoHFTag"))                ),
    metPuppiTag_(       consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metsPuppiTag"))               ),
    metFilterBitsTag_(	consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("metFilterBitsTag"))		),
    prunedGenTag_(	consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("prunedTag"))	),
    puInfoTag_(         consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puInfoTag"))     ),
    genInfoTag_(        consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfoTag"))                ),
    packedGenTag_(	consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedTag"))	),
    genjetTag_(		consumes<edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJetsTag"))		),
    triggerBits_(	consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))			),
    triggerObjects_(	consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
    triggerPrescales_(	consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))		),
    DoubleMuTrigs_(	iConfig.getParameter<std::vector<std::string> >("DoubleMuTrigs")				),
    DoubleEleTrigs_(    iConfig.getParameter<std::vector<std::string> >("DoubleEleTrigs")                               ),
    SingleMuTrigs_(	iConfig.getParameter<std::vector<std::string> >("SingleMuTrigs")				),
    SingleEleTrigs_(    iConfig.getParameter<std::vector<std::string> >("SingleEleTrigs")				),
    MuEGTrigs_(		iConfig.getParameter<std::vector<std::string> >("MuEGTrigs")					),
//    DoubleTauTrigs_(     iConfig.getParameter<std::vector<std::string> >("DoubleTauTrigs")                                ),
    controlHistos_(	iConfig.getParameter<std::string>("dtag")							),
    isPythia8_(		iConfig.getParameter<bool>("isPythia8")								),
    isMC_(		iConfig.getParameter<bool>("isMC")								),
    verbose_(		iConfig.getParameter<bool>("verbose")								),
    //rhoAllTag_(			consumes<double>(iConfig.getParameter<edm::InputTag>("rhoAll"))				),
    rhoFastjetAllTag_(  	consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAll")) 			),
    //rhoFastjetAllCaloTag_( 	consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetAllCalo")) 		),
    //rhoFastjetCentralCaloTag_(  consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralCalo")) 		),
    //rhoFastjetCentralChargedPileUpTag_(  consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralChargedPileUp")) ),
    //rhoFastjetCentralNeutralTag_(  	consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastjetCentralNeutral"))	 ),
    eleMediumIdMapTokenTrig_(	consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMapTrig"))	),
    eleTightIdMapTokenTrig_(	consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMapTrig"))	),
    mvaValuesMapTokenTrig_(	consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMapTrig"))	),
    mvaCategoriesMapTokenTrig_(	consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMapTrig"))	),
    curAvgInstLumi_(0),
    curIntegLumi_(0)


{

//    std::string objs[]= {"triggerPaths"};
//    for(size_t iobj=0; iobj<sizeof(objs)/sizeof(string); iobj++)
//        objConfig_[ objs[iobj] ] = iConfig.getParameter<edm::ParameterSet>( objs[iobj] );

    consumesMany<LHEEventProduct>();

    //now do what ever initialization is needed
    edm::Service<TFileService> fs;
    summaryHandler_.initTree(  fs->make<TTree>("data","Event Summary") );
    TFileDirectory baseDir=fs->mkdir(iConfig.getParameter<std::string>("dtag"));

    //cut flow histograms for book keeping
//    TString selFilters[]= {"Reco","no scrap","#geq 1 vertex","HB/HE noise","No beam halo"};
//    const size_t nselFilters=sizeof(selFilters)/sizeof(TString);
//    controlHistos_.addHistogram("cutflow", ";Steps; Events", nselFilters, 0.,nselFilters);
//    TH1 *h = controlHistos_.getHisto("cutflow");
//    for(size_t istep=0; istep<nselFilters; istep++) h->GetXaxis()->SetBinLabel(istep+1,selFilters[istep]);

    controlHistos_.addHistogram("nevents",";nevents; nevents",1,-0.5,0.5);
    controlHistos_.addHistogram("n_negevents",";n_negevents; n_negevents",1,-0.5,0.5);
    controlHistos_.addHistogram("n_posevents",";n_posevents; n_posevents",1,-0.5,0.5);
    controlHistos_.addHistogram("pileup", ";Pileup; Events",100,-0.5,99.5);
    //controlHistos_.addHistogram("integlumi", ";Integrated luminosity ; Events",100,0,1e5);
    //controlHistos_.addHistogram("instlumi", ";Max average inst. luminosity; Events",100,0,1e5);
    controlHistos_.addHistogram("pileuptrue", ";True pileup; Events",100,-0.5,99.5);


    //set up electron MVA ID
    /*
        //twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        //twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2
        //twiki.cern.ch/twiki/bin/viewauth/CMS/HEEPElectronIdentificationRun2
        std::vector<std::string> myTrigWeights;
    //    myTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/TrigIDMVA_25ns_EB_BDT.weights.xml").fullPath().c_str());
    //    myTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/TrigIDMVA_25ns_EE_BDT.weights.xml").fullPath().c_str());
        myTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/TrigIDMVA_50ns_EB_BDT.weights.xml").fullPath().c_str());
        myTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/TrigIDMVA_50ns_EE_BDT.weights.xml").fullPath().c_str());


        myMVATrig = new EGammaMvaEleEstimatorCSA14();
        myMVATrig->initialize("BDT",
                              EGammaMvaEleEstimatorCSA14::kTrig,
                              true,
                              myTrigWeights);

        std::vector<std::string> myNonTrigWeights;
    //    myNonTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/EIDmva_EB_5_25ns_BDT.weights.xml").fullPath().c_str());
    //    myNonTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/EIDmva_EE_5_25ns_BDT.weights.xml").fullPath().c_str());
    //    myNonTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/EIDmva_EB_10_25ns_BDT.weights.xml").fullPath().c_str());
    //    myNonTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/EIDmva_EE_10_25ns_BDT.weights.xml").fullPath().c_str());
        myNonTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/EIDmva_EB_5_50ns_BDT.weights.xml").fullPath().c_str());
        myNonTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/EIDmva_EE_5_50ns_BDT.weights.xml").fullPath().c_str());
        myNonTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/EIDmva_EB_10_50ns_BDT.weights.xml").fullPath().c_str());
        myNonTrigWeights.push_back(edm::FileInPath("llvvAnalysis/DMAnalysis/data/CSA14/EIDmva_EE_10_50ns_BDT.weights.xml").fullPath().c_str());


        myMVANonTrig = new EGammaMvaEleEstimatorCSA14();
        myMVANonTrig->initialize("BDT",
                                 EGammaMvaEleEstimatorCSA14::kNonTrig,
                                 true,
                                 myNonTrigWeights);

    */


}


MainAnalyzer::~MainAnalyzer()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MainAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& iSetup)
{
    controlHistos_.fillHisto("nevents","all",0); //increment event count

    summaryHandler_.resetStruct();
    //event summary to be filled
    DataEvtSummary_t &ev=summaryHandler_.getEvent();

    //event header
    ev.run    = event.id().run();
    ev.lumi   = event.luminosityBlock();
    ev.event  = event.id().event();
    ev.curAvgInstLumi=curAvgInstLumi_;
    ev.curIntegLumi=curIntegLumi_;

    // Pruned particles are the one containing "important" stuff
    edm::Handle<edm::View<reco::GenParticle> > pruned;
    //event.getByToken(prunedGenTag_,pruned);


    //MC truth
    if(isMC_) {
        event.getByToken(prunedGenTag_,pruned);
        getMCtruth(event, iSetup);
    }

    //
    // trigger
    //

    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

    event.getByToken(triggerBits_, triggerBits);
    event.getByToken(triggerObjects_, triggerObjects);
    event.getByToken(triggerPrescales_, triggerPrescales);

    bool hasDoubleMuTrigs(false);
    bool hasDoubleEleTrigs(false);
    bool hasSingleMuTrigs(false);
    bool hasSingleEleTrigs(false);
    bool hasMuEGTrigs(false);
//    bool hasDoubleTauTrigs(false);

    const edm::TriggerNames &names = event.triggerNames(*triggerBits);

    if (verbose_) {
        std::cout << "\n === TRIGGER PATHS === " << std::endl;
        for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
            std::cout << "Trigger " << names.triggerName(i) <<
                      ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
                      ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
                      << std::endl;
        }
    }

    for(size_t it=0; it<DoubleMuTrigs_.size(); it++) {
        hasDoubleMuTrigs |= checkIfTriggerFired(triggerBits, names, DoubleMuTrigs_[it]);
    }
    for(size_t it=0; it<DoubleEleTrigs_.size(); it++) {
        hasDoubleEleTrigs |= checkIfTriggerFired(triggerBits, names, DoubleEleTrigs_[it]);
    }
    for(size_t it=0; it<SingleMuTrigs_.size(); it++) {
        hasSingleMuTrigs |= checkIfTriggerFired(triggerBits, names, SingleMuTrigs_[it]);
    }
    for(size_t it=0; it<SingleEleTrigs_.size(); it++) {
        hasSingleEleTrigs |= checkIfTriggerFired(triggerBits, names, SingleEleTrigs_[it]);
    }
    for(size_t it=0; it<MuEGTrigs_.size(); it++) {
        hasMuEGTrigs |= checkIfTriggerFired(triggerBits, names, MuEGTrigs_[it]);
    }
//    for(size_t it=0; it<DoubleTauTrigs_.size(); it++) {
//        hasDoubleTauTrigs |= checkIfTriggerFired(triggerBits, names, DoubleTauTrigs_[it]);
//    }

    ev.hasTrigger = (hasDoubleMuTrigs || hasDoubleEleTrigs || hasSingleMuTrigs || hasSingleEleTrigs || hasMuEGTrigs);//|| hasDoubleTauTrigs);

    ev.triggerType = ( hasDoubleMuTrigs  << 0 )
                     | ( hasSingleMuTrigs  << 1 )
                     | ( hasDoubleEleTrigs << 2 )
                     | ( hasSingleEleTrigs << 3 )
                     | ( hasMuEGTrigs	 << 4 );
    //| ( hasDoubleTauTrigs << 5 );

    //if(!isMC_ && !ev.hasTrigger) return; // skip the event if no trigger, only for Data
    if(!ev.hasTrigger) return; // skip the event if no trigger, for both Data and MC



    //
    // vertex and beam spot
    //

    edm::Handle<reco::BeamSpot> beamSpot;
    event.getByToken(beamSpotTag_, beamSpot);

    edm::Handle<reco::VertexCollection> vertices;
    event.getByToken(vtxTag_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

    ev.vtx_x = PV.x();
    ev.vtx_y = PV.y();
    ev.vtx_z = PV.z();

    ev.nvtx = 0;
    //select good vertices
    for(unsigned int i = 0; i < vertices->size(); i++) {
        if(vertices->at(i).isValid() && !vertices->at(i).isFake()) ev.nvtx++;
    }
    if(ev.nvtx == 0) return;


    //edm::Handle<double> rhoAll;
    edm::Handle<double> rhoFastjetAll;
    //edm::Handle<double> rhoFastjetAllCalo;
    //edm::Handle<double> rhoFastjetCentralCalo;
    //edm::Handle<double> rhoFastjetCentralChargedPileUp;
    //edm::Handle<double> rhoFastjetCentralNeutral;

    //event.getByToken(rhoAllTag_,rhoAll);
    event.getByToken(rhoFastjetAllTag_,rhoFastjetAll);
    //event.getByToken(rhoFastjetAllCaloTag_,rhoFastjetAllCalo);
    //event.getByToken(rhoFastjetCentralCaloTag_,rhoFastjetCentralCalo);
    //event.getByToken(rhoFastjetCentralChargedPileUpTag_,rhoFastjetCentralChargedPileUp);
    //event.getByToken(rhoFastjetCentralNeutralTag_,rhoFastjetCentralNeutral);

    //get rho
    //ev.fixedGridRhoAll = *rhoAll;
    float fixedGridRhoFastjetAll = *rhoFastjetAll;
    //ev.fixedGridRhoFastjetAllCalo = *rhoFastjetAllCalo;
    //ev.fixedGridRhoFastjetCentralCalo = *rhoFastjetCentralCalo;
    //ev.fixedGridRhoFastjetCentralChargedPileUp = *rhoFastjetCentralChargedPileUp;
    //ev.fixedGridRhoFastjetCentralNeutral = *rhoFastjetCentralNeutral;



    //met filters
    edm::Handle<edm::TriggerResults> metFilterBits;
    event.getByToken(metFilterBitsTag_, metFilterBits);
    const edm::TriggerNames &metNames = event.triggerNames(*metFilterBits);
    bool passMETFilters(true);
    for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i) {
        if(strcmp(metNames.triggerName(i).c_str(), "Flag_goodVertices") == 0)
            passMETFilters &= metFilterBits->accept(i);
        else if(strcmp(metNames.triggerName(i).c_str(), "Flag_eeBadScFilter") == 0)
            passMETFilters &= metFilterBits->accept(i);
    }
    if(!passMETFilters) return;



    //
    // muon selection
    //

    edm::Handle<pat::MuonCollection> muons;
    event.getByToken(muonTag_, muons);
    ev.mn=0;
    for(const pat::Muon &mu : *muons) {
        if(mu.pt() < 3) continue;
        ev.mn_px[ev.mn] = mu.px();
        ev.mn_py[ev.mn] = mu.py();
        ev.mn_pz[ev.mn] = mu.pz();
        ev.mn_en[ev.mn] = mu.energy();
        ev.mn_id[ev.mn] = 13*mu.charge();

        ev.mn_d0[ev.mn] = -mu.muonBestTrack()->dxy(PV.position());
        ev.mn_dZ[ev.mn] = mu.muonBestTrack()->dz(PV.position());
        ev.mn_ip3d[ev.mn] = mu.dB(pat::Muon::PV3D);
        ev.mn_ip3dsig[ev.mn] = mu.dB(pat::Muon::PV3D)/mu.edB(pat::Muon::PV3D);

        ev.mn_IsLoose[ev.mn] = mu.isLooseMuon();
        ev.mn_IsMedium[ev.mn] = mu.isMediumMuon();
        ev.mn_IsTight[ev.mn] = mu.isTightMuon(PV);
        ev.mn_IsSoft[ev.mn] = mu.isSoftMuon(PV);
        ev.mn_IsHighPt[ev.mn] = mu.isHighPtMuon(PV);

        ev.mn_pileupIsoR03[ev.mn] = mu.pfIsolationR03().sumPUPt;
        ev.mn_chargedIsoR03[ev.mn] = mu.pfIsolationR03().sumChargedHadronPt;
        ev.mn_photonIsoR03[ev.mn] = mu.pfIsolationR03().sumPhotonEt;
        ev.mn_neutralHadIsoR03[ev.mn] = mu.pfIsolationR03().sumNeutralHadronEt;

        ev.mn_pileupIsoR04[ev.mn] = mu.pfIsolationR04().sumPUPt;
        ev.mn_chargedIsoR04[ev.mn] = mu.pfIsolationR04().sumChargedHadronPt;
        ev.mn_photonIsoR04[ev.mn] = mu.pfIsolationR04().sumPhotonEt;
        ev.mn_neutralHadIsoR04[ev.mn] = mu.pfIsolationR04().sumNeutralHadronEt;

        //ev.mn_nMatches[ev.mn]                   = mu.numberOfMatches();
        //ev.mn_nMatchedStations[ev.mn]           = mu.numberOfMatchedStations();
        //ev.mn_validMuonHits[ev.mn]              = mu.isGlobalMuon() ? mu.globalTrack().hitPattern().numberOfValidMuonHits() : 0.;
        //ev.mn_innerTrackChi2[ev.mn]             = mu.isTrackerMuon() ? mu.innerTrack().normalizedChi2() : 0.;
        //ev.mn_trkLayersWithMeasurement[ev.mn]   = mu.track().hitPattern().trackerLayersWithMeasurement();
        //ev.mn_pixelLayersWithMeasurement[ev.mn] = mu.isTrackerMuon() ? mu.innerTrack().hitPattern().pixelLayersWithMeasurement() : 0.;

        ev.mn_type[ev.mn]   = (	mu.isMuon() 		<< 0)
                              | (	mu.isGlobalMuon() 	<< 1)
                              | (	mu.isTrackerMuon() 	<< 2)
                              | (	mu.isStandAloneMuon()	<< 3)
                              | (	mu.isCaloMuon() 	<< 4)
                              | (	mu.isPFMuon() 		<< 5)
                              | (	mu.isRPCMuon()		<< 6);
        ev.mn++;
    }


    //
    // electron selection
    //

    //edm::Handle<pat::ElectronCollection> electrons;
    Handle<edm::View<pat::Electron> > electrons;
    event.getByToken(electronTag_, electrons);

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Electron_ID_Working_Points_WP_de
    // Get the electron ID data from the event stream.
    // Note: this implies that the VID ID modules have been run upstream.
    // If you need more info, check with the EGM group.
    edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
    edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
    edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
    edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
    edm::Handle<edm::ValueMap<bool> > heep_id_decisions;

    event.getByToken(electronVetoIdTag_,veto_id_decisions);
    event.getByToken(electronLooseIdTag_,loose_id_decisions);
    event.getByToken(electronMediumIdTag_,medium_id_decisions);
    event.getByToken(electronTightIdTag_,tight_id_decisions);
    event.getByToken(electronHEEPIdTag_,heep_id_decisions);


    // Get the electron ID data from the event stream.
    // Note: this implies that the VID ID modules have been run upstream.
    // If you need more info, check with the EGM group.
    edm::Handle<edm::ValueMap<bool> > mvaTrig_medium_id_decisions;
    edm::Handle<edm::ValueMap<bool> > mvaTrig_tight_id_decisions;
    event.getByToken(eleMediumIdMapTokenTrig_,mvaTrig_medium_id_decisions);
    event.getByToken(eleTightIdMapTokenTrig_,mvaTrig_tight_id_decisions);

    // Get MVA values and categories (optional)
    edm::Handle<edm::ValueMap<float> > mvaTrigValues;
    edm::Handle<edm::ValueMap<int> > mvaTrigCategories;
    event.getByToken(mvaValuesMapTokenTrig_,mvaTrigValues);
    event.getByToken(mvaCategoriesMapTokenTrig_,mvaTrigCategories);


    ev.en=0;
    //for (const pat::Electron &el : *electrons) {
    for( View<pat::Electron>::const_iterator el = electrons->begin(); el != electrons->end(); el++ ) {
        float pt_ = el->pt();
        if (pt_ < 5) continue;

        // Kinematics
        ev.en_px[ev.en] = el->px();
        ev.en_py[ev.en] = el->py();
        ev.en_pz[ev.en] = el->pz();
        ev.en_en[ev.en] = el->energy();
        ev.en_id[ev.en] = 11*el->charge();

        /*
                ev.en_EtaSC[ev.en] = el->superCluster()->eta();
                ev.en_PhiSC[ev.en] = el->superCluster()->phi();
                ev.en_EnSC[ev.en] = el->superCluster()->energy();

                // ID and matching
                ev.en_dEtaIn[ev.en] = el->deltaEtaSuperClusterTrackAtVtx();
                ev.en_dPhiIn[ev.en] = el->deltaPhiSuperClusterTrackAtVtx();
                ev.en_hOverE[ev.en] = el->hcalOverEcal();
                ev.en_R9[ev.en] = el->r9();
                ev.en_sigmaIetaIeta[ev.en] = el->sigmaIetaIeta();
                ev.en_sigmaIetaIeta5x5[ev.en] = el->full5x5_sigmaIetaIeta();

                // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
                // The if protects against ecalEnergy == inf or zero (always
                // the case for electrons below 5 GeV in miniAOD)
                if( el->ecalEnergy() == 0 ) {
                    //printf("Electron energy is zero!\n");
                    ev.en_ooEmooP[ev.en] = 1e30;
                } else if( !std::isfinite(el->ecalEnergy())) {
                    //printf("Electron energy is not finite!\n");
                    ev.en_ooEmooP[ev.en] = 1e30;
                } else {
                    ev.en_ooEmooP[ev.en] = fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() );
                }

                // Impact parameter
                ev.en_d0[ev.en] = (-1) * el->gsfTrack()->dxy(PV.position());
                ev.en_dZ[ev.en] = el->gsfTrack()->dz(PV.position());
        */
        //Isolation
        GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
        ev.en_pileupIso[ev.en] = pfIso.sumPUPt;
        ev.en_chargedIso[ev.en] = pfIso.sumChargedHadronPt;
        ev.en_photonIso[ev.en] = pfIso.sumPhotonEt;
        ev.en_neutralHadIso[ev.en] = pfIso.sumNeutralHadronEt;

        // Compute isolation with effective area correction for PU
        // Find eta bin first. If eta>2.5, the last eta bin is used.
        int etaBin = 0;
        float etaSC_ = el->superCluster()->eta();
        while ( etaBin < EffectiveAreas::nEtaBins-1 && abs(etaSC_) > EffectiveAreas::etaBinLimits[etaBin+1] ) {
            ++etaBin;
        };

        double area = EffectiveAreas::effectiveAreaValues[etaBin];
        double rho_ = fixedGridRhoFastjetAll;
        ev.en_relIsoWithEA[ev.en] = ( pfIso.sumChargedHadronPt + max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho_ * area ) )/pt_;

        // Compute isolation with delta beta correction for PU
        float absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
        ev.en_relIsoWithDBeta[ev.en] = absiso/pt_;

        // this part should move to DMPhysicsEvent.h

        // Conversion rejection
        ev.en_MissingHits[ev.en] = el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
        ev.en_passConversionVeto[ev.en] = el->passConversionVeto();


        // Look up the ID decision for this electron in
        // the ValueMap object and store it. We need a Ptr object as the key.
        const Ptr<pat::Electron> elPtr(electrons, el - electrons->begin() );
        ev.en_passVeto[ev.en]  = (*veto_id_decisions)[ elPtr ];
        ev.en_passLoose[ev.en] = (*loose_id_decisions)[ elPtr ];
        ev.en_passMedium[ev.en]= (*medium_id_decisions)[ elPtr ];
        ev.en_passTight[ev.en] = (*tight_id_decisions)[ elPtr ];
        ev.en_passHEEP[ev.en]  = (*heep_id_decisions) [ elPtr ];

        //ID MVA
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2
        ev.en_passMVATrigMedium[ev.en] 	= (*mvaTrig_medium_id_decisions)[ elPtr ];
        ev.en_passMVATrigTight[ev.en] 	= (*mvaTrig_tight_id_decisions) [ elPtr ];
        ev.en_IDMVATrigValue[ev.en] 	= (*mvaTrigValues)	        [ elPtr ];
        ev.en_IDMVATrigCategory[ev.en] 	= (*mvaTrigCategories)		[ elPtr ];

        //ev.en_IDMVATrig[ev.en] = myMVATrig->mvaValue(*el,false);
        //ev.en_IDMVANonTrig[ev.en] = myMVANonTrig->mvaValue(*el,false);

        //
        // Explicit loop over gen candidates method
        //
        if(isMC_) ev.en_istrue[ev.en] = matchToTruth( *el, pruned);
        else ev.en_istrue[ev.en] = 0;

        ev.en++;
    }


    //
    // tau selection
    //
    edm::Handle<pat::TauCollection> taus;
    event.getByToken(tauTag_, taus);
    ev.ta=0;
    for (const pat::Tau &tau : *taus) {
        if(tau.pt() < 20) continue;

        ev.ta_px[ev.ta] = tau.px();
        ev.ta_py[ev.ta] = tau.py();
        ev.ta_pz[ev.ta] = tau.pz();
        ev.ta_en[ev.ta] = tau.energy();
        ev.ta_id[ev.ta] = 15*tau.charge();

        //Decay Mode Reconstruction
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
        //ev.ta_dm[ev.ta] = bool(tau.tauID("decayModeFindingOldDMs"));
        ev.ta_dm[ev.ta] = bool(tau.tauID("decayModeFinding"));
        ev.ta_newdm[ev.ta] = bool(tau.tauID("decayModeFindingNewDMs"));

        //Isolation
        ev.ta_IsLooseIso[ev.ta] = bool(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
        ev.ta_IsMediumIso[ev.ta] = bool(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
        ev.ta_IsTightIso[ev.ta] = bool(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
        ev.ta_combIsoDBeta3Hits[ev.ta] = tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        ev.ta_chargedIso[ev.ta] = tau.tauID("chargedIsoPtSum");
        ev.ta_neutralIso[ev.ta] = tau.tauID("neutralIsoPtSum");
        ev.ta_pileupIso[ev.ta]  = tau.tauID("puCorrPtSum");

        //Electron Rejection
        ev.ta_passEleVetoLoose[ev.ta] = bool(tau.tauID("againstElectronLooseMVA5"));
        ev.ta_passEleVetoMedium[ev.ta] = bool(tau.tauID("againstElectronMediumMVA5"));
        ev.ta_passEleVetoTight[ev.ta] = bool(tau.tauID("againstElectronTightMVA5"));

        //Muon Rejection
        ev.ta_passMuVetoLoose3[ev.ta] = bool(tau.tauID("againstMuonLoose3"));
        ev.ta_passMuVetoTight3[ev.ta] = bool(tau.tauID("againstMuonTight3") );

        ev.ta++;
    }


    //
    // jet selection (ak4PFJetsCHS)
    //
    edm::Handle<pat::JetCollection> jets;
    event.getByToken(jetTag_, jets);
    ev.jet=0;
    PFJetIDSelectionFunctor looseJetIdSelector(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::LOOSE);
    PFJetIDSelectionFunctor tightJetIdSelector(PFJetIDSelectionFunctor::FIRSTDATA,PFJetIDSelectionFunctor::TIGHT);
    pat::strbitset hasLooseId = looseJetIdSelector.getBitTemplate();
    pat::strbitset hasTightId = tightJetIdSelector.getBitTemplate();

    for (const pat::Jet &j : *jets) {
        if(j.pt() < 20) continue;

        //jet id
        hasLooseId.set(false);
        hasTightId.set(false);
        bool passLooseId(looseJetIdSelector( j, hasLooseId ));
        bool passTightId(tightJetIdSelector( j, hasTightId ));
        ev.jet_PFLoose[ev.jet] = passLooseId;
        ev.jet_PFTight[ev.jet] = passTightId;

        ev.jet_px[ev.jet] = j.px(); //correctedP4(0).px();
        ev.jet_py[ev.jet] = j.py(); //correctedP4(0).py();
        ev.jet_pz[ev.jet] = j.pz(); //correctedP4(0).pz();
        ev.jet_en[ev.jet] = j.energy(); //correctedP4(0).energy();

        ev.jet_btag0[ev.jet] = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        ev.jet_btag1[ev.jet] = j.bDiscriminator("pfJetBProbabilityBJetTags");
        ev.jet_btag2[ev.jet] = j.bDiscriminator("pfJetProbabilityBJetTags");
        ev.jet_btag3[ev.jet] = j.bDiscriminator("pfTrackCountingHighPurBJetTags");
        ev.jet_btag4[ev.jet] = j.bDiscriminator("pfTrackCountingHighEffBJetTags");
        ev.jet_btag5[ev.jet] = j.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags");
        ev.jet_btag6[ev.jet] = j.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags");
        ev.jet_btag7[ev.jet] = j.bDiscriminator("combinedSecondaryVertexBJetTags");
        ev.jet_btag8[ev.jet] = j.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
        ev.jet_btag9[ev.jet] = j.bDiscriminator("pfCombinedSecondaryVertexSoftLeptonBJetTags");
        ev.jet_btag10[ev.jet] = j.bDiscriminator("pfCombinedMVABJetTags");

        ev.jet_mass[ev.jet] = j.mass(); //correctedP4(0).M();
        ev.jet_area[ev.jet] = j.jetArea();
        ev.jet_pu[ev.jet] = j.pileup();
        ev.jet_puId[ev.jet] = j.userFloat("pileupJetId:fullDiscriminant");
        ev.jet_partonFlavour[ev.jet] = j.partonFlavour();
        ev.jet_hadronFlavour[ev.jet] = j.hadronFlavour();
        const reco::GenJet *gJet=j.genJet();
        if(gJet) ev.jet_genpt[ev.jet] = gJet->pt();
        else     ev.jet_genpt[ev.jet] = 0;

        ev.jet++;
    }


    //
    // slimmedJetsPuppi
    //
    edm::Handle<pat::JetCollection> puppijets;
    event.getByToken(jetPuppiTag_, puppijets);

    ev.pjet=0;
    if(puppijets.isValid()) {
        for (const pat::Jet &j : *puppijets) {
            if(j.pt() < 20) continue;
            ev.pjet_px[ev.pjet] = j.px(); //correctedP4(0).px();
            ev.pjet_py[ev.pjet] = j.px(); //correctedP4(0).py();
            ev.pjet_pz[ev.pjet] = j.pz(); //correctedP4(0).pz();
            ev.pjet_en[ev.pjet] = j.energy(); //correctedP4(0).energy();

            const reco::GenJet *gJet=j.genJet();
            if(gJet) ev.pjet_genpt[ev.pjet] = gJet->pt();
            else     ev.pjet_genpt[ev.pjet] = 0;

            ev.pjet_btag0[ev.pjet] = j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
            ev.pjet_btag1[ev.pjet] = j.bDiscriminator("pfJetBProbabilityBJetTags");
            ev.pjet_btag2[ev.pjet] = j.bDiscriminator("pfJetProbabilityBJetTags");
            ev.pjet_btag3[ev.pjet] = j.bDiscriminator("pfTrackCountingHighPurBJetTags");
            ev.pjet_btag4[ev.pjet] = j.bDiscriminator("pfTrackCountingHighEffBJetTags");
            ev.pjet_btag5[ev.pjet] = j.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags");
            ev.pjet_btag6[ev.pjet] = j.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags");
            ev.pjet_btag7[ev.pjet] = j.bDiscriminator("combinedSecondaryVertexBJetTags");
            ev.pjet_btag8[ev.pjet] = j.bDiscriminator("pfCombinedSecondaryVertexV2BJetTags");
            ev.pjet_btag9[ev.pjet] = j.bDiscriminator("pfCombinedSecondaryVertexSoftLeptonBJetTags");
            ev.pjet_btag10[ev.pjet] = j.bDiscriminator("pfCombinedMVABJetTags");
            ev.pjet++;
        }
    }

    //
    // fat jet selection (ak8PFJetsCHS)
    //
    /*
        edm::Handle<pat::JetCollection> fatjets;
        event.getByToken(fatjetTag_, fatjets);
        ev.fjet=0;
        for (const pat::Jet &j : *fatjets) {
            ev.fjet_px[ev.fjet] = j.correctedP4(0).px();
            ev.fjet_py[ev.fjet] = j.correctedP4(0).py();
            ev.fjet_pz[ev.fjet] = j.correctedP4(0).pz();
            ev.fjet_en[ev.fjet] = j.correctedP4(0).energy();

            const reco::GenJet *gJet=j.genJet();
            if(gJet) ev.fjet_genpt[ev.fjet] = gJet->pt();
            else     ev.fjet_genpt[ev.fjet] = 0;


            ev.fjet_prunedM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSPrunedLinks");
            ev.fjet_trimmedM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSTrimmedLinks");
            ev.fjet_filteredM[ev.fjet] = (float) j.userFloat("ak8PFJetsCHSFilteredLinks");
            ev.fjet_tau1[ev.fjet] =  (float) j.userFloat("NjettinessAK8:tau1");
            ev.fjet_tau2[ev.fjet] =  (float) j.userFloat("NjettinessAK8:tau2");
            ev.fjet_tau3[ev.fjet] =  (float) j.userFloat("NjettinessAK8:tau3");

            ev.fjet++;
        }
    */



    //
    // met selection
    //
    edm::Handle<pat::METCollection> mets;
    event.getByToken(metTag_, mets);
    const pat::MET &met = mets->front();

    //PF type-1 ETmiss
    ev.met_pt = met.pt();
    ev.met_phi = met.phi();
    ev.met_sumMET = met.sumEt();

    // raw PF ETmiss
    ev.rawpfmet_pt = met.uncorPt();
    ev.rawpfmet_phi = met.uncorPhi();
    ev.rawpfmet_sumMET = met.uncorSumEt();

    // raw calo ETmiss
    ev.rawcalomet_pt = met.caloMETPt();
    ev.rawcalomet_phi = met.caloMETPhi();
    ev.rawcalomet_sumMET = met.caloMETSumEt();

    // type1 PF MET but excluding HF
    edm::Handle<pat::METCollection> metsNoHF;
    event.getByToken(metNoHFTag_, metsNoHF);
    if(metsNoHF.isValid()) {
        const pat::MET &metNoHF = metsNoHF->front();
        ev.metNoHF_pt = metNoHF.pt();
        ev.metNoHF_phi = metNoHF.phi();
        ev.metNoHF_sumMET = metNoHF.sumEt();
    }

    // puppi-corrected MET
    edm::Handle<pat::METCollection> metsPuppi;
    event.getByToken(metPuppiTag_, metsPuppi);
    if(metsPuppi.isValid()) {
        const pat::MET &metPuppi = metsPuppi->front();
        ev.metPuppi_pt = metPuppi.pt();
        ev.metPuppi_phi = metPuppi.phi();
        ev.metPuppi_sumMET = metPuppi.sumEt();
    }



    /*
        //met filters
        edm::Handle<edm::TriggerResults> metFilterBits;
        event.getByToken(metFilterBitsTag_, metFilterBits);
        const edm::TriggerNames &metNames = event.triggerNames(*metFilterBits);
        for(unsigned int i = 0, n = metFilterBits->size(); i < n; ++i) {
            if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseFilter") == 0)
                ev.flag_HBHENoiseFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseIsoFilter") == 0)
                ev.flag_HBHENoiseIsoFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_CSCTightHaloFilter") == 0)
                ev.flag_CSCTightHaloFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_hcalLaserEventFilter") == 0)
                ev.flag_hcalLaserEventFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellTriggerPrimitiveFilter") == 0)
                ev.flag_EcalDeadCellTriggerPrimitiveFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellBoundaryEnergyFilter") == 0)
                ev.flag_EcalDeadCellBoundaryEnergyFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_goodVertices") == 0)
                ev.flag_goodVertices = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trackingFailureFilter") == 0)
                ev.flag_trackingFailureFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_eeBadScFilter") == 0)
                ev.flag_eeBadScFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_ecalLaserCorrFilter") == 0)
                ev.flag_ecalLaserCorrFilter = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOGFilters") == 0)
                ev.flag_trkPOGFilters = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_manystripclus53X") == 0)
                ev.flag_trkPOG_manystripclus53X = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_toomanystripclus53X") == 0)
                ev.flag_trkPOG_toomanystripclus53X = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_trkPOG_logErrorTooManyClusters") == 0)
                ev.flag_trkPOG_logErrorTooManyClusters = metFilterBits->accept(i);
            else if(strcmp(metNames.triggerName(i).c_str(), "Flag_METFilters") == 0)
                ev.flag_METFilters = metFilterBits->accept(i);
        }
    */

    /*
        edm::Handle<pat::PhotonCollection> photons;
        event.getByToken(photonTag_, photons);
        for (const pat::Photon &pho : *photons) {
            if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
            printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
                   pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
        }


    */

    summaryHandler_.fillTree();
}


//
const reco::Candidate* MainAnalyzer::findFirstMotherWithDifferentID(const reco::Candidate *particle)
{

    if( particle == 0 ) {
        printf("ERROR! null candidate pointer, this should never happen\n");
        return 0;
    }

    // Is this the first parent with a different ID? If yes, return, otherwise
    // go deeper into recursion
    if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
        if (particle->pdgId() == particle->mother(0)->pdgId()) {
            return findFirstMotherWithDifferentID(particle->mother(0));
        } else {
            return particle->mother(0);
        }
    }

    return 0;
}


void
MainAnalyzer::getMCtruth(const edm::Event& event, const edm::EventSetup& iSetup)
{
    DataEvtSummary_t &ev=summaryHandler_.getEvent();
    ev.nmcparticles = 0;
    //if(event.isRealData()) return;

    edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
    //event.getByLabel("addPileupInfo", puInfoH);
    //event.getByLabel("slimmedAddPileupInfo", puInfoH);
    event.getByToken(puInfoTag_,puInfoH);
    int npuOOT(0),npuIT(0),npuOOTm1(0);
    float truePU(0);
    if(puInfoH.isValid()) {
        for(std::vector<PileupSummaryInfo>::const_iterator it = puInfoH->begin(); it != puInfoH->end(); it++) {
            if(it->getBunchCrossing()==0) {
                npuIT += it->getPU_NumInteractions();
                truePU = it->getTrueNumInteractions();
            } else                          npuOOT += it->getPU_NumInteractions();
            if(it->getBunchCrossing()<0)  npuOOTm1 += it->getPU_NumInteractions();

        }
    }
    ev.ngenITpu=npuIT;
    ev.ngenOOTpu=npuOOT;
    ev.ngenOOTpum1=npuOOTm1;
    ev.ngenTruepu=truePU;
    controlHistos_.fillHisto("pileup","all",ev.ngenITpu);
    controlHistos_.fillHisto("pileuptrue","all",truePU);



    //retrieve pdf info
    edm::Handle<GenEventInfoProduct> genEventInfoProd;
    //event.getByLabel("generator", genEventInfoProd);
    event.getByToken(genInfoTag_, genEventInfoProd);
    ev.genWeight = genEventInfoProd->weight();
    ev.qscale = genEventInfoProd->qScale();
    if(genEventInfoProd->pdf()) {
        ev.x1  = genEventInfoProd->pdf()->x.first;
        ev.x2  = genEventInfoProd->pdf()->x.second;
        ev.id1 = genEventInfoProd->pdf()->id.first;
        ev.id2 = genEventInfoProd->pdf()->id.second;
    }
    if(genEventInfoProd->binningValues().size()>0) ev.pthat = genEventInfoProd->binningValues()[0];

    if(ev.genWeight<0) controlHistos_.fillHisto("n_negevents","all",0); //increment negative event count
    if(ev.genWeight>0) controlHistos_.fillHisto("n_posevents","all",0); //increment positive event count


    //scale variations
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW
    vector<edm::Handle<LHEEventProduct> > EvtHandles;
    event.getManyByType(EvtHandles);
    if(EvtHandles.size()>0) {
        edm::Handle<LHEEventProduct> EvtHandle = EvtHandles.front();
        if(EvtHandle.isValid() && EvtHandle->weights().size()>=9) {
            ev.weight_QCDscale_muR1_muF1 		= EvtHandle->weights()[0].wgt/EvtHandle->originalXWGTUP();
            ev.weight_QCDscale_muR1_muF2 		= EvtHandle->weights()[1].wgt/EvtHandle->originalXWGTUP();
            ev.weight_QCDscale_muR1_muF0p5 		= EvtHandle->weights()[2].wgt/EvtHandle->originalXWGTUP();
            ev.weight_QCDscale_muR2_muF1 		= EvtHandle->weights()[3].wgt/EvtHandle->originalXWGTUP();
            ev.weight_QCDscale_muR2_muF2 		= EvtHandle->weights()[4].wgt/EvtHandle->originalXWGTUP();
            ev.weight_QCDscale_muR2_muF0p5 		= EvtHandle->weights()[5].wgt/EvtHandle->originalXWGTUP();
            ev.weight_QCDscale_muR0p5_muF1 		= EvtHandle->weights()[6].wgt/EvtHandle->originalXWGTUP();
            ev.weight_QCDscale_muR0p5_muF2 		= EvtHandle->weights()[7].wgt/EvtHandle->originalXWGTUP();
            ev.weight_QCDscale_muR0p5_muF0p5 		= EvtHandle->weights()[8].wgt/EvtHandle->originalXWGTUP();
        }
    }






    //
    // gen particles
    //
    // Pruned particles are the one containing "important" stuff
    edm::Handle<edm::View<reco::GenParticle> > pruned;
    event.getByToken(prunedGenTag_,pruned);

    // Packed particles are all the status 1, so usable to remake jets
    // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    edm::Handle<edm::View<pat::PackedGenParticle> > packed;
    event.getByToken(packedGenTag_,packed);


    std::vector<TLorentzVector> chLeptons;

    if(isPythia8_) {

        std::vector<const reco::Candidate*> prunedV;
        for(size_t i=0; i<packed->size(); i++) {
            const Candidate * genParticle = &(*packed)[i];
            int pid=genParticle->pdgId();
            if(abs(pid) == 2000012 || abs(pid) == 5000039) {
                prunedV.push_back(genParticle);
            }
        }

        //Allows easier comparison for mother finding
        //std::vector<const reco::Candidate*> prunedV;
        for(size_t i=0; i<pruned->size(); i++) {
            const Candidate * genParticle = &(*pruned)[i];
            int status=genParticle->status();
            int pid=genParticle->pdgId();
            if(	//( abs(pid) >= 1  && abs(pid) <= 6 && status < 30 )
                ( abs(pid) >= 11 && abs(pid) <= 16 )
                //|| ( abs(pid) == 21 && status < 30 )
                || ( abs(pid) >= 23 && abs(pid) <= 25 && status < 30 )
                || ( abs(pid) == 1008 )
                || ( abs(pid) >= 1  && abs(pid) <= 6 && status < 30 )
                //|| ( abs(pid) >= 32 && abs(pid) <= 42 )
            ) {
//		cout << "pid: " << pid << " status: " << status << endl;
                prunedV.push_back(genParticle);
            }
        }

        //Look for mother particle and Fill gen variables
        for(unsigned int i = 0; i < prunedV.size(); i++) {
            if(prunedV[i]->numberOfMothers() > 0) {
                //find the ID of the first mother that has a different ID than the particle itself
                const reco::Candidate* mom = findFirstMotherWithDifferentID(prunedV[i]);
                if(mom) {
                    int pid    = prunedV[i]->pdgId();
                    int mompid = mom->pdgId();
                    int status = prunedV[i]->status();

                    if( !(abs(pid)>=1 && abs(pid)<=6) && abs(pid)!=23 && abs(pid)!=2000012 && abs(pid)!=5000039 && abs(mompid)!=23 && abs(mompid)!=24 && abs(mompid)!=25) continue;
                    if( status!=1 && status!=2 && status!=22  && status!=23 ) continue;

                    //cout << "saving --> pid: " << pid << " mom: " << mompid << " status: " << status << endl;

                    ev.mc_px[ev.nmcparticles] = prunedV[i]->px();
                    ev.mc_py[ev.nmcparticles] = prunedV[i]->py();
                    ev.mc_pz[ev.nmcparticles] = prunedV[i]->pz();
                    ev.mc_en[ev.nmcparticles] = prunedV[i]->energy();
                    ev.mc_id[ev.nmcparticles] = prunedV[i]->pdgId();
                    ev.mc_mom[ev.nmcparticles] = mom->pdgId();
                    ev.mc_status[ev.nmcparticles] = prunedV[i]->status();

                    TLorentzVector p4( prunedV[i]->px(), prunedV[i]->py(), prunedV[i]->pz(), prunedV[i]->energy() );
                    if(abs(pid)==11 || abs(pid)==13 || abs(pid)==15) {
                        chLeptons.push_back(p4);
                    }
                    ev.nmcparticles++;
                }
            }
        }

    } else {

        // for PYTHIA 6 based
        for(size_t i=0; i<pruned->size(); i++) {
            const Candidate * genParticle = &(*pruned)[i];
            int status=genParticle->status();
            if(status!=3) continue; // PYTHIA 6 base
            int pid=genParticle->pdgId();

            ev.mc_px[ev.nmcparticles] = genParticle->px();
            ev.mc_py[ev.nmcparticles] = genParticle->py();
            ev.mc_pz[ev.nmcparticles] = genParticle->pz();
            ev.mc_en[ev.nmcparticles] = genParticle->energy();
            ev.mc_id[ev.nmcparticles] = genParticle->pdgId();
            ev.mc_mom[ev.nmcparticles] = 0;
            ev.mc_status[ev.nmcparticles] = genParticle->status();

            TLorentzVector p4( genParticle->px(), genParticle->py(), genParticle->pz(), genParticle->energy() );
            if(abs(pid)==11 || abs(pid)==13 || abs(pid)==15) {
                chLeptons.push_back(p4);

            }

            // check mother particle id just for leptons
            if((abs(pid)>10 && abs(pid)<17) || abs(pid)==23 || abs(pid)==24 || abs(pid)==25) {
                size_t nGenMothers = genParticle->numberOfMothers();

                // Madgraph, Pythia6, POWHEG, not for Sherpa
                if(nGenMothers==1) {
                    const reco::GenParticle* genMother = dynamic_cast<const reco::GenParticle*>(genParticle->mother(0));
                    ev.mc_mom[ev.nmcparticles] = genMother->pdgId();
                }
            }
            ev.nmcparticles++;
        }
    }

//
// gen jets
//
    edm::Handle<edm::View<reco::GenJet> > genJets;
    event.getByToken(genjetTag_, genJets);
    if(!genJets.isValid())     cerr << "  WARNING: genJets is not valid! " << endl;
    std::vector<TLorentzVector> jets;
    for(size_t j=0; j<genJets->size(); j++) {
        const Candidate * genJet = &(*genJets)[j];
        TLorentzVector p4( genJet->px(), genJet->py(), genJet->pz(), genJet->energy() );
        if(p4.Pt()<10 || fabs(p4.Eta())>2.5) continue;

        bool matchesLepton(false);
        for(size_t i=0; i<chLeptons.size(); i++) {
            float dR=p4.DeltaR(chLeptons[i]);
            if(dR>0.4) continue;
            matchesLepton=true;
            break;
        }
        if(matchesLepton) continue;

        jets.push_back(p4);
        ev.mc_px[ev.nmcparticles]=genJet->px();
        ev.mc_py[ev.nmcparticles]=genJet->py();
        ev.mc_pz[ev.nmcparticles]=genJet->pz();
        ev.mc_en[ev.nmcparticles]=genJet->energy();
        ev.mc_status[ev.nmcparticles]=0; // special for genjet
        ev.mc_id[ev.nmcparticles]=1; // special for genjet
        ev.mc_mom[ev.nmcparticles] = 0;
        ev.nmcparticles++;
    }


}


//
bool
MainAnalyzer::checkIfTriggerFired(edm::Handle<edm::TriggerResults> &allTriggerBits, const edm::TriggerNames &triggerNames, std::string triggerPath)
{
    for (size_t itrig = 0; itrig != allTriggerBits->size(); ++itrig) {
        std::string trigName = triggerNames.triggerName(itrig);
        if( !allTriggerBits->wasrun(itrig) ) continue;
        if( allTriggerBits->error(itrig)   ) continue;
        if( !allTriggerBits->accept(itrig) ) continue;
        if(trigName.find(triggerPath)!= std::string::npos) return true;
    }
    return false;
}


//
// The function that uses algorith from Josh Bendavid with
// an explicit loop over gen particles.
int MainAnalyzer::matchToTruth(const pat::Electron &el,
                               const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles)
{

    //
    // Explicit loop and geometric matching method (advised by Josh Bendavid)
    //

    // Find the closest status 1 gen electron to the reco electron
    double dR = 999;
    const reco::Candidate *closestElectron = 0;
    for(size_t i=0; i<prunedGenParticles->size(); i++) {
        const reco::Candidate *particle = &(*prunedGenParticles)[i];
        // Drop everything that is not electron or not status 1
        if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
            continue;
        //
        double dRtmp = ROOT::Math::VectorUtil::DeltaR( el.p4(), particle->p4() );
        if( dRtmp < dR ) {
            dR = dRtmp;
            closestElectron = particle;
        }
    }
    // See if the closest electron (if it exists) is close enough.
    // If not, no match found.
    if( !(closestElectron != 0 && dR < 0.1) ) {
        return UNMATCHED;
    }

    //
    int ancestorPID = -999;
    int ancestorStatus = -999;
    findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

    if( ancestorPID == -999 && ancestorStatus == -999 ) {
        // No non-electron parent??? This should never happen.
        // Complain.
        printf("MainAnalyzer: ERROR! Electron does not apper to have a non-electron parent\n");
        return UNMATCHED;
    }

    if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
        return TRUE_NON_PROMPT_ELECTRON;

    if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
        return TRUE_ELECTRON_FROM_TAU;

    // What remains is true prompt electrons
    return TRUE_PROMPT_ELECTRON;
}


void MainAnalyzer::findFirstNonElectronMother(const reco::Candidate *particle,
        int &ancestorPID, int &ancestorStatus)
{

    if( particle == 0 ) {
        printf("MainAnalyzer: ERROR! null candidate pointer, this should never happen\n");
        return;
    }

    // Is this the first non-electron parent? If yes, return, otherwise
    // go deeper into recursion
    if( abs(particle->pdgId()) == 11 ) {
        findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
    } else {
        ancestorPID = particle->pdgId();
        ancestorStatus = particle->status();
    }

    return;
}








// ------------ method called once each job just before starting event loop  ------------
void
MainAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
MainAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
MainAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
MainAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MainAnalyzer::beginLuminosityBlock(edm::LuminosityBlock &iLumi, edm::EventSetup &iSetup)
{
  edm::Handle<LumiSummary> l;
  iLumi.getByLabel("lumiProducer", l);
  if (!l.isValid())  return;
  curAvgInstLumi_ = l->avgInsDelLumi();
  curIntegLumi_   = l->intgDelLumi();
  controlHistos_.fillHisto("instlumi","all",curAvgInstLumi_);
  controlHistos_.fillHisto("integlumi","all",curIntegLumi_);
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MainAnalyzer::endLuminosityBlock(edm::LuminosityBlock &iLumi, edm::EventSetup &iSetup)
{
    TString filterCtrs[]= {"startCounter","noScrapCounter","goodVertexCounter","noHBHEnoiseCounter","nobeamHaloCounter"};
    const size_t nselFilters=sizeof(filterCtrs)/sizeof(TString);
    for(size_t istep=0; istep<nselFilters; istep++) {
        std::string fname(filterCtrs[istep].Data());
        try {
            edm::Handle<edm::MergeableCounter> ctrHandle;
            iLumi.getByLabel(fname, ctrHandle);
            if(!ctrHandle.isValid()) continue;
            controlHistos_.fillHisto("cutflow","all",istep,ctrHandle->value);
        } catch(std::exception) {
            controlHistos_.fillHisto("cutflow","all",istep);
        }
    }
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MainAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MainAnalyzer);
