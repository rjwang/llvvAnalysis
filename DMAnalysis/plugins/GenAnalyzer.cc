// -*- C++ -*-
//
// Package:    GenAnalyzer
// Class:      GenAnalyzer
//
/**\class GenAnalyzer GenAnalyzer.cc llvvAnalysis/DMAnalysis/plugins/GenAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Renjie Wang,32 4-C21,,
//         Created:  Mon Jan 19 18:18:24 CET 2015
// $Id$
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
// class declaration
//

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/METReco/interface/GenMET.h"


#include "llvvAnalysis/DMAnalysis/interface/GenEvtSummaryHandler.h"
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>


using namespace std;
using namespace edm;
using namespace reco;

#define MAXMCPARTICLES 250

class GenAnalyzer : public edm::EDAnalyzer {
public:
    GenAnalyzer(const edm::ParameterSet&);
    ~GenAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    edm::InputTag _GenJets_;
    edm::InputTag _GenParticles_;
    bool isPythia8_;
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

    const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);

    // --------- tree ----------------------------
    TTree *tree_;
    GenEvtSummaryHandler summaryHandler_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
{
    //now do what ever initialization is needed
    _GenJets_ =  iConfig.getUntrackedParameter<edm::InputTag>("GenJets");
    _GenParticles_ = iConfig.getUntrackedParameter<edm::InputTag>("GenParticles");
    isPythia8_ = iConfig.getParameter<bool>("isPythia8");

}

GenAnalyzer::~GenAnalyzer()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


    //event summary to be filled
    GenEvtSummary_t &ev = summaryHandler_.getEvent();
    ev.nmcparticles = 0;

    //event header
    ev.run    = iEvent.id().run();
    ev.lumi   = iEvent.luminosityBlock();
    ev.event  = iEvent.id().event();

    //
    // gen particles
    //
    edm::Handle< std::vector<reco::GenParticle> > genParticles;
    iEvent.getByLabel(_GenParticles_, genParticles);
    if(!genParticles.isValid())     cerr << "  WARNING: genParticles is not valid! " << endl;
    std::vector<TLorentzVector> chLeptons;


    if(isPythia8_) {
        //Allows easier comparison for mother finding
        std::vector<const reco::Candidate*> prunedV;
        for(size_t i=0; i<genParticles->size(); i++) {
            const Candidate * genParticle = &(*genParticles)[i];
            int status=genParticle->status();
            int pid=genParticle->pdgId();
            if(	//( abs(pid) >= 1  && abs(pid) <= 6 && status < 30 )
                ( abs(pid) >= 11 && abs(pid) <= 16 )
                //|| ( abs(pid) == 21 && status < 30 )
                || ( abs(pid) >= 23 && abs(pid) <= 25 && status < 30 )
                || ( abs(pid) == 2000012 ) //Dark Matter
                || ( abs(pid) == 5000039 ) //Graviton/Unparticle
            ) {
		//cout << "status: " << status << " PID: " << pid << endl;
                prunedV.push_back(genParticle);
            }
        }

        //Look for mother particle and Fill gen variables
        for(unsigned int i = 0; i < prunedV.size(); i++) {
            if(prunedV[i]->numberOfMothers() > 0) {
                //find the ID of the first mother that has a different ID than the particle itself
                const reco::Candidate* mom = findFirstMotherWithDifferentID(prunedV[i]);
                if(mom) {
                    int pid=prunedV[i]->pdgId();
                    int mompid = mom->pdgId();
                    if(abs(pid) != 2000012 && abs(pid) != 5000039 && abs(mompid)!=23 && abs(mompid)!=24 && abs(mompid)!=25) continue;
                    if(prunedV[i]->status()!=1 && prunedV[i]->status()!=2) continue;


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

        for(size_t i=0; i<genParticles->size(); i++) {
            const Candidate * genParticle = &(*genParticles)[i];
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
    edm::Handle< std::vector<reco::GenJet> > genJets;
    //iEvent.getByLabel("ak5GenJets", genJets);
    //iEvent.getByLabel("ak5GenJetsNoMuNoNu", genJets);
    iEvent.getByLabel(_GenJets_, genJets);
    if(!genJets.isValid())     cerr << "  WARNING: genJets is not valid! " << endl;
    std::vector<TLorentzVector> jets;
    for(std::vector<reco::GenJet>::const_iterator genJet=genJets->begin(); genJet!=genJets->end(); genJet++) {
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
        ev.mc_status[ev.nmcparticles]=0; // special for jet
        ev.mc_id[ev.nmcparticles]=1; // special for jet
        ev.nmcparticles++;
    }

    tree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
GenAnalyzer::beginJob()
{
    //book the histograms
    edm::Service<TFileService> fs;
    tree_              = fs->make<TTree>("data","Event Summary");
    summaryHandler_.initTree(tree_);

}

const reco::Candidate* GenAnalyzer::findFirstMotherWithDifferentID(const reco::Candidate *particle)
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


// ------------ method called once each job just after ending the event loop  ------------
void
GenAnalyzer::endJob()
{
}

// ------------ method called when starting to processes a run  ------------
void
GenAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
GenAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
GenAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
GenAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
