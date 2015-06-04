#ifndef dataevtsummaryhandler_h
#define dataevtsummaryhandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>

#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include "TTree.h"
#include "TLorentzVector.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;

#define MAXPARTICLES 50
#define MAXMCPARTICLES 250

struct DataEvtSummary_t {

    Int_t run,lumi,event;
    Float_t curAvgInstLumi,curIntegLumi;
    Bool_t hasTrigger;
    Int_t triggerType;

    //primary vertex
    Int_t nvtx;
    Float_t vtx_x, vtx_y, vtx_z;
    Float_t fixedGridRhoAll,fixedGridRhoFastjetAll,fixedGridRhoFastjetAllCalo,fixedGridRhoFastjetCentralCalo,fixedGridRhoFastjetCentralChargedPileUp,fixedGridRhoFastjetCentralNeutral;

    //generator info
    Int_t ngenITpu,ngenOOTpu,ngenOOTpum1, ngenTruepu;
    Float_t puWeight, hptWeights[5];
    Float_t pthat,genWeight, qscale, x1,x2;
    Int_t id1,id2;

    //gen level event
    Int_t nmcparticles;
    Float_t mc_px[MAXMCPARTICLES],mc_py[MAXMCPARTICLES],mc_pz[MAXMCPARTICLES],mc_en[MAXMCPARTICLES];
    Int_t mc_id[MAXMCPARTICLES], mc_status[MAXMCPARTICLES], mc_mom[MAXMCPARTICLES];

    //muon
    Int_t mn;
    Float_t mn_px[MAXPARTICLES],mn_py[MAXPARTICLES],mn_pz[MAXPARTICLES],mn_en[MAXPARTICLES];
    Int_t mn_id[MAXPARTICLES], mn_type[MAXPARTICLES];
    Float_t mn_d0[MAXPARTICLES],mn_dZ[MAXPARTICLES],mn_ip3d[MAXPARTICLES],mn_ip3dsig[MAXPARTICLES];
    Bool_t mn_IsLoose[MAXPARTICLES],mn_IsTight[MAXPARTICLES],mn_IsSoft[MAXPARTICLES],mn_IsHighPt[MAXPARTICLES];
    Float_t mn_pileupIsoR03[MAXPARTICLES],mn_chargedIsoR03[MAXPARTICLES],mn_photonIsoR03[MAXPARTICLES],mn_neutralHadIsoR03[MAXPARTICLES];
    Float_t mn_pileupIsoR04[MAXPARTICLES],mn_chargedIsoR04[MAXPARTICLES],mn_photonIsoR04[MAXPARTICLES],mn_neutralHadIsoR04[MAXPARTICLES];

    //Float_t mn_nMatches[MAXPARTICLES],mn_nMatchedStations[MAXPARTICLES],mn_validMuonHits[MAXPARTICLES],mn_innerTrackChi2[MAXPARTICLES],mn_trkLayersWithMeasurement[MAXPARTICLES],mn_pixelLayersWithMeasurement[MAXPARTICLES];

    //electron
    Int_t en;
    Float_t en_px[MAXPARTICLES],en_py[MAXPARTICLES],en_pz[MAXPARTICLES],en_en[MAXPARTICLES];
    Int_t en_id[MAXPARTICLES];
    Float_t en_d0[MAXPARTICLES],en_dZ[MAXPARTICLES];
    Float_t en_EtaSC[MAXPARTICLES],en_PhiSC[MAXPARTICLES],en_EnSC[MAXPARTICLES];
    Float_t en_dEtaIn[MAXPARTICLES],en_dPhiIn[MAXPARTICLES],en_hOverE[MAXPARTICLES],en_R9[MAXPARTICLES],en_sigmaIetaIeta[MAXPARTICLES],en_sigmaIetaIeta5x5[MAXPARTICLES],en_ooEmooP[MAXPARTICLES];
    Float_t en_pileupIso[MAXPARTICLES],en_chargedIso[MAXPARTICLES],en_photonIso[MAXPARTICLES],en_neutralHadIso[MAXPARTICLES];
    Float_t en_relIsoWithEA[MAXPARTICLES],en_relIsoWithDBeta[MAXPARTICLES],en_MissingHits[MAXPARTICLES],en_passConversionVeto[MAXPARTICLES];
    Bool_t en_passVeto[MAXPARTICLES],en_passLoose[MAXPARTICLES],en_passMedium[MAXPARTICLES],en_passTight[MAXPARTICLES];
    Float_t en_IDMVATrig[MAXPARTICLES],en_IDMVANonTrig[MAXPARTICLES];
    Int_t en_istrue[MAXPARTICLES];


    //tau
    Int_t ta;
    Float_t ta_px[MAXPARTICLES],ta_py[MAXPARTICLES],ta_pz[MAXPARTICLES],ta_en[MAXPARTICLES];
    Int_t ta_id[MAXPARTICLES];
    Bool_t ta_dm[MAXPARTICLES],ta_newdm[MAXPARTICLES];
    Bool_t ta_IsLooseIso[MAXPARTICLES],ta_IsMediumIso[MAXPARTICLES],ta_IsTightIso[MAXPARTICLES];
    Float_t ta_combIsoDBeta3Hits[MAXPARTICLES],ta_chargedIso[MAXPARTICLES],ta_neutralIso[MAXPARTICLES],ta_pileupIso[MAXPARTICLES];
    Bool_t ta_passEleVetoLoose[MAXPARTICLES],ta_passEleVetoMedium[MAXPARTICLES],ta_passEleVetoTight[MAXPARTICLES],ta_passMuVetoLoose3[MAXPARTICLES],ta_passMuVetoTight3[MAXPARTICLES];

    //jet (ak4PFJetsCHS)
    Int_t jet;
    Float_t jet_px[MAXPARTICLES],jet_py[MAXPARTICLES],jet_pz[MAXPARTICLES],jet_en[MAXPARTICLES];
    Float_t jet_btag0[MAXPARTICLES],jet_btag1[MAXPARTICLES],jet_btag2[MAXPARTICLES],jet_btag3[MAXPARTICLES];
    Float_t jet_btag4[MAXPARTICLES],jet_btag5[MAXPARTICLES],jet_btag6[MAXPARTICLES],jet_btag7[MAXPARTICLES];
    Float_t jet_mass[MAXPARTICLES],jet_area[MAXPARTICLES],jet_pu[MAXPARTICLES],jet_puId[MAXPARTICLES];
    Bool_t jet_PFLoose[MAXPARTICLES], jet_PFTight[MAXPARTICLES];
    Int_t jet_partonFlavour[MAXPARTICLES], jet_hadronFlavour[MAXPARTICLES];


    //fjet (ak8PFJetsCHS)
    Int_t fjet;
    Float_t fjet_px[MAXPARTICLES],fjet_py[MAXPARTICLES],fjet_pz[MAXPARTICLES],fjet_en[MAXPARTICLES];
    Float_t fjet_prunedM[MAXPARTICLES],fjet_trimmedM[MAXPARTICLES],fjet_filteredM[MAXPARTICLES];
    Float_t fjet_tau1[MAXPARTICLES],fjet_tau2[MAXPARTICLES],fjet_tau3[MAXPARTICLES];


    //met
    Float_t met_pt,met_phi,met_sumMET;
    Bool_t flag_HBHENoiseFilter,flag_CSCTightHaloFilter,flag_hcalLaserEventFilter,flag_EcalDeadCellTriggerPrimitiveFilter,flag_goodVertices;
    Bool_t flag_trackingFailureFilter,flag_eeBadScFilter,flag_ecalLaserCorrFilter,flag_trkPOGFilters,flag_trkPOG_manystripclus53X,flag_trkPOG_toomanystripclus53X;
    Bool_t flag_trkPOG_logErrorTooManyClusters,flag_METFilters;

};

class DataEvtSummaryHandler {
public:
    //
    DataEvtSummaryHandler();
    ~DataEvtSummaryHandler();

    //current event
    DataEvtSummary_t evSummary_;
    DataEvtSummary_t &getEvent() {
        return evSummary_;
    }

    //write mode
    bool initTree(TTree *t);
    void fillTree();

    //read mode
    bool attachToTree(TTree *t);
    int getEntries() { return (t_ ? t_->GetEntriesFast() : 0); }
    void getEntry(int ientry) {
    	resetStruct();
    	if(t_) t_->GetEntry(ientry);
    }

    void resetStruct();

private:
    //the tree
    TTree *t_;
};

#endif
