#include "llvvAnalysis/DMAnalysis/interface/DataEvtSummaryHandler.h"

using namespace std;

//
DataEvtSummaryHandler::DataEvtSummaryHandler()
{
}

//
bool DataEvtSummaryHandler::initTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;

    //event info
    t_->Branch("run",        	&evSummary_.run,            "run/I");
    t_->Branch("lumi",       	&evSummary_.lumi,           "lumi/I");
    t_->Branch("event",      	&evSummary_.event,          "event/I");
    t_->Branch("curAvgInstLumi",&evSummary_.curAvgInstLumi, "curAvgInstLumi/F");
    t_->Branch("curIntegLumi", 	&evSummary_.curIntegLumi,   "curIntegLumi/F");

    t_->Branch("hasTrigger",   	&evSummary_.hasTrigger,     "hasTrigger/O");
    t_->Branch("triggerType",	&evSummary_.triggerType,    "triggerType/I");

    //primary vertex
    t_->Branch("nvtx",      	&evSummary_.nvtx,           "nvtx/I");
    t_->Branch("vtx_x",         &evSummary_.vtx_x,          "vtx_x/F");
    t_->Branch("vtx_y",         &evSummary_.vtx_y,          "vtx_y/F");
    t_->Branch("vtx_z",         &evSummary_.vtx_z,          "vtx_z/F");

    t_->Branch("fixedGridRhoAll",                           &evSummary_.fixedGridRhoAll,                        "fixedGridRhoAll/F");
    t_->Branch("fixedGridRhoFastjetAll",                    &evSummary_.fixedGridRhoFastjetAll,                 "fixedGridRhoFastjetAll/F");
    t_->Branch("fixedGridRhoFastjetAllCalo",                &evSummary_.fixedGridRhoFastjetAllCalo,             "fixedGridRhoFastjetAllCalo/F");
    t_->Branch("fixedGridRhoFastjetCentralCalo",            &evSummary_.fixedGridRhoFastjetCentralCalo,         "fixedGridRhoFastjetCentralCalo/F");
    t_->Branch("fixedGridRhoFastjetCentralChargedPileUp",   &evSummary_.fixedGridRhoFastjetCentralChargedPileUp,"fixedGridRhoFastjetCentralChargedPileUp/F");
    t_->Branch("fixedGridRhoFastjetCentralNeutral",         &evSummary_.fixedGridRhoFastjetCentralNeutral,      "fixedGridRhoFastjetCentralNeutral/F");

    //generator level info
    t_->Branch("ngenITpu",   	&evSummary_.ngenITpu,       "ngenITpu/I");
    t_->Branch("ngenOOTpu",  	&evSummary_.ngenOOTpu,      "ngenOOTpu/I");
    t_->Branch("ngenOOTpum1",  	&evSummary_.ngenOOTpum1,    "ngenOOTpum1/I");
    t_->Branch("ngenTruepu",  	&evSummary_.ngenTruepu,     "ngenTruepu/I");
    t_->Branch("puweight",     	&evSummary_.puWeight,       "puweight/F");
    t_->Branch("hptWeights", 	evSummary_.hptWeights,      "hptWeights[5]/F");
    t_->Branch("pthat",      	&evSummary_.pthat,          "pthat/F");
    t_->Branch("genWeight",  	&evSummary_.genWeight,      "genWeight/F");
    t_->Branch("qscale",     	&evSummary_.qscale,         "qscale/F");
    t_->Branch("x1",         	&evSummary_.x1,             "x1/F");
    t_->Branch("x2",         	&evSummary_.x2,             "x2/F");
    t_->Branch("id1",        	&evSummary_.id1,            "id1/I");
    t_->Branch("id2",        	&evSummary_.id2,            "id2/I");


    //mc truth
    t_->Branch("nmcparticles", 	&evSummary_.nmcparticles,   "nmcparticles/I");
    t_->Branch("mc_px",         evSummary_.mc_px,           "mc_px[nmcparticles]/F");
    t_->Branch("mc_py",         evSummary_.mc_py,           "mc_py[nmcparticles]/F");
    t_->Branch("mc_pz",         evSummary_.mc_pz,           "mc_pz[nmcparticles]/F");
    t_->Branch("mc_en",         evSummary_.mc_en,           "mc_en[nmcparticles]/F");
    t_->Branch("mc_id",         evSummary_.mc_id,           "mc_id[nmcparticles]/I");
    t_->Branch("mc_status", 	evSummary_.mc_status,       "mc_status[nmcparticles]/I");
    t_->Branch("mc_mom",        evSummary_.mc_mom,          "mc_mom[nmcparticles]/I");


    //muon
    t_->Branch("mn",                    &evSummary_.mn,                     "mn/I");
    t_->Branch("mn_px",                 evSummary_.mn_px,                   "mn_px[mn]/F");
    t_->Branch("mn_py",                 evSummary_.mn_py,                   "mn_py[mn]/F");
    t_->Branch("mn_pz",                 evSummary_.mn_pz,                   "mn_pz[mn]/F");
    t_->Branch("mn_en",                 evSummary_.mn_en,                   "mn_en[mn]/F");
    t_->Branch("mn_id",                 evSummary_.mn_id,                   "mn_id[mn]/I");
    t_->Branch("mn_type",               evSummary_.mn_type,                 "mn_type[mn]/I");
    t_->Branch("mn_d0",                 evSummary_.mn_d0,                   "mn_d0[mn]/F");
    t_->Branch("mn_dZ",                 evSummary_.mn_dZ,                   "mn_dZ[mn]/F");
    t_->Branch("mn_ip3d",               evSummary_.mn_ip3d,                 "mn_ip3d[mn]/F");
    t_->Branch("mn_ip3dsig",            evSummary_.mn_ip3dsig,              "mn_ip3dsig[mn]/F");
    t_->Branch("mn_IsLoose",            evSummary_.mn_IsLoose,              "mn_IsLoose[mn]/O");
    t_->Branch("mn_IsTight",            evSummary_.mn_IsTight,              "mn_IsTight[mn]/O");
    t_->Branch("mn_IsSoft",            	evSummary_.mn_IsSoft,               "mn_IsSoft[mn]/O");
    t_->Branch("mn_IsHighPt",           evSummary_.mn_IsHighPt,             "mn_IsHighPt[mn]/O");
    t_->Branch("mn_pileupIsoR03",       evSummary_.mn_pileupIsoR03,         "mn_pileupIsoR03[mn]/F");
    t_->Branch("mn_chargedIsoR03",      evSummary_.mn_chargedIsoR03,        "mn_chargedIsoR03[mn]/F");
    t_->Branch("mn_photonIsoR03",       evSummary_.mn_photonIsoR03,         "mn_photonIsoR03[mn]/F");
    t_->Branch("mn_neutralHadIsoR03",   evSummary_.mn_neutralHadIsoR03,     "mn_neutralHadIsoR03[mn]/F");
    t_->Branch("mn_pileupIsoR04",       evSummary_.mn_pileupIsoR04,         "mn_pileupIsoR04[mn]/F");
    t_->Branch("mn_chargedIsoR04",      evSummary_.mn_chargedIsoR04,        "mn_chargedIsoR04[mn]/F");
    t_->Branch("mn_photonIsoR04",       evSummary_.mn_photonIsoR04,         "mn_photonIsoR04[mn]/F");
    t_->Branch("mn_neutralHadIsoR04",   evSummary_.mn_neutralHadIsoR04,     "mn_neutralHadIsoR04[mn]/F");

    /*
        t_->Branch("mn_nMatches",                   evSummary_.mn_nMatches,                     "mn_nMatches[mn]/F");
        t_->Branch("mn_nMatchedStations",           evSummary_.mn_nMatchedStations,             "mn_nMatchedStations[mn]/F");
        t_->Branch("mn_validMuonHits",              evSummary_.mn_validMuonHits,                "mn_validMuonHits[mn]/F");
        t_->Branch("mn_innerTrackChi2",             evSummary_.mn_innerTrackChi2,               "mn_innerTrackChi2[mn]/F");
        t_->Branch("mn_trkLayersWithMeasurement",   evSummary_.mn_trkLayersWithMeasurement,     "mn_trkLayersWithMeasurement[mn]/F");
        t_->Branch("mn_pixelLayersWithMeasurement", evSummary_.mn_pixelLayersWithMeasurement,   "mn_pixelLayersWithMeasurement[mn]/F");
    */

    //electron
    t_->Branch("en",                    &evSummary_.en,                     "en/I");
    t_->Branch("en_px",                 evSummary_.en_px,                   "en_px[en]/F");
    t_->Branch("en_py",                 evSummary_.en_py,                   "en_py[en]/F");
    t_->Branch("en_pz",                 evSummary_.en_pz,                   "en_pz[en]/F");
    t_->Branch("en_en",                 evSummary_.en_en,                   "en_en[en]/F");
    t_->Branch("en_id",                 evSummary_.en_id,                   "en_id[en]/I");
    t_->Branch("en_d0",                 evSummary_.en_d0,                   "en_d0[en]/F");
    t_->Branch("en_dZ",                 evSummary_.en_dZ,                   "en_dZ[en]/F");
    t_->Branch("en_EtaSC",              evSummary_.en_EtaSC,                "en_EtaSC[en]/F");
    t_->Branch("en_PhiSC",              evSummary_.en_PhiSC,                "en_PhiSC[en]/F");
    t_->Branch("en_EnSC",               evSummary_.en_EnSC,                 "en_EnSC[en]/F");
    t_->Branch("en_dEtaIn",             evSummary_.en_dEtaIn,               "en_dEtaIn[en]/F");
    t_->Branch("en_dPhiIn",             evSummary_.en_dPhiIn,               "en_dPhiIn[en]/F");
    t_->Branch("en_hOverE",             evSummary_.en_hOverE,               "en_hOverE[en]/F");
    t_->Branch("en_R9",                 evSummary_.en_R9,                   "en_R9[en]/F");
    t_->Branch("en_sigmaIetaIeta",      evSummary_.en_sigmaIetaIeta,        "en_sigmaIetaIeta[en]/F");
    t_->Branch("en_sigmaIetaIeta5x5",   evSummary_.en_sigmaIetaIeta5x5,     "en_sigmaIetaIeta5x5[en]/F");
    t_->Branch("en_ooEmooP",          	evSummary_.en_ooEmooP,              "en_ooEmooP[en]/F");
    t_->Branch("en_pileupIso",          evSummary_.en_pileupIso,            "en_pileupIso[en]/F");
    t_->Branch("en_chargedIso",         evSummary_.en_chargedIso,           "en_chargedIso[en]/F");
    t_->Branch("en_photonIso",          evSummary_.en_photonIso,            "en_photonIso[en]/F");
    t_->Branch("en_neutralHadIso",      evSummary_.en_neutralHadIso,        "en_neutralHadIso[en]/F");
    t_->Branch("en_relIsoWithEA",       evSummary_.en_relIsoWithEA, 	    "en_relIsoWithEA[en]/F");
    t_->Branch("en_relIsoWithDBeta",    evSummary_.en_relIsoWithDBeta,	    "en_relIsoWithDBeta[en]/F");
    t_->Branch("en_MissingHits",        evSummary_.en_MissingHits,          "en_MissingHits[en]/F");
    t_->Branch("en_passConversionVeto", evSummary_.en_passConversionVeto,   "en_passConversionVeto[en]/F");
    t_->Branch("en_passVeto",           evSummary_.en_passVeto,             "en_passVeto[en]/O");
    t_->Branch("en_passLoose",          evSummary_.en_passLoose,            "en_passLoose[en]/O");
    t_->Branch("en_passMedium",         evSummary_.en_passMedium,           "en_passMedium[en]/O");
    t_->Branch("en_passTight",          evSummary_.en_passTight,            "en_passTight[en]/O");
    t_->Branch("en_passHEEP",           evSummary_.en_passHEEP,             "en_passHEEP[en]/O");
    t_->Branch("en_passMVATrigMedium",  evSummary_.en_passMVATrigMedium,             "en_passMVATrigMedium[en]/O");
    t_->Branch("en_passMVATrigTight",   evSummary_.en_passMVATrigTight,             "en_passMVATrigTight[en]/O");
    t_->Branch("en_IDMVATrigValue",     evSummary_.en_IDMVATrigValue,            "en_IDMVATrigValue[en]/F");
    t_->Branch("en_IDMVATrigCategory",  evSummary_.en_IDMVATrigCategory,         "en_IDMVATrigCategory[en]/I");
    t_->Branch("en_istrue",             evSummary_.en_istrue,               "en_istrue[en]/I");



    //tau
    t_->Branch("ta",                    &evSummary_.ta,                     "ta/I");
    t_->Branch("ta_px",                 evSummary_.ta_px,                   "ta_px[ta]/F");
    t_->Branch("ta_py",                 evSummary_.ta_py,                   "ta_py[ta]/F");
    t_->Branch("ta_pz",                 evSummary_.ta_pz,                   "ta_pz[ta]/F");
    t_->Branch("ta_en",                 evSummary_.ta_en,                   "ta_en[ta]/F");
    t_->Branch("ta_id",                 evSummary_.ta_id,                   "ta_id[ta]/I");
    t_->Branch("ta_dm",                 evSummary_.ta_dm,                   "ta_dm[ta]/O");
    t_->Branch("ta_newdm",              evSummary_.ta_newdm,                "ta_newdm[ta]/O");
    t_->Branch("ta_IsLooseIso",         evSummary_.ta_IsLooseIso,           "ta_IsLooseIso[ta]/O");
    t_->Branch("ta_IsMediumIso",        evSummary_.ta_IsMediumIso,          "ta_IsMediumIso[ta]/O");
    t_->Branch("ta_IsTightIso",         evSummary_.ta_IsTightIso,           "ta_IsTightIso[ta]/O");
    t_->Branch("ta_combIsoDBeta3Hits",  evSummary_.ta_combIsoDBeta3Hits,    "ta_combIsoDBeta3Hits[ta]/F");
    t_->Branch("ta_chargedIso",         evSummary_.ta_chargedIso,           "ta_chargedIso[ta]/F");
    t_->Branch("ta_neutralIso",         evSummary_.ta_neutralIso,           "ta_neutralIso[ta]/F");
    t_->Branch("ta_pileupIso",          evSummary_.ta_pileupIso,            "ta_pileupIso[ta]/F");
    t_->Branch("ta_passEleVetoLoose",   evSummary_.ta_passEleVetoLoose,     "ta_passEleVetoLoose[ta]/O");
    t_->Branch("ta_passEleVetoMedium",  evSummary_.ta_passEleVetoMedium,    "ta_passEleVetoMedium[ta]/O");
    t_->Branch("ta_passEleVetoTight",   evSummary_.ta_passEleVetoTight,     "ta_passEleVetoTight[ta]/O");
    t_->Branch("ta_passMuVetoLoose3",   evSummary_.ta_passMuVetoLoose3,     "ta_passMuVetoLoose3[ta]/O");
    t_->Branch("ta_passMuVetoTight3",   evSummary_.ta_passMuVetoTight3,     "ta_passMuVetoTight3[ta]/O");



    //jet (ak4PFJetsCHS)
    t_->Branch("jet",                   &evSummary_.jet,                    "jet/I");
    t_->Branch("jet_px",                evSummary_.jet_px,                  "jet_px[jet]/F");
    t_->Branch("jet_py",                evSummary_.jet_py,                  "jet_py[jet]/F");
    t_->Branch("jet_pz",                evSummary_.jet_pz,                  "jet_pz[jet]/F");
    t_->Branch("jet_en",                evSummary_.jet_en,                  "jet_en[jet]/F");
    t_->Branch("jet_btag0",             evSummary_.jet_btag0,               "jet_btag0[jet]/F");
    t_->Branch("jet_btag1",             evSummary_.jet_btag1,               "jet_btag1[jet]/F");
    t_->Branch("jet_btag2",             evSummary_.jet_btag2,               "jet_btag2[jet]/F");
    t_->Branch("jet_btag3",             evSummary_.jet_btag3,               "jet_btag3[jet]/F");
    t_->Branch("jet_btag4",             evSummary_.jet_btag4,               "jet_btag4[jet]/F");
    t_->Branch("jet_btag5",             evSummary_.jet_btag5,               "jet_btag5[jet]/F");
    t_->Branch("jet_btag6",             evSummary_.jet_btag6,               "jet_btag6[jet]/F");
    t_->Branch("jet_btag7",             evSummary_.jet_btag7,               "jet_btag7[jet]/F");
    t_->Branch("jet_btag8",             evSummary_.jet_btag8,               "jet_btag8[jet]/F");
    t_->Branch("jet_btag9",             evSummary_.jet_btag9,               "jet_btag9[jet]/F");
    t_->Branch("jet_btag10",            evSummary_.jet_btag10,              "jet_btag10[jet]/F");
    t_->Branch("jet_mass",              evSummary_.jet_mass,                "jet_mass[jet]/F");
    t_->Branch("jet_area",              evSummary_.jet_area,                "jet_area[jet]/F");
    t_->Branch("jet_pu",                evSummary_.jet_pu,                  "jet_pu[jet]/F");
    t_->Branch("jet_puId",              evSummary_.jet_puId,                "jet_puId[jet]/F");
    t_->Branch("jet_PFLoose",           evSummary_.jet_PFLoose,             "jet_PFLoose[jet]/O");
    t_->Branch("jet_PFTight",           evSummary_.jet_PFTight,             "jet_PFTight[jet]/O");
    t_->Branch("jet_partonFlavour",     evSummary_.jet_partonFlavour,       "jet_partonFlavour[jet]/I");
    t_->Branch("jet_hadronFlavour",     evSummary_.jet_hadronFlavour,       "jet_hadronFlavour[jet]/I");
    t_->Branch("jet_genpt",             evSummary_.jet_genpt,               "jet_genpt[jet]/F");

    //jet (slimmedJetsPuppi)
    t_->Branch("pjet",                   &evSummary_.pjet,                    "pjet/I");
    t_->Branch("pjet_px",                evSummary_.pjet_px,                  "pjet_px[pjet]/F");
    t_->Branch("pjet_py",                evSummary_.pjet_py,                  "pjet_py[pjet]/F");
    t_->Branch("pjet_pz",                evSummary_.pjet_pz,                  "pjet_pz[pjet]/F");
    t_->Branch("pjet_en",                evSummary_.pjet_en,                  "pjet_en[pjet]/F");
    t_->Branch("pjet_genpt",             evSummary_.pjet_genpt,               "pjet_genpt[pjet]/F");
    t_->Branch("pjet_btag0",             evSummary_.pjet_btag0,               "pjet_btag0[pjet]/F");
    t_->Branch("pjet_btag1",             evSummary_.pjet_btag1,               "pjet_btag1[pjet]/F");
    t_->Branch("pjet_btag2",             evSummary_.pjet_btag2,               "pjet_btag2[pjet]/F");
    t_->Branch("pjet_btag3",             evSummary_.pjet_btag3,               "pjet_btag3[pjet]/F");
    t_->Branch("pjet_btag4",             evSummary_.pjet_btag4,               "pjet_btag4[pjet]/F");
    t_->Branch("pjet_btag5",             evSummary_.pjet_btag5,               "pjet_btag5[pjet]/F");
    t_->Branch("pjet_btag6",             evSummary_.pjet_btag6,               "pjet_btag6[pjet]/F");
    t_->Branch("pjet_btag7",             evSummary_.pjet_btag7,               "pjet_btag7[pjet]/F");
    t_->Branch("pjet_btag8",             evSummary_.pjet_btag8,               "pjet_btag8[pjet]/F");
    t_->Branch("pjet_btag9",             evSummary_.pjet_btag9,               "pjet_btag9[pjet]/F");
    t_->Branch("pjet_btag10",            evSummary_.pjet_btag10,              "pjet_btag10[pjet]/F");


    //fjet (ak8PFJetsCHS)
    t_->Branch("fjet",                  &evSummary_.fjet,                   "fjet/I");
    t_->Branch("fjet_px",               evSummary_.fjet_px,                 "fjet_px[fjet]/F");
    t_->Branch("fjet_py",               evSummary_.fjet_py,                 "fjet_py[fjet]/F");
    t_->Branch("fjet_pz",               evSummary_.fjet_pz,                 "fjet_pz[fjet]/F");
    t_->Branch("fjet_en",               evSummary_.fjet_en,                 "fjet_en[fjet]/F");
    t_->Branch("fjet_genpt",            evSummary_.fjet_genpt,              "fjet_genpt[fjet]/F");
    t_->Branch("fjet_prunedM",          evSummary_.fjet_prunedM,            "fjet_prunedM[fjet]/F");
    t_->Branch("fjet_trimmedM",         evSummary_.fjet_trimmedM,           "fjet_trimmedM[fjet]/F");
    t_->Branch("fjet_filteredM",        evSummary_.fjet_filteredM,          "fjet_filteredM[fjet]/F");
    t_->Branch("fjet_tau1",             evSummary_.fjet_tau1,               "fjet_tau1[fjet]/F");
    t_->Branch("fjet_tau2",             evSummary_.fjet_tau2,               "fjet_tau2[fjet]/F");
    t_->Branch("fjet_tau3",             evSummary_.fjet_tau3,               "fjet_tau3[fjet]/F");


    //met
    t_->Branch("met_pt",               	&evSummary_.met_pt,                 "met_pt/F");
    t_->Branch("met_phi",               &evSummary_.met_phi,                "met_phi/F");
    t_->Branch("met_sumMET",            &evSummary_.met_sumMET,             "met_sumMET/F");
    t_->Branch("rawpfmet_pt",               	&evSummary_.rawpfmet_pt,                 "rawpfmet_pt/F");
    t_->Branch("rawpfmet_phi",               &evSummary_.rawpfmet_phi,                "rawpfmet_phi/F");
    t_->Branch("rawpfmet_sumMET",            &evSummary_.rawpfmet_sumMET,             "rawpfmet_sumMET/F");
    t_->Branch("rawcalomet_pt",               	&evSummary_.rawcalomet_pt,                 "rawcalomet_pt/F");
    t_->Branch("rawcalomet_phi",               &evSummary_.rawcalomet_phi,                "rawcalomet_phi/F");
    t_->Branch("rawcalomet_sumMET",            &evSummary_.rawcalomet_sumMET,             "rawcalomet_sumMET/F");
    t_->Branch("metNoHF_pt",                &evSummary_.metNoHF_pt,                 "metNoHF_pt/F");
    t_->Branch("metNoHF_phi",               &evSummary_.metNoHF_phi,                "metNoHF_phi/F");
    t_->Branch("metNoHF_sumMET",            &evSummary_.metNoHF_sumMET,             "metNoHF_sumMET/F");
    t_->Branch("metPuppi_pt",                &evSummary_.metPuppi_pt,                 "metPuppi_pt/F");
    t_->Branch("metPuppi_phi",               &evSummary_.metPuppi_phi,                "metPuppi_phi/F");
    t_->Branch("metPuppi_sumMET",            &evSummary_.metPuppi_sumMET,             "metPuppi_sumMET/F");

    /*
        t_->Branch("flag_HBHENoiseFilter",                      &evSummary_.flag_HBHENoiseFilter,                       "flag_HBHENoiseFilter/O");
        t_->Branch("flag_HBHENoiseIsoFilter",                   &evSummary_.flag_HBHENoiseIsoFilter,                    "flag_HBHENoiseIsoFilter/O");
        t_->Branch("flag_EcalDeadCellBoundaryEnergyFilter",     &evSummary_.flag_EcalDeadCellBoundaryEnergyFilter,      "flag_EcalDeadCellBoundaryEnergyFilter/O");
        t_->Branch("flag_CSCTightHaloFilter",                   &evSummary_.flag_CSCTightHaloFilter,                    "flag_CSCTightHaloFilter/O");
        t_->Branch("flag_hcalLaserEventFilter",                 &evSummary_.flag_hcalLaserEventFilter,                  "flag_hcalLaserEventFilter/O");
        t_->Branch("flag_EcalDeadCellTriggerPrimitiveFilter",   &evSummary_.flag_EcalDeadCellTriggerPrimitiveFilter,    "flag_EcalDeadCellTriggerPrimitiveFilter/O");
        t_->Branch("flag_goodVertices",                         &evSummary_.flag_goodVertices,                          "flag_goodVertices/O");
        t_->Branch("flag_trackingFailureFilter",                &evSummary_.flag_trackingFailureFilter,                 "flag_trackingFailureFilter/O");
        t_->Branch("flag_eeBadScFilter",                        &evSummary_.flag_eeBadScFilter,                         "flag_eeBadScFilter/O");
        t_->Branch("flag_ecalLaserCorrFilter",                  &evSummary_.flag_ecalLaserCorrFilter,                   "flag_ecalLaserCorrFilter/O");
        t_->Branch("flag_trkPOGFilters",                        &evSummary_.flag_trkPOGFilters,                         "flag_trkPOGFilters/O");
        t_->Branch("flag_trkPOG_manystripclus53X",              &evSummary_.flag_trkPOG_manystripclus53X,               "flag_trkPOG_manystripclus53X/O");
        t_->Branch("flag_trkPOG_toomanystripclus53X",           &evSummary_.flag_trkPOG_toomanystripclus53X,            "flag_trkPOG_toomanystripclus53X/O");
        t_->Branch("flag_trkPOG_logErrorTooManyClusters",       &evSummary_.flag_trkPOG_logErrorTooManyClusters,        "flag_trkPOG_logErrorTooManyClusters/O");
        t_->Branch("flag_METFilters",                           &evSummary_.flag_METFilters,                            "flag_METFilters/O");
    */

    return true;
}

//
bool DataEvtSummaryHandler::attachToTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    //event info
    t_->SetBranchAddress("run",             &evSummary_.run);
    t_->SetBranchAddress("lumi",            &evSummary_.lumi);
    t_->SetBranchAddress("event",           &evSummary_.event);
    t_->SetBranchAddress("curAvgInstLumi",  &evSummary_.curAvgInstLumi);
    t_->SetBranchAddress("curIntegLumi", 	&evSummary_.curIntegLumi);

    t_->SetBranchAddress("hasTrigger",   	&evSummary_.hasTrigger);
    t_->SetBranchAddress("triggerType",     &evSummary_.triggerType);

    //primary vertex
    t_->SetBranchAddress("nvtx",            &evSummary_.nvtx);
    t_->SetBranchAddress("vtx_x",           &evSummary_.vtx_x);
    t_->SetBranchAddress("vtx_y",           &evSummary_.vtx_y);
    t_->SetBranchAddress("vtx_z",           &evSummary_.vtx_z);

    t_->SetBranchAddress("fixedGridRhoAll",                           &evSummary_.fixedGridRhoAll);
    t_->SetBranchAddress("fixedGridRhoFastjetAll",                    &evSummary_.fixedGridRhoFastjetAll);
    t_->SetBranchAddress("fixedGridRhoFastjetAllCalo",                &evSummary_.fixedGridRhoFastjetAllCalo);
    t_->SetBranchAddress("fixedGridRhoFastjetCentralCalo",            &evSummary_.fixedGridRhoFastjetCentralCalo);
    t_->SetBranchAddress("fixedGridRhoFastjetCentralChargedPileUp",   &evSummary_.fixedGridRhoFastjetCentralChargedPileUp);
    t_->SetBranchAddress("fixedGridRhoFastjetCentralNeutral",         &evSummary_.fixedGridRhoFastjetCentralNeutral);

    //generator level info
    t_->SetBranchAddress("ngenITpu",        &evSummary_.ngenITpu);
    t_->SetBranchAddress("ngenOOTpu",       &evSummary_.ngenOOTpu);
    t_->SetBranchAddress("ngenOOTpum1",  	&evSummary_.ngenOOTpum1);
    t_->SetBranchAddress("ngenTruepu",  	&evSummary_.ngenTruepu);
    t_->SetBranchAddress("puweight",     	&evSummary_.puWeight);
    t_->SetBranchAddress("hptWeights",      evSummary_.hptWeights);
    t_->SetBranchAddress("pthat",           &evSummary_.pthat);
    t_->SetBranchAddress("genWeight",       &evSummary_.genWeight);
    t_->SetBranchAddress("qscale",          &evSummary_.qscale);
    t_->SetBranchAddress("x1",              &evSummary_.x1);
    t_->SetBranchAddress("x2",              &evSummary_.x2);
    t_->SetBranchAddress("id1",             &evSummary_.id1);
    t_->SetBranchAddress("id2",             &evSummary_.id2);


    //mc truth
    t_->SetBranchAddress("nmcparticles", 	&evSummary_.nmcparticles);
    t_->SetBranchAddress("mc_px",           evSummary_.mc_px);
    t_->SetBranchAddress("mc_py",           evSummary_.mc_py);
    t_->SetBranchAddress("mc_pz",           evSummary_.mc_pz);
    t_->SetBranchAddress("mc_en",           evSummary_.mc_en);
    t_->SetBranchAddress("mc_id",           evSummary_.mc_id);
    t_->SetBranchAddress("mc_status",       evSummary_.mc_status);
    t_->SetBranchAddress("mc_mom",          evSummary_.mc_mom);


    //muon
    t_->SetBranchAddress("mn",                      &evSummary_.mn);
    t_->SetBranchAddress("mn_px",                   evSummary_.mn_px);
    t_->SetBranchAddress("mn_py",                   evSummary_.mn_py);
    t_->SetBranchAddress("mn_pz",                   evSummary_.mn_pz);
    t_->SetBranchAddress("mn_en",                   evSummary_.mn_en);
    t_->SetBranchAddress("mn_id",                   evSummary_.mn_id);
    t_->SetBranchAddress("mn_type",                 evSummary_.mn_type);
    t_->SetBranchAddress("mn_d0",                   evSummary_.mn_d0);
    t_->SetBranchAddress("mn_dZ",                   evSummary_.mn_dZ);
    t_->SetBranchAddress("mn_ip3d",                 evSummary_.mn_ip3d);
    t_->SetBranchAddress("mn_ip3dsig",              evSummary_.mn_ip3dsig);
    t_->SetBranchAddress("mn_IsLoose",              evSummary_.mn_IsLoose);
    t_->SetBranchAddress("mn_IsTight",              evSummary_.mn_IsTight);
    t_->SetBranchAddress("mn_IsSoft",               evSummary_.mn_IsSoft);
    t_->SetBranchAddress("mn_IsHighPt",             evSummary_.mn_IsHighPt);
    t_->SetBranchAddress("mn_pileupIsoR03",         evSummary_.mn_pileupIsoR03);
    t_->SetBranchAddress("mn_chargedIsoR03",        evSummary_.mn_chargedIsoR03);
    t_->SetBranchAddress("mn_photonIsoR03",         evSummary_.mn_photonIsoR03);
    t_->SetBranchAddress("mn_neutralHadIsoR03",     evSummary_.mn_neutralHadIsoR03);
    t_->SetBranchAddress("mn_pileupIsoR04",         evSummary_.mn_pileupIsoR04);
    t_->SetBranchAddress("mn_chargedIsoR04",        evSummary_.mn_chargedIsoR04);
    t_->SetBranchAddress("mn_photonIsoR04",         evSummary_.mn_photonIsoR04);
    t_->SetBranchAddress("mn_neutralHadIsoR04",     evSummary_.mn_neutralHadIsoR04);

    /*
     t_->SetBranchAddress("mn_nMatches",                   evSummary_.mn_nMatches);
     t_->SetBranchAddress("mn_nMatchedStations",           evSummary_.mn_nMatchedStations);
     t_->SetBranchAddress("mn_validMuonHits",              evSummary_.mn_validMuonHits);
     t_->SetBranchAddress("mn_innerTrackChi2",             evSummary_.mn_innerTrackChi2);
     t_->SetBranchAddress("mn_trkLayersWithMeasurement",   evSummary_.mn_trkLayersWithMeasurement);
     t_->SetBranchAddress("mn_pixelLayersWithMeasurement", evSummary_.mn_pixelLayersWithMeasurement);
     */

    //electron
    t_->SetBranchAddress("en",                      &evSummary_.en);
    t_->SetBranchAddress("en_px",                   evSummary_.en_px);
    t_->SetBranchAddress("en_py",                   evSummary_.en_py);
    t_->SetBranchAddress("en_pz",                   evSummary_.en_pz);
    t_->SetBranchAddress("en_en",                   evSummary_.en_en);
    t_->SetBranchAddress("en_id",                   evSummary_.en_id);
    t_->SetBranchAddress("en_d0",                   evSummary_.en_d0);
    t_->SetBranchAddress("en_dZ",                   evSummary_.en_dZ);
    t_->SetBranchAddress("en_EtaSC",                evSummary_.en_EtaSC);
    t_->SetBranchAddress("en_PhiSC",                evSummary_.en_PhiSC);
    t_->SetBranchAddress("en_EnSC",                 evSummary_.en_EnSC);
    t_->SetBranchAddress("en_dEtaIn",               evSummary_.en_dEtaIn);
    t_->SetBranchAddress("en_dPhiIn",               evSummary_.en_dPhiIn);
    t_->SetBranchAddress("en_hOverE",               evSummary_.en_hOverE);
    t_->SetBranchAddress("en_R9",                   evSummary_.en_R9);
    t_->SetBranchAddress("en_sigmaIetaIeta",        evSummary_.en_sigmaIetaIeta);
    t_->SetBranchAddress("en_sigmaIetaIeta5x5",     evSummary_.en_sigmaIetaIeta5x5);
    t_->SetBranchAddress("en_ooEmooP",          	evSummary_.en_ooEmooP);
    t_->SetBranchAddress("en_pileupIso",            evSummary_.en_pileupIso);
    t_->SetBranchAddress("en_chargedIso",           evSummary_.en_chargedIso);
    t_->SetBranchAddress("en_photonIso",            evSummary_.en_photonIso);
    t_->SetBranchAddress("en_neutralHadIso",        evSummary_.en_neutralHadIso);
    t_->SetBranchAddress("en_relIsoWithEA",         evSummary_.en_relIsoWithEA);
    t_->SetBranchAddress("en_relIsoWithDBeta",      evSummary_.en_relIsoWithDBeta);
    t_->SetBranchAddress("en_MissingHits",          evSummary_.en_MissingHits);
    t_->SetBranchAddress("en_passConversionVeto",   evSummary_.en_passConversionVeto);
    t_->SetBranchAddress("en_passVeto",             evSummary_.en_passVeto);
    t_->SetBranchAddress("en_passLoose",            evSummary_.en_passLoose);
    t_->SetBranchAddress("en_passMedium",           evSummary_.en_passMedium);
    t_->SetBranchAddress("en_passTight",            evSummary_.en_passTight);
    t_->SetBranchAddress("en_passHEEP",             evSummary_.en_passHEEP);
    t_->SetBranchAddress("en_passMVATrigMedium",    evSummary_.en_passMVATrigMedium);
    t_->SetBranchAddress("en_passMVATrigTight",     evSummary_.en_passMVATrigTight);
    t_->SetBranchAddress("en_IDMVATrigValue",       evSummary_.en_IDMVATrigValue);
    t_->SetBranchAddress("en_IDMVATrigCategory",    evSummary_.en_IDMVATrigCategory);
    t_->SetBranchAddress("en_istrue",               evSummary_.en_istrue);



    //tau
    t_->SetBranchAddress("ta",                      &evSummary_.ta);
    t_->SetBranchAddress("ta_px",                   evSummary_.ta_px);
    t_->SetBranchAddress("ta_py",                   evSummary_.ta_py);
    t_->SetBranchAddress("ta_pz",                   evSummary_.ta_pz);
    t_->SetBranchAddress("ta_en",                   evSummary_.ta_en);
    t_->SetBranchAddress("ta_id",                   evSummary_.ta_id);
    t_->SetBranchAddress("ta_dm",                   evSummary_.ta_dm);
    t_->SetBranchAddress("ta_newdm",                evSummary_.ta_newdm);
    t_->SetBranchAddress("ta_IsLooseIso",           evSummary_.ta_IsLooseIso);
    t_->SetBranchAddress("ta_IsMediumIso",          evSummary_.ta_IsMediumIso);
    t_->SetBranchAddress("ta_IsTightIso",           evSummary_.ta_IsTightIso);
    t_->SetBranchAddress("ta_combIsoDBeta3Hits",    evSummary_.ta_combIsoDBeta3Hits);
    t_->SetBranchAddress("ta_chargedIso",           evSummary_.ta_chargedIso);
    t_->SetBranchAddress("ta_neutralIso",           evSummary_.ta_neutralIso);
    t_->SetBranchAddress("ta_pileupIso",            evSummary_.ta_pileupIso);
    t_->SetBranchAddress("ta_passEleVetoLoose",     evSummary_.ta_passEleVetoLoose);
    t_->SetBranchAddress("ta_passEleVetoMedium",    evSummary_.ta_passEleVetoMedium);
    t_->SetBranchAddress("ta_passEleVetoTight",     evSummary_.ta_passEleVetoTight);
    t_->SetBranchAddress("ta_passMuVetoLoose3",     evSummary_.ta_passMuVetoLoose3);
    t_->SetBranchAddress("ta_passMuVetoTight3",     evSummary_.ta_passMuVetoTight3);



    //jet (ak4PFJetsCHS)
    t_->SetBranchAddress("jet",                     &evSummary_.jet);
    t_->SetBranchAddress("jet_px",                  evSummary_.jet_px);
    t_->SetBranchAddress("jet_py",                  evSummary_.jet_py);
    t_->SetBranchAddress("jet_pz",                  evSummary_.jet_pz);
    t_->SetBranchAddress("jet_en",                  evSummary_.jet_en);
    t_->SetBranchAddress("jet_btag0",               evSummary_.jet_btag0);
    t_->SetBranchAddress("jet_btag1",               evSummary_.jet_btag1);
    t_->SetBranchAddress("jet_btag2",               evSummary_.jet_btag2);
    t_->SetBranchAddress("jet_btag3",               evSummary_.jet_btag3);
    t_->SetBranchAddress("jet_btag4",               evSummary_.jet_btag4);
    t_->SetBranchAddress("jet_btag5",               evSummary_.jet_btag5);
    t_->SetBranchAddress("jet_btag6",               evSummary_.jet_btag6);
    t_->SetBranchAddress("jet_btag7",               evSummary_.jet_btag7);
    t_->SetBranchAddress("jet_btag8",               evSummary_.jet_btag8);
    t_->SetBranchAddress("jet_btag9",               evSummary_.jet_btag9);
    t_->SetBranchAddress("jet_btag10",               evSummary_.jet_btag10);
    t_->SetBranchAddress("jet_mass",                evSummary_.jet_mass);
    t_->SetBranchAddress("jet_area",                evSummary_.jet_area);
    t_->SetBranchAddress("jet_pu",                  evSummary_.jet_pu);
    t_->SetBranchAddress("jet_puId",                evSummary_.jet_puId);
    t_->SetBranchAddress("jet_PFLoose",             evSummary_.jet_PFLoose);
    t_->SetBranchAddress("jet_PFTight",             evSummary_.jet_PFTight);
    t_->SetBranchAddress("jet_partonFlavour",       evSummary_.jet_partonFlavour);
    t_->SetBranchAddress("jet_hadronFlavour",       evSummary_.jet_hadronFlavour);
    t_->SetBranchAddress("jet_genpt",               evSummary_.jet_genpt);


    //pjet: slimmedJetsPuppi
    t_->SetBranchAddress("pjet",                     &evSummary_.pjet);
    t_->SetBranchAddress("pjet_px",                  evSummary_.pjet_px);
    t_->SetBranchAddress("pjet_py",                  evSummary_.pjet_py);
    t_->SetBranchAddress("pjet_pz",                  evSummary_.pjet_pz);
    t_->SetBranchAddress("pjet_en",                  evSummary_.pjet_en);
    t_->SetBranchAddress("pjet_genpt",               evSummary_.pjet_genpt);
    t_->SetBranchAddress("pjet_btag0",               evSummary_.pjet_btag0);
    t_->SetBranchAddress("pjet_btag1",               evSummary_.pjet_btag1);
    t_->SetBranchAddress("pjet_btag2",               evSummary_.pjet_btag2);
    t_->SetBranchAddress("pjet_btag3",               evSummary_.pjet_btag3);
    t_->SetBranchAddress("pjet_btag4",               evSummary_.pjet_btag4);
    t_->SetBranchAddress("pjet_btag5",               evSummary_.pjet_btag5);
    t_->SetBranchAddress("pjet_btag6",               evSummary_.pjet_btag6);
    t_->SetBranchAddress("pjet_btag7",               evSummary_.pjet_btag7);
    t_->SetBranchAddress("pjet_btag8",               evSummary_.pjet_btag8);
    t_->SetBranchAddress("pjet_btag9",               evSummary_.pjet_btag9);
    t_->SetBranchAddress("pjet_btag10",               evSummary_.pjet_btag10);



    //fjet (ak8PFJetsCHS)
    t_->SetBranchAddress("fjet",                    &evSummary_.fjet);
    t_->SetBranchAddress("fjet_px",                 evSummary_.fjet_px);
    t_->SetBranchAddress("fjet_py",                 evSummary_.fjet_py);
    t_->SetBranchAddress("fjet_pz",                 evSummary_.fjet_pz);
    t_->SetBranchAddress("fjet_en",                 evSummary_.fjet_en);
    t_->SetBranchAddress("fjet_genpt",              evSummary_.fjet_genpt);
    t_->SetBranchAddress("fjet_prunedM",            evSummary_.fjet_prunedM);
    t_->SetBranchAddress("fjet_trimmedM",           evSummary_.fjet_trimmedM);
    t_->SetBranchAddress("fjet_filteredM",          evSummary_.fjet_filteredM);
    t_->SetBranchAddress("fjet_tau1",               evSummary_.fjet_tau1);
    t_->SetBranchAddress("fjet_tau2",               evSummary_.fjet_tau2);
    t_->SetBranchAddress("fjet_tau3",               evSummary_.fjet_tau3);


    //met
    t_->SetBranchAddress("met_pt",                  &evSummary_.met_pt);
    t_->SetBranchAddress("met_phi",                 &evSummary_.met_phi);
    t_->SetBranchAddress("met_sumMET",              &evSummary_.met_sumMET);
    t_->SetBranchAddress("rawpfmet_pt",                  &evSummary_.rawpfmet_pt);
    t_->SetBranchAddress("rawpfmet_phi",                 &evSummary_.rawpfmet_phi);
    t_->SetBranchAddress("rawpfmet_sumMET",              &evSummary_.rawpfmet_sumMET);
    t_->SetBranchAddress("rawcalomet_pt",                  &evSummary_.rawcalomet_pt);
    t_->SetBranchAddress("rawcalomet_phi",                 &evSummary_.rawcalomet_phi);
    t_->SetBranchAddress("rawcalomet_sumMET",              &evSummary_.rawcalomet_sumMET);
    t_->SetBranchAddress("metNoHF_pt",                  &evSummary_.metNoHF_pt);
    t_->SetBranchAddress("metNoHF_phi",                 &evSummary_.metNoHF_phi);
    t_->SetBranchAddress("metNoHF_sumMET",              &evSummary_.metNoHF_sumMET);
    t_->SetBranchAddress("metPuppi_pt",                  &evSummary_.metPuppi_pt);
    t_->SetBranchAddress("metPuppi_phi",                 &evSummary_.metPuppi_phi);
    t_->SetBranchAddress("metPuppi_sumMET",              &evSummary_.metPuppi_sumMET);


    /*
        t_->SetBranchAddress("flag_HBHENoiseFilter",                      &evSummary_.flag_HBHENoiseFilter);
        t_->SetBranchAddress("flag_HBHENoiseIsoFilter",                   &evSummary_.flag_HBHENoiseIsoFilter);
        t_->SetBranchAddress("flag_EcalDeadCellBoundaryEnergyFilter",     &evSummary_.flag_EcalDeadCellBoundaryEnergyFilter);
        t_->SetBranchAddress("flag_CSCTightHaloFilter",                   &evSummary_.flag_CSCTightHaloFilter);
        t_->SetBranchAddress("flag_hcalLaserEventFilter",                 &evSummary_.flag_hcalLaserEventFilter);
        t_->SetBranchAddress("flag_EcalDeadCellTriggerPrimitiveFilter",   &evSummary_.flag_EcalDeadCellTriggerPrimitiveFilter);
        t_->SetBranchAddress("flag_goodVertices",                         &evSummary_.flag_goodVertices);
        t_->SetBranchAddress("flag_trackingFailureFilter",                &evSummary_.flag_trackingFailureFilter);
        t_->SetBranchAddress("flag_eeBadScFilter",                        &evSummary_.flag_eeBadScFilter);
        t_->SetBranchAddress("flag_ecalLaserCorrFilter",                  &evSummary_.flag_ecalLaserCorrFilter);
        t_->SetBranchAddress("flag_trkPOGFilters",                        &evSummary_.flag_trkPOGFilters);
        t_->SetBranchAddress("flag_trkPOG_manystripclus53X",              &evSummary_.flag_trkPOG_manystripclus53X);
        t_->SetBranchAddress("flag_trkPOG_toomanystripclus53X",           &evSummary_.flag_trkPOG_toomanystripclus53X);
        t_->SetBranchAddress("flag_trkPOG_logErrorTooManyClusters",       &evSummary_.flag_trkPOG_logErrorTooManyClusters);
        t_->SetBranchAddress("flag_METFilters",                           &evSummary_.flag_METFilters);
    */


    return true;
}


//
void DataEvtSummaryHandler::resetStruct()
{
    evSummary_.nmcparticles=0;
    evSummary_.run=0;
    evSummary_.lumi=0;
    evSummary_.event=0;
    evSummary_.mn=0;
    evSummary_.en=0;
    evSummary_.ta=0;
    evSummary_.jet=0;
    evSummary_.pjet=0;
    evSummary_.fjet=0;
}

//
void DataEvtSummaryHandler::fillTree()
{
    if(t_) t_->Fill();
}

//
DataEvtSummaryHandler::~DataEvtSummaryHandler()
{
}
