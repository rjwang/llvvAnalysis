import FWCore.ParameterSet.Config as cms

process = cms.Process("DMAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )



process.mainAnalyzer = cms.EDAnalyzer('MainAnalyzer',
    dtag = cms.string('llvv'),
    isMC = cms.bool(True),
    isPythia8 = cms.bool(True),

    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    verticesTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

    rhoAll = cms.InputTag("fixedGridRhoAll"),
    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),

    muonsTag = cms.InputTag("slimmedMuons"),

    electronsTag = cms.InputTag("slimmedElectrons"),
    	electronVetoIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
	electronLooseIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
	electronMediumIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    	electronTightIdTag = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),

    tausTag = cms.InputTag("slimmedTaus"),


    photonsTag = cms.InputTag("slimmedPhotons"),
    jetsTag = cms.InputTag("slimmedJets"),
    fatjetsTag = cms.InputTag("slimmedJetsAK8"),
    metsTag = cms.InputTag("slimmedMETs"),
    metFilterBitsTag = cms.InputTag("TriggerResults"),
    packedTag = cms.InputTag("packedGenParticles"),
    prunedTag = cms.InputTag("prunedGenParticles"),
    genJetsTag = cms.InputTag("slimmedGenJets"),


    ##trigger
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),

    	DoubleMuTrigs = cms.vstring("HLT_Mu17_Mu8_v",
					"HLT_Mu17_TkMu8_v",
					"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",
					"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"
				   ),

    	DoubleEleTrigs = cms.vstring("HLT_Ele23_Ele12_CaloId_TrackId_Iso_v",
					"HLT_Ele20WP60_Ele8_Mass55_v",
					"HLT_Ele25WP60_SC4_Mass55_v",
					"HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v",
					"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v",
					"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v"),

    	SingleMuTrigs = cms.vstring("HLT_IsoMu24_IterTrk02_v",
					"HLT_IsoTkMu24_IterTrk02_v",
					"HLT_IsoMu20_eta2p1_IterTrk02_v",
					"HLT_IsoTkMu20_eta2p1_IterTrk02_v",
					"HLT_Mu40_v"),

    	SingleEleTrigs = cms.vstring("HLT_Ele27_eta2p1_WP85_Gsf_v",
					"HLT_Ele95_CaloIdVT_GsfTrkIdT_v"),

   	MuEGTrigs = cms.vstring("HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v",
					"HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v"),

	#DoubleTauTrigs = cms.vstring("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v")
)
