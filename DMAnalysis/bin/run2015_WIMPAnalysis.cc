#include <iostream>

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "llvvAnalysis/DMAnalysis/interface/MacroUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/DataEvtSummaryHandler.h"
#include "llvvAnalysis/DMAnalysis/interface/DMPhysicsEvent.h"
#include "llvvAnalysis/DMAnalysis/interface/SmartSelectionMonitor.h"
#include "llvvAnalysis/DMAnalysis/interface/METUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/BTagUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/WIMPReweighting.h"
#include "llvvAnalysis/DMAnalysis/interface/EventCategory.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"


#include "llvvAnalysis/DMAnalysis/interface/PDFInfo.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
#include "llvvAnalysis/DMAnalysis/interface/BTagCalibrationStandalone.h"

using namespace std;

int main(int argc, char* argv[])
{
    //##################################################################################
    //##########################    GLOBAL INITIALIZATION     ##########################
    //##################################################################################

    // check arguments
    if(argc<2) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        exit(0);
    }

    // load framework libraries
    gSystem->Load( "libFWCoreFWLite" );
    AutoLibraryLoader::enable();

    // configure the process
    const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");

    bool isMC       = runProcess.getParameter<bool>("isMC");
    int mctruthmode = runProcess.getParameter<int>("mctruthmode");
    bool doWIMPreweighting = runProcess.getParameter<bool>("doWIMPreweighting");
    bool usemetNoHF = runProcess.getParameter<bool>("usemetNoHF");

    TString url=runProcess.getParameter<std::string>("input");
    TString outFileUrl(gSystem->BaseName(url));
    outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }
    TString outdir=runProcess.getParameter<std::string>("outdir");
    TString outUrl( outdir );
    gSystem->Exec("mkdir -p " + outUrl);


    TString outTxtUrl_final= outUrl + "/" + outFileUrl + "_FinalList.txt";
    FILE* outTxtFile_final = NULL;
    outTxtFile_final = fopen(outTxtUrl_final.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl_final.Data());
    fprintf(outTxtFile_final,"run lumi event\n");

    int fType(0);
    if(url.Contains("DoubleEG")) fType=EE;
    if(url.Contains("DoubleMuon"))  fType=MUMU;
    if(url.Contains("MuonEG"))      fType=EMU;
    if(url.Contains("SingleMuon"))  fType=MUMU;
    if(url.Contains("SingleElectron")) fType=EE;
    bool isSingleMuPD(!isMC && url.Contains("SingleMuon"));
    bool isDoubleMuPD(!isMC && url.Contains("DoubleMuon"));
    bool isSingleElePD(!isMC && url.Contains("SingleElectron"));
    bool isDoubleElePD(!isMC && url.Contains("DoubleEG"));

    bool isMC_ZZ  = isMC && ( string(url.Data()).find("TeV_ZZTo2L2Nu")  != string::npos);
//  bool isMC_WZ  = isMC && ( string(url.Data()).find("TeV_WZ_")  != string::npos);
//  bool isMC_WW  = isMC && ( string(url.Data()).find("TeV_WW_")  != string::npos);
    bool isMC_ttbar = isMC && (string(url.Data()).find("TeV_TT")  != string::npos);
    bool isMC_stop  = isMC && (string(url.Data()).find("TeV_SingleT")  != string::npos);
    bool isMC_WIMP  = isMC && (string(url.Data()).find("TeV_DM_V_Mx") != string::npos
                               || string(url.Data()).find("TeV_DM_A_Mx") != string::npos
                               || string(url.Data()).find("TeV_EWKDM_S_Mx") != string::npos
                               || string(url.Data()).find("TeV_EWKDM_P_Mx") != string::npos);
    bool isMC_ADD  = isMC && (string(url.Data()).find("TeV_ADD_D") != string::npos);
    bool isMC_Unpart = isMC && (string(url.Data()).find("TeV_Unpart") != string::npos);


    bool isSignal = (isMC_WIMP || isMC_ADD || isMC_Unpart);


    BTagUtils myBtagUtils(runProcess);


    // https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
    // setup calibration readers
    TString BtagSF = runProcess.getParameter<std::string>("BtagSF");
    gSystem->ExpandPathName(BtagSF);
    cout << "Loading btag sacle factor: " << BtagSF << endl;
    BTagCalibration btagcalib("csvv2", BtagSF.Data());
    BTagCalibrationReader btag_reader(&btagcalib, BTagEntry::OP_LOOSE, "mujets", "central");
    BTagCalibrationReader btag_reader_up(&btagcalib, BTagEntry::OP_LOOSE, "mujets", "up");  // sys up
    BTagCalibrationReader btag_reader_down(&btagcalib, BTagEntry::OP_LOOSE, "mujets", "down");  // sys down



    WIMPReweighting myWIMPweights(runProcess);


    //systematics
    bool runSystematics                        = runProcess.getParameter<bool>("runSystematics");
    std::vector<TString> varNames(1,"");
    if(runSystematics) {
        cout << "Systematics will be computed for this analysis" << endl;
        varNames.push_back("_jerup"); 	//1
        varNames.push_back("_jerdown"); //2
        varNames.push_back("_jesup"); 	//3
        varNames.push_back("_jesdown"); //4
        varNames.push_back("_umetup"); 	//5
        varNames.push_back("_umetdown");//6
        varNames.push_back("_lesup"); 	//7
        varNames.push_back("_lesdown"); //8
        varNames.push_back("_puup"); 	//9
        varNames.push_back("_pudown"); 	//10
        varNames.push_back("_btagup"); 	//11
        varNames.push_back("_btagdown");//12
        if(isSignal)            {
            varNames.push_back("_pdfup");
            varNames.push_back("_pdfdown");
        }
    }
    size_t nvarsToInclude=varNames.size();



    //tree info
    int evStart     = runProcess.getParameter<int>("evStart");
    int evEnd       = runProcess.getParameter<int>("evEnd");
    TString dirname = runProcess.getParameter<std::string>("dirName");

    //jet energy scale uncertainties
    TString uncFile = runProcess.getParameter<std::string>("jesUncFileName");
    gSystem->ExpandPathName(uncFile);
    cout << "Loading jet energy scale uncertainties: " << uncFile << endl;
    JetCorrectionUncertainty jecUnc(uncFile.Data());



    //pdf info
    PDFInfo *mPDFInfo=0;
    if(isSignal) {
        TString pdfUrl = runProcess.getParameter<std::string>("pdfInput");
        std::string Url = runProcess.getParameter<std::string>("input");
        if(doWIMPreweighting) {
            if(Url.find("TeV_DM_V_Mx") != string::npos) Url = "root://eoscms//eos/cms/store/group/phys_exotica/monoZ/llvvNtuple_25ns_29Oct2015/MC13TeV_DM_V_Mx1Mv200.root";
            if(Url.find("TeV_DM_A_Mx") != string::npos) Url = "root://eoscms//eos/cms/store/group/phys_exotica/monoZ/llvvNtuple_25ns_29Oct2015/MC13TeV_DM_A_Mx1Mv200.root";

        }

        std::size_t found = Url.find_last_of("/\\");
        pdfUrl += '/';
        pdfUrl += Url.substr(found+1);


        if(pdfUrl.Contains("K1_0.1_K2_1")) pdfUrl.ReplaceAll("K1_0.1_K2_1","K1_1_K2_1");
        if(pdfUrl.Contains("K1_0.2_K2_1")) pdfUrl.ReplaceAll("K1_0.2_K2_1","K1_1_K2_1");
        if(pdfUrl.Contains("K1_0.3_K2_1")) pdfUrl.ReplaceAll("K1_0.3_K2_1","K1_1_K2_1");
        if(pdfUrl.Contains("K1_0.5_K2_1")) pdfUrl.ReplaceAll("K1_0.5_K2_1","K1_1_K2_1");
        if(pdfUrl.Contains("K1_2_K2_1"))   pdfUrl.ReplaceAll("K1_2_K2_1","K1_1_K2_1");
        if(pdfUrl.Contains("K1_3_K2_1"))   pdfUrl.ReplaceAll("K1_3_K2_1","K1_1_K2_1");
        if(pdfUrl.Contains("K1_5_K2_1"))   pdfUrl.ReplaceAll("K1_5_K2_1","K1_1_K2_1");
        if(pdfUrl.Contains("K1_10_K2_1"))  pdfUrl.ReplaceAll("K1_10_K2_1","K1_1_K2_1");
        pdfUrl.ReplaceAll(".root","_pdf.root");

        mPDFInfo=new PDFInfo(pdfUrl,"cteq66.LHgrid");
        cout << "Readout " << mPDFInfo->numberPDFs() << " pdf variations: " << pdfUrl << endl;
    }


    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;


    TH1F *h=(TH1F*) mon.addHistogram( new TH1F ("eventflow", ";;Events", 10,0,10) );
    h->GetXaxis()->SetBinLabel(1,"Trigger && 2 leptons");
    h->GetXaxis()->SetBinLabel(2,"|#it{m}_{ll}-#it{m}_{Z}|<15");
    h->GetXaxis()->SetBinLabel(3,"#it{p}_{T}^{ll}>50");
    h->GetXaxis()->SetBinLabel(4,"3^{rd}-lepton veto");
    h->GetXaxis()->SetBinLabel(5,"b-veto");
    h->GetXaxis()->SetBinLabel(6,"#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss})>2.7");
    h->GetXaxis()->SetBinLabel(7,"|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T}<0.2");
    h->GetXaxis()->SetBinLabel(8,"E_{T}^{miss}>80");


    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

    mon.addHistogram( new TH1F( "nvtx_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxwgt_raw",	";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "zpt_raw",      ";#it{p}_{T}^{ll} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "pfmet_raw",    ";E_{T}^{miss} [GeV];Events", 50,0,500) );
    mon.addHistogram( new TH1F( "mt_raw",       ";#it{m}_{T} [GeV];Events", 100,0,2000) );
    double MTBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000,2000};
    const int nBinsMT = sizeof(MTBins)/sizeof(double) - 1;
    mon.addHistogram( new TH1F( "mt2_raw",       ";#it{m}_{T} [GeV];Events", nBinsMT,MTBins) );
    mon.addHistogram( new TH1F( "zmass_raw",    ";#it{m}_{ll} [GeV];Events", 100,40,250) );

    mon.addHistogram( new TH2F( "ptlep1vs2_raw",";#it{p}_{T}^{l1} [GeV];#it{p}_{T}^{l2} [GeV];Events",250,0,500, 250,0,500) );

    mon.addHistogram( new TH1F( "leadlep_pt_raw", ";Leading lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "leadlep_eta_raw",";Leading lepton #eta^{l};Events", 50,-2.6,2.6) );
    mon.addHistogram( new TH1F( "trailep_pt_raw", ";Trailing lepton #it{p}_{T}^{l};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "trailep_eta_raw",";Trailing lepton #eta^{l};Events", 50,-2.6,2.6) );


    mon.addHistogram( new TH1F( "jet_pt_raw", ";all jet #it{p}_{T}^{j};Events", 50,0,500) );
    mon.addHistogram( new TH1F( "jet_eta_raw",";all jet #eta^{j};Events", 50,-2.6,2.6) );


    TH1F *h1 = (TH1F*) mon.addHistogram( new TH1F( "nleptons_raw", ";Lepton multiplicity;Events", 3,2,5) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        label +="=";
        label += (ibin+1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }

    TH1F *h2 = (TH1F *)mon.addHistogram( new TH1F("njets_raw",  ";Jet multiplicity (#it{p}_{T}>30 GeV);Events",5,0,5) );
    for(int ibin=1; ibin<=h2->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h2->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h2->GetXaxis()->SetBinLabel(ibin,label);
    }

    TH1F *h3 = (TH1F *)mon.addHistogram( new TH1F("nbjets_raw",    ";b-tag Jet multiplicity;Events",5,0,5) );
    for(int ibin=1; ibin<=h3->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        if(ibin==h3->GetXaxis()->GetNbins()) label +="#geq";
        else                                label +="=";
        label += (ibin-1);
        h3->GetXaxis()->SetBinLabel(ibin,label);
    }


    // preselection plots
    double METBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500};
    const int nBinsMET = sizeof(METBins)/sizeof(double) - 1;

    double METBins2[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,500,1000};
    const int nBinsMET2 = sizeof(METBins2)/sizeof(double) - 1;

    mon.addHistogram( new TH1F( "pfmet_presel",      ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET, METBins));
    mon.addHistogram( new TH1F( "pfmet2_presel",     ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET2, METBins2));
    mon.addHistogram( new TH1F( "dphiZMET_presel",   ";#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss});Events", 10,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "balancedif_presel", ";|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T};Events", 5,0,1.0) );
    mon.addHistogram( new TH1F( "mt_presel",         ";#it{m}_{T} [GeV];Events", 12,0,1200) );
    mon.addHistogram( new TH1F( "mt2_presel",         ";#it{m}_{T} [GeV];Events", nBinsMT,MTBins) );
    mon.addHistogram( new TH1F( "axialpfmet_presel", ";Axial E_{T}^{miss} [GeV];Events", 200,-150,150) );

    //adding N-1 plots


    //MET X-Y shift correction
    mon.addHistogram( new TH2F( "pfmetx_vs_nvtx_presel",";Vertices;E_{X}^{miss} [GeV];Events",50,0,50, 200,-75,75) );
    mon.addHistogram( new TH2F( "pfmety_vs_nvtx_presel",";Vertices;E_{Y}^{miss} [GeV];Events",50,0,50, 200,-75,75) );
    mon.addHistogram( new TH1F( "pfmetphi_wocorr_presel",";#it{#phi}(E_{T}^{miss});Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "pfmetphi_wicorr_presel",";#it{#phi}(E_{T}^{miss});Events", 50,-1.*TMath::Pi(),TMath::Pi()) );

    mon.addHistogram( new TH1F( "pfmet_wicorr_presel",      ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET, METBins));
    mon.addHistogram( new TH1F( "pfmet2_wicorr_presel",     ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET2, METBins2));


    // generator level plots
    mon.addHistogram( new TH1F( "pileup", ";pileup;Events", 50,0,50) );
    mon.addHistogram( new TH1F( "met_Gen", ";#it{p}_{T}(#bar{#chi}#chi) [GeV];Events", nBinsMET, METBins) );
    mon.addHistogram( new TH1F( "met2_Gen", ";#it{p}_{T}(#bar{#chi}#chi) [GeV];Events", 500, 0,1000) );
    mon.addHistogram( new TH1F( "zpt_Gen", ";#it{p}_{T}(Z) [GeV];Events", 800,0,800) );
    mon.addHistogram( new TH1F( "dphi_Gen", ";#Delta#phi(Z,#bar{#chi}#chi) [rad];Events", 100,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "zmass_Gen", ";#it{m}_{ll} [GeV] [GeV];Events", 250,0,250) );
    mon.addHistogram( new TH2F( "ptlep1vs2_Gen",";#it{p}_{T}^{l1} [GeV];#it{p}_{T}^{l2} [GeV];Events",250,0,500, 250,0,500) );

    h=(TH1F *)mon.addHistogram( new TH1F ("acceptance", ";;Events", 2,0,2) );
    h->GetXaxis()->SetBinLabel(1,"Gen");
    h->GetXaxis()->SetBinLabel(2,"Gen Acc");


    // btaging efficiency
    std::vector<TString> CSVkey;
    CSVkey.push_back("CSVL");
    CSVkey.push_back("CSVM");
    CSVkey.push_back("CSVT");
    for(size_t csvtag=0; csvtag<CSVkey.size(); csvtag++) {
        mon.addHistogram( new TH1F( TString("beff_Denom_")+CSVkey[csvtag],      "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("ceff_Denom_")+CSVkey[csvtag],      "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("udsgeff_Denom_")+CSVkey[csvtag],   "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("beff_Num_")+CSVkey[csvtag],        "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("ceff_Num_")+CSVkey[csvtag],        "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
        mon.addHistogram( new TH1F( TString("udsgeff_Num_")+CSVkey[csvtag],     "; Jet #it{p}_{T} [GeV];Events", 20,20,500) );
    }


    mon.addHistogram( new TH1F( "mt_final",             ";#it{m}_{T} [GeV];Events", 12,0,1200) );
    mon.addHistogram( new TH1F( "mt_final120",             ";#it{m}_{T} [GeV];Events", 12,0,1200) );
    mon.addHistogram( new TH1F( "pfmet_final",      ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET, METBins));
    mon.addHistogram( new TH1F( "pfmet2_final",     ";E_{T}^{miss} [GeV];Events / 1 GeV", nBinsMET2, METBins2));




    //#################################################
    //############# CONTROL PLOTS #####################
    //#################################################
    // WW control plots, for k-method (for emu channel)
    mon.addHistogram( new TH1F( "zpt_wwctrl_raw",   ";#it{p}_{T}^{ll} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "zmass_wwctrl_raw", ";#it{m}_{ll} [GeV];Events", 100,20,300) );
    mon.addHistogram( new TH1F( "pfmet_wwctrl_raw", ";E_{T}^{miss} [GeV];Events", 50,0,300) );
    mon.addHistogram( new TH1F( "mt_wwctrl_raw",";#it{m}_{T}(#it{ll}, E_{T}^{miss}) [GeV];Events", 50,0,300) );


    //WZ control
    mon.addHistogram( new TH1F( "pfmet_WZctrl",         ";E_{T}^{miss} [GeV];Events", 20,0,500));
    mon.addHistogram( new TH1F( "balancedif_WZctrl",    ";|E_{T}^{miss}-#it{q}_{T}|/#it{q}_{T};Events", 20,0,1.0) );
    mon.addHistogram( new TH1F( "DphiZMET_WZctrl",      ";#Delta#it{#phi}(#it{l^{+}l^{-}},E_{T}^{miss});Events", 20,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "zpt_WZctrl",           ";#it{p}_{T}^{ll} [GeV];Events", 20,0,300) );

    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET40",    ";E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET40",       ";#it{m}_{T} [GeV];Events", 20,0,800) );
    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET50",    ";E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET50",       ";#it{m}_{T} [GeV];Events", 20,0,800) );
    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET60",    ";E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET60",       ";#it{m}_{T} [GeV];Events", 20,0,800) );
    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET70",    ";E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET70",       ";#it{m}_{T} [GeV];Events", 20,0,800) );
    mon.addHistogram( new TH1F( "pfmet_WZctrl_ZZlike_MET80",    ";E_{T}^{miss} [GeV];Events", 20,0,400));
    mon.addHistogram( new TH1F( "mt_WZctrl_ZZlike_MET80",       ";#it{m}_{T} [GeV];Events", 20,0,800) );






    //##################################################################################
    //########################## STUFF FOR CUTS OPTIMIZATION  ##########################
    //##################################################################################
    //optimization
    std::vector<double> optim_Cuts1_MET;
    std::vector<double> optim_Cuts1_Balance;
    std::vector<double> optim_Cuts1_DphiZMET;

    bool runOptimization = runProcess.getParameter<bool>("runOptimization");
    if(runOptimization) {
        // for optimization
        cout << "Optimization will be performed for this analysis" << endl;
        for(double met=80; met<=200; met+=10) {
            for(double balance=0.2; balance<=0.2; balance+=0.05) {
                for(double dphi=2.7; dphi<2.8; dphi+=0.1) {
                    optim_Cuts1_MET     .push_back(met);
                    optim_Cuts1_Balance .push_back(balance);
                    optim_Cuts1_DphiZMET.push_back(dphi);
                }
            }
        }
    }

    size_t nOptims = optim_Cuts1_MET.size();


    //make it as a TProfile so hadd does not change the value
    TProfile* Hoptim_cuts1_MET      = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_MET",";cut index;met",nOptims,0,nOptims) );
    TProfile* Hoptim_cuts1_Balance  = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_Balance",";cut index;balance",nOptims,0,nOptims) );
    TProfile* Hoptim_cuts1_DphiZMET = (TProfile*) mon.addHistogram( new TProfile ("optim_cut1_DphiZMET",";cut index;dphi",nOptims,0,nOptims) );

    for(size_t index=0; index<nOptims; index++) {
        Hoptim_cuts1_MET        ->Fill(index, optim_Cuts1_MET[index]);
        Hoptim_cuts1_Balance    ->Fill(index, optim_Cuts1_Balance[index]);
        Hoptim_cuts1_DphiZMET   ->Fill(index, optim_Cuts1_DphiZMET[index]);
    }

    TH1F* Hoptim_systs  =  (TH1F*) mon.addHistogram( new TH1F ("optim_systs"    , ";syst;", nvarsToInclude,0,nvarsToInclude) ) ;



    for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
        Hoptim_systs->GetXaxis()->SetBinLabel(ivar+1, varNames[ivar]);

        //1D shapes for limit setting
        mon.addHistogram( new TH2F (TString("mt_shapes")+varNames[ivar],";cut index; #it{m}_{T} [GeV];#Events (/100GeV)",nOptims,0,nOptims,12,0,1200) );
        //mon.addHistogram( new TH2F (TString("pfmet_shapes")+varNames[ivar],";cut index; E_{T}^{miss} [GeV];#Events",nOptims,0,nOptims,nBinPFMET,xbinsPFMET) );

        //2D shapes for limit setting
        //
        //
    }

    //##################################################################################
    //#############         GET READY FOR THE EVENT LOOP           #####################
    //##################################################################################

    //open the file and get events tree
    DataEvtSummaryHandler summaryHandler_;
    if(doWIMPreweighting) {
        if(url.Contains("TeV_DM_V_Mx")) url = "root://eoscms//eos/cms/store/group/phys_exotica/monoZ/llvvNtuple_25ns_29Oct2015/MC13TeV_DM_V_Mx1Mv200.root";
        if(url.Contains("TeV_DM_A_Mx")) url = "root://eoscms//eos/cms/store/group/phys_exotica/monoZ/llvvNtuple_25ns_29Oct2015/MC13TeV_DM_A_Mx1Mv200.root";

        if(url.Contains("K1_0.1_K2_1")) url.ReplaceAll("K1_0.1_K2_1","K1_1_K2_1");
        if(url.Contains("K1_0.2_K2_1")) url.ReplaceAll("K1_0.2_K2_1","K1_1_K2_1");
        if(url.Contains("K1_0.3_K2_1")) url.ReplaceAll("K1_0.3_K2_1","K1_1_K2_1");
        if(url.Contains("K1_0.5_K2_1")) url.ReplaceAll("K1_0.5_K2_1","K1_1_K2_1");
        if(url.Contains("K1_2_K2_1"))   url.ReplaceAll("K1_2_K2_1","K1_1_K2_1");
        if(url.Contains("K1_3_K2_1"))   url.ReplaceAll("K1_3_K2_1","K1_1_K2_1");
        if(url.Contains("K1_5_K2_1"))   url.ReplaceAll("K1_5_K2_1","K1_1_K2_1");
        if(url.Contains("K1_10_K2_1"))  url.ReplaceAll("K1_10_K2_1","K1_1_K2_1");

    }
    TFile *file = TFile::Open(url);
    printf("Looping on %s\n",url.Data());
    if(file==0) return -1;
    if(file->IsZombie()) return -1;
    if( !summaryHandler_.attachToTree( (TTree *)file->Get(dirname) ) ) {
        file->Close();
        return -1;
    }


    //check run range to compute scale factor (if not all entries are used)
    const Int_t totalEntries= summaryHandler_.getEntries();
    float rescaleFactor( evEnd>0 ?  float(totalEntries)/float(evEnd-evStart) : -1 );
    if(evEnd<0 || evEnd>summaryHandler_.getEntries() ) evEnd=totalEntries;
    if(evStart > evEnd ) {
        file->Close();
        return -1;
    }

    //MC normalization (to 1/pb)
    float cnorm=1.0;
    if(isMC) {
        //TH1F* cutflowH = (TH1F *) file->Get("mainAnalyzer/llvv/nevents");
        //if(cutflowH) cnorm=cutflowH->GetBinContent(1);
        TH1F* posH = (TH1F *) file->Get("mainAnalyzer/llvv/n_posevents");
        TH1F* negH = (TH1F *) file->Get("mainAnalyzer/llvv/n_negevents");
        if(posH && negH) cnorm = posH->GetBinContent(1) + negH->GetBinContent(1);
        if(rescaleFactor>0) cnorm /= rescaleFactor;
        printf("cnorm = %f\n",cnorm);
    }
    Hcutflow->SetBinContent(1,cnorm);



    //pileup weighting: based on vtx for now...
    std::vector<float> dataPileupDistribution;
    TString PU_Central = runProcess.getParameter<std::string>("PU_Central");
    gSystem->ExpandPathName(PU_Central);
    cout << "Loading PU weights: " << PU_Central << endl;
    TFile *PU_Central_File = TFile::Open(PU_Central);
    TH1F* histo = (TH1F *) PU_Central_File->Get("pileup");
    if(!histo)std::cout<<"pileup histogram is null!!!\n";
    for(int i=1; i<=histo->GetNbinsX(); i++) {
        dataPileupDistribution.push_back(histo->GetBinContent(i)/histo->Integral());
    }
    delete histo;
    PU_Central_File->Close();

    std::vector<float> mcPileupDistribution;
    bool useObservedPU(false);
    if(isMC) {
        TString puDist("mainAnalyzer/llvv/pileuptrue");
        if(useObservedPU) puDist="mainAnalyzer/llvv/pileup";
        TH1F* histo = (TH1F *) file->Get(puDist);
        if(!histo)std::cout<<"pileup histogram is null!!!\n";
        for(int i=1; i<=histo->GetNbinsX(); i++) {
            mcPileupDistribution.push_back(histo->GetBinContent(i)/histo->Integral());
        }
        delete histo;
        if(dataPileupDistribution.size()==0) dataPileupDistribution=mcPileupDistribution;
    }

    while(mcPileupDistribution.size()<dataPileupDistribution.size())  mcPileupDistribution.push_back(0.0);
    while(mcPileupDistribution.size()>dataPileupDistribution.size())  dataPileupDistribution.push_back(0.0);

    gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE
    edm::LumiReWeighting *LumiWeights=0;
    PuShifter_t PuShifters;
    if(isMC) {
        LumiWeights = new edm::LumiReWeighting(mcPileupDistribution,dataPileupDistribution);
        PuShifters=getPUshifters(dataPileupDistribution,0.05);
    }



    // event categorizer
    EventCategory eventCategoryInst(1);   //jet(0,1,>=2) binning


    //####################################################################################################################
    //###########################################           EVENT LOOP         ###########################################
    //####################################################################################################################

    // loop on all the events
    printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    printf("Scanning the ntuple :");
    int treeStep = (evEnd-evStart)/50;
    if(treeStep==0)treeStep=1;
    DuplicatesChecker duplicatesChecker;
    int nDuplicates(0);
    for( int iev=evStart; iev<evEnd; iev++) {
        if((iev-evStart)%treeStep==0) {
            printf("#");
            fflush(stdout);
        }

        //##############################################   EVENT LOOP STARTS   ##############################################
        //load the event content from tree
        summaryHandler_.getEntry(iev);
        DataEvtSummary_t &ev=summaryHandler_.getEvent();
        if(!isMC && duplicatesChecker.isDuplicate( ev.run, ev.lumi, ev.event) ) {
            nDuplicates++;
            cout << "nDuplicates: " << nDuplicates << endl;
            continue;
        }


        //prepare the tag's vectors for histo filling
        std::vector<TString> tags(1,"all");

        //genWeight
        float genWeight = 1.0;
        if(isMC && ev.genWeight<0) genWeight = -1.0;

        //pileup weight
        float weight = 1.0;
        if(isMC) weight *= genWeight;
        double TotalWeight_plus = 1.0;
        double TotalWeight_minus = 1.0;

        if(isMC) mon.fillHisto("pileup", "all", ev.ngenTruepu, weight);

        if(isMC) {
            weight            *= LumiWeights->weight(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
            TotalWeight_plus  *= PuShifters[PUUP]->Eval(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
            TotalWeight_minus *= PuShifters[PUDOWN]->Eval(useObservedPU ? ev.ngenITpu : ev.ngenTruepu);
        }

        Hcutflow->Fill(1,1);
        Hcutflow->Fill(2,weight);
        Hcutflow->Fill(3,weight*TotalWeight_minus);
        Hcutflow->Fill(4,weight*TotalWeight_plus);


        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);

        // FIXME need to have a function: loop all leptons, find a Z candidate,
        // can have input, ev.mn, ev.en
        // assign ee,mm,emu channel
        // check if channel name is consistent with trigger
        // store dilepton candidate in lep1 lep2, and other leptons in 3rdleps


        bool hasMMtrigger = ev.triggerType & 0x1;
        bool hasMtrigger  = (ev.triggerType >> 1 ) & 0x1;
        bool hasEEtrigger = (ev.triggerType >> 2 ) & 0x1;
        bool hasEtrigger  = (ev.triggerType >> 3 ) & 0x1;
        bool hasEMtrigger = (ev.triggerType >> 4 ) & 0x1;



        //#########################################################################
        //####################  Generator Level Reweighting  ######################
        //#########################################################################


        //for Wimps
        if(isMC_WIMP || isMC_ADD || isMC_Unpart || isMC_ZZ) {
            if(phys.genleptons.size()!=2) continue;
            if(phys.genGravitons.size()!=1 && phys.genWIMPs.size()!=2 && phys.genneutrinos.size()!=2) continue;

            LorentzVector genmet(0,0,0,0);
            if(phys.genWIMPs.size()==2) 	  genmet = phys.genWIMPs[0]+phys.genWIMPs[1];
            else if(phys.genGravitons.size()==1)  genmet = phys.genGravitons[0];
            else if(phys.genneutrinos.size()==2)  genmet = phys.genneutrinos[0]+phys.genneutrinos[1];

            LorentzVector dilep = phys.genleptons[0]+phys.genleptons[1];
            double dphizmet = fabs(deltaPhi(dilep.phi(),genmet.phi()));

            //reweighting
            if(doWIMPreweighting) {
                if(url.Contains("TeV_DM_V_Mx") || url.Contains("TeV_DM_A_Mx")) weight *= myWIMPweights.get1DWeights(genmet.pt(),"genmet");
                if(url.Contains("TeV_EWKDM_S_Mx")) weight *= myWIMPweights.get1DWeights(genmet.pt(),"pt_chichi");
            }
            //if(doWIMPreweighting) weight *= myWIMPweights.get2DWeights(genmet.pt(),dphizmet,"dphi_vs_met");

            mon.fillHisto("met_Gen", tags, genmet.pt(), weight, true);
            mon.fillHisto("met2_Gen", tags, genmet.pt(), weight);
            mon.fillHisto("zpt_Gen", tags, dilep.pt(), weight);
            mon.fillHisto("dphi_Gen", tags, dphizmet, weight);
            mon.fillHisto("zmass_Gen", tags, dilep.mass(), weight);
            if(phys.genleptons[0].pt() > phys.genleptons[1].pt()) mon.fillHisto("ptlep1vs2_Gen", tags, phys.genleptons[0].pt(), phys.genleptons[1].pt(), weight);
            else mon.fillHisto("ptlep1vs2_Gen", tags, phys.genleptons[1].pt(), phys.genleptons[0].pt(), weight);


            bool isInGenAcceptance = (phys.genleptons[0].pt()>20 && phys.genleptons[1].pt()>20);
            isInGenAcceptance &= ((fabs(dilep.mass()-91)<15) && (genmet.pt()>80) && (dilep.pt()>50) && (dphizmet>2.7));
            isInGenAcceptance &= (fabs(genmet.pt()-dilep.pt())/dilep.pt() < 0.2);

            mon.fillHisto("acceptance",tags,0,1.0);
            if(isInGenAcceptance) mon.fillHisto("acceptance",tags,1,1.0);

        }





        //#########################################################################
        //#####################      Objects Selection       ######################
        //#########################################################################

        //
        // MET ANALYSIS
        //
        //apply Jet Energy Resolution corrections to jets (and compute associated variations on the MET variable)
        std::vector<PhysicsObjectJetCollection> variedJets;
        LorentzVectorCollection variedMET;

        METUtils::computeVariation(phys.jets, phys.leptons, (usemetNoHF ? phys.metNoHF : phys.met), variedJets, variedMET, &jecUnc);

        LorentzVector metP4=variedMET[0];
        PhysicsObjectJetCollection &corrJets = variedJets[0];

        //
        // LEPTON ANALYSIS
        //

        // looping leptons (electrons + muons)
        int nGoodLeptons(0);
        std::vector<std::pair<int,LorentzVector> > goodLeptons;
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            LorentzVector lep=phys.leptons[ilep];
            int lepid = phys.leptons[ilep].id;
            if(lep.pt()<20) continue;
            if(abs(lepid)==13 && fabs(lep.eta())> 2.4) continue;
            if(abs(lepid)==11 && fabs(lep.eta())> 2.5) continue;

            bool hasTightIdandIso(true);
            if(abs(lepid)==13) { //muon
                hasTightIdandIso &= phys.leptons[ilep].isTightMu;
                //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2?sortcol=1;table=7;up=0#Muon_Isolation
                hasTightIdandIso &= ( phys.leptons[ilep].m_pfRelIsoDbeta() < 0.15 );
            } else if(abs(lepid)==11) { //electron
                hasTightIdandIso &= phys.leptons[ilep].isElpassMedium;
            } else continue;


            if(!hasTightIdandIso) continue;
            nGoodLeptons++;
            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(lepid,lep);
            goodLeptons.push_back(goodlep);

        }

        if(nGoodLeptons<2) continue; // 2 tight leptons

        float _MASSDIF_(999.);
        int id1(0),id2(0);
        LorentzVector lep1(0,0,0,0),lep2(0,0,0,0);
        for(size_t ilep=0; ilep<goodLeptons.size(); ilep++) {
            int id1_ = goodLeptons[ilep].first;
            LorentzVector lep1_ = goodLeptons[ilep].second;

            for(size_t jlep=ilep+1; jlep<goodLeptons.size(); jlep++) {
                int id2_ = goodLeptons[jlep].first;
                LorentzVector lep2_ = goodLeptons[jlep].second;
                if(id1_*id2_>0) continue; // opposite charge

                LorentzVector dilepton=lep1_+lep2_;
                double massdif = fabs(dilepton.mass()-91.);
                if(massdif < _MASSDIF_) {
                    lep1.SetPxPyPzE(lep1_.px(),lep1_.py(),lep1_.pz(),lep1_.energy());
                    lep2.SetPxPyPzE(lep2_.px(),lep2_.py(),lep2_.pz(),lep2_.energy());
                    id1 = id1_;
                    id2 = id2_;
                }
            }
        }



        if(id1*id2==0) continue;
        LorentzVector zll(lep1+lep2);
        bool passZmass(fabs(zll.mass()-91)<10);
        bool passZpt(zll.pt()>50);


        TString tag_cat;
        int evcat = getDileptonId(abs(id1),abs(id2));
        switch(evcat) {
        case MUMU :
            tag_cat = "mumu";
            break;
        case EE   :
            tag_cat = "ee";
            break;
        case EMU  :
            tag_cat = "emu";
            break;
        default   :
            continue;
        }


        bool hasTrigger(false);

        if(!isMC) {
            if(evcat!=fType) continue;

            if(evcat==EE   && !(hasEEtrigger||hasEtrigger) ) continue;
            if(evcat==MUMU && !(hasMMtrigger||hasMtrigger) ) continue;
            if(evcat==EMU  && !hasEMtrigger && !(hasEtrigger && hasMtrigger) ) continue;

            //this is a safety veto for the single mu PD
            if(isSingleMuPD) {
                if(!hasMtrigger) continue;
                if(hasMtrigger && hasMMtrigger) continue;
            }
            if(isDoubleMuPD) {
                if(!hasMMtrigger) continue;

            }

            //this is a safety veto for the single Ele PD
            if(isSingleElePD) {
                if(!hasEtrigger) continue;
                if(hasEtrigger && hasEEtrigger) continue;
            }
            if(isDoubleElePD) {
                if(!hasEEtrigger) continue;
            }

            hasTrigger=true;

        } else {
            if(evcat==EE   && (hasEEtrigger || hasEtrigger) ) hasTrigger=true;
            if(evcat==MUMU && (hasMMtrigger || hasMtrigger) ) hasTrigger=true;
            if(evcat==EMU  && (hasEMtrigger || (hasEtrigger && hasMtrigger)) ) hasTrigger=true;
            if(!hasTrigger) continue;
        }

        tags.push_back(tag_cat); //add ee, mumu, emu category


        // pielup reweightiing
        mon.fillHisto("nvtx_raw",   tags, phys.nvtx,      1.0);
        //if(isMC) weight *= myWIMPweights.get1DWeights(phys.nvtx,"pileup_weights");
        mon.fillHisto("nvtxwgt_raw",tags, phys.nvtx,      weight);


        mon.fillHisto("eventflow",tags,0,weight);
        mon.fillHisto("nleptons_raw",tags, nGoodLeptons, weight);


        if(lep1.pt()>lep2.pt()) {
            mon.fillHisto("leadlep_pt_raw",   tags, lep1.pt(), weight);
            mon.fillHisto("leadlep_eta_raw",  tags, lep1.eta(), weight);
            mon.fillHisto("trailep_pt_raw",   tags, lep2.pt(), weight);
            mon.fillHisto("trailep_eta_raw",  tags, lep2.eta(), weight);
        } else {
            mon.fillHisto("leadlep_pt_raw",   tags, lep2.pt(), weight);
            mon.fillHisto("leadlep_eta_raw",  tags, lep2.eta(), weight);
            mon.fillHisto("trailep_pt_raw",   tags, lep1.pt(), weight);
            mon.fillHisto("trailep_eta_raw",  tags, lep1.eta(), weight);
        }


        //
        // 3rd LEPTON ANALYSIS
        //

        //loop over all lepton again, check deltaR with dilepton,
        bool pass3dLeptonVeto(true);
        bool hasTight3dLepton(false);
        int n3rdLeptons(0), nTight3rdLeptons(0);
        vector<LorentzVector> allLeptons;
        allLeptons.push_back(lep1);
        allLeptons.push_back(lep2);
        std::vector<LorentzVector> extraTight10Leptons;
        for(size_t ilep=0; ilep<phys.leptons.size(); ilep++) {
            LorentzVector lep=phys.leptons[ilep];
            int lepid = phys.leptons[ilep].id;
            if(abs(lepid)==13 && fabs(lep.eta())> 2.4) continue;
            if(abs(lepid)==11 && fabs(lep.eta())> 2.5) continue;

            bool isMatched(false);
            isMatched |= (deltaR(lep1,lep) < 0.01);
            isMatched |= (deltaR(lep2,lep) < 0.01);
            if(isMatched) continue;

            bool hasLooseIdandIso(true);
            bool hasTightIdandIso(true);
            if(abs(lepid)==13) { //muon
                hasLooseIdandIso &= ( phys.leptons[ilep].isLooseMu && phys.leptons[ilep].m_pfRelIsoDbeta()<0.25 && phys.leptons[ilep].pt()>10 );
                hasLooseIdandIso |= ( phys.leptons[ilep].isSoftMu  && phys.leptons[ilep].pt()>3 );
                //
                hasTightIdandIso &= ( phys.leptons[ilep].isTightMu && phys.leptons[ilep].m_pfRelIsoDbeta()<0.15 && phys.leptons[ilep].pt()>10 );
            } else if(abs(lepid)==11) { //electron
                hasLooseIdandIso &= ( phys.leptons[ilep].isElpassVeto && phys.leptons[ilep].pt()>10 );
                //
                hasTightIdandIso &= ( phys.leptons[ilep].isElpassMedium && phys.leptons[ilep].pt()>10 );

            } else continue;


            if(!hasLooseIdandIso) continue;
            allLeptons.push_back(lep);
            n3rdLeptons++;


            if(hasTightIdandIso) {
                extraTight10Leptons.push_back(lep);
                nTight3rdLeptons++;
            }
        }

        pass3dLeptonVeto=(n3rdLeptons==0);
        hasTight3dLepton=(nTight3rdLeptons==1);




        //
        //JET AND BTAGING ANALYSIS
        //
        PhysicsObjectJetCollection GoodIdJets;
        bool passBveto(true);
        int nJetsGood30(0);
        int nCSVLtags(0),nCSVMtags(0),nCSVTtags(0);
        double BTagScaleFactor(1.0);
        for(size_t ijet=0; ijet<corrJets.size(); ijet++) {

            if(corrJets[ijet].pt()<20) continue;
            if(fabs(corrJets[ijet].eta())>5.0) continue;

            //jet ID
            if(!corrJets[ijet].isPFLoose) continue;
            //if(corrJets[ijet].pumva<0.5) continue;

            //check overlaps with selected leptons
            double minDR(999.);
            for(vector<LorentzVector>::iterator lIt = allLeptons.begin(); lIt != allLeptons.end(); lIt++) {
                double dR = deltaR( corrJets[ijet], *lIt );
                if(dR > minDR) continue;
                minDR = dR;
            }
            if(minDR < 0.4) continue;

            GoodIdJets.push_back(corrJets[ijet]);
            if(corrJets[ijet].pt()>30) nJetsGood30++;


            //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation74X
            if(corrJets[ijet].pt()>30 && fabs(corrJets[ijet].eta())<2.4)  {

                nCSVLtags += (corrJets[ijet].btag0>0.605);
                nCSVMtags += (corrJets[ijet].btag0>0.89);
                nCSVTtags += (corrJets[ijet].btag0>0.97);

                bool isCSVLtagged(corrJets[ijet].btag0>0.605);
                if(abs(corrJets[ijet].flavid)==5) {
                    //BTagScaleFactor *= myBtagUtils.getBTagWeight(isCSVLtagged,corrJets[ijet].pt(),corrJets[ijet].eta(),abs(corrJets[ijet].flavid),"CSVL","CSVL/b_eff").first;
                    if(isCSVLtagged) BTagScaleFactor *= btag_reader.eval( BTagEntry::FLAV_B, corrJets[ijet].eta(), (corrJets[ijet].pt()<670. ? corrJets[ijet].pt() : 670.) );
                } else if(abs(corrJets[ijet].flavid)==4) {
                    //BTagScaleFactor *= myBtagUtils.getBTagWeight(isCSVLtagged,corrJets[ijet].pt(),corrJets[ijet].eta(),abs(corrJets[ijet].flavid),"CSVL","CSVL/c_eff").first;
                    if(isCSVLtagged) BTagScaleFactor *= btag_reader.eval( BTagEntry::FLAV_C, corrJets[ijet].eta(), (corrJets[ijet].pt()<670. ? corrJets[ijet].pt() : 670.) );
                } else {
                    //BTagScaleFactor *= myBtagUtils.getBTagWeight(isCSVLtagged,corrJets[ijet].pt(),corrJets[ijet].eta(),abs(corrJets[ijet].flavid),"CSVL","CSVL/udsg_eff").first;
                    //no SF available yet for light flavor jet (CSVv2.csv)
                    //BTagScaleFactor *= btag_reader.eval( BTagEntry::FLAV_UDSG, corrJets[ijet].eta(), (corrJets[ijet].pt()<670. ? corrJets[ijet].pt() : 670.) );
                }


                // begin btagging efficiency
                if(isMC_ttbar||isMC_stop) {
                    int flavid = abs(corrJets[ijet].flavid);
                    for(size_t csvtag=0; csvtag<CSVkey.size(); csvtag++) {
                        bool isBTag(false);
                        if     (CSVkey[csvtag]=="CSVL" && (corrJets[ijet].btag0>0.605)) isBTag = true;
                        else if(CSVkey[csvtag]=="CSVM" && (corrJets[ijet].btag0>0.89)) isBTag = true;
                        else if(CSVkey[csvtag]=="CSVT" && (corrJets[ijet].btag0>0.97)) isBTag = true;

                        if(flavid==5) {
                            mon.fillHisto(TString("beff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                            if(isBTag) mon.fillHisto(TString("beff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        } else if(flavid==4) {
                            mon.fillHisto(TString("ceff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                            if(isBTag) mon.fillHisto(TString("ceff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        } else {
                            mon.fillHisto(TString("udsgeff_Denom_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                            if(isBTag) mon.fillHisto(TString("udsgeff_Num_")+CSVkey[csvtag],tags,corrJets[ijet].pt(),weight);
                        }

                    }
                } //end btagging efficiency


            }


        }

        //using CSV Loose WP
        passBveto=(nCSVLtags==0);

        for(size_t ij=0; ij<GoodIdJets.size(); ij++) {
            mon.fillHisto("jet_pt_raw",   tags, GoodIdJets[ij].pt(),weight);
            mon.fillHisto("jet_eta_raw",  tags, GoodIdJets[ij].eta(),weight);
        }

        double dphiZMET=fabs(deltaPhi(zll.phi(),metP4.phi()));
        bool passDphiZMETcut(dphiZMET>2.7);

        TVector2 dil2(zll.px(),zll.py());
        TVector2 met2(metP4.px(),metP4.py());
        double axialmet=dil2*met2;
        axialmet /= -zll.pt();

        //missing ET
        bool passMETcut=(metP4.pt()>80);
        bool passMETcut120=(metP4.pt()>120);


        double response = METUtils::response(zll,metP4);
        bool passResponseCut = (response>-1 && response<1);

        //missing ET balance
        bool passBalanceCut=(metP4.pt()/zll.pt()>0.80 && metP4.pt()/zll.pt()<1.20);
        double balanceDif = fabs(1-metP4.pt()/zll.pt());

        //transverse mass
        double MT_massless = METUtils::transverseMass(zll,metP4,false);

        //#########################################################
        //####  RUN PRESELECTION AND CONTROL REGION PLOTS  ########
        //#########################################################

        //event category
        int eventSubCat  = eventCategoryInst.Get(phys,&GoodIdJets);
        TString tag_subcat = eventCategoryInst.GetLabel(eventSubCat);

        tags.push_back(tag_cat+tag_subcat); // add jet binning category
        if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") tags.push_back(tag_cat+"lesq1jets");
        if(tag_cat=="mumu" || tag_cat=="ee") {
            tags.push_back("ll"+tag_subcat);
            if(tag_subcat=="eq0jets" || tag_subcat=="eq1jets") tags.push_back("lllesq1jets");
        }


        //apply weights
        if(isMC) weight *= BTagScaleFactor;


        mon.fillHisto("zpt_raw"                         ,tags, zll.pt(),   weight);
        mon.fillHisto("pfmet_raw"                       ,tags, metP4.pt(), weight);

        mon.fillHisto("mt_raw"                          ,tags, MT_massless, weight);
        mon.fillHisto("mt2_raw"                          ,tags, MT_massless, weight,true);
        mon.fillHisto("zmass_raw"                       ,tags, zll.mass(), weight);
        mon.fillHisto("njets_raw"                       ,tags, nJetsGood30, weight);
        mon.fillHisto("nbjets_raw"                      ,tags, nCSVLtags, weight);
        if(lep1.pt()>lep2.pt()) mon.fillHisto("ptlep1vs2_raw"                   ,tags, lep1.pt(), lep2.pt(), weight);
        else 			mon.fillHisto("ptlep1vs2_raw"                   ,tags, lep2.pt(), lep1.pt(), weight);


        // WW/ttbar/Wt/tautau control
        mon.fillHisto("zpt_wwctrl_raw"                      ,tags, zll.pt(),   weight);
        mon.fillHisto("zmass_wwctrl_raw"                    ,tags, zll.mass(), weight);
        mon.fillHisto("pfmet_wwctrl_raw"                    ,tags, metP4.pt(), weight);
        mon.fillHisto("mt_wwctrl_raw"              	    ,tags, MT_massless, weight);

        //
        // WZ control region
        //
        if(hasTight3dLepton && passZmass && passBveto) {
            mon.fillHisto("pfmet_WZctrl"                ,tags, metP4.pt(), weight);
            mon.fillHisto("balancedif_WZctrl"           ,tags, balanceDif, weight);
            mon.fillHisto("DphiZMET_WZctrl"             ,tags, dphiZMET, weight);
            mon.fillHisto("zpt_WZctrl"                  ,tags, zll.pt(),   weight);

            if(metP4.pt()>40) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET40"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET40"     ,tags, fakeMT, weight);
            }

            if(metP4.pt()>50) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET50"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET50"     ,tags, fakeMT, weight);
            }

            if(metP4.pt()>60) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET60"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET60"     ,tags, fakeMT, weight);
            }

            if(metP4.pt()>70) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET70"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET70"     ,tags, fakeMT, weight);
            }

            if(metP4.pt()>80) {
                LorentzVector fakeMET = extraTight10Leptons[0]+metP4;
                double fakeMT=METUtils::transverseMass(zll,fakeMET,false);
                mon.fillHisto("pfmet_WZctrl_ZZlike_MET80"     ,tags, fakeMET.pt(), weight);
                mon.fillHisto("mt_WZctrl_ZZlike_MET80"     ,tags, fakeMT, weight);
            }


        }

        //##############################################
        //########  Main Event Selection        ########
        //##############################################

        //for MET X-Y shift correction
        mon.fillHisto("pfmetx_vs_nvtx_presel",tags,phys.nvtx,metP4.px(), weight);
        mon.fillHisto("pfmety_vs_nvtx_presel",tags,phys.nvtx,metP4.py(), weight);
        LorentzVector metP4_XYCorr = METUtils::applyMETXYCorr(metP4,isMC,phys.nvtx);
        mon.fillHisto("pfmetphi_wocorr_presel",tags, metP4.phi(), weight);
        mon.fillHisto("pfmetphi_wicorr_presel",tags, metP4_XYCorr.phi(), weight);

        mon.fillHisto("pfmet_wicorr_presel",tags, metP4_XYCorr.pt(), weight, true);
        mon.fillHisto("pfmet2_wicorr_presel",tags, metP4_XYCorr.pt(), weight, true);


        if(passZmass) {
            mon.fillHisto("eventflow",  tags, 1, weight);

            if(passZpt) {
                mon.fillHisto("eventflow",  tags, 2, weight);

                if(pass3dLeptonVeto) {
                    mon.fillHisto("eventflow",  tags, 3, weight);

                    if(passBveto) {
                        mon.fillHisto("eventflow",  tags, 4, weight);

                        //preselection plots
                        if(!isMC) {
                            //data blinded
                            if(metP4.pt()<80) {
                                mon.fillHisto("pfmet_presel",tags, metP4.pt(), weight, true);
                                mon.fillHisto("pfmet2_presel",tags, metP4.pt(), weight, true);
                                mon.fillHisto("mt_presel",   tags, MT_massless, weight);
                                mon.fillHisto("mt2_presel",   tags, MT_massless, weight, true);
                            }
                        } else {
                            mon.fillHisto("pfmet_presel",tags, metP4.pt(), weight, true);
                            mon.fillHisto("pfmet2_presel",tags, metP4.pt(), weight, true);
                            mon.fillHisto("mt_presel",   tags, MT_massless, weight);
                            mon.fillHisto("mt2_presel",   tags, MT_massless, weight, true);
                        }

                        mon.fillHisto("dphiZMET_presel",tags, dphiZMET, weight);
                        mon.fillHisto("balancedif_presel",tags, balanceDif, weight);

                        mon.fillHisto("axialpfmet_presel",    tags,  axialmet,  weight);


                        if(passDphiZMETcut && passResponseCut) {
                            mon.fillHisto("eventflow",  tags, 5, weight);

                            if(passBalanceCut) {
                                mon.fillHisto("eventflow",  tags, 6, weight);

                                if(passMETcut) {
                                    mon.fillHisto("eventflow",  tags, 7, weight);
                                    mon.fillHisto("mt_final",   tags, MT_massless, weight);
                                    mon.fillHisto("pfmet_final",tags, metP4.pt(), weight, true);
                                    mon.fillHisto("pfmet2_final",tags, metP4.pt(), weight, true);

                                    if(passMETcut120) mon.fillHisto("mt_final120",   tags, MT_massless, weight);

                                    if(!isMC) fprintf(outTxtFile_final,"%d | %d | %d | pfmet: %f | mt: %f | mass: %f \n",ev.run,ev.lumi,ev.event,metP4.pt(), MT_massless,zll.mass());

                                } //passMETcut

                            } //passBalanceCut

                        } //passDphiZMETcut

                    } //passBveto

                } //pass3dLeptonVeto

            } //passZpt

        } //passZmass











        //##############################################################################
        //### HISTOS FOR STATISTICAL ANALYSIS (include systematic variations)
        //##############################################################################


        //Fill histogram for posterior optimization, or for control regions
        for(size_t ivar=0; ivar<nvarsToInclude; ivar++) {
            float iweight = weight;                                               //nominal
            if(varNames[ivar]=="_puup")        iweight *=TotalWeight_plus;        //pu up
            if(varNames[ivar]=="_pudown")      iweight *=TotalWeight_minus;       //pu down

            if(isSignal && (varNames[ivar]=="_pdfup" || varNames[ivar]=="_pdfdown")) {
                if(mPDFInfo) {
                    float PDFWeight_plus(1.0), PDFWeight_down(1.0);
                    std::vector<float> wgts=mPDFInfo->getWeights(iev);
                    for(size_t ipw=0; ipw<wgts.size(); ipw++) {
                        PDFWeight_plus = TMath::Max(PDFWeight_plus,wgts[ipw]);
                        PDFWeight_down = TMath::Min(PDFWeight_down,wgts[ipw]);
                    }
                    if(varNames[ivar]=="_pdfup")    iweight *= PDFWeight_plus;
                    if(varNames[ivar]=="_pdfdown")  iweight *= PDFWeight_down;
                }
            }

            //##############################################
            // recompute MET/MT if JES/JER was varied
            //##############################################
            //LorentzVector vMET = variedMET[ivar>8 ? 0 : ivar];
            //PhysicsObjectJetCollection &vJets = ( ivar<=4 ? variedJets[ivar] : variedJets[0] );


            LorentzVector vMET = variedMET[0];
            if(varNames[ivar]=="_jerup" || varNames[ivar]=="_jerdown" || varNames[ivar]=="_jesup" || varNames[ivar]=="_jesdown" ||
                    varNames[ivar]=="_umetup" || varNames[ivar]=="_umetdown" || varNames[ivar]=="_lesup" || varNames[ivar]=="_lesdown") {
                vMET = variedMET[ivar];
            }

            PhysicsObjectJetCollection &vJets = variedJets[0];
            if(varNames[ivar]=="_jerup" || varNames[ivar]=="_jerdown" || varNames[ivar]=="_jesup" || varNames[ivar]=="_jesdown") {
                vJets = variedJets[ivar];
            }


            bool passLocalBveto(true);
            for(size_t ijet=0; ijet<vJets.size(); ijet++) {

                if(vJets[ijet].pt()<20) continue;
                if(fabs(vJets[ijet].eta())>5.0) continue;

                //jet ID
                if(!vJets[ijet].isPFLoose) continue;
                //if(vJets[ijet].pumva<0.5) continue;

                //check overlaps with selected leptons
                double minDR(999.);
                for(vector<LorentzVector>::iterator lIt = allLeptons.begin(); lIt != allLeptons.end(); lIt++) {
                    double dR = deltaR( vJets[ijet], *lIt );
                    if(dR > minDR) continue;
                    minDR = dR;
                }
                if(minDR < 0.4) continue;


                if(vJets[ijet].pt()>30 && fabs(vJets[ijet].eta())<2.4) {

                    passLocalBveto &= (vJets[ijet].btag0<0.605);
                    bool isLocalCSVLtagged(vJets[ijet].btag0>0.605);

                    //double val=1., valerr=0.;
		    double BTagWeights_Up=1., BTagWeights_Down=1.;
                    if(abs(vJets[ijet].flavid)==5) {
                        //val = myBtagUtils.getBTagWeight(isLocalCSVLtagged,vJets[ijet].pt(),vJets[ijet].eta(),abs(vJets[ijet].flavid),"CSVL","CSVL/b_eff").first;
                        //valerr = myBtagUtils.getBTagWeight(isLocalCSVLtagged,vJets[ijet].pt(),vJets[ijet].eta(),abs(vJets[ijet].flavid),"CSVL","CSVL/b_eff").second;
			if(isLocalCSVLtagged) BTagWeights_Up *= btag_reader_up.eval( BTagEntry::FLAV_B, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
			if(isLocalCSVLtagged) BTagWeights_Down *= btag_reader_down.eval( BTagEntry::FLAV_B, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );

                    } else if(abs(vJets[ijet].flavid)==4) {
                        //val = myBtagUtils.getBTagWeight(isLocalCSVLtagged,vJets[ijet].pt(),vJets[ijet].eta(),abs(vJets[ijet].flavid),"CSVL","CSVL/c_eff").first;
                        //valerr = myBtagUtils.getBTagWeight(isLocalCSVLtagged,vJets[ijet].pt(),vJets[ijet].eta(),abs(vJets[ijet].flavid),"CSVL","CSVL/c_eff").second;
                        if(isLocalCSVLtagged) BTagWeights_Up *= btag_reader_up.eval( BTagEntry::FLAV_C, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );
                        if(isLocalCSVLtagged) BTagWeights_Down *= btag_reader_down.eval( BTagEntry::FLAV_C, vJets[ijet].eta(), (vJets[ijet].pt()<670. ? vJets[ijet].pt() : 670.) );

                    } else {
                        //val = myBtagUtils.getBTagWeight(isLocalCSVLtagged,vJets[ijet].pt(),vJets[ijet].eta(),abs(vJets[ijet].flavid),"CSVL","CSVL/udsg_eff").first;
                        //valerr= myBtagUtils.getBTagWeight(isLocalCSVLtagged,vJets[ijet].pt(),vJets[ijet].eta(),abs(vJets[ijet].flavid),"CSVL","CSVL/udsg_eff").second;
                    }
                    //double BTagWeights_Up = (val+valerr)/val;
                    //double BTagWeights_Down = (val-valerr)/val;
                    if(varNames[ivar]=="_btagup") 	iweight *= BTagWeights_Up;
                    if(varNames[ivar]=="_btagdown")	iweight *= BTagWeights_Down;
                }

            }


            double responseVal = METUtils::response(zll,vMET);
            bool passResponseValCut = (fabs(responseVal)<1);

            bool passBaseSelection( passZmass && passZpt && pass3dLeptonVeto && passLocalBveto && passResponseValCut);

            double mt_massless = METUtils::transverseMass(zll,vMET,false); //massless mt
            double LocalDphiZMET=fabs(deltaPhi(zll.phi(),vMET.phi()));

            //############
            //optimization
            //############
            for(unsigned int index=0; index<nOptims; index++) {


                double minMET = optim_Cuts1_MET[index];
                double minBalance = optim_Cuts1_Balance[index];
                double minDphi = optim_Cuts1_DphiZMET[index];

                bool passLocalMETcut(vMET.pt()>minMET);
                bool passLocalBalanceCut=(vMET.pt()/zll.pt()>(1.-minBalance) && vMET.pt()/zll.pt()<(1.+minBalance));
                bool passLocalDphiZMETcut(LocalDphiZMET>minDphi);

                bool passOptimSelection(passBaseSelection && passLocalMETcut && passLocalBalanceCut && passLocalDphiZMETcut);



                // fill shapes for limit setting
                if( passOptimSelection ) {
                    mon.fillHisto(TString("mt_shapes")+varNames[ivar],tags,index, mt_massless, iweight);
                }


            }//all optimization END






        }//Systematic variation END

















    } // loop on all events END




    printf("\n");
    file->Close();

    //##############################################
    //########     SAVING HISTO TO FILE     ########
    //##############################################
    //save control plots to file
    outUrl += "/";
    outUrl += outFileUrl + ".root";
    printf("Results saved in %s\n", outUrl.Data());

    //save all to the file
    TFile *ofile=TFile::Open(outUrl, "recreate");
    mon.Write();

    ofile->Close();

    if(outTxtFile_final)fclose(outTxtFile_final);
}
