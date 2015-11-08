#include <iostream>
#include <boost/shared_ptr.hpp>
#include "Math/GenVector/Boost.h"

#include <sstream>
#include "llvvAnalysis/DMAnalysis/src/tdrstyle.C"
#include "llvvAnalysis/DMAnalysis/src/JSONWrapper.cc"
#include "llvvAnalysis/DMAnalysis/interface/setStyle.h"
#include "llvvAnalysis/DMAnalysis/interface/MacroUtils.h"
#include "llvvAnalysis/DMAnalysis/interface/plotter.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"

#include<iostream>
#include<fstream>
#include<map>
#include<algorithm>
#include<vector>
#include<set>

using namespace std;
double NonResonnantSyst = 0.25; //0.1;//0.25;

double WWtopSyst_ee0jet = 0.;
double WWtopSyst_ee1jet = 0.;
double WWtopSyst_eelesq1jet = 0.;
double WWtopSyst_mm0jet = 0.;
double WWtopSyst_mm1jet = 0.;
double WWtopSyst_mmlesq1jet = 0.;

double GammaJetSyst = 0.7;//1.0; //0.5;//0.5, 1.0;
double ZjetsExtropSyst = 0.6;
double WjetsSyst_ee = 0.146;
double WjetsSyst_mm = 0.229;

struct YIELDS_T {
    double WZ;
    double ZZ;
    double Zjets;
    double WWtop;
    double Wjets;
    double Data;
    double Sig;
    double totBkg;

    double WZ_StatErr;
    double ZZ_StatErr;
    double Zjets_StatErr;
    double WWtop_StatErr;
    double Wjets_StatErr;
    double Data_StatErr;
    double Sig_StatErr;
    double totBkg_StatErr;

    YIELDS_T():WZ(0.),ZZ(0.),Zjets(0.),WWtop(0.),Wjets(0.),Data(0.),Sig(0.),totBkg(0.),
        WZ_StatErr(0.),ZZ_StatErr(0.),Zjets_StatErr(0.),WWtop_StatErr(0.),
        Wjets_StatErr(0.),Data_StatErr(0.),Sig_StatErr(0.),totBkg_StatErr(0.) { }
};




TFile *fOut_binbybinNuisance   = TFile::Open("BinbyBinNuisance.root","RECREATE");
TH1F  *nuisance_h = new TH1F("nuisance_h","nuisance_h",100,1,1.1);

// Map for non-shape systematics (lnN)
std::map<TString, double> normSysts;
void initNormalizationSysts();

//wrapper for a projected shape for a given set of cuts
class Shape_t {
public:
    TH1* data, *totalBckg;
    std::vector<TH1 *> bckg, signal, bckgInSignal;
    //the key corresponds to the proc name
    //the key is the name of the variation: e.g. jesup, jesdown, etc.
    std::map<TString,std::vector<std::pair<TString, TH1*> > > bckgVars, signalVars, bckgInSignalVars;

    std::map<TString, double> xsections;
    std::map<TString, double> BRs;

    Shape_t() {}
    ~Shape_t() {}


    void clear() {
        std::cout<<"shape is destructed...";
        if(data)delete data;
        if(totalBckg)totalBckg->Delete();
        for(unsigned int i=0; i<bckg.  size(); i++) {
            delete bckg  [i];
        }
        bckg  .clear();
        for(unsigned int i=0; i<signal.size(); i++) {
            delete signal[i];
        }
        signal.clear();
        for(std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it=bckgVars  .begin(); it!=bckgVars  .end(); it++) {
            for(unsigned int i=0; i<(*it).second.size(); i++) {
                delete (*it).second[i].second;
            }
        }
        bckgVars  .clear();
        for(std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it=signalVars.begin(); it!=signalVars.end(); it++) {
            for(unsigned int i=0; i<(*it).second.size(); i++) {
                delete (*it).second[i].second;
            }
        }
        signalVars.clear();
        std::cout<<"done\n";

    }


};

typedef std::pair<TString,TString> RateKey_t;
struct DataCardInputs {
    TString shapesFile;
    std::vector<TString> ch;
    std::vector<TString> procs;
    std::map<RateKey_t, Double_t> obs;
    std::map<RateKey_t, Double_t> rates;
    std::map<TString, std::map<RateKey_t,Double_t> > systs;
    int nsignalproc;
};


void printHelp();
Shape_t getShapeFromFile(TFile* inF, TString ch, TString shapeName, int cutBin,JSONWrapper::Object &Root,double minCut=0, double maxCut=9999, bool onlyData=false);
void showShape(std::vector<TString>& selCh ,map<TString, Shape_t>& allShapes, TString mainHisto, TString SaveName);
void Draw1DHistogram(TH1* mc, THStack *stack, TH1 *mcPlusRelUnc, TGraphErrors *errors, std::vector<TH1 *>& spimpose, TH1 *data, TLegend *legA, bool noLog, TString finstate, TString AnalysisBins);
void Draw2DHistogram(std::map<TString, TH1*>& mapbkgs, std::vector<TH1 *>& spimpose, TH1 *data, TString finstate, TString AnalysisBins);

void getYieldsFromShape(std::vector<TString> ch, const map<TString, Shape_t> &allShapes, TString shName, bool isdataBlinded);
void getEffFromShape(std::vector<TString> ch, const map<TString, Shape_t> &allShapes, TString shName);




void convertHistosForLimits_core(DataCardInputs& dci, TString& proc, TString& bin, TString& ch, std::vector<TString>& systs, std::vector<TH1*>& hshapes);
DataCardInputs convertHistosForLimits(Int_t mass,TString histo="finalmt",TString url="plotter.root",TString Json="");
DataCardInputs convertHistosForLimits(TString atgcpar,Int_t mass,TString histo="finalmt",TString url="plotter.root",TString Json="");
std::vector<TString> buildDataCard(Int_t mass, TString histo="finalmt", TString url="plotter.root",TString Json="");
std::vector<TString> buildDataCard(TString atgcpar,Int_t mass, TString histo="finalmt", TString url="plotter.root",TString Json="");
void doBackgroundSubtraction(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t> &allShapes, TString mainHisto, TString sideBandHisto, TString url, JSONWrapper::Object &Root, bool isMCclosureTest);
void dodataDrivenWWtW(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t> &allShapes, TString mainHisto,bool isMCclosureTest);
void doWjetsBackground(std::vector<TString>& selCh,map<TString, Shape_t> &allShapes, TString mainHisto);
void doDYextrapolation(std::vector<TString>& selCh,map<TString, Shape_t> &allShapes,TString mainHisto,TString DY_EXTRAPOL_SHAPES,float DY_EXTRAPOL,bool isleftCR);
void doQCDBackground(std::vector<TString>& selCh,map<TString, Shape_t> &allShapes, TString mainHisto);
void doDYReplacement(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString metHistoForRescale);
void doWZSubtraction(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto);
void BlindData(std::vector<TString>& selCh, map<TString, Shape_t>& allShapes, TString mainHisto, bool addSignal);

void RescaleForInterference(std::vector<TString>& selCh,map<TString, Shape_t>& allShapes, TString mainHisto);


bool subNRB2011 = false;
bool subNRB2012 = false;
bool MCclosureTest = false;

bool mergeWWandZZ = false;
bool skipWW = true;
std::vector<TString> Channels;
std::vector<TString> AnalysisBins;
bool fast = false;
bool skipGGH = false;
bool skipQQH = false;
bool subDY = false;
bool subWZ = false;
double DDRescale = 1.0;
double MCRescale = 1.0;
bool blindData = false;
bool blindWithSignal = false;
TString DYFile ="";
TString inFileUrl(""),jsonFile(""), histo("");
TString postfix="";
TString systpostfix="";
double shapeMin = 0;
double shapeMax = 9999;
double shapeMinVBF = 0;
double shapeMaxVBF = 9999;
bool doInterf = false;

float DYMET_EXTRAPOL = 9999;
float DYRESP_EXTRAPOL = 9999;
float DYDPHI_EXTRAPOL = 9999;

int indexvbf = -1;
int indexcut   = -1, indexcutL=-1, indexcutR=-1;
int mass=-1, massL=-1, massR=-1;
int MV=-1, MA=-1;
TString atgcpar="", atgcpar2="";
bool runSystematics = false;
bool shape = false;
float sysSherpa=1.;

void initNormalizationSysts()
{
    normSysts["lumi_7TeV"] = 0.022;
    normSysts["lumi_8TeV"] = 0.026;
    normSysts["lumi_13TeV"] = 0.026;
    normSysts["accept_7TeV"] = 0.;//0.02;//0.003; //RJ
    normSysts["accept_8TeV"] = 0.;//0.02;//0.018; //RJ
    normSysts["sherpa_kin_syst"] = sysSherpa-1.0;
    normSysts["CMS_eff_e"] = 0.03;
    normSysts["CMS_eff_m"] = 0.04;
    normSysts["CMS_scale_e"] = 0.01; // do we need it? There's shape uncertainty "les"...
    normSysts["CMS_scale_m"] = 0.01; // do we need it? There's shape uncertainty "les"...
    normSysts["QCDscale_VV_zz_7Tev"] = 0.07021;
    normSysts["QCDscale_VV_wz_7Tev"] = 0.059;
    normSysts["QCDscale_VV_zz_8Tev"] = 0.0944;
    normSysts["QCDscale_VV_wz_8Tev"] = 0.054;
    normSysts["QCDscale_VV1in_zz_7Tev"] = 0.91657-1.0; // s'ha da fa' cosi'...
    normSysts["QCDscale_VV1in_zz_8Tev"] = 0.92617-1.0; // s'ha da fa' cosi'...
    normSysts["pdf_VV_zz_7TeV"] = 0.0115;
    normSysts["pdf_VV_wz_7TeV"] = 0.0116;
    normSysts["pdf_VV_zz_8TeV"] = 0.0112;
    normSysts["pdf_VV_wz_8TeV"] = 0.0120;
    //normSysts["sys_zlldata_7TeV"] = GammaJetSyst;
    //normSysts["sys_zlldata_8TeV"] = GammaJetSyst;
    //normSysts["sys_topwwwjetsdata_8TeV"] = NonResonnantSyst;
    //normSysts["sys_topwwwjetsdata_7TeV"] = NonResonnantSyst;
    normSysts["EM_7TeV"] = 1.0;
    normSysts["EM_8TeV"] = 1.0;
    //
    normSysts["CMS_zllwimps_mumueq0jets_leptonVeto"] = 0.01;
    normSysts["CMS_zllwimps_eeeq0jets_leptonVeto"] = 0.01;
    normSysts["CMS_zllwimps_mumueq1jets_leptonVeto"] = 0.013;
    normSysts["CMS_zllwimps_eeeq1jets_leptonVeto"] = 0.013;

    //
    normSysts["CMS_zllwimps_pdf"] = 0.;

    //unparticle
    normSysts["QCDscale_UnPart1p01"]=1.027561608;
    normSysts["QCDscale_UnPart1p02"]=1.027560594;
    normSysts["QCDscale_UnPart1p04"]=1.027107579;
    normSysts["QCDscale_UnPart1p06"]=1.026726974;
    normSysts["QCDscale_UnPart1p09"]=1.026031746;
    normSysts["QCDscale_UnPart1p10"]=1.025737817;
    normSysts["QCDscale_UnPart1p20"]=1.023378111;
    normSysts["QCDscale_UnPart1p30"]=1.019920319;
    normSysts["QCDscale_UnPart1p40"]=1.016645015;
    normSysts["QCDscale_UnPart1p50"]=1.012751678;
    normSysts["QCDscale_UnPart1p60"]=1.009233738;
    normSysts["QCDscale_UnPart1p70"]=1.004749134;
    normSysts["QCDscale_UnPart1p80"]=1.001779359;
    normSysts["QCDscale_UnPart1p90"]=1.006003132;
    normSysts["QCDscale_UnPart2p00"]=1.008754209;
    normSysts["QCDscale_UnPart2p20"]=1.017029328;

}

void printHelp()
{
    printf("Options\n");
    printf("--in        --> input file with from plotter\n");
    printf("--json      --> json file with the sample descriptor\n");
    printf("--histo     --> name of histogram to be used\n");
    printf("--shapeMin  --> left cut to apply on the shape histogram\n");
    printf("--shapeMax  --> right cut to apply on the shape histogram\n");
    printf("--shapeMinVBF  --> left cut to apply on the shape histogram for Vbf bin\n");
    printf("--shapeMaxVBF  --> right cut to apply on the shape histogram for Vbf bin\n");
    printf("--indexvbf  --> index of selection to be used for the vbf bin (if unspecified same as --index)\n");
    printf("--index     --> index of selection to be used (Xbin in histogram to be used)\n");
    printf("--indexL    --> index of selection to be used (Xbin in histogram to be used) used for interpolation\n");
    printf("--indexR    --> index of selection to be used (Xbin in histogram to be used) used for interpolation\n");
    printf("--m         --> higgs mass to be considered\n");
    printf("--mL        --> higgs mass on the left  of the mass to be considered (used for interpollation\n");
    printf("--mR        --> higgs mass on the right of the mass to be considered (used for interpollation\n");
    printf("--atgc      --> aTGC parameter (ex. string format: \"f4z=-0.01\")\n");
    printf("--syst      --> use this flag if you want to run systematics, default is no systematics\n");
    printf("--shape     --> use this flag if you want to run shapeBased analysis, default is cut&count\n");
    printf("--subNRB    --> use this flag if you want to subtract non-resonant-backgounds similarly to what was done in 2011 (will also remove H->WW)\n");
    printf("--subNRB12  --> use this flag if you want to subtract non-resonant-backgounds using a new technique that keep H->WW\n");
    printf("--subDY     --> histogram that contains the Z+Jets background estimated from Gamma+Jets)\n");
    printf("--subWZ     --> use this flag if you want to subtract WZ background by the 3rd lepton SB)\n");
    printf("--DDRescale --> factor to be used in order to multiply/rescale datadriven estimations\n");
    printf("--closure   --> use this flag if you want to perform a MC closure test (use only MC simulation)\n");
    printf("--bins      --> list of bins to be used (they must be comma separated without space)\n");
    printf("--HWW       --> use this flag to consider HWW signal)\n");
    printf("--skipGGH   --> use this flag to skip GGH signal)\n");
    printf("--skipQQH   --> use this flag to skip GGH signal)\n");
    printf("--blind     --> use this flag to replace observed data by total predicted background)\n");
    printf("--blindWithSignal --> use this flag to replace observed data by total predicted background+signal)\n");
    printf("--fast      --> use this flag to only do assymptotic prediction (very fast but inaccurate))\n");
    printf("--postfix    --> use this to specify a postfix that will be added to the process names)\n");
    printf("--systpostfix    --> use this to specify a syst postfix that will be added to the process names)\n");
    printf("--MCRescale    --> use this to specify a syst postfix that will be added to the process names)\n");
    printf("--interf     --> use this to rescale xsection according to WW interferences)\n");
    printf("--aTGC_Syst    --> use this to specify that you want to add and extra systematic in the shape)\n");
}

//
int main(int argc, char* argv[])
{
    setTDRStyle();
    gStyle->SetPadTopMargin   (0.06);
    gStyle->SetPadBottomMargin(0.12);
    //gStyle->SetPadRightMargin (0.16);
    gStyle->SetPadRightMargin (0.06);
    gStyle->SetPadLeftMargin  (0.14);
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetTitleXOffset(1.1);
    gStyle->SetTitleYOffset(1.45);
    gStyle->SetPalette(1);
    gStyle->SetNdivisions(505);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);


    //get input arguments
    for(int i=1; i<argc; i++) {
        string arg(argv[i]);
        if(arg.find("--help")          !=string::npos) {
            printHelp();
            return -1;
        } else if(arg.find("--subNRB12") !=string::npos) {
            subNRB2012=true;
            skipWW=false;
            printf("subNRB2012 = True\n");
        } else if(arg.find("--subNRB")   !=string::npos) {
            subNRB2011=true;
            skipWW=true;
            printf("subNRB2011 = True\n");
        } else if(arg.find("--subDY")    !=string::npos) {
            subDY=true;
            DYFile=argv[i+1];
            i++;
            printf("Z+Jets will be replaced by %s\n",DYFile.Data());
        } else if(arg.find("--subWZ")    !=string::npos) {
            subWZ=true;
            printf("WZ will be estimated from 3rd lepton SB\n");
        } else if(arg.find("--DDRescale")!=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&DDRescale);
            i++;
        } else if(arg.find("--MCRescale")!=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&MCRescale);
            i++;
        } else if(arg.find("--metEXTRAPOL") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%f",&DYMET_EXTRAPOL);
            i++;
            printf("DYMET_EXTRAPOL = %f\n", DYMET_EXTRAPOL);
        } else if(arg.find("--respEXTRAPOL") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%f",&DYRESP_EXTRAPOL);
            i++;
            printf("DYRESP_EXTRAPOL = %f\n", DYRESP_EXTRAPOL);
        } else if(arg.find("--dphiEXTRAPOL") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%f",&DYDPHI_EXTRAPOL);
            i++;
            printf("DYDPHI_EXTRAPOL = %f\n", DYDPHI_EXTRAPOL);
        } else if(arg.find("--HWW")      !=string::npos) {
            skipWW=false;
            printf("HWW = True\n");
        } else if(arg.find("--skipGGH")  !=string::npos) {
            skipGGH=true;
            printf("skipGGH = True\n");
        } else if(arg.find("--skipQQH")  !=string::npos) {
            skipQQH=true;
            printf("skipQQH = True\n");
        } else if(arg.find("--blindWithSignal")  !=string::npos) {
            blindData=true;
            blindWithSignal=true;
            printf("blindData = True; blindWithSignal = True\n");
        } else if(arg.find("--blind")    !=string::npos) {
            blindData=true;
            printf("blindData = True\n");
        } else if(arg.find("--closure")  !=string::npos) {
            MCclosureTest=true;
            printf("MCclosureTest = True\n");
        } else if(arg.find("--shapeMinVBF") !=string::npos && i+1<argc) {
            sscanf(argv[i+1],"%lf",&shapeMinVBF);
            i++;
            printf("Min cut on shape for VBF = %f\n", shapeMinVBF);
        } else if(arg.find("--shapeMaxVBF") !=string::npos && i+1<argc) {
            sscanf(argv[i+1],"%lf",&shapeMaxVBF);
            i++;
            printf("Max cut on shape for VBF = %f\n", shapeMaxVBF);
        } else if(arg.find("--shapeMin") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&shapeMin);
            i++;
            printf("Min cut on shape = %f\n", shapeMin);
        } else if(arg.find("--shapeMax") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%lf",&shapeMax);
            i++;
            printf("Max cut on shape = %f\n", shapeMax);
        } else if(arg.find("--interf")   !=string::npos) {
            doInterf=true;
            printf("doInterf = True\n");
        } else if(arg.find("--indexvbf") !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexvbf);
            i++;
            printf("indexVBF = %i\n", indexvbf);
        } else if(arg.find("--indexL")   !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcutL);
            i++;
            printf("indexL = %i\n", indexcutL);
        } else if(arg.find("--indexR")   !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcutR);
            i++;
            printf("indexR = %i\n", indexcutR);
        } else if(arg.find("--index")    !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&indexcut);
            i++;
            printf("index = %i\n", indexcut);
        } else if(arg.find("--in")       !=string::npos && i+1<argc)  {
            inFileUrl = argv[i+1];
            i++;
            printf("in = %s\n", inFileUrl.Data());
        } else if(arg.find("--json")     !=string::npos && i+1<argc)  {
            jsonFile  = argv[i+1];
            i++;
            printf("json = %s\n", jsonFile.Data());
        } else if(arg.find("--histo")    !=string::npos && i+1<argc)  {
            histo     = argv[i+1];
            i++;
            printf("histo = %s\n", histo.Data());
        } else if(arg.find("--mL")       !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&massL );
            i++;
            printf("massL = %i\n", massL);
        } else if(arg.find("--mR")       !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&massR );
            i++;
            printf("massR = %i\n", massR);
        } else if(arg.find("--mV")       !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&MV );
            i++;
            printf("MV = %i\n", MV);
        } else if(arg.find("--mA")       !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&MA );
            i++;
            printf("MA = %i\n", MA);
        } else if(arg.find("--m")        !=string::npos && i+1<argc)  {
            sscanf(argv[i+1],"%i",&mass );
            i++;
            printf("mass = %i\n", mass);
        } else if(arg.find("--atgc")     !=string::npos && i+1<argc)  {
            atgcpar = argv[i+1];
            i++;
            printf("aTGC parameter: %s\n", atgcpar.Data());
        } else if(arg.find("--bins")     !=string::npos && i+1<argc)  {
            char* pch = strtok(argv[i+1],",");
            printf("bins are : ");
            while (pch!=NULL) {
                printf(" %s ",pch);
                AnalysisBins.push_back(pch);
                pch = strtok(NULL,",");
            }
            printf("\n");
            i++;
        } else if(arg.find("--channels") !=string::npos && i+1<argc)  {
            char* pch = strtok(argv[i+1],",");
            printf("channels are : ");
            while (pch!=NULL) {
                printf(" %s ",pch);
                Channels.push_back(pch);
                pch = strtok(NULL,",");
            }
            printf("\n");
            i++;
        } else if(arg.find("--fast")     !=string::npos) {
            fast=true;
            printf("fast = True\n");
        } else if(arg.find("--postfix")  !=string::npos && i+1<argc)  {
            postfix = argv[i+1];
            systpostfix = argv[i+1];
            i++;
            printf("postfix '%s' will be used\n", postfix.Data());
        } else if(arg.find("--systpostfix") !=string::npos && i+1<argc)  {
            systpostfix = argv[i+1];
            i++;
            printf("systpostfix '%s' will be used\n", systpostfix.Data());
        } else if(arg.find("--syst")     !=string::npos) {
            runSystematics=true;
            printf("syst = True\n");
        } else if(arg.find("--shape")    !=string::npos) {
            shape=true;
            printf("shapeBased = True\n");
        } else if(arg.find("--aTGC_Syst")!=string::npos) {
            sscanf(argv[i+1],"%f",&sysSherpa);
            i++;
            printf("Additional systematic on the shape setted: %f\n", sysSherpa);
        }
    }

    if(jsonFile.IsNull() || inFileUrl.IsNull() || histo.IsNull() || indexcut == -1 || mass==-1) {
        printHelp();
        return -1;
    }
    if(AnalysisBins.size()==0)AnalysisBins.push_back("");
    if(Channels.size()==0) {
        Channels.push_back("ee");
        Channels.push_back("mumu");
        //Channels.push_back("ll"); //RENJIE, add all
    }

    atgcpar2 = atgcpar;
    atgcpar2.ReplaceAll("-", "M");
    atgcpar2.ReplaceAll("+", "P");
    atgcpar2.ReplaceAll("=0.", "=P0.");
    atgcpar2.ReplaceAll(".", "p");
    atgcpar2.ReplaceAll("=", "_");

    //if(systpostfix.Contains('8')) {NonResonnantSyst = 0.25; GammaJetSyst = 0.40;}
    if(atgcpar.Length()>0) GammaJetSyst = 0.0;

    initNormalizationSysts();

    //build the datacard for this mass point
    //std::vector<TString> dcUrls = buildDataCard(mass,histo,inFileUrl, jsonFile);
    std::vector<TString> dcUrls = buildDataCard(atgcpar,mass,histo,inFileUrl, jsonFile);

    fOut_binbybinNuisance->cd();
    nuisance_h->Write("Bin by Bin Nuisance");
    fOut_binbybinNuisance->Close();
}




//
Shape_t getShapeFromFile(TFile* inF, TString ch, TString shapeName, int cutBin, JSONWrapper::Object &Root, double minCut, double maxCut, bool onlyData)
{
    gROOT->cd();  //THIS LINE IS NEEDED TO MAKE SURE THAT HISTOGRAM INTERNALLY PRODUCED IN LumiReWeighting ARE NOT DESTROYED WHEN CLOSING THE FILE

    Shape_t shape;
    shape.totalBckg=NULL;
    shape.data=NULL;

    std::vector<TString> BackgroundsInSignal;

    //iterate over the processes required
    std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
    for(unsigned int i=0; i<Process.size(); i++) {
        TString procCtr("");
        procCtr+=i;
        TString proc=(Process[i])["tag"].toString();
        TDirectory *pdir = (TDirectory *)inF->Get(proc);
        if(pdir==0) {
            /*printf("Skip Proc=%s because its directory is missing in root file\n", proc.Data());*/ continue;
        }

        bool isData(Process[i]["isdata"].toBool());
        if(onlyData && !isData)continue; //just here to speedup the NRB prediction

        bool isSignal(Process[i].isTag("issignal") && Process[i]["issignal"].toBool());
        if(Process[i]["spimpose"].toBool() && (proc.Contains("ggH") || proc.Contains("qqH")))isSignal=true;
        bool isInSignal(Process[i].isTag("isinsignal") && Process[i]["isinsignal"].toBool());
        int color(1);
        if(Process[i].isTag("color" ) ) color  = (int)Process[i]["color" ].toInt();
        int lcolor(color);
        if(Process[i].isTag("lcolor") ) lcolor = (int)Process[i]["lcolor"].toInt();
        int mcolor(color);
        if(Process[i].isTag("mcolor") ) mcolor = (int)Process[i]["mcolor"].toInt();
        //int fcolor(color);
        //if(Process[i].isTag("fcolor") ) fcolor = (int)Process[i]["fcolor"].toInt();
        int lwidth(1);
        if(Process[i].isTag("lwidth") ) lwidth = (int)Process[i]["lwidth"].toInt();
        int lstyle(1);
        if(Process[i].isTag("lstyle") ) lstyle = (int)Process[i]["lstyle"].toInt();
        int fill(1001);
        if(Process[i].isTag("fill"  ) ) fill   = (int)Process[i]["fill"  ].toInt();
        int marker(20);
        if(Process[i].isTag("marker") ) marker = (int)Process[i]["marker"].toInt();

        TH1* syst = (TH1*)pdir->Get("optim_systs");
        if(syst==NULL) {
            cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>> NULL" << endl;
            syst =new TH1F("optim_systs","optim_systs",1,0,1);
            syst->GetXaxis()->SetBinLabel(1,"");
        }
        for(int ivar = 1; ivar<=syst->GetNbinsX(); ivar++) {
            TH1D* hshape   = NULL;

            TString varName = syst->GetXaxis()->GetBinLabel(ivar);
            TString histoName = ch+"_"+shapeName+varName ;
            TH2* hshape2D = (TH2*)pdir->Get(histoName );
            if(!hshape2D) {
//            if(varName==""){
                //replace by empty histogram (take inclusive histo to make sure it has same binning)
                hshape2D = (TH2*)pdir->Get(shapeName+varName);
                if(hshape2D)hshape2D->Reset();
//            }else{
//               continue;
//            }
            }

            if(hshape2D) {
                histoName.ReplaceAll(ch,ch+"_proj"+procCtr);
                hshape   = hshape2D->ProjectionY(histoName,cutBin,cutBin);
                if(hshape->Integral()<=0 && varName=="" && !isData) {
                    hshape->Reset();
                    hshape->SetBinContent(1, 1E-10);
                }

                if(isnan((float)hshape->Integral())) {
                    hshape->Reset();
                }
                hshape->SetDirectory(0);
                hshape->SetTitle(proc);
                fixExtremities(hshape,true,true);
                hshape->SetFillColor(color);
                hshape->SetLineColor(lcolor);
                hshape->SetMarkerColor(mcolor);
                hshape->SetFillStyle(fill);
                hshape->SetLineWidth(lwidth);
                hshape->SetMarkerStyle(marker);
                hshape->SetLineStyle(lstyle);
            } else {
                if(!histoName.Contains("FR_WjetCtrl") && !histoName.Contains("FR_QCDCtrl")
                        && !histoName.Contains("minus_shapes") && !histoName.Contains("_NRBctrl") && !histoName.Contains("_NRBsyst"))
                    printf("Histo %s does not exist for syst:%s\n", histoName.Data(), varName.Data());
                continue;
            }


            //if current shape is the one to cut on, then apply the cuts
            if(shapeName == histo) {
                for(int x=0; x<=hshape->GetXaxis()->GetNbins()+1; x++) {
                    if(hshape->GetXaxis()->GetBinCenter(x)<=minCut || hshape->GetXaxis()->GetBinCenter(x)>=maxCut) {
                        hshape->SetBinContent(x,0);
                        hshape->SetBinError(x,0);
                    }
                }
                //hshape->Rebin(2);
                hshape->GetYaxis()->SetTitle("Entries (/25GeV)");
            }

            hshape->Scale(MCRescale);



            //save in structure
            if(isData) {
                if(varName=="")  shape.data=hshape;
                else continue;
            } else if(isSignal) {
                if(skipGGH && proc.Contains("ggH"))continue;
                if(skipQQH && proc.Contains("qqH"))continue;

                if(skipWW && string(proc.Data()).find("WW")!=string::npos )continue;
                if(!skipWW && mergeWWandZZ) {
                    proc.ReplaceAll("WW","VV");
                    proc.ReplaceAll("ZZ","VV");
                }

                if(varName=="") {
                    if(Process[i]["data"].daughters()[0].isTag("xsec"))shape.xsections[proc] = Process[i]["data"].daughters()[0]["xsec"].toDouble();
                    if(Process[i]["data"].daughters()[0].isTag("br")) {
                        std::vector<JSONWrapper::Object> BRs = Process[i]["data"].daughters()[0]["br"].daughters();
                        double totalBR=1.0;
                        for(size_t ipbr=0; ipbr<BRs.size(); ipbr++) {
                            totalBR*=BRs[ipbr].toDouble();
                        }
                        shape.BRs[proc] = totalBR;
                    }

                    int procIndex = -1;
                    for(unsigned int i=0; i<shape.signal.size(); i++) {
                        if(string(proc.Data())==shape.signal[i]->GetTitle() ) {
                            procIndex=i;
                            break;
                        }
                    }
                    if(procIndex>=0) shape.signal[procIndex]->Add(hshape);
                    else             {
                        hshape->SetTitle(proc);
                        shape.signal.push_back(hshape);
                    }

                    //printf("Adding signal %s\n",proc.Data());
                } else {
                    std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator it = shape.signalVars.find(proc);

                    bool newVar = true;
                    if(it!=shape.signalVars.end()) {
                        for(unsigned int i=0; i<it->second.size(); i++) {
                            if( string(it->second[i].first.Data()) == varName ) {
                                it->second[i].second->Add(hshape);
                                newVar=false;
                                break;
                            }
                        }
                    }

                    if(newVar) {
                        shape.signalVars[proc].push_back( std::pair<TString,TH1*>(varName,hshape) );
                    }
                }
            } else {
                if(isInSignal) {
                    if(varName=="")  BackgroundsInSignal.push_back(proc);
                    if(varName=="")  shape.bckgInSignal.push_back(hshape);
                    else             shape.bckgInSignalVars[proc].push_back( std::pair<TString,TH1*>(varName,hshape) );
                } else {
                    if(varName=="")  shape.bckg.push_back(hshape);
                    else             shape.bckgVars[proc].push_back( std::pair<TString,TH1*>(varName,hshape) );
                }

                //printf("histoName = B %i -- %i  -- %s - %s --> %s\n", i, int(varName==""), proc.Data(), histoName.Data(), hshape->GetTitle());
            }
        }
        //delete syst;
    }

    //compute the total
    for(size_t i=0; i<shape.bckg.size(); i++) {
        if(i==0) {
            shape.totalBckg = (TH1 *)shape.bckg[i]->Clone(ch+"_"+shapeName+"_total");
            shape.totalBckg->SetDirectory(0);
        } else     {
            shape.totalBckg->Add(shape.bckg[i]);
        }
    }

    if(MCclosureTest) {
        if(shape.totalBckg) {
            if(!shape.data) {
                shape.data=(TH1F*)shape.totalBckg->Clone("data");
                shape.data->SetDirectory(0);
                shape.data->SetTitle("data");
            } else {
                shape.data->Reset();
                shape.data->Add(shape.totalBckg, 1);
            }
        }
    }


    //subtract background that are included to signal sample
    //if(shape.signal.size()>1 && BackgroundsInSignal.size()>0)printf("YOU ARE TRYING TO SUBTRACT BACKGROUNDS FROM SIGNAL SAMPLE WHILE YOU HAVE MORE THAN ONE SIGNAL GIVEN!  ASSUME THAT FIRST SIGNAL SAMPLE IS THE ONE CONTAINING THE BACKGROUND\n");
    if(shape.signal.size()>1 && BackgroundsInSignal.size()>0) printf("YOU ARE TRYING TO SUBTRACT BACKGROUNDS FROM SIGNAL SAMPLE WHILE YOU HAVE MORE THAN ONE SIGNAL GIVEN! ASSUME THAT ALL SIGNAL SAMPLES CONTAIN THE BACKGROUND, AND SUBTRACT THEM ALL!\n");
    //for(size_t i=0;i<BackgroundsInSignal.size()  && shape.signal.size()>0  ;i++){
    //deal with central value
    //for(size_t j=0;j<shape.bckg.size();j++){
    for(size_t j=0; j<shape.bckgInSignal.size(); j++) {
        //printf("Compare %s and %s\n",shape.bckg[j]->GetTitle(), BackgroundsInSignal[i].Data());
        //if(BackgroundsInSignal[i]!=shape.bckg[j]->GetTitle())continue;
        //printf("Subtract %s to %s\n",shape.signal[0]->GetTitle(), BackgroundsInSignal[i].Data());
        //shape.signal[0]->Add(shape.bckg[j],-1);
        std::map<TString, std::vector<double> > atgcShiftMap;
        for(size_t k=0; k<shape.signal.size(); ++k) {
            printf("Subtract %s to %s\n", shape.signal[k]->GetTitle(), shape.bckgInSignal[j]->GetTitle());
            shape.signal[k]->Add(shape.bckgInSignal[j],-1);
            // Get bin values for interpolated points
            TString atgcCorrIdx = ch+"_"+shape.signal[k]->GetTitle();
            if(systpostfix.Contains("7"))      atgcCorrIdx.Prepend("7TeV_");
            else if(systpostfix.Contains("8")) atgcCorrIdx.Prepend("8TeV_");
            std::vector<double> atgcCorrVect; //, atgcErrVect;
            //if(extrvalpoints.count(atgcCorrIdx.Data())>0) atgcCorrVect = extrvalpoints[atgcCorrIdx.Data()];
            //else std::cout << "No ATGC correction found with label " << atgcCorrIdx.Data() << "!" << std::endl;
            //if(extrerrpoints.count(atgcCorrIdx.Data())>0) atgcErrVect = extrerrpoints[atgcCorrIdx.Data()];
            if(atgcCorrVect.size()!=((unsigned int)shape.signal[k]->GetNbinsX())) {
                std::cout << "ATGC correction has different size than signal shape!" << std::endl;
                std::cout << "           atgcCorrVect.size() = " << atgcCorrVect.size() << "  (" << atgcCorrIdx.Data() << ")" << std::endl;
                std::cout << "  shape.signal[k]->GetNbinsX() = " << shape.signal[k]->GetNbinsX() << "  (" << shape.signal[k]->GetTitle() << ")" << std::endl;
                atgcCorrVect.clear(); //atgcErrVect.clear();
            }
            std::vector<double> atgcShiftVec(atgcCorrVect.size(), 0.);
            for(int x=0; x<=shape.signal[k]->GetNbinsX()+1; x++) {
                if(shape.signal[k]->GetBinContent(x)<0)shape.signal[k]->SetBinContent(x,0.0);
                if(atgcCorrVect.size()==0) continue;
                if(x==0 || x==shape.signal[k]->GetNbinsX()+1) continue;
                if(atgcCorrVect[x-1]<0.) continue;
                atgcShiftVec[x-1] = atgcCorrVect[x-1] - shape.signal[k]->GetBinContent(x);
                // before changing, save the difference (the syst plots will be shifted by the same amount)

                shape.signal[k]->SetBinContent(x, atgcCorrVect[x-1]);
                // replace the content of the bin with the expected value for the extrapolated point
            }
            if(atgcShiftVec.size()>0) atgcShiftMap[atgcCorrIdx] = atgcShiftVec;
        }
        //}

        //deal with systematics
        //const std::vector<std::pair<TString, TH1*> >& bckgSysts = shape.bckgVars  [BackgroundsInSignal[i] ];
        const std::vector<std::pair<TString, TH1*> >& bckgSysts = shape.bckgInSignalVars[shape.bckgInSignal[j]->GetTitle()];
        for(size_t k=0; k<shape.signal.size(); ++k) {
            const std::vector<std::pair<TString, TH1*> >& signSysts = shape.signalVars[shape.signal[k]->GetTitle() ];
            if(bckgSysts.size()!=signSysts.size())printf("Problem the two vectors have different size!\n");
            TString atgcCorrIdx = ch+"_"+shape.signal[k]->GetTitle();
            if(systpostfix.Contains("7"))      atgcCorrIdx.Prepend("7TeV_");
            else if(systpostfix.Contains("8")) atgcCorrIdx.Prepend("8TeV_");
            for(size_t jsys=0; jsys<bckgSysts.size(); jsys++) {
                signSysts[jsys].second->Add(bckgSysts[jsys].second, -1);
                for(int x=0; x<=signSysts[jsys].second->GetNbinsX()+1; x++) {
                    if(signSysts[jsys].second->GetBinContent(x)<0)signSysts[jsys].second->SetBinContent(x,0.0);
                    if(atgcShiftMap.count(atgcCorrIdx)==0) continue;
                    if(x==0 || x==shape.signal[k]->GetNbinsX()+1) continue;
                    signSysts[jsys].second->SetBinContent( x, signSysts[jsys].second->GetBinContent(x) + atgcShiftMap[atgcCorrIdx][x-1] );
                    // apply to the syst plots the same shift used for the main distribution
                }
            }
        }
    }



    //all done
    return shape;
}

//
void showShape(std::vector<TString>& selCh ,map<TString, Shape_t>& allShapes, TString mainHisto, TString SaveName)
{
    cout << "########################## showShape ##########################" << endl;
    std::map<TString, TH1*> CCHistos;

    TCanvas* c1 = new TCanvas("c1","c1",300*AnalysisBins.size(),300*selCh.size());
    c1->SetTopMargin(0.00);
    c1->SetRightMargin(0.00);
    c1->SetBottomMargin(0.00);
    c1->SetLeftMargin(0.00);
    TPad* t1 = new TPad("t1","t1", 0.03, 0.03, 1.0, 1.00);
    t1->Draw();
    t1->cd();
    t1->Divide(AnalysisBins.size(), selCh.size(), 0, 0);

    t1->cd(selCh.size()*AnalysisBins.size());
    TLegend* legA  = new TLegend(0.6,0.5,0.99,0.85, "NDC");
    legA->SetNColumns(1); //RJ
    t1->cd(0);

    bool haserrorsleg(false);

    for(size_t s=0; s<selCh.size(); s++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            TVirtualPad* pad = t1->cd(1+s*AnalysisBins.size()+b);
            pad->SetTopMargin(0.06);
            pad->SetRightMargin(0.03);
            pad->SetBottomMargin(0.07);
            pad->SetLeftMargin(0.06);

            TH1* allbkg=NULL;
            std::map<TString, TH1*> mapbkg;
            std::map<TString, TH1*> mapsig;
            TH1* alldata=NULL;

            Shape_t& shape = allShapes.find(selCh[s]+AnalysisBins[b]+mainHisto)->second;
            cout << " Limit setting shapes: " << selCh[s]+AnalysisBins[b]+"_"+mainHisto << endl;

            if(!allbkg) {
                allbkg=(TH1*)shape.totalBckg->Clone("mc");
            } else {
                allbkg->Add(shape.totalBckg);
            }
            if(!alldata) {
                alldata=(TH1*)shape.data->Clone("data");
            } else {
                alldata->Add(shape.data);
            }

            for(size_t i=0; i<shape.bckg.size(); i++) {
                if(shape.bckg[i]->Integral()<=1E-6) continue;
                cout << " Backgrounds: " << shape.bckg[i]->GetTitle() << endl;
                if(mapbkg.find(shape.bckg[i]->GetTitle())!=mapbkg.end()) {
                    mapbkg[shape.bckg[i]->GetTitle()]->Add(shape.bckg[i]);
                } else {
                    mapbkg[shape.bckg[i]->GetTitle()]=(TH1*)shape.bckg[i]->Clone(shape.bckg[i]->GetTitle());
                }

                //cut & count plot
                double VAL, ERR;
                VAL = shape.bckg[i]->IntegralAndError(1,shape.bckg[i]->GetXaxis()->GetNbins(),ERR);
                if(CCHistos.find(shape.bckg[i]->GetTitle())==CCHistos.end()) {
                    CCHistos[shape.bckg[i]->GetTitle()] = new TH1D(TString(shape.bckg[i]->GetTitle())+"CC", "", 25,0,25);
                    CCHistos[shape.bckg[i]->GetTitle()]->SetLineColor(shape.bckg[i]->GetLineColor());
                    CCHistos[shape.bckg[i]->GetTitle()]->SetLineStyle(shape.bckg[i]->GetLineStyle());
                    CCHistos[shape.bckg[i]->GetTitle()]->SetLineWidth(shape.bckg[i]->GetLineWidth());
                    CCHistos[shape.bckg[i]->GetTitle()]->SetFillColor(shape.bckg[i]->GetFillColor());
                    CCHistos[shape.bckg[i]->GetTitle()]->SetFillStyle(shape.bckg[i]->GetFillStyle());
                }
                CCHistos[shape.bckg[i]->GetTitle()]->SetBinContent(1+s*AnalysisBins.size()+b,VAL);
                CCHistos[shape.bckg[i]->GetTitle()]->SetBinError(1+s*AnalysisBins.size()+b,VAL);
            }

            TString massStr("");
            if(mass>0) massStr += mass;
            if(atgcpar.Length()>0) massStr = atgcpar;
            TString massStr2 = atgcpar2;
            for(size_t ip=0; ip<shape.signal.size(); ip++) {
                TString proc(shape.signal[ip]->GetTitle());
                if(mass>0 && !proc.Contains(massStr))continue;
                if(mass>0 && proc.Contains("ggH") && proc.Contains("ZZ"))proc = "ggHZZ2l2v";
                else if(mass>0 && proc.Contains("D1")			     )proc = "D1_"+massStr+"GeV";
                else if(mass>0 && proc.Contains("D4")                        )proc = "D4_"+massStr+"GeV";
                else if(mass>0 && proc.Contains("D5")                        )proc = "D5_"+massStr+"GeV";
                else if(mass>0 && proc.Contains("D8")                        )proc = "D8_"+massStr+"GeV";
                else if(mass>0 && proc.Contains("D9")                        )proc = "D9_"+massStr+"GeV";
                else if(mass>0 && proc.Contains("C3")                        )proc = "C3_"+massStr+"GeV";
                else if(mass>0 && proc.Contains("qqH") && proc.Contains("ZZ"))proc = "qqHZZ2l2v";
                else if(mass>0 && proc.Contains("ggH") && proc.Contains("WW"))proc = "ggHWW2l2v";
                else if(mass>0 && proc.Contains("qqH") && proc.Contains("WW"))proc = "qqHWW2l2v";
                else if(mass>0 && proc.Contains("ZH")                        )proc = "ZH"+massStr+"2lMET";
                else continue;

                if(proc=="qqHZZ") {
                    shape.signal[ip]->SetLineStyle(2);
                }

                //if(atgcpar.Length()>0 && !proc.Contains(massStr)) continue;
                if(atgcpar.Length()>0 && !proc.EndsWith(massStr)) continue;
                if(atgcpar.Length()>0) proc = "ZZ2l2nu_"+massStr2;

                if(mapsig.find(shape.signal[ip]->GetTitle())!=mapsig.end()) {
                    mapsig[shape.signal[ip]->GetTitle()]->Add(shape.signal[ip]);
                } else {
                    mapsig[shape.signal[ip]->GetTitle()]=(TH1*)shape.signal[ip]->Clone(shape.signal[ip]->GetTitle());
                }

                //cut & count plot
                double VAL, ERR;
                VAL = shape.signal[ip]->IntegralAndError(1,shape.signal[ip]->GetXaxis()->GetNbins(),ERR);
                if(CCHistos.find(shape.signal[ip]->GetTitle())==CCHistos.end()) {
                    CCHistos[shape.signal[ip]->GetTitle()] = new TH1D(TString(shape.signal[ip]->GetTitle())+"CC", "", 25,0,25);
                    CCHistos[shape.signal[ip]->GetTitle()]->SetLineColor(shape.signal[ip]->GetLineColor());
                    CCHistos[shape.signal[ip]->GetTitle()]->SetLineStyle(shape.signal[ip]->GetLineStyle());
                    CCHistos[shape.signal[ip]->GetTitle()]->SetLineWidth(shape.signal[ip]->GetLineWidth());
                    CCHistos[shape.signal[ip]->GetTitle()]->SetFillColor(shape.signal[ip]->GetFillColor());
                    CCHistos[shape.signal[ip]->GetTitle()]->SetFillStyle(shape.signal[ip]->GetFillStyle());
                }
                CCHistos[shape.signal[ip]->GetTitle()]->SetBinContent(1+s*AnalysisBins.size()+b,VAL);
                CCHistos[shape.signal[ip]->GetTitle()]->SetBinError(1+s*AnalysisBins.size()+b,ERR);
            }

            THStack *stack=0;
            TH1* mc=allbkg;
            TH1* mcPlusRelUnc=0;
            TGraphErrors* errors=0;

            if(allbkg) {
                mcPlusRelUnc = (TH1 *) allbkg->Clone("totalmcwithunc");
                mcPlusRelUnc->SetDirectory(0);
                mcPlusRelUnc->Reset();
                stack = new THStack("stack","stack");
                for(std::map<TString, TH1*>::iterator it=mapbkg.begin(); it!=mapbkg.end(); ++it) {
                    it->second->SetLineColor( it->second->GetFillColor());
                    stack->Add(it->second,"HIST");
                    if(s==0 && b==0)legA->AddEntry(it->second,it->second->GetTitle(),"F");

                    double baseRelUnc = it->second->GetBinError(0)/it->second->Integral();
                    //if(TString(it->second->GetTitle()).Contains("rightarrow ll (data)")){printf("replace uncertainty %g",baseRelUnc); baseRelUnc=1.0;}
                    for(int ibin=1; ibin<=mcPlusRelUnc->GetXaxis()->GetNbins(); ibin++) {
                        double val = it->second->GetBinContent(ibin);
                        double err = it->second->GetBinError(ibin);
                        double value = mcPlusRelUnc->GetBinContent(ibin) + val;
                        double error = sqrt(pow(mcPlusRelUnc->GetBinError(ibin),2) + pow(err,2) + pow(val*baseRelUnc,2));

                        //Syst Err from Shape
                        std::vector<float> ErrShape_tot;
                        ErrShape_tot.clear();
                        // Loop on proc
                        //std::cout << " ************** " << it->first.Data() << " ************** " << std::endl;
                        for(std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator jt=shape.bckgVars.begin(); jt!=shape.bckgVars.end(); ++jt) {
                            if(it->first.Contains("data")) continue;
                            if(it->first.CompareTo(jt->first)!=0) continue;
                            //std::cout << " ************** " << jt->first.Data() << " ************** " << std::endl;
                            //float Err(0);
                            //cout<<"proc: "<<jt->first<<endl; //Z#rightarrow ll
                            std::vector<float> ErrShape_up;
                            ErrShape_up.clear();
                            std::vector<float> ErrShape_down;
                            ErrShape_down.clear();
                            // Loop on syst
                            for(int isel=0; isel<12 ; isel++) {
                                if((jt->second)[isel].first.Contains("up")) {
                                    //cout << "isel: " << isel << ":  " << (jt->second)[isel].first <<endl;
                                    TH1F *h1 = (TH1F*) (jt->second)[isel].second;
                                    float err_tmp = h1->GetBinContent(ibin);
                                    err_tmp = (err_tmp-value)/value;
                                    //cout<<"Syst Name: "<<(jt->second)[isel].first<<" : "<<h1->GetBinContent(ibin)<<" "<<value<<" -> "<<err_tmp<<endl;
                                    if(err_tmp!=-1 && value!=0) ErrShape_up.push_back(err_tmp);
                                    else                        ErrShape_up.push_back(0.);
                                } else if((jt->second)[isel].first.Contains("down")) {
                                    //cout << "isel: " << isel << ":  " << (jt->second)[isel].first <<endl;
                                    TH1F *h1 = (TH1F*) (jt->second)[isel].second;
                                    float err_tmp = h1->GetBinContent(ibin);
                                    err_tmp = (err_tmp-value)/value;
                                    //cout<<"Syst Name: "<<(jt->second)[isel].first<<" : "<<h1->GetBinContent(ibin)<<" "<<value<<" -> "<<err_tmp<<endl;
                                    if(err_tmp!=-1 && value!=0) ErrShape_down.push_back(err_tmp);
                                    else                        ErrShape_down.push_back(0.);
                                } else {
                                    cout<<"WARNING: Syst not computed properly: "<<endl;
                                }
                            }//syst
                            std::vector<float> ErrShape;
                            ErrShape.clear();
                            for(unsigned int isel=0; isel<ErrShape_down.size() ; isel++) {
                                if(ErrShape_down.size() != ErrShape_up.size()) cout<<"WARNING: Syst not computed properly"<<endl;
                                ErrShape.push_back( std::max(fabs(ErrShape_down[isel]), fabs(ErrShape_up[isel])) );
                            }
                            // Non-shape uncertainties ("lnN")
                            for(std::map<TString, double>::const_iterator kt=normSysts.begin(); kt!=normSysts.end(); ++kt) {
                                if(kt->first.Contains("7TeV") && !systpostfix.Contains("7TeV")) continue;
                                if(kt->first.Contains("8TeV") && !systpostfix.Contains("8TeV")) continue;
                                if(kt->first.EndsWith("_e") && !selCh[s].Contains("ee")) continue;
                                if(kt->first.EndsWith("_m") && !selCh[s].Contains("mumu")) continue;
                                if(kt->first.Contains("_zz_") && !it->first.Contains("ZZ")) continue;
                                if(kt->first.Contains("_wz_") && !it->first.Contains("WZ")) continue;
                                if(kt->first.Contains("_zlldata_") && !it->first.Contains("Z#rightarrow ll (data)")) continue;
                                if(kt->first.Contains("_topwwztautaudata_") && !it->first.Contains("Top/WW/Ztautau (data)")) continue;

                                ErrShape.push_back(kt->second);
                            }
                            float TotErrShape(0);
                            for(unsigned int iEr=0; iEr<ErrShape.size(); iEr++) TotErrShape+=pow(ErrShape[iEr],2);
                            TotErrShape=sqrt(TotErrShape);
                            //cout<<"Sing process Unc:"<<TotErrShape<<endl;
                            ErrShape_tot.push_back(TotErrShape);
                        }//process
                        float TotErrShape_tot(0.);
                        for(unsigned int iEr=0; iEr<ErrShape_tot.size(); iEr++) TotErrShape_tot+=pow(ErrShape_tot[iEr],2);
                        TotErrShape_tot=sqrt(TotErrShape_tot);
                        //cout<<"BIN unc:"<<TotErrShape_tot<<" , the other is: "<<error<<" summed are: "<<sqrt(TotErrShape_tot*TotErrShape_tot+error*error)<<endl;

                        //Syst Err not form datacard (harcoded)

                        mcPlusRelUnc->SetBinContent(ibin,value);
                        mcPlusRelUnc->SetBinError(ibin,error);
                    }
                }


                double mcErorMax=0;
                errors = new TGraphErrors(mcPlusRelUnc->GetXaxis()->GetNbins());
                int icutg=0;
                for(int ibin=1; ibin<=mcPlusRelUnc->GetXaxis()->GetNbins(); ibin++) {
                    if(mcPlusRelUnc->GetBinContent(ibin)>0)
                        errors->SetPoint(icutg,mcPlusRelUnc->GetXaxis()->GetBinCenter(ibin), mcPlusRelUnc->GetBinContent(ibin));
                    errors->SetPointError(icutg,mcPlusRelUnc->GetXaxis()->GetBinWidth(ibin)/2.0, mcPlusRelUnc->GetBinError(ibin));
                    icutg++;
                    mcErorMax = std::max(mcErorMax, mcPlusRelUnc->GetBinContent(ibin) + mcPlusRelUnc->GetBinError(ibin));
                }

                TH1* axis = (TH1*)allbkg->Clone("axis");
                axis->Reset();
                axis->GetXaxis()->SetTitle(mc->GetXaxis()->GetTitle());
                axis->GetYaxis()->SetTitle(b==0?mc->GetYaxis()->GetTitle():"");
                //axis->SetMinimum(mc->GetMinimum());
                //axis->SetMaximum(1.1*std::max(mcErorMax, alldata->GetMaximum()));
                axis->SetMaximum(1.9*std::max(mcErorMax, alldata->GetMaximum()));
                //axis->GetXaxis()->SetRangeUser(150,700);//
                axis->GetXaxis()->SetRangeUser(0,400);//
                axis->Draw();
                stack->Draw("same");


                errors->Set(icutg);
                //errors->SetFillStyle(3427);
                //errors->SetFillStyle(3002);
                errors->SetFillStyle(3254); //RJ
                errors->SetFillColor(kGray+3);
                errors->SetLineStyle(1);
                errors->SetLineColor(0);
                errors->Draw("2 same");
                if(!haserrorsleg) {
                    legA->AddEntry(errors,"Syst.+Stat. Unc.","F"); //RJ
                    haserrorsleg=true;
                }
            }


            std::vector<TH1 *> sigvec;
            for(std::map<TString, TH1*>::iterator it=mapsig.begin(); it!=mapsig.end(); it++) {
                it->second->SetFillStyle(0);
                it->second->Draw("hist same");
                if(s==0 && b==0)legA->AddEntry(it->second,it->second->GetTitle(),"L");
                sigvec.push_back(it->second);
            }



            if(alldata) {
                alldata->Draw("E1 same");
                if(s==0 && b==0)legA->AddEntry(alldata,alldata->GetTitle(),"PL");
            }


            TPaveText* Label = new TPaveText(0.1,0.81,0.94,0.89, "NDC");
            Label->SetFillColor(0);
            Label->SetFillStyle(0);
            Label->SetLineColor(0);
            Label->SetBorderSize(0);
            Label->SetTextAlign(31);
            TString LabelText = selCh[s]+"  -  "+AnalysisBins[b];
            LabelText.ReplaceAll("mumu","#mu#mu");
            LabelText.ReplaceAll("geq2jets","#geq2jets");
            LabelText.ReplaceAll("eq0jets","0jet");
            LabelText.ReplaceAll("eq1jets","1jet");
            Label->AddText(LabelText);
            Label->Draw();



            if(s==selCh.size()-1 && b==AnalysisBins.size()-1) {
                legA->SetFillColor(0);
                legA->SetFillStyle(0);
                legA->SetLineColor(0);
                legA->SetBorderSize();
                legA->SetHeader("");
                legA->Draw("same");
                legA->SetTextFont(42);
            }

            //Draw1DHistogram(mc, stack, mcPlusRelUnc, errors, sigvec, alldata, legA, true, selCh[s], AnalysisBins[b]);
            //if(histo.Contains("dphiLL_mt_shapes"))
            //    Draw2DHistogram(mapbkg, sigvec, alldata, selCh[s], AnalysisBins[b]);


            /*
              TH1 *ratio=0;
              if(allbkg && alldata){
                  c1->cd();
                  TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.2);     t2->Draw();
                  t2->cd();
                  t2->SetGridy(true);
                  t2->SetTopMargin(0);   t2->SetBottomMargin(0.5);
                  float yscale = (1.0-0.2)/(0.18-0);
                  TH1 *ratio = (TH1*)alldata->Clone("RatioHistogram");
                  ratio->SetDirectory(0);
                  ratio->Divide(allbkg);
                  ratio->GetYaxis()->SetTitle("Obs/Ref");
                  ratio->GetXaxis()->SetTitle("");
                  ratio->SetMinimum(0);
                  ratio->SetMaximum(2.2);
                  ratio->GetXaxis()->SetTitleOffset(1.3);
                  ratio->GetXaxis()->SetLabelSize(0.033*yscale);
                  ratio->GetXaxis()->SetTitleSize(0.036*yscale);
                  ratio->GetXaxis()->SetTickLength(0.03*yscale);
                  ratio->GetYaxis()->SetTitleOffset(0.3);
                  ratio->GetYaxis()->SetNdivisions(5);
                  ratio->GetYaxis()->SetLabelSize(0.033*yscale);
                  ratio->GetYaxis()->SetTitleSize(0.036*yscale);
                  ratio->GetXaxis()->SetRangeUser(175,450);
                  ratio->Draw("E1");
              }
            */
            pad->SetLogy();
            pad->Update();
        }
    }

    t1->cd();
    t1->Update();
    c1->cd();
    TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetBorderSize(0);
    T->SetTextAlign(22);
    if(systpostfix.Contains('8')) {
        T->AddText("CMS preliminary, ZH #rightarrow ll+MET, #sqrt{s}=8.0 TeV, #scale[0.5]{#int} L=19.6  fb^{-1}");
    } else {
        T->AddText("CMS preliminary, ZH #rightarrow ll+MET, #sqrt{s}=7.0 TeV, #scale[0.5]{#int} L=5.1  fb^{-1}");
    }
    T->Draw();
    c1->Update();
    //c1->SaveAs(SaveName+"_Shape.eps");
    //c1->SaveAs(SaveName+"_Shape.png");
    //c1->SaveAs(SaveName+"_Shape.pdf");
    //c1->SaveAs(SaveName+"_Shape.C");
    delete c1;
    cout << "########################## showShape ##########################" << endl;
}


void Draw1DHistogram(TH1* mc, THStack *stack, TH1 *mcPlusRelUnc, TGraphErrors *errors, std::vector<TH1 *>& spimpose, TH1 *data, TLegend *legA, bool noLog, TString finstate, TString AnalysisBins)
{

    bool showsigvsBkg(true);
    bool isDataBlind = true;

    TCanvas* c1d = new TCanvas("c1d","c1d",750,800);//,800,800);

    TPad* t1d = new TPad();
    if(!isDataBlind) t1d = new TPad("t1d","t1d", 0.0, 0.3, 1.0, 1.0);
    else t1d = new TPad("t1d","t1d", 0.0, 0.0, 1.0, 1.0);

    t1d->Draw();
    t1d->cd();
    if(!isDataBlind) t1d->SetBottomMargin(0); //RJ
    TString name_denRelUncH;
    if(!noLog) t1d->SetLogy(true);
    if(histo.Contains("redMet_rebin_shapes") ||
            histo.Contains("zpt_rebin_shapes") ) t1d->SetLogx(true);
    float maximumFound(noLog);

    if(mc &&  maximumFound<mc->GetMaximum()) maximumFound=mc->GetMaximum() * (noLog ? (systpostfix.Contains('8') ? 1.2 : 1.4) : 1.1);
    if(data && maximumFound<data->GetMaximum()) maximumFound=data->GetMaximum() * (noLog ? (systpostfix.Contains('8') ? 1.2 : 1.4) : 1.1);

    bool canvasIsFilled(false);
    if(stack && stack->GetStack() && stack->GetStack()->GetEntriesFast()>0) {
        stack->Draw("");
        TH1 *hist=(TH1*)stack->GetStack()->At(0);
        //Set YTitle, RJ
        float binsize = hist->GetBinWidth(1);
        std::ostringstream strs;
        strs << binsize;
        TString binSize = strs.str();
        if(stack->GetXaxis()) {
            stack->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
            if(histo.Contains("new_mt_shapes")) stack->GetXaxis()->SetTitle("M_{T}(ll,#slash{E}_{T}) [GeV]");
            if(isDataBlind) stack->GetYaxis()->SetTitle(hist->GetYaxis()->GetTitle());
            //stack->GetYaxis()->SetTitle("Entries");
            TString xtitle = hist->GetXaxis()->GetTitle();
            if(xtitle.Contains("GeV") )
                stack->GetYaxis()->SetTitle("Events / "+binSize+" GeV");
            else	stack->GetYaxis()->SetTitle("Events / "+binSize);

            stack->GetXaxis()->SetTitleSize(0.05);
            stack->GetYaxis()->SetTitleSize(0.05);
            stack->GetXaxis()->SetMoreLogLabels(true);
            name_denRelUncH = hist->GetXaxis()->GetTitle(); //RJ
            if(histo.Contains("redMet_rebin_shapes")) stack->GetXaxis()->SetLimits(40., stack->GetXaxis()->GetXmax());
            if(histo.Contains("zpt_rebin_shapes")) stack->GetXaxis()->SetLimits(40., stack->GetXaxis()->GetXmax());
            if(histo.Contains("redMet_shapes")) {
                stack->GetXaxis()->SetLimits(40., 400.);
            }
            if(histo.Contains("zpt_shapes")) {
                stack->GetXaxis()->SetLimits(40., 400.);
            }

            stack->SetMinimum(hist->GetMinimum());
            stack->SetMaximum(maximumFound);
            if(noLog) {
                stack->SetMaximum(1.2*maximumFound);
                //stack->GetXaxis()->SetRangeUser(hist->GetMinimum(),maximumFound);
            }
        }
        canvasIsFilled = true;
    }

    stack->GetYaxis()->SetTitleOffset(1.3);
    stack->GetYaxis()->SetLabelSize(0.04);
    stack->GetYaxis()->SetTitleSize(0.04);

    stack->GetXaxis()->SetTitleOffset(1.3);
    stack->GetXaxis()->SetLabelSize(0.04);
    stack->GetXaxis()->SetTitleSize(0.04);



    if(errors) {
        errors->Draw("2 same");
    }

    if(data) {
        //if(!isDataBlind)
        //data->Draw(canvasIsFilled ? "E1 same" : "E1"); //RJ, blind data point for Higgs-Exotic Meeting
        canvasIsFilled = true;
    }

    for(size_t ip=0; ip<spimpose.size(); ip++) {
        //double mrksz = spimpose[ip]->GetMarkerSize(); spimpose[ip]->SetMarkerSize(mrksz*0.1); spimpose[ip]->Draw((canvasIsFilled ? "e1same": "e1") ); canvasIsFilled = true; //spimpose[ip]->SetMarkerSize(mrksz);
        TString opt = "hist";
        spimpose[ip]->Draw( (opt + (canvasIsFilled ? "same": "")).Data() );
        canvasIsFilled = true;
    }

    TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    char Buffer[1024];
    //bool isSim = data ? false : true;
    double iEcm  = systpostfix.Contains('8') ? 8.0 : 7.0;
    double iLumi = systpostfix.Contains('8') ? 19577 : 5051;


    T = new TPaveText(0.18,0.9,0.33,0.8, "NDC");
    sprintf(Buffer, "#splitline{#bf{CMS}}{Preliminary}");
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.35,0.9,0.5,0.8, "NDC");
    if(finstate.Contains("mumu"))   sprintf(Buffer, "#splitline{#it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}}{#it{#mu#mu channel}}");
    if(finstate.Contains("ee"))     sprintf(Buffer, "#splitline{#it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}}{#it{ee channel}}");
    if(finstate.Contains("ll"))     sprintf(Buffer, "#splitline{#it{ZH #rightarrow l^{+}l^{-}+#slash{E}_{T}}}{#it{ll channel}}");
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.35,0.8,0.5,0.75, "NDC");
    sprintf(Buffer, "#sqrt{s} = %.1f TeV", iEcm);
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);

    T = new TPaveText(0.35,0.75,0.52,0.7, "NDC");
    sprintf(Buffer, "#scale[1.1]{#font[12]{#int}}Ldt = %.1f fb^{-1}",iLumi/1000);
    T->AddText(Buffer);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);


    legA->SetX1(0.55);
    legA->SetY1(0.65);
    legA->SetX2(0.95);
    legA->SetY2(0.92);
    legA->SetFillColor(0);
    legA->SetFillStyle(0);
    legA->SetLineColor(0);
    legA->SetHeader("");
    legA->Draw("same");
    legA->SetTextFont(42);

    std::vector<TH1 *> compDists;
    if(data)                   compDists.push_back(data);
    else if(spimpose.size()>0) compDists = spimpose;

    if(mc && compDists.size() && !isDataBlind) {
        c1d->cd();
        TPad* t2d = new TPad("t2d","t2d", 0.0, 0.0, 1.0, 0.3);
        t2d->Draw();
        t2d->cd();
        t2d->SetGridy(true);
        t2d->SetTopMargin(0);
        t2d->SetBottomMargin(0.5);
        if(histo.Contains("redMet_rebin_shapes") ||
                histo.Contains("zpt_rebin_shapes") ) t2d->SetLogx(true);

        // MC stats + syst
        TH1D *denRelUncH=0;
        if(mcPlusRelUnc) denRelUncH=(TH1D *) mcPlusRelUnc->Clone("mcrelunc");
        else             denRelUncH=(TH1D *) mc->Clone("mcrelunc");
        for(int xbin=1; xbin<=denRelUncH->GetXaxis()->GetNbins(); xbin++) {
            if(denRelUncH->GetBinContent(xbin)==0) continue;
            Double_t err=denRelUncH->GetBinError(xbin)/denRelUncH->GetBinContent(xbin);
            denRelUncH->SetBinContent(xbin,1);
            denRelUncH->SetBinError(xbin,err);
        }
        TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
        denRelUnc->SetLineColor(1);
        denRelUnc->SetFillStyle(3001);
        denRelUnc->SetFillColor(kGray);
        denRelUnc->SetMarkerColor(1);
        denRelUnc->SetMarkerStyle(1);
        denRelUncH->Reset("ICE");
        denRelUncH->Draw();
        denRelUnc->Draw("3");
        float yscale = (1.0-0.2)/(0.18-0);
        denRelUncH->GetYaxis()->SetTitle("Data/#Sigma MC");
        denRelUncH->SetMinimum(0.1);//0.4);
        denRelUncH->SetMaximum(1.9);//1.6);
        denRelUncH->GetXaxis()->SetTitle("");
        //denRelUncH->SetMinimum(0);
        //denRelUncH->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.10);
        denRelUncH->GetXaxis()->SetTitleOffset(1.3);
        denRelUncH->GetXaxis()->SetLabelSize(0.033*yscale);
        denRelUncH->GetXaxis()->SetTitleSize(0.036*yscale);
        denRelUncH->GetXaxis()->SetTickLength(0.03*yscale);
        denRelUncH->GetXaxis()->SetTitle(name_denRelUncH); //RJ
        denRelUncH->GetYaxis()->SetTitleOffset(0.3);
        denRelUncH->GetYaxis()->SetNdivisions(5);
        denRelUncH->GetYaxis()->SetLabelSize(0.033*yscale);
        denRelUncH->GetYaxis()->SetTitleSize(0.036*yscale);
        denRelUncH->GetXaxis()->SetMoreLogLabels(true);
        //denRelUncH->GetXaxis()->SetNdivisions(9);
        if(histo.Contains("redMet_rebin_shapes")) denRelUncH->GetXaxis()->SetRangeUser(40., denRelUncH->GetXaxis()->GetXmax());
        if(histo.Contains("zpt_rebin_shapes")) denRelUncH->GetXaxis()->SetRangeUser(40, denRelUncH->GetXaxis()->GetXmax());
        if(histo.Contains("redMet_shapes")) denRelUncH->GetXaxis()->SetRangeUser(40., 400.);
        if(histo.Contains("zpt_shapes")) denRelUncH->GetXaxis()->SetRangeUser(40, 400.);

        // Add comparisons
        for(size_t icd=0; icd<compDists.size(); icd++) {
            TString name("CompHistogram");
            name+=icd;
            TH1D *dataToObsH = (TH1D*)compDists[icd]->Clone(name);
            dataToObsH->Divide(mc);
            if(!isDataBlind) dataToObsH->Draw("same"); //RJ, Blind data point
        }
    } else {
        //if not comparison resize the canvas
        //c1d->SetWindowSize(600,400);
        //c1d->SetCanvasSize(600,400);
        t1d->SetPad(0,0,1,1);
    }

    c1d->Modified();
    c1d->Update();
    c1d->cd();
    //cout << "AnalysisBins: " << AnalysisBins << endl;

    string SavePath = (finstate + AnalysisBins + "_" + histo + systpostfix + ".eps").Data();
    while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
    while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
    while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
    while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
    while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
    while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
    while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
    while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
    system(string(("rm -f ") + SavePath).c_str());
    c1d->SaveAs(SavePath.c_str());
    while(SavePath.find(".eps")!=std::string::npos)SavePath.replace(SavePath.find(".eps"),4,".png");
    c1d->SaveAs(SavePath.c_str());
    while(SavePath.find(".png")!=std::string::npos)SavePath.replace(SavePath.find(".png"),4,".pdf");
    c1d->SaveAs(SavePath.c_str());
    delete c1d;
    delete T;


    /******************/
    if(showsigvsBkg) {
        TCanvas* cc = new TCanvas("cc","cc",800,600);//,800,800);
        //cc->cd();
        //TPad* t2 = new TPad("t2d","t2d", 0.0, 0.0, 1.0, 0.3);
        //t2d->Draw();
        //t2d->cd();
        //t2d->SetGridy(true);
        //t2d->SetTopMargin(0);
        //t2d->SetBottomMargin(0.5);

        // Add comparisons
        for(size_t ip=0; ip<spimpose.size(); ip++) {
            TH1D *sigTobkg = (TH1D*)spimpose[ip]->Clone();
            sigTobkg->Divide(mc);
            sigTobkg->Draw();

            float binsize = sigTobkg->GetBinWidth(1);
            std::ostringstream strs;
            strs << binsize;
            TString binSize = strs.str();

            sigTobkg->GetYaxis()->SetTitle("Sig/Bkg / "+binSize);
        }


        //TPaveText* T = new TPaveText(0.1,0.995,0.84,0.95, "NDC");
        TPaveText* T = new TPaveText(0.1,0.995,0.9,0.95, "NDC");
        T->SetFillColor(0);
        T->SetFillStyle(0);
        T->SetLineColor(0);
        T->SetTextAlign(22);
        char Buffer[1024];
        bool isSim = data ? false : true;
        double iEcm  = systpostfix.Contains('8') ? 8.0 : 7.0;
        double iLumi = systpostfix.Contains('8') ? 19577 : 5051;
        //isSim = true; //RJ
        if(isSim) sprintf(Buffer, "CMS simulation, ZH #rightarrow ll+MET, #sqrt{s}=%.1f TeV, #scale[0.7]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
        else      sprintf(Buffer, "CMS preliminary, ZH #rightarrow ll+MET, #sqrt{s}=%.1f TeV, #scale[0.7]{#int} L=%.1f fb^{-1}", iEcm, iLumi/1000);
        T->AddText(Buffer);
        T->Draw("same");
        T->SetBorderSize(0);

        //RJ adding ee or mumu channel number!
        TPaveText* T2 = new TPaveText(0.2,0.91,0.3,0.81, "NDC");
        T2->SetFillColor(0);
        T2->SetFillStyle(0);
        T2->SetLineColor(0);
        T2->SetTextAlign(22);
        char Buffer2[1024];
        if(finstate.Contains("ee"))       sprintf(Buffer2,"(ee)");
        if(finstate.Contains("mumu"))     sprintf(Buffer2,"(#mu#mu)");
        if(finstate.Contains("lleq"))	  sprintf(Buffer2,"(ee+#mu#mu)");
        T2->AddText(Buffer2);
        T2->Draw("same");

        string mysavePath_eps = (finstate + AnalysisBins + "_" + histo + systpostfix + "_sigToBkg.eps").Data();
        string mysavePath_png = (finstate + AnalysisBins + "_" + histo + systpostfix + "_sigToBkg.png").Data();
        string mysavePath_pdf = (finstate + AnalysisBins + "_" + histo + systpostfix + "_sigToBkg.pdf").Data();
        cc->SaveAs(mysavePath_eps.c_str());
        cc->SaveAs(mysavePath_png.c_str());
        cc->SaveAs(mysavePath_pdf.c_str());
        delete cc;
    }
    /******************/
}


void Draw2DHistogram(std::map<TString, TH1*>& mapbkgs, std::vector<TH1 *>& spimpose, TH1 *data, TString finstate, TString AnalysisBins)
{

    double DphiLLXaxis[4];//BalanceXaxis[4], DphiLLXaxis[4], CosThetaXaxis[4];
    //BalanceXaxis[0] = 0.8;
    //BalanceXaxis[1] = 0.9;
    //BalanceXaxis[2] = 1.1;
    //BalanceXaxis[3] = 1.2;
    DphiLLXaxis[0] = 0.;
    DphiLLXaxis[1] = 0.6;
    DphiLLXaxis[2] = 1.4;
    DphiLLXaxis[3] = TMath::Pi();
    //CosThetaXaxis[0] = -1.;
    //CosThetaXaxis[1] = -0.5;
    //CosThetaXaxis[2] = 0.;
    //CosThetaXaxis[3] = 1.;

    std::map<TString, TH1*> toDraw;
    TH1* totBkg=NULL;
    for(std::map<TString, TH1*>::iterator it=mapbkgs.begin(); it!=mapbkgs.end(); ++it) {
        toDraw[it->first]=it->second ;
        if(!totBkg) totBkg=(TH1*)it->second->Clone();
        else	totBkg->Add(it->second);
    }
    for(size_t ip=0; ip<spimpose.size(); ip++) {
        toDraw[spimpose[ip]->GetTitle()]=spimpose[ip];
    }
    toDraw["Data"]=data;
    toDraw["Background"]=totBkg;

    for(std::map<TString, TH1*>::iterator it=toDraw.begin(); it!=toDraw.end(); ++it) {
        TCanvas* c_1 = new TCanvas("c_1","c_1",800,800);
        c_1->SetTopMargin(0.08);
        c_1->SetRightMargin(0.15);
        //c_1->SetBottomMargin(0.00);
        //c_1->SetLeftMargin(0.00);
        //TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.0);
        //t1->Draw();
        //t1->cd();

        //t1->SetPadTopMargin   (0.06);
        //t1->SetPadBottomMargin(0.12);
        //t1->SetPadRightMargin (0.1);
        //t1->SetPadLeftMargin  (0.14);

        TH2F *hist = new TH2F("",";M_{T} [GeV];#Delta#phi_{ll} [rad];Events",12,0,1200,3,DphiLLXaxis);
        hist->SetTitleOffset(1.2,"yz");
        TH1* hist_id = it->second;
        TString iTitle = it->first;
        for(int ybin=1; ybin<=3; ybin++) {
            for(int xbin=1; xbin<=12; xbin++) {
                double value = hist_id->GetBinContent(xbin+12*(ybin-1));
                //printf("val: \t%f\n",value);
                hist->SetBinContent(xbin,ybin,value);
            }
        }
        gStyle->SetPalette(1);
        hist->SetMinimum(-1.0e-10);
        hist->Draw("COLZ");


        TPaveText* T = new TPaveText(0.4,0.98,0.86,0.90, "NDC");
        T->SetFillColor(0);
        T->SetFillStyle(0);
        T->SetLineColor(0);
        T->SetTextAlign(22);
        T->SetTextFont(42);
        char Buffer[1024];
        bool isSim = data ? false : true;
        double iEcm  = systpostfix.Contains('8') ? 8.0 : 7.0;
        double iLumi = systpostfix.Contains('8') ? 19577 : 5051;
        if(finstate.Contains("mumu")) sprintf(Buffer, "#it{ZH #rightarrow ll+#slash{E}_{T}}, #it{#mu#mu channel}, #sqrt{s} = %.1f TeV, L = %.1f fb^{-1}", iEcm, iLumi/1000);
        if(finstate.Contains("ee")) sprintf(Buffer, "#it{ZH #rightarrow ll+#slash{E}_{T}}, #it{ee channel}, #sqrt{s} = %.1f TeV, L = %.1f fb^{-1}", iEcm, iLumi/1000);
        T->AddText(Buffer);
        T->Draw("same");
        T->SetBorderSize(0);

        T = new TPaveText(0.58,0.995,0.86,0.95, "NDC");
        isSim = true;
        if(isSim) sprintf(Buffer, "#bf{CMS Preliminary}");
        T->AddText(Buffer);
        T->Draw("same");
        T->SetBorderSize(0);

        T = new TPaveText(0.12,0.98,0.3,0.90, "NDC");
        if(iTitle.Contains("ZH(105)"))		sprintf(Buffer, "SM M_{H}=105GeV");
        else if(iTitle.Contains("ZH(115)"))	sprintf(Buffer, "SM M_{H}=115GeV");
        else if(iTitle.Contains("ZH(125)"))	sprintf(Buffer, "SM M_{H}=125GeV");
        else if(iTitle.Contains("ZH(135)"))     sprintf(Buffer, "SM M_{H}=135GeV");
        else if(iTitle.Contains("ZH(145)"))     sprintf(Buffer, "SM M_{H}=145GeV");
        else if(iTitle.Contains("Data")) 	sprintf(Buffer, "  Data         ");
        else if(iTitle.Contains("WZ"))  	sprintf(Buffer, " WZ#rightarrow 3l#nu     ");
        else if(iTitle.Contains("ZZ")) 		sprintf(Buffer, " ZZ#rightarrow 2l2#nu    ");
        else if(iTitle.Contains("WW")) {
            T = new TPaveText(0.12,0.98,0.37,0.90, "NDC");
            sprintf(Buffer, iTitle);
        } else sprintf(Buffer, iTitle);
        T->AddText(Buffer);
        T->Draw("same");
        T->SetBorderSize(0);

        string SavePath = (finstate + AnalysisBins + "_" + histo + systpostfix + "_" + iTitle + ".eps").Data();
        while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
        while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
        while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
        while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
        while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
        while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
        while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
        while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"_");
        while(SavePath.find(" ")!=std::string::npos)SavePath.replace(SavePath.find(" "),1,"_");
        c_1->SaveAs(SavePath.c_str());
        while(SavePath.find(".eps")!=std::string::npos)SavePath.replace(SavePath.find(".eps"),4,".png");
        c_1->SaveAs(SavePath.c_str());
        while(SavePath.find(".png")!=std::string::npos)SavePath.replace(SavePath.find(".png"),4,".pdf");
        c_1->SaveAs(SavePath.c_str());
        delete c_1;
    }

}



//
void getYieldsFromShape(std::vector<TString> ch, const map<TString, Shape_t> &allShapes, TString shName, bool isdataBlinded)
{
    cout << endl;
    cout << "########################## getYieldsFromShape ##########################" << endl;

    TString massStr("");
    if(mass>0)massStr += mass;
    if(atgcpar.Length()>0) massStr = atgcpar;
    TString massStr2 = atgcpar2;

    TString MVStr("");
    if(MV>0)MVStr += MV;

    TString MAStr("");
    if(MA>0)MAStr += MA;


    YIELDS_T ee0jet_Yields;
    YIELDS_T ee1jet_Yields;
    YIELDS_T eelesq1jet_Yields;
    YIELDS_T mm0jet_Yields;
    YIELDS_T mm1jet_Yields;
    YIELDS_T mmlesq1jet_Yields;

    TH1* h;
    Double_t valerr, val;// syst;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t ich=0; ich<ch.size(); ich++) {
            YIELDS_T fortableYields;
            TString YieldsLabel = ch[ich]+AnalysisBins[b];
            cout << "Channels: " << YieldsLabel << endl;

            //nbckgs
            size_t nbckg=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.bckg.size();
            double sum_allbkgs(0.);
            double err_allbkgs(0.);

            for(size_t ibckg=0; ibckg<nbckg; ibckg++) {
                TH1* h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.bckg[ibckg];
                TString procTitle(h->GetTitle());
                cout << "backgrounds:  " << procTitle << endl;
                //if(procTitle.Contains("QCD"))continue;
                //if(procTitle.Contains("Z#rightarrow #tau#tau"))continue;

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                //syst = h->GetBinError(0)<=0 ? -1 : h->GetBinError(0);


                if(procTitle.Contains("ZZ#rightarrow 2l2#nu")) {
                    fortableYields.ZZ = val;
                    fortableYields.ZZ_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                } else if(procTitle.Contains("WZ#rightarrow 3l#nu")) {
                    fortableYields.WZ = val;
                    fortableYields.WZ_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                } else if(procTitle.Contains("Z+jets (data)")) {
                    fortableYields.Zjets = val;
                    fortableYields.Zjets_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                } else if(procTitle.Contains("Top/WW/Ztautau (data)")) {
                    fortableYields.WWtop = val;
                    fortableYields.WWtop_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                } else if(procTitle.Contains("Wjets (data)")) {
                    fortableYields.Wjets = val;
                    fortableYields.Wjets_StatErr = valerr;
                    sum_allbkgs += val;
                    err_allbkgs += valerr*valerr;
                }
            } //nbckgs END!


            //total bckg
            ////h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.totalBckg;
            ////val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
            //syst = h->GetBinError(0)<=0 ? -1 : h->GetBinError(0);

            fortableYields.totBkg = sum_allbkgs;
            fortableYields.totBkg_StatErr = sqrt(err_allbkgs);

            //signal
            size_t nsig=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal.size();
            for(size_t isig=0; isig<nsig; isig++) {
                h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal[isig];
                TString procTitle(h->GetTitle());
                procTitle.ReplaceAll("#","\\");

                if(mass>0 && !procTitle.Contains(massStr))continue;
                //MonoZ
                if(procTitle.Contains("D1")) {
                    if(mass>0 && !procTitle.Contains("D1("+massStr+"GeV)"))continue;
                }
                if(procTitle.Contains("D4")) {
                    if(mass>0 && !procTitle.Contains("D4("+massStr+"GeV)"))continue;
                }
                if(procTitle.Contains("D5")) {
                    if(mass>0 && !procTitle.Contains("D5("+massStr+"GeV)"))continue;
                }
                if(procTitle.Contains("D8")) {
                    if(mass>0 && !procTitle.Contains("D8("+massStr+"GeV)"))continue;
                }
                if(procTitle.Contains("D9")) {
                    if(mass>0 && !procTitle.Contains("D9("+massStr+"GeV)"))continue;
                }
                if(procTitle.Contains("C3")) {
                    if(mass>0 && !procTitle.Contains("C3("+massStr+"GeV)"))continue;
                }
                if(procTitle.Contains("Unpart")) {
                    if(mass>0 && !procTitle.Contains("Unpart("+massStr+")"))continue;
                }



                if(procTitle.Contains("DM") && procTitle.Contains("MV")) {
                    if(mass<0 || MV<0) continue;
                    if(!procTitle.Contains("DM("+massStr+")MV("+MVStr+")")) continue;
                }

                if(procTitle.Contains("DM") && procTitle.Contains("MA")) {
                    if(mass<0 || MA<0) continue;
                    if(!procTitle.Contains("DM("+massStr+")MA("+MAStr+")")) continue;
                }


                cout << "Signals >>>>>>> " << procTitle << endl;


                if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("ZZ"))procTitle = "ggH("+massStr+")";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("ZZ"))procTitle = "qqH("+massStr+")";
                else if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("WW"))procTitle = "ggH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("WW"))procTitle = "qqH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("ZH")                             )procTitle = "ZH("+massStr+")2lMET";

                //if(atgcpar.Length()>0 && !procTitle.Contains(massStr)) continue;
                if(atgcpar.Length()>0 && !procTitle.EndsWith(massStr)) continue;
                if(atgcpar.Length()>0) procTitle = "ZZ("+massStr+")";

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);

                fortableYields.Sig = val;
                fortableYields.Sig_StatErr = valerr;

            }// signal END!


            //data
            h=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.data;
            val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);

            fortableYields.Data = val;
            fortableYields.Data_StatErr = valerr;

            if(YieldsLabel.Contains("eeeq0jets")) 	ee0jet_Yields = fortableYields;
            if(YieldsLabel.Contains("eeeq1jets")) 	ee1jet_Yields = fortableYields;
            if(YieldsLabel.Contains("eelesq1jets")) 	eelesq1jet_Yields = fortableYields;
            if(YieldsLabel.Contains("mumueq0jets")) 	mm0jet_Yields = fortableYields;
            if(YieldsLabel.Contains("mumueq1jets")) 	mm1jet_Yields = fortableYields;
            if(YieldsLabel.Contains("mumulesq1jets")) 	mmlesq1jet_Yields = fortableYields;

        }//end channel
    }// end analysis bin

    //print out yields table

    FILE* pFile = fopen("Yields.tex","w");
    fprintf(pFile,"\n");
    fprintf(pFile,"\\begin{table}[h!]\n");
    fprintf(pFile,"\\begin{center}\n");
    fprintf(pFile,"\\begin{tabular}{c|cc|cc}\n");
    fprintf(pFile,"\\hline\\hline\n");
    fprintf(pFile,"Process & $ee$ & $\\mu\\mu$ & $ee$ & $\\mu\\mu$ \\\\\n");
    fprintf(pFile,"        & \\multicolumn{2}{c|}{0 jet selection} & \\multicolumn{2}{c}{1 jet selection} \\\\\n");
    fprintf(pFile,"\\hline\\hline\n");
    fprintf(pFile,"%40s & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f \\\\\n"
            ,"Signal"
            ,ee0jet_Yields.Sig, ee0jet_Yields.Sig_StatErr
            ,mm0jet_Yields.Sig, mm0jet_Yields.Sig_StatErr
            ,ee1jet_Yields.Sig, ee1jet_Yields.Sig_StatErr
            ,mm1jet_Yields.Sig, mm1jet_Yields.Sig_StatErr
           );

    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$Z/\\gamma^*\\rightarrow\\ell^+\\ell^-$"
            ,ee0jet_Yields.Zjets, ee0jet_Yields.Zjets_StatErr
            ,mm0jet_Yields.Zjets, mm0jet_Yields.Zjets_StatErr
            ,ee1jet_Yields.Zjets, ee1jet_Yields.Zjets_StatErr
            ,mm1jet_Yields.Zjets, mm1jet_Yields.Zjets_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$WZ\\rightarrow 3\\ell\\nu$"
            ,ee0jet_Yields.WZ, ee0jet_Yields.WZ_StatErr
            ,mm0jet_Yields.WZ, mm0jet_Yields.WZ_StatErr
            ,ee1jet_Yields.WZ, ee1jet_Yields.WZ_StatErr
            ,mm1jet_Yields.WZ, mm1jet_Yields.WZ_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$ZZ\\rightarrow 2\\ell2\\nu$"
            ,ee0jet_Yields.ZZ, ee0jet_Yields.ZZ_StatErr
            ,mm0jet_Yields.ZZ, mm0jet_Yields.ZZ_StatErr
            ,ee1jet_Yields.ZZ, ee1jet_Yields.ZZ_StatErr
            ,mm1jet_Yields.ZZ, mm1jet_Yields.ZZ_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"Top/WW/$Z\\to\\tau^+\\tau^-$"
            ,ee0jet_Yields.WWtop, ee0jet_Yields.WWtop_StatErr
            ,mm0jet_Yields.WWtop, mm0jet_Yields.WWtop_StatErr
            ,ee1jet_Yields.WWtop, ee1jet_Yields.WWtop_StatErr
            ,mm1jet_Yields.WWtop, mm1jet_Yields.WWtop_StatErr
           );


    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"W+jets"
            ,ee0jet_Yields.Wjets, ee0jet_Yields.Wjets_StatErr
            ,mm0jet_Yields.Wjets, mm0jet_Yields.Wjets_StatErr
            ,ee1jet_Yields.Wjets, ee1jet_Yields.Wjets_StatErr
            ,mm1jet_Yields.Wjets, mm1jet_Yields.Wjets_StatErr
           );
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"total bkg."
            ,ee0jet_Yields.totBkg, ee0jet_Yields.totBkg_StatErr
            ,mm0jet_Yields.totBkg, mm0jet_Yields.totBkg_StatErr
            ,ee1jet_Yields.totBkg, ee1jet_Yields.totBkg_StatErr
            ,mm1jet_Yields.totBkg, mm1jet_Yields.totBkg_StatErr
           );

    fprintf(pFile,"\\hline\n");
    if(isdataBlinded) {
        fprintf(pFile,"%40s & %s & %s & %s & %s \\\\\n","Data","-","-","-","-");
    } else {
        fprintf(pFile,"%40s & %.0f & %.0f & %.0f & %.0f \\\\\n","Data",ee0jet_Yields.Data,mm0jet_Yields.Data,ee1jet_Yields.Data,mm1jet_Yields.Data);
    }
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\hline\n");

    fprintf(pFile,"\\end{tabular}\n");
    fprintf(pFile,"\\caption{Observed(blinded) yields, pre-fit background estimates and signal predictions.}\n");
    fprintf(pFile,"\\label{tab:seltable}\n");
    fprintf(pFile,"\\end{center}\n");
    fprintf(pFile,"\\end{table}\n");
    fprintf(pFile,"\n");
    fclose(pFile);


    // printout <= 0,1 jet
    pFile = fopen("Yields_lesq1jet.tex","w");
    fprintf(pFile,"\n");
    fprintf(pFile,"\\begin{table}[h!]\n");
    fprintf(pFile,"\\begin{center}\n");
    fprintf(pFile,"\\begin{tabular}{c|c|c}\n");
    fprintf(pFile,"\\hline\\hline\n");
    fprintf(pFile,"Process & $ee$ & $\\mu\\mu$ \\\\\n");
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.3f $\\pm$ %.3f & %.3f $\\pm$ %.3f \\\\\n"
            ,"Signal"
            ,eelesq1jet_Yields.Sig, eelesq1jet_Yields.Sig_StatErr
            ,mmlesq1jet_Yields.Sig, mmlesq1jet_Yields.Sig_StatErr
           );

    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$Z/\\gamma^*\\rightarrow\\ell^+\\ell^-$"
            ,eelesq1jet_Yields.Zjets, eelesq1jet_Yields.Zjets_StatErr
            ,mmlesq1jet_Yields.Zjets, mmlesq1jet_Yields.Zjets_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$WZ\\rightarrow 3\\ell\\nu$"
            ,eelesq1jet_Yields.WZ, eelesq1jet_Yields.WZ_StatErr
            ,mmlesq1jet_Yields.WZ, mmlesq1jet_Yields.WZ_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"$ZZ\\rightarrow 2\\ell2\\nu$"
            ,eelesq1jet_Yields.ZZ, eelesq1jet_Yields.ZZ_StatErr
            ,mmlesq1jet_Yields.ZZ, mmlesq1jet_Yields.ZZ_StatErr
           );
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"Top/WW/$Z\\to\\tau^+\\tau^-$"
            ,eelesq1jet_Yields.WWtop, eelesq1jet_Yields.WWtop_StatErr
            ,mmlesq1jet_Yields.WWtop, mmlesq1jet_Yields.WWtop_StatErr
           );


    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"W+jets"
            ,eelesq1jet_Yields.Wjets, eelesq1jet_Yields.Wjets_StatErr
            ,mmlesq1jet_Yields.Wjets, mmlesq1jet_Yields.Wjets_StatErr
           );
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"%40s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f \\\\\n"
            ,"total bkg."
            ,eelesq1jet_Yields.totBkg, eelesq1jet_Yields.totBkg_StatErr
            ,mmlesq1jet_Yields.totBkg, mmlesq1jet_Yields.totBkg_StatErr
           );

    fprintf(pFile,"\\hline\n");
    if(isdataBlinded) {
        fprintf(pFile,"%40s & %s & %s \\\\\n","Data","-","-");
    } else {
        fprintf(pFile,"%40s & %.0f & %.0f \\\\\n","Data",eelesq1jet_Yields.Data,mmlesq1jet_Yields.Data);
    }
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\hline\n");

    fprintf(pFile,"\\end{tabular}\n");
    fprintf(pFile,"\\caption{Observed(blinded) yields, pre-fit background estimates and signal predictions.}\n");
    fprintf(pFile,"\\label{tab:seltable}\n");
    fprintf(pFile,"\\end{center}\n");
    fprintf(pFile,"\\end{table}\n");
    fprintf(pFile,"\n");
    fclose(pFile);



    cout << "########################## getYieldsFromShape ##########################" << endl;
    cout << endl;
}


void getEffFromShape(std::vector<TString> ch, map<TString, Shape_t> &allShapes, TString shName)
{
    FILE* pFile = fopen("Efficiency.tex","w");
//  fprintf(pFile,"\\begin{sidewaystable}[htp]\n\\begin{center}\n\\caption{Event yields expected for background and signal processes and observed in data.}\n\\label{tab:table}\n");

    string Ccol   = "\\begin{tabular}{|c|";
    string Cname  = "channel";
    string Cval   = "";

    TString massStr("");
    if(mass>0)massStr += mass;
    if(atgcpar.Length()>0) massStr = atgcpar;
    TString massStr2 = atgcpar2;


    TH1* h;
    Double_t valerr, val;//syst;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t ich=0; ich<ch.size(); ich++) {
            TString icol(ch[ich]+"-"+AnalysisBins[b]);
            icol.ReplaceAll("mu","\\mu");
            icol.ReplaceAll("_"," ");
            Cval = "$ "+icol+" $";

            //signal
            size_t nsig=allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second.signal.size();
            for(size_t isig=0; isig<nsig; isig++) {
                Shape_t& shape  = allShapes.find(ch[ich]+AnalysisBins[b]+shName)->second;
                h=shape.signal[isig];
                TString procTitle(h->GetTitle());
                procTitle.ReplaceAll("#","\\");

                if(mass>0 && !procTitle.Contains(massStr))continue;
                if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("ZZ"))procTitle = "ggH("+massStr+")";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("ZZ"))procTitle = "qqH("+massStr+")";
                else if(mass>0 && procTitle.Contains("ggH") && procTitle.Contains("WW"))procTitle = "ggH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("qqH") && procTitle.Contains("WW"))procTitle = "qqH("+massStr+")WW";
                else if(mass>0 && procTitle.Contains("ZH")                             )procTitle = "ZH("+massStr+")2lMET";

                //if(atgcpar.Length()>0 && !procTitle.Contains(massStr)) continue;
                if(atgcpar.Length()>0 && !procTitle.EndsWith(massStr)) continue;
                if(atgcpar.Length()>0) procTitle = "ZZ("+massStr+")";

                if(b==0&&ich==0)Ccol  += "c|";
                if(b==0&&ich==0)Cname += "&$" + procTitle+"$";

//       printf("%s --> Xsec=%E x %E\n",h->GetTitle(), shape.xsections[h->GetTitle()], shape.BRs[h->GetTitle()]);
                double xsecXbr = shape.xsections[h->GetTitle()] * shape.BRs[h->GetTitle()];

                val = h->IntegralAndError(1,h->GetXaxis()->GetNbins(),valerr);
                if(val<1E-6) {
                    val=0.0;
                    valerr=0.0;
                }
                Cval += "&" + toLatexRounded(val/xsecXbr,valerr/xsecXbr);

                fprintf(pFile,"%30s %30s %4.0f %6.2E %6.2E %6.2E %6.2E\n",icol.Data(), procTitle.Data(), (double)mass, shape.xsections[h->GetTitle()], shape.BRs[h->GetTitle()], val/xsecXbr, valerr/xsecXbr);
            }

            //endline
//    if(b==0&&ich==0)fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
//    if(b==0&&ich==0)fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
//    fprintf(pFile,"%s\\\\\n", Cval.c_str());
        }
    }
//  fprintf(pFile,"\\hline\n");
//  fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{sidewaystable}\n");
//  fprintf(pFile,"\n\n\n\n");
    fclose(pFile);
}



std::vector<TString>  buildDataCard(Int_t mass, TString histo, TString url, TString Json)
{
    return buildDataCard("", mass, histo, url, Json);
}

std::vector<TString>  buildDataCard(TString atgcpar, Int_t mass, TString histo, TString url, TString Json)
{
    std::vector<TString> dcUrls;

    //get the datacard inputs
    DataCardInputs dci = convertHistosForLimits(atgcpar,mass,histo,url,Json);

    TString eecard = "";
    TString mumucard = "";
    TString combinedcard = "";
    TString combinedcardLL = "";

    //build the datacard separately for each channel
    for(size_t i=1; i<=dci.ch.size(); i++) {
        TString dcName=dci.shapesFile;
        dcName.ReplaceAll(".root","_"+dci.ch[i-1]+".dat");
        FILE* pFile = fopen(dcName.Data(),"w");

        if(!dci.ch[i-1].Contains("lleq")) combinedcard += dci.ch[i-1]+"="+dcName+" ";
        if(dci.ch[i-1].Contains("lleq")) combinedcardLL += dci.ch[i-1]+"="+dcName+" ";
        if(dci.ch[i-1].Contains("ee"))eecard += dci.ch[i-1]+"="+dcName+" ";
        if(dci.ch[i-1].Contains("mumu"))mumucard += dci.ch[i-1]+"="+dcName+" ";

        //header
        fprintf(pFile, "imax 1\n");
        fprintf(pFile, "jmax *\n");
        fprintf(pFile, "kmax *\n");
        fprintf(pFile, "-------------------------------\n");
        if(shape) {
            fprintf(pFile, "shapes * * %s %s/$PROCESS %s/$PROCESS_$SYSTEMATIC\n",dci.shapesFile.Data(), dci.ch[i-1].Data(), dci.ch[i-1].Data());
            fprintf(pFile, "-------------------------------\n");
        }
        //observations
        fprintf(pFile, "bin 1\n");
        fprintf(pFile, "Observation %f\n",dci.obs[RateKey_t("obs",dci.ch[i-1])]);
        fprintf(pFile, "-------------------------------\n");

        fprintf(pFile,"%55s ", "bin");
        for(size_t j=1; j<=dci.procs.size(); j++) {
            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            fprintf(pFile,"%6i ", 1);
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%55s ", "process");
        for(size_t j=1; j<=dci.procs.size(); j++) {
            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            fprintf(pFile,"%6s ", (postfix+dci.procs[j-1]).Data());
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%55s ", "process");
        int procCtr(1-dci.nsignalproc);
        for(size_t j=1; j<=dci.procs.size(); j++) {
            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            fprintf(pFile,"%6i ", procCtr );
            procCtr++;
        }
        fprintf(pFile,"\n");

        fprintf(pFile,"%55s ", "rate");
        for(size_t j=1; j<=dci.procs.size(); j++) {
            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            fprintf(pFile,"%6f ", dci.rates[RateKey_t(dci.procs[j-1],dci.ch[i-1])] );
            //fprintf(pFile,"%e ", dci.rates[RateKey_t(dci.procs[j-1],dci.ch[i-1])] );
        }
        fprintf(pFile,"\n");
        fprintf(pFile, "-------------------------------\n");


        //systematics
        char sFile[2048];
        bool isSyst;
        if(runSystematics) {
            if(systpostfix.Contains("13")) {
                fprintf(pFile,"%45s %10s ", "lumi_13TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("EM") && !dci.procs[j-1].Contains("Zjets") && !dci.procs[j-1].Contains("Wjets")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["lumi_13TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } else if(systpostfix.Contains('8')) {
                fprintf(pFile,"%45s %10s ", "lumi_8TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("EM") && !dci.procs[j-1].Contains("Zjets") && !dci.procs[j-1].Contains("Wjets")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["lumi_8TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } else {
                fprintf(pFile,"%45s %10s ", "lumi_7TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("EM") && !dci.procs[j-1].Contains("Zjets")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["lumi_7TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }/*
            if(systpostfix.Contains('8')) {
                fprintf(pFile,"%45s %10s ", "accept_8TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("EM") && !dci.procs[j-1].Contains("Zjets")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["accept_8TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } else {
                fprintf(pFile,"%45s %10s ", "accept_7TeV", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("EM") && !dci.procs[j-1].Contains("Zjets")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["accept_7TeV"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } */
            if(sysSherpa != 1.) {
                fprintf(pFile,"%45s %10s ", "sherpa_kin_syst", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("zz2l2nu") && dci.procs[j-1].Contains("_") ) {
                        fprintf(pFile,"%6f ",sysSherpa);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }

//         if(mass>0){
//         fprintf(pFile,"%45s %10s ", "theoryUncXS_HighMH", "lnN");
//         for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
//            if((int)j<=dci.nsignalproc){fprintf(pFile,"%6f ",std::min(1.0+1.5*pow((mass/1000.0),3),2.0));}else{fprintf(pFile,"%6s ","-");}
//         }fprintf(pFile,"\n");
//         }

            //Id+Trigger efficiencies combined
            if(dci.ch[i-1].Contains("ee")) {
                fprintf(pFile,"%45s %10s ", "CMS_eff_e", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("EM") && !dci.procs[j-1].Contains("Zjets") && !dci.procs[j-1].Contains("Wjets")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_eff_e"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            } else {
                fprintf(pFile,"%45s %10s ", "CMS_eff_m", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(!dci.procs[j-1].Contains("EM") && !dci.procs[j-1].Contains("Zjets") && !dci.procs[j-1].Contains("Wjets")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_eff_m"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }


//##############################
//# temporily take values from Guillelmo
            /*
                        if(systpostfix.Contains('8')) {
                            fprintf(pFile,"%45s %10s ", "CMS_zllwimps_EM_8TeV", "lnN");
                            for(size_t j=1; j<=dci.procs.size(); j++) {
                                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                                if(dci.procs[j-1].Contains("EM")) {
            			if(dci.ch[i-1].Contains("eeeq0jets")) fprintf(pFile,"%6f ",2.109);
            			else if(dci.ch[i-1].Contains("eeeq1jets")) fprintf(pFile,"%6f ",2.459);
            			else if(dci.ch[i-1].Contains("mumueq0jets")) fprintf(pFile,"%6f ",2.018);
            			else if(dci.ch[i-1].Contains("mumueq1jets")) fprintf(pFile,"%6f ",1.73);
            			else fprintf(pFile,"%6s ","-");
                                } else {
                                    fprintf(pFile,"%6s ","-");
                                }
                            }
                            fprintf(pFile,"\n");
                        } else {
                            fprintf(pFile,"%45s %10s ", "CMS_zllwimps_EM_7TeV", "lnN");
                            for(size_t j=1; j<=dci.procs.size(); j++) {
                                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                                if(dci.procs[j-1].Contains("EM")) {
                                    if(dci.ch[i-1].Contains("eeeq0jets")) fprintf(pFile,"%6f ",4.16);
                                    else if(dci.ch[i-1].Contains("eeeq1jets")) fprintf(pFile,"%6f ",2.018);
                                    else if(dci.ch[i-1].Contains("mumueq0jets")) fprintf(pFile,"%6f ",3.453);
                                    else if(dci.ch[i-1].Contains("mumueq1jets")) fprintf(pFile,"%6f ",2.328);
            			else fprintf(pFile,"%6s ","-");
                                } else {
                                    fprintf(pFile,"%6s ","-");
                                }
                            }
                            fprintf(pFile,"\n");
                        }
            */
//##############################
            /*

            	    //k-factor  Nov 6,2014
            	    fprintf(pFile,"%45s %10s ", "CMS_zllwimps_WZ_kfactor", "lnN");
                            for(size_t j=1; j<=dci.procs.size(); j++) {
                                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                                if(dci.procs[j-1].Contains("WZ")) {
                                    fprintf(pFile,"%6f ",1.25);
                                } else {
                                    fprintf(pFile,"%6s ","-");
                                }
                            }
                         fprintf(pFile,"\n");

            	    fprintf(pFile,"%45s %10s ", "CMS_zllwimps_ZZ_kfactor", "lnN");
                            for(size_t j=1; j<=dci.procs.size(); j++) {
                                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                                if(dci.procs[j-1].Contains("ZZ")) {
                                    fprintf(pFile,"%6f ",1.25);
                                } else {
                                    fprintf(pFile,"%6s ","-");
                                }
                            }
                         fprintf(pFile,"\n");

            */



            fprintf(pFile,"%45s %10s ", "CMS_zllwimps_ZZ_generator", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                if(dci.procs[j-1].Contains("ZZ")) {
                    if(DYMET_EXTRAPOL<=80)        fprintf(pFile,"%6f ",1.1424);
                    else if (DYMET_EXTRAPOL==90 ) fprintf(pFile,"%6f ",1.1437);
                    else if (DYMET_EXTRAPOL==100) fprintf(pFile,"%6f ",1.1507);
                    else if (DYMET_EXTRAPOL==110) fprintf(pFile,"%6f ",1.1644);
                    else if (DYMET_EXTRAPOL==120) fprintf(pFile,"%6f ",1.1636);
                    else if (DYMET_EXTRAPOL==130) fprintf(pFile,"%6f ",1.1809);
                    else if (DYMET_EXTRAPOL==140) fprintf(pFile,"%6f ",1.1866);
                    else if (DYMET_EXTRAPOL==150) fprintf(pFile,"%6f ",1.1944);
                    else if (DYMET_EXTRAPOL==170) fprintf(pFile,"%6f ",1.2331);
                    else if (DYMET_EXTRAPOL==190) fprintf(pFile,"%6f ",1.2449);
                    else if (DYMET_EXTRAPOL==210) fprintf(pFile,"%6f ",1.2573);
                    else if (DYMET_EXTRAPOL==230) fprintf(pFile,"%6f ",1.3149);
                    else if (DYMET_EXTRAPOL==250) fprintf(pFile,"%6f ",1.3596);
                    else if (DYMET_EXTRAPOL>=260) fprintf(pFile,"%6f ",1.3797);
                    else fprintf(pFile,"%6s ","-");
                } else {
                    fprintf(pFile,"%6s ","-");
                }
            }
            fprintf(pFile,"\n");




            //Oct 30, 2013

            if(dci.ch[i-1].Contains("mumueq0jets")) {
                fprintf(pFile,"%45s %10s ", "CMS_zllwimps_WZ3l", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("WZ")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_zllwimps_mumueq0jets_leptonVeto"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }

            if(dci.ch[i-1].Contains("mumueq1jets") || dci.ch[i-1].Contains("mumulesq1jets")) {
                fprintf(pFile,"%45s %10s ", "CMS_zllwimps_WZ3l", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("WZ")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_zllwimps_mumueq1jets_leptonVeto"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }




            if(dci.ch[i-1].Contains("eeeq0jets")) {
                fprintf(pFile,"%45s %10s ", "CMS_zllwimps_WZ3l", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("WZ")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_zllwimps_eeeq0jets_leptonVeto"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }

            if(dci.ch[i-1].Contains("eeeq1jets") || dci.ch[i-1].Contains("eelesq1jets")) {
                fprintf(pFile,"%45s %10s ", "CMS_zllwimps_WZ3l", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].Contains("WZ")) {
                        fprintf(pFile,"%6f ",1.0+normSysts["CMS_zllwimps_eeeq1jets_leptonVeto"]);
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            }



            // Scale-e and scale-mu replaced by shape-unncertainty "les"
//          if(dci.ch[i-1].Contains("ee")){
//             fprintf(pFile,"%45s %10s ", "CMS_scale_e", "lnN");
//             for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
//                if(!dci.procs[j-1].Contains("data")){fprintf(pFile,"%6f ",1.0+normSysts["CMS_scale_e"]);}else{fprintf(pFile,"%6s ","-");}
//             }fprintf(pFile,"\n");
//          }else{
//             fprintf(pFile,"%45s %10s ", "CMS_scale_m", "lnN");
//             for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
//                if(!dci.procs[j-1].Contains("data")){fprintf(pFile,"%6f ",1.0+normSysts["CMS_scale_m"]);}else{fprintf(pFile,"%6s ","-");}
//             }fprintf(pFile,"\n");
//          }

            /*
                        if(mass>0) {

                            fprintf(pFile,"%45s %10s ", "UEPS", "lnN");
                            for(size_t j=1; j<=dci.procs.size(); j++) {
                                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                                if(dci.procs[j-1].Contains("ZH")) {
                                    fprintf(pFile,"%6f ",1.03);
                                }else {
                                    fprintf(pFile,"%6s ","-");
                                }
                            }
                            fprintf(pFile,"\n");

                        } //mass>0
            */

            ///////////////////////////////////////////////
            // RJ, for count and cut
            ///////////////////////////////////////////////

            fprintf(pFile,"%45s %10s ", "QCDscale_VV", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                if(dci.procs[j-1].BeginsWith("ZZ")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.067);
                    } else {
                        fprintf(pFile,"%6f ",1.067);
                    }
                } else if(dci.procs[j-1].BeginsWith("WZ")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.077);
                    } else {
                        fprintf(pFile,"%6f ",1.077);
                    }
                } else {
                    fprintf(pFile,"%6s ","-");
                }
            }
            fprintf(pFile,"\n");

            /*
                        fprintf(pFile,"%45s %10s ", "QCDscale_DM", "lnN");
                        for(size_t j=1; j<=dci.procs.size(); j++) {
                            if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                            if(dci.procs[j-1].BeginsWith("D5") || dci.procs[j-1].BeginsWith("D8")
            			|| dci.procs[j-1].BeginsWith("C3") || dci.procs[j-1].BeginsWith("D9") ) {
            		    cout << "QCDscale_DM " << dci.procs[j-1] << " mass: " << mass << endl;
                                fprintf(pFile,"%6f ",1.054);
                            } else {
                                fprintf(pFile,"%6s ","-");
                            }
                        }
                        fprintf(pFile,"\n");
            */

            fprintf(pFile,"%45s %10s ", "QCDscale_DM", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                if(dci.procs[j-1].BeginsWith("D5")) {
                    if(mass==1) 	fprintf(pFile,"%6f ",1.0520);
                    if(mass==10) 	fprintf(pFile,"%6f ",1.0525);
                    if(mass==100) 	fprintf(pFile,"%6f ",1.0532);
                    if(mass==200) 	fprintf(pFile,"%6f ",1.0537);
                    if(mass==300) 	fprintf(pFile,"%6f ",1.0543);
                    if(mass==500) 	fprintf(pFile,"%6f ",1.0552);
                    if(mass==1000) 	fprintf(pFile,"%6f ",1.0570);
                } else if(dci.procs[j-1].BeginsWith("D8")) {
                    if(mass==1)         fprintf(pFile,"%6f ",1.0525);
                    if(mass==10)        fprintf(pFile,"%6f ",1.0527);
                    if(mass==100)       fprintf(pFile,"%6f ",1.0528);
                    if(mass==200)       fprintf(pFile,"%6f ",1.0543);
                    if(mass==300)       fprintf(pFile,"%6f ",1.0546);
                    if(mass==500)       fprintf(pFile,"%6f ",1.0557);
                    if(mass==1000)      fprintf(pFile,"%6f ",1.0564);
                } else if(dci.procs[j-1].BeginsWith("D9")) {
                    if(mass==1)         fprintf(pFile,"%6f ",1.0551);
                    if(mass==10)        fprintf(pFile,"%6f ",1.0555);
                    if(mass==100)       fprintf(pFile,"%6f ",1.0558);
                    if(mass==200)       fprintf(pFile,"%6f ",1.0559);
                    if(mass==300)       fprintf(pFile,"%6f ",1.0567);
                    if(mass==500)       fprintf(pFile,"%6f ",1.0573);
                    if(mass==1000)      fprintf(pFile,"%6f ",1.0572);
                } else if(dci.procs[j-1].BeginsWith("C3")) {
                    if(mass==1)         fprintf(pFile,"%6f ",1.0516);
                    if(mass==10)        fprintf(pFile,"%6f ",1.0526);
                    if(mass==100)       fprintf(pFile,"%6f ",1.0535);
                    if(mass==200)       fprintf(pFile,"%6f ",1.0542);
                    if(mass==300)       fprintf(pFile,"%6f ",1.0552);
                    if(mass==500)       fprintf(pFile,"%6f ",1.0561);
                    if(mass==1000)      fprintf(pFile,"%6f ",1.0561);
                } else {
                    fprintf(pFile,"%6s ","-");
                }
            }
            fprintf(pFile,"\n");





            fprintf(pFile,"%45s %10s ", "QCDscale_UnPart", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                if     (dci.procs[j-1].BeginsWith("UnPart1p01") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p01"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p02") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p02"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p04") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p04"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p06") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p06"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p09") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p09"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p10") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p10"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p20") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p20"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p30") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p30"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p40") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p40"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p50") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p50"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p60") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p60"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p70") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p70"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p80") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p80"]);
                else if(dci.procs[j-1].BeginsWith("UnPart1p90") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart1p90"]);
                else if(dci.procs[j-1].BeginsWith("UnPart2p00") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart2p00"]);
                else if(dci.procs[j-1].BeginsWith("UnPart2p20") ) fprintf(pFile,"%6f ",normSysts["QCDscale_UnPart2p20"]);
                else {
                    fprintf(pFile,"%6s ","-");
                }
            }
            fprintf(pFile,"\n");





            /*
                fprintf(pFile,"%45s %10s ", "QCDscale_VH", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++) {
                    if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                    if(dci.procs[j-1].BeginsWith("ZH")) {
                        if(dci.ch[i-1].Contains("eq0jets")) {
                            if(systpostfix.Contains('8')) {
                                fprintf(pFile,"%6f ",1.072);
                            } else {
                                fprintf(pFile,"%6f ",1.072);
                            }
                        } else if(dci.ch[i-1].Contains("eq1jets")) {
                            if(systpostfix.Contains('8')) {
                                fprintf(pFile,"%6f ",1.105);
                            } else {
                                fprintf(pFile,"%6f ",1.105);
                            }
                        }
                    } else {
                        fprintf(pFile,"%6s ","-");
                    }
                }
                fprintf(pFile,"\n");
            */

            /*
            	 if(dci.ch[i-1].Contains("ee"))
                     	fprintf(pFile,"%45s %10s ", "scale_e", "lnN");
            	 else fprintf(pFile,"%45s %10s ", "scale_mu", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++) { if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            	          if(dci.procs[j-1].BeginsWith("ZH")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.01);}else{fprintf(pFile,"%6f ",1.01);}}
            		  else if(dci.procs[j-1].BeginsWith("zz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.01);}else{fprintf(pFile,"%6f ",1.01);}}
            		  else if(dci.procs[j-1].BeginsWith("zll") && !dci.procs[j-1].BeginsWith("zlldata")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.01);}else{fprintf(pFile,"%6f ",1.01);}}
            		  else if(dci.procs[j-1].BeginsWith("wz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.01);}else{fprintf(pFile,"%6f ",1.01);}}
            		  else {fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");
            */

//RJ May 17
            /*
                     fprintf(pFile,"%45s %10s ", "jes", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++) { if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            	          if(dci.procs[j-1].BeginsWith("ZH")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.02);}else{fprintf(pFile,"%6f ",1.02);}}
            		  else if(dci.procs[j-1].BeginsWith("zz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.02);}else{fprintf(pFile,"%6f ",1.02);}}
            		  else if(dci.procs[j-1].BeginsWith("wz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.03);}else{fprintf(pFile,"%6f ",1.03);}}
            		  else {fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");
            */

            /*
            	//remove it if use data-driven DY sample
                     fprintf(pFile,"%45s %10s ", "err_zll", "lnN");
                     for(size_t j=1; j<=dci.procs.size(); j++) { if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
            	          if(dci.procs[j-1].BeginsWith("zll")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",2.0);}else{fprintf(pFile,"%6f ",2.0);}}
            		  else {fprintf(pFile,"%6s ","-");}
                     }fprintf(pFile,"\n");
            */
            ///////////////////////////////////////////////


            ///////////////////////////////////////////////
            /*
                fprintf(pFile,"%45s %10s ", "XSec_sys_WW", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                   if(dci.procs[j-1].BeginsWith("ww")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.097);}else{fprintf(pFile,"%6f ",1.097);}
                   }else{fprintf(pFile,"%6s ","-");}
                }fprintf(pFile,"\n");

                fprintf(pFile,"%45s %10s ", "XSec_sys_WZ", "lnN");
                for(size_t j=1; j<=dci.procs.size(); j++){ if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                   if(dci.procs[j-1].BeginsWith("wz")){if(systpostfix.Contains('8')){fprintf(pFile,"%6f ",1.056);}else{fprintf(pFile,"%6f ",1.056);}
                   }else{fprintf(pFile,"%6s ","-");}
                }fprintf(pFile,"\n");
            */

            ///////////////////////////////////////////////


            for(std::map<TString, std::map<RateKey_t,Double_t> >::iterator it=dci.systs.begin(); it!=dci.systs.end(); it++) {
                if(!runSystematics && string(it->first.Data()).find("stat")>0 )continue;

                isSyst=false;
                if(it->first.Contains("_sys_") || it->first.Contains("_interpol_")) {
                    sprintf(sFile,"%45s %10s ", it->first.Data(), "lnN");
                    for(size_t j=1; j<=dci.procs.size(); j++) {
                        if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                        if(it->second.find(RateKey_t(dci.procs[j-1],dci.ch[i-1])) != it->second.end()) {
                            Double_t systUnc = it->second[RateKey_t(dci.procs[j-1],dci.ch[i-1])];
                            if(systUnc<=0) {
                                sprintf(sFile,"%s%6s ",sFile,"-");
                            } else {
                                sprintf(sFile,"%s%6f ",sFile,(1.0+ (systUnc / dci.rates[RateKey_t(dci.procs[j-1],dci.ch[i-1])]) ));
                                isSyst=true;
                            }
                        } else {
                            sprintf(sFile,"%s%6s ",sFile,"-");
                        }
                    }
                    if(isSyst)fprintf(pFile,"%s\n",sFile);

                } else {
                    if(shape) {
                        if(it->first.Contains("CMS_zllwimps_les") && dci.ch[i-1].Contains("ee")) {
                            sprintf(sFile,"%45s %10s ", "CMS_p_scale_e", "shape");
                        } else
                            sprintf(sFile,"%45s %10s ", it->first.Data(), "shape");
                    } else {
                        sprintf(sFile,"%45s %10s ", it->first.Data(), "lnN");
                    }
                    for(size_t j=1; j<=dci.procs.size(); j++) {
                        if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                        if(it->first.Contains("sherpa")) continue; //RJ
                        //if(it->first.Contains("CMS_scale_j")) continue; //RJ
                        //if(it->first.Contains("CMS_res_j")) continue; //RJ
                        if(it->first.Contains("CMS_zllwimps_les") && dci.ch[i-1].Contains("mumu")) continue; //RJ

                        if(it->second.find(RateKey_t(dci.procs[j-1],dci.ch[i-1])) != it->second.end()) {
                            if(it->first.Contains("sherpa") && (!dci.procs[j-1].Contains("zz2l2nu"))) {
                                sprintf(sFile,"%s%6s ",sFile,"-");
                            } else {
                                if(it->first.Contains("CMS_zllwimps_pdf")) {
                                    normSysts["CMS_zllwimps_pdf"]=it->second[RateKey_t(dci.procs[j-1],dci.ch[i-1])];
                                    continue;
                                }
                                sprintf(sFile,"%s%6.3f ",sFile,it->second[RateKey_t(dci.procs[j-1],dci.ch[i-1])]);
                                isSyst=true;
                            }
                        } else {
                            sprintf(sFile,"%s%6s ",sFile,"-");
                        }
                    }
                    if(isSyst)fprintf(pFile,"%s\n",sFile);
                }
            }//end of shape uncertainty




            fprintf(pFile,"%45s %10s ", "pdf_qqbar", "lnN");
            for(size_t j=1; j<=dci.procs.size(); j++) {
                //cout << "pdf_qqbar: dci.procs[j-1] " << dci.procs[j-1] << endl;
                if(dci.rates.find(RateKey_t(dci.procs[j-1],dci.ch[i-1]))==dci.rates.end()) continue;
                if(dci.procs[j-1].BeginsWith("ZH")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.055);
                    } else {
                        fprintf(pFile,"%6f ",1.055);
                    }
                } else if(dci.procs[j-1].BeginsWith("ZZ")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.057);
                    } else {
                        fprintf(pFile,"%6f ",1.057);
                    }
                } else if(dci.procs[j-1].BeginsWith("WZ")) {
                    if(systpostfix.Contains('8')) {
                        fprintf(pFile,"%6f ",1.048);
                    } else {
                        fprintf(pFile,"%6f ",1.048);
                    }
                } else if(dci.procs[j-1].BeginsWith("D5") || dci.procs[j-1].BeginsWith("D8")
                          || dci.procs[j-1].BeginsWith("C3") || dci.procs[j-1].BeginsWith("D9")
                          || dci.procs[j-1].BeginsWith("UnPart")
                          || dci.procs[j-1].BeginsWith("dm") ) {
                    fprintf(pFile,"%6f ",normSysts["CMS_zllwimps_pdf"]); //convert shape -> normalization
                } else {
                    fprintf(pFile,"%6s ","-");
                }
            }
            fprintf(pFile,"\n");





        }


        fclose(pFile);
        cout << "Data card for " << dci.shapesFile << " and " << dci.ch[i-1] << " channel @ " << dcName << endl;
        dcUrls.push_back(dcName);
    }

    FILE* pFile = fopen("combineCards.sh","w");
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + combinedcard + " > " + "card_combined.dat").Data());
    //fprintf(pFile,"%s;\n",(TString("combineCards.py ") + combinedcardLL + " > " + "card_combinedLL.dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + eecard       + " > " + "card_ee.dat").Data());
    fprintf(pFile,"%s;\n",(TString("combineCards.py ") + mumucard     + " > " + "card_mumu.dat").Data());
    fclose(pFile);

//  gSystem->Exec(TString("combineCards.py ") + combinedcard + " > " + outDir+"/card_combined.dat");
//  gSystem->Exec(TString("combineCards.py ") + eecard       + " > " + outDir+"/card_ee.dat");
//  gSystem->Exec(TString("combineCards.py ") + mumucard     + " > " + outDir+"/card_mumu.dat");

    //all done
    return dcUrls;
}

//
DataCardInputs convertHistosForLimits(Int_t mass,TString histo,TString url,TString Json)
{
    return convertHistosForLimits("",mass,histo,url,Json);
}

DataCardInputs convertHistosForLimits(TString atgcpar,Int_t mass,TString histo,TString url,TString Json)
{
    DataCardInputs dci;

    //init the json wrapper
    JSONWrapper::Object Root(Json.Data(), true);

    //init globalVariables
    TString massStr("");
    if(mass>0)massStr += mass;
    if(atgcpar.Length()>0) massStr = atgcpar;
    TString massStr2 = atgcpar2;
    if(massStr2.CompareTo("")==0) massStr2 = massStr;
    //std::set<TString> allCh,allProcs;
    std::vector<TString> allCh,allProcs;

    TString MVStr("");
    if(MV>0)MVStr += MV;

    TString MAStr("");
    if(MA>0)MAStr += MA;

    //open input file
    TFile* inF = TFile::Open(url);
    if( !inF || inF->IsZombie() ) {
        cout << "Invalid file name : " << url << endl;
        return dci;
    }

    //get the shapes for each channel
    map<TString, Shape_t> allShapes;
    map<TString, Shape_t> allShapesL;
    map<TString, Shape_t> allShapesR;
    TString ch[]= {"mumu","ee","emu","ll"};
    const size_t nch=sizeof(ch)/sizeof(TString);
    std::vector<TString> sh;
    sh.push_back(histo);
    //if(subNRB2011 || subNRB2012)sh.push_back("nonresbckg_ctrl");
    if(subNRB2011 || subNRB2012)sh.push_back(histo+"_NRBctrl");
    if(subNRB2012)		sh.push_back(histo+"BTagSB");
    if(subWZ)			sh.push_back(histo+"_3rdLepton");
    sh.push_back("FR_WjetCtrl_"+histo);
    sh.push_back(histo+"_NRBsyst");
    sh.push_back("FR_QCDCtrl_mt_shapes");
    sh.push_back("pfmet_minus_shapes");
    sh.push_back("balancedif_minus_shapes");
    sh.push_back("dphizmet_minus_shapes");
    const size_t nsh=sh.size();
    for(size_t i=0; i<nch; i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            int indexcut_ = indexcut;
            double cutMin=shapeMin;
            double cutMax=shapeMax;
            if(indexvbf>=0 && AnalysisBins[b].Contains("vbf")) {
                printf("use vbf index(%i) for bin %s\n", indexvbf, AnalysisBins[b].Data());
                indexcut_ = indexvbf;
                cutMin=shapeMinVBF;
                cutMax=shapeMaxVBF;
            }
            for(size_t j=0; j<nsh; j++) {
                //printf("i=%i b=%i j=%i\n",(int)i,(int)b,(int)j);
                allShapes[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexcut_,Root,cutMin, cutMax);
                cout << "Loading shapes: " << ch[i]+AnalysisBins[b]+sh[j] << endl;
                if(indexcutL>=0 && indexcutR>=0) {
                    if(indexvbf>=0 && AnalysisBins[b].Contains("vbf")) {
                        allShapesL[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexvbf,Root,cutMin, cutMax);
                        allShapesR[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexvbf,Root,cutMin, cutMax);
                    } else {
                        allShapesL[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexcutL,Root,cutMin, cutMax);
                        allShapesR[ch[i]+AnalysisBins[b]+sh[j]]=getShapeFromFile(inF, ch[i]+AnalysisBins[b],sh[j],indexcutR,Root,cutMin, cutMax);
                    }
                }
            }
        }
    }

    //all done with input file
    inF->Close();

    //printf("done loading all shapes\n");

    //define vector for search
    std::vector<TString>& selCh = Channels;
    //selCh.push_back("ee"); selCh.push_back("mumu");

    //non-resonant background estimation
    //estimateNonResonantBackground(selCh,"emu",allShapes,"nonresbckg_ctrl");
    /*
        // estimate jet induced background
        //>>>>>>>>>>>>>>>>>>>>
        doWjetsBackground(selCh,allShapes,histo); // W+jets -> Wjets (data)
        //<<<<<<<<<<<<<<<<<<<<
        //doQCDBackground(selCh,allShapes,histo);
    */


    // Mono-Z analysis (new)   Top/WW/Ztautau (data)
    //MC closure test, adding systematics
    dodataDrivenWWtW(selCh,"emu",allShapes,histo,true);
    //new data-driven WW/tW/Ztautau background
    dodataDrivenWWtW(selCh,"emu",allShapes,histo,false);

    /*

        //remove the non-resonant background from data
        //if(subNRB2011 || subNRB2012) doBackgroundSubtraction(selCh,"emu",allShapes,histo,histo+"_NRBctrl", url, Root, false);
        //MC closure test
        //if(subNRB2011 || subNRB2012) doBackgroundSubtraction(selCh,"emu",allShapes,histo,histo+"_NRBctrl", url, Root, true);
        //replace WZ by its estimate from 3rd Lepton SB
        //if(subWZ)doWZSubtraction(selCh,"emu",allShapes,histo,histo+"_3rdLepton");


        // MC Z+jets -> DYExtrapolation
        doDYextrapolation(selCh,allShapes,histo,"pfmet_minus_shapes", DYMET_EXTRAPOL, true);
        //doDYextrapolation(selCh,allShapes,histo,"balancedif_minus_shapes", DYRESP_EXTRAPOL, false);
        //doDYextrapolation(selCh,allShapes,histo,"dphizmet_minus_shapes", DYDPHI_EXTRAPOL, true);

    */
    //replace Z+Jet background by Gamma+Jet estimates
    ////if(subDY)doDYReplacement(selCh,"gamma",allShapes,histo,"met_met");
    /////
    //if(subDY)doDYReplacement(selCh,"gamma",allShapes,histo,"met_redMet");

    //replace data by total MC background
    if(blindData)BlindData(selCh,allShapes,histo, blindWithSignal);

    //interpollate signal sample if desired mass point is not available
    //SignalInterpolation(selCh,allShapesL, allShapes, allShapesR, histo);


    if(doInterf)RescaleForInterference(selCh,allShapes, histo);



    //print event yields from the mt shapes
    //if(!fast)getYieldsFromShape(selCh,allShapes,histo,true); //blind data
    if(!fast)getYieldsFromShape(selCh,allShapes,histo,false); //unblind data

    if(!fast)getEffFromShape(selCh,allShapes,histo);

    if(!fast)showShape(selCh,allShapes,histo,"plot");


    //prepare the output
    dci.shapesFile="zllwimps_"+massStr2+systpostfix+".root";
    TFile *fout=TFile::Open(dci.shapesFile,"recreate");

    //loop on channel/proc/systematics
    for(size_t ich=0; ich<selCh.size(); ich++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            TString chbin = selCh[ich]+AnalysisBins[b];
            fout->mkdir(chbin);
            fout->cd(chbin);
            allCh.push_back(chbin);
            Shape_t shapeSt = allShapes.find(chbin+histo)->second;

            //signals
            dci.nsignalproc = 0;
            size_t nsignal=allShapes.find(chbin+histo)->second.signal.size();
            for(size_t isignal=0; isignal<nsignal; isignal++) {
                TH1* h=shapeSt.signal[isignal];

                TString proc(h->GetTitle());
                if(mass>0 && !proc.Contains(massStr))continue;

                if(proc.Contains("D1")) {
                    if(mass>0 && !proc.Contains("D1("+massStr+"GeV)"))continue;
                }
                if(proc.Contains("D4")) {
                    if(mass>0 && !proc.Contains("D4("+massStr+"GeV)"))continue;
                }
                if(proc.Contains("D5")) {
                    if(mass>0 && !proc.Contains("D5("+massStr+"GeV)"))continue;
                }
                if(proc.Contains("D8")) {
                    if(mass>0 && !proc.Contains("D8("+massStr+"GeV)"))continue;
                }
                if(proc.Contains("D9")) {
                    if(mass>0 && !proc.Contains("D9("+massStr+"GeV)"))continue;
                }
                if(proc.Contains("C3")) {
                    if(mass>0 && !proc.Contains("C3("+massStr+"GeV)"))continue;
                }
                if(proc.Contains("Unpart")) {
                    if(mass>0 && !proc.Contains("Unpart("+massStr+")"))continue;
                }

                if(proc.Contains("DM") && proc.Contains("MV")) {
                    if(mass<0 || MV<0) continue;
                    if(!proc.Contains("DM("+massStr+")MV("+MVStr+")")) continue;
                }

                if(proc.Contains("DM") && proc.Contains("MA")) {
                    if(mass<0 || MA<0) continue;
                    if(!proc.Contains("DM("+massStr+")MA("+MAStr+")")) continue;
                }



                if(mass>0 && proc.Contains("ZH")                        )proc = "ZH"+massStr+"2lMET";

                cout << "############## Signal: " << proc << "##############" << endl;

                //if(atgcpar.Length()>0 && !proc.Contains(massStr)) continue;
                if(atgcpar.Length()>0 && !proc.EndsWith(massStr)) continue;
                if(atgcpar.Length()>0) proc = "ZZ2l2nu_"+massStr2;

                std::vector<std::pair<TString, TH1*> > vars = shapeSt.signalVars[h->GetTitle()];
                std::vector<TString> systs;
                std::vector<TH1*>    hshapes;
                systs.push_back("");
                hshapes.push_back(shapeSt.signal[isignal]);
                for(size_t v=0; v<vars.size(); v++) {
                    //printf("SYSTEMATIC FOR SIGNAL %s : %s\n",h->GetTitle(), vars[v].first.Data());
                    systs.push_back(vars[v].first);
                    hshapes.push_back(vars[v].second);
                }

                convertHistosForLimits_core(dci, proc, AnalysisBins[b], chbin, systs, hshapes);
                if(ich==0 && b==0)allProcs.push_back(proc);
                dci.nsignalproc++;
            }

            //backgrounds
            size_t nbckg=allShapes.find(chbin+histo)->second.bckg.size();
            size_t nNonNullBckg=0;
            for(size_t ibckg=0; ibckg<nbckg; ibckg++) {
                TH1* h=shapeSt.bckg[ibckg];
                if(h->Integral()<=1e-10) continue;
                cout << "\n" << h->GetTitle() << " has rate: " << h->Integral() << endl;


                //check if any bin content has a negative value
                double tot_integralval = h->Integral();
                for(int bin=0; bin<h->GetXaxis()->GetNbins()+1; bin++) {
                    //cout << "bin: " << bin << " val: " << h->GetBinContent(bin) << endl;
                    if(h->GetBinContent(bin)<0) h->SetBinContent(bin,0);
                }
                h->Scale(tot_integralval/h->Integral());
                //cout << "---------------------------------" << endl;
                //for(int bin=0; bin<h->GetXaxis()->GetNbins()+1; bin++) {
                //cout << "bin: " << bin << " val: " << h->GetBinContent(bin) << endl;
                //}


                std::vector<std::pair<TString, TH1*> > vars = shapeSt.bckgVars[h->GetTitle()];

                std::vector<TString> systs;
                std::vector<TH1*>    hshapes;
                systs.push_back("");
                hshapes.push_back(shapeSt.bckg[ibckg]);
                for(size_t v=0; v<vars.size(); v++) {
                    printf("SYSTEMATIC FOR BACKGROUND %s : %s\n",h->GetTitle(), vars[v].first.Data());

                    TH1* h= vars[v].second;

                    //check if any bin content has a negative value
                    double tot_integralval = h->Integral();
                    for(int bin=0; bin<h->GetXaxis()->GetNbins()+1; bin++) {
                        //cout << "bin: " << bin << " val: " << h->GetBinContent(bin) << endl;
                        if(h->GetBinContent(bin)<0) h->SetBinContent(bin,0);
                    }
                    h->Scale(tot_integralval/h->Integral());
                    //cout << "---------------------------------" << endl;
                    //for(int bin=0; bin<h->GetXaxis()->GetNbins()+1; bin++) {
                    //        cout << "bin: " << bin << " val: " << h->GetBinContent(bin) << endl;
                    //}


                    systs.push_back(vars[v].first);
                    hshapes.push_back(h/*vars[v].second*/);
                }

                TString proc(h->GetTitle());
                convertHistosForLimits_core(dci, proc, AnalysisBins[b], chbin, systs, hshapes);
                if(ich==0 && b==0)allProcs.push_back(proc);

                //remove backgrounds with rate=0 (but keep at least one background)
                std::map<RateKey_t, Double_t>::iterator it = dci.rates.find(RateKey_t(proc,chbin));
                if(it==dci.rates.end()) {
                    printf("proc=%s not found --> THIS SHOULD NEVER HAPPENS.  PLEASE CHECK THE CODE\n",proc.Data());
                } else {
                    if(it->second>0) {
                        nNonNullBckg++;
                    } else if(ibckg<nbckg-1 || nNonNullBckg>0) {
                        dci.rates.erase(dci.rates.find(RateKey_t(proc,chbin)));
                    }
                }

            }





            //data
            TH1* h=shapeSt.data;
            std::vector<TString> systs;
            std::vector<TH1*>    hshapes;
            systs.push_back("");
            hshapes.push_back(h);
            TString proc(h->GetTitle());
            convertHistosForLimits_core(dci, proc, AnalysisBins[b], chbin, systs, hshapes);

            //return to parent dir
            fout->cd("..");
        }
    }

    dci.ch.resize(allCh.size());
    std::copy(allCh.begin(), allCh.end(),dci.ch.begin());
    dci.procs.resize(allProcs.size());
    std::copy(allProcs.begin(), allProcs.end(),dci.procs.begin());



    //all done
    fout->Close();

    return dci;
}


//print
//bin-by-bin uncertainty
void convertHistosForLimits_core(DataCardInputs& dci, TString& proc, TString& bin, TString& ch, std::vector<TString>& systs, std::vector<TH1*>& hshapes)
{
    proc.ReplaceAll("#bar{t}","tbar");
    proc.ReplaceAll("Z-#gamma^{*}+jets#rightarrow ll","dy");
    proc.ReplaceAll("#rightarrow","");
    proc.ReplaceAll("(","");
    proc.ReplaceAll(")","");
    proc.ReplaceAll("+","");
    proc.ReplaceAll(" ","");
    proc.ReplaceAll("/","");
    proc.ReplaceAll("#","");
    proc.ReplaceAll("=","");
    proc.ReplaceAll(".","");
    proc.ReplaceAll("^","");
    proc.ReplaceAll("}","");
    proc.ReplaceAll("{","");
    proc.ToLower();

    //cout << "Process: " << proc << endl;
    proc.ReplaceAll("zh1052lmet","ZH_hinv");
    proc.ReplaceAll("zh1152lmet","ZH_hinv");
    proc.ReplaceAll("zh1252lmet","ZH_hinv");
    proc.ReplaceAll("zh1352lmet","ZH_hinv");
    proc.ReplaceAll("zh1452lmet","ZH_hinv");
    proc.ReplaceAll("zh1752lmet","ZH_hinv");
    proc.ReplaceAll("zh2002lmet","ZH_hinv");
    proc.ReplaceAll("zh3002lmet","ZH_hinv");


    proc.ReplaceAll("c31gev","C3");
    proc.ReplaceAll("d11gev","D1");
    proc.ReplaceAll("d41gev","D4");
    proc.ReplaceAll("d51gev","D5");
    proc.ReplaceAll("d81gev","D8");
    proc.ReplaceAll("d91gev","D9");

    proc.ReplaceAll("c310gev","C3");
    proc.ReplaceAll("d110gev","D1");
    proc.ReplaceAll("d410gev","D4");
    proc.ReplaceAll("d510gev","D5");
    proc.ReplaceAll("d810gev","D8");
    proc.ReplaceAll("d910gev","D9");

    proc.ReplaceAll("c3100gev","C3");
    proc.ReplaceAll("d1100gev","D1");
    proc.ReplaceAll("d4100gev","D4");
    proc.ReplaceAll("d5100gev","D5");
    proc.ReplaceAll("d8100gev","D8");
    proc.ReplaceAll("d9100gev","D9");

    proc.ReplaceAll("c3200gev","C3");
    proc.ReplaceAll("d1200gev","D1");
    proc.ReplaceAll("d4200gev","D4");
    proc.ReplaceAll("d5200gev","D5");
    proc.ReplaceAll("d8200gev","D8");
    proc.ReplaceAll("d9200gev","D9");

    proc.ReplaceAll("c3300gev","C3");
    proc.ReplaceAll("d1300gev","D1");
    proc.ReplaceAll("d4300gev","D4");
    proc.ReplaceAll("d5300gev","D5");
    proc.ReplaceAll("d8300gev","D8");
    proc.ReplaceAll("d9300gev","D9");

    proc.ReplaceAll("c3500gev","C3");
    proc.ReplaceAll("d1500gev","D1");
    proc.ReplaceAll("d4500gev","D4");
    proc.ReplaceAll("d5500gev","D5");
    proc.ReplaceAll("d8500gev","D8");
    proc.ReplaceAll("d9500gev","D9");

    proc.ReplaceAll("c31000gev","C3");
    proc.ReplaceAll("d11000gev","D1");
    proc.ReplaceAll("d41000gev","D4");
    proc.ReplaceAll("d51000gev","D5");
    proc.ReplaceAll("d81000gev","D8");
    proc.ReplaceAll("d91000gev","D9");

    proc.ReplaceAll("unpart101","UnPart1p01");
    proc.ReplaceAll("unpart102","UnPart1p02");
    proc.ReplaceAll("unpart104","UnPart1p04");
    proc.ReplaceAll("unpart106","UnPart1p06");
    proc.ReplaceAll("unpart109","UnPart1p09");
    proc.ReplaceAll("unpart110","UnPart1p10");
    proc.ReplaceAll("unpart120","UnPart1p20");
    proc.ReplaceAll("unpart130","UnPart1p30");
    proc.ReplaceAll("unpart140","UnPart1p40");
    proc.ReplaceAll("unpart150","UnPart1p50");
    proc.ReplaceAll("unpart160","UnPart1p60");
    proc.ReplaceAll("unpart170","UnPart1p70");
    proc.ReplaceAll("unpart180","UnPart1p80");
    proc.ReplaceAll("unpart190","UnPart1p90");
    proc.ReplaceAll("unpart200","UnPart2p00");
    proc.ReplaceAll("unpart220","UnPart2p20");

    proc.ReplaceAll("msdmvg10mx50m100","MSDMVg1p0mx50");
    proc.ReplaceAll("msdmvg10mx50m200","MSDMVg1p0mx50");
    proc.ReplaceAll("msdmvg10mx50m300","MSDMVg1p0mx50");
    proc.ReplaceAll("msdmvg10mx50m500","MSDMVg1p0mx50");
    proc.ReplaceAll("msdmvg10mx50m1000","MSDMVg1p0mx50");
    proc.ReplaceAll("msdmvg10mx50m5000","MSDMVg1p0mx50");


    proc.ReplaceAll("zz2l2nu","ZZ");
    proc.ReplaceAll("wz3lnu","WZ");
    proc.ReplaceAll("zjetsdata","Zjets");
    proc.ReplaceAll("topwwztautaudata","EM");
    proc.ReplaceAll("wjetsdata","Wjets");




    for(unsigned int i=0; i<systs.size(); i++) {
        TString syst   = systs[i];
        TH1*    hshape = hshapes[i];
        hshape->SetDirectory(0);

        //Do Renaming and cleaning
        syst.ReplaceAll("down","Down");
        syst.ReplaceAll("up","Up");

        if(syst=="") {
            syst="";
        } else if(syst.BeginsWith("_jes")) {
            syst.ReplaceAll("_jes","_CMS_scale_j");
        } else if(syst.BeginsWith("_jer")) {
            syst.ReplaceAll("_jer","_CMS_res_j"); // continue;//skip res for now
        } else if(syst.BeginsWith("_btag")) {
            syst.ReplaceAll("_btag","_CMS_eff_b");
        } else if(syst.BeginsWith("_pu" )) {
            syst.ReplaceAll("_pu", "_CMS_zllwimps_pu");
        } else if(syst.BeginsWith("_les" )) {
            if(ch.Contains("ee")) syst.ReplaceAll("_les", "_CMS_scale_e");
            if(ch.Contains("mumu")) syst.ReplaceAll("_les", "_CMS_scale_m");
        } else if(syst.BeginsWith("_umet" )) {
            syst.ReplaceAll("_umet", "_CMS_scale_met");
        } else if(syst.BeginsWith("_ren" )) {
            continue;   //already accounted for in QCD scales
        } else if(syst.BeginsWith("_fact" )) {
            continue; //skip this one
        } else if(syst.BeginsWith("_interf" )) {
        } else {
            syst="_CMS_zllwimps"+syst;
        }

        double systUncertainty = hshape->GetBinError(0);
        hshape->SetBinError(0,0);

        //If cut&count keep only 1 bin in the histo
        if(!shape) {
            //hshape = hshape->Rebin(hshape->GetXaxis()->GetNbins(), TString(hshape->GetName())+"_Rebin");
            hshape = hshape->Rebin(hshape->GetXaxis()->GetNbins());
            //make sure to also count the underflow and overflow
            double bin  = hshape->GetBinContent(0) + hshape->GetBinContent(1) + hshape->GetBinContent(2);
            double bine = sqrt(hshape->GetBinError(0)*hshape->GetBinError(0) + hshape->GetBinError(1)*hshape->GetBinError(1) + hshape->GetBinError(2)*hshape->GetBinError(2));
            hshape->SetBinContent(0,0);
            hshape->SetBinError  (0,0);
            hshape->SetBinContent(1,bin);
            hshape->SetBinError  (1,bine);
            hshape->SetBinContent(2,0);
            hshape->SetBinError  (2,0);
        }

        if(syst=="") {
            //central shape (for data call it data_obs)
            hshape->SetName(proc);
            if(proc=="data") {
                hshape->Write("data_obs");
            } else {
                hshape->Write(proc+postfix);

                TString zjetsch="";
                if(ch.Contains("eq0jets")) zjetsch="_eq0jets";
                if(ch.Contains("eq1jets")) zjetsch="_eq1jets";
                if(ch.Contains("lesq1jets")) zjetsch="_lesq1jets";

                if(hshape->Integral()>0) {
                    hshape->SetName(proc+syst);
                    TH1* statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Up");
                    TH1* statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Down");
                    //RENJIE
                    if(proc.Contains("EM")) {
                        statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Up");
                        statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"Down");
                    }
                    if(proc.Contains("Zjets")) {
                        statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+zjetsch+proc+"Up");
                        statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+zjetsch+proc+"Down");
                    }

                    for(int ibin=1; ibin<=statup->GetXaxis()->GetNbins(); ibin++) {
                        statup  ->SetBinContent(ibin,std::min(2*hshape->GetBinContent(ibin), std::max(0.01*hshape->GetBinContent(ibin), statup  ->GetBinContent(ibin) + statup  ->GetBinError(ibin))));
                        statdown->SetBinContent(ibin,std::min(2*hshape->GetBinContent(ibin), std::max(0.01*hshape->GetBinContent(ibin), statdown->GetBinContent(ibin) - statdown->GetBinError(ibin))));
                    }

                    //############################################################################
                    //### bin-by-bin uncertainty
                    //############################################################################
                    bool useBinbyBin(false);
                    useBinbyBin &= (proc.Contains("ZH") || proc.Contains("ZZ") || proc.Contains("WZ") || proc.Contains("EM") ||
                                    proc.Contains("D1") || proc.Contains("D4") || proc.Contains("D5") ||
                                    proc.Contains("D8") ||  proc.Contains("D9") || proc.Contains("C3")
                                   );
                    useBinbyBin &= (shape);


                    if(!useBinbyBin) {

                        if(proc.Contains("EM")) {
                            statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"Up");
                            statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"Down");
                        } else if(proc.Contains("Zjets")) {
                            statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"Up");
                            statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"Down");
                        } else {
                            statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"Up");
                            statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"Down");
                        }

                        if(shape) { //RENJIE  Jun15
                            if(proc.Contains("EM"))
                                //cout << "Skiping " << proc <<" stat" << endl;
                                //cout << "removing EM stat" << endl;
                                dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=1.0;
                            else if(proc.Contains("Zjets"))
                                dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=1.0;
                            else
                                dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=1.0;
                        } else {
                            if(proc.Contains("EM"))
                                //cout << "removing EM stat" << endl;
                                dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                            else if(proc.Contains("Zjets"))
                                dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                            else
                                dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                        }
                    }

                    //############################################################################
                    //### bin-by-bin uncertainty
                    //############################################################################
                    //bool useBinbyBin(true);
                    if(useBinbyBin) {

                        hshape->SetName(proc+syst);
                        vector<int> saveBins;
                        for(int iBin=1; iBin<=hshape->GetXaxis()->GetNbins(); iBin++) {

                            if(hshape->GetBinContent(iBin)<=0) continue;

                            stringstream Bin_t;
                            Bin_t << iBin;
                            string sbin = Bin_t.str();
                            TString Bin = sbin;

                            TH1* statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"_Bin"+Bin+"Up");
                            TH1* statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"_Bin"+Bin+"Down");


                            for(int jbin=1; jbin<=statup->GetXaxis()->GetNbins(); jbin++) {
                                if(jbin!=iBin) continue;
                                statup  ->SetBinContent(jbin,std::min(2*hshape->GetBinContent(jbin), std::max(0.01*hshape->GetBinContent(jbin), statup  ->GetBinContent(jbin) + statup  ->GetBinError(jbin))));
                                statdown->SetBinContent(jbin,std::min(2*hshape->GetBinContent(jbin), std::max(0.01*hshape->GetBinContent(jbin), statdown->GetBinContent(jbin) - statdown->GetBinError(jbin))));
                                //double variance = (statup->Integral()/hshapes[0]->Integral());
                            }

                            //#########Debug###########
                            double variance = (statup->Integral()/hshapes[0]->Integral());
                            TString name_variance = proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin+"Up";
                            if(proc.Contains("EM")) name_variance = proc+postfix+"_CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin+"Up";
                            cout << name_variance << ": " << variance << endl;
                            nuisance_h->Fill(variance);
                            //#########Debug###########
                            //if((proc=="zz2l2nu") && variance<1.0001) {saveBins.push_back(iBin); continue;}
                            if(ch.Contains("eq0jets")) {
                                if(proc.Contains("ZZ") && variance < 1.005) {
                                    saveBins.push_back(iBin);
                                    continue;
                                }
                                if(!proc.Contains("ZZ") && variance < 1.008) {
                                    saveBins.push_back(iBin);
                                    continue;
                                }
                            } else if(ch.Contains("eq1jets")) {
                                if(proc.Contains("ZZ") && variance < 1.01) {
                                    saveBins.push_back(iBin);
                                    continue;
                                }
                                if(!proc.Contains("ZZ") && variance < 1.02) {
                                    saveBins.push_back(iBin);    //save bins
                                    continue;
                                }
                            } else continue;


//                            statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin+"Up");
//                            statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin+"Down");

                            if(proc.Contains("EM")) {
                                statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin+"Up");
                                statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin+"Down");
                            } else {
                                statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin+"Up");
                                statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin+"Down");
                            }

                            if(shape) {
                                if(proc.Contains("EM")) dci.systs["CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin][RateKey_t(proc,ch)]=1.0;
                                else dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin][RateKey_t(proc,ch)]=1.0;

                            } else {
                                if(proc.Contains("EM")) dci.systs["CMS_zllwimps_stat_"+proc+systpostfix+"_Bin"+Bin][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                                else dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_Bin"+Bin][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                            }

                        }

                        //#################### Bin-by-Bin Rest bins #####################
                        /*
                        			for(unsigned int ij=0; ij<saveBins.size(); ij++){
                        			    int iBin=saveBins[ij];
                        			    if(hshape->GetBinContent(iBin)<=0) continue;
                                                    TH1* statup=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"_BinRestUp");
                                                    TH1* statdown=(TH1 *)hshape->Clone(proc+"_stat"+ch+proc+"_BinRestDown");

                                                    statup  ->SetBinContent(iBin,std::min(2*hshape->GetBinContent(iBin), std::max(0.01*hshape->GetBinContent(iBin), statup  ->GetBinContent(iBin) + statup  ->GetBinError(iBin))));
                                                    statdown->SetBinContent(iBin,std::min(2*hshape->GetBinContent(iBin), std::max(0.01*hshape->GetBinContent(iBin), statdown->GetBinContent(iBin) - statdown->GetBinError(iBin))));

                        			}


                        			if(saveBins.size()>0){
                                                    if(proc=="topwwztautaudata"){
                                                        statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+proc+systpostfix+"_BinRestUp");
                                                        statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+proc+systpostfix+"_BinRestDown");
                                                    }else{
                                                        statup  ->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_BinRestUp");
                                                        statdown->Write(proc+postfix+"_CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_BinRestDown");
                                                    }

                                                    if(shape) {
                                                        if(proc=="topwwztautaudata") dci.systs["CMS_zllwimps_stat_"+proc+systpostfix+"_BinRest"][RateKey_t(proc,ch)]=1.0;
                                                        else dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_BinRest"][RateKey_t(proc,ch)]=1.0;

                                                    }else{
                                                        if(proc=="topwwztautaudata") dci.systs["CMS_zllwimps_stat_"+proc+systpostfix+"_BinRest"][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                                                        else dci.systs["CMS_zllwimps_stat_"+ch+"_"+proc+systpostfix+"_BinRest"][RateKey_t(proc,ch)]=(statup->Integral()/hshapes[0]->Integral());
                                                    }
                        			}
                        */




                    }
                    //############################################################################
                    //############################################################################



                    //data driven Systematics
                    if(systUncertainty>0) {
                        printf("SYST in %s - %s - %s = %f\n",bin.Data(), ch.Data(), proc.Data(), systUncertainty);
                        //makesure that syst+stat error is never bigger than 100%
                        //double valerr, val  = hshape->IntegralAndError(1,hshape->GetXaxis()->GetNbins(),valerr);
                        //if(sqrt(pow(valerr,2)+pow(systUncertainty,2))>val){systUncertainty = sqrt(std::max(0.0, pow(val,2) - pow(valerr,2)));}

                        //add syst uncertainty as bin dependent or not
                        //dci.systs["CMS_zllwimps_sys_"+bin+"_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;

                        TString wjetsch="";
                        if(ch.Contains("ee")) wjetsch="ee";
                        if(ch.Contains("mumu")) wjetsch="mumu";

                        if(proc.Contains("EM") || proc.Contains("Zjets")) dci.systs["CMS_zllwimps_sys_"+ch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                        else if(proc.Contains("Wjets")) dci.systs["CMS_zllwimps_sys_"+wjetsch+"_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                        else dci.systs["CMS_zllwimps_sys_"+proc+systpostfix][RateKey_t(proc,ch)]=systUncertainty;
                    }
                }
            }
        } else if(runSystematics && proc!="data" && (syst.Contains("Up") || syst.Contains("Down"))) {
            //if empty histogram --> no variation is applied
            if(hshape->Integral()<hshapes[0]->Integral()*0.01 || isnan((float)hshape->Integral())) {
                hshape->Reset();
                hshape->Add(hshapes[0],1);
            }

            //write variation to file
            hshape->SetName(proc+syst);
            hshape->Write(proc+postfix+syst);
        } else if(runSystematics) {
            //for one sided systematics the down variation mirrors the difference bin by bin
            hshape->SetName(proc+syst);
            hshape->Write(proc+postfix+syst+"Up");
            TH1 *hmirrorshape=(TH1 *)hshape->Clone(proc+syst+"Down");
            for(int ibin=1; ibin<=hmirrorshape->GetXaxis()->GetNbins(); ibin++) {
                //double bin = hmirrorshape->GetBinContent(ibin);
                double bin = 2*hshapes[0]->GetBinContent(ibin)-hmirrorshape->GetBinContent(ibin);
                if(bin<0)bin=0;
                hmirrorshape->SetBinContent(ibin,bin);
            }
            if(hmirrorshape->Integral()<=0)hmirrorshape->SetBinContent(1, 1E-10);
            hmirrorshape->Write(proc+postfix+syst+"Down");
        }

        if(runSystematics && syst!="") {
            TString systName(syst);
            systName.ReplaceAll("Up","");
            systName.ReplaceAll("Down","");//  systName.ReplaceAll("_","");
            if(systName.First("_")==0)systName.Remove(0,1);


            TH1 *temp=(TH1*) hshape->Clone();
            temp->Add(hshapes[0],-1);
            if(temp->Integral()!=0 || syst.BeginsWith("_interf" )) {
                if(shape && !systName.Contains("_pdf")) { //convert pdf shape uncertaint to normalization one
                    dci.systs[systName][RateKey_t(proc,ch)]=1.0;
                } else {
                    double Unc = 1 + fabs(temp->Integral()/hshapes[0]->Integral());
                    if(dci.systs.find(systName)==dci.systs.end() || dci.systs[systName].find(RateKey_t(proc,ch))==dci.systs[systName].end() ) {
                        dci.systs[systName][RateKey_t(proc,ch)]=Unc;
                    } else {
                        dci.systs[systName][RateKey_t(proc,ch)]=(dci.systs[systName][RateKey_t(proc,ch)] + Unc)/2.0;
                    }
                }
            }


            delete temp;
        } else if(proc!="data" && syst=="") {
            dci.rates[RateKey_t(proc,ch)]= hshape->Integral();
        } else if(proc=="data" && syst=="") {
            dci.obs[RateKey_t("obs",ch)]=hshape->Integral();
        }
    }

}//convertHistosForLimits_core END

void doDYextrapolation(std::vector<TString>& selCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString DY_EXTRAPOL_SHAPES, float DY_EXTRAPOL, bool isleftCR)
{
    cout << "\n#################### doDYextrapolation #######################\n" << endl;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t i=0; i<selCh.size(); i++) {

            THStack *stack = new THStack("stack","stack");

            Shape_t& shapeChan = allShapes.find(selCh[i]+AnalysisBins[b]+DY_EXTRAPOL_SHAPES)->second;
            TH1* hChan_data=shapeChan.data;
            cout << "DY extrapolation shapes: " << selCh[i]+AnalysisBins[b]+"_"+DY_EXTRAPOL_SHAPES << endl;
            cout << "Bins: " << hChan_data->GetXaxis()->GetNbins() << endl;

            TH1* hChan_MCnonDY = (TH1*)shapeChan.totalBckg->Clone("hChan_totMCnonDY");
            TH1* hChan_MCDY=NULL;//= (TH1*)shapeChan.totalBckg->Clone("hChan_MCDY");
            hChan_MCnonDY->Reset();
            //hChan_MCDY->Reset();
            for(size_t ibckg=0; ibckg<shapeChan.bckg.size(); ibckg++) {
                TString proc(shapeChan.bckg[ibckg]->GetTitle());
                //cout << "proc: " << proc << endl;
                if(proc.Contains("Z+jets")) {
                    hChan_MCDY=(TH1*)shapeChan.bckg[ibckg]->Clone("hChan_MCDY");
                    hChan_MCDY->Reset();
                    hChan_MCDY->Add(shapeChan.bckg[ibckg], 1);
                } else {
                    hChan_MCnonDY->Add(shapeChan.bckg[ibckg], 1);
                    stack->Add(shapeChan.bckg[ibckg],"HIST");
                }

                //stack->Add(shapeChan.bckg[ibckg],"HIST");
            }

            int thre_bin(0);
            int total_bin=hChan_MCnonDY->GetXaxis()->GetNbins();
            for(int bin=0; bin<=hChan_MCnonDY->GetXaxis()->GetNbins()+1; bin++) {
                double lowedge = hChan_MCnonDY->GetBinLowEdge(bin);
                if(lowedge >= DY_EXTRAPOL) {
                    thre_bin = bin;
                    break;
                }
            }

            double data_CR(0.);
            double subt_CR(0.);
            double zjet_SR(0.),zjet_CR(0.);
            double data_err2_CR(0.);
            double subt_err2_CR(0.);
            double zjet_err2_SR(0.),zjet_err2_CR(0.);

            if(isleftCR) {
                data_CR = hChan_data	->IntegralAndError(0,thre_bin-1,data_err2_CR);
                subt_CR = hChan_MCnonDY	->IntegralAndError(0,thre_bin-1,subt_err2_CR);
                zjet_CR = hChan_MCDY	->IntegralAndError(0,thre_bin-1,zjet_err2_CR);
                zjet_SR = hChan_MCDY	->IntegralAndError(thre_bin,total_bin+1,zjet_err2_SR);
            } else {
                data_CR = hChan_data	->IntegralAndError(thre_bin,total_bin+1,data_err2_CR);
                subt_CR = hChan_MCnonDY	->IntegralAndError(thre_bin,total_bin+1,subt_err2_CR);
                zjet_CR = hChan_MCDY	->IntegralAndError(thre_bin,total_bin+1,zjet_err2_CR);
                zjet_SR = hChan_MCDY	->IntegralAndError(0,thre_bin-1,zjet_err2_SR);
            }

            double sf=(data_CR-subt_CR)/zjet_CR;
            if(sf<0) sf=1;
            double zjet_Est = sf*zjet_SR;
            cout << "sf: " << sf << "\t" << "Est: " << zjet_Est << "\t" << "MC: " << zjet_SR << endl;

            hChan_MCDY->Scale(sf);
            stack->Add(hChan_MCDY,"HIST");


            TCanvas *c = new TCanvas("c", "c", 700, 700);
            TPad* t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.00);
            t1->Draw();
            t1->SetLogy(true);
            t1->cd();
            //t1->SetBottomMargin(0.3);
            t1->SetRightMargin(0.05);
            //c->Divide(1,2);

            stack->Draw("");
            hChan_data->Draw("E1 same");

            stack->SetMinimum(5e-3);
            stack->GetXaxis()->SetTitle(hChan_data->GetXaxis()->GetTitle());
            stack->GetYaxis()->SetTitle("Events");

            c->SaveAs(selCh[i]+AnalysisBins[b]+"_"+DY_EXTRAPOL_SHAPES+"_DYEXTRAPOL.pdf");
            c->SaveAs(selCh[i]+AnalysisBins[b]+"_"+DY_EXTRAPOL_SHAPES+"_DYEXTRAPOL.png");
            delete c;


            //add data-driven backgrounds
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
            TH1* hChan_DATADY = NULL;

            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                if(proc.Contains("Z+jets")) {
                    hChan_DATADY=(TH1*)shapeChan_SI.bckg[ibckg]->Clone("hChan_DATADY");
                    hChan_DATADY->SetTitle("Z+jets (data)");
                    //hChan_DATADY->Reset();
                    //remove MC Z+jets
                    shapeChan_SI.bckg.erase(shapeChan_SI.bckg.begin()+ibckg);
                    ibckg--;
                }
            }


            //update  Feb 24
            TString tag_dy = selCh[i]+AnalysisBins[b];
            if(tag_dy=="eelesq1jets") sf = 1.25;//1.49205;
            else if(tag_dy=="mumulesq1jets") sf = 1.16;//1.46399;
            else sf = 1.;

            cout << "Scale Factor: " << sf << " DY MC: " << hChan_DATADY->Integral() << endl;
            hChan_DATADY->Scale(sf);

            //assign Systematics
            if(tag_dy=="eelesq1jets") hChan_DATADY->SetBinError(0, 0.10/*0.193361*/*hChan_DATADY->Integral());
            if(tag_dy=="mumulesq1jets") hChan_DATADY->SetBinError(0, 0.11/*0.48201*/*hChan_DATADY->Integral());

            cout << "DY Extrapolation: " << hChan_DATADY->Integral() << endl;
            if(hChan_DATADY->Integral() > 1E-6) shapeChan_SI.bckg.push_back(hChan_DATADY);

            cout << "\n";


        }
    }
    cout << "\n#################### doDYextrapolation #######################\n" << endl;
}


void doWjetsBackground(std::vector<TString>& selCh,map<TString, Shape_t>& allShapes, TString mainHisto)
{
    cout << "\n#################### doWjetsBackground #######################\n" << endl;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t i=0; i<selCh.size(); i++) {
            TString label_Syst = selCh[i]+AnalysisBins[b];

            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
            Shape_t& shapeChan_Wjet = allShapes.find(selCh[i]+AnalysisBins[b]+"FR_WjetCtrl_"+mainHisto)->second;
            TH1* hChan_data=shapeChan_Wjet.data;

            double valerr_hChan_data, val_hChan_data;
            val_hChan_data = hChan_data->IntegralAndError(1,hChan_data->GetXaxis()->GetNbins(),valerr_hChan_data);

            cout << "Wjets Shapes: " << selCh[i]+AnalysisBins[b]+"_FR_WjetCtrl_"+mainHisto << "\t"
                 << val_hChan_data << endl;

            TH1* h_Wjet = NULL;
            h_Wjet = (TH1*)hChan_data->Clone("Wjets (data)");
            h_Wjet->SetTitle("Wjets (data)");

            //remove mc Wjet estimation
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                if(proc.Contains("W+jets")) {
                    shapeChan_SI.bckg.erase(shapeChan_SI.bckg.begin()+ibckg);
                    ibckg--;
                }
            }

            //add the background estimate
            if(val_hChan_data <= 0) continue;

            double minVal_bin(999.);
            for(int b=1; b<=h_Wjet->GetXaxis()->GetNbins()+1; b++) {
                double val = h_Wjet->GetBinContent(b);
                double err = h_Wjet->GetBinError(b);

                cout << "bin: " << b << "\t" << val << " +/- " << err << endl;
                if(val<minVal_bin) minVal_bin = val;
            }

            if(minVal_bin<0) {
                for(int b=1; b<=h_Wjet->GetXaxis()->GetNbins()+1; b++) {
                    double val = h_Wjet->GetBinContent(b);
                    if(val!=0) val -= minVal_bin;
                    //cout << "bin: " << b << "\t" << val << endl;
                    h_Wjet->SetBinContent(b, val);
                }
            }

            double newsum = h_Wjet->Integral();
            h_Wjet->Scale(val_hChan_data/newsum);

            cout << ">>>>>>>>>> remove negative bins <<<<<<<<<<" << endl;
            for(int b=1; b<=h_Wjet->GetXaxis()->GetNbins()+1; b++) {
                double val = h_Wjet->GetBinContent(b);
                double err = h_Wjet->GetBinError(b);
                cout << "bin: " << b << "\t" << val << " +/- " << err << endl;
            }


            double sysErr(0.);
            if(label_Syst.Contains("ee"))   sysErr = val_hChan_data*WjetsSyst_ee;
            else if(label_Syst.Contains("mumu")) sysErr = val_hChan_data*WjetsSyst_mm;
            h_Wjet->SetBinError(0,sysErr);



            shapeChan_SI.bckg.push_back(h_Wjet);
        }
    }
    cout << "\n#################### doWjetsBackground #######################\n" << endl;
}



void doQCDBackground(std::vector<TString>& selCh,map<TString, Shape_t>& allShapes, TString mainHisto)
{
    cout << "\n#################### doQCDBackground #######################\n" << endl;
    for(size_t b=0; b<AnalysisBins.size(); b++) {
        for(size_t i=0; i<selCh.size(); i++) {


            cout << "QCD Shapes: " << selCh[i]+AnalysisBins[b]+mainHisto << endl;

            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
            TH1* hChan_data=shapeChan_SI.data;

            double valerr_hChan_data, val_hChan_data;
            val_hChan_data = hChan_data->IntegralAndError(1,hChan_data->GetXaxis()->GetNbins(),valerr_hChan_data);

            cout << "QCD Shapes: " << selCh[i]+AnalysisBins[b]+mainHisto << "\t"
                 << val_hChan_data << endl;

        }
    }
    cout << "\n#################### doQCDBackground #######################\n" << endl;
}








void dodataDrivenWWtW(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, bool isMCclosureTest)
{

    if(isMCclosureTest) cout << "\n############ dodataDrivenWWtW (MCclosureTest) ###############\n" << endl;
    else 		cout << "\n#################### dodataDrivenWWtW #######################\n" << endl;
    FILE* pFile = NULL;
    if(isMCclosureTest) pFile = fopen("DataDrivenWWtW_MCclosureTest.tex","w");
    else 		pFile = fopen("DataDrivenWWtW.tex","w");
    fprintf(pFile,"\n\n\n");
    fprintf(pFile,"\\setlength{\\tabcolsep}{2pt}\n");
    fprintf(pFile,"\\begin{table}[ht!]\n");
    //fprintf(pFile,"\\begin{sidewaystable}[htp]\n");
    fprintf(pFile,"\\begin{center}\n");
    fprintf(pFile,"\\centering\n");
    fprintf(pFile,"\\scriptsize\\begin{tabular}{c|ccccc|cc|cc}\n");
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\hline\n");
    if(isMCclosureTest) {
        fprintf(pFile,"Bin & $N_{ee}^{mc}$ & $N_{\\mu\\mu}^{mc}$ & $N_{e\\mu}^{mc}$ & $N_{e\\mu}^{mc,subtr}$ & $N_{e\\mu}^{corr}$ \n");
        fprintf(pFile,"& $N_{bkg,ee}^{closure}$ & $N_{bkg,ee}^{mc}$ & $N_{bkg,\\mu\\mu}^{closure}$ & $N_{bkg,\\mu\\mu}^{mc}$ \\\\\n");
    } else {
        fprintf(pFile,"Bin & $N_{ee}^{data}$ & $N_{\\mu\\mu}^{data}$ & $N_{e\\mu}^{data}$ & $N_{e\\mu}^{mc,subtr}$ & $N_{e\\mu}^{corr}$ \n");
        fprintf(pFile,"& $N_{bkg,ee}^{est}$ & $N_{bkg,ee}^{mc}$ & $N_{bkg,\\mu\\mu}^{est}$ & $N_{bkg,\\mu\\mu}^{mc}$ \\\\\n");
    }
    fprintf(pFile,"\\hline\n");

    for(size_t b=0; b<AnalysisBins.size(); b++) {


        double N_data_ee(0.), N_data_mm(0.), N_data_em(0.);
        double N_mc_ee(0.), N_mc_mm(0.);
        double N_mcsubtr_em(0.);

        double Err_data_ee(0.), Err_data_mm(0.), Err_data_em(0.);
        double Err_mc_ee(0.), Err_mc_mm(0.);
        double Err_mcsubtr_em(0.);

        for(size_t i=0; i<selCh.size(); i++) {

            cout << ctrlCh+AnalysisBins[b]+mainHisto << " : " << selCh[i]+AnalysisBins[b]+mainHisto << endl;

            Shape_t& shapeCtrl_SI = allShapes.find(ctrlCh+AnalysisBins[b]+mainHisto)->second;
            TH1* hCtrl_data=shapeCtrl_SI.data;
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
            TH1* hChan_data=shapeChan_SI.data;

            if(isMCclosureTest) {
                hCtrl_data=shapeCtrl_SI.totalBckg;
                hChan_data=shapeChan_SI.totalBckg;
            }


            TH1* hCtrl_MCnonNRB = (TH1*)shapeCtrl_SI.totalBckg->Clone("hCtrl_MCnonNRB");
            TH1* hCtrl_MCNRB = (TH1*)shapeCtrl_SI.totalBckg->Clone("hCtrl_MCNRB");
            TH1* hChan_MCnonNRB = (TH1*)shapeChan_SI.totalBckg->Clone("hChan_MCnonNRB");
            TH1* hChan_MCNRB = (TH1*)shapeChan_SI.totalBckg->Clone("hChan_MCNRB");
            hCtrl_MCnonNRB->Reset();
            hCtrl_MCNRB->Reset();
            hChan_MCnonNRB->Reset();
            hChan_MCNRB->Reset();


            //emu channel
            for(size_t ibckg=0; ibckg<shapeCtrl_SI.bckg.size(); ibckg++) {
                TString proc(shapeCtrl_SI.bckg[ibckg]->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu")
                        || proc.Contains("WZ#rightarrow 3l#nu")
                        || proc.Contains("Z+jets")
                        || proc.Contains("W+jets")
                  ) {
                    hCtrl_MCnonNRB->Add(shapeCtrl_SI.bckg[ibckg], 1);
                }
                if(proc.Contains("t#bar{t}")
                        || proc.Contains("Single top")
                        || proc.Contains("WW#rightarrow l#nul#nu")
                        || proc.Contains("Z#rightarrow #tau#tau")
                  ) {
                    hCtrl_MCNRB->Add(shapeCtrl_SI.bckg[ibckg], 1);
                }
            }

            //ee, mumu channel
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu")
                        || proc.Contains("WZ#rightarrow 3l#nu")
                        || proc.Contains("Z+jets")
                        || proc.Contains("W+jets")
                  ) {
                    hChan_MCnonNRB->Add(shapeChan_SI.bckg[ibckg], 1);
                }
                if(proc.Contains("t#bar{t}")
                        || proc.Contains("Single top")
                        || proc.Contains("WW#rightarrow l#nul#nu")
                        || proc.Contains("Z#rightarrow #tau#tau")
                  ) {
                    hChan_MCNRB->Add(shapeChan_SI.bckg[ibckg], 1);
                }
            }


            double valerr_hCtrl_data, val_hCtrl_data;
            val_hCtrl_data = hCtrl_data->IntegralAndError(1,hCtrl_data->GetXaxis()->GetNbins(),valerr_hCtrl_data);

            double valerr_hChan_data, val_hChan_data;
            val_hChan_data = hChan_data->IntegralAndError(1,hChan_data->GetXaxis()->GetNbins(),valerr_hChan_data);


            double valerr_hCtrl_MCnonNRB, val_hCtrl_MCnonNRB;
            double valerr_hCtrl_MCNRB, val_hCtrl_MCNRB;
            val_hCtrl_MCnonNRB = hCtrl_MCnonNRB->IntegralAndError(1,hCtrl_MCnonNRB->GetXaxis()->GetNbins(),valerr_hCtrl_MCnonNRB);
            val_hCtrl_MCNRB = hCtrl_MCNRB->IntegralAndError(1,hCtrl_MCNRB->GetXaxis()->GetNbins(),valerr_hCtrl_MCNRB);

            double valerr_hChan_MCnonNRB, val_hChan_MCnonNRB;
            double valerr_hChan_MCNRB, val_hChan_MCNRB;
            val_hChan_MCnonNRB = hChan_MCnonNRB->IntegralAndError(1,hChan_MCnonNRB->GetXaxis()->GetNbins(),valerr_hChan_MCnonNRB);
            val_hChan_MCNRB = hChan_MCNRB->IntegralAndError(1,hChan_MCNRB->GetXaxis()->GetNbins(),valerr_hChan_MCNRB);


            cout << "val_hCtrl_data: " << val_hCtrl_data << "\t"
                 << "val_hCtrl_MCnonNRB: " << val_hCtrl_MCnonNRB << "\t"
                 << "val_hCtrl_MCNRB:	" << val_hCtrl_MCNRB << endl;

            cout << "val_hChan_data: " << val_hChan_data << "\t"
                 << "val_hChan_MCnonNRB: " << val_hChan_MCnonNRB << "\t"
                 << "val_hChan_MCNRB:	" << val_hChan_MCNRB << endl;


            if(selCh[i].Contains("ee")) {
                N_data_ee = val_hChan_data;
                N_mc_ee = val_hChan_MCNRB;
                Err_data_ee = valerr_hChan_data;
                Err_mc_ee = valerr_hChan_MCNRB;
            }
            if(selCh[i].Contains("mumu")) {
                N_data_mm = val_hChan_data;
                N_mc_mm = val_hChan_MCNRB;
                Err_data_mm = valerr_hChan_data;
                Err_mc_mm = valerr_hChan_MCNRB;
            }
            if(ctrlCh.Contains("emu")) {
                N_data_em = val_hCtrl_data;
                N_mcsubtr_em = val_hCtrl_MCnonNRB;
                Err_data_em = valerr_hCtrl_data;
                Err_mcsubtr_em = valerr_hCtrl_MCnonNRB;
            }


        } //ee, mumu, emu END


        cout << endl;
        cout << "N_data_ee: " << N_data_ee << "\t"
             << "N_data_mm: " << N_data_mm << "\t"
             << "N_data_em: " << N_data_em << "\t"
             << "N_mcsubtr_em: " << N_mcsubtr_em << endl;

        double k_ee = 0.5*sqrt(N_data_ee/N_data_mm);
        double k_mm = 0.5*sqrt(N_data_mm/N_data_ee);
        double N_data_corr = N_data_em-N_mcsubtr_em;
        double Err_data_corr = sqrt(pow(Err_data_em,2)+pow(Err_mcsubtr_em,2));
        double N_est_ee = N_data_corr*k_ee;
        double N_est_mm = N_data_corr*k_mm;


        double Err_est_ee = sqrt( (pow(Err_data_em,2)+pow(Err_mcsubtr_em,2))*N_data_ee/(4*N_data_mm)
                                  + pow(N_data_em-N_mcsubtr_em,2)*(pow(Err_data_ee/N_data_ee,2)+pow(Err_data_mm/N_data_mm,2))*(N_data_ee/(16*N_data_mm)) );

        double Err_est_mm = sqrt( (pow(Err_data_em,2)+pow(Err_mcsubtr_em,2))*N_data_mm/(4*N_data_ee)
                                  + pow(N_data_em-N_mcsubtr_em,2)*(pow(Err_data_ee/N_data_ee,2)+pow(Err_data_mm/N_data_mm,2))*(N_data_mm/(16*N_data_ee)) );



        cout << "N_est_ee: " << N_est_ee << "\t"
             << "N_est_mm: " << N_est_mm << endl;

        if(isMCclosureTest) {
            double sysEE = fabs(N_est_ee-N_mc_ee)/N_mc_ee;
            double sysMM = fabs(N_est_mm-N_mc_mm)/N_mc_mm;

            cout << "AnalysisBins: " << AnalysisBins[b] << " sysEE: " << sysEE << "\t" << "sysMM: " << sysMM << endl;

            if(AnalysisBins[b] == "eq0jets") {
                WWtopSyst_ee0jet = sysEE;
                WWtopSyst_mm0jet = sysMM;
            }
            if(AnalysisBins[b] == "eq1jets") {
                WWtopSyst_ee1jet = sysEE;
                WWtopSyst_mm1jet = sysMM;
            }
            if(AnalysisBins[b] == "lesq1jets") {
                WWtopSyst_eelesq1jet = sysEE;
                WWtopSyst_mmlesq1jet = sysMM;
            }
        }





        if(isMCclosureTest) {
            fprintf(pFile,"%s & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f  \\\\ \n",
                    (AnalysisBins[b]).Data(),
                    N_data_ee,Err_data_ee,
                    N_data_mm,Err_data_mm,
                    N_data_em,Err_data_em,
                    N_mcsubtr_em,Err_mcsubtr_em,
                    N_data_corr,Err_data_corr,
                    N_est_ee,Err_est_ee, N_mc_ee,Err_mc_ee,
                    N_est_mm,Err_est_mm, N_mc_mm,Err_mc_mm
                   );
        } else {
            fprintf(pFile,"%s & %.0f $\\pm$ %.2f & %.0f $\\pm$ %.2f & %.0f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f & %.2f $\\pm$ %.2f  \\\\ \n",
                    (AnalysisBins[b]).Data(),
                    N_data_ee,Err_data_ee,
                    N_data_mm,Err_data_mm,
                    N_data_em,Err_data_em,
                    N_mcsubtr_em,Err_mcsubtr_em,
                    N_data_corr,Err_data_corr,
                    N_est_ee,Err_est_ee, N_mc_ee,Err_mc_ee,
                    N_est_mm,Err_est_mm, N_mc_mm,Err_mc_mm
                   );
        }





        if(!isMCclosureTest) {
            //add data-driven backgrounds
            for(size_t i=0; i<selCh.size(); i++) {
                TString labelChan = selCh[i]+AnalysisBins[b];
                Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
                TH1* hChan_MCNRB = (TH1*)shapeChan_SI.totalBckg->Clone("hChan_MCNRB");
                hChan_MCNRB->Reset();

                //ee, mumu channel
                for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                    TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                    if(proc.Contains("t#bar{t}")
                            || proc.Contains("Single top")
                            || proc.Contains("WW#rightarrow l#nul#nu")
                            || proc.Contains("Z#rightarrow #tau#tau")
                      ) {
                        hChan_MCNRB->Add(shapeChan_SI.bckg[ibckg], 1);
                        //remove the separate parts
                        shapeChan_SI.bckg.erase(shapeChan_SI.bckg.begin()+ibckg);
                        ibckg--;

                    }
                }

                TH1* NonResonant = NULL;
                NonResonant = (TH1*)hChan_MCNRB->Clone("Top/WW/Ztautau (data)");
                NonResonant->SetTitle("Top/WW/Ztautau (data)");
                //add the background estimate

                double dataDrivenScale(1.0);
                double dataDarvenScale_err(1.0);
                if(selCh[i].Contains("ee")) {
                    dataDrivenScale = N_est_ee/N_mc_ee;
                    dataDarvenScale_err = dataDrivenScale*sqrt(pow(Err_est_ee/N_est_ee,2)+pow(Err_mc_ee/N_mc_ee,2));
                }
                if(selCh[i].Contains("mumu")) {
                    dataDrivenScale = N_est_mm/N_mc_mm;
                    dataDarvenScale_err = dataDrivenScale*sqrt(pow(Err_est_mm/N_est_mm,2)+pow(Err_mc_mm/N_mc_mm,2));
                }

                cout << "DATA/MC scale factors: " << dataDrivenScale << endl;

                for(int b=1; b<=NonResonant->GetXaxis()->GetNbins()+1; b++) {
                    double val = NonResonant->GetBinContent(b);
                    double newval = val*dataDrivenScale;
                    double err = NonResonant->GetBinError(b);
                    double newerr = newval*sqrt(pow(err/val,2)+pow(dataDarvenScale_err/dataDrivenScale,2));

                    NonResonant->SetBinContent(b, newval );
                    if(newval!=0)NonResonant->SetBinError(b, newerr );
                }

                double Syst_WWTop(0.);
                if(labelChan.Contains("eeeq0jets")) 	Syst_WWTop = WWtopSyst_ee0jet;
                if(labelChan.Contains("mumueq0jets")) 	Syst_WWTop = WWtopSyst_mm0jet;
                if(labelChan.Contains("eeeq1jets")) 	Syst_WWTop = WWtopSyst_ee1jet;
                if(labelChan.Contains("mumueq1jets")) 	Syst_WWTop = WWtopSyst_mm1jet;
                if(labelChan.Contains("eelesq1jets")) 	Syst_WWTop = WWtopSyst_eelesq1jet;
                if(labelChan.Contains("mumulesq1jets")) 	Syst_WWTop = WWtopSyst_mmlesq1jet;


                cout << "labelchan: " << labelChan << " >>> Adding syst: " << Syst_WWTop << endl;
                NonResonant->SetBinError(0, Syst_WWTop*NonResonant->Integral());


                shapeChan_SI.bckg.push_back(NonResonant);
            }
        }

        cout << "\n\n";


    } //0jet, 1jet END


    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\hline\n");
    fprintf(pFile,"\\end{tabular}\n");
    fprintf(pFile,"\\end{center}\n");
    //fprintf(pFile,"\\end{sidewaystable}\n");
    fprintf(pFile,"\\end{table}\n");
    fprintf(pFile,"\n\n\n");
    fclose(pFile);
    if(isMCclosureTest) cout << "############ dodataDrivenWWtW (MCclosureTest) ###############\n" << endl;
    else 		cout << "#################### dodataDrivenWWtW #######################\n" << endl;
}

//
void doBackgroundSubtraction(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto, TString url, JSONWrapper::Object &Root, bool isMCclosureTest)
{

    cout << "########################## doBackgroundSubtraction ##########################" << endl;
    string Lcol   = "\\begin{tabular}{|c";
    string Lchan  = "channel";
    string Lalph1 = "$\\alpha$ measured";
    string Lalph2 = "$\\alpha$ used";
    string Lyield   = "yield data";
    string LyieldMC = "yield mc";

    string Ccol   = "\\begin{tabular}{|c|c|c|c|c|c|";
    string Cname  = "channel & $\\alpha$ measured & $\\alpha$ used & yields $e\\mu$ & yield data & yield mc";
    string Cval   = "";
    FILE* pFile = NULL;
    if(!fast) {
        pFile = fopen("NonResonnant.tex","w");
        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
    }

    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Lcol += " |c";
            Lchan += string(" &")+selCh[i]+string(" - ")+AnalysisBins[b];
            Cval   = selCh[i]+string(" - ")+AnalysisBins[b];


            Shape_t& shapeCtrl_SB = allShapes.find(ctrlCh+AnalysisBins[b]+sideBandHisto)->second;
            TH1* hCtrl_SB=shapeCtrl_SB.data;

            Shape_t& shapeCtrl_SI = allShapes.find(ctrlCh+AnalysisBins[b]+mainHisto)->second;
            TH1* hCtrl_SI=shapeCtrl_SI.data;

            Shape_t& shapeChan_SB = allShapes.find(selCh[i]+AnalysisBins[b]+sideBandHisto)->second;
            TH1* hChan_SB=shapeChan_SB.data;

            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;


            if(isMCclosureTest) {
                hCtrl_SB = shapeCtrl_SB.totalBckg;
                hCtrl_SI = shapeCtrl_SI.totalBckg;
                hChan_SB = shapeChan_SB.totalBckg;
            }

            cout << ">>>>>>>>>>>>>>>>> 2" << endl;

            //to subtract non WW/Top background from data
            TH1* tosubtrChan_SB = (TH1*)shapeChan_SB.totalBckg->Clone("tosubtrChan_SB");
            TH1* tosubtrCtrl_SB = (TH1*)shapeCtrl_SB.totalBckg->Clone("tosubtrCtrl_SB");
            tosubtrChan_SB->Reset();
            tosubtrCtrl_SB->Reset();


            for(size_t ibckg=0; ibckg<shapeChan_SB.bckg.size(); ibckg++) {
                TString proc(shapeChan_SB.bckg[ibckg]->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu")
                        || proc.Contains("WZ#rightarrow 3l#nu")
                        || proc.Contains("Z+jets")
                        || proc.Contains("W+jets")
                  ) {
                    tosubtrChan_SB->Add(shapeChan_SB.bckg[ibckg], 1);
                }
            }
            for(size_t ibckg=0; ibckg<shapeCtrl_SB.bckg.size(); ibckg++) {
                TString proc(shapeCtrl_SB.bckg[ibckg]->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu")
                        || proc.Contains("WZ#rightarrow 3l#nu")
                        || proc.Contains("Z+jets")
                        || proc.Contains("W+jets")
                  ) {
                    tosubtrCtrl_SB->Add(shapeCtrl_SB.bckg[ibckg], 1);
                }
            }

            double alpha=0 ,alpha_err=0;


            // >= 1 b-jets
            if(hCtrl_SB->GetBinContent(5)>0) {
                //new for monoZ analysis
                alpha = (hChan_SB->GetBinContent(5)-tosubtrChan_SB->GetBinContent(5))/(hCtrl_SB->GetBinContent(5)-tosubtrCtrl_SB->GetBinContent(5));
                double err1 = sqrt(pow(hChan_SB->GetBinError(5),2)+pow(tosubtrChan_SB->GetBinError(5),2));
                double err2 = sqrt(pow(hCtrl_SB->GetBinError(5),2)+pow(tosubtrCtrl_SB->GetBinError(5),2));
                alpha_err = alpha * sqrt(pow(err1/(hChan_SB->GetBinContent(5)-tosubtrChan_SB->GetBinContent(5)),2)
                                         + pow(err2/(hCtrl_SB->GetBinContent(5)-tosubtrCtrl_SB->GetBinContent(5)),2));
                cout << "alpha: " << alpha << " +/- " << alpha_err << endl;
            }

            // =0 b-jets
            int BinPosition = 2;
            double alpha_bveto=0 ,alpha_bveto_err=0;
            if(hCtrl_SB->GetBinContent(BinPosition)>0) {
                alpha_bveto = (hChan_SB->GetBinContent(BinPosition)-tosubtrChan_SB->GetBinContent(BinPosition))
                              /(hCtrl_SB->GetBinContent(BinPosition)-tosubtrCtrl_SB->GetBinContent(BinPosition));
                double err1 = sqrt(pow(hChan_SB->GetBinError(BinPosition),2)+pow(tosubtrChan_SB->GetBinError(BinPosition),2));
                double err2 = sqrt(pow(hCtrl_SB->GetBinError(BinPosition),2)+pow(tosubtrCtrl_SB->GetBinError(BinPosition),2));
                alpha_bveto_err = alpha_bveto * sqrt(pow(err1/(hChan_SB->GetBinContent(BinPosition)-tosubtrChan_SB->GetBinContent(BinPosition)),2)
                                                     + pow(err2/(hCtrl_SB->GetBinContent(BinPosition)-tosubtrCtrl_SB->GetBinContent(BinPosition)),2));
                cout << "alpha_bveto: " << alpha_bveto << " +/- " << alpha_bveto_err << endl;
            }

            cout << "<--------- Alpha systematics --------->" << endl;
            double syst_alpha = fabs(alpha_bveto-alpha)/alpha;
            double syst_alpha_err = pow(alpha_bveto_err/alpha_bveto,2);
            syst_alpha_err += pow(alpha_err/alpha,2);
            syst_alpha_err = sqrt(syst_alpha_err);
            cout << "syst_alpha: " << syst_alpha*100. << " +/- " << syst_alpha*100.*syst_alpha_err << " (%)"<< endl;
            cout << "<--------- Alpha systematics --------->" << endl;



            Lalph1 += string(" &") + toLatexRounded(alpha,alpha_err);
            Cval   += string(" &") + toLatexRounded(alpha,alpha_err);

            Lalph2 += string(" &") + toLatexRounded(alpha,alpha_err);
            Cval   += string(" &") + toLatexRounded(alpha,alpha_err);


            TH1* NonResonant = NULL;
            if(subNRB2011) {
                NonResonant = (TH1*)hCtrl_SI->Clone("Top/WW/W+Jets (data)");
                NonResonant->SetTitle("Top/WW/W+Jets (data)");
            } else if(subNRB2012) {
                Shape_t& shapeChan_BTag = allShapes.find(selCh[i]+mainHisto+"BTagSB")->second;
                NonResonant = (TH1*)shapeChan_BTag.data->Clone("Top (data)");
                NonResonant->SetTitle("Top (data)");
            } else {
                return;
            }



            TH1* emu_mcnonNRB = (TH1*)shapeCtrl_SI.totalBckg->Clone("emu_mcnonNRB");
            emu_mcnonNRB->Reset();

            for(size_t ibckg=0; ibckg<shapeCtrl_SI.bckg.size(); ibckg++) {
                TString proc(shapeCtrl_SI.bckg[ibckg]->GetTitle());
                if(proc.Contains("ZZ#rightarrow 2l2#nu") ||
                        proc.Contains("WZ#rightarrow 3l#nu") ||
                        proc.Contains("Z+jets") ||
                        proc.Contains("W+jets")
                  ) {
                    emu_mcnonNRB->Add(shapeCtrl_SI.bckg[ibckg], 1);
                }
            }


            double valerr_MCnonNRB_emu;
            double val_MCnonNRB_emu = emu_mcnonNRB->IntegralAndError(1,emu_mcnonNRB->GetXaxis()->GetNbins(),valerr_MCnonNRB_emu);
            cout << "Total MCnonNRB in EMU: " << val_MCnonNRB_emu << endl;


            double valvalerr, valval;
            valval = NonResonant->IntegralAndError(1,NonResonant->GetXaxis()->GetNbins(),valvalerr);

            Cval   += string(" &") + toLatexRounded(valval,valvalerr);

            cout << ctrlCh+AnalysisBins[b]+mainHisto << ": " << valval << endl;
            for(int b=1; b<=NonResonant->GetXaxis()->GetNbins()+1; b++) {
                double val = NonResonant->GetBinContent(b);
                double subtr = emu_mcnonNRB->GetBinContent(b);
                double err = NonResonant->GetBinError(b);
                double newval = (val-subtr)*alpha;
                double newerr = sqrt(pow(err*alpha,2) + pow(val*alpha_err,2));
                NonResonant->SetBinContent(b, newval );
                NonResonant->SetBinError  (b, newerr );
            }
            NonResonant->Scale(DDRescale);

            Double_t valerr;
            Double_t val = NonResonant->IntegralAndError(1,NonResonant->GetXaxis()->GetNbins(),valerr);
            Double_t systError = val*NonResonnantSyst;
            NonResonant->SetBinError(0,systError);//save syst error in underflow bin that is always empty
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
                systError=-1;
            }
            Lyield += string(" &") + toLatexRounded(val,valerr,systError);
            Cval   += string(" &") + toLatexRounded(val,valerr,systError);

            //Clean background collection
            TH1* MCNRB = (TH1*)shapeChan_SI.totalBckg->Clone("MCNRB");
            MCNRB->Reset();

            TH1* MCnonNRB = (TH1*)shapeChan_SI.totalBckg->Clone("MCnonNRB");
            MCnonNRB->Reset();

            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                cout << "proc: " << proc << endl;
                if(proc.Contains("t#bar{t}")
                        || proc.Contains("Single top")
                        || proc.Contains("WW#rightarrow l#nul#nu")
                        || proc.Contains("Z#rightarrow #tau#tau")
                  ) {
                    MCNRB->Add(shapeChan_SI.bckg[ibckg], 1);
                    NonResonant->SetFillColor(shapeChan_SI.bckg[ibckg]->GetFillColor());
                    NonResonant->SetLineColor(shapeChan_SI.bckg[ibckg]->GetLineColor());
                    shapeChan_SI.bckg.erase(shapeChan_SI.bckg.begin()+ibckg);
                    ibckg--;
                } else {
                    MCnonNRB->Add(shapeChan_SI.bckg[ibckg], 1);
                }
            }


            NonResonant->SetFillStyle(1001);
            NonResonant->SetFillColor(393);//592);
            NonResonant->SetLineColor(393);//592);


            //add the background estimate
            shapeChan_SI.bckg.push_back(NonResonant);


            //recompute total background
            shapeChan_SI.totalBckg->Reset();
            for(size_t i=0; i<shapeChan_SI.bckg.size(); i++) {
                shapeChan_SI.totalBckg->Add(shapeChan_SI.bckg[i]);
            }

            val = MCNRB->IntegralAndError(1,MCNRB->GetXaxis()->GetNbins(),valerr);
            cout << "MCNRB: " << val << endl;
            if(val<1E-6) {
                val=0.0;
                valerr=0.0;
            }
            LyieldMC += string(" &") + toLatexRounded(val,valerr);
            Cval     += string(" &") + toLatexRounded(val,valerr);

            if(pFile) {
                fprintf(pFile,"%s\\\\\n", Cval.c_str());
            }
        }
    }

    if(pFile) {
        fprintf(pFile,"\\hline\n");
        fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
        fprintf(pFile,"\n\n\n\n");

        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s|}\\hline\n", Lcol.c_str());
        fprintf(pFile,"%s\\\\\n", Lchan.c_str());
        fprintf(pFile,"%s\\\\\n", Lalph1.c_str());
        fprintf(pFile,"%s\\\\\n", Lalph2.c_str());
        fprintf(pFile,"%s\\\\\n", Lyield.c_str());
        fprintf(pFile,"%s\\\\\n", LyieldMC.c_str());
        fprintf(pFile,"\\hline\n");
        fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
        fclose(pFile);
    }


    cout << "########################## doBackgroundSubtraction ##########################" << endl;
}



void doDYReplacement(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString metHistoForRescale)
{
    cout << "########################## doDYReplacement ##########################" << endl;
    TString DYProcName = "Z+jets";
    TString GammaJetProcName = "Instr. background (data)";
    std::map<TString, double> LowMetIntegral;


    string Ccol   = "\\begin{tabular}{|l|c|c|c|";
    string Cname  = "channel & rescale & yield data & yield mc";
    string Cval   = "";
    FILE* pFile = NULL;
    if(!fast) {
        pFile = fopen("GammaJets.tex","w");
        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Instrumental background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
    }

    FILE* gjFile = NULL;
    gjFile = fopen("DataDY_Yields.log","w");


    //open input file
    TFile* inF = TFile::Open(inFileUrl);
    if( !inF || inF->IsZombie() ) {
        cout << "Invalid file name : " << inFileUrl << endl;
        return;
    }

    TDirectory *pdir = (TDirectory *)inF->Get(DYProcName);
    if(!pdir) {
        printf("Skip Z+Jet estimation because %s directory is missing in root file\n", DYProcName.Data());
        return;
    }

    //all done with input file
    inF->Close();


    //open gamma+jet file
    inF = TFile::Open(DYFile);
    if( !inF || inF->IsZombie() ) {
        cout << "Invalid file name : " << DYFile << endl;
        return;
    }

    pdir = (TDirectory *)inF->Get(GammaJetProcName);
    if(!pdir) {
        printf("Skip Z+Jet estimation because %s directory is missing in Gamma+Jet file\n", GammaJetProcName.Data());
        return;
    }

    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Cval   = selCh[i]+string(" - ")+AnalysisBins[b];
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;

            //find DY background
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                TString proc(shapeChan_SI.bckg[ibckg]->GetTitle());
                if( proc.Contains(DYProcName) ) {

                    double rescaleFactor = 1.0;
                    char buffer[255];
                    sprintf(buffer,"%6.3f",rescaleFactor);
                    Cval   += string(" &") + buffer;

                    Double_t valerr, val = shapeChan_SI.bckg[ibckg]->IntegralAndError(1,shapeChan_SI.bckg[ibckg]->GetXaxis()->GetNbins(),valerr);
                    if(val<1E-6) {
                        val=0.0;
                        valerr=0.0;
                    }
                    string MCTruth = toLatexRounded(val,valerr);
                    double systError = val;

                    //replace DY histogram by G+Jets data
                    int indexcut_ = indexcut;
                    double cutMin=shapeMin;
                    double cutMax=shapeMax;
                    if(indexvbf>=0 && AnalysisBins[b].Contains("vbf")) {
                        printf("use vbf index(%i) for bin %s\n", indexvbf, AnalysisBins[b].Data());
                        indexcut_ = indexvbf;
                        cutMin=shapeMinVBF;
                        cutMax=shapeMaxVBF;
                    }

                    std::cout<<">>>>>  " << selCh[i]+AnalysisBins[b]+"_"+mainHisto <<"\n";
                    fprintf(gjFile,"[%-10.10s - %s]\n", (selCh[i]+AnalysisBins[b]).Data(), mainHisto.Data());

                    TH2* gjets2Dshape  = (TH2*)pdir->Get(selCh[i]+AnalysisBins[b]+"_"+mainHisto);
                    if(!gjets2Dshape)printf("Can't find histo: %s in g+jets template\n",(selCh[i]+AnalysisBins[b]+"_"+mainHisto).Data());

                    TH1* gjets1Dshape  = gjets2Dshape->ProjectionY("tmpName",indexcut_,indexcut_);

                    //apply the cuts
                    for(int x=0; x<=gjets1Dshape->GetXaxis()->GetNbins()+1; x++) {
                        if(gjets1Dshape->GetXaxis()->GetBinCenter(x)<=cutMin || gjets1Dshape->GetXaxis()->GetBinCenter(x)>=cutMax) {
                            gjets1Dshape->SetBinContent(x,0);
                            gjets1Dshape->SetBinError(x,0);
                        }
                    }
                    //gjets1Dshape->Rebin(2);

                    shapeChan_SI.bckg[ibckg]->SetTitle(DYProcName + " (data)");
                    for(int i=1; i<=shapeChan_SI.bckg[ibckg]->GetNbinsX(); i++) {
                        double val = gjets1Dshape->GetBinContent(i);
                        cout << "i: " << i << " \tval: " << val << endl;
                        fprintf(gjFile,"bin: %d, \tval: %f\n", i, val);
                        double err = gjets1Dshape->GetBinError(i);
                        //if(val<0) val=0;
                        //if(AnalysisBins[b]=="eq1jets" && val<0) val=0;
                        double newval = rescaleFactor*val;
                        double newerr = rescaleFactor*err;
                        shapeChan_SI.bckg[ibckg]->SetBinContent(i, newval);
                        shapeChan_SI.bckg[ibckg]->SetBinError  (i, newerr);
                    }
                    delete gjets1Dshape;

                    shapeChan_SI.bckg[ibckg]->Scale(DDRescale);

                    val  = shapeChan_SI.bckg[ibckg]->IntegralAndError(1,shapeChan_SI.bckg[ibckg]->GetXaxis()->GetNbins(),valerr);
                    cout << "Normalization: " << val << " $pm$ " << valerr << "(stat.)" << endl;
                    //systError = std::min(GammaJetSyst * val, fabs(systError - val) ); //syst is difference between G+Jet estimate and MC
                    systError = GammaJetSyst * val; //syst completely coming from G+Jet estimate
                    shapeChan_SI.bckg[ibckg]->SetBinError(0,systError);//save syst error in underflow bin that is always empty
                    if(val<1E-6) {
                        val=0.0;
                        valerr=0.0;
                        systError=-1;
                    }
                    Cval+= string(" &") + toLatexRounded(val,valerr,systError) +" & " + MCTruth;

                    //erase Systematic relatated to DY background
                    std::map<TString,std::vector<std::pair<TString, TH1*> > >::iterator dyvars = shapeChan_SI.bckgVars.find(DYProcName);
                    if(dyvars!=shapeChan_SI.bckgVars.end()) {
                        shapeChan_SI.bckgVars.erase(dyvars);
                    }

                    //recompute total background
                    shapeChan_SI.totalBckg->Reset();
                    for(size_t i=0; i<shapeChan_SI.bckg.size(); i++) {
                        shapeChan_SI.totalBckg->Add(shapeChan_SI.bckg[i]);
                    }

                    if(pFile) {
                        fprintf(pFile,"%s\\\\\n", Cval.c_str());
                    }
                }
            }


        }
    }

    if(pFile) {
        fprintf(pFile,"\\hline\n");
        fprintf(pFile,"\\end{tabular}\n\\end{center}\n\\end{table}\n");
        fprintf(pFile,"\n\n\n\n");
        fclose(pFile);
    }

    //all done with gamma+jet file
    inF->Close();
    cout << "########################## doDYReplacement ##########################" << endl;
}



void doWZSubtraction(std::vector<TString>& selCh,TString ctrlCh,map<TString, Shape_t>& allShapes, TString mainHisto, TString sideBandHisto)
{
    string Ccol   = "\\begin{tabular}{|l|c|c|c|c|";
    string Cname  = "channel & $\\alpha$ measured & $\\alpha$ used & yield data & yield mc";
    string Cval   = "";
    FILE* pFile = NULL;
    if(!fast) {
        pFile = fopen("WZ.tex","w");
        fprintf(pFile,"\\begin{table}[htp]\n\\begin{center}\n\\caption{Non resonant background estimation.}\n\\label{tab:table}\n");
        fprintf(pFile,"%s}\\hline\n", Ccol.c_str());
        fprintf(pFile,"%s\\\\\\hline\n", Cname.c_str());
    }

    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Cval   = selCh[i]+string(" - ")+AnalysisBins[b];

            Shape_t& shapeCtrl_SB = allShapes.find(ctrlCh+AnalysisBins[b]+sideBandHisto)->second;
            Shape_t& shapeCtrl_SI = allShapes.find(ctrlCh+AnalysisBins[b]+mainHisto)->second;
            Shape_t& shapeChan_SB = allShapes.find(selCh[i]+AnalysisBins[b]+sideBandHisto)->second;
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;

            fprintf(pFile,"#############%s:\n",(string(" &")+selCh[i]+string(" - ")+AnalysisBins[b]).Data());
            fprintf(pFile,"MC: em 3leptons=%6.2E  em 2leptons=%6.2E  ll 3leptons=%6.2E  ll 2leptons=%6.2E\n",shapeCtrl_SB.totalBckg->Integral(), shapeCtrl_SI.totalBckg->Integral(), shapeChan_SB.totalBckg->Integral(), shapeChan_SI.totalBckg->Integral());

            TH1* histo1=NULL, *histo2=NULL, *histo3=NULL, *histo4=NULL;
            for(size_t ibckg=0; ibckg<shapeCtrl_SB.bckg.size(); ibckg++) {
                if(TString(shapeCtrl_SB.bckg[ibckg]->GetTitle()).Contains("WZ"))histo1=shapeCtrl_SB.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeCtrl_SI.bckg.size(); ibckg++) {
                if(TString(shapeCtrl_SI.bckg[ibckg]->GetTitle()).Contains("WZ"))histo2=shapeCtrl_SI.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeChan_SB.bckg.size(); ibckg++) {
                if(TString(shapeChan_SB.bckg[ibckg]->GetTitle()).Contains("WZ"))histo3=shapeChan_SB.bckg[ibckg];
            }
            for(size_t ibckg=0; ibckg<shapeChan_SI.bckg.size(); ibckg++) {
                if(TString(shapeChan_SI.bckg[ibckg]->GetTitle()).Contains("WZ"))histo4=shapeChan_SI.bckg[ibckg];
            }
            fprintf(pFile,"WZ: em 3leptons=%6.2E  em 2leptons=%6.2E  ll 3leptons=%6.2E  ll 2leptons=%6.2E\n",histo1->Integral(), histo2->Integral(), histo3->Integral(), histo4->Integral());

            double Num, Denom, NumError, DenomError;
            Num = histo4->IntegralAndError(1,histo4->GetXaxis()->GetNbins(),NumError);
            Denom = histo3->IntegralAndError(1,histo3->GetXaxis()->GetNbins(),DenomError);
            double ratio = Num/Denom;
            double ratio_err  = sqrt(pow(Num*DenomError,2) + pow(Denom*NumError,2))/ pow(Denom,2);
            double ratio_syst = fabs(histo3->Integral() - shapeChan_SB.totalBckg->Integral())/shapeChan_SB.totalBckg->Integral();
            fprintf(pFile,"Ratio = %s\n",toLatexRounded(ratio,ratio_err,ratio_syst).c_str());

            Double_t valerr;
            Double_t val = histo4->IntegralAndError(1,histo4->GetXaxis()->GetNbins(),valerr);

            Double_t valerr2, valsyst2;
            Double_t val2 = shapeChan_SB.totalBckg->IntegralAndError(1,shapeChan_SB.data->GetXaxis()->GetNbins(),valerr2);
            valerr2= sqrt(pow(valerr2*ratio,2) + pow(val2*ratio_err,2) );
            valsyst2 = val2*ratio_syst;
            val2=val2*ratio;
            fprintf(pFile,"WZ (MC closure test): %s --> %s\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val2,valerr2,valsyst2).c_str());

            bool noDataObserved=false;
            Double_t valerr3, valsyst3;
            Double_t val3 = shapeChan_SB.data->IntegralAndError(1,shapeChan_SB.data->GetXaxis()->GetNbins(),valerr3);
            if(val3<=0) {
                noDataObserved=true;
                val3=1.0;
                valerr3=1.0;
            }
            valerr3= sqrt(pow(valerr3*ratio,2) + pow(val3*ratio_err,2) );
            valsyst3 = val3*ratio_syst;
            val3=val3*ratio;
            if(!noDataObserved) {
                fprintf(pFile,"WZ (from data)      : %s --> %s\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val3,valerr3,valsyst3).c_str());
            } else {
                fprintf(pFile,"WZ (from data)      : %s --> smaller than %s (because no data was observed in 3dlepton SideBand--> assume 1+-1 observed data for rescale)\n",toLatexRounded(val,valerr).c_str(), toLatexRounded(val3,valerr3,valsyst3).c_str());
            }


        }
    }

    if(pFile) {
        fclose(pFile);
    }
}

void BlindData(std::vector<TString>& selCh, map<TString, Shape_t>& allShapes, TString mainHisto, bool addSignal)
{
    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Shape_t& shapeChan_SI = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
            shapeChan_SI.data->Reset();
            shapeChan_SI.data->Add(shapeChan_SI.totalBckg,1);
            if(addSignal) {
                for(unsigned int s=0; s<shapeChan_SI.signal.size(); s++) {
                    shapeChan_SI.data->Add(shapeChan_SI.signal[s], 1);
                }
            }
        }
    }
}




void RescaleForInterference(std::vector<TString>& selCh,map<TString, Shape_t>& allShapes, TString mainHisto)
{
    TString massStr("");
    if(mass>0)massStr += mass;
    for(size_t i=0; i<selCh.size(); i++) {
        for(size_t b=0; b<AnalysisBins.size(); b++) {
            Shape_t& shape  = allShapes.find(selCh[i]+AnalysisBins[b]+mainHisto)->second;
            //signals
            size_t nsignal=shape.signal.size();
            for(size_t isignal=0; isignal<nsignal; isignal++) {
                TString proc(((TH1D*)shape.signal[isignal])->GetTitle());
                if(mass>0 && !proc.Contains(massStr))continue;
                if(!proc.Contains("ggH"))continue;
                if(mass<400)continue;
                double scaleFactor = 1.45 - 0.00218 * mass + 0.000002625 * mass * mass;
                double scaleFactorDown = 1.0;
                double scaleFactorUp = 1 + (scaleFactor-1)*2;
                printf("Scale Factor : %f [%f,%f] applied on %s\n",scaleFactor, scaleFactorDown, scaleFactorUp, proc.Data());

                ((TH1D*)shape.signal[isignal])->Scale(scaleFactor);
                std::vector<std::pair<TString, TH1*> >& vars = shape.signalVars[proc];
                for(size_t v=0; v<vars.size(); v++) {
                    ((TH1D*)vars[v].second)->Scale(scaleFactor);
                }

                //((TH1D*)shape.signal[isignal])->SetBinContent(0,scaleFactor);

                TH1* down = (TH1D*)shape.signal[isignal]->Clone(proc+"interf_ggHDown");
                down->Scale(scaleFactorDown);
                TH1* up   = (TH1D*)shape.signal[isignal]->Clone(proc+"interf_ggHUp"  );
                up  ->Scale(scaleFactorUp  );
                vars.push_back(std::make_pair("_interf_ggHDown", down) );
                vars.push_back(std::make_pair("_interf_ggHUp"  , up  ) );
            }
        }
    }
}




