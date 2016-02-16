//
//  BTagUtils.cpp
//
//
//  Created by RENJIE WANG on 11/24/14.
//
//

#include "llvvAnalysis/DMAnalysis/interface/BTagUtils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom.h"

using namespace std;

//
BTagUtils::BTagUtils(const edm::ParameterSet &runProcess)
{
    std::vector<std::string> BtagEffFiles = runProcess.getParameter<std::vector<std::string> >("BtagEffFiles");
    for(size_t ifile=0; ifile<BtagEffFiles.size(); ifile++) {

        TString File(BtagEffFiles[ifile].c_str());
        gSystem->ExpandPathName(File);
        TFile *btagFile=TFile::Open(File);

        if(btagFile) {
            cout << "[BTagUtils] retrieving b-tagging efficiency from: " << File << endl;
            std::vector<TString> BtagKeys;
            BtagKeys.push_back("CSVL");
            BtagKeys.push_back("CSVM");
            BtagKeys.push_back("CSVT");

            std::vector<TString> typeseff;
            typeseff.push_back("b_eff");
            typeseff.push_back("c_eff");
            typeseff.push_back("udsg_eff");

            for(size_t itag=0; itag<BtagKeys.size(); itag++) {
                for(size_t ieff=0; ieff<typeseff.size(); ieff++) {
                    TString key = BtagKeys[itag]+"/"+typeseff[ieff];
                    cout << "key: " << key << endl;
                    TH1F *h = (TH1F *) btagFile->Get(key);
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    btageff1DH_[key] = h;
                }//ieff
            }//itag

        }
        btagFile->Close();
        cout << "[BTagUtils] close file: " << File << endl;

    }
}


//
double BTagUtils::getBTagEff(double pt, TString key)
{
    double btag_eff = 1.;

    TH1F* beff_h = btageff1DH_[key];
    if(beff_h==0) cout << "cannot find hist: " << key << endl;
    if(pt < 20) pt=20;
    if(pt > 500) pt=500;
    int binx = beff_h->GetXaxis()->FindBin(pt);
    btag_eff = beff_h->GetBinContent(binx);

    return btag_eff;
}


//
int BTagUtils::getptbin_for_btag(float pt)
{
    if(pt<30) return 0;
    else if(pt<40) return 1;
    else if(pt<50) return 2;
    else if(pt<60) return 3;
    else if(pt<70) return 4;
    else if(pt<80) return 5;
    else if(pt<100) return 6;
    else if(pt<120) return 7;
    else if(pt<160) return 8;
    else if(pt<210) return 9;
    else if(pt<260) return 10;
    else if(pt<320) return 11;
    else if(pt<400) return 12;
    else if(pt<500) return 13;
    else if(pt<600) return 14;
    else return 15;

}

//
int BTagUtils::get_eta_bin_jet(float eta)
{
    eta = fabs(eta);
    if(eta<0.9) return 0;
    else if(eta<1.2) return 1;
    else if(eta<2.1) return 2;
    else if(eta<2.4) return 3;
    else return -1;
}


//
std::pair<double,double> BTagUtils::getBTagSF(double pt, double eta, int flavorJet, TString btag)
{
    double x = pt;
    eta = fabs(eta);

    double SF=1., SFerr=0.;

    if(x<20) x=20;
    if(x>800) x= 800;

    if(abs(flavorJet)==5 || abs(flavorJet) == 4) { //for b or c.
        if(btag=="CSVL")      SF = 0.997942*((1.+(0.00923753*x))/(1.+(0.0096119*x)));
        else if(btag=="CSVM") SF = (0.938887+(0.00017124*x))+(-2.76366e-07*(x*x));
        else if(btag=="CSVT") SF = (0.927563+(1.55479e-05*x))+(-1.90666e-07*(x*x));

        int ptbin = getptbin_for_btag(pt);
        if(btag=="CSVL")      SFerr = SFb_error_CSVL[ptbin];
        else if(btag=="CSVM") SFerr = SFb_error_CSVM[ptbin];
        else if(btag=="CSVT") SFerr = SFb_error_CSVT[ptbin];

        if(x>800 || x<20 ) SFerr *= 2;
        if(abs(flavorJet) == 4 ) SFerr *= 2;

    } else { ///SFlight
        float max=999.;
        float min=-999.;

        if(btag=="CSVL") {
            if(eta<=0.5) {
                SF = ((1.01177+(0.0023066*x))+(-4.56052e-06*(x*x)))+(2.57917e-09*(x*(x*x)));
                max = ((1.04582+(0.00290226*x))+(-5.89124e-06*(x*x)))+(3.37128e-09*(x*(x*x)));
                min = ((0.977761+(0.00170704*x))+(-3.2197e-06*(x*x)))+(1.78139e-09*(x*(x*x)));
            } else if(eta<=1.0) {
                SF = ((0.975966+(0.00196354*x))+(-3.83768e-06*(x*x)))+(2.17466e-09*(x*(x*x)));
                max = ((1.00683+(0.00246404*x))+(-4.96729e-06*(x*x)))+(2.85697e-09*(x*(x*x)));
                min = ((0.945135+(0.00146006*x))+(-2.70048e-06*(x*x)))+(1.4883e-09*(x*(x*x)));
            } else if(eta<=1.5) {
                SF = ((0.93821+(0.00180935*x))+(-3.86937e-06*(x*x)))+(2.43222e-09*(x*(x*x)));
                max = ((0.964787+(0.00219574*x))+(-4.85552e-06*(x*x)))+(3.09457e-09*(x*(x*x)));
                min = ((0.911657+(0.00142008*x))+(-2.87569e-06*(x*x)))+(1.76619e-09*(x*(x*x)));
            } else if(eta<=2.4) {
                SF = ((1.00022+(0.0010998*x))+(-3.10672e-06*(x*x)))+(2.35006e-09*(x*(x*x)));
                max = ((1.03039+(0.0013358*x))+(-3.89284e-06*(x*x)))+(3.01155e-09*(x*(x*x)));
                min = ((0.970045+(0.000862284*x))+(-2.31714e-06*(x*x)))+(1.68866e-09*(x*(x*x)));
            }
        } else if(btag=="CSVM") {
            if(eta<=0.8) {
                SF = ((1.07541+(0.00231827*x))+(-4.74249e-06*(x*x)))+(2.70862e-09*(x*(x*x)));
                max = ((1.18638+(0.00314148*x))+(-6.68993e-06*(x*x)))+(3.89288e-09*(x*(x*x)));
                min = ((0.964527+(0.00149055*x))+(-2.78338e-06*(x*x)))+(1.51771e-09*(x*(x*x)));
            } else if(eta<=1.6) {
                SF = ((1.05613+(0.00114031*x))+(-2.56066e-06*(x*x)))+(1.67792e-09*(x*(x*x)));
                max = ((1.16624+(0.00151884*x))+(-3.59041e-06*(x*x)))+(2.38681e-09*(x*(x*x)));
                min = ((0.946051+(0.000759584*x))+(-1.52491e-06*(x*x)))+(9.65822e-10*(x*(x*x)));
            }  else if(eta<=2.4) {
                SF = ((1.05625+(0.000487231*x))+(-2.22792e-06*(x*x)))+(1.70262e-09*(x*(x*x)));
                max = ((1.15575+(0.000693344*x))+(-3.02661e-06*(x*x)))+(2.39752e-09*(x*(x*x)));
                min = ((0.956736+(0.000280197*x))+(-1.42739e-06*(x*x)))+(1.0085e-09*(x*(x*x)));
            }
        } else if(btag=="CSVT") {
            if(eta<=2.4) {
                SF = ((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
                max = ((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)));
                min = ((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)));
            }
        }

        SFerr = fabs(max-SF)>fabs(min-SF)? fabs(max-SF):fabs(min-SF);
    }

    std::pair<float,float> SF_SFerr(SF,SFerr);
    return SF_SFerr;
}


//
std::pair<double,double> BTagUtils::getBTagWeight(bool istag, double pt, double eta, int flavorJet, TString btag, TString key)
{
    float mcTag = 1.;
    float mcNoTag = 1.;
    float dataTag = 1.;
    float dataNoTag = 1.;

    float err1 = 0;
    float err2 = 0;
    float err3 = 0;
    float err4 = 0;

    double SF = getBTagSF(pt,eta,flavorJet,btag).first;
    double SFerr = getBTagSF(pt,eta,flavorJet,btag).second;
    double eff = getBTagEff(pt,key);


    if(istag) {
        mcTag *= eff;
        dataTag *= eff*SF;

        if(abs(flavorJet)==5 || abs(flavorJet)==4)  err1 += SFerr/SF; ///correlated for b/c
        else err3 += SFerr/SF; //correlated for light


    } else {
        mcNoTag *= (1- eff);
        dataNoTag *= (1- eff*SF);

        if(abs(flavorJet)==5 || abs(flavorJet)==4 ) err2 += (-eff*SFerr)/(1-eff*SF); /// /correlated for b/c
        else err4 +=  (-eff*SFerr)/(1-eff*SF);  ////correlated for light

    }

    double wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag );
    double wtbtagErr = sqrt( pow(err1+err2,2) + pow( err3 + err4,2)) * wtbtag;  ///un-correlated for b/c and light


    std::pair<double,double> BTagWeights(wtbtag,wtbtagErr);
    return BTagWeights;
}



//
double BTagUtils::getNewBTagWeight(bool istag, double pt, double SF, TString btag, TString key)
{
    float mcTag = 1.;
    float mcNoTag = 1.;
    float dataTag = 1.;
    float dataNoTag = 1.;

    double eff = getBTagEff(pt,key);

    if(istag) {
        mcTag *= eff;
        dataTag *= eff*SF;
    } else {
        mcNoTag *= (1- eff);
        dataNoTag *= (1- eff*SF);
    }

    double wtbtag = (dataNoTag * dataTag ) / ( mcNoTag * mcTag );

    return wtbtag;
}


//
BTagUtils::~BTagUtils()
{
}






