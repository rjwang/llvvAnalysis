//
//  BTagUtils.h
//
//
//  Created by RENJIE WANG on 11/24/14.
//
//

#ifndef ____BTagUtils__
#define ____BTagUtils__

#include <iostream>
#include <TString.h>
#include <map>

#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"

#include "Math/LorentzVector.h"
#include "TVector2.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class BTagUtils {

public:

    BTagUtils(const edm::ParameterSet &runProcess);
    double getBTagEff(double pt, TString key);
    std::pair<double,double> getBTagSF(double pt, double eta, int flavorJet, TString btag);
    std::pair<double,double> getBTagWeight(bool istag, double pt, double eta, int flavorJet, TString btag, TString key);

    int getptbin_for_btag(float pt);
    int get_eta_bin_jet(float eta);

    ~BTagUtils();

private:

    std::map<TString, TH1F *> btageff1DH_;


    double SFb_error_CSVL[100] = {
        0.033299,
        0.0146768,
        0.013803,
        0.0170145,
        0.0166976,
        0.0137879,
        0.0149072,
        0.0153068,
        0.0133077,
        0.0123737,
        0.0157152,
        0.0175161,
        0.0209241,
        0.0278605,
        0.0346928,
        0.0350099
    };

    double SFb_error_CSVM[100] = {
        0.0415707,
        0.0204209,
        0.0223227,
        0.0206655,
        0.0199325,
        0.0174121,
        0.0202332,
        0.0182446,
        0.0159777,
        0.0218531,
        0.0204688,
        0.0265191,
        0.0313175,
        0.0415417,
        0.0740446,
        0.0596716
    };

    double SFb_error_CSVT[100] = {
        0.0515703,
        0.0264008,
        0.0272757,
        0.0275565,
        0.0248745,
        0.0218456,
        0.0253845,
        0.0239588,
        0.0271791,
        0.0273912,
        0.0379822,
        0.0411624,
        0.0786307,
        0.0866832,
        0.0942053,
        0.102403
    };

};


#endif /* defined(____BTagUtils__) */
