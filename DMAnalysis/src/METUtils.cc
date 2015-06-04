//
//  METUtils.cpp
//
//
//  Created by RENJIE WANG on 2/23/15.
//
//

#include "llvvAnalysis/DMAnalysis/interface/METUtils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom.h"

using namespace std;


namespace METUtils {

double transverseMass(LorentzVector &visible, LorentzVector &invisible, bool assumeSameMass)
{
    if(assumeSameMass) {
        LorentzVector sum=visible+invisible;
        double tMass = TMath::Power(TMath::Sqrt(TMath::Power(visible.pt(),2)+pow(visible.mass(),2))+TMath::Sqrt(TMath::Power(invisible.pt(),2)+pow(visible.mass(),2)),2);
        tMass-=TMath::Power(sum.pt(),2);
        return TMath::Sqrt(tMass);
    } else {
        double dphi=fabs(deltaPhi(invisible.phi(),visible.phi()));
        return TMath::Sqrt(2*invisible.pt()*visible.pt()*(1-TMath::Cos(dphi)));
    }
    return -1;
}

}
