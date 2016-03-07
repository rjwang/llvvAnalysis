//
//  DMPhysicsEvent.cpp
//
//
//  Created by RENJIE WANG on 2/17/15.
//
//

#include "llvvAnalysis/DMAnalysis/interface/DMPhysicsEvent.h"

using namespace std;


//
PhysicsEvent_t getPhysicsEventFrom(DataEvtSummary_t &ev)
{
    PhysicsEvent_t phys;

    phys.run=ev.run;
    phys.event=ev.event;
    phys.lumi=ev.lumi;
    phys.nvtx = ev.nvtx;

    // Leptons
    size_t nlep(0);
    for(Int_t i=0; i<ev.mn; i++) {
        LorentzVector P4(ev.mn_px[i],ev.mn_py[i],ev.mn_pz[i],ev.mn_en[i]);
        if(P4.pt()>0) {
            phys.leptons.push_back(PhysicsObject_Lepton(P4, ev.mn_id[i]));
            phys.leptons[nlep].setLeptonIDInfo(ev.mn_IsLoose[i], ev.mn_IsTight[i], ev.mn_IsMedium[i], ev.mn_IsSoft[i],ev.mn_IsHighPt[i],
                                               ev.en_passVeto[i], ev.en_passLoose[i], ev.en_passMedium[i], ev.en_passTight[i],
                                               ev.ta_dm[i]
                                              );
            phys.leptons[nlep].setLeptonIsoInfo(ev.mn_pileupIsoR03[i],ev.mn_chargedIsoR03[i],ev.mn_photonIsoR03[i],ev.mn_neutralHadIsoR03[i],
                                                ev.en_pileupIso[i],ev.en_chargedIso[i],ev.en_photonIso[i],ev.en_neutralHadIso[i],ev.en_relIsoWithEA[i],
                                                ev.ta_IsLooseIso[i], ev.ta_IsMediumIso[i], ev.ta_IsTightIso[i]
                                               );

            nlep++;
        }
    }

    for(Int_t i=0; i<ev.en; i++) {
        LorentzVector P4(ev.en_px[i],ev.en_py[i],ev.en_pz[i],ev.en_en[i]);
        if(P4.pt()>0) {
            phys.leptons.push_back(PhysicsObject_Lepton(P4, ev.en_id[i]));
            phys.leptons[nlep].setLeptonIDInfo(ev.mn_IsLoose[i], ev.mn_IsTight[i], ev.mn_IsMedium[i], ev.mn_IsSoft[i],ev.mn_IsHighPt[i],
                                               ev.en_passVeto[i], ev.en_passLoose[i], ev.en_passMedium[i], ev.en_passTight[i],
                                               ev.ta_dm[i]
                                              );
            phys.leptons[nlep].setLeptonIsoInfo(ev.mn_pileupIsoR03[i],ev.mn_chargedIsoR03[i],ev.mn_photonIsoR03[i],ev.mn_neutralHadIsoR03[i],
                                                ev.en_pileupIso[i],ev.en_chargedIso[i],ev.en_photonIso[i],ev.en_neutralHadIso[i],ev.en_relIsoWithEA[i],
                                                ev.ta_IsLooseIso[i], ev.ta_IsMediumIso[i], ev.ta_IsTightIso[i]
                                               );


            nlep++;
        }
    }


    for(Int_t i=0; i<ev.ta; i++) {
        LorentzVector P4(ev.ta_px[i],ev.ta_py[i],ev.ta_pz[i],ev.ta_en[i]);
        if(P4.pt()>0) {
            phys.leptons.push_back(PhysicsObject_Lepton(P4, ev.ta_id[i]));
            phys.leptons[nlep].setLeptonIDInfo(ev.mn_IsLoose[i], ev.mn_IsTight[i], ev.mn_IsMedium[i], ev.mn_IsSoft[i],ev.mn_IsHighPt[i],
                                               ev.en_passVeto[i], ev.en_passLoose[i], ev.en_passMedium[i], ev.en_passTight[i],
                                               ev.ta_dm[i]
                                              );
            phys.leptons[nlep].setLeptonIsoInfo(ev.mn_pileupIsoR03[i],ev.mn_chargedIsoR03[i],ev.mn_photonIsoR03[i],ev.mn_neutralHadIsoR03[i],
                                                ev.en_pileupIso[i],ev.en_chargedIso[i],ev.en_photonIso[i],ev.en_neutralHadIso[i],ev.en_relIsoWithEA[i],
                                                ev.ta_IsLooseIso[i], ev.ta_IsMediumIso[i], ev.ta_IsTightIso[i]
                                               );

            nlep++;
        }
    }





    // MET
    phys.met 	 = LorentzVector( ev.met_pt*cos(ev.met_phi), ev.met_pt*sin(ev.met_phi), 0, ev.met_pt );
    phys.metNoHF = LorentzVector( ev.metNoHF_pt*cos(ev.metNoHF_phi), ev.metNoHF_pt*sin(ev.metNoHF_phi), 0, ev.metNoHF_pt );




    // Jet
    size_t njet(0);
    for(Int_t i=0; i<ev.jet; i++) {
        LorentzVector P4( ev.jet_px[i],ev.jet_py[i],ev.jet_pz[i],ev.jet_en[i] );
        if(P4.pt()>0) {
            phys.jets.push_back( PhysicsObject_Jet( P4,ev.jet_puId[i],ev.jet_PFLoose[i],ev.jet_PFTight[i] ) );
            phys.jets[i].setBtagInfo(ev.jet_btag0[i],ev.jet_btag1[i],ev.jet_btag2[i],ev.jet_btag3[i],ev.jet_btag4[i],ev.jet_btag5[i],ev.jet_btag6[i],ev.jet_btag7[i]);
            phys.jets[i].setGenInfo(ev.jet_partonFlavour[i], ev.jet_genpt[i]);

            njet++;
        }
    }


    //generator level particles
    for(Int_t ipart=0; ipart<ev.nmcparticles; ipart++) {
        LorentzVector p4(ev.mc_px[ipart],ev.mc_py[ipart],ev.mc_pz[ipart],ev.mc_en[ipart]);
        if(ev.mc_status[ipart]==2 && abs(ev.mc_id[ipart])==15) phys.genleptons.push_back( PhysicsObject(p4,ev.mc_id[ipart]) );
        if(ev.mc_status[ipart]!=1) continue;
        switch( ev.mc_id[ipart] ) {
        case 12:
        case -12:
        case 14:
        case -14:
        case 16:
        case -16: {
            phys.genneutrinos.push_back( PhysicsObject(p4,ev.mc_id[ipart]) );
        }
        break;
        case 2000012:
        case -2000012: {
            phys.genWIMPs.push_back( PhysicsObject(p4,ev.mc_id[ipart]) );
        }
        break;
        case 5000039:
        case -5000039: {
            phys.genGravitons.push_back( PhysicsObject(p4,ev.mc_id[ipart]) );
        }
        break;
        case 11:
        case -11:
        case 13:
        case -13:
            //case 15:
            //case -15:
        {
            phys.genleptons.push_back( PhysicsObject(p4,ev.mc_id[ipart]) );
        }
        break;
        }
    }

    //add generator jets separately
    for(Int_t ipart=0; ipart<ev.nmcparticles; ipart++) {
        if(ev.mc_id[ipart]!=1 || ev.mc_status[ipart]!=0) continue; //genjet id set
        LorentzVector p4(ev.mc_px[ipart],ev.mc_py[ipart],ev.mc_pz[ipart],ev.mc_en[ipart]);

        bool overlap(false);
        //check overlap with previous jets
        for(size_t igj=0; igj<phys.genjets.size(); igj++) {
            if(deltaR(p4,phys.genjets[igj])<0.1) {
                overlap=true;
                break;
            }
        }
        if(overlap) continue;

        //check overlap with any other gen particle
        for(Int_t jpart=0; jpart<ev.nmcparticles; jpart++) {
            if(ev.mc_id[jpart]==1 && ev.mc_status[jpart]==0) continue; // not genjet
            if(ev.mc_status[jpart]!=1) continue; // no intermediator
            LorentzVector jp4(ev.mc_px[jpart],ev.mc_py[jpart],ev.mc_pz[jpart],ev.mc_en[jpart]);
            if(deltaR(p4,jp4)<0.1) {
                overlap=true;
                break;
            }
        }
        if(overlap) continue;

        phys.genjets.push_back( PhysicsObject(p4,ev.mc_id[ipart]) );
    }



    return phys;
}


//
int getDileptonId(int id1, int id2)
{
    if( id1==MUON     && id2==MUON)     return MUMU;
    if( id1==ELECTRON && id2==ELECTRON) return EE;
    if( ( id1==MUON && id2==ELECTRON) || ( id1==ELECTRON && id2==MUON)  ) return EMU;
    return UNKNOWN;
}


bool isDYToLL(int id1, int id2)
{
    if(id1==11 || id2==-11) return true;
    else if(id1==-11 || id2==11) return true;
    else if(id1==13 || id2==-13) return true;
    else if(id1==-13 || id2==13) return true;
    else return false;
}

bool isDYToTauTau(int id1, int id2)
{
    if(id1==15 || id2==-15) return true;
    else if(id1==-15 || id2==15) return true;
    else return false;
}


float getSFfrom2DHist(double xval, double yval, TH2F* h_)
{

    if(h_==NULL) {
        cout << "[getSFfrom2DHist]: empty hist! " << endl;
        return 1;
    }
    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

    int ybins = h_->GetYaxis()->GetNbins();
    if(yval > h_->GetYaxis()->GetBinUpEdge(ybins)    ) yval = h_->GetYaxis()->GetBinUpEdge(ybins);
    if(yval < h_->GetYaxis()->GetBinLowEdge(1)       ) yval = h_->GetYaxis()->GetBinLowEdge(1);

    int binx = h_->GetXaxis()->FindBin(xval);
    int biny = h_->GetYaxis()->FindBin(yval);
    float sf_ = h_->GetBinContent(binx,biny);

    return sf_;
}

float getNLOEWKZZWeight(double trailing_pt)
{
    double ewk_w = 1.0;

    if(     trailing_pt<60.)  ewk_w = 1.-( 4.0/100.);
    else if(trailing_pt<80.)  ewk_w = 1.-( 5.0/100.);
    else if(trailing_pt<100.) ewk_w = 1.-( 6.3/100.);
    else if(trailing_pt<120.) ewk_w = 1.-( 7.6/100.);
    else if(trailing_pt<140.) ewk_w = 1.-( 9.2/100.);
    else if(trailing_pt<160.) ewk_w = 1.-(10.0/100.);
    else if(trailing_pt<180.) ewk_w = 1.-(11.4/100.);
    else if(trailing_pt<200.) ewk_w = 1.-(12.8/100.);
    else if(trailing_pt<220.) ewk_w = 1.-(14.2/100.);
    else if(trailing_pt<240.) ewk_w = 1.-(15.6/100.);
    else if(trailing_pt<260.) ewk_w = 1.-(17.0/100.);
    else if(trailing_pt<280.) ewk_w = 1.-(18.4/100.);
    else if(trailing_pt<300.) ewk_w = 1.-(20.0/100.);
    else if(trailing_pt<320.) ewk_w = 1.-(21.2/100.);
    else if(trailing_pt<340.) ewk_w = 1.-(22.4/100.);
    else if(trailing_pt<360.) ewk_w = 1.-(23.6/100.);
    else if(trailing_pt<380.) ewk_w = 1.-(24.8/100.);
    else                      ewk_w = 1.-(26.0/100.);

    return ewk_w;
}


float kfactor_qqZZ_qcd_dPhi(float GENdPhiZZ)
{

    float k = 0.0;
    k += 1.513834489150 * (abs(GENdPhiZZ)>0.0 && abs(GENdPhiZZ)<=0.1);
    k += 1.541738780180 * (abs(GENdPhiZZ)>0.1 && abs(GENdPhiZZ)<=0.2);
    k += 1.497829632510 * (abs(GENdPhiZZ)>0.2 && abs(GENdPhiZZ)<=0.3);
    k += 1.534956782920 * (abs(GENdPhiZZ)>0.3 && abs(GENdPhiZZ)<=0.4);
    k += 1.478217033060 * (abs(GENdPhiZZ)>0.4 && abs(GENdPhiZZ)<=0.5);
    k += 1.504330859290 * (abs(GENdPhiZZ)>0.5 && abs(GENdPhiZZ)<=0.6);
    k += 1.520626246850 * (abs(GENdPhiZZ)>0.6 && abs(GENdPhiZZ)<=0.7);
    k += 1.507013090030 * (abs(GENdPhiZZ)>0.7 && abs(GENdPhiZZ)<=0.8);
    k += 1.494243156250 * (abs(GENdPhiZZ)>0.8 && abs(GENdPhiZZ)<=0.9);
    k += 1.450536096150 * (abs(GENdPhiZZ)>0.9 && abs(GENdPhiZZ)<=1.0);
    k += 1.460812521660 * (abs(GENdPhiZZ)>1.0 && abs(GENdPhiZZ)<=1.1);
    k += 1.471603622200 * (abs(GENdPhiZZ)>1.1 && abs(GENdPhiZZ)<=1.2);
    k += 1.467700038200 * (abs(GENdPhiZZ)>1.2 && abs(GENdPhiZZ)<=1.3);
    k += 1.422408690640 * (abs(GENdPhiZZ)>1.3 && abs(GENdPhiZZ)<=1.4);
    k += 1.397184022730 * (abs(GENdPhiZZ)>1.4 && abs(GENdPhiZZ)<=1.5);
    k += 1.375593447520 * (abs(GENdPhiZZ)>1.5 && abs(GENdPhiZZ)<=1.6);
    k += 1.391901318370 * (abs(GENdPhiZZ)>1.6 && abs(GENdPhiZZ)<=1.7);
    k += 1.368564350560 * (abs(GENdPhiZZ)>1.7 && abs(GENdPhiZZ)<=1.8);
    k += 1.317884804290 * (abs(GENdPhiZZ)>1.8 && abs(GENdPhiZZ)<=1.9);
    k += 1.314019950800 * (abs(GENdPhiZZ)>1.9 && abs(GENdPhiZZ)<=2.0);
    k += 1.274641749910 * (abs(GENdPhiZZ)>2.0 && abs(GENdPhiZZ)<=2.1);
    k += 1.242346606820 * (abs(GENdPhiZZ)>2.1 && abs(GENdPhiZZ)<=2.2);
    k += 1.244727403840 * (abs(GENdPhiZZ)>2.2 && abs(GENdPhiZZ)<=2.3);
    k += 1.146259351670 * (abs(GENdPhiZZ)>2.3 && abs(GENdPhiZZ)<=2.4);
    k += 1.107804993520 * (abs(GENdPhiZZ)>2.4 && abs(GENdPhiZZ)<=2.5);
    k += 1.042053646740 * (abs(GENdPhiZZ)>2.5 && abs(GENdPhiZZ)<=2.6);
    k += 0.973608545141 * (abs(GENdPhiZZ)>2.6 && abs(GENdPhiZZ)<=2.7);
    k += 0.872169942668 * (abs(GENdPhiZZ)>2.7 && abs(GENdPhiZZ)<=2.8);
    k += 0.734505279177 * (abs(GENdPhiZZ)>2.8 && abs(GENdPhiZZ)<=2.9);
    k += 1.163152837230 * (abs(GENdPhiZZ)>2.9 && abs(GENdPhiZZ)<=3.1416);

    if (k==0.0) return 1.1; // if something goes wrong return inclusive k-factor
    else return k;
}
