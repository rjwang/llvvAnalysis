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


float getSFfrom2DHist(double xval, double yval, TH2F* h_){

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

