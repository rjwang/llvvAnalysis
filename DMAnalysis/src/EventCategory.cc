//
//  EventCategory.cpp
//
//
//  Created by RENJIE WANG on 3/3/15.
//
//

#include "llvvAnalysis/DMAnalysis/interface/EventCategory.h"


using namespace std;

EventCategory::EventCategory(int mode_)
{
    mode = mode_;
    if(mode==1) {
        NStates = 3;
        EvtCategoryLabel = new TString[NStates];
        EvtCategoryLabel[0] = "eq0jets";
        EvtCategoryLabel[1] = "eq1jets";
        EvtCategoryLabel[2] = "geq2jets";
    } else if(mode==2) {
        NStates = 2;
        EvtCategoryLabel = new TString[NStates];
        EvtCategoryLabel[0] = "lesq1jets";
        EvtCategoryLabel[1] = "geq2jets";
    } else {
        mode = 0;
        NStates = 1;
        EvtCategoryLabel = new TString[NStates];
        EvtCategoryLabel[0] = "";
    }
}

EventCategory::~EventCategory() {}

//
int EventCategory::Get(const PhysicsEvent_t& phys, PhysicsObjectJetCollection* JetsP4)
{
    if(!JetsP4) {
        printf("NoJetCollection given to EventCategory::Get\n");
        exit(0);
    }

    NJets = 0;
    PhysicsObjectJetCollection jets = *JetsP4;

    for(size_t ijet=0; ijet<jets.size(); ijet++) {
        if(jets[ijet].pt()<=30)continue;
        NJets++;
    }

    switch(mode) {
    case 1: {
        if(NJets>=2) return 2;
        if(NJets==1) return 1;
        return 0;
    }
    break;
    case 2: {
        if(NJets>=2) return 1;
        return 0;
    }
    break;
    case 0:
    default: {
        return 0;
    }
    break;
    }
}

//
TString EventCategory::GetLabel(int CategoryType)
{
    if(mode==0) {
        if(CategoryType<=2) return EvtCategoryLabel[0];
        return EvtCategoryLabel[1];
    } else {
        return EvtCategoryLabel[CategoryType];
    }
}

//
TString EventCategory::GetLabel(const PhysicsEvent_t& phys)
{
    return GetLabel(Get(phys));
}
