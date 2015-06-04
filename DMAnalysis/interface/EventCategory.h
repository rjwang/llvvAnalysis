//
//  EventCategory.h
//
//
//  Created by RENJIE WANG on 3/3/15.
//
//

#ifndef ____EventCategory__
#define ____EventCategory__

#include <iostream>
#include "llvvAnalysis/DMAnalysis/interface/DMPhysicsEvent.h"
#include "llvvAnalysis/DMAnalysis/interface/DataEvtSummaryHandler.h"

class EventCategory {
public:
    int mode;
    int NStates;
    TString* EvtCategoryLabel;


    /// Constructor
    EventCategory(int mode=0);

    /// Destructor
    virtual ~EventCategory();

    int GetLabelSize() {
        return NStates;
    }

    int Get(const PhysicsEvent_t& phys, PhysicsObjectJetCollection* JetsP4=NULL);
    TString GetLabel(int CategoryType);
    TString GetLabel(const PhysicsEvent_t& phys);


private:
    int NJets;

};


#endif /* defined(____EventCategory__) */
