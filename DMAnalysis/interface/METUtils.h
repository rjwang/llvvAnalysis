//
//  METUtils.h
//
//
//  Created by RENJIE WANG on 2/23/15.
//
//

#ifndef ____METUtils__
#define ____METUtils__



#include <utility>
#include <vector>

#include "Math/LorentzVector.h"
#include "TVector2.h"

#include "CondFormats/JetMETObjects/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "llvvAnalysis/DMAnalysis/interface/DMPhysicsEvent.h"
#include "llvvAnalysis/DMAnalysis/interface/DataEvtSummaryHandler.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;




namespace METUtils {

  double transverseMass(LorentzVector &visible, LorentzVector &invisible, bool assumeSameMass);
  double response(LorentzVector &Z, LorentzVector &MET);
  PhysicsObject_Jet smearedJet(const LorentzVector &origJet, double genJetPt, int mode=0);


  enum UncertaintyVariations { JER, JER_UP, JER_DOWN, JES_UP, JES_DOWN,UMET_UP,UMET_DOWN,LES_UP,LES_DOWN};
  void computeVariation(PhysicsObjectJetCollection& jets,
                        PhysicsObjectLeptonCollection &leptons,
                        LorentzVector& met,
                        std::vector<PhysicsObjectJetCollection>& jetsVar,
                        LorentzVectorCollection& metsVar,
                        JetCorrectionUncertainty *jecUnc);
  LorentzVector applyMETXYCorr(LorentzVector met, bool isMC, int nvtx);
}
#endif /* defined(____METUtils__) */
