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


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;




namespace METUtils {

  double transverseMass(LorentzVector &visible, LorentzVector &invisible, bool assumeSameMass);

}
#endif /* defined(____METUtils__) */
