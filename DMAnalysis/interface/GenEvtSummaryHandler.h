#ifndef genevtsummaryhandler_h
#define genevtsummaryhandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>

#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include "TTree.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;

#define MAXPARTICLES 50
#define MAXMCPARTICLES 250

struct GenEvtSummary_t
{

  Int_t run,lumi,event;

  //gen level event
  int nmcparticles;
  float mc_px[MAXMCPARTICLES],mc_py[MAXMCPARTICLES],mc_pz[MAXMCPARTICLES],mc_en[MAXMCPARTICLES];
  int mc_id[MAXPARTICLES], mc_status[MAXPARTICLES], mc_mom[MAXPARTICLES];

};

class GenEvtSummaryHandler{
 public:
  //
  GenEvtSummaryHandler();
  ~GenEvtSummaryHandler();

  //current event
  GenEvtSummary_t evSummary_;
  GenEvtSummary_t &getEvent() { return evSummary_; }

  //write mode
  bool initTree(TTree *t);
  void resetStruct();
  //void fillTree();

 private:
  //the tree
  TTree *t_;
};

#endif
