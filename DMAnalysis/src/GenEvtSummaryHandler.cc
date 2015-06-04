#include "llvvAnalysis/DMAnalysis/interface/GenEvtSummaryHandler.h"

using namespace std;

//
GenEvtSummaryHandler::GenEvtSummaryHandler()
{
}

//
bool GenEvtSummaryHandler::initTree(TTree *t)
{
  if(t==0) return false;
  t_ = t;

  //event info
  t_->Branch("run",        &evSummary_.run,        "run/I");
  t_->Branch("lumi",       &evSummary_.lumi,       "lumi/I");
  t_->Branch("event",      &evSummary_.event,      "event/I");

  //mc truth
  t_->Branch("nmcparticles", &evSummary_.nmcparticles, "nmcparticles/I");
  t_->Branch("mc_px", evSummary_.mc_px, "mc_px[nmcparticles]/F");
  t_->Branch("mc_py", evSummary_.mc_py, "mc_py[nmcparticles]/F");
  t_->Branch("mc_pz", evSummary_.mc_pz, "mc_pz[nmcparticles]/F");
  t_->Branch("mc_en", evSummary_.mc_en, "mc_en[nmcparticles]/F");
  t_->Branch("mc_id", evSummary_.mc_id, "mc_id[nmcparticles]/I");
  t_->Branch("mc_status", evSummary_.mc_status, "mc_status[nmcparticles]/I");
  t_->Branch("mc_mom", evSummary_.mc_mom, "mc_mom[nmcparticles]/I");

  return true;
}

//
void GenEvtSummaryHandler::resetStruct()
{
  evSummary_.nmcparticles=0;
  evSummary_.run=0;    evSummary_.lumi=0;   evSummary_.event=0;
}

//
/*
void GenEvtSummaryHandler::fillTree()
{
  if(t_) t_->Fill();
}
*/

//
GenEvtSummaryHandler::~GenEvtSummaryHandler()
{
}
