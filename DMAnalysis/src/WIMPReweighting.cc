//
//  WIMPReweighting.cpp
//
//
//  Created by RENJIE WANG on 3/8/15.
//
//

#include "llvvAnalysis/DMAnalysis/interface/WIMPReweighting.h"

using namespace std;

//
WIMPReweighting::WIMPReweighting(const edm::ParameterSet &runProcess)
{
    std::vector<std::string> WIMPFiles = runProcess.getParameter<std::vector<std::string> >("genwimpweights");
    std::vector<std::string> MCFiles = runProcess.getParameter<std::vector<std::string> >("MCweights");

    for(size_t ifile=0; ifile<WIMPFiles.size(); ifile++) {

        TString File(WIMPFiles[ifile].c_str());
	if(!File.Contains(".root")) continue;
        gSystem->ExpandPathName(File);
        TFile *wimpFile=TFile::Open(File);

        if(wimpFile) {
            cout << "[WIMPReweighting] retrieving WIMP weights from: " << File << endl;
            std::vector<TString> WimpKeys1D;
            WimpKeys1D.push_back("genmet");
	    WimpKeys1D.push_back("pt_chichi");
	    //WimpKeys1D.push_back("pileup_weights");
            for(size_t itag=0; itag<WimpKeys1D.size(); itag++) {
                    TString key = WimpKeys1D[itag];
                    cout << "key: " << key << endl;
                    TH1F *h = (TH1F *) wimpFile->Get(key);
		    if(h==0) continue;
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    wimpWeights1DH_[key] = h;
            }//itag

            std::vector<TString> WimpKeys2D;
            WimpKeys2D.push_back("dphi_vs_met");

            for(size_t itag=0; itag<WimpKeys2D.size(); itag++) {
                    TString key = WimpKeys2D[itag];
                    cout << "key: " << key << endl;
                    TH2F *h = (TH2F *) wimpFile->Get(key);
		    if(h==0) continue;
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    wimpWeights2DH_[key] = h;
            }//itag


	    cout << "[WIMPReweighting] close file: " << File << endl;
        }
        wimpFile->Close();
    }


    for(size_t ifile=0; ifile<MCFiles.size(); ifile++) {

        TString File(MCFiles[ifile].c_str());
	if(!File.Contains(".root")) continue;
        gSystem->ExpandPathName(File);
        TFile *wimpFile=TFile::Open(File);

        if(wimpFile) {
            cout << "[WIMPReweighting] retrieving MC weights from: " << File << endl;
            std::vector<TString> WimpKeys1D;
	    WimpKeys1D.push_back("pileup_weights");
            for(size_t itag=0; itag<WimpKeys1D.size(); itag++) {
                    TString key = WimpKeys1D[itag];
                    cout << "key: " << key << endl;
                    TH1F *h = (TH1F *) wimpFile->Get(key);
		    if(h==0) continue;
                    h->SetDirectory(0); //THIS IS IMPORTANT FOR TH1 Weight File!!!
                    wimpWeights1DH_[key] = h;
            }//itag

	    cout << "[WIMPReweighting] close file: " << File << endl;
        }
        wimpFile->Close();
    }







}


//
double WIMPReweighting::get1DWeights(double xval, TString key)
{
    double weight_ = 1.;
    TH1F* h_ = wimpWeights1DH_[key];
    if(h_==0) {
	//cout << "cannot find hist: " << key << " weight will be return as 1" << endl;
	return weight_;
    }

    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

    int binx = h_->GetXaxis()->FindBin(xval);
    weight_ = h_->GetBinContent(binx);

    return weight_;
}


//
double WIMPReweighting::get2DWeights(double xval, double yval, TString key)
{
    double weight_ = 1.;
    TH2F* h_ = wimpWeights2DH_[key];
    if(h_==0) {
	//cout << "cannot find hist: " << key << " weight will be return as 1" << endl;
	return weight_;
    }

    int xbins = h_->GetXaxis()->GetNbins();
    if(xval > h_->GetXaxis()->GetBinUpEdge(xbins)    ) xval = h_->GetXaxis()->GetBinUpEdge(xbins);
    if(xval < h_->GetXaxis()->GetBinLowEdge(1)       ) xval = h_->GetXaxis()->GetBinLowEdge(1);

    int ybins = h_->GetYaxis()->GetNbins();
    if(yval > h_->GetYaxis()->GetBinUpEdge(ybins)    ) yval = h_->GetYaxis()->GetBinUpEdge(ybins);
    if(yval < h_->GetYaxis()->GetBinLowEdge(1)       ) yval = h_->GetYaxis()->GetBinLowEdge(1);

    int binx = h_->GetXaxis()->FindBin(xval);
    int biny = h_->GetYaxis()->FindBin(yval);
    weight_ = h_->GetBinContent(binx,biny);

    return weight_;
}


//
WIMPReweighting::~WIMPReweighting()
{
}

