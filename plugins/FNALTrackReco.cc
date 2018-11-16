#include "FNALTrackReco.h"

//----------Begin-------------------------------------------------------------------------
bool FNALTrackReco::Begin(CfgManager& opts, uint64* index)
{
    //---create a position tree
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;
  
    //---create a position tree
    string treeName = opts.OptExist(instanceName_+".treeName") ?
        opts.GetOpt<string>(instanceName_+".treeName") : "track_tree";

    RegisterSharedData(new TTree(treeName.c_str(), treeName.c_str()), treeName.c_str(), storeTree);
  
    trackTree_ = new TrackTree(index, (TTree*)data_.back().obj);
    trackTree_->Init();

    return true;
}

bool FNALTrackReco::BeginLoop(int iLoop, CfgManager& opts)
{
    return true;
}


//----------ProcessEvent------------------------------------------------------------------
bool FNALTrackReco::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---fill output tree
    trackTree_->Clear();

    trackTree_->n_tracks=event.nTracks;

    //get track par
    std::vector<float> trackPar(4,0.);
    trackPar[0]=event.trackX;
    trackPar[1]=event.trackY;
    trackPar[2]=event.trackXSlope;
    trackPar[4]=event.trackYSlope;
    std::vector<float> covariance(10,-999.);

    //assign track par
    TrackPar par;
    par.value.assign(trackPar.begin(),trackPar.begin()+4);
    par.covariance.assign(covariance.begin(),covariance.begin()+10);

    //fill Track Tree
    trackTree_->fitResult.push_back(par);
    trackTree_->fitStatus.push_back(0);
    trackTree_->trackPattern.push_back(0);
    trackTree_->trackHits.push_back(0);
    trackTree_->trackChi2.push_back(event.trackChi2);

    trackTree_->Fill();
    return true;
}
