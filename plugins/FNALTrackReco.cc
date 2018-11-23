#include "FNALTrackReco.h"

static_assert(sizeof(unsigned int) == sizeof(float));

unsigned int castFloatToUInt(float f) {
  union { float f; unsigned int i; } u;
  u.f = f;
  return u.i;
}

float castUIntToFloat(unsigned int i) {
  union { float f; unsigned int i; } u;
  u.i = i;
  return u.f;
}

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

    std::vector<float> values(4,0.);
    std::vector<float> covariances(10,999.);
    for(unsigned int iADC=0; iADC<event.nAdcChannels; ++iADC)
    {
      if ( event.adcBoard[iADC] !=1 )
	continue;
      if ( event.adcChannel[iADC] == 1)
	trackTree_->n_tracks=event.adcData[iADC];
      if ( event.adcChannel[iADC] == 2)
	trackTree_->trackChi2.push_back(castUIntToFloat(event.adcData[iADC]));
      if ( event.adcChannel[iADC] > 2 &&  event.adcChannel[iADC] < 7)
	values[ event.adcChannel[iADC] - 3 ] = castUIntToFloat(event.adcData[iADC]);
    }
        
    
    TrackPar par;
    par.value.assign(values.begin(),values.begin()+4);
    par.covariance.assign(covariances.begin(),covariances.begin()+10);
    trackTree_->fitResult.push_back(par);
    trackTree_->fitStatus.push_back(0);
    trackTree_->trackPattern.push_back(0);
    trackTree_->trackHits.push_back(0);

    trackTree_->Fill();
    
    return true;
}
