#include "WireChamberReco.h"

//**********Utils*************************************************************************
//----------Begin-------------------------------------------------------------------------
bool WireReco::Begin(CfgManager& opts, int* index)
{
    //---get wire chamber mapping
    chXl_ = opts.GetOpt<int>(instanceName_+".chXleft");
    chXr_ = opts.GetOpt<int>(instanceName_+".chXright");
    chYu_ = opts.GetOpt<int>(instanceName_+".chYup");
    chYd_ = opts.GetOpt<int>(instanceName_+".chYdown");
    
    //---create a position tree
    RegisterSharedData(new TTree("wire", "wire_tree"), "wire_tree", true);
    wireTree_ = PositionTree(index, (TTree*)data_.back().obj, 1);
    wireTree_.Init();
    
    return true;
}

//----------ProcessEvent------------------------------------------------------------------
bool WireReco::ProcessEvent(const H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---search the first (in time) hit for each channel
    vector<float> timeL, timeR, timeU, timeD;
    for(int iCh=0; iCh<h4Tree.nTdcChannels; ++iCh)
    {
        if(h4Tree.tdcChannel[iCh]==chXl_)
            timeL.push_back(h4Tree.tdcData[iCh]);
        if(h4Tree.tdcChannel[iCh]==chXr_)
            timeR.push_back(h4Tree.tdcData[iCh]);
        if(h4Tree.tdcChannel[iCh]==chYu_)
            timeU.push_back(h4Tree.tdcData[iCh]);
        if(h4Tree.tdcChannel[iCh]==chYd_)
            timeD.push_back(h4Tree.tdcData[iCh]);
    }

    //---compute X and Y from channels times
    if(timeR.size()!=0 && timeL.size()!=0)
        wireTree_.X[0] = (*min_element(timeR.begin(), timeR.begin()+timeR.size()) -
                  *min_element(timeL.begin(), timeL.begin()+timeL.size()))*0.005;
    else
        wireTree_.X[0] = -1000;
    if(timeU.size()!=0 && timeD.size()!=0)
        wireTree_.Y[0] = (*min_element(timeU.begin(), timeU.begin()+timeU.size()) -
                  *min_element(timeD.begin(), timeD.begin()+timeD.size()))*0.005;
    else
        wireTree_.Y[0] = -1000;

    //---fill output tree
    wireTree_.Fill();
    
    return true;
}