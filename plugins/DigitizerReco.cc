#include "DigitizerReco.h"

//----------Utils-------------------------------------------------------------------------
bool DigitizerReco::Begin(CfgManager& opts, uint64* index)
{
    //---inputs---
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");

    for(auto& channel : channelsNames_)
    {
        nSamples_[channel] = opts.GetOpt<int>(channel+".nSamples");
        auto tUnit = opts.GetOpt<float>(channel+".tUnit");
        if(opts.OptExist(channel+".type"))
        {
            if(opts.GetOpt<string>(channel+".type") == "NINO")
                WFs[channel] = new WFClassNINO(opts.GetOpt<int>(channel+".polarity"), tUnit);
            else if(opts.GetOpt<string>(channel+".type") == "Clock")
                WFs[channel] = new WFClassClock(tUnit);
        }
        else
            WFs[channel] = new WFClass(opts.GetOpt<int>(channel+".polarity"), tUnit);
        RegisterSharedData(WFs[channel], channel, false);
    }
    
    return true;
}

bool DigitizerReco::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{        
    //---read the digitizer
    //---set time reference from digitized trigger
    int trigRef=0;
    // FIXME
    // for(int iSample=nSamples_[channel]*8; iSample<nSamples_[channel]*9; ++iSample)
    // {
    //     if(event.digiSampleValue[iSample] < 1000)
    //     {
    //         trigRef = iSample-nSamples_[channel]*8;
    //         break;
    //     }
    // }
    
    //---user channels
    bool evtStatus = true;
    for(auto& channel : channelsNames_)
    {
        //---reset and read new WFs
        WFs[channel]->Reset();
        unsigned int digiBd = opts.GetOpt<unsigned int>(channel+".digiBoard");
        unsigned int digiGr = opts.GetOpt<unsigned int>(channel+".digiGroup");
        unsigned int digiCh = opts.GetOpt<unsigned int>(channel+".digiChannel");
        int offset = event.digiMap.at(make_tuple(digiBd, digiGr, digiCh));
        for(int iSample=offset; iSample<offset+nSamples_[channel]; ++iSample)
        {
            //---H4DAQ bug: sometimes ADC value is out of bound.
            //---skip everything if one channel is bad
            if(event.digiSampleValue[iSample] > 1e6)
            {
                evtStatus = false;
                WFs[channel]->AddSample(4095,999.);
            }
            else
	      WFs[channel]->AddSample(event.digiSampleValue[iSample],event.digiSampleTime[iSample]);
        }
        if(opts.OptExist(channel+".useTrigRef") && opts.GetOpt<bool>(channel+".useTrigRef"))
            WFs[channel]->SetTrigRef(trigRef);
    }

    if(!evtStatus)
        cout << ">>>DigiReco WARNING: bad amplitude detected" << endl;

    return evtStatus;
}
