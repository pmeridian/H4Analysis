#include "WFAnalyzer.h"

//----------Utils-------------------------------------------------------------------------
bool WFAnalyzer::Begin(CfgManager& opts, uint64* index)
{
    //---inputs---
    if(!opts.OptExist(instanceName_+".srcInstanceName"))
    {
        cout << ">>> WFAnalyzer ERROR: no source plugin specified" << endl;
        return false;
    }
    srcInstance_ = opts.GetOpt<string>(instanceName_+".srcInstanceName");
    channelsNames_ = opts.GetOpt<vector<string> >(instanceName_+".channelsNames");
    timeRecoTypes_ = opts.GetOpt<vector<string> >(instanceName_+".timeRecoTypes");

    //---channels setup
    string templateTag="";
    if(opts.OptExist(instanceName_+".templateTags"))
        for(auto& tag : opts.GetOpt<vector<string> >(instanceName_+".templateTags"))
            for(auto& run : opts.GetOpt<vector<string> >(tag+".runList"))
                if(run == opts.GetOpt<string>("h4reco.run"))
                    templateTag = tag;

    int nSamples = 0;
    for(auto& channel : channelsNames_)
    {        
        nSamples += opts.GetOpt<int>(channel+".nSamples");
        if(opts.OptExist(channel+".templateFit.file"))
        {            
            TFile* templateFile = TFile::Open(opts.GetOpt<string>(channel+".templateFit.file", 0).c_str(), ".READ");
            TH1* wfTemplate=(TH1*)templateFile->Get((opts.GetOpt<string>(channel+".templateFit.file", 1)+templateTag).c_str());
            templates_[channel] = (TH1F*) wfTemplate->Clone();
            templates_[channel] -> SetDirectory(0);
            templateFile->Close();
        }
        //---keep track of all the possible time reco method requested
        for(auto type_name : timeRecoTypes_)
        {
            if(opts.OptExist(channel+"."+type_name))
                timeOpts_[channel+"."+type_name] = opts.GetOpt<vector<float> >(channel+"."+type_name);
        }
    }
    
    //---outputs---
    string digiTreeName = opts.OptExist(instanceName_+".digiTreeName") ?
        opts.GetOpt<string>(instanceName_+".digiTreeName") : "digi";
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;
    RegisterSharedData(new TTree(digiTreeName.c_str(), "digi_tree"), "digi_tree", storeTree);
    digiTree_ = DigiTree(index, (TTree*)data_.back().obj);
    digiTree_.Init(channelsNames_, timeRecoTypes_);
    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
    {
        string wfTreeName = opts.OptExist(instanceName_+".wfTreeName") ?
            opts.GetOpt<string>(instanceName_+".wfTreeName") : "wf";
        RegisterSharedData(new TTree(wfTreeName.c_str(), "wf_tree"), "wf_tree", true);
        outWFTree_ = WFTree(nSamples, index, (TTree*)data_.back().obj);
        outWFTree_.Init();
    }

    return true;
}

bool WFAnalyzer::ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts)
{
    //---setup output event 
    int outCh=0;
    bool fillWFtree=false;
    if(opts.GetOpt<int>(instanceName_+".fillWFtree"))
        fillWFtree = *digiTree_.index % opts.GetOpt<int>(instanceName_+".WFtreePrescale") == 0;

    //---load WFs from source instance shared data
    for(auto& channel : channelsNames_)
    {
        auto shared_data = plugins[srcInstance_]->GetSharedData(srcInstance_+"_"+channel, "", false);
        if(shared_data.size() != 0)
            WFs_[channel] = (WFClass*)shared_data.at(0).obj;
        else
            cout << "[WFAnalizer::" << instanceName_ << "]: channels samples not found check DigiReco step" << endl; 
    }

    //---compute reco variables
    for(auto& channel : channelsNames_)
    {
        //---skip dead channels
        if(WFs_.find(channel) == WFs_.end())
        {
            ++outCh;
            continue;
        }

        //---subtract a specified channel if requested
        if(opts.OptExist(channel+".subtractChannel") && WFs_.find(opts.GetOpt<string>(channel+".subtractChannel")) != WFs_.end())
            *WFs_[channel] -= *WFs_[opts.GetOpt<string>(channel+".subtractChannel")];        
        WFs_[channel]->SetBaselineWindow(opts.GetOpt<int>(channel+".baselineWin", 0), 
                                         opts.GetOpt<int>(channel+".baselineWin", 1));
        WFs_[channel]->SetBaselineIntegralWindow(opts.GetOpt<int>(channel+".baselineInt", 0),
                                                 opts.GetOpt<int>(channel+".baselineInt", 1));
        WFs_[channel]->SetSignalWindow(opts.GetOpt<int>(channel+".signalWin", 0), 
                                       opts.GetOpt<int>(channel+".signalWin", 1));
        WFs_[channel]->SetSignalIntegralWindow(opts.GetOpt<int>(channel+".signalInt", 0),
                                               opts.GetOpt<int>(channel+".signalInt", 1));
        WFBaseline baselineInfo = WFs_[channel]->SubtractBaseline();
        //FIXME MAREMMA MAIALA
        WFFitResults interpolAmpMax;
        if(opts.OptExist(channel+".signalWin", 4))
        {
            string max_function = opts.OptExist(channel+".signalWin", 4) ? opts.GetOpt<string>(channel+".signalWin", 4) : "pol2";
            int nParams = opts.OptExist(channel+".signalWin", 5) ? opts.GetOpt<int>(channel+".signalWin", 5) : 0;
            std::vector<float> max_params;
            for(int ipar = 0; ipar < nParams; ++ipar)
                max_params.push_back( opts.GetOpt<float>(channel+".signalWin", 6+ipar) );
            interpolAmpMax = WFs_[channel]->GetInterpolatedAmpMax(-1,-1,
                                                                  opts.GetOpt<int>(channel+".signalWin", 2),
                                                                  opts.GetOpt<int>(channel+".signalWin", 3), max_function, max_params);
        }
        else
        {
            string max_function = opts.OptExist(channel+".signalWin", 3) ? opts.GetOpt<string>(channel+".signalWin", 3) : "pol2";
            interpolAmpMax = WFs_[channel]->GetInterpolatedAmpMax(-1,-1,
                                                                  opts.GetOpt<int>(channel+".signalWin", 2)/2,
                                                                  opts.GetOpt<int>(channel+".signalWin", 2)/2,
                                                                  max_function);
        }
        digiTree_.pedestal[outCh] = baselineInfo.baseline;
        digiTree_.b_charge[outCh] = WFs_[channel]->GetIntegral(opts.GetOpt<int>(channel+".baselineInt", 0), 
                                                               opts.GetOpt<int>(channel+".baselineInt", 1));        
        digiTree_.b_slope[outCh] = baselineInfo.slope;
        digiTree_.b_rms[outCh] = baselineInfo.rms;
        digiTree_.maximum[outCh] = WFs_[channel]->GetAmpMax();
        digiTree_.time_maximum[outCh] = WFs_[channel]->GetTimeCF(1).time;
        digiTree_.amp_max[outCh] = interpolAmpMax.ampl;
        digiTree_.time_max[outCh] = interpolAmpMax.time;
        digiTree_.chi2_max[outCh] = interpolAmpMax.chi2;
        digiTree_.charge_tot[outCh] = WFs_[channel]->GetModIntegral(opts.GetOpt<int>(channel+".baselineInt", 1), 
                                                                    WFs_[channel]->GetNSample());
        digiTree_.charge_sig[outCh] = WFs_[channel]->GetSignalIntegral(opts.GetOpt<int>(channel+".signalInt", 0), 
                                                                       opts.GetOpt<int>(channel+".signalInt", 1));
        //---compute time with all the requested time reconstruction method
        for(unsigned int iT=0; iT<timeRecoTypes_.size(); ++iT)
        {
            //---compute time with selected method or store default value (-99)
            if(timeOpts_.find(channel+"."+timeRecoTypes_[iT]) != timeOpts_.end())
            {
                WFFitResults timeInfo = WFs_[channel]->GetTime(timeRecoTypes_[iT], timeOpts_[channel+"."+timeRecoTypes_[iT]]);
                digiTree_.time[outCh+iT*channelsNames_.size()] = timeInfo.time;
                digiTree_.time_chi2[outCh+iT*channelsNames_.size()] = timeInfo.chi2;
                digiTree_.time_slope[outCh+iT*channelsNames_.size()] = timeInfo.slope;
            }
            else
            {
                digiTree_.time[outCh+iT*channelsNames_.size()] = -99;
                digiTree_.time_chi2[outCh+iT*channelsNames_.size()] = -99;
                digiTree_.time_slope[outCh+iT*channelsNames_.size()] = -99;
            }
        }

        //---template fit (only specified channels)
        WFFitResults fitResults{-1, -1000, -1};
        if(opts.OptExist(channel+".templateFit.file"))
        {
            WFs_[channel]->SetTemplate(templates_[channel]);
            fitResults = WFs_[channel]->TemplateFit(opts.GetOpt<float>(channel+".templateFit.fitWin", 0),
                                                    opts.GetOpt<int>(channel+".templateFit.fitWin", 1),
                                                    opts.GetOpt<int>(channel+".templateFit.fitWin", 2));
            digiTree_.fit_ampl[outCh] = fitResults.ampl;
            digiTree_.fit_time[outCh] = fitResults.time;
            digiTree_.fit_chi2[outCh] = fitResults.chi2;
        }            
        //---calibration constant for each channel if needed
        if(opts.OptExist(channel+".calibration.calibrationConst"))
            digiTree_.calibration[outCh]=opts.GetOpt<float>(channel+".calibration.calibrationConst");
        else
            digiTree_.calibration[outCh]=1;
	
        //---WFs---
        if(fillWFtree)
        {
            auto analizedWF = WFs_[channel]->GetSamples();
            auto times = WFs_[channel]->GetTimes();
	    //            float tUnit = WFs_[channel]->GetTUnit();
            for(unsigned int jSample=0; jSample<analizedWF->size(); ++jSample)
            {
                outWFTree_.WF_ch.push_back(outCh);
                outWFTree_.WF_time.push_back(times->at(jSample));
                outWFTree_.WF_val.push_back(analizedWF->at(jSample));
            }
        }
        //---increase output tree channel counter
        ++outCh;
    }

    //---fill the output trees 
    //---reco var
    digiTree_.Fill();
    //---WFs
    if(fillWFtree)
        outWFTree_.Fill();

    return true;
}
