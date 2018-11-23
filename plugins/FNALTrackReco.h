#ifndef __FNALTRACK_RECO__
#define __FNALTRACK_RECO__

#include "interface/PluginBase.h"
#include "interface/TrackTree.h"
#include "interface/Track.h"

using namespace std;


class FNALTrackReco: public PluginBase
{
public:

    //---ctors---
    FNALTrackReco() {};
  
    //---dtor---
    ~FNALTrackReco() {};
   
    //---utils---
    bool Begin(CfgManager& opts, uint64* index);
    bool BeginLoop(int iLoop, CfgManager& opts);
    bool ProcessEvent(H4Tree& event, map<string, PluginBase*>& plugins, CfgManager& opts);
    
private:

    TrackTree* trackTree_;
};

DEFINE_PLUGIN(FNALTrackReco);

#endif
