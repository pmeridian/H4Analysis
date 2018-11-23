#ifndef __H4_TREE__
#define __H4_TREE__

#include <iostream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <string>

#include "DynamicTTree/interface/DynamicTTreeBase.h"

using namespace std;

#define MAX_ADC_CHANNELS 500000
#define MAX_TDC_CHANNELS 200


typedef unsigned long int uint32;
typedef unsigned long long int uint64;
 
//****************************************************************************************
//----------Helper functions--------------------------------------------------------------
//---board+group+channel map key
typedef std::tuple<unsigned int, unsigned int, unsigned int> bgc_key_t;

//---hash functions
struct key_hash : public std::unary_function<bgc_key_t, size_t>
{
    size_t operator()(const bgc_key_t& k) const
        {
            return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
        }
};

struct key_equal : public std::binary_function<bgc_key_t, bgc_key_t, bool>
{
    bool operator()(const bgc_key_t& v0, const bgc_key_t& v1) const
        {
            return (
                std::get<0>(v0) == std::get<0>(v1) &&
                std::get<1>(v0) == std::get<1>(v1) &&
                std::get<2>(v0) == std::get<2>(v1)
                );
        }
};

//---map type
typedef std::unordered_map<const bgc_key_t, int, key_hash, key_equal> bgc_map_t;

//****************************************************************************************
//----------Tree reader class-------------------------------------------------------------
#define DYNAMIC_TREE_NAME H4TreeBase

#define DATA_TABLE                              \
    DATA(unsigned int,  evtTimeStart)           \
    DATA(unsigned int,  runNumber)              \
    DATA(unsigned int,  spillNumber)            \
    DATA(unsigned int,  evtNumber)              \
    DATA(unsigned int,  nEvtTimes)              \
    DATA(unsigned int,  nAdcChannels)           \
    DATA(unsigned int,  nTdcChannels)           \
    DATA(unsigned int,  nPatterns)              \
    DATA(unsigned int,  nDigiSamples)

#define DATA_VECT_TABLE                                 \
    DATA(unsigned int, evtTimeBoard, nEvtTimes)         \
    DATA(uint64,       evtTime, nEvtTimes)              \
    DATA(unsigned int, adcBoard, MAX_ADC_CHANNELS)      \
    DATA(unsigned int, adcChannel, MAX_ADC_CHANNELS)    \
    DATA(unsigned int, adcData, MAX_ADC_CHANNELS)       \
    DATA(unsigned int, tdcChannel, MAX_TDC_CHANNELS)    \
    DATA(unsigned int, tdcData, MAX_TDC_CHANNELS)       \
    DATA(unsigned int, pattern, nPatterns)              \
    DATA(unsigned int, patternBoard, nPatterns)         \
    DATA(unsigned int, patternChannel, nPatterns)       \
    DATA(unsigned int, digiBoard, nDigiSamples)         \
    DATA(unsigned int, digiGroup, nDigiSamples)         \
    DATA(unsigned int, digiChannel, nDigiSamples)       \
    DATA(float,     digiSampleTime, nDigiSamples)       \
    DATA(float,     digiSampleValue, nDigiSamples)     

#include "DynamicTTree/interface/DynamicTTreeInterface.h"

class H4Tree : public H4TreeBase
{
public:
    //---ctors
    H4Tree(TChain* t):
        H4TreeBase(t)
        {
            Init();
        }
    H4Tree(TTree* t):
        H4TreeBase(t)
        {
            Init();
        }    
    //---dtor
    ~H4Tree();
    
    //---members
    void Init();
    uint64 GetEntries(){ return tree_->GetEntriesFast(); };
    
    bgc_map_t digiMap;
};
   
#endif 
