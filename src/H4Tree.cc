#include "interface/H4Tree.h"

void H4Tree::Init()
{
    //---fill map< <board, group, channel>, pointer to first sample>
    //tree_->GetEntry(0);
    NextEntry();
    unsigned int currentDigiBoard=-1, currentDigiGroup=-1, currentDigiChannel=-1;
    for(unsigned int iSample=0; iSample<nDigiSamples; ++iSample)
    {
        if(digiChannel[iSample] != currentDigiChannel ||
           digiGroup[iSample] != currentDigiGroup ||
           digiBoard[iSample] != currentDigiBoard)
        {
            currentDigiChannel = digiChannel[iSample];
            currentDigiGroup = digiGroup[iSample];
            currentDigiBoard = digiBoard[iSample];            
            digiMap[make_tuple(currentDigiBoard, currentDigiGroup, currentDigiChannel)] = iSample;
        }
    }

    NextEntry(GetTTreePtr()->GetEntriesFast()+1);
}

H4Tree::~H4Tree()
{
    delete[] evtTimeBoard;
    delete[] evtTime;
    delete[] adcBoard;
    delete[] adcChannel;
    delete[] adcData;
    delete[] tdcChannel;
    delete[] tdcData;
    delete[] pattern;
    delete[] patternBoard;
    delete[] patternChannel;
    delete[] digiBoard;
    delete[] digiGroup;
    delete[] digiChannel;
    delete[] digiSampleValue;
    delete[] digiSampleTime;
}
