#include "FitpixReco.h"

//**********Utils*************************************************************************
//----------Begin*************************************************************************
bool FitpixReco::Begin(CfgManager& opts, uint64* index)
{
    //---create a position tree
    bool storeTree = opts.OptExist(instanceName_+".storeTree") ?
        opts.GetOpt<bool>(instanceName_+".storeTree") : true;

    //---create a position tree
    string treeName = opts.OptExist(instanceName_+".treeName") ?
        opts.GetOpt<string>(instanceName_+".treeName") : "fitpix_tree";

    RegisterSharedData(new TTree(treeName.c_str(), treeName.c_str()), treeName.c_str(), storeTree);

    fitpixTree_ = FitpixTree(index, (TTree*)data_.back().obj);
    fitpixTree_.Init();

    boardId_ = opts.GetOpt<int>(instanceName_+".boardId"); 

    return true;
}

bool FitpixReco::ProcessEvent(const H4Tree& h4Tree, map<string, PluginBase*>& plugins, CfgManager& opts)
{
  hits_.clear();
  clusters_.clear();
  fitpixTree_.Clear();

  std::vector<FPHit*> clustered_hits;

  for(unsigned int i=0; i<h4Tree.nAdcChannels; ++i)
    {
      if(h4Tree.adcBoard[i] != boardId_)
	continue;

      int x = h4Tree.adcChannel[i] % FITPIX_PIXELS_X;
      int y = floor(h4Tree.adcChannel[i]/FITPIX_PIXELS_Y);
      hits_.push_back(FPHit(x,y,h4Tree.adcData[i]));
      clustered_hits.push_back(&hits_.back());
    }

  while( clustered_hits.size()>0 ) {
    std::vector<FPHit*>::iterator i_hit = clustered_hits.begin();
    FPCluster this_cluster;
    this_cluster.add_hit( *(*i_hit) ); // the seed
    clustered_hits.erase( i_hit ); // erase takes care of moving the pointer fwd so no need for ++

    
    std::vector<FPHit*>::iterator j_hit = clustered_hits.begin();

    while( j_hit != clustered_hits.end() ) 
      { // add all adjacent hits to the cluster
	if( this_cluster.isAdjacent(*(*j_hit)) )
	  {
	    this_cluster.add_hit(*(*j_hit));
	    clustered_hits.erase(j_hit); // erase takes care of moving the pointer fwd so no need for ++
	  } 
	else
	  j_hit++;
      } // while
    clusters_.push_back(this_cluster);
  }  // while i

  //  cout << "[FITPIX]:: Hits " << hits_.size() << " Clusters " << clusters_.size() << endl;
  fitpixTree_.n_hits=hits_.size();
  for (int i=0;i<hits_.size();++i)
    {
      fitpixTree_.hitX.push_back(hits_[i].x_);
      fitpixTree_.hitY.push_back(hits_[i].y_);
      fitpixTree_.hitCharge.push_back(hits_[i].c_);
      //      cout << "i_hit\t" << i << "\t" << hits_[i].x_ << "\t" << hits_[i].y_  << "\t" << hits_[i].c_ << endl;
    }
  fitpixTree_.n_clusters=clusters_.size();
  for (int i=0;i<clusters_.size();++i)
    {
      fitpixTree_.clusterX.push_back(clusters_[i].x());
      fitpixTree_.clusterY.push_back(clusters_[i].y());
      fitpixTree_.clusterCharge.push_back(clusters_[i].charge());
      fitpixTree_.clusterSize.push_back(clusters_[i].nhits());
      //      cout << "i_cluster\t" << i << "\t" << clusters_[i].x() << "\t" << clusters_[i].y()  << "\t" << clusters_[i].charge() << endl; 
    }
  //---fill output tree
  fitpixTree_.Fill();
  
  return true;
}
