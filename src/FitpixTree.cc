#include "interface/FitpixTree.h"

FitpixTree::FitpixTree(uint64* idx, TTree* tree)
{
    tree_ = tree ? tree : new TTree();

    n_hits = 0;
    n_clusters = 0;

    index=idx;
}

void FitpixTree::Init()
{
    //---global branches
    tree_->Branch("index", index, "index/l");

    tree_->Branch("n_hits", &n_hits, "n_planes/I");
    tree_->Branch("n_clusters", &n_clusters, "n_clusters/I");
    tree_->Branch("hitX", &hitX);
    tree_->Branch("hitY", &hitY);
    tree_->Branch("hitCharge", &hitCharge);

    tree_->Branch("clusterX", &clusterX);
    tree_->Branch("clusterY", &clusterY);
    tree_->Branch("clusterCharge", &clusterCharge);
    tree_->Branch("clusterSize", &clusterSize);
}
