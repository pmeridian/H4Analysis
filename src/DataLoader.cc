#include "interface/DataLoader.h"

//**********Contructor********************************************************************
DataLoader::DataLoader(CfgManager& opts):
    currentFile_(NULL), inTree_(NULL), firstEventInSpill_(false)
{
    opts_ = opts.GetSubCfg("h4reco");
    ReadInputFiles();
}

//**********Utils*************************************************************************
//----------Read input files from disk----------------------------------------------------
//---Get list of files from run folder and store the list
bool DataLoader::ReadInputFiles()
{
    nFiles_=0;
    int firstSpill=opts_.GetOpt<int>("h4reco.firstSpill");
    string ls_command;
    string file;
    string path=opts_.GetOpt<string>("h4reco.path2data");
    string run=opts_.GetOpt<string>("h4reco.run");    

    //---Get file list searching in specified path (eos or locally)
    if(path.find("/eos/cms") != string::npos)
    {
        // if ( getMachineDomain() != "cern.ch" )
        //     ls_command = string("gfal-ls root://eoscms/"+path+run+" | grep 'root' > /tmp/"+run+".list");
        // else
        ls_command = string("ls "+path+run+" | grep 'root' > /tmp/"+run+".list");
    }
    else if(path.find("srm://") != string::npos)
        ls_command = string("echo "+path+run+"/`gfal-ls "+path+run+
                            "` | sed -e 's:^.*\\/cms\\/:root\\:\\/\\/xrootd-cms.infn.it\\/\\/:g' | grep 'root' > /tmp/"+run+".list");
    else
        ls_command = string("ls "+path+run+" | grep 'root' > /tmp/"+run+".list");
    system(ls_command.c_str());

    ifstream waveList(string("/tmp/"+run+".list").c_str(), ios::in);
    while(waveList >> file && (opts_.GetOpt<int>("h4reco.maxFiles")<0 || nFiles_<opts_.GetOpt<int>("h4reco.maxFiles")) )
    {
        //---skip files before specified spill
        auto currentSpill = std::stoi(file.substr(0, file.size()-4));
        if(firstSpill == -1 || currentSpill == firstSpill)
        {
            if(path.find("/eos/cms") != string::npos)
            {
                std::cout << "+++ Adding file " << (path+run+"/"+file).c_str() << std::endl;
                fileList_.push_back((path+run+"/"+file).c_str());
            }
            else if(path.find("srm://") != string::npos)
            {
                std::cout << "+++ Adding file " << file << std::endl;
                fileList_.push_back((file).c_str());
            }
            else
            {
                std::cout << "+++ Adding file " << (path+run+"/"+file).c_str() << std::endl;
                fileList_.push_back((path+run+"/"+file).c_str());
            }
            ++nFiles_;
        }
    }
    std::cout << "+++ Added " << nFiles_ << std::endl;
    
    //---reverse list to use pop_back later on
    reverse(fileList_.begin(), fileList_.end());

    return true;
}
    
//----------Load next event---------------------------------------------------------------
bool DataLoader::NextEvent()
{
    //---try to load next entry from current tree/file, otherwise try to load new tree/file
    if(!currentFile_ || !inTree_->NextEntry())
    {
        if(LoadNextFile())
        {
            firstEventInSpill_ = true;
            inTree_->NextEntry();

            return true;
        }
        else
            return false;
    }
    else
    {
        firstEventInSpill_ &= false;
        return true;
    }
}

//----------Load next tree/file-----------------------------------------------------------
bool DataLoader::LoadNextFile()
{
    //---clean data from previous spill
    if(inTree_)
    {
        inTree_->GetTTreePtr()->Delete();            
        delete inTree_;
    }
    if(currentFile_ && currentFile_->IsOpen())
        currentFile_->Close();
        
    //---get new file and data tree
    if(fileList_.size() > 0)
    {
        currentFile_ = TFile::Open(fileList_.back().c_str(), "READ");
        if(currentFile_)
        {
	  std::cout << "Opened " << fileList_.back().c_str() << std::endl;
	  inTree_ = new H4Tree((TTree*)currentFile_->Get("H4tree"));
	  inTree_->Init();
	  fileList_.pop_back();
	  return true;
        }
        else
            return false;
    }
    else
        return false;
    
}
