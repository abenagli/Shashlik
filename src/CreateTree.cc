#include "CreateTree.hh"



CreateTree* CreateTree::fInstance = NULL ;


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::CreateTree (TString name)
{
  if ( fInstance )
  {
    return ;
  }

  this->fInstance = this ;
  this->fname     = name ;
  this->ftree     = new TTree (name,name) ;
  
  this->GetTree ()->Branch ("Event",                  &this->Event,                  "Event/I") ;
  this->GetTree ()->Branch ("totalPhLengthInChamfer", &this->totalPhLengthInChamfer, "totalPhLengthInChamfer[4]/F") ;
  this->GetTree ()->Branch ("numPhLengthInChamfer",   &this->numPhLengthInChamfer,   "numPhLengthInChamfer[4]/I") ;
  
  this->Clear () ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


CreateTree::~CreateTree ()
{}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


int CreateTree::Fill () 
{ 
  for (std::map <int, std::pair<int, float> >::const_iterator iMap = fsingleGammaInfo.begin () ;
       iMap != fsingleGammaInfo.end () ;
       ++iMap)
    {
      ++numPhLengthInChamfer[iMap->second.first] ;
    }
  return this->GetTree ()->Fill () ; 
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


bool CreateTree::Write ()
{
  TString filename = this->GetName () ;
  filename+=".root" ;
  TFile* file = new TFile (filename,"RECREATE") ;
  this->GetTree ()->Write () ;
  file->Write () ;
  file->Close () ;
  return true ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void CreateTree::addPhoton (int trackId, float length, int chamferId)
{
  if (fsingleGammaInfo.find (trackId) == fsingleGammaInfo.end ())
    fsingleGammaInfo[trackId] = std::pair<int, float> (0, length) ;
  else  
    fsingleGammaInfo[trackId].second += length ;
  return ;
}


// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void CreateTree::Clear ()
{
  Event	= 0 ;
  for (int i = 0 ; i < 4 ; ++i) 
    {
      totalPhLengthInChamfer[i] = 0. ;
      numPhLengthInChamfer[i] = 0. ;
    }
  fsingleGammaInfo.clear () ;
}
